#!/usr/bin/env julia

# Diagnostic-only whole-ROI mold assembly. This file does not enter fitting,
# N selection, calibration, or production unit assignment.

if !isdefined(Main, :JointProxyRegistry)
    include(joinpath(@__DIR__, "lib", "joint_proxy", "proxy_registry.jl"))
end

const WHOLE_ROI_REGISTRY = Main.JointProxyRegistry

function _whole_roi_template(entry, typ::Int, parity::Int, mirror::Int)
    matches = [template for template in entry.templates
        if template.type == typ && template.parity == parity && template.mirror == mirror]
    length(matches) == 1 || throw(ArgumentError(
        "expected one template for type=$typ parity=$parity mirror=$mirror"))
    return only(matches).pixels
end

function _whole_roi_parity(lobe::Int, n::Int, direction::Int, phase::Int)
    direction in (0, 1) || throw(ArgumentError("direction must be 0 or 1"))
    phase in (0, 1) || throw(ArgumentError("phase must be 0 or 1"))
    index = direction == 0 ? lobe : n - lobe + 1
    return mod(index - 1 + phase, 2)
end

function _template_matrix(pixels::AbstractVector{<:Real}, grid_n::Int)
    length(pixels) == grid_n^2 || throw(ArgumentError("template/grid size mismatch"))
    # Pixels are ordered with u outer and t inner. Julia matrices are column-major.
    return permutedims(reshape(Float64.(pixels), grid_n, grid_n))
end

function _sample_template(matrix::Matrix{Float64}, t::Float64, u::Float64,
                          half_nm::Float64, step_nm::Float64)
    tolerance = 32eps(Float64) * max(1.0, half_nm)
    (t < -half_nm - tolerance || t > half_nm + tolerance ||
     u < -half_nm - tolerance || u > half_nm + tolerance) && return 0.0
    n = size(matrix, 1)
    size(matrix, 2) == n || throw(ArgumentError("template matrix must be square"))
    qt = clamp((t + half_nm) / step_nm + 1.0, 1.0, Float64(n))
    qu = clamp((u + half_nm) / step_nm + 1.0, 1.0, Float64(n))
    it = min(floor(Int, qt), n - 1)
    iu = min(floor(Int, qu), n - 1)
    ft = qt - it
    fu = qu - iu
    return (1 - ft) * (1 - fu) * matrix[iu, it] +
           ft * (1 - fu) * matrix[iu, it + 1] +
           (1 - ft) * fu * matrix[iu + 1, it] +
           ft * fu * matrix[iu + 1, it + 1]
end

"""Superpose every selected typed mold on one common image grid.

Overlap is represented by ordinary addition in image space. Values outside each
tracked mold's finite support are zero; this makes truncation explicit rather
than silently extrapolating edge pixels.
"""
function assemble_whole_roi(entry, sequence::AbstractVector{<:Integer},
                            centers::AbstractVector, xs::Vector{Float64},
                            ys::Vector{Float64}; half_nm::Float64,
                            step_nm::Float64, theta::Float64=0.0,
                            direction::Int=0, phase::Int=0, mirror::Int=0)
    half_nm > 0 || throw(ArgumentError("half_nm must be positive"))
    step_nm > 0 || throw(ArgumentError("step_nm must be positive"))
    length(sequence) == length(centers) || throw(ArgumentError("sequence/center length mismatch"))
    all(typ -> typ in (0, 1), sequence) || throw(ArgumentError("sequence values must be 0 or 1"))
    mirror in (0, 1) || throw(ArgumentError("mirror must be 0 or 1"))
    grid_n = length(collect(-half_nm:step_nm:half_nm))
    grid_n >= 2 || throw(ArgumentError("mold grid must contain at least two points"))
    cos_theta, sin_theta = cos(theta), sin(theta)
    assembled = zeros(Float64, length(ys), length(xs))
    n = length(sequence)
    for lobe in eachindex(sequence)
        parity = _whole_roi_parity(lobe, n, direction, phase)
        pixels = _whole_roi_template(entry, Int(sequence[lobe]), parity, mirror)
        template = _template_matrix(pixels, grid_n)
        center_x, center_y = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx = xs[ix] - center_x
            dy = ys[iy] - center_y
            t = dx * cos_theta + dy * sin_theta
            u = -dx * sin_theta + dy * cos_theta
            assembled[iy, ix] += _sample_template(template, t, u, half_nm, step_nm)
        end
    end
    return assembled
end

function _whole_roi_candidate(observed, predictors::Vector{Matrix{Float64}}, active)
    good = isfinite.(observed)
    for predictor in predictors
        good .&= isfinite.(predictor)
    end
    count(good) >= length(active) + 1 || throw(ArgumentError("insufficient finite ROI pixels"))
    y = Float64.(observed[good])
    design = ones(Float64, length(y), length(active) + 1)
    for (column, predictor_index) in enumerate(active)
        design[:, column + 1] .= predictors[predictor_index][good]
    end
    coefficients = design \ y
    any(coefficients[2:end] .< -sqrt(eps(Float64))) && return nothing
    coefficients[2:end] .= max.(coefficients[2:end], 0.0)
    residual = y - design * coefficients
    full = zeros(Float64, 2)
    for (column, predictor_index) in enumerate(active)
        full[predictor_index] = coefficients[column + 1]
    end
    return (; sse=sum(abs2, residual), background=coefficients[1], coefficients=full,
        n_pixels=length(y))
end

"""Fit background + nonnegative backbone + nonnegative assembled-mold terms."""
function fit_whole_roi_nuisance(observed::Matrix{Float64}, backbone::Matrix{Float64},
                                mold::Matrix{Float64})
    size(observed) == size(backbone) == size(mold) || throw(ArgumentError("ROI shape mismatch"))
    predictors = [backbone, mold]
    candidates = filter(!isnothing, [_whole_roi_candidate(observed, predictors, active)
        for active in (Int[], [1], [2], [1, 2])])
    isempty(candidates) && throw(ArgumentError("whole-ROI nuisance fit failed"))
    best = candidates[argmin(getproperty.(candidates, :sse))]
    return (; sse=best.sse, background=best.background,
        backbone_coefficient=best.coefficients[1], mold_coefficient=best.coefficients[2],
        n_pixels=best.n_pixels)
end

"""Exhaustive short-chain diagnostic search using one global ROI score."""
function search_whole_roi_sequences(observed::Matrix{Float64}, backbone::Matrix{Float64},
                                    entry, centers::AbstractVector,
                                    xs::Vector{Float64}, ys::Vector{Float64};
                                    half_nm::Float64, step_nm::Float64,
                                    theta::Float64=0.0, directions=(0, 1),
                                    phases=(0, 1), mirrors=(0, 1), max_n::Int=12)
    1 <= max_n <= 12 || throw(ArgumentError("max_n must be between 1 and 12"))
    n = length(centers)
    n <= max_n || throw(ArgumentError("N=$n exceeds exhaustive diagnostic max_n=$max_n"))
    n >= 1 || throw(ArgumentError("whole-ROI search requires at least one lobe"))
    best = nothing
    for bits in 0:(2^n - 1)
        sequence = [Int((bits >> (lobe - 1)) & 1) for lobe in 1:n]
        for direction in directions, phase in phases, mirror in mirrors
            mold = assemble_whole_roi(entry, sequence, centers, xs, ys;
                half_nm=half_nm, step_nm=step_nm, theta=theta,
                direction=direction, phase=phase, mirror=mirror)
            fit = fit_whole_roi_nuisance(observed, backbone, mold)
            candidate = (; sequence=copy(sequence), direction, phase, mirror, fit, mold)
            if best === nothing || candidate.fit.sse < best.fit.sse
                best = candidate
            end
        end
    end
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Diagnostic helper loaded. Run test/joint_proxy/test_whole_roi_mold_model.jl for the synthetic gate.")
end
