#!/usr/bin/env julia

# Diagnostic-only integration seam: synthesize an observation, discard truth,
# register only the chemistry-independent common mold, freeze that result, then
# score type contrast without reopening geometry.

if !isdefined(Main, :FrozenRegistration)
    include(joinpath(@__DIR__, "whole_roi_common_registration.jl"))
end
if !isdefined(Main, :FrozenContrastResult)
    include(joinpath(@__DIR__, "whole_roi_frozen_contrast.jl"))
end

using Random

"""Truth-free inputs accepted by the common-registration/scoring pipeline."""
struct WholeRoiObservedCase
    observed::Matrix{Float64}
    backbone::Matrix{Float64}
    centers::Vector{Tuple{Float64,Float64}}
    xs::Vector{Float64}
    ys::Vector{Float64}
    theta_base::Float64
end

"""Frozen registration, scoring evidence, and explicit search-wall status."""
struct WholeRoiDiagnosticResult
    registration::FrozenRegistration
    scoring::FrozenContrastResult
    boundary_parameters::Vector{Symbol}
end

function _e2e_validate_geometry(backbone, centers, xs, ys)
    size(backbone) == (length(ys), length(xs)) || throw(ArgumentError(
        "backbone shape $(size(backbone)) does not match grid ($(length(ys)), $(length(xs)))"))
    length(xs) >= 2 || throw(ArgumentError("xs requires at least two grid points"))
    length(ys) >= 2 || throw(ArgumentError("ys requires at least two grid points"))
    isempty(centers) && throw(ArgumentError("at least one lobe center is required"))
    return nothing
end

@inline function _e2e_sample(template::Matrix{Float64}, t::Float64, u::Float64,
                             half_nm::Float64, step_nm::Float64)
    ft = (t + half_nm) / step_nm + 1.0
    fu = (u + half_nm) / step_nm + 1.0
    n_u, n_t = size(template)
    (1.0 <= ft <= n_t && 1.0 <= fu <= n_u) || return 0.0
    t0 = clamp(floor(Int, ft), 1, n_t)
    u0 = clamp(floor(Int, fu), 1, n_u)
    t1 = min(t0 + 1, n_t)
    u1 = min(u0 + 1, n_u)
    wt = ft - t0
    wu = fu - u0
    return (1 - wu) * ((1 - wt) * template[u0, t0] + wt * template[u0, t1]) +
           wu * ((1 - wt) * template[u1, t0] + wt * template[u1, t1])
end

function _e2e_forward_components(set, sequence, centers, xs, ys;
        half_nm, step_nm, theta_total, direction, phase, mirror,
        shift_t, shift_u)
    common = zeros(Float64, length(ys), length(xs))
    contrast = zeros(Float64, length(ys), length(xs))
    n = length(centers)
    cos_theta, sin_theta = cos(theta_total), sin(theta_total)
    for lobe in eachindex(centers)
        index = direction == 0 ? lobe : n - lobe + 1
        parity = mod(index - 1 + phase, 2)
        mold = common_mold(set, parity, mirror)
        sign = sequence[lobe] == 1 ? 1.0 : -1.0
        center_x, center_y = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx = xs[ix] - center_x
            dy = ys[iy] - center_y
            t = dx * cos_theta + dy * sin_theta + shift_t
            u = -dx * sin_theta + dy * cos_theta + shift_u
            common[iy, ix] += _e2e_sample(
                mold.common, t, u, half_nm, step_nm)
            contrast[iy, ix] += sign * _e2e_sample(
                mold.contrast, t, u, half_nm, step_nm)
        end
    end
    return common, contrast
end

function _e2e_blur(image::Matrix{Float64}, sigma_px::Float64)
    sigma_px <= 0 && return copy(image)
    radius = max(1, ceil(Int, 3sigma_px))
    kernel = [exp(-0.5 * (offset / sigma_px)^2) for offset in -radius:radius]
    kernel ./= sum(kernel)
    ny, nx = size(image)
    horizontal = similar(image)
    @inbounds for y in 1:ny, x in 1:nx
        horizontal[y, x] = sum(kernel[k] * image[y, clamp(x + k - radius - 1, 1, nx)]
                               for k in eachindex(kernel))
    end
    output = similar(image)
    @inbounds for y in 1:ny, x in 1:nx
        output[y, x] = sum(kernel[k] * horizontal[clamp(y + k - radius - 1, 1, ny), x]
                           for k in eachindex(kernel))
    end
    return output
end

"""Generate data with a known sequence, but return no truth-bearing field.

The common and signed-contrast images are independently assembled under the
same injected rigid nuisance and blurred separately before combination. The
sequence is consumed only in this forward generator and is not retained in the
returned `WholeRoiObservedCase`.
"""
function synthesize_whole_roi_observation(
        set::CommonContrastSet, sequence::AbstractVector{<:Integer},
        backbone::AbstractMatrix{<:Real}, centers::AbstractVector,
        xs::Vector{Float64}, ys::Vector{Float64};
        half_nm::Real, step_nm::Real, theta_base::Real=0.0,
        direction::Int=0, phase::Int=0, mirror::Int=0,
        shift_t_nm::Real=0.0, shift_u_nm::Real=0.0,
        rotation_deg::Real=0.0, blur_sigma_nm::Real=0.0,
        background::Real=0.0, backbone_amp::Real=1.0,
        common_amp::Real=1.0, contrast_amp::Real=1.0,
        noise_sigma::Real=0.0, seed::Integer=0)
    _e2e_validate_geometry(backbone, centers, xs, ys)
    length(sequence) == length(centers) ||
        throw(ArgumentError("sequence/center length mismatch"))
    all(bit -> bit in (0, 1), sequence) ||
        throw(ArgumentError("sequence values must be 0 or 1"))
    blur_sigma_nm >= 0 || throw(ArgumentError("blur_sigma_nm must be nonnegative"))
    noise_sigma >= 0 || throw(ArgumentError("noise_sigma must be nonnegative"))

    half_f = Float64(half_nm)
    step_f = Float64(step_nm)
    theta_total = Float64(theta_base) + deg2rad(Float64(rotation_deg))
    common, contrast = _e2e_forward_components(set, sequence, centers, xs, ys;
        half_nm=half_f, step_nm=step_f, theta_total=theta_total,
        direction=direction, phase=phase, mirror=mirror,
        shift_t=Float64(shift_t_nm), shift_u=Float64(shift_u_nm))

    pixel_step = xs[2] - xs[1]
    pixel_step > 0 || throw(ArgumentError("pixel spacing must be positive"))
    blur_px = Float64(blur_sigma_nm) / pixel_step
    if blur_px > 0
        common = _e2e_blur(common, blur_px)
        contrast = _e2e_blur(contrast, blur_px)
    end

    backbone_copy = Matrix{Float64}(backbone)
    observed = Float64(background) .+
        Float64(backbone_amp) .* backbone_copy .+
        Float64(common_amp) .* common .+
        Float64(contrast_amp) .* contrast
    if noise_sigma > 0
        observed .+= Float64(noise_sigma) .* randn(MersenneTwister(seed), size(observed))
    end
    center_copy = Tuple{Float64,Float64}[(Float64(c[1]), Float64(c[2])) for c in centers]
    return WholeRoiObservedCase(observed, backbone_copy, center_copy,
        copy(xs), copy(ys), Float64(theta_base))
end

function _e2e_boundary_parameters(reg::FrozenRegistration, config)
    hits = Symbol[]
    for (name, value, bounds) in (
            (:shift_t_nm, reg.shift_t_nm, config.shift_t_nm),
            (:shift_u_nm, reg.shift_u_nm, config.shift_u_nm),
            (:rotation_deg, reg.rotation_deg, config.rotation_deg))
        values = Main.WholeRoiDiagnosticConfig.search_values(bounds.coarse)
        if length(values) > 1 &&
           (isapprox(value, first(values); atol=eps(Float64), rtol=0) ||
            isapprox(value, last(values); atol=eps(Float64), rtol=0))
            push!(hits, name)
        end
    end
    blur_values = Main.WholeRoiDiagnosticConfig.search_values(config.blur_sigma_nm.coarse)
    if length(blur_values) > 1 &&
       isapprox(reg.blur_sigma_nm, last(blur_values); atol=eps(Float64), rtol=0)
        push!(hits, :blur_sigma_nm)
    end
    return hits
end

"""Run chemistry-blind registration followed by frozen contrast scoring."""
function run_whole_roi_common_contrast(case::WholeRoiObservedCase,
                                       set::CommonContrastSet;
                                       config, half_nm::Real, step_nm::Real,
                                       max_n::Int=16)
    reg = register_common(set, case.observed, case.backbone, case.centers,
        case.xs, case.ys; config=config, half_nm=Float64(half_nm),
        step_nm=Float64(step_nm), theta_base=case.theta_base)
    frozen = frozen_scorer_input(reg; theta_base=case.theta_base)
    scoring = score_frozen_contrast(frozen, set, case.observed, case.backbone,
        case.centers, case.xs, case.ys; half_nm=half_nm, step_nm=step_nm,
        config=config, max_n=max_n)
    return WholeRoiDiagnosticResult(reg, scoring,
        _e2e_boundary_parameters(reg, config))
end

"""Canonical registration serialization used by swap/determinism controls."""
function freeze_registration_tsv(reg::FrozenRegistration)
    fields = (reg.direction, reg.phase, reg.mirror, reg.shift_t_nm,
        reg.shift_u_nm, reg.rotation_deg, reg.blur_sigma_nm,
        reg.common_image_hash, reg.grid_n, reg.fit.sse, reg.fit.background,
        join(reg.fit.coefficients, ','), join(reg.fit.active, ','), reg.tie_info)
    return join(fields, '\t') * "\n"
end
