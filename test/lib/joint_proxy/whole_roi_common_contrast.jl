#!/usr/bin/env julia

# Diagnostic-only common/contrast decomposition and generalized exhaustive
# active-set NNLS for the whole-ROI mold model.
#
# This file extends test/diagnose_joint_proxy_whole_roi.jl WITHOUT duplicating
# its parity formula, template-matrix reshape, or bilinear sampler: it includes
# that helper and reuses _whole_roi_parity, _template_matrix,
# _whole_roi_template, _sample_template, assemble_whole_roi, and
# fit_whole_roi_nuisance. There is exactly one parity implementation and one
# pixel-ordering convention (u outer, t inner) shared with the reference.
#
# Scientific invariants (AGENTS.md): no benchmark label, expected N, NKNNKN,
# class count, or composition is read here. The common mold is always the exact
# average of the two type templates; it is never substituted with the type-0
# template ("never all-type-0 common"). Contrast is the exact signed
# half-difference, so M0 = common - contrast and M1 = common + contrast hold
# to floating-point precision for every (parity, mirror).

if !isdefined(Main, :JointProxyRegistry)
    include(joinpath(@__DIR__, "proxy_registry.jl"))
end

if !isdefined(Main, :_whole_roi_parity)
    include(joinpath(@__DIR__, "..", "..", "diagnose_joint_proxy_whole_roi.jl"))
end

using LinearAlgebra

const _CC_REG = Main.JointProxyRegistry

# --- typed results -----------------------------------------------------------

"""Common and signed half-difference (contrast) template matrices for one
(parity, mirror) pair, in the [iu, it] convention used by _template_matrix.

Satisfies, exactly: M0 = common - contrast, M1 = common + contrast, where M0/M1
are the type 0/type 1 template matrices for this (parity, mirror)."""
struct CommonContrastMold
    parity::Int
    mirror::Int
    common::Matrix{Float64}    # (M0 + M1) / 2
    contrast::Matrix{Float64}  # (M1 - M0) / 2
    grid_n::Int
end

"""All four (parity, mirror) common/contrast pairs for one proxy entry."""
struct CommonContrastSet
    molds::Dict{Tuple{Int,Int}, CommonContrastMold}
    grid_n::Int
end

"""Free-intercept least-squares fit with nonnegative predictor coefficients,
selected by exhaustive active-set enumeration. `coefficients` has length
`n_predictors` with all entries >= 0; `active` lists the chosen active set."""
struct ActiveSetFit
    sse::Float64
    background::Float64
    coefficients::Vector{Float64}
    active::Vector{Int}
    n_observations::Int
end

# --- decomposition -----------------------------------------------------------

"""Derive common and contrast templates for every (parity, mirror).

Validates that exactly one template exists for each of the eight
(type, parity, mirror) states — _whole_roi_template raises ArgumentError on a
missing or duplicate state, so a naive all-type-0 fallback cannot occur. The
common template is the exact average of the two type templates."""
function derive_common_contrast(entry)
    # Probe one state to derive grid_n; _whole_roi_template errors if missing.
    m0_00 = Main._whole_roi_template(entry, 0, 0, 0)
    npix = length(m0_00)
    grid_n = isqrt(npix)
    grid_n^2 == npix || throw(ArgumentError(
        "template pixel length $npix is not a perfect square"))
    grid_n >= 2 || throw(ArgumentError("mold grid must contain at least two points"))
    molds = Dict{Tuple{Int,Int}, CommonContrastMold}()
    for parity in (0, 1), mirror in (0, 1)
        m0_pixels = Main._whole_roi_template(entry, 0, parity, mirror)
        m1_pixels = Main._whole_roi_template(entry, 1, parity, mirror)
        length(m0_pixels) == npix == length(m1_pixels) || throw(ArgumentError(
            "inconsistent template pixel length for (type=*, parity=$parity, mirror=$mirror)"))
        M0 = Main._template_matrix(m0_pixels, grid_n)
        M1 = Main._template_matrix(m1_pixels, grid_n)
        common = (M0 .+ M1) ./ 2.0
        contrast = (M1 .- M0) ./ 2.0
        molds[(parity, mirror)] = CommonContrastMold(parity, mirror, common, contrast, grid_n)
    end
    return CommonContrastSet(molds, grid_n)
end

"""Fetch the common/contrast pair for (parity, mirror).

`parity` must be 0 or 1 (ArgumentError otherwise); `mirror` is looked up
directly and raises KeyError if the pair was never built."""
function common_mold(set::CommonContrastSet, parity::Int, mirror::Int)
    parity in (0, 1) || throw(ArgumentError("parity must be 0 or 1, got $parity"))
    return set.molds[(parity, mirror)]
end

# --- assembly (chemistry-blind common, signed typed contrast) ----------------

function _cc_check_grid(set::CommonContrastSet, half_nm::Float64, step_nm::Float64)
    half_nm > 0 || throw(ArgumentError("half_nm must be positive"))
    step_nm > 0 || throw(ArgumentError("step_nm must be positive"))
    expected = length(collect(-half_nm:step_nm:half_nm))
    expected == set.grid_n || throw(ArgumentError(
        "stale grid: half_nm=$half_nm step_nm=$step_nm imply grid_n=$expected but molds have $(set.grid_n)"))
    return nothing
end

"""Assemble the chemistry-blind common image: superpose the common template for
each lobe using its parity. Sequence-independent. Pixel ordering, parity
formula, and bilinear sampling are identical to assemble_whole_roi."""
function assemble_common_image(set::CommonContrastSet, centers::AbstractVector,
                               xs::Vector{Float64}, ys::Vector{Float64};
                               half_nm::Float64, step_nm::Float64,
                               theta::Float64=0.0, direction::Int=0,
                               phase::Int=0, mirror::Int=0)
    _cc_check_grid(set, half_nm, step_nm)
    direction in (0, 1) || throw(ArgumentError("direction must be 0 or 1"))
    phase in (0, 1) || throw(ArgumentError("phase must be 0 or 1"))
    mirror in (0, 1) || throw(ArgumentError("mirror must be 0 or 1"))
    n = length(centers)
    n >= 1 || throw(ArgumentError("whole-ROI assembly requires at least one lobe"))
    cos_theta, sin_theta = cos(theta), sin(theta)
    assembled = zeros(Float64, length(ys), length(xs))
    for lobe in eachindex(centers)
        parity = Main._whole_roi_parity(lobe, n, direction, phase)
        mold = common_mold(set, parity, mirror)
        center_x, center_y = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx = xs[ix] - center_x
            dy = ys[iy] - center_y
            t = dx * cos_theta + dy * sin_theta
            u = -dx * sin_theta + dy * cos_theta
            assembled[iy, ix] += Main._sample_template(mold.common, t, u, half_nm, step_nm)
        end
    end
    return assembled
end

"""Assemble the signed typed-contrast image: for each lobe add +contrast for
type 1 and -contrast for type 0. Added to the chemistry-blind common image it
exactly reproduces assemble_whole_roi for the given sequence."""
function assemble_contrast_image(set::CommonContrastSet,
                                sequence::AbstractVector{<:Integer},
                                centers::AbstractVector, xs::Vector{Float64},
                                ys::Vector{Float64};
                                half_nm::Float64, step_nm::Float64,
                                theta::Float64=0.0, direction::Int=0,
                                phase::Int=0, mirror::Int=0)
    _cc_check_grid(set, half_nm, step_nm)
    direction in (0, 1) || throw(ArgumentError("direction must be 0 or 1"))
    phase in (0, 1) || throw(ArgumentError("phase must be 0 or 1"))
    mirror in (0, 1) || throw(ArgumentError("mirror must be 0 or 1"))
    length(sequence) == length(centers) || throw(ArgumentError("sequence/center length mismatch"))
    all(typ -> typ in (0, 1), sequence) || throw(ArgumentError("sequence values must be 0 or 1"))
    n = length(centers)
    cos_theta, sin_theta = cos(theta), sin(theta)
    assembled = zeros(Float64, length(ys), length(xs))
    for lobe in eachindex(sequence)
        parity = Main._whole_roi_parity(lobe, n, direction, phase)
        mold = common_mold(set, parity, mirror)
        sign = Int(sequence[lobe]) == 1 ? 1.0 : -1.0
        center_x, center_y = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx = xs[ix] - center_x
            dy = ys[iy] - center_y
            t = dx * cos_theta + dy * sin_theta
            u = -dx * sin_theta + dy * cos_theta
            assembled[iy, ix] += sign * Main._sample_template(mold.contrast, t, u, half_nm, step_nm)
        end
    end
    return assembled
end

# --- exhaustive active-set NNLS with a free intercept ------------------------

const _NNLS_NEG_TOL = sqrt(eps(Float64))
const _NNLS_MAX_PREDICTORS = 16  # 2^16 subsets; diagnostic scope never exceeds this.
const _NNLS_DEFAULT_RTOL = sqrt(eps(Float64))  # scale-aware tie tolerance by default.

"""Solve y ≈ background + sum_i coef_i * predictors[i] with a free intercept and
nonnegative predictor coefficients by enumerating every active set (predictor
subset). Three feasibility rules govern candidate admissibility:

  * rank-deficient designs (a predictor collinear with the intercept or another
    predictor) are skipped — their slope is unidentifiable, so the active set is
    infeasible rather than solved with an arbitrary least-norm pick;
  * active sets whose unrestricted solution strongly wants a negative predictor
    slope (< -sqrt(eps)) are rejected;
  * small numerical negatives are clamped to zero.

Among surviving candidates the minimum-SSE fit wins, with ties (within
atol + rtol*||y||^2) resolved in favor of the smallest active set enumerated
first (popcount, then lexicographic) — deterministic and parsimony-favoring.
The default rtol is sqrt(eps) so floating-point noise does not overturn the
parsimony rule; callers needing exact strict improvement pass rtol=0. The empty
active set (intercept only) is always full-rank and admissible, so a fit is
always returned for non-empty y. Rank deficiency is handled by skipping the
candidate (a documented feasibility decision), not by silently catching errors."""
function active_set_nnls(y::AbstractVector{<:Real},
                         predictors::AbstractVector{<:AbstractVector{<:Real}};
                         rtol::Real=_NNLS_DEFAULT_RTOL, atol::Real=0.0)
    rtol ≥ 0 || throw(ArgumentError("rtol must be non-negative"))
    atol ≥ 0 || throw(ArgumentError("atol must be non-negative"))
    n = length(y)
    n ≥ 1 || throw(ArgumentError("observations must be non-empty"))
    n_pred = length(predictors)
    n_pred ≤ _NNLS_MAX_PREDICTORS || throw(ArgumentError(
        "active-set enumeration capped at $_NNLS_MAX_PREDICTORS predictors, got $n_pred"))
    for p in predictors
        length(p) == n || throw(ArgumentError("predictor length mismatch"))
    end
    yv = Float64.(y)
    tss = sum(abs2, yv)                          # problem-scale reference for ties
    pv = [Float64.(p) for p in predictors]
    order = sort(collect(0:(2^n_pred - 1)); by=mask -> (count_ones(mask), mask))
    best = nothing
    for mask in order
        active = Int[i for i in 1:n_pred if (mask >> (i - 1)) & 1 == 1]
        n ≥ length(active) + 1 || continue        # skip underdetermined designs
        design = ones(Float64, n, length(active) + 1)
        for (col, idx) in enumerate(active)
            design[:, col + 1] .= pv[idx]
        end
        ncols = size(design, 2)
        # Skip rank-deficient (collinear) designs: the predictor slope is
        # unidentifiable, so the active set is infeasible. This is a documented
        # feasibility decision, not a swallowed error.
        rank(design) == ncols || continue
        coef = design \ yv
        if length(active) > 0 && any(c -> c < -_NNLS_NEG_TOL, coef[2:end])
            continue   # strongly negative predictor -> reject this active set
        end
        if length(active) > 0
            coef[2:end] .= max.(coef[2:end], 0.0)  # clamp small numerical negatives
        end
        residual = yv .- design * coef
        sse = sum(abs2, residual)
        full = zeros(Float64, n_pred)
        for (col, idx) in enumerate(active)
            full[idx] = coef[col + 1]
        end
        candidate = (sse=sse, background=coef[1], coefficients=full, active=active)
        if best === nothing
            best = candidate
        else
            tie_tol = Float64(atol) + Float64(rtol) * max(tss, 1.0)
            candidate.sse < best.sse - tie_tol && (best = candidate)
        end
    end
    best === nothing && throw(ArgumentError(
        "active-set NNLS failed: no admissible design"))
    return ActiveSetFit(best.sse, best.background, best.coefficients, best.active, n)
end

"""Matrix form: applies the finite-pixel good mask (matching
fit_whole_roi_nuisance) and delegates to the vector solver."""
function active_set_nnls(observed::AbstractMatrix{<:Real},
                         predictors::AbstractVector{<:AbstractMatrix{<:Real}};
                         rtol::Real=0.0, atol::Real=0.0)
    for p in predictors
        size(p) == size(observed) || throw(ArgumentError("predictor shape mismatch"))
    end
    good = isfinite.(observed)
    for p in predictors
        good .&= isfinite.(p)
    end
    ngood = count(good)
    needed = max(length(predictors) + 1, 1)
    ngood ≥ needed || throw(ArgumentError("insufficient finite ROI pixels"))
    yv = Float64.(observed[good])
    pv = [Float64.(p[good]) for p in predictors]
    return active_set_nnls(yv, pv; rtol=rtol, atol=atol)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Diagnostic common/contrast helper loaded. Run ",
            "test/joint_proxy/test_whole_roi_common_contrast.jl for the focused gate.")
end
