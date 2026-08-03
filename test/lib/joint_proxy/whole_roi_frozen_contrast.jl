#!/usr/bin/env julia

# STUB for red-first TDD. Real structs + adapter so tests can construct inputs;
# _verify_frozen and score_frozen_contrast throw so the focused test run is RED.
# (Replaced by the real implementation after red evidence is captured.)
#
# Diagnostic-only frozen type-contrast sequence scoring.
# See .omo/plans/common-template-frozen-contrast-v2.md (Todo 4) and AGENTS.md.

if !isdefined(Main, :JointProxyRegistry)
    include(joinpath(@__DIR__, "proxy_registry.jl"))
end
if !isdefined(Main, :_whole_roi_parity)
    include(joinpath(@__DIR__, "..", "..", "diagnose_joint_proxy_whole_roi.jl"))
end
if !isdefined(Main, :CommonContrastSet)
    include(joinpath(@__DIR__, "whole_roi_common_contrast.jl"))
end
if !isdefined(Main, :WholeRoiDiagnosticConfig)
    include(joinpath(@__DIR__, "whole_roi_diagnostic_config.jl"))
end

if !isdefined(Main, :_gaussian_blur)
    function _gaussian_blur(z::Matrix{Float64}, sigma_px::Float64)
        sigma_px <= 0.0 && return copy(z)
        radius = max(1, ceil(Int, 3.0 * sigma_px))
        ny, nx = size(z)
        kernel = Float64[exp(-0.5 * (i / sigma_px)^2) for i in -radius:radius]
        kernel ./= sum(kernel)
        nk = length(kernel)
        tmp = similar(z)
        @inbounds for iy in 1:ny, ix in 1:nx
            s = 0.0
            for j in 1:nk
                idx = clamp(ix + j - 1 - radius, 1, nx)
                s += kernel[j] * z[iy, idx]
            end
            tmp[iy, ix] = s
        end
        out = similar(z)
        @inbounds for iy in 1:ny, ix in 1:nx
            s = 0.0
            for j in 1:nk
                idx = clamp(iy + j - 1 - radius, 1, ny)
                s += kernel[j] * tmp[idx, ix]
            end
            out[iy, ix] = s
        end
        return out
    end
end

const _FC_DC = Main.WholeRoiDiagnosticConfig

# --- frozen input contract ---------------------------------------------------

"""Narrow immutable frozen input contract for the type-contrast scorer.

Compatible with Todo 3's `FrozenRegistration` via `frozen_scorer_input(reg)`.
Carries exactly what the scorer needs to reassemble signed contrast under the
already-chosen global state and geometric transform, plus the chemistry-blind
common-only baseline SSE from registration. The scorer never re-registers and
never re-derives direction/phase/mirror/geometry: it verifies the frozen hashes
and grid identity, then enumerates sequences under the frozen state only."""
struct FrozenScorerInput
    direction::Int          # frozen global state (0 or 1)
    phase::Int              # frozen global state (0 or 1)
    mirror::Int             # frozen global state (0 or 1)
    theta_total::Float64    # effective total rotation (rad) used at registration
    shift_t_nm::Float64     # frozen chain-axis shift
    shift_u_nm::Float64     # frozen transverse shift
    blur_sigma_nm::Float64  # frozen blur sigma
    common_image::Matrix{Float64}   # frozen chemistry-blind common predictor
    common_image_hash::UInt64       # integrity hash of common_image
    mask::Matrix{Bool}              # frozen finite-pixel good mask
    mask_hash::UInt64               # integrity hash of mask
    grid_n::Int                     # frozen mold grid size
    common_only_sse::Float64        # baseline SSE from registration (bg+backbone+common)
end

"""Adapt any record with the Todo 3 `FrozenRegistration` field names into the
scorer's narrow contract. Duck-typed by name so the scorer does not import
Todo 3's source. `theta_base` defaults to 0.0 to match the registration default;
pass the same `theta_base` you gave `register_common` for a nonzero base angle."""
function frozen_scorer_input(reg; theta_base::Float64=0.0)
    theta_total = Float64(theta_base) + deg2rad(Float64(getproperty(reg, :rotation_deg)))
    common_image = getproperty(reg, :common_image)
    mask = getproperty(reg, :mask)
    return FrozenScorerInput(
        Int(getproperty(reg, :direction)),
        Int(getproperty(reg, :phase)),
        Int(getproperty(reg, :mirror)),
        theta_total,
        Float64(getproperty(reg, :shift_t_nm)),
        Float64(getproperty(reg, :shift_u_nm)),
        Float64(getproperty(reg, :blur_sigma_nm)),
        common_image,
        UInt64(getproperty(reg, :common_image_hash)),
        mask,
        Base.hash(mask),
        Int(getproperty(reg, :grid_n)),
        Float64(getproperty(reg, :fit).sse))
end

# --- typed results ------------------------------------------------------------

"""One sequence's frozen-contrast fit. `bits` is the deterministic integer
encoding of `sequence` (lobe 1 is the least-significant bit) used only for
tie-break ordering; it carries no chemical meaning."""
struct SequenceScore
    sequence::Vector{Int}
    bits::Int
    sse::Float64
    background::Float64
    backbone_coef::Float64
    common_coef::Float64
    contrast_coef::Float64
    n_pixels::Int
end

"""Full frozen type-contrast scoring result with ranking, margins, and named
abstention reasons. `abstain` is true iff `abstention_reasons` is non-empty;
reasons are structural only (zero contrast gain, zero contrast coefficient,
nonfinite evidence, best/runner tie, complement tie) — there is no arbitrary
chemistry-confidence floor anywhere in this struct or its construction."""
struct FrozenContrastResult
    best::SequenceScore
    runner_up::Union{Nothing,SequenceScore}
    complement::SequenceScore
    ranking::Vector{SequenceScore}
    common_only_sse::Float64
    incremental_contrast_gain::Float64
    best_runner_margin::Float64
    complement_margin::Float64
    abstain::Bool
    abstention_reasons::Vector{String}
    tie_rtol::Float64
    tie_atol::Float64
    n_sequence::Int
    direction::Int
    phase::Int
    mirror::Int
end

# --- local helpers (pure coordinate/template/blur; NOT registration) ----------
# These mirror the frozen-transform application used by registration's
# `_assemble_common_shifted`, but for the SIGNED contrast. They read the frozen
# transform parameters from the input record and apply them via the shared
# parity formula, bilinear sampler, and Gaussian blur. They never re-register,
# never enumerate states, and never consult a type sequence beyond the sign
# convention below.

function _check_grid(set::CommonContrastSet, half_nm::Real, step_nm::Real)
    hf = Float64(half_nm); st = Float64(step_nm)
    hf > 0 || throw(ArgumentError("half_nm must be positive"))
    st > 0 || throw(ArgumentError("step_nm must be positive"))
    expected = length(collect(-hf:st:hf))
    expected == set.grid_n || throw(ArgumentError(
        "stale grid: half_nm=$hf step_nm=$st imply grid_n=$expected but molds have $(set.grid_n)"))
    return nothing
end

"""Assemble the signed typed-contrast image under the FROZEN global state and
geometric transform. For each lobe, add +contrast for type 1 and -contrast for
type 0, sampled at the frozen (direction, phase, mirror) parity and the frozen
(theta_total, shift_t, shift_u) transform. Pixel ordering and parity formula are
identical to `assemble_contrast_image` and `_assemble_common_shifted`."""
function _assemble_contrast_frozen(set::CommonContrastSet,
                                   sequence::AbstractVector{<:Integer},
                                   centers::AbstractVector,
                                   xs::Vector{Float64}, ys::Vector{Float64};
                                   half_nm::Float64, step_nm::Float64,
                                   theta_total::Float64,
                                   direction::Int, phase::Int, mirror::Int,
                                   shift_t::Float64, shift_u::Float64)
    n = length(centers)
    cos_theta, sin_theta = cos(theta_total), sin(theta_total)
    assembled = zeros(Float64, length(ys), length(xs))
    for lobe in eachindex(sequence)
        parity = Main._whole_roi_parity(lobe, n, direction, phase)
        mold = common_mold(set, parity, mirror)
        sign = Int(sequence[lobe]) == 1 ? 1.0 : -1.0
        center_x, center_y = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx = xs[ix] - center_x
            dy = ys[iy] - center_y
            t = dx * cos_theta + dy * sin_theta + shift_t
            u = -dx * sin_theta + dy * cos_theta + shift_u
            assembled[iy, ix] += sign * Main._sample_template(
                mold.contrast, t, u, half_nm, step_nm)
        end
    end
    return assembled
end

"""Apply the frozen Gaussian blur to a signed image. Linear, so blurring the
signed contrast is identical to (blur(common+contrast) - blur(common)); the
frozen common image already carries its blur, so the contrast must be blurred
the same way for the decomposition to remain consistent."""
function _frozen_blur(img::Matrix{Float64}, blur_sigma_nm::Float64,
                      pixel_step::Float64)
    blur_px = blur_sigma_nm / pixel_step
    return blur_px > 0.0 ? Main._gaussian_blur(img, blur_px) : img
end

# --- frozen integrity verification -------------------------------------------

"""Verify the frozen record BEFORE scoring: common_image and mask hashes must
match the carried integrity hashes (catches post-adaptation tampering of either
the pixels or the hash), the grid identity must match the mold set, the global
state must be 0/1, the ROI shapes must agree, and the observation's finite-pixel
mask must match the frozen mask (catches a drifted observed/backbone passed
against a record frozen against different data)."""
function _verify_frozen(input::FrozenScorerInput, set::CommonContrastSet,
                        observed, backbone)
    Base.hash(input.common_image) == input.common_image_hash ||
        throw(ArgumentError(
            "frozen common_image hash mismatch: the common image or its stored hash was changed after freezing"))
    Base.hash(input.mask) == input.mask_hash ||
        throw(ArgumentError(
            "frozen mask hash mismatch: the mask or its stored hash was changed after freezing"))
    input.grid_n == set.grid_n ||
        throw(ArgumentError(
            "frozen grid identity mismatch: record grid_n=$(input.grid_n) but mold grid_n=$(set.grid_n)"))
    input.direction in (0, 1) ||
        throw(ArgumentError("frozen direction must be 0 or 1, got $(input.direction)"))
    input.phase in (0, 1) ||
        throw(ArgumentError("frozen phase must be 0 or 1, got $(input.phase)"))
    input.mirror in (0, 1) ||
        throw(ArgumentError("frozen mirror must be 0 or 1, got $(input.mirror)"))
    size(observed) == size(input.common_image) ||
        throw(ArgumentError("observed shape $(size(observed)) != frozen common shape $(size(input.common_image))"))
    size(backbone) == size(input.common_image) ||
        throw(ArgumentError("backbone shape $(size(backbone)) != frozen common shape $(size(input.common_image))"))
    size(input.mask) == size(input.common_image) ||
        throw(ArgumentError("mask shape $(size(input.mask)) != frozen common shape $(size(input.common_image))"))
    # Observation-context consistency: the finite-pixel mask implied by the
    # passed observed/backbone/common must equal the frozen mask, otherwise the
    # frozen common-only baseline SSE is invalid for this observation.
    auto_mask = isfinite.(observed) .& isfinite.(input.common_image) .& isfinite.(backbone)
    auto_mask == input.mask ||
        throw(ArgumentError(
            "observation mask drifted from the frozen mask: the observed/backbone finite-pixel set no longer matches the record frozen against the registration data"))
    return nothing
end

# --- frozen type-contrast sequence scoring -----------------------------------

const _FROZEN_MAX_N = 16  # 2^16 sequences; diagnostic scope never reaches this.

"""Score every binary sequence under the frozen global state and geometric
transform, fitting `background + backbone + frozen common + signed contrast`
with nonnegative predictor coefficients, and return the full ranking plus
explicit best / runner-up / complement metrics and named abstention reasons.

The scorer never re-registers and never re-derives direction/phase/mirror or the
geometric transform: all of these are read from the immutable `input` and applied
verbatim. There is no composition prior (every binary sequence is enumerated and
ranked by SSE alone), no expected-sequence assertion, and no arbitrary
chemistry-confidence floor: abstention fires ONLY for the named structural
conditions (zero contrast gain, zero contrast coefficient, nonfinite evidence, a
best/runner-up tie, or a complement tie) under the configured scale-aware
numerical tolerance."""
function score_frozen_contrast(input::FrozenScorerInput, set::CommonContrastSet,
                               observed, backbone, centers, xs, ys;
                               half_nm::Real, step_nm::Real,
                               config::_FC_DC.DiagnosticConfig, max_n::Int=_FROZEN_MAX_N)
    _verify_frozen(input, set, observed, backbone)
    _check_grid(set, half_nm, step_nm)
    n = length(centers)
    1 ≤ n ||
        throw(ArgumentError("frozen scoring requires n >= 1 lobes, got $n"))
    n ≤ max_n ||
        throw(ArgumentError("n=$n exceeds exhaustive diagnostic max_n=$max_n"))
    1 ≤ max_n ≤ _FROZEN_MAX_N ||
        throw(ArgumentError("max_n must be in [1, $_FROZEN_MAX_N], got $max_n"))

    pixel_step = xs[2] - xs[1]
    pixel_step > 0 ||
        throw(ArgumentError("pixel spacing must be positive, got $pixel_step"))

    half_f = Float64(half_nm)
    step_f = Float64(step_nm)
    rtol = Float64(config.tie_rtol)
    atol = Float64(config.tie_atol)

    scores = Vector{SequenceScore}(undef, 2^n)
    @inbounds for bits in 0:(2^n - 1)
        sequence = Int[(bits >> (i - 1)) & 1 for i in 1:n]
        contrast_img = _assemble_contrast_frozen(set, sequence, centers, xs, ys;
            half_nm=half_f, step_nm=step_f, theta_total=input.theta_total,
            direction=input.direction, phase=input.phase, mirror=input.mirror,
            shift_t=input.shift_t_nm, shift_u=input.shift_u_nm)
        contrast_img = _frozen_blur(contrast_img, input.blur_sigma_nm, pixel_step)
        fit = active_set_nnls(observed, [backbone, input.common_image, contrast_img])
        scores[bits + 1] = SequenceScore(sequence, bits, fit.sse, fit.background,
            fit.coefficients[1], fit.coefficients[2], fit.coefficients[3],
            fit.n_observations)
    end

    # Deterministic ranking: SSE ascending, then integer encoding ascending.
    sort!(scores; by=s -> (s.sse, s.bits))
    best = scores[1]
    runner_up = length(scores) ≥ 2 ? scores[2] : nothing
    complement_bits = xor(best.bits, 2^n - 1)
    complement = scores[findfirst(s -> s.bits == complement_bits, scores)]

    common_only_sse = input.common_only_sse
    gain = common_only_sse - best.sse
    best_runner_margin = runner_up === nothing ? Inf : (runner_up.sse - best.sse)
    complement_margin = complement.sse - best.sse

    # Named structural abstention reasons. Each uses the configured scale-aware
    # numerical tie tolerance (isapprox with rtol/atol). There is deliberately
    # no absolute confidence floor anywhere here.
    reasons = String[]
    nonfinite = !isfinite(best.sse) || !isfinite(gain) ||
        !isfinite(best_runner_margin) || !isfinite(complement_margin) ||
        !isfinite(best.contrast_coef) || !isfinite(common_only_sse)
    if nonfinite
        push!(reasons, "nonfinite_evidence")
    end
    if isfinite(common_only_sse) && isfinite(best.sse) &&
       isapprox(common_only_sse, best.sse; rtol=rtol, atol=atol)
        push!(reasons, "zero_contrast_gain")
    end
    if isfinite(best.contrast_coef) && best.contrast_coef == 0.0
        push!(reasons, "zero_contrast_coefficient")
    end
    if runner_up !== nothing && isfinite(best.sse) && isfinite(runner_up.sse) &&
       isapprox(best.sse, runner_up.sse; rtol=rtol, atol=atol)
        push!(reasons, "best_runner_tie")
    end
    if isfinite(best.sse) && isfinite(complement.sse) &&
       isapprox(best.sse, complement.sse; rtol=rtol, atol=atol)
        push!(reasons, "complement_tie")
    end

    return FrozenContrastResult(best, runner_up, complement, scores,
        common_only_sse, gain, best_runner_margin, complement_margin,
        !isempty(reasons), reasons, rtol, atol, n,
        input.direction, input.phase, input.mirror)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Frozen type-contrast helper loaded. Run ",
            "test/joint_proxy/test_whole_roi_frozen_contrast.jl for the focused gate.")
end
