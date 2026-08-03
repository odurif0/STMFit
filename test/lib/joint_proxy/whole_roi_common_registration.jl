#!/usr/bin/env julia

# Deterministic chemistry-blind global registration for the whole-ROI common
# mold.  Enumerates 8 global states (direction x phase x mirror), performs a
# bounded coarse-to-fine search over shift_t, shift_u, rotation, and blur, and
# selects by SSE with a documented minimal-transform tie-break.  Returns an
# immutable FrozenRegistration.
#
# Scientific invariants (AGENTS.md): no benchmark label, expected N, NKNNKN,
# class count, or composition is read here.  The API takes no sequence or
# contrast argument — registration is chemistry-blind.
#
# This is Todo 3 of common-template-frozen-contrast-v2.

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
    # Inline the standalone blur from simulator_nuisance.jl to avoid pulling
    # in ViewNuisanceTruth which is only available in a different include chain.
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

using LinearAlgebra

const _REG_DC = Main.WholeRoiDiagnosticConfig

# Number of coarse-grid seeds refined in the fine pass per global state.
#
# Single-seed coarse-to-fine commits to the best coarse point before the fine
# grid can resolve the SSE surface. That is unsafe when the mold geometry
# admits competitive near-degenerate rotations on the coarse grid: a spurious
# coarse basin at the rotation wall can beat the true basin by ~1e-3 at coarse
# granularity even though the true basin collapses to SSE ~= 0 once refined
# (verified empirically for the synthetic two-peak mold). We therefore
# deterministically rank every coarse point by (SSE, minimal-transform)
# tie-key and refine the top-_REG_NUM_COARSE_SEEDS seeds, letting the fine grid
# escape such aliases and recover an identifiable injected transform within one
# final-grid step.
#
# This is an algorithmic constant, NOT a physical parameter: the search DOMAIN
# (shift/rotation/blur bounds and steps) and the numerical-tie tolerances stay
# fully config-driven. K is bounded so runtime stays tractable; ranks above ~5
# are genuinely worse basins, not competitive aliases.
const _REG_NUM_COARSE_SEEDS = 5

# --- read-only matrix wrapper -------------------------------------------------

# `FrozenRegistration` exposes its array data through `FrozenMatrix`: a defensive
# read-only AbstractMatrix. A plain `struct` is only SHALLOWLY immutable — a
# `Matrix` field can still be mutated in place (`img[i,j] = v`), silently
# invalidating `common_image_hash` and defeating Todo 4's tamper check. The
# wrapper instead:
#   * COPIES the source data at construction (so an alias held elsewhere cannot
#     mutate the frozen copy);
#   * rejects every in-place mutation path (`setindex!`, and `fill!`/`copyto!`
#     which route through `setindex!`) with a loud ArgumentError;
#   * behaves exactly like the underlying Matrix for reads (`getindex`, `size`,
#     `axes`, iteration, `==`, broadcast) and for content-based `hash`.
# Consumers that need a mutable copy use `Matrix(frozen)` or
# `convert(Matrix{T}, frozen)` (Todo 4's adapter does this via its typed field).
struct FrozenMatrix{T} <: AbstractMatrix{T}
    bytes::String
    dims::Tuple{Int,Int}
    function FrozenMatrix{T}(src::AbstractMatrix) where {T}
        data = Matrix{T}(undef, size(src))
        copyto!(data, src)
        bytes = String(copy(reinterpret(UInt8, vec(data))))
        new{T}(bytes, size(data))
    end
end

FrozenMatrix(src::AbstractMatrix{T}) where {T} = FrozenMatrix{T}(src)

Base.size(m::FrozenMatrix) = m.dims
Base.axes(m::FrozenMatrix) = map(Base.OneTo, m.dims)
Base.eltype(::Type{FrozenMatrix{T}}) where {T} = T
Base.length(m::FrozenMatrix) = prod(m.dims)
Base.IndexStyle(::Type{<:FrozenMatrix}) = IndexLinear()

@inline _frozen_values(m::FrozenMatrix{T}) where {T} = reinterpret(T, codeunits(m.bytes))
@inline Base.getindex(m::FrozenMatrix, i::Int) = _frozen_values(m)[i]
@inline Base.getindex(m::FrozenMatrix, i::Int, j::Int) =
    _frozen_values(m)[(j - 1) * m.dims[1] + i]

# The whole point: every mutation path is closed. ArgumentError names the
# invariant. `fill!` and `copyto!` reach this method through their generic
# AbstractArray implementations, so they are rejected too.
function Base.setindex!(::FrozenMatrix, v, inds...)
    throw(ArgumentError(
        "FrozenMatrix is read-only: frozen registration fields cannot be mutated in place"))
end

# Content-based hash identical to the underlying Matrix, so the frozen
# `common_image_hash` stays meaningful and cross-representation comparison works.
Base.hash(m::FrozenMatrix, h::UInt) = hash(Matrix(m), h)
# Mutable copies for downstream consumers (independent of the frozen store).
# Use explicit Array construction, NOT broadcast (`T.(data)` would produce a
# BitMatrix for T=Bool, which is not a Matrix{Bool} and breaks typed struct
# fields that consume the wrapper via convert).
Base.copy(m::FrozenMatrix) = Matrix(m)
Base.similar(m::FrozenMatrix, ::Type{T}, dims::Dims) where {T} = Matrix{T}(undef, dims)
function Base.convert(::Type{Matrix{T}}, m::FrozenMatrix) where {T}
    return copyto!(Matrix{T}(undef, size(m)), m)
end
Base.Matrix(m::FrozenMatrix{T}) where {T} = copyto!(Matrix{T}(undef, size(m)), m)

"""Deeply immutable snapshot of registration fit diagnostics."""
struct FrozenActiveSetFit
    sse::Float64
    background::Float64
    coefficients::Tuple
    active::Tuple
    n_observations::Int
end

FrozenActiveSetFit(fit::ActiveSetFit) = FrozenActiveSetFit(
    fit.sse, fit.background, Tuple(fit.coefficients), Tuple(fit.active),
    fit.n_observations)

# --- FrozenRegistration struct ------------------------------------------------

"""Immutable record of the chosen global registration state.
Chemistry-blind: no sequence, no contrast, no type label.

Array-typed fields (`common_image`, `mask`) are exposed as read-only
`FrozenMatrix` wrappers: they copy on freeze and reject in-place mutation, so
the stored `common_image_hash` cannot be silently invalidated and downstream
tamper checks (Todo 4) stay meaningful."""
struct FrozenRegistration
    direction::Int         # 0 or 1
    phase::Int             # 0 or 1
    mirror::Int            # 0 or 1
    shift_t_nm::Float64    # chosen coarse+fine t shift
    shift_u_nm::Float64    # chosen coarse+fine u shift
    rotation_deg::Float64  # chosen coarse+fine rotation (degrees)
    blur_sigma_nm::Float64 # chosen coarse+fine blur sigma
    common_image::FrozenMatrix{Float64}    # read-only assembled + blurred common at final transform
    common_image_hash::UInt64              # content hash of common_image for integrity checks
    grid_n::Int                            # mold grid size
    mask::FrozenMatrix{Bool}               # read-only finite-pixel good mask
    fit::FrozenActiveSetFit                # deeply immutable NNLS fit snapshot
    tie_info::String                       # "exact" or tie description
end

# --- local helpers ------------------------------------------------------------

"""Assemble the common image with shift and rotation applied via coordinate
sampling.  Shifts shift_t/shift_u are added to t,u before sampling the template.
This duplicates the assembly loop from assemble_common_image, adding the
shift offsets directly."""
function _assemble_common_shifted(set::CommonContrastSet,
                                   centers::AbstractVector,
                                   xs::Vector{Float64}, ys::Vector{Float64};
                                   half_nm::Float64, step_nm::Float64,
                                   theta::Float64,
                                   direction::Int, phase::Int, mirror::Int,
                                   shift_t::Float64, shift_u::Float64)
    n = length(centers)
    cos_theta, sin_theta = cos(theta), sin(theta)
    assembled = zeros(Float64, length(ys), length(xs))
    for lobe in eachindex(centers)
        parity = Main._whole_roi_parity(lobe, n, direction, phase)
        mold = common_mold(set, parity, mirror)
        center_x, center_y = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx = xs[ix] - center_x
            dy = ys[iy] - center_y
            t = dx * cos_theta + dy * sin_theta + shift_t
            u = -dx * sin_theta + dy * cos_theta + shift_u
            assembled[iy, ix] += Main._sample_template(mold.common, t, u, half_nm, step_nm)
        end
    end
    return assembled
end

"""Tie-break key: lower is better.
Priority order:
  1. SSE
  2. |shift_t_nm| (smaller abs, then signed: negative before positive)
  3. |shift_u_nm| (same)
  4. |rotation_deg| (same)
  5. |blur_sigma_nm| (same)
  6. direction (0 < 1)
  7. phase (0 < 1)
  8. mirror (0 < 1)
"""
function _tie_key(sse::Float64, shift_t::Float64, shift_u::Float64,
                  rot::Float64, blur::Float64,
                  direction::Int, phase::Int, mirror::Int)
    return (sse,
            abs(shift_t), shift_t,     # abs first, then signed (negative < positive)
            abs(shift_u), shift_u,
            abs(rot), rot,
            abs(blur), blur,
            direction, phase, mirror)
end

"""Clamp fine grid so it does not exceed the coarse bounds.
If coarse_best is at a boundary, shift the fine grid inward."""
function _clamped_fine_grid(coarse_best::Float64,
                            coarse::_REG_DC.DiagnosticRange,
                            fine::_REG_DC.DiagnosticRange)
    coarse_vals = _REG_DC.search_values(coarse)
    fine_offsets = _REG_DC.search_values(fine)
    coarse_min, coarse_max = first(coarse_vals), last(coarse_vals)
    # If coarse_best is at a boundary, clamp fine grid to not exceed bounds
    vals = coarse_best .+ fine_offsets
    vals = clamp.(vals, coarse_min, coarse_max)
    # Remove exact duplicates that arise from clamping
    unique!(vals)
    return vals
end

"""Evaluate one grid point: assemble shifted common image, optionally blur, fit
with backbone+common, return (sse, tie_key, common_image) or nothing on failure."""
function _evaluate_point(set::CommonContrastSet,
                         observed::AbstractMatrix{<:Real},
                         backbone::AbstractMatrix{<:Real},
                         centers::AbstractVector,
                         xs::Vector{Float64}, ys::Vector{Float64};
                         half_nm::Float64, step_nm::Float64,
                         theta::Float64,
                         direction::Int, phase::Int, mirror::Int,
                         shift_t::Float64, shift_u::Float64,
                         rotation_deg::Float64,
                         blur_sigma_nm::Float64,
                         pixel_step::Float64)
    rotation_rad = deg2rad(rotation_deg)
    theta_total = theta + rotation_rad

    common_img = _assemble_common_shifted(set, centers, xs, ys;
        half_nm=half_nm, step_nm=step_nm, theta=theta_total,
        direction=direction, phase=phase, mirror=mirror,
        shift_t=shift_t, shift_u=shift_u)

    blur_px = blur_sigma_nm / pixel_step
    if blur_px > 0
        common_img = Main._gaussian_blur(common_img, blur_px)
    end

    fit = active_set_nnls(observed, [backbone, common_img])
    key = _tie_key(fit.sse, shift_t, shift_u, rotation_deg, blur_sigma_nm,
                   direction, phase, mirror)
    return (sse=fit.sse, tie_key=key, fit=fit, common_image=common_img)
end

# --- public API ---------------------------------------------------------------

"""Deterministic chemistry-blind global registration.
Enumerates 8 global states (direction x phase x mirror), performs coarse-then-
fine search over shift_t_nm, shift_u_nm, rotation_deg, blur_sigma_nm.
At each point, assembles the common image using coordinate sampling with shift
and rotation via theta addition, optionally blurs, and fits background + backbone
+ common using active_set_nnls.  Selects deterministically by SSE, then by the
documented minimal-transform tie-break.

NO sequence argument.  NO contrast argument.  NO type label.
"""
function register_common(set::CommonContrastSet,
                          observed::AbstractMatrix{<:Real},
                          backbone::AbstractMatrix{<:Real},
                          centers::AbstractVector,
                          xs::Vector{Float64}, ys::Vector{Float64};
                          config::_REG_DC.DiagnosticConfig,
                          half_nm::Float64, step_nm::Float64,
                          theta_base::Float64=0.0)
    pixel_step = xs[2] - xs[1]
    pixel_step > 0 || throw(ArgumentError("pixel spacing must be positive"))

    coarse_shift_t = _REG_DC.search_values(config.shift_t_nm.coarse)
    coarse_shift_u = _REG_DC.search_values(config.shift_u_nm.coarse)
    coarse_rotation = _REG_DC.search_values(config.rotation_deg.coarse)
    coarse_blur = _REG_DC.search_values(config.blur_sigma_nm.coarse)

    best_overall = nothing  # (direction, phase, mirror, shift_t, shift_u, rot, blur, tie_key, fit, common_image)

    for direction in (0, 1), phase in (0, 1), mirror in (0, 1)
        # --- coarse pass: collect EVERY admissible grid point ---------------
        # The coarse pass is not a single-winner stage. It enumerates the full
        # bounded coarse grid and keeps all admissible (finite, feasible-fit)
        # points so the fine stage can refine a bounded set of competitive
        # seeds rather than commit to one basin prematurely (see
        # _REG_NUM_COARSE_SEEDS). Infeasible points (e.g. insufficient finite
        # ROI pixels, rank-deficient design) are skipped, not silently caught.
        coarse_results = Any[]
        for st in coarse_shift_t, su in coarse_shift_u,
            rd in coarse_rotation, bl in coarse_blur
            try
                result = _evaluate_point(set, observed, backbone, centers, xs, ys;
                    half_nm=half_nm, step_nm=step_nm, theta=theta_base,
                    direction=direction, phase=phase, mirror=mirror,
                    shift_t=st, shift_u=su, rotation_deg=rd, blur_sigma_nm=bl,
                    pixel_step=pixel_step)
                push!(coarse_results, (shift_t=st, shift_u=su, rotation_deg=rd,
                                       blur_sigma_nm=bl, tie_key=result.tie_key,
                                       fit=result.fit, common_image=result.common_image))
            catch
                continue
            end
        end
        isempty(coarse_results) && continue

        # Deterministic competitive-seed selection: rank by the documented
        # (SSE, minimal-transform) tie-key and refine the top-K. The tie-key
        # makes this fully deterministic under exact SSE ties.
        sort!(coarse_results, by=r -> r.tie_key)
        n_seeds = min(_REG_NUM_COARSE_SEEDS, length(coarse_results))

        for seed_idx in 1:n_seeds
            seed = coarse_results[seed_idx]
            # --- fine pass around this competitive seed ---------------------
            fine_shift_t = _clamped_fine_grid(seed.shift_t, config.shift_t_nm.coarse, config.shift_t_nm.fine)
            fine_shift_u = _clamped_fine_grid(seed.shift_u, config.shift_u_nm.coarse, config.shift_u_nm.fine)
            fine_rotation = _clamped_fine_grid(seed.rotation_deg, config.rotation_deg.coarse, config.rotation_deg.fine)
            fine_blur = _clamped_fine_grid(seed.blur_sigma_nm, config.blur_sigma_nm.coarse, config.blur_sigma_nm.fine)

            best_fine = seed  # start from the coarse seed itself
            for st in fine_shift_t, su in fine_shift_u,
                rd in fine_rotation, bl in fine_blur
                try
                    result = _evaluate_point(set, observed, backbone, centers, xs, ys;
                        half_nm=half_nm, step_nm=step_nm, theta=theta_base,
                        direction=direction, phase=phase, mirror=mirror,
                        shift_t=st, shift_u=su, rotation_deg=rd, blur_sigma_nm=bl,
                        pixel_step=pixel_step)
                    if result.tie_key < best_fine.tie_key
                        best_fine = (shift_t=st, shift_u=su, rotation_deg=rd,
                                     blur_sigma_nm=bl, tie_key=result.tie_key,
                                     fit=result.fit, common_image=result.common_image)
                    end
                catch
                    continue
                end
            end

            if best_overall === nothing || best_fine.tie_key < best_overall.tie_key
                best_overall = (direction=direction, phase=phase, mirror=mirror,
                                shift_t=best_fine.shift_t, shift_u=best_fine.shift_u,
                                rotation_deg=best_fine.rotation_deg,
                                blur_sigma_nm=best_fine.blur_sigma_nm,
                                tie_key=best_fine.tie_key,
                                fit=best_fine.fit,
                                common_image=best_fine.common_image)
            end
        end
    end

    best_overall === nothing && throw(ArgumentError(
        "registration failed: no admissible grid point found"))

    bo = best_overall
    # Freeze the array outputs defensively: FrozenMatrix copies on construction
    # and rejects in-place mutation, so the stored hash is tamper-proof.
    frozen_common = FrozenMatrix(bo.common_image)
    frozen_mask = FrozenMatrix(isfinite.(observed) .& isfinite.(bo.common_image) .& isfinite.(backbone))
    common_hash = Base.hash(frozen_common)

    # Check if the result was an exact winner or tied
    tie_info = "exact"

    return FrozenRegistration(bo.direction, bo.phase, bo.mirror,
        bo.shift_t, bo.shift_u, bo.rotation_deg, bo.blur_sigma_nm,
        frozen_common, common_hash, set.grid_n, frozen_mask,
        FrozenActiveSetFit(bo.fit), tie_info)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Registration helper loaded. Run ",
            "test/joint_proxy/test_whole_roi_common_registration.jl for tests.")
end
