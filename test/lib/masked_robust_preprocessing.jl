"""
Synthetic-first, opt-in background estimation. No production imports or input I/O.
Coordinates and heights are in nm; storage is z[y,x]. See `masked_background`.
"""
module MaskedRobustPreprocessing

using LinearAlgebra
using Statistics

export read_options, masked_background

struct Options
    huber_delta::Float64
    maxiter::Int
    coefficient_rtol::Float64
    weight_atol::Float64
    scale_floor_nm::Float64
    rank_rtol::Float64
    min_background_pixels::Int
    min_background_fraction::Float64
    min_row_background_pixels::Int
end

const OPTION_KEYS = String.(fieldnames(Options))
const NORMAL_MEDIAN_ABS = 0.6744897501960817 # Φ⁻¹(3/4); not a noise calibration

function _validate(o::Options)
    for key in (:huber_delta, :scale_floor_nm)
        v = getfield(o, key)
        isfinite(v) && v > 0 || throw(ArgumentError("$key must be finite and positive"))
    end
    for key in (:coefficient_rtol, :weight_atol)
        v = getfield(o, key)
        isfinite(v) && 0 <= v < 1 || throw(ArgumentError("$key must be in [0,1)"))
    end
    isfinite(o.rank_rtol) && 0 < o.rank_rtol < 1 ||
        throw(ArgumentError("rank_rtol must be in (0,1)"))
    isfinite(o.min_background_fraction) && 0 <= o.min_background_fraction <= 1 ||
        throw(ArgumentError("min_background_fraction must be in [0,1]"))
    for key in (:maxiter, :min_background_pixels, :min_row_background_pixels)
        getfield(o, key) >= 1 || throw(ArgumentError("$key must be positive"))
    end
    return o
end

"""
    read_options(config::AbstractDict)

Read all nine required keys from the explicit `[masked_robust_preprocessing]`
section of a parsed TOML dictionary. Other sections are ignored. Unknown keys
inside this section, missing keys, booleans as numbers, and invalid values are
errors. No physical or selection parameter has a hidden default.
"""
function read_options(config::AbstractDict)
    section = get(config, "masked_robust_preprocessing", nothing)
    section isa AbstractDict || throw(ArgumentError("missing [masked_robust_preprocessing] section"))
    Set(keys(section)) == Set(OPTION_KEYS) ||
        throw(ArgumentError("masked_robust_preprocessing requires exactly: " * join(OPTION_KEYS, ", ")))
    values = map(fieldnames(Options)) do key
        value = section[String(key)]
        if key in (:maxiter, :min_background_pixels, :min_row_background_pixels)
            value isa Integer && !(value isa Bool) && typemin(Int) <= value <= typemax(Int) ||
                throw(ArgumentError("$key must be an integer representable as Int"))
            return Int(value)
        end
        value isa Real && !(value isa Bool) || throw(ArgumentError("$key must be a real number"))
        return Float64(value)
    end
    return _validate(Options(values...))
end

function _coordinates(values, name)
    values isa AbstractVector{<:Real} || throw(ArgumentError("$name must be a real vector"))
    isempty(values) && throw(ArgumentError("$name must not be empty"))
    converted = Float64.(values)
    all(isfinite, converted) || throw(ArgumentError("$name must be finite in Float64"))
    all(converted[i] < converted[i+1] for i in 1:length(converted)-1) ||
        throw(ArgumentError("$name must be strictly increasing"))
    return converted
end

# Background-only extrema determine the numerical coordinate transformation.
function _normalize(values, used)
    lo, hi = extrema(used)
    center = lo / 2 + hi / 2
    scale = max(abs(lo - center), abs(hi - center))
    normalized = scale == 0 ? zeros(length(values)) : (values .- center) ./ scale
    return normalized, center, scale
end

# SVD of the weighted design, never normal equations or an inverse of A'A.
function _solve(design, response, weights, rank_rtol)
    rootweights = sqrt.(weights)
    weighted = design .* rootweights
    rhs = response .* rootweights
    if !all(isfinite, weighted) || !all(isfinite, rhs)
        return nothing, 0, :nonfinite_arithmetic
    end
    factor = svd(weighted; full=false)
    rank = count(s -> s > rank_rtol * maximum(factor.S), factor.S)
    rank == 3 || return nothing, rank, :rank_deficient
    beta = factor.V * ((factor.U' * rhs) ./ factor.S)
    all(isfinite, beta) || return nothing, rank, :nonfinite_arithmetic
    return beta, rank, :ok
end

function _weights(residuals, scale, delta)
    cutoff = delta * scale
    return [abs(r) <= cutoff ? 1.0 : cutoff / abs(r) for r in residuals]
end

# Huber rho in nm², with a FIXED cutoff; OLS uses half RSS on the same support.
function _objective(residuals, scale, delta, loss)
    loss == :ols && return sum(abs2, residuals) / 2
    cutoff = delta * scale
    return sum(abs(r) <= cutoff ? r^2 / 2 : cutoff * (abs(r) - cutoff / 2) for r in residuals)
end

"""
    masked_background(xs, ys, z_nm, background_mask, options; loss=:huber)

Fit a plane to **only** `background_mask .& isfinite.(z_nm)`, then subtract the
median plane residual of each sufficiently supported row. There is no mask
selection, filling, interpolation, smoothing, registration, or fit-engine call.
Inputs are not modified. `z_nm` and the Boolean mask must have size
`(length(ys), length(xs))`. Coordinates must be finite and strictly increasing.
Support fraction uses the entire input rectangle as its denominator, including
missing pixels. Support counts and rank are checked before estimating a plane.

The Huber variant fixes its scale once: normalized MAD of background-only OLS
plane residuals, floored by the explicit `scale_floor_nm`. This is NOT calibrated
noise. The OLS initializer can inflate that scale under contamination; fixing it
avoids claiming convergence for an unanalysed parameter-dependent scale update.
IRLS starts at OLS and uses weighted SVD with a relative singular-value rank
check on every step. Convergence needs both (a) relative coefficient change in
the centered/scaled design, using max(old norm, new norm, robust scale), and
(b) maximum absolute weight change. Heights are median-centered before solving.
`iterations` counts IRLS weighted solves (OLS reports one solve). Reaching the
limit without BOTH criteria is `:NONCONVERGED/:maxiter`, even with finite values.

The sequential plane and arbitrary row offsets are not a unique decomposition
of generating y-tilt and row shifts. The plane is estimated first, independently
of offsets. Each usable row's observed background corrected median is zero; no
global nonnegative shift is applied. Unequal within-row x support can correlate
row offsets with x and bias the plane's x slope even with a correct exclusion
footprint. Zero corrected row medians do not certify a correct total background.
This convention also removes row-constant signal that was not excluded by the
caller. Huber handles background outliers, not arbitrary foreground omitted
from the exclusion footprint.

Returned fields:
- `corrected`, `background`: corrected observations and estimated total
  background, respectively. Both are NaN at raw missing pixels or unusable rows.
- `plane`: estimated plane on the coordinate grid, including unobserved pixels;
  these values are model predictions, NEVER imputed measurements.
- `row_offsets`, `row_status`, `row_background_counts`: row medians, availability
  reasons, and counts on the supplied finite background support.
- `observed`: raw finite mask; `valid`: corrected availability mask;
  `background_mask`: exact caller mask intersected with raw finite support (also
  retained if estimation fails; availability/status must be checked separately).
- `coefficients=(intercept_nm, slope_x, slope_y)`: physical plane coefficients.
- `status`, `reason`, `converged`, `iterations`, `robust_scale`,
  `coefficient_change`, `weight_change`, `scale_floor_active`, `design_rank`,
  `loss`: numerical diagnostics, not chemical or statistical confidence.
- `initial_objective_nm2`, `objective_nm2`: initial/final fixed-scale Huber
  objective (OLS: half RSS), on the supplied finite background only.
  `stationarity_inf_nm` is max|A' ψ(r)| in the centered/scaled design, with
  ψ(r)=clamp(r, -delta*scale, delta*scale) for Huber, ψ(r)=r for OLS. It is a
  raw gradient diagnostic, not a calibrated score or an extra convergence gate.

`:OK` means the plane converged and every row is supported; it does not imply
complete raw observations. `:PARTIAL` retains unsupported rows as unavailable.
Global plane failure returns `:UNAVAILABLE`. Nonconvergence retains its final
plane/coefficients for diagnosis but no corrected/background observations or
row offsets. Global failures do not fall back to unmasked or filled data.
"""
function masked_background(xs, ys, z_nm, background_mask, options::Options; loss=:huber)
    _validate(options)
    loss in (:ols, :huber) || throw(ArgumentError("loss must be :ols or :huber"))
    x = _coordinates(xs, "xs")
    y = _coordinates(ys, "ys")
    z_nm isa AbstractMatrix{<:Real} || throw(ArgumentError("z_nm must be a real matrix"))
    background_mask isa AbstractMatrix{Bool} || throw(ArgumentError("background_mask must be Boolean"))
    size(z_nm) == (length(y), length(x)) || throw(DimensionMismatch("z_nm must be z[y,x]"))
    size(background_mask) == size(z_nm) || throw(DimensionMismatch("background_mask must match z_nm"))
    Base.require_one_based_indexing(xs, ys, z_nm, background_mask)

    # No values from excluded pixels enter the response or its scale estimate.
    observed = isfinite.(z_nm)
    used = BitMatrix(background_mask .& observed)
    row_background_counts = vec(sum(used; dims=2))
    corrected = fill(NaN, size(z_nm))
    background = fill(NaN, size(z_nm))
    plane = fill(NaN, size(z_nm))
    row_offsets = fill(NaN, length(y))
    row_status = fill(:plane_unavailable, length(y))
    valid = falses(size(z_nm))
    coefficients = (intercept_nm=NaN, slope_x=NaN, slope_y=NaN)
    converged, iterations, design_rank = false, 0, 0
    robust_scale, coefficient_change, weight_change = NaN, NaN, NaN
    initial_objective_nm2, objective_nm2, stationarity_inf_nm = NaN, NaN, NaN
    scale_floor_active = false
    finish(status, reason) = (; corrected, background, plane, row_offsets, valid,
        background_mask=used, status, reason, converged, iterations, coefficients,
        robust_scale, row_status, row_background_counts, observed, coefficient_change,
        weight_change, scale_floor_active, design_rank, loss,
        initial_objective_nm2, objective_nm2, stationarity_inf_nm)

    indices = findall(used)
    length(indices) >= options.min_background_pixels ||
        return finish(:UNAVAILABLE, :insufficient_background_pixels)
    length(indices) / length(z_nm) >= options.min_background_fraction ||
        return finish(:UNAVAILABLE, :insufficient_background_fraction)

    xn, xc, xscale = _normalize(x, [x[i[2]] for i in indices])
    yn, yc, yscale = _normalize(y, [y[i[1]] for i in indices])
    all(isfinite, xn) && all(isfinite, yn) || return finish(:UNAVAILABLE, :nonfinite_arithmetic)
    design = [ones(length(indices)) [xn[i[2]] for i in indices] [yn[i[1]] for i in indices]]
    heights = [Float64(z_nm[i]) for i in indices]
    zcenter = median(heights)
    response = heights .- zcenter
    beta, design_rank, solve_reason = _solve(design, response, ones(length(indices)), options.rank_rtol)
    solve_reason == :ok || return finish(:UNAVAILABLE, solve_reason)
    residuals = response - design * beta
    all(isfinite, residuals) || return finish(:UNAVAILABLE, :nonfinite_arithmetic)
    mad_scale = median(abs.(residuals .- median(residuals))) / NORMAL_MEDIAN_ABS
    isfinite(mad_scale) || return finish(:UNAVAILABLE, :nonfinite_arithmetic)
    scale_floor_active = mad_scale <= options.scale_floor_nm
    robust_scale = max(mad_scale, options.scale_floor_nm)
    initial_objective_nm2 = _objective(residuals, robust_scale, options.huber_delta, loss)

    if loss == :huber
        weights = _weights(residuals, robust_scale, options.huber_delta)
        for step in 1:options.maxiter
            iterations = step
            updated, design_rank, solve_reason = _solve(design, response, weights, options.rank_rtol)
            solve_reason == :ok || return finish(:UNAVAILABLE,
                solve_reason == :rank_deficient ? :weighted_rank_deficient : solve_reason)
            residuals = response - design * updated
            all(isfinite, residuals) || return finish(:UNAVAILABLE, :nonfinite_arithmetic)
            newweights = _weights(residuals, robust_scale, options.huber_delta)
            coefficient_change = norm(updated - beta) / max(norm(beta), norm(updated), robust_scale)
            weight_change = maximum(abs.(newweights - weights))
            beta, weights = updated, newweights
            if coefficient_change <= options.coefficient_rtol && weight_change <= options.weight_atol
                converged = true
                break
            end
        end
    else
        converged, iterations = true, 1
        coefficient_change, weight_change = 0.0, 0.0
    end

    objective_nm2 = _objective(residuals, robust_scale, options.huber_delta, loss)
    psi = loss == :huber ? clamp.(residuals, -options.huber_delta * robust_scale,
                                options.huber_delta * robust_scale) : residuals
    stationarity_inf_nm = maximum(abs.(design' * psi))

    # Positive xscale/yscale follow from the full-rank design.
    slope_x, slope_y = beta[2] / xscale, beta[3] / yscale
    coefficients = (intercept_nm=zcenter + beta[1] - slope_x * xc - slope_y * yc,
                    slope_x=slope_x, slope_y=slope_y)
    all(isfinite, coefficients) || return finish(:UNAVAILABLE, :nonfinite_coefficients)
    for j in eachindex(y), i in eachindex(x)
        plane[j,i] = zcenter + beta[1] + beta[2] * xn[i] + beta[3] * yn[j]
    end
    if !all(isfinite, plane)
        fill!(plane, NaN)
        return finish(:UNAVAILABLE, :nonfinite_plane)
    end
    converged || return finish(:NONCONVERGED, :maxiter)

    for j in eachindex(y)
        if row_background_counts[j] < options.min_row_background_pixels
            row_status[j] = :insufficient_background
            continue
        end
        residual_row = [Float64(z_nm[j,i]) - plane[j,i] for i in eachindex(x) if used[j,i]]
        offset = median(residual_row)
        observed_columns = findall(view(observed, j, :))
        corrected_row = [Float64(z_nm[j,i]) - plane[j,i] - offset for i in observed_columns]
        background_row = [plane[j,i] + offset for i in observed_columns]
        if !isfinite(offset) || !all(isfinite, corrected_row) || !all(isfinite, background_row)
            row_status[j] = :nonfinite_row_arithmetic
            continue
        end
        row_offsets[j] = offset
        row_status[j] = :ok
        for (k, i) in enumerate(observed_columns)
            corrected[j,i], background[j,i] = corrected_row[k], background_row[k]
            valid[j,i] = true
        end
    end
    n_supported = count(==(:ok), row_status)
    n_supported == 0 && return finish(:UNAVAILABLE, :no_supported_rows)
    any(==(:nonfinite_row_arithmetic), row_status) && return finish(:PARTIAL, :nonfinite_row_arithmetic)
    n_supported < length(y) && return finish(:PARTIAL, :insufficient_row_background)
    return finish(:OK, :ok)
end

end # module
