module STMFitCore

export FWHM_TO_SIGMA,
       sigma_from_fwhm, fwhm_from_sigma,
       PhysicalChainConstraints,
       effective_spacing_min, support_n_bounds, chain_can_fit_support,
       endpoint_overrun, best_valid_or_best,
       overlap_condition_number, kappa_penalty, adjacent_kappa_max,
       durbin_watson, runs_test, ResidualDiagnostics, compute_residual_diagnostics

const FWHM_TO_SIGMA = 2 * sqrt(2 * log(2))

sigma_from_fwhm(fwhm::Real) = Float64(fwhm) / FWHM_TO_SIGMA
fwhm_from_sigma(sigma::Real) = Float64(sigma) * FWHM_TO_SIGMA

Base.@kwdef struct PhysicalChainConstraints
    spacing_min_nm::Float64 = 0.35
    spacing_max_nm::Float64 = 0.75
    fwhm_min_nm::Float64 = 0.45
    fwhm_max_nm::Float64 = 1.20
    max_overlap::Float64 = 0.60
    endpoint_tolerance_nm::Float64 = 1e-6
end

function effective_spacing_min(spacing_min_nm::Real, spacing_max_nm::Real,
                               sigma_max_nm::Real, max_overlap::Real)
    if !(0.0 < max_overlap < 1.0)
        return Float64(spacing_min_nm)
    end
    overlap_spacing = sqrt(-2.0 * log(Float64(max_overlap))) * Float64(sigma_max_nm)
    return min(max(Float64(spacing_min_nm), overlap_spacing), Float64(spacing_max_nm))
end

function effective_spacing_min(c::PhysicalChainConstraints)
    return effective_spacing_min(c.spacing_min_nm, c.spacing_max_nm,
                                 sigma_from_fwhm(c.fwhm_max_nm), c.max_overlap)
end

function support_n_bounds(support_length_nm::Real, spacing_min_eff_nm::Real, spacing_max_nm::Real;
                          n_min_config::Int=2, n_max_config::Int=typemax(Int))
    support = max(Float64(support_length_nm), 0.0)
    n_max_data = max(1, Int(floor(support / max(Float64(spacing_min_eff_nm), eps(Float64)))) + 1)
    n_min_data = max(2, Int(floor(support / max(Float64(spacing_max_nm), eps(Float64)))))
    n_min_eff = max(n_min_config, n_min_data)
    n_max_eff = min(n_max_config, n_max_data)
    if n_min_eff > n_max_eff
        n_min_eff = n_max_eff
    end
    return n_min_eff, n_max_eff
end

function support_n_bounds(support_length_nm::Real, c::PhysicalChainConstraints;
                          n_min_config::Int=2, n_max_config::Int=typemax(Int))
    return support_n_bounds(support_length_nm, effective_spacing_min(c), c.spacing_max_nm;
                            n_min_config=n_min_config, n_max_config=n_max_config)
end

chain_can_fit_support(n::Integer, support_length_nm::Real, spacing_min_eff_nm::Real) =
    max(n - 1, 0) * Float64(spacing_min_eff_nm) <= Float64(support_length_nm) + 1e-9

function chain_can_fit_support(n::Integer, support_length_nm::Real, c::PhysicalChainConstraints)
    return chain_can_fit_support(n, support_length_nm, effective_spacing_min(c))
end

function endpoint_overrun(ts::AbstractVector{<:Real}, tmin::Real, tmax::Real)
    isempty(ts) && return 0.0
    return max(Float64(tmin) - minimum(ts), 0.0) + max(maximum(ts) - Float64(tmax), 0.0)
end

function best_valid_or_best(results; score = r -> getproperty(r, :bic))
    valid = [r for r in results if getproperty(r, :success) && getproperty(r, :valid) && isfinite(score(r))]
    if !isempty(valid)
        return sort(valid; by=score)[1]
    end
    ok = [r for r in results if getproperty(r, :success) && isfinite(score(r))]
    isempty(ok) && return nothing
    return sort(ok; by=score)[1]
end

# ===========================================================================
# Overlap condition number κ
# ===========================================================================

"""
    overlap_condition_number(d, sigma) -> Float64

Condition number κ = (1+ρ)/(1-ρ) of the 2×2 Gram matrix for two unit-amplitude
Gaussians separated by `d` with width `sigma`. ρ = exp(-d²/(4σ²)) is their
normalized correlation. Returns 1.0 (well-separated) to Inf (coincident).
"""
function overlap_condition_number(d::Real, sigma::Real)
    sigma <= 0 && return Inf
    d <= 0 && return Inf
    rho = exp(-0.25 * (Float64(d) / Float64(sigma))^2)
    rho >= 1.0 && return Inf
    return (1.0 + rho) / (1.0 - rho)
end

"""
    kappa_penalty(kappa; kappa_max=25.0, weight=1.0) -> Float64

Progressive penalty on condition number: zero for κ ≤ κ_max,
then `weight × ((κ - κ_max)/κ_max)²`.
"""
function kappa_penalty(kappa::Real; kappa_max::Real=25.0, weight::Real=1.0)
    Float64(kappa_max) <= 0 && return 0.0
    Float64(kappa) <= Float64(kappa_max) && return 0.0
    return Float64(weight) * ((Float64(kappa) - Float64(kappa_max)) / Float64(kappa_max))^2
end

"""
    adjacent_kappa_max(deltas, sigmas) -> Float64

Maximum condition number across all adjacent peak pairs.
Uses σ = max(σ_i, σ_{i+1}) for each pair.
"""
function adjacent_kappa_max(deltas::AbstractVector{<:Real}, sigmas::AbstractVector{<:Real})
    isempty(deltas) && return 1.0
    length(sigmas) < length(deltas) + 1 && return Inf
    kmax = 1.0
    for i in eachindex(deltas)
        sigma_local = max(Float64(sigmas[i]), Float64(sigmas[i + 1]))
        k = overlap_condition_number(deltas[i], sigma_local)
        k > kmax && (kmax = k)
    end
    return kmax
end

# ===========================================================================
# Residual diagnostics
# ===========================================================================

# Normal CDF via Abramowitz & Stegun approximation (no external dependency)
function _norm_cdf(x::Float64)
    x < -8.0 && return 0.0
    x > 8.0 && return 1.0
    t = 1.0 / (1.0 + 0.2316419 * abs(x))
    d = 0.3989422804014327  # 1/√(2π)
    p = d * exp(-0.5 * x * x) *
        ((((1.330274429 * t - 1.821255978) * t + 1.781477937) * t -
          0.356563782) * t + 0.319381530) * t
    return x > 0 ? 1.0 - p : p
end

"""
    durbin_watson(residuals) -> (dw_stat, p_value)

Durbin-Watson statistic for residual autocorrelation.
DW ≈ 2 → no autocorrelation. DW < 1.5 → positive autocorrelation (missed structure).
Returns (NaN, NaN) for fewer than 5 points.
"""
function durbin_watson(residuals::AbstractVector{<:Real})
    n = length(residuals)
    n < 5 && return (NaN, NaN)
    rss = sum(abs2, residuals)
    rss <= 0 && return (NaN, NaN)
    dw = sum(abs2, diff(residuals)) / rss
    z = (dw - 2.0) / (2.0 / sqrt(n))
    p = 2.0 * min(1.0, _norm_cdf(-abs(z)))
    return (Float64(dw), p)
end

"""
    runs_test(residuals) -> (n_runs, expected, p_value)

Wald-Wolfowitz runs test on residual signs.
Too few runs → systematic bias. Too many → oscillation/overfitting.
"""
function runs_test(residuals::AbstractVector{<:Real})
    n = length(residuals)
    n < 5 && return (0, NaN, NaN)
    signs = residuals .> 0
    n_pos = count(signs)
    n_neg = n - n_pos
    (n_pos == 0 || n_neg == 0) && return (1, NaN, NaN)
    n_runs = 1
    for i in 2:n
        signs[i] != signs[i-1] && (n_runs += 1)
    end
    expected = 1.0 + 2.0 * n_pos * n_neg / n
    var_runs = 2.0 * n_pos * n_neg * (2*n_pos*n_neg - n) / (n*n*(n-1))
    z = var_runs > 0 ? (n_runs - expected) / sqrt(var_runs) : 0.0
    p = 2.0 * min(1.0, _norm_cdf(-abs(z)))
    return (n_runs, expected, p)
end

Base.@kwdef struct ResidualDiagnostics
    durbin_watson::Float64     = NaN
    durbin_watson_p::Float64   = NaN
    runs_n::Int                = 0
    runs_expected::Float64     = NaN
    runs_p::Float64            = NaN
    residual_rms::Float64      = NaN
    residual_max::Float64      = NaN
end

"""
    compute_residual_diagnostics(residuals) -> ResidualDiagnostics

Compute all residual diagnostics: Durbin-Watson, runs test, RMS, max.
"""
function compute_residual_diagnostics(residuals::AbstractVector{<:Real})
    dw, dw_p = durbin_watson(residuals)
    rn, re, rp = runs_test(residuals)
    return ResidualDiagnostics(
        durbin_watson=dw, durbin_watson_p=dw_p,
        runs_n=rn, runs_expected=re, runs_p=rp,
        residual_rms=sqrt(sum(abs2, residuals) / length(residuals)),
        residual_max=maximum(abs, residuals),
    )
end

end # module
