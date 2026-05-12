module STMFitCore

export FWHM_TO_SIGMA,
       sigma_from_fwhm, fwhm_from_sigma,
       PhysicalChainConstraints,
       effective_spacing_min, support_n_bounds, chain_can_fit_support,
       endpoint_overrun, best_valid_or_best,
       overlap_condition_number, kappa_penalty, adjacent_kappa_max

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

end # module
