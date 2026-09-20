"""Final hard-cluster covariance only; no experimental labels or composition prior.

Ledoit-Wolf shrinkage uses centered p-by-n observations, the maximum-likelihood
sample covariance S, and the spherical target tr(S)/p I. The closed-form mixing
coefficient is computed from the same observations, never selected by a grade.
This is a numerical regularizer, not an iid noise or chemical-confidence claim:
lobes from a scan are dependent and cluster memberships are learned from data.
"""
module AssignmentCovariance

using LinearAlgebra, TOML
export load_covariance_config, final_covariance

function load_covariance_config(config::AbstractDict)
    model = get(config, "model", Dict())
    mode = get(model, "gmm_final_covariance", nothing)
    mode in ("ridge", "ledoit_wolf") || throw(ArgumentError("explicit gmm_final_covariance required"))
    ridge = get(model, "gmm_covariance_ridge", nothing)
    ridge isa Real && !(ridge isa Bool) && isfinite(ridge) && ridge > 0 ||
        throw(ArgumentError("explicit positive finite gmm_covariance_ridge required"))
    return (mode=String(mode), ridge=Float64(ridge))
end
load_covariance_config(path::AbstractString) = load_covariance_config(TOML.parsefile(path))

function final_covariance(centered::AbstractMatrix; mode::AbstractString, ridge::Real)
    p, n = size(centered)
    p > 0 && n > 0 || throw(ArgumentError("nonempty covariance observations required"))
    all(isfinite, centered) || throw(ArgumentError("nonfinite covariance input"))
    mode in ("ridge", "ledoit_wolf") || throw(ArgumentError("unknown covariance mode"))
    !(ridge isa Bool) && isfinite(ridge) && ridge > 0 || throw(ArgumentError("positive finite ridge required"))
    sample = centered * transpose(centered) / n
    all(isfinite, sample) || throw(ArgumentError("covariance overflow"))
    shrinkage = 0.0
    selected = sample
    if mode == "ledoit_wolf"
        target_scale = tr(sample) / p
        target = target_scale * Matrix{Float64}(I, p, p)
        dispersion = sum(abs2, sample - target)
        # Equivalent to sum_i ||x_i*x_i' - S||_F^2 / n^2; use sufficient
        # statistics rather than allocating n outer products.
        fourth = sum(sum(abs2, x)^2 for x in eachcol(centered))
        variance = max(0.0, (fourth / n - sum(abs2, sample)) / n)
        isfinite(dispersion) && isfinite(variance) || throw(ArgumentError("shrinkage overflow"))
        # If S already equals the spherical target, all mixing coefficients
        # give the same estimate. Report zero, including the one-point case.
        shrinkage = dispersion > 0 ? clamp(variance / dispersion, 0.0, 1.0) : 0.0
        selected = (1 - shrinkage) * sample + shrinkage * target
    end
    covariance = selected + ridge * I
    return (covariance=covariance, sample=sample, shrinkage=shrinkage)
end

end # module
