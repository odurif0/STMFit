# Shared two-component diagonal Gaussian emissions and one-component
# companion, with FIXED class priors (0.5, 0.5), covariance floors,
# monotone-likelihood checks, and deterministic multi-start EM.

const cov_floor_default = 1e-4
const default_n_starts = 5

function _positive_integer_parameter(name::AbstractString, value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer (Bool is not allowed)"))
    0 < value <= typemax(Int) ||
        throw(ArgumentError("$name must be a positive integer representable as Int"))
    return Int(value)
end

function _nonnegative_integer_parameter(name::AbstractString, value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a non-negative integer (Bool is not allowed)"))
    0 <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be a non-negative integer representable as Int"))
    return Int(value)
end

function _positive_finite_parameter(name::AbstractString, value)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive finite real (Bool is not allowed)"))
    converted = try
        Float64(value)
    catch
        throw(ArgumentError("$name must be a positive finite real representable as Float64"))
    end
    isfinite(converted) && converted > 0.0 ||
        throw(ArgumentError("$name must be a positive finite real representable as Float64"))
    return converted
end

struct TwoComponentFit
    means::Matrix{Float64}      # K × p
    vars::Matrix{Float64}       # K × p
    weights::Vector{Float64}    # K, always [0.5, 0.5]
    loglik::Float64             # final total log-likelihood
    loglik_trace::Vector{Float64}
    converged::Bool
    monotone::Bool
end

function _m_step(Z::AbstractMatrix, resp::Matrix{Float64}, cov_floor::Float64)
    n, p = size(Z)
    K = size(resp, 2)
    means = zeros(K, p)
    vars = zeros(K, p)
    for k in 1:K
        Nk = sum(view(resp, :, k))
        if Nk <= 1e-12
            means[k, :] .= 0.0
            vars[k, :] .= cov_floor
        else
            for j in 1:p
                zj = view(Z, :, j)
                rk = view(resp, :, k)
                μ = sum(rk[i] * zj[i] for i in 1:n) / Nk
                means[k, j] = μ
                v = sum(rk[i] * (zj[i] - μ)^2 for i in 1:n) / Nk
                vars[k, j] = max(v, cov_floor)
            end
        end
    end
    return means, vars
end

function _run_em(Z::AbstractMatrix, means0::Matrix{Float64},
                 vars0::Matrix{Float64}, weights::Vector{Float64},
                 cov_floor::Float64, max_iter::Int, tol::Float64)
    n, p = size(Z)
    K = length(weights)
    means = copy(means0)
    vars = copy(vars0)
    accepted_means = copy(means)
    accepted_vars = copy(vars)
    trace = Float64[]
    converged = false
    monotone = true
    prev_ll = _total_log_likelihood(Z, means, vars, weights)
    push!(trace, prev_ll)
    for _ in 1:max_iter
        lr = _log_resp_matrix(Z, means, vars, weights)
        _normalize_log_resp!(lr)
        resp = exp.(lr)
        means, vars = _m_step(Z, resp, cov_floor)
        ll = _total_log_likelihood(Z, means, vars, weights)
        if ll + 1e-9 < prev_ll
            monotone = false
            break
        end
        accepted_means = copy(means)
        accepted_vars = copy(vars)
        push!(trace, ll)
        if abs(ll - prev_ll) <= tol * (1.0 + abs(prev_ll))
            converged = true
            prev_ll = ll
            break
        end
        prev_ll = ll
    end
    return accepted_means, accepted_vars, trace, converged, monotone
end

# Deterministic initialization grid. Init candidates span quantile pairs of
# increasing extremity plus a robust min/max start, so heavily unbalanced
# clusters (e.g. 50 vs 4 lobes) still seed the two means in distinct modes.
function _init_candidates(Z::AbstractMatrix, K::Int, n_starts::Int)
    n, p = size(Z)
    cols = [filter(isfinite, Z[:, j]) for j in 1:p]
    starts = Tuple{Matrix{Float64},Matrix{Float64}}[]
    base_var = 1.0
    quantile_deltas = [0.25, 0.40, 0.45, 0.475, 0.49, 0.495]
    n_quantile = min(n_starts, length(quantile_deltas))
    for s in 1:n_quantile
        delta = quantile_deltas[s]
        lo_q = max(0.005, 0.5 - delta)
        hi_q = min(0.995, 0.5 + delta)
        means = zeros(K, p)
        vars = fill(base_var, K, p)
        for j in 1:p
            col = cols[j]
            isempty(col) && continue
            means[1, j] = quantile(col, lo_q)
            means[2, j] = quantile(col, hi_q)
            if means[1, j] > means[2, j]
                means[1, j], means[2, j] = means[2, j], means[1, j]
            end
        end
        push!(starts, (copy(means), copy(vars)))
    end
    if n_starts > n_quantile
        means = zeros(K, p)
        vars = fill(base_var, K, p)
        for j in 1:p
            col = cols[j]
            isempty(col) && continue
            lo = minimum(col)
            hi = maximum(col)
            means[1, j] = lo
            means[2, j] = hi
            if lo == hi
                means[1, j] = lo - 1.0
                means[2, j] = hi + 1.0
            end
        end
        push!(starts, (copy(means), copy(vars)))
    end
    while length(starts) < n_starts
        push!(starts, deepcopy(starts[1]))
    end
    return starts
end

# Public entry point: 1D or 2D input. seed selects the first start; the
# remaining starts follow the deterministic grid.
function fit_em_two_component(X, first_seed;
                              n_starts=default_n_starts,
                              cov_floor=cov_floor_default,
                              max_iter::Int=200,
                              tol=1e-6)
    n_starts = _positive_integer_parameter("n_starts", n_starts)
    first_seed = _nonnegative_integer_parameter("first_seed", first_seed)
    first_seed < n_starts ||
        throw(ArgumentError("first_seed must be less than n_starts"))
    tol = _positive_finite_parameter("tol", tol)
    cov_floor = _positive_finite_parameter("cov_floor", cov_floor)
    Z = X isa AbstractVector ? reshape(Float64.(X), :, 1) : Float64.(X)
    n, p = size(Z)
    K = 2
    weights = copy(EQUAL_PRIOR_WEIGHTS)
    valid_rows = [all(isfinite, view(Z, i, :)) for i in 1:n]
    Zv = Z[valid_rows, :]
    size(Zv, 1) >= K ||
        error("fit_em_two_component needs at least $K finite rows, got $(size(Zv,1))")
    starts = _init_candidates(Zv, K, n_starts)
    starts = starts[first_seed + 1:end]
    isempty(starts) && (starts = _init_candidates(Zv, K, 1))
    best = nothing
    for (means0, vars0) in starts
        means, vars, trace, converged, monotone =
            _run_em(Zv, means0, vars0, weights, cov_floor, max_iter, tol)
        ll = isempty(trace) ? -Inf : trace[end]
        cand = (means = means, vars = vars, trace = trace,
                converged = converged, monotone = monotone, ll = ll)
        if best === nothing || cand.ll > best.ll
            best = cand
        end
    end
    return TwoComponentFit(best.means, best.vars, weights, best.ll,
                           best.trace, best.converged, best.monotone)
end

function responsibilities(fit::TwoComponentFit, X)
    Z = X isa AbstractVector ? reshape(Float64.(X), :, 1) : Float64.(X)
    n = size(Z, 1)
    K = length(fit.weights)
    resp = fill(NaN, n, K)
    for i in 1:n
        all(isfinite, view(Z, i, :)) || continue
        lr = Vector{Float64}(undef, K)
        for k in 1:K
            lr[k] = log(fit.weights[k]) +
                    _log_gaussian_diag(view(Z, i, :), view(fit.means, k, :), view(fit.vars, k, :))
        end
        lse = _logsumexp_row(lr)
        if isfinite(lse)
            for k in 1:K
                resp[i, k] = exp(lr[k] - lse)
            end
        end
    end
    return resp
end

function log_likelihood_under(fit::TwoComponentFit, X)
    Z = X isa AbstractVector ? reshape(Float64.(X), :, 1) : Float64.(X)
    return _total_log_likelihood(Z, fit.means, fit.vars, fit.weights)
end

# One-component companion likelihood (single diagonal Gaussian).
struct OneComponentFit
    mean::Vector{Float64}
    var::Vector{Float64}
    loglik::Float64
end

function fit_one_component(X; cov_floor::Float64=cov_floor_default)
    Z = X isa AbstractVector ? reshape(Float64.(X), :, 1) : Float64.(X)
    n, p = size(Z)
    valid_rows = [all(isfinite, view(Z, i, :)) for i in 1:n]
    Zv = Z[valid_rows, :]
    nv = size(Zv, 1)
    mean = [nv >= 1 ? sum(view(Zv, :, j)) / nv : 0.0 for j in 1:p]
    var = [max(sum((view(Zv, :, j) .- mean[j]).^2) / max(nv, 1), cov_floor) for j in 1:p]
    weights = [1.0]
    ll = _total_log_likelihood(Zv, reshape(mean, 1, p), reshape(var, 1, p), weights)
    return OneComponentFit(mean, var, ll)
end
