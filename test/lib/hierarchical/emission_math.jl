# Numerical helpers shared by the two-component and one-component emissions.

const LOG2PI = log(2π)

# Log density of a diagonal Gaussian at row x_row given mean/var vectors.
function _log_gaussian_diag(x_row::AbstractVector, mean::AbstractVector,
                            var::AbstractVector)
    p = length(x_row)
    ll = 0.0
    @inbounds for j in 1:p
        v = max(Float64(var[j]), 1e-300)
        d = Float64(x_row[j]) - Float64(mean[j])
        ll += -0.5 * (LOG2PI + log(v) + d * d / v)
    end
    return ll
end

# Row log-likelihoods under each component: returns an n × K matrix of
# log(weight_k * N(x_i | mu_k, var_k)).
function _log_resp_matrix(Z::AbstractMatrix, means::AbstractMatrix,
                          vars::AbstractMatrix, weights::Vector{Float64})
    n, p = size(Z)
    K = length(weights)
    lr = zeros(n, K)
    @inbounds for i in 1:n, k in 1:K
        lr[i, k] = log(weights[k]) +
                   _log_gaussian_diag(view(Z, i, :), view(means, k, :), view(vars, k, :))
    end
    return lr
end

function _logsumexp_row(v::AbstractVector)
    m = maximum(v)
    isfinite(m) || return -Inf
    return m + log(sum(x -> exp(x - m), v))
end

function _normalize_log_resp!(lr::Matrix{Float64})
    n, K = size(lr)
    for i in 1:n
        lse = _logsumexp_row(view(lr, i, :))
        if isfinite(lse)
            for k in 1:K
                lr[i, k] = lr[i, k] - lse
            end
        else
            for k in 1:K
                lr[i, k] = log(1.0 / K)
            end
        end
    end
    return lr
end

function _total_log_likelihood(Z::AbstractMatrix, means::AbstractMatrix,
                               vars::AbstractMatrix, weights::Vector{Float64})
    lr = _log_resp_matrix(Z, means, vars, weights)
    ll = 0.0
    for i in 1:size(lr, 1)
        ll += _logsumexp_row(view(lr, i, :))
    end
    return ll
end
