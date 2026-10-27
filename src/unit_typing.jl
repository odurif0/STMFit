# Label-free typing of units in two families from their repeated measurements
# (campaign T, docs/src/journal.md 2026-10-27): the reproducible subspace of unit
# means over two halves of the scans, Gaussian mixtures with full covariance
# compared by BIC, and the agreement of the types between the halves.
# No label, position or composition prior enters.
module UnitTyping

using LinearAlgebra, Random, Statistics

export Mixture, fit_mixture, mixture_bic, posterior, reproducible_subspace, cohen_kappa, robust_z,
       psd_shrink, contrastive_direction, winsorize,
       fit_two_types, fit_continuum, type_posterior

"Gaussian mixture with full covariances: weights, means (k × p) and covariances."
struct Mixture
    weights::Vector{Float64}
    means::Matrix{Float64}
    covs::Vector{Matrix{Float64}}
    loglik::Float64
end

function _logdens(X, m, S)
    p = size(X, 2); C = cholesky(Symmetric(S))
    D = (X .- m') / C.U
    return -0.5 .* (vec(sum(abs2, D; dims=2)) .+ 2sum(log, diag(C.U)) .+ p * log(2π))
end

function _logsumexp_rows(L)
    m = maximum(L; dims=2)
    return vec(m) .+ log.(vec(sum(exp.(L .- m); dims=2)))
end

"""
    fit_mixture(X, k; starts=20, iters=500, tol=1e-9, seed=20261027, ridge=1e-6) -> Mixture

Maximum-likelihood Gaussian mixture of `k` full-covariance components on the
rows of `X` (expectation–maximisation from `starts` random starts, best
likelihood kept; `ridge` is added to every covariance).
"""
function fit_mixture(X::AbstractMatrix, k::Integer; starts::Integer=20, iters::Integer=500, tol::Real=1e-9,
                     seed::Integer=20261027, ridge::Real=1e-6)
    n, p = size(X); n > k || throw(ArgumentError("need more rows than components"))
    rng = MersenneTwister(seed); best = nothing
    S0 = cov(X) + ridge * I
    for _ in 1:(k == 1 ? 1 : starts)
        m = X[randperm(rng, n)[1:k], :]; S = [Matrix(S0) for _ in 1:k]; w = fill(1 / k, k); ll = -Inf
        for _ in 1:iters
            L = hcat([log(w[j]) .+ _logdens(X, m[j, :], S[j]) for j in 1:k]...)
            lse = _logsumexp_rows(L); newll = sum(lse)
            G = exp.(L .- lse)
            w = vec(mean(G; dims=1)); w = max.(w, 1e-12); w ./= sum(w)
            for j in 1:k
                g = G[:, j]; sg = sum(g)
                sg > 1e-9 || continue
                m[j, :] = vec(sum(g .* X; dims=1)) ./ sg
                D = X .- m[j, :]'
                S[j] = (D' * (g .* D)) ./ sg + ridge * I
            end
            abs(newll - ll) < tol * max(1.0, abs(newll)) && (ll = newll; break)
            ll = newll
        end
        (best === nothing || ll > best.loglik) && (best = Mixture(copy(w), copy(m), deepcopy(S), ll))
    end
    return best
end

"BIC of a mixture fitted on `n` rows of dimension `p`."
function mixture_bic(M::Mixture, n::Integer)
    k = length(M.weights); p = size(M.means, 2)
    params = k * p + k * p * (p + 1) ÷ 2 + (k - 1)
    return -2M.loglik + params * log(n)
end

"Posterior probabilities (rows of `X` × components) under a mixture."
function posterior(M::Mixture, X::AbstractMatrix)
    L = hcat([log(M.weights[j]) .+ _logdens(X, M.means[j, :], M.covs[j]) for j in eachindex(M.weights)]...)
    return exp.(L .- _logsumexp_rows(L))
end

"""
    reproducible_subspace(A, B; share=0.8, kmax=4) -> (V, eigenvalues)

Leading eigenvectors (columns of V) of the symmetrised cross-covariance of two
matrices of the same units measured twice (rows = units, columns = features,
already centred and scaled): the directions of variation that repeat. Kept up
to `share` of the positive eigenvalues, at most `kmax`, at least one.
"""
function reproducible_subspace(A::AbstractMatrix, B::AbstractMatrix; share::Real=0.8, kmax::Integer=4)
    size(A) == size(B) || throw(ArgumentError("the two measurements must have the same shape"))
    C = (A' * B + B' * A) ./ (2size(A, 1))
    e = eigen(Symmetric(C)); vals = reverse(e.values); vecs = e.vectors[:, end:-1:1]
    pos = sum(v for v in vals if v > 0; init=0.0)
    k = 1
    while k < min(kmax, length(vals)) && vals[k+1] > 0 && sum(vals[1:k]) < share * pos
        k += 1
    end
    return vecs[:, 1:k], vals
end

"Symmetric matrix made positive semi-definite (negative eigenvalues set to 0) and shrunk by `shrink` × its mean eigenvalue on the diagonal."
function psd_shrink(C::AbstractMatrix; shrink::Real=0.1)
    e = eigen(Symmetric(Matrix(C))); λ = max.(e.values, 0.0)
    P = e.vectors * Diagonal(λ) * e.vectors'
    return Symmetric(P + shrink * max(mean(λ), eps()) * I)
end

"""
    contrastive_direction(Ca, Cc; shrink=0.1) -> (v, λ)

Leading direction of excess variation of `Ca` over `Cc` (largest λ of
Ca v = λ Cc v, with `Cc` made positive semi-definite and shrunk), scaled to
unit length.
"""
function contrastive_direction(Ca::AbstractMatrix, Cc::AbstractMatrix; shrink::Real=0.1)
    B = psd_shrink(Cc; shrink)
    e = eigen(Symmetric(Matrix((Ca + Ca') / 2)), B)
    i = argmax(e.values); v = e.vectors[:, i]
    return v ./ norm(v), e.values[i]
end

"Values clamped to the `lo` and `hi` quantiles of the finite ones."
function winsorize(v::AbstractVector{<:Real}, lo::Real=0.025, hi::Real=0.975)
    ok = filter(isfinite, v); a, b = quantile(ok, lo), quantile(ok, hi)
    return clamp.(v, a, b)
end

_normlog(x, m, v) = -0.5 * (log(2π * v) + (x - m)^2 / v)

"""
    fit_two_types(s, v; starts=20, iters=1000, seed=20261027) -> (π, μ1, μ2, loglik)

Two Gaussian types on scores `s`, every score with its own known variance
`v[i]` (the same for both types: type width plus measurement noise), weight π
of type 1; expectation–maximisation from `starts` starts.
"""
function fit_two_types(s::AbstractVector, v::AbstractVector; starts::Integer=20, iters::Integer=1000, seed::Integer=20261027)
    rng = MersenneTwister(seed); best = (π=NaN, μ1=NaN, μ2=NaN, loglik=-Inf)
    for _ in 1:starts
        μ1, μ2 = quantile(s, rand(rng) * 0.4 + 0.05), quantile(s, 0.95 - rand(rng) * 0.4); π = 0.5; ll = -Inf
        for _ in 1:iters
            a = [log(π) + _normlog(s[i], μ1, v[i]) for i in eachindex(s)]
            b = [log(1 - π) + _normlog(s[i], μ2, v[i]) for i in eachindex(s)]
            m = max.(a, b); lse = m .+ log.(exp.(a .- m) .+ exp.(b .- m)); newll = sum(lse)
            r = exp.(a .- lse)
            π = clamp(mean(r), 1e-6, 1 - 1e-6)
            w1 = r ./ v; w2 = (1 .- r) ./ v
            μ1 = sum(w1 .* s) / sum(w1); μ2 = sum(w2 .* s) / sum(w2)
            abs(newll - ll) < 1e-10 * max(1.0, abs(newll)) && (ll = newll; break)
            ll = newll
        end
        ll > best.loglik && (best = (π=π, μ1=μ1, μ2=μ2, loglik=ll))
    end
    return best
end

"Posterior probability of type 1 for scores `s` with variances `v` under `fit_two_types` parameters."
function type_posterior(s::AbstractVector, v::AbstractVector, π, μ1, μ2)
    return [1 / (1 + exp(log(1 - π) + _normlog(s[i], μ2, v[i]) - log(π) - _normlog(s[i], μ1, v[i]))) for i in eachindex(s)]
end

"""
    fit_continuum(s, noise) -> (μ, τ2, loglik)

One Gaussian on scores `s` whose variance is τ2 + noise[i], τ2 ≥ 0 free
(golden-section search on τ2; μ the weighted mean for each τ2).
"""
function fit_continuum(s::AbstractVector, noise::AbstractVector)
    function ll(t2)
        v = t2 .+ noise; w = 1 ./ v; μ = sum(w .* s) / sum(w)
        return sum(_normlog(s[i], μ, v[i]) for i in eachindex(s)), μ
    end
    lo, hi = 0.0, 10var(s) + maximum(noise); g = (sqrt(5) - 1) / 2
    a, b = hi - g * (hi - lo), lo + g * (hi - lo)
    for _ in 1:200
        ll(a)[1] > ll(b)[1] ? (hi = b) : (lo = a)
        a, b = hi - g * (hi - lo), lo + g * (hi - lo)
    end
    t2 = (lo + hi) / 2; L, μ = ll(t2)
    L0, μ0 = ll(0.0)
    return L0 > L ? (μ=μ0, τ2=0.0, loglik=L0) : (μ=μ, τ2=t2, loglik=L)
end

"Cohen's κ of two Boolean vectors."
function cohen_kappa(a::AbstractVector{Bool}, b::AbstractVector{Bool})
    length(a) == length(b) && !isempty(a) || throw(ArgumentError("two non-empty vectors of the same length"))
    po = mean(a .== b); pa, pb = mean(a), mean(b); pe = pa * pb + (1 - pa) * (1 - pb)
    return pe ≈ 1 ? NaN : (po - pe) / (1 - pe)
end

"Robust z of the finite values of `v` (minus median, over 1.4826 × MAD; scale 1 when the MAD is 0); NaN stays NaN."
function robust_z(v::AbstractVector{<:Real})
    ok = filter(isfinite, v)
    isempty(ok) && return fill(NaN, length(v))
    med = median(ok); s = 1.4826 * median(abs.(ok .- med))
    s > 0 || (s = 1.0)
    return [isfinite(x) ? (x - med) / s : NaN for x in v]
end

end
