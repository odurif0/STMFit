module AssignmentMixtures

using LinearAlgebra, Statistics

export MixtureFitError, initialize_state, expectation, update_state, fit_mixture, student_component_score

"A failed numerical seed is unavailable, never silently replaced by a Gaussian."
struct MixtureFitError <: Exception
    reason::String
end
Base.showerror(io::IO, e::MixtureFitError) = print(io, e.reason)

function checked_cholesky(A)
    all(isfinite, A) || throw(MixtureFitError("nonfinite_matrix"))
    try
        return cholesky(Symmetric(A))
    catch err
        err isa PosDefException || rethrow()
        throw(MixtureFitError("nonpositive_matrix"))
    end
end

"Student log mass-density, omitting only the common fixed-dimension/fixed-df constant."
function student_component_score(x, mean, scale, weight; df, guard)
    df > 0 && isfinite(df) && weight > 0 && isfinite(weight) ||
        throw(ArgumentError("positive finite Student df and mass required"))
    chol = checked_cholesky(scale + guard * I)
    distance = sum(abs2, chol.L \ (x - mean))
    return log(weight) - sum(log, diag(chol.L)) - ((df + length(x)) / 2) * log1p(distance / df)
end

function initialize_state(means, covs, weights, settings; ridge)
    p, k = size(means)
    k == 2 && length(covs) == k && length(weights) == k || throw(ArgumentError("two components required"))
    all(isfinite, means) && all(isfinite, weights) && all(>(0), weights) ||
        throw(MixtureFitError("invalid_initial_parameters"))
    loadings, noise = Matrix{Float64}[], Vector{Float64}[]
    if settings.family == "factor_analyzer"
        q = settings.factor_rank
        1 <= q < p || throw(ArgumentError("factor rank must be strictly below feature dimension"))
        for C in covs
            size(C) == (p, p) || throw(ArgumentError("invalid covariance dimensions"))
            all(isfinite, C) || throw(MixtureFitError("nonfinite_initial_covariance"))
            eig = eigen(Symmetric(C - ridge * I))
            residual = max(mean(eig.values[1:p-q]), ridge)
            L = eig.vectors[:, p-q+1:p] * Diagonal(sqrt.(max.(eig.values[p-q+1:p] .- residual, 0)))
            push!(loadings, L)
            push!(noise, fill(residual, p))
        end
        covs = [L * L' + Diagonal(psi) for (L, psi) in zip(loadings, noise)]
    elseif settings.family != "student_t"
        throw(ArgumentError("unknown alternative learning family"))
    end
    return (means=copy(means), covs=deepcopy(covs), weights=copy(weights), loadings, noise)
end

"Component kernels omit only a constant shared by both components at fixed dimension/df."
function expectation(X, state, settings)
    p, n = size(X)
    log_r, distances = zeros(n, 2), zeros(n, 2)
    for c in 1:2
        chol = checked_cholesky(state.covs[c] + settings.cholesky_guard * I)
        centered = X .- state.means[:, c]
        distances[:, c] = vec(sum(abs2, chol.L \ centered; dims=1))
        logdet_half = sum(log, diag(chol.L))
        if settings.family == "student_t"
            log_r[:, c] = log(state.weights[c]) .- logdet_half .-
                ((settings.student_df + p) / 2) .* log1p.(distances[:, c] ./ settings.student_df)
        else
            log_r[:, c] = log(state.weights[c]) .- logdet_half .- distances[:, c] ./ 2
        end
    end
    all(isfinite, log_r) || throw(MixtureFitError("nonfinite_component_scores"))
    normalizers = maximum(log_r; dims=2)
    normalizers .+= log.(sum(exp.(log_r .- normalizers); dims=2))
    responsibilities = exp.(log_r .- normalizers)
    return (responsibilities, distances, log_kernel_sum=sum(normalizers))
end

"One exact latent-variable M-step, with diagonal noise floor / ridge regularization."
function update_state(X, state, responsibilities, distances, settings; ridge)
    p, n = size(X)
    size(responsibilities) == (n, 2) && size(distances) == (n, 2) ||
        throw(ArgumentError("invalid update dimensions"))
    all(isfinite, responsibilities) && all(>=(0), responsibilities) &&
        all(s -> isapprox(s, 1; atol=1e-12, rtol=1e-12), vec(sum(responsibilities; dims=2))) ||
        throw(ArgumentError("invalid responsibilities"))
    masses = vec(sum(responsibilities; dims=1))
    all(>(settings.min_mass), masses) || throw(MixtureFitError("collapsed_component"))
    means = similar(state.means)
    covs, loadings, noise = Matrix{Float64}[], Matrix{Float64}[], Vector{Float64}[]
    for c in 1:2
        r, nk = responsibilities[:, c], masses[c]
        if settings.family == "student_t"
            u = (settings.student_df + p) ./ (settings.student_df .+ distances[:, c])
            ru = r .* u
            sum(ru) > settings.min_mass || throw(MixtureFitError("collapsed_precision_mass"))
            means[:, c] = X * ru / sum(ru)
            centered = X .- means[:, c]
            # This is a Student scale, not its nu/(nu-2) covariance. Denominator is Nk, not sum(r*u).
            push!(covs, (centered .* reshape(ru, 1, :)) * centered' / nk + ridge * I)
        else
            L, psi = state.loadings[c], state.noise[c] .+ settings.cholesky_guard
            q = size(L, 2)
            scaled = L ./ psi
            V = Matrix(checked_cholesky(Matrix{Float64}(I, q, q) + L' * scaled) \ Matrix{Float64}(I, q, q))
            Z = (V * scaled') * (X .- state.means[:, c])
            rz = Z * r
            T = [nk * V + (Z .* reshape(r, 1, :)) * Z' rz; rz' nk]
            B = hcat((X .* reshape(r, 1, :)) * Z', X * r)
            # Augment latent z by a constant: fit loadings and mean jointly (Ghahramani/Hinton Eq. 15).
            W = Matrix((checked_cholesky(T) \ B')')
            next_L, next_mu = W[:, 1:q], W[:, q+1]
            residual = vec(sum(abs2.(X) .* reshape(r, 1, :); dims=2)) / nk .-
                vec(sum(W .* B; dims=2)) / nk
            # Stored noise excludes the fixed numerical guard included in the E-step.
            next_psi = max.(residual .- settings.cholesky_guard, ridge)
            means[:, c] = next_mu
            push!(loadings, next_L); push!(noise, next_psi)
            push!(covs, next_L * next_L' + Diagonal(next_psi))
        end
    end
    all(isfinite, means) && all(C -> all(isfinite, C), covs) || throw(MixtureFitError("nonfinite_update"))
    return (means, covs, weights=masses / n, loadings, noise)
end

"Learn one fixed family, then retain it through the explicitly configured hard rule."
function fit_mixture(X, initial, settings; ridge, hard_iterations, trace=nothing)
    all(isfinite, X) || throw(ArgumentError("nonfinite training observations"))
    hard_iterations >= 0 || throw(ArgumentError("negative hard iteration count"))
    state = initialize_state(initial..., settings; ridge)
    previous = -Inf
    converged, updates = false, 0
    for iteration in 1:settings.maxiter
        e = expectation(X, state, settings)
        if isfinite(previous) && abs(e.log_kernel_sum - previous) < settings.tolerance * (1 + abs(previous))
            converged = true
            break
        end
        previous = e.log_kernel_sum
        state = update_state(X, state, e.responsibilities, e.distances, settings; ridge)
        updates += 1
        trace === nothing || push!(trace, (stage="em", iteration, state=deepcopy(state)))
    end
    for iteration in 1:hard_iterations
        e = expectation(X, state, settings)
        assignments = settings.hard_assignment == "student_density" ?
            [argmax(view(e.responsibilities, i, :)) for i in axes(X, 2)] :
            [argmin(view(e.distances, i, :)) for i in axes(X, 2)]
        r = Float64.([assignments[i] == c for i in axes(X, 2), c in 1:2])
        state = update_state(X, state, r, e.distances, settings; ridge)
        trace === nothing || push!(trace, (stage="hard", iteration, state=deepcopy(state)))
    end
    e = expectation(X, state, settings)
    return (state, updates, converged, log_kernel_sum=e.log_kernel_sum,
        precision_range=settings.family == "student_t" ?
            extrema((settings.student_df + size(X, 1)) ./ (settings.student_df .+ e.distances)) : (NaN, NaN))
end

end
