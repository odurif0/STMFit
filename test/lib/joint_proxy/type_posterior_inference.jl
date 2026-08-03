struct BinaryChainPosterior
    log_evidence::Float64
    state_marginals::Matrix{Float64}
    log_state_marginals::Matrix{Float64}
    forward::Matrix{Float64}
    backward::Matrix{Float64}
end

function binary_chain_forward_backward(unary::AbstractMatrix{<:Real};
                                       log_transition::AbstractMatrix{<:Real}=zeros(2, 2))
    size(unary, 2) == 2 || throw(ArgumentError("binary_chain_forward_backward expects 2 states"))
    size(log_transition) == (2, 2) || throw(ArgumentError("transition matrix must be 2x2"))
    all(isfinite, unary) || throw(ArgumentError("non-finite unary evidence"))
    all(isfinite, log_transition) || throw(ArgumentError("non-finite transition evidence"))
    n = size(unary, 1)
    fwd = fill(-Inf, n, 2)
    bwd = fill(-Inf, n, 2)
    fwd[1, :] .= unary[1, :]
    for i in 2:n, s in 1:2
        fwd[i, s] = unary[i, s] + _logsumexp(fwd[i - 1, :] .+ log_transition[:, s])
    end
    bwd[n, :] .= 0.0
    for i in (n - 1):-1:1, s in 1:2
        bwd[i, s] = _logsumexp(log_transition[s, :] .+ unary[i + 1, :] .+ bwd[i + 1, :])
    end
    log_evidence = _logsumexp(@view fwd[n, :])
    log_post = fwd .+ bwd .- log_evidence
    marginals = exp.(log_post)
    return BinaryChainPosterior(log_evidence, marginals, log_post, fwd, bwd)
end

function binary_chain_viterbi(unary::AbstractMatrix{<:Real};
                              log_transition::AbstractMatrix{<:Real}=zeros(2, 2))
    size(unary, 2) == 2 || throw(ArgumentError("binary_chain_viterbi expects 2 states"))
    size(log_transition) == (2, 2) || throw(ArgumentError("transition matrix must be 2x2"))
    all(isfinite, unary) || throw(ArgumentError("non-finite unary evidence"))
    all(isfinite, log_transition) || throw(ArgumentError("non-finite transition evidence"))
    n = size(unary, 1)
    dp = fill(-Inf, n, 2)
    prev = fill(0, n, 2)
    dp[1, :] .= unary[1, :]
    for i in 2:n, s in 1:2
        vals = dp[i - 1, :] .+ log_transition[:, s]
        best = argmax(vals)
        dp[i, s] = unary[i, s] + vals[best]
        prev[i, s] = best
    end
    last = argmax(@view dp[n, :])
    seq = fill(0, n)
    seq[n] = last - 1
    state = last
    for i in (n - 1):-1:1
        state = prev[i + 1, state]
        seq[i] = state - 1
    end
    return seq, dp[n, last]
end

function _family_sensitivity(states::Vector{TypePosteriorGlobalState})
    by_family = Dict{String,Float64}()
    for s in states
        by_family[s.family] = get(by_family, s.family, 0.0) + s.posterior
    end
    return sort(collect(by_family); by=first)
end

function infer_type_posterior(lobes::AbstractVector{TypePosteriorLobeEvidence}, ensemble::ProxyEnsemble;
                              priors::TypePosteriorPriors=TypePosteriorPriors(),
                              rho::Real=NaN,
                              effective_factor::Union{Nothing,Real}=nothing,
                              pair_bond_evidence=nothing)
    pair_bond_evidence === nothing || throw(ArgumentError(
        "concatenated unary bond evidence is rejected in v1"))
    scale = effective_view_scale(rho, effective_factor)
    _validate_patch_views(lobes)
    pri = _normalized_priors(priors, length(ensemble.entries))
    type_log = reshape(log.(pri.type), 1, 2)
    state_logw = Float64[]
    states = TypePosteriorGlobalState[]
    state_posts = Vector{Matrix{Float64}}()
    best_logw = -Inf
    best_seq = Int[]
    for dir in 0:1, phase in 0:1, mirror in 0:1, (sidx, entry) in enumerate(ensemble.entries)
        unary = _unary_evidence_for_state(lobes, entry, dir, phase, mirror, scale)
        unary .+= type_log
        fb = binary_chain_forward_backward(unary)
        seq, vll = binary_chain_viterbi(unary)
        srcw = pri.source === nothing ? entry.source.weight : pri.source[sidx]
        lp = log(pri.direction[dir + 1]) + log(pri.phase[phase + 1]) +
             log(pri.mirror[mirror + 1]) + log(srcw)
        lw = lp + fb.log_evidence
        push!(state_logw, lw)
        push!(state_posts, fb.state_marginals)
        push!(states, TypePosteriorGlobalState(dir, phase, mirror, sidx, entry.source.name,
                                               entry.source.family, lp, fb.log_evidence, 0.0))
        if isempty(best_seq) || lw > best_logw
            best_logw = lw
            best_seq = seq
        end
        _ = vll
    end
    state_probs, log_evidence = _normalize_logweights(state_logw)
    states = [TypePosteriorGlobalState(s.direction, s.phase, s.mirror, s.source_index,
                                       s.source_name, s.family, s.log_prior, s.log_evidence,
                                       state_probs[i]) for (i, s) in enumerate(states)]
    lobe_marginals = zeros(Float64, length(lobes), 2)
    for (p, marg) in zip(state_probs, state_posts)
        lobe_marginals .+= p .* marg
    end
    best_seq = isempty(best_seq) ? fill(0, length(lobes)) : best_seq
    return TypePosteriorResult(lobe_marginals, best_seq, log_evidence, states,
                               _family_sensitivity(states))
end
