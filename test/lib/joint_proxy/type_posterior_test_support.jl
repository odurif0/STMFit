_basis(i, n) = [j == i ? 1.0 : 0.0 for j in 1:n]
_logsumexp(v) = (m = maximum(v); isfinite(m) ? m + log(sum(exp.(v .- m))) : m)

function _template_vec(type::Int, parity::Int, mirror::Int; source_bias::Float64=0.0)
    v = _basis(1 + type + 2 * parity + 4 * mirror, 8)
    source_bias != 0.0 && (v[8] += source_bias)
    return v
end

function _ensemble(; swap::Bool=false)
    entries = ProxyEntry[]
    for (family, name, bias, weight) in (("geom", "geom", 0.0, 0.6), ("stm", "stm_h050", 0.15, 0.4))
        tmpls = ProxyTemplate[]
        for typ in 0:1, parity in 0:1, mirror in 0:1
            push!(tmpls, ProxyTemplate(typ, parity, mirror,
                                       _template_vec(swap ? 1 - typ : typ, parity, mirror; source_bias=bias)))
        end
        push!(entries, ProxyEntry(ProxySource(family, name, "$(name).tsv", "sha_$name", 0.5, weight, true), tmpls))
    end
    return ProxyEnsemble(entries, 0.32, 0.08, 9, 81, "payload", String[])
end

function _geometric_templates()
    path = joinpath(@__DIR__, "..", "..", "..", "templates", "chitosan_geometric_sites.tsv")
    return render_geometric_templates(path; half_nm=0.32, step_nm=0.08,
                                      parity_flip="t", mirror_flip="u", normalize="zscore")
end

function _actual_ensemble(; shuffle::Bool=false, swap_metadata::Bool=false)
    tmpls = _geometric_templates()
    swap_metadata && (tmpls = [ProxyTemplate(1 - t.type, t.parity, t.mirror, copy(t.pixels)) for t in tmpls])
    shuffle && (tmpls = tmpls[[8, 1, 6, 3, 4, 2, 7, 5]])
    return ProxyEnsemble([ProxyEntry(ProxySource("geometric", "geometric", "templates/chitosan_geometric_sites.tsv",
                                                 "sha_geometric", 0.0, 1.0, true), tmpls)],
                         0.32, 0.08, 9, 81, "payload", String[])
end

function _patch_vec(true_t::Vector{Float64}, false_t::Vector{Float64}, α::Float64)
    v = true_t .+ α .* false_t
    return v ./ norm(v)
end

function _evidence(n::Int; duplicate::Bool=false, reversed::Bool=false)
    true_types = isodd(n) ? [mod1(i, 2) - 1 for i in 1:n] : [0 for _ in 1:n]
    reversed && (true_types = reverse(true_types))
    lobes = TypePosteriorLobeEvidence[]
    for (i, typ) in enumerate(true_types)
        parity = mod(i - 1, 2)
        true_t = _template_vec(typ, parity, 0)
        false_t = _template_vec(1 - typ, parity, 0)
        fwd = _patch_vec(true_t, false_t, 0.25 + 0.02 * i)
        bwd = duplicate ? copy(fwd) : _patch_vec(true_t, false_t, 0.10 + 0.01 * i)
        push!(lobes, TypePosteriorLobeEvidence(Dict("fwd" => fwd, "bwd" => bwd)))
    end
    return lobes, true_types, duplicate ? 1.0 : 0.0
end

function _symmetric_evidence(n::Int; reversed::Bool=false)
    true_types = [mod1(i, 2) - 1 for i in 1:n]
    reversed && (true_types = reverse(true_types))
    lobes = TypePosteriorLobeEvidence[]
    for (i, typ) in enumerate(true_types)
        parity = mod(i - 1, 2)
        true_t = _template_vec(typ, parity, 0)
        false_t = _template_vec(1 - typ, parity, 0)
        push!(lobes, TypePosteriorLobeEvidence(Dict("fwd" => _patch_vec(true_t, false_t, 0.2),
                                                    "bwd" => _patch_vec(true_t, false_t, 0.2))))
    end
    return lobes, 0.0
end

function _geometric_evidence(n::Int; duplicate::Bool=false, reversed::Bool=false)
    entry = _actual_ensemble().entries[1]
    true_types = [mod1(i, 2) - 1 for i in 1:n]
    reversed && (true_types = reverse(true_types))
    lobes = TypePosteriorLobeEvidence[]
    for (i, typ) in enumerate(true_types)
        parity = mod(i - 1, 2)
        true_t = JointProxyTypePosterior._template_for(entry, typ, parity, 0)
        false_t = JointProxyTypePosterior._template_for(entry, 1 - typ, parity, 0)
        fwd = _patch_vec(true_t, false_t, 0.25 + 0.01 * i)
        bwd = duplicate ? copy(fwd) : _patch_vec(true_t, false_t, 0.08 + 0.01 * i)
        push!(lobes, TypePosteriorLobeEvidence(Dict("fwd" => fwd, "bwd" => bwd)))
    end
    return lobes, duplicate ? 1.0 : 0.0
end

function _bruteforce(lobes::Vector{TypePosteriorLobeEvidence}, ensemble::ProxyEnsemble,
                     priors::TypePosteriorPriors; rho::Real=NaN, factor::Union{Nothing,Real}=nothing)
    n = length(lobes)
    weights = Float64[]
    seqs = Vector{Int}[]
    states = NamedTuple[]
    type_prior = priors.type ./ sum(priors.type)
    dir_prior = priors.direction ./ sum(priors.direction)
    phase_prior = priors.phase ./ sum(priors.phase)
    mirror_prior = priors.mirror ./ sum(priors.mirror)
    source_prior = priors.source === nothing ? [e.source.weight for e in ensemble.entries] : priors.source
    source_prior = source_prior ./ sum(source_prior)
    scale = isnothing(factor) ? effective_view_factor(rho) : Float64(factor)
    for dir in 0:1, phase in 0:1, mirror in 0:1, (sidx, entry) in enumerate(ensemble.entries)
        state_log = log(dir_prior[dir + 1]) + log(phase_prior[phase + 1]) +
                    log(mirror_prior[mirror + 1]) + log(source_prior[sidx])
        unary = zeros(Float64, n, 2)
        for (i, lobe) in enumerate(lobes)
            parity = mod((dir == 0 ? i : (n - i + 1)) - 1 + phase, 2)
            for typ in 0:1
                templ = JointProxyTypePosterior._template_for(entry, typ, parity, mirror)
                score = sum(-score_ncc_cost(lobe.view_patches[v], templ) for v in sort(collect(keys(lobe.view_patches))))
                unary[i, typ + 1] = scale * score + log(type_prior[typ + 1])
            end
        end
        fb = binary_chain_forward_backward(unary)
        vseq, vlog = binary_chain_viterbi(unary)
        for seq in Iterators.product(ntuple(_ -> 0:1, n)...)
            ll = state_log
            for (i, typ) in enumerate(seq)
                ll += unary[i, typ + 1]
            end
            push!(weights, ll)
            push!(seqs, collect(seq))
            push!(states, (fb=fb, vseq=vseq, vlog=vlog))
        end
    end
    logz = _logsumexp(weights)
    probs = exp.(weights .- logz)
    return logz, probs, seqs, states
end
