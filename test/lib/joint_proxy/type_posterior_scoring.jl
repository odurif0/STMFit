function score_ncc_cost(patch::AbstractVector{<:Real}, template::AbstractVector{<:Real})
    length(patch) == length(template) || throw(ArgumentError("patch/template length mismatch"))
    n_good = 0
    dot_pt = 0.0
    norm_p2 = 0.0
    norm_t2 = 0.0
    for i in eachindex(patch)
        p = Float64(patch[i])
        t = Float64(template[i])
        if isfinite(p) && isfinite(t)
            n_good += 1
            dot_pt += p * t
            norm_p2 += p * p
            norm_t2 += t * t
        end
    end
    n_good < max(5, cld(length(patch), 2)) && return Inf
    denom = sqrt(norm_p2 * norm_t2)
    denom <= eps(Float64) && return Inf
    return -dot_pt / denom
end

log_ncc_score(patch::AbstractVector{<:Real}, template::AbstractVector{<:Real}) = -score_ncc_cost(patch, template)

effective_view_factor(rho::Real) = isnan(rho) ? 1.0 : 1.0 / (1.0 + clamp(Float64(rho), 0.0, 0.95))

function effective_view_scale(rho::Real, factor::Union{Nothing,Real}=nothing)
    if factor === nothing
        return effective_view_factor(rho)
    end
    isfinite(Float64(factor)) && Float64(factor) > 0 ? Float64(factor) :
        throw(ArgumentError("effective view factor must be positive and finite"))
end

function _view_labels(lobe::TypePosteriorLobeEvidence)
    labels = sort(collect(keys(lobe.view_patches)))
    isempty(labels) && throw(ArgumentError("lobe has no view patches"))
    return labels
end

function _validate_patch_views(lobes::AbstractVector{TypePosteriorLobeEvidence})
    isempty(lobes) && throw(ArgumentError("candidate has no lobes"))
    labels = _view_labels(lobes[1])
    for lobe in lobes[2:end]
        _view_labels(lobe) == labels || throw(ArgumentError("lobe view-label mismatch"))
    end
    return labels
end

function _parity_for_lobe(lobe::Int, n::Int, direction::Int, phase::Int)
    idx = direction == 0 ? lobe : (n - lobe + 1)
    return mod(idx - 1 + phase, 2)
end

function _template_for(entry::ProxyEntry, typ::Int, parity::Int, mirror::Int)
    found = ProxyTemplate[]
    for t in entry.templates
        if t.type == typ && t.parity == parity && t.mirror == mirror
            push!(found, t)
        end
    end
    isempty(found) && throw(ArgumentError(
        "missing ProxyTemplate for type=$typ parity=$parity mirror=$mirror in $(entry.source.name)"))
    length(found) == 1 || throw(ArgumentError(
        "duplicate ProxyTemplate for type=$typ parity=$parity mirror=$mirror in $(entry.source.name)"))
    return found[1].pixels
end

function _unary_evidence_for_state(lobes::AbstractVector{TypePosteriorLobeEvidence}, entry::ProxyEntry,
                                   direction::Int, phase::Int, mirror::Int,
                                   scale::Float64)
    labels = _validate_patch_views(lobes)
    n = length(lobes)
    unary = zeros(Float64, n, 2)
    for (i, lobe) in enumerate(lobes)
        parity = _parity_for_lobe(i, n, direction, phase)
        for typ in 0:1
            templ = _template_for(entry, typ, parity, mirror)
            total = 0.0
            for view in labels
                patch = lobe.view_patches[view]
                total += log_ncc_score(patch, templ)
            end
            unary[i, typ + 1] = scale * total
        end
    end
    all(isfinite, unary) || throw(ArgumentError("non-finite unary evidence"))
    return unary
end

function _logsumexp(v::AbstractVector{<:Real})
    m = maximum(v)
    isfinite(m) || return Float64(m)
    return m + log(sum(exp.(Float64.(v) .- m)))
end

function _normalize_logweights(logw::AbstractVector{<:Real})
    z = _logsumexp(logw)
    isfinite(z) || throw(ArgumentError("log weights are not normalizable"))
    return exp.(Float64.(logw) .- z), z
end
