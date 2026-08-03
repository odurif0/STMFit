struct TypePosteriorPriors
    type::Vector{Float64}
    direction::Vector{Float64}
    phase::Vector{Float64}
    mirror::Vector{Float64}
    source::Union{Nothing,Vector{Float64}}
end

TypePosteriorPriors(; type=[1.0, 1.0], direction=[1.0, 1.0], phase=[1.0, 1.0],
                    mirror=[1.0, 1.0], source=nothing) =
    TypePosteriorPriors(Float64.(type), Float64.(direction), Float64.(phase),
                         Float64.(mirror), source === nothing ? nothing : Float64.(source))

struct TypePosteriorGlobalState
    direction::Int
    phase::Int
    mirror::Int
    source_index::Int
    source_name::String
    family::String
    log_prior::Float64
    log_evidence::Float64
    posterior::Float64
end

struct TypePosteriorLobeEvidence
    view_patches::Dict{String,Vector{Float64}}
end

struct TypePosteriorResult
    lobe_marginals::Matrix{Float64}
    map_sequence::Vector{Int}
    log_evidence::Float64
    global_state_posterior::Vector{TypePosteriorGlobalState}
    proxy_family_sensitivity::Vector{Pair{String,Float64}}
end

function _reject_bad_weights(label::String, weights::AbstractVector{<:Real})
    isempty(weights) && throw(ArgumentError("$label prior vector is empty"))
    all(isfinite, weights) || throw(ArgumentError("$label prior contains non-finite values"))
    all(>(0), weights) || throw(ArgumentError("$label prior must be strictly positive"))
    s = sum(weights)
    isfinite(s) && s > 0 || throw(ArgumentError("$label prior sum must be positive and finite"))
    return Float64.(weights) ./ s
end

function _normalized_priors(priors::TypePosteriorPriors, n_sources::Int)
    type = _reject_bad_weights("type", priors.type)
    length(type) == 2 || throw(ArgumentError("type prior must have length 2"))
    direction = _reject_bad_weights("direction", priors.direction)
    phase = _reject_bad_weights("phase", priors.phase)
    mirror = _reject_bad_weights("mirror", priors.mirror)
    source = priors.source === nothing ? nothing : _reject_bad_weights("source", priors.source)
    if source !== nothing && length(source) != n_sources
        throw(ArgumentError("source prior length mismatch"))
    end
    return (; type, direction, phase, mirror, source)
end

function _logweights(p::AbstractVector{<:Real})
    all(isfinite, p) || throw(ArgumentError("posterior weights contain non-finite values"))
    s = sum(p)
    isfinite(s) && s > 0 || throw(ArgumentError("posterior weights must sum to a positive finite value"))
    return log.(Float64.(p) ./ s)
end

function _type_prior_logrow(priors::TypePosteriorPriors)
    type = _reject_bad_weights("type", priors.type)
    return log.(type)'
end
