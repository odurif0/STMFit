struct TypeScoringConfig
    temperature::Float64
    confidence_threshold::Float64
    config_hash::String
    manifest_hash::String
    source_hash::String
    calibrated::Bool
end

_clamp_probability(p::Real) = clamp(Float64(p), eps(Float64), 1.0 - eps(Float64))

function _logit(p::Real)
    q = _clamp_probability(p)
    return log(q / (1.0 - q))
end

function _family_scaled_probabilities(result::TypePosteriorResult, temperature::Float64)
    temperature > 0 && isfinite(temperature) || return nothing
    mats = result.lobe_marginals
    size(mats, 2) == 2 || return nothing
    all(isfinite, mats) || return nothing
    scaled = Matrix{Float64}(undef, size(mats, 1), 2)
    for i in 1:size(mats, 1)
        p1 = _clamp_probability(mats[i, 2])
        q1 = 1.0 / (1.0 + exp(-_logit(p1) / temperature))
        scaled[i, 1] = 1.0 - q1
        scaled[i, 2] = q1
    end
    return scaled
end

function _all_half_probability(family_probs::Vector{Matrix{Float64}})
    for probs in family_probs
        all(abs.(probs .- 0.5) .<= 1e-9) || return false
    end
    return true
end

function _family_votes(family_probs::Vector{Matrix{Float64}}, lobe::Int)
    votes = Int[]
    for probs in family_probs
        push!(votes, probs[lobe, 2] >= probs[lobe, 1] ? 1 : 0)
    end
    return votes
end

_lobe_family_conflict(family_probs::Vector{Matrix{Float64}}, lobe::Int) = length(unique(_family_votes(family_probs, lobe))) > 1

function _combine_family_probs(family_probs::Vector{Matrix{Float64}})
    n = size(first(family_probs), 1)
    combined = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        p1 = mean(p[i, 2] for p in family_probs)
        combined[i, 2] = p1
        combined[i, 1] = 1.0 - p1
    end
    return combined
end

function _base_rows(combined::Matrix{Float64})
    rows = TypePosteriorRow[]
    for i in 1:size(combined, 1)
        p0 = combined[i, 1]
        p1 = combined[i, 2]
        confidence = max(p0, p1)
        push!(rows, TypePosteriorRow(i, p0, p1, confidence >= 0.5 ? (p1 >= p0 ? 1 : 0) : missing,
            confidence, p1 >= p0 ? 1 : 2))
    end
    return rows
end

function _score_type_case(case::TypeCalibrationCase, cfg::TypeScoringConfig)
    cfg.temperature > 0 && isfinite(cfg.temperature) || throw(TypeCalibrationError("temperature must be finite and positive"))
    isempty(case.family_reports) && throw(TypeCalibrationError("type report missing family evidence"))

    family_probs = Matrix{Float64}[]
    for (_, result) in case.family_reports
        probs = _family_scaled_probabilities(result, cfg.temperature)
        probs === nothing && return _abstained_case(case, "ood_likelihood")
        size(probs, 1) == length(case.true_types) || return _abstained_case(case, "malformed_report")
        push!(family_probs, probs)
    end
    all(size(p, 1) == size(first(family_probs), 1) for p in family_probs) || return _abstained_case(case, "malformed_report")

    combined = _combine_family_probs(family_probs)
    rows = _base_rows(combined)
    confidences = [r.confidence for r in rows]
    flags = String[]

    if !cfg.calibrated
        push!(flags, "uncalibrated")
    end
    stale_config = cfg.config_hash != "" && case.config_hash != cfg.config_hash
    if stale_config
        push!(flags, "stale_config")
    end
    stale_hash = false
    if cfg.manifest_hash != "" && case.manifest_hash != cfg.manifest_hash
        push!(flags, "stale_hash")
        stale_hash = true
    end
    if cfg.source_hash != "" && case.source_hash != cfg.source_hash
        push!(flags, "stale_hash")
        stale_hash = true
    end
    null_control = case.control in (:null, :identical, :identical_null) || _all_half_probability(family_probs)
    if null_control
        push!(flags, "null_control")
    end

    has_family_conflict = any(i -> _lobe_family_conflict(family_probs, i), 1:length(rows))
    has_family_conflict && push!(flags, "family_conflict")

    control_conflict = case.control in (:family_conflict, :conflict)
    global_hard_abstain = !cfg.calibrated || null_control || control_conflict || stale_hash || stale_config ||
        any(f -> !isfinite(f), confidences) || any(r -> !isfinite(r.probability_0) || !isfinite(r.probability_1), rows)

    predicted = Vector{Union{Int,Missing}}(undef, length(rows))
    abstained = Vector{Bool}(undef, length(rows))
    for (i, row) in enumerate(rows)
        lobe_conflict = !control_conflict && _lobe_family_conflict(family_probs, i)
        abstain = global_hard_abstain || lobe_conflict || row.confidence < cfg.confidence_threshold
        abstained[i] = abstain
        predicted[i] = abstain ? missing : row.predicted_type
        if lobe_conflict
            push!(flags, "lobe_$(i)_family_conflict")
        elseif !global_hard_abstain && row.confidence < cfg.confidence_threshold
            push!(flags, "lobe_$(i)_low_confidence")
        end
    end

    nll = mean(-log(_clamp_probability(combined[i, case.true_types[i] + 1])) for i in eachindex(case.true_types))
    return TypeCasePosterior(case.case_id, copy(case.true_types), rows, predicted, confidences, abstained,
        nll, isempty(flags) ? "ok" : join(unique(flags), ';'))
end

function _abstained_case(case::TypeCalibrationCase, flag::String)
    rows = TypePosteriorRow[]
    predicted = Union{Int,Missing}[]
    confidences = Float64[]
    abstained = Bool[]
    for (i, _) in enumerate(case.true_types)
        push!(rows, TypePosteriorRow(i, 0.5, 0.5, missing, 0.5, 1))
        push!(predicted, missing)
        push!(confidences, 0.5)
        push!(abstained, true)
    end
    return TypeCasePosterior(case.case_id, copy(case.true_types), rows, predicted, confidences, abstained,
        NaN, flag)
end

function score_type_case(case::TypeCalibrationCase, temperature::Real; threshold::Real=0.5)
    cfg = TypeScoringConfig(Float64(temperature), Float64(threshold), "", "", "", true)
    return _score_type_case(case, cfg)
end

function score_type_case(case::TypeCalibrationCase, model::TypeCalibrationModel)
    cfg = TypeScoringConfig(model.temperature, model.confidence_threshold, model.config_hash, model.manifest_hash, model.source_hash, true)
    if case.config_hash != cfg.config_hash && cfg.config_hash != ""
        return _abstained_case(case, "stale_config")
    end
    return _score_type_case(case, cfg)
end

function predict_type(case::TypeCalibrationCase, model::Union{Nothing,TypeCalibrationModel}=nothing)
    if model === nothing
        cfg = TypeScoringConfig(1.0, 1.1, "", "", "", false)
        return _score_type_case(case, cfg)
    end
    if model.config_hash != "" && case.config_hash != model.config_hash
        return _abstained_case(case, "stale_config")
    end
    return score_type_case(case, model)
end
