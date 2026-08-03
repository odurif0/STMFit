function _case_nll(case::TypeCalibrationCase, temperature::Float64)
    scored = _score_type_case(case, TypeScoringConfig(temperature, 0.0, "", "", "", true))
    return scored.nll
end

function _grid_search_temperature(cases::Vector{TypeCalibrationCase}, temps::Vector{Float64})
    nlls = Float64[]
    for t in temps
        push!(nlls, mean(_case_nll(c, t) for c in cases))
    end
    best_idx = argmin(nlls)
    return temps, nlls, best_idx
end

function _row_samples(cases::Vector{TypeCalibrationCase}, scored::Vector{TypeCasePosterior})
    samples = NamedTuple[]
    for (case, score) in zip(cases, scored)
        for row in score.rows
            row.predicted_type === missing && continue
            true_type = case.true_types[row.lobe]
            p_true = true_type == 0 ? row.probability_0 : row.probability_1
            push!(samples, (confidence=row.confidence, correct=row.predicted_type == true_type,
                nll=-log(clamp(p_true, eps(Float64), 1.0 - eps(Float64)))))
        end
    end
    return samples
end

function _threshold_and_bins(cases::Vector{TypeCalibrationCase}, scored::Vector{TypeCasePosterior}; bin_count::Int=5)
    samples = _row_samples(cases, scored)
    if isempty(samples)
        return 1.0, "all_abstained", TypeReliabilityBin[]
    end

    wrong_confs = [s.confidence for s in samples if !s.correct]
    right_confs = [s.confidence for s in samples if s.correct]
    null_confs = [maximum(score.confidences) for (case, score) in zip(cases, scored) if case.control in (:null, :identical, :identical_null)]

    source = "heldout_correct_min"
    threshold = isempty(right_confs) ? 1.0 : minimum(right_confs)
    if !isempty(wrong_confs) || !isempty(null_confs)
        threshold = min(1.0, max(maximum(vcat(wrong_confs, null_confs...)), threshold) + 1e-12)
        source = isempty(null_confs) ? "heldout_error_max" : "heldout_error_null_max"
    end

    edges = collect(range(0.0, 1.0; length=bin_count + 1))
    bins = TypeReliabilityBin[]
    for i in 1:bin_count
        lo, hi = edges[i], edges[i + 1]
        members = [s for s in samples if (i < bin_count ? (s.confidence >= lo && s.confidence < hi) : (s.confidence >= lo && s.confidence <= hi))]
        isempty(members) && push!(bins, TypeReliabilityBin(lo, hi, 0, 0, NaN, NaN))
        if !isempty(members)
            push!(bins, TypeReliabilityBin(lo, hi, length(members), count(s -> s.correct, members),
                mean(s.confidence for s in members), mean(s.nll for s in members)))
        end
    end
    return threshold, source, bins
end

function calibrate_type_temperature(cases::AbstractVector{TypeCalibrationCase};
                                    seed::Integer, config_hash::AbstractString,
                                    manifest_hash::AbstractString, source_hash::AbstractString,
                                    grid_min::Float64=0.25, grid_max::Float64=8.0,
                                    grid_points::Int=17, refine_rounds::Int=3,
                                    bin_count::Int=5)
    cal, held = _validate_calibration_cases(cases)
    _ensure_same_seed_config_manifest_source(cases, seed, config_hash, manifest_hash, source_hash)
    grid_min > 0 && grid_max > grid_min || throw(TypeCalibrationError("temperature grid must be positive and ordered"))
    grid_points >= 3 || throw(TypeCalibrationError("temperature grid must have at least 3 points"))

    lo, hi = grid_min, grid_max
    temps = Float64[]
    nlls = Float64[]
    best_idx = 1
    for _ in 1:max(refine_rounds, 1)
        temps = collect(exp.(range(log(lo), log(hi); length=grid_points)))
        temps, nlls, best_idx = _grid_search_temperature(cal, temps)
        best_temp = temps[best_idx]
        span = sqrt(3.0)
        lo = max(best_temp / span, eps(Float64))
        hi = best_temp * span
    end

    temperature = temps[best_idx]
    calibration_nll = nlls[best_idx]
    baseline_nll = mean(_case_nll(c, 1.0) for c in held)
    held_scored = [_score_type_case(c, TypeScoringConfig(temperature, 0.0, config_hash, manifest_hash, source_hash, true)) for c in held]
    heldout_nll = mean(s.nll for s in held_scored)
    threshold, source, bins = _threshold_and_bins(held, held_scored; bin_count=bin_count)
    diag = TypeCalibrationDiagnostics(copy(temps), copy(nlls), calibration_nll, heldout_nll,
        baseline_nll, best_idx, source, length(held), length(cal))
    return TypeCalibrationModel(temperature, threshold, calibration_nll, heldout_nll,
        baseline_nll, diag, bins, config_hash, manifest_hash, source_hash)
end
