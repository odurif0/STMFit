function _case_nll(case::CountCalibrationCase, temperature::Float64)
    scored = score_count_report(case.report, temperature; threshold=0.0, true_n=case.true_n)
    return scored.nll
end

function _threshold_and_bins(scored::Vector{CountCasePosterior}; bin_count::Int=5)
    emitted = [s for s in scored if !s.abstained]
    confs = [s.confidence for s in emitted]
    correct = [s.predicted_n == s.true_n for s in emitted]
    if isempty(emitted)
        return 1.0, "all_abstained", CountReliabilityBin[]
    end
    wrong_confs = [s.confidence for s in emitted if s.predicted_n != s.true_n]
    right_confs = [s.confidence for s in emitted if s.predicted_n == s.true_n]
    if isempty(wrong_confs)
        threshold = isempty(right_confs) ? 1.0 : minimum(right_confs)
        source = "heldout_correct_min"
    else
        threshold = min(1.0, maximum(wrong_confs) + 1e-12)
        source = "heldout_error_max"
    end
    edges = collect(range(0.0, 1.0; length=bin_count + 1))
    bins = CountReliabilityBin[]
    for i in 1:bin_count
        lo, hi = edges[i], edges[i + 1]
        members = [s for s in scored if (i < bin_count ? (s.confidence >= lo && s.confidence < hi) : (s.confidence >= lo && s.confidence <= hi))]
        isempty(members) && push!(bins, CountReliabilityBin(lo, hi, 0, 0, NaN, NaN))
        if !isempty(members)
            push!(bins, CountReliabilityBin(lo, hi, length(members), count(s -> s.predicted_n == s.true_n, members),
                mean(s.confidence for s in members), mean(s.nll for s in members)))
        end
    end
    return threshold, source, bins
end

function _grid_search_temperature(cases::Vector{CountCalibrationCase}, temps::Vector{Float64})
    nlls = Float64[]
    for t in temps
        push!(nlls, mean(_case_nll(c, t) for c in cases))
    end
    best_idx = argmin(nlls)
    return temps, nlls, best_idx
end

function calibrate_count_temperature(cases::AbstractVector{CountCalibrationCase};
                                     seed::Integer, config_hash::AbstractString,
                                     grid_min::Float64=0.25, grid_max::Float64=8.0,
                                     grid_points::Int=17, refine_rounds::Int=3,
                                     bin_count::Int=5)
    isempty(cases) && throw(CountCalibrationError("zero cases: synthetic calibration requires at least one case"))
    _ensure_same_seed_and_config(cases, seed, config_hash)
    cal = [c for c in cases if c.split == :calibration]
    held = [c for c in cases if c.split == :heldout]
    isempty(cal) && throw(CountCalibrationError("zero cases: no calibration split provided"))
    isempty(held) && throw(CountCalibrationError("zero cases: no heldout split provided"))
    any(length(_finite_candidates(c.report)) < 2 for c in cal) && throw(CountCalibrationError(
        "insufficient one-candidate calibration: at least one calibration case has only one valid candidate"))

    lo, hi = grid_min, grid_max
    temps = Float64[]
    nlls = Float64[]
    best_idx = 1
    for _ in 1:max(refine_rounds, 1)
        lo > 0 && hi > lo || throw(CountCalibrationError("temperature grid must be positive and ordered"))
        temps = collect(exp.(range(log(lo), log(hi); length=grid_points)))
        temps, nlls, best_idx = _grid_search_temperature(cal, temps)
        best_temp = temps[best_idx]
        span = sqrt(3.0)
        lo = max(best_temp / span, eps(Float64))
        hi = best_temp * span
    end
    temperature = temps[best_idx]
    calibration_nll = nlls[best_idx]
    held_scored = [score_count_report(c.report, temperature; threshold=0.0, true_n=c.true_n) for c in held]
    baseline_nll = mean(_case_nll(c, 1.0) for c in held)
    heldout_nll = mean(s.nll for s in held_scored)
    threshold, source, bins = _threshold_and_bins(held_scored; bin_count=bin_count)
    diag = CountCalibrationDiagnostics(copy(temps), copy(nlls), calibration_nll, heldout_nll,
        baseline_nll, best_idx, source, length(held), length(cal))
    return CountCalibrationModel(temperature, threshold, calibration_nll, heldout_nll,
        baseline_nll, diag, bins)
end
