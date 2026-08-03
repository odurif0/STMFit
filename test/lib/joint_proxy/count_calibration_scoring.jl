function _report_candidates(report)
    hasproperty(report, :candidates) || throw(CountCalibrationError("report missing candidates field"))
    return getproperty(report, :candidates)
end

function _stable_softmax(logits::AbstractVector{<:Real})
    isempty(logits) && throw(CountCalibrationError("no logits provided"))
    m = maximum(logits)
    vals = exp.(Float64.(logits) .- m)
    z = sum(vals)
    z > 0 || throw(CountCalibrationError("numerically unstable softmax: all weights vanished"))
    return vals ./ z
end

function _candidate_rows(report, temperature::Float64)
    temperature > 0 && isfinite(temperature) || throw(CountCalibrationError(
        "temperature must be finite and positive"))
    cands = sort(_finite_candidates(report); by=c -> (c.joint_gcv, c.n))
    min_gcv = minimum(c.joint_gcv for c in cands)
    scale = max(abs(min_gcv), eps(Float64))
    deltas = [(c.joint_gcv - min_gcv) / scale for c in cands]
    logits = -Float64.(deltas) ./ temperature
    probs = _stable_softmax(logits)
    rows = CountPosteriorRow[]
    tol = 1e-12 * max(abs(min_gcv), 1.0)
    for (i, c) in pairs(cands)
        rank = 1 + count(d -> d < deltas[i] - tol, deltas)
        push!(rows, CountPosteriorRow(c.n, c.joint_gcv, deltas[i], probs[i], rank))
    end
    return rows
end

function _pick_prediction(rows::Vector{CountPosteriorRow}, threshold::Float64)
    isempty(rows) && throw(CountCalibrationError("no rows to predict from"))
    threshold >= 0 && threshold <= 1 || throw(CountCalibrationError(
        "confidence threshold must lie in [0, 1]"))
    best = maximum(r.probability for r in rows)
    winners = [r for r in rows if isapprox(r.probability, best; atol=1e-12, rtol=0.0)]
    if best < threshold || length(winners) != 1
        return missing, best, true
    end
    return winners[1].n, best, false
end

function score_count_report(report, temperature::Float64;
                            threshold::Float64=0.5, true_n::Union{Nothing,Int}=nothing)
    rows = _candidate_rows(report, temperature)
    pred, conf, abstained = _pick_prediction(rows, threshold)
    true_rank = 0
    nll = NaN
    if true_n !== nothing
        idx = findfirst(r -> r.n == true_n, rows)
        idx === nothing && throw(CountCalibrationError(
            "missing true candidate: N=$true_n not present among valid candidates"))
        true_rank = rows[idx].rank
        nll = -log(rows[idx].probability)
    end
    return CountCasePosterior("", true_n === nothing ? -1 : true_n, rows, pred, conf, abstained, true_rank, nll)
end

function predict_count(report, model::CountCalibrationModel; true_n::Union{Nothing,Int}=nothing)
    return score_count_report(report, model.temperature; threshold=model.confidence_threshold, true_n=true_n)
end
