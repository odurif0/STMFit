module UnitAssignmentReport

include(joinpath(@__DIR__, "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv

export load_truth, write_position_rows

_pct(num::Integer, den::Integer) = den > 0 ? string(round(100 * num / den; digits=1), "%") : "NA"
_basename_file(s::AbstractString) = basename(strip(s))
_flip_pred(x) = ismissing(x) ? missing : 1 - x

function _parse_seq(s::AbstractString)
    return [parse(Int, c) for c in strip(s) if c in "01"]
end

function _parse_prediction(s::AbstractString)
    t = strip(s)
    (isempty(t) || t in ("?", "NA", "NaN", "nan")) && return missing
    parsed = tryparse(Int, t)
    parsed === nothing && error("Invalid predicted label: $s")
    parsed in (0, 1) || error("Predicted label must be 0, 1, or ?: $s")
    return parsed
end

function load_truth(path::String, primary_only::Bool)
    _, rows = _read_tsv(path)
    truth = Dict{String,Vector{Int}}()
    for row in rows
        quality = get(row, "quality", "clean")
        primary_only && quality == "poor_quality" && continue
        seq = _parse_seq(row["sequence"])
        isempty(seq) && continue
        truth[_basename_file(row["file"])] = seq
    end
    return truth
end

function _load_predictions(path::String)
    _, rows = _read_tsv(path)
    by_file = Dict{String,Vector{Dict{String,String}}}()
    for row in rows
        file = _basename_file(row["file"])
        push!(get!(by_file, file, Dict{String,String}[]), row)
    end
    for rows_for_file in values(by_file)
        sort!(rows_for_file, by=r -> parse(Int, r["lobe"]))
    end
    return rows, by_file
end

function _physical_remap(rows::Vector{Dict{String,String}})
    amps = Dict{Int,Vector{Float64}}()
    for row in rows
        label = _parse_prediction(row["predicted"])
        ismissing(label) && continue
        amp = tryparse(Float64, strip(get(row, "amplitude", "")))
        amp === nothing && continue
        isfinite(amp) || continue
        push!(get!(amps, label, Float64[]), amp)
    end
    means = Dict(k => sum(v) / length(v) for (k, v) in amps if !isempty(v))
    length(keys(means)) < 2 && return identity
    high_label = first(sort(collect(keys(means)); by=k -> means[k], rev=true))
    return high_label == 1 ? identity : _flip_pred
end

function _alignment_score(pred, truth::Vector{Int})
    n = min(length(pred), length(truth))
    correct = classified = 0
    for i in 1:n
        ismissing(pred[i]) && continue
        classified += 1
        correct += pred[i] == truth[i] ? 1 : 0
    end
    return correct, classified
end

function _best_physical_alignment(pred, truth::Vector{Int})
    candidates = (pred, reverse(pred))
    scores = [_alignment_score(c, truth) for c in candidates]
    best = _accuracy(scores[2]) > _accuracy(scores[1]) ? 2 : 1
    return candidates[best]
end

_accuracy(score) = score[2] > 0 ? score[1] / score[2] : 0.0

function _reject_symlink_leaf(path::String)
    islink(path) && error("Refusing to overwrite symlink: $path")
    return nothing
end

function write_position_rows(profile_name::String, prediction_path::String, truth, out_path::String;
                             include_missing_truth_positions::Bool=false)
    pred_rows, by_file = _load_predictions(prediction_path)
    remap = _physical_remap(pred_rows)
    stats = Dict{Int,Dict{String,Int}}()
    for file in sort(collect(keys(truth)))
        rows = get(by_file, file, Dict{String,String}[])
        isempty(rows) && !include_missing_truth_positions && continue
        pred = [remap(_parse_prediction(r["predicted"])) for r in rows]
        aligned = _best_physical_alignment(pred, truth[file])
        limit = include_missing_truth_positions ? length(truth[file]) : min(length(aligned), length(truth[file]))
        for i in 1:limit
            d = get!(stats, i, Dict("possible" => 0, "classified" => 0, "correct" => 0, "wrong" => 0, "uncertain" => 0))
            d["possible"] += 1
            pred_i = i <= length(aligned) ? aligned[i] : missing
            if ismissing(pred_i)
                d["uncertain"] += 1
            else
                d["classified"] += 1
                d[pred_i == truth[file][i] ? "correct" : "wrong"] += 1
            end
        end
    end
    _ensure_parent(out_path)
    _reject_symlink_leaf(out_path)
    open(out_path, "a") do io
        for pos in sort(collect(keys(stats)))
            d = stats[pos]
            println(io, join([profile_name, pos, d["possible"], d["classified"], d["correct"], d["wrong"],
                              d["uncertain"], _pct(d["correct"], d["classified"]), _pct(d["correct"], d["possible"])], '\t'))
        end
    end
end

end
