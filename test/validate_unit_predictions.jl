#!/usr/bin/env julia

# Validate unknown-production 0/1/? unit-assignment prediction artifacts.
# This script is label-free: it checks TSV integrity and optional feature-key
# coverage only, never truth labels, benchmark manifests, or expected counts.

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _read_tsv

const FORBIDDEN_FLAGS = Set(["--truth", "--manifest", "--full145", "--control-sequence", "--expected-N", "--expected-n", "--control"])
const ALLOWED_PREDICTION_COLUMNS = Set([
    "file", "lobe", "predicted", "confidence", "amplitude", "probability_1",
    "views_used", "invalid_reason", "model", "model_version", "provenance_sha256",
])
const FORBIDDEN_PATH_TOKENS = ("benchmark", "truth", "manifest", "report", "grade", "grading")

struct Options
    predictions::String
    features::String
end

mutable struct Stats
    rows::Int
    files::Int
    predicted_0::Int
    predicted_1::Int
    uncertain::Int
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

function _flag_name(arg::AbstractString)
    return split(String(arg), '='; limit=2)[1]
end

function _parse_cli(args)
    predictions = ""
    features = ""
    i = 1
    while i <= length(args)
        arg = args[i]
        flag = _flag_name(arg)
        if flag in FORBIDDEN_FLAGS
            error("Forbidden benchmark-only flag in unknown prediction validator: $flag")
        elseif arg == "--predictions"
            predictions = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--predictions=")
            predictions = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--features"
            features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--features=")
            features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/validate_unit_predictions.jl --predictions PATH [--features PATH]

            Checks required columns, duplicate (file,lobe) keys, positive and
            contiguous per-file lobe indices, prediction values in 0/1/?,
            confidence/probability ranges, the explicit base/hierarchical
            prediction-column allowlist, forbidden path-valued cells, and
            optional feature-key coverage.

            Label-free constraint: validates artifact structure only; no external
            labels, benchmark metadata, control motifs, or expected counts.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(predictions) && error("--predictions is required")
    isfile(predictions) || error("Prediction TSV not found: $predictions")
    !isempty(features) && !isfile(features) && error("Feature TSV not found: $features")
    return Options(predictions, features)
end

function _require_columns(header::Vector{String}, required, source::AbstractString)
    for col in required
        col in header || error("$source missing required column: $col")
    end
    return nothing
end

function _column_key(column::AbstractString)
    return replace(lowercase(strip(String(column))), r"[^a-z0-9]+" => "_")
end

function _forbidden_column(column::AbstractString)
    key = strip(_column_key(column), '_')
    compact = replace(key, "_" => "")
    compact in ("sequence", "controlsequence", "expectedn", "targetn",
                "benchmarklabel", "benchmarktruth", "truth", "truthlabel",
                "grade", "gradepath", "gradingreport", "report", "reportpath",
                "manifest", "benchmarkmanifest", "denominatormanifest") && return true
    return occursin("control_sequence", key) || occursin("expected_n", key) ||
           occursin("target_n", key) || occursin("benchmark", key) ||
           occursin("truth", key) || startswith(key, "grade_") ||
           endswith(key, "_grade") || startswith(key, "report_") ||
           endswith(key, "_report") || endswith(key, "_sequence")
end

function _reject_forbidden_columns(header::Vector{String}, source::AbstractString)
    hits = sort(filter(_forbidden_column, header))
    isempty(hits) || error("$source has forbidden benchmark/truth/control/grade/report columns: $(join(hits, ", "))")
    return nothing
end

function _reject_unlisted_prediction_columns(header::Vector{String})
    extras = sort(collect(setdiff(Set(header), ALLOWED_PREDICTION_COLUMNS)))
    isempty(extras) || error("Prediction TSV has columns outside the explicit allowlist: $(join(extras, ", "))")
    return nothing
end


function _forbidden_path_value(value::AbstractString)
    text = lowercase(strip(String(value)))
    isempty(text) && return false
    uri_like = occursin(r"^[a-z][a-z0-9+.-]*:", text)
    uri_like && return true
    endswith(text, ".sxm") && !occursin('/', text) && !occursin('\\', text) && return false
    path_like = occursin('/', text) || occursin('\\', text) ||
                occursin(r"\.(tsv|csv|toml|json|jl|md)$", text)
    return path_like && any(token -> occursin(token, text), FORBIDDEN_PATH_TOKENS)
end

function _scan_forbidden_paths(rows, source::AbstractString)
    for (row_index, row) in enumerate(rows), (column, value) in row
        _forbidden_path_value(value) &&
            error("$source row $row_index column $column contains a forbidden benchmark truth/manifest/report/grade path")
    end
    return nothing
end

function _parse_lobe(row, source::AbstractString)
    lobe = tryparse(Int, strip(row["lobe"]))
    lobe === nothing && error("$source has invalid lobe index: $(row["lobe"])")
    lobe >= 1 || error("$source lobe index must be >= 1, got $lobe")
    return lobe
end

function _parse_unit_prediction(value::AbstractString, source::AbstractString)
    pred = strip(value)
    pred in ("0", "1", "?") || error("$source has invalid prediction: $pred")
    return pred
end

function _check_unit_interval(value::AbstractString, label::AbstractString, source::AbstractString)
    text = strip(value)
    text in ("", "NA", "NaN", "nan", "?") && return nothing
    parsed = tryparse(Float64, text)
    parsed === nothing && error("$source has invalid $label value: $value")
    0.0 <= parsed <= 1.0 || error("$source $label must be in [0, 1], got $parsed")
    return nothing
end

function _key(row, source::AbstractString)
    file = basename(strip(row["file"]))
    isempty(file) && error("$source has empty file value")
    return (file, _parse_lobe(row, source))
end

function _load_feature_keys(path::String)
    isempty(path) && return Set{Tuple{String,Int}}()
    header, rows = _read_tsv(path)
    columns = String.(header)
    _require_columns(columns, ("file", "lobe"), "Feature TSV")
    _reject_forbidden_columns(columns, "Feature TSV")
    _scan_forbidden_paths(rows, "Feature TSV")
    return Set(_key(row, "Feature TSV") for row in rows)
end

function _validate_contiguous(by_file::Dict{String,Set{Int}})
    bad = String[]
    for file in sort(collect(keys(by_file)))
        lobes = by_file[file]
        max_lobe = maximum(lobes)
        if length(lobes) != max_lobe || !all(in(lobes), 1:max_lobe)
            push!(bad, "$file($(join(sort(collect(lobes)), ',')))")
        end
    end
    isempty(bad) || error("Noncontiguous lobe indices: $(join(bad[1:min(end, 12)], "; "))")
    return nothing
end

function _validate_predictions(opt::Options)
    header_raw, rows = _read_tsv(opt.predictions)
    header = String.(header_raw)
    _require_columns(header, ("file", "lobe", "predicted"), "Prediction TSV")
    _reject_forbidden_columns(header, "Prediction TSV")
    _reject_unlisted_prediction_columns(header)
    _scan_forbidden_paths(rows, "Prediction TSV")
    isempty(rows) && error("Prediction TSV has no rows: $(opt.predictions)")

    seen = Set{Tuple{String,Int}}()
    by_file = Dict{String,Set{Int}}()
    counts = Dict("0" => 0, "1" => 0, "?" => 0)
    for row in rows
        key = _key(row, "Prediction TSV")
        key in seen && error("Duplicate prediction row for $(key[1]) lobe $(key[2])")
        push!(seen, key)
        push!(get!(by_file, key[1], Set{Int}()), key[2])
        pred = _parse_unit_prediction(row["predicted"], "Prediction TSV $(key[1]) lobe $(key[2])")
        counts[pred] += 1
        "confidence" in header && _check_unit_interval(row["confidence"], "confidence", "Prediction TSV $(key[1]) lobe $(key[2])")
        "probability_1" in header && _check_unit_interval(row["probability_1"], "probability_1", "Prediction TSV $(key[1]) lobe $(key[2])")
    end
    _validate_contiguous(by_file)

    feature_keys = _load_feature_keys(opt.features)
    if !isempty(feature_keys)
        missing = sort(collect(setdiff(feature_keys, seen)))
        extra = sort(collect(setdiff(seen, feature_keys)))
        isempty(missing) || error("Predictions missing feature keys: $(length(missing))")
        isempty(extra) || error("Predictions contain keys absent from features: $(length(extra))")
    end
    return Stats(length(rows), length(by_file), counts["0"], counts["1"], counts["?"])
end

function main(args=ARGS)
    opt = _parse_cli(args)
    stats = _validate_predictions(opt)
    println("Unit prediction validation")
    println("  predictions: ", opt.predictions)
    println("  features:    ", isempty(opt.features) ? "none" : opt.features)
    println("  files:       ", stats.files)
    println("  rows:        ", stats.rows)
    println("  predicted_0: ", stats.predicted_0)
    println("  predicted_1: ", stats.predicted_1)
    println("  uncertain:   ", stats.uncertain)
    println("  status:      ok")
    return nothing
end

main()
