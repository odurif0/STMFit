#!/usr/bin/env julia

# Build abstention variants from an existing 0/1 unit-assignment prediction TSV.
#
# This script is label-free: it does not read the benchmark truth. It keeps the
# base prediction when confidence and optional auxiliary-model agreement pass the
# requested gates; otherwise it writes "?". Grade the resulting TSV separately
# with grade_unit_assignment.jl.

using Printf

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv

const DEFAULT_BASE = "results/unit_assignment/best_labelfree_ensemble3_forced_predictions.tsv"
const DEFAULT_OUT = "results/unit_assignment/best_labelfree_ensemble3_abstain_variant_predictions.tsv"

struct Options
    base::String
    out_tsv::String
    confidence_min::Float64
    agree_with::Vector{String}
    min_agree::Int
end

function _parse_cli(args)
    base = DEFAULT_BASE
    out_tsv = DEFAULT_OUT
    confidence_min = 0.80
    agree_with = String[]
    min_agree = -1

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--base"
            base = args[i+1]; i += 2
        elseif startswith(arg, "--base=")
            base = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--confidence-min"
            confidence_min = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--confidence-min=")
            confidence_min = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--agree-with"
            push!(agree_with, args[i+1]); i += 2
        elseif startswith(arg, "--agree-with=")
            push!(agree_with, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--min-agree"
            min_agree = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--min-agree=")
            min_agree = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/build_unit_abstention_variants.jl [options]

            Options:
              --base PATH             Base forced predictions TSV [$(DEFAULT_BASE)]
              --out PATH              Output abstention TSV [$(DEFAULT_OUT)]
              --confidence-min FLOAT  Keep only base rows with confidence >= FLOAT [0.80]
              --agree-with PATH       Auxiliary prediction TSV to require agreement with. Repeatable.
              --min-agree INT         Minimum number of auxiliary predictions that must agree.
                                      Default: all --agree-with files.

            Input TSVs must contain file, lobe, predicted. The base TSV must also
            contain confidence. The output preserves base columns, rewrites
            predicted to ? on rejected rows, and appends forced_predicted,
            agreement_count, agreement_available, and abstention_reason.

            Label-free constraint: this script never reads truth labels. Use
            test/grade_unit_assignment.jl afterwards for post-hoc grading.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    0.0 <= confidence_min <= 1.0 || error("--confidence-min must be in [0, 1]")
    min_agree = min_agree < 0 ? length(agree_with) : min_agree
    0 <= min_agree <= length(agree_with) || error("--min-agree must be between 0 and number of --agree-with files")
    isfile(base) || error("Base predictions not found: $base")
    for path in agree_with
        isfile(path) || error("Auxiliary predictions not found: $path")
    end
    return Options(base, out_tsv, confidence_min, agree_with, min_agree)
end

_key(row) = (basename(strip(row["file"])), parse(Int, strip(row["lobe"])))

function _parse_prediction(s::AbstractString)
    t = strip(s)
    t == "?" && return missing
    parsed = tryparse(Int, t)
    parsed === nothing && error("Invalid predicted label: $s")
    parsed in (0, 1) || error("Predicted label must be 0, 1, or ?: $s")
    return parsed
end

function _load_aux(path::String)
    _, rows = _read_tsv(path)
    aux = Dict{Tuple{String,Int},Union{Missing,Int}}()
    for row in rows
        haskey(row, "file") || error("Missing file column in $path")
        haskey(row, "lobe") || error("Missing lobe column in $path")
        haskey(row, "predicted") || error("Missing predicted column in $path")
        aux[_key(row)] = _parse_prediction(row["predicted"])
    end
    return aux
end

function _row_values(header::AbstractVector{<:AbstractString}, row::Dict{String,String})
    return [get(row, h, "") for h in header]
end

function main(args=ARGS)
    opt = _parse_cli(args)
    header_raw, rows = _read_tsv(opt.base)
    header = String.(header_raw)
    isempty(rows) && error("No prediction rows in $(opt.base)")
    for col in ("file", "lobe", "predicted", "confidence")
        col in header || error("Base predictions missing column: $col")
    end

    auxs = [_load_aux(path) for path in opt.agree_with]
    out_header = copy(header)
    for col in ("forced_predicted", "agreement_count", "agreement_available", "abstention_reason")
        col in out_header || push!(out_header, col)
    end

    n_total = n_kept = n_low_conf = n_disagree = 0
    _ensure_parent(opt.out_tsv)
    open(opt.out_tsv, "w") do io
        println(io, join(out_header, '\t'))
        for row in rows
            n_total += 1
            forced = _parse_prediction(row["predicted"])
            ismissing(forced) && error("Base predictions must be forced 0/1, got ?: $(row["file"]) lobe $(row["lobe"])")
            confidence = parse(Float64, strip(row["confidence"]))
            agreement_count = 0
            agreement_available = 0
            for aux in auxs
                k = _key(row)
                haskey(aux, k) || continue
                aux_pred = aux[k]
                ismissing(aux_pred) && continue
                agreement_available += 1
                agreement_count += aux_pred == forced ? 1 : 0
            end

            reason = "accepted"
            keep = true
            if confidence < opt.confidence_min
                keep = false
                reason = "low_confidence"
                n_low_conf += 1
            elseif agreement_count < opt.min_agree
                keep = false
                reason = "aux_disagreement"
                n_disagree += 1
            else
                n_kept += 1
            end

            out_row = copy(row)
            out_row["forced_predicted"] = string(forced)
            out_row["predicted"] = keep ? string(forced) : "?"
            out_row["agreement_count"] = string(agreement_count)
            out_row["agreement_available"] = string(agreement_available)
            out_row["abstention_reason"] = reason
            println(io, join(_row_values(out_header, out_row), '\t'))
        end
    end

    println("Unit abstention variant")
    println("  base:           ", opt.base)
    println("  out:            ", opt.out_tsv)
    println("  confidence_min: ", @sprintf("%.3f", opt.confidence_min))
    println("  agree_with:     ", isempty(opt.agree_with) ? "none" : join(opt.agree_with, ", "))
    println("  min_agree:      ", opt.min_agree)
    println("  kept:           ", n_kept, "/", n_total)
    println("  abstained:      ", n_total - n_kept, " (low_confidence=", n_low_conf,
            ", aux_disagreement=", n_disagree, ")")
end

main()
