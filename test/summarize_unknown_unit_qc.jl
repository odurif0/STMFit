#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv
using Statistics

struct Options
    predictions::String
    out_tsv::String
    profile::String
    plots_dir::String
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

function _parse_cli(args)
    predictions = ""
    out_tsv = ""
    profile = "auto"
    plots_dir = ""
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--predictions"
            predictions = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--predictions=")
            predictions = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--profile"
            profile = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--profile=")
            profile = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--plots-dir"
            plots_dir = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--plots-dir=")
            plots_dir = split(arg, "="; limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/summarize_unknown_unit_qc.jl --predictions PATH --out PATH [--profile NAME] [--plots-dir DIR]

            Default label-free review flags:
              high_uncertain_fraction: uncertain_fraction >= 0.50
              low_mean_confidence: finite mean_confidence < 0.60
              missing_optional_view: views_used lower than profile view count
              noncontiguous_lobes: per-file lobe index gap
              n_outlier: >=8 files and N outside Q1 - 1.5*IQR, Q3 + 1.5*IQR
              plot_missing: expected standalone plot absent/empty when --plots-dir is supplied
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(predictions) && error("--predictions is required")
    isempty(out_tsv) && error("--out is required")
    isfile(predictions) || error("Prediction TSV not found: $predictions")
    profile in ("auto", "base_bwd_consensus", "base_split_log_skew", "base") || error("Unsupported --profile: $profile")
    return Options(predictions, out_tsv, profile, plots_dir)
end

function _required_cols(header::Vector{String}, cols)
    for col in cols
        col in header || error("Predictions missing required column: $col")
    end
    return nothing
end

function _profile_view_count(opt::Options)
    opt.profile == "base_bwd_consensus" && return 3
    opt.profile == "base_split_log_skew" && return 1
    opt.profile == "base" && return 1
    name = lowercase(basename(opt.predictions))
    occursin("base_bwd_consensus", name) && return 3
    occursin("base_split_log_skew", name) && return 1
    return 0
end

function _parse_conf(row)
    haskey(row, "confidence") || return NaN
    text = strip(row["confidence"])
    text in ("", "NA", "NaN", "nan", "?") && return NaN
    return parse(Float64, text)
end

function _parse_views(row)
    haskey(row, "views_used") || return typemax(Int)
    parsed = tryparse(Int, strip(row["views_used"]))
    parsed === nothing && error("Invalid views_used: $(row["views_used"])")
    return parsed
end

function _prediction(row)
    pred = strip(row["predicted"])
    pred in ("0", "1", "?") || error("Invalid prediction: $pred")
    return pred
end

function _plot_path(plots_dir::String, file::String)
    isempty(plots_dir) && return ""
    return joinpath(plots_dir, "standalone", replace(file, ".sxm" => "_chain.png"))
end

function _has_noncontiguous(lobes::Vector{Int})
    isempty(lobes) && return true
    uniq = sort(unique(lobes))
    return length(uniq) != maximum(uniq) || !all(in(uniq), 1:maximum(uniq))
end

function _n_fences(ns::Vector{Int})
    length(ns) < 8 && return (-Inf, Inf)
    q1 = quantile(ns, 0.25)
    q3 = quantile(ns, 0.75)
    iqr = q3 - q1
    return (q1 - 1.5 * iqr, q3 + 1.5 * iqr)
end

function _summarize(opt::Options)
    header_raw, rows = _read_tsv(opt.predictions)
    header = String.(header_raw)
    _required_cols(header, ("file", "lobe", "predicted"))
    by_file = Dict{String,Vector{Dict{String,String}}}()
    for row in rows
        file = basename(strip(row["file"]))
        isempty(file) && error("Empty file value")
        push!(get!(by_file, file, Dict{String,String}[]), row)
    end
    ns = [length(v) for v in values(by_file)]
    low_n, high_n = _n_fences(ns)
    expected_views = _profile_view_count(opt)
    return sort(collect(keys(by_file))), by_file, low_n, high_n, expected_views
end

function _write_queue(opt::Options)
    files, by_file, low_n, high_n, expected_views = _summarize(opt)
    _ensure_parent(opt.out_tsv)
    islink(opt.out_tsv) && error("Refusing to overwrite symlink: $(opt.out_tsv)")
    open(opt.out_tsv, "w") do io
        println(io, join(["file", "N_predicted_lobes", "uncertain_fraction", "mean_confidence", "review_status", "review_reasons"], '\t'))
        for file in files
            rows = by_file[file]
            n = length(rows)
            uncertain = count(row -> _prediction(row) == "?", rows)
            confidences = [_parse_conf(row) for row in rows]
            finite_conf = filter(isfinite, confidences)
            mean_conf = isempty(finite_conf) ? NaN : mean(finite_conf)
            lobes = [parse(Int, strip(row["lobe"])) for row in rows]
            reasons = String[]
            uncertain / n >= 0.50 && push!(reasons, "high_uncertain_fraction")
            isfinite(mean_conf) && mean_conf < 0.60 && push!(reasons, "low_mean_confidence")
            expected_views > 0 && any(row -> _parse_views(row) < expected_views, rows) && push!(reasons, "missing_optional_view")
            _has_noncontiguous(lobes) && push!(reasons, "noncontiguous_lobes")
            (n < low_n || n > high_n) && push!(reasons, "n_outlier")
            plot = _plot_path(opt.plots_dir, file)
            !isempty(plot) && (!isfile(plot) || filesize(plot) == 0) && push!(reasons, "plot_missing")
            status = isempty(reasons) ? "ok" : "review"
            reason_text = isempty(reasons) ? "ok" : join(reasons, ",")
            conf_text = isfinite(mean_conf) ? string(round(mean_conf; digits=8)) : "NA"
            println(io, join([file, string(n), string(round(uncertain / n; digits=8)), conf_text, status, reason_text], '\t'))
        end
    end
    return length(files)
end

function main(args=ARGS)
    opt = _parse_cli(args)
    files = _write_queue(opt)
    println("Unknown unit QC review queue")
    println("  predictions: ", opt.predictions)
    println("  out:         ", opt.out_tsv)
    println("  files:       ", files)
end

main()
