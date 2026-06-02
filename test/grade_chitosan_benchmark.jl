#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# grade_chitosan_benchmark.jl — External grading of chitosan benchmark results
#
# This script evaluates already-produced result TSVs against the benchmark
# manifest.  It is intentionally outside the fitting/model-selection path: it
# may use human validation targets to grade candidate methods, but those labels
# must not be used by the fitter itself.
#
# Examples:
#   julia --project=. test/grade_chitosan_benchmark.jl \
#       --results results/best_plots/summary_overlap060_hard.tsv \
#       --column N_eff \
#       --out results/benchmark_grades/baseline_N_eff.tsv
#
#   julia --project=. test/grade_chitosan_benchmark.jl \
#       --results results/resolved_lobes_audit/resolved_lobes.tsv \
#       --file-column file \
#       --column N_resolved \
#       --selected-column is_selected \
#       --out results/benchmark_grades/resolved_selected.tsv
# ──────────────────────────────────────────────────────────────────────────────

using Printf, TOML

const DEFAULT_MANIFEST = "benchmarks/chitosan_240817.toml"
const DEFAULT_RESULTS = "results/best_plots/summary_overlap060_hard.tsv"
const DEFAULT_OUT = "results/benchmark_grades/chitosan_grade.tsv"

struct Options
    manifest::String
    results::String
    out_tsv::String
    file_column::String
    result_column::String
    selected_column::Union{Nothing,String}
    target_n::Int
end

function _parse_cli(args)
    manifest = DEFAULT_MANIFEST
    results = DEFAULT_RESULTS
    out_tsv = DEFAULT_OUT
    file_column = "filepath"
    result_column = "N_eff"
    selected_column::Union{Nothing,String} = nothing
    target_n = 6

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--manifest"
            i < length(args) || error("--manifest requires a path")
            manifest = args[i+1]; i += 2
        elseif startswith(arg, "--manifest=")
            manifest = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--results"
            i < length(args) || error("--results requires a TSV path")
            results = args[i+1]; i += 2
        elseif startswith(arg, "--results=")
            results = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            i < length(args) || error("--out requires a TSV path")
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--file-column"
            i < length(args) || error("--file-column requires a name")
            file_column = args[i+1]; i += 2
        elseif startswith(arg, "--file-column=")
            file_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--column"
            i < length(args) || error("--column requires a name")
            result_column = args[i+1]; i += 2
        elseif startswith(arg, "--column=")
            result_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--selected-column"
            i < length(args) || error("--selected-column requires a name")
            selected_column = args[i+1]; i += 2
        elseif startswith(arg, "--selected-column=")
            selected_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--target-N"
            i < length(args) || error("--target-N requires an integer")
            target_n = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--target-N=")
            target_n = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/grade_chitosan_benchmark.jl [options]

            Options:
              --manifest PATH          Benchmark manifest [$(DEFAULT_MANIFEST)]
              --results PATH           Result TSV to grade [$(DEFAULT_RESULTS)]
              --out PATH               Per-file grade TSV [$(DEFAULT_OUT)]
              --file-column NAME       File column in result TSV [filepath]
              --column NAME            Integer result column to grade [N_eff]
              --selected-column NAME   Optional boolean column; if present, only selected rows are graded
              --target-N INT           Default human-validation target for primary files [6]
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    return Options(manifest, results, out_tsv, file_column, result_column,
                   selected_column, target_n)
end

function _read_tsv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("Empty TSV: $path")
    header = split(lines[1], '\t'; keepempty=true)
    rows = Vector{Dict{String,String}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (j, name) in enumerate(header)
            row[name] = j <= length(vals) ? vals[j] : ""
        end
        push!(rows, row)
    end
    return header, rows
end

_as_bool(s::AbstractString) = lowercase(strip(s)) in ("true", "t", "1", "yes", "y")

function _parse_int_or_missing(s::AbstractString)
    t = strip(s)
    if isempty(t) || t in ("NA", "NaN", "missing")
        return missing
    end
    try
        return parse(Int, t)
    catch
        try
            return round(Int, parse(Float64, t))
        catch
            return missing
        end
    end
end

function _basename_file(s::AbstractString)
    b = basename(strip(s))
    if isempty(b)
        return b
    end
    m = match(r"(\d{6}[^/\\]*\.sxm)", b)
    return m === nothing ? b : m.captures[1]
end

function _quality_info(manifest::Dict{String,Any}, file::String)
    files = get(manifest, "files", Dict{String,Any}())
    info = get(files, file, Dict{String,Any}())
    quality = String(get(info, "quality", "clean"))
    exclude_primary = Bool(get(info, "exclude_primary", false))
    target = get(info, "target_N_for_human_validation", missing)
    target = target === missing ? missing : Int(target)
    return info, quality, exclude_primary, target
end

function _result_rows_by_file(rows, opt::Options)
    by_file = Dict{String,Dict{String,String}}()
    for row in rows
        haskey(row, opt.file_column) || error("Missing file column $(opt.file_column)")
        haskey(row, opt.result_column) || error("Missing result column $(opt.result_column)")
        if opt.selected_column !== nothing
            haskey(row, opt.selected_column) || error("Missing selected column $(opt.selected_column)")
            _as_bool(row[opt.selected_column]) || continue
        end
        file = _basename_file(row[opt.file_column])
        by_file[file] = row
    end
    return by_file
end

function _pct(n::Int, d::Int)
    d == 0 && return "NA"
    return @sprintf("%.1f%%", 100n / d)
end

function _write_report(opt::Options, manifest, by_file)
    eval_cfg = get(manifest, "evaluation", Dict{String,Any}())
    manifest_files = get(manifest, "files", Dict{String,Any}())
    grade_only_manifest_files = Bool(get(eval_cfg, "grade_only_manifest_files", false))
    primary_exclude_quality = Set(String.(get(eval_cfg, "primary_exclude_quality", ["poor_quality", "excluded"])))
    stress_quality = Set(String.(get(eval_cfg, "stress_quality", ["poor_quality"])))

    files = sort(collect(keys(by_file)))
    mkpath(dirname(opt.out_tsv))

    primary_total = primary_ok = primary_pm1 = primary_over = primary_under = 0
    target_total = target_ok = target_pm1 = target_over = target_under = 0
    stress_total = 0
    missing_result = String[]
    primary_failures = String[]
    target_failures = String[]

    open(opt.out_tsv, "w") do io
        println(io, join([
            "file", "quality", "set", "target_N", "observed_N", "pass",
            "pass_pm1", "source_column", "status", "classification", "ambiguous_eff",
            "runnerup_N_eff", "delta_GCV_rel_eff"
        ], '\t'))

        for file in files
            grade_only_manifest_files && !haskey(manifest_files, file) && continue
            row = by_file[file]
            _, quality, exclude_primary, manifest_target = _quality_info(manifest, file)
            is_stress = quality in stress_quality
            is_primary = !(exclude_primary || quality in primary_exclude_quality)
            is_target = quality == "clean_target"
            set_name = is_target ? "target" : (is_primary ? "primary" : (is_stress ? "stress" : "excluded"))
            target_n = manifest_target === missing ? opt.target_n : manifest_target
            observed = _parse_int_or_missing(row[opt.result_column])
            passed = observed !== missing && observed == target_n
            passed_pm1 = observed !== missing && abs(observed - target_n) <= 1

            if observed === missing && (is_primary || is_target || is_stress)
                push!(missing_result, file)
            end
            if is_primary
                primary_total += 1
                if passed
                    primary_ok += 1
                else
                    push!(primary_failures, file)
                end
                passed_pm1 && (primary_pm1 += 1)
                if observed !== missing
                    observed > target_n && (primary_over += 1)
                    observed < target_n && (primary_under += 1)
                end
            end
            if is_target
                target_total += 1
                if passed
                    target_ok += 1
                else
                    push!(target_failures, file)
                end
                passed_pm1 && (target_pm1 += 1)
                if observed !== missing
                    observed > target_n && (target_over += 1)
                    observed < target_n && (target_under += 1)
                end
            end
            if is_stress
                stress_total += 1
            end

            println(io, join([
                file,
                quality,
                set_name,
                string(target_n),
                observed === missing ? "NA" : string(observed),
                string(passed),
                string(passed_pm1),
                opt.result_column,
                get(row, "status", ""),
                get(row, "classification", ""),
                get(row, "ambiguous_eff", ""),
                get(row, "runnerup_N_eff", ""),
                get(row, "delta_GCV_rel_eff", ""),
            ], '\t'))
        end
    end

    println("Benchmark grade")
    println("  manifest:       ", opt.manifest)
    println("  results:        ", opt.results)
    println("  result column:  ", opt.result_column)
    println("  output:         ", opt.out_tsv)
    println("  primary score:  ", primary_ok, "/", primary_total, " (", _pct(primary_ok, primary_total), ")")
    println("  primary ±1:     ", primary_pm1, "/", primary_total, " (", _pct(primary_pm1, primary_total), ")")
    println("  primary over/under: ", primary_over, "/", primary_under)
    println("  target score:   ", target_ok, "/", target_total, " (", _pct(target_ok, target_total), ")")
    println("  target ±1:      ", target_pm1, "/", target_total, " (", _pct(target_pm1, target_total), ")")
    println("  target over/under:  ", target_over, "/", target_under)
    println("  stress rows:    ", stress_total)
    !isempty(primary_failures) && println("  primary fails:  ", join(primary_failures, ", "))
    !isempty(target_failures) && println("  target fails:   ", join(target_failures, ", "))
    !isempty(missing_result) && println("  missing result: ", join(missing_result, ", "))
end

function main(args=ARGS)
    opt = _parse_cli(args)
    manifest = TOML.parsefile(opt.manifest)
    _, rows = _read_tsv(opt.results)
    by_file = _result_rows_by_file(rows, opt)
    isempty(by_file) && error("No result rows to grade after filtering")
    _write_report(opt, manifest, by_file)
end

main()
