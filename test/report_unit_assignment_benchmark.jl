#!/usr/bin/env julia

# Reproducible report for chitosan 0/1/? unit-assignment benchmark grades.
# This is post-hoc grading/reporting only: truth labels never enter prediction.

using Printf, TOML

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv
include(joinpath(@__DIR__, "lib", "unit_assignment_report.jl"))
using .UnitAssignmentReport: load_truth, write_position_rows

const DEFAULT_TRUTH = "benchmarks/chitosan_240817_unit_sequences.tsv"
const DEFAULT_MANIFEST = "benchmarks/chitosan_6mer_counting_confirmed.toml"
const DEFAULT_OUTDIR = "results/unit_assignment/benchmark_report"
const DEFAULT_CONTROL_SEQUENCE = "010010"
const DEFAULT_PROFILES = [
    ("forced_ensemble3", "results/unit_assignment/best_labelfree_ensemble3_forced_predictions.tsv",
     "Frozen full-coverage three-view label-free ensemble."),
    ("agreebase65", "results/unit_assignment/best_labelfree_ensemble3_abstain_agreebase65_predictions.tsv",
     "Frozen 0/1/? profile: confidence >= 0.65 plus base-model agreement."),
    ("err05", "results/unit_assignment/best_labelfree_ensemble3_abstain_err05_predictions.tsv",
     "Frozen strict emitted-error profile, not a newly tuned threshold."),
]

struct ProfileSpec
    name::String
    path::String
    note::String
end

struct Options
    truth::String
    manifest::String
    outdir::String
    profiles::Vector{ProfileSpec}
    primary_only::Bool
    full145::Bool
    full145_own_n::Bool
    control_sequence::String
end

_pct(num::Integer, den::Integer) = den > 0 ? @sprintf("%.1f%%", 100 * num / den) : "NA"

function _default_profiles()
    return [ProfileSpec(name, path, note) for (name, path, note) in DEFAULT_PROFILES]
end

function _parse_profile(s::AbstractString)
    parts = split(String(s), '='; limit=3)
    length(parts) >= 2 || error("--profile must be NAME=PATH or NAME=PATH=NOTE")
    name, path = strip(parts[1]), strip(parts[2])
    isempty(name) && error("--profile name is empty")
    occursin(r"^[A-Za-z0-9._-]+$", name) || error("--profile name may contain only letters, digits, '.', '_', and '-': $name")
    isempty(path) && error("--profile path is empty")
    note = length(parts) == 3 ? strip(parts[3]) : "Custom frozen benchmark profile."
    return ProfileSpec(name, path, note)
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

function _parse_cli(args)
    truth = DEFAULT_TRUTH
    manifest = DEFAULT_MANIFEST
    outdir = DEFAULT_OUTDIR
    profiles = ProfileSpec[]
    primary_only = true
    full145 = false
    full145_own_n = false
    control_sequence = DEFAULT_CONTROL_SEQUENCE

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--truth"
            truth = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--truth=")
            truth = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--manifest"
            manifest = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--manifest=")
            manifest = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--outdir"
            outdir = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--full145"
            full145 = true; i += 1
        elseif arg == "--full145-own-n"
            full145 = true; full145_own_n = true; i += 1
        elseif arg == "--control-sequence"
            control_sequence = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--control-sequence=")
            control_sequence = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--profile"
            push!(profiles, _parse_profile(_arg_value(args, i, arg))); i += 2
        elseif startswith(arg, "--profile=")
            push!(profiles, _parse_profile(split(arg, "="; limit=2)[2])); i += 1
        elseif arg == "--include-stress"
            primary_only = false; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/report_unit_assignment_benchmark.jl [options]

            Options:
              --truth PATH          Unit truth TSV for subset mode [$(DEFAULT_TRUTH)]
              --full145             Grade the full 145-file benchmark using manifest control NKNNKN.
              --full145-own-n       Grade full145 predictions that use each file's own N_selected.
                                    Missing control positions count uncertain; extra predicted lobes are reported.
              --manifest PATH       Full benchmark manifest [$(DEFAULT_MANIFEST)]
              --control-sequence S  External control sequence for --full145 [$(DEFAULT_CONTROL_SEQUENCE)]
              --outdir PATH         Report output directory [$(DEFAULT_OUTDIR)]
              --profile NAME=PATH   Profile to grade. Repeatable. Defaults to the three frozen profiles.
                                    Optional form: NAME=PATH=NOTE
              --include-stress      Include poor_quality truth rows instead of strict primary-only grading.

            Outputs:
              summary.tsv, lobe_position_errors.tsv, report.md, and per-profile grade TSVs.

            Label-free constraint: this script is post-hoc benchmark reporting only.
            It reads truth for grading, never to tune predictions or thresholds.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    profiles = isempty(profiles) ? _default_profiles() : profiles
    full145 ? (isfile(manifest) || error("Benchmark manifest not found: $manifest")) : (isfile(truth) || error("Truth TSV not found: $truth"))
    occursin(r"^[01]+$", control_sequence) || error("--control-sequence must contain only 0/1 digits")
    length(control_sequence) == 6 || error("--control-sequence must have length 6 for the 6mer benchmark")
    for p in profiles
        isfile(p.path) || error("Prediction TSV not found for $(p.name): $(p.path)")
    end
    return Options(truth, manifest, outdir, profiles, primary_only, full145, full145_own_n, control_sequence)
end

function _reject_symlink_leaf(path::String)
    islink(path) && error("Refusing to overwrite symlink: $path")
    return nothing
end

function _grade_path(outdir::String, profile_name::String)
    grade_dir = abspath(joinpath(outdir, "grades"))
    path = abspath(joinpath(grade_dir, profile_name * ".tsv"))
    startswith(path, grade_dir * Base.Filesystem.path_separator) || error("Grade path escaped output directory: $path")
    _reject_symlink_leaf(path)
    return path
end

function _write_full145_truth(opt::Options)
    path = abspath(joinpath(opt.outdir, "control_full145_truth.tsv"))
    _reject_symlink_leaf(path)
    manifest = TOML.parsefile(opt.manifest)
    files = manifest["files"]
    _ensure_parent(path)
    open(path, "w") do io
        println(io, join(["file", "sequence", "quality", "target_N", "notes"], '\t'))
        for file in sort(collect(keys(files)))
            row = files[file]
            quality = String(get(row, "quality", "clean_target"))
            target_n = string(get(row, "target_N_for_human_validation", 6))
            relative_path = String(get(row, "relative_path", file))
            note = "full145 external control NKNNKN; grader-only, label-free method must not read"
            println(io, join([relative_path, opt.control_sequence, quality, target_n, note], '\t'))
        end
    end
    return path
end

function _run_grader(profile::ProfileSpec, opt::Options, truth_path::String)
    grade_path = _grade_path(opt.outdir, profile.name)
    _ensure_parent(grade_path)
    script = joinpath(@__DIR__, "grade_unit_assignment.jl")
    cmd = opt.primary_only ?
        `$(Base.julia_cmd()) $script --predictions $(profile.path) --truth $truth_path --out $grade_path --primary-only` :
        `$(Base.julia_cmd()) $script --predictions $(profile.path) --truth $truth_path --out $grade_path`
    run(cmd)
    return grade_path
end

function _summarize_grade(profile::ProfileSpec, grade_path::String, opt::Options, truth_path::String)
    _, rows = _read_tsv(grade_path)
    possible = sum(parse(Int, r["N_truth"]) for r in rows)
    prediction_lobes = sum(parse(Int, r["N_pred"]) for r in rows)
    missing_truth_lobes = sum(max(parse(Int, r["N_truth"]) - parse(Int, r["N_pred"]), 0) for r in rows)
    extra_prediction_lobes = sum(max(parse(Int, r["N_pred"]) - parse(Int, r["N_truth"]), 0) for r in rows)
    files_short_n = count(r -> parse(Int, r["N_pred"]) < parse(Int, r["N_truth"]), rows)
    files_extra_n = count(r -> parse(Int, r["N_pred"]) > parse(Int, r["N_truth"]), rows)
    classified = sum(parse(Int, r["N_classified"]) for r in rows)
    correct = sum(parse(Int, r["phys_correct"]) for r in rows)
    exact = count(r -> strip(r["phys_seq_match"]) == "true", rows)
    oracle_total = sum(parse(Int, r["N_classified"]) for r in rows)
    oracle_correct = sum(parse(Int, r["oracle_correct"]) for r in rows)
    errors = classified - correct
    if opt.full145 && (length(rows) != 145 || possible != 870)
        error("Profile $(profile.name) did not grade the full145 denominator: got $(length(rows)) files / $possible lobes")
    elseif !opt.full145 && opt.primary_only && truth_path == DEFAULT_TRUTH && (length(rows) != 35 || possible != 210)
        error("Profile $(profile.name) did not grade the strict default denominator: got $(length(rows)) files / $possible lobes")
    end
    return Dict(
        "profile" => profile.name,
        "prediction_path" => profile.path,
        "note" => profile.note,
        "files" => string(length(rows)),
        "lobes_possible" => string(possible),
        "prediction_lobes" => string(prediction_lobes),
        "missing_truth_lobes" => string(missing_truth_lobes),
        "extra_prediction_lobes" => string(extra_prediction_lobes),
        "files_short_N" => string(files_short_n),
        "files_extra_N" => string(files_extra_n),
        "lobes_classified" => string(classified),
        "coverage" => _pct(classified, possible),
        "physical_correct" => string(correct),
        "physical_errors_emitted" => string(errors),
        "physical_accuracy_classified" => _pct(correct, classified),
        "honest_correct" => string(correct),
        "honest_uncertain" => string(possible - correct),
        "honest_correct_frac" => _pct(correct, possible),
        "sequence_exact" => string(exact),
        "oracle_correct" => string(oracle_correct),
        "oracle_accuracy_classified" => _pct(oracle_correct, oracle_total),
        "grade_tsv" => grade_path,
    )
end

function _format_examples(values; limit::Int=12)
    items = collect(values)
    isempty(items) && return "none"
    shown = join(items[1:min(end, limit)], ", ")
    length(items) > limit && return shown * ", ..."
    return shown
end

function _preflight_full145_profile(profile::ProfileSpec, truth; allow_variable_n::Bool=false)
    _, rows = _read_tsv(profile.path)
    by_file = Dict{String,Set{Int}}()
    duplicate_lobes = String[]
    invalid_lobes = String[]
    for row in rows
        haskey(row, "file") || error("Profile $(profile.name) is missing required column: file")
        haskey(row, "lobe") || error("Profile $(profile.name) is missing required column: lobe")
        file = basename(strip(row["file"]))
        lobe = tryparse(Int, strip(row["lobe"]))
        lobe === nothing && error("Profile $(profile.name) has invalid lobe index for $(row["file"]): $(row["lobe"])")
        if lobe < 1
            push!(invalid_lobes, "$file:$lobe")
            continue
        end
        lobes = get!(by_file, file, Set{Int}())
        lobe in lobes && push!(duplicate_lobes, "$file:$lobe")
        push!(lobes, lobe)
    end

    expected_files = Set(keys(truth))
    observed_files = Set(keys(by_file))
    missing_files = sort(collect(setdiff(expected_files, observed_files)))
    extra_files = sort(collect(setdiff(observed_files, expected_files)))
    incomplete = String[]
    noncontiguous = String[]
    short_n = String[]
    extra_n = String[]
    expected_lobes = sum(length(seq) for seq in values(truth))
    observed_lobes = length(rows)
    for file in sort(collect(intersect(expected_files, observed_files)))
        observed = by_file[file]
        expected_count = length(truth[file])
        if allow_variable_n
            observed_count = isempty(observed) ? 0 : maximum(observed)
            if !isempty(observed) && !(length(observed) == observed_count && all(in(1:observed_count), observed))
                push!(noncontiguous, "$file($(join(sort(collect(observed)), ',')))")
            end
            length(observed) < expected_count && push!(short_n, "$file($(length(observed))/$expected_count)")
            length(observed) > expected_count && push!(extra_n, "$file($(length(observed))/$expected_count)")
        else
            length(observed) == expected_count && all(in(1:expected_count), observed) && continue
            push!(incomplete, "$file($(length(observed))/$expected_count)")
        end
    end

    if !isempty(missing_files) || !isempty(extra_files) || !isempty(incomplete) ||
       !isempty(noncontiguous) || !isempty(duplicate_lobes) || !isempty(invalid_lobes)
        error("Profile $(profile.name) is not full145-ready: expected $(length(expected_files)) files / $expected_lobes lobes, found $(length(observed_files)) files / $observed_lobes lobe rows. Missing files ($(length(missing_files))): $(_format_examples(missing_files)). Extra files ($(length(extra_files))): $(_format_examples(extra_files)). Incomplete files ($(length(incomplete))): $(_format_examples(incomplete)). Noncontiguous files ($(length(noncontiguous))): $(_format_examples(noncontiguous)). Duplicate lobes ($(length(duplicate_lobes))): $(_format_examples(duplicate_lobes)). Invalid lobes ($(length(invalid_lobes))): $(_format_examples(invalid_lobes)).")
    end
    if allow_variable_n && (!isempty(short_n) || !isempty(extra_n))
        println("Profile $(profile.name) own-N audit: short-N files ($(length(short_n))): $(_format_examples(short_n)); extra-N files ($(length(extra_n))): $(_format_examples(extra_n)).")
    end
    return nothing
end

function _write_summary(rows, path::String)
    fields = ["profile", "files", "lobes_possible", "lobes_classified", "coverage",
              "prediction_lobes", "missing_truth_lobes", "extra_prediction_lobes",
              "files_short_N", "files_extra_N",
              "physical_correct", "physical_errors_emitted", "physical_accuracy_classified",
              "honest_correct", "honest_uncertain", "honest_correct_frac", "sequence_exact",
              "oracle_correct", "oracle_accuracy_classified", "prediction_path", "grade_tsv", "note"]
    _ensure_parent(path)
    _reject_symlink_leaf(path)
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            println(io, join([row[f] for f in fields], '\t'))
        end
    end
end

function _write_markdown(rows, truth, path::String, opt::Options)
    _ensure_parent(path)
    _reject_symlink_leaf(path)
    total_lobes = sum(length(seq) for seq in values(truth))
    open(path, "w") do io
        println(io, "# Chitosan 0/1/? Unit-Assignment Benchmark Report\n")
        println(io, "This report is post-hoc benchmark grading only. Truth sequences are read only by the grader/report, never by the prediction rules or fitting path. The listed profiles are frozen operating points; this script does not sweep thresholds or choose a new best profile.\n")
        if opt.full145
            println(io, "Full benchmark denominator: **$(length(truth)) files / $(total_lobes) control positions**. External control sequence: `$(opt.control_sequence)` (`NKNNKN` convention-dependent encoding).\n")
            if opt.full145_own_n
                println(io, "Own-N mode is enabled: each profile may contain the lobes produced by that file's label-free `N_selected`. Missing control positions count as uncertain in the honest view; extra predicted lobes are reported but not aligned to the 6-position external control.\n")
            end
        else
            println(io, "Historical subset denominator: **$(length(truth)) files / $(total_lobes) lobes** (`clean` + `clean_target`; stress rows excluded).\n")
        end
        println(io, "| profile | coverage | N audit | classified accuracy | honest correct | honest uncertain | emitted errors | exact chains |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in rows
            n_audit = "pred $(row["prediction_lobes"]); missing $(row["missing_truth_lobes"]); extra $(row["extra_prediction_lobes"])"
            println(io, "| $(row["profile"]) | $(row["lobes_classified"])/$(row["lobes_possible"]) ($(row["coverage"])) | $n_audit | $(row["physical_correct"])/$(row["lobes_classified"]) ($(row["physical_accuracy_classified"])) | $(row["honest_correct"])/$(row["lobes_possible"]) ($(row["honest_correct_frac"])) | $(row["honest_uncertain"])/$(row["lobes_possible"]) | $(row["physical_errors_emitted"]) | $(row["sequence_exact"])/$(row["files"]) |")
        end
        println(io, "\nCompanion files: `summary.tsv`, `lobe_position_errors.tsv`, and per-profile grade TSVs under `grades/`.")
    end
end

function main(args=ARGS)
    opt = _parse_cli(args)
    mkpath(opt.outdir)
    truth_path = opt.full145 ? _write_full145_truth(opt) : opt.truth
    truth = load_truth(truth_path, opt.primary_only)
    total_lobes = sum(length(seq) for seq in values(truth))
    if opt.full145 && (length(truth) != 145 || total_lobes != 870)
        error("Full145 denominator changed: expected 145 files / 870 lobes, got $(length(truth)) / $total_lobes")
    elseif !opt.full145 && opt.primary_only && truth_path == DEFAULT_TRUTH && (length(truth) != 35 || total_lobes != 210)
        error("Strict default denominator changed: expected 35 files / 210 lobes, got $(length(truth)) / $total_lobes")
    end

    pos_path = joinpath(opt.outdir, "lobe_position_errors.tsv")
    _reject_symlink_leaf(pos_path)
    open(pos_path, "w") do io
        println(io, join(["profile", "lobe", "possible", "classified", "correct", "wrong", "uncertain",
                          "classified_accuracy", "honest_correct_frac"], '\t'))
    end

    summaries = Dict{String,String}[]
    for profile in opt.profiles
        opt.full145 && _preflight_full145_profile(profile, truth; allow_variable_n=opt.full145_own_n)
        grade_path = _run_grader(profile, opt, truth_path)
        push!(summaries, _summarize_grade(profile, grade_path, opt, truth_path))
        write_position_rows(profile.name, profile.path, truth, pos_path;
                            include_missing_truth_positions=opt.full145_own_n)
    end

    _write_summary(summaries, joinpath(opt.outdir, "summary.tsv"))
    _write_markdown(summaries, truth, joinpath(opt.outdir, "report.md"), opt)
    println("\nUnit-assignment benchmark report")
    println("  outdir:     ", opt.outdir)
    scope = opt.full145 ? (opt.full145_own_n ? "full145_own_n" : "full145") : "historical_subset"
    println("  scope:      ", scope)
    println("  truth:      ", truth_path)
    println("  denominator:", length(truth), " files / ", total_lobes, " lobes")
    println("  summary:    ", joinpath(opt.outdir, "summary.tsv"))
    println("  report:     ", joinpath(opt.outdir, "report.md"))
end

main()
