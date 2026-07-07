#!/usr/bin/env julia

# Production wrapper for unknown-chitosan 0/1/? unit-assignment predictions.
#
# This script is deliberately label-free: it accepts selected-N feature artifacts
# and writes fixed prediction profiles plus a production summary. Benchmark truth,
# control motifs, benchmark manifests, and grading outputs are rejected here.

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv

const DEFAULT_PROFILE = "default"
const FORBIDDEN_FLAGS = Set(["--truth", "--control-sequence", "--manifest", "--benchmark-manifest", "--expected-N", "--expected-n", "--full145", "--control"])
const FORBIDDEN_COLUMNS = Set(["sequence", "expected_N", "target_N", "control_sequence"])
const LOCAL_BASE = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel"]
const GAUSSIAN_BASE = ["amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated"]

struct Options
    features::String
    split_features::String
    patches::String
    outdir::String
    profile::String
    seeds::Int
    first_seed::Int
    interactions::Bool
end

struct ProfileRun
    name::String
    out_tsv::String
    views::Vector{Pair{String,Vector{String}}}
    validation_log::String
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

function _forbidden_flag(arg::AbstractString)
    name = split(String(arg), '='; limit=2)[1]
    return name in FORBIDDEN_FLAGS
end

function _show_help()
    println("""
    Usage: julia --project=. test/run_unknown_unit_assignment.jl [options]

    Required: --features PATH --outdir PATH
    Optional: --split-features PATH, --patches PATH, --profile NAME,
              --first-seed INT, --seeds INT, --interactions
    Profiles: base_bwd_consensus, base_split_log_skew, all, default [$(DEFAULT_PROFILE)]
    Outputs: predictions_base_bwd_consensus.tsv,
             predictions_base_split_log_skew.tsv, summary.tsv

    Label-free/no-truth constraint: rejects benchmark-only flags/columns such as
    --truth, --manifest, --expected-N, --full145, sequence, expected_N,
    target_N, and control_sequence. Grade only after outputs are frozen.
    """)
end

function _parse_cli(args)
    features = ""
    split_features = ""
    patches = ""
    outdir = ""
    profile = DEFAULT_PROFILE
    seeds = 20
    first_seed = 0
    interactions = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if _forbidden_flag(arg)
            error("Forbidden benchmark-only flag in unknown production runner: $(split(arg, '='; limit=2)[1])")
        elseif arg == "--features"
            features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--features=")
            features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--split-features"
            split_features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--split-features=")
            split_features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--patches"
            patches = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--patches=")
            patches = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--outdir"
            outdir = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--profile"
            profile = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--profile=")
            profile = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--first-seed"
            first_seed = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--first-seed=")
            first_seed = parse(Int, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--seeds"
            seeds = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--seeds=")
            seeds = parse(Int, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--interactions"
            interactions = true; i += 1
        elseif arg in ("-h", "--help")
            _show_help()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isempty(features) && error("--features is required")
    isempty(outdir) && error("--outdir is required")
    profile in ("base_bwd_consensus", "base_split_log_skew", "all", "default") ||
        error("--profile must be one of base_bwd_consensus, base_split_log_skew, all, default")
    seeds > 0 || error("--seeds must be positive")
    isfile(features) || error("Feature TSV not found: $features")
    !isempty(split_features) && !isfile(split_features) && error("Split feature TSV not found: $split_features")
    !isempty(patches) && !isfile(patches) && error("Patch TSV not found: $patches")
    return Options(features, split_features, patches, outdir, profile, seeds, first_seed, interactions)
end

function _check_forbidden_columns(path::String, label::String)
    header, _ = _read_tsv(path)
    hits = sort(collect(intersect(Set(String.(header)), FORBIDDEN_COLUMNS)))
    isempty(hits) || error("$label contains benchmark-control columns forbidden in unknown production: $(join(hits, ", "))")
    return String.(header)
end

function _base_features(header::Vector{String})
    available = Set(header)
    if all(in(available), LOCAL_BASE)
        return LOCAL_BASE
    elseif all(in(available), GAUSSIAN_BASE)
        return GAUSSIAN_BASE
    end
    error("Feature TSV has neither local base features ($(join(LOCAL_BASE, ", "))) nor gaussian base features ($(join(GAUSSIAN_BASE, ", ")))")
end

function _profile_runs(opt::Options, base::Vector{String})
    runs = ProfileRun[]
    selected = opt.profile in ("default", "all") ? ["base_bwd_consensus", "base_split_log_skew"] : [opt.profile]
    for name in selected
        if name == "base_bwd_consensus"
            isempty(opt.patches) && error("Profile base_bwd_consensus requires --patches")
            views = [
                "base" => base,
                "base_bwd_neg_com_t" => vcat(base, ["bwd_neg_com_t"]),
                "base_bwd_neg_diag45" => vcat(base, ["bwd_neg_diag45"]),
            ]
            push!(runs, ProfileRun(name, joinpath(opt.outdir, "predictions_base_bwd_consensus.tsv"), views, joinpath(opt.outdir, "validation_base_bwd_consensus.log")))
        elseif name == "base_split_log_skew"
            isempty(opt.split_features) && error("Profile base_split_log_skew requires --split-features")
            views = ["base_split_log_skew" => vcat(base, ["split_log_skew"])]
            push!(runs, ProfileRun(name, joinpath(opt.outdir, "predictions_base_split_log_skew.tsv"), views, joinpath(opt.outdir, "validation_base_split_log_skew.log")))
        else
            error("Unknown fixed profile: $name")
        end
    end
    return runs
end

function _prediction_command(opt::Options, run_spec::ProfileRun)
    script = joinpath(@__DIR__, "build_labelfree_unit_predictions.jl")
    project_dir = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --project=$project_dir $script --features $(opt.features) --out $(run_spec.out_tsv) --first-seed $(opt.first_seed) --seeds $(opt.seeds)`
    if !isempty(opt.split_features)
        cmd = `$cmd --split-features $(opt.split_features)`
    end
    if !isempty(opt.patches)
        cmd = `$cmd --patches $(opt.patches)`
    end
    opt.interactions && (cmd = `$cmd --interactions`)
    for view in run_spec.views
        cmd = `$cmd --view $(first(view))=$(join(last(view), ','))`
    end
    return cmd
end

function _validation_command(opt::Options, run_spec::ProfileRun)
    script = joinpath(@__DIR__, "validate_unit_predictions.jl")
    project_dir = dirname(Base.active_project())
    return `$(Base.julia_cmd()) --project=$project_dir $script --predictions $(run_spec.out_tsv) --features $(opt.features)`
end

function _read_prediction_summary(path::String)
    _, rows = _read_tsv(path)
    files = Set{String}()
    predicted_0 = 0
    predicted_1 = 0
    uncertain = 0
    conf_sum = 0.0
    conf_count = 0
    for row in rows
        haskey(row, "file") || error("Prediction TSV missing file column: $path")
        haskey(row, "predicted") || error("Prediction TSV missing predicted column: $path")
        push!(files, basename(strip(row["file"])))
        pred = strip(row["predicted"])
        if pred == "0"
            predicted_0 += 1
        elseif pred == "1"
            predicted_1 += 1
        elseif pred == "?"
            uncertain += 1
        else
            error("Invalid predicted value in $path: $pred")
        end
        if haskey(row, "confidence")
            conf = tryparse(Float64, strip(row["confidence"]))
            if conf !== nothing && isfinite(conf)
                conf_sum += conf
                conf_count += 1
            end
        end
    end
    mean_confidence = conf_count > 0 ? conf_sum / conf_count : NaN
    return (files=length(files), lobes=length(rows), predicted_0=predicted_0,
            predicted_1=predicted_1, uncertain=uncertain,
            mean_confidence=mean_confidence)
end

function _write_summary(path::String, opt::Options, runs::Vector{ProfileRun})
    fields = ["profile", "prediction_path", "files", "lobes", "predicted_0", "predicted_1", "uncertain", "mean_confidence", "features", "split_features", "patches", "seeds", "first_seed", "interactions"]
    _ensure_parent(path)
    islink(path) && error("Refusing to overwrite symlink: $path")
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for run_spec in runs
            s = _read_prediction_summary(run_spec.out_tsv)
            row = [
                run_spec.name,
                run_spec.out_tsv,
                string(s.files),
                string(s.lobes),
                string(s.predicted_0),
                string(s.predicted_1),
                string(s.uncertain),
                isfinite(s.mean_confidence) ? string(round(s.mean_confidence; digits=8)) : "NA",
                opt.features,
                isempty(opt.split_features) ? "none" : opt.split_features,
                isempty(opt.patches) ? "none" : opt.patches,
                string(opt.seeds),
                string(opt.first_seed),
                string(opt.interactions),
            ]
            println(io, join(row, '\t'))
        end
    end
    return nothing
end

function _write_manifest(path::String, opt::Options, runs::Vector{ProfileRun})
    fields = ["profile", "prediction_path", "validation_log", "validation_status", "files", "lobes", "input_path_features", "input_path_split_features", "input_path_patches", "split_available", "patches_available", "seeds", "first_seed", "interactions", "views"]
    _ensure_parent(path)
    islink(path) && error("Refusing to overwrite symlink: $path")
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for run_spec in runs
            s = _read_prediction_summary(run_spec.out_tsv)
            row = [run_spec.name, run_spec.out_tsv, run_spec.validation_log, "ok", string(s.files), string(s.lobes), opt.features, isempty(opt.split_features) ? "none" : opt.split_features, isempty(opt.patches) ? "none" : opt.patches, string(!isempty(opt.split_features)), string(!isempty(opt.patches)), string(opt.seeds), string(opt.first_seed), string(opt.interactions), join(first.(run_spec.views), ",")]
            println(io, join(row, '\t'))
        end
    end
    return nothing
end

function main(args=ARGS)
    opt = _parse_cli(args)
    feature_header = _check_forbidden_columns(opt.features, "Feature TSV")
    !isempty(opt.split_features) && _check_forbidden_columns(opt.split_features, "Split feature TSV")
    !isempty(opt.patches) && _check_forbidden_columns(opt.patches, "Patch TSV")
    base = _base_features(feature_header)
    runs = _profile_runs(opt, base)

    mkpath(opt.outdir)
    for run_spec in runs
        cmd = _prediction_command(opt, run_spec)
        println("running profile=", run_spec.name)
        println("  out=", run_spec.out_tsv)
        run(cmd)
        _ensure_parent(run_spec.validation_log)
        islink(run_spec.validation_log) && error("Refusing to overwrite symlink: $(run_spec.validation_log)")
        open(run_spec.validation_log, "w") do io
            run(pipeline(_validation_command(opt, run_spec), stdout=io, stderr=io))
        end
    end
    summary_path = joinpath(opt.outdir, "summary.tsv")
    manifest_path = joinpath(opt.outdir, "manifest.tsv")
    _write_summary(summary_path, opt, runs)
    _write_manifest(manifest_path, opt, runs)

    println("\nUnknown unit-assignment production run")
    println("  features: ", opt.features)
    println("  split:    ", isempty(opt.split_features) ? "none" : opt.split_features)
    println("  patches:  ", isempty(opt.patches) ? "none" : opt.patches)
    println("  outdir:   ", opt.outdir)
    println("  profiles: ", join([r.name for r in runs], ", "))
    println("  summary:  ", summary_path)
    println("  manifest: ", manifest_path)
    return nothing
end

main()
