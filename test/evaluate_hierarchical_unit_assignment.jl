#!/usr/bin/env julia

module HierarchicalUnitAssignmentEvaluation

using Dates
using Printf
using Random
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl"))
using .HierarchicalUnitAssignment

export EvaluationError, DateFold, BOOTSTRAP_SEEDS, SHARD_COLUMNS,
       parse_leading_date, build_date_folds, whole_scan_bootstrap,
       bootstrap_training_units, evaluate_fold, run_cli

struct EvaluationError <: Exception
    msg::String
end
Base.showerror(io::IO, err::EvaluationError) = print(io, "EvaluationError: ", err.msg)
_fail(msg) = throw(EvaluationError(String(msg)))

const BOOTSTRAP_SEEDS = collect(0:499)
const RESAMPLE_UNIT = "whole_scan_within_training_dates"
const SHARD_COLUMNS = String[
    "schema_version", "fold_date", "bootstrap_kind", "bootstrap_seed",
    "bootstrap_block", "scan_key", "scan_date", "training_scan_count",
    "heldout_scan_count", "lobe_count", "ll_one", "ll_two", "improvement",
    "improvement_per_lobe", "forward_backward_agreement",
    "perturbation_min_agreement", "occupancy_1", "config_sha256",
    "input_sha256", "source_sha256", "resample_unit", "row_complete",
    "shard_complete", "expected_rows", "expected_seed_count",
]

struct DateFold
    fold_date::String
    training_scans::Vector{String}
    heldout_scans::Vector{String}
end
Base.:(==)(a::DateFold, b::DateFold) =
    a.fold_date == b.fold_date && a.training_scans == b.training_scans &&
    a.heldout_scans == b.heldout_scans

function parse_leading_date(path::AbstractString)::String
    normalized = replace(strip(String(path)), '\\' => '/')
    isempty(normalized) && _fail("scan path is empty")
    tokens = String[]
    for component in split(normalized, '/'; keepempty=false)
        m = match(r"^(\d{8})(?:$|[_.-])", component)
        m === nothing || push!(tokens, m.captures[1])
    end
    length(tokens) == 1 || _fail(
        "scan path must contain exactly one leading YYYYMMDD token; found $(length(tokens)): $path")
    token = only(tokens)
    try
        Date(token, dateformat"yyyymmdd")
    catch
        _fail("invalid leading YYYYMMDD token $token in $path")
    end
    return token
end

function build_date_folds(scan_paths)::Vector{DateFold}
    paths = sort!(unique(String.(scan_paths)))
    isempty(paths) && _fail("no scan paths supplied")
    dates = Dict(path => parse_leading_date(path) for path in paths)
    unique_dates = sort!(unique(collect(values(dates))))
    length(unique_dates) >= 2 || _fail("leave-date-out evaluation needs at least two dates")
    return [DateFold(date,
                     [p for p in paths if dates[p] != date],
                     [p for p in paths if dates[p] == date])
            for date in unique_dates]
end

function whole_scan_bootstrap(training_scans, seed::Integer)::Vector{String}
    scans = sort!(unique(String.(training_scans)))
    isempty(scans) && _fail("cannot bootstrap an empty training set")
    by_date = Dict{String,Vector{String}}()
    for scan in scans
        push!(get!(by_date, parse_leading_date(scan), String[]), scan)
    end
    rng = MersenneTwister(Int(seed))
    sampled = String[]
    for date in sort!(collect(keys(by_date)))
        group = sort!(by_date[date])
        append!(sampled, [group[rand(rng, eachindex(group))] for _ in eachindex(group)])
    end
    return sampled
end

function bootstrap_training_units(training_scans, seed::Integer;
                                  resample_unit::AbstractString=RESAMPLE_UNIT)
    resample_unit == RESAMPLE_UNIT ||
        _fail("per-lobe bootstrap is forbidden; resample_unit must be $RESAMPLE_UNIT")
    seed in BOOTSTRAP_SEEDS || _fail("bootstrap seed must be in 0:499, got $seed")
    return whole_scan_bootstrap(training_scans, seed)
end

_sha256_file(path) = bytes2hex(SHA.sha256(read(path)))

function _source_sha256()
    paths = [abspath(@__FILE__), joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl")]
    append!(paths, sort!(filter(isfile, readdir(joinpath(@__DIR__, "lib", "hierarchical"); join=true))))
    io = IOBuffer()
    for path in paths
        println(io, relpath(path, dirname(@__DIR__)), '\t', _sha256_file(path))
    end
    return bytes2hex(SHA.sha256(take!(io)))
end

function _raw_scan_paths(features_path::AbstractString)
    header, rows = HierarchicalUnitAssignment._read_tsv(features_path)
    "file" in header || _fail("feature TSV lacks file column: $features_path")
    paths = sort!(unique(strip(row["file"]) for row in rows))
    any(isempty, paths) && _fail("feature TSV contains an empty scan path")
    by_base = Dict{String,String}()
    for path in paths
        base = basename(path)
        haskey(by_base, base) && by_base[base] != path &&
            _fail("duplicate basename maps to multiple relative scan paths: $base")
        by_base[base] = path
        parse_leading_date(path)
    end
    return paths, by_base
end

function _default_views(records)
    available = Set{String}()
    for rec in records
        union!(available, keys(rec.features))
    end
    base_local = ["amp_prominence", "amp_neighbor_ratio",
                  "integrated_prominence", "amp_rel"]
    base_gaussian = ["amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated"]
    base = all(in(available), base_local) ? base_local :
           all(in(available), base_gaussian) ? base_gaussian :
           _fail("no default hierarchical feature base is available")
    specs = Pair{String,Vector{String}}["base_local" => copy(base)]
    for (name, extra) in (("base_local+bwd_neg_com_t", "bwd_neg_com_t"),
                          ("base_local+bwd_neg_diag45", "bwd_neg_diag45"),
                          ("base_local+split_log_skew", "split_log_skew"))
        extra in available && push!(specs, name => vcat(base, [extra]))
    end
    return specs
end

function _parse_view(value::AbstractString)
    parts = split(String(value), '='; limit=2)
    length(parts) == 2 || _fail("--view requires NAME=feature1,feature2,...")
    name = strip(parts[1])
    features = [strip(x) for x in split(parts[2], ',') if !isempty(strip(x))]
    isempty(name) && _fail("--view name is empty")
    isempty(features) && _fail("--view $name has no features")
    any(HierarchicalUnitAssignment.is_forbidden_feature, features) &&
        _fail("--view $name contains a forbidden feature")
    return name => features
end

function _one_ll(fit, Z)
    return HierarchicalUnitAssignment._total_log_likelihood(
        Z, reshape(fit.mean, 1, :), reshape(fit.var, 1, :), [1.0])
end

function _mapped_labels(fit, high_cluster::Int, Z)
    _require_stable_fit(fit, "held-out prediction")
    resp = HierarchicalUnitAssignment.responsibilities(fit, Z)
    return [argmax(view(resp, i, :)) == high_cluster ? 1 : 0 for i in 1:size(Z, 1)]
end

function _require_stable_fit(fit, context::AbstractString)
    fit.converged || _fail("$context produced an unstable_model fit")
    fit.monotone || _fail("$context produced a nonmonotone_model fit")
    return fit
end

function _fit_emission(Ztrain, train_records; first_seed=0)
    fit2 = HierarchicalUnitAssignment.fit_em_two_component(Ztrain, first_seed)
    _require_stable_fit(fit2, "training emission")
    fit1 = HierarchicalUnitAssignment.fit_one_component(Ztrain)
    resp = HierarchicalUnitAssignment.responsibilities(fit2, Ztrain)
    hard = [argmax(view(resp, i, :)) for i in 1:size(resp, 1)]
    high = HierarchicalUnitAssignment.physical_high_amplitude_cluster(train_records, hard)
    return fit1, fit2, high
end

_agreement(a, b) = isempty(a) || length(a) != length(b) ? NaN :
                    count(identity, a .== b) / length(a)

function _validate_candidate_config(path::AbstractString)
    cfg = TOML.parsefile(path)
    get(cfg["bootstrap"], "count", nothing) == 500 || _fail("config bootstrap count must be 500")
    get(cfg["bootstrap"], "seeds_range_start", nothing) == 0 || _fail("config seed start must be 0")
    get(cfg["bootstrap"], "seeds_range_stop", nothing) == 499 || _fail("config seed stop must be 499")
    get(cfg["bootstrap"], "resample_unit", nothing) == RESAMPLE_UNIT || _fail("config resample unit mismatch")
    get(cfg["features"], "policy", nothing) == "equal_priors" || _fail("config must retain equal priors")
    return cfg
end

function evaluate_fold(features_path::AbstractString, config_path::AbstractString,
                       fold_date::AbstractString, block::AbstractString;
                       view_specs=Pair{String,Vector{String}}[])
    _validate_candidate_config(config_path)
    records = HierarchicalUnitAssignment.load_records(features_path)
    scan_paths, by_base = _raw_scan_paths(features_path)
    folds = build_date_folds(scan_paths)
    fold_idx = findfirst(f -> f.fold_date == fold_date, folds)
    fold_idx === nothing && _fail("unknown fold date $fold_date")
    fold = folds[fold_idx]
    specs = isempty(view_specs) ? _default_views(records) : view_specs
    base_to_path = Dict(base => path for (base, path) in by_base)
    path_to_indices = Dict(path => findall(r -> base_to_path[r.file] == path, records)
                           for path in scan_paths)
    runs = if block == "unbootstrapped"
        [(kind="unbootstrapped", seed=-1, scans=fold.training_scans)]
    else
        m = match(r"^(\d+):(\d+)$", String(block))
        m === nothing && _fail("bootstrap block must be unbootstrapped or START:STOP")
        first_seed, last_seed = parse.(Int, m.captures)
        first_seed <= last_seed || _fail("bootstrap block start exceeds stop")
        all(in(BOOTSTRAP_SEEDS), first_seed:last_seed) || _fail("bootstrap block must stay within 0:499")
        [(kind="bootstrap", seed=seed,
          scans=bootstrap_training_units(fold.training_scans, seed))
         for seed in first_seed:last_seed]
    end

    config_hash = _sha256_file(config_path)
    input_hash = _sha256_file(features_path)
    source_hash = _source_sha256()
    output = Dict{String,Any}[]
    for run in runs
        per_scan = Dict(scan => (ll1=Float64[], ll2=Float64[], labels=Vector{Vector{Int}}(),
                                 names=String[], perturb=Float64[]) for scan in fold.heldout_scans)
        for (view_name, features) in specs
            X, _ = HierarchicalUnitAssignment.feature_matrix(records, features)
            Z = HierarchicalUnitAssignment.normalize_per_scan(records, X, features)
            train_idx = reduce(vcat, (path_to_indices[path] for path in run.scans); init=Int[])
            train_idx = [i for i in train_idx if all(isfinite, view(Z, i, :))]
            length(train_idx) >= 2 || _fail("fold $fold_date view $view_name has fewer than two finite training lobes")
            Ztrain = Z[train_idx, :]
            fit1, fit2, high = _fit_emission(Ztrain, records[train_idx])
            _, perturbed, perturbed_high = _fit_emission(Ztrain, records[train_idx]; first_seed=1)
            for scan in fold.heldout_scans
                idx = [i for i in path_to_indices[scan] if all(isfinite, view(Z, i, :))]
                isempty(idx) && continue
                Zheld = Z[idx, :]
                labels = _mapped_labels(fit2, high, Zheld)
                alt_labels = _mapped_labels(perturbed, perturbed_high, Zheld)
                acc = per_scan[scan]
                push!(acc.ll1, _one_ll(fit1, Zheld) / length(idx))
                push!(acc.ll2, HierarchicalUnitAssignment.log_likelihood_under(fit2, Zheld) / length(idx))
                push!(acc.labels, labels)
                push!(acc.names, view_name)
                push!(acc.perturb, _agreement(labels, alt_labels))
            end
        end
        for scan in fold.heldout_scans
            acc = per_scan[scan]
            isempty(acc.ll1) && _fail("fold $fold_date held-out scan $scan has no finite view")
            n_lobes = length(path_to_indices[scan])
            ll1 = mean(acc.ll1) * n_lobes
            ll2 = mean(acc.ll2) * n_lobes
            base_i = findfirst(name -> !occursin("bwd", lowercase(name)), acc.names)
            bwd_i = findall(name -> occursin("bwd", lowercase(name)), acc.names)
            fb = base_i === nothing || isempty(bwd_i) ? NaN :
                 mean(_agreement(acc.labels[base_i], acc.labels[i]) for i in bwd_i)
            base_labels = base_i === nothing ? first(acc.labels) : acc.labels[base_i]
            occupancy = count(==(1), base_labels) / length(base_labels)
            improvement = ll2 - ll1
            push!(output, Dict(
                "schema_version" => 1, "fold_date" => fold.fold_date,
                "bootstrap_kind" => run.kind, "bootstrap_seed" => run.seed,
                "bootstrap_block" => String(block), "scan_key" => scan,
                "scan_date" => parse_leading_date(scan),
                "training_scan_count" => length(fold.training_scans),
                "heldout_scan_count" => length(fold.heldout_scans),
                "lobe_count" => n_lobes, "ll_one" => ll1, "ll_two" => ll2,
                "improvement" => improvement,
                "improvement_per_lobe" => improvement / n_lobes,
                "forward_backward_agreement" => fb,
                "perturbation_min_agreement" => minimum(acc.perturb),
                "occupancy_1" => occupancy, "config_sha256" => config_hash,
                "input_sha256" => input_hash, "source_sha256" => source_hash,
                "resample_unit" => RESAMPLE_UNIT, "row_complete" => true,
                "shard_complete" => true, "expected_rows" => 0,
                "expected_seed_count" => 500))
        end
    end
    for row in output
        row["expected_rows"] = length(output)
    end
    sort!(output; by=row -> (row["fold_date"], row["bootstrap_seed"], row["scan_key"]))
    return output
end

_format_value(x::AbstractFloat) = isfinite(x) ? @sprintf("%.17g", x) : "NA"
_format_value(x) = string(x)

function _atomic_write(writer::Function, path::AbstractString)
    destination = abspath(path)
    parent = dirname(destination)
    mkpath(parent)
    temp_path, io = mktemp(parent; cleanup=false)
    committed = false
    try
        writer(io)
        flush(io)
        close(io)
        Base.Filesystem.rename(temp_path, destination)
        committed = true
    finally
        isopen(io) && close(io)
        !committed && ispath(temp_path) && rm(temp_path; force=true)
    end
    return nothing
end

function _write_rows(path::AbstractString, rows)
    _atomic_write(path) do io
        println(io, join(SHARD_COLUMNS, '\t'))
        for row in rows
            println(io, join((_format_value(row[col]) for col in SHARD_COLUMNS), '\t'))
        end
    end
end

function _show_help()
    println("""
Usage: julia --project=. test/evaluate_hierarchical_unit_assignment.jl \\
  --features PATH --config PATH --fold YYYYMMDD \\
  --bootstrap-block unbootstrapped|START:STOP --out PATH [--view NAME=LIST]

Deterministic label-free leave-one-date-out evaluator. Bootstrap blocks use
exact seeds 0:499 and resample whole scans within each training date. The
held-out date is scored only after fitting on training-date scans. Optional
--config-sha256, --input-sha256, and --source-sha256 values are verified before
writing a strict shard.
""")
end

function _value(args, i, flag)
    i < length(args) || _fail("$flag requires a value")
    return String(args[i + 1])
end

function run_cli(args=ARGS)
    any(a -> a in ("-h", "--help"), args) && (_show_help(); return nothing)
    values = Dict("--config" => "config/unit_assignment_candidate.toml")
    views = Pair{String,Vector{String}}[]
    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--view"
            push!(views, _parse_view(_value(args, i, arg))); i += 2
        elseif startswith(arg, "--view=")
            push!(views, _parse_view(split(arg, '='; limit=2)[2])); i += 1
        elseif arg in ("--features", "--config", "--fold", "--bootstrap-block", "--out",
                       "--config-sha256", "--input-sha256", "--source-sha256")
            values[arg] = _value(args, i, arg); i += 2
        elseif startswith(arg, "--") && occursin('=', arg)
            key, value = split(arg, '='; limit=2)
            key in ("--features", "--config", "--fold", "--bootstrap-block", "--out",
                    "--config-sha256", "--input-sha256", "--source-sha256") ||
                _fail("unknown argument: $arg")
            values[key] = value; i += 1
        else
            _fail("unknown argument: $arg")
        end
    end
    for key in ("--features", "--fold", "--bootstrap-block", "--out")
        haskey(values, key) || _fail("$key is required")
    end
    isfile(values["--features"]) || _fail("feature TSV not found: $(values["--features"])")
    isfile(values["--config"]) || _fail("config not found: $(values["--config"])")
    rows = evaluate_fold(values["--features"], values["--config"], values["--fold"],
                         values["--bootstrap-block"]; view_specs=views)
    observed = Dict("--config-sha256" => rows[1]["config_sha256"],
                    "--input-sha256" => rows[1]["input_sha256"],
                    "--source-sha256" => rows[1]["source_sha256"])
    for (flag, actual) in observed
        haskey(values, flag) && lowercase(values[flag]) != actual &&
            _fail("$flag mismatch: expected $(values[flag]), observed $actual")
    end
    _write_rows(values["--out"], rows)
    println("hierarchical evaluation shard complete")
    println("  fold: ", values["--fold"], "  block: ", values["--bootstrap-block"])
    println("  rows: ", length(rows), "  out: ", values["--out"])
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        run_cli(ARGS)
    catch err
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end

end
