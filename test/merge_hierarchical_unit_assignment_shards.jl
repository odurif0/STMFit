#!/usr/bin/env julia

module HierarchicalUnitAssignmentMerge

using Printf
using SHA
using Statistics

export MergeError, MergeContract, GateResult, MergeTransactionInterruption, SHARD_COLUMNS,
       evaluate_gate, merge_shards, run_cli

struct MergeError <: Exception
    msg::String
end
Base.showerror(io::IO, err::MergeError) = print(io, "MergeError: ", err.msg)
_fail(msg) = throw(MergeError(String(msg)))

const SHARD_COLUMNS = String[
    "schema_version", "fold_date", "bootstrap_kind", "bootstrap_seed",
    "bootstrap_block", "scan_key", "scan_date", "training_scan_count",
    "heldout_scan_count", "lobe_count", "ll_one", "ll_two", "improvement",
    "improvement_per_lobe", "forward_backward_agreement",
    "perturbation_min_agreement", "occupancy_1", "config_sha256",
    "input_sha256", "source_sha256", "resample_unit", "row_complete",
    "shard_complete", "expected_rows", "expected_seed_count",
]
const EXPECTED_SEEDS = Set(0:499)
const RESAMPLE_UNIT = "whole_scan_within_training_dates"

struct MergeContract
    folds::Vector{String}
    config_sha256::String
    input_sha256::String
    source_sha256::String
    function MergeContract(folds, config_sha256, input_sha256, source_sha256)
        normalized_folds = sort!(unique(String.(folds)))
        isempty(normalized_folds) && _fail("at least one expected fold is required")
        hashes = lowercase.([String(config_sha256), String(input_sha256), String(source_sha256)])
        all(h -> occursin(r"^[0-9a-f]{64}$", h), hashes) ||
            _fail("config/input/source hashes must be 64 lowercase hexadecimal characters")
        new(normalized_folds, hashes...)
    end
end

struct GateResult
    passed::Bool
    every_fold_positive::Bool
    lower_bound::Float64
    fold_improvements::Vector{Pair{String,Float64}}
    bootstrap_seed_means::Vector{Float64}
end
Base.:(==)(a::GateResult, b::GateResult) =
    a.passed == b.passed && a.every_fold_positive == b.every_fold_positive &&
    a.lower_bound == b.lower_bound && a.fold_improvements == b.fold_improvements &&
    a.bootstrap_seed_means == b.bootstrap_seed_means

function _parse_bool(value, field, path)
    text = lowercase(strip(String(value)))
    text == "true" && return true
    text == "false" && return false
    _fail("$path has invalid $field boolean: $value")
end

function _parse_int(value, field, path)
    parsed = tryparse(Int, strip(String(value)))
    parsed === nothing && _fail("$path has invalid $field integer: $value")
    return parsed
end

function _parse_float(value, field, path; allow_na=false)
    text = strip(String(value))
    allow_na && text == "NA" && return NaN
    parsed = tryparse(Float64, text)
    parsed === nothing && _fail("$path has invalid $field number: $value")
    isfinite(parsed) || _fail("$path has non-finite $field: $value")
    return parsed
end

function _read_shard(path::AbstractString)
    isfile(path) || _fail("shard not found: $path")
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    isempty(lines) && _fail("empty shard: $path")
    header = String.(split(lines[1], '\t'; keepempty=true))
    header == SHARD_COLUMNS || _fail("malformed shard schema in $path")
    rows = Dict{String,Any}[]
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        values = split(line, '\t'; keepempty=true)
        length(values) == length(header) ||
            _fail("$path:$line_number has $(length(values)) fields, expected $(length(header))")
        raw = Dict(header[i] => values[i] for i in eachindex(header))
        row = Dict{String,Any}(
            "schema_version" => _parse_int(raw["schema_version"], "schema_version", path),
            "fold_date" => raw["fold_date"], "bootstrap_kind" => raw["bootstrap_kind"],
            "bootstrap_seed" => _parse_int(raw["bootstrap_seed"], "bootstrap_seed", path),
            "bootstrap_block" => raw["bootstrap_block"], "scan_key" => raw["scan_key"],
            "scan_date" => raw["scan_date"],
            "training_scan_count" => _parse_int(raw["training_scan_count"], "training_scan_count", path),
            "heldout_scan_count" => _parse_int(raw["heldout_scan_count"], "heldout_scan_count", path),
            "lobe_count" => _parse_int(raw["lobe_count"], "lobe_count", path),
            "ll_one" => _parse_float(raw["ll_one"], "ll_one", path),
            "ll_two" => _parse_float(raw["ll_two"], "ll_two", path),
            "improvement" => _parse_float(raw["improvement"], "improvement", path),
            "improvement_per_lobe" => _parse_float(raw["improvement_per_lobe"], "improvement_per_lobe", path),
            "forward_backward_agreement" => _parse_float(raw["forward_backward_agreement"], "forward_backward_agreement", path; allow_na=true),
            "perturbation_min_agreement" => _parse_float(raw["perturbation_min_agreement"], "perturbation_min_agreement", path),
            "occupancy_1" => _parse_float(raw["occupancy_1"], "occupancy_1", path),
            "config_sha256" => lowercase(raw["config_sha256"]),
            "input_sha256" => lowercase(raw["input_sha256"]),
            "source_sha256" => lowercase(raw["source_sha256"]),
            "resample_unit" => raw["resample_unit"],
            "row_complete" => _parse_bool(raw["row_complete"], "row_complete", path),
            "shard_complete" => _parse_bool(raw["shard_complete"], "shard_complete", path),
            "expected_rows" => _parse_int(raw["expected_rows"], "expected_rows", path),
            "expected_seed_count" => _parse_int(raw["expected_seed_count"], "expected_seed_count", path),
        )
        push!(rows, row)
    end
    isempty(rows) && _fail("shard has no data rows: $path")
    expected_rows = unique(row["expected_rows"] for row in rows)
    length(expected_rows) == 1 && only(expected_rows) == length(rows) ||
        _fail("shard completeness mismatch in $path")
    all(row -> row["row_complete"] && row["shard_complete"], rows) ||
        _fail("partial/incomplete shard: $path")
    return rows
end

function _validate_row(row, contract::MergeContract)
    row["schema_version"] == 1 || _fail("unsupported schema version")
    row["fold_date"] in contract.folds || _fail("wrong fold: $(row["fold_date"])")
    row["scan_date"] == row["fold_date"] || _fail("held-out scan date does not match fold")
    row["config_sha256"] == contract.config_sha256 || _fail("wrong config hash")
    row["input_sha256"] == contract.input_sha256 || _fail("wrong input hash")
    row["source_sha256"] == contract.source_sha256 || _fail("wrong source hash")
    row["resample_unit"] == RESAMPLE_UNIT || _fail("per-lobe/unknown bootstrap semantics rejected")
    row["expected_seed_count"] == 500 || _fail("expected seed count must be exactly 500")
    row["training_scan_count"] > 0 || _fail("training scan count must be positive")
    row["heldout_scan_count"] > 0 || _fail("held-out scan count must be positive")
    row["lobe_count"] > 0 || _fail("lobe count must be positive")
    isapprox(row["improvement"], row["ll_two"] - row["ll_one"]; rtol=1e-10, atol=1e-10) ||
        _fail("improvement does not equal ll_two - ll_one")
    isapprox(row["improvement_per_lobe"], row["improvement"] / row["lobe_count"];
             rtol=1e-10, atol=1e-10) || _fail("per-lobe improvement mismatch")
    0.0 <= row["occupancy_1"] <= 1.0 || _fail("descriptive occupancy is outside [0,1]")
    if row["bootstrap_kind"] == "unbootstrapped"
        row["bootstrap_seed"] == -1 || _fail("unbootstrapped row has a bootstrap seed")
        row["bootstrap_block"] == "unbootstrapped" || _fail("unbootstrapped row has wrong block")
    elseif row["bootstrap_kind"] == "bootstrap"
        seed = row["bootstrap_seed"]
        seed in EXPECTED_SEEDS || _fail("wrong bootstrap seed: $seed")
        m = match(r"^(\d+):(\d+)$", row["bootstrap_block"])
        m === nothing && _fail("malformed bootstrap block")
        lo, hi = parse.(Int, m.captures)
        lo <= seed <= hi || _fail("seed $seed is outside declared bootstrap block")
    else
        _fail("unknown bootstrap kind: $(row["bootstrap_kind"])")
    end
end

function _validate_complete(rows, contract::MergeContract)
    seen = Set{Tuple{String,String,Int,String}}()
    for row in rows
        _validate_row(row, contract)
        key = (row["fold_date"], row["bootstrap_kind"], row["bootstrap_seed"], row["scan_key"])
        key in seen && _fail("duplicate fold/kind/seed/scan row: $key")
        push!(seen, key)
    end
    for fold in contract.folds
        unboot = [row for row in rows if row["fold_date"] == fold && row["bootstrap_kind"] == "unbootstrapped"]
        isempty(unboot) && _fail("missing unbootstrapped fold $fold")
        scans = sort!(unique(String(row["scan_key"]) for row in unboot))
        length(scans) == length(unboot) || _fail("duplicate unbootstrapped scan in fold $fold")
        all(row -> row["heldout_scan_count"] == length(scans), unboot) ||
            _fail("wrong held-out scan count in fold $fold")
        boot = [row for row in rows if row["fold_date"] == fold && row["bootstrap_kind"] == "bootstrap"]
        seeds = Set(Int(row["bootstrap_seed"]) for row in boot)
        seeds == EXPECTED_SEEDS || _fail("partial or wrong seed set in fold $fold")
        for seed in 0:499
            seed_scans = sort!(String[row["scan_key"] for row in boot if row["bootstrap_seed"] == seed])
            seed_scans == scans || _fail("missing/wrong held-out scans for fold $fold seed $seed")
        end
    end
    return nothing
end

function evaluate_gate(rows, folds)::GateResult
    expected_folds = sort!(unique(String.(folds)))
    fold_values = Pair{String,Float64}[]
    for fold in expected_folds
        values = Float64[row["improvement_per_lobe"] for row in rows
                         if row["fold_date"] == fold && row["bootstrap_kind"] == "unbootstrapped"]
        isempty(values) && _fail("missing unbootstrapped rows for gate fold $fold")
        push!(fold_values, fold => mean(values))
    end
    every_positive = all(pair -> last(pair) > 0.0, fold_values)
    seed_means = Float64[]
    for seed in 0:499
        values = Float64[row["improvement_per_lobe"] for row in rows
                         if row["bootstrap_kind"] == "bootstrap" && row["bootstrap_seed"] == seed]
        isempty(values) && _fail("missing bootstrap seed $seed for gate")
        push!(seed_means, mean(values))
    end
    lower_bound = quantile(seed_means, 0.025)
    return GateResult(every_positive && lower_bound > 0.0, every_positive,
                      lower_bound, fold_values, seed_means)
end

_format(x::AbstractFloat) = isfinite(x) ? @sprintf("%.17g", x) : "NA"
_format(x) = string(x)

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

function _write_rows(path, rows)
    _atomic_write(path) do io
        println(io, join(SHARD_COLUMNS, '\t'))
        for row in rows
            println(io, join((_format(row[col]) for col in SHARD_COLUMNS), '\t'))
        end
    end
end

function _write_gate(path, gate::GateResult)
    _atomic_write(path) do io
        println(io, "terminal_status\tevery_fold_positive\tlower_bound_95\tbootstrap_seeds")
        println(io, join((gate.passed ? "PASS" : "FAIL", gate.every_fold_positive,
                          _format(gate.lower_bound), length(gate.bootstrap_seed_means)), '\t'))
        println(io, "fold_date\tunbootstrapped_improvement_per_lobe")
        for pair in gate.fold_improvements
            println(io, first(pair), '\t', _format(last(pair)))
        end
    end
end

struct MergeTransactionInterruption <: Exception
    point::Symbol
end

Base.showerror(io::IO, err::MergeTransactionInterruption) =
    print(io, "simulated merge transaction interruption at ", err.point)

_txn_present(path::AbstractString) = ispath(path) || islink(path)
_txn_marker(gate::AbstractString) = string(gate, ".stmfit-txn")
_txn_signature(destinations) = bytes2hex(sha256(join(abspath.(destinations), '\0')))
_txn_stage(path, id) = string(path, ".", id, ".stmfit-txn.new", splitext(path)[2])
_txn_backup(path, id) = string(path, ".", id, ".stmfit-txn.old")
_txn_absent(path, id) = string(path, ".", id, ".stmfit-txn.absent")

function _txn_rename(source::AbstractString, destination::AbstractString)
    rc = ccall(:rename, Cint, (Cstring, Cstring), source, destination)
    rc == 0 || Base.systemerror("rename $source to $destination", true)
    return nothing
end

function _txn_fsync_file(path::AbstractString)
    open(path, "r") do io
        rc = ccall(:fsync, Cint, (Cint,), fd(io))
        rc == 0 || Base.systemerror("fsync $path", true)
    end
    return nothing
end

function _txn_fsync_parent(path::AbstractString)
    directory = dirname(abspath(path))
    directory_fd = ccall(:open, Cint, (Cstring, Cint), directory, 0)
    directory_fd >= 0 || Base.systemerror("open directory $directory", true)
    try
        rc = ccall(:fsync, Cint, (Cint,), directory_fd)
        rc == 0 || Base.systemerror("fsync directory $directory", true)
    finally
        ccall(:close, Cint, (Cint,), directory_fd)
    end
    return nothing
end

function _txn_write_marker(marker, id, phase, signature)
    temp_path, io = mktemp(dirname(marker); cleanup=false)
    try
        println(io, "stmfit-output-transaction-v1")
        println(io, "id\t", id)
        println(io, "phase\t", phase)
        println(io, "signature\t", signature)
        flush(io)
        rc = ccall(:fsync, Cint, (Cint,), fd(io))
        rc == 0 || Base.systemerror("fsync merge transaction marker", true)
        close(io)
        _txn_rename(temp_path, marker)
        _txn_fsync_parent(marker)
    finally
        isopen(io) && close(io)
        rm(temp_path; force=true)
    end
    return nothing
end

function _txn_read_marker(marker, signature)
    islink(marker) && _fail("refusing symlinked merge transaction marker: $marker")
    lines = readlines(marker)
    length(lines) == 4 && lines[1] == "stmfit-output-transaction-v1" ||
        _fail("malformed merge transaction marker: $marker")
    fields = Dict{String,String}()
    for line in lines[2:end]
        parts = split(line, '\t'; limit=2)
        length(parts) == 2 || _fail("malformed merge transaction marker: $marker")
        fields[parts[1]] = parts[2]
    end
    get(fields, "signature", "") == signature ||
        _fail("merge transaction destination set does not match marker: $marker")
    phase = get(fields, "phase", "")
    phase in ("prepared", "committed") ||
        _fail("invalid merge transaction phase in $marker")
    id = get(fields, "id", "")
    occursin(r"^[0-9a-f]{16}$", id) || _fail("invalid merge transaction id in $marker")
    return id, phase
end

function _txn_cleanup(destinations, gate, id)
    for destination in destinations
        rm(_txn_stage(destination, id); force=true)
        rm(_txn_backup(destination, id); force=true)
        rm(_txn_absent(destination, id); force=true)
    end
    rm(_txn_marker(gate); force=true)
    _txn_fsync_parent(gate)
    return nothing
end

function _txn_restore_one(destination, id)
    backup = _txn_backup(destination, id)
    absent = _txn_absent(destination, id)
    if _txn_present(backup)
        _txn_present(destination) && rm(destination; force=true)
        _txn_rename(backup, destination)
    elseif isfile(absent)
        _txn_present(destination) && rm(destination; force=true)
        rm(absent; force=true)
    end
    return nothing
end

function _txn_rollback(destinations, gate, id)
    gate_processed = _txn_present(_txn_backup(gate, id)) || isfile(_txn_absent(gate, id))
    gate_processed && _txn_present(gate) && rm(gate; force=true)
    for destination in destinations
        destination == gate || _txn_restore_one(destination, id)
    end
    _txn_restore_one(gate, id)
    _txn_fsync_parent(gate)
    _txn_cleanup(destinations, gate, id)
    return nothing
end

function _txn_recover(destinations, gate)
    marker = _txn_marker(gate)
    _txn_present(marker) || return nothing
    id, phase = _txn_read_marker(marker, _txn_signature(destinations))
    phase == "committed" ? _txn_cleanup(destinations, gate, id) :
        _txn_rollback(destinations, gate, id)
    return nothing
end

function _txn_inject(failpoint, interruptpoint, point)
    interruptpoint == point && throw(MergeTransactionInterruption(point))
    failpoint == point && error("injected merge transaction failure at $point")
    return nothing
end

function _with_merge_transaction(writer::Function, destinations;
        gate=last(destinations), failpoint=nothing, interruptpoint=nothing)
    ordered = unique(abspath.(String.(destinations)))
    isempty(ordered) && throw(ArgumentError("merge transaction requires destinations"))
    gate = abspath(String(gate))
    gate in ordered || throw(ArgumentError("merge transaction gate must be a destination"))
    length(unique(dirname.(ordered))) == 1 ||
        throw(ArgumentError("merge transaction destinations must share one directory"))
    for destination in ordered
        isdir(destination) && !islink(destination) &&
            throw(ArgumentError("merge transaction destination is a directory: $destination"))
    end
    _txn_recover(ordered, gate)

    id = bytes2hex(sha256(string(time_ns(), ':', getpid(), ':', objectid(writer))))[1:16]
    staged = Dict(destination => _txn_stage(destination, id) for destination in ordered)
    marker = _txn_marker(gate)
    signature = _txn_signature(ordered)
    for destination in ordered
        for side_path in (staged[destination], _txn_backup(destination, id),
                          _txn_absent(destination, id))
            _txn_present(side_path) && error("merge transaction side path already exists: $side_path")
        end
    end
    try
        writer(staged)
        all(destination -> isfile(staged[destination]) && !islink(staged[destination]), ordered) ||
            error("one or more staged merge outputs are missing or invalid")
        foreach(_txn_fsync_file, values(staged))
        _txn_write_marker(marker, id, "prepared", signature)
        _txn_inject(failpoint, interruptpoint, :after_marker)

        backup_order = vcat([gate], filter(!=(gate), ordered))
        for (index, destination) in enumerate(backup_order)
            if _txn_present(destination)
                _txn_rename(destination, _txn_backup(destination, id))
            else
                write(_txn_absent(destination, id), "absent\n")
                _txn_fsync_file(_txn_absent(destination, id))
            end
            _txn_fsync_parent(destination)
            _txn_inject(failpoint, interruptpoint, Symbol("after_backup_$index"))
        end

        install_order = vcat(filter(!=(gate), ordered), [gate])
        for (index, destination) in enumerate(install_order)
            point = destination == gate ? :before_gate : Symbol("before_install_$index")
            _txn_inject(failpoint, interruptpoint, point)
            _txn_rename(staged[destination], destination)
            _txn_fsync_parent(destination)
            point = destination == gate ? :after_gate : Symbol("after_install_$index")
            _txn_inject(failpoint, interruptpoint, point)
        end
        _txn_write_marker(marker, id, "committed", signature)
        _txn_inject(failpoint, interruptpoint, :after_commit)
        _txn_cleanup(ordered, gate, id)
    catch err
        err isa MergeTransactionInterruption && rethrow()
        if _txn_present(marker)
            try
                marker_id, phase = _txn_read_marker(marker, signature)
                marker_id == id || error("merge transaction marker changed during commit")
                phase == "committed" ? _txn_cleanup(ordered, gate, id) :
                    _txn_rollback(ordered, gate, id)
            catch rollback_error
                error("merge transaction failed and recovery failed: $(sprint(showerror, err)); " *
                      "recovery: $(sprint(showerror, rollback_error))")
            end
        else
            foreach(path -> rm(path; force=true), values(staged))
        end
        rethrow()
    end
    return nothing
end

function merge_shards(paths, out_path, contract::MergeContract;
                      failpoint=nothing, interruptpoint=nothing)
    isempty(paths) && _fail("at least one shard is required")
    rows = Dict{String,Any}[]
    for path in sort!(unique(String.(paths)))
        append!(rows, _read_shard(path))
    end
    _validate_complete(rows, contract)
    sort!(rows; by=row -> (row["fold_date"], row["bootstrap_kind"] == "unbootstrapped" ? -1 : row["bootstrap_seed"], row["scan_key"]))
    gate = evaluate_gate(rows, contract.folds)
    gate_path = string(out_path, ".gate.tsv")
    destinations = [abspath(String(out_path)), abspath(gate_path)]
    _with_merge_transaction(destinations; gate=destinations[2], failpoint, interruptpoint) do staged
        _write_rows(staged[destinations[1]], rows)
        _write_gate(staged[destinations[2]], gate)
    end
    return (rows=rows, gate=gate, out=String(out_path), gate_out=gate_path)
end

function _show_help()
    println("""
Usage: julia --project=. test/merge_hierarchical_unit_assignment_shards.jl \\
  --shard PATH [--shard PATH ...] --out PATH --expected-folds D1,D2 \\
  --config-sha256 HEX --input-sha256 HEX --source-sha256 HEX

Strictly validates and deterministically merges complete leave-date-out shards.
Every fold must contain one unbootstrapped result per held-out scan and exactly
500 whole-scan bootstrap seeds (0:499). Prints terminal PASS or FAIL without
silently skipping or repairing corrupt shards.
""")
end

function _value(args, i, flag)
    i < length(args) || _fail("$flag requires a value")
    return String(args[i + 1])
end

function run_cli(args=ARGS)
    any(a -> a in ("-h", "--help"), args) && (_show_help(); return nothing)
    shards = String[]
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--shard"
            push!(shards, _value(args, i, arg)); i += 2
        elseif startswith(arg, "--shard=")
            push!(shards, split(arg, '='; limit=2)[2]); i += 1
        elseif arg in ("--out", "--expected-folds", "--config-sha256", "--input-sha256", "--source-sha256")
            values[arg] = _value(args, i, arg); i += 2
        elseif startswith(arg, "--") && occursin('=', arg)
            key, value = split(arg, '='; limit=2)
            key in ("--out", "--expected-folds", "--config-sha256", "--input-sha256", "--source-sha256") ||
                _fail("unknown argument: $arg")
            values[key] = value; i += 1
        else
            _fail("unknown argument: $arg")
        end
    end
    isempty(shards) && _fail("at least one --shard is required")
    for key in ("--out", "--expected-folds", "--config-sha256", "--input-sha256", "--source-sha256")
        haskey(values, key) || _fail("$key is required")
    end
    folds = [strip(x) for x in split(values["--expected-folds"], ',') if !isempty(strip(x))]
    contract = MergeContract(folds, values["--config-sha256"], values["--input-sha256"], values["--source-sha256"])
    result = merge_shards(shards, values["--out"], contract)
    println("TERMINAL ", result.gate.passed ? "PASS" : "FAIL")
    println("  rows: ", length(result.rows), "  seeds: ", length(result.gate.bootstrap_seed_means))
    println("  lower_bound_95: ", _format(result.gate.lower_bound))
    println("  merged: ", result.out)
    println("  gate: ", result.gate_out)
    return result
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
