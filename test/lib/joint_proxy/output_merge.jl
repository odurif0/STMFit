module JointProxyShardMerge

include("output_validation.jl")
using .JointProxyOutputValidation

include(joinpath(@__DIR__, "..", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv

using Printf
using TOML

export MergeError, MergeOptions, parse_cli, merge_artifacts, run_cli

struct MergeError <: Exception
    msg::String
end

Base.showerror(io::IO, e::MergeError) = print(io, "MergeError: ", e.msg)

Base.@kwdef struct MergeOptions
    shards::Vector{String}
    outdir::String
end

function _fail(msg::AbstractString)
    throw(MergeError(String(msg)))
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || _fail("$flag requires a value")
    return String(args[i + 1])
end

function parse_cli(args::AbstractVector{<:AbstractString})::MergeOptions
    shards = String[]
    outdir = ""
    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--shard"
            push!(shards, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--shard=")
            push!(shards, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--shards"
            append!(shards, [strip(x) for x in split(_arg_value(args, i, arg), ',') if !isempty(strip(x))]); i += 2
        elseif startswith(arg, "--shards=")
            append!(shards, [strip(x) for x in split(split(arg, "="; limit=2)[2], ',') if !isempty(strip(x))]); i += 1
        elseif arg == "--outdir"
            outdir = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "="; limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/merge_joint_proxy_shards.jl --shard DIR [--shard DIR ...] --outdir DIR

            Merges disjoint joint-proxy artifact shards into canonical TSV/TOML outputs.
            """)
            exit(0)
        else
            _fail("Unknown argument: $arg")
        end
    end
    isempty(shards) && _fail("at least one --shard is required")
    isempty(outdir) && _fail("--outdir is required")
    return MergeOptions(shards=shards, outdir=outdir)
end

function _read_rows(dir::AbstractString)
    return (
        candidate_n = _read_tsv(joinpath(dir, "candidate_n.tsv"))[2],
        candidate_lobes = _read_tsv(joinpath(dir, "candidate_lobes.tsv"))[2],
        predictions = _read_tsv(joinpath(dir, "predictions.tsv"))[2],
        summary = _read_tsv(joinpath(dir, "chain_summary.tsv"))[2],
        manifest = TOML.parsefile(joinpath(dir, "run_manifest.toml")),
        calibration = TOML.parsefile(joinpath(dir, "joint_proxy_calibration.toml")),
    )
end

function _contract(manifest::AbstractDict)
    prov = manifest["provenance"]
    inputs = manifest["inputs"]
    return (
        config_sha256 = lowercase(String(prov["config_sha256"])),
        source_sha256 = lowercase(String(prov["source_sha256"])),
        payload_sha256 = lowercase(String(prov["payload_sha256"])),
        uncalibrated = Bool(prov["uncalibrated"]),
        config = String(inputs["config"]),
        data_dir = String(inputs["data_dir"]),
    )
end

function _sort_component(row, field::String)
    value = string(get(row, field, ""))
    return field in ("n", "lobe") ? (0, parse(Int, value)) : (1, value)
end

function _stable_sort!(rows, fields::Vector{String})
    sort!(rows; by = row -> Tuple(_sort_component(row, field) for field in fields))
    return rows
end

function _write_tsv(path::AbstractString, fields::Vector{String}, rows)
    _ensure_parent(path)
    open(path, "w") do io
        println(io, join(fields, '\t'))
        for row in rows
            println(io, join((string(get(row, f, "")) for f in fields), '\t'))
        end
    end
end

_toml_q(s::AbstractString) = '"' * replace(String(s), '"' => "\\\"") * '"'
_toml_fmt(x::AbstractString) = _toml_q(x)
_toml_fmt(x::Bool) = x ? "true" : "false"
_toml_fmt(x::Integer) = string(x)
_toml_fmt(x::Real) = isfinite(Float64(x)) ? @sprintf("%.15g", Float64(x)) : "nan"
_toml_fmt(x::AbstractVector) = "[" * join((_toml_fmt(v) for v in x), ", ") * "]"
_toml_fmt(x::Tuple) = _toml_fmt(collect(x))
_toml_fmt(x) = x === nothing ? "\"\"" : _toml_q(string(x))

function _write_toml(io::IO, name::String, value)
    println(io, "[$name]")
    if value isa AbstractDict
        for (k, v) in pairs(value)
            if v isa AbstractDict || v isa NamedTuple
                _write_toml(io, "$name.$(String(k))", v)
            elseif v isa AbstractVector && !isempty(v) && (first(v) isa AbstractDict || first(v) isa NamedTuple)
                for item in v
                    println(io, "[[$name.$(String(k))]]")
                    for (ik, iv) in pairs(item)
                        println(io, "$(String(ik)) = $(_toml_fmt(iv))")
                    end
                end
            else
                println(io, "$(String(k)) = $(_toml_fmt(v))")
            end
        end
    else
        for (k, v) in pairs(value)
            println(io, "$(String(k)) = $(_toml_fmt(v))")
        end
    end
    println(io)
end

function _validate_shards(shards)
    first_contract = nothing
    files = String[]
    candidate_n = Any[]
    candidate_lobes = Any[]
    predictions = Any[]
    summaries = Any[]
    seen_files = Set{String}()
    for shard in shards
        JointProxyOutputValidation.validate_artifacts(ValidationOptions(artifacts=shard, calibration=joinpath(shard, "joint_proxy_calibration.toml")))
        rows = _read_rows(shard)
        contract = _contract(rows.manifest)
        first_contract === nothing && (first_contract = contract)
        contract == first_contract || _fail("manifest/hash/schema mismatch across shards")
        shard_files = sort(String.(rows.manifest["inputs"]["files"]))
        for file in shard_files
            file in seen_files && _fail("overlapping shard: $file")
            push!(seen_files, file)
            push!(files, file)
        end
        append!(candidate_n, rows.candidate_n)
        append!(candidate_lobes, rows.candidate_lobes)
        append!(predictions, rows.predictions)
        append!(summaries, rows.summary)
    end
    return first_contract, sort(files), candidate_n, candidate_lobes, predictions, summaries
end

function _manifest_from(contract, files, counts, outdir::AbstractString)
    return Dict(
        "provenance" => Dict(
            "config_sha256" => contract.config_sha256,
            "source_sha256" => contract.source_sha256,
            "payload_sha256" => contract.payload_sha256,
            "calibration_path" => "joint_proxy_calibration.toml",
            "uncalibrated" => contract.uncalibrated,
        ),
        "inputs" => Dict(
            "config" => contract.config,
            "data_dir" => contract.data_dir,
            "files" => files,
            "chunk" => "none",
        ),
        "outputs" => Dict(
            "candidate_n_rows" => counts[1],
            "candidate_lobes_rows" => counts[2],
            "predictions_rows" => counts[3],
            "chain_summary_rows" => counts[4],
        ),
    )
end

function merge_artifacts(opt::MergeOptions)
    contract, files, candidate_n, candidate_lobes, predictions, summaries = _validate_shards(opt.shards)
    _stable_sort!(candidate_n, JointProxyOutputValidation.CANDIDATE_N_FIELDS)
    _stable_sort!(candidate_lobes, JointProxyOutputValidation.CANDIDATE_LOBE_FIELDS)
    _stable_sort!(predictions, JointProxyOutputValidation.PREDICTION_FIELDS)
    _stable_sort!(summaries, JointProxyOutputValidation.SUMMARY_FIELDS)

    outdir = abspath(opt.outdir)
    mkpath(outdir)
    _write_tsv(joinpath(outdir, "candidate_n.tsv"), JointProxyOutputValidation.CANDIDATE_N_FIELDS, candidate_n)
    _write_tsv(joinpath(outdir, "candidate_lobes.tsv"), JointProxyOutputValidation.CANDIDATE_LOBE_FIELDS, candidate_lobes)
    _write_tsv(joinpath(outdir, "predictions.tsv"), JointProxyOutputValidation.PREDICTION_FIELDS, predictions)
    _write_tsv(joinpath(outdir, "chain_summary.tsv"), JointProxyOutputValidation.SUMMARY_FIELDS, summaries)
    mani = _manifest_from(contract, files, (length(candidate_n), length(candidate_lobes), length(predictions), length(summaries)), outdir)
    open(joinpath(outdir, "run_manifest.toml"), "w") do io
        for (name, value) in pairs(mani)
            _write_toml(io, String(name), value)
        end
    end
    cp(joinpath(first(opt.shards), "joint_proxy_calibration.toml"), joinpath(outdir, "joint_proxy_calibration.toml"); force=true)
    return outdir
end

function run_cli(args::AbstractVector{<:AbstractString}=String[])
    opt = parse_cli(args)
    outdir = merge_artifacts(opt)
    println("joint proxy merge ok")
    println("  outdir: ", outdir)
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cli(ARGS)
end

end # module
