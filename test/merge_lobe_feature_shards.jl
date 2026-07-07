#!/usr/bin/env julia

# Merge resumable per-lobe feature shards into one deterministic TSV.
#
# The reference TSV defines the allowed (file, lobe) keys and output row order.
# This is useful for split-width refits that are run in many small batches: old
# partial shards can be combined with resumed shards while preserving the exact
# selected-N denominator and rejecting missing or duplicate rows.

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent

const DEFAULT_REFERENCE = "results/unit_separability/lobe_features_selectedN.tsv"
const DEFAULT_OUT = "results/unit_separability/lobe_features_selectedN_merged.tsv"

struct Options
    reference::String
    shards::Vector{String}
    out_tsv::String
    ignore_extra::Bool
    dry_run::Bool
end

struct TsvData
    header::Vector{String}
    rows::Vector{Vector{String}}
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

function _append_shards!(shards::Vector{String}, value::AbstractString)
    for path in split(String(value), ','; keepempty=false)
        p = strip(path)
        isempty(p) || push!(shards, p)
    end
    return nothing
end

function _parse_cli(args)
    reference = DEFAULT_REFERENCE
    shards = String[]
    out_tsv = DEFAULT_OUT
    ignore_extra = false
    dry_run = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--reference"
            reference = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--reference=")
            reference = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--shard"
            push!(shards, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--shard=")
            push!(shards, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--shards"
            _append_shards!(shards, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--shards=")
            _append_shards!(shards, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--out"
            out_tsv = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--ignore-extra"
            ignore_extra = true; i += 1
        elseif arg == "--dry-run"
            dry_run = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/merge_lobe_feature_shards.jl [options]

            Options:
              --reference PATH    Reference feature TSV defining row keys/order [$(DEFAULT_REFERENCE)]
              --shard PATH        Shard TSV to merge. Repeatable.
              --shards LIST       Comma-separated shard TSV paths.
              --out PATH          Output merged TSV [$(DEFAULT_OUT)]
              --ignore-extra      Ignore shard rows whose (file,lobe) key is not in the reference.
                                  Without this flag, extra rows are an error.
              --dry-run           Validate and print counts without writing --out.

            The output header is copied from the first shard and rows are ordered exactly
            as in --reference. Missing reference keys and duplicate shard keys are errors.

            Label-free constraint: this script reads no truth sequence, expected N,
            control motif, or composition count. It only merges already generated feature TSVs.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isfile(reference) || error("Reference TSV not found: $reference")
    isempty(shards) && error("At least one --shard or --shards path is required")
    for shard in shards
        isfile(shard) || error("Shard TSV not found: $shard")
    end
    return Options(reference, shards, out_tsv, ignore_extra, dry_run)
end

function _read_raw_tsv(path::AbstractString)
    lines = readlines(path)
    data = filter(l -> !isempty(strip(l)) && !startswith(strip(l), '#'), lines)
    isempty(data) && error("No TSV data in $path")
    header = split(data[1], '\t'; keepempty=true)
    rows = [split(line, '\t'; keepempty=true) for line in data[2:end]]
    return TsvData(header, rows)
end

function _column_index(header::Vector{String}, name::String, source::AbstractString)
    idx = findfirst(==(name), header)
    idx === nothing && error("$source missing $name column")
    return idx
end

function _key(row::Vector{String}, file_idx::Int, lobe_idx::Int, source::AbstractString)
    length(row) >= max(file_idx, lobe_idx) || error("Short TSV row in $source")
    file = basename(strip(row[file_idx]))
    isempty(file) && error("Empty file value in $source")
    lobe = tryparse(Int, strip(row[lobe_idx]))
    lobe === nothing && error("Invalid lobe index in $source: $(row[lobe_idx])")
    return (file, lobe)
end

function _reference_keys(reference::TsvData, source::AbstractString)
    file_idx = _column_index(reference.header, "file", source)
    lobe_idx = _column_index(reference.header, "lobe", source)
    keys = Tuple{String,Int}[]
    seen = Set{Tuple{String,Int}}()
    for row in reference.rows
        key = _key(row, file_idx, lobe_idx, source)
        key in seen && error("Duplicate reference key: $(key[1]) lobe $(key[2])")
        push!(seen, key)
        push!(keys, key)
    end
    return keys
end

function _merge_shards(opt::Options)
    reference = _read_raw_tsv(opt.reference)
    keys = _reference_keys(reference, opt.reference)
    key_set = Set(keys)

    merged = Dict{Tuple{String,Int},Vector{String}}()
    merged_header = String[]
    duplicates = Tuple{String,Int}[]
    extras = Tuple{String,Int}[]

    for shard in opt.shards
        data = _read_raw_tsv(shard)
        if isempty(merged_header)
            merged_header = data.header
        elseif data.header != merged_header
            error("Header mismatch in shard: $shard")
        end

        file_idx = _column_index(data.header, "file", shard)
        lobe_idx = _column_index(data.header, "lobe", shard)
        for row in data.rows
            key = _key(row, file_idx, lobe_idx, shard)
            if !(key in key_set)
                push!(extras, key)
                continue
            end
            if haskey(merged, key)
                push!(duplicates, key)
                continue
            end
            merged[key] = row
        end
    end

    missing = [key for key in keys if !haskey(merged, key)]
    return reference, keys, merged_header, merged, missing, duplicates, extras
end

function _print_examples(label::String, keys::Vector{Tuple{String,Int}})
    isempty(keys) && return nothing
    println(label, ":")
    for key in keys[1:min(end, 20)]
        println("  ", key[1], "\tlobe=", key[2])
    end
    length(keys) > 20 && println("  ... ", length(keys) - 20, " more")
    return nothing
end

function main()
    opt = _parse_cli(ARGS)
    _, keys, header, merged, missing, duplicates, extras = _merge_shards(opt)

    reference_files = Set(key[1] for key in keys)
    merged_files = Set(key[1] for key in keys if haskey(merged, key))
    extra_files = Set(key[1] for key in extras)

    println("reference_files=", length(reference_files))
    println("reference_rows=", length(keys))
    println("merged_files=", length(merged_files))
    println("merged_rows=", length(merged))
    println("missing_rows=", length(missing))
    println("duplicate_rows=", length(duplicates))
    println("extra_rows=", length(extras))
    isempty(extras) || println("extra_files=", join(sort(collect(extra_files)), ","))

    _print_examples("Missing rows", missing)
    _print_examples("Duplicate rows", duplicates)
    (!opt.ignore_extra) && _print_examples("Extra rows", extras)

    isempty(missing) || error("Cannot write merged TSV with missing reference rows")
    isempty(duplicates) || error("Cannot write merged TSV with duplicate shard rows")
    opt.ignore_extra || isempty(extras) || error("Shard rows outside the reference were found; pass --ignore-extra to skip them")

    if opt.dry_run
        println("dry_run=true")
        return nothing
    end

    _ensure_parent(opt.out_tsv)
    open(opt.out_tsv, "w") do io
        println(io, join(header, '\t'))
        for key in keys
            println(io, join(merged[key], '\t'))
        end
    end
    println("wrote=", opt.out_tsv)
    return nothing
end

main()
