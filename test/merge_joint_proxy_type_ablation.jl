#!/usr/bin/env julia

const LOBE_NAME = "real_type_ablation.tsv"
const SUMMARY_NAME = "real_type_ablation_summary.tsv"

function parse_args(args)
    shards = String[]
    outdir = ""
    i = 1
    while i <= length(args)
        flag = String(args[i])
        i < length(args) || error("$flag requires a value")
        if flag == "--shard"
            push!(shards, String(args[i + 1]))
        elseif flag == "--outdir"
            outdir = String(args[i + 1])
        else
            error("unknown argument: $flag")
        end
        i += 2
    end
    isempty(shards) && error("at least one --shard is required")
    isempty(outdir) && error("--outdir is required")
    return shards, outdir
end

function read_table(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("empty table: $path")
    header = split(lines[1], '\t')
    rows = [split(line, '\t') for line in lines[2:end] if !isempty(strip(line))]
    all(row -> length(row) == length(header), rows) || error("row width mismatch: $path")
    return header, rows
end

function merge_table(shards, name::String, key_fields::Vector{String}, output::String)
    expected_header = nothing
    merged = Vector{Vector{SubString{String}}}()
    seen = Set{Tuple}()
    for shard in shards
        header, rows = read_table(joinpath(shard, name))
        expected_header === nothing ? (expected_header = header) : header == expected_header ||
            error("header mismatch in $(joinpath(shard, name))")
        indices = [something(findfirst(==(field), header), 0) for field in key_fields]
        all(>(0), indices) || error("missing merge key in $name")
        for row in rows
            key = Tuple(row[index] for index in indices)
            key in seen && error("duplicate key in $name: $key")
            push!(seen, key)
            push!(merged, row)
        end
    end
    sort!(merged; by=row -> Tuple(begin
        value = row[something(findfirst(==(field), expected_header))]
        field in ("n", "lobe") ? lpad(value, 12, '0') : value
    end for field in key_fields))
    mkpath(dirname(output))
    open(output, "w") do io
        println(io, join(expected_header, '\t'))
        foreach(row -> println(io, join(row, '\t')), merged)
    end
    return length(merged)
end

function main(args)
    shards, outdir = parse_args(args)
    lobe_count = merge_table(shards, LOBE_NAME, ["file", "n", "lobe", "patch_kind", "ablation"],
        joinpath(outdir, LOBE_NAME))
    summary_count = merge_table(shards, SUMMARY_NAME, ["file", "patch_kind", "ablation"],
        joinpath(outdir, SUMMARY_NAME))
    println("type-ablation merge ok")
    println("  lobe rows: $lobe_count")
    println("  summary rows: $summary_count")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
