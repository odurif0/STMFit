#!/usr/bin/env julia
# Collect label-free repeat-scan candidates of target molecules from a raw tree.
#
# Usage:
#   julia --project=. test/collect_repeat_scans.jl --targets targets.txt \
#       --raw-root RAW_TREE --config config/molecule_consensus.toml --outdir NEW_DIR
#
# For every target scan, same-session scans within collect_max_gap_min of its
# acquisition time, whose frame centre lies within collect_max_offset_nm and
# whose range is at most collect_max_range_nm, are candidates. Only SXM headers
# are read; byte-identical duplicates are dropped. The output directory holds
# symlinks to targets plus candidates and manifest.tsv. Whether a candidate is
# the same unchanged molecule is decided later by image registration.
using STMSXMIO, TOML, Dates, SHA, Printf
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: write_table

function header_only(path)
    io = open(path); buf = UInt8[]
    while !eof(io)
        append!(buf, read(io, 65536))
        idx = findfirst(Vector{UInt8}(":SCANIT_END:"), buf)
        if idx !== nothing
            close(io); return STMSXMIO._parse_header(String(buf[1:first(idx)-1]))
        end
    end
    close(io); error("No SCANIT_END in $path")
end

function frame_meta(path)
    h = header_only(path)
    t = DateTime(strip(h["REC_DATE"]) * " " * strip(h["REC_TIME"]), dateformat"dd.mm.yyyy HH:MM:SS")
    off = parse.(Float64, split(strip(h["SCAN_OFFSET"]))) .* 1e9
    rng = parse.(Float64, split(strip(h["SCAN_RANGE"]))) .* 1e9
    return (time=t, session=strip(h["REC_DATE"]), offset=off, range=rng)
end

function main(args=ARGS)
    ("--help" in args || isempty(args)) && return println("collect_repeat_scans.jl --targets FILE --raw-root DIR --config TOML --outdir NEW_DIR")
    o = Dict(args[i] => args[i+1] for i in 1:2:length(args))
    for k in ("--targets", "--raw-root", "--config", "--outdir"); haskey(o, k) || error("Missing $k"); end
    ispath(o["--outdir"]) && error("Output exists")
    sel = TOML.parsefile(o["--config"])["selection"]
    gap = Float64(sel["collect_max_gap_min"]); dmax = Float64(sel["collect_max_offset_nm"]); rmax = Float64(sel["collect_max_range_nm"])
    targets = [strip(l) for l in readlines(o["--targets"]) if !isempty(strip(l))]
    paths = Dict{String,String}(); ambiguous = Set{String}()
    for (d, _, names) in walkdir(o["--raw-root"]), n in names
        endswith(lowercase(n), ".sxm") || continue
        p = joinpath(d, n)
        if haskey(paths, n)
            # Identical copies in two folders are one scan; differing contents are ambiguous.
            read(paths[n]) == read(p) || push!(ambiguous, n)
        else
            paths[n] = p
        end
    end
    for n in ambiguous; delete!(paths, n); end
    all(t -> haskey(paths, t), targets) || error("A target is missing or ambiguous in the raw tree")
    meta = Dict{String,Any}()
    for (n, p) in paths
        try meta[n] = frame_meta(p) catch; end
    end
    chosen = Dict{String,Dict{String,String}}(); hashes = Dict{String,String}()
    for t in targets
        mt = meta[t]; hashes[bytes2hex(sha256(read(paths[t])))] = t
        chosen[t] = Dict("file" => t, "target" => t, "role" => "target", "dt_min" => "0", "offset_nm" => "0")
    end
    for t in targets, (n, m) in meta
        n in keys(chosen) && continue
        mt = meta[t]
        m.session == mt.session || continue
        dt = abs(Dates.value(m.time - mt.time)) / 60000
        dist = hypot((m.offset .- mt.offset)...)
        dt <= gap && dist <= dmax && m.range[1] <= rmax && m.range[2] <= rmax || continue
        h = bytes2hex(sha256(read(paths[n])))
        haskey(hashes, h) && continue  # byte-identical duplicate of a chosen scan
        hashes[h] = n
        chosen[n] = Dict("file" => n, "target" => t, "role" => "repeat_candidate",
                         "dt_min" => @sprintf("%.2f", dt), "offset_nm" => @sprintf("%.3f", dist))
    end
    out = abspath(o["--outdir"]); raw = joinpath(out, "raw"); mkpath(raw)
    for n in sort(collect(keys(chosen))); symlink(abspath(paths[n]), joinpath(raw, n)); end
    write(joinpath(out, "ambiguous_basenames.txt"), join(sort(collect(ambiguous)), "\n"))
    write_table(joinpath(out, "manifest.tsv"), ["file", "target", "role", "dt_min", "offset_nm"],
                [chosen[n] for n in sort(collect(keys(chosen)))])
    println("Collected $(length(targets)) targets and $(length(chosen) - length(targets)) repeat candidates -> $raw")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
