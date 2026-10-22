#!/usr/bin/env julia
# Data index and working sets (grading side: reads benchmark membership).
#
# Usage:
#   julia -t 4 --project=. scripts/grading/report_data_sets.jl --data data \
#       --table benchmarks/chitosan_6mer_cu100_triage.tsv --table benchmarks/chitosan_10_20mer_targets.tsv \
#       --table benchmarks/chitosan_10_20mer_cu100_triage.tsv --table benchmarks/chitosan_10_20mer_curated.tsv
#
# DATA/raw/ holds one copy of every raw file, by sample and session. Each
# table (columns file, set, motif, source, note, sha256, raw; one row per
# decided scan) places scans by sha256: a non-empty set puts the scan in
# DATA/sets/<set>/<file> as a relative symlink to raw/<raw>; an empty set means
# discarded, for the reason in motif. Rebuilds DATA/sets/ (symlinks only) and
# DATA/INDEX.tsv: one row per raw .sxm with where it is, what it is (sample,
# substrate, session, acquisition date and time, frame, pixels, bias, sha256)
# and what was decided (set, discarded or not triaged; motif, who, note).
# A raw file with the same bytes as another (for example a "- Copy" file) is
# listed as duplicate of the copy whose name has no "copy".
# Nothing here feeds the method.
using STMFit.ReconstructedUnitAssignment
using SHA, STMSXMIO, Printf

const SET_NAME = r"^[A-Za-z0-9_]+$"

function data_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    data = ""; tables = String[]; i = 1
    while i <= length(args)
        k = args[i]
        k in ("--data", "--table") && i < length(args) || error("Unknown option or missing value: $k")
        k == "--data" ? (data = args[i+1]) : push!(tables, args[i+1]); i += 2
    end
    (isempty(data) || isempty(tables)) && error("Required: --data DIR --table TSV [--table TSV ...]")
    isdir(joinpath(data, "raw")) || error("No raw/ folder in $data")
    return data, tables
end

"One row per raw .sxm: relative path, sha256 and header fields (one read per file)."
function raw_scans(raw)
    paths = String[]
    for (d, _, names) in walkdir(raw), n in names
        endswith(lowercase(n), ".sxm") && push!(paths, joinpath(d, n))
    end
    sort!(paths)
    rows = Vector{Dict{String,String}}(undef, length(paths))
    Threads.@threads for i in eachindex(paths)
        bytes = read(paths[i]); rel = relpath(paths[i], raw)
        idx = findfirst(Vector{UInt8}(":SCANIT_END:"), bytes)
        h = idx === nothing ? Dict{String,String}() : STMSXMIO._parse_header(String(bytes[1:first(idx)-1]))
        g(k) = strip(get(h, k, ""))
        nums(s) = [x for x in tryparse.(Float64, split(s)) if x !== nothing]
        r = nums(g("SCAN_RANGE")) .* 1e9
        session = basename(dirname(rel))
        rows[i] = Dict("raw_path" => rel, "sample" => first(splitpath(rel)), "session" => session,
            "file" => basename(rel), "sha256" => bytes2hex(sha256(bytes)),
            "substrate" => occursin("Ag111", session) ? "Ag(111)" : occursin("Cu100", session) ? "Cu(100)" : "",
            "date" => g("REC_DATE"), "time" => g("REC_TIME"), "bias_v" => g("BIAS"),
            "frame_nm" => length(r) >= 2 ? @sprintf("%.2fx%.2f", r[1], r[2]) : "",
            "pixels" => join(split(g("SCAN_PIXELS")), "x"))
    end
    return rows
end

function data_main(args=ARGS)
    parsed = data_options(args)
    parsed === nothing && return println("report_data_sets.jl --data DIR --table TSV [--table TSV ...]")
    data, tables = parsed
    raw = joinpath(data, "raw"); scans = raw_scans(raw)
    canon(s) = (occursin("copy", lowercase(s["file"])), length(s["raw_path"]), s["raw_path"])
    by_sha = Dict{String,Dict{String,String}}(); paths = Dict(s["raw_path"] => s for s in scans)
    for s in sort(scans; by=canon)
        haskey(by_sha, s["sha256"]) ? (s["duplicate_of"] = by_sha[s["sha256"]]["raw_path"]) : (by_sha[s["sha256"]] = s)
    end
    decided = Dict{String,Dict{String,String}}(); problems = String[]
    for t in tables, r in last(read_table(t))
        sha = r["sha256"]
        haskey(decided, sha) && push!(problems, "$(r["file"]): decided twice")
        haskey(by_sha, sha) || (push!(problems, "$(r["file"]) ($t): not found in raw/"); continue)
        haskey(paths, r["raw"]) && paths[r["raw"]]["sha256"] == sha ||
            push!(problems, "$(r["file"]): raw/$(r["raw"]) missing or with other bytes (found at $(by_sha[sha]["raw_path"]))")
        set = r["set"]
        isempty(set) ? (isempty(r["motif"]) && push!(problems, "$(r["file"]): discarded without motif")) :
            (occursin(SET_NAME, set) || push!(problems, "$(r["file"]): bad set name $set"))
        decided[sha] = r
    end
    names = Dict{Tuple{String,String},Int}()
    for r in values(decided); isempty(r["set"]) || (names[(r["set"], r["file"])] = get(names, (r["set"], r["file"]), 0) + 1); end
    for ((set, f), n) in names; n > 1 && push!(problems, "$f twice in set $set"); end
    isempty(problems) || error("Tables do not match raw/:\n  " * join(first(sort(problems), 20), "\n  "))
    sets = joinpath(data, "sets")
    if ispath(sets)
        for (d, _, fs) in walkdir(sets), f in fs
            islink(joinpath(d, f)) || error("sets/ holds a real file, not removed: $(joinpath(d, f))")
        end
        rm(sets; recursive=true)
    end
    for r in values(decided)
        isempty(r["set"]) && continue
        dir = mkpath(joinpath(sets, r["set"]))
        symlink(relpath(joinpath(raw, r["raw"]), dir), joinpath(dir, r["file"]))
    end
    cols = ["raw_path", "sample", "substrate", "session", "file", "date", "time", "frame_nm", "pixels", "bias_v",
            "sha256", "status", "motif", "decided_by", "note"]
    rows = map(scans) do s
        haskey(s, "duplicate_of") && return merge(s, Dict("status" => "duplicate", "motif" => "", "decided_by" => "auto",
                                                          "note" => "same bytes as $(s["duplicate_of"])"))
        r = get(decided, s["sha256"], nothing)
        status = r === nothing ? "not_triaged" : isempty(r["set"]) ? "discarded" : r["set"]
        merge(s, Dict("status" => status, "motif" => r === nothing ? "" : r["motif"],
                      "decided_by" => r === nothing ? "" : r["source"], "note" => r === nothing ? "" : r["note"]))
    end
    index = joinpath(data, "INDEX.tsv"); isfile(index) && rm(index)
    write_table(index, cols, rows)
    tally = Dict{String,Int}(); for r in rows; tally[r["status"]] = get(tally, r["status"], 0) + 1; end
    println("INDEX.tsv: $(length(rows)) raw scans; ", join(["$k $v" for (k, v) in sort(collect(tally))], ", "))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && data_main()
