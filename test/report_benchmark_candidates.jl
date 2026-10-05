#!/usr/bin/env julia
# Review sheets for new 6-mer benchmark candidates (grading side, human decision).
#
# Usage:
#   julia -t 4 --project=. test/report_benchmark_candidates.jl --catalog CATALOG_TSV \
#       --count-config config/chitosan.toml --outdir NEW_DIR \
#       [--manifest benchmarks/chitosan_6mer_counting_confirmed.toml] \
#       [--review benchmarks/chitosan_6mer_preassignment_review.tsv]
#
# Splits the usable scans of a catalog_sxm_candidates.jl catalog into
# in_benchmark (manifest files), reviewed_before (a decision in the
# pre-assignment review) and never_reviewed, matched on session/file. Writes
# candidates.tsv and contact sheets of the never-reviewed scans in
# sheets/no_flags/ and sheets/flagged/. Benchmark membership stays a human
# decision; nothing here feeds the method.
include(joinpath(@__DIR__, "catalog_sxm_candidates.jl"))

const DEFAULT_MANIFEST = joinpath(@__DIR__, "..", "benchmarks", "chitosan_6mer_counting_confirmed.toml")
const DEFAULT_REVIEW = joinpath(@__DIR__, "..", "benchmarks", "chitosan_6mer_preassignment_review.tsv")
const NOT_REVIEWED = "unclassified_no_batch_plot"

"Grading-side TSV reader (the method's read_table refuses label columns)."
function grading_rows(path)
    lines = filter(l -> !isempty(strip(l)) && !startswith(l, '#'), readlines(path))
    h = split(lines[1], '\t')
    return [Dict(zip(h, split(l, '\t'; keepempty=true))) for l in lines[2:end]]
end

function candidate_groups(rows, manifest, review)
    bench = Set(v["relative_path"] for v in values(manifest["files"]))
    reviewed = Dict(r["relative_path"] => r for r in review if r["pre_assignment"] != NOT_REVIEWED)
    out = Dict{String,String}[]
    for r in rows
        r["status"] == "usable" || continue
        key = string(r["session"], "/", r["file"])
        group = key in bench ? "in_benchmark" : haskey(reviewed, key) ? "reviewed_before" : "never_reviewed"
        note = group == "reviewed_before" ? string(reviewed[key]["pre_assignment"], ": ", reviewed[key]["reason"]) : ""
        push!(out, merge(r, Dict("group" => group, "review" => note)))
    end
    return out
end

function report_main(args=ARGS)
    ("--help" in args || "-h" in args || isempty(args)) &&
        return println("report_benchmark_candidates.jl --catalog TSV --count-config TOML --outdir NEW_DIR [--manifest TOML] [--review TSV]")
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--catalog", "--count-config", "--outdir", "--manifest", "--review") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--catalog", "--count-config", "--outdir")) || error("Required: --catalog --count-config --outdir")
    out = o["--outdir"]; (ispath(out) || islink(out)) && error("Output exists: $out")
    rows = grading_rows(o["--catalog"])
    cands = candidate_groups(rows, TOML.parsefile(get(o, "--manifest", DEFAULT_MANIFEST)),
                             grading_rows(get(o, "--review", DEFAULT_REVIEW)))
    mkpath(out)
    cols = ["group", "file", "session", "central_name", "range_nm", "flags", "review", "path"]
    write_table(joinpath(out, "candidates.tsv"), cols,
        [Dict(c => replace(get(r, c, ""), '\t' => ' ', '\n' => ' ') for c in cols)
         for r in sort(cands; by=r -> (r["group"], r["session"], r["file"]))])
    pcfg = count_preprocessing(o["--count-config"])
    new = [r for r in cands if r["group"] == "never_reviewed"]
    views = Vector{Any}(undef, length(new))
    Threads.@threads for k in eachindex(new)
        views[k] = last(inspect_scan(new[k]["path"], pcfg, Inf))
    end
    vd = Dict(new[k]["path"] => views[k] for k in eachindex(new))
    clean = [r for r in new if isempty(r["flags"])]; flagged = [r for r in new if !isempty(r["flags"])]
    p1 = contact_sheets(clean, vd, joinpath(out, "sheets", "no_flags"); title="Never reviewed, no flag")
    p2 = contact_sheets(flagged, vd, joinpath(out, "sheets", "flagged"); title="Never reviewed, flagged")
    tally = Dict{String,Int}(); for r in cands; tally[r["group"]] = get(tally, r["group"], 0) + 1; end
    println("Usable scans: ", join(["$k $v" for (k, v) in sort(collect(tally))], ", "),
            "; never reviewed: $(length(clean)) without flag ($p1 pages), $(length(flagged)) flagged ($p2 pages)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && report_main()
