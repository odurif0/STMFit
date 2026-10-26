#!/usr/bin/env julia
# Label-free quality page of one run of run_molecule_consensus_chitosan.jl: the
# benchmark-free checks gathered in one place (process item 3, 2026-10-26;
# docs/src/journal.md). Reads the run's own tables, SXM headers and, with
# --control-run, a run of a fixed-sequence control; no label, no grade. The
# numbers are diagnostics: nothing is selected or changed from them.
#
# Usage:
#   julia --project=. scripts/pipeline/report_run_qc.jl --run RUN_DIR --data-dir RAW_DIR \
#       --outdir NEW_DIR [--control-run CONTROL_RUN_DIR]
#
# Sections: cohort; counts (consensus rules, count range across the scans of a
# molecule, split-half agreement of the molecule count); track integrity
# (scripts/data/check_track_integrity.jl); positions (position_repeatability);
# unit calls (fusion parameters, uncertain calls); with --control-run, the
# between-unit variance of brightness and shape against the control
# (scripts/data/report_unit_reliability.jl, scripts/data/report_unit_shape.jl).
using STMFit.ReconstructedUnitAssignment, STMFit.MoleculeConsensus
using STMSXMIO: read_sxm
using Printf, Statistics

const ROOT = dirname(dirname(@__DIR__))

function qc_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--data-dir", "--outdir", "--control-run") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    for k in ("--run", "--data-dir", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

function raw_index(dir)
    index = Dict{String,String}()
    for (d, _, names) in walkdir(dir; follow_symlinks=true), n in names
        endswith(lowercase(n), ".sxm") && (index[n] = joinpath(d, n))
    end
    return index
end

"Strict-majority count, else the lower median."
function vote(ns)
    t = Dict{Int,Int}()
    for x in ns
        t[x] = get(t, x, 0) + 1
    end
    n, k = first(sort(collect(t); by=p -> (-last(p), first(p))))
    return 2k > length(ns) ? n : sort(ns)[(length(ns) + 1) ÷ 2]
end

"Split-half differences of the molecule count: per track of at least 4 scans in time order, |vote(odd) − vote(even)|."
function split_half(members, n_scan)
    return [abs(vote([n_scan[f] for f in fs[1:2:end]]) - vote([n_scan[f] for f in fs[2:2:end]])) for fs in values(members) if length(fs) >= 4]
end

count_range(members, n_scan; min_scans=3) = [maximum(n_scan[f] for f in fs) - minimum(n_scan[f] for f in fs) for fs in values(members) if length(fs) >= min_scans]

run_script(path, args) = Base.run(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(ROOT, path)) $args`)

pct(a, b) = b == 0 ? "NA" : @sprintf("%.0f %%", 100a / b)

function qc_main(args=ARGS)
    o = qc_options(args)
    o === nothing && return println("report_run_qc.jl --run RUN_DIR --data-dir RAW_DIR --outdir NEW_DIR [--control-run CONTROL_RUN_DIR]")
    run, out = o["--run"], o["--outdir"]
    _, cons = read_table(joinpath(run, "consensus", "consensus.tsv"))
    n_scan = Dict(r["file"] => parse(Int, r["N_scan"]) for r in cons)
    n_final = Dict(r["file"] => parse(Int, r["N_final"]) for r in cons)
    members = Dict{String,Vector{String}}()
    for r in cons
        push!(get!(members, r["track"], String[]), r["file"])
    end
    raw = raw_index(o["--data-dir"])
    for fs in values(members)
        length(fs) >= 4 || continue
        times = Dict(f => haskey(raw, f) ? scan_frame(read_sxm(raw[f])).acquired : error("Raw scan missing: $f") for f in fs)
        sort!(fs; by=f -> (times[f], f))
    end
    mkpath(out)
    md = IOBuffer()
    println(md, "# Run quality (label-free): $(basename(abspath(run)))\n")
    println(md, "Diagnostics only: nothing is selected or changed from them. A molecule is a consensus track.\n")
    sizes = [length(fs) for fs in values(members)]
    println(md, "## Cohort\n")
    println(md, "| scans | molecules | molecules of ≥ 3 scans | scans in them | single-scan molecules | largest molecule |\n|---:|---:|---:|---:|---:|---:|")
    println(md, "| $(length(cons)) | $(length(members)) | $(count(>=(3), sizes)) | $(sum(s for s in sizes if s >= 3; init=0)) | $(count(==(1), sizes)) | $(maximum(sizes)) scans |\n")
    nf = [n_final[f] for f in keys(n_final)]
    println(md, "Final counts: median $(median(nf)), range $(minimum(nf))–$(maximum(nf)).\n")
    println(md, "## Counts\n")
    rules = Dict{String,Int}()
    for r in cons
        rules[r["rule"]] = get(rules, r["rule"], 0) + 1
    end
    println(md, "Consensus rules: " * join(["$k $v" for (k, v) in sort(collect(rules))], ", ") * ".\n")
    cr = count_range(members, n_scan); sh = split_half(members, n_scan)
    println(md, "| check | molecules | value |\n|---|---:|---|")
    println(md, "| per-scan counts all equal | $(length(cr)) | $(count(==(0), cr)) ($(pct(count(==(0), cr), length(cr)))) |")
    println(md, "| per-scan counts within 1 | $(length(cr)) | $(count(<=(1), cr)) ($(pct(count(<=(1), cr), length(cr)))) |")
    println(md, "| per-scan counts more than 2 apart | $(length(cr)) | $(count(>(2), cr)) ($(pct(count(>(2), cr), length(cr)))) |")
    println(md, "| molecule count, odd and even scans agree | $(length(sh)) | $(count(==(0), sh)) ($(pct(count(==(0), sh), length(sh)))); mean difference $(isempty(sh) ? "NA" : @sprintf("%.2f", mean(sh))) |\n")
    println(md, "## Track integrity (SXM headers and tip positions)\n")
    integ = joinpath(out, "integrity")
    try
        run_script("scripts/data/check_track_integrity.jl", ["--run", run, "--data-dir", o["--data-dir"], "--outdir", integ])
        lines = readlines(joinpath(integ, "report.md"))
        i = findfirst(l -> startswith(l, "| check"), lines)
        println(md, join(lines[i:min(end, i + 7)], "\n"), "\n")
        l = findfirst(l -> occursin("flagged", l), lines)
        l === nothing || println(md, lines[l], " Details: `integrity/`.\n")
    catch e
        println(md, "Not available: ", sprint(showerror, e), "\n")
    end
    println(md, "## Positions over repeated scans\n")
    pos = joinpath(run, "position_repeatability", "summary.tsv")
    if isfile(pos)
        _, rows = read_table(pos)
        println(md, "| frame | position | lobes | molecules | sd along (nm) | sd across (nm) |\n|---|---|---:|---:|---:|---:|")
        for r in rows
            println(md, "| $(r["mode"]) | $(r["position"]) | $(r["lobes"]) | $(r["tracks"]) | $(r["sd_along_nm"]) | $(r["sd_across_nm"]) |")
        end
        println(md)
    else
        println(md, "No position_repeatability/summary.tsv.\n")
    end
    println(md, "## Unit calls\n")
    fp = joinpath(run, "fusion", "fusion_params.tsv")
    if isfile(fp)
        _, f = read_table(fp); f = f[1]
        println(md, "Latent-class fusion: class-1 share π $(f["pi"]), false-call rate θ0 $(f["theta0"]), detection rate θ1 $(f["theta1"]), $(f["physical_lobes_fused"]) physical lobes fused. A high θ0 means per-scan calls disagree across repeat scans.\n")
    end
    cr_path = joinpath(run, "chain_report.tsv")
    if isfile(cr_path)
        _, ch = read_table(cr_path)
        units = sum(length(r["assignment"]) for r in ch); unc = sum(parse(Int, r["uncertain"]) for r in ch)
        conf = [parse(Float64, r["mean_confidence"]) for r in ch if r["mean_confidence"] != "NaN"]
        println(md, "Uncertain unit calls: $unc of $units ($(pct(unc, units))); median of the per-scan mean confidence $(isempty(conf) ? "NA" : @sprintf("%.2f", median(conf))).\n")
    end
    println(md, "The classes are physical (1 = brighter); their chemical identity is not established.\n")
    if haskey(o, "--control-run")
        println(md, "## Against a fixed-sequence control (`$(o["--control-run"])`)\n")
        try
            run_script("scripts/data/report_unit_reliability.jl", ["--run", "run=$run", "--run", "control=$(o["--control-run"])", "--control", "control", "--outdir", joinpath(out, "reliability")])
            run_script("scripts/data/report_unit_shape.jl", ["--run", "control=$(o["--control-run"])", "--run", "run=$run", "--control", "control", "--outdir", joinpath(out, "shape")])
            for (dir, title) in (("reliability", "Brightness and fitted widths"), ("shape", "Shape of the units (patches)"))
                lines = readlines(joinpath(out, dir, "report.md"))
                i = findfirst(l -> occursin("ratio", l) && startswith(l, "| run"), lines)
                i === nothing && continue
                j = findnext(l -> !startswith(l, "|"), lines, i)
                println(md, "**$title**, between-unit variance of this run over the control (details in `$dir/`):\n")
                println(md, join(lines[i:(j === nothing ? length(lines) : j - 1)], "\n"), "\n")
            end
        catch e
            println(md, "Not available: ", sprint(showerror, e), "\n")
        end
    end
    write(joinpath(out, "report.md"), take!(md))
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && qc_main()
