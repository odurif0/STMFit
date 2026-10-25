#!/usr/bin/env julia
# Chain lengths of one or more runs against length readings (10-20mers).
# Grading side: reads the readings, which the method never sees.
#
# Usage:
#   julia --project=. scripts/grading/report_length_readings.jl --lengths TSV \
#       --run NAME=RUN_DIR [--run NAME=RUN_DIR ...] --outdir NEW_DIR
#
# --lengths: columns file, length, source (for example experimenters or
# third_party). Each run's chain_report.tsv gives N_scan (per scan) and
# N_selected (after the molecule consensus). For every run, source and count
# (scan, final) the summary gives the readings found in the run, exact, within
# one unit and the mean |N - length|; per_scan.tsv lists every comparison.
using STMFit.ReconstructedUnitAssignment: read_table, write_table
using Printf, Statistics

function options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,Any}("--run" => Pair{String,String}[]); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--lengths", "--run", "--outdir") && i < length(args) || error("Unknown option or missing value: $k")
        if k == "--run"
            name, dir = split(args[i+1], "="; limit=2)
            push!(o[k], String(name) => String(dir))
        else
            o[k] = args[i+1]
        end
        i += 2
    end
    haskey(o, "--lengths") && haskey(o, "--outdir") && !isempty(o["--run"]) || error("Required: --lengths --run --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

function report(o)
    _, readings = read_table(o["--lengths"])
    rows = Dict{String,String}[]
    for (name, dir) in o["--run"]
        _, chains = read_table(joinpath(dir, "chain_report.tsv"))
        byfile = Dict(basename(r["file"]) => r for r in chains)
        for r in readings
            haskey(byfile, r["file"]) || continue
            c = byfile[r["file"]]; len = parse(Int, r["length"])
            push!(rows, Dict("run" => name, "file" => r["file"], "source" => r["source"], "length" => r["length"],
                "N_scan" => c["N_scan"], "N_final" => c["N_selected"], "count_rule" => c["count_rule"],
                "d_scan" => string(parse(Int, c["N_scan"]) - len), "d_final" => string(parse(Int, c["N_selected"]) - len)))
        end
    end
    out = mkpath(o["--outdir"])
    write_table(joinpath(out, "per_scan.tsv"),
        ["run", "file", "source", "length", "N_scan", "N_final", "count_rule", "d_scan", "d_final"], rows)
    summary = Dict{String,String}[]
    for (name, _) in o["--run"], source in sort(unique(r["source"] for r in readings)), kind in ("scan", "final")
        d = [abs(parse(Int, r["d_$kind"])) for r in rows if r["run"] == name && r["source"] == source]
        isempty(d) && continue
        push!(summary, Dict("run" => name, "source" => source, "count" => kind, "readings" => string(length(d)),
            "exact" => string(count(==(0), d)), "within_one" => string(count(<=(1), d)),
            "mean_abs_diff" => @sprintf("%.3f", mean(d))))
    end
    write_table(joinpath(out, "summary.tsv"), ["run", "source", "count", "readings", "exact", "within_one", "mean_abs_diff"], summary)
    for s in summary
        println(join((s["run"], s["source"], s["count"], s["readings"], s["exact"], s["within_one"], s["mean_abs_diff"]), '\t'))
    end
end

function main(args=ARGS)
    o = options(args)
    o === nothing ? println("report_length_readings.jl --lengths TSV --run NAME=RUN_DIR [--run ...] --outdir NEW_DIR") : report(o)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
