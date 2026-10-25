#!/usr/bin/env julia
# Frozen image-quality score per scan, from a catalog of catalog_sxm_candidates.jl.
# Label-free: image values only, no fit, no count, no method output.
#
# Usage:
#   julia --project=. scripts/data/score_scan_quality.jl --catalog CATALOG_TSV \
#       --config config/scan_quality.toml --out SCORE_TSV
#   julia --project=. scripts/data/score_scan_quality.jl --catalog CATALOG_TSV --derive-reference
#
# Seven axes, each turned into a percentile rank against a frozen reference
# stored in the config ([reference.AXIS]: 21 quantile knots at 0, 5, ..., 100%
# of the reference scans). rank = share of reference knots below the value
# (ties count half), reversed for the axes where lower is better. The score is
# the mean rank, so a scan is scored the same whatever set it is in.
# --derive-reference prints the [reference] tables of the technically usable
# scans of the catalog, to paste into the config.
using STMFit.ReconstructedUnitAssignment: write_table
using Printf, Statistics, TOML

const HIGHER = ("contrast_to_noise", "fwd_bwd_corr", "roi_border_margin_nm")
const LOWER = ("line_discontinuity", "stripe_periodicity", "fwd_bwd_nrmse", "second_object_ratio")
const AXES = (HIGHER..., LOWER...)

function score_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        if k == "--derive-reference"; o[k] = "true"; i += 1; continue; end
        k in ("--catalog", "--config", "--out") && i < length(args) || error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    haskey(o, "--catalog") || error("Required: --catalog")
    haskey(o, "--derive-reference") || (haskey(o, "--config") && haskey(o, "--out")) || error("Required: --config --out (or --derive-reference)")
    haskey(o, "--out") && (ispath(o["--out"]) || islink(o["--out"])) && error("Output exists: $(o["--out"])")
    return o
end

function read_catalog(path)
    lines = [l for l in readlines(path) if !isempty(strip(l))]
    header = String.(split(lines[1], '\t'; keepempty=true))
    return [Dict(zip(header, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

value(r, k) = something(tryparse(Float64, get(r, k, "")), NaN)

"Share of the reference knots below x, ties counting half (the empirical CDF at x)."
cdf(knots, x) = mean(k < x ? 1.0 : (k == x ? 0.5 : 0.0) for k in knots)

function axis_rank(knots, axis, x)
    isfinite(x) || return NaN
    c = cdf(knots, x)
    return axis in HIGHER ? c : 1 - c
end

function main(args=ARGS)
    o = score_options(args)
    o === nothing && return println("score_scan_quality.jl --catalog TSV --config TOML --out TSV | --catalog TSV --derive-reference")
    rows = read_catalog(o["--catalog"])
    if haskey(o, "--derive-reference")
        usable = [r for r in rows if r["status"] == "usable"]
        println("# reference knots (0, 5, ..., 100%) over $(length(usable)) usable scans of $(o["--catalog"])")
        for a in AXES
            v = filter(isfinite, [value(r, a) for r in usable])
            println("[reference.$a]"); println("knots = [", join((@sprintf("%.6g", quantile(v, q)) for q in 0:0.05:1), ", "), "]")
        end
        return
    end
    ref = TOML.parsefile(o["--config"])["reference"]
    out = Dict{String,String}[]
    for r in rows
        r["status"] == "usable" || continue
        ranks = [axis_rank(Float64.(ref[a]["knots"]), a, value(r, a)) for a in AXES]
        score = any(isnan, ranks) ? NaN : mean(ranks)
        push!(out, merge(Dict("file" => r["file"], "sha256" => get(r, "sha256", ""), "session" => get(r, "session", ""),
                              "score" => isnan(score) ? "NaN" : @sprintf("%.4f", score)),
                         Dict(a => @sprintf("%.3f", x) for (a, x) in zip(AXES, ranks))))
    end
    write_table(o["--out"], vcat(["file", "sha256", "session", "score"], collect(AXES)), out)
    println("Scored $(length(out)) usable scans; median score $(round(median(filter(isfinite, parse.(Float64, [r["score"] for r in out]))), digits=3))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
