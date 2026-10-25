#!/usr/bin/env julia
# Automatic image-quality decision per scan, from a catalog of
# catalog_sxm_candidates.jl. Label-free: no fit, no count, no method output.
#
# Usage:
#   julia --project=. scripts/data/decide_scan_quality.jl --catalog CATALOG_TSV \
#       --config config/scan_quality.toml --out QUALITY_TSV [--size-check]
#   julia --project=. scripts/data/decide_scan_quality.jl --catalog CATALOG_TSV --derive
#
# A scan is accepted when its technical status is usable, its molecule region
# touches neither the border nor the unscanned area (a scan stopped early is
# fine if the molecule is whole), and every value passes its bound in
# [selection]. --size-check (6-mer benchmark only, never for unknown lengths)
# also rejects a ROI length outside size_window x the median ROI length of the
# otherwise accepted scans; the decision without it is kept in its own column.
# --derive prints the robust bounds that the config stores as fixed numbers:
# median +- 3 x 1.4826 x MAD of log(value) over the scans whose status is
# usable (the values are positive and strongly right-skewed).
using STMFit.ReconstructedUnitAssignment: write_table
using Printf, Statistics, TOML

const UPPER = ("line_discontinuity", "stripe_periodicity", "fwd_bwd_nrmse")
const LOWER = ("contrast_to_noise",)

function quality_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        if k in ("--derive", "--size-check")
            o[k] = "true"; i += 1; continue
        end
        k in ("--catalog", "--config", "--out") && i < length(args) || error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    haskey(o, "--catalog") || error("Required: --catalog")
    haskey(o, "--derive") || (haskey(o, "--config") && haskey(o, "--out")) || error("Required: --config --out (or --derive)")
    haskey(o, "--out") && (ispath(o["--out"]) || islink(o["--out"])) && error("Output exists: $(o["--out"])")
    return o
end

function read_catalog(path)
    lines = [l for l in readlines(path) if !isempty(strip(l))]
    header = String.(split(lines[1], '\t'; keepempty=true))
    return header, [Dict(zip(header, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

value(r, k) = something(tryparse(Float64, get(r, k, "")), NaN)

"Robust bound on the log scale: exp(median +- 3 x 1.4826 x MAD of log(v))."
function robust_bound(v, side)
    x = log.(filter(t -> isfinite(t) && t > 0, v))
    m = median(x); s = 1.4826 * median(abs.(x .- m))
    return exp(side == :upper ? m + 3s : m - 3s)
end

function decide(rows, sel; size_check=false)
    out = Dict{String,String}[]
    for r in rows
        why = r["status"] == "usable" ? String[] : [r["status"]]
        if isempty(why)
            flags = split(get(r, "flags", ""), ","; keepempty=false)
            any(in(("roi_touches_border", "roi_touches_unscanned")), flags) && push!(why, "roi_at_edge")
            value(r, "second_object_ratio") < sel["max_second_object_ratio"] || push!(why, "several_objects")
            value(r, "fwd_bwd_corr") >= sel["min_fwd_bwd_corr"] || push!(why, "fwd_bwd_mismatch")
            for k in UPPER; value(r, k) <= sel["max_" * k] || push!(why, k); end
            for k in LOWER; value(r, k) >= sel["min_" * k] || push!(why, "low_" * k); end
        end
        push!(out, Dict("file" => r["file"], "sha256" => get(r, "sha256", ""), "session" => get(r, "session", ""),
                        "reasons_without_size" => join(why, ","), "decision_without_size" => isempty(why) ? "accept" : "reject",
                        "roi_length_nm" => get(r, "roi_length_nm", "")))
    end
    accepted = [value(Dict("x" => o["roi_length_nm"]), "x") for o in out if o["decision_without_size"] == "accept"]
    median_len = isempty(accepted) ? NaN : median(filter(isfinite, accepted))
    lo, hi = sel["size_window"]
    for o in out
        why = split(o["reasons_without_size"], ","; keepempty=false)
        if size_check && o["decision_without_size"] == "accept"
            len = value(Dict("x" => o["roi_length_nm"]), "x")
            lo * median_len <= len <= hi * median_len || (why = vcat(why, "size"))
        end
        o["reasons"] = join(why, ","); o["decision"] = isempty(why) ? "accept" : "reject"
    end
    return out, median_len
end

function main(args=ARGS)
    o = quality_options(args)
    o === nothing && return println("decide_scan_quality.jl --catalog TSV --config TOML --out TSV [--size-check] | --catalog TSV --derive")
    header, rows = read_catalog(o["--catalog"])
    usable = [r for r in rows if r["status"] == "usable"]
    if haskey(o, "--derive")
        println("# robust bounds over $(length(usable)) usable scans of $(o["--catalog"])")
        for k in UPPER; @printf("max_%s = %.4f\n", k, robust_bound([value(r, k) for r in usable], :upper)); end
        for k in LOWER; @printf("min_%s = %.4f\n", k, robust_bound([value(r, k) for r in usable], :lower)); end
        return
    end
    sel = TOML.parsefile(o["--config"])["selection"]
    out, median_len = decide(rows, sel; size_check=haskey(o, "--size-check"))
    metrics = ["status", "flags", "finite_fraction", "second_object_ratio", "fwd_bwd_corr", "fwd_bwd_nrmse",
               "line_discontinuity", "stripe_periodicity", "contrast_to_noise"]
    byfile = Dict(r["file"] * "|" * get(r, "sha256", "") => r for r in rows)
    for q in out, k in metrics; q[k] = get(byfile[q["file"] * "|" * q["sha256"]], k, ""); end
    write_table(o["--out"], vcat(["file", "sha256", "session", "decision", "reasons", "decision_without_size",
                                  "reasons_without_size", "roi_length_nm"], metrics), out)
    tally(k) = join(["$v $(count(q -> q[k] == v, out))" for v in ("accept", "reject")], ", ")
    println("Scans: $(length(out)); with size check: $(tally("decision")); without: $(tally("decision_without_size"))",
            haskey(o, "--size-check") ? @sprintf("; median ROI length of accepted scans %.2f nm", median_len) : "")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
