#!/usr/bin/env julia
# Compare the automatic quality decision (decide_scan_quality.jl) with the
# earlier human triage decisions on the same scans. Grading side: reads the
# triage table; no benchmark grade is involved.
#
# Usage:
#   julia --project=. scripts/grading/report_scan_quality.jl --quality QUALITY_TSV \
#       --triage benchmarks/chitosan_6mer_cu100_triage.tsv --outdir NEW_DIR [--ignore-motif vue_large]
#
# Scans are matched by sha256. Human kept = non-empty set; human rejected =
# empty set, with its motif. Rows whose motif is ignored (default vue_large, an
# automatic overview rejection) are left out. Writes report.md, confusion.tsv,
# by_motif.tsv (share of each human rejection reason the rule rejects) and
# kept_rejected.tsv (human-kept scans the rule rejects, with its reasons).
using STMFit.ReconstructedUnitAssignment: write_table
using Printf

function options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--ignore-motif" => "vue_large"); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--quality", "--triage", "--outdir", "--ignore-motif") && i < length(args) || error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--quality", "--triage", "--outdir")) || error("Required: --quality --triage --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

function tsv(path)
    lines = [l for l in readlines(path) if !isempty(strip(l)) && !startswith(l, "#")]
    header = String.(split(lines[1], '\t'; keepempty=true))
    return [Dict(zip(header, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

pct(a, b) = b == 0 ? "NA" : @sprintf("%.1f%%", 100a / b)

function report(o)
    quality = Dict(q["sha256"] => q for q in tsv(o["--quality"]))
    ignored = Set(split(o["--ignore-motif"], ","; keepempty=false))
    rows = Dict{String,String}[]
    for t in tsv(o["--triage"])
        t["motif"] in ignored && isempty(t["set"]) && continue
        q = get(quality, t["sha256"], nothing)
        q === nothing && error("No quality row for $(t["file"]) ($(t["sha256"]))")
        push!(rows, Dict("file" => t["file"], "human" => isempty(t["set"]) ? "rejected" : "kept",
                         "motif" => isempty(t["set"]) ? t["motif"] : "", "auto" => q["decision"],
                         "auto_without_size" => q["decision_without_size"], "reasons" => q["reasons"]))
    end
    out = mkpath(o["--outdir"])
    conf = Dict{String,String}[]
    for col in ("auto", "auto_without_size"), h in ("kept", "rejected")
        n = count(r -> r["human"] == h, rows)
        a = count(r -> r["human"] == h && r[col] == "accept", rows)
        push!(conf, Dict("decision" => col, "human" => h, "scans" => string(n), "auto_accept" => string(a),
                         "auto_reject" => string(n - a)))
    end
    write_table(joinpath(out, "confusion.tsv"), ["decision", "human", "scans", "auto_accept", "auto_reject"], conf)
    motifs = sort(unique(r["motif"] for r in rows if r["human"] == "rejected"))
    bym = [Dict("motif" => m, "scans" => string(count(r -> r["motif"] == m, rows)),
                "auto_reject" => string(count(r -> r["motif"] == m && r["auto"] == "reject", rows)),
                "auto_reject_without_size" => string(count(r -> r["motif"] == m && r["auto_without_size"] == "reject", rows)))
           for m in motifs]
    write_table(joinpath(out, "by_motif.tsv"), ["motif", "scans", "auto_reject", "auto_reject_without_size"], bym)
    kr = sort([r for r in rows if r["human"] == "kept" && r["auto"] == "reject"]; by=r -> r["file"])
    write_table(joinpath(out, "kept_rejected.tsv"), ["file", "reasons"], kr)
    reasons = Dict{String,Int}()
    for r in kr, x in split(r["reasons"], ","; keepempty=false); reasons[x] = get(reasons, x, 0) + 1; end
    agree(col) = count(r -> (r["human"] == "kept") == (r[col] == "accept"), rows)
    open(joinpath(out, "report.md"), "w") do io
        println(io, "# Automatic quality decision against the human triage\n")
        println(io, "$(length(rows)) scans matched by sha256 (ignored motif: $(o["--ignore-motif"])).\n")
        println(io, "| decision | human | scans | auto accept | auto reject |\n|---|---|---:|---:|---:|")
        for c in conf; println(io, "| $(c["decision"]) | $(c["human"]) | $(c["scans"]) | $(c["auto_accept"]) | $(c["auto_reject"]) |"); end
        println(io, "\nAgreement: $(agree("auto")) / $(length(rows)) ($(pct(agree("auto"), length(rows)))) with the size check, ",
                "$(agree("auto_without_size")) / $(length(rows)) without.\n")
        println(io, "## Human rejection reasons caught by the rule\n\n| motif | scans | auto reject | without size |\n|---|---:|---:|---:|")
        for m in bym
            println(io, "| $(m["motif"]) | $(m["scans"]) | $(m["auto_reject"]) ($(pct(parse(Int, m["auto_reject"]), parse(Int, m["scans"])))) | $(m["auto_reject_without_size"]) |")
        end
        println(io, "\n## Human-kept scans the rule rejects: $(length(kr))\n\n| reason | scans |\n|---|---:|")
        for (k, v) in sort(collect(reasons); by=last, rev=true); println(io, "| $k | $v |"); end
    end
    print(read(joinpath(out, "report.md"), String))
end

function main(args=ARGS)
    o = options(args)
    o === nothing ? println("report_scan_quality.jl --quality TSV --triage TSV --outdir NEW_DIR [--ignore-motif vue_large]") : report(o)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
