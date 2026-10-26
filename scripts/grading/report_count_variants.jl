#!/usr/bin/env julia
# Exact-N grade of count variants by molecule, with paired bootstrap intervals
# (grading side; campaign K, docs/src/journal.md 2026-10-26). Reads the
# variant_counts.tsv of scripts/data/replay_count_policy.jl and the benchmark
# manifest; molecule = track, session = first folder of the manifest path.
#
# Usage:
#   julia --project=. scripts/grading/report_count_variants.jl --counts variant_counts.tsv \
#       --manifest TOML --outdir NEW_DIR [--reference V0_consensus] [--draws 2000] [--seed 20261026]
#
# Every column ending in `_consensus`, and N_final, is graded on the manifest
# files: exact N by molecule (primary) and by scan, per session, and the paired
# difference against --reference.
include(joinpath(@__DIR__, "report_grade_uncertainty.jl"))

function variant_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--reference" => "V0_consensus", "--draws" => "2000", "--seed" => "20261026"); i = 1
    while i <= length(args)
        args[i] in ("--counts", "--manifest", "--outdir", "--reference", "--draws", "--seed") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    for k in ("--counts", "--manifest", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Scan grades of one count column: only exact N is filled (the other fields are neutral)."
function count_grades(rows, col, target, sessions)
    return [ScanGrade(r["file"], r["track"], get(sessions, r["file"], ""), parse(Int, r[col]) == target[r["file"]], 0, 1, 0, false)
            for r in rows if haskey(target, r["file"])]
end

function variant_main(args=ARGS)
    o = variant_options(args)
    o === nothing && return println("report_count_variants.jl --counts variant_counts.tsv --manifest TOML --outdir NEW_DIR [--reference V0_consensus]")
    draws = parse(Int, o["--draws"]); rng = MersenneTwister(parse(Int, o["--seed"]))
    files = TOML.parsefile(o["--manifest"])["files"]
    target = Dict(f => Int(v["target_N_for_human_validation"]) for (f, v) in files if haskey(v, "target_N_for_human_validation"))
    sessions = manifest_sessions(o["--manifest"])
    rows = read_tsv(o["--counts"])
    cols = vcat(["N_final"], sort([c for c in keys(rows[1]) if endswith(c, "_consensus")]))
    o["--reference"] in cols || error("Reference column $(o["--reference"]) not found")
    grades = Dict(c => count_grades(rows, c, target, sessions) for c in cols)
    mkpath(o["--outdir"])
    md = IOBuffer()
    println(md, "# Exact N of count variants, by molecule\n")
    println(md, "$(length(grades[cols[1]])) graded scans; a molecule counts once; 95% intervals, bootstrap over molecules ($draws draws).\n")
    println(md, "| count | exact N by molecule (%) | 95% interval | exact N by scan (%) | 95% interval |\n|---|---:|---|---:|---|")
    for c in cols
        _, groups = by_molecule(grades[c])
        mo = molecule_rates(groups); ci = bootstrap(molecule_rates, groups; draws, rng)
        sc = pooled(groups); cs = bootstrap(pooled, groups; draws, rng)
        println(md, "| $c | $(pct(mo.exact_n)) | $(interval_pct(ci.exact_n)) | $(pct(sc.exact_n)) | $(interval_pct(cs.exact_n)) |")
    end
    println(md, "\n## Paired differences against `$(o["--reference"])`\n")
    println(md, "| count | molecules | difference by molecule (points) | 95% interval | difference by scan | 95% interval |\n|---|---:|---:|---|---:|---|")
    ref = grades[o["--reference"]]
    for c in cols
        c == o["--reference"] && continue
        p = paired_groups(ref, grades[c])
        dm = paired_stat(p; weight=:molecule); ds = paired_stat(p)
        cim = bootstrap(q -> paired_stat(q; weight=:molecule), p; draws, rng); cis = bootstrap(q -> paired_stat(q), p; draws, rng)
        println(md, "| $c | $(length(p)) | $(pct(dm.exact_n)) | $(interval_pct(cim.exact_n)) | $(pct(ds.exact_n)) | $(interval_pct(cis.exact_n)) |")
    end
    println(md, "\n## By session (exact N by molecule, %)\n")
    println(md, "| session | molecules | scans | " * join(cols, " | ") * " |\n|---|---:|---:|" * repeat("---:|", length(cols)))
    for ses in sort(unique(x.session for x in grades[cols[1]]))
        cells = String[]
        nm = 0; ns = 0
        for c in cols
            _, groups = by_molecule([x for x in grades[c] if x.session == ses])
            nm, ns = length(groups), sum(length, groups)
            push!(cells, pct(molecule_rates(groups).exact_n))
        end
        println(md, "| $(isempty(ses) ? "(none)" : ses) | $nm | $ns | " * join(cells, " | ") * " |")
    end
    write(joinpath(o["--outdir"], "report.md"), take!(md))
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && variant_main()
