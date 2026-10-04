#!/usr/bin/env julia
# Verify, then externally grade, one run of test/run_molecule_consensus_chitosan.jl.
#
# Usage:
#   julia --project=. test/grade_consensus_run.jl --run RUN_DIR --outdir NEW_DIR \
#       [--data-dir RAW_DIR] [--verify-only] [--profile NAME=PREDICTIONS.tsv ...]
#
# 1. Label-free integrity checks on the run (cohort, per-stage counts,
#    consensus provenance, fusion integrity; raw SXM hashes when --data-dir is
#    given). Any failure stops before grading.
# 2. External grading, outside the method: benchmark membership filtering,
#    report_unit_assignment_benchmark.jl --full145-own-n for the final (fused)
#    and per-scan calls plus any --profile, and grade_chitosan_benchmark.jl for
#    the per-scan and final counts. Labels are read only here.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using SHA, TOML

const ROOT = dirname(@__DIR__)
const MANIFEST = joinpath(ROOT, "benchmarks", "chitosan_6mer_counting_confirmed.toml")
const KNOWN_RULES = Set(["agrees", "consensus_applied", "track_too_short", "no_strict_majority", "footprint_outside_frame"])
const CHANGED_RULES = Set(["consensus_applied", "consensus_registered_roi"])

function grade_options(args)
    ("--help" in args || "-h" in args) && return nothing
    o = Dict{String,String}(); profiles = Pair{String,String}[]; i = 1
    while i <= length(args)
        k = args[i]
        if k == "--verify-only"; o[k] = "true"; i += 1; continue; end
        k in ("--run", "--outdir", "--data-dir", "--profile") && i < length(args) || error("Unknown option or missing value: $k")
        if k == "--profile"
            name, path = split(args[i+1], '='; limit=2)
            push!(profiles, String(name) => String(path))
        else
            haskey(o, k) && error("Repeated option $k")
            o[k] = args[i+1]
        end
        i += 2
    end
    haskey(o, "--run") && isdir(o["--run"]) || error("--run RUN_DIR is required")
    haskey(o, "--verify-only") || haskey(o, "--outdir") || error("--outdir NEW_DIR is required")
    haskey(o, "--outdir") && (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o, profiles
end

rows_per_file(rows) = (d = Dict{String,Int}(); for r in rows; f = basename(r["file"]); d[f] = get(d, f, 0) + 1; end; d)

"""Label-free integrity checks. Returns true when every check passes."""
function verify_run(run; data_dir="")
    ok = true
    check(c, msg) = (println(c ? "PASS " : "FAIL ", msg); ok &= c)
    tab(p) = last(read_table(joinpath(run, p)))
    counting = tab("counting_summary.tsv")
    nscan = Dict(basename(r["filepath"]) => parse(Int, r["N_selected"]) for r in counting)
    files = Set(keys(nscan))
    check(length(files) == length(counting), "counting has one row per scan ($(length(files)) scans)")
    check(all(r["status"] == "ok" for r in counting), "all counting rows ok")
    check(rows_per_file(tab("scan_geometry.tsv")) == nscan, "scan geometry rows equal per-scan N")
    cons = tab("consensus/consensus.tsv")
    final = tab("consensus/consensus_summary_final.tsv")
    nfinal = Dict(basename(r["filepath"]) => parse(Int, r["N_selected"]) for r in final)
    rule = Dict(basename(r["filepath"]) => r["count_rule"] for r in final)
    check(Set(r["file"] for r in cons) == files && Set(keys(nfinal)) == files, "consensus covers the cohort")
    check(all(parse(Int, r["N_scan"]) == nscan[r["file"]] for r in cons), "consensus N_scan equals counting N")
    check(all(r["rule"] in KNOWN_RULES for r in cons), "known consensus rules only")
    check(all(nfinal[f] == nscan[f] || rule[f] in CHANGED_RULES for f in files), "only consensus rows change N")
    check(all(nfinal[f] == nscan[f] for f in files if rule[f] == "consensus_fit_failed"), "failed refits keep the per-scan N")
    scan_pred = tab("assignment/predictions.tsv")
    check(rows_per_file(scan_pred) == nfinal, "per-scan predictions rows equal final N")
    check(rows_per_file(tab("assignment/features.tsv")) == nfinal, "assignment geometry rows equal final N")
    check(all(r["predicted"] in ("0", "1", "?") for r in scan_pred), "per-scan calls in {0,1,?}")
    check(!isfile(joinpath(run, "failures.tsv")), "no failures.tsv")
    if isdir(joinpath(run, "fusion"))
        fused = tab("fusion/predictions_fused.tsv")
        check(read(joinpath(run, "predictions.tsv")) == read(joinpath(run, "fusion", "predictions_fused.tsv")),
              "final predictions.tsv equals the fusion output")
        orig = Dict((r["file"], r["lobe"]) => r["predicted"] for r in scan_pred)
        check(all(r["scan_predicted"] == orig[(r["file"], r["lobe"])] for r in fused), "scan_predicted equals the per-scan call")
        isfused(r) = endswith(r["model"], "+latent_class_fusion")
        p = only(tab("fusion/fusion_params.tsv"))
        check(count(isfused, fused) == parse(Int, p["rows_fused"]), "fused row count matches fusion_params.tsv")
        check(parse(Float64, p["theta1"]) > parse(Float64, p["theta0"]), "detection rate exceeds false-call rate")
        check(all(r["predicted"] == r["scan_predicted"] for r in fused if !isfused(r)), "unfused rows keep per-scan calls")
        println("fusion: pi=$(p["pi"]) theta0=$(p["theta0"]) theta1=$(p["theta1"]) rows_fused=$(p["rows_fused"]) ",
                "calls_changed=$(count(r -> isfused(r) && r["predicted"] != r["scan_predicted"], fused))")
    end
    if isfile(joinpath(run, "positions.tsv"))
        pos = tab("positions.tsv")
        check(rows_per_file(pos) == nfinal, "positions rows equal final N")
        check(all(r["sd_source"] in ("track", "cohort", "none") for r in pos), "known position-uncertainty sources")
    end
    if !isempty(data_dir)
        raw = Dict{String,String}()
        for (d, _, names) in walkdir(data_dir), n in names
            endswith(lowercase(n), ".sxm") && (raw[n] = joinpath(d, n))
        end
        hashes = tab("raw_hashes.tsv")
        check(length(hashes) == length(files) &&
              all(haskey(raw, r["file"]) && bytes2hex(open(sha256, raw[r["file"]])) == r["sha256"] for r in hashes),
              "raw SXM hashes match $data_dir")
    end
    tally = Dict{String,Int}(); for f in files; tally[rule[f]] = get(tally, rule[f], 0) + 1; end
    println("final count rules: ", tally)
    return ok
end

function membership_filter(predictions, out, files)
    header, rows = read_table(predictions)
    present = Set(basename(r["file"]) for r in rows)
    isempty(setdiff(files, present)) || error("$predictions misses benchmark files; no partial grade")
    write_table(out, header, [r for r in rows if basename(r["file"]) in files])
end

function grade_run(run, outdir, profiles)
    files = Set(String.(keys(TOML.parsefile(MANIFEST)["files"])))
    inputs = joinpath(outdir, "inputs"); mkpath(inputs)
    sources = vcat(["final" => joinpath(run, "predictions.tsv"), "per_scan" => joinpath(run, "assignment", "predictions.tsv")], profiles)
    args = ["--full145-own-n", "--outdir", joinpath(outdir, "units")]
    for (name, path) in sources
        filtered = joinpath(inputs, name * ".tsv")
        membership_filter(path, filtered, files)
        append!(args, ["--profile", name * "=" * filtered])
    end
    julia = Base.julia_cmd()
    cd(() -> run_cmd(`$julia --project=$ROOT $(joinpath(ROOT, "test", "report_unit_assignment_benchmark.jl")) $args`), ROOT)
    for (name, summary) in (("final", joinpath(run, "consensus", "consensus_summary_final.tsv")),
                            ("per_scan", joinpath(run, "counting_summary.tsv")))
        cd(() -> run_cmd(`$julia --project=$ROOT $(joinpath(ROOT, "test", "grade_chitosan_benchmark.jl")) --manifest $MANIFEST
                          --results $summary --column N_selected --out $(joinpath(outdir, "counts_" * name * ".tsv"))`), ROOT)
    end
end
run_cmd(c) = Base.run(c)

function main(args=ARGS)
    parsed = grade_options(args)
    if parsed === nothing
        println("grade_consensus_run.jl --run RUN_DIR --outdir NEW_DIR [--data-dir RAW_DIR] [--verify-only] [--profile NAME=TSV ...]")
        return
    end
    o, profiles = parsed
    verify_run(o["--run"]; data_dir=get(o, "--data-dir", "")) || error("Integrity checks failed; no grade emitted")
    haskey(o, "--verify-only") && return println("ALL PASS")
    grade_run(abspath(o["--run"]), abspath(o["--outdir"]), profiles)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
