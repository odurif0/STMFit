#!/usr/bin/env julia
# Raw SXM -> per-scan GCV counts -> repeated-scan molecule consensus -> 0/1/?.
#
# Usage:
#   julia -t 4 --project=. scripts/pipeline/run_molecule_consensus_chitosan.jl \
#       --data-dir RAW_DIR --count-config config/chitosan_decoupled_support.toml \
#       --config config/unit_assignment_corroborated_training.toml \
#       --consensus-config config/molecule_consensus.toml \
#       --templates dft/molds/chitosan_cc_molds_native_v1.tsv --outdir NEW_DIR [--selected-summary TSV] \
#       [--count-shards N] [--dry-run]
#
# --count-shards N runs the per-scan counting as N single-thread processes
# (default: min(scans, 4)), for a full cluster node. Every scan is counted
# independently, so the shard count does not change any result.
# If the counting config sets model.geometry_support_threshold_rule to another
# rule than support_threshold_rule, N is counted on the counting support and
# every later fixed-N refit uses the geometry support, or the counting support
# for scans whose fit fails there (geometry_summary.tsv records the choice).
# --selected-summary reuses a label-free per-scan counting summary instead of
# counting again (development only). Every other stage is recomputed from raw
# pixels. No benchmark label, expected count, sequence, grade or reference
# prediction is an input. External grading remains a separate step.
module MoleculeConsensusRun
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
using STMFit.MoleculeConsensus: load_consensus_settings
using SHA, Statistics

const REQUIRED = Set(["--data-dir", "--count-config", "--config", "--consensus-config", "--templates", "--outdir"])
const OPTIONAL = Set(["--selected-summary", "--count-shards"])

function run_options(args)
    ("--help" in args || "-h" in args) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        haskey(o, k) && error("Repeated option $k")
        if k == "--dry-run"; o[k] = "true"; i += 1; continue; end
        (k in REQUIRED || k in OPTIONAL) && i < length(args) && !startswith(args[i+1], "--") ||
            error("Missing/forbidden option $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in REQUIRED) || error("Required: $(join(sort(collect(REQUIRED)), ", "))")
    isdir(o["--data-dir"]) || error("Missing raw directory")
    for k in ("--count-config", "--config", "--consensus-config", "--templates", "--selected-summary")
        haskey(o, k) && !isfile(o[k]) && error("Missing input $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    haskey(o, "--count-shards") && something(tryparse(Int, o["--count-shards"]), 0) < 1 && error("--count-shards must be a positive integer")
    return o
end

const COUNT_POLICIES = ("support_midpoint_hybrid", "adaptive_support_rescue", "support_midpoint")

function merge_count_shards(paths, files, output; policy="support_midpoint_hybrid")
    header = String[]; rows = Dict{String,String}[]; seen = Set{String}()
    for path in paths
        h, rs = read_table(path)
        isempty(header) && (header = h)
        h == header || error("Counting shard schemas differ")
        counts = selected_counts(path)
        isempty(intersect(seen, keys(counts))) || error("Duplicate scan across counting shards")
        union!(seen, keys(counts))
        all(get(r, "selection_policy", "") == policy for r in rs) || error("Wrong counting policy")
        policy == "support_midpoint_hybrid" && any(startswith(get(r, k, ""), "adaptive_support") for r in rs for k in
            ("refined_policy", "refined_source", "selection_source")) && error("Adaptive support is not part of the hybrid policy")
        append!(rows, rs)
    end
    seen == Set(files) || error("Counting omitted or added scans")
    sort!(rows; by=r -> basename(r["filepath"]))
    write_table(output, header, rows)
    return selected_counts(output)
end

const SUPPORT_RULES = ("profile_quantile_noise", "half_maximum_cap")

"Counting and geometry support rules of a counting config (geometry defaults to counting)."
function support_rules(model)
    count_rule = String(get(model, "support_threshold_rule", "profile_quantile_noise"))
    geometry_rule = String(get(model, "geometry_support_threshold_rule", count_rule))
    count_rule in SUPPORT_RULES && geometry_rule in SUPPORT_RULES || error("Unknown support threshold rule")
    return count_rule, geometry_rule
end

"Copy of a counting summary with a per-scan `support_threshold_rule` column."
function write_support_rules(summary, output, rule_of)
    header, rows = read_table(summary)
    "support_threshold_rule" in header || push!(header, "support_threshold_rule")
    for r in rows; r["support_threshold_rule"] = rule_of(basename(r["filepath"])); end
    return write_table(output, header, rows)
end

"""
Decoupled supports (`geometry_support_threshold_rule`): `initial` holds every
scan refit at its counted N on the geometry support. A scan whose fit there is
missing or flagged is refit on the counting support instead. The returned
`geometry_summary.tsv` records each scan's rule, which every later fixed-N refit
replays. Only fit success decides. Returns (summary, geometry, moved scans).
"""
function geometry_fallback(out, raw, count_config, selected, initial, counts, rules, runner)
    count_rule, geometry_rule = rules
    _, lobes = lobe_table(initial)
    fitted = Dict{String,Int}(); flagged = Set{String}()
    for ((f, _), r) in lobes
        fitted[f] = get(fitted, f, 0) + 1
        endswith(r["source"], "_flagged") && push!(flagged, f)
    end
    moved = sort([f for (f, n) in counts if get(fitted, f, 0) != n || f in flagged])
    summary = write_support_rules(selected, joinpath(out, "geometry_summary.tsv"),
                                  f -> f in moved ? count_rule : geometry_rule)
    geometry = joinpath(out, "scan_geometry.tsv")
    isempty(moved) && return summary, cp(initial, geometry), moved
    part = joinpath(out, "scan_geometry_count_support.tsv")
    runner(out, "scan_geometry_count_support", "extract_lobe_features.jl", ["--config", count_config, "--data-dir", raw,
        "--selected-summary", summary, "--files", join(moved, ","), "--out", part]; threads=1)
    header, kept = read_table(initial)
    part_header, refit = read_table(part)
    part_header == header || error("Fallback geometry schema differs")
    rows = vcat([r for r in kept if !(r["file"] in moved)], refit)
    sort!(rows; by=r -> (r["file"], parse(Int, r["lobe"])))
    return summary, write_table(geometry, header, rows), moved
end

"""
A consensus count is kept only if the unchanged fixed-N fitter can fit that
scan at it: first on the scan's own support, then once on the registered
molecule region (`registered_roi`, rule `consensus_registered_roi`).
Otherwise the scan keeps its own GCV count (`consensus_fit_failed`).
Only changed scans are refit here; no label or grade is involved.
"""
function check_consensus_refits(out, raw, count_config, proposed, runner)
    header, rows = read_table(proposed)
    if !("refit_support" in header)
        push!(header, "refit_support"); for r in rows; r["refit_support"] = "scan"; end
    end
    dir = dirname(proposed)
    target(r) = parse(Int, r["N_selected"])
    function refit(files, summary, name)
        fitted = Dict{String,Int}()
        isempty(files) && return fitted
        check = joinpath(dir, name * ".tsv")
        runner(out, name, "extract_lobe_features.jl", ["--config", count_config, "--data-dir", raw,
            "--selected-summary", summary, "--files", join(sort(files), ","), "--out", check]; threads=1)
        _, fr = read_table(check)
        for r in fr; fitted[r["file"]] = get(fitted, r["file"], 0) + 1; end
        return fitted
    end
    byfile = Dict(basename(r["filepath"]) => r for r in rows)
    changed = [f for (f, r) in byfile if r["N_selected"] != r["N_scan"]]
    ok = refit(changed, proposed, "consensus_refit_check")
    failed = [f for f in changed if get(ok, f, 0) != target(byfile[f])]
    if !isempty(failed)
        for f in failed; byfile[f]["refit_support"] = "registered_roi"; end
        retry = joinpath(dir, "consensus_summary_registered_roi.tsv")
        write_table(retry, header, rows)
        ok2 = refit(failed, retry, "consensus_refit_registered_roi")
        for f in failed
            r = byfile[f]
            if get(ok2, f, 0) == target(r)
                r["count_rule"] = "consensus_registered_roi"
            else
                r["N_selected"] = r["N_scan"]; r["count_rule"] = "consensus_fit_failed"; r["refit_support"] = "scan"
            end
        end
    end
    final = joinpath(dir, "consensus_summary_final.tsv")
    write_table(final, header, rows)
    return final
end

function write_chain_report(consensus_summary, predictions, output)
    _, cons = read_table(consensus_summary)
    _, pred = lobe_table(predictions; required=["predicted", "confidence"])
    rows = Dict{String,String}[]
    for c in sort(cons; by=r -> r["filepath"])
        f = c["filepath"]; keys_f = sort([k for k in keys(pred) if first(k) == f]; by=last)
        labels = [pred[k]["predicted"] for k in keys_f]
        conf = [parse(Float64, pred[k]["confidence"]) for k in keys_f if pred[k]["predicted"] != "?"]
        push!(rows, Dict("file" => f, "N_scan" => c["N_scan"], "N_selected" => c["N_selected"],
            "count_rule" => c["count_rule"], "track" => c["track"], "track_size" => c["track_size"],
            "count_agreement" => c["count_agreement"], "assignment" => join(labels),
            "uncertain" => string(count(==("?"), labels)),
            "min_confidence" => isempty(conf) ? "NaN" : @sprintf("%.4f", minimum(conf)),
            "mean_confidence" => isempty(conf) ? "NaN" : @sprintf("%.4f", mean(conf))))
    end
    write_table(output, ["file", "N_scan", "N_selected", "count_rule", "track", "track_size", "count_agreement",
                         "assignment", "uncertain", "min_confidence", "mean_confidence"], rows)
end

function execute_run(o; runner=run_stage, thread_budget=min(Threads.nthreads(), 4))
    VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
    thread_budget > 0 || error("Positive thread budget required")
    cfg = TOML.parsefile(o["--count-config"])
    cfg["model"]["selection_criterion"] == "gcv" || error("GCV required")
    policy = cfg["model"]["selection_policy"]
    policy in COUNT_POLICIES || error("Counting policy must be one of $(join(COUNT_POLICIES, ", "))")
    rules = support_rules(cfg["model"]); decoupled = rules[2] != rules[1]
    assignment = load_config(o["--config"])
    assignment["preprocessing"]["patch_residual_filter"] == "smooth_residual" || error("Symmetric residuals required")
    load_consensus_settings(TOML.parsefile(o["--consensus-config"]))
    rawfiles = raw_index(o["--data-dir"]); files = sort(collect(keys(rawfiles)))
    nchunks = haskey(o, "--count-shards") ? min(length(files), parse(Int, o["--count-shards"])) : min(length(files), thread_budget, 4)
    tag = @sprintf("%03d", round(Int, 100cfg["model"]["max_overlap"]))
    reuse = haskey(o, "--selected-summary")
    println("$(length(files)) raw scans; ", reuse ? "reused per-scan counts" : "$nchunks counting shards",
            " -> scan geometry -> molecule consensus -> assignment -> molecule fusion (if configured) -> positions -> plots")
    decoupled && println("Counting support: $(rules[1]); geometry support: $(rules[2]), counting support where that fit fails")
    haskey(o, "--dry-run") && return println("Dry-run: no pixels read, no output written; no labels or grades are inputs")
    out = abspath(o["--outdir"]); mkpath(joinpath(out, "logs")); stage = "inputs"
    try
        raw = prepare_raw(o["--data-dir"], out, Set(files))
        inputs = [k for k in ("--count-config", "--config", "--consensus-config", "--templates", "--selected-summary") if haskey(o, k)]
        write_table(joinpath(out, "input_hashes.tsv"), ["input", "path", "sha256"],
            [Dict("input" => k, "path" => abspath(o[k]), "sha256" => bytes2hex(sha256(read(o[k])))) for k in inputs])
        write_table(joinpath(out, "raw_hashes.tsv"), ["file", "sha256"],
            [Dict("file" => f, "sha256" => bytes2hex(sha256(read(rawfiles[f])))) for f in files])
        selected = joinpath(out, "counting_summary.tsv")
        if reuse
            cp(o["--selected-summary"], selected)
            Set(keys(selected_counts(selected))) == Set(files) || error("Reused counts do not cover the raw cohort")
        else
            stage = "counting"
            paths = String[]
            for i in 1:nchunks
                name = nchunks == 1 ? "summary_overlap$(tag)_hard.tsv" :
                    @sprintf("summary_overlap%s_hard_chunk%02dof%02d.tsv", tag, i, nchunks)
                push!(paths, joinpath(out, "counting_chunk$i", name))
            end
            @sync for i in 1:nchunks
                Threads.@spawn runner(out, "counting_chunk$i", "batch_full.jl",
                    [string(length(files)), "--config", abspath(o["--count-config"]), "--data-dir", raw,
                     "--outdir", dirname(paths[i]), "--tsv", joinpath(out, "no_triage_input.tsv"),
                     "--skip-1d", "--chunk", "$i/$nchunks"]; threads=1)
            end
            merge_count_shards(paths, files, selected; policy)
        end
        counts = selected_counts(selected)
        stage = "scan_geometry"
        support_summary = selected
        if decoupled
            initial_summary = write_support_rules(selected, joinpath(out, "geometry_summary_initial.tsv"), _ -> rules[2])
            initial = export_features(out, stage, "extract_lobe_features.jl", ["--config", abspath(o["--count-config"]),
                "--data-dir", raw, "--selected-summary", initial_summary], joinpath(out, "scan_geometry_initial.tsv");
                nfiles=length(files), thread_budget, runner)
            support_summary, geometry, moved = geometry_fallback(out, raw, abspath(o["--count-config"]), selected,
                initial, counts, rules, runner)
            println("Geometry on the counting support for $(length(moved)) of $(length(files)) scans")
        else
            geometry = joinpath(out, "scan_geometry.tsv")
            export_features(out, stage, "extract_lobe_features.jl", ["--config", abspath(o["--count-config"]),
                "--data-dir", raw, "--selected-summary", selected], geometry; nfiles=length(files),
                thread_budget, runner)
        end
        check_counts(geometry, counts)
        stage = "consensus"
        consensus_dir = joinpath(out, "consensus")
        runner(out, stage, "build_molecule_consensus.jl", ["--data-dir", raw, "--selected-summary", support_summary,
            "--features", geometry, "--config", abspath(o["--consensus-config"]), "--outdir", consensus_dir]; threads=1)
        proposed = joinpath(consensus_dir, "consensus_summary.tsv")
        stage = "consensus_refit_check"
        consensus = check_consensus_refits(out, raw, abspath(o["--count-config"]), proposed, runner)
        final_counts = selected_counts(consensus)
        Set(keys(final_counts)) == Set(files) || error("Consensus omitted or added scans")
        sel = assignment["selection"]
        stage = "assignment"
        dest = joinpath(out, "assignment")
        runner(out, stage, "run_reconstructed_chitosan.jl",
            ["--data-dir", raw, "--count-config", abspath(o["--count-config"]),
             "--config", abspath(o["--config"]), "--templates", abspath(o["--templates"]),
             "--selected-summary", consensus, "--outdir", dest,
             (get(sel, "assignment_training_scans", "all") == "corroborated_counts" ? ["--training-scans", consensus] : String[])...];
            threads=min(thread_budget, 4))
        for name in ("features.tsv", "features_split.tsv", "predictions.tsv")
            check_counts(joinpath(dest, name), final_counts)
        end
        final_predictions = joinpath(dest, "predictions.tsv")
        if get(TOML.parsefile(o["--consensus-config"])["selection"], "fusion", "none") == "latent_class"
            stage = "fusion"
            fusion_dir = joinpath(out, "fusion")
            runner(out, stage, "build_molecule_fusion.jl", ["--data-dir", raw, "--consensus", consensus_dir,
                "--features", joinpath(dest, "features.tsv"), "--predictions", final_predictions,
                "--config", abspath(o["--consensus-config"]), "--outdir", fusion_dir]; threads=1)
            final_predictions = joinpath(fusion_dir, "predictions_fused.tsv")
            check_counts(final_predictions, final_counts)
        end
        cp(final_predictions, joinpath(out, "predictions.tsv"))
        stage = "report"
        write_chain_report(consensus, final_predictions, joinpath(out, "chain_report.tsv"))
        # Label-free position uncertainty from repeated scans (does not change any earlier output).
        stage = "positions"
        runner(out, stage, "report_position_repeatability.jl", ["--run", out, "--data-dir", raw,
            "--outdir", joinpath(out, "position_repeatability")]; threads=1)
        cp(joinpath(out, "position_repeatability", "positions.tsv"), joinpath(out, "positions.tsv"))
        # One label-free plot per molecule scan (final calls, positions, consensus outcome).
        stage = "plots"
        runner(out, stage, "plot_final_molecules.jl", ["--run", out, "--data-dir", raw,
            "--config", abspath(o["--count-config"]), "--outdir", joinpath(out, "plots")]; threads=1)
    catch e
        write_table(joinpath(out, "failures.tsv"), ["stage", "reason"],
            [Dict("stage" => stage, "reason" => replace(sprint(showerror, e), '\n' => ' ', '\t' => ' '))])
        rethrow()
    end
    println("Raw-to-prediction molecule-consensus run complete: $out. External grading remains separate.")
end

function main(args=ARGS)
    o = run_options(args)
    o === nothing ? println("run_molecule_consensus_chitosan.jl --data-dir DIR --count-config TOML --config TOML --consensus-config TOML --templates TSV --outdir NEW_DIR [--selected-summary TSV] [--count-shards N] [--dry-run]") :
        execute_run(o)
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && MoleculeConsensusRun.main()
