#!/usr/bin/env julia
# Leave-one-session-out application of the frozen unit assignment to a run of
# run_molecule_consensus_chitosan.jl (experiment F, docs/src/journal.md
# 2026-10-27). Label-free: sessions come from the SXM header date; no label,
# grade or benchmark membership is read.
#
# Usage:
#   julia --project=. scripts/pipeline/apply_frozen_assignment_loso.jl --run RUN_DIR --data-dir RAW_DIR \
#       --config config/unit_assignment_frozen_labelfree.toml --outdir NEW_DIR
#
# For each session, a frozen model is trained on the other sessions as
# build_frozen_assignment.jl does on a whole run (reference: their corroborated
# scans, at most reference_max_scans_per_track per molecule; emission rates on
# all their scans), then applied to the session's scans with fusion inside each
# molecule track. Writes NEW_DIR/predictions.tsv (the columns of
# apply_frozen_assignment.jl plus session) and NEW_DIR/sessions.tsv.
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment, STMFit.MoleculeConsensus
using STMSXMIO: read_sxm
using LinearAlgebra, Printf

function loso_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--data-dir", "--config", "--outdir") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--data-dir", "--config", "--outdir")) || error("Required: --run --data-dir --config --outdir")
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

function loso_main(args=ARGS)
    o = loso_options(args)
    o === nothing && return println("apply_frozen_assignment_loso.jl --run RUN_DIR --data-dir RAW_DIR --config TOML --outdir NEW_DIR")
    BLAS.set_num_threads(1)
    cfg = load_config(o["--config"]); st = load_frozen_settings(cfg)
    R = load_run_assignment(o["--run"], cfg)
    raw = raw_index(o["--data-dir"])
    session = Dict(f => scan_frame(read_sxm(raw[f])).session for f in unique(R.data.file))
    fused = fill(-1, length(R.data.file)); scan = fill(-1, length(R.data.file)); prob = fill(NaN, length(R.data.file))
    log = Dict{String,Any}[]
    for s in sort(unique(values(session)))
        held = findall(f -> session[f] == s, R.data.file)
        others = [r for r in R.consensus if session[basename(r["filepath"])] != s]
        ref = reference_files(others, st.max_per_track)
        ref_rows = findall(in(ref), R.data.file)
        theta_rows = findall(f -> session[f] != s, R.data.file)
        entry = Dict{String,Any}("session" => s, "scans" => length(unique(R.data.file[held])), "reference_scans" => length(ref))
        try
            model = train_reference_model(R.data, R.fisher, ref_rows, theta_rows, R.physical, cfg)
            f, p = apply_reference_model(model, R.data, R.fisher, held, R.physical, cfg)
            for (k, i) in enumerate(held)
                fused[i] = f[k]; scan[i] = p.call[k]; prob[i] = p.p_vote[k]
            end
            merge!(entry, Dict("theta0" => model.theta0, "theta1" => model.theta1, "status" => "ok"))
        catch e
            merge!(entry, Dict("theta0" => NaN, "theta1" => NaN, "status" => replace(sprint(showerror, e), '\t' => ' ', '\n' => ' ')))
        end
        push!(log, entry)
    end
    c(x) = x < 0 ? "?" : string(x)
    name = cfg["model"]["name"] * "+leave_one_session_out"
    rows = [Dict("file" => R.data.file[i], "lobe" => string(R.data.lobe[i]), "predicted" => c(fused[i]),
                 "scan_predicted" => c(scan[i]), "scan_probability_1" => isfinite(prob[i]) ? @sprintf("%.8f", prob[i]) : "NA",
                 "physical_lobe" => R.physical[i], "changed_by_fusion" => string(fused[i] != scan[i]),
                 "session" => session[R.data.file[i]], "model" => name) for i in eachindex(R.data.file)]
    mkpath(o["--outdir"])
    write_table(joinpath(o["--outdir"], "predictions.tsv"),
                ["file", "lobe", "predicted", "scan_predicted", "scan_probability_1", "physical_lobe", "changed_by_fusion", "session", "model"], rows)
    write_table(joinpath(o["--outdir"], "sessions.tsv"), ["session", "scans", "reference_scans", "theta0", "theta1", "status"],
                [Dict(k => string(v) for (k, v) in e) for e in log])
    @printf("Leave-one-session-out frozen assignment: %d sessions, %d failed, %d of %d lobes unclassified -> %s\n",
            length(log), count(e -> e["status"] != "ok", log), count(<(0), fused), length(fused), o["--outdir"])
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && loso_main()
