#!/usr/bin/env julia
# Apply a frozen unit-assignment model to a run of
# run_molecule_consensus_chitosan.jl (label-free).
#
# Usage:
#   julia --project=. scripts/pipeline/apply_frozen_assignment.jl --run RUN_DIR --model MODEL.tsv \
#       --config config/unit_assignment_frozen_labelfree.toml --outdir NEW_DIR
#
# Recomputes emp_fisher with the model's parity models, scores every lobe of
# the run with the frozen clusters, and fuses the per-scan calls within each
# molecule track with the model's emission rates and prior. Nothing is learned
# from the run. Writes NEW_DIR/predictions.tsv: file, lobe, predicted (final),
# scan_predicted, scan_probability_1, physical_lobe, changed_by_fusion, model.
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment
using LinearAlgebra, Printf

function apply_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--model", "--config", "--outdir") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        haskey(o, args[i]) && error("Repeated option $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--model", "--config", "--outdir")) || error("Required: --run --model --config --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

function apply_main(args=ARGS)
    o = apply_options(args)
    o === nothing && return println("apply_frozen_assignment.jl --run RUN_DIR --model MODEL.tsv --config TOML --outdir NEW_DIR")
    BLAS.set_num_threads(1)
    cfg = load_config(o["--config"]); load_frozen_settings(cfg)
    R = load_run_assignment(o["--run"], cfg)
    model = read_frozen_model(o["--model"])
    rows = collect(eachindex(R.data.file))
    fused, p = apply_reference_model(model, R.data, R.fisher, rows, R.physical, cfg)
    c(x) = x < 0 ? "?" : string(x)
    name = cfg["model"]["name"]
    out = [Dict("file" => R.data.file[i], "lobe" => string(R.data.lobe[i]), "predicted" => c(fused[i]),
                "scan_predicted" => c(p.call[i]), "physical_lobe" => R.physical[i], "changed_by_fusion" => string(fused[i] != p.call[i]),
                "scan_probability_1" => isfinite(p.p_vote[i]) ? @sprintf("%.8f", p.p_vote[i]) : "NA", "model" => name) for i in rows]
    mkpath(o["--outdir"])
    write_table(joinpath(o["--outdir"], "predictions.tsv"),
                ["file", "lobe", "predicted", "scan_predicted", "scan_probability_1", "physical_lobe", "changed_by_fusion", "model"], out)
    @printf("Frozen assignment: %d lobes of %d scans, %d per-scan calls changed by fusion, %d unclassified -> %s\n",
            length(rows), length(unique(R.data.file)), count(i -> fused[i] != p.call[i], rows), count(<(0), fused),
            joinpath(o["--outdir"], "predictions.tsv"))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && apply_main()
