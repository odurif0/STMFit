#!/usr/bin/env julia
# Build a fully frozen unit-assignment model from a run of
# run_molecule_consensus_chitosan.jl (label-free).
#
# Usage:
#   julia --project=. scripts/pipeline/build_frozen_assignment.jl --run RUN_DIR \
#       --config config/unit_assignment_frozen_labelfree.toml --out MODEL.tsv
#
# The reference cohort follows the config's label-free rule (corroborated
# scans, at most reference_max_scans_per_track per molecule track, first in
# file-name order). The empirical Fisher per lobe parity, k-means, GMM and
# physical naming are learned on it; the fusion emission rates on the per-scan
# calls of every scan of the run. Writes MODEL.tsv (write_frozen_model) and
# MODEL.reference.tsv (file, track of every reference scan).
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment
using LinearAlgebra

function build_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--config", "--out") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        haskey(o, args[i]) && error("Repeated option $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--config", "--out")) || error("Required: --run --config --out")
    for p in (o["--out"], o["--out"] * ".reference.tsv")
        (ispath(p) || islink(p)) && error("Output exists: $p")
    end
    return o
end

function build_main(args=ARGS)
    o = build_options(args)
    o === nothing && return println("build_frozen_assignment.jl --run RUN_DIR --config TOML --out MODEL.tsv")
    BLAS.set_num_threads(1)
    cfg = load_config(o["--config"]); st = load_frozen_settings(cfg)
    R = load_run_assignment(o["--run"], cfg)
    ref = reference_files(R.consensus, st.max_per_track)
    ref_rows = findall(in(ref), R.data.file)
    model = train_reference_model(R.data, R.fisher, ref_rows, collect(eachindex(R.data.file)), R.physical, cfg)
    write_frozen_model(o["--out"], model)
    track = Dict(basename(r["filepath"]) => r["track"] for r in R.consensus)
    write_table(o["--out"] * ".reference.tsv", ["file", "track"],
                [Dict("file" => f, "track" => track[f]) for f in sort(collect(ref))])
    println("Frozen model: $(length(ref)) reference scans of $(length(unique(R.data.file))); theta0=$(model.theta0) ",
            "theta1=$(model.theta1) pi=$(model.pi) -> $(o["--out"])")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && build_main()
