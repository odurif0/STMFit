#!/usr/bin/env julia
# Latent-class fusion of per-scan unit calls over repeated scans of one molecule.
#
# Usage:
#   julia --project=. test/build_molecule_fusion.jl --data-dir RAW_DIR \
#       --consensus CONSENSUS_DIR --features assignment/features.tsv \
#       --predictions assignment/predictions.tsv \
#       --config config/molecule_consensus.toml --outdir NEW_DIR
#
# Outputs predictions_fused.tsv (every lobe; fused lobes carry the posterior
# as probability_1 and the per-scan call in scan_predicted), physical_lobes.tsv
# and fusion_params.tsv. Label-free: no benchmark label, sequence or count.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
include(joinpath(@__DIR__, "lib", "molecule_consensus.jl"))
include(joinpath(@__DIR__, "lib", "molecule_fusion.jl"))
using .ReconstructedUnitAssignment, .MoleculeConsensus, .MoleculeFusion
using STMSXMIO, TOML, Printf

const FUSION_OPTIONS = Set(["--data-dir", "--consensus", "--features", "--predictions", "--config", "--outdir"])

function fusion_options(args)
    ("--help" in args || "-h" in args) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        k in FUSION_OPTIONS && i < length(args) || error("Unknown option or missing value: $k")
        haskey(o, k) && error("Repeated option $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in FUSION_OPTIONS) || error("Required: $(join(sort(collect(FUSION_OPTIONS)), ", "))")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

function build_fusion(o)
    st = load_fusion_settings(TOML.parsefile(o["--config"]))
    _, cons = read_table(joinpath(o["--consensus"], "consensus.tsv"))
    _, final = read_table(joinpath(o["--consensus"], "consensus_summary_final.tsv"))
    tracks = Dict(r["file"] => parse(Int, r["track"]) for r in cons)
    counts = Dict(basename(r["filepath"]) => parse(Int, r["N_selected"]) for r in final)
    Set(keys(tracks)) == Set(keys(counts)) || error("Consensus tables disagree on the cohort")
    _, geo = lobe_table(o["--features"]; required=["x_nm", "y_nm"])
    _, pred = lobe_table(o["--predictions"]; required=["predicted", "probability_1", "confidence"])
    require_same_keys(geo, pred, "predictions")
    raw = Dict{String,String}()
    for (d, _, names) in walkdir(o["--data-dir"]), n in names
        endswith(lowercase(n), ".sxm") && (raw[n] = joinpath(d, n))
    end
    frames = Dict(f => scan_frame(read_sxm(raw[f])) for f in keys(counts))
    lobes_abs = Dict{String,Vector{NTuple{2,Float64}}}()
    for key in sort(collect(keys(geo)))
        r = geo[key]
        push!(get!(lobes_abs, key[1], NTuple{2,Float64}[]),
              absolute_xy(frames[key[1]], parse(Float64, r["x_nm"]), parse(Float64, r["y_nm"])))
    end
    all(length(lobes_abs[f]) == counts[f] for f in keys(counts)) || error("Geometry rows differ from final counts")
    order = sort(collect(keys(counts)); by=f -> (frames[f].acquired, f))
    phys = physical_lobes(order, tracks, counts, lobes_abs)
    calls = Dict{String,Vector{Int}}()
    for (key, r) in pred
        haskey(phys, key) || continue
        r["predicted"] in ("0", "1") || continue
        push!(get!(calls, phys[key], Int[]), r["predicted"] == "1" ? 1 : 0)
    end
    fused_keys = sort([k for (k, v) in calls if length(v) >= st.min_scans])
    isempty(fused_keys) && error("No physical lobe observed in at least $(st.min_scans) scans")
    fit = fit_latent_class([sum(calls[k]) for k in fused_keys], [length(calls[k]) for k in fused_keys], st)
    post = Dict(zip(fused_keys, fit.posterior))
    out = abspath(o["--outdir"]); mkpath(out)
    rows = Dict{String,String}[]; nfused = 0
    for key in sort(collect(keys(pred)))
        r = copy(pred[key]); pk = get(phys, key, "")
        r["scan_predicted"] = r["predicted"]; r["scan_probability_1"] = r["probability_1"]
        r["physical_lobe"] = pk
        r["fusion_scans"] = haskey(calls, pk) ? string(length(calls[pk])) : "0"
        r["fusion_calls_1"] = haskey(calls, pk) ? string(sum(calls[pk])) : "0"
        if haskey(post, pk)
            p = post[pk]; nfused += 1
            r["predicted"] = p >= 0.5 ? "1" : "0"
            r["probability_1"] = @sprintf("%.8f", p)
            r["confidence"] = @sprintf("%.8f", abs(2p - 1))
            r["scan_predicted"] == "?" && (r["invalid_reason"] = "fused_from_repeated_scans")
            r["model"] = r["model"] * "+latent_class_fusion"
        end
        push!(rows, r)
    end
    write_table(joinpath(out, "predictions_fused.tsv"),
        ["file", "lobe", "predicted", "confidence", "probability_1", "invalid_reason", "model",
         "scan_predicted", "scan_probability_1", "physical_lobe", "fusion_scans", "fusion_calls_1"], rows)
    write_table(joinpath(out, "physical_lobes.tsv"), ["physical_lobe", "scans", "calls_1", "posterior_1"],
        [Dict("physical_lobe" => k, "scans" => string(length(calls[k])), "calls_1" => string(sum(calls[k])),
              "posterior_1" => haskey(post, k) ? @sprintf("%.8f", post[k]) : "NA") for k in sort(collect(keys(calls)))])
    write_table(joinpath(out, "fusion_params.tsv"), ["pi", "theta0", "theta1", "iterations", "physical_lobes_fused", "rows_fused"],
        [Dict("pi" => @sprintf("%.8f", fit.pi), "theta0" => @sprintf("%.8f", fit.theta0), "theta1" => @sprintf("%.8f", fit.theta1),
              "iterations" => string(fit.iterations), "physical_lobes_fused" => string(length(fused_keys)), "rows_fused" => string(nfused))])
    @printf("Molecule fusion: %d physical lobes, %d rows fused; pi=%.3f theta0=%.3f theta1=%.3f\n",
            length(fused_keys), nfused, fit.pi, fit.theta0, fit.theta1)
    return joinpath(out, "predictions_fused.tsv")
end

function main(args=ARGS)
    o = fusion_options(args)
    o === nothing ? println("build_molecule_fusion.jl --data-dir DIR --consensus DIR --features TSV --predictions TSV --config TOML --outdir NEW_DIR") :
        build_fusion(o)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
