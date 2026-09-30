#!/usr/bin/env julia
# Research harness (2026-09-30 end-lobe robustness campaign, see journal):
# rerun only the classifier stages (GMM, k-means, vote, fusion) of a finished
# run_molecule_consensus_chitosan.jl run on its cached features, with the
# feature views of a declared variant. No fitting; labels are never read.
#   julia --project=. test/assignment_view_variants.jl RUN_DIR VARIANT NEW_OUT_DIR RAW_DIR
# VARIANT: v0_current (reproduces the run byte for byte), v1_interior, v2_shape.
# OUT_DIR is run-like (symlinks to the run's counting/consensus files), so
# test/grade_consensus_run.jl verifies and grades it unchanged.
const REPO = dirname(@__DIR__)
include(joinpath(REPO, "test", "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using Statistics, Printf

const CONFIG = joinpath(REPO, "config", "unit_assignment_corroborated_training.toml")
const BASE4 = "amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
const INT4 = "amp_int,integ_int,amp_prom_int,integ_prom_int"

"Neighbour contrasts that never use a terminal lobe as reference (interior-referenced BASE4)."
function add_interior_features!(header, rows)
    byfile = Dict{String,Vector{Dict{String,String}}}()
    for r in rows; push!(get!(byfile, r["file"], Dict{String,String}[]), r); end
    for rs in values(byfile)
        sort!(rs; by=r -> parse(Int, r["lobe"]))
        n = length(rs)
        amp = [parse(Float64, r["amplitude"]) for r in rs]
        integ = [parse(Float64, r["amplitude"]) * parse(Float64, r["sigma_parallel_nm"]) * parse(Float64, r["sigma_perp_nm"]) for r in rs]
        interior = n >= 3 ? (2:n-1) : (1:n)
        ma, mi = median(amp[interior]), median(integ[interior])
        for i in 1:n
            # neighbours excluding terminal lobes (a terminal lobe's own neighbour is kept)
            nb = [j for j in (i-1, i+1) if 1 <= j <= n && (j in interior || !(i in interior))]
            isempty(nb) && (nb = [j for j in (i-1, i+1) if 1 <= j <= n])
            rs[i]["amp_int"] = @sprintf("%.10g", log(amp[i] / ma))
            rs[i]["integ_int"] = @sprintf("%.10g", log(integ[i] / mi))
            rs[i]["amp_prom_int"] = @sprintf("%.10g", log(amp[i] / mean(amp[nb])))
            rs[i]["integ_prom_int"] = @sprintf("%.10g", log(integ[i] / mean(integ[nb])))
        end
    end
    append!(header, ["amp_int", "integ_int", "amp_prom_int", "integ_prom_int"])
end

function variant_views(variant)
    variant == "v0_current" && return (gmm="v_cc=$BASE4,patch_u_asym_reconstructed,mold_cc_fwd,mold_cc_bwd,emp_fisher",
        km=["v_base=$BASE4", "v_split=$BASE4,split_log_skew", "v_comt=$BASE4,bwd_neg_com_t", "v_diag45=$BASE4,bwd_neg_diag45"])
    variant == "v1_interior" && return (gmm="v_cc=$INT4,patch_u_asym_reconstructed,mold_cc_fwd,mold_cc_bwd,emp_fisher",
        km=["v_base=$INT4", "v_split=$INT4,split_log_skew", "v_comt=$INT4,bwd_neg_com_t", "v_diag45=$INT4,bwd_neg_diag45"])
    variant == "v2_shape" && return (gmm="v_cc=patch_u_asym_reconstructed,mold_cc_fwd,mold_cc_bwd,emp_fisher",
        km=["v_split=split_log_skew,patch_u_asym_reconstructed", "v_comt=bwd_neg_com_t,patch_u_asym_reconstructed",
            "v_diag45=bwd_neg_diag45,patch_u_asym_reconstructed", "v_mold=mold_cc_fwd,mold_cc_bwd,emp_fisher"])
    error("unknown variant $variant")
end

function run_variant(run, variant, out)
    (ispath(out) || islink(out)) && error("Output exists: $out")
    mkpath(joinpath(out, "assignment"))
    for f in ("counting_summary.tsv", "scan_geometry.tsv", "consensus", "raw_hashes.tsv", "input_hashes.tsv")
        symlink(joinpath(run, f), joinpath(out, f))
    end
    src = joinpath(run, "assignment")
    symlink(joinpath(src, "features.tsv"), joinpath(out, "assignment", "features.tsv"))
    header, rows = read_table(joinpath(src, "features_predictor.tsv"))
    variant == "v1_interior" && add_interior_features!(header, rows)
    table = joinpath(out, "assignment", "features_predictor.tsv")
    write_table(table, header, rows)
    v = variant_views(variant)
    consensus = joinpath(run, "consensus", "consensus_summary_final.tsv")
    julia = Base.julia_cmd()
    common = ["--features", table, "--first-seed", "0", "--interactions", "--training-scans", consensus]
    gmm = joinpath(out, "assignment", "pred_gmm.tsv"); km = joinpath(out, "assignment", "pred_kmeans.tsv")
    Base.run(pipeline(`$julia --project=$REPO $(joinpath(REPO, "test", "build_labelfree_gmm_predictions.jl")) $common --config $CONFIG --out $gmm --view $(v.gmm) --seeds 10 --selftrain 2`,
                      stdout=joinpath(out, "gmm.log"), stderr=joinpath(out, "gmm.log")))
    kmviews = vcat([["--view", s] for s in v.km]...)
    Base.run(pipeline(`$julia --project=$REPO $(joinpath(REPO, "test", "build_labelfree_unit_predictions.jl")) $common --out $km --patches $(joinpath(src, "patches_bwd17.tsv")) $kmviews --seeds 20`,
                      stdout=joinpath(out, "kmeans.log"), stderr=joinpath(out, "kmeans.log")))
    pred = joinpath(out, "assignment", "predictions.tsv")
    write_soft_vote(table, km, gmm, pred, CONFIG)
end

function fuse(run, out, rawdir)
    julia = Base.julia_cmd()
    Base.run(pipeline(`$julia --project=$REPO $(joinpath(REPO, "test", "build_molecule_fusion.jl")) --data-dir $rawdir --consensus $(joinpath(run, "consensus")) --features $(joinpath(out, "assignment", "features.tsv")) --predictions $(joinpath(out, "assignment", "predictions.tsv")) --config $(joinpath(REPO, "config", "molecule_consensus.toml")) --outdir $(joinpath(out, "fusion"))`,
                      stdout=joinpath(out, "fusion.log"), stderr=joinpath(out, "fusion.log")))
    cp(joinpath(out, "fusion", "predictions_fused.tsv"), joinpath(out, "predictions.tsv"))
end

run, variant, out, rawdir = abspath(ARGS[1]), ARGS[2], abspath(ARGS[3]), abspath(ARGS[4])
run_variant(run, variant, out)
fuse(run, out, rawdir)
println("done ", variant, " -> ", out)
