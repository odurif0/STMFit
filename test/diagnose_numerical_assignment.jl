#!/usr/bin/env julia
# Fixed diagnostics for the two independent numerical arms. No benchmark input.
module NumericalAssignmentDiagnostics
using LinearAlgebra, Statistics
include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .EmpiricalFisherNative, .ReconstructedUnitAssignment
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end

const ROOT = dirname(@__DIR__)
const FEATURES = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
    "patch_u_asym_reconstructed", "mold_cc_fwd", "mold_cc_bwd", "emp_fisher"]
config(name) = joinpath(ROOT, "config", "unit_assignment_" * name * ".toml")
key(r) = (r.file, r.lobe)
number(s) = something(tryparse(Float64, s), NaN)

"Two complementary whole-scan halves in sorted filename order, never lobe halves."
function scan_partitions(records)
    files = sort!(unique(r.file for r in records))
    length(files) >= 2 || error("At least two scans required")
    membership = Dict(f => mod1(i, 2) for (i, f) in enumerate(files))
    return [findall(r -> membership[r.file] == half, records) for half in 1:2]
end

"Describe final covariances and verify identical pre-shrink groups, means and weights."
function covariance_pairs(base, candidate, partition)
    length(base) == length(candidate) || error("Seed coverage differs")
    rows = Dict{String,Any}[]
    for (b, c) in zip(base, candidate)
        b.seed == c.seed && b.indices == c.indices || error("Seed/input mismatch")
        length(b.clusters) == length(c.clusters) || error("Cluster coverage differs")
        for (bc, cc) in zip(b.clusters, c.clusters)
            (bc.component == cc.component && bc.members == cc.members &&
             bc.mean == cc.mean && bc.weight == cc.weight && bc.sample == cc.sample) ||
                error("Shrinkage changed a pre-final-score fit")
            for (mode, cluster) in (("ridge", bc), ("ledoit_wolf", cc))
                eigenvalues = eigvals(Symmetric(cluster.covariance))
                push!(rows, Dict("partition"=>partition, "seed"=>b.seed,
                    "component"=>cluster.component, "mode"=>mode,
                    "members"=>length(cluster.members), "weight"=>cluster.weight,
                    "shrinkage"=>cluster.shrinkage, "eigen_min"=>minimum(eigenvalues),
                    "eigen_max"=>maximum(eigenvalues),
                    "condition"=>maximum(eigenvalues)/minimum(eigenvalues)))
            end
        end
    end
    return rows
end

function fisher_offsets(patches, options)
    rows = Dict{String,Any}[]
    models = Dict{Int,FisherModel}()
    options.patch_projection == "none" && options.score_center == "legacy_centered_mean" ||
        error("This diagnostic requires the unprojected legacy control")
    for parity in (0, 1)
        train = findall(i -> isempty(patches.invalid_reasons[i]) &&
            mod(patches.keys[i][2], 2) == parity, eachindex(patches.keys))
        model = fit_fisher(patches.X[train,:], patches.amplitudes[train], options)
        offset = dot(vec(mean(patches.X[train,:]; dims=1)) - model.mid, model.w_p)
        models[parity] = model
        push!(rows, Dict("training_parity"=>parity, "training_rows"=>length(train),
            "score_offset"=>offset, "training_mean_norm"=>norm(mean(patches.X[train,:]; dims=1)),
            "legacy_mid_norm"=>norm(model.mid), "converged"=>model.gmm.converged))
    end
    return rows, models
end

function run_diagnostics(control, centered, shrunk, outdir)
    (ispath(outdir) || islink(outdir)) && error("Output already exists: $outdir")
    baseconfig, centerconfig, shrinkconfig = config.(("affine_residual", "centered_fisher", "shrunk_gmm"))
    # Strict readers reject label columns; the predictor view is fixed above.
    tables = [last(lobe_table(joinpath(dir, "features_predictor.tsv"); required=FEATURES))
        for dir in (control, centered, shrunk)]
    for t in tables[2:end]; require_same_keys(tables[1], t, "numerical arm"); end
    read(joinpath(control,"features_predictor.tsv")) == read(joinpath(shrunk,"features_predictor.tsv")) ||
        error("Shrinkage arm changed predictor inputs")
    options = load_fisher_config(baseconfig)
    patches = load_patches(joinpath(control,"patches_fwd17.tsv"), "res", options)
    foldrows, models = fisher_offsets(patches, options)
    fscores = [last(lobe_table(joinpath(dir,"fisher_cv.tsv"); required=["score","invalid_reason"]))
        for dir in (control, centered)]
    for t in fscores; require_same_keys(tables[1],t,"Fisher scores"); end
    Set(patches.keys) == Set(keys(tables[1])) || error("Patch coverage mismatch")
    for (i,k) in enumerate(patches.keys)
        b, c = fscores[1][k], fscores[2][k]
        b["invalid_reason"] == c["invalid_reason"] || error("Fisher validity changed")
        if isempty(b["invalid_reason"])
            parity = 1-mod(k[2],2)
            raw = maxmirror_score(patches.X[i,:],models[parity],patches.grid)
            offset = foldrows[parity+1]["score_offset"]
            # Exported margins have six decimals. This is a serialization check,
            # not a learned or scientific acceptance threshold.
            isapprox(number(b["score"]),raw;atol=6e-7,rtol=1e-10) || error("Control Fisher replay mismatch")
            isapprox(number(c["score"]),raw-offset;atol=6e-7,rtol=1e-10) || error("Centered Fisher offset mismatch")
        end
    end
    records = GMM._load_records(joinpath(control,"features_predictor.tsv"))
    crecords = GMM._load_records(joinpath(centered,"features_predictor.tsv"))
    key.(records) == key.(crecords) || error("Predictor order differs")
    base_scaling = load_gmm_normalization(load_config(baseconfig))
    center_scaling = load_gmm_normalization(load_config(centerconfig))
    X, valid = GMM._standardized_matrix(records,FEATURES;interactions=true,
        normalization=base_scaling.mode,scale_fallback=base_scaling.scale_fallback)
    C, cvalid = GMM._standardized_matrix(crecords,FEATURES;interactions=true,
        normalization=center_scaling.mode,scale_fallback=center_scaling.scale_fallback)
    valid == cvalid || error("Standardization validity changed")
    zrows = [Dict("file"=>r.file,"lobe"=>r.lobe,"valid"=>valid[i],
        "fisher_z_control"=>X[i,8],"fisher_z_centered"=>C[i,8],
        "expanded_feature_distance"=>norm(X[i,:]-C[i,:])) for (i,r) in enumerate(records)]
    opts = [GMM._parse_cli(["--features",joinpath(control,"features_predictor.tsv"),
        "--config",cfg,"--view","v_cc="*join(FEATURES,','),"--seeds",string(load_config(cfg)["selection"]["gmm_seeds"]),
        "--first-seed",string(load_config(cfg)["selection"]["first_seed"]),"--interactions",
        "--selftrain",string(load_config(cfg)["selection"]["gmm_selftrain"])]) for cfg in (baseconfig,shrinkconfig)]
    diagnostics = [[],[]]
    full = [GMM._view_probability(records,FEATURES,opts[m];diagnostics=diagnostics[m]) for m in 1:2]
    covrows = covariance_pairs(diagnostics[1],diagnostics[2],"full")
    for (m,dir) in enumerate((control,shrunk))
        saved = last(lobe_table(joinpath(dir,"pred_gmm.tsv");required=["probability_1"]))
        require_same_keys(tables[1],saved,"GMM saved scores")
        for (i,r) in enumerate(records)
            p=number(saved[key(r)]["probability_1"])
            (isnan(p) && isnan(full[m][i])) || isapprox(p,full[m][i];atol=1e-8,rtol=0) ||
                error("GMM replay mismatch")
        end
    end
    stability = Dict{String,Any}[]
    for (half,idxs) in enumerate(scan_partitions(records))
        diag = [[],[]]
        probs = [GMM._view_probability(records[idxs],FEATURES,opts[m];diagnostics=diag[m]) for m in 1:2]
        append!(covrows,covariance_pairs(diag[1],diag[2],"scan_half_$half"))
        for (j,i) in enumerate(idxs), (m,mode) in enumerate(("ridge","ledoit_wolf"))
            p, q = full[m][i], probs[m][j]
            finite = isfinite(p) && isfinite(q)
            push!(stability,Dict("file"=>records[i].file,"lobe"=>records[i].lobe,
                "mode"=>mode,"scan_half"=>half,"full_vote"=>p,"half_vote"=>q,
                "both_available"=>finite,"absolute_vote_change"=>abs(p-q),
                "decision_changed"=>finite ? string((p>=0.5)!=(q>=0.5)) : "NA"))
        end
    end
    write_table(joinpath(outdir,"fisher_offsets.tsv"),["training_parity","training_rows","score_offset",
        "training_mean_norm","legacy_mid_norm","converged"],foldrows)
    write_table(joinpath(outdir,"standardized_features.tsv"),["file","lobe","valid","fisher_z_control",
        "fisher_z_centered","expanded_feature_distance"],zrows)
    write_table(joinpath(outdir,"covariances.tsv"),["partition","seed","component","mode","members",
        "weight","shrinkage","eigen_min","eigen_max","condition"],covrows)
    write_table(joinpath(outdir,"scan_withdrawal.tsv"),["file","lobe","mode","scan_half","full_vote",
        "half_vote","both_available","absolute_vote_change","decision_changed"],stability)
    println("Numerical diagnostics: $(length(records)) keys, $(length(covrows)) covariance rows.")
    println("Whole-scan withdrawal is a training perturbation, NOT held-out recognition or calibrated confidence.")
    return nothing
end

function main(args=ARGS)
    if args in (["--help"],["-h"])
        println("Usage: diagnose_numerical_assignment.jl --control DIR --centered DIR --shrunk DIR --outdir NEW")
        println("Fixed numerical diagnostics and two sorted whole-scan halves; no labels, grade or parameter search.")
        return 0
    end
    allowed=("--control","--centered","--shrunk","--outdir")
    length(args)==8 || error("Exactly four option/value pairs required; use --help")
    opts=Dict{String,String}()
    for i in 1:2:8
        args[i] in allowed && !haskey(opts,args[i]) || error("Unknown or duplicate option: $(args[i])")
        opts[args[i]]=args[i+1]
    end
    run_diagnostics((opts[k] for k in allowed)...)
    return 0
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(NumericalAssignmentDiagnostics.main())
