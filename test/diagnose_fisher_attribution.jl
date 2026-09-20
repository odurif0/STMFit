#!/usr/bin/env julia
"""
Explicit fixed native Fisher replay with newly fitted fold coefficients.
This is NOT the no-refit saved-patch representation diagnostic. Real cohorts
run only in the approved Viper job; local tests use synthetic patches only.
"""
module FisherAttributionDiagnostic
using LinearAlgebra, Statistics, Printf, SHA
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
const RD = ReconstructedRepresentationDiagnostics
const EF = RD.EmpiricalFisherNative
const RU = RD.ReconstructedUnitAssignment
export replay_folds, compare_saved, main

const USAGE = """
Usage: julia --startup-file=no --threads=1 --project=. test/diagnose_fisher_attribution.jl \\
  --patches PATH --prefix res --config PATH --saved-scores PATH --outdir NEW_DIR

Explicit fixed native PCA/GMM/Fisher opposite-parity fold REPLAY, not recovery of
saved weights. All options are required. Original producer uses FORWARD 17x17
patches, prefix res. Only frozen config settings are used; no tuning or search.
Compare every key, six-decimal score text and native reason with saved fisher_cv
BEFORE attributing any saved pipeline behavior. Differences are retained, never
repaired. Exact exported-score agreement does not establish unique/original weights.
All invalid keys remain. Output must not exist, including empty dirs or symlinks.
Real cohorts: approved Viper job only. Local use: synthetic tests only.
"""

"Native cv_scores loop with models retained, no changed fitting/scoring arithmetic."
function replay_folds(patches::EF.PatchTable, options::EF.FisherOptions)
    options.patch_projection == "none" ||
        throw(ArgumentError("Raw-patch attribution requires fisher_patch_projection=none"))
    n = length(patches.keys)
    size(patches.X,1) == n && length(patches.amplitudes) == n && length(patches.invalid_reasons) == n ||
        throw(DimensionMismatch("patch table rows must match keys"))
    grid = EF.fisher_grid(options)
    grid.disk_indices == patches.grid.disk_indices && grid.coords == patches.grid.coords ||
        throw(ArgumentError("patch table grid does not match Fisher config"))
    reasons = copy(patches.invalid_reasons)
    for i in 1:n
        if isempty(reasons[i]) && !(all(isfinite,patches.X[i,:]) && isfinite(patches.amplitudes[i]))
            reasons[i] = "nonfinite_patch_input"
        end
    end
    margins = fill(NaN,n)
    valid = isempty.(reasons)
    models = Dict{Int,EF.FisherModel}()
    folds = Dict{String,Any}[]
    for parity in (0,1)
        train = [i for i in 1:n if valid[i] && mod(patches.keys[i][2],2) == parity]
        held = [i for i in 1:n if valid[i] && mod(patches.keys[i][2],2) != parity]
        fold = parity == 0 ? "even" : "odd"
        info = Dict{String,Any}("training_parity"=>fold,"training_rows"=>length(train),
            "heldout_rows"=>length(held),"status"=>"not_fit_no_heldout_rows", "fit_reason"=>"",
            "gmm_converged"=>"NA","gmm_iterations"=>"NA","gmm_lower_bound"=>NaN,
            "shared_centering_offset"=>NaN,"weight_l2_norm"=>NaN,"mid_l2_norm"=>NaN)
        push!(folds,info)
        isempty(held) && continue
        model = try
            EF.fit_fisher(patches.X[train,:],patches.amplitudes[train],options)
        catch error
            error isa EF.FisherFitError || rethrow()
            info["status"], info["fit_reason"] = "invalid_fit", error.reason
            for i in held
                reasons[i] = "training_$fold:" * error.reason
            end
            continue
        end
        models[parity] = model
        info["status"] = "fitted_new_realization"
        info["gmm_converged"], info["gmm_iterations"] = model.gmm.converged, model.gmm.iterations
        info["gmm_lower_bound"] = model.gmm.lower_bound
        info["shared_centering_offset"] = -dot(model.mid,model.w_p)
        info["weight_l2_norm"], info["mid_l2_norm"] = norm(model.w_p), norm(model.mid)
        model.gmm.converged || @warn "Fisher GMM reached configured iteration limit" fold iterations=model.gmm.iterations
        for i in held
            value = EF.maxmirror_score(patches.X[i,:],model,grid)
            if isfinite(value)
                margins[i] = value
            else
                reasons[i] = "nonfinite_fisher_score"
            end
        end
    end
    scores = [EF.FisherScore(key[1],key[2],margins[i],reasons[i]) for (i,key) in enumerate(patches.keys)]
    return (; scores, models, folds)
end

"Compare complete exported score/reason/key text, not rounded residuals or correlations."
function compare_saved(saved_path, replay_path)
    _, saved = RU.lobe_table(saved_path; required=["score","invalid_reason"], contiguous=false)
    _, replay = RU.lobe_table(replay_path; required=["score","invalid_reason"], contiguous=false)
    keys = sort(collect(union(Set(Base.keys(saved)),Set(Base.keys(replay)))))
    rows = Dict{String,Any}[]
    for key in keys
        s, r = get(saved,key,nothing), get(replay,key,nothing)
        flags = String[]
        s === nothing && push!(flags,"missing_saved_key")
        r === nothing && push!(flags,"missing_replay_key")
        if s !== nothing && r !== nothing
            s["score"] == r["score"] || push!(flags,"score_text_difference")
            s["invalid_reason"] == r["invalid_reason"] || push!(flags,"reason_difference")
        end
        sv, rv = s === nothing ? "NA" : s["score"], r === nothing ? "NA" : r["score"]
        push!(rows, Dict("file"=>key[1],"lobe"=>key[2],"saved_score"=>sv,"replayed_score"=>rv,
            "saved_reason"=>s===nothing ? "missing_key" : s["invalid_reason"],
            "replayed_reason"=>r===nothing ? "missing_key" : r["invalid_reason"],
            "replayed_minus_saved_serialized"=>RD.number(rv)-RD.number(sv),
            "status"=>isempty(flags) ? "match" : join(flags,",")))
    end
    all_match = all(r["status"] == "match" for r in rows)
    return (; rows, all_match, saved_keys=length(saved), replay_keys=length(replay),
            matching=count(r->r["status"] == "match",rows))
end

function attribution_rows(patches, replay, comparison)
    status = comparison.all_match ? "new_realization_matches_all_saved_exports_not_original_weights" :
                                    "new_realization_only_saved_exports_differ"
    rows = Dict{String,Any}[]
    for (i,key) in enumerate(patches.keys)
        parity = 1-mod(key[2],2)
        score = replay.scores[i]
        row = Dict{String,Any}("file"=>key[1],"lobe"=>key[2],
            "training_parity"=>parity==0 ? "even" : "odd", "attribution_scope"=>status,
            "native_invalid_reason"=>score.invalid_reason,
            "original_response"=>NaN,"reflected_response"=>NaN,"even_uncentered"=>NaN,
            "shared_centering_offset"=>NaN,"even_centered"=>NaN,"odd_response"=>NaN,
            "absolute_odd_bonus"=>NaN,"maximum_response"=>NaN,"analytic_identity_error"=>NaN)
        if isempty(score.invalid_reason) && isfinite(score.score) && haskey(replay.models,parity)
            model = replay.models[parity]
            parts = RD.response_decomposition(patches.X[i,:],model.w_p,model.mid,patches.grid)
            row["original_response"], row["reflected_response"] = parts.original, parts.mirrored
            row["even_uncentered"], row["shared_centering_offset"] = parts.even_uncentered, parts.shared_centering_offset
            row["even_centered"], row["odd_response"] = parts.even_centered, parts.odd
            row["absolute_odd_bonus"], row["maximum_response"] = parts.reflection_bonus, parts.direct_max
            row["analytic_identity_error"] = parts.analytic_max-parts.direct_max
        end
        push!(rows,row)
    end
    return rows
end

function run_replay(patch_path,prefix,config,saved_path,outdir)
    VERSION.major == 1 && VERSION.minor == 13 || error("Requires Julia 1.13")
    Threads.nthreads() == 1 || error("Use --threads=1")
    prefix == "res" || error("Original Fisher uses forward patches with prefix res, not backward patches")
    (ispath(outdir) || islink(outdir)) && error("Output already exists: $outdir")
    paths = Dict("patches"=>patch_path,"config"=>config,"saved_scores"=>saved_path)
    hashes = Dict(k=>bytes2hex(sha256(read(path))) for (k,path) in paths)
    options = EF.load_fisher_config(config)
    # Fail before fitting if the supplied comparison is not a native score table.
    header,_ = RU.lobe_table(saved_path; required=["score","invalid_reason"],contiguous=false)
    any(c->lowercase(c) in ("predicted","label","probability_1","confidence"),header) &&
        error("Saved comparison must contain native Fisher scores, not chemical predictions")
    patches = EF.load_patches(patch_path,prefix,options)
    isempty(patches.keys) && error("No patch rows")
    replay = replay_folds(patches,options)
    mkpath(dirname(abspath(outdir)))
    mkdir(outdir)
    EF.write_scores(joinpath(outdir,"replayed_fisher_cv.tsv"),replay.scores)
    # Comparison occurs before constructing any attribution of these new weights.
    comparison = compare_saved(saved_path,joinpath(outdir,"replayed_fisher_cv.tsv"))
    RD.write_rows(joinpath(outdir,"saved_score_comparison.tsv"),
        ["file","lobe","saved_score","replayed_score","saved_reason","replayed_reason",
         "replayed_minus_saved_serialized","status"],comparison.rows)
    rows = attribution_rows(patches,replay,comparison)
    RD.write_rows(joinpath(outdir,"replayed_response_attribution.tsv"),
        ["file","lobe","training_parity","attribution_scope","native_invalid_reason",
         "original_response","reflected_response","even_uncentered","shared_centering_offset",
         "even_centered","odd_response","absolute_odd_bonus","maximum_response",
         "analytic_identity_error"],rows)
    weights = Dict{String,Any}[]
    for parity in sort(collect(keys(replay.models)))
        model = replay.models[parity]
        for (disk_index, full_index) in enumerate(patches.grid.disk_indices)
            u_index,t_index = divrem(full_index-1,patches.grid.side)
            push!(weights,Dict("training_parity"=>parity==0 ? "even" : "odd", "disk_index"=>disk_index,
                "full_patch_index"=>full_index,"t_nm"=>patches.grid.coords[t_index+1],
                "u_nm"=>patches.grid.coords[u_index+1],"new_weight"=>model.w_p[disk_index],
                "new_centered_mid"=>model.mid[disk_index],"new_centered_low_mean"=>model.g0[disk_index],
                "new_centered_high_mean"=>model.g1[disk_index]))
        end
    end
    RD.write_rows(joinpath(outdir,"new_fold_weights_mid.tsv"),["training_parity","disk_index","full_patch_index",
        "t_nm","u_nm","new_weight","new_centered_mid","new_centered_low_mean","new_centered_high_mean"],weights)
    RD.write_rows(joinpath(outdir,"new_fold_fit_status.tsv"),["training_parity","training_rows","heldout_rows",
        "status","fit_reason","gmm_converged","gmm_iterations","gmm_lower_bound",
        "shared_centering_offset","weight_l2_norm","mid_l2_norm"],replay.folds)
    all(bytes2hex(sha256(read(paths[k])))==h for (k,h) in hashes) || error("Input changed during replay")
    RD.write_rows(joinpath(outdir,"inputs.tsv"),["input","path","sha256"],
        [Dict("input"=>k,"path"=>abspath(paths[k]),"sha256"=>hashes[k]) for k in sort(collect(keys(paths)))])
    open(joinpath(outdir,"report.md"),"w") do io
        println(io,"# Explicit fixed native Fisher replay\n")
        println(io,"Julia $(VERSION), one thread. **Newly fitted realization**, not recovered historical weights. This is distinct from the earlier no-refit saved-patch audit. Frozen config, original forward patches/prefix res, opposite-parity folds, native initialization and scoring were used unchanged. No downstream classifier, grade, alternative settings, seed search or correction ran.\n")
        println(io,"## Complete saved-export comparison (before attribution)\n")
        println(io,"Saved keys: $(comparison.saved_keys); replay keys: $(comparison.replay_keys); exact key/score-text/reason matches: $(comparison.matching)/$(length(comparison.rows)). Full six-decimal export agreement: **$(comparison.all_match)**. All differences and invalids remain in saved_score_comparison.tsv.\n")
        if comparison.all_match
            println(io,"The newly fitted folds reproduce every saved six-decimal score and reason. This permits a response-attribution diagnostic for a realization consistent with those exports. It does NOT prove identical or unique original coefficients; rounding hides smaller score differences.\n")
        else
            println(io,"**Do not attribute these new coefficients or response parts to the saved pipeline realization.** They describe only this fixed replay. Differences are retained without tuning or retries to obtain a match.\n")
        end
        println(io,"## Exact response identity\n")
        println(io,"For R the unchanged legacy disk reflection, e=w⋅(x+R*x)/2, o=w⋅(x-R*x)/2, c=-w⋅mid. Native max-reflection response = e+c+abs(o). The shared offset is not mirrored; mid is the native mean of centered training patches, not the original training mean. Legacy reflection reverses physical t, not u. No convention was corrected.\n")
        println(io,"New weights and midpoints are in new_fold_weights_mid.tsv; fit termination and training/held-out support are in new_fold_fit_status.tsv. Invalid rows have NA attribution and native reasons. All numeric attribution values are explicitly NEW replay responses, not new chemical predictions.\n")
        println(io,"| Quantity | Finite | Median | Minimum | Maximum |\n|---|---:|---:|---:|---:|")
        for metric in ("original_response","reflected_response","even_centered","shared_centering_offset",
                       "absolute_odd_bonus","maximum_response","analytic_identity_error")
            s = RD.finite_summary([row[metric] for row in rows])
            println(io,"| $metric | $(s.n) | $(RD.text(s.median)) | $(RD.text(s.min)) | $(RD.text(s.max)) |")
        end
        println(io,"\nThe nonnegative abs(odd) term can increase the maximum relative to the centered even response. Neither its size nor response positivity proves chemical meaning, a preprocessing fault, or unique acquisition causation. Scores in opposite held-out parities use different fitted weights; they are not independently calibrated chemical probabilities.")
    end
    println("Fixed Fisher replay: $(length(patches.keys)) keys; saved six-decimal/reason matches $(comparison.matching)/$(length(comparison.rows)); all_match=$(comparison.all_match)")
    println("New realization only; output: ",outdir)
    return comparison
end

function main(args=ARGS)
    args in (["--help"],["-h"]) && (print(USAGE);return 0)
    required=["--patches","--prefix","--config","--saved-scores","--outdir"]
    iseven(length(args)) || error("Expected option/value pairs; use --help")
    opts=Dict{String,String}()
    for i in 1:2:length(args)
        k,v=args[i],args[i+1]
        k in required || error("Unknown option: $k")
        haskey(opts,k) && error("Duplicate option: $k")
        isempty(v) && error("Empty value: $k")
        opts[k]=v
    end
    all(haskey(opts,k) for k in required) || error("Missing required option; use --help")
    run_replay(opts["--patches"],opts["--prefix"],opts["--config"],opts["--saved-scores"],opts["--outdir"])
    return 0
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(FisherAttributionDiagnostic.main())
    catch err
        err isa InterruptException && rethrow()
        print(stderr,"Fixed Fisher attribution replay: ");showerror(stderr,err);println(stderr)
        exit(1)
    end
end
