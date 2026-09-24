#!/usr/bin/env julia
# Real-data learning runs on Viper only. No refit, N selection or benchmark labels.
module ScanExclusionDiagnostic
using LinearAlgebra, TOML, Printf
include(joinpath(@__DIR__, "diagnose_frozen_learning.jl"))
const D = FrozenLearningDiagnostic
const RU = D.RU
const STAGES = ("pred_gmm", "pred_kmeans", "predictions")

function settings(path)
    s=TOML.parsefile(path)
    s==Dict("model"=>Dict("method"=>"native_saved_input_leave_one_scan_out", "repetitions"=>2,
            "assignment_config"=>"unit_assignment_patch_support.toml"),
        "selection"=>Dict("input_policy"=>"first_observed_fit_execution_all_keys",
            "arms"=>["control","observed"], "exclude_from"=>["fisher","gmm","kmeans","cluster_naming"],
            "fisher_inner_folds"=>"unchanged_lobe_parity", "require_exact_native_reapplication"=>true,
            "grade_policy"=>"external_after_complete_verification"),
        "preprocessing"=>Dict("local_inputs"=>"saved_unmodified",
            "normalization"=>"each_scan_own_unlabelled_features", "missing_features"=>"native_abstention_no_imputation")) ||
        error("Only the declared scan-exclusion experiment is supported")
    s
end

function parse_cli(args)
    args==["--help"] && return nothing
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in ("--input","--settings","--outdir","--arm") && i<length(args) && !startswith(args[i+1],"--") ||
            error("Unknown, missing or forbidden option: $k")
        o[k]=k=="--arm" ? args[i+1] : abspath(args[i+1]); i+=2
    end
    all(haskey(o,k) for k in ("--input","--settings","--outdir","--arm")) || error("All options required")
    o["--arm"] in ("control","observed") || error("Unknown local input arm")
    isdir(o["--input"]) || error("Input directory missing")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    settings(o["--settings"]); o
end

function training_patches(p, excluded)
    idx=findall(k->k[1]!=excluded,p.keys)
    isempty(idx) && error("Empty training cohort")
    isempty(excluded) || any(k->k[1]==excluded,p.keys) || error("Excluded file not present")
    D.EF.PatchTable(p.keys[idx],p.X[idx,:],p.amplitudes[idx],p.invalid_reasons[idx],p.grid)
end

function split_records(rec, excluded)
    train=(g=filter(r->r.file!=excluded,rec.g),k=filter(r->r.file!=excluded,rec.k))
    target=isempty(excluded) ? rec : (g=filter(r->r.file==excluded,rec.g),k=filter(r->r.file==excluded,rec.k))
    (isempty(train.g) || isempty(target.g)) && error("Empty training or target records")
    isempty(excluded) || isempty(intersect(D.key.(train.g),D.key.(target.g))) || error("Training/target overlap")
    train,target
end

"Target-only preprocessing; no fitted mixture, centers or class names."
function local_scaling(records,opt)
    g=[]
    D.Gaussian._standardized_matrix(records.g,D.GFEATURES;interactions=opt.g.interactions,
        normalization=opt.g.normalization,scale_fallback=opt.g.scale_fallback,normalization_state=g)
    k=Dict()
    for (_,features) in D.KVIEWS
        states=[]
        D.KMeans._standardized_matrix(records.k,features;interactions=opt.k.interactions,normalization_state=states)
        k[join(features,',')]=(;scales=states)
    end
    (;g=(;scales=g),k)
end

function scaling_rows(bank)
    rows=Dict{String,String}[]
    for (view,h) in vcat(["gmm"=>bank.g],[name=>bank.k[join(f,',')] for (name,f) in D.KVIEWS]), s in h.scales
        push!(rows,Dict("view"=>view,"file"=>s.file,"feature"=>s.feature,
            "center"=>string(s.center),"scale"=>string(s.scale)))
    end
    rows
end

function write_target(dir,records,scales,bank,opt,cfgpath)
    p=D.predict_heads(records,scales,bank,opt)
    D.Gaussian._write_predictions(joinpath(dir,"pred_gmm.tsv"),records.g,p.g,Int.(isfinite.(p.g)))
    D.KMeans._write_predictions(joinpath(dir,"pred_kmeans.tsv"),records.k,p.k,p.kviews)
    RU.write_soft_vote(joinpath(dir,"target_features.tsv"),joinpath(dir,"pred_kmeans.tsv"),
        joinpath(dir,"pred_gmm.tsv"),joinpath(dir,"predictions.tsv"),cfgpath)
end

function fit_fold(root,arm,dir,excluded,patches,opts,data)
    mkdir(dir)
    # Subset BEFORE PCA/GMM/Fisher. Indices in saved Fisher states refer to this
    # sorted training-only PatchTable, not the full cohort. No mask-policy change.
    training=training_patches(patches,excluded)
    fisher=D.capture_fisher(training,opts)
    all_scores=D.frozen_fisher(patches,fisher.models)
    lookup=Dict((r.file,r.lobe)=>r for r in all_scores)
    all(isequal(lookup[(r.file,r.lobe)].score,r.score) &&
        lookup[(r.file,r.lobe)].invalid_reason==r.invalid_reason for r in fisher.native) || error("Fisher subset reapplication differs")
    D.EF.write_scores(joinpath(dir,"fisher_cv.tsv"),all_scores)
    ft=joinpath(dir,"features_predictor.tsv"); D.joined(root,arm,joinpath(dir,"fisher_cv.tsv"),ft,data.cfg)
    pp=joinpath(root,arm,"patches_bwd17.tsv"); rec=D.head_records(ft,pp)
    train,target=split_records(rec,excluded)
    opt=D.options(ft,pp,joinpath(dir,"unused.tsv"),data.cfgpath,data.cfg)
    g=D.capture_head("gmm",train.g,D.GFEATURES,opt.g)
    k=Dict(join(features,',')=>D.capture_head("kmeans",train.k,features,opt.k) for (_,features) in D.KVIEWS)
    bank=(;g,k); D.save_bank(joinpath(dir,"models.toml"),bank,fisher)
    saved=D.load_bank(joinpath(dir,"models.toml"))
    isequal(saved.fisher,fisher.models) || error("Fisher serialization differs")
    for h in (bank.g,values(bank.k)...)
        restored=h.kind=="gmm" ? saved.bank.g : saved.bank.k[join(h.features,',')]
        isequal(h.models,restored.models) && isequal(h.scales,restored.scales) || error("Head serialization differs")
    end
    scales=local_scaling(target,opt)
    RU.write_table(joinpath(dir,"target_scales.tsv"),["view","file","feature","center","scale"],scaling_rows(scales))
    header,rows=RU.lobe_table(ft)
    RU.write_table(joinpath(dir,"target_features.tsv"),header,[rows[D.key(r)] for r in target.g])
    write_target(dir,target,scales,saved.bank,opt,data.cfgpath)
    targetcounts=isempty(excluded) ? data.counts : Dict(excluded=>data.counts[excluded])
    for stage in STAGES; D.Pipeline.check_counts(joinpath(dir,stage*".tsv"),targetcounts); end
    if isempty(excluded)
        for stage in D.REFERENCES
            read(joinpath(dir,stage*".tsv"))==read(joinpath(root,arm,stage*".tsv")) || error("Native reference differs: $stage")
        end
    end
    nothing
end

function merge_predictions(out,folds,counts)
    for stage in STAGES
        header=nothing; rows=Dict{String,String}[]
        for fold in folds
            h,r=RU.read_table(joinpath(out,fold["fold"],stage*".tsv"))
            header===nothing ? (header=h) : (header==h || error("Fold schema mismatch"))
            all(row["file"]==fold["file"] for row in r) || error("Unexpected held-out file")
            append!(rows,r)
        end
        RU.write_table(joinpath(out,stage*".tsv"),header,rows)
        D.Pipeline.check_counts(joinpath(out,stage*".tsv"),counts)
    end
end

function execute(root,settings_path,out,arm)
    v"1.13"<=VERSION<v"1.14" || error("Julia 1.13 required")
    Threads.nthreads()==1 || error("One Julia thread per arm/repetition required")
    BLAS.set_num_threads(1)
    settings(settings_path); arm in ("control","observed") || error("Unknown arm")
    data=D.inputs(root); length(data.counts)>=3 || error("At least three scans required")
    (ispath(out)||islink(out)) && error("Output exists")
    mkpath(out); cp(settings_path,joinpath(out,"settings.toml")); cp(data.cfgpath,joinpath(out,"assignment.toml"))
    cp(joinpath(root,"selected_summary.tsv"),joinpath(out,"selected_summary.tsv"))
    RU.write_table(joinpath(out,"inputs.tsv"),["path","sha256"],
        [Dict("path"=>p,"sha256"=>h) for (p,h) in sort(collect(data.hashes))])
    RU.write_table(joinpath(out,"groups.tsv"),["file","group"],
        [Dict("file"=>f,"group"=>data.groups[f]) for f in sort(collect(keys(data.groups)))])
    RU.write_table(joinpath(out,"arm.tsv"),["arm"],[Dict("arm"=>arm)])
    opts=D.EF.load_fisher_config(data.cfg)
    patches=D.EF.load_patches(joinpath(root,arm,"patches_fwd17.tsv"),"res",opts)
    println("Native exact reference: $arm"); flush(stdout)
    fit_fold(root,arm,joinpath(out,"native"),"",patches,opts,data)
    files=sort(collect(keys(data.counts)))
    folds=[Dict("fold"=>@sprintf("fold%04d",i),"file"=>f,"training_scans"=>string(length(files)-1),
        "training_lobes"=>string(sum(values(data.counts))-data.counts[f]),"target_lobes"=>string(data.counts[f])) for (i,f) in enumerate(files)]
    RU.write_table(joinpath(out,"folds.tsv"),["fold","file","training_scans","training_lobes","target_lobes"],folds)
    for (i,fold) in enumerate(folds)
        println("$arm $i/$(length(folds)): excluding $(fold["file"])"); flush(stdout)
        fit_fold(root,arm,joinpath(out,fold["fold"]),fold["file"],patches,opts,data)
    end
    merge_predictions(out,folds,data.counts)
    all(D.sha(p)==h for (p,h) in data.hashes) || error("Inputs changed")
    println("Complete: $arm, $(length(files)) excluded-scan folds, $(sum(values(data.counts))) target lobes. Grading remains external.")
end

function main(args=ARGS)
    o=parse_cli(args)
    o===nothing && return println("diagnose_scan_exclusion.jl --input OBSERVED_RUN --settings TOML --outdir NEW_DIR --arm control|observed [--dry-run]")
    if haskey(o,"--dry-run")
        d=D.inputs(o["--input"])
        println("Dry run: $(length(d.counts)) scans / $(sum(values(d.counts))) lobes, arm $(o["--arm"]); no learning or output")
        return
    end
    execute(o["--input"],o["--settings"],o["--outdir"],o["--arm"])
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ScanExclusionDiagnostic.main()
