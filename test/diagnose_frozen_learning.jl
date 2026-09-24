#!/usr/bin/env julia
# Saved local inputs -> newly captured native learners -> fixed counterfactuals.
# Real cohorts run on Viper. No labels, image fit, threshold search or promotion.
module FrozenLearningDiagnostic
using LinearAlgebra, Statistics, TOML, SHA, Printf
module Pipeline
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
end
module Gaussian
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
module KMeans
include(joinpath(@__DIR__, "build_labelfree_unit_predictions.jl"))
end
include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
const EF = EmpiricalFisherNative
const RU = Pipeline.ReconstructedUnitAssignment
const ROOT = dirname(@__DIR__)
const BASE = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel"]
const GFEATURES = vcat(BASE, ["patch_u_asym_reconstructed", "mold_cc_fwd", "mold_cc_bwd", "emp_fisher"])
const KVIEWS = ["base"=>BASE, "split"=>vcat(BASE,["split_log_skew"]),
    "comt"=>vcat(BASE,["bwd_neg_com_t"]), "diag45"=>vcat(BASE,["bwd_neg_diag45"])]
# Ordered path: local observations -> Fisher -> scaling -> classifier parameters/names.
# The reverse frozen-input check exposes dependence on the choice of source bank.
const ARMS = [
    (name="control", input="control", fisher="control", scales="control", heads="control"),
    (name="local_only", input="observed", fisher="control", scales="control", heads="control"),
    (name="fisher", input="observed", fisher="observed", scales="control", heads="control"),
    (name="normalization", input="observed", fisher="observed", scales="observed", heads="control"),
    (name="observed", input="observed", fisher="observed", scales="observed", heads="observed"),
    (name="reverse_local_only", input="control", fisher="observed", scales="observed", heads="observed")]
const INPUT_TABLES = ("features_descriptor", "features_split", "score_fwd", "score_bwd",
    "patches_fwd17", "patches_bwd17")
const REFERENCES = ("fisher_cv", "features_predictor", "pred_gmm", "pred_kmeans", "predictions")
sha(path) = open(io->bytes2hex(sha256(io)), path)
key(r) = (r.file,r.lobe)

function settings(path)
    s=TOML.parsefile(path)
    s==Dict("model"=>Dict("method"=>"native_saved_input_learning_counterfactual", "repetitions"=>2,
            "assignment_config"=>"unit_assignment_patch_support.toml"),
        "selection"=>Dict("input_policy"=>"first_observed_fit_execution_all_keys",
            "path"=>[a.name for a in ARMS[1:5]], "reverse_control"=>"reverse_local_only",
            "require_exact_native_reapplication"=>true, "grade_policy"=>"external_after_complete_verification"),
        "preprocessing"=>Dict("local_inputs"=>"saved_unmodified", "normalization"=>"captured_per_scan_per_feature",
            "missing_features"=>"native_abstention_no_imputation")) || error("Only the declared attribution is supported")
    s
end

function parse_cli(args)
    args==["--help"] && return nothing
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option $k")
        if k=="--dry-run"; opts[k]="true"; i+=1; continue; end
        k in ("--input","--settings","--outdir") && i<length(args) && !startswith(args[i+1],"--") ||
            error("Unknown, missing or forbidden option: $k")
        opts[k]=abspath(args[i+1]); i+=2
    end
    all(haskey(opts,k) for k in ("--input","--settings","--outdir")) || error("All paths required")
    isdir(opts["--input"]) || error("Input directory missing")
    (ispath(opts["--outdir"]) || islink(opts["--outdir"])) && error("Output exists")
    settings(opts["--settings"])
    opts
end

function inputs(root)
    counts=Pipeline.selected_counts(joinpath(root,"selected_summary.tsv"))
    cfgpath=joinpath(root,"assignment.toml"); cfg=RU.load_config(cfgpath)
    cfg==RU.load_config(joinpath(ROOT,"config","unit_assignment_patch_support.toml")) ||
        error("Assignment settings differ from the declared native reference")
    files=[joinpath(root,"selected_summary.tsv"),cfgpath]
    for arm in ("control","observed"), stage in (INPUT_TABLES...,REFERENCES...)
        p=joinpath(root,arm,stage*".tsv")
        Pipeline.check_counts(p,counts)
        push!(files,p)
    end
    audits=[p for p in readdir(root;join=true) if startswith(basename(p),"refit_chunk") && endswith(p,".fits.tsv")]
    isempty(audits) && error("Missing observation/reuse audit")
    rows=vcat([last(RU.read_table(p)) for p in audits]...)
    slots=[(r["file"],r["arm"],r["profile"]) for r in rows]
    length(unique(slots))==length(slots)==4length(counts) || error("Incomplete/duplicate fit slots")
    Set(slots)==Set((f,a,p) for f in keys(counts) for a in ("control","observed") for p in ("gaussian","split")) || error("Wrong fit cohort")
    all(r["valid"]=="true" for r in rows) || error("Invalid input geometry")
    groups=Dict{String,String}()
    for f in keys(counts)
        rs=filter(r->r["file"]==f && r["arm"]=="observed",rows)
        reuse=all(r["reused_control"]=="true" for r in rs)
        all((r["reused_control"]=="true")==reuse for r in rs) || error("Mixed reuse")
        all((r["observed_pixels"]==r["total_pixels"])==reuse for r in rs) || error("Mask/reuse mismatch")
        groups[f]=reuse ? "fully_observed" : "partial"
    end
    # Check the local-input control, before any learning or external grading.
    for stage in INPUT_TABLES
        c=last(RU.lobe_table(joinpath(root,"control",stage*".tsv")))
        o=last(RU.lobe_table(joinpath(root,"observed",stage*".tsv")))
        all(groups[k[1]]!="fully_observed" || c[k]==o[k] for k in keys(c)) || error("Local control changed: $stage")
    end
    append!(files,audits)
    (;counts,cfg,cfgpath,groups,hashes=Dict(p=>sha(p) for p in files))
end

"Apply captured statistics without recomputing any target-image statistic."
function frozen_matrix(records, features, states; interactions)
    stat=Dict((r.file,r.feature)=>(r.center,r.scale) for r in states)
    length(stat)==length(states) || error("Duplicate normalization state")
    expected=Set((r.file,f) for r in records for f in features)
    Set(keys(stat))==expected || error("Normalization keys differ")
    z=fill(NaN,length(records),length(features))
    for (i,r) in enumerate(records), (j,f) in enumerate(features)
        center,scale=stat[(r.file,f)]; x=get(r.features,f,NaN)
        (isnan(center) && isnan(scale)) && continue
        isfinite(center) && isfinite(scale) && scale>0 || error("Invalid frozen statistic")
        isfinite(x) && (z[i,j]=(x-center)/scale)
    end
    if interactions && length(features)>1
        extra=hcat([z[:,a].*z[:,b] for a in 1:length(features)-1 for b in a+1:length(features)]...)
        z=hcat(z,extra)
    end
    z,[all(isfinite,r) for r in eachrow(z)]
end

"Score fixed centers; no refit and no physical renaming on the target."
function frozen_kmeans(X, valid, models)
    p=fill(NaN,size(X,1)); votes=zeros(length(p)); n=0
    for m in models
        m.high_cluster==0 && continue
        m.high_cluster in (1,2) || error("Invalid frozen cluster name")
        size(m.centers)==(size(X,2),2) || error("Frozen center shape differs")
        for i in findall(valid)
            assigned=argmin([sum(abs2,X[i,:]-m.centers[:,c]) for c in 1:2])
            votes[i]+=assigned==m.high_cluster
        end
        n+=1
    end
    n>0 && (p[valid]=votes[valid]./n)
    p
end

"Native final scores with fixed means, covariances, weights and physical names."
function frozen_gmm(X, valid, models, opt)
    p=fill(NaN,size(X,1)); votes=zeros(length(p)); n=0
    for m in models
        m.high_cluster==0 && continue
        m.high_cluster in (1,2) || error("Invalid frozen cluster name")
        size(m.means)==(size(X,2),2) || error("Frozen mixture shape differs")
        for i in findall(valid)
            lr=[Gaussian._final_component_score(X[i,:],m.means[:,c],m.covariances[c],m.weights[c],opt) for c in 1:2]
            r=exp.(lr .- maximum(lr)); r./=sum(r)
            all(isfinite,r) || error("Invalid frozen mixture scores")
            votes[i]+=Gaussian._seed_contribution(r,m.high_cluster,opt.seed_aggregation)
        end
        n+=1
    end
    n>0 && (p[valid]=votes[valid]./n)
    p
end

function capture_head(kind,records,features,opt)
    mod=kind=="gmm" ? Gaussian : KMeans
    scales=[]; models=[]
    native=mod._view_probability(records,features,opt;normalization_state=scales,model_state=models)
    X,valid=frozen_matrix(records,features,scales;interactions=opt.interactions)
    expected=kind=="gmm" ? mod._standardized_matrix(records,features;interactions=opt.interactions,
        normalization=opt.normalization,scale_fallback=opt.scale_fallback) :
        mod._standardized_matrix(records,features;interactions=opt.interactions)
    isequal((X,valid),expected) || error("Captured scaling does not reproduce native arithmetic")
    for m in models
        assigned=[kind=="gmm" ? argmax([Gaussian._final_component_score(X[i,:],m.means[:,c],m.covariances[c],m.weights[c],opt) for c in 1:2]) :
            argmin([sum(abs2,X[i,:]-m.centers[:,c]) for c in 1:2]) for i in m.indices]
        assigned==m.assignments || error("Frozen $kind seed $(m.seed) does not reproduce native memberships")
    end
    reapplied=kind=="gmm" ? frozen_gmm(X,valid,models,opt) : frozen_kmeans(X,valid,models)
    isequal(native,reapplied) || error("Captured $kind does not exactly reproduce its native probabilities")
    (;kind,features=copy(features),scales,models,native)
end

function capture_fisher(patches,options)
    options.patch_projection=="none" && options.cv_scheme=="lobe_parity" && options.training_support=="all_admissible" ||
        error("Frozen attribution only supports the declared Fisher path")
    diagnostics=[]; native=EF.cv_scores(patches,options;diagnostics)
    models=Dict{Int,Any}()
    for (parity,d) in enumerate(diagnostics)
        d.fold==(parity==1 ? "even" : "odd") || error("Unexpected fold order")
        d.model===nothing && error("Cannot freeze an unavailable Fisher fold: $(d.reason)")
        models[parity-1]=(w_p=copy(d.model.w_p),mid=copy(d.model.mid),train=copy(d.train),held=copy(d.held),
            converged=d.model.gmm.converged,iterations=d.model.gmm.iterations)
    end
    reapplied=frozen_fisher(patches,models)
    all(isequal(a.score,b.score) && a.invalid_reason==b.invalid_reason for (a,b) in zip(native,reapplied)) ||
        error("Frozen Fisher does not reproduce native scores")
    (;models,native)
end

function frozen_fisher(patches,models)
    Set(keys(models))==Set((0,1)) || error("Both Fisher folds required")
    rows=EF.FisherScore[]
    for (i,k) in enumerate(patches.keys)
        reason=patches.invalid_reasons[i]; value=NaN
        if isempty(reason)
            x=patches.X[i,:]
            if all(isfinite,x) && isfinite(patches.amplitudes[i])
                m=models[1-mod(k[2],2)]
                value=max(dot(x-m.mid,m.w_p),dot(EF.flip_u_disk(x,patches.grid)-m.mid,m.w_p))
                isfinite(value) || (reason="nonfinite_fisher_score")
            else
                reason="nonfinite_patch_input"
            end
        end
        push!(rows,EF.FisherScore(k[1],k[2],value,reason))
    end
    rows
end

function options(features,patches,out,cfgpath,cfg)
    s=cfg["selection"]
    common=["--features",features,"--out",out,"--first-seed",string(s["first_seed"])]
    s["interactions"] && push!(common,"--interactions")
    g=Gaussian._parse_cli(vcat(common,["--config",cfgpath,"--view","v_cc="*join(GFEATURES,','),
        "--seeds",string(s["gmm_seeds"]),"--selftrain",string(s["gmm_selftrain"])]))
    k=KMeans._parse_cli(vcat(common,["--patches",patches,"--seeds",string(s["kmeans_seeds"])]))
    (;g,k)
end

function head_records(features,patches)
    g=Gaussian._load_records(features); k=KMeans._load_records(features)
    KMeans._merge_patches!(k,patches)
    key.(g)==key.(k) || error("Head keys differ")
    (;g,k)
end

function predict_heads(records,scales,heads,opt)
    x,v=frozen_matrix(records.g,GFEATURES,scales.g.scales;interactions=opt.g.interactions)
    g=frozen_gmm(x,v,heads.g.models,opt.g)
    total=zeros(length(g)); counts=zeros(Int,length(g))
    for (_,features) in KVIEWS
        name=join(features,','); s=scales.k[name]; h=heads.k[name]
        x,v=frozen_matrix(records.k,features,s.scales;interactions=opt.k.interactions)
        p=frozen_kmeans(x,v,h.models)
        for i in eachindex(p)
            isfinite(p[i]) || continue
            total[i]+=p[i]; counts[i]+=1
        end
    end
    k=[counts[i]>0 ? total[i]/counts[i] : NaN for i in eachindex(total)]
    (;g,k,kviews=counts)
end

# Human-readable, full Float64 states. No deserialization of executable objects.
function save_bank(path,bank,fisher)
    function pack(h)
        models=[h.kind=="gmm" ? Dict("seed"=>m.seed,"dimension"=>size(m.means,1),"means"=>vec(m.means),
            "covariances"=>vec.(m.covariances),"weights"=>m.weights,"high_cluster"=>m.high_cluster,
            "indices"=>m.indices,"assignments"=>m.assignments) :
            Dict("seed"=>m.seed,"dimension"=>size(m.centers,1),"centers"=>vec(m.centers),
            "high_cluster"=>m.high_cluster,"indices"=>m.indices,"assignments"=>m.assignments,
            "converged"=>m.converged) for m in h.models]
        Dict("kind"=>h.kind,"features"=>h.features,"models"=>models,
            "scales"=>[Dict(string(k)=>v for (k,v) in pairs(s)) for s in h.scales])
    end
    data=Dict("gmm"=>pack(bank.g),"kmeans"=>[pack(bank.k[join(f,',')]) for (_,f) in KVIEWS],
        "fisher"=>[merge(Dict("parity"=>p),Dict(string(k)=>v for (k,v) in pairs(fisher.models[p]))) for p in 0:1])
    open(io->TOML.print(io,data;sorted=true),path,"w")
end

function load_bank(path)
    d=TOML.parsefile(path)
    function unpack(h)
        models=[h["kind"]=="gmm" ? (seed=m["seed"],means=reshape(Float64.(m["means"]),m["dimension"],2),
            covariances=[reshape(Float64.(v),m["dimension"],m["dimension"]) for v in m["covariances"]],
            weights=Float64.(m["weights"]),high_cluster=m["high_cluster"],indices=Int.(m["indices"]),assignments=Int.(m["assignments"])) :
            (seed=m["seed"],centers=reshape(Float64.(m["centers"]),m["dimension"],2),high_cluster=m["high_cluster"],
            indices=Int.(m["indices"]),assignments=Int.(m["assignments"]),converged=m["converged"]) for m in h["models"]]
        scales=[(file=s["file"],feature=s["feature"],center=Float64(s["center"]),scale=Float64(s["scale"])) for s in h["scales"]]
        (;kind=h["kind"],features=String.(h["features"]),scales,models)
    end
    g=unpack(d["gmm"]); k=Dict(join(h["features"],',')=>unpack(h) for h in d["kmeans"])
    fisher=Dict(m["parity"] => (w_p=Float64.(m["w_p"]),mid=Float64.(m["mid"]),train=Int.(m["train"]),
        held=Int.(m["held"]),converged=m["converged"],iterations=m["iterations"]) for m in d["fisher"])
    (;bank=(;g,k),fisher)
end

function joined(root,input,fisher_path,out,cfg)
    dir=joinpath(root,input)
    Pipeline.join_predictor_features(joinpath(dir,"features_descriptor.tsv"),joinpath(dir,"features_split.tsv"),
        joinpath(dir,"score_fwd.tsv"),joinpath(dir,"score_bwd.tsv"),fisher_path,out;
        margin_mode=cfg["model"]["mold_margin_mode"])
end

function compare_tables(a,b,columns)
    ar=last(RU.lobe_table(a;required=columns)); br=last(RU.lobe_table(b;required=columns))
    RU.require_same_keys(ar,br,"saved comparison")
    sum(any(ar[k][c]!=br[k][c] for c in columns) for k in keys(ar))
end

function execute(root,settings_path,out)
    v"1.13"<=VERSION<v"1.14" || error("Julia 1.13 required")
    Threads.nthreads()==1 || error("One Julia thread per independent repetition required")
    BLAS.set_num_threads(1)
    settings(settings_path); data=inputs(root)
    (ispath(out)||islink(out)) && error("Output exists")
    mkpath(out); cp(settings_path,joinpath(out,"settings.toml")); cp(data.cfgpath,joinpath(out,"assignment.toml"))
    cp(joinpath(root,"selected_summary.tsv"),joinpath(out,"selected_summary.tsv"))
    RU.write_table(joinpath(out,"groups.tsv"),["file","group"],
        [Dict("file"=>f,"group"=>data.groups[f]) for f in sort(collect(keys(data.groups)))])
    RU.write_table(joinpath(out,"inputs.tsv"),["path","sha256"],
        [Dict("path"=>p,"sha256"=>h) for (p,h) in sort(collect(data.hashes))])
    opts=EF.load_fisher_config(data.cfg); patches=Dict(); banks=Dict(); fishers=Dict(); native=Dict()
    for arm in ("control","observed")
        println("Capturing native learning: $arm"); flush(stdout)
        dir=joinpath(out,"native_"*arm); mkdir(dir)
        patches[arm]=EF.load_patches(joinpath(root,arm,"patches_fwd17.tsv"),"res",opts)
        f=capture_fisher(patches[arm],opts)
        EF.write_scores(joinpath(dir,"fisher_cv.tsv"),f.native)
        ft=joinpath(dir,"features_predictor.tsv"); joined(root,arm,joinpath(dir,"fisher_cv.tsv"),ft,data.cfg)
        patchpath=joinpath(root,arm,"patches_bwd17.tsv"); rec=head_records(ft,patchpath)
        opt=options(ft,patchpath,joinpath(dir,"unused.tsv"),data.cfgpath,data.cfg)
        g=capture_head("gmm",rec.g,GFEATURES,opt.g)
        k=Dict(join(features,',')=>capture_head("kmeans",rec.k,features,opt.k) for (_,features) in KVIEWS)
        bank=(;g,k); save_bank(joinpath(dir,"models.toml"),bank,f)
        restored=load_bank(joinpath(dir,"models.toml"))
        for h in (bank.g,values(bank.k)...)
            r=h.kind=="gmm" ? restored.bank.g : restored.bank.k[join(h.features,',')]
            isequal(h.scales,r.scales) && isequal(h.models,r.models) || error("Model serialization changed values")
        end
        isequal(f.models,restored.fisher) || error("Fisher serialization changed values")
        banks[arm]=restored.bank; fishers[arm]=restored.fisher
        p=predict_heads(rec,bank,bank,opt)
        isequal(p.g,g.native) || error("Native GMM replay changed")
        native[arm]=p
        Gaussian._write_predictions(joinpath(dir,"pred_gmm.tsv"),rec.g,p.g,Int.(isfinite.(p.g)))
        KMeans._write_predictions(joinpath(dir,"pred_kmeans.tsv"),rec.k,p.k,p.kviews)
        RU.write_soft_vote(ft,joinpath(dir,"pred_kmeans.tsv"),joinpath(dir,"pred_gmm.tsv"),joinpath(dir,"predictions.tsv"),data.cfgpath)
    end
    for a in ARMS
        println("Frozen application: $(a.name)"); flush(stdout)
        dir=joinpath(out,a.name); mkdir(dir)
        EF.write_scores(joinpath(dir,"fisher_cv.tsv"),frozen_fisher(patches[a.input],fishers[a.fisher]))
        ft=joinpath(dir,"features_predictor.tsv"); joined(root,a.input,joinpath(dir,"fisher_cv.tsv"),ft,data.cfg)
        patchpath=joinpath(root,a.input,"patches_bwd17.tsv"); rec=head_records(ft,patchpath)
        opt=options(ft,patchpath,joinpath(dir,"unused.tsv"),data.cfgpath,data.cfg)
        p=predict_heads(rec,banks[a.scales],banks[a.heads],opt)
        Gaussian._write_predictions(joinpath(dir,"pred_gmm.tsv"),rec.g,p.g,Int.(isfinite.(p.g)))
        KMeans._write_predictions(joinpath(dir,"pred_kmeans.tsv"),rec.k,p.k,p.kviews)
        RU.write_soft_vote(ft,joinpath(dir,"pred_kmeans.tsv"),joinpath(dir,"pred_gmm.tsv"),joinpath(dir,"predictions.tsv"),data.cfgpath)
        for stage in REFERENCES; Pipeline.check_counts(joinpath(dir,stage*".tsv"),data.counts); end
        if a.name in ("control","observed")
            for stage in REFERENCES
                read(joinpath(out,"native_"*a.name,stage*".tsv"))==read(joinpath(dir,stage*".tsv")) || error("Native reapplication mismatch: $(a.name)/$stage")
            end
        end
    end
    # Native same-input preservation checks are properties of frozen learning,
    # not empirical benchmark criteria. Nothing is excluded from grading.
    for (a,b) in (("control","local_only"),("observed","reverse_local_only"))
        c=last(RU.lobe_table(joinpath(out,a,"predictions.tsv"))); o=last(RU.lobe_table(joinpath(out,b,"predictions.tsv")))
        all(data.groups[k[1]]!="fully_observed" || c[k]==o[k] for k in keys(c)) || error("Frozen unchanged-input predictions differ")
    end
    # Saved outputs are compared only after every counterfactual has been emitted.
    comparison=Dict{String,String}[]
    for arm in ("control","observed"), stage in REFERENCES
        cols=stage=="fisher_cv" ? ["score","invalid_reason"] : stage=="features_predictor" ?
            first(RU.lobe_table(joinpath(root,arm,stage*".tsv"))) : ["predicted","probability_1","confidence","invalid_reason"]
        n=compare_tables(joinpath(root,arm,stage*".tsv"),joinpath(out,arm,stage*".tsv"),cols)
        push!(comparison,Dict("arm"=>arm,"stage"=>stage,"changed_rows"=>string(n)))
    end
    RU.write_table(joinpath(out,"saved_replay_comparison.tsv"),["arm","stage","changed_rows"],comparison)
    all(sha(p)==h for (p,h) in data.hashes) || error("Inputs changed during diagnostic")
    RU.write_table(joinpath(out,"arms.tsv"),["name","input","fisher","scales","heads"],
        [Dict(string(k)=>v for (k,v) in pairs(a)) for a in ARMS])
    println("Complete: $(length(data.counts)) scans, $(sum(values(data.counts))) lobes, six frozen counterfactuals. External grading remains separate.")
end

function main(args=ARGS)
    o=parse_cli(args)
    o===nothing && return println("diagnose_frozen_learning.jl --input OBSERVED_RUN --settings TOML --outdir NEW_DIR [--dry-run]")
    if haskey(o,"--dry-run")
        d=inputs(o["--input"])
        println("Dry run: $(length(d.counts)) scans / $(sum(values(d.counts))) lobes; six declared counterfactuals; no fitting or output")
        return
    end
    execute(o["--input"],o["--settings"],o["--outdir"])
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && FrozenLearningDiagnostic.main()
