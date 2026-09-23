#!/usr/bin/env julia
# Saved-cost decoder ablations: no image refit, external grade or label inputs.
module MoldStateComparison
module Pipeline
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
end
include(joinpath(@__DIR__,"lib","mold_state_decode.jl"))
using .Pipeline.ReconstructedUnitAssignment, SHA, LinearAlgebra
const M = MoldStateDecode
const INPUTS = ("features.tsv","features_split.tsv","features_local.tsv","features_descriptor.tsv",
    "patches_fwd17.tsv","patches_bwd17.tsv","patches_bwd9.tsv","fisher_cv.tsv","training_support.tsv")
const REPLAY_FILES = ("score_fwd.tsv","score_bwd.tsv","features_predictor.tsv",
    "pred_gmm.tsv","pred_kmeans.tsv","predictions.tsv","summary.tsv")
const OPTIONS = Set(["--reference-dir","--config","--settings","--outdir"])

function parse_cli(args)
    if "--help" in args
        println("run_mold_state_comparison.jl --reference-dir SAVED_TANGENT_DIR --config TOML --settings TOML --outdir NEW_DIR [--dry-run]")
        return nothing
    end
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option: $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in OPTIONS && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden option: $k")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in OPTIONS) || error("All inputs required")
    for k in ("--config","--settings"); isfile(o[k]) || error("Missing $k"); end
    isdir(o["--reference-dir"]) || error("Missing saved reference")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output already exists")
    return o
end

function load_inputs(o)
    ref=o["--reference-dir"]; opt=M.settings(o["--settings"])
    cfg=load_config(o["--config"])
    cfg["model"]["mold_margin_mode"]=="absolute_cost_margin" || error("Frozen absolute margins required")
    cfg["preprocessing"]["patch_residual_filter"]=="smooth_residual" || error("Matched residual input required")
    _,base=lobe_table(joinpath(ref,"features.tsv"))
    for f in INPUTS
        _,rows=lobe_table(joinpath(ref,f)); require_same_keys(base,rows,f)
    end
    paths=[M.audit_paths(ref,v) for v in ("fwd","bwd")]
    views=[M.load_costs(p,base) for p in paths]
    files=sort(unique(first.(keys(base))))
    for f in files; M.records_for(base,f); end
    return (;ref,opt,cfg,base,paths,views,files)
end

function write_scores(dir, data, mode)
    ios=[IOBuffer(),IOBuffer()]; diagnostics=Dict{String,String}[]
    for io in ios; println(io,join(M.Connected.SCORE_HEADER,'\t')); end
    for f in data.files
        records=M.records_for(data.base,f)
        decoded=M.decode(records,data.views,data.opt;mode)
        for j in 1:2
            b=decoded[j]; M.Connected.write_decoded(ios[j],f,records,b)
            push!(diagnostics,Dict("file"=>f,"view"=>j==1 ? "fwd" : "bwd",
                "phase"=>string(mod(b.phase+b.direction*(length(records)-1),2)),
                "mirror"=>string(b.mirror),"file_cost"=>string(b.total),
                "available_lobes"=>string(count(i->all(isfinite,b.costs[i,:]),axes(b.costs,1)))))
        end
    end
    for (v,io) in zip(("fwd","bwd"),ios)
        path=joinpath(dir,"score_"*v*".tsv"); ispath(path) && error("Score output exists")
        write(path,take!(io))
    end
    write_table(joinpath(dir,"state_selection.tsv"),["file","view","phase","mirror","file_cost","available_lobes"],diagnostics)
end

function classify(dir, cfg, config; runner=Pipeline.run_stage)
    table=joinpath(dir,"features_predictor.tsv")
    Pipeline.join_predictor_features(joinpath(dir,"features_descriptor.tsv"),joinpath(dir,"features_split.tsv"),
        joinpath(dir,"score_fwd.tsv"),joinpath(dir,"score_bwd.tsv"),joinpath(dir,"fisher_cv.tsv"),table;
        margin_mode=cfg["model"]["mold_margin_mode"])
    sel=cfg["selection"]
    common=["--features",table,"--first-seed",string(sel["first_seed"])]
    sel["interactions"] && push!(common,"--interactions")
    training=sel["assignment_training_support"]=="complete_patches" ?
        ["--training-support",joinpath(dir,"training_support.tsv")] : String[]
    bootstrap=load_gmm_resampling(cfg).mode=="whole_scans" ?
        ["--bootstrap-audit",joinpath(dir,"gmm_scan_bootstrap.tsv")] : String[]
    gmm=joinpath(dir,"pred_gmm.tsv"); km=joinpath(dir,"pred_kmeans.tsv"); b=Pipeline.BASE4
    runner(dir,"gmm","build_labelfree_gmm_predictions.jl",vcat(common,training,bootstrap,
        ["--config",config,"--out",gmm,"--view","v_cc=$b,patch_u_asym_reconstructed,mold_cc_fwd,mold_cc_bwd,emp_fisher",
         "--seeds",string(sel["gmm_seeds"]),"--selftrain",string(sel["gmm_selftrain"])]))
    runner(dir,"kmeans","build_labelfree_unit_predictions.jl",vcat(common,
        ["--out",km,"--patches",joinpath(dir,"patches_bwd17.tsv"),"--view","v_base=$b",
         "--view","v_split=$b,split_log_skew","--view","v_comt=$b,bwd_neg_com_t",
         "--view","v_diag45=$b,bwd_neg_diag45","--seeds",string(sel["kmeans_seeds"])]))
    predictions=joinpath(dir,"predictions.tsv")
    write_soft_vote(table,km,gmm,predictions,config)
    runner(dir,"validation","validate_unit_predictions.jl",["--predictions",predictions,"--features",joinpath(dir,"features.tsv")])
    Pipeline.write_summary(predictions,joinpath(dir,"summary.tsv"))
end

function execute(o; classifier=classify)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    BLAS.set_num_threads(1)
    data=load_inputs(o)
    haskey(o,"--dry-run") && return println("Three decoder arms: $(length(data.files)) scans / $(length(data.base)) keys; no inference or output")
    out=abspath(o["--outdir"]); mkpath(out)
    inputs=vcat([joinpath(data.ref,f) for f in (INPUTS...,REPLAY_FILES...)],data.paths...,[o["--config"],o["--settings"]])
    write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
        [Dict("input"=>abspath(p),"sha256"=>bytes2hex(sha256(read(p)))) for p in inputs])
    cp(o["--config"],joinpath(out,"assignment_config.toml")); cp(o["--settings"],joinpath(out,"decoder_settings.toml"))
    stage="initialization"
    try
        for mode in (:reference,:finite,:shared)
            stage=string(mode); dir=joinpath(out,stage); mkpath(joinpath(dir,"logs"))
            for f in INPUTS; cp(joinpath(data.ref,f),joinpath(dir,f)); end
            write_scores(dir,data,mode)
            if mode==:reference
                for f in ("score_fwd.tsv","score_bwd.tsv")
                    read(joinpath(dir,f))==read(joinpath(data.ref,f)) || error("Reference score replay differs: $f")
                end
            end
            classifier(dir,data.cfg,abspath(o["--config"]))
            _,pred=lobe_table(joinpath(dir,"predictions.tsv")); require_same_keys(data.base,pred,"complete prediction cohort")
            if mode==:reference
                for f in ("features_predictor.tsv","pred_gmm.tsv","pred_kmeans.tsv","predictions.tsv","summary.tsv")
                    read(joinpath(dir,f))==read(joinpath(data.ref,f)) || error("Reference replay differs: $f")
                end
            end
        end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Three complete decoder arms; benchmark grading remains external")
end
function main(args=ARGS)
    o=parse_cli(args); o===nothing || execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && MoldStateComparison.main()
