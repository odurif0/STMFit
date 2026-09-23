#!/usr/bin/env julia
# One saved scan: bounded multistart, then zero/free translation with spatial holdouts.
module PairedShiftDiagnostic
include(joinpath(@__DIR__,"diagnose_paired_convergence.jl"))
include(joinpath(@__DIR__,"lib","paired_shift.jl"))
using .PairedConvergenceDiagnostic: table, value, write_records
using .PairedShift, TOML, SHA, LinearAlgebra
const D=PairedConvergenceDiagnostic
const H=PairedShift
const CASES=[(profile=p,family=f) for p in ("gaussian","split") for f in ("circ","ell")]
const OPTIONS=union(D.OPTIONS,Set(["--solver-settings","--shift-settings"]))

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); old=String[]; i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option")
        if k=="--dry-run"; opts[k]="true"; push!(old,k); i+=1; continue; end
        k in union(OPTIONS,Set(["--chunk"])) && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden pilot option")
        opts[k]=args[i+1]
        k in ("--solver-settings","--shift-settings") || append!(old,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All pilot inputs required")
    D.parse_cli(old); H.S.settings(opts["--solver-settings"]); H.settings(opts["--shift-settings"])
    return opts
end

function inputs(opts)
    input=D.inputs(opts); settings=H.settings(opts["--shift-settings"])
    settings.file==input.file || error("Pilot file differs")
    solver=merge(H.S.settings(opts["--solver-settings"]),(slsqp_maxeval=settings.slsqp_maxeval,
        max_time_s=settings.max_time_s,evaluation_checkpoints=[1,settings.slsqp_maxeval]))
    step=H.grid_step(input.pixels)
    rows=parse.(Int,getindex.(input.pixels,"row")); cols=parse.(Int,getindex.(input.pixels,"column"))
    folds=H.folds(rows,cols,settings.holdout_block_px,settings.holdout_buffer_px)
    all(!isempty(f.train) && !isempty(f.test) for f in folds) || error("Empty spatial fold")
    return merge(input,(;settings,solver,step,folds))
end

function record(out,input,case,ctx,p,fit,audit,mode,fold,seed,scale; extra=NamedTuple())
    id="$(case.profile).$(case.family).$mode.$fold.$seed"
    test=fold==0 ? Int[] : input.folds[fold].test
    pred=fit.predictions; d=ctx.data
    heldout=isempty(test) ? NaN : sum(abs2,ctx.fwd[test].-pred.fwd[test])+sum(abs2,ctx.bwd[test].-pred.bwd[test])
    row=Dict("id"=>id,"file"=>input.file,"profile"=>case.profile,"family"=>case.family,"mode"=>mode,"fold"=>string(fold),"seed"=>string(seed),
        "N"=>string(ctx.initial.n),"valid"=>string(fit.valid),"mean_valid"=>string(fit.result.valid),"mean_reason"=>fit.result.reason,
        "mean_rss"=>string(fit.result.rss),"mean_peak_snr"=>string(fit.result.residual_peak_snr),
        "noise"=>string(d.noise),"threshold"=>string(ctx.cfg.residual_peak_snr_threshold),
        "heldout_rss"=>string(heldout),"heldout_n"=>string(2length(test)),"objective_scale"=>string(scale),
        "shift_x_px"=>mode=="shifted" ? string(fit.params[end-1]) : "0.0",
        "shift_y_px"=>mode=="shifted" ? string(fit.params[end]) : "0.0", "step_x_nm"=>string(input.step[1]),"step_y_nm"=>string(input.step[2]))
    for k in (:rss,:gcv,:nd,:np,:peak_fwd,:peak_bwd,:rss_fwd,:rss_bwd,:initial_rss,:elapsed_s,:setup_s,:status,:evaluations,:gradient_evaluations,:model_evaluations)
        row[string(k)]=string(getproperty(fit,k))
    end
    for k in (:agreement,:passed,:half_step_projected_gradient,:gradient_error,:gradient_tolerance,:rank,:condition)
        row["stationarity_"*string(k)]=string(getproperty(audit,k))
    end
    for (key,value) in pairs(extra); row[string(key)]=string(value); end
    write_records(joinpath(out,"$id.fit.tsv"),[row])
    write_records(joinpath(out,"$id.parameters.tsv"),[Dict("parameter"=>string(j),"value"=>string(fit.params[j]),
        "initial"=>string(fit.diagnostic.initial[j]),"lower"=>string(fit.lower[j]),"upper"=>string(fit.upper[j]),
        "gradient"=>string(audit.gradient[j]),"projected"=>string(audit.projected[j])) for j in eachindex(fit.params)])
    train=Set(p.indices); held=Set(test)
    write_records(joinpath(out,"$id.residuals.tsv"),[Dict("row"=>input.pixels[j]["row"],"column"=>input.pixels[j]["column"],
        "training"=>string(j in train),"heldout"=>string(j in held),
        "mean_prediction"=>string(pred.mean[j]),"fwd_prediction"=>string(pred.fwd[j]),"bwd_prediction"=>string(pred.bwd[j]),
        "mean_residual"=>string(d.z[j]-pred.mean[j]),"fwd_residual"=>string(ctx.fwd[j]-pred.fwd[j]),
        "bwd_residual"=>string(ctx.bwd[j]-pred.bwd[j])) for j in eachindex(d.z)])
    write_records(joinpath(out,"$id.trace.tsv"),[Dict(string(k)=>string(getproperty(r,k)) for k in propertynames(r)) for r in fit.history])
    return row
end

function execute(opts)
    BLAS.set_num_threads(1); input=inputs(opts)
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/')); cases=CASES[chunk[1]:chunk[2]:end]
    for case in cases; D.context(input,case); end
    haskey(opts,"--dry-run") && return println("Shift pilot dry run: $(input.file), $(length(cases)) families, $(length(input.pixels)) unchanged pixels, steps=$(input.step); no fits/output.")
    out=opts["--outdir"]; mkpath(out)
    hashes=Dict{String,String}[]
    for key in ("--physical-config","--assignment-config","--model-settings","--settings","--solver-settings","--shift-settings")
        file=replace(key,"--"=>"")*".toml"; cp(opts[key],joinpath(out,file))
        push!(hashes,Dict("input"=>key,"sha256"=>bytes2hex(sha256(read(opts[key])))))
    end
    for (tag,key,files) in (("native","--native-dir",["fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv"]),
        ("reference","--reference-dir",["fits.tsv","parameters.tsv","native_hashes.tsv","paired_settings.toml"]))
        for file in files
            push!(hashes,Dict("input"=>"$tag/$file","sha256"=>bytes2hex(sha256(read(joinpath(opts[key],file))))))
        end
    end
    write_records(joinpath(out,"input_hashes.tsv"),hashes)
    fits=Dict{String,String}[]; failures=Dict{String,String}[]
    for case in cases
        ctx=D.context(input,case)
        # A fixed schedule, not a retry/selection based on the observed outcome.
        schedule=vcat([(mode=m,fold=0) for m in ("fused","paired","shifted")],
            [(mode=m,fold=f) for f in 1:2 for m in ("paired","shifted")])
        for spec in schedule
            mode,fold=spec.mode,spec.fold
            ix=fold==0 ? eachindex(ctx.data.z) : input.folds[fold].train
            p=H.problem(ctx,mode,input.model_options,input.settings,input.step;indices=ix)
            # Zero-shift and free-shift share exactly this native-start objective scale.
            scale=max(sum(abs2,p.predict(p.q0).-p.target),input.options.objective_scale_floor)
            length(p.target)>length(p.q0) || error("Too few training observations")
            for seed in input.settings.start_seeds
                id="$(case.profile).$(case.family).$mode.$fold.$seed"
                try
                    println("START $id"); flush(stdout)
                    q0=H.start(p,seed,input.settings.start_unit_box_radius)
                    fit=H.finish(H.S.solve(merge(p,(;q0)),input.solver,input.options;objective_scale=scale),p,ctx,mode)
                    audit=H.C.stationarity(merge(fit,(initial_rss=scale,)),input.options)
                    push!(fits,record(out,input,case,ctx,p,fit,audit,mode,fold,seed,scale))
                    println("DONE $id rss=$(fit.rss) valid=$(fit.valid) pg=$(audit.half_step_projected_gradient)"); flush(stdout)
                catch err
                    reason=replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
                    push!(failures,Dict("id"=>id,"reason"=>reason)); println(stderr,"FAILED $id $reason"); flush(stderr)
                end
            end
        end
    end
    write_records(joinpath(out,"fits.tsv"),fits); write_records(joinpath(out,"failures.tsv"),failures)
    isempty(failures) || error("Pilot failures; no fallback, retry or partial grade")
end

function summaries(fits,settings)
    winners=Dict{String,String}[]
    for case in CASES,mode in ("fused","paired","shifted"),fold in (mode=="fused" ? (0,) : (0,1,2))
        rows=filter(r->r["profile"]==case.profile && r["family"]==case.family && r["mode"]==mode && parse(Int,r["fold"])==fold,fits)
        Set(parse.(Int,getindex.(rows,"seed")))==Set(settings.start_seeds) && length(rows)==length(settings.start_seeds) || error("Incomplete start group")
        push!(winners,first(sort(rows;by=r->(value(r,"rss"),parse(Int,r["seed"])))))
    end
    pilots=Dict{String,String}[]
    for case in CASES
        getrow(m,f)=only(r for r in winners if r["profile"]==case.profile && r["family"]==case.family && r["mode"]==m && parse(Int,r["fold"])==f)
        c=getrow("paired",0); s=getrow("shifted",0)
        cv=[getrow("paired",f) for f in 1:2]; sv=[getrow("shifted",f) for f in 1:2]
        gains=[value(cv[i],"heldout_rss")-value(sv[i],"heldout_rss") for i in 1:2]
        interior=all(max(abs(value(r,"shift_x_px")),abs(value(r,"shift_y_px")))<settings.shift_max_px-settings.shift_interior_margin_px for r in vcat([s],sv))
        spread=maximum(abs(value(sv[1],k)-value(sv[2],k)) for k in ("shift_x_px","shift_y_px"))
        eligible=s["valid"]=="true" && value(s,"gcv")<value(c,"gcv") && all(>(0),gains) && interior && spread<=settings.fold_shift_agreement_px
        push!(pilots,Dict("profile"=>case.profile,"family"=>case.family,"control_id"=>c["id"],"shifted_id"=>s["id"],
            "full_valid"=>s["valid"],"full_gcv_gain"=>string(value(c,"gcv")-value(s,"gcv")),"fold1_rss_gain"=>string(gains[1]),
            "fold2_rss_gain"=>string(gains[2]),"shift_interior"=>string(interior),"fold_shift_spread_px"=>string(spread),"eligible"=>string(eligible)))
    end
    return (;winners,pilots)
end

function merge_chunks(out)
    dirs=[joinpath(out,"chunk$i") for i in 1:4]; fits=Dict{String,String}[]
    names=["input_hashes.tsv","physical-config.toml","assignment-config.toml","model-settings.toml","settings.toml","solver-settings.toml","shift-settings.toml"]
    for dir in dirs
        isfile(joinpath(dir,"failures.tsv")) && error("Shard failed")
        append!(fits,table(joinpath(dir,"fits.tsv")))
        all(read(joinpath(dir,n))==read(joinpath(first(dirs),n)) for n in names) || error("Shard inputs differ")
    end
    settings=H.settings(joinpath(first(dirs),"shift-settings.toml"))
    length(fits)==length(unique(getindex.(fits,"id")))==28length(settings.start_seeds) || error("Incomplete pilot")
    summary=summaries(fits,settings)
    for name in names; cp(joinpath(first(dirs),name),joinpath(out,name)); end
    write_records(joinpath(out,"fits.tsv"),fits)
    write_records(joinpath(out,"selected_starts.tsv"),summary.winners)
    write_records(joinpath(out,"pilot.tsv"),summary.pilots)
    eligible=all(any(r["profile"]==p && r["eligible"]=="true" for r in summary.pilots) for p in ("gaussian","split"))
    println("All $(length(fits)) endpoints accounted for. Complete Gaussian/split pilot eligible: $eligible. No benchmark grade.")
end
function main(args=ARGS)
    length(args)==2 && first(args)=="--merge" && return merge_chunks(last(args))
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_paired_shift.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--chunk I/N] [--dry-run]; or --merge RUN_DIR")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && PairedShiftDiagnostic.main()
