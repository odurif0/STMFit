#!/usr/bin/env julia
# Fixed saved pixels and starts; no raw preprocessing, labels, N search or grade.
module PairedConvergenceDiagnostic
include(joinpath(@__DIR__,"diagnose_registered_fit_failure.jl"))
include(joinpath(@__DIR__,"lib","paired_convergence.jl"))
using .RegisteredFailureDiagnostic: G, F, read_table, write_table
using .PairedConvergence
using TOML, SHA, LinearAlgebra
const C=PairedConvergence
const P=C.P
const OPTIONS=Set(["--native-dir","--reference-dir","--physical-config","--assignment-config","--model-settings","--settings","--outdir"])
const CASES=[(profile=profile,family=family,mode=mode) for profile in ("gaussian","split") for family in ("circ","ell") for mode in ("fused","paired")]
table(path)=last(read_table(path))
value(row,k)=parse(Float64,row[k])

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option")
        if k=="--dry-run"; opts[k]="true"; i+=1; continue; end
        k in union(OPTIONS,Set(["--chunk"])) && i<length(args) || error("Missing/forbidden convergence option")
        opts[k]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All convergence inputs required")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    for k in ("--physical-config","--assignment-config","--model-settings","--settings")
        isfile(opts[k]) || error("Missing $k")
    end
    for (k,files) in (("--native-dir",["fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv"]),
            ("--reference-dir",["fits.tsv","parameters.tsv","native_hashes.tsv","paired_settings.toml"]))
        all(isfile(joinpath(opts[k],f)) for f in files) || error("Incomplete $k")
    end
    (ispath(opts["--outdir"]) || islink(opts["--outdir"])) && error("Output exists")
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2]<=4 || error("Invalid chunk")
    C.settings(opts["--settings"]); P.settings(opts["--model-settings"])
    return opts
end

function inputs(opts)
    native,reference=opts["--native-dir"],opts["--reference-dir"]
    raw=TOML.parsefile(opts["--physical-config"]); assignment=TOML.parsefile(opts["--assignment-config"])
    hashes=table(joinpath(native,"input_hashes.tsv"))
    for (old,new) in (("--config","--physical-config"),("--assignment-config","--assignment-config"))
        row=only(r for r in hashes if r["input"]==old)
        row["sha256"]==bytes2hex(sha256(read(opts[new]))) || error("Native configuration differs")
    end
    for row in table(joinpath(reference,"native_hashes.tsv"))
        row["sha256"]==bytes2hex(sha256(read(joinpath(native,row["input"])))) || error("Native/reference binding differs")
    end
    read(opts["--model-settings"])==read(joinpath(reference,"paired_settings.toml")) || error("Paired physical settings differ")
    options=C.settings(opts["--settings"]); model_options=P.settings(opts["--model-settings"])
    options.control_maxiter==model_options.local_maxiter || error("Control budget must replay previous model")
    fits=table(joinpath(native,"fits.tsv")); params=table(joinpath(native,"parameters.tsv"))
    previous=table(joinpath(reference,"fits.tsv")); previous_params=table(joinpath(reference,"parameters.tsv"))
    file=TOML.parsefile(opts["--settings"])["selection"]["file"]
    length(fits)==length(previous)==8 && all(r["file"]==file for r in vcat(fits,previous)) || error("One declared complete scan required")
    pixels=table(joinpath(native,"registered.pixels.tsv"))
    length(Set((r["row"],r["column"]) for r in pixels))==length(pixels) || error("Duplicate saved pixel")
    return (;raw,assignment,options,model_options,fits,params,previous,previous_params,pixels,file)
end

function context(input,case)
    profile,family=case.profile,case.family
    rows=filter(r->r["arm"]=="registered" && r["profile"]==profile && r["family"]==family,input.fits)
    row=only(rows); n=parse(Int,row["N"]); row["success"]=="true" || error("Missing numerical start")
    model=deepcopy(input.raw["model"]); model["peak_profile"]=profile
    profile=="split" && (model["skew_ratio_max"]=input.assignment["model"]["split_skew_ratio_max"])
    _,ec,cc=F.Extractor._configs(model,input.raw["preprocessing"],"unused"); cfg=family=="circ" ? cc : ec
    cfg.cv_method=="gcv" || error("GCV required; saved pixels cannot run kfold")
    params=sort(filter(r->r["arm"]=="registered" && r["profile"]==profile && r["family"]==family && r["stage"]=="final",input.params);by=r->parse(Int,r["parameter"]))
    length(params)==G._chain_nparams(n,cfg) && all(r["start"]=="1" for r in params) || error("Invalid saved parameter layout")
    p=value.(params,"value")
    initial=G.ChainModelResult(n=n,params=p,success=true,amp_min=value(row,"amp_min"),amp_range=value(row,"amp_range"))
    ax,ay=value(row,"axis_x"),value(row,"axis_y")
    axis=(origin=(value(row,"origin_x_nm"),value(row,"origin_y_nm")),axis=[ax,ay],perp=[-ay,ax],
        tmin=value(row,"support_tmin"),tmax=value(row,"support_tmax"))
    x=value.(input.pixels,"x_nm"); y=value.(input.pixels,"y_nm"); z=value.(input.pixels,"z_nm")
    # GCV finalization uses only observed x/y/z, noise and the frozen axis.
    # Grid fields are deliberately empty: neither an image nor missing pixels
    # are reconstructed, and no initialization/ROI operation is called.
    data=(x=x,y=y,z=z,zfull=z,noise=value(row,"noise"),axisctx=axis,xs=Float64[],ys=Float64[],zimg=zeros(0,0))
    length(z)==parse(Int,row["n_pixels"]) || error("Saved pixel count differs")
    fwd=value.(input.pixels,"fwd_nm").-value(row,"offset"); bwd=value.(input.pixels,"bwd_nm").-value(row,"offset")
    all(isfinite,z) && all(isfinite,fwd) && all(isfinite,bwd) || error("Nonfinite saved objective pixel")
    isapprox((fwd.+bwd)./2,z;rtol=1e-12,atol=1e-14) || error("Saved mean differs")
    return (;initial,data,cfg,fwd,bwd,row)
end

function write_records(path,rows)
    isempty(rows) || write_table(path,sort(collect(keys(first(rows)))),rows)
end

function record_case(out,input,case,pair,ctx)
    profile,family,mode=case.profile,case.family,case.mode
    prefix="$profile.$family.$mode"
    summaries=Dict{String,String}[]; parameters=Dict{String,String}[]; gradients=Dict{String,String}[]
    traces=Dict{String,String}[]; checkpoints=Dict{String,String}[]
    for stage in ("control","extended")
        fit=getproperty(pair,Symbol(stage)); audit=getproperty(pair,Symbol(stage*"_stationarity")); trace=getproperty(pair,Symbol(stage*"_trace"))
        base=Dict("file"=>input.file,"profile"=>profile,"family"=>family,"mode"=>mode,"stage"=>stage,"N"=>string(fit.result.n))
        summary=copy(base)
        for k in (:valid,:reason,:rss,:gcv,:nd,:np,:rss_fwd,:rss_bwd,:peak_fwd,:peak_bwd,:initial_rss,:elapsed_s,:converged,:iterations)
            summary[string(k)]=string(getproperty(fit,k))
        end
        for k in (:scale,:gradient_error,:gradient_tolerance,:agreement,:passed,:rank,:condition,:projected_gradient,
                :half_step_projected_gradient,:solver_projected_gradient,:raw_gradient_norm)
            summary["stationarity_"*string(k)]=string(getproperty(audit,k))
        end
        summary["mean_peak_snr"]=F.fmt(fit.result.residual_peak_snr); summary["mean_rss"]=F.fmt(fit.result.rss)
        summary["native_valid"]=string(fit.result.valid); summary["threshold"]=F.fmt(ctx.cfg.residual_peak_snr_threshold)
        summary["amp_min"]=F.fmt(fit.result.amp_min); summary["amp_range"]=F.fmt(fit.result.amp_range)
        summary["native_convergence_reason"]=fit.converged ? "native_small_step_or_gradient" :
            fit.iterations==(stage=="control" ? input.options.control_maxiter : input.options.extended_maxiter) ? "iteration_cap" : "time_cap"
        previous=only(r for r in input.previous if r["profile"]==profile && r["family"]==family && r["mode"]==mode)
        summary["previous_300_rss"]=previous["rss"]
        push!(summaries,summary)
        for j in eachindex(fit.params)
            push!(parameters,merge(base,Dict("parameter"=>string(j),"value"=>F.fmt(fit.params[j]),
                "initial"=>F.fmt(fit.diagnostic.initial[j]),"lower"=>F.fmt(fit.lower[j]),"upper"=>F.fmt(fit.upper[j]))))
            push!(gradients,merge(base,Dict("parameter"=>string(j),"gradient_step"=>F.fmt(audit.first_gradient[j]),
                "gradient_half_step"=>F.fmt(audit.gradient[j]),"projected_half_step"=>F.fmt(audit.projected[j]))))
        end
        for r in trace.rows
            push!(traces,merge(base,Dict(string(k)=>string(getproperty(r,k)) for k in propertynames(r))))
        end
        for iteration in sort(collect(keys(trace.points)))
            q=trace.points[iteration]
            for j in eachindex(q)
                push!(checkpoints,merge(base,Dict("iteration"=>string(iteration),"parameter"=>string(j),"value"=>F.fmt(q[j]))))
            end
        end
        residuals=[Dict("row"=>input.pixels[j]["row"],"column"=>input.pixels[j]["column"],
            "mean_prediction"=>F.fmt(fit.predictions.mean[j]),"fwd_prediction"=>F.fmt(fit.predictions.fwd[j]),
            "bwd_prediction"=>F.fmt(fit.predictions.bwd[j]),"mean_residual"=>F.fmt(ctx.data.z[j]-fit.predictions.mean[j]),
            "fwd_residual"=>F.fmt(ctx.fwd[j]-fit.predictions.fwd[j]),"bwd_residual"=>F.fmt(ctx.bwd[j]-fit.predictions.bwd[j])) for j in eachindex(ctx.data.z)]
        write_records(joinpath(out,"$prefix.$stage.residuals.tsv"),residuals)
        write_records(joinpath(out,"$prefix.$stage.singular_values.tsv"),[Dict("index"=>string(i),"value"=>F.fmt(v)) for (i,v) in enumerate(audit.singular_values)])
    end
    for (name,rows) in (("fits",summaries),("parameters",parameters),("gradients",gradients),("trace",traces),("checkpoints",checkpoints))
        write_records(joinpath(out,"$prefix.$name.tsv"),rows)
    end
    return summaries
end

function execute(opts)
    BLAS.set_num_threads(1)
    input=inputs(opts)
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    cases=CASES[chunk[1]:chunk[2]:end]
    # Metadata and saved-pixel validation is cheap; no optimizer during dry run.
    for case in cases; context(input,case); end
    haskey(opts,"--dry-run") && return println("Convergence dry run: $(input.file), $(length(cases)) cases, $(input.options.control_maxiter)/$(input.options.extended_maxiter) iterations; no optimization/output.")
    out=opts["--outdir"]; mkpath(out)
    cp(opts["--settings"],joinpath(out,"settings.toml")); cp(opts["--model-settings"],joinpath(out,"model_settings.toml"))
    hashes=Dict{String,String}[]
    for k in ("--physical-config","--assignment-config","--model-settings","--settings")
        push!(hashes,Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))))
    end
    for (tag,dir,files) in (("native",opts["--native-dir"],["fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv"]),
            ("reference",opts["--reference-dir"],["fits.tsv","parameters.tsv","native_hashes.tsv","paired_settings.toml"]))
        for file in files; push!(hashes,Dict("input"=>"$tag/$file","sha256"=>bytes2hex(sha256(read(joinpath(dir,file)))))); end
    end
    write_records(joinpath(out,"input_hashes.tsv"),hashes)
    summaries=Dict{String,String}[]; failures=Dict{String,String}[]
    for case in cases
        prefix="$(case.profile).$(case.family).$(case.mode)"
        try
            ctx=context(input,case)
            println("START $prefix"); flush(stdout)
            pair=C.fit_pair(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,case.mode,input.model_options,input.options)
            append!(summaries,record_case(out,input,case,pair,ctx))
            a=pair.extended_stationarity; r=pair.extended
            println("DONE $prefix iterations=$(r.iterations) LM=$(r.converged) stationarity=$(a.passed) pg=$(a.half_step_projected_gradient) valid=$(r.valid) mean_peak=$(r.result.residual_peak_snr)"); flush(stdout)
        catch err
            reason=replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
            push!(failures,Dict("case"=>prefix,"reason"=>reason)); println(stderr,"FAILED $prefix: $reason"); flush(stderr)
        end
    end
    write_records(joinpath(out,"fits.tsv"),summaries); write_records(joinpath(out,"failures.tsv"),failures)
    isempty(failures) || error("Explicit diagnostic failures; no retry or partial conclusion")
end

function main(args=ARGS)
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_paired_convergence.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--chunk I/N] [--dry-run]")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && PairedConvergenceDiagnostic.main()
