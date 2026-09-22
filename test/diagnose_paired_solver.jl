#!/usr/bin/env julia
module PairedSolverDiagnostic
include(joinpath(@__DIR__,"diagnose_paired_convergence.jl"))
include(joinpath(@__DIR__,"lib","paired_solver.jl"))
using TOML, SHA, LinearAlgebra
const D=PairedConvergenceDiagnostic
const S=PairedSolver
const C=S.C
const P=S.P
const OPTIONS=union(D.OPTIONS,Set(["--solver-settings","--comparison-dir"]))
table=D.table
value=D.value
write_records=D.write_records

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); oldargs=String[]; i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option")
        if k=="--dry-run"; opts[k]="true"; push!(oldargs,k); i+=1; continue; end
        k in union(OPTIONS,Set(["--chunk"])) && i<length(args) || error("Forbidden/missing solver option")
        opts[k]=args[i+1]
        k in ("--solver-settings","--comparison-dir") || append!(oldargs,[k,args[i+1]])
        i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All solver inputs required")
    D.parse_cli(oldargs); S.settings(opts["--solver-settings"])
    all(isfile(joinpath(opts["--comparison-dir"],f)) for f in ("fits.tsv","input_hashes.tsv","settings.toml","model_settings.toml")) || error("Incomplete previous comparison")
    return opts
end

function inputs(opts)
    input=D.inputs(opts); solver_options=S.settings(opts["--solver-settings"])
    solver_options.lm_maxiter==input.options.extended_maxiter || error("Matched LM budget differs")
    solver_options.max_time_s==input.options.max_time_s || error("Common wall limit differs")
    comparison=opts["--comparison-dir"]
    for (file,key) in (("settings.toml","--settings"),("model_settings.toml","--model-settings"))
        read(joinpath(comparison,file))==read(opts[key]) || error("Previous settings differ")
    end
    for r in table(joinpath(comparison,"input_hashes.tsv"))
        k=r["input"]
        path=startswith(k,"--") ? opts[k] : joinpath(opts[startswith(k,"native/") ? "--native-dir" : "--reference-dir"],split(k,'/')[2])
        r["sha256"]==bytes2hex(sha256(read(path))) || error("Previous input binding differs")
    end
    oldfits=table(joinpath(comparison,"fits.tsv"))
    length(oldfits)==16 && all(r["file"]==input.file for r in oldfits) || error("Incomplete previous cases")
    oldparams=Dict{String,String}[]
    for (i,case) in enumerate(D.CASES)
        name="$(case.profile).$(case.family).$(case.mode).parameters.tsv"
        append!(oldparams,table(joinpath(comparison,"chunk$(mod1(i,4))",name)))
    end
    return merge(input,(;solver_options,oldfits,oldparams))
end

function input_paths(opts)
    paths=Dict(k=>opts[k] for k in ("--physical-config","--assignment-config","--model-settings","--settings","--solver-settings"))
    for (tag,key,files) in (("native","--native-dir",["fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv"]),
            ("reference","--reference-dir",["fits.tsv","parameters.tsv","native_hashes.tsv","paired_settings.toml"]),
            ("comparison","--comparison-dir",["fits.tsv","input_hashes.tsv","settings.toml","model_settings.toml"]))
        for file in files; paths["$tag/$file"]=joinpath(opts[key],file); end
    end
    for (i,case) in enumerate(D.CASES)
        file="chunk$(mod1(i,4))/$(case.profile).$(case.family).$(case.mode).parameters.tsv"
        paths["comparison/$file"]=joinpath(opts["--comparison-dir"],file)
    end
    return paths
end

function record(out,input,case,ctx,fit,audit,solver,trace,points)
    prefix="$(case.profile).$(case.family).$(case.mode).$solver"
    base=Dict("file"=>input.file,"profile"=>case.profile,"family"=>case.family,"mode"=>case.mode,"solver"=>solver,"N"=>string(fit.result.n))
    row=copy(base)
    for k in (:valid,:reason,:rss,:gcv,:nd,:np,:rss_fwd,:rss_bwd,:peak_fwd,:peak_bwd,:initial_rss,:elapsed_s)
        row[string(k)]=string(getproperty(fit,k))
    end
    for k in (:scale,:gradient_error,:gradient_tolerance,:agreement,:passed,:rank,:condition,:projected_gradient,
            :half_step_projected_gradient,:solver_projected_gradient,:raw_gradient_norm)
        row["stationarity_"*string(k)]=string(getproperty(audit,k))
    end
    row["mean_peak_snr"]=D.F.fmt(fit.result.residual_peak_snr); row["mean_rss"]=D.F.fmt(fit.result.rss)
    row["native_valid"]=string(fit.result.valid); row["threshold"]=D.F.fmt(ctx.cfg.residual_peak_snr_threshold)
    row["amp_min"]=D.F.fmt(fit.result.amp_min); row["amp_range"]=D.F.fmt(fit.result.amp_range)
    row["lm_converged"]=solver=="lm" ? string(fit.converged) : "NA"
    row["iterations"]=solver=="lm" ? string(fit.iterations) : "NA"
    row["solver_return"]=solver=="slsqp" ? fit.status : fit.converged ? "native_small_step_or_gradient" :
        fit.iterations==input.solver_options.lm_maxiter ? "iteration_cap" : "time_cap"
    for k in (:evaluations,:gradient_evaluations,:model_evaluations,:setup_s)
        row[string(k)]=solver=="slsqp" ? string(getproperty(fit,k)) : "NA"
    end
    old=only(r for r in input.oldfits if r["profile"]==case.profile && r["family"]==case.family && r["mode"]==case.mode && r["stage"]=="extended")
    oldp=sort(filter(r->r["profile"]==case.profile && r["family"]==case.family && r["mode"]==case.mode && r["stage"]=="extended",input.oldparams);by=r->parse(Int,r["parameter"]))
    row["previous_lm_rss"]=old["rss"]
    row["previous_lm_parameter_difference"]=D.F.fmt(maximum(abs.(fit.params.-value.(oldp,"value"))))
    parameters=[merge(base,Dict("parameter"=>string(j),"value"=>D.F.fmt(fit.params[j]),"initial"=>D.F.fmt(fit.diagnostic.initial[j]),
        "lower"=>D.F.fmt(fit.lower[j]),"upper"=>D.F.fmt(fit.upper[j]))) for j in eachindex(fit.params)]
    gradients=[merge(base,Dict("parameter"=>string(j),"gradient_step"=>D.F.fmt(audit.first_gradient[j]),
        "gradient_half_step"=>D.F.fmt(audit.gradient[j]),"projected_half_step"=>D.F.fmt(audit.projected[j]))) for j in eachindex(fit.params)]
    residuals=[Dict("row"=>input.pixels[j]["row"],"column"=>input.pixels[j]["column"],
        "mean_prediction"=>D.F.fmt(fit.predictions.mean[j]),"fwd_prediction"=>D.F.fmt(fit.predictions.fwd[j]),"bwd_prediction"=>D.F.fmt(fit.predictions.bwd[j]),
        "mean_residual"=>D.F.fmt(ctx.data.z[j]-fit.predictions.mean[j]),"fwd_residual"=>D.F.fmt(ctx.fwd[j]-fit.predictions.fwd[j]),
        "bwd_residual"=>D.F.fmt(ctx.bwd[j]-fit.predictions.bwd[j])) for j in eachindex(ctx.data.z)]
    checkpoints=[merge(base,Dict("index"=>string(index),"parameter"=>string(j),"value"=>D.F.fmt(points[index][j]))) for index in sort(collect(keys(points))) for j in eachindex(points[index])]
    for (name,rows) in (("fit",[row]),("parameters",parameters),("gradients",gradients),("residuals",residuals),("checkpoints",checkpoints))
        write_records(joinpath(out,"$prefix.$name.tsv"),rows)
    end
    write_records(joinpath(out,"$prefix.trace.tsv"),[Dict(string(k)=>string(getproperty(r,k)) for k in propertynames(r)) for r in trace])
    write_records(joinpath(out,"$prefix.singular_values.tsv"),[Dict("index"=>string(i),"value"=>D.F.fmt(v)) for (i,v) in enumerate(audit.singular_values)])
    return row
end

function execute(opts)
    BLAS.set_num_threads(1); input=inputs(opts)
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/')); cases=D.CASES[chunk[1]:chunk[2]:end]
    for case in cases; D.context(input,case); end
    haskey(opts,"--dry-run") && return println("Solver dry run: $(input.file), $(length(cases)) LM/SLSQP pairs; no optimization/output.")
    out=opts["--outdir"]; mkpath(out)
    for (key,file) in (("--settings","settings.toml"),("--solver-settings","solver_settings.toml"),("--model-settings","model_settings.toml"))
        cp(opts[key],joinpath(out,file))
    end
    write_records(joinpath(out,"input_hashes.tsv"),[Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(p)))) for (k,p) in sort(collect(input_paths(opts));by=first)])
    fits=Dict{String,String}[]; failures=Dict{String,String}[]
    for case in cases
        name="$(case.profile).$(case.family).$(case.mode)"
        try
            ctx=D.context(input,case); println("START $name LM"); flush(stdout)
            lm=P.fit(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,case.mode,
                merge(input.model_options,(local_maxiter=input.solver_options.lm_maxiter,));
                optimizer_options=(maxTime=input.solver_options.max_time_s,x_tol=input.options.x_tol,g_tol=input.options.g_tol),diagnostic=true)
            lt=C.replay_trace(lm,input.options); la=C.stationarity(lm,input.options)
            push!(fits,record(out,input,case,ctx,lm,la,"lm",lt.rows,lt.points))
            println("START $name SLSQP"); flush(stdout)
            problem=S.problem(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,case.mode,input.model_options)
            problem.q0==lm.diagnostic.initial || error("Starts differ")
            sl=S.finish(S.solve(problem,input.solver_options,input.options),ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,case.mode)
            sa=C.stationarity(sl,input.options)
            push!(fits,record(out,input,case,ctx,sl,sa,"slsqp",sl.history,sl.checkpoints))
            println("DONE $name SLSQP=$(sl.status) evaluations=$(sl.evaluations) stationary=$(sa.passed) pg=$(sa.half_step_projected_gradient) valid=$(sl.valid) mean_peak=$(sl.result.residual_peak_snr)"); flush(stdout)
        catch err
            reason=replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
            push!(failures,Dict("case"=>name,"reason"=>reason)); println(stderr,"FAILED $name: $reason"); flush(stderr)
        end
    end
    write_records(joinpath(out,"fits.tsv"),fits); write_records(joinpath(out,"failures.tsv"),failures)
    isempty(failures) || error("Diagnostic failure; no retry, fallback or partial grade")
end

function merge_chunks(out)
    dirs=[joinpath(out,"chunk$i") for i in 1:4]; fits=Dict{String,String}[]
    for dir in dirs
        isfile(joinpath(dir,"failures.tsv")) && error("A shard failed")
        append!(fits,table(joinpath(dir,"fits.tsv")))
    end
    keys=[(r["profile"],r["family"],r["mode"],r["solver"]) for r in fits]
    expected=Set((c.profile,c.family,c.mode,s) for c in D.CASES for s in ("lm","slsqp"))
    length(keys)==length(Set(keys))==16 && Set(keys)==expected || error("Incomplete/duplicate comparison")
    for name in ("input_hashes.tsv","settings.toml","solver_settings.toml","model_settings.toml")
        all(read(joinpath(dir,name))==read(joinpath(first(dirs),name)) for dir in dirs) || error("Shard inputs differ")
        cp(joinpath(first(dirs),name),joinpath(out,name))
    end
    sort!(fits;by=r->(r["profile"],r["family"],r["mode"],r["solver"])); write_records(joinpath(out,"fits.tsv"),fits)
    rows=Dict{String,String}[]
    for solver in ("lm","slsqp"),mode in ("fused","paired"),profile in ("gaussian","split")
        valid=sort(filter(r->r["solver"]==solver && r["mode"]==mode && r["profile"]==profile && r["valid"]=="true",fits);by=r->value(r,"gcv"))
        push!(rows,Dict("solver"=>solver,"mode"=>mode,"profile"=>profile,"valid_families"=>string(length(valid)),
            "selected_family"=>isempty(valid) ? "none" : first(valid)["family"],"selected_stationary"=>isempty(valid) ? "NA" : first(valid)["stationarity_passed"]))
    end
    write_records(joinpath(out,"eligibility.tsv"),rows)
    println("All sixteen solver endpoints accounted for; no recognition grade.")
end
function main(args=ARGS)
    length(args)==2 && first(args)=="--merge" && return merge_chunks(last(args))
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_paired_solver.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--chunk I/N] [--dry-run]; or --merge RUN_DIR")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && PairedSolverDiagnostic.main()
