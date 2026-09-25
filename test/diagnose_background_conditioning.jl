#!/usr/bin/env julia
# Continuation at saved predicted N, not a count search or a benchmark grade.
module BackgroundConditioningDiagnostic
include(joinpath(@__DIR__,"lib","background_conditioning.jl"))
include(joinpath(@__DIR__,"diagnose_local_sigma_counting.jl"))
using .BackgroundConditioning, LinearAlgebra, Statistics, TOML, SHA
const B=BackgroundConditioning
const C=B.C
const G=B.G
const L=LocalSigmaCounting
const VP=B.VP
const REQUIRED=Set(["--input-root","--settings","--outdir"])
rows(filename)=last(L.P.read_table(filename))
number(r,k)=parse(Float64,r[k])
integer(r,k)=parse(Int,r[k])
digest(filename)=bytes2hex(sha256(read(filename)))
str(x)=replace(string(x),'\n'=>' ','\r'=>' ','\t'=>' ')
record(x)=Dict(string(k)=>str(v) for (k,v) in pairs(x))
function write_rows(filename,rs)
    isempty(rs) && return
    columns=sort(collect(union((Set(keys(r)) for r in rs)...)))
    L.P.write_table(filename,columns,[Dict(k=>str(get(r,k,"NA")) for k in columns) for r in rs])
end

function options(args)
    args==["--help"] && return nothing
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in union(REQUIRED,Set(["--chunk"])) && i<length(args) && !startswith(args[i+1],"--") || error("Forbidden/missing option")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in REQUIRED) || error("All inputs required")
    isdir(o["--input-root"]) && isfile(o["--settings"]) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2] && chunk[2] in (1,4) || error("One or four shards required")
    B.settings(o["--settings"])
    o
end

"Select four fixed ranks using ONLY repeatability of saved unlabelled GCV fits."
function inputs(o)
    settings=B.settings(o["--settings"])
    roots=[joinpath(o["--input-root"],"chunk$i") for i in 1:4]
    selected=Dict{String,String}[]; candidates=Dict{String,String}[]; contexts=Dict{String,String}[]
    sources=Dict{String,String}(); hashes=Dict("settings"=>digest(o["--settings"]))
    cohort=String[]
    for (i,root) in enumerate(roots)
        isfile(joinpath(root,"failures.tsv")) && error("Saved source has failures")
        for name in ("selected.tsv","candidates.tsv","context.tsv","cohort.tsv","count.toml","settings.toml","global_sigma_max.ell.toml")
            hashes["chunk$i/$name"]=digest(joinpath(root,name))
        end
        for r in rows(joinpath(root,"cohort.tsv"))
            r["file"] in cohort && error("Duplicate source file")
            push!(cohort,r["file"]); sources[r["file"]]=root
        end
        append!(selected,rows(joinpath(root,"selected.tsv")))
        append!(candidates,rows(joinpath(root,"candidates.tsv")))
        append!(contexts,rows(joinpath(root,"context.tsv")))
    end
    length(cohort)==settings.raw["selection"]["cohort_size"] || error("Incomplete inference cohort")
    ranking=Dict{String,String}[]
    for file in sort(cohort)
        pair=[only(r for r in selected if r["file"]==file && r["arm"]=="global_sigma_max" && integer(r,"repetition")==rep) for rep in 1:2]
        all(r->r["status"]=="ok" && r["source"]=="ell",pair) || error("Unavailable/nonelliptical saved selection")
        pair[1]["N_selected"]==pair[2]["N_selected"] || error("Saved N differs between repeats")
        a,b=number.(pair,"gcv")
        all(isfinite,(a,b)) && min(a,b)>0 || error("Invalid saved GCV")
        push!(ranking,Dict("file"=>file,"N"=>pair[1]["N_selected"],"family"=>"ell",
            "repeat1_gcv"=>str(a),"repeat2_gcv"=>str(b),"relative_difference"=>str(abs(a-b)/max(a,b))))
    end
    sort!(ranking;by=r->(number(r,"relative_difference"),r["file"]))
    for (i,r) in enumerate(ranking); r["rank"]=str(i); end
    chosen=ranking[settings.raw["selection"]["case_ranks"]]
    cases=NamedTuple[]
    for (i,r) in enumerate(chosen)
        file=r["file"]; root=sources[file]; n=integer(r,"N")
        candidate=only(c for c in candidates if c["file"]==file && c["arm"]=="global_sigma_max" &&
            integer(c,"repetition")==1 && integer(c,"N")==n && c["family"]=="ell")
        candidate["success"]==candidate["valid"]=="true" || error("Invalid saved start")
        number(candidate,"gcv")==number(r,"repeat1_gcv") || error("Selection/candidate mismatch")
        context=only(c for c in contexts if c["file"]==file)
        cfg=L.configs(TOML.parsefile(joinpath(root,"count.toml")),L.settings(joinpath(root,"settings.toml")),"global_sigma_max").ell
        snapshot=Dict(string(k)=>(getfield(cfg,k) isa Symbol ? string(getfield(cfg,k)) : getfield(cfg,k)) for k in fieldnames(typeof(cfg)))
        isequal(snapshot,TOML.parsefile(joinpath(root,"global_sigma_max.ell.toml"))) || error("Physical configuration changed")
        dataname="fit_data/$file.tsv"
        hashes["$(basename(root))/$dataname"]=digest(joinpath(root,dataname))
        push!(cases,(;index=i,file,n,root,candidate,context,cfg,ranking=r))
    end
    (;settings,ranking,cases,hashes)
end

function context(case)
    r=case.context
    pixels=rows(joinpath(case.root,"fit_data",case.file*".tsv"))
    x=number.(pixels,"x_nm"); y=number.(pixels,"y_nm"); z=number.(pixels,"z_nm")
    length(z)==integer(r,"n_data") && all(isfinite,vcat(x,y,z)) || error("Invalid saved pixels")
    axis=(origin=(number(r,"origin_x_nm"),number(r,"origin_y_nm")),axis=[number(r,"axis_x"),number(r,"axis_y")],
        perp=[-number(r,"axis_y"),number(r,"axis_x")],tmin=number(r,"support_tmin_nm"),tmax=number(r,"support_tmax_nm"))
    # GCV finalization does not access these unused grid slots or refit anything.
    data=(;x,y,z,zfull=z,xs=Float64[],ys=Float64[],zimg=zeros(0,0),noise=number(r,"noise"),axisctx=axis)
    p0=parse.(Float64,split(case.candidate["parameters"],';'))
    lo,hi=VP.native_raw_bounds(case.n,case.cfg)
    length(p0)==length(lo) && all(lo.<=p0.<=hi) || error("Invalid saved parameter vector")
    amin=number(case.candidate,"amp_min"); arange=number(case.candidate,"amp_range")
    predict=p->G._chain_model_values(x,y,p,case.n,axis,case.cfg;amp_min=amin,amp_range=arange)
    problem=(q0=p0,lo,hi,predict,target=z)
    initial=G.ChainModelResult(n=case.n,params=copy(p0),success=true,amp_min=amin,amp_range=arange)
    VP.finalize_native!(initial,data,case.cfg)
    isapprox(initial.rss,number(case.candidate,"rss");rtol=1e-10,atol=1e-12) && initial.valid || error("Saved starting fit did not replay")
    (;problem,data,initial)
end

function save_fit(out,case,ctx,fit,audit,arm,repetition)
    prefix="$(case.file).repeat$repetition.$arm"
    base=Dict("file"=>case.file,"arm"=>arm,"repetition"=>str(repetition),"N"=>str(case.n),"family"=>"ell")
    native=G.ChainModelResult(n=case.n,params=copy(fit.params),success=true,amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
    VP.finalize_native!(native,ctx.data,case.cfg)
    summary=copy(base)
    for k in (:rss,:initial_rss,:elapsed_s,:setup_s,:solve_total_s,:audit_s,:status,:converged,:iterations,:evaluations,:gradient_evaluations,
            :model_evaluations,:infeasible_trials,:raw_endpoint_violation,:endpoint_roundoff)
        summary[string(k)]=str(getproperty(fit,k))
    end
    for k in (:valid,:reason,:gcv,:overlap,:residual_peak_snr,:endpoint_overrun_nm)
        summary[string(k)]=str(getproperty(native,k))
    end
    for k in (:scale,:gradient_error,:gradient_tolerance,:agreement,:passed,:rank,:condition,
            :projected_gradient,:half_step_projected_gradient,:solver_projected_gradient,:raw_gradient_norm)
        summary["stationarity_"*string(k)]=str(getproperty(audit,k))
    end
    J=C.finite_jacobian(ctx.problem.predict,fit.params,fit.lower,fit.upper,
        case.options.derivative_relative_step*case.options.finite_difference_half_step_factor)
    unit=B.coordinates(ctx.problem,ctx.data.x,ctx.data.y,"unit_box_slsqp",case.options)
    centered=B.coordinates(ctx.problem,ctx.data.x,ctx.data.y,"centered_plane_slsqp",case.options)
    for (name,map) in (("unit",unit),("centered",centered))
        s=svdvals(J*map.M)
        summary[name*"_jacobian_condition"]=str(first(s)/last(s))
    end
    pars=[merge(base,Dict("index"=>str(j),"value"=>str(fit.params[j]),"initial"=>str(ctx.problem.q0[j]),
        "lower"=>str(fit.lower[j]),"upper"=>str(fit.upper[j]),"gradient_step"=>str(audit.first_gradient[j]),
        "gradient_half_step"=>str(audit.gradient[j]),"projected_half_step"=>str(audit.projected[j]))) for j in eachindex(fit.params)]
    pred=ctx.problem.predict(fit.params)
    pixels=[Dict("index"=>str(j),"prediction"=>str(pred[j]),"residual"=>str(ctx.data.z[j]-pred[j])) for j in eachindex(pred)]
    points=[Dict("checkpoint"=>str(k),"index"=>str(j),"value"=>str(v[j])) for (k,v) in sort(collect(fit.checkpoints);by=first) for j in eachindex(v)]
    _,lobes,ts,us,sp,sq=G._decode_chain(fit.params,case.n,ctx.data.axisctx,case.cfg;amp_min=native.amp_min,amp_range=native.amp_range)
    ls=[Dict("lobe"=>str(j),"amplitude"=>str(lobes[j].amplitude),"x_nm"=>str(lobes[j].x_nm),"y_nm"=>str(lobes[j].y_nm),
        "t_nm"=>str(ts[j]),"u_nm"=>str(us[j]),"sigma_parallel_nm"=>str(sp[j]),"sigma_perp_nm"=>str(sq[j])) for j in 1:case.n]
    for (name,rs) in (("fit",[summary]),("parameters",pars),("pixels",pixels),("checkpoints",points),("lobes",ls))
        write_rows(joinpath(out,"$prefix.$name.tsv"),rs)
    end
    write_rows(joinpath(out,"$prefix.trace.tsv"),record.(fit.history))
    summary
end

function initial_diagnostics(ctx,options,audit_options)
    p=ctx.problem.q0; lo=ctx.problem.lo; hi=ctx.problem.hi; predict=ctx.problem.predict
    J=C.finite_jacobian(predict,p,lo,hi,options.derivative_relative_step)
    rss=sum(abs2,predict(p).-ctx.data.z)
    a=C.stationarity((params=p,lower=lo,upper=hi,initial_rss=rss,
        diagnostic=(;predict,target=ctx.data.z,jacobian=J)),audit_options)
    result=Dict("rss"=>str(rss),"stationarity_passed"=>str(a.passed),"stationarity_agreement"=>str(a.agreement),
        "projected_gradient"=>str(a.half_step_projected_gradient))
    Jhalf=C.finite_jacobian(predict,p,lo,hi,options.derivative_relative_step*options.finite_difference_half_step_factor)
    for (name,mode) in (("unit","unit_box_slsqp"),("centered","centered_plane_slsqp"))
        map=B.coordinates(ctx.problem,ctx.data.x,ctx.data.y,mode,options)
        singular=svdvals(Jhalf*map.M)
        result[name*"_jacobian_condition"]=str(first(singular)/last(singular))
        result[name*"_plane_condition"]=str(sqrt(cond(map.plane_gram)))
    end
    result
end

function execute(o)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    Threads.nthreads()==1 || error("One Julia thread per shard required")
    BLAS.set_num_threads(1)
    input=inputs(o); chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    cases=input.cases[chunk[1]:chunk[2]:end]
    if haskey(o,"--dry-run")
        for c in cases; println("Metadata only: case=$(c.index) file=$(c.file) saved_N=$(c.n) rank=$(c.ranking["rank"]) drift=$(c.ranking["relative_difference"])"); end
        return
    end
    out=o["--outdir"]; (ispath(out)||islink(out)) && error("Output exists"); mkpath(out)
    cp(o["--settings"],joinpath(out,"settings.toml"))
    write_rows(joinpath(out,"input_hashes.tsv"),[Dict("input"=>k,"sha256"=>v) for (k,v) in sort(collect(input.hashes);by=first)])
    write_rows(joinpath(out,"ranking.tsv"),input.ranking)
    write_rows(joinpath(out,"cases.tsv"),[merge(copy(c.ranking),Dict("case"=>str(c.index))) for c in input.cases])
    fits=Dict{String,String}[]; failures=Dict{String,String}[]
    for item in cases
        case=merge(item,(options=input.settings.options,)); ctx=context(case)
        # Keep the original physical snapshot and fit pixels with each case.
        cp(joinpath(case.root,"global_sigma_max.ell.toml"),joinpath(out,case.file*".physical.toml"))
        cp(joinpath(case.root,"fit_data",case.file*".tsv"),joinpath(out,case.file*".input_pixels.tsv"))
        write_rows(joinpath(out,case.file*".context.tsv"),[case.context])
        write_rows(joinpath(out,case.file*".initial.tsv"),[initial_diagnostics(ctx,case.options,input.settings.audit)])
        for rep in 1:case.options.repetitions, arm in (isodd(rep) ? B.ARMS : reverse(B.ARMS))
            println("START $(case.file) repeat=$rep $arm"); flush(stdout)
            try
                started=time_ns()
                fit=B.solve(ctx.problem,ctx.data.x,ctx.data.y,arm,case.options,input.settings.audit)
                solve_total_s=(time_ns()-started)/1e9; audit_started=time_ns()
                audit=C.stationarity(fit,input.settings.audit)
                fit=merge(fit,(;solve_total_s,audit_s=(time_ns()-audit_started)/1e9))
                row=save_fit(out,case,ctx,fit,audit,arm,rep); push!(fits,row)
                println("DONE $(case.file) repeat=$rep $arm status=$(fit.status) rss=$(fit.rss) stationary=$(audit.passed) pg=$(audit.half_step_projected_gradient) valid=$(row["valid"])"); flush(stdout)
            catch err
                push!(failures,Dict("file"=>case.file,"repetition"=>str(rep),"arm"=>arm,"error"=>str(sprint(showerror,err))))
                println(stderr,"FAILED $(case.file) repeat=$rep $arm: ",sprint(showerror,err)); flush(stderr)
            end
        end
    end
    write_rows(joinpath(out,"fits.tsv"),fits); write_rows(joinpath(out,"failures.tsv"),failures)
    isempty(failures) || error("Numerical comparison has failures; keep outputs, no retry or grade")
end

function main(args=ARGS)
    o=options(args)
    o===nothing ? println("diagnose_background_conditioning.jl --input-root SAVED_COUNT_RUN --settings TOML --outdir NEW [--chunk I/4] [--dry-run]\nFour numerical-repeatability ranks; fixed predicted N/family, no labels or grade.") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && BackgroundConditioningDiagnostic.main()
