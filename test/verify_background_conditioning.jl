#!/usr/bin/env julia
# Saved-output checks only; no optimizer, labels, count selection or grade.
module BackgroundConditioningVerification
include(joinpath(@__DIR__,"diagnose_background_conditioning.jl"))
using .BackgroundConditioningDiagnostic, Test, LinearAlgebra, Statistics, TOML
const D=BackgroundConditioningDiagnostic
const B=D.B
const C=B.C
n=D.number
i=D.integer
rows=D.rows
close(a,b,opts)=isapprox(a,b;atol=opts.replay_atol,rtol=opts.replay_rtol)

"Independent free-gap, globally oriented elliptical Gaussian forward model."
function independent(p,case,ctx)
    cfg=case.cfg; N=case.n; axis=ctx.data.axisctx
    cfg.chain_spacing_model=="free" && cfg.chain_peak_orientation=="global" &&
        cfg.overlap_constraint=="global_sigma_max" && !cfg.chain_circular_sigmas &&
        cfg.shared_sigma_types==0 && cfg.peak_profile==:gaussian || error("Unsupported independent model")
    logistic=v->1/(1+exp(-v))
    amplitudes=ctx.initial.amp_min .+ ctx.initial.amp_range .* logistic.(p[4:3+N])
    start_raw=p[4+N]; k=5+N
    floor_gap=max(cfg.spacing_min_nm,sqrt(-2log(cfg.max_overlap))*cfg.sigma_parallel_max_nm)
    support=axis.tmax-axis.tmin; gaps=Float64[]; used=0.
    for j in 1:N-1
        maximum_gap=max(floor_gap,min(cfg.spacing_max_nm,support-used-(N-1-j)*floor_gap))
        gap=floor_gap+(maximum_gap-floor_gap)*logistic(p[k]); push!(gaps,gap); used+=gap; k+=1
    end
    first_t=axis.tmin+max(axis.tmax-used-axis.tmin,D.G.EPS)*logistic(start_raw)
    ts=vcat(first_t,first_t.+cumsum(gaps))
    us=cfg.lateral_max_nm.*tanh.(p[k:k+N-1]); k+=N
    sp=cfg.sigma_parallel_min_nm .+ (cfg.sigma_parallel_max_nm-cfg.sigma_parallel_min_nm).*logistic.(p[k:k+N-1]); k+=N
    sq=cfg.sigma_perp_min_nm .+ (cfg.sigma_perp_max_nm-cfg.sigma_perp_min_nm).*logistic.(p[k:k+N-1]); k+=N
    k==length(p)+1 || error("Parameter accounting differs")
    xs=axis.origin[1].+axis.axis[1].*ts.+axis.perp[1].*us
    ys=axis.origin[2].+axis.axis[2].*ts.+axis.perp[2].*us
    pred=p[1].+p[2].*ctx.data.x.+p[3].*ctx.data.y
    for j in 1:N
        dt=(ctx.data.x.-xs[j]).*axis.axis[1].+(ctx.data.y.-ys[j]).*axis.axis[2]
        du=(ctx.data.x.-xs[j]).*axis.perp[1].+(ctx.data.y.-ys[j]).*axis.perp[2]
        pred .+= amplitudes[j] .* exp.(-0.5 .* ((dt ./ sp[j]).^2 .+ (du ./ sq[j]).^2))
    end
    (;pred,amplitudes,xs,ys,ts,us,sp,sq)
end

function verify(source,out; settings=joinpath(@__DIR__,"..","config","background_conditioning.toml"))
    BLAS.set_num_threads(1)
    input=D.inputs(Dict("--input-root"=>source,"--settings"=>settings))
    options=input.settings.options; audit_options=input.settings.audit
    allfits=Dict{String,String}[]; repeatability=Dict{String,String}[]
    @testset "Exact affine problem, all endpoints, independent physics and stationarity" begin
        for item in input.cases
            case=merge(item,(;options)); ctx=D.context(case); dir=joinpath(out,"chunk$(case.index)")
            @test !isfile(joinpath(dir,"failures.tsv"))
            @test read(joinpath(dir,"settings.toml"))==read(settings)
            @test rows(joinpath(dir,"ranking.tsv"))==[Dict(k=>string(v) for (k,v) in r) for r in input.ranking]
            hashes=Dict(r["input"]=>r["sha256"] for r in rows(joinpath(dir,"input_hashes.tsv")))
            @test hashes==input.hashes
            @test read(joinpath(dir,case.file*".physical.toml"))==read(joinpath(case.root,"global_sigma_max.ell.toml"))
            @test read(joinpath(dir,case.file*".input_pixels.tsv"))==read(joinpath(case.root,"fit_data",case.file*".tsv"))
            initial=only(rows(joinpath(dir,case.file*".initial.tsv")))
            recalculated=D.initial_diagnostics(ctx,options,audit_options)
            @test Set(keys(initial))==Set(keys(recalculated))
            for k in keys(initial)
                @test k in ("stationarity_passed","stationarity_agreement") ? initial[k]==recalculated[k] : close(n(initial,k),n(recalculated,k),options)
            end
            fits=rows(joinpath(dir,"fits.tsv"))
            @test length(fits)==6 && all(r["file"]==case.file for r in fits)
            @test Set((r["arm"],i(r,"repetition")) for r in fits)==Set((arm,rep) for arm in B.ARMS for rep in 1:2)
            saved=Dict{Tuple{String,Int},Vector{Float64}}()
            for r in fits
                arm=r["arm"]; rep=i(r,"repetition"); prefix="$(case.file).repeat$rep.$arm"
                ps=sort(rows(joinpath(dir,"$prefix.parameters.tsv"));by=r->i(r,"index"))
                @test i.(ps,"index")==collect(eachindex(ctx.problem.q0))
                p=n.(ps,"value"); lo=n.(ps,"lower"); hi=n.(ps,"upper")
                @test n.(ps,"initial")==ctx.problem.q0 && lo==ctx.problem.lo && hi==ctx.problem.hi
                @test all(lo.<=p.<=hi) && i(r,"N")==case.n && r["family"]=="ell"
                @test n(r,"raw_endpoint_violation")<=options.coordinate_roundoff_tolerance
                saved[(arm,rep)]=p
                physical=independent(p,case,ctx); pred=physical.pred
                pixels=sort(rows(joinpath(dir,"$prefix.pixels.tsv"));by=r->i(r,"index"))
                @test i.(pixels,"index")==collect(eachindex(pred))
                @test all(close(a,b,options) for (a,b) in zip(pred,n.(pixels,"prediction")))
                @test all(close(a,b,options) for (a,b) in zip(ctx.data.z.-pred,n.(pixels,"residual")))
                rss=sum(abs2,pred.-ctx.data.z); pixel_count=length(pred); np=3+5case.n
                @test close(rss,n(r,"rss"),options) && close(pixel_count/(pixel_count-np)^2*rss,n(r,"gcv"),options)
                @test close(sum(abs2,independent(ctx.problem.q0,case,ctx).pred.-ctx.data.z),n(r,"initial_rss"),options)
                ls=sort(rows(joinpath(dir,"$prefix.lobes.tsv"));by=r->i(r,"lobe"))
                @test i.(ls,"lobe")==collect(1:case.n)
                for (field,key) in ((:amplitudes,"amplitude"),(:xs,"x_nm"),(:ys,"y_nm"),(:ts,"t_nm"),(:us,"u_nm"),(:sp,"sigma_parallel_nm"),(:sq,"sigma_perp_nm"))
                    @test all(close(a,b,options) for (a,b) in zip(getproperty(physical,field),n.(ls,key)))
                end
                ov=maximum(exp(-0.5*(hypot(physical.xs[a]-physical.xs[b],physical.ys[a]-physical.ys[b])/max(mean(physical.sp),mean(physical.sq)))^2) for a in 1:case.n-1 for b in a+1:case.n)
                peak=maximum(abs,pred.-ctx.data.z)/ctx.data.noise
                @test close(ov,n(r,"overlap"),options) && close(peak,n(r,"residual_peak_snr"),options)
                @test (r["valid"]=="true")==(ov<=case.cfg.max_overlap && peak<=case.cfg.residual_peak_snr_threshold && n(r,"endpoint_overrun_nm")<=1e-6)
                predict=q->independent(q,case,ctx).pred
                J=C.finite_jacobian(predict,p,lo,hi,options.derivative_relative_step)
                fit=(params=p,lower=lo,upper=hi,initial_rss=n(r,"initial_rss"),diagnostic=(;predict,target=ctx.data.z,jacobian=J))
                a=C.stationarity(fit,audit_options)
                tol=options.gradient_agreement_atol+options.gradient_agreement_rtol*max(norm(a.gradient,Inf),norm(n.(ps,"gradient_half_step"),Inf))
                @test norm(a.gradient.-n.(ps,"gradient_half_step"),Inf)<=tol
                @test norm(a.projected.-n.(ps,"projected_half_step"),Inf)<=tol
                @test a.agreement==(r["stationarity_agreement"]=="true") && a.passed==(r["stationarity_passed"]=="true")
                for mode in B.ARMS[2:3]
                    map=B.coordinates(ctx.problem,ctx.data.x,ctx.data.y,mode,options)
                    v=B.to_coordinates(p,ctx.problem,map)
                    @test all(map.lower.-1e-10 .<= v .<= map.upper.+1e-10)
                    @test maximum(map.A*v.-map.bconstraint)<=options.coordinate_roundoff_tolerance
                    @test all(close(x,y,options) for (x,y) in zip(B.from_coordinates(v,ctx.problem,map),p))
                    mode=="centered_plane_slsqp" && @test isapprox(map.plane_gram,Matrix{Float64}(I,3,3)*map.scale^2;rtol=1e-10,atol=1e-14)
                end
                trace=rows(joinpath(dir,"$prefix.trace.tsv")); checkpoints=rows(joinpath(dir,"$prefix.checkpoints.tsv"))
                @test !isempty(trace) && !isempty(checkpoints)
                if arm=="native_lm"
                    @test i.(trace,"iteration")==collect(0:i(r,"iterations"))
                    @test i(r,"iterations")<=options.lm_maxiter && all(diff(n.(trace,"rss")).<=0)
                else
                    @test i.(trace,"evaluation")==collect(1:i(r,"evaluations"))
                    @test i(r,"evaluations")<=options.slsqp_maxeval && all(diff(n.(trace,"best_rss")).<=0)
                    @test count(q->n(q,"bound_violation")>options.coordinate_roundoff_tolerance,trace)==i(r,"infeasible_trials")
                end
                for index in unique(i.(checkpoints,"checkpoint"))
                    point=sort(filter(q->i(q,"checkpoint")==index,checkpoints);by=q->i(q,"index"))
                    @test i.(point,"index")==collect(eachindex(p))
                    q=n.(point,"value")
                    expected=index==0 ? n(r,"initial_rss") : n(only(t for t in trace if i(t,arm=="native_lm" ? "iteration" : "evaluation")==index),"rss")
                    @test close(sum(abs2,predict(q).-ctx.data.z),expected,options)
                end
            end
            for arm in B.ARMS
                a,b=[only(r for r in fits if r["arm"]==arm && i(r,"repetition")==rep) for rep in 1:2]
                push!(repeatability,Dict("file"=>case.file,"arm"=>arm,
                    "max_parameter_difference"=>D.str(maximum(abs,saved[(arm,1)].-saved[(arm,2)])),
                    "relative_rss_difference"=>D.str(abs(n(a,"rss")-n(b,"rss"))/max(n(a,"rss"),n(b,"rss"))),
                    "same_validity"=>D.str(a["valid"]==b["valid"]),"same_stationarity"=>D.str(a["stationarity_passed"]==b["stationarity_passed"])))
            end
            append!(allfits,fits)
        end
        @test length(allfits)==24 && length(repeatability)==12
    end
    D.write_rows(joinpath(out,"fits.tsv"),allfits); D.write_rows(joinpath(out,"repeatability.tsv"),repeatability)
    println("Complete: 4 saved predicted counts, 3 numerical methods, 2 repeats, 24 endpoints. No grade.")
    (;fits=allfits,repeatability)
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("verify_background_conditioning.jl SAVED_COUNT_RUN OUTPUT_RUN")
    BackgroundConditioningVerification.verify(ARGS...)
end
