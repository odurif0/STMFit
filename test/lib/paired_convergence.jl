module PairedConvergence
using TOML, LinearAlgebra
include(joinpath(@__DIR__,"paired_acquisition_fit.jl"))
const P=PairedAcquisitionFit

function settings(path)
    raw=TOML.parsefile(path)
    Set(keys(raw))==Set(["model","selection","preprocessing"]) || error("Unexpected convergence section")
    m=raw["model"]
    keys_expected=Set(["method","control_maxiter","extended_maxiter","max_time_s","x_tol","g_tol","checkpoints",
        "finite_difference_relative_step","finite_difference_half_step_factor","gradient_agreement_atol",
        "gradient_agreement_rtol","stationarity_tolerance","objective_scale_floor","jacobian_rank_rtol","replay_atol","replay_rtol"])
    Set(keys(m))==keys_expected || error("Unexpected convergence model key")
    m["method"]=="native_lsqfit_finite_difference" || error("Only budget extension is authorized")
    for key in setdiff(keys_expected,Set(["method","checkpoints"]))
        v=m[key]; v isa Real && !(v isa Bool) && isfinite(v) && v>0 || error("Positive finite $key required")
    end
    m["control_maxiter"] isa Integer && m["extended_maxiter"] isa Integer || error("Integer iteration budgets required")
    m["extended_maxiter"]>m["control_maxiter"] || error("Extended budget must be larger")
    m["x_tol"]==1e-8 && m["g_tol"]==1e-12 || error("Do not loosen native stopping tolerances")
    c=m["checkpoints"]
    c isa Vector && all(x->x isa Integer,c) && issorted(c) && length(unique(c))==length(c) &&
        first(c)==0 && last(c)==m["extended_maxiter"] && m["control_maxiter"] in c || error("Invalid checkpoints")
    m["finite_difference_half_step_factor"]==0.5 || error("Derivative cross-check must halve the step")
    0<m["jacobian_rank_rtol"]<1 || error("Invalid numerical rank threshold")
    raw["selection"]==Dict("file"=>"240817_006.sxm","count_policy"=>"reuse_native_diagnostic_N",
        "family_policy"=>"valid_minimum_full_parameter_gcv",
        "convergence_policy"=>"report_lm_stop_and_independent_box_stationarity_separately",
        "failure_policy"=>"explicit_no_retry_no_partial_grade") || error("Unexpected convergence selection setting")
    raw["preprocessing"]==Dict("input"=>"verified_saved_registered_fit_pixels","support"=>"reuse_native_frame_and_bounds",
        "noise"=>"reuse_native_noise") || error("Unexpected convergence preprocessing setting")
    return (; (Symbol(k)=>v for (k,v) in m)...)
end

"Second-order, bound-respecting finite differences, independent of the solver Jacobian."
function finite_jacobian(predict,q,lo,hi,relative_step)
    all(lo.<hi) && all(lo.<=q.<=hi) || error("Invalid box point")
    f0=predict(q); J=Matrix{Float64}(undef,length(f0),length(q))
    for j in eachindex(q)
        h=min(relative_step*max(1.,abs(q[j])),(hi[j]-lo[j])/4)
        a=copy(q); b=copy(q)
        if q[j]-h>=lo[j] && q[j]+h<=hi[j]
            a[j]+=h; b[j]-=h
            J[:,j]=(predict(a).-predict(b))./(2h)
        elseif q[j]+2h<=hi[j]
            a[j]+=h; b[j]+=2h
            J[:,j]=(-3 .* f0 .+ 4 .* predict(a) .- predict(b))./(2h)
        else
            a[j]-=h; b[j]-=2h
            J[:,j]=(3 .* f0 .- 4 .* predict(a) .+ predict(b))./(2h)
        end
    end
    all(isfinite,J) || error("Nonfinite independent Jacobian")
    return J
end

"Projected gradient of RSS/initial_RSS on the unit parameter box, not a validity gate."
function box_gradient(q,half_gradient,lo,hi,scale)
    width=hi.-lo; u=(q.-lo)./width
    gradient=2 .* width .* half_gradient ./ scale
    projected=u.-clamp.(u.-gradient,0.,1.)
    return (;gradient,projected,norm=norm(projected,Inf))
end

function stationarity(fit,options)
    diag=fit.diagnostic; diag===nothing && error("Diagnostics required")
    q,lo,hi=fit.params,fit.lower,fit.upper
    residual=diag.predict(q).-diag.target
    scale=max(fit.initial_rss,options.objective_scale_floor)
    J1=finite_jacobian(diag.predict,q,lo,hi,options.finite_difference_relative_step)
    J2=finite_jacobian(diag.predict,q,lo,hi,options.finite_difference_relative_step*options.finite_difference_half_step_factor)
    b1=box_gradient(q,J1'*residual,lo,hi,scale)
    b2=box_gradient(q,J2'*residual,lo,hi,scale)
    solver=box_gradient(q,diag.jacobian'*residual,lo,hi,scale)
    error=norm(b1.gradient.-b2.gradient,Inf)
    tolerance=options.gradient_agreement_atol+options.gradient_agreement_rtol*max(norm(b1.gradient,Inf),norm(b2.gradient,Inf))
    agreement=error<=tolerance
    s=svdvals(J2.*reshape(hi.-lo,1,:))
    rank=isempty(s) || first(s)==0 ? 0 : count(>(options.jacobian_rank_rtol*first(s)),s)
    condition=isempty(s) || last(s)==0 ? Inf : first(s)/last(s)
    passed=agreement && max(b1.norm,b2.norm)<=options.stationarity_tolerance
    return (;scale,gradient_error=error,gradient_tolerance=tolerance,agreement,passed,rank,condition,
        projected_gradient=b1.norm,half_step_projected_gradient=b2.norm,solver_projected_gradient=solver.norm,
        raw_gradient_norm=norm(J2'*residual,Inf),singular_values=s,
        gradient=b2.gradient,projected=b2.projected,first_gradient=b1.gradient)
end

"Replay accepted LM updates from the stored trace; rejected trial steps are not applied."
function replay_trace(fit,options)
    trace=fit.diagnostic.trace; p=copy(fit.diagnostic.initial)
    !isempty(trace) && first(trace).iteration==0 || error("Initial trace state missing")
    isapprox(first(trace).value,fit.initial_rss;atol=options.replay_atol,rtol=options.replay_rtol) || error("Trace initial RSS mismatch")
    rows=NamedTuple[]; points=Dict(0=>copy(p)); previous=first(trace).value
    for (i,state) in enumerate(trace)
        state.iteration==i-1 || error("Trace iterations not contiguous")
        state.value<=previous || error("LM RSS increased")
        accepted=i>1 && state.value<previous
        step=i==1 ? zeros(length(p)) : state.metadata["dx"]
        accepted && (p.+=step)
        all(fit.lower.<=p.<=fit.upper) || error("Replayed step outside bounds")
        push!(rows,(iteration=state.iteration,rss=state.value,reported_gradient=state.g_norm,
            trial_step_norm=norm(step),lambda=state.metadata["lambda"],accepted=accepted))
        (state.iteration in options.checkpoints || i==length(trace)) && (points[state.iteration]=copy(p))
        previous=state.value
    end
    all(isapprox.(p,fit.params;atol=options.replay_atol,rtol=options.replay_rtol)) || error("Trace did not reconstruct final parameters")
    return (;rows,points)
end

function fit_pair(initial,data,cfg,fwd,bwd,mode,model_options,options)
    kw=(maxTime=Float64(options.max_time_s),x_tol=Float64(options.x_tol),g_tol=Float64(options.g_tol))
    control=P.fit(initial,data,cfg,fwd,bwd,mode,merge(model_options,(local_maxiter=options.control_maxiter,));optimizer_options=kw,diagnostic=true)
    extended=P.fit(initial,data,cfg,fwd,bwd,mode,merge(model_options,(local_maxiter=options.extended_maxiter,));optimizer_options=kw,diagnostic=true)
    ct=replay_trace(control,options); et=replay_trace(extended,options)
    length(et.rows)>=length(ct.rows) || error("Extended trace ended before control")
    for (a,b) in zip(ct.rows,et.rows)
        for k in (:rss,:reported_gradient,:trial_step_norm,:lambda)
            x,y=getproperty(a,k),getproperty(b,k)
            (isnan(x) && isnan(y)) || isapprox(x,y;atol=options.replay_atol,rtol=options.replay_rtol) || error("Trace prefixes differ: $k")
        end
        a.accepted==b.accepted || error("Accepted trace steps differ")
    end
    if haskey(et.points,control.iterations)
        all(isapprox.(et.points[control.iterations],control.params;atol=options.replay_atol,rtol=options.replay_rtol)) || error("Common-iteration parameters differ")
    end
    return (;control,extended,control_trace=ct,extended_trace=et,
        control_stationarity=stationarity(control,options),extended_stationarity=stationarity(extended,options))
end
end
