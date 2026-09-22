module PairedSolver
include(joinpath(@__DIR__,"paired_convergence.jl"))
using .PairedConvergence, TOML, LinearAlgebra
const C=PairedConvergence
const P=C.P
const G=P.G
const NL=G.OptimizationNLopt.NLopt

function settings(path)
    t=TOML.parsefile(path)
    Set(keys(t))==Set(["model","selection","preprocessing"]) || error("Unexpected solver section")
    m=t["model"]
    Set(keys(m))==Set(["method","lm_maxiter","slsqp_maxeval","max_time_s","slsqp_ftol_rel","slsqp_xtol_rel",
        "derivative_relative_step","evaluation_checkpoints","coordinate_roundoff_tolerance"]) || error("Unexpected solver key")
    m["method"]=="unit_box_slsqp" || error("Only the declared SLSQP comparison is allowed")
    for k in ("lm_maxiter","slsqp_maxeval","max_time_s","slsqp_ftol_rel","slsqp_xtol_rel","derivative_relative_step","coordinate_roundoff_tolerance")
        v=m[k]; v isa Real && !(v isa Bool) && isfinite(v) && v>0 || error("Positive finite $k required")
    end
    m["lm_maxiter"] isa Integer && m["slsqp_maxeval"] isa Integer || error("Integer budgets required")
    c=m["evaluation_checkpoints"]
    !isempty(c) && all(x->x isa Integer && x>0,c) && issorted(c) && length(unique(c))==length(c) &&
        first(c)==1 && last(c)==m["slsqp_maxeval"] || error("Invalid evaluation checkpoints")
    t["selection"]==Dict("file"=>"240817_006.sxm","count_policy"=>"reuse_native_diagnostic_N",
        "family_policy"=>"valid_minimum_full_parameter_gcv","stationarity"=>"reuse_paired_convergence_audit_unchanged",
        "failure_policy"=>"explicit_no_retry_no_partial_grade") || error("Unexpected solver selection")
    t["preprocessing"]==Dict("input"=>"verified_saved_registered_fit_pixels","support"=>"reuse_native_frame_and_bounds",
        "noise"=>"reuse_native_noise") || error("Unexpected solver preprocessing")
    return (;(Symbol(k)=>v for (k,v) in m)...)
end

"The same native/paired least-squares problem; no optimizer or initialization."
function problem(initial,data,cfg,fwd,bwd,mode,model_options)
    mode in ("fused","paired") || error("Unknown mode")
    length(fwd)==length(bwd)==length(data.z) && all(isfinite,fwd) && all(isfinite,bwd) || error("Invalid observed samples")
    isapprox((fwd.+bwd)./2,data.z;rtol=1e-12,atol=1e-14) || error("Saved mean differs")
    initial.success && all(isfinite,initial.params) || error("Invalid native start")
    lo,hi=P.VP.native_raw_bounds(initial.n,cfg); q0=copy(initial.params)
    if mode=="paired"
        b=[model_options.gain_delta_max,model_options.background_delta_max_nm,model_options.tilt_delta_max,model_options.tilt_delta_max]
        lo=vcat(lo,-b); hi=vcat(hi,b); q0=vcat(q0,zeros(4))
    end
    all(lo.<=q0.<=hi) || error("Start outside native box")
    predict=q->mode=="fused" ? G._chain_model_values(data.x,data.y,q,initial.n,data.axisctx,cfg;
        amp_min=initial.amp_min,amp_range=initial.amp_range) : begin
        p=P.predictions(q,initial.n,data.x,data.y,data.axisctx,cfg,initial.amp_min,initial.amp_range)
        vcat(p.fwd,p.bwd)
    end
    target=mode=="fused" ? data.z : vcat(fwd,bwd)
    return (;q0,lo,hi,predict,target)
end

"Affine unit-box coordinates anchored at the exact original starting vector."
function from_unit(u,q0,lo,hi,tolerance)
    all(isfinite,u) && all(-tolerance.<=u.<=1+tolerance) || error("Solver left its unit box")
    width=hi.-lo; u0=(q0.-lo)./width
    # Clipping only absorbs affine endpoint roundoff, never an infeasible trial.
    q=q0.+width.*(u.-u0)
    all(lo.-tolerance.*width.<=q.<=hi.+tolerance.*width) || error("Affine map left the native box")
    return clamp.(q,lo,hi)
end

"SLSQP callback counts are objective evaluations, not LM iterations."
function solve(p,options,audit_options)
    (;q0,lo,hi,predict,target)=p
    scale=max(sum(abs2,predict(q0).-target),audit_options.objective_scale_floor)
    width=hi.-lo; u0=(q0.-lo)./width
    from_unit(u0,q0,lo,hi,options.coordinate_roundoff_tolerance)==q0 || error("Initial vector changed")
    evaluations=Ref(0); gradient_evaluations=Ref(0); model_evaluations=Ref(0)
    counted=q->begin model_evaluations[]+=1; predict(q) end
    history=NamedTuple[]; checkpoints=Dict(0=>copy(q0)); best=Ref(Inf); last_q=Ref(copy(q0))
    callback=function(u,gradient)
        q=from_unit(u,q0,lo,hi,options.coordinate_roundoff_tolerance)
        r=counted(q).-target; rss=sum(abs2,r); pg=NaN
        if !isempty(gradient)
            J=C.finite_jacobian(counted,q,lo,hi,options.derivative_relative_step)
            gradient.=2 .* width .* (J'*r)./scale
            gradient_evaluations[]+=1
            pg=norm(u.-clamp.(u.-gradient,0,1),Inf)
        end
        evaluations[]+=1; best[]=min(best[],rss); last_q[]=copy(q)
        push!(history,(evaluation=evaluations[],rss=rss,best_rss=best[],gradient_requested=!isempty(gradient),projected_gradient=pg))
        evaluations[] in options.evaluation_checkpoints && (checkpoints[evaluations[]]=copy(q))
        return rss/scale
    end
    # Compile the objective/Jacobian before the optimizer timer; do not take a step.
    setup_started=time_ns()
    warm_gradient=zeros(length(q0)); callback(u0,warm_gradient)
    setup_s=(time_ns()-setup_started)/1e9
    empty!(history); empty!(checkpoints); checkpoints[0]=copy(q0)
    evaluations[]=0; gradient_evaluations[]=0; model_evaluations[]=0; best[]=Inf
    opt=NL.Opt(:LD_SLSQP,length(q0)); opt.lower_bounds=zeros(length(q0)); opt.upper_bounds=ones(length(q0))
    opt.ftol_rel=options.slsqp_ftol_rel; opt.xtol_rel=options.slsqp_xtol_rel
    opt.maxeval=options.slsqp_maxeval; opt.maxtime=options.max_time_s; opt.min_objective=callback
    started=time_ns(); objective,u,status=NL.optimize(opt,u0); elapsed_s=(time_ns()-started)/1e9
    status in (:SUCCESS,:STOPVAL_REACHED,:FTOL_REACHED,:XTOL_REACHED,:MAXEVAL_REACHED,:MAXTIME_REACHED,:ROUNDOFF_LIMITED) || error("SLSQP failed: $status")
    q=from_unit(u,q0,lo,hi,options.coordinate_roundoff_tolerance)
    checkpoints[evaluations[]]=last_q[]
    rss=sum(abs2,predict(q).-target)
    isfinite(rss) && isapprox(rss/scale,objective;atol=audit_options.replay_atol,rtol=audit_options.replay_rtol) || error("Returned objective differs")
    evaluations[]==NL.numevals(opt) || error("Objective evaluation count differs")
    J=C.finite_jacobian(predict,q,lo,hi,options.derivative_relative_step)
    return (;params=q,lower=lo,upper=hi,initial_rss=sum(abs2,predict(q0).-target),rss,elapsed_s,setup_s,
        status=string(status),evaluations=evaluations[],gradient_evaluations=gradient_evaluations[],model_evaluations=model_evaluations[],
        history,checkpoints,diagnostic=(initial=copy(q0),target=copy(target),jacobian=J,predict=predict))
end

"Evaluate the returned endpoint with unchanged native validity and full GCV."
function finish(solved,initial,data,cfg,fwd,bwd,mode)
    q=solved.params; k=G._chain_nparams(initial.n,cfg)
    r=G.ChainModelResult(n=initial.n,params=copy(q[1:k]),success=true,amp_min=initial.amp_min,amp_range=initial.amp_range)
    P.VP.finalize_native!(r,data,cfg)
    both=P.predictions(mode=="paired" ? q : vcat(q,zeros(4)),r.n,data.x,data.y,data.axisctx,cfg,r.amp_min,r.amp_range)
    rss_fwd=sum(abs2,fwd.-both.fwd); rss_bwd=sum(abs2,bwd.-both.bwd)
    peak_fwd=maximum(abs,fwd.-both.fwd)/data.noise; peak_bwd=maximum(abs,bwd.-both.bwd)/data.noise
    nd=length(solved.diagnostic.target); np=length(q); gcv=nd>np ? nd/(nd-np)^2*solved.rss : Inf
    valid=r.valid && (mode=="fused" || max(peak_fwd,peak_bwd)<=cfg.residual_peak_snr_threshold)
    reason=valid ? "ok" : join(filter(!isempty,[r.valid ? "" : r.reason,
        mode=="paired" && max(peak_fwd,peak_bwd)>cfg.residual_peak_snr_threshold ? "view residual high" : ""]),",")
    return merge(solved,(;result=r,predictions=both,rss_fwd,rss_bwd,peak_fwd,peak_bwd,nd,np,gcv,valid,reason,mode))
end
end
