"Numerical diagnostic only; never loaded by production."
module BackgroundConditioning
include(joinpath(@__DIR__, "paired_convergence.jl"))
using .PairedConvergence, LinearAlgebra, Statistics, TOML
const C=PairedConvergence
const G=C.P.G
const VP=C.P.VP
const NL=G.OptimizationNLopt.NLopt
const ARMS=("native_lm","unit_box_slsqp","centered_plane_slsqp")

function settings(filename)
    raw=TOML.parsefile(filename)
    Set(keys(raw))==Set(["model","selection","preprocessing"]) || error("Unexpected settings section")
    m=raw["model"]
    expected=Set(["method","arms","repetitions","lm_maxiter","slsqp_maxeval","max_time_s","x_tol","g_tol",
        "slsqp_ftol_rel","slsqp_xtol_rel","constraint_tolerance","coordinate_roundoff_tolerance",
        "derivative_relative_step","finite_difference_half_step_factor","gradient_agreement_atol",
        "gradient_agreement_rtol","stationarity_tolerance","objective_scale_floor","jacobian_rank_rtol",
        "replay_atol","replay_rtol","checkpoints","evaluation_checkpoints","background_scaling"])
    Set(keys(m))==expected || error("Unexpected numerical setting")
    m["method"]=="exact_affine_background_conditioning" && m["arms"]==collect(ARMS) || error("Unexpected methods")
    m["background_scaling"]=="observed_target_population_standard_deviation" || error("Unexpected plane scale")
    m["repetitions"]==2 && m["x_tol"]==1e-8 && m["g_tol"]==1e-12 || error("Frozen repeats/native tolerances required")
    m["finite_difference_half_step_factor"]==0.5 || error("Half-step check required")
    for k in setdiff(expected,Set(["method","arms","checkpoints","evaluation_checkpoints","background_scaling"]))
        v=m[k]; v isa Real && !(v isa Bool) && isfinite(v) && v>0 || error("Positive finite $k required")
    end
    for k in ("lm_maxiter","slsqp_maxeval","repetitions"); m[k] isa Integer || error("Integer $k required"); end
    for (k,last_index,first_index) in (("checkpoints",m["lm_maxiter"],0),("evaluation_checkpoints",m["slsqp_maxeval"],1))
        v=m[k]; v isa Vector && !isempty(v) && all(x->x isa Integer,v) && issorted(v) &&
            length(unique(v))==length(v) && first(v)==first_index && last(v)==last_index || error("Invalid $k")
    end
    raw["selection"]==Dict("input_arm"=>"global_sigma_max","input_repetition"=>1,"cohort_size"=>146,
        "case_ranks"=>[1,49,98,146],"ranking"=>"ascending_relative_repeat_gcv_difference_then_filename",
        "count_and_family"=>"reuse_saved_selected_no_reselection",
        "start"=>"same_saved_selected_endpoint_for_all_arms_and_repetitions",
        "failure_policy"=>"retain_all_no_retry_no_grade") || error("Unexpected case selection")
    raw["preprocessing"]==Dict("input"=>"unchanged_saved_fused_fit_pixels","physical_config"=>"unchanged_saved_count_config",
        "support"=>"unchanged_saved_axis_support_noise") || error("Unexpected preprocessing")
    opts=(;(Symbol(k)=>v for (k,v) in m)...)
    audit=merge(opts,(finite_difference_relative_step=opts.derivative_relative_step,))
    return (;raw,options=opts,audit)
end

"Exact affine native-parameter map. The outer box is NOT the feasible set."
function coordinates(problem,x,y,mode,options)
    mode in ARMS[2:3] || error("Only the two declared SLSQP coordinates")
    (;q0,lo,hi,target)=problem
    length(q0)>=3 && all(lo.<hi) && all(lo.<=q0.<=hi) || error("Invalid native box")
    length(x)==length(y)==length(target)>3 || error("Insufficient shared pixels")
    M=Matrix(Diagonal(hi.-lo)); center=[mean(x),mean(y)]
    scale=sqrt(mean(abs2,target.-mean(target)))
    isfinite(scale) && scale>0 || error("Unidentifiable constant target scale; no fallback")
    X=hcat(ones(length(x)),x,y)
    Xc=hcat(ones(length(x)),x.-center[1],y.-center[2])
    singular=svdvals(Xc)
    last(singular)>options.jacobian_rank_rtol*first(singular) || error("Unidentifiable background plane; no fallback")
    if mode=="centered_plane_slsqp"
        # beta = T*p[1:3] expresses the SAME plane in centered coordinates.
        # Xc/R has orthonormal columns. Scale to a target-derived RMS height.
        T=[1. center[1] center[2]; 0. 1. 0.; 0. 0. 1.]
        R=Matrix(qr(Xc).R)
        M[1:3,1:3]=T\(R\Matrix{Float64}(I,3,3)) .* (sqrt(length(x))*scale)
    end
    inverse=inv(M)
    a=lo.-q0; b=hi.-q0
    lower=[sum(min(inverse[i,j]*a[j],inverse[i,j]*b[j]) for j in eachindex(q0)) for i in eachindex(q0)]
    upper=[sum(max(inverse[i,j]*a[j],inverse[i,j]*b[j]) for j in eachindex(q0)) for i in eachindex(q0)]
    # Six exact, normalized inequalities retain the original background box.
    # The other parameters have a diagonal map and exact scalar bounds.
    A=vcat(-M[1:3,:]./(hi.-lo)[1:3],M[1:3,:]./(hi.-lo)[1:3])
    bconstraint=vcat((q0.-lo)[1:3]./(hi.-lo)[1:3],(hi.-q0)[1:3]./(hi.-lo)[1:3])
    return (;M,inverse,lower,upper,A,bconstraint,center,scale,mode,
        plane_gram=(X*M[1:3,1:3])'*(X*M[1:3,1:3])/length(x))
end

from_coordinates(v,problem,map)=problem.q0.+map.M*v
to_coordinates(p,problem,map)=map.inverse*(p.-problem.q0)
violation(p,lo,hi)=max(0.,maximum((lo.-p)./(hi.-lo)),maximum((p.-hi)./(hi.-lo)))

"Same native derivative calculation in both SLSQP arms, plus exact plane columns."
function jacobian(predict,p,lo,hi,x,y,step)
    # SLSQP can visit an infeasible linear-plane trial. Molecular derivatives do
    # not depend on the additive plane, so evaluate them at its box projection.
    q=copy(p); q[1:3]=clamp.(q[1:3],lo[1:3],hi[1:3])
    J=C.finite_jacobian(predict,q,lo,hi,step)
    J[:,1].=1.; J[:,2].=x; J[:,3].=y
    J
end

"All methods start from exactly the same saved physical vector."
function solve(problem,x,y,arm,options,audit)
    (;q0,lo,hi,predict,target)=problem
    all(isfinite,q0) && all(lo.<=q0.<=hi) || error("Invalid shared start")
    initial_rss=sum(abs2,predict(q0).-target)
    scale=max(initial_rss,options.objective_scale_floor)
    if arm=="native_lm"
        model=(unused,p)->predict(p)
        setup_started=time_ns(); model(nothing,q0); setup_s=(time_ns()-setup_started)/1e9
        started=time_ns()
        fit=G.LsqFit.curve_fit(model,zeros(length(target)),target,copy(q0);lower=lo,upper=hi,
            maxIter=options.lm_maxiter,maxTime=options.max_time_s,x_tol=options.x_tol,g_tol=options.g_tol,
            autodiff=:finite,store_trace=true)
        elapsed_s=(time_ns()-started)/1e9
        params=copy(fit.param); iterations=last(fit.trace).iteration
        result=(;params,lower=lo,upper=hi,initial_rss,rss=sum(abs2,predict(params).-target),elapsed_s,setup_s,
            status=fit.converged ? "native_small_step_or_gradient" : iterations==options.lm_maxiter ? "iteration_cap" : "time_cap",
            converged=fit.converged,iterations,evaluations=-1,gradient_evaluations=-1,model_evaluations=-1,
            infeasible_trials=0,raw_endpoint_violation=violation(params,lo,hi),endpoint_roundoff=0.,
            diagnostic=(initial=copy(q0),predict,target,jacobian=copy(fit.jacobian),trace=deepcopy(fit.trace)))
        replay=C.replay_trace(result,audit)
        return merge(result,(history=replay.rows,checkpoints=replay.points,map=nothing))
    end
    map=coordinates(problem,x,y,arm,options)
    from_coordinates(zeros(length(q0)),problem,map)==q0 || error("Start changed")
    evaluations=Ref(0); gradients=Ref(0); model_evaluations=Ref(0); infeasible=Ref(0)
    counted=p->begin model_evaluations[]+=1; predict(p) end
    history=NamedTuple[]; checkpoints=Dict(0=>copy(q0)); best=Ref(Inf); last_p=Ref(copy(q0))
    callback=function(v,gradient)
        p=from_coordinates(v,problem,map)
        # Only floating-point roundoff in diagonal molecular bounds is clipped.
        width=hi.-lo; tol=options.coordinate_roundoff_tolerance
        all(lo[4:end].-tol.*width[4:end].<=p[4:end].<=hi[4:end].+tol.*width[4:end]) || error("Molecular trial outside bounds")
        p[4:end]=clamp.(p[4:end],lo[4:end],hi[4:end])
        err=violation(p,lo,hi); infeasible[]+=err>tol
        residual=counted(p).-target; rss=sum(abs2,residual)
        if !isempty(gradient)
            J=jacobian(counted,p,lo,hi,x,y,options.derivative_relative_step)
            gradient.=2 .* (map.M'*(J'*residual))./scale
            gradients[]+=1
        end
        evaluations[]+=1; best[]=min(best[],rss); last_p[]=copy(p)
        push!(history,(evaluation=evaluations[],rss,best_rss=best[],bound_violation=err,gradient_requested=!isempty(gradient)))
        evaluations[] in options.evaluation_checkpoints && (checkpoints[evaluations[]]=copy(p))
        rss/scale
    end
    setup_started=time_ns(); callback(zeros(length(q0)),zeros(length(q0))); setup_s=(time_ns()-setup_started)/1e9
    empty!(history); empty!(checkpoints); checkpoints[0]=copy(q0)
    evaluations[]=gradients[]=model_evaluations[]=infeasible[]=0; best[]=Inf
    opt=NL.Opt(:LD_SLSQP,length(q0)); opt.lower_bounds=map.lower; opt.upper_bounds=map.upper
    opt.ftol_rel=options.slsqp_ftol_rel; opt.xtol_rel=options.slsqp_xtol_rel
    opt.maxeval=options.slsqp_maxeval; opt.maxtime=options.max_time_s; opt.min_objective=callback
    if arm=="centered_plane_slsqp"
        for i in axes(map.A,1)
            ai=copy(map.A[i,:]); bi=map.bconstraint[i]
            NL.inequality_constraint!(opt,(v,g)->begin
                isempty(g) || (g.=ai)
                dot(ai,v)-bi
            end,options.constraint_tolerance)
        end
    end
    started=time_ns(); objective,v,status=NL.optimize(opt,zeros(length(q0))); elapsed_s=(time_ns()-started)/1e9
    status in (:SUCCESS,:STOPVAL_REACHED,:FTOL_REACHED,:XTOL_REACHED,:MAXEVAL_REACHED,:MAXTIME_REACHED,:ROUNDOFF_LIMITED) || error("SLSQP failed: $status")
    raw=from_coordinates(v,problem,map); err=violation(raw,lo,hi)
    err<=options.coordinate_roundoff_tolerance || error("Infeasible returned endpoint: $err; no fallback")
    params=clamp.(raw,lo,hi); rss=sum(abs2,predict(params).-target)
    isapprox(sum(abs2,predict(raw).-target)/scale,objective;atol=audit.replay_atol,rtol=audit.replay_rtol) || error("Returned objective differs")
    evaluations[]==NL.numevals(opt) || error("Callback count differs")
    checkpoints[evaluations[]]=last_p[]
    J=jacobian(predict,params,lo,hi,x,y,options.derivative_relative_step)
    return (;params,lower=lo,upper=hi,initial_rss,rss,elapsed_s,setup_s,status=string(status),converged=false,iterations=-1,
        evaluations=evaluations[],gradient_evaluations=gradients[],model_evaluations=model_evaluations[],
        infeasible_trials=infeasible[],raw_endpoint_violation=err,endpoint_roundoff=maximum(abs,params.-raw),
        history,checkpoints,map,diagnostic=(initial=copy(q0),predict,target,jacobian=J))
end
end
