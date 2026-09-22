module PairedAcquisitionFit
using GaussianFit2D, LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__,"counting_variable_projection.jl"))
const G=GaussianFit2D
const VP=CountingVariableProjection

function settings(path)
    t=TOML.parsefile(path)
    Set(keys(t))==Set(["model","selection","preprocessing"]) || error("Unexpected paired configuration section")
    m=t["model"]
    Set(keys(m))==Set(["method","profiles","gain_delta_max","background_delta_max_nm","tilt_delta_max","local_maxiter"]) || error("Unexpected paired model setting")
    m["method"]=="shared_shape_antisymmetric_gain_plane" && m["profiles"]==["gaussian","split"] || error("Unsupported paired model")
    for k in ("gain_delta_max","background_delta_max_nm","tilt_delta_max","local_maxiter")
        m[k] isa Real && !(m[k] isa Bool) && isfinite(m[k]) && m[k]>0 || error("Positive finite $k required")
    end
    m["gain_delta_max"]<1 || error("Both gains must stay strictly positive")
    m["local_maxiter"] isa Integer || error("Integer iteration budget required")
    t["selection"]==Dict("count_policy"=>"reuse_saved_N","family_policy"=>"valid_minimum_full_parameter_gcv",
        "validity"=>"native_mean_and_both_view_residuals","failure_policy"=>"report_and_stop_no_partial_grade") || error("Unexpected paired selection setting")
    t["preprocessing"]==Dict("support"=>"original_native_intersect_both_observed",
        "registration"=>"reuse_accepted_integer_shifts","reference_frame"=>"forward") || error("Unexpected paired preprocessing setting")
    return (gain_delta_max=Float64(m["gain_delta_max"]),background_delta_max_nm=Float64(m["background_delta_max_nm"]),
        tilt_delta_max=Float64(m["tilt_delta_max"]),local_maxiter=Int(m["local_maxiter"]))
end

"Identifiable mean model: f=M+g*L+d; b=M-g*L-d, where L excludes the mean plane."
function predictions(q,n,x,y,axis,cfg,amp_min,amp_range)
    k=G._chain_nparams(n,cfg)
    length(q)==k+4 || throw(DimensionMismatch("Four acquisition coefficients required"))
    cfg.chain_tilted_baseline || error("Tilted native mean required")
    p=view(q,1:k)
    mean=G._chain_model_values(x,y,p,n,axis,cfg;amp_min,amp_range)
    lobes=mean.-(p[1].+p[2].*x.+p[3].*y)
    delta=q[k+1].*lobes.+q[k+2].+q[k+3].*x.+q[k+4].*y
    return (mean=mean,fwd=mean.+delta,bwd=mean.-delta,lobes=lobes)
end

"Matched LM continuation from the same native candidate; never a saved-result fallback."
function fit(initial,data,cfg,fwd,bwd,mode,options)
    mode in ("fused","paired") || error("Unknown acquisition fit mode")
    length(fwd)==length(bwd)==length(data.z) || throw(DimensionMismatch("Unequal supports"))
    all(isfinite,fwd) && all(isfinite,bwd) || error("Unobserved objective pixel")
    initial.success && all(isfinite,initial.params) || error("No numerical initial candidate")
    isapprox((fwd.+bwd)./2,data.z;rtol=1e-12,atol=1e-14) || error("Views and native mean differ")
    lo,hi=VP.native_raw_bounds(initial.n,cfg)
    p0=copy(initial.params)
    if mode=="paired"
        p0=vcat(p0,zeros(4))
        bounds=[options.gain_delta_max,options.background_delta_max_nm,options.tilt_delta_max,options.tilt_delta_max]
        lo=vcat(lo,-bounds); hi=vcat(hi,bounds)
    end
    all(lo.<=p0.<=hi) || error("Initial parameters outside bounds")
    xy=vcat(reshape(data.x,1,:),reshape(data.y,1,:))
    target=mode=="fused" ? data.z : vcat(fwd,bwd)
    model=function(xy,q)
        if mode=="fused"
            return G._chain_model_values(view(xy,1,:),view(xy,2,:),q,initial.n,data.axisctx,cfg;
                amp_min=initial.amp_min,amp_range=initial.amp_range)
        end
        pred=predictions(q,initial.n,view(xy,1,:),view(xy,2,:),data.axisctx,cfg,initial.amp_min,initial.amp_range)
        return vcat(pred.fwd,pred.bwd)
    end
    # Compile the objective before measuring numerical fitting time.
    initial_rss=sum(abs2,target.-model(xy,p0))
    started=time_ns()
    fitted=G.LsqFit.curve_fit(model,xy,target,p0;lower=lo,upper=hi,maxIter=options.local_maxiter,
        autodiff=:finite,store_trace=true)
    elapsed_s=(time_ns()-started)/1e9
    q=copy(fitted.param); k=G._chain_nparams(initial.n,cfg)
    p=copy(q[1:k])
    r=G.ChainModelResult(n=initial.n,params=p,success=true,amp_min=initial.amp_min,amp_range=initial.amp_range)
    VP.finalize_native!(r,data,cfg)
    both=mode=="paired" ? predictions(q,r.n,data.x,data.y,data.axisctx,cfg,r.amp_min,r.amp_range) :
        predictions(vcat(q,zeros(4)),r.n,data.x,data.y,data.axisctx,cfg,r.amp_min,r.amp_range)
    rss_fwd=sum(abs2,fwd.-both.fwd); rss_bwd=sum(abs2,bwd.-both.bwd)
    peak_fwd=maximum(abs,fwd.-both.fwd)/data.noise
    peak_bwd=maximum(abs,bwd.-both.bwd)/data.noise
    np=length(q); nd=length(target); rss=sum(abs2,target.-model(xy,q))
    gcv=nd>np ? nd/(nd-np)^2*rss : Inf
    # These are the existing maximum-residual guard and physical checks, not
    # a likelihood calibration or a claim of independent acquisition noise.
    valid=r.valid && (mode=="fused" || max(peak_fwd,peak_bwd)<=cfg.residual_peak_snr_threshold)
    reason=valid ? "ok" : join(filter(!isempty,[r.valid ? "" : r.reason,
        mode=="paired" && max(peak_fwd,peak_bwd)>cfg.residual_peak_snr_threshold ? "view residual high" : ""]),",")
    return (;result=r,params=q,lower=lo,upper=hi,valid,reason,rss,gcv,nd,np,rss_fwd,rss_bwd,peak_fwd,peak_bwd,
        initial_rss,elapsed_s,converged=fitted.converged,iterations=isempty(fitted.trace) ? 0 : last(fitted.trace).iteration,
        predictions=both,mode)
end
end
