module PairedShift
include(joinpath(@__DIR__,"paired_solver.jl"))
using .PairedSolver, TOML, LinearAlgebra, Random
const S=PairedSolver
const C=S.C
const P=S.P
const G=S.G

function settings(path)
    t=TOML.parsefile(path)
    Set(keys(t))==Set(["model","selection","preprocessing"]) || error("Unexpected shift section")
    m,s,p=t["model"],t["selection"],t["preprocessing"]
    Set(keys(m))==Set(["method","start_seeds","start_unit_box_radius","shift_max_px","slsqp_maxeval","max_time_s"]) || error("Unexpected shift model key")
    m["method"]=="shared_shape_backward_translation" || error("Unknown shift model")
    seeds=m["start_seeds"]
    length(seeds)>=2 && first(seeds)==0 && length(unique(seeds))==length(seeds) &&
        all(v->v isa Integer && !(v isa Bool) && v>=0,seeds) || error("Distinct deterministic starts required")
    for k in ("start_unit_box_radius","shift_max_px","slsqp_maxeval","max_time_s")
        v=m[k]; v isa Real && !(v isa Bool) && isfinite(v) && v>0 || error("Invalid $k")
    end
    m["start_unit_box_radius"]<0.5 && m["slsqp_maxeval"] isa Integer || error("Invalid start radius/budget")
    Set(keys(s))==Set(["file","count_policy","family_policy","start_policy","holdout_policy","shift_interior_margin_px","fold_shift_agreement_px","failure_policy"]) || error("Unexpected shift selection key")
    for (k,v) in (("file","240817_006.sxm"),("count_policy","reuse_native_diagnostic_N"),
        ("family_policy","valid_minimum_full_parameter_gcv"),("start_policy","minimum_training_rss_then_existing_validity"),
        ("holdout_policy","both_folds_improve_paired_rss"),("failure_policy","explicit_no_retry_no_partial_grade"))
        s[k]==v || error("Unexpected $k")
    end
    for k in ("shift_interior_margin_px","fold_shift_agreement_px")
        v=s[k]; v isa Real && !(v isa Bool) && isfinite(v) && 0<v<m["shift_max_px"] || error("Invalid $k")
    end
    Set(keys(p))==Set(["input","support","reference_frame","background_coordinates","noise","holdout_block_px","holdout_buffer_px"]) || error("Unexpected shift preprocessing key")
    for (k,v) in (("input","verified_saved_registered_fit_pixels"),("support","reuse_all_observed_pixels_without_resampling"),
        ("reference_frame","forward"),("background_coordinates","observed_grid_in_both_views"),("noise","reuse_native_noise"))
        p[k]==v || error("Unexpected $k")
    end
    b,g=p["holdout_block_px"],p["holdout_buffer_px"]
    b isa Integer && g isa Integer && !(b isa Bool) && !(g isa Bool) && b>=4 && 0<=g<b/2 || error("Invalid spatial blocks")
    return (;(Symbol(k)=>v for d in (m,s,p) for (k,v) in d)...)
end

"One backward translation of the molecular signal; planes stay on observed coordinates."
function predictions(q,n,x,y,axis,cfg,amin,arange,step)
    k=G._chain_nparams(n,cfg)
    length(q)==k+6 || throw(DimensionMismatch("Four view coefficients and two shifts required"))
    base=P.predictions(view(q,1:k+4),n,x,y,axis,cfg,amin,arange)
    dx,dy=q[k+5]*step[1],q[k+6]*step[2]
    if dx==0 && dy==0
        return base # Exact nesting, including floating-point arithmetic.
    end
    xx=x.-dx; yy=y.-dy; p=view(q,1:k)
    shifted=G._chain_model_values(xx,yy,p,n,axis,cfg;amp_min=amin,amp_range=arange)
    lobes=shifted.-(p[1].+p[2].*xx.+p[3].*yy)
    plane=(p[1]-q[k+2]).+(p[2]-q[k+3]).*x.+(p[3]-q[k+4]).*y
    backward=(1-q[k+1]).*lobes.+plane
    return (mean=(base.fwd.+backward)./2,fwd=base.fwd,bwd=backward,lobes=base.lobes)
end

function problem(ctx,mode,model,options,step; indices=eachindex(ctx.data.z))
    mode in ("fused","paired","shifted") || error("Unknown fit mode")
    base=S.problem(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,mode=="fused" ? "fused" : "paired",model)
    q0,lo,hi=copy(base.q0),copy(base.lo),copy(base.hi)
    if mode=="shifted"
        append!(q0,zeros(2)); append!(lo,fill(-options.shift_max_px,2)); append!(hi,fill(options.shift_max_px,2))
    end
    ix=collect(indices); length(unique(ix))==length(ix)>0 || error("Invalid objective support")
    x,y=ctx.data.x[ix],ctx.data.y[ix]; k=length(ctx.initial.params)
    allpred=q-> mode=="shifted" ? predictions(q,ctx.initial.n,ctx.data.x,ctx.data.y,ctx.data.axisctx,ctx.cfg,
        ctx.initial.amp_min,ctx.initial.amp_range,step) : P.predictions(mode=="fused" ? vcat(q,zeros(4)) : q,
        ctx.initial.n,ctx.data.x,ctx.data.y,ctx.data.axisctx,ctx.cfg,ctx.initial.amp_min,ctx.initial.amp_range)
    predict=q->begin
        if mode=="fused"
            G._chain_model_values(x,y,q,ctx.initial.n,ctx.data.axisctx,ctx.cfg;amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
        else
            a=mode=="shifted" ? predictions(q,ctx.initial.n,x,y,ctx.data.axisctx,ctx.cfg,ctx.initial.amp_min,ctx.initial.amp_range,step) :
                P.predictions(q,ctx.initial.n,x,y,ctx.data.axisctx,ctx.cfg,ctx.initial.amp_min,ctx.initial.amp_range)
            vcat(a.fwd,a.bwd)
        end
    end
    target=mode=="fused" ? ctx.data.z[ix] : vcat(ctx.fwd[ix],ctx.bwd[ix])
    return (;q0,lo,hi,predict,target,allpred,indices=ix,molecular_count=k)
end

"Only raw molecular coordinates change: no background/gain/shift seed or data-driven retry."
function start(p,seed,radius)
    seed==0 && return copy(p.q0)
    q=copy(p.q0); rng=Xoshiro(seed)
    for j in 4:p.molecular_count
        u=(q[j]-p.lo[j])/(p.hi[j]-p.lo[j])
        u=clamp(u+radius*(2rand(rng)-1),0.,1.)
        q[j]=clamp(p.lo[j]+u*(p.hi[j]-p.lo[j]),p.lo[j],p.hi[j])
    end
    return q
end

"Full observed-grid mean and both view guards, with unchanged physical geometry checks."
function finish(solved,p,ctx,mode)
    q=solved.params; pred=p.allpred(q); n=ctx.initial.n; k=length(ctx.initial.params); d=ctx.data
    r=G.ChainModelResult(n=n,params=copy(q[1:k]),success=true,amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
    # The mean prediction now includes the translation. No pseudo-data are created.
    G._finalize_chain_result!(r,d.z,pred.mean,d.noise,n,max(10,length(d.z)÷9),d.zfull,d.xs,d.ys,d.zimg,d.x,d.y,d.axisctx,ctx.cfg)
    rf=ctx.fwd.-pred.fwd; rb=ctx.bwd.-pred.bwd
    peak_fwd=maximum(abs,rf)/d.noise; peak_bwd=maximum(abs,rb)/d.noise
    valid=r.valid && (mode=="fused" || max(peak_fwd,peak_bwd)<=ctx.cfg.residual_peak_snr_threshold)
    nd=length(p.target); np=length(q); gcv=nd>np ? nd/(nd-np)^2*solved.rss : Inf
    # Native mean IC/GCV do not count acquisition terms and are not exported or used.
    return merge(solved,(;result=r,predictions=pred,valid,nd,np,gcv,peak_fwd,peak_bwd,rss_fwd=sum(abs2,rf),rss_bwd=sum(abs2,rb)))
end

"Reserve whole 16x16 pixel blocks; a training buffer separates neighboring test blocks."
function folds(rows,cols,block,buffer)
    length(rows)==length(cols) && !isempty(rows) || error("Invalid pixel keys")
    parity=[mod(fld(r-1,block)+fld(c-1,block),2) for (r,c) in zip(rows,cols)]
    interior=[buffer<=mod(r-1,block)<block-buffer && buffer<=mod(c-1,block)<block-buffer for (r,c) in zip(rows,cols)]
    return [(train=findall((parity.!=f).&interior),test=findall(parity.==f)) for f in 0:1]
end

function grid_step(pixels)
    out=Float64[]
    for (coord,key) in (("x_nm","column"),("y_nm","row"))
        v=parse.(Float64,getindex.(pixels,coord)); i=parse.(Int,getindex.(pixels,key))
        a,b=argmin(i),argmax(i); i[b]>i[a] || error("Degenerate image grid")
        step=(v[b]-v[a])/(i[b]-i[a]); isfinite(step) && step>0 || error("Invalid pixel step")
        maximum(abs.(v.-(v[a].+(i.-i[a]).*step)))<=1e-10 || error("Saved coordinates do not match grid")
        push!(out,step)
    end
    return Tuple(out)
end
end
