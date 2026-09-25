module VacuumSurfaceTransfer
using LinearAlgebra, Statistics, TOML
using STMSXMIO
include(joinpath(@__DIR__, "acquisition_registration.jl"))
const A = AcquisitionRegistration

function settings(path)
    c=TOML.parsefile(path)
    Set(keys(c))==Set(["model","selection","preprocessing"]) || error("Unexpected sections")
    fields=Dict(
        "model"=>["observable","height_reference","interpolation","surface_run","glcn_surface_sha256","glcnac_surface_sha256",
            "glcn_xml_sha256","glcnac_xml_sha256","glcn_frame_sha256","glcnac_frame_sha256","map_half_nm",
            "map_step_nm","intervals","sample_bias_v","bias_atol_v","rotation_step_deg","translations_nm",
            "height_gain","photometry","arms"],
        "selection"=>["policy","patch_radius_nm","anchor_bright_quantile","anchor_hf_mad_multiplier","anchor_plateau",
            "min_patch_pixels","bands","min_band_pixels","min_band_rows","min_row_pixels","min_correlation",
            "min_distant_peak_gap","peak_neighborhood_nm","max_band_disagreement_nm","unresolved_policy"],
        "preprocessing"=>["channel","stride","fit_flatten","anchor_flatten","anchor_smooth_radius_px",
            "plane_rank_rtol","calibration_block_rows","holdout_buffer_rows","max_lag_nm",
            "max_lag_width_fraction","missing_policy","target_correction"])
    for (section,expected) in fields
        Set(keys(c[section]))==Set(expected) || error("Unknown/missing setting in $section")
    end
    m,q,p=c["model"],c["selection"],c["preprocessing"]
    m["observable"]=="fixed_continuous_surface_local_height" || error("Wrong observable")
    m["height_reference"]=="lowest_copper_atom" && m["height_gain"]==1 || error("Changed height convention")
    m["interpolation"]=="bilinear_certified_grid_heights" || error("Wrong lateral interpolation")
    m["photometry"]=="one_common_offset_per_source_scan" || error("Wrong photometry")
    m["arms"]==["common","chemical"] && m["intervals"]==[1,2,3] || error("All arms required")
    q["policy"]=="report_all_isovalues_no_grade_no_promotion" || error("Wrong policy")
    q["anchor_plateau"]=="connected_eight_neighbor_centroid" || error("Wrong plateau policy")
    p["missing_policy"]=="observed_only_no_imputation" || error("Wrong missing policy")
    p["target_correction"]=="calibration_rows_only_difference_plane" || error("Wrong target correction")
    p["channel"]=="Z" && p["stride"]==1 && p["fit_flatten"]=="plane" &&
        p["anchor_flatten"]=="plane+rows" || error("Unexpected preprocessing")
    for (d,ks) in ((m,["map_half_nm","map_step_nm","bias_atol_v","rotation_step_deg"]),
                   (q,["patch_radius_nm","anchor_hf_mad_multiplier"]),
                   (p,["plane_rank_rtol","max_lag_nm","max_lag_width_fraction"]))
        all(k->d[k] isa Real && !(d[k] isa Bool) && isfinite(d[k]) && d[k]>0,ks) || error("Invalid positive setting")
    end
    0<q["anchor_bright_quantile"]<1 && 0<p["plane_rank_rtol"]<1 || error("Invalid quantile/tolerance")
    isinteger(360/m["rotation_step_deg"]) && isinteger(2m["map_half_nm"]/m["map_step_nm"]) || error("Invalid grid")
    d=m["translations_nm"]; !isempty(d) && issorted(d) && all(isfinite,d) && 0. in d || error("Invalid translations")
    q["patch_radius_nm"]+sqrt(2)*maximum(abs,d)<m["map_half_nm"] || error("Transform exceeds map support")
    for k in ("calibration_block_rows","anchor_smooth_radius_px","holdout_buffer_rows")
        p[k] isa Int && p[k]>=0 || error("Invalid pixel integer")
    end
    2p["holdout_buffer_rows"]<p["calibration_block_rows"] || error("Empty row blocks")
    q["min_patch_pixels"] isa Int && q["min_patch_pixels"]>=3 || error("Too few patch pixels")
    for k in ("bands","min_band_pixels","min_band_rows","min_row_pixels")
        q[k] isa Int && q[k]>=2 || error("Invalid registration integer")
    end
    0<q["min_correlation"]<=1 && 0<q["min_distant_peak_gap"]<1 &&
        0<q["peak_neighborhood_nm"]<p["max_lag_nm"] && q["max_band_disagreement_nm"]>0 &&
        p["max_lag_width_fraction"]<.5 && q["unresolved_policy"]=="zero_shift_keep_all_files" || error("Invalid registration settings")
    m["sample_bias_v"]==-.300 || error("Bias differs from the pinned surface calculation")
    for species in ("glcn","glcnac"), kind in ("surface","xml","frame")
        occursin(r"^[0-9a-f]{64}$",m[species*"_"*kind*"_sha256"]) || error("Invalid input hash")
    end
    c
end

"Equal-curvature parabolas n*(b-mu)^2+v: lower envelope in increasing b."
function envelope(mu,v,n)
    length(mu)==length(v)>0 && n>0 && all(isfinite,mu) && all(x->isfinite(x)&&x>=0,v) || error("Invalid costs")
    order=sortperm(eachindex(mu);by=i->(mu[i],v[i],i))
    ids=Int[]; starts=Float64[]
    for id in order
        !isempty(ids) && mu[id]==mu[last(ids)] && continue
        start=-Inf
        while !isempty(ids)
            j=last(ids)
            start=(mu[id]+mu[j])/2+(v[id]-v[j])/(2n*(mu[id]-mu[j]))
            start>last(starts) && break
            pop!(ids); pop!(starts)
        end
        isempty(ids) && (start=-Inf)
        push!(ids,id); push!(starts,start)
    end
    (;ids,starts)
end

"Global finite-dictionary minimum, with one common height offset, no class counts."
function common_offset(costs)
    isempty(costs) && error("No patches")
    env=[envelope(c.mu,c.v,c.n) for c in costs]
    events=sort([(e.starts[k],j,k) for (j,e) in enumerate(env) for k in 2:length(e.ids)])
    states=[first(e.ids) for e in env]; totaln=sum(c.n for c in costs)
    summean=sum(c.n*c.mu[s] for (c,s) in zip(costs,states))
    bestloss=Inf; bestb=NaN; lower=-Inf; event=1
    while true
        upper=event<=length(events) ? events[event][1] : Inf
        b=clamp(summean/totaln,lower,upper)
        loss=sum(c.v[s]+c.n*(b-c.mu[s])^2 for (c,s) in zip(costs,states))
        if loss<bestloss; bestloss=loss; bestb=b; end
        isinf(upper) && break
        while event<=length(events) && events[event][1]==upper
            _,j,k=events[event]; next=env[j].ids[k]; c=costs[j]
            summean+=c.n*(c.mu[next]-c.mu[states[j]])
            states[j]=next; event+=1
        end
        lower=upper
    end
    choices=[argmin(c.v.+c.n.*(bestb.-c.mu).^2) for c in costs]
    bestb=sum(c.n*c.mu[k] for (c,k) in zip(costs,choices))/totaln
    loss=sum(c.v[k]+c.n*(bestb-c.mu[k])^2 for (c,k) in zip(costs,choices))
    (;offset=bestb,loss,choices,segments=length(events)+1)
end

"The two sets are disjoint before looking at target heights."
function row_masks(n,c)
    p=c["preprocessing"]; block=p["calibration_block_rows"]; g=p["holdout_buffer_rows"]
    interior=[g<=mod(y-1,block)<block-g for y in 1:n]
    cal=BitVector([interior[y] && iseven(fld(y-1,block)) for y in 1:n])
    test=BitVector([interior[y] && isodd(fld(y-1,block)) for y in 1:n])
    (;cal,test)
end

function registration_settings(c)
    q,p=c["selection"],c["preprocessing"]
    s=Dict{String,Any}(k=>q[k] for k in ["bands","min_band_pixels","min_band_rows","min_row_pixels",
        "min_correlation","min_distant_peak_gap","peak_neighborhood_nm","max_band_disagreement_nm","unresolved_policy"])
    merge!(s,Dict(k=>p[k] for k in ["max_lag_nm","max_lag_width_fraction"]))
    s
end

"Only calibration rows enter registration and the target/source difference plane."
function calibrate_pair(source,target,xs,ys,c)
    size(source)==size(target)==(length(ys),length(xs)) || error("Pair grids differ")
    length(xs)>=2 && length(ys)>=2 || error("Insufficient coordinates")
    masks=row_masks(length(ys),c)
    a=copy(source); b=copy(target)
    a[.!masks.cal,:].=NaN; b[.!masks.cal,:].=NaN
    reg=A.scan(a,b,xs[2]-xs[1],registration_settings(c))
    lag=reg.applied_dx_px; inds=findall(reg.support)
    length(inds)>=3 || error("No calibration support")
    xv=[xs[i[2]] for i in inds]; yv=[ys[i[1]] for i in inds]
    dx=xv.-mean(xv); dy=yv.-mean(yv)
    xx=sum(abs2,dx); yy=sum(abs2,dy); xy=sum(dx.*dy)
    xx>0 && yy>0 && xx*yy-xy^2>c["preprocessing"]["plane_rank_rtol"]*xx*yy || error("Calibration plane unidentifiable")
    delta=[target[i[1],i[2]+lag]-source[i] for i in inds]
    X=hcat(ones(length(inds)),dx,dy); coeff=X\delta
    coeff[1]-=coeff[2]*mean(xv)+coeff[3]*mean(yv)
    all(isfinite,coeff) || error("Nonfinite calibration plane")
    (;reg,coeff,masks,calibration_pixels=length(inds))
end

madscale(v)=isempty(v) ? NaN : 1.4826*median(abs.(v.-median(v)))

"Source-only image extrema. Anchors are diagnostics, not selected molecular N."
function anchors(xs,ys,raw,z,c)
    q=c["selection"]; radius=q["patch_radius_nm"]
    good=isfinite.(z); any(good) || return NamedTuple[]
    threshold=quantile(z[good],q["anchor_bright_quantile"])
    hf=Float64[]
    for y in axes(raw,1), x in 2:size(raw,2)
        isfinite(raw[y,x]) && isfinite(raw[y,x-1]) && push!(hf,raw[y,x]-raw[y,x-1])
    end
    noise=madscale(hf)/sqrt(2)
    isfinite(noise) || return NamedTuple[]
    rx=ceil(Int,radius/(xs[2]-xs[1])); ry=ceil(Int,radius/(ys[2]-ys[1]))
    candidates=NamedTuple[]; visited=falses(size(z))
    for x in 2:length(xs)-1,y in 2:length(ys)-1
        v=z[y,x]; !visited[y,x] && isfinite(v) && v>threshold || continue
        component=[CartesianIndex(y,x)]; visited[y,x]=true; regional=true; observed=true; next=1
        while next<=length(component)
            iy,ix=Tuple(component[next]); next+=1
            for yy in iy-1:iy+1,xx in ix-1:ix+1
                yy==iy && xx==ix && continue
                if !(1<=yy<=length(ys) && 1<=xx<=length(xs)) || !isfinite(z[yy,xx])
                    observed=false; continue
                end
                z[yy,xx]>v && (regional=false)
                if z[yy,xx]==v && !visited[yy,xx]
                    visited[yy,xx]=true; push!(component,CartesianIndex(yy,xx))
                end
            end
        end
        regional && observed || continue
        cx=mean(xs[i[2]] for i in component); cy=mean(ys[i[1]] for i in component)
        center=first(sort(component;by=i->((xs[i[2]]-cx)^2+(ys[i[1]]-cy)^2,i[1],i[2])))
        iy,ix=Tuple(center)
        localvals=filter(isfinite,vec(z[max(1,iy-ry):min(end,iy+ry),max(1,ix-rx):min(end,ix+rx)]))
        prominence=v-minimum(localvals)
        prominence>0 && prominence>=q["anchor_hf_mad_multiplier"]*noise || continue
        push!(candidates,(;x=xs[ix],y=ys[iy],ix,iy,height=v,prominence,noise,threshold,plateau_size=length(component)))
    end
    sort!(candidates;by=a->(-a.height,a.iy,a.ix))
    selected=NamedTuple[]
    for a in candidates
        all(b->hypot(a.x-b.x,a.y-b.y)>2radius,selected) && push!(selected,a)
    end
    sort!(selected;by=a->(a.iy,a.ix))
end

function patches(xs,ys,z,aa,c)
    result=NamedTuple[]; radius=c["selection"]["patch_radius_nm"]
    for (j,a) in enumerate(aa)
        xr=searchsortedfirst(xs,a.x-radius):searchsortedlast(xs,a.x+radius)
        yr=searchsortedfirst(ys,a.y-radius):searchsortedlast(ys,a.y+radius)
        inds=[CartesianIndex(y,x) for x in xr for y in yr if
            hypot(xs[x]-a.x,ys[y]-a.y)<=radius && isfinite(z[y,x])]
        length(inds)>=c["selection"]["min_patch_pixels"] || continue
        dx=[xs[i[2]]-a.x for i in inds]; dy=[ys[i[1]]-a.y for i in inds]
        push!(result,(;anchor=j,inds,dx,dy,values=z[inds]))
    end
    result
end

@inline function surface_value(map,t,u,half,step)
    a=(t+half)/step+1; b=(u+half)/step+1
    1<=a<=size(map,2) && 1<=b<=size(map,1) || error("Outside map; no extrapolation")
    ix=min(floor(Int,a),size(map,2)-1); iy=min(floor(Int,b),size(map,1)-1)
    fx=a-ix; fy=b-iy
    (1-fy)*((1-fx)*map[iy,ix]+fx*map[iy,ix+1])+fy*((1-fx)*map[iy+1,ix]+fx*map[iy+1,ix+1])
end

function states(c,arm)
    m=c["model"]; types=arm=="common" ? [-1] : arm=="chemical" ? [0,1] : error("Unknown arm")
    [(;type=k,angle=Float64(a),tx=Float64(x),ty=Float64(y)) for k in types
        for a in 0:m["rotation_step_deg"]:(360-m["rotation_step_deg"])
        for y in m["translations_nm"] for x in m["translations_nm"]]
end

function prediction(p,map,state,c)
    m=c["model"]; co=cosd(state.angle); si=sind(state.angle)
    [surface_value(map,co*(x-state.tx)+si*(y-state.ty),-si*(x-state.tx)+co*(y-state.ty),
        m["map_half_nm"],m["map_step_nm"]) for (x,y) in zip(p.dx,p.dy)]
end

function patch_costs(p,maps,ss,c)
    mu=Float64[]; v=Float64[]
    for state in ss
        residual=p.values-prediction(p,maps[state.type],state,c)
        meanr=mean(residual); push!(mu,meanr); push!(v,sum(abs2,residual.-meanr))
    end
    (;n=length(p.values),mu,v)
end

function fit_dictionary(pp,maps,c,arm)
    ss=states(c,arm); costs=[patch_costs(p,maps,ss,c) for p in pp]
    fit=common_offset(costs)
    pred=[fit.offset .+ prediction(p,maps[ss[k].type],ss[k],c) for (p,k) in zip(pp,fit.choices)]
    (;ss,costs,fit,pred)
end

end
