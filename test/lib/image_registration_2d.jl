module ImageRegistration2D

using TOML, Statistics

function settings(path)
    t=TOML.parsefile(path)
    fields=Dict("model"=>["method","reference_frame","interpolation","refinement_step_px","refinement_radius_px"],
        "selection"=>["file","bands","min_band_pixels","min_band_rows","min_row_pixels","min_correlation",
            "min_distant_peak_gap","peak_neighborhood_nm","max_band_disagreement_nm","max_fold_disagreement_px","unresolved_policy"],
        "preprocessing"=>["residual_window_nm","max_lag_width_fraction","missing_policy","support","holdout_block_px","holdout_buffer_px"])
    Set(keys(t))==Set(keys(fields)) || error("Unexpected registration section")
    for (section,keys_) in fields
        Set(keys(t[section]))==Set(keys_) || error("Unknown/missing registration field")
    end
    s=merge(values(t)...)
    for (k,v) in (("method","signed_row_centered_2d_translation"),("reference_frame","forward"),
        ("interpolation","bilinear_observed_corners"),("unresolved_policy","report_no_molecular_fit_no_retry"),
        ("missing_policy","restore_raw_mask_before_smoothing"),("support","fixed_for_all_displacements_and_interpolation_corners"))
        s[k]==v || error("Unsupported $k")
    end
    s["file"]==basename(s["file"]) && endswith(s["file"],".sxm") || error("One SXM basename required")
    for k in ("bands","min_band_pixels","min_band_rows","min_row_pixels","refinement_radius_px","holdout_block_px")
        s[k] isa Integer && !(s[k] isa Bool) && s[k]>=1 || error("Invalid $k")
    end
    for k in ("refinement_step_px","min_correlation","min_distant_peak_gap","peak_neighborhood_nm",
        "max_band_disagreement_nm","max_fold_disagreement_px","residual_window_nm","max_lag_width_fraction")
        v=s[k]; v isa Real && !(v isa Bool) && isfinite(v) && v>0 || error("Invalid $k")
    end
    s["min_correlation"]<=1 && s["min_distant_peak_gap"]<1 || error("Invalid correlation threshold")
    s["max_lag_width_fraction"]<.5 && s["residual_window_nm"]>s["peak_neighborhood_nm"] || error("Invalid search window")
    h=s["refinement_step_px"]; h<=1 && isinteger(1/h) || error("Refinement must divide a pixel")
    b,g=s["holdout_block_px"],s["holdout_buffer_px"]
    g isa Integer && !(g isa Bool) && 0<=g<b/2 || error("Invalid block buffer")
    return (;(Symbol(k)=>v for (k,v) in s)...)
end

"Both views use this parity in the already integer-registered frame."
function training_mask(shape,fold,block,buffer)
    fold in 0:2 || error("Unknown fold")
    fold==0 && return trues(shape)
    return BitMatrix([mod(fld(y-1,block)+fld(x-1,block),2)!=fold-1 &&
        buffer<=mod(y-1,block)<block-buffer && buffer<=mod(x-1,block)<block-buffer
        for y in 1:shape[1],x in 1:shape[2]])
end

"Fixed support: every candidate backward source, including fractional corners, is observed and allowed."
function common_support(a,b,base,limits,min_row_pixels; allowed=trues(size(a)))
    size(a)==size(b)==size(allowed) || throw(DimensionMismatch("Image/mask shapes differ"))
    all(>=(0),limits) || error("Negative search limits")
    ny,nx=size(a); lx,ly=limits; bx,by=base
    bad=ones(Int,ny,nx)
    for x in 1:nx,y in 1:ny
        yy,xx=y-by,x-bx
        bad[y,x]=!(isfinite(b[y,x]) && 1<=yy<=ny && 1<=xx<=nx && allowed[yy,xx])
    end
    prefix=zeros(Int,ny+1,nx+1)
    prefix[2:end,2:end].=cumsum(cumsum(bad;dims=1);dims=2)
    mask=falses(ny,nx)
    for x in 1:nx,y in 1:ny
        isfinite(a[y,x]) && allowed[y,x] || continue
        y1,y2=y+by-ly,y+by+ly; x1,x2=x+bx-lx,x+bx+lx
        1<=y1<=y2<=ny && 1<=x1<=x2<=nx || continue
        mask[y,x]=prefix[y2+1,x2+1]-prefix[y1,x2+1]-prefix[y2+1,x1]+prefix[y1,x1]==0
    end
    for y in 1:ny
        count(@view mask[y,:])>=min_row_pixels || (mask[y,:].=false)
    end
    return mask
end

"No wrap/padding; integer positions do not touch zero-weight neighboring pixels."
@inline function sample(b,y,x)
    iy,ix=floor(Int,y),floor(Int,x); fy,fx=y-iy,x-ix
    v0=fx==0 ? b[iy,ix] : (1-fx)*b[iy,ix]+fx*b[iy,ix+1]
    fy==0 && return v0
    v1=fx==0 ? b[iy+1,ix] : (1-fx)*b[iy+1,ix]+fx*b[iy+1,ix+1]
    return (1-fy)*v0+fy*v1
end

function row_data(a,mask,bands)
    out=NamedTuple[]; ny=size(a,1)
    for y in axes(a,1)
        xs=findall(@view mask[y,:]); isempty(xs) && continue
        av=a[y,xs]; av.-=mean(av)
        push!(out,(;y,xs,av,aa=sum(abs2,av),band=min(bands,1+(y-1)*bands÷ny)))
    end
    return out
end

"Signed, within-row centered NCC; offsets/gains are not chemical features."
function correlations(b,rows,base,shift,bands)
    aa=zeros(bands); bb=zeros(bands); ab=zeros(bands)
    dx,dy=base[1]+shift[1],base[2]+shift[2]
    for r in rows
        n=length(r.xs); sb=0.; sb2=0.; sab=0.
        # Stable two-pass centering avoids cancellation on large row offsets.
        for x in r.xs; sb+=sample(b,r.y+dy,x+dx); end
        mb=sb/n
        for (j,x) in enumerate(r.xs)
            v=sample(b,r.y+dy,x+dx)-mb; sb2+=v*v; sab+=r.av[j]*v
        end
        aa[r.band]+=r.aa; bb[r.band]+=sb2; ab[r.band]+=sab
    end
    result=fill(NaN,bands+1)
    for i in 0:bands
        va=i==0 ? sum(aa) : aa[i]; vb=i==0 ? sum(bb) : bb[i]; cab=i==0 ? sum(ab) : ab[i]
        va>0 && vb>0 && (result[i+1]=clamp(cab/sqrt(va*vb),-1.,1.))
    end
    return result
end

best_index(shifts,values)=begin
    ix=findall(isfinite,values)
    isempty(ix) ? nothing : first(sort(ix;by=j->(-values[j],sum(abs2,shifts[j]),shifts[j]...)))
end

function peak(shifts,values,coarse,coarsevalues,rows,region,step,limits,s)
    i=best_index(shifts,values); xy=i===nothing ? (0.,0.) : shifts[i]
    score=i===nothing ? NaN : values[i]
    far=[j for j in eachindex(coarse) if isfinite(coarsevalues[j]) &&
        maximum(abs.((coarse[j].-xy).*step))>s.peak_neighborhood_nm]
    distant=isempty(far) ? nothing : first(sort(far;by=j->(-coarsevalues[j],sum(abs2,coarse[j]),coarse[j]...)))
    gap=distant===nothing ? NaN : score-coarsevalues[distant]
    rr=region==0 ? rows : filter(r->r.band==region,rows)
    np=sum(length(r.xs) for r in rr;init=0); nr=length(rr)
    boundary=any(abs(xy[j])>=limits[j] for j in 1:2)
    reasons=String[]
    np>=s.min_band_pixels && nr>=s.min_band_rows || push!(reasons,"insufficient_support")
    isfinite(score) || push!(reasons,"constant_signal")
    isfinite(score) && score<s.min_correlation && push!(reasons,"weak_or_negative_correlation")
    isfinite(gap) && gap>=s.min_distant_peak_gap || push!(reasons,"ambiguous_peak")
    boundary && push!(reasons,"boundary_peak")
    return (dx=Float64(xy[1]),dy=Float64(xy[2]),correlation=score,distant_peak_gap=gap,
        distant_dx=distant===nothing ? NaN : Float64(coarse[distant][1]),
        distant_dy=distant===nothing ? NaN : Float64(coarse[distant][2]),
        pixels=np,rows=nr,boundary,usable=isempty(reasons),status=isempty(reasons) ? "identified" : join(reasons,";"))
end

function scan(a,b,step,base,s;fold=0)
    size(a)==size(b) || throw(DimensionMismatch("Images differ"))
    all(v->isfinite(v)&&v>0,step) || error("Invalid pixel scale")
    limits=Tuple(min(floor(Int,s.residual_window_nm/step[j]),floor(Int,s.max_lag_width_fraction*size(a,3-j))) for j in 1:2)
    allowed=training_mask(size(a),fold,s.holdout_block_px,s.holdout_buffer_px)
    mask=common_support(a,b,base,limits,s.min_row_pixels;allowed)
    rows=row_data(a,mask,s.bands)
    coarse=[(Float64(x),Float64(y)) for y in -limits[2]:limits[2] for x in -limits[1]:limits[1]]
    scores=hcat([correlations(b,rows,base,xy,s.bands) for xy in coarse]...)
    fine=NamedTuple[]; peaks=NamedTuple[]
    for region in 0:s.bands
        cv=vec(scores[region+1,:]); i=best_index(coarse,cv)
        center=i===nothing ? (0.,0.) : coarse[i]
        ranges=[range(max(-limits[j],center[j]-s.refinement_radius_px),
            min(limits[j],center[j]+s.refinement_radius_px);step=s.refinement_step_px) for j in 1:2]
        xy=[(x,y) for y in ranges[2] for x in ranges[1]]
        rr=region==0 ? rows : filter(r->r.band==region,rows)
        v=[correlations(b,rr,base,q,s.bands)[region+1] for q in xy]
        push!(fine,(;region,shifts=xy,values=v))
        push!(peaks,peak(xy,v,coarse,cv,rows,region,step,limits,s))
    end
    g=first(peaks); reasons=String[]
    all(limits[j]*step[j]>s.peak_neighborhood_nm+step[j] for j in 1:2) || push!(reasons,"insufficient_lag_span")
    for (i,p) in enumerate(peaks)
        p.usable || push!(reasons,"region$(i-1):"*p.status)
        maximum(abs.((p.dx-g.dx,p.dy-g.dy).*step))<=s.max_band_disagreement_nm || push!(reasons,"region$(i-1)_disagreement")
    end
    return (;fold,limits,mask,coarse,scores,fine,peaks,accepted=isempty(reasons),status=isempty(reasons) ? "identified" : join(reasons,";"))
end

function eligibility(runs,s)
    length(runs)==3 && getproperty.(runs,:fold)==[0,1,2] || error("Full image and both folds required")
    shifts=[(first(r.peaks).dx,first(r.peaks).dy) for r in runs]
    spread=maximum(maximum(abs.(a.-b)) for a in shifts,b in shifts)
    return (eligible=all(r.accepted for r in runs) && spread<=s.max_fold_disagreement_px,spread_px=spread)
end
end
