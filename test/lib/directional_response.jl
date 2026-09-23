module DirectionalResponse
using TOML, Statistics

function settings(path)
    t=TOML.parsefile(path)
    Set(keys(t))==Set(["model","selection","preprocessing"]) || error("Unexpected response section")
    m,s,p=t["model"],t["selection"],t["preprocessing"]
    Set(keys(m))==Set(["method","mean_lag_px","kernel_tail_l1","residual_shift_max_px",
        "residual_shift_step_px","gain_min","gain_max","negative_control"]) || error("Unexpected response model key")
    Set(keys(s))==Set(["panel_quantiles","panel_policy","candidate_policy","min_pixels","min_rows",
        "min_pixels_per_row","min_heldout_gain_fraction","min_reverse_advantage_fraction",
        "lag_agreement_px","shift_agreement_px","min_supporting_scans","failure_policy"]) || error("Unexpected response selection key")
    Set(keys(p))==Set(["input","registration","support","contrast","holdout_band_rows",
        "holdout_buffer_rows","noise"]) || error("Unexpected response preprocessing key")
    for (d,k,v) in ((m,"method","mirrored_first_order_phase_transfer"),
        (m,"negative_control","reversed_scan_direction"),
        (s,"panel_policy","sorted_acquisition_filenames_without_labels"),
        (s,"candidate_policy","minimum_training_bidirectional_contrast_rss"),
        (s,"failure_policy","report_all_scans_no_retuning_no_partial_grade"),
        (p,"input","unsmoothed_native_z_actual_observation_masks"),
        (p,"registration","reuse_saved_integer_x_shifts"),
        (p,"support","common_observed_sources_for_all_candidates_both_directions"),
        (p,"contrast","within_row_centered_on_common_support"),
        (p,"noise","no_noise_calibration_or_independence_claim"))
        d[k]==v || error("Unexpected $k")
    end
    positive(x)=x isa Real && !(x isa Bool) && isfinite(x) && x>0
    lags=m["mean_lag_px"]
    length(lags)>=3 && first(lags)==0 && all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>=0 && isinteger(2x),lags) &&
        issorted(lags) && length(unique(lags))==length(lags) || error("Increasing nonnegative half-pixel lags required")
    for k in ("kernel_tail_l1","residual_shift_max_px","residual_shift_step_px","gain_min","gain_max")
        positive(m[k]) || error("Invalid $k")
    end
    0<m["kernel_tail_l1"]<1e-6 && m["gain_min"]<1<m["gain_max"] || error("Invalid kernel/gain bounds")
    isinteger(m["residual_shift_max_px"]/m["residual_shift_step_px"]) || error("Shift grid must contain zero")
    q=s["panel_quantiles"]
    length(q)>=3 && first(q)==0 && last(q)==1 && issorted(q) && length(unique(q))==length(q) &&
        all(x->x isa Real && !(x isa Bool) && isfinite(x) && 0<=x<=1,q) || error("Invalid acquisition quantiles")
    for (d,keys) in ((s,("min_pixels","min_rows","min_pixels_per_row","min_supporting_scans")),
                    (p,("holdout_band_rows","holdout_buffer_rows")))
        all(k->d[k] isa Integer && !(d[k] isa Bool) && d[k]>0,keys) || error("Positive integer settings required")
    end
    for k in ("min_heldout_gain_fraction","min_reverse_advantage_fraction","lag_agreement_px","shift_agreement_px")
        positive(s[k]) || error("Invalid $k")
    end
    s["min_heldout_gain_fraction"]<1 && s["min_reverse_advantage_fraction"]<1 || error("Invalid fractional gain")
    s["min_supporting_scans"]<=length(q) || error("Panel support exceeds panel")
    p["holdout_buffer_rows"]>=ceil(Int,m["residual_shift_max_px"])+1 &&
        2p["holdout_buffer_rows"]<p["holdout_band_rows"] || error("Unsafe spatial fold buffer")
    return (;(Symbol(k)=>v for d in (m,s,p) for (k,v) in d)...)
end

function panel(files,s)
    f=sort(collect(files))
    length(unique(f))==length(f)>=length(s.panel_quantiles) || error("Insufficient/distinct acquisition filenames")
    all(x->!isempty(x) && x==basename(x) && endswith(x,".sxm"),f) || error("Invalid scan filename")
    ix=round.(Int,1 .+(length(f)-1).*s.panel_quantiles)
    length(unique(ix))==length(ix) || error("Repeated panel index")
    return [(file=f[i],index=i,cohort_size=length(f)) for i in ix]
end

"""
Relative transfer of two mirrored, normalized first-order scan responses.
For a=mu/(1+mu), Hf=(1-a)/(1-a exp(-iw)), Hb=conj(Hf).
Hb/Hf has coefficients -a at -1, (1-a^2)a^k at k>=0 (sample offsets).
The integer 2mu recentering removes its mean advance, separating shape from
the nuisance translation. The infinite transfer has UNIT frequency magnitude:
unlike cross-smoothing, increasing mu cannot win simply by erasing contrasts.
The L1 omitted tail is explicitly bounded. No periodic boundary or padding.
"""
function kernel(mu::Real,sign::Int,tail::Real)
    isfinite(mu) && mu>=0 && isinteger(2mu) && sign in (-1,1) && 0<tail<1 || error("Invalid response kernel")
    mu==0 && return (offsets=[0],weights=[1.0],a=0.0,tail=0.0)
    a=Float64(mu/(1+mu)); kmax=max(0,ceil(Int,log(tail/(1+a))/log(a))-1)
    offsets=sign.*(vcat(-1,collect(0:kmax)).-round(Int,2mu))
    weights=vcat(-a,(1-a*a).*a.^(0:kmax))
    return (;offsets,weights,a,tail=(1+a)*a^(kmax+1))
end

function transfer(z,k)
    k.offsets==[0] && k.weights==[1.0] && return copy(z)
    ny,nx=size(z); out=fill(NaN,ny,nx)
    left=max(1,1-minimum(k.offsets)); right=min(nx,nx-maximum(k.offsets))
    left>right && return out
    out[:,left:right].=0.0
    # NaNs propagate: an unavailable source invalidates the output. There is
    # no substituted/imputed observation, including at truncated kernel edges.
    for (offset,weight) in zip(k.offsets,k.weights), x in left:right, y in 1:ny
        @inbounds out[y,x]+=weight*z[y,x+offset]
    end
    return out
end

"Conservative common rectangle of actual sources for ALL shifts/kernels/views."
function common_mask(a,b,s)
    size(a)==size(b) || throw(DimensionMismatch("View grids differ"))
    ny,nx=size(a)
    radius=maximum(maximum(abs,kernel(mu,1,s.kernel_tail_l1).offsets) for mu in s.mean_lag_px)
    rx=radius+ceil(Int,s.residual_shift_max_px)+1
    ry=ceil(Int,s.residual_shift_max_px)+1
    bad=.!isfinite.(a).|.!isfinite.(b)
    prefix=zeros(Int,ny+1,nx+1)
    for x in 1:nx,y in 1:ny
        prefix[y+1,x+1]=Int(bad[y,x])+prefix[y,x+1]+prefix[y+1,x]-prefix[y,x]
    end
    mask=falses(ny,nx); fold=zeros(Int,ny)
    for y in 1+ry:ny-ry
        t=mod(y-1,s.holdout_band_rows)
        s.holdout_buffer_rows<=t<s.holdout_band_rows-s.holdout_buffer_rows || continue
        for x in 1+rx:nx-rx
            mask[y,x]=prefix[y+ry+1,x+rx+1]-prefix[y-ry,x+rx+1]-prefix[y+ry+1,x-rx]+prefix[y-ry,x-rx]==0
        end
        if count(view(mask,y,:))<s.min_pixels_per_row
            mask[y,:].=false
        else
            fold[y]=1+mod(fld(y-1,s.holdout_band_rows),2)
        end
    end
    return (;mask,fold,rx,ry)
end

@inline function sample(z,y,x)
    yi=floor(Int,y); xi=floor(Int,x); fy=y-yi; fx=x-xi
    # Never access a zero-weight interpolation corner (not an observation).
    @inbounds a=fx==0 ? z[yi,xi] : (1-fx)*z[yi,xi]+fx*z[yi,xi+1]
    fy==0 && return a
    @inbounds b=fx==0 ? z[yi+1,xi] : (1-fx)*z[yi+1,xi]+fx*z[yi+1,xi+1]
    return (1-fy)*a+fy*b
end

"Fold sufficient statistics of observed within-row contrasts, no chemical inputs."
function moments(a,b,afrom_b,bfrom_a,dx,dy,support)
    sums=zeros(6,2); counts=zeros(Int,2); rows=zeros(Int,2)
    for y in axes(a,1)
        fold=support.fold[y]; fold==0 && continue
        n=0; sa=sb=spa=spb=aa=bb=pa2=pb2=apa=bpb=0.0
        for x in axes(a,2)
            support.mask[y,x] || continue
            @inbounds va=a[y,x]; @inbounds vb=b[y,x]
            pa=sample(afrom_b,y+dy,x+dx); pb=sample(bfrom_a,y-dy,x-dx)
            all(isfinite,(va,vb,pa,pb)) || error("Unobserved response objective source")
            n+=1; sa+=va; sb+=vb; spa+=pa; spb+=pb
            aa+=va*va; bb+=vb*vb; pa2+=pa*pa; pb2+=pb*pb; apa+=va*pa; bpb+=vb*pb
        end
        n==0 && continue
        sums[:,fold].+=(pa2-spa*spa/n,aa-sa*sa/n,apa-sa*spa/n,
                       pb2-spb*spb/n,bb-sb*sb/n,bpb-sb*spb/n)
        counts[fold]+=n; rows[fold]+=1
    end
    return (;sums,counts,rows)
end

function gains(v,s)
    v[1]>0 && v[4]>0 && v[2]>0 && v[5]>0 || return (NaN,NaN)
    return (clamp(v[3]/v[1],s.gain_min,s.gain_max),clamp(v[6]/v[4],s.gain_min,s.gain_max))
end
loss(v,g)=(max(0.,v[2]+g[1]^2*v[1]-2g[1]*v[3]),max(0.,v[5]+g[2]^2*v[4]-2g[2]*v[6]))

function scores(m,fold,s)
    train=fold==0 ? vec(sum(m.sums;dims=2)) : m.sums[:,fold]
    test=fold==0 ? train : m.sums[:,3-fold]
    g=gains(train,s)
    if !all(isfinite,g)
        return (;gain_fwd=NaN,gain_bwd=NaN,train_fwd=Inf,train_bwd=Inf,test_fwd=Inf,test_bwd=Inf,train_rss=Inf,test_rss=Inf)
    end
    tr=loss(train,g); te=loss(test,g)
    return (;gain_fwd=g[1],gain_bwd=g[2],train_fwd=tr[1],train_bwd=tr[2],test_fwd=te[1],test_bwd=te[2],train_rss=sum(tr),test_rss=sum(te))
end

relative_gain(control,candidate)=isfinite(control) && control>0 && isfinite(candidate) ? (control-candidate)/control : NaN

function selections(records,s)
    chosen=NamedTuple[]
    for fold in 0:2,mode in ("translation","physical","reversed")
        candidates=filter(r->r.fold==fold && (r.lag==0 || (mode=="physical" && r.sign==1) || (mode=="reversed" && r.sign == -1)),records)
        mode=="translation" && filter!(r->r.lag==0,candidates)
        sort!(candidates;by=r->(r.train_rss,r.lag,abs(r.dx)+abs(r.dy),r.dx,r.dy))
        isfinite(first(candidates).train_rss) || error("No finite within-row contrast candidate")
        push!(chosen,merge(first(candidates),(;mode)))
    end
    physical=filter(r->r.mode=="physical",chosen)
    checks=NamedTuple[]
    for f in 1:2
        c=only(r for r in chosen if r.fold==f && r.mode=="translation")
        p=only(r for r in chosen if r.fold==f && r.mode=="physical")
        r=only(r for r in chosen if r.fold==f && r.mode=="reversed")
        push!(checks,(fold=f,gain_fwd=relative_gain(c.test_fwd,p.test_fwd),gain_bwd=relative_gain(c.test_bwd,p.test_bwd),
            reverse_advantage=relative_gain(r.test_rss,p.test_rss)))
    end
    interior=all(r->0<r.lag<last(s.mean_lag_px) && abs(r.dx)<s.residual_shift_max_px && abs(r.dy)<s.residual_shift_max_px,physical)
    gain_interior=all(r->s.gain_min<r.gain_fwd<s.gain_max && s.gain_min<r.gain_bwd<s.gain_max,physical)
    lag_spread=maximum(r.lag for r in physical)-minimum(r.lag for r in physical)
    shift_spread=maximum(maximum(getproperty(r,key) for r in physical)-minimum(getproperty(r,key) for r in physical) for key in (:dx,:dy))
    predictive=all(r->isfinite(r.gain_fwd) && isfinite(r.gain_bwd) && min(r.gain_fwd,r.gain_bwd)>=s.min_heldout_gain_fraction,checks)
    directional=all(r->isfinite(r.reverse_advantage) && r.reverse_advantage>=s.min_reverse_advantage_fraction,checks)
    supported=interior && gain_interior && predictive && directional && lag_spread<=s.lag_agreement_px && shift_spread<=s.shift_agreement_px
    return (;chosen,checks,summary=(;supported,interior,gain_interior,predictive,directional,lag_spread,shift_spread))
end

function scan(a,b,s; progress=(_...)->nothing)
    support=common_mask(a,b,s)
    for f in 1:2
        rows=findall(==(f),support.fold)
        length(rows)>=s.min_rows && sum(support.mask[rows,:])>=s.min_pixels || error("Insufficient common observed fold $f")
    end
    shift=collect(-s.residual_shift_max_px:s.residual_shift_step_px:s.residual_shift_max_px)
    records=NamedTuple[]
    for lag in s.mean_lag_px,sign in (lag==0 ? (1,) : (1,-1))
        bf=transfer(a,kernel(lag,sign,s.kernel_tail_l1))
        af=transfer(b,kernel(lag,-sign,s.kernel_tail_l1))
        for dy in shift,dx in shift
            m=moments(a,b,af,bf,dx,dy,support)
            for fold in 0:2
                push!(records,merge((;lag=Float64(lag),sign,dx,dy,fold),scores(m,fold,s)))
            end
        end
        progress(lag,sign)
    end
    return merge((;records,support),selections(records,s))
end
end
