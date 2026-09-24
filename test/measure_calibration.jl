#!/usr/bin/env julia
# Observability audit, not an automatic production-calibration writer.
# Legacy placeholders are reported as placeholders; no expected N is accepted.
module CalibrationMeasurements
using STMSXMIO, GaussianFit2D
using Statistics, LinearAlgebra, TOML, SHA, Printf

const ROOT = dirname(@__DIR__)
const SETTINGS = joinpath(ROOT,"config/calibration_measurements.toml")
const METHODS = ("legacy_bootstrap", "full_profile_prominence")
const DIRECTIONS = ("fwd", "bwd")
const SCHEMA = Dict(
    "model"=>["axis_bright_quantile","strip_halfwidth_pixels","bin_width_pixels",
        "axis_degeneracy_rtol","width_quantiles","legacy_weight_floor",
        "legacy_fwhm_fallback_nm","legacy_spacing_fallback_nm"],
    "selection"=>["prominence_hf_mad_multiplier","min_width_pixels",
        "legacy_height_std_multiplier","legacy_width_window_bins","legacy_max_width_nm"],
    "preprocessing"=>["channel","stride","flatten","smooth_radius_px","plane_rank_rtol"])

function settings(path)
    c=TOML.parsefile(path)
    Set(keys(c))==Set(keys(SCHEMA)) || error("Measurement settings require only model, selection, preprocessing")
    for (section,keys_expected) in SCHEMA
        Set(keys(c[section]))==Set(keys_expected) || error("Unknown/missing setting in $section")
    end
    m=c["model"]; s=c["selection"]; p=c["preprocessing"]
    for (d,ks) in ((m,["strip_halfwidth_pixels","bin_width_pixels","axis_degeneracy_rtol",
                       "legacy_weight_floor","legacy_spacing_fallback_nm"]),
                  (s,["min_width_pixels","legacy_max_width_nm"]))
        all(k->d[k] isa Real && isfinite(d[k]) && d[k]>0,ks) || error("Positive finite settings required")
    end
    0<m["axis_bright_quantile"]<1 || error("Invalid axis quantile")
    0<m["axis_degeneracy_rtol"]<1 || error("Invalid numerical tolerance")
    for key in ("width_quantiles","legacy_fwhm_fallback_nm")
        v=m[key]; length(v)==2 && all(isfinite,v) && v[1]<v[2] || error("Invalid $key")
    end
    0<=m["width_quantiles"][1]<m["width_quantiles"][2]<=1 || error("Invalid width quantiles")
    m["legacy_fwhm_fallback_nm"][1]>0 || error("Invalid legacy widths")
    all(k->s[k] isa Real && isfinite(s[k]) && s[k]>=0,
        ["prominence_hf_mad_multiplier","legacy_height_std_multiplier"]) || error("Invalid contrast filter")
    s["legacy_width_window_bins"] isa Integer && s["legacy_width_window_bins"]>=1 || error("Invalid legacy window")
    p["stride"] isa Integer && p["stride"]>=1 || error("Invalid stride")
    p["smooth_radius_px"] isa Integer && p["smooth_radius_px"]>=0 || error("Invalid smoothing radius")
    0<p["plane_rank_rtol"]<1 || error("Invalid plane rank tolerance")
    p["flatten"] in ("none","plane","rows","plane+rows") || error("Unknown flattening")
    p["channel"] isa String && !isempty(p["channel"]) || error("Missing channel")
    c
end

mad_scale(v)=1.4826*median(abs.(v .- median(v)))
qsummary(v,q)=isempty(v) ? (NaN,NaN,NaN) : (quantile(v,q[1]),median(v),quantile(v,q[2]))

"Axis weights enter the covariance once, not squared as in the legacy SVD."
function axis_frame(xs,ys,z,c; legacy=false)
    valid=isfinite.(z); count(valid)>=2 || return nothing
    m=c["model"]; threshold=quantile(z[valid],m["axis_bright_quantile"])
    mask=valid .& (z .>= threshold)
    inds=[CartesianIndex(iy,ix) for iy in axes(z,1) for ix in axes(z,2) if mask[iy,ix]]
    length(inds)>=2 || return nothing
    xb=[xs[i[2]] for i in inds]; yb=[ys[i[1]] for i in inds]; zb=z[inds]
    w=zb .- minimum(zb)
    if legacy
        w .+= m["legacy_weight_floor"]
    elseif maximum(w)<=0
        return nothing
    end
    sw=sum(w); sw>0 || return nothing
    origin=[sum(xb.*w)/sw,sum(yb.*w)/sw]
    X=hcat(xb .- origin[1],yb .- origin[2])
    if legacy
        _,sing,V=svd(Diagonal(w ./ maximum(w))*X)
        axis=collect(V[:,1]); eigenvalues=reverse(sing.^2)
    else
        cov=Symmetric(transpose(X)*(X .* reshape(w,:,1))/sw)
        eig=eigen(cov); eigenvalues=eig.values
        eigenvalues[2]>0 || return nothing
        eigenvalues[2]-eigenvalues[1]>m["axis_degeneracy_rtol"]*eigenvalues[2] || return nothing
        axis=collect(eig.vectors[:,2])
    end
    axis ./= norm(axis)
    (axis[2]<0 || (axis[2]==0 && axis[1]<0)) && (axis .*= -1)
    (;origin,axis,threshold,mask,eigenvalues)
end

"Full strip keeps low pixels; unobserved bins are NA, never manufactured zeros."
function strip_profile(xs,ys,z,frame,c; legacy=false)
    m=c["model"]
    dx=abs(xs[2]-xs[1]); dy=abs(ys[2]-ys[1]); px=legacy ? dx : min(dx,dy)
    half=m["strip_halfwidth_pixels"]*px; binw=m["bin_width_pixels"]*px
    ts=Float64[]; zs=Float64[]
    for iy in eachindex(ys), ix in eachindex(xs)
        legacy && !frame.mask[iy,ix] && continue
        isfinite(z[iy,ix]) || continue
        x=xs[ix]-frame.origin[1]; y=ys[iy]-frame.origin[2]
        u=-x*frame.axis[2]+y*frame.axis[1]
        abs(u)<=half || continue
        push!(ts,x*frame.axis[1]+y*frame.axis[2]); push!(zs,z[iy,ix])
    end
    isempty(ts) && return (;t=Float64[],z=Float64[],n=Int[],binw,px)
    lo,hi=extrema(ts)
    edges=legacy ? collect(lo:binw:hi) : Float64[]
    nb=legacy ? length(edges)-1 : max(1,ceil(Int,(hi-lo)/binw))
    nb>0 || return (;t=Float64[],z=Float64[],n=Int[],binw,px)
    sums=zeros(nb); cnt=zeros(Int,nb)
    for k in eachindex(ts)
        b=clamp((legacy ? Int(fld(ts[k]-lo,binw)) : floor(Int,(ts[k]-lo)/binw))+1,1,nb)
        sums[b]+=zs[k]; cnt[b]+=1
    end
    prof=[cnt[k]>0 ? sums[k]/cnt[k] : (legacy ? 0.0 : NaN) for k in 1:nb]
    cents=legacy ? edges[1:nb] .+ binw/2 : lo .+ ((0:nb-1) .+ .5).*binw
    (;t=cents,z=prof,n=cnt,binw,px)
end

"Plateau centres, excluding endpoints and gaps. No minimum peak count or spacing."
function maxima(y)
    peaks=Int[]; i=2
    while i<length(y)
        if isfinite(y[i-1]) && isfinite(y[i]) && y[i]>y[i-1]
            j=i
            while j<length(y) && y[j+1]==y[i]; j+=1; end
            j<length(y) && isfinite(y[j+1]) && y[j]>y[j+1] && push!(peaks,fld(i+j,2))
            i=j+1
        else
            i+=1
        end
    end
    peaks
end

# Algorithm definition (not a dependency on SciPy):
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_widths.html
"Topographic prominence and linearly interpolated half-prominence width."
function peak_measurement(t,y,k)
    1<k<length(y) || error("Peak must be interior")
    isfinite(y[k]) || error("Nonfinite peak")
    left=k; right=k; lb=k; rb=k
    while left>1 && isfinite(y[left-1]) && y[left-1]<=y[k]
        left-=1
        y[left]<y[lb] && (lb=left)
    end
    while right<length(y) && isfinite(y[right+1]) && y[right+1]<=y[k]
        right+=1
        y[right]<y[rb] && (rb=right)
    end
    baseline=max(y[lb],y[rb]); prominence=y[k]-baseline
    level=y[k]-prominence/2
    l=k; r=k
    while l>lb && y[l]>level; l-=1; end
    while r<rb && y[r]>level; r+=1; end
    lx=t[l]; rx=t[r]
    if l<k && y[l]<level
        lx+=(level-y[l])/(y[l+1]-y[l])*(t[l+1]-t[l])
    end
    if r>k && y[r]<level
        rx-=(level-y[r])/(y[r-1]-y[r])*(t[r]-t[r-1])
    end
    (;index=k,t_nm=t[k],height_nm=y[k],prominence_nm=prominence,
      left_base_nm=t[lb],right_base_nm=t[rb],left_crossing_nm=lx,
      right_crossing_nm=rx,width_nm=rx-lx,contour_nm=baseline,
      evaluation_height_nm=level,segment_left=left,segment_right=right)
end

function observations(profile,c,hf_mad; legacy=false)
    t=profile.t; y=profile.z; s=c["selection"]; q=c["model"]["width_quantiles"]
    peaks=NamedTuple[]; widths=Float64[]; positions=Int[]
    if legacy && !isempty(y)
        threshold=median(y)+s["legacy_height_std_multiplier"]*std(y)
        for k in 2:length(y)-1
            y[k]>threshold && y[k]>=y[k-1] && y[k]>=y[k+1] || continue
            push!(positions,k)
            win=max(1,k-s["legacy_width_window_bins"]):min(length(y),k+s["legacy_width_window_bins"])
            amp=y[k]-minimum(y[win]); half=minimum(y[win])+amp/2
            above=win[y[win] .>= half]
            width=isempty(above) ? NaN : t[min(last(above)+1,length(t))]-t[first(above)]
            valid=amp>0 && width>s["min_width_pixels"]*profile.px && width<s["legacy_max_width_nm"]
            valid && push!(widths,width)
            push!(peaks,(;index=k,t_nm=t[k],height_nm=y[k],prominence_nm=amp,
                width_nm=width,accepted=valid,reason=valid ? "legacy_window_estimate" : "legacy_width_rejected"))
        end
    elseif !legacy
        for k in maxima(y)
            r=peak_measurement(t,y,k)
            contrast=r.prominence_nm>0 && r.prominence_nm>=s["prominence_hf_mad_multiplier"]*hf_mad
            resolved=r.width_nm>s["min_width_pixels"]*profile.px
            valid=contrast && resolved
            if valid; push!(widths,r.width_nm); push!(positions,k); end
            reason=!contrast ? "low_prominence" : (!resolved ? "sampling_limited" : "apparent_only")
            push!(peaks,merge(r,(;accepted=valid,reason)))
        end
    end
    # Never bridge a missing profile interval to assert an observed spacing.
    gaps=Float64[]
    for j in 2:length(positions)
        a,b=positions[j-1],positions[j]
        all(isfinite,@view y[a:b]) && push!(gaps,t[b]-t[a])
    end
    lo,med,hi=qsummary(widths,q)
    (;peaks,n_peaks=length(positions),n_widths=length(widths),width_low_nm=lo,
      width_median_nm=med,width_high_nm=hi,n_spacings=length(gaps),
      spacing_median_nm=isempty(gaps) ? NaN : median(gaps),
      status=isempty(widths) && isempty(gaps) ? "unavailable" : "apparent_only")
end

const SUMMARY_HEADER=["file","direction","method","status","reason","raw_nonfinite",
    "observed_pixels","smoothed_observed_pixels","preprocessing_mode",
    "pixel_x_nm","pixel_y_nm","pipeline_dispersion_nm","hf_mad_nm","axis_angle_deg",
    "profile_bins","empty_bins","observed_peaks","width_samples","width_low_nm",
    "width_median_nm","width_high_nm","spacing_samples","spacing_median_nm",
    "legacy_width_fallback_would_apply","legacy_spacing_fallback_would_apply",
    "legacy_reported_width_low_nm","legacy_reported_width_high_nm","legacy_reported_spacing_nm",
    "production_calibration_emitted"]
const PROFILE_HEADER=["file","direction","method","bin","t_nm","height_nm","pixels"]
const PEAK_HEADER=["file","direction","method","index","t_nm","height_nm","prominence_nm",
    "width_nm","accepted","reason","left_base_nm","right_base_nm","left_crossing_nm",
    "right_crossing_nm","contour_nm","evaluation_height_nm"]

function unavailable(file,direction,method,reason; raw_nonfinite=0)
    r=Dict{String,Any}(k=>NaN for k in SUMMARY_HEADER)
    merge!(r,Dict("file"=>file,"direction"=>direction,"method"=>method,"status"=>"unavailable",
        "reason"=>reason,"raw_nonfinite"=>raw_nonfinite,"profile_bins"=>0,"empty_bins"=>0,
        "observed_pixels"=>0,"smoothed_observed_pixels"=>0,"preprocessing_mode"=>"unavailable",
        "observed_peaks"=>0,"width_samples"=>0,"spacing_samples"=>0,
        "legacy_width_fallback_would_apply"=>false,"legacy_spacing_fallback_would_apply"=>false,
        "production_calibration_emitted"=>false))
    r
end

function measure_image(img,c)
    rows=Dict{String,Any}[]; profiles=Dict{String,Any}[]; peaks=Dict{String,Any}[]
    file=basename(img.filepath); p=c["preprocessing"]; m=c["model"]
    for direction in DIRECTIONS
        channels=filter(ch->lowercase(ch.name)==lowercase(p["channel"]) && lowercase(ch.direction)==direction,img.channels)
        reason=length(channels)!=1 ? "missing_or_ambiguous_direction" : ""
        bad=isempty(reason) ? count(!isfinite,only(channels).data) : 0
        isempty(reason) && bad==length(only(channels).data) && (reason="no_observed_pixels")
        min(img.width,img.height)<2 && (reason="insufficient_pixel_grid")
        if !isempty(reason)
            append!(rows,[unavailable(file,direction,method,reason;raw_nonfinite=bad) for method in METHODS]); continue
        end
        cfg=GaussianFit2D.PatternConfig(filepath=img.filepath,channel=p["channel"],direction=direction,
            stride=p["stride"],flatten=p["flatten"],smooth_radius_px=p["smooth_radius_px"],no_plot=true)
        xs,ys,raw,z,zs,unit,disp=GaussianFit2D.preprocess_channel(img,only(channels),cfg)
        if unit!="nm" || min(length(xs),length(ys))<2 || !all(isfinite,zs)
            append!(rows,[unavailable(file,direction,method,"invalid_units_or_preprocessed_grid") for method in METHODS]); continue
        end
        masked=STMSXMIO.preprocess_observed_channel(img,only(channels);
            stride=p["stride"],flatten=p["flatten"],smooth_radius_px=p["smooth_radius_px"],
            plane_rank_rtol=p["plane_rank_rtol"])
        dx=abs(xs[2]-xs[1]); dy=abs(ys[2]-ys[1])
        for method in METHODS
            legacy=method=="legacy_bootstrap"
            row=unavailable(file,direction,method,"axis_unavailable";raw_nonfinite=bad)
            field=legacy ? zs : masked.z_smooth
            eligible=isfinite.(field)
            count_raw=count(masked.observed); count_smooth=count(masked.smoothed_observed)
            mode=legacy ? "legacy_pipeline_median_imputation" : "observed_only_strict_box"
            merge!(row,Dict("observed_pixels"=>count_raw,"smoothed_observed_pixels"=>count_smooth,
                "preprocessing_mode"=>mode))
            if !legacy && (masked.status!="observed_only" || !any(eligible))
                row["reason"]=masked.status=="observed_only" ? "no_complete_smoothing_footprint" : masked.status
                push!(rows,row); continue
            end
            hf=legacy ? mad_scale(vec(z .- zs)) : mad_scale((masked.z.-field)[eligible])
            merge!(row,Dict("pixel_x_nm"=>dx,"pixel_y_nm"=>dy,"pipeline_dispersion_nm"=>disp,"hf_mad_nm"=>hf))
            frame=axis_frame(xs,ys,field,c;legacy)
            if frame===nothing; push!(rows,row); continue; end
            profile=strip_profile(xs,ys,field,frame,c;legacy)
            ob=observations(profile,c,hf;legacy)
            for i in eachindex(profile.t)
                push!(profiles,Dict("file"=>file,"direction"=>direction,"method"=>method,"bin"=>i,
                    "t_nm"=>profile.t[i],"height_nm"=>profile.z[i],"pixels"=>profile.n[i]))
            end
            for peak in ob.peaks
                pr=Dict{String,Any}(k=>NaN for k in PEAK_HEADER)
                for (k,v) in pairs(peak); haskey(pr,string(k)) && (pr[string(k)]=v); end
                merge!(pr,Dict("file"=>file,"direction"=>direction,"method"=>method)); push!(peaks,pr)
            end
            wf=legacy && ob.n_widths==0; sf=legacy && ob.n_spacings==0
            status=legacy && bad>0 ? "legacy_imputed_"*ob.status : ob.status
            merge!(row,Dict("status"=>status,"reason"=>legacy ? "legacy_local_window_not_physical_FWHM" : "half_prominence_not_physical_FWHM",
                "axis_angle_deg"=>rad2deg(atan(frame.axis[2],frame.axis[1])),
                "profile_bins"=>length(profile.t),"empty_bins"=>count(==(0),profile.n),
                "observed_peaks"=>ob.n_peaks,"width_samples"=>ob.n_widths,
                "width_low_nm"=>ob.width_low_nm,"width_median_nm"=>ob.width_median_nm,"width_high_nm"=>ob.width_high_nm,
                "spacing_samples"=>ob.n_spacings,"spacing_median_nm"=>ob.spacing_median_nm,
                "legacy_width_fallback_would_apply"=>wf,"legacy_spacing_fallback_would_apply"=>sf,
                "legacy_reported_width_low_nm"=>legacy ? (wf ? m["legacy_fwhm_fallback_nm"][1] : ob.width_low_nm) : NaN,
                "legacy_reported_width_high_nm"=>legacy ? (wf ? m["legacy_fwhm_fallback_nm"][2] : ob.width_high_nm) : NaN,
                "legacy_reported_spacing_nm"=>legacy ? (sf ? m["legacy_spacing_fallback_nm"] : ob.spacing_median_nm) : NaN))
            push!(rows,row)
        end
    end
    (;rows,profiles,peaks)
end

function raw_files(dir)
    files=Dict{String,String}()
    for (root,_,names) in walkdir(dir), name in names
        endswith(lowercase(name),".sxm") || continue
        haskey(files,name) && error("Duplicate raw basename $name")
        files[name]=abspath(joinpath(root,name))
    end
    isempty(files) && error("No SXM inputs")
    [files[k] for k in sort(collect(keys(files)))]
end

function options(args)
    args==["--help"] && return nothing
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]
        if !startswith(k,"--") && !haskey(o,"--file") && i==1
            o["--file"]=k; i+=1; continue
        end
        haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in ("--file","--data-dir","--config","--outdir","--chunk") && i<length(args) && !startswith(args[i+1],"--") || error("Forbidden/missing option $k")
        o[k]=args[i+1]; i+=2
    end
    xor(haskey(o,"--file"),haskey(o,"--data-dir")) || error("Exactly one raw file or directory required")
    get!(o,"--config",SETTINGS); settings(o["--config"])
    files=haskey(o,"--file") ? [abspath(o["--file"])] : raw_files(o["--data-dir"])
    all(isfile,files) || error("Missing raw input")
    all(f->endswith(lowercase(f),".sxm"),files) || error("Only SXM input accepted")
    get!(o,"--chunk","1/1")
    ij=try parse.(Int,split(o["--chunk"],'/')) catch; Int[] end
    length(ij)==2 && 1<=ij[1]<=ij[2]<=length(files) || error("Invalid shard")
    if !haskey(o,"--outdir")
        haskey(o,"--file") || error("Directory input requires --outdir")
        o["--outdir"]=joinpath(ROOT,"results/calibration_measurements",splitext(basename(only(files)))[1])
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    o
end

format_value(x)=x isa AbstractFloat ? (isnan(x) ? "NA" : @sprintf("%.17g",x)) : string(x)
function table(path,header,rows)
    open(path,"w") do io
        println(io,join(header,'\t'))
        for r in rows
            values=[format_value(r[k]) for k in header]
            any(v->occursin(r"[\t\r\n]",v),values) && error("Invalid table value")
            println(io,join(values,'\t'))
        end
    end
end

function execute(o; reader=read_sxm)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    c=settings(o["--config"])
    files=haskey(o,"--file") ? [abspath(o["--file"])] : raw_files(o["--data-dir"])
    i,n=parse.(Int,split(o["--chunk"],'/')); selected=files[i:n:end]
    haskey(o,"--dry-run") && return println("Metadata only: ",length(files)," raw scans, shard ",i,"/",n,": ",length(selected),"; no fit, labels or calibration emitted")
    out=abspath(o["--outdir"])
    (ispath(out) || islink(out)) && error("Output exists")
    mkpath(out)
    rows=Dict{String,Any}[]; profiles=Dict{String,Any}[]; peaks=Dict{String,Any}[]
    hashes=Dict{String,Any}[]
    for file in selected
        push!(hashes,Dict("file"=>basename(file),"sha256"=>bytes2hex(sha256(read(file)))))
        r=measure_image(reader(file),c)
        append!(rows,r.rows); append!(profiles,r.profiles); append!(peaks,r.peaks)
        println(basename(file),": ",join([x["direction"]*"/"*x["method"]*"="*x["status"] for x in r.rows],", ")); flush(stdout)
    end
    table(joinpath(out,"measurements.tsv"),SUMMARY_HEADER,rows)
    table(joinpath(out,"profiles.tsv"),PROFILE_HEADER,profiles)
    table(joinpath(out,"peaks.tsv"),PEAK_HEADER,peaks)
    table(joinpath(out,"raw_hashes.tsv"),["file","sha256"],hashes)
    # Exact measurement settings, never a TOML with production fit parameters.
    cp(o["--config"],joinpath(out,"measurement_settings.toml"))
    println("Completed ",length(selected)," raw scans. Apparent measurements only; no production calibration, N selection, assignment or grade.")
    rows
end

function main(args=ARGS)
    o=options(args)
    o===nothing ? println("measure_calibration.jl FILE.sxm [--config SETTINGS] [--outdir NEW] [--dry-run]\n  or --data-dir RAW --outdir NEW [--chunk I/N]\nDiagnostic apparent widths/spacings and legacy fallback audit; no production TOML or expected N.") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && CalibrationMeasurements.main()
