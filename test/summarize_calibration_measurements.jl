#!/usr/bin/env julia
# Post-measurement, label-free comparison; never writes a production config.
module CalibrationMeasurementReport
using Statistics, TOML, SHA, Printf

function read_table(path)
    lines=readlines(path); isempty(lines) && error("Empty table: $path")
    header=split(first(lines),'\t')
    length(unique(header))==length(header) || error("Duplicate columns: $path")
    rows=Dict{String,String}[]
    for line in lines[2:end]
        fields=split(line,'\t';keepempty=true)
        length(fields)==length(header) || error("Malformed row: $path")
        push!(rows,Dict(zip(header,fields)))
    end
    rows
end
number(r,k)=r[k]=="NA" ? NaN : parse(Float64,r[k])
integer(r,k)=parse(Int,r[k])
flag(r,k)=parse(Bool,r[k])
key(r)=(r["file"],r["direction"],r["method"])
function grouped(rows)
    out=Dict{Tuple{String,String,String},Vector{Dict{String,String}}}()
    for r in rows; push!(get!(out,key(r),Dict{String,String}[]),r); end
    out
end
check(ok,message)=ok || error(message)
same(a,b)=isnan(a) && isnan(b) || isapprox(a,b;rtol=1e-12,atol=1e-12)
fmt(v)=isempty(v) ? "NA" : @sprintf("%.6g [%.6g, %.6g]",median(v),minimum(v),maximum(v))
validvalues(rows,k)=filter(isfinite,[number(r,k) for r in rows])

"Verify serialized measurements and compare actual views without any labels."
function inspect(root,rawdir,config)
    chunks=sort(filter(p->isdir(joinpath(root,p)) && occursin(r"^chunk[1-9][0-9]*$",p),readdir(root)))
    isempty(chunks) && error("No measurement chunks")
    raw=Dict{String,String}()
    for (dir,_,names) in walkdir(rawdir), name in names
        endswith(lowercase(name),".sxm") || continue
        haskey(raw,name) && error("Duplicate raw basename")
        raw[name]=joinpath(dir,name)
    end
    isempty(raw) && error("No raw files")
    rows=Dict{String,String}[]; profiles=Dict{String,String}[]; peaks=Dict{String,String}[]
    hashes=Dict{String,String}[]
    for chunk in chunks
        dir=joinpath(root,chunk)
        check(Set(readdir(dir))==Set(["measurements.tsv","profiles.tsv","peaks.tsv","raw_hashes.tsv","measurement_settings.toml"]),"Unexpected output in $chunk")
        check(read(joinpath(dir,"measurement_settings.toml"))==read(config),"Settings changed in $chunk")
        append!(rows,read_table(joinpath(dir,"measurements.tsv")))
        append!(profiles,read_table(joinpath(dir,"profiles.tsv")))
        append!(peaks,read_table(joinpath(dir,"peaks.tsv")))
        append!(hashes,read_table(joinpath(dir,"raw_hashes.tsv")))
    end
    check(length(hashes)==length(raw),"Incomplete or duplicated input hashes")
    check(Set(h["file"] for h in hashes)==Set(keys(raw)),"Raw cohort mismatch")
    for h in hashes
        check(bytes2hex(sha256(read(raw[h["file"]])))==h["sha256"],"Raw bytes changed: $(h["file"])")
    end
    methods=("legacy_bootstrap","full_profile_prominence"); directions=("fwd","bwd")
    expected=Set((file,d,m) for file in keys(raw) for d in directions for m in methods)
    check(length(rows)==length(expected) && Set(key.(rows))==expected,"Incomplete or duplicated measurement rows")
    pg=grouped(profiles); kg=grouped(peaks)
    check(issubset(Set(keys(pg)),expected) && issubset(Set(keys(kg)),expected),"Unknown profile or peak")
    settings=TOML.parsefile(config); quantiles=settings["model"]["width_quantiles"]
    for r in rows
        id=key(r); legacy=r["method"]==first(methods)
        check(!flag(r,"production_calibration_emitted"),"Unexpected production calibration")
        ps=sort(get(pg,id,Dict{String,String}[]);by=p->integer(p,"bin"))
        ks=sort(get(kg,id,Dict{String,String}[]);by=p->integer(p,"index"))
        check(length(ps)==integer(r,"profile_bins"),"Profile size mismatch: $id")
        check([integer(p,"bin") for p in ps]==collect(1:length(ps)),"Profile order mismatch: $id")
        check(count(p->integer(p,"pixels")==0,ps)==integer(r,"empty_bins"),"Empty-bin mismatch: $id")
        for p in ps
            check(integer(p,"pixels")>=0,"Negative bin support")
            z=number(p,"height_nm")
            check(integer(p,"pixels")>0 ? isfinite(z) : (legacy ? z==0 : isnan(z)),"Filled or lost profile bin: $id")
        end
        ts=[number(p,"t_nm") for p in ps]; ys=[number(p,"height_nm") for p in ps]
        check(all(>(0),diff(ts)),"Non-increasing profile axis")
        check(length(unique(integer(p,"index") for p in ks))==length(ks),"Duplicate peak")
        widths=Float64[]; indices=Int[]
        for p in ks
            i=integer(p,"index")
            check(1<i<length(ps),"Non-interior peak")
            check(same(ts[i],number(p,"t_nm")) && same(ys[i],number(p,"height_nm")),"Peak/profile mismatch")
            if legacy || flag(p,"accepted"); push!(indices,i); end
            flag(p,"accepted") && push!(widths,number(p,"width_nm"))
            if !legacy
                lx=number(p,"left_crossing_nm"); rx=number(p,"right_crossing_nm")
                check(number(p,"left_base_nm")<=lx<=ts[i]<=rx<=number(p,"right_base_nm"),"Crossing outside bases")
                check(same(rx-lx,number(p,"width_nm")),"Width/crossing mismatch")
                level=number(p,"evaluation_height_nm"); prominence=number(p,"prominence_nm")
                check(same(ys[i]-number(p,"contour_nm"),prominence) && same(ys[i]-prominence/2,level),"Prominence contour mismatch")
                for x in (lx,rx)
                    j=clamp(searchsortedlast(ts,x),1,length(ts)-1)
                    z=ys[j]+(ys[j+1]-ys[j])*(x-ts[j])/(ts[j+1]-ts[j])
                    check(same(z,level),"Interpolated crossing mismatch")
                end
                selected=prominence>0 && prominence>=settings["selection"]["prominence_hf_mad_multiplier"]*number(r,"hf_mad_nm") &&
                    number(p,"width_nm")>settings["selection"]["min_width_pixels"]*min(number(r,"pixel_x_nm"),number(r,"pixel_y_nm"))
                check(selected==flag(p,"accepted"),"Changed descriptive filter")
            end
        end
        check(length(widths)==integer(r,"width_samples") && length(indices)==integer(r,"observed_peaks"),"Peak counts mismatch")
        expected_widths=isempty(widths) ? (NaN,NaN,NaN) : (quantile(widths,quantiles[1]),median(widths),quantile(widths,quantiles[2]))
        for (col,value) in zip(("width_low_nm","width_median_nm","width_high_nm"),expected_widths)
            check(same(number(r,col),value),"Width summary mismatch")
        end
        gaps=Float64[]
        for j in 2:length(indices)
            a,b=indices[j-1],indices[j]
            all(isfinite,ys[a:b]) && push!(gaps,ts[b]-ts[a])
        end
        check(length(gaps)==integer(r,"spacing_samples"),"Spacing count mismatch")
        check(same(number(r,"spacing_median_nm"),isempty(gaps) ? NaN : median(gaps)),"Spacing summary mismatch")
        if legacy && !isempty(ps)
            check(flag(r,"legacy_width_fallback_would_apply")==isempty(widths),"Width fallback mismatch")
            check(flag(r,"legacy_spacing_fallback_would_apply")==isempty(gaps),"Spacing fallback mismatch")
            wl,wh=settings["model"]["legacy_fwhm_fallback_nm"]
            check(same(number(r,"legacy_reported_width_low_nm"),isempty(widths) ? wl : expected_widths[1]),"Legacy low width mismatch")
            check(same(number(r,"legacy_reported_width_high_nm"),isempty(widths) ? wh : expected_widths[3]),"Legacy high width mismatch")
            check(same(number(r,"legacy_reported_spacing_nm"),isempty(gaps) ? settings["model"]["legacy_spacing_fallback_nm"] : median(gaps)),"Legacy spacing mismatch")
        end
    end
    byid=Dict(key(r)=>r for r in rows)
    pairs=Dict{String,Any}[]
    for file in sort(collect(keys(raw))), method in methods
        f=byid[(file,"fwd",method)]; b=byid[(file,"bwd",method)]
        p=Dict{String,Any}("file"=>file,"method"=>method)
        for measure in ("width_median_nm","spacing_median_nm")
            x,y=number(f,measure),number(b,measure)
            p[measure*"_fwd"]=x; p[measure*"_bwd"]=y
            p[measure*"_relative_difference"]=isfinite(x)&&isfinite(y)&&x+y>0 ? 2abs(x-y)/(x+y) : NaN
        end
        x,y=number(f,"axis_angle_deg"),number(b,"axis_angle_deg")
        p["axis_difference_deg"]=isfinite(x)&&isfinite(y) ? min(abs(x-y),180-abs(x-y)) : NaN
        p["peak_count_equal"]=integer(f,"observed_peaks")==integer(b,"observed_peaks")
        push!(pairs,p)
    end
    (;rows,pairs,scans=length(raw),profiles=length(profiles),peaks=length(peaks))
end

function report(io,r)
    println(io,"# Apparent calibration measurements: two-view comparison\n")
    println(io,"$(r.scans) raw scans; $(length(r.rows)) measurement rows; $(r.profiles) profile bins; $(r.peaks) peak records. Raw hashes, settings, complete cohort, profile support, serialized width crossings, filters and summary arithmetic verified. No benchmark label or saved fit/prediction is read.\n")
    println(io,"These are apparent image measurements, not monomer FWHM, physical support bounds, a noise significance test or a recognition score. Methods differ in preprocessing, axis and profile construction, so this is not an isolated causal ablation.\n")
    println(io,"All distribution summaries below are **median [minimum, maximum]**. Paired relative difference is 2|fwd−bwd|/(fwd+bwd); no acceptance threshold is chosen. Unavailable pairs are not filled.\n")
    for method in ("legacy_bootstrap","full_profile_prominence")
        rows=filter(x->x["method"]==method,r.rows); pairs=filter(x->x["method"]==method,r.pairs)
        println(io,"## $method\n")
        for direction in ("fwd","bwd")
            rs=filter(x->x["direction"]==direction,rows)
            println(io,"- $direction: width available $(count(x->isfinite(number(x,"width_median_nm")),rs))/$(r.scans); spacing available $(count(x->isfinite(number(x,"spacing_median_nm")),rs))/$(r.scans).")
        end
        println(io,"- Apparent median width (nm): ",fmt(validvalues(rows,"width_median_nm")))
        println(io,"- Apparent median spacing (nm): ",fmt(validvalues(rows,"spacing_median_nm")))
        for k in ("width_median_nm_relative_difference","spacing_median_nm_relative_difference","axis_difference_deg")
            v=filter(isfinite,[p[k] for p in pairs])
            println(io,"- $k: ",fmt(v),"; $(length(v))/$(r.scans) available pairs.")
        end
        println(io,"- Equal descriptive peak counts: $(count(p->p["peak_count_equal"],pairs))/$(r.scans), including unavailable/empty views; not a count-accuracy grade.")
        println(io,"- Legacy width/spacing fallback flags: $(count(x->flag(x,"legacy_width_fallback_would_apply"),rows)) / $(count(x->flag(x,"legacy_spacing_fallback_would_apply"),rows)).")
        println(io,"- Rows with nonfinite raw pixels: $(count(x->integer(x,"raw_nonfinite")>0,rows))/$(length(rows)).")
        println(io,"- Observed pixels: ",fmt([Float64(integer(x,"observed_pixels")) for x in rows]),"; complete smoothing footprints: ",fmt([Float64(integer(x,"smoothed_observed_pixels")) for x in rows]),".")
        println(io,"\nStatuses:\n")
        for status in sort(unique(x["status"] for x in rows))
            println(io,"- $status: ",count(x->x["status"]==status,rows))
        end
        println(io)
    end
    println(io,"Legacy imputed observations are explicitly flagged; their fallback columns never enter the measured-width summaries. The observed-only method retains holes and rejects incomplete smoothing footprints. The legacy pipeline-dispersion column remains a reference, not a mask-aware noise estimate.\n")
    println(io,"No production calibration, N, unit prediction, benchmark grade or champion changes follow from this report.")
end

function main(args=ARGS)
    length(args)==3 || error("Usage: summarize_calibration_measurements.jl OUTPUT_ROOT RAW_DIR SETTINGS.toml")
    root,raw,config=args
    out=joinpath(root,"comparison.md")
    paired=joinpath(root,"paired_measurements.tsv")
    (ispath(out) || ispath(paired)) && error("Report exists")
    r=inspect(root,raw,config)
    header=["file","method","width_median_nm_fwd","width_median_nm_bwd",
        "width_median_nm_relative_difference","spacing_median_nm_fwd",
        "spacing_median_nm_bwd","spacing_median_nm_relative_difference",
        "axis_difference_deg","peak_count_equal"]
    open(paired,"w") do io
        println(io,join(header,'\t'))
        for p in r.pairs
            println(io,join([p[k] isa AbstractFloat && isnan(p[k]) ? "NA" : string(p[k]) for k in header],'\t'))
        end
    end
    open(io->report(io,r),out,"w")
    report(stdout,r)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && CalibrationMeasurementReport.main()
