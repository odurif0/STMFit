#!/usr/bin/env julia
# Saved-output checks and descriptive aggregation only: no image/model refit.
module VerifyImageForeground
using Test, Statistics, TOML
include(joinpath(@__DIR__,"diagnose_image_foreground.jl"))
const D=DiagnoseImageForeground
const R=D.R
const t=D.t
const n=D.n
const i=D.i
b(r,k)=parse(Bool,r[k])
close(a,b)=isequal(a,b) || isapprox(a,b;rtol=1e-9,atol=1e-12)

function load_direction(dir)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml")); ny,nx=meta["rows"],meta["columns"]
    arrays=Dict(r["file"]=>D.get_array(dir,r["file"],endswith(r["file"],"f64") ? Float64 : UInt8,ny,nx)
        for r in t(joinpath(dir,"arrays.tsv")))
    (;meta,arrays)
end

"Independent integral-image guard, not the detector's separable implementation."
function guard_reference(mask,rx,ry)
    ny,nx=size(mask); acc=zeros(Int,ny+1,nx+1)
    acc[2:end,2:end]=cumsum(cumsum(Int.(mask);dims=1);dims=2)
    [begin
        lo,hi=max(1,y-ry),min(ny,y+ry); left,right=max(1,x-rx),min(nx,x+rx)
        acc[hi+1,right+1]-acc[lo,right+1]-acc[hi+1,left]+acc[lo,left]>0
    end for y in 1:ny,x in 1:nx]
end

function verify_direction(dir,c)
    a=load_direction(dir); m,z=a.meta,a.arrays; ny,nx=m["rows"],m["columns"]
    expected=Set(["raw_centered.f64","initial_smooth.f64","joint_smooth.f64","initial_core.u8","joint_core.u8","guard.u8"])
    @test Set(keys(z))==expected
    hashes=t(joinpath(dir,"arrays.tsv")); @test length(hashes)==6
    @test all(R.sha(joinpath(dir,r["file"]))==r["sha256"] && filesize(joinpath(dir,r["file"]))==i(r,"bytes") for r in hashes)
    raw=z["raw_centered.f64"]; xs=m["xs_nm"]; ys=m["ys_nm"]; xn=(xs.-m["xcenter_nm"])./m["xscale_nm"]
    @test length(xs)==nx && length(ys)==ny && all(diff(xs).>0) && all(diff(ys).>0)
    @test all(q in (0,1) for name in ("initial_core.u8","joint_core.u8","guard.u8") for q in z[name])
    differences=[raw[y,x]-raw[y,x-1] for y in 1:ny for x in 2:nx if isfinite(raw[y,x]) && isfinite(raw[y,x-1])]
    @test close(m["source_hf_scale_nm"],1.4826median(abs.(differences.-median(differences)))/sqrt(2))
    @test abs(median(raw[isfinite.(raw)]))<=1e-12
    guard=z["guard.u8"].==1
    @test m["guard_rx"]==ceil(Int,c["model"]["guard_padding_nm"]/(xs[2]-xs[1]))
    @test m["guard_ry"]==ceil(Int,c["model"]["guard_padding_nm"]/(ys[2]-ys[1]))
    @test guard==guard_reference(z["initial_core.u8"].==1,m["guard_rx"],m["guard_ry"])
    rows=t(joinpath(dir,"background_rows.tsv")); @test i.(rows,"row")==collect(1:ny)
    corrected=fill(NaN,ny,nx); used=0; xx=0.; xz=0.; xnormal=0.
    for row in rows
        y=i(row,"row"); inds=findall(isfinite.(raw[y,:]).&.!guard[y,:])
        @test length(inds)==i(row,"background_pixels")
        eligible=length(inds)>=c["selection"]["min_row_background_pixels"]
        @test b(row,"usable")== (eligible && m["background_status"]=="ok")
        eligible || continue
        used+=length(inds); xm=mean(xn[inds]); zm=mean(raw[y,inds])
        @test close(xm,n(row,"xmean")) && close(zm,n(row,"zmean_centered_nm"))
        xx+=sum(abs2,xn[inds].-xm); xz+=sum((xn[inds].-xm).*(raw[y,inds].-zm))
        if b(row,"usable")
            level=n(row,"row_level_centered_nm"); slope=m["slope_normalized_nm"]
            @test close(level,zm-slope*xm)
            for x in 1:nx
                isfinite(raw[y,x]) && (corrected[y,x]=raw[y,x]-level-slope*xn[x])
            end
            residual=corrected[y,inds]
            @test abs(sum(residual))<=1e-10length(inds)
            xnormal+=sum(residual.*xn[inds])
        end
    end
    @test used==m["background_pixels"] && close(xx,m["within_xx"]) && close(xz,m["within_xz_nm"])
    status=used<max(c["selection"]["min_background_pixels"],c["selection"]["min_background_fraction"]*length(raw)) ? "insufficient_background" :
        !(isfinite(xx) && isfinite(xz) && xx>c["preprocessing"]["rank_rtol"]*used) ? "rank_failure" : "ok"
    @test status==m["background_status"]
    @test status!="ok" || (close(m["slope_normalized_nm"],xz/xx) && abs(xnormal)<=1e-10used)
    radius=c["preprocessing"]["smooth_radius_px"]
    smooth=[mean(corrected[max(1,y-radius):min(ny,y+radius),max(1,x-radius):min(nx,x+radius)]) for y in 1:ny,x in 1:nx]
    @test all(close(smooth[k],z["joint_smooth.f64"][k]) for k in eachindex(smooth))
    summaries=t(joinpath(dir,"summary.tsv")); @test getindex.(summaries,"arm")==["initial","joint_background"]
    components=t(joinpath(dir,"components.tsv")); labels=Dict{String,Matrix{Int}}()
    for (arm,prefix) in (("initial","initial"),("joint_background","joint"))
        q=only(r for r in summaries if r["arm"]==arm); zz=z[prefix*"_smooth.f64"]; mask=z[prefix*"_core.u8"].==1
        values=sort(zz[isfinite.(zz)]); available=length(values)
        @test i(q,"raw_pixels")==count(isfinite,raw) && i(q,"available_pixels")==available
        @test i(q,"foreground_pixels")==count(mask) && mask==(isfinite.(zz).&(zz.>n(q,"threshold_nm")))
        @test all(!mask[k] && !isfinite(zz[k]) for k in eachindex(raw) if !isfinite(raw[k]))
        if available>0
            med=median(values); floor=med+max(c["selection"]["noise_multiplier"]*m["source_hf_scale_nm"],c["selection"]["numeric_height_floor_nm"])
            @test close(n(q,"median_nm"),med) && close(n(q,"floor_threshold_nm"),floor)
            @test close(n(q,"threshold_nm"),max(floor,n(q,"otsu_threshold_nm")))
            if first(values)<last(values)
                # Alternate variance formula, centered at the median. No segmentation call.
                centered=values.-med; sums=cumsum(centered); total=sums[end]; N=length(values)
                scores=[values[j]<values[j+1] ? (sums[j]-j*total/N)^2/(j*(N-j)) : -Inf for j in 1:N-1]
                split=i(q,"otsu_split")
                @test 1<=split<N && values[split]<values[split+1]
                @test close(scores[split],maximum(scores)) && close(scores[split],n(q,"otsu_score_nm2"))
                @test close(n(q,"otsu_threshold_nm"),(values[split]+values[split+1])/2)
            else
                @test i(q,"otsu_split")==0 && n(q,"otsu_threshold_nm")==first(values) && !any(mask)
            end
        else
            @test q["detector_status"]=="unavailable" && !any(mask)
        end
        # Replay the separately synthetic-tested topology, not selection or fitting.
        cc=D.F.components(mask,zz,xs,ys); labels[prefix]=cc.labels
        saved=filter(r->r["arm"]==arm,components)
        @test length(saved)==length(cc.rows)==i(q,"components")
        for (savedrow,row) in zip(saved,cc.rows), key in keys(row)
            value=getproperty(row,key); textkey=string(key)
            @test value isa Bool ? b(savedrow,textkey)==value : close(n(savedrow,textkey),value)
        end
    end
    path=joinpath(dir,"anchor_membership.tsv")
    if isfile(path)
        for row in t(path)
            y,x=i(row,"row"),i(row,"column")
            @test n(row,"x_nm")==xs[x] && n(row,"y_nm")==ys[y]
            for prefix in ("initial","joint")
                @test b(row,prefix*"_available")==isfinite(z[prefix*"_smooth.f64"][y,x])
                @test b(row,prefix*"_core")==Bool(z[prefix*"_core.u8"][y,x])
                @test i(row,prefix*"_component")==labels[prefix][y,x]
            end
        end
    end
    a
end

function verify_pair(dir,a,other)
    status=TOML.parsefile(joinpath(dir,"pair_status.toml"))
    status["status"]=="ok" || return NamedTuple[]
    rows=t(joinpath(dir,"pair_overlap.tsv")); @test length(rows)==4
    ny,nx=size(a.arrays["raw_centered.f64"]); lag=status["applied_dx_px"]
    for row in rows
        @test i(row,"applied_dx_px")==lag && b(row,"identified")==status["identified"]
        prefix=row["arm"]=="initial" ? "initial" : "joint"
        counts=zeros(Int,4)
        for y in 1:ny,x in 1:nx
            xx=x+lag; 1<=xx<=nx || continue
            prefixes=row["support"]=="common_methods" ? ("initial","joint") : (prefix,)
            all(isfinite(a.arrays[p*"_smooth.f64"][y,x]) && isfinite(other.arrays[p*"_smooth.f64"][y,xx]) for p in prefixes) || continue
            aa=a.arrays[prefix*"_core.u8"][y,x]==1; bb=other.arrays[prefix*"_core.u8"][y,xx]==1
            counts .+= (1,aa,bb,aa&&bb)
        end
        @test counts==[i(row,k) for k in ("pixels","source_foreground","target_foreground","intersection")]
        union_pixels=counts[2]+counts[3]-counts[4]
        @test union_pixels==i(row,"union_pixels")
        @test close(n(row,"iou"),union_pixels>0 ? counts[4]/union_pixels : NaN)
        @test close(n(row,"dice"),counts[2]+counts[3]>0 ? 2counts[4]/(counts[2]+counts[3]) : NaN)
    end
    rows
end

function predictive_rows(dir,old,controls,file,view)
    isfile(joinpath(old,"patches.tsv")) || return (NamedTuple[],NamedTuple[])
    anchors=t(joinpath(dir,"anchor_membership.tsv")); byanchor=Dict(i(r,"anchor")=>r for r in anchors)
    original=t(joinpath(old,"anchors.tsv")); @test length(original)==length(byanchor)==length(anchors)
    @test only(t(joinpath(dir,"upstream_anchor_hash.tsv")))["sha256"]==R.sha(joinpath(old,"anchors.tsv"))
    @test all(haskey(byanchor,i(r,"anchor")) && b(r,"used")==b(byanchor[i(r,"anchor")],"used") for r in original)
    identified=TOML.parsefile(joinpath(old,"calibration.toml"))["identified"]
    upstream=t(joinpath(old,"patches.tsv")); extra=t(joinpath(controls,"patches.tsv")); rows=NamedTuple[]
    arms=[("original_"*arm,filter(r->r["arm"]==arm,upstream),"test_sse_nm2") for arm in ("common","chemical")]
    append!(arms,[(arm,filter(r->r["arm"]==arm,extra),"test_sse_nm2") for arm in ("height_contrast","local_constant","local_plane","plane_common","plane_chemical")])
    push!(arms,("source_copy",filter(r->r["arm"]=="common",upstream),"source_copy_sse_nm2"))
    for interval in 1:3
        reference=filter(r->i(r,"interval")==interval && r["arm"]=="common",upstream)
        keys=Set((i(r,"patch"),i(r,"anchor"),i(r,"test_pixels")) for r in reference)
        for (arm,patches,ssekey) in arms
            selected=filter(r->i(r,"interval")==interval,patches)
            @test length(selected)==length(keys) && Set((i(r,"patch"),i(r,"anchor"),i(r,"test_pixels")) for r in selected)==keys
            for (mask,prefix) in (("initial","initial"),("joint_background","joint"))
                groups=Dict(g=>Dict{String,String}[] for g in ("retained","rejected","unavailable"))
                for patch in selected
                    a=byanchor[i(patch,"anchor")]; @test b(a,"used")
                    group=!b(a,prefix*"_available") ? "unavailable" : b(a,prefix*"_core") ? "retained" : "rejected"
                    push!(groups[group],patch)
                end
                @test sum(length,values(groups))==length(selected)
                groups["all"]=selected
                for group in ("all","retained","rejected","unavailable")
                    rr=groups[group]; pixels=sum((i(r,"test_pixels") for r in rr);init=0)
                    sse=sum((n(r,ssekey) for r in rr);init=0.)
                    push!(rows,(;file,view,identified,interval,mask,group,arm,patches=length(rr),test_pixels=pixels,test_sse_nm2=sse,
                        rms_pm=pixels>0 ? 1000sqrt(sse/pixels) : NaN))
                end
            end
        end
    end
    hashes=[(;file,view,input=p,sha256=R.sha(p)) for p in (joinpath(old,"patches.tsv"),joinpath(controls,"patches.tsv"))]
    rows,hashes
end

function aggregate_predictions(rows)
    out=NamedTuple[]
    for subset in ("all_scored","identified"),interval in 1:3,mask in ("initial","joint_background"),
        group in ("all","retained","rejected","unavailable"),arm in sort(unique(getproperty.(rows,:arm)))
        rr=filter(r->(subset=="all_scored" || r.identified) && r.interval==interval && r.mask==mask && r.group==group && r.arm==arm,rows)
        pixels=sum((r.test_pixels for r in rr);init=0); sse=sum((r.test_sse_nm2 for r in rr);init=0.)
        push!(out,(;subset,interval,mask,group,arm,files=length(unique(r.file for r in rr if r.test_pixels>0)),
            views=count(r->r.test_pixels>0,rr),patches=sum((r.patches for r in rr);init=0),test_pixels=pixels,test_sse_nm2=sse,
            rms_pm=pixels>0 ? 1000sqrt(sse/pixels) : NaN))
    end
    out
end

function main(args=ARGS)
    length(args)==4 || error("verify_image_foreground.jl REPO SAVED_RUN RAW_DIR NEW_REPORT_DIR")
    root,run,rawdir,out=args; c=D.F.settings(joinpath(run,"settings.toml")); R.S.newdir(out)
    old=joinpath(root,c["model"]["transfer_run"]); controls=joinpath(root,c["model"]["controls_run"])
    names=getindex.(t(joinpath(run,"files.tsv")),"file"); expected=getindex.(t(joinpath(old,"files.tsv")),"file")
    summaries=NamedTuple[]; pairs=NamedTuple[]; predictions=NamedTuple[]; hashes=NamedTuple[]; statusrows=NamedTuple[]
    @testset "Saved source-only support and unchanged-prediction aggregation" begin
        @test R.sha(joinpath(old,"settings.toml"))==c["model"]["transfer_settings_sha256"]
        @test R.sha(joinpath(old,"files.tsv"))==c["model"]["files_sha256"]
        @test names==expected || names==[first(expected)]
        for file in names
            stem=splitext(file)[1]; dir=joinpath(run,stem); input=joinpath(old,stem)
            meta=TOML.parsefile(joinpath(dir,"input.toml")); upstream=TOML.parsefile(joinpath(input,"input.toml"))
            @test meta["file"]==file && meta["raw_sha256"]==upstream["raw_sha256"]==R.sha(joinpath(rawdir,file))
            statuses=t(joinpath(dir,"status.tsv")); @test getindex.(statuses,"direction")==["fwd","bwd"]
            loaded=Dict{String,Any}()
            for status in statuses
                view=status["direction"]; path=joinpath(dir,view)
                push!(statusrows,(;file,view,status=status["status"]))
                startswith(status["status"],"failed:") && continue
                loaded[view]=verify_direction(path,c)
                @test loaded[view].meta["background_status"]==status["status"]
                append!(summaries,[(;file,view,(Symbol(k)=>v for (k,v) in row)...) for row in t(joinpath(path,"summary.tsv"))])
                pp,hh=predictive_rows(path,joinpath(input,view),joinpath(controls,stem,view),file,view)
                append!(predictions,pp); append!(hashes,hh)
            end
            if length(loaded)==2
                for (view,other) in (("fwd","bwd"),("bwd","fwd"))
                    rr=verify_pair(joinpath(dir,view),loaded[view],loaded[other])
                    append!(pairs,[(;file,view,(Symbol(k)=>v for (k,v) in row)...) for row in rr])
                    previous=joinpath(input,view,"calibration.toml")
                    if isfile(previous) && !isempty(rr)
                        cal=TOML.parsefile(previous); now=TOML.parsefile(joinpath(dir,view,"pair_status.toml"))
                        if cal["status"]=="ok"
                            @test cal["identified"]==now["identified"] && cal["applied_dx_px"]==now["applied_dx_px"]
                            @test all(close.(cal["difference_plane"],now["difference_plane"]))
                        end
                    end
                end
            end
            println("VERIFIED ",file); flush(stdout)
        end
    end
    R.table(joinpath(out,"statuses.tsv"),statusrows); R.table(joinpath(out,"support.tsv"),summaries)
    R.table(joinpath(out,"pair_overlap.tsv"),pairs); R.table(joinpath(out,"prediction_cases.tsv"),predictions)
    R.table(joinpath(out,"prediction_summary.tsv"),aggregate_predictions(predictions))
    R.table(joinpath(out,"prediction_input_hashes.tsv"),hashes)
    println("Verified ",length(names)," scans. Conditional saved errors only; no new fit, grade or promotion.")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyImageForeground.main()
