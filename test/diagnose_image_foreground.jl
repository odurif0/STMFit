#!/usr/bin/env julia
# Source-only support first; old anchors and paired agreement are reports only.
module DiagnoseImageForeground
using STMSXMIO, Statistics, TOML, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_vacuum_surface_transfer.jl"))
include(joinpath(@__DIR__,"lib/image_foreground.jl"))
const R=DiagnoseVacuumSurfaceTransfer
const F=ImageForeground
n(r,k)=R.num(r,k)
i(r,k)=R.integer(r,k)
t(path)=R.V.tsv(path)

function put_array(path,a)
    Base.ENDIAN_BOM==0x04030201 || error("Little-endian output required")
    open(io->write(io,a),path,"w")
end

function get_array(dir,name,T,ny,nx)
    path=joinpath(dir,name); filesize(path)==ny*nx*sizeof(T) || error("Wrong array size")
    a=Matrix{T}(undef,ny,nx); open(io->read!(io,a),path); a
end

function direction_data(img,view,c)
    p=c["preprocessing"]
    ch=only(filter(q->lowercase(q.name)==lowercase(p["channel"]) && q.direction==view,img.channels))
    obs=STMSXMIO.preprocess_observed_channel(img,ch;stride=p["stride"],flatten=p["initial_flatten"],
        smooth_radius_px=p["smooth_radius_px"],plane_rank_rtol=p["rank_rtol"])
    obs.unit=="nm" && obs.status=="observed_only" || error("Unavailable initial physical image: "*obs.status)
    r=F.analyse(obs.raw,obs.z_smooth,obs.xs,obs.ys,c)
    (;obs,r)
end

function save_direction(dir,source,c)
    obs,r=source.obs,source.r; ny,nx=size(obs.raw)
    arrays=[("raw_centered.f64",r.centered),("initial_smooth.f64",obs.z_smooth),
        ("joint_smooth.f64",r.final_smooth),("initial_core.u8",UInt8.(r.initial.mask)),
        ("guard.u8",UInt8.(r.guard)),("joint_core.u8",UInt8.(r.final.mask))]
    for (name,a) in arrays; put_array(joinpath(dir,name),a); end
    b=r.background
    meta=Dict("rows"=>ny,"columns"=>nx,"xs_nm"=>obs.xs,"ys_nm"=>obs.ys,
        "source_reference_nm"=>r.reference,"source_hf_scale_nm"=>r.noise,
        "guard_rx"=>r.rx,"guard_ry"=>r.ry,"background_status"=>b.status,
        "slope_normalized_nm"=>b.slope,"xcenter_nm"=>b.xcenter,"xscale_nm"=>b.xscale,
        "within_xx"=>b.within_xx,"within_xz_nm"=>b.within_xz,"background_pixels"=>b.background_pixels)
    open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    R.table(joinpath(dir,"background_rows.tsv"),b.rows)
    summaries=NamedTuple[]; components=NamedTuple[]
    for (arm,seg,cc,z) in (("initial",r.initial,r.initial_components,obs.z_smooth),
                          ("joint_background",r.final,r.final_components,r.final_smooth))
        available=count(isfinite,z); pixels=count(seg.mask)
        push!(summaries,(;arm,status=arm=="initial" ? "ok" : b.status,detector_status=seg.status,
            raw_pixels=count(isfinite,obs.raw),available_pixels=available,foreground_pixels=pixels,
            foreground_fraction=available>0 ? pixels/available : NaN,components=length(cc.rows),
            edge_components=count(q->q.touches_edge,cc.rows),missing_components=count(q->q.touches_missing,cc.rows),
            threshold_nm=seg.threshold,otsu_threshold_nm=seg.otsu_threshold,otsu_score_nm2=seg.otsu_score,
            otsu_split=seg.otsu_split,median_nm=seg.median_nm,floor_threshold_nm=seg.floor_threshold_nm))
        append!(components,[(;arm,row...) for row in cc.rows])
    end
    R.table(joinpath(dir,"summary.tsv"),summaries)
    R.table(joinpath(dir,"components.tsv"),components;header=["arm","component","pixels","first_row","first_column",
        "centroid_x_nm","centroid_y_nm","mean_height_nm","max_height_nm","touches_edge","touches_missing"])
    R.table(joinpath(dir,"arrays.tsv"),[(;file=name,bytes=filesize(joinpath(dir,name)),sha256=R.sha(joinpath(dir,name))) for (name,_) in arrays])
end

"Only after masks are final: annotate every old anchor, never choose a mask from it."
function membership(upstream,dir,source)
    path=joinpath(upstream,"anchors.tsv")
    isfile(path) || return
    old=t(path); rows=NamedTuple[]; obs,r=source.obs,source.r
    for a in old
        y,x=i(a,"iy"),i(a,"ix")
        obs.xs[x]==n(a,"x") && obs.ys[y]==n(a,"y") || error("Changed anchor coordinates")
        push!(rows,(;anchor=i(a,"anchor"),row=y,column=x,x_nm=obs.xs[x],y_nm=obs.ys[y],used=a["used"]=="true",
            initial_available=isfinite(obs.z_smooth[y,x]),joint_available=isfinite(r.final_smooth[y,x]),
            initial_core=r.initial.mask[y,x],joint_core=r.final.mask[y,x],
            initial_component=r.initial_components.labels[y,x],joint_component=r.final_components.labels[y,x]))
    end
    R.table(joinpath(dir,"anchor_membership.tsv"),rows;header=["anchor","row","column","x_nm","y_nm","used",
        "initial_available","joint_available","initial_core","joint_core","initial_component","joint_component"])
    R.table(joinpath(dir,"upstream_anchor_hash.tsv"),[(;file="anchors.tsv",sha256=R.sha(path))])
end

"Descriptive cross-view mask overlap; the target never changes a source mask."
function compare_pair(source,target,base)
    pair=R.F.calibrate_pair(source.obs.raw,target.obs.raw,source.obs.xs,source.obs.ys,base)
    lag=pair.reg.applied_dx_px; ny,nx=size(source.obs.raw)
    avail_source=(isfinite.(source.obs.z_smooth),isfinite.(source.r.final_smooth))
    avail_target=(isfinite.(target.obs.z_smooth),isfinite.(target.r.final_smooth))
    source_masks=(source.r.initial.mask,source.r.final.mask)
    target_masks=(target.r.initial.mask,target.r.final.mask)
    rows=NamedTuple[]
    for (k,arm) in enumerate(("initial","joint_background")), support in ("own","common_methods")
        total=0; a=0; b=0; intersection=0
        for y in 1:ny,x in 1:nx
            xx=x+lag; 1<=xx<=nx || continue
            good=avail_source[k][y,x] && avail_target[k][y,xx]
            if support=="common_methods"
                good &= all(avail_source[j][y,x] && avail_target[j][y,xx] for j in 1:2)
            end
            good || continue
            total+=1; aa=source_masks[k][y,x]; bb=target_masks[k][y,xx]
            a+=aa; b+=bb; intersection+=aa && bb
        end
        union_pixels=a+b-intersection
        push!(rows,(;arm,support,identified=pair.reg.accepted,registration_status=pair.reg.status,
            applied_dx_px=lag,pixels=total,source_foreground=a,target_foreground=b,intersection,
            union_pixels,iou=union_pixels>0 ? intersection/union_pixels : NaN,dice=a+b>0 ? 2intersection/(a+b) : NaN))
    end
    (;rows,pair)
end

function process_file(file,upstream,c,base,out)
    R.S.newdir(out); hash=R.sha(file); img=read_sxm(file)
    open(io->TOML.print(io,Dict("file"=>basename(file),"raw_sha256"=>hash,"bias_v"=>R.bias(img.header))),joinpath(out,"input.toml"),"w")
    sources=Dict{String,Any}(); statuses=NamedTuple[]
    for view in ("fwd","bwd")
        dir=joinpath(out,view); R.S.newdir(dir)
        try
            source=direction_data(img,view,c)
            save_direction(dir,source,c)
            membership(joinpath(upstream,view),dir,source)
            sources[view]=source
            push!(statuses,(;direction=view,status=source.r.background.status))
        catch e
            e isa InterruptException && rethrow()
            push!(statuses,(;direction=view,status="failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' ')))
        end
    end
    if length(sources)==2
        for (view,other) in (("fwd","bwd"),("bwd","fwd"))
            try
                result=compare_pair(sources[view],sources[other],base)
                R.table(joinpath(out,view,"pair_overlap.tsv"),result.rows)
                open(io->TOML.print(io,Dict("status"=>"ok","identified"=>result.pair.reg.accepted,
                    "registration_status"=>result.pair.reg.status,"applied_dx_px"=>result.pair.reg.applied_dx_px,
                    "difference_plane"=>result.pair.coeff)),joinpath(out,view,"pair_status.toml"),"w")
            catch e
                e isa InterruptException && rethrow()
                open(io->TOML.print(io,Dict("status"=>"failed:"*sprint(showerror,e))),joinpath(out,view,"pair_status.toml"),"w")
            end
        end
    end
    R.sha(file)==hash || error("Raw input changed")
    R.table(joinpath(out,"status.tsv"),statuses)
end

function main(args=ARGS)
    "--help" in args && return println("diagnose_image_foreground.jl --root REPO --data-dir RAW --config TOML --outdir NEW_DIR [--file ONE_SXM] [--dry-run]")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    opts=Dict{String,String}(); j=1
    while j<=length(args)
        k=args[j]; haskey(opts,k) && error("Duplicate option")
        if k=="--dry-run"; opts[k]="true"; j+=1; continue; end
        k in ("--root","--data-dir","--config","--outdir","--file") && j<length(args) || error("Unexpected option")
        opts[k]=args[j+1]; j+=2
    end
    all(haskey(opts,k) for k in ("--root","--data-dir","--config","--outdir")) || error("Missing option")
    c=F.settings(opts["--config"]); upstream=joinpath(opts["--root"],c["model"]["transfer_run"])
    R.sha(joinpath(upstream,"settings.toml"))==c["model"]["transfer_settings_sha256"] || error("Changed upstream settings")
    R.sha(joinpath(upstream,"files.tsv"))==c["model"]["files_sha256"] || error("Changed input list")
    base=R.F.settings(joinpath(upstream,"settings.toml")); names=getindex.(t(joinpath(upstream,"files.tsv")),"file")
    actual=sort(filter(f->endswith(lowercase(f),".sxm") && isfile(joinpath(opts["--data-dir"],f)),readdir(opts["--data-dir"])))
    names==actual || error("Raw cohort differs")
    if haskey(opts,"--file"); opts["--file"] in names || error("Unknown file"); names=[opts["--file"]]; end
    if haskey(opts,"--dry-run")
        for file in names
            meta=TOML.parsefile(joinpath(upstream,splitext(file)[1],"input.toml"))
            R.sha(joinpath(opts["--data-dir"],file))==meta["raw_sha256"] || error("Changed raw image")
        end
        println("Metadata/hash only: ",length(names)," raw scans, both directions, all biases, no image decoding or output.")
        return
    end
    out=opts["--outdir"]; R.S.newdir(out); cp(opts["--config"],joinpath(out,"settings.toml"))
    R.table(joinpath(out,"files.tsv"),[(;file=f) for f in names]); BLAS.set_num_threads(1)
    Threads.@threads :static for k in eachindex(names)
        file=names[k]; started=time_ns()
        process_file(joinpath(opts["--data-dir"],file),joinpath(upstream,splitext(file)[1]),c,base,joinpath(out,splitext(file)[1]))
        println("FOREGROUND_FILE_COMPLETE ",file," elapsed_s=",(time_ns()-started)/1e9); flush(stdout)
    end
    println("FOREGROUND_COMPLETE; diagnostic bright support, no counts, labels or promotion")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseImageForeground.main()
