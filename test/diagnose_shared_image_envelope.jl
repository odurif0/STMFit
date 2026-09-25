#!/usr/bin/env julia
# Whole native image, source-only GCV; cross-view scores and DFT probes afterwards.
module DiagnoseSharedImageEnvelope
using LinearAlgebra, SparseArrays, Statistics, TOML
include(joinpath(@__DIR__,"diagnose_image_foreground.jl"))
include(joinpath(@__DIR__,"lib/shared_image_envelope.jl"))
const D=DiagnoseImageForeground
const R=D.R
const F=SharedImageEnvelope
const t=D.t
const n=D.n
const i=D.i
const PROGRESS_LOCK=ReentrantLock()

function load_source(dir)
    m=TOML.parsefile(joinpath(dir,"metadata.toml")); ny,nx=m["rows"],m["columns"]
    arrays=t(joinpath(dir,"arrays.tsv"))
    all(R.sha(joinpath(dir,r["file"]))==r["sha256"] for r in arrays) || error("Changed saved image array")
    raw=D.get_array(dir,"raw_centered.f64",Float64,ny,nx)
    guard=D.get_array(dir,"guard.u8",UInt8,ny,nx).==1
    mask=D.get_array(dir,"joint_core.u8",UInt8,ny,nx).==1
    rows=t(joinpath(dir,"background_rows.tsv")); i.(rows,"row")==collect(1:ny) || error("Changed row order")
    xn=(m["xs_nm"].-m["xcenter_nm"])./m["xscale_nm"]
    level=n.(rows,"row_level_centered_nm"); valid=isfinite.(raw)
    for row in rows; row["usable"]=="true" || (valid[i(row,"row"),:].=false); end
    b=F.background_operator(valid,guard,xn)
    bg=[level[q[1]]+m["slope_normalized_nm"]*xn[q[2]] for q in b.inds]
    all(isfinite,bg) || error("Missing background")
    isapprox(F.background(b,raw[valid]),bg;rtol=1e-9,atol=1e-12) || error("Saved background does not replay")
    z=raw[valid]-bg
    (;meta=m,raw,guard,mask,rows,xn,valid,b,bg,z)
end

"Target values only score; subtract the SOURCE background, never target-fitted rows."
function target_vector(a,target,cal,base)
    m=a.meta; ny,nx=size(a.raw); lag=cal["applied_dx_px"]; plane=cal["difference_plane"]
    size(target.raw)==(ny,nx) && target.meta["xs_nm"]==m["xs_nm"] && target.meta["ys_nm"]==m["ys_nm"] || error("Different target grid")
    holdout=R.F.row_masks(ny,base).test
    v=fill(NaN,length(a.b.inds))
    for (k,q) in enumerate(a.b.inds)
        y,x=Tuple(q); xx=x+lag
        holdout[y] && 1<=xx<=nx && isfinite(target.raw[y,xx]) || continue
        v[k]=target.raw[y,xx]+target.meta["source_reference_nm"]-m["source_reference_nm"]-
            (plane[1]+plane[2]*m["xs_nm"][x]+plane[3]*m["ys_nm"][y])-a.bg[k]
    end
    v
end

function score(pred,target,mask,groups)
    rows=NamedTuple[]; finite=isfinite.(target)
    for group in groups
        good=finite .& (group=="foreground" ? mask : group=="background" ? .!mask : trues(length(mask)))
        push!(rows,(;group,test_pixels=count(good),test_sse_nm2=sum(abs2,target[good]-pred[good])))
    end
    rows
end

function probes(a,fit,transfer,controls,foreground,maps,base,out)
    headers=["interval","patch","anchor","status","pixels","input_sse_nm2","background_sse_nm2",
        "corrected_sse_nm2","envelope_sse_nm2","residual_sse_nm2","remaining_fraction"]
    if !isfile(joinpath(controls,"patches.tsv"))
        R.table(joinpath(out,"probes.tsv"),NamedTuple[];header=headers)
        return (;status="upstream_no_supported_physical_patches",count=0)
    end
    membership=Dict(i(r,"anchor")=>r for r in t(joinpath(foreground,"anchor_membership.tsv")))
    observations=t(joinpath(transfer,"observations.tsv")); groups=Dict{Int,Vector{eltype(observations)}}()
    for r in observations; push!(get!(groups,i(r,"patch"),eltype(observations)[]),r); end
    common=filter(r->r["arm"]=="plane_common",t(joinpath(controls,"patches.tsv")))
    index=zeros(Int,size(a.raw)); index[a.valid]=collect(eachindex(a.z))
    result=NamedTuple[]; coefficients=NamedTuple[]
    for saved in common
        patch=i(saved,"patch"); anchor=i(saved,"anchor"); interval=i(saved,"interval")
        membership[anchor]["joint_core"]=="true" || continue
        rr=groups[patch]; inds=[index[i(r,"row"),i(r,"column")] for r in rr]
        if any(iszero,inds)
            push!(result,(;interval,patch,anchor,status="missing_source_support",pixels=length(inds),
                input_sse_nm2=NaN,background_sse_nm2=NaN,corrected_sse_nm2=NaN,envelope_sse_nm2=NaN,
                residual_sse_nm2=NaN,remaining_fraction=NaN)); continue
        end
        p=(;dx=n.(rr,"dx_nm"),dy=n.(rr,"dy_nm"))
        state=(;angle=n(saved,"angle"),tx=n(saved,"tx"),ty=n(saved,"ty"))
        h=maps.heights[interval]
        values=R.F.prediction(p,h[1],state,base)-R.F.prediction(p,h[0],state,base)
        delta=zeros(length(a.z)); delta[inds]=values
        probe=F.probe(delta,a.b,fit)
        for (j,v) in enumerate(probe.beta); push!(coefficients,(;interval,patch,coefficient=j,value_nm=v)); end
        push!(result,(;interval,patch,anchor,status="ok",pixels=length(inds),input_sse_nm2=probe.input_sse,
            background_sse_nm2=probe.background_sse,corrected_sse_nm2=probe.corrected_sse,
            envelope_sse_nm2=probe.envelope_sse,residual_sse_nm2=probe.residual_sse,
            remaining_fraction=probe.input_sse>0 ? probe.residual_sse/probe.input_sse : NaN))
    end
    R.table(joinpath(out,"probes.tsv"),result;header=headers)
    R.table(joinpath(out,"probe_coefficients.tsv"),coefficients;header=["interval","patch","coefficient","value_nm"])
    (;status="ok",count=length(result))
end

function one_direction(a,target,paths,maps,base,c,out)
    R.S.newdir(out); consumed=vcat([("foreground",name) for name in
        ("metadata.toml","arrays.tsv","background_rows.tsv","pair_status.toml")],
        [("target",name) for name in ("metadata.toml","arrays.tsv")])
    for (role,name) in (("foreground","anchor_membership.tsv"),("transfer","observations.tsv"),("controls","patches.tsv"))
        isfile(joinpath(paths[role],name)) && push!(consumed,(role,name))
    end
    hashes=[(;role,file=name,sha256=R.sha(joinpath(paths[role],name))) for (role,name) in consumed]
    result=F.fit_all(a.z,a.valid,a.meta["xs_nm"],a.meta["ys_nm"],a.b,c)
    # Everything above is independent of target values and chemical templates.
    cal=TOML.parsefile(joinpath(paths["foreground"],"pair_status.toml")); cal["status"]=="ok" || error("Missing frozen calibration")
    target_values=target_vector(a,target,cal,base)
    fits=result.fits; chosen=result.chosen; models=NamedTuple[]; scores=NamedTuple[]; coeffs=NamedTuple[]; arrays=NamedTuple[]
    for (k,fit) in enumerate(fits)
        arm=k==1 ? "background_only" : "spline_$(k-1)"
        if fit===nothing
            push!(models,(;arm,chosen=false,status=result.statuses[k],width_nm=c["model"]["minimum_knot_intervals_nm"][k-1],
                pixels=length(a.z),parameters=0,background_rank=a.b.rank,trace_pb=NaN,degrees_of_freedom=NaN,
                source_sse_nm2=NaN,gcv=Inf,normal=NaN,segments_x=0,segments_y=0,step_x_nm=NaN,step_y_nm=NaN)); continue
        end
        push!(models,(;arm,chosen=k==chosen,status="ok",width_nm=fit.width,pixels=length(a.z),parameters=fit.parameters,
            background_rank=a.b.rank,trace_pb=fit.trace_pb,degrees_of_freedom=fit.df,source_sse_nm2=fit.rss,gcv=fit.gcv,
            normal=fit.normal,segments_x=fit.segments_x,segments_y=fit.segments_y,step_x_nm=fit.step_x,step_y_nm=fit.step_y))
        append!(coeffs,[(;arm,coefficient=j,basis_column=fit.active[j],value_nm=v) for (j,v) in enumerate(fit.beta)])
        pred=fill(NaN,size(a.raw)); pred[a.valid]=fit.pred; name=arm*".f64"; D.put_array(joinpath(out,name),pred)
        push!(arrays,(;file=name,bytes=filesize(joinpath(out,name)),sha256=R.sha(joinpath(out,name))))
        for row in score(fit.pred,target_values,a.mask[a.valid],c["selection"]["groups"])
            push!(scores,(;arm,identified=cal["identified"],chosen=k==chosen,row...))
        end
    end
    for row in score(a.z,target_values,a.mask[a.valid],c["selection"]["groups"])
        push!(scores,(;arm="source_copy",identified=cal["identified"],chosen=false,row...))
    end
    for (name,values) in (("source_corrected.f64",a.z),("target_corrected.f64",target_values))
        z=fill(NaN,size(a.raw)); z[a.valid]=values; D.put_array(joinpath(out,name),z)
        push!(arrays,(;file=name,bytes=filesize(joinpath(out,name)),sha256=R.sha(joinpath(out,name))))
    end
    probe=probes(a,fits[chosen],paths["transfer"],paths["controls"],paths["foreground"],maps,base,out)
    all(R.sha(joinpath(paths[r.role],r.file))==r.sha256 for r in hashes) || error("Changed input during fit")
    R.table(joinpath(out,"input_hashes.tsv"),hashes); R.table(joinpath(out,"arrays.tsv"),arrays)
    R.table(joinpath(out,"models.tsv"),models); R.table(joinpath(out,"scores.tsv"),scores)
    R.table(joinpath(out,"coefficients.tsv"),coeffs;header=["arm","coefficient","basis_column","value_nm"])
    open(io->TOML.print(io,Dict("chosen"=>models[chosen].arm,"identified"=>cal["identified"],
        "probe_status"=>probe.status,"probes"=>probe.count)),joinpath(out,"selection.toml"),"w")
    cp(joinpath(paths["foreground"],"pair_status.toml"),joinpath(out,"calibration.toml"))
    (;status="ok",pixels=length(a.z),probes=probe.count)
end

function process_file(file,paths,maps,base,c,out)
    R.S.newdir(out); stem=splitext(file)[1]; statuses=NamedTuple[]; sources=Dict{String,Any}()
    for view in ("fwd","bwd")
        try; sources[view]=load_source(joinpath(paths["foreground"],stem,view))
        catch e
            e isa InterruptException && rethrow()
            push!(statuses,(;direction=view,status="load_failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' '),pixels=0,probes=0))
        end
    end
    for (view,other) in (("fwd","bwd"),("bwd","fwd"))
        haskey(sources,view) || continue
        try
            p=Dict(k=>joinpath(paths[k],stem,view) for k in keys(paths)); p["target"]=joinpath(paths["foreground"],stem,other)
            result=one_direction(sources[view],sources[other],p,maps,base,c,joinpath(out,view))
            push!(statuses,(;direction=view,result...))
        catch e
            e isa InterruptException && rethrow()
            push!(statuses,(;direction=view,status="failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' '),pixels=0,probes=0))
        end
    end
    R.table(joinpath(out,"status.tsv"),statuses)
end

function main(args=ARGS)
    "--help" in args && return println("diagnose_shared_image_envelope.jl --root REPO --config TOML --outdir NEW_DIR [--file ONE_SXM] [--dry-run]")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    opts=Dict{String,String}(); j=1
    while j<=length(args)
        key=args[j]; haskey(opts,key) && error("Duplicate option")
        if key=="--dry-run"; opts[key]="true"; j+=1; continue; end
        key in ("--root","--config","--outdir","--file") && j<length(args) || error("Unexpected option")
        opts[key]=args[j+1]; j+=2
    end
    all(haskey(opts,k) for k in ("--root","--config","--outdir")) || error("Missing option")
    c=F.settings(opts["--config"]); paths=Dict(k=>joinpath(opts["--root"],c["model"][k*"_run"]) for k in ("foreground","transfer","controls"))
    for (k,dir) in paths
        R.sha(joinpath(dir,"settings.toml"))==c["model"][k*"_settings_sha256"] || error("Changed $k settings")
        R.sha(joinpath(dir,"files.tsv"))==c["model"]["files_sha256"] || error("Changed input list")
    end
    names=getindex.(t(joinpath(paths["foreground"],"files.tsv")),"file")
    if haskey(opts,"--file"); opts["--file"] in names || error("Unknown input"); names=[opts["--file"]]; end
    base=R.F.settings(joinpath(paths["transfer"],"settings.toml")); maps=R.load_maps(opts["--root"],base)
    if haskey(opts,"--dry-run")
        for file in names,view in ("fwd","bwd")
            dir=joinpath(paths["foreground"],splitext(file)[1],view)
            for name in ("metadata.toml","arrays.tsv","background_rows.tsv","pair_status.toml"); isfile(joinpath(dir,name)) || error("Missing input"); end
        end
        println("Metadata/maps only: ",length(names)," files, both directions. No fitting or writes."); return
    end
    out=opts["--outdir"]; R.S.newdir(out); cp(opts["--config"],joinpath(out,"settings.toml"))
    R.table(joinpath(out,"files.tsv"),[(;file=f) for f in names]); R.table(joinpath(out,"map_inputs.tsv"),maps.refs)
    BLAS.set_num_threads(1)
    Threads.@threads :static for k in eachindex(names)
        file=names[k]; started=time_ns(); process_file(file,paths,maps,base,c,joinpath(out,splitext(file)[1]))
        lock(PROGRESS_LOCK) do
            println("SHARED_ENVELOPE_FILE_COMPLETE ",file," elapsed_s=",(time_ns()-started)/1e9); flush(stdout)
        end
    end
    println("SHARED_ENVELOPE_COMPLETE; source-only image diagnostic, no units or promotion")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseSharedImageEnvelope.main()
