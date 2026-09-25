#!/usr/bin/env julia
# Reconstruct saved source-fitted candidates; no fitting/selection and no labels.
module DiagnoseFixedShapeReplay
using TOML, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_vacuum_surface_controls.jl"))
const L=DiagnoseVacuumSurfaceControls
const R=L.R
const t=L.t
const n=L.n
const i=L.i
const PROGRESS_LOCK=ReentrantLock()

function progress(io,file,seconds)
    lock(PROGRESS_LOCK) do
        println(io,"FIXED_SHAPE_FILE_COMPLETE ",file," elapsed_s=",seconds)
        flush(io)
    end
end

function settings(path)
    c=TOML.parsefile(path)
    keysby=Dict("model"=>["transfer_run","controls_run","foreground_run","transfer_settings_sha256",
        "controls_settings_sha256","foreground_settings_sha256","files_sha256","arms","geometry","height_gain"],
        "selection"=>["policy","mask_use","masks","groups"],
        "preprocessing"=>["observations","coefficients","target_use","replay_rtol","replay_atol_nm","replay_atol_nm2"])
    Set(keys(c))==Set(keys(keysby)) || error("Changed sections")
    for (section,fields) in keysby; Set(keys(c[section]))==Set(fields) || error("Unknown/missing setting"); end
    m,s,p=c["model"],c["selection"],c["preprocessing"]
    m["arms"]==["fixed_glcn","fixed_glcnac","source_selected","local_plane","plane_common"] &&
        m["geometry"]=="saved_source_plane_common_argmin" && !(m["height_gain"] isa Bool) && m["height_gain"]==1 || error("Changed replay model")
    s["policy"]=="all_saved_candidates_no_new_selection_no_grade" &&
        s["mask_use"]=="post_scoring_center_membership_only" && s["masks"]==["initial","joint_background"] &&
        s["groups"]==["all","retained","rejected","unavailable"] || error("Changed report population")
    p["observations"]=="saved_native_pixels_no_refiltering" && p["coefficients"]=="saved_source_plane_no_refit" &&
        p["target_use"]=="score_only_never_fit_or_choose" || error("Changed observation boundary")
    for key in ("transfer_settings_sha256","controls_settings_sha256","foreground_settings_sha256","files_sha256")
        occursin(r"^[0-9a-f]{64}$",m[key]) || error("Invalid hash")
    end
    for key in ("replay_rtol","replay_atol_nm","replay_atol_nm2")
        v=p[key]; v isa Real && !(v isa Bool) && isfinite(v) && 0<v<1 || error("Invalid replay tolerance")
    end
    c
end

state(row)=(;type=i(row,"type"),angle=n(row,"angle"),tx=n(row,"tx"),ty=n(row,"ty"))
beta(row)=[n(row,"beta0_nm"),n(row,"betax"),n(row,"betay")]
close(a,b,c;coefficient=false)=isapprox(a,b;rtol=c["preprocessing"]["replay_rtol"],
    atol=c["preprocessing"][coefficient ? "replay_atol_nm" : "replay_atol_nm2"])

"No target, mask, type choice or optimization in this reconstruction."
function prediction(p,map,geometry,coeff,base)
    length(coeff)==3 && all(isfinite,coeff) || error("Invalid saved coefficients")
    R.F.prediction(p,map,geometry,base) .+ coeff[1] .+ coeff[2].*p.dx .+ coeff[3].*p.dy
end

function replay_patch(p,target,maps,geometry,candidates,selected,oldplane,oldcommon,base,c)
    length(target)==length(p.values)==length(p.dx)==length(p.dy) || error("Changed pixel support")
    all(isfinite,p.values) || error("Nonfinite source")
    length(candidates)==2 && Set(i(r,"type") for r in candidates)==Set([0,1]) || error("Both candidates required")
    chosen=only(filter(r->r["chosen"]=="true",candidates))
    i(chosen,"type")==i(selected,"type") && i(chosen,"state")==i(selected,"choice") || error("Changed saved choice")
    all(n(selected,k)==getproperty(geometry,Symbol(k)) for k in ("angle","tx","ty")) || error("Changed chemical geometry")
    all(close(beta(chosen)[j],beta(selected)[j],c;coefficient=true) for j in 1:3) || error("Changed selected plane")
    predictions=Dict{String,Vector{Float64}}(); coefficients=Dict{String,Vector{Float64}}(); types=Dict{String,Int}()
    for (type,arm) in ((0,"fixed_glcn"),(1,"fixed_glcnac"))
        row=only(filter(r->i(r,"type")==type,candidates)); coeff=beta(row)
        i(row,"state")==type+1 || error("Changed candidate identity")
        pred=prediction(p,maps[type],geometry,coeff,base)
        close(sum(abs2,p.values-pred),n(row,"training_sse_nm2"),c) || error("Saved source candidate replay differs")
        predictions[arm]=pred; coefficients[arm]=coeff; types[arm]=type
    end
    arm=i(selected,"type")==0 ? "fixed_glcn" : "fixed_glcnac"
    predictions["source_selected"]=predictions[arm]; coefficients["source_selected"]=coefficients[arm]; types["source_selected"]=types[arm]
    coefficients["local_plane"]=beta(oldplane); types["local_plane"]=-1
    coeff=coefficients["local_plane"]
    predictions["local_plane"]=coeff[1] .+ coeff[2].*p.dx .+ coeff[3].*p.dy
    coefficients["plane_common"]=beta(oldcommon); types["plane_common"]=-1
    predictions["plane_common"]=prediction(p,maps[-1],geometry,beta(oldcommon),base)
    good=isfinite.(target); rows=NamedTuple[]
    for arm in c["model"]["arms"]
        pred=predictions[arm]; coeff=coefficients[arm]
        train=sum(abs2,p.values-pred); test=sum(abs2,target[good]-pred[good])
        push!(rows,(;arm,type=types[arm],angle=geometry.angle,tx=geometry.tx,ty=geometry.ty,
            beta0_nm=coeff[1],betax=coeff[2],betay=coeff[3],training_pixels=length(p.values),training_sse_nm2=train,
            test_pixels=count(good),test_sse_nm2=test))
        if arm in ("source_selected","local_plane","plane_common")
            old=arm=="source_selected" ? selected : arm=="local_plane" ? oldplane : oldcommon
            i(old,"training_pixels")==length(p.values) && i(old,"test_pixels")==count(good) || error("Changed old population")
            close(train,n(old,"training_sse_nm2"),c) && close(test,n(old,"test_sse_nm2"),c) || error("Saved selected/control error changed")
        end
    end
    (;rows,predictions)
end

function one_direction(transfer,controls,foreground,maps,base,c,out)
    R.S.newdir(out)
    inputs=vcat([("transfer",f) for f in ("observations.tsv","calibration.toml")],
        [("controls",f) for f in ("patches.tsv","input_hashes.tsv",("costs_$(j)_plane_chemical.tsv" for j in base["model"]["intervals"])...)],
        [("foreground","anchor_membership.tsv")])
    dirs=Dict("transfer"=>transfer,"controls"=>controls,"foreground"=>foreground)
    hashes=[(;role,file,sha256=R.sha(joinpath(dirs[role],file))) for (role,file) in inputs]
    previous=t(joinpath(controls,"input_hashes.tsv"))
    for file in ("observations.tsv","calibration.toml")
        R.sha(joinpath(transfer,file))==only(filter(r->r["file"]==file,previous))["sha256"] || error("Old control input changed")
    end
    obs=L.read_observations(joinpath(transfer,"observations.tsv")); old=t(joinpath(controls,"patches.tsv"))
    indexed=Dict((i(r,"interval"),r["arm"],i(r,"patch"))=>r for r in old)
    length(indexed)==length(old)==15length(obs.patches) || error("Changed controls population")
    anchors=t(joinpath(foreground,"anchor_membership.tsv")); byanchor=Dict(i(r,"anchor")=>r for r in anchors)
    length(byanchor)==length(anchors) || error("Duplicate anchor membership")
    cal=TOML.parsefile(joinpath(transfer,"calibration.toml")); rows=NamedTuple[]
    for interval in base["model"]["intervals"]
        costs=t(joinpath(controls,"costs_$(interval)_plane_chemical.tsv"))
        length(costs)==2length(obs.patches) && Set(i.(costs,"patch"))==Set(eachindex(obs.patches)) || error("Incomplete candidates")
        for (j,(patch,target)) in enumerate(zip(obs.patches,obs.targets))
            selected=indexed[(interval,"plane_chemical",j)]; common=indexed[(interval,"plane_common",j)]
            plain=indexed[(interval,"local_plane",j)]
            all(i(r,"anchor")==patch.anchor for r in (selected,common,plain)) || error("Changed anchor")
            a=byanchor[patch.anchor]; a["used"]=="true" || error("Unexpected unscored anchor")
            result=replay_patch(patch,target,maps.heights[interval],state(common),filter(r->i(r,"patch")==j,costs),selected,plain,common,base,c)
            # Mask annotations are read only after the five predictions are complete.
            for row in result.rows
                push!(rows,(;interval,patch=j,anchor=patch.anchor,identified=cal["identified"],row...,
                    initial_available=parse(Bool,a["initial_available"]),initial_core=parse(Bool,a["initial_core"]),
                    joint_available=parse(Bool,a["joint_available"]),joint_core=parse(Bool,a["joint_core"])))
            end
        end
    end
    all(R.sha(joinpath(dirs[r.role],r.file))==r.sha256 for r in hashes) || error("Input changed during replay")
    R.table(joinpath(out,"input_hashes.tsv"),hashes); R.table(joinpath(out,"patches.tsv"),rows)
    cp(joinpath(transfer,"calibration.toml"),joinpath(out,"calibration.toml"))
    (;patches=length(obs.patches),status="ok")
end

function process_file(file,paths,maps,base,c,out)
    R.S.newdir(out); stem=splitext(file)[1]; statuses=t(joinpath(paths["controls"],stem,"status.tsv")); rows=NamedTuple[]
    cp(joinpath(paths["transfer"],stem,"input.toml"),joinpath(out,"upstream_input.toml"))
    cp(joinpath(paths["controls"],stem,"status.tsv"),joinpath(out,"upstream_status.tsv"))
    for row in statuses
        view=row["direction"]
        if row["status"]!="ok"
            push!(rows,(;direction=view,status=row["status"],patches=i(row,"patches"))); continue
        end
        try
            result=one_direction(joinpath(paths["transfer"],stem,view),joinpath(paths["controls"],stem,view),
                joinpath(paths["foreground"],stem,view),maps,base,c,joinpath(out,view))
            result.patches==i(row,"patches") || error("Lost patches")
            push!(rows,(;direction=view,result...))
        catch e
            e isa InterruptException && rethrow()
            push!(rows,(;direction=view,status="failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' '),patches=i(row,"patches")))
        end
    end
    R.table(joinpath(out,"status.tsv"),rows)
end

function main(args=ARGS)
    "--help" in args && return println("diagnose_fixed_shape_replay.jl --root REPO --config TOML --outdir NEW_DIR [--file FIRST_SXM] [--dry-run]")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    opts=Dict{String,String}(); j=1
    while j<=length(args)
        key=args[j]; haskey(opts,key) && error("Duplicate option")
        if key=="--dry-run"; opts[key]="true"; j+=1; continue; end
        key in ("--root","--config","--outdir","--file") && j<length(args) || error("Unexpected option")
        opts[key]=args[j+1]; j+=2
    end
    all(haskey(opts,k) for k in ("--root","--config","--outdir")) || error("Missing option")
    c=settings(opts["--config"]); paths=Dict(k=>joinpath(opts["--root"],c["model"][k*"_run"]) for k in ("transfer","controls","foreground"))
    for k in keys(paths)
        R.sha(joinpath(paths[k],"settings.toml"))==c["model"][k*"_settings_sha256"] || error("Changed $k settings")
        R.sha(joinpath(paths[k],"files.tsv"))==c["model"]["files_sha256"] || error("Changed $k cohort")
    end
    base=R.F.settings(joinpath(paths["transfer"],"settings.toml")); maps=R.load_maps(opts["--root"],base)
    names=getindex.(t(joinpath(paths["controls"],"files.tsv")),"file")
    if haskey(opts,"--file"); opts["--file"]==first(names) || error("Only the fixed first file is a local witness"); names=[first(names)]; end
    if haskey(opts,"--dry-run")
        statuses=[r for file in names for r in t(joinpath(paths["controls"],splitext(file)[1],"status.tsv"))]
        println("Metadata/maps only: ",length(names)," scans; ",count(r->r["status"]=="ok",statuses)," supported views. No patch prediction or output.")
        return
    end
    out=opts["--outdir"]; R.S.newdir(out); cp(opts["--config"],joinpath(out,"settings.toml"))
    R.table(joinpath(out,"files.tsv"),[(;file) for file in names]); R.table(joinpath(out,"surface_hashes.tsv"),maps.refs)
    BLAS.set_num_threads(1)
    Threads.@threads :static for k in eachindex(names)
        file=names[k]; start=time_ns(); process_file(file,paths,maps,base,c,joinpath(out,splitext(file)[1]))
        progress(stdout,file,(time_ns()-start)/1e9)
    end
    println("FIXED_SHAPE_COMPLETE; saved coefficients only, no fit, labels or promotion")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseFixedShapeReplay.main()
