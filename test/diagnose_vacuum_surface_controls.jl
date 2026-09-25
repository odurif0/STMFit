#!/usr/bin/env julia
# Source-only control fits on frozen native observations; never reads benchmark labels.
module DiagnoseVacuumSurfaceControls
using LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__, "diagnose_vacuum_surface_transfer.jl"))
include(joinpath(@__DIR__, "lib/vacuum_surface_controls.jl"))
const R = DiagnoseVacuumSurfaceTransfer
const C = VacuumSurfaceControls
n(r,k) = R.num(r,k)
i(r,k) = R.integer(r,k)
t(path) = R.V.tsv(path)

function read_observations(path)
    rows = t(path); groups = Dict{Int,Vector{eltype(rows)}}()
    for row in rows
        push!(get!(groups, i(row,"patch"), eltype(rows)[]), row)
    end
    ids = sort(collect(keys(groups)))
    ids == collect(1:length(ids)) || error("Noncontiguous patches")
    patches = NamedTuple[]; targets = Vector{Float64}[]
    for j in ids
        rr = groups[j]
        i.(rr,"pixel") == collect(1:length(rr)) || error("Changed pixel order")
        length(unique(i.(rr,"anchor"))) == 1 || error("Mixed anchor")
        push!(patches, (; anchor=i(first(rr),"anchor"), dx=n.(rr,"dx_nm"), dy=n.(rr,"dy_nm"), values=n.(rr,"source_nm")))
        push!(targets, n.(rr,"target_nm"))
    end
    (; patches, targets, rows)
end

function one_direction(dir, maps, base, c, out)
    R.S.newdir(out)
    consumed = ["observations.tsv", "models.tsv", "calibration.toml"]
    hashes = [(; file=f, sha256=R.sha(joinpath(dir,f))) for f in consumed]
    obs = read_observations(joinpath(dir,"observations.tsv"))
    old = t(joinpath(dir,"models.tsv")); cal = TOML.parsefile(joinpath(dir,"calibration.toml"))
    allrows = NamedTuple[]; models = NamedTuple[]; replayrows = NamedTuple[]
    for interval in base["model"]["intervals"]
        result = C.fit_patches(obs.patches, maps.heights[interval], base, c)
        reference = only(r for r in old if i(r,"interval")==interval && r["arm"]=="common")
        p = c["preprocessing"]
        isapprox(result.replay.offset,n(reference,"offset_nm");rtol=p["replay_rtol"],atol=p["replay_atol_nm"]) &&
            isapprox(result.replay.loss,n(reference,"training_sse_nm2");rtol=p["replay_rtol"],atol=p["replay_atol_nm2"]) || error("Common source replay changed")
        push!(replayrows,(; interval, offset_nm=result.replay.offset, training_sse_nm2=result.replay.loss,
            mean_contrast_nm=result.delta, offset_segments=result.height.segments))
        R.table(joinpath(out,"costs_$(interval)_height_contrast.tsv"),[(;patch=j,state=k,s...,pixels=cost.n,
            residual_mean_nm=cost.mu[k],centered_sse_nm2=cost.v[k],chosen=k==result.height.choices[j])
            for (j,cost) in enumerate(result.height_costs) for (k,s) in enumerate(result.hs)])
        for arm in ("plane_common","plane_chemical")
            costs = arm=="plane_common" ? result.plane_costs : result.chemical_costs
            R.table(joinpath(out,"costs_$(interval)_$(arm).tsv"),[(;patch=j,state=k,
                type=arm=="plane_common" ? -1 : k-1,
                training_sse_nm2=cost.losses[k],beta0_nm=cost.betas[k][1],betax=cost.betas[k][2],betay=cost.betas[k][3],
                chosen=k==result.fits[arm][j].choice)
                for (j,cost) in enumerate(costs) for k in eachindex(cost.losses)])
        end
        for arm in c["model"]["arms"]
            train=0.; test=0.; np=0
            for (j,(patch,target,fit)) in enumerate(zip(obs.patches,obs.targets,result.fits[arm]))
                good=isfinite.(target); countgood=count(good)
                training=sum(abs2,patch.values-fit.pred); testing=sum(abs2,target[good]-fit.pred[good])
                push!(allrows,(;interval,arm,patch=j,anchor=patch.anchor,fit.state...,choice=fit.choice,
                    beta0_nm=fit.beta[1],betax=fit.beta[2],betay=fit.beta[3],
                    training_pixels=length(patch.values),training_sse_nm2=training,
                    test_pixels=countgood,test_sse_nm2=testing))
                train+=training; test+=testing; np+=countgood
            end
            np==i(reference,"test_pixels") || error("Target support changed")
            push!(models,(;interval,isovalue=maps.isos[interval],arm,identified=cal["identified"],
                mean_contrast_nm=result.delta,training_pixels=length(obs.rows),training_sse_nm2=train,
                test_pixels=np,test_sse_nm2=test))
        end
    end
    all(R.sha(joinpath(dir,r.file))==r.sha256 for r in hashes) || error("Input changed")
    R.table(joinpath(out,"input_hashes.tsv"),hashes)
    R.table(joinpath(out,"models.tsv"),models); R.table(joinpath(out,"patches.tsv"),allrows)
    R.table(joinpath(out,"common_replay.tsv"),replayrows)
    (;patches=length(obs.patches),status="ok")
end

function process_file(file, input, maps, base, c, out)
    R.S.newdir(out); dir=joinpath(input,splitext(file)[1]); statuses=NamedTuple[]
    for name in ("input.toml","status.tsv"); cp(joinpath(dir,name),joinpath(out,"upstream_"*name)); end
    for row in t(joinpath(dir,"status.tsv"))
        direction=row["direction"]; reason=row["status"]
        if reason != "ok"
            push!(statuses,(;direction,status="upstream:"*reason,patches=i(row,"patches"))); continue
        end
        try
            result=one_direction(joinpath(dir,direction),maps,base,c,joinpath(out,direction))
            result.patches==i(row,"patches") || error("Lost source patches")
            push!(statuses,(;direction,result...))
        catch e
            e isa InterruptException && rethrow()
            push!(statuses,(;direction,status="failed:"*replace(sprint(showerror,e),'\n'=>' ','\t'=>' '),patches=i(row,"patches")))
        end
    end
    R.table(joinpath(out,"status.tsv"),statuses)
end

function main(args=ARGS)
    "--help" in args && return println("diagnose_vacuum_surface_controls.jl --root REPO --config TOML --outdir NEW_DIR [--dry-run]")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    opts=Dict{String,String}(); j=1
    while j<=length(args)
        k=args[j]; haskey(opts,k) && error("Duplicate option")
        if k=="--dry-run"; opts[k]="true"; j+=1; continue; end
        k in ("--root","--config","--outdir") && j<length(args) || error("Unexpected option")
        opts[k]=args[j+1]; j+=2
    end
    all(haskey(opts,k) for k in ("--root","--config","--outdir")) || error("Missing option")
    c=C.settings(opts["--config"]); input=joinpath(opts["--root"],c["model"]["input_run"])
    R.sha(joinpath(input,"settings.toml"))==c["model"]["input_settings_sha256"] || error("Changed upstream settings")
    R.sha(joinpath(input,"files.tsv"))==c["model"]["input_files_sha256"] || error("Changed cohort")
    base=R.F.settings(joinpath(input,"settings.toml")); maps=R.load_maps(opts["--root"],base)
    names=getindex.(t(joinpath(input,"files.tsv")),"file")
    names==sort(unique(names)) && all(f->basename(f)==f && endswith(f,".sxm"),names) || error("Invalid cohort")
    statuses=[r for f in names for r in t(joinpath(input,splitext(f)[1],"status.tsv"))]
    if haskey(opts,"--dry-run")
        println("Metadata-only: ",length(names)," upstream scans, ",count(r->r["status"]=="ok",statuses)," usable views; all failures retained. No pixel fit/output.")
        return
    end
    out=opts["--outdir"]; R.S.newdir(out)
    cp(opts["--config"],joinpath(out,"settings.toml")); cp(joinpath(input,"files.tsv"),joinpath(out,"files.tsv"))
    cp(joinpath(input,"settings.toml"),joinpath(out,"upstream_settings.toml"))
    R.table(joinpath(out,"surface_hashes.tsv"),maps.refs); BLAS.set_num_threads(1)
    Threads.@threads :static for k in eachindex(names)
        file=names[k]; started=time_ns()
        process_file(file,input,maps,base,c,joinpath(out,splitext(file)[1]))
        println("CONTROLS_FILE_COMPLETE ",file," elapsed_s=",(time_ns()-started)/1e9); flush(stdout)
    end
    println("SURFACE_CONTROLS_COMPLETE; no benchmark, production change or chemical validation")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseVacuumSurfaceControls.main()
