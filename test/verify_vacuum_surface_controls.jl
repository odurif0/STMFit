#!/usr/bin/env julia
# Read-only saved prediction and linear normal-equation checks, not a new model fit.
module VerifyVacuumSurfaceControls
using Test, TOML, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "diagnose_vacuum_surface_controls.jl"))
include(joinpath(@__DIR__, "verify_vacuum_surface_transfer.jl"))
const L = DiagnoseVacuumSurfaceControls
const V = VerifyVacuumSurfaceTransfer
const R = L.R
n(r,k)=L.n(r,k)
i(r,k)=L.i(r,k)
t(path)=L.t(path)
close(a,b)=isapprox(a,b;rtol=1e-10,atol=1e-18)
beta(r)=[n(r,"beta0_nm"),n(r,"betax"),n(r,"betay")]

function prediction(row, patch, arm, maps, base, delta)
    physical = if arm=="height_contrast"
        V.at(maps[-1],row,patch,base)+(2i(patch,"type")-1)*delta
    elseif arm=="plane_common"
        V.at(maps[-1],row,patch,base)
    elseif arm=="plane_chemical"
        V.at(maps[i(patch,"type")],row,patch,base)
    elseif arm in ("local_constant","local_plane")
        0.
    else
        error("Unknown control arm")
    end
    physical+n(patch,"beta0_nm")+n(patch,"betax")*n(row,"dx_nm")+n(patch,"betay")*n(row,"dy_nm")
end

function verify_direction(dir, input, maps, base, c)
    obs=t(joinpath(input,"observations.tsv")); grouped=V.by_patch(obs)
    models=t(joinpath(dir,"models.tsv")); patches=t(joinpath(dir,"patches.tsv"))
    modelby=Dict((i(r,"interval"),r["arm"])=>r for r in models)
    patchby=Dict((i(r,"interval"),r["arm"],i(r,"patch"))=>r for r in patches)
    old=t(joinpath(input,"models.tsv")); cal=TOML.parsefile(joinpath(input,"calibration.toml"))
    hashes=t(joinpath(dir,"input_hashes.tsv"))
    @test Set(r["file"] for r in hashes)==Set(["observations.tsv","models.tsv","calibration.toml"])
    @test length(hashes)==3 && all(R.sha(joinpath(input,r["file"]))==r["sha256"] for r in hashes)
    @test length(models)==length(modelby)==15
    @test Set(keys(modelby))==Set((j,a) for j in 1:3 for a in c["model"]["arms"])
    @test length(patches)==length(patchby)==15length(grouped)
    @test Set(keys(patchby))==Set((j,a,k) for j in 1:3 for a in c["model"]["arms"] for k in keys(grouped))
    replay=t(joinpath(dir,"common_replay.tsv"))
    @test i.(replay,"interval")==[1,2,3]
    ss=R.F.states(base,"common"); hs=R.F.states(base,"chemical")
    for interval in 1:3
        map=maps.heights[interval]; delta=(mean(map[1])-mean(map[0]))/2
        re=only(r for r in replay if i(r,"interval")==interval)
        ref=only(r for r in old if i(r,"interval")==interval && r["arm"]=="common")
        @test close(n(re,"offset_nm"),n(ref,"offset_nm")) && close(n(re,"training_sse_nm2"),n(ref,"training_sse_nm2"))
        @test isapprox(delta,n(re,"mean_contrast_nm");rtol=1e-12,atol=1e-14)
        costs=Dict(a=>V.by_patch(t(joinpath(dir,"costs_$(interval)_$(a).tsv"))) for a in ("height_contrast","plane_common","plane_chemical"))
        for groups in values(costs); @test Set(keys(groups))==Set(keys(grouped)); end
        for arm in c["model"]["arms"]
            model=modelby[(interval,arm)]
            @test parse(Bool,model["identified"])==cal["identified"]
            @test isapprox(n(model,"mean_contrast_nm"),delta;rtol=1e-12,atol=1e-14)
            train=0.; test=0.; np=0; residualsum=0.; offsets=Float64[]
            for (j,rows) in grouped
                patch=patchby[(interval,arm,j)]
                @test i(patch,"anchor")==i(first(rows),"anchor")
                source=n.(rows,"source_nm"); target=n.(rows,"target_nm"); good=isfinite.(target)
                X=hcat(ones(length(rows)),n.(rows,"dx_nm"),n.(rows,"dy_nm"))
                pred=[prediction(r,patch,arm,map,base,delta) for r in rows]
                @test all(isfinite,pred)
                training=sum(abs2,source-pred); testing=sum(abs2,target[good]-pred[good])
                @test close(training,n(patch,"training_sse_nm2")) && close(testing,n(patch,"test_sse_nm2"))
                @test i(patch,"training_pixels")==length(rows) && i(patch,"test_pixels")==count(good)
                if arm=="height_contrast"
                    @test n(patch,"betax")==n(patch,"betay")==0
                    push!(offsets,n(patch,"beta0_nm")); residualsum+=sum(source-pred)
                else
                    @test abs(sum(source-pred))<1e-10length(rows)
                    if arm=="local_constant"
                        @test n(patch,"betax")==n(patch,"betay")==0
                        @test close(n(patch,"beta0_nm"),mean(source))
                    else
                        @test norm(X'*(source-pred),Inf)<1e-10length(rows)
                    end
                end
                if arm in ("height_contrast","plane_common","plane_chemical")
                    cc=costs[arm][j]; chosen=only(r for r in cc if r["chosen"]=="true")
                    @test i(chosen,"state")==i(patch,"choice")
                    if arm=="height_contrast"
                        @test length(cc)==length(hs)==432
                        state=hs[i(patch,"choice")]
                        @test all(n(patch,string(k))==getproperty(state,k) for k in keys(state))
                        @test all(i(r,"pixels")==length(rows) for r in cc)
                        offset=n(patch,"beta0_nm")
                        objectives=[n(r,"centered_sse_nm2")+length(rows)*(offset-n(r,"residual_mean_nm"))^2 for r in cc]
                        @test close(minimum(objectives),training)
                        residue=source-(pred.-offset)
                        @test close(n(chosen,"residual_mean_nm"),mean(residue))
                        @test close(n(chosen,"centered_sse_nm2"),sum(abs2,residue.-mean(residue)))
                    else
                        @test length(cc)==(arm=="plane_common" ? 216 : 2)
                        @test beta(chosen)==beta(patch)
                        @test close(n(chosen,"training_sse_nm2"),training)
                        @test close(minimum(n.(cc,"training_sse_nm2")),training)
                        if arm=="plane_common"
                            state=ss[i(patch,"choice")]
                            @test all(n(patch,string(k))==getproperty(state,k) for k in keys(state))
                        else
                            common=patchby[(interval,"plane_common",j)]
                            @test all(patch[k]==common[k] for k in ("angle","tx","ty"))
                            @test i(patch,"type")==i(patch,"choice")-1
                            # Verify both chemical alternatives at exactly the frozen geometry.
                            for candidate in cc
                                physical=[V.at(map[i(candidate,"type")],r,common,base) for r in rows]
                                residue=source-physical-X*beta(candidate)
                                @test norm(X'*residue,Inf)<1e-10length(rows)
                                @test close(sum(abs2,residue),n(candidate,"training_sse_nm2"))
                            end
                        end
                    end
                end
                train+=training; test+=testing; np+=count(good)
            end
            @test i(model,"training_pixels")==length(obs) && i(model,"test_pixels")==np==i(ref,"test_pixels")
            @test close(n(model,"training_sse_nm2"),train) && close(n(model,"test_sse_nm2"),test)
            if arm=="height_contrast"
                @test length(unique(offsets))==1
                @test abs(residualsum)<1e-10length(obs)
            end
        end
    end
    for arm in ("local_constant","local_plane"), j in keys(grouped), interval in (2,3)
        a=patchby[(1,arm,j)]; b=patchby[(interval,arm,j)]
        @test all(a[k]==b[k] for k in keys(a) if k!="interval")
    end
    (;models,cal,patches=length(grouped))
end

function verify(root,run)
    c=L.C.settings(joinpath(run,"settings.toml")); input=joinpath(root,c["model"]["input_run"])
    base=R.F.settings(joinpath(input,"settings.toml")); maps=R.load_maps(root,base)
    names=getindex.(t(joinpath(input,"files.tsv")),"file")
    cases=NamedTuple[]; statuses=Dict{String,Int}(); sourcepatches=0
    @testset "Saved surface controls: accounting, prediction, geometry and nuisance" begin
        @test R.sha(joinpath(input,"settings.toml"))==c["model"]["input_settings_sha256"]
        @test R.sha(joinpath(input,"files.tsv"))==c["model"]["input_files_sha256"]
        @test read(joinpath(run,"files.tsv"))==read(joinpath(input,"files.tsv"))
        @test read(joinpath(run,"upstream_settings.toml"))==read(joinpath(input,"settings.toml"))
        refs=t(joinpath(run,"surface_hashes.tsv"))
        @test length(refs)==6 && Set((r["species"],r["kind"],r["sha256"]) for r in refs)==Set((r.species,r.kind,r.sha256) for r in maps.refs)
        for file in names
            dir=joinpath(run,splitext(file)[1]); upstream=joinpath(input,splitext(file)[1])
            for name in ("input.toml","status.tsv")
                @test read(joinpath(dir,"upstream_"*name))==read(joinpath(upstream,name))
            end
            old=t(joinpath(upstream,"status.tsv")); now=t(joinpath(dir,"status.tsv"))
            @test getindex.(now,"direction")==getindex.(old,"direction")
            for (before,after) in zip(old,now)
                reason=after["status"]; view=after["direction"]; statuses[reason]=get(statuses,reason,0)+1
                @test i(before,"patches")==i(after,"patches")
                if before["status"]!="ok"
                    @test reason=="upstream:"*before["status"]; continue
                end
                reason=="ok" || (@test startswith(reason,"failed:"); continue)
                got=verify_direction(joinpath(dir,view),joinpath(upstream,view),maps,base,c)
                @test got.patches==i(after,"patches"); sourcepatches+=got.patches
                bm=t(joinpath(upstream,view,"models.tsv"))
                for interval in 1:3
                    common=only(r for r in bm if i(r,"interval")==interval && r["arm"]=="common")
                    chemical=only(r for r in bm if i(r,"interval")==interval && r["arm"]=="chemical")
                    arms=Dict(r["arm"]=>n(r,"test_sse_nm2") for r in got.models if i(r,"interval")==interval)
                    push!(cases,(;file,view,interval,identified=got.cal["identified"],pixels=i(common,"test_pixels"),
                        common=n(common,"test_sse_nm2"),chemical=n(chemical,"test_sse_nm2"),source_copy=n(common,"source_copy_sse_nm2"),
                        height_contrast=arms["height_contrast"],local_constant=arms["local_constant"],local_plane=arms["local_plane"],
                        plane_common=arms["plane_common"],plane_chemical=arms["plane_chemical"]))
                end
            end
            println("CONTROLS_CHECKED ",file); flush(stdout)
        end
    end
    println((;files=length(names),statuses,sourcepatches))
    for group in ("all_scored","identified_only"), interval in 1:3
        rows=filter(r->r.interval==interval && r.pixels>0 && (group=="all_scored" || r.identified),cases)
        isempty(rows) && continue
        np=sum(r.pixels for r in rows)
        losses=Dict(a=>sum(getproperty(r,a) for r in rows) for a in (:common,:chemical,:source_copy,:height_contrast,:local_constant,:local_plane,:plane_common,:plane_chemical))
        println((;group,interval,views=length(rows),files=length(unique(r.file for r in rows)),pixels=np,
            rms_pm=Dict(a=>1000sqrt(v/np) for (a,v) in losses),
            chemical_gain_vs_height=1-losses[:chemical]/losses[:height_contrast],
            plane_chemical_gain=1-losses[:plane_chemical]/losses[:plane_common],
            plane_view_wins=count(r->r.plane_chemical<r.plane_common,rows)))
    end
    (;cases,statuses,sourcepatches)
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    ARGS==["--help"] ? println("verify_vacuum_surface_controls.jl ROOT SAVED_RUN (read-only)") :
        (length(ARGS)==2 ? VerifyVacuumSurfaceControls.verify(ARGS...) : error("Expected root and saved run"))
end
