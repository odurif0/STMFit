#!/usr/bin/env julia
# Read-only arithmetic on saved observations/decisions; no new fitting or grading.
module VerifyVacuumSurfaceTransfer
using Test, TOML, Statistics
include(joinpath(@__DIR__,"diagnose_vacuum_surface_transfer.jl"))
const R=DiagnoseVacuumSurfaceTransfer
const F=R.F
n(r,k)=R.num(r,k)
i(r,k)=R.integer(r,k)
t(path)=R.V.tsv(path)
close(a,b)=isapprox(a,b;rtol=1e-10,atol=1e-18)

"Independent four-weight interpolation of the saved physical height grid."
function at(map,row,state,c)
    co,si=cosd(n(state,"angle")),sind(n(state,"angle"))
    x=n(row,"dx_nm")-n(state,"tx"); y=n(row,"dy_nm")-n(state,"ty")
    half,step=c["model"]["map_half_nm"],c["model"]["map_step_nm"]
    t=1+(co*x+si*y+half)/step; u=1+(-si*x+co*y+half)/step
    ix=min(floor(Int,t),size(map,2)-1); iy=min(floor(Int,u),size(map,1)-1)
    fx,fy=t-ix,u-iy
    sum(map[iy+dy,ix+dx]*(dx==0 ? 1-fx : fx)*(dy==0 ? 1-fy : fy) for dy in 0:1,dx in 0:1)
end

function verify_direction(dir,maps,c)
    models=t(joinpath(dir,"models.tsv")); rows=t(joinpath(dir,"patches.tsv")); obs=t(joinpath(dir,"observations.tsv"))
    cal=TOML.parsefile(joinpath(dir,"calibration.toml"))
    grouped=Dict(j=>filter(r->i(r,"patch")==j,obs) for j in unique(i.(obs,"patch")))
    @test Set((i(r,"interval"),r["arm"]) for r in models)==Set((j,a) for j in 1:3 for a in ("common","chemical"))
    @test length(models)==6
    pixels=[(i(r,"row"),i(r,"column")) for r in obs]
    @test length(unique(pixels))==length(pixels)
    @test all(isfinite(n(r,"source_nm")) for r in obs)
    block=c["preprocessing"]["calibration_block_rows"]; buffer=c["preprocessing"]["holdout_buffer_rows"]
    for r in obs
        if isfinite(n(r,"target_nm"))
            @test isodd(fld(i(r,"row")-1,block))
            @test buffer<=mod(i(r,"row")-1,block)<block-buffer
        end
    end
    for model in models
        interval=i(model,"interval"); arm=model["arm"]; offset=n(model,"offset_nm")
        @test close(n(model,"flat_height_nm"),mean(n.(obs,"source_nm")))
        pp=filter(r->i(r,"interval")==interval && r["arm"]==arm,rows)
        @test Set(i.(pp,"patch"))==Set(keys(grouped)) && length(pp)==length(grouped)
        costs=t(joinpath(dir,"costs_$(interval)_$(arm).tsv"))
        training=0.; testing=0.; flatloss=0.; copyloss=0.; testpixels=0; residualsum=0.
        for patch in pp
            j=i(patch,"patch"); observations=grouped[j]
            @test n(patch,"offset_nm")==offset
            prototype=maps.heights[interval][i(patch,"type")]
            pred=[offset+at(prototype,r,patch,c) for r in observations]
            source=n.(observations,"source_nm"); target=n.(observations,"target_nm"); good=isfinite.(target)
            train=sum(abs2,source-pred); test=sum(abs2,target[good]-pred[good])
            flat=sum(abs2,target[good].-n(model,"flat_height_nm")); copy=sum(abs2,target[good]-source[good])
            @test close(train,n(patch,"training_sse_nm2")) && close(test,n(patch,"test_sse_nm2"))
            @test close(flat,n(patch,"flat_sse_nm2")) && close(copy,n(patch,"source_copy_sse_nm2"))
            @test length(source)==i(patch,"training_pixels") && count(good)==i(patch,"test_pixels")
            cc=filter(r->i(r,"patch")==j,costs); chosen=filter(r->r["chosen"]=="true",cc)
            @test length(cc)==(arm=="common" ? 216 : 432) && length(chosen)==1
            selected=only(chosen)
            @test all(selected[k]==patch[k] for k in ("type","angle","tx","ty"))
            @test all(i(r,"pixels")==length(source) for r in cc)
            mu=mean(source-(pred.-offset)); centered=sum(abs2,source-(pred.-offset).-mu)
            @test close(mu,n(selected,"residual_mean_nm")) && close(centered,n(selected,"centered_sse_nm2"))
            objectives=[n(r,"centered_sse_nm2")+i(r,"pixels")*(offset-n(r,"residual_mean_nm"))^2 for r in cc]
            @test close(minimum(objectives),train)
            residualsum+=sum(source-pred)
            training+=train; testing+=test; flatloss+=flat; copyloss+=copy; testpixels+=count(good)
        end
        @test abs(residualsum)<=1e-10*max(1,length(obs))
        @test close(training,n(model,"training_sse_nm2")) && close(training,n(model,"optimizer_sse_nm2"))
        @test close(testing,n(model,"test_sse_nm2")) && testpixels==i(model,"test_pixels")
        @test close(flatloss,n(model,"flat_sse_nm2")) && close(copyloss,n(model,"source_copy_sse_nm2"))
        @test i(model,"training_pixels")==length(obs)
    end
    for j in 1:3
        common=only(r for r in models if i(r,"interval")==j && r["arm"]=="common")
        chemical=only(r for r in models if i(r,"interval")==j && r["arm"]=="chemical")
        @test all(common[k]==chemical[k] for k in ("test_pixels","training_pixels","flat_height_nm","flat_sse_nm2","source_copy_sse_nm2"))
    end
    (;models,cal,patches=length(grouped))
end

function verify(root,run,rawdir)
    c=F.settings(joinpath(run,"settings.toml")); maps=R.load_maps(root,c)
    files=getindex.(t(joinpath(run,"files.tsv")),"file")
    actual=sort(filter(f->endswith(lowercase(f),".sxm") && isfile(joinpath(rawdir,f)),readdir(rawdir)))
    cases=NamedTuple[]; statuses=Dict{String,Int}(); unsupported=String[]; patches=0
    @testset "Saved surface transfer: cohort, disjoint target rows and physical predictions" begin
        @test files==actual
        references=t(joinpath(run,"surface_hashes.tsv"))
        @test length(references)==length(maps.refs)==6
        @test Set((r["species"],r["kind"],r["sha256"]) for r in references)==
            Set((r.species,r.kind,r.sha256) for r in maps.refs)
        for file in files
            dir=joinpath(run,splitext(file)[1]); meta=TOML.parsefile(joinpath(dir,"input.toml"))
            @test meta["file"]==file && meta["raw_sha256"]==R.sha(joinpath(rawdir,file))
            status=t(joinpath(dir,"status.tsv"))
            if !isapprox(meta["bias_v"],c["model"]["sample_bias_v"];rtol=0,atol=c["model"]["bias_atol_v"])
                @test length(status)==1 && only(status)["status"]=="unsupported_bias"
                push!(unsupported,file); continue
            end
            @test getindex.(status,"direction")==["fwd","bwd"]
            for row in status
                reason=row["status"]; statuses[reason]=get(statuses,reason,0)+1
                view=row["direction"]; path=joinpath(dir,view)
                isfile(joinpath(path,"models.tsv")) || continue
                result=verify_direction(path,maps,c); patches+=result.patches
                @test result.patches==i(row,"patches")
                for interval in 1:3
                    common=only(r for r in result.models if i(r,"interval")==interval && r["arm"]=="common")
                    chemical=only(r for r in result.models if i(r,"interval")==interval && r["arm"]=="chemical")
                    push!(cases,(;file,view,interval,identified=result.cal["identified"],pixels=i(common,"test_pixels"),
                        common=n(common,"test_sse_nm2"),chemical=n(chemical,"test_sse_nm2"),
                        flat=n(common,"flat_sse_nm2"),source_copy=n(common,"source_copy_sse_nm2")))
                end
            end
        end
    end
    println((;raw_files=length(files),unsupported_bias=unsupported,statuses,source_patches=patches))
    for group in ("all_scored","identified_only"), interval in 1:3
        rows=filter(r->r.interval==interval && r.pixels>0 && (group=="all_scored" || r.identified),cases)
        isempty(rows) && (println((;group,interval,cases=0)); continue)
        np=sum(r.pixels for r in rows)
        losses=Dict(key=>sum(getproperty(r,key) for r in rows) for key in (:common,:chemical,:flat,:source_copy))
        wins=count(r->r.chemical<r.common,rows); ties=count(r->r.chemical==r.common,rows)
        perfile=Dict(f=>sum(r.chemical-r.common for r in rows if r.file==f) for f in unique(getproperty.(rows,:file)))
        println((;group,interval,views=length(rows),files=length(perfile),pixels=np,
            rms_pm=Dict(k=>1000sqrt(v/np) for (k,v) in losses),
            chemical_relative_gain=1-losses[:chemical]/losses[:common],
            view_wins=wins,view_ties=ties,view_losses=length(rows)-wins-ties,
            file_wins=count(<(0),values(perfile)),file_losses=count(>(0),values(perfile))))
    end
    (;cases,statuses,unsupported,patches)
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    ARGS==["--help"] ? println("verify_vacuum_surface_transfer.jl ROOT SAVED_RUN RAW_DIR (read-only; no new fitting)") :
        (length(ARGS)==3 ? VerifyVacuumSurfaceTransfer.verify(ARGS...) : error("Expected root, run and raw directory"))
end
