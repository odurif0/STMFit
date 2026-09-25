#!/usr/bin/env julia
# Independent interpolation and saved-coefficient checks; no fitting or grading.
module VerifyFixedShapeReplay
using Test, TOML, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_fixed_shape_replay.jl"))
include(joinpath(@__DIR__,"verify_vacuum_surface_transfer.jl"))
include(joinpath(@__DIR__,"verify_image_foreground.jl"))
const P=DiagnoseFixedShapeReplay
const V=VerifyVacuumSurfaceTransfer
const G=VerifyImageForeground
const R=P.R
const t=P.t
const n=P.n
const i=P.i
b(r,k)=parse(Bool,r[k])

function verify_direction(dir,transfer,controls,foreground,maps,base,c)
    dirs=Dict("transfer"=>transfer,"controls"=>controls,"foreground"=>foreground)
    hashes=t(joinpath(dir,"input_hashes.tsv"))
    @test length(hashes)==8 && length(unique((r["role"],r["file"]) for r in hashes))==8
    @test all(R.sha(joinpath(dirs[r["role"]],r["file"]))==r["sha256"] for r in hashes)
    @test read(joinpath(dir,"calibration.toml"))==read(joinpath(transfer,"calibration.toml"))
    cal=TOML.parsefile(joinpath(dir,"calibration.toml"))
    observations=V.by_patch(t(joinpath(transfer,"observations.tsv"))); rows=t(joinpath(dir,"patches.tsv"))
    old=t(joinpath(controls,"patches.tsv")); oldby=Dict((i(r,"interval"),r["arm"],i(r,"patch"))=>r for r in old)
    by=Dict((i(r,"interval"),r["arm"],i(r,"patch"))=>r for r in rows)
    anchors=Dict(i(r,"anchor")=>r for r in t(joinpath(foreground,"anchor_membership.tsv")))
    @test length(by)==length(rows)==15length(observations)
    @test Set(keys(by))==Set((interval,arm,j) for interval in 1:3 for arm in c["model"]["arms"] for j in keys(observations))
    for interval in 1:3
        costs=V.by_patch(t(joinpath(controls,"costs_$(interval)_plane_chemical.tsv")))
        @test Set(keys(costs))==Set(keys(observations))
        for j in sort(collect(keys(observations)))
            pixels=observations[j]; source=n.(pixels,"source_nm"); target=n.(pixels,"target_nm"); good=isfinite.(target)
            X=hcat(ones(length(pixels)),n.(pixels,"dx_nm"),n.(pixels,"dy_nm"))
            common=oldby[(interval,"plane_common",j)]; selected=oldby[(interval,"plane_chemical",j)]
            anchor=anchors[i(common,"anchor")]; candidates=sort(costs[j];by=r->i(r,"type"))
            @test i.(candidates,"type")==[0,1] && i.(candidates,"state")==[1,2]
            @test i(selected,"type")+1==argmin(n.(candidates,"training_sse_nm2"))
            @test count(r->b(r,"chosen"),candidates)==1 && b(candidates[i(selected,"type")+1],"chosen")
            for arm in c["model"]["arms"]
                row=by[(interval,arm,j)]; type=i(row,"type")
                @test b(row,"identified")==cal["identified"]
                @test i(row,"anchor")==i(common,"anchor")==i(first(pixels),"anchor")
                @test all(n(row,k)==n(common,k) for k in ("angle","tx","ty"))
                for key in ("initial_available","initial_core","joint_available","joint_core")
                    @test row[key]==anchor[key]
                end
                reference=if arm=="fixed_glcn"
                    @test type==0; candidates[1]
                elseif arm=="fixed_glcnac"
                    @test type==1; candidates[2]
                else
                    oldby[(interval,arm=="source_selected" ? "plane_chemical" : arm,j)]
                end
                @test P.beta(row)==P.beta(reference)
                @test type==i(reference,"type")
                physical=arm=="local_plane" ? zeros(length(pixels)) :
                    [V.at(maps.heights[interval][type],pixel,row,base) for pixel in pixels]
                pred=physical+X*P.beta(row); train=sum(abs2,source-pred); test=sum(abs2,target[good]-pred[good])
                @test P.close(train,n(row,"training_sse_nm2"),c) && P.close(test,n(row,"test_sse_nm2"),c)
                @test i(row,"training_pixels")==length(source) && i(row,"test_pixels")==count(good)
                @test norm(X'*(source-pred),Inf)<=c["preprocessing"]["replay_rtol"]*length(source)
                @test P.close(train,n(reference,"training_sse_nm2"),c)
                if arm in ("source_selected","local_plane","plane_common")
                    @test P.close(test,n(reference,"test_sse_nm2"),c)
                end
            end
            chosen=by[(interval,"source_selected",j)]
            fixed=by[(interval,i(chosen,"type")==0 ? "fixed_glcn" : "fixed_glcnac",j)]
            @test all(chosen[k]==fixed[k] for k in keys(chosen) if k!="arm")
        end
    end
    (;rows,cal,patches=length(observations))
end

function case_rows(rows,file,view,c)
    result=NamedTuple[]
    for interval in 1:3,arm in c["model"]["arms"],(mask,prefix) in (("initial","initial"),("joint_background","joint"))
        rr=filter(r->i(r,"interval")==interval && r["arm"]==arm,rows)
        groups=Dict(g=>eltype(rr)[] for g in ("retained","rejected","unavailable"))
        for r in rr
            group=!b(r,prefix*"_available") ? "unavailable" : b(r,prefix*"_core") ? "retained" : "rejected"
            push!(groups[group],r)
        end
        @test sum(length,values(groups))==length(rr)
        groups["all"]=rr
        for group in c["selection"]["groups"]
            pp=groups[group]; pixels=sum((i(r,"test_pixels") for r in pp);init=0); sse=sum((n(r,"test_sse_nm2") for r in pp);init=0.)
            push!(result,(;file,view,identified=b(first(rr),"identified"),interval,mask,group,arm,patches=length(pp),
                test_pixels=pixels,test_sse_nm2=sse,rms_pm=pixels>0 ? 1000sqrt(sse/pixels) : NaN))
        end
    end
    result
end

function comparisons(cases)
    result=NamedTuple[]
    for subset in ("all_scored","identified"),mask in ("initial","joint_background"),interval in 1:3,
        group in ("all","retained","rejected","unavailable"),reference in ("fixed_glcn","fixed_glcnac","local_plane","plane_common")
        rr=filter(r->r.mask==mask && r.interval==interval && r.group==group && (subset=="all_scored" || r.identified),cases)
        aa=Dict((r.file,r.view)=>r for r in rr if r.arm=="source_selected")
        bb=Dict((r.file,r.view)=>r for r in rr if r.arm==reference)
        @test Set(keys(aa))==Set(keys(bb))
        @test all(aa[k].test_pixels==bb[k].test_pixels && aa[k].patches==bb[k].patches for k in keys(aa))
        kk=sort([k for k in keys(aa) if aa[k].test_pixels>0]); differences=[aa[k].test_sse_nm2-bb[k].test_sse_nm2 for k in kk]
        byfile=Dict{String,Float64}()
        for (k,delta) in zip(kk,differences); byfile[k[1]]=get(byfile,k[1],0.)+delta; end
        reference_sse=sum((bb[k].test_sse_nm2 for k in kk);init=0.); delta=sum(differences;init=0.)
        push!(result,(;subset,mask,interval,group,reference,views=length(kk),files=length(byfile),
            selector_better_views=count(<(0),differences),selector_worse_views=count(>(0),differences),tied_views=count(==(0),differences),
            selector_better_files=count(<(0),values(byfile)),selector_worse_files=count(>(0),values(byfile)),
            test_pixels=sum((aa[k].test_pixels for k in kk);init=0),selector_minus_reference_sse_nm2=delta,
            selector_mse_gain=reference_sse>0 ? -delta/reference_sse : NaN))
    end
    result
end

function main(args=ARGS)
    length(args)==3 || error("verify_fixed_shape_replay.jl REPO SAVED_RUN NEW_REPORT_DIR")
    root,run,out=args; c=P.settings(joinpath(run,"settings.toml")); R.S.newdir(out)
    paths=Dict(k=>joinpath(root,c["model"][k*"_run"]) for k in ("transfer","controls","foreground"))
    base=R.F.settings(joinpath(paths["transfer"],"settings.toml")); maps=R.load_maps(root,base)
    files=getindex.(t(joinpath(run,"files.tsv")),"file"); expected=getindex.(t(joinpath(paths["controls"],"files.tsv")),"file")
    cases=NamedTuple[]; statuses=NamedTuple[]; detail=NamedTuple[]
    @testset "Independent fixed-shape replay and unchanged populations" begin
        @test files==expected || files==[first(expected)]
        for role in keys(paths)
            @test R.sha(joinpath(paths[role],"settings.toml"))==c["model"][role*"_settings_sha256"]
            @test R.sha(joinpath(paths[role],"files.tsv"))==c["model"]["files_sha256"]
        end
        @test Set((r["species"],r["kind"],r["sha256"]) for r in t(joinpath(run,"surface_hashes.tsv")))==Set((r.species,r.kind,r.sha256) for r in maps.refs)
        for file in files
            stem=splitext(file)[1]; dir=joinpath(run,stem); control=joinpath(paths["controls"],stem)
            @test read(joinpath(dir,"upstream_input.toml"))==read(joinpath(paths["transfer"],stem,"input.toml"))
            @test read(joinpath(dir,"upstream_status.tsv"))==read(joinpath(control,"status.tsv"))
            previous=t(joinpath(control,"status.tsv")); now=t(joinpath(dir,"status.tsv"))
            @test getindex.(previous,"direction")==getindex.(now,"direction")
            for (old,row) in zip(previous,now)
                view=row["direction"]; reason=row["status"]
                push!(statuses,(;file,view,status=reason,patches=i(row,"patches")))
                @test i(row,"patches")==i(old,"patches")
                if old["status"]!="ok"; @test reason==old["status"]; continue; end
                @test reason=="ok"
                reason=="ok" || continue
                got=verify_direction(joinpath(dir,view),joinpath(paths["transfer"],stem,view),joinpath(control,view),
                    joinpath(paths["foreground"],stem,view),maps,base,c)
                @test got.patches==i(row,"patches")
                append!(cases,case_rows(got.rows,file,view,c))
                append!(detail,[(;file,view,(Symbol(k)=>v for (k,v) in r)...) for r in got.rows])
            end
            println("FIXED_SHAPE_VERIFIED ",file); flush(stdout)
        end
        compared=comparisons(cases)
        R.table(joinpath(out,"comparisons.tsv"),compared)
    end
    R.table(joinpath(out,"statuses.tsv"),statuses); R.table(joinpath(out,"cases.tsv"),cases)
    R.table(joinpath(out,"summary.tsv"),G.aggregate_predictions(cases)); R.table(joinpath(out,"patches.tsv"),detail)
    println("Verified ",length(files)," scans; fixed-shape counterfactuals only, no chemical validation or grade.")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyFixedShapeReplay.main()
