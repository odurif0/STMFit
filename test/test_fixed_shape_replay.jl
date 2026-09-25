# Synthetic only, including upstream control regression. No benchmark inputs.
include(joinpath(@__DIR__,"test_vacuum_surface_controls.jl"))
include(joinpath(@__DIR__,"diagnose_fixed_shape_replay.jl"))
include(joinpath(@__DIR__,"verify_fixed_shape_replay.jl"))
const P=DiagnoseFixedShapeReplay
const W=VerifyFixedShapeReplay
const REPLAY_PATH=joinpath(@__DIR__,"../config/fixed_shape_replay.toml")
const REPLAY=P.settings(REPLAY_PATH)
strings(row)=Dict(string(k)=>string(getproperty(row,k)) for k in keys(row))

@testset "Concurrent progress records do not interleave" begin
    io=IOBuffer()
    @sync for k in 1:200
        Threads.@spawn P.progress(io,"synthetic_$(k).sxm",Float64(k))
    end
    lines=split(chomp(String(take!(io))),'\n')
    @test length(lines)==200
    @test Set(lines)==Set("FIXED_SHAPE_FILE_COMPLETE synthetic_$(k).sxm elapsed_s=$(Float64(k))" for k in 1:200)
end

function old_fixture(f,j;target=f.pp[j].values.+.002)
    result=U.fit_patches(f.pp,f.maps,f.small,CONTROL); p=f.pp[j]
    function oldrow(arm)
        fit=result.fits[arm][j]; good=isfinite.(target)
        strings((;fit.state...,choice=fit.choice,beta0_nm=fit.beta[1],betax=fit.beta[2],betay=fit.beta[3],
            training_pixels=length(p.values),test_pixels=count(good),training_sse_nm2=sum(abs2,p.values-fit.pred),
            test_sse_nm2=sum(abs2,target[good]-fit.pred[good])))
    end
    costs=result.chemical_costs[j]; chosen=result.fits["plane_chemical"][j].choice
    candidates=[strings((;type=k-1,state=k,training_sse_nm2=costs.losses[k],
        beta0_nm=costs.betas[k][1],betax=costs.betas[k][2],betay=costs.betas[k][3],chosen=k==chosen)) for k in 1:2]
    geometry=result.fits["plane_common"][j].state
    (;p,target,candidates,geometry,selected=oldrow("plane_chemical"),plain=oldrow("local_plane"),common=oldrow("plane_common"))
end

function evaluate(f,old)
    P.replay_patch(old.p,old.target,f.maps,old.geometry,old.candidates,old.selected,old.plain,old.common,f.small,REPLAY)
end

@testset "Both fixed shapes, saved selector and source/target boundary" begin
    f=control_fixture()
    for j in eachindex(f.pp)
        old=old_fixture(f,j); got=evaluate(f,old)
        @test getproperty.(got.rows,:arm)==REPLAY["model"]["arms"]
        choice=j==1 ? "fixed_glcn" : "fixed_glcnac"
        @test got.predictions["source_selected"]===got.predictions[choice]
        @test got.predictions[choice]≈old.p.values
        wrong=choice=="fixed_glcn" ? "fixed_glcnac" : "fixed_glcn"
        @test norm(got.predictions[wrong]-old.p.values)>1e-4
        for row in got.rows
            @test row.training_pixels==row.test_pixels==length(old.p.values)
            @test row.test_sse_nm2≈sum(abs2,old.target-got.predictions[row.arm])
        end
        altered=copy(old.target); altered[1:3].=NaN; altered[4:end].+=8.
        changed=old_fixture(f,j;target=altered); other=evaluate(f,changed)
        @test isequal(got.predictions,other.predictions)
        for (a,b) in zip(got.rows,other.rows)
            @test b.test_pixels==a.test_pixels-3 && b.test_sse_nm2>a.test_sse_nm2
            @test a.training_sse_nm2==b.training_sse_nm2
            @test (a.type,a.angle,a.tx,a.ty,a.beta0_nm,a.betax,a.betay)==(b.type,b.angle,b.tx,b.ty,b.beta0_nm,b.betax,b.betay)
        end
        no_target=old_fixture(f,j;target=fill(NaN,length(old.target))); none=evaluate(f,no_target)
        @test all(r.test_pixels==0 && r.test_sse_nm2==0 for r in none.rows)
        @test isequal(got.predictions,none.predictions)
        bad=deepcopy(old); bad.candidates[1]["beta0_nm"]=string(P.n(bad.candidates[1],"beta0_nm")+.1)
        @test_throws ErrorException evaluate(f,bad)
        bad=deepcopy(old); bad.selected["tx"]="0.1"
        @test_throws ErrorException evaluate(f,bad)
        @test_throws ErrorException P.replay_patch(old.p,old.target,f.maps,old.geometry,old.candidates[1:1],old.selected,old.plain,old.common,f.small,REPLAY)
    end
    old=old_fixture(f,1); original=evaluate(f,old)
    swapped=merge(f,(maps=Dict(-1=>f.maps[-1],0=>f.maps[1],1=>f.maps[0]),))
    reverse=evaluate(swapped,old_fixture(swapped,1))
    @test reverse.predictions["source_selected"]≈original.predictions["source_selected"]
    @test reverse.predictions["fixed_glcn"]≈original.predictions["fixed_glcnac"]
    @test reverse.predictions["fixed_glcnac"]≈original.predictions["fixed_glcn"]
    null=merge(f,(maps=Dict(-1=>f.maps[-1],0=>f.maps[-1],1=>f.maps[-1]),))
    nullresult=evaluate(null,old_fixture(null,1))
    @test nullresult.predictions["fixed_glcn"]==nullresult.predictions["fixed_glcnac"]
end

@testset "Saved synthetic workflow and mask-only reporting" begin
    f=control_fixture(); maps=(heights=Dict(j=>f.maps for j in 1:3),isos=Dict(j=>Float64(j) for j in 1:3))
    mktempdir() do dir
        upstream=joinpath(dir,"upstream"); controls=joinpath(dir,"controls"); foreground=joinpath(dir,"foreground"); mkdir(foreground)
        R.one_direction(toy_image(),"fwd",maps,C,upstream); L.one_direction(upstream,maps,C,CONTROL,controls)
        anchors=P.t(joinpath(upstream,"anchors.tsv"))
        member=[(;anchor=P.i(a,"anchor"),used=parse(Bool,a["used"]),initial_available=true,initial_core=true,joint_available=true,joint_core=false) for a in anchors]
        R.table(joinpath(foreground,"anchor_membership.tsv"),member)
        out=joinpath(dir,"replay"); result=P.one_direction(upstream,controls,foreground,maps,C,REPLAY,out)
        rows=P.t(joinpath(out,"patches.tsv"))
        @test result.status=="ok" && length(rows)==15result.patches
        @test all(r["initial_core"]=="true" && r["joint_core"]=="false" for r in rows)
        verified=W.verify_direction(out,upstream,controls,foreground,maps,C,REPLAY)
        @test verified.patches==result.patches
        reported=W.case_rows(verified.rows,"synthetic.sxm","fwd",REPLAY)
        compared=W.comparisons(reported)
        @test length(reported)==120 && length(compared)==192
        @test all(r.patches==result.patches for r in reported if r.mask=="initial" && r.group=="retained")
        @test all(r.patches==0 for r in reported if r.mask=="joint_background" && r.group=="retained")
        R.table(joinpath(foreground,"anchor_membership.tsv"),[merge(r,(initial_core=false,joint_core=true)) for r in member])
        changed=joinpath(dir,"membership_changed"); P.one_direction(upstream,controls,foreground,maps,C,REPLAY,changed)
        other=P.t(joinpath(changed,"patches.tsv"))
        @test all(a[k]==b[k] for (a,b) in zip(rows,other) for k in keys(a) if !(k in ("initial_core","joint_core")))
        @test_throws ErrorException P.one_direction(upstream,controls,foreground,maps,C,REPLAY,out)
    end
    @test_throws ErrorException P.main(["--expected-N","6"])
    mktempdir() do dir
        for mutation in (c->(c["model"]["composition"]=.5),c->(c["model"]["height_gain"]=2),
            c->(c["selection"]["masks"]=["joint_background"]),c->(c["preprocessing"]["target_use"]="choose"))
            bad=deepcopy(REPLAY); mutation(bad); path=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException P.settings(path)
        end
    end
end
