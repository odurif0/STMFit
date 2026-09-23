#!/usr/bin/env julia
# Synthetic-only tests; no real images or grading inputs.
using Test, TOML, LinearAlgebra, Random
include(joinpath(@__DIR__,"profile_cross_view_residuals.jl"))
const C=CrossViewResidual
const F=C.F
const G=C.G
const ROOT=normpath(joinpath(@__DIR__,".."))
const SETTINGS=joinpath(ROOT,"config","cross_view_residual.toml")
const OPTIONS=C.read_settings(SETTINGS)
const RAW=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
BLAS.set_num_threads(1)

function cross_fixture()
    _,cfg,_=F.Extractor._configs(RAW["model"],RAW["preprocessing"],"unused")
    axis=(origin=(0.0,0.0),axis=[1.0,0.0],perp=[0.0,1.0],tmin=-1.0,tmax=1.0)
    rows=[Dict("file"=>"sample.sxm","N"=>"2","lobe"=>string(i),"source"=>"ell",
        "amplitude"=>"0.6","x_nm"=>string(cx),"y_nm"=>"0.0","t_nm"=>string(cx),"u_nm"=>"0.0",
        "sigma_parallel_nm"=>"0.18","sigma_perp_nm"=>"0.15","skew_ratio"=>"1.0",
        "axis_x"=>"1.00000000","axis_y"=>"0.00000000","origin_x_nm"=>"0.000000","origin_y_nm"=>"0.000000",
        "baseline"=>"0.01","tilt_x"=>"0.0","tilt_y"=>"0.0","amp_rel"=>"1.0","gcv"=>"0.01") for (i,cx) in enumerate((-0.35,0.35))]
    rng=MersenneTwister(923)
    x=2rand(rng,201).-1; y=rand(rng,201).-0.5
    A=F.frozen_design(x,y,rows,true)
    lo,hi=F.coefficient_bounds(2,cfg,axis,[1.0;;])
    zf=A*[0.01,0.02,-0.01,0.5,0.7]
    zb=A*[0.02,-0.01,0.02,0.7,0.5]
    return (;rows,cfg,axis,x,y,A,lo,hi,zf,zb)
end

module Driver
include(joinpath(@__DIR__,"run_cross_view_residual_comparison.jl"))
end
@testset "Three matched arms; no grading or incomplete cohort" begin
    mktempdir() do tmp
        t=cross_fixture(); header=sort(collect(keys(first(t.rows))))
        raw=joinpath(tmp,"raw"); mkdir(raw); touch(joinpath(raw,"sample.sxm"))
        features=F.write_table(joinpath(tmp,"features.tsv"),header,t.rows)
        split=F.write_table(joinpath(tmp,"split.tsv"),header,t.rows)
        templates=joinpath(tmp,"templates.tsv"); touch(templates)
        out=joinpath(tmp,"out")
        args=["--data-dir",raw,"--features",features,"--split-features",split,"--templates",templates,
            "--config",joinpath(ROOT,"config/unit_assignment_patch_support.toml"),"--count-config",joinpath(ROOT,"config/chitosan.toml"),
            "--settings",SETTINGS,"--outdir",out]
        Driver.cross_view_comparison(vcat(args,["--dry-run"]))
        @test !ispath(out)
        calls=[]
        fake_export=function(out,stage,script,as,output;nfiles)
            @test stage=="profile" && script=="profile_cross_view_residuals.jl" && nfiles==1
            for (path,view) in ((output,"fwd"),(output*".bwd.tsv","bwd"))
                F.write_table(path,vcat(header,["profile_view"]),[merge(r,Dict("profile_view"=>view)) for r in t.rows])
            end
        end
        fake_runner=function(out,stage,script,as)
            push!(calls,(stage,script,as)); arm=as[findfirst(==("--outdir"),as)+1]; mkpath(arm)
            F.write_table(joinpath(arm,"predictions.tsv"),header,t.rows)
        end
        Driver.cross_view_comparison(args;runner=fake_runner,exporter=fake_export)
        @test first.(calls)==["reference","same_view","cross_view"]
        @test all(c[2]=="run_reconstructed_chitosan.jl" for c in calls)
        @test !("--residual-features-fwd" in calls[1][3])
        @test calls[2][3][findfirst(==("--residual-features-fwd"),calls[2][3])+1]==joinpath(out,"profile_fwd.tsv")
        @test calls[3][3][findfirst(==("--residual-features-fwd"),calls[3][3])+1]==joinpath(out,"profile_bwd.tsv")
        for arm in ("reference","same_view","cross_view")
            @test read(joinpath(out,arm,"features.tsv"))==read(features)
            @test read(joinpath(out,arm,"features_split.tsv"))==read(split)
        end
        @test_throws ErrorException Driver.cross_view_comparison(args)
        failed=joinpath(tmp,"failed")
        @test_throws ErrorException Driver.cross_view_comparison(replace.(args,out=>failed);exporter=(a...;kw...)->error("intentional"))
        @test isfile(joinpath(failed,"failures.tsv")) && !ispath(joinpath(failed,"reference"))
        baseopts=Dict("--config"=>joinpath(ROOT,"config/unit_assignment_patch_support.toml"),"--outdir"=>joinpath(tmp,"no_output"),
            "--features"=>features,"--split-features"=>split,"--residual-features-fwd"=>joinpath(out,"profile_fwd.tsv"),
            "--residual-features-bwd"=>joinpath(out,"profile_bwd.tsv"),"--dry-run"=>"true")
        Driver.execute_pipeline(baseopts)
        @test !ispath(baseopts["--outdir"])
        for k in ("--features","--split-features","--residual-features-bwd")
            bad=copy(baseopts); delete!(bad,k); @test_throws ErrorException Driver.execute_pipeline(bad)
        end
        for k in ("--patches-fwd","--patches-bwd","--descriptor-patches","--patch-frames","--acquisition-shifts")
            @test_throws ErrorException Driver.execute_pipeline(merge(baseopts,Dict(k=>"bad.tsv")))
        end
        @test !ispath(baseopts["--outdir"])
        touch(joinpath(raw,"extra.sxm"))
        @test_throws ErrorException Driver.cross_view_comparison(vcat(replace.(args,out=>joinpath(tmp,"other")),["--dry-run"]))
    end
    shell=joinpath(ROOT,"hpc/compare_cross_view_residuals.sbatch")
    @test success(`bash -n $shell`)
    @test occursin("#SBATCH --time=02:00:00",read(shell,String))
    @test occursin("SLURM_JOB_ID:?",read(shell,String))
end

@testset "Two conditional source-only linear profiles" begin
    t=cross_fixture()
    pair=C.profile_pair(t.rows,t.x,t.y,t.zf,t.zb,t.lo,t.hi,t.cfg,OPTIONS)
    @test pair.fwd.result.x≈[0.01,0.02,-0.01,0.5,0.7] atol=1e-10
    @test pair.bwd.result.x≈[0.02,-0.01,0.02,0.7,0.5] atol=1e-10
    changed=C.profile_pair(t.rows,t.x,t.y,t.zf.+0.001,t.zb,t.lo,t.hi,t.cfg,OPTIONS)
    @test changed.bwd.result.x==pair.bwd.result.x
    @test changed.fwd.result.x!=pair.fwd.result.x
    swapped=C.profile_pair(t.rows,t.x,t.y,t.zb,t.zf,t.lo,t.hi,t.cfg,OPTIONS)
    @test swapped.fwd.result.x==pair.bwd.result.x && swapped.bwd.result.x==pair.fwd.result.x
    d=C.mean_diagnostics(t.rows,t.x,t.y,(t.zf+t.zb)/2,0.01,t.axis,t.cfg,pair)
    @test d.valid && d.snr<1e-10
    @test d.pfull==G._chain_nparams(2,t.cfg)>5
    @test d.gcv≈201/(201-d.pfull)^2*d.rss
    @test d.overlap≈exp(-0.5*(0.7/0.18)^2)
    @test d.overrun==0
    for r in (pair.fwd,pair.bwd)
        @test r.result.converged && r.result.feasible
        @test r.result.kkt_violation<=r.result.kkt_tolerance
        @test r.result.rss<=r.initial_rss
        for (a,b) in zip(t.rows,r.rows),k in setdiff(keys(a),F.MUTABLE_COLUMNS)
            @test a[k]==b[k]
        end
    end
    badz=(t.zf+t.zb)/2; badz[1]+=1
    @test !C.mean_diagnostics(t.rows,t.x,t.y,badz,0.01,t.axis,t.cfg,pair).valid
    narrow=merge(t.axis,(tmin=-0.2,tmax=0.2))
    @test !C.mean_diagnostics(t.rows,t.x,t.y,(t.zf+t.zb)/2,0.01,narrow,t.cfg,pair).valid
    tight=deepcopy(t.cfg); tight.max_overlap=0.0
    @test !C.mean_diagnostics(t.rows,t.x,t.y,(t.zf+t.zb)/2,0.01,t.axis,tight,pair).valid
    @test_throws ErrorException C.parse_cli(["--labels","no"])
    @test_throws ErrorException C.parse_cli(["--expected-N","6"])
    source=read(joinpath(@__DIR__,"profile_cross_view_residuals.jl"),String)
    @test !occursin("_fit_chain_n(",source) && !occursin("refine_profile(",source)
end

@testset "Strict settings and metadata before any output" begin
    mktempdir() do tmp
        t=cross_fixture(); header=sort(collect(keys(first(t.rows))))
        features=joinpath(tmp,"features.tsv"); F.write_table(features,header,t.rows)
        raw=joinpath(tmp,"raw"); mkdir(raw); touch(joinpath(raw,"sample.sxm"))
        out=joinpath(tmp,"profile.tsv")
        args=["--features",features,"--data-dir",raw,"--config",joinpath(ROOT,"config","chitosan.toml"),"--settings",SETTINGS,"--out",out]
        C.main(vcat(args,["--dry-run"]))
        @test !ispath(out) && !ispath(out*".fit_data")
        touch(out*".bwd.tsv")
        @test_throws ErrorException C.parse_cli(args)
        settings=TOML.parsefile(SETTINGS); settings["cross_view_residual"]["validity"]="none"
        bad=joinpath(tmp,"bad.toml"); open(io->TOML.print(io,settings),bad,"w")
        @test_throws ErrorException C.read_settings(bad)
        settings["cross_view_residual"]["validity"]="native_fused_mean"; settings["model"]["N"]=6
        open(io->TOML.print(io,settings),bad,"w")
        @test_throws ErrorException C.read_settings(bad)
    end
end
