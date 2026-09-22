#!/usr/bin/env julia
# Synthetic and metadata checks only. No cohort fitting or external labels.
using Test, TOML, LinearAlgebra, Printf, Random
include(joinpath(@__DIR__,"run_frozen_amplitude_comparison.jl"))
const F=FrozenAmplitudeProfile
const VP=F.VP
const G=F.G
BLAS.set_num_threads(1)
const SETTINGS=TOML.parsefile(joinpath(ROOT,"config","label_free_exploration.toml"))
const OPT=VP.read_options(SETTINGS["counting_variable_projection"])
const RAW=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))

function fixture(n=3;circular=false)
    _,cfg,_=F.Extractor._configs(RAW["model"],RAW["preprocessing"],"unused")
    cfg.chain_circular_sigmas=circular
    axis=(origin=(0.2,-0.1),axis=[0.8,0.6],perp=[-0.6,0.8],tmin=-1.4,tmax=1.5)
    p=zeros(G._chain_nparams(n,cfg)); p[1:3].=[0.01,0.005,-0.002]
    _,feats,ts,us,sp,sq=G._decode_chain(p,n,axis,cfg;amp_min=0.3,amp_range=0.7)
    rows=[Dict("file"=>"sample.sxm","N"=>string(n),"lobe"=>string(i),"source"=>(circular ? "circ" : "ell"),
        "amplitude"=>@sprintf("%.8e",feats[i].amplitude),"x_nm"=>@sprintf("%.6f",feats[i].x_nm),
        "y_nm"=>@sprintf("%.6f",feats[i].y_nm),"t_nm"=>@sprintf("%.6f",ts[i]),"u_nm"=>@sprintf("%.6f",us[i]),
        "sigma_parallel_nm"=>@sprintf("%.6f",sp[i]),"sigma_perp_nm"=>@sprintf("%.6f",sq[i]),
        "skew_ratio"=>"1.000000","axis_x"=>"0.80000000","axis_y"=>"0.60000000",
        "origin_x_nm"=>"0.200000","origin_y_nm"=>"-0.100000","baseline"=>@sprintf("%.8e",p[1]),
        "tilt_x"=>@sprintf("%.8e",p[2]),"tilt_y"=>@sprintf("%.8e",p[3]),"amp_rel"=>"1.000000",
        "spacing_prev_nm"=>"NA","gcv"=>"0.01") for i in 1:n]
    return rows,cfg,axis,p
end

@testset "Frozen decoded basis, native endpoints, full complexity" begin
    rng=MersenneTwister(9622)
    for n in (1,3,8),circular in (false,true)
        rows,cfg,axis,p=fixture(n;circular)
        x,y=randn(rng,151),randn(rng,151)
        A=F.frozen_design(x,y,rows,true)
        for (i,r) in enumerate(rows)
            dx=x.-parse(Float64,r["x_nm"]); dy=y.-parse(Float64,r["y_nm"])
            dt=0.8dx+0.6dy; du=-0.6dx+0.8dy
            expected=exp.(-0.5.*((dt./parse(Float64,r["sigma_parallel_nm"])).^2+(du./parse(Float64,r["sigma_perp_nm"])).^2))
            @test A[:,3+i]≈expected atol=1e-15
        end
        @test A[:,1]==ones(151) && A[:,2]==x && A[:,3]==y
        lo,hi=F.coefficient_bounds(n,cfg,axis,[0.1 1.0])
        @test lo[1:3]==[-5,-1,-1] && hi[1:3]==[5,1,1]
        @test all(isapprox.(lo[4:end],0.3+0.7/(1+exp(5))))
        @test all(isapprox.(hi[4:end],0.3+0.7/(1+exp(-5))))
        target=vcat([0.013,0.002,-0.001],n==1 ? [0.6] : collect(range(0.4,0.8;length=n)))
        z=A*target
        r=F.profile_rows(rows,x,y,z,lo,hi,cfg,OPT)
        @test r.result.converged && r.result.feasible
        @test r.result.x≈target atol=1e-8
        @test r.result.rss<=r.initial_rss
        @test r.pfull==G._chain_nparams(n,cfg)>n+3
        @test r.gcv≈151/(151-r.pfull)^2*r.result.rss
        for (a,b) in zip(rows,r.rows),c in setdiff(keys(a),F.MUTABLE_COLUMNS)
            @test a[c]==b[c]
        end
        recovered=vcat([parse(Float64,r.rows[1][c]) for c in ("baseline","tilt_x","tilt_y")],
            [parse(Float64,row["amplitude"]) for row in r.rows])
        @test recovered==r.result.x
    end
end

@testset "Bounded solution, serialization and explicit failure" begin
    rows,cfg,axis,_=fixture(2)
    x=collect(range(-1.5,1.5;length=101)); y=sin.(x)
    A=F.frozen_design(x,y,rows,true); lo,hi=F.coefficient_bounds(2,cfg,axis,[1.0;;])
    z=A*[0.02,0.0,0.0,0.05,1.8]
    r=F.profile_rows(rows,x,y,z,lo,hi,cfg,OPT)
    @test r.result.converged && r.result.feasible
    @test !isempty(r.result.active_lower) || !isempty(r.result.active_upper)
    @test r.result.kkt_violation<=r.result.kkt_tolerance
    @test F.decimal_radius("3.12345678e-02")≈5e-11
    @test F.decimal_radius("0.123456")≈5e-7
    @test_throws ErrorException F.validate_chain([merge(rows[1],Dict("N"=>"3")),rows[2]])
    for (key,val) in (("skew_ratio","1.1"),("amplitude","NaN"),("sigma_parallel_nm","-0.2"),("source","unknown"))
        broken=deepcopy(rows); broken[1][key]=val
        @test_throws ErrorException F.validate_chain(broken)
    end
    broken=deepcopy(rows); broken[1]["amplitude"]="3.0"
    @test_throws ErrorException F.profile_rows(broken,x,y,z,lo,hi,cfg,OPT)
    tiny=VP.read_options(merge(SETTINGS["counting_variable_projection"],Dict("linear_maxiter"=>1)))
    @test_throws ErrorException F.profile_rows(rows,x,y,z,lo,hi,cfg,tiny)
    @test_throws ErrorException F.parse_cli(["--expected-N","6"])
    @test_throws ErrorException F.parse_cli(["--benchmark-manifest","labels.toml"])
end

@testset "Metadata, two frozen arms and no nonlinear fitter" begin
    mktempdir() do tmp
        raw=joinpath(tmp,"raw"); mkdir(raw); touch(joinpath(raw,"sample.sxm"))
        rows,_,_,_=fixture(2); header=sort(collect(keys(rows[1])))
        features=joinpath(tmp,"features.tsv"); split=joinpath(tmp,"split.tsv")
        write_table(features,header,rows); write_table(split,header,rows)
        templates=joinpath(tmp,"templates.tsv"); touch(templates)
        args=["--data-dir",raw,"--count-config",joinpath(ROOT,"config","chitosan.toml"),
            "--config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",joinpath(ROOT,"config","label_free_exploration.toml"),"--features",features,
            "--split-features",split,"--templates",templates,"--outdir",joinpath(tmp,"out")]
        profile_comparison(vcat(args,["--dry-run"]))
        @test !ispath(joinpath(tmp,"out"))
        calls=[]
        fake_runner=function(out,stage,script,as)
            push!(calls,(stage,script,as)); arm=as[findfirst(==("--outdir"),as)+1]; mkpath(arm)
            write_table(joinpath(arm,"predictions.tsv"),header,rows)
        end
        fake_export=function(out,stage,script,as,output; nfiles)
            @test stage=="profile" && script=="profile_frozen_amplitudes.jl" && nfiles==1
            write_table(output,header,rows)
        end
        profile_comparison(args;runner=fake_runner,exporter=fake_export)
        @test first.(calls)==["control","profiled"]
        @test all(c[2]=="run_reconstructed_chitosan.jl" for c in calls)
        @test all("--features" in c[3] && "--split-features" in c[3] for c in calls)
        @test read(joinpath(tmp,"out","control","features.tsv"))==read(features)
        @test read(joinpath(tmp,"out","profiled","features_split.tsv"))==read(split)
        @test_throws ErrorException profile_comparison(args)
        out2=joinpath(tmp,"failed"); failed=replace.(args,joinpath(tmp,"out")=>out2)
        @test_throws ErrorException profile_comparison(failed;runner=(args...)->error("intentional"))
        @test isfile(joinpath(out2,"failures.tsv"))
        touch(joinpath(raw,"extra.sxm"))
        other=replace.(args,joinpath(tmp,"out")=>joinpath(tmp,"other"))
        @test_throws ErrorException profile_comparison(vcat(other,["--dry-run"]))
    end
    source=read(joinpath(@__DIR__,"profile_frozen_amplitudes.jl"),String)
    @test !occursin("refine_profile(",source) && !occursin("_fit_chain_n(",source) && !occursin("chain_gaussian_sweep(",source)
    shell=read(joinpath(ROOT,"hpc","compare_frozen_amplitudes.sbatch"),String)
    @test occursin("#SBATCH --time=02:00:00",shell)
    @test occursin("SLURM_JOB_ID:?",shell) && occursin("SLURM_CPUS_PER_TASK",shell)
    @test occursin("config/unit_assignment_patch_support.toml",shell)
    @test success(`bash -n $(joinpath(ROOT,"hpc","compare_frozen_amplitudes.sbatch"))`)
end
