using Test, TOML
include(joinpath(@__DIR__,"verify_local_sigma_counting.jl"))
const V=LocalSigmaVerification
const C=V.C
const G=C.G
const P=C.P
const ROOT=dirname(@__DIR__)
const COUNT=joinpath(ROOT,"config/chitosan.toml")
const SETTINGS=joinpath(ROOT,"config/local_sigma_counting.toml")
const PHYS=TOML.parsefile(COUNT)

function image_fixture(;partial=false)
    f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
    b=copy(f)
    partial && (f[end-2:end,:].=NaN;b[3,4]=Inf)
    G.SXMImage("synthetic.sxm",Dict{String,String}(),65,41,(5.,2.),(0.,0.),
        [G.SXMChannel("Z","nm","fwd",f),G.SXMChannel("Z","nm","bwd",b)])
end

function fixture(dir)
    rawdir=joinpath(dir,"raw");mkdir(rawdir);touch(joinpath(rawdir,"a.sxm"))
    raw=deepcopy(PHYS)
    raw["model"]["global_maxtime"]=.05;raw["model"]["global_maxiter"]=10;raw["model"]["max_iter"]=10
    raw["preprocessing"]["flatten"]="none";raw["preprocessing"]["smooth_radius_px"]=0
    config=joinpath(dir,"count.toml");open(io->TOML.print(io,raw),config,"w")
    Dict("--data-dir"=>rawdir,"--count-config"=>config,"--settings"=>SETTINGS,"--outdir"=>joinpath(dir,"out"))
end

@testset "Counting scope, exhaustive range and valid GCV only" begin
    s=C.settings(SETTINGS)
    for arm in C.ARMS
        cfg=C.configs(PHYS,s,arm)
        @test cfg.circ.overlap_constraint==cfg.ell.overlap_constraint==arm
        @test !cfg.circ.intelligent_sweep && !cfg.ell.intelligent_sweep
        @test cfg.circ.n_min==2 && cfg.circ.n_max==14 && cfg.ell.skip_global && cfg.ell.max_iter==50
        @test isempty(cfg.circ.init_centers_t) && isempty(cfg.ell.init_centers_t)
        @test cfg.circ.spacing_min_nm==PHYS["model"]["spacing_min_nm"]
        @test cfg.circ.max_overlap==PHYS["model"]["max_overlap"]
        for partial in (false,true)
            img=image_fixture(;partial);d=C.data_bundle(img,cfg.pcfg,cfg.circ)
            other=C.data_bundle(img,cfg.pcfg,cfg.ell)
            @test isequal(d,other) && all(isfinite,d.z)
            @test d.observed==65*41-(partial ? 196 : 0)
        end
    end
    for partial in (false,true)
        a=C.configs(PHYS,s,C.ARMS[1]);b=C.configs(PHYS,s,C.ARMS[2]);img=image_fixture(;partial)
        @test isequal(C.data_bundle(img,a.pcfg,a.circ),C.data_bundle(img,b.pcfg,b.circ))
    end
    d=(axisctx=(tmin=0.,tmax=1.3),)
    @test collect(C.candidate_ns(d,C.configs(PHYS,s,C.ARMS[1]).circ))==[2,3]
    @test collect(C.candidate_ns(d,C.configs(PHYS,s,C.ARMS[2]).circ))==[2,3,4]
    @test isempty(C.candidate_ns((axisctx=(tmin=0.,tmax=.1),),C.configs(PHYS,s,C.ARMS[1]).circ))
    r(n,g;valid=true,success=true)=G.ChainModelResult(n=n,gcv=g,valid=valid,success=success)
    q=[(result=r(7,0.;valid=false),family="ell"),(result=r(4,1.),family="circ"),
       (result=r(3,1.),family="circ"),(result=r(3,1.),family="ell"),(result=r(2,Inf),family="ell")]
    @test C.select_candidate(q).result.n==3 && C.select_candidate(q).family=="ell"
    @test C.select_candidate(q[1:1])===nothing && C.select_candidate(NamedTuple[])===nothing
    for k in ("--manifest","--expected-N","--sequence","--truth","--selected-summary","--features","--selection-policy")
        @test_throws ErrorException C.options([k,"bad"])
    end
    @test C.options(["--help"])===nothing
    mktempdir() do dir
        o=fixture(dir);args=vcat([[k,v] for (k,v) in o]...)
        @test C.options(args)==o
        @test_throws ErrorException C.options(vcat(args,["--chunk","0/4"]))
        @test_throws ErrorException C.options(vcat(args,["--settings",SETTINGS]))
        C.execute(merge(o,Dict("--dry-run"=>"true"));reader=p->error("Pixels read during dry-run"))
        @test !ispath(o["--outdir"])
        bad=deepcopy(s);bad["model"]["n_max"]=6
        path=joinpath(dir,"bad.toml");open(io->TOML.print(io,bad),path,"w")
        @test_throws ErrorException C.settings(path)
    end
end

@testset "Small synthetic raw-to-saved-state comparison and failure accounting" begin
    mktempdir() do dir
        o=fixture(dir)
        C.execute(o;reader=p->image_fixture(;partial=true))
        checked=V.verify(o["--outdir"];rawdir=o["--data-dir"])
        @test length(checked.selected)==4 && checked.files==Set(["a.sxm"])
        @test all(r["status"]=="ok" for r in checked.selected)
        @test all(only(filter(q->q["file"]==r["file"] && q["repetition"]==r["repetition"] && q["arm"]==r["arm"] &&
            q["N"]==r["N_selected"] && q["family"]==r["source"],checked.fits))["valid"]=="true" for r in checked.selected)
        @test_throws ErrorException C.options(vcat([[k,v] for (k,v) in o]...))
        failed=merge(o,Dict("--outdir"=>joinpath(dir,"failed")))
        @test_throws ErrorException C.execute(failed;reader=p->image_fixture(),fitter=(args...;kwargs...)->error("synthetic failure"))
        rs=V.rows(joinpath(failed["--outdir"],"selected.tsv"))
        @test length(rs)==4 && all(r["status"]=="unavailable" for r in rs)
        fs=V.rows(joinpath(failed["--outdir"],"candidates.tsv"))
        @test all(r["success"]==r["valid"]=="false" for r in fs)
    end
end

@testset "Simple four-shard Slurm boundary" begin
    script=joinpath(ROOT,"hpc/diagnose_local_sigma_counting.sbatch")
    @test success(`bash -n $script`)
    text=read(script,String)
    @test occursin("--time=04:00:00",text) && occursin("--cpus-per-task=4",text)
    @test occursin("SLURM_JOB_ID",text) && occursin("--dry-run",text)
    @test occursin("verify_local_sigma_counting.jl",text)
    @test !occursin("grade_chitosan",text) && !occursin("sbatch ",replace(text,"sbatch --export"=>""))
end
