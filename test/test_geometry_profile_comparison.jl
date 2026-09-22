using Test, TOML, LinearAlgebra, Random
include(joinpath(@__DIR__, "run_geometry_profile_comparison.jl"))
const R = FixedNGeometryRefinement
const G = R.G
const VP = R.VP
const SETTINGS_PATH = joinpath(ROOT, "config", "geometry_profile_comparison.toml")
const PHYS = TOML.parsefile(joinpath(ROOT, "config", "chitosan.toml"))
const RAWSET = TOML.parsefile(SETTINGS_PATH)
BLAS.set_num_threads(1)

function fixture(n=3; circular=false)
    _, cfg, _ = R.F.Extractor._configs(PHYS["model"],PHYS["preprocessing"],"unused")
    cfg.chain_circular_sigmas = circular
    axis = (origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
    xs = collect(range(-1.8,1.8;length=31)); ys = collect(range(-.5,.5;length=11))
    x = repeat(xs;inner=length(ys)); y = repeat(ys;outer=length(xs))
    p = zeros(G._chain_nparams(n,cfg)); p[1:3] .= [.01,.001,-.002]
    z = G._chain_model_values(x,y,p,n,axis,cfg;amp_min=.03,amp_range=.07)
    zimg = reshape(z,length(ys),length(xs))
    data = (xs=xs,ys=ys,zimg=zimg,x=x,y=y,z=z,zfull=z,noise=.01,axisctx=axis)
    p0 = copy(p); p0[4:3+n] .+= .1; p0[end] += .15
    initial = R.evaluate(p0,n,data,cfg,.03,.07).result
    @assert initial.valid
    return initial,data,cfg,p
end

@testset "Matched options, objective and native diagnostics" begin
    options, algorithm = R.settings(SETTINGS_PATH)
    @test algorithm == :LN_BOBYQA
    @test options.outer_maxeval == 800 && options.outer_maxtime_s == 30
    old = TOML.parsefile(joinpath(ROOT,"config","label_free_exploration.toml"))
    @test RAWSET["counting_variable_projection"] == old["counting_variable_projection"]
    for n in (1,3), circular in (true,false)
        initial,data,cfg,p = fixture(n;circular)
        e = R.evaluate(p,n,data,cfg,.03,.07)
        @test e.valid && e.result.rss <= 1e-28
        _,_,ts,_,sp,sq = G._decode_chain(initial.params,n,data.axisctx,cfg;amp_min=.03,amp_range=.07)
        kappa = n == 1 ? 1. : G.adjacent_kappa_max(diff(ts),max.(sp,sq))
        expected = initial.rss * (1 + (n == 1 ? 0. : G.kappa_penalty(kappa;kappa_max=cfg.kappa_max,weight=cfg.kappa_weight)))
        @test R.evaluate(initial.params,n,data,cfg,.03,.07).objective ≈ expected
        @test initial.gcv ≈ length(data.z)/(length(data.z)-G._chain_nparams(n,cfg))^2 * initial.rss
        rows = R.feature_rows("sample.sxm",initial,data,cfg,circular ? "circ" : "ell")
        @test R.F.validate_chain(rows) == n
        @test all(parse(Int,r["N"]) == n for r in rows)
    end
    # Exercise an active penalty without changing the physical production config.
    initial,data,cfg,_ = fixture(3)
    cfg.kappa_max = 1.01
    penalized = R.evaluate(initial.params,3,data,cfg,.03,.07)
    factor = 1 + G.kappa_penalty(penalized.result.kappa_max_adj;
        kappa_max=cfg.kappa_max,weight=cfg.kappa_weight)
    @test factor > 1
    @test penalized.objective ≈ penalized.result.rss * factor
    prof = VP.fixed_geometry_profile(initial.params,3,data.x,data.y,data.z,data.axisctx,cfg;
        amp_min=.03,amp_range=.07,options)
    after = R.evaluate(prof.params,3,data,cfg,.03,.07)
    @test prof.success && after.valid
    @test after.result.kappa_max_adj == penalized.result.kappa_max_adj
    @test after.objective ≈ after.result.rss * factor
    @test after.objective <= penalized.objective
    mktempdir() do dir
        path=joinpath(dir,"settings.toml")
        raw=deepcopy(RAWSET); raw["model"]["max_overlap"]=1
        open(io->TOML.print(io,raw),path,"w")
        @test_throws ErrorException R.settings(path)
    end
end

@testset "Same-start paired optimization, saved trace and budget stops" begin
    options = VP.read_options(merge(RAWSET["counting_variable_projection"],Dict("outer_maxeval"=>120)))
    for circular in (true,false)
        initial,data,cfg,_ = fixture(3;circular)
        p0 = copy(initial.params); alltraces = Dict()
        results = Dict()
        for mode in ("joint","profiled")
            trace=[]
            r = R.refine(initial,data,cfg,mode,options,:LN_BOBYQA;trace=x->push!(trace,x))
            results[mode] = r; alltraces[mode] = trace
            @test initial.params == p0
            @test r.best.valid && r.best.objective <= r.initial_objective
            @test r.evaluations == length(trace) <= options.outer_maxeval
            @test [x.evaluation for x in trace] == collect(1:length(trace))
            @test r.rejected == count(x->!x.valid,trace)
            @test r.best.result.n == initial.n
            @test r.best.result.gcv ≈ VP.full_gcv(r.best.result.rss,length(data.z),initial.n,cfg)
            @test r.status in ("FTOL_REACHED","XTOL_REACHED","MAXEVAL_REACHED","MAXTIME_REACHED","SUCCESS","ROUNDOFF_LIMITED")
            @test r.converged == (r.status in ("FTOL_REACHED","XTOL_REACHED","SUCCESS"))
            lo,hi=VP.native_raw_bounds(initial.n,cfg)
            @test all(lo .<= r.best.result.params .<= hi)
            if mode == "profiled"
                @test r.profile !== nothing && r.profile.success
                @test r.profile.linear.kkt_violation <= r.profile.linear.kkt_tolerance
                @test all(x.linear_kkt <= x.linear_tolerance for x in trace)
                @test r.best.result.params == r.profile.params
            end
        end
        @test results["joint"].initial_objective == results["profiled"].initial_objective
        @test results["joint"].initial_rss == results["profiled"].initial_rss
        @test_throws ErrorException R.refine(initial,data,cfg,"unknown",options,:LN_BOBYQA)
        # Deterministic synthetic replay, far from the time ceiling.
        replay = R.refine(initial,data,cfg,"profiled",options,:LN_BOBYQA)
        @test replay.best.result.params == results["profiled"].best.result.params
        @test replay.evaluations == results["profiled"].evaluations
        invalidcfg=deepcopy(cfg); invalidcfg.max_overlap=0.
        @test !R.evaluate(p0,initial.n,data,invalidcfg,.03,.07).valid
        @test_throws ErrorException R.refine(initial,data,invalidcfg,"joint",options,:LN_BOBYQA)
    end
    @test_throws ErrorException R.parse_cli(["--expected-N","6"])
    @test_throws ErrorException R.parse_cli(["--benchmark-manifest","labels.toml"])
end

@testset "Three complete arms, shared split cache, no selectors or labels" begin
    mktempdir() do tmp
        raw=joinpath(tmp,"raw"); mkdir(raw); touch(joinpath(raw,"sample.sxm"))
        initial,data,cfg,_=fixture(3)
        rows=R.feature_rows("sample.sxm",initial,data,cfg,"ell"); header=sort(collect(keys(first(rows))))
        base=joinpath(tmp,"base.tsv"); split=joinpath(tmp,"split.tsv")
        write_table(base,header,rows); write_table(split,header,rows)
        templates=joinpath(tmp,"templates.tsv"); touch(templates)
        args=["--data-dir",raw,"--count-config",joinpath(ROOT,"config","chitosan.toml"),
            "--config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",SETTINGS_PATH,"--features",base,"--split-features",split,
            "--templates",templates,"--outdir",joinpath(tmp,"out")]
        geometry_comparison(vcat(args,["--dry-run"]))
        @test !ispath(joinpath(tmp,"out"))
        calls=[]
        runner=function(out,stage,script,as)
            push!(calls,(stage,script,as)); target=as[findfirst(==("--outdir"),as)+1]; mkpath(target)
            write_table(joinpath(target,"predictions.tsv"),header,rows)
        end
        exporter=function(out,stage,script,as,path;nfiles)
            @test nfiles==1 && script=="refine_fixed_n_geometry.jl"
            write_table(path,header,rows)
            for mode in ("native","joint","profiled"); write_table(path*".$mode.tsv",header,rows); end
        end
        geometry_comparison(args;runner,exporter)
        @test first.(calls)==["reference","joint","profiled"]
        @test all(c[2]=="run_reconstructed_chitosan.jl" for c in calls)
        for arm in ("reference","joint","profiled")
            @test read(joinpath(tmp,"out",arm,"features_split.tsv"))==read(split)
        end
        @test read(joinpath(tmp,"out","reference","features.tsv"))==read(base)
        @test_throws ErrorException geometry_comparison(args)
        failargs=replace.(args,joinpath(tmp,"out")=>joinpath(tmp,"failed"))
        @test_throws ErrorException geometry_comparison(failargs;runner=(args...)->error("intentional"))
        @test isfile(joinpath(tmp,"failed","failures.tsv"))
        touch(joinpath(raw,"extra.sxm"))
        @test_throws ErrorException geometry_comparison(vcat(replace.(args,joinpath(tmp,"out")=>joinpath(tmp,"other")),["--dry-run"]))
    end
    script=read(joinpath(@__DIR__,"refine_fixed_n_geometry.jl"),String)
    @test !occursin("chain_gaussian_sweep(",script)
    @test !occursin("_select_chain_model(",script)
    shell=read(joinpath(ROOT,"hpc","compare_geometry_profile.sbatch"),String)
    @test occursin("#SBATCH --time=04:00:00",shell)
    @test occursin("SLURM_JOB_ID:?",shell)
    @test success(`bash -n $(joinpath(ROOT,"hpc","compare_geometry_profile.sbatch"))`)
end
