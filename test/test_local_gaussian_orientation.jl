using Test, TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__,"run_local_gaussian_comparison.jl"))
const L = LocalGaussianOrientation
const G = L.G
const VP = L.VP
const SETTINGS = joinpath(ROOT,"config/local_gaussian_orientation.toml")
const COUNT = joinpath(ROOT,"config/chitosan.toml")
const ASSIGNMENT = joinpath(ROOT,"config/unit_assignment_patch_support.toml")
BLAS.set_num_threads(1)

function fixture(n=3)
    raw=TOML.parsefile(COUNT)
    _,ell,circ=L.F.Extractor._configs(raw["model"],raw["preprocessing"],"unused")
    axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
    xs=collect(range(-1.8,1.8;length=33)); ys=collect(range(-.7,.7;length=17))
    x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
    p=zeros(G._chain_nparams(n,ell)); p[1:3].=[.01,.001,-.002]
    k=3+n+G._chain_spacing_param_count(n,ell)
    p[k+1:k+n].=collect(range(.2,.6;length=n)); p[k+div(n,2)+1]=-.5
    p[k+n+1:k+2n].=-.8; p[k+2n+1:k+3n].=.7
    localcfg=deepcopy(ell); localcfg.chain_peak_orientation="local_tangent"
    z=G._chain_model_values(x,y,p,n,axis,localcfg;amp_min=.03,amp_range=.07)
    data=(xs=xs,ys=ys,zimg=reshape(z,length(ys),length(xs)),x=x,y=y,z=z,zfull=z,noise=.02,axisctx=axis)
    return data,ell,circ,localcfg,p
end

@testset "Settings and model-basis consistency" begin
    opt=L.settings(SETTINGS)
    @test opt==(degree=2,initial_maxiter=50,continuation_maxiter=300)
    data,ell,circ,localcfg,p=fixture()
    for cfg in (ell,localcfg)
        basis=VP.design_matrix(data.x,data.y,p,3,data.axisctx,cfg;amp_min=.03,amp_range=.07)
        coeff=VP.linear_coefficients(p,3,data.axisctx,cfg;amp_min=.03,amp_range=.07)
        @test basis*coeff≈G._chain_model_values(data.x,data.y,p,3,data.axisctx,cfg;amp_min=.03,amp_range=.07) atol=2e-16
        r=G.ChainModelResult(n=3,params=p,success=true,amp_min=.03,amp_range=.07)
        VP.finalize_native!(r,data,cfg)
        @test r.gcv≈length(data.z)/(length(data.z)-G._chain_nparams(3,cfg))^2*r.rss
        c=(fit=(result=r,),cfg)
        rows=L.export_rows("sample.sxm",c,data)
        @test length(PatchFrames.read_model_axes(rows))==3
        @test all(row["model_orientation"]==cfg.chain_peak_orientation for row in rows)
        @test all(parse(Int,row["N"])==3 for row in rows)
    end
    @test VP.native_raw_bounds(3,ell)==VP.native_raw_bounds(3,localcfg)
    mktempdir() do dir
        raw=TOML.parsefile(SETTINGS)
        for (section,key,value) in (("model","chain_tangent_degree",true),("model","chain_tangent_degree",3),
            ("model","expected_N",6),("selection","criterion","truth"),("selection","continuation_maxiter",0),
            ("selection","tie_breaker","local"),("preprocessing","patch_orientation","local_tangent"))
            bad=deepcopy(raw); bad[section][key]=value
            path=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException L.settings(path)
        end
    end
    @test_throws ErrorException L.parse_cli(["--expected-N","6"])
    @test_throws ErrorException local_gaussian_comparison(["--truth","forbidden"])
    @test_throws ErrorException local_gaussian_comparison(["--patch-frames","forbidden"])
end

@testset "Matched starts and valid minimum GCV without hidden local fallback" begin
    data,ell,circ,localcfg,p=fixture()
    opt=L.settings(SETTINGS)
    calls=[]
    function fake(n,d,cfg;warm_start=nothing)
        push!(calls,(cfg=deepcopy(cfg),start=deepcopy(warm_start)))
        j=length(calls)
        params=warm_start===nothing ? zeros(G._chain_nparams(n,cfg)) : warm_start .+ .01
        r=G.ChainModelResult(n=n,params=params,success=true,valid=true,rss=1.,gcv=[3.,2.,1.,.5][j],amp_min=.03,amp_range=.07)
        return (result=r,diagnostic=nothing,elapsed_s=0.)
    end
    result=L.fit_candidates(3,data,ell,circ,opt;fitter=fake)
    @test length(calls)==4
    @test result.control.name=="global_elliptical"
    @test result.selected.name=="local_tangent_elliptical"
    @test calls[3].start==calls[4].start==result.candidates[2].fit.result.params
    @test calls[3].cfg.max_iter==calls[4].cfg.max_iter==300
    @test calls[3].cfg.skip_global && calls[4].cfg.skip_global
    @test Set(k for k in fieldnames(typeof(ell)) if !isequal(getfield(calls[3].cfg,k),getfield(calls[4].cfg,k)))==Set([:chain_peak_orientation])
    result.candidates[4].fit.result.gcv=1.
    @test L.best_valid(result.candidates).name=="global_elliptical" # exact ties global
    result.candidates[4].fit.result.gcv=.1; result.candidates[4].fit.result.valid=false
    @test L.best_valid(result.candidates).name=="global_elliptical"
    result.candidates[3].fit.result.valid=false
    @test L.best_valid(result.candidates).name=="initial_elliptical"
    for c in result.candidates; c.fit.result.valid=false; end
    @test_throws ErrorException L.best_valid(result.candidates)
end

@testset "Actual synthetic LM fit and serialization" begin
    data,ell,circ,localcfg,p=fixture()
    for cfg in (ell,localcfg); cfg.skip_global=true; cfg.max_iter=300; end
    a=L.fit_endpoint(3,data,ell;warm_start=copy(p))
    b=L.fit_endpoint(3,data,localcfg;warm_start=copy(p))
    @test a.diagnostic.initial==b.diagnostic.initial
    @test isempty(a.diagnostic.lm_error) && isempty(b.diagnostic.lm_error)
    @test b.result.success && b.result.valid
    @test b.result.rss<a.result.rss
    @test b.result.rss<1e-8
    rows=L.export_rows("sample.sxm",(cfg=localcfg,fit=b),data)
    axes=PatchFrames.read_model_axes(rows)
    @test length(axes)==3
    @test any(abs(v[2])>.05 for v in values(axes))
    @test b.diagnostic.lm_iterations<=300
    @test b.diagnostic.lm_converged == (b.diagnostic.lm_status=="converged")
    println("Synthetic global/local RSS: ",a.result.rss," / ",b.result.rss,
        "; local stopping: ",b.diagnostic.lm_status)
end

@testset "Three complete arms, fresh residuals and strict failures" begin
    mktempdir() do dir
        data,ell,circ,localcfg,p=fixture()
        r=G.ChainModelResult(n=3,params=p,success=true,valid=true,rss=1.,gcv=.01,amp_min=.03,amp_range=.07)
        rows=L.R.feature_rows("sample.sxm",r,data,ell,"ell"); header=sort(collect(keys(rows[1])))
        features=write_table(joinpath(dir,"features.tsv"),header,rows)
        raw=joinpath(dir,"raw"); mkdir(raw); write(joinpath(raw,"sample.sxm"),"fixture, never read as SXM")
        templates=joinpath(dir,"templates.tsv"); write(templates,"fixture")
        out=joinpath(dir,"out")
        args=["--features",features,"--split-features",features,"--data-dir",raw,"--count-config",COUNT,
            "--settings",SETTINGS,"--config",ASSIGNMENT,"--templates",templates,"--outdir",out]
        local_gaussian_comparison(vcat(args,["--dry-run"]))
        @test !ispath(out)
        calls=String[]
        function fake_runner(outdir,name,script,cli)
            @test script=="run_reconstructed_chitosan.jl"
            @test !("--patch-frames" in cli)
            opts=parse_options(cli); push!(calls,name)
            if name!="reference"
                @test_throws ErrorException execute_pipeline(merge(opts,Dict("--config"=>joinpath(ROOT,"config/unit_assignment_reconstructed.toml"),"--patches-fwd"=>features)))
                @test !ispath(opts["--outdir"])
            end
            mkdir(opts["--outdir"]); cp(features,joinpath(opts["--outdir"],"predictions.tsv"))
        end
        function fake_export(outdir,name,script,cli,output;nfiles)
            @test name=="refit" && script=="refit_local_gaussian_orientation.jl" && nfiles==1
            for (mode,cfg) in (("global",ell),("gcv",localcfg))
                rs=L.export_rows("sample.sxm",(fit=(result=r,),cfg),data)
                write_table(output*".$mode.tsv",vcat(header,L.AXIS_FIELDS),rs)
            end
            cp(output*".gcv.tsv",output)
        end
        local_gaussian_comparison(args;runner=fake_runner,exporter=fake_export)
        @test calls==["reference","global_refit","gcv_orientation"]
        for name in calls
            @test read(joinpath(out,name,"features_split.tsv"))==read(features)
        end
        @test !isfile(joinpath(out,"failures.tsv"))
        @test_throws ErrorException local_gaussian_comparison(args;runner=fake_runner,exporter=fake_export)
        bad=copy(args); bad[end]=joinpath(dir,"failed")
        @test_throws ErrorException local_gaussian_comparison(bad;runner=fake_runner,exporter=(a...;kw...)->error("fixture failure"))
        _,failures=read_table(joinpath(dir,"failed/failures.tsv"))
        @test only(failures)["stage"]=="refit"
    end
    @test success(`bash -n $(joinpath(ROOT,"hpc/compare_local_gaussian.sbatch"))`)
end
