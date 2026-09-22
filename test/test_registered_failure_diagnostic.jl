using Test, TOML, LinearAlgebra
using STMSXMIO: SXMImage, SXMChannel
include(joinpath(@__DIR__,"diagnose_registered_fit_failure.jl"))
const D=RegisteredFailureDiagnostic
const G=D.G
const ROOT=dirname(@__DIR__)
BLAS.set_num_threads(1)

@testset "Native optimizer observer changes no fixed-start fit or validity" begin
    raw=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
    for profile in (:gaussian,:split), circular in (false,true)
        _,cfg,_=D.F.Extractor._configs(raw["model"],raw["preprocessing"],"unused")
        cfg.peak_profile=profile; cfg.chain_circular_sigmas=circular
        cfg.skip_global=true; cfg.max_iter=5
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
        xs=collect(range(-2.,2.;length=25)); ys=collect(range(-.5,.5;length=9))
        x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
        p=zeros(G._chain_nparams(3,cfg)); p[1]=.01
        z=G._chain_model_values(x,y,p,3,axis,cfg;amp_min=.03,amp_range=.07)
        zimg=reshape(z,length(ys),length(xs))
        data=(;xs,ys,x,y,z,zimg,zfull=z,noise=.01,axisctx=axis)
        plain=G._fit_chain_n(xs,ys,zimg,x,y,z,.01,3,axis,cfg;warm_start=p)
        rows=NamedTuple[]
        observed=G._fit_chain_n(xs,ys,zimg,x,y,z,.01,3,axis,cfg;warm_start=p,diagnostics=r->push!(rows,r))
        @test length(rows)==1
        @test plain.params==observed.params==rows[1].params
        @test isequal(plain.param_perr,observed.param_perr)
        @test plain.amp_min==observed.amp_min && plain.amp_range==observed.amp_range
        @test rows[1].global_status=="skipped" && rows[1].global_evaluations==0
        @test isempty(rows[1].lm_error)
        @test rows[1].lm_status in ("converged","not_converged")
        @test 0<=rows[1].lm_iterations<=cfg.max_iter
        D.VP.finalize_native!(plain,data,cfg); D.VP.finalize_native!(observed,data,cfg)
        @test plain.valid==observed.valid && plain.reason==observed.reason && plain.rss==observed.rss && plain.gcv==observed.gcv
        # Even a mutating observer cannot alter the result.
        mutated=G._fit_chain_n(xs,ys,zimg,x,y,z,.01,3,axis,cfg;warm_start=p,
            diagnostics=r->(r.params.=NaN; r.initial.=NaN; r.global_params.=NaN))
        @test mutated.params==plain.params
    end
end

@testset "One-file diagnostic writes every stage once, then complete summaries" begin
    mktempdir() do dir
        raw=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
        raw["model"]["global_maxtime"]=0.1; raw["model"]["global_maxiter"]=10; raw["model"]["max_iter"]=10
        raw["preprocessing"]["flatten"]="none"; raw["preprocessing"]["smooth_radius_px"]=0
        f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
        img=SXMImage("synthetic.sxm",Dict{String,String}(),65,41,(5.,2.),(0.,0.),
            [SXMChannel("Z","nm","fwd",f),SXMChannel("Z","nm","bwd",copy(f))])
        pcfg,ec,_=D.F.Extractor._configs(raw["model"],raw["preprocessing"],dir)
        original=D.RR.original_support(img,pcfg,ec)
        p=zeros(G._chain_nparams(3,ec)); p[1]=.01
        r=G.ChainModelResult(n=3,params=p,amp_min=.03,amp_range=.07,gcv=.1)
        rows=D.R.feature_rows("synthetic.sxm",r,(axisctx=original.axis,),ec,"ell")
        base=joinpath(dir,"base.tsv"); D.write_table(base,sort(collect(keys(first(rows)))),rows)
        shifts=joinpath(dir,"shifts.tsv"); D.write_table(shifts,["file","bwd_sample_dx_px"],
            [Dict("file"=>"synthetic.sxm","bwd_sample_dx_px"=>"0")])
        cfgpath=joinpath(dir,"physical.toml"); open(io->TOML.print(io,raw),cfgpath,"w")
        touch(joinpath(dir,"synthetic.sxm"))
        opts=D.parse_cli(["--file","synthetic.sxm","--features",base,"--shifts",shifts,"--data-dir",dir,
            "--config",cfgpath,"--assignment-config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",joinpath(ROOT,"config","registered_refit_original_support.toml"),"--outdir",joinpath(dir,"out")])
        D.execute(opts;reader=path->img)
        _,fits=D.read_table(joinpath(dir,"out","fits.tsv"))
        _,optimizers=D.read_table(joinpath(dir,"out","optimizers.tsv"))
        @test length(fits)==length(optimizers)==8
        @test all(r["N"]=="3" for r in fits)
        for arm in ("control","registered"),profile in ("gaussian","split"),family in ("circ","ell")
            @test isfile(joinpath(dir,"out","$arm.$profile.$family.fits.tsv"))
            @test isfile(joinpath(dir,"out","$arm.$profile.$family.parameters.tsv"))
            @test isfile(joinpath(dir,"out","$arm.$profile.$family.residuals.tsv"))
        end
        @test_throws ErrorException D.parse_cli(["--file","synthetic.sxm","--features",base,"--shifts",shifts,"--data-dir",dir,
            "--config",cfgpath,"--assignment-config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",joinpath(ROOT,"config","registered_refit_original_support.toml"),"--outdir",joinpath(dir,"out")])
    end
end

@testset "Diagnostic CLI rejects labels and missing inputs" begin
    @test D.parse_cli(["--help"])===nothing
    for arg in ("--expected-n","--labels","--sequence","--threshold","--retry")
        @test_throws ErrorException D.parse_cli([arg,"6"])
    end
end
