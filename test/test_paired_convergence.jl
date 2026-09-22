using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_paired_convergence.jl"))
module Merger
include(joinpath(@__DIR__,"merge_paired_convergence.jl"))
end
include(joinpath(@__DIR__,"verify_paired_convergence.jl"))
const D=PairedConvergenceDiagnostic
const C=D.C
const P=D.P
const G=D.G
const ROOT=dirname(@__DIR__)
const SETTINGS=joinpath(ROOT,"config","paired_convergence.toml")
BLAS.set_num_threads(1)

@testset "Independent derivatives and box stationarity at interior and active bounds" begin
    f=q->[exp(q[1]),q[1]*q[2],q[2]^2]
    lo=[-1.,0.]; hi=[2.,3.]
    for q in ([-1.,3.],[.2,.7],[2.,0.]),step in (1e-5,5e-6)
        J=C.finite_jacobian(f,q,lo,hi,step)
        @test J≈[exp(q[1]) 0.; q[2] q[1]; 0. 2q[2]] rtol=1e-8 atol=1e-8
    end
    @test C.box_gradient([-1.,3.],[1.,-1.],lo,hi,1.).norm==0
    @test C.box_gradient([-1.,3.],[-1.,1.],lo,hi,1.).norm>0
    @test C.box_gradient([.2,.7],zeros(2),lo,hi,1.).norm==0
    @test_throws ErrorException C.finite_jacobian(f,[-2.,3.],lo,hi,1e-5)
end

@testset "Strict convergence settings: budgets not physics or label controls" begin
    raw=TOML.parsefile(SETTINGS)
    @test C.settings(SETTINGS).extended_maxiter==10000
    mktempdir() do dir
        path=joinpath(dir,"bad.toml")
        for (section,key,value) in (("selection","expected_N",6),("selection","count_policy","reselect"),
                ("model","x_tol",1e-4),("model","g_tol",1e-4),("model","extended_maxiter",300),
                ("model","finite_difference_half_step_factor",1.),("model","checkpoints",[0,10000]))
            bad=deepcopy(raw); bad[section][key]=value
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException C.settings(path)
        end
    end
    @test_throws ErrorException D.parse_cli(["--labels","truth"])
end

@testset "Saved-pixel runner: complete cases, identical trace prefixes, no raw scan" begin
    mktempdir() do dir
        native=joinpath(dir,"native"); reference=joinpath(dir,"reference"); mkpath(native); mkpath(reference)
        physical=joinpath(ROOT,"config","chitosan.toml"); assign=joinpath(ROOT,"config","unit_assignment_patch_support.toml")
        modelsettings=TOML.parsefile(joinpath(ROOT,"config","paired_acquisition.toml")); modelsettings["model"]["local_maxiter"]=3
        modelpath=joinpath(dir,"model.toml"); open(io->TOML.print(io,modelsettings),modelpath,"w")
        cp(modelpath,joinpath(reference,"paired_settings.toml"))
        settings=TOML.parsefile(SETTINGS); settings["model"]["control_maxiter"]=3; settings["model"]["extended_maxiter"]=6
        settings["model"]["checkpoints"]=[0,3,6]
        configpath=joinpath(dir,"convergence.toml"); open(io->TOML.print(io,settings),configpath,"w")
        raw=TOML.parsefile(physical); file=settings["selection"]["file"]; n=3
        _,ec,_=D.F.Extractor._configs(raw["model"],raw["preprocessing"],dir)
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
        xs=collect(range(-2.,2.;length=21)); ys=collect(range(-.5,.5;length=7))
        x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
        p=zeros(G._chain_nparams(n,ec)); p[1]=.01
        z=G._chain_model_values(x,y,p,n,axis,ec;amp_min=.03,amp_range=.07).+0.0001 .* sin.(eachindex(x))
        pixels=[Dict("row"=>string(mod1(j,length(ys))),"column"=>string(cld(j,length(ys))),
            "x_nm"=>D.F.fmt(x[j]),"y_nm"=>D.F.fmt(y[j]),"z_nm"=>D.F.fmt(z[j]),
            "fwd_nm"=>D.F.fmt(z[j]+.001),"bwd_nm"=>D.F.fmt(z[j]-.001)) for j in eachindex(x)]
        D.write_records(joinpath(native,"registered.pixels.tsv"),pixels)
        fits=Dict{String,String}[]; parameters=Dict{String,String}[]; prev=Dict{String,String}[]; prevp=Dict{String,String}[]
        for profile in ("gaussian","split"),family in ("circ","ell")
            cfg=deepcopy(ec); cfg.peak_profile=Symbol(profile); cfg.chain_circular_sigmas=family=="circ"
            pp=zeros(G._chain_nparams(n,cfg)); pp[1]=.01
            for arm in ("control","registered")
                push!(fits,Dict("file"=>file,"N"=>string(n),"arm"=>arm,"profile"=>profile,"family"=>family,"success"=>"true",
                    "axis_x"=>"1.0","axis_y"=>"0.0","origin_x_nm"=>"0.0","origin_y_nm"=>"0.0",
                    "support_tmin"=>"-1.5","support_tmax"=>"1.5","amp_min"=>"0.03","amp_range"=>"0.07",
                    "noise"=>"0.01","offset"=>"0.0","n_pixels"=>string(length(z))))
                for j in eachindex(pp)
                    push!(parameters,Dict("arm"=>arm,"profile"=>profile,"family"=>family,"stage"=>"final","start"=>"1","parameter"=>string(j),"value"=>D.F.fmt(pp[j])))
                end
            end
            for mode in ("fused","paired")
                push!(prev,Dict("file"=>file,"profile"=>profile,"family"=>family,"mode"=>mode,"rss"=>"0.1"))
                for j in 1:(length(pp)+(mode=="paired" ? 4 : 0))
                    push!(prevp,Dict("profile"=>profile,"family"=>family,"mode"=>mode,"parameter"=>string(j),"value"=>D.F.fmt(j<=length(pp) ? pp[j] : 0.)))
                end
            end
        end
        for (name,rows) in (("fits",fits),("parameters",parameters)); D.write_records(joinpath(native,"$name.tsv"),rows); end
        D.write_records(joinpath(native,"input_hashes.tsv"),[Dict("input"=>key,"sha256"=>bytes2hex(sha256(read(path)))) for (key,path) in (("--config",physical),("--assignment-config",assign))])
        D.write_records(joinpath(reference,"native_hashes.tsv"),[Dict("input"=>f,"sha256"=>bytes2hex(sha256(read(joinpath(native,f))))) for f in ("fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv")])
        D.write_records(joinpath(reference,"fits.tsv"),prev); D.write_records(joinpath(reference,"parameters.tsv"),prevp)
        args=["--native-dir",native,"--reference-dir",reference,"--physical-config",physical,"--assignment-config",assign,
            "--model-settings",modelpath,"--settings",configpath,"--outdir",joinpath(dir,"out")]
        dry=D.parse_cli(vcat(args,["--dry-run"])); D.execute(dry)
        @test !ispath(joinpath(dir,"out"))
        opts=D.parse_cli(args); D.execute(opts)
        actual=D.table(joinpath(dir,"out","fits.tsv"))
        @test length(actual)==16
        @test !isfile(joinpath(dir,"out","failures.tsv"))
        @test all(r["N"]=="3" && r["threshold"]=="3.5" for r in actual)
        for case in D.CASES
            prefix="$(case.profile).$(case.family).$(case.mode)"
            tr=D.table(joinpath(dir,"out",prefix*".trace.tsv"))
            ct=filter(r->r["stage"]=="control",tr); et=filter(r->r["stage"]=="extended",tr)
            @test length(et)>=length(ct)
            @test getindex.(ct,"rss")==getindex.(et[1:length(ct)],"rss")
            @test isfile(joinpath(dir,"out",prefix*".gradients.tsv"))
            @test isfile(joinpath(dir,"out",prefix*".extended.residuals.tsv"))
        end
        @test_throws ErrorException D.parse_cli(args)
        # Merely adding diagnostic exports and explicit native tolerances must
        # not change the short-run numerical result.
        input=D.inputs(opts); ctx=D.context(input,first(D.CASES))
        plain=P.fit(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,"fused",input.model_options)
        observed=P.fit(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,"fused",input.model_options;
            optimizer_options=(maxTime=600.,x_tol=1e-8,g_tol=1e-12),diagnostic=true)
        @test plain.params==observed.params && plain.rss==observed.rss
        @test plain.diagnostic===nothing
        @test_throws ErrorException P.fit(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,"fused",input.model_options;optimizer_options=(expected_N=3,))
        merged=joinpath(dir,"merged"); mkpath(merged)
        for i in 1:4
            shard=joinpath(merged,"chunk$i"); mkpath(shard)
            cases=Set((c.profile,c.family,c.mode) for c in D.CASES[i:4:end])
            D.write_records(joinpath(shard,"fits.tsv"),filter(r->(r["profile"],r["family"],r["mode"]) in cases,actual))
            for file in ("input_hashes.tsv","settings.toml","model_settings.toml")
                cp(joinpath(dir,"out",file),joinpath(shard,file))
            end
            for case in D.CASES[i:4:end]
                prefix="$(case.profile).$(case.family).$(case.mode)."
                for file in filter(f->startswith(f,prefix),readdir(joinpath(dir,"out")))
                    cp(joinpath(dir,"out",file),joinpath(shard,file))
                end
            end
        end
        Merger.main([merged])
        @test length(D.table(joinpath(merged,"fits.tsv")))==16
        @test length(D.table(joinpath(merged,"eligibility.tsv")))==8
        @test_throws Exception Merger.main([merged])
        VerifyPairedConvergence.verify(native,reference,merged)
        launcher=joinpath(ROOT,"hpc","diagnose_paired_convergence.sbatch")
        @test success(`bash -n $launcher`)
        @test occursin("SLURM_JOB_ID:?",read(launcher,String))
    end
end
