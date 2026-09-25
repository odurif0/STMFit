using Test, TOML, LinearAlgebra, Statistics, Random
include(joinpath(@__DIR__,"verify_background_conditioning.jl"))
const V=BackgroundConditioningVerification
const D=V.D
const B=D.B
const G=B.G
const ROOT=dirname(@__DIR__)
const SETTINGS=joinpath(ROOT,"config","background_conditioning.toml")
const OPTIONS=B.settings(SETTINGS)
BLAS.set_num_threads(1)

@testset "Strict scope and exact transformed domain" begin
    @test OPTIONS.options.arms==collect(B.ARMS)
    for bad in ("--expected-N","--labels","--manifest","--truth","--family","--grade")
        @test_throws ErrorException D.options([bad,"anything"])
    end
    mktempdir() do dir
        for (section,key,val) in (("selection","expected_N",6),("model","profile","split"),
                ("selection","case_ranks",[1,2,3,4]),("model","max_time_s",0.),
                ("model","g_tol",1e-3),("preprocessing","support","rebuild"))
            t=deepcopy(OPTIONS.raw); t[section][key]=val; filename=joinpath(dir,"invalid.toml")
            open(io->TOML.print(io,t),filename,"w")
            @test_throws ErrorException B.settings(filename)
        end
    end
    x=repeat(collect(range(2.,4.;length=11));inner=7)
    y=0.7.*x.+repeat(collect(range(4.,4.2;length=7));outer=11)
    z=sin.(x).+cos.(y); p0=[.013,-.003,.002,.7]
    lo=[-5.,-1.,-1.,-2.]; hi=-lo
    problem=(q0=p0,lo,hi,target=z,predict=p->p[1].+p[2].*x.+p[3].*y.+p[4].*sin.(x))
    rng=MersenneTwister(4321)
    for mode in B.ARMS[2:3]
        map=B.coordinates(problem,x,y,mode,OPTIONS.options)
        @test B.from_coordinates(zeros(4),problem,map)==p0
        @test inv(map.M)≈map.inverse
        for code in 0:15
            p=[isodd(code÷2^(j-1)) ? hi[j] : lo[j] for j in 1:4]
            v=B.to_coordinates(p,problem,map)
            @test B.from_coordinates(v,problem,map)≈p atol=1e-12
            @test all(map.lower.-1e-10 .<= v .<= map.upper.+1e-10)
            @test maximum(map.A*v.-map.bconstraint)<=1e-12
        end
        for _ in 1:100
            v=map.lower.+rand(rng,4).*(map.upper.-map.lower)
            p=B.from_coordinates(v,problem,map)
            @test (maximum(map.A*v.-map.bconstraint)<=1e-12)==(B.violation(p,lo,hi)<=1e-12)
        end
        if mode=="centered_plane_slsqp"
            @test map.plane_gram≈Matrix{Float64}(I,3,3)*map.scale^2 rtol=1e-12 atol=1e-14
            # The enclosing coordinate box alone really does enlarge the domain.
            @test any(B.violation(B.from_coordinates([map.upper[1],map.lower[2],map.upper[3],0.],problem,map),lo,hi)>1e-8 for _ in 1:1)
        end
        J=B.jacobian(problem.predict,p0,lo,hi,x,y,OPTIONS.options.derivative_relative_step)
        @test J≈hcat(ones(length(x)),x,y,sin.(x)) rtol=1e-8 atol=1e-8
        Jv=D.C.finite_jacobian(v->problem.predict(B.from_coordinates(v,problem,map)),zeros(4),map.lower,map.upper,
            OPTIONS.options.derivative_relative_step)
        @test J*map.M≈Jv rtol=1e-7 atol=1e-7
    end
    @test_throws ErrorException B.coordinates(problem,x,2x,"centered_plane_slsqp",OPTIONS.options)
    @test_throws ErrorException B.coordinates(merge(problem,(target=ones(length(x)),)),x,y,"centered_plane_slsqp",OPTIONS.options)
end

@testset "Same constrained optimum, not a larger background box" begin
    x=repeat(collect(range(-1.,1.;length=11));inner=7)
    y=repeat(collect(range(-.5,.5;length=7));outer=11)
    A=hcat(ones(length(x)),x,y,sin.(3x).*cos.(y))
    optimum=[5.,.2,-.1,.3]
    problem=(q0=[0.,0.,0.,0.],lo=[-5.,-1.,-1.,-2.],hi=[5.,1.,1.,2.],predict=p->A*p,target=A*optimum.+1.)
    opts=merge(OPTIONS.options,(slsqp_maxeval=1000,max_time_s=10.,evaluation_checkpoints=[1,100,1000]))
    for arm in B.ARMS[2:3]
        fit=B.solve(problem,x,y,arm,opts,OPTIONS.audit)
        @test fit.params≈optimum atol=1e-6
        @test fit.rss≈length(x) atol=1e-8
        @test fit.raw_endpoint_violation<=opts.coordinate_roundoff_tolerance
        @test D.C.stationarity(fit,OPTIONS.audit).passed
        @test fit.diagnostic.initial==problem.q0
        @test fit.evaluations==length(fit.history)
        @test fit.model_evaluations>fit.evaluations && fit.gradient_evaluations>0
        @test fit.rss<=fit.initial_rss
        stopped=B.solve(problem,x,y,arm,merge(opts,(slsqp_maxeval=1,evaluation_checkpoints=[1])),OPTIONS.audit)
        @test stopped.status=="MAXEVAL_REACHED" && !D.C.stationarity(stopped,OPTIONS.audit).passed
    end
end

function fixture(dir)
    source=joinpath(dir,"source"); out=joinpath(dir,"run"); mkpath(source); mkdir(out)
    physical=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
    count_settings=D.L.settings(joinpath(ROOT,"config","local_sigma_counting.toml"))
    cfg=D.L.configs(physical,count_settings,"global_sigma_max").ell
    N=3; axis=(origin=(3.,4.),axis=[.8,.6],perp=[-.6,.8],tmin=-1.5,tmax=1.5)
    t=repeat(collect(range(-1.7,1.7;length=21));inner=7)
    u=repeat(collect(range(-.15,.15;length=7));outer=21)
    x=3 .+.8t.-.6u; y=4 .+.6t.+.8u
    p=zeros(3+5N); p[1:3]=[.012,.002,-.003]
    amin=.03; arange=.1
    pred=G._chain_model_values(x,y,p,N,axis,cfg;amp_min=amin,amp_range=arange)
    z=pred.+.0001sin.(collect(eachindex(x)))
    rss=sum(abs2,z.-pred); gcv=length(z)/(length(z)-length(p))^2*rss
    for chunk in 1:4
        folder=joinpath(source,"chunk$chunk"); mkpath(joinpath(folder,"fit_data"))
        open(io->TOML.print(io,physical),joinpath(folder,"count.toml"),"w")
        open(io->TOML.print(io,count_settings),joinpath(folder,"settings.toml"),"w")
        snapshot=Dict(string(k)=>(getfield(cfg,k) isa Symbol ? string(getfield(cfg,k)) : getfield(cfg,k)) for k in fieldnames(typeof(cfg)))
        open(io->TOML.print(io,snapshot),joinpath(folder,"global_sigma_max.ell.toml"),"w")
        cohort=Dict{String,String}[]; selected=Dict{String,String}[]; candidates=Dict{String,String}[]; contexts=Dict{String,String}[]
        for j in chunk:4:146
            file="synthetic_$(lpad(j,3,'0')).sxm"
            push!(cohort,Dict("file"=>file,"raw_sha256"=>repeat("0",64)))
            push!(contexts,Dict("file"=>file,"raw_sha256"=>repeat("0",64),"n_data"=>string(length(z)),"noise"=>"0.01",
                "origin_x_nm"=>"3.0","origin_y_nm"=>"4.0","axis_x"=>"0.8","axis_y"=>"0.6","support_tmin_nm"=>"-1.5","support_tmax_nm"=>"1.5"))
            for rep in 1:2
                push!(selected,Dict("file"=>file,"repetition"=>string(rep),"arm"=>"global_sigma_max","status"=>"ok","source"=>"ell",
                    "N_selected"=>string(N),"gcv"=>string(gcv*(rep==1 ? 1. : 1. + j*1e-6))))
            end
            push!(candidates,Dict("file"=>file,"repetition"=>"1","arm"=>"global_sigma_max","N"=>string(N),"family"=>"ell",
                "success"=>"true","valid"=>"true","gcv"=>string(gcv),"rss"=>string(rss),"amp_min"=>string(amin),"amp_range"=>string(arange),"parameters"=>join(p,';')))
            D.write_rows(joinpath(folder,"fit_data",file*".tsv"),[Dict("x_nm"=>string(x[k]),"y_nm"=>string(y[k]),"z_nm"=>string(z[k])) for k in eachindex(z)])
        end
        for (name,rs) in (("cohort",cohort),("selected",selected),("candidates",candidates),("context",contexts))
            D.write_rows(joinpath(folder,"$name.tsv"),rs)
        end
    end
    settings=deepcopy(OPTIONS.raw)
    settings["model"]["lm_maxiter"]=3; settings["model"]["slsqp_maxeval"]=5
    settings["model"]["max_time_s"]=10.; settings["model"]["checkpoints"]=[0,3]; settings["model"]["evaluation_checkpoints"]=[1,5]
    filename=joinpath(dir,"test_settings.toml"); open(io->TOML.print(io,settings),filename,"w")
    (;source,out,settings=filename)
end

@testset "Synthetic full saved-input workflow and independent verifier" begin
    mktempdir() do dir
        f=fixture(dir)
        opts=Dict("--input-root"=>f.source,"--settings"=>f.settings,"--outdir"=>joinpath(f.out,"unused"))
        input=D.inputs(opts)
        @test [D.integer(c.ranking,"rank") for c in input.cases]==[1,49,98,146]
        @test length(input.ranking)==146
        D.execute(merge(opts,Dict("--dry-run"=>"true")))
        @test !ispath(opts["--outdir"])
        ctx=D.context(first(input.cases)); case=first(input.cases)
        for p in (ctx.problem.q0,ctx.problem.lo,ctx.problem.hi,(ctx.problem.lo.+ctx.problem.hi)./3)
            @test V.independent(p,case,ctx).pred≈ctx.problem.predict(p) rtol=1e-12 atol=1e-12
        end
        for chunk in 1:4
            D.execute(merge(opts,Dict("--outdir"=>joinpath(f.out,"chunk$chunk"),"--chunk"=>"$chunk/4")))
        end
        result=V.verify(f.source,f.out;settings=f.settings)
        @test length(result.fits)==24 && length(result.repeatability)==12
        @test all(r["same_validity"]=="true" for r in result.repeatability)
        @test all(D.number(r,"rss")<=D.number(r,"initial_rss")+1e-12 for r in result.fits)
        @test_throws ErrorException D.options(["--input-root",f.source,"--settings",f.settings,"--outdir",f.out])
    end
end
