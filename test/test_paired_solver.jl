using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_paired_solver.jl"))
const D=PairedSolverDiagnostic
const S=D.S
const C=S.C
const P=S.P
const G=P.G
const ROOT=dirname(@__DIR__)
const CONFIG=joinpath(ROOT,"config","paired_solver.toml")
const AUDIT=joinpath(ROOT,"config","paired_convergence.toml")
BLAS.set_num_threads(1)

@testset "Strict label-free solver settings and affine box" begin
    raw=TOML.parsefile(CONFIG); opts=S.settings(CONFIG)
    @test opts.method=="unit_box_slsqp" && opts.lm_maxiter==10000
    mktempdir() do dir
        path=joinpath(dir,"bad.toml")
        for (section,key,val) in (("selection","expected_N",6),("model","labels","NKNNKN"),("model","method","try_every_solver"),
                ("model","slsqp_maxeval",2.5),("model","max_time_s",0.),("model","evaluation_checkpoints",[1,100,50]))
            bad=deepcopy(raw); bad[section][key]=val
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException S.settings(path)
        end
    end
    lo=[-5.,-.1]; hi=[5.,.3]; q0=[.123456789,.05]; u0=(q0.-lo)./(hi.-lo)
    @test S.from_unit(u0,q0,lo,hi,opts.coordinate_roundoff_tolerance)==q0
    @test S.from_unit(zeros(2),q0,lo,hi,opts.coordinate_roundoff_tolerance)≈lo atol=1e-15
    @test S.from_unit(ones(2),q0,lo,hi,opts.coordinate_roundoff_tolerance)≈hi atol=1e-15
    @test_throws ErrorException S.from_unit([-0.1,.5],q0,lo,hi,opts.coordinate_roundoff_tolerance)
    @test_throws ErrorException D.parse_cli(["--labels","truth"])
end

@testset "Coupled constrained least squares: active optimum and independent stationarity" begin
    opts=merge(S.settings(CONFIG),(slsqp_maxeval=500,evaluation_checkpoints=[1,100,500]))
    audit=C.settings(AUDIT)
    p=(q0=[1.,1.],lo=[-1.,-2.],hi=[2.,3.],predict=q->[4q[1]+q[2],q[2]],target=[10.,-1.])
    before=copy(p.q0); fit=S.solve(p,opts,audit); a=C.stationarity(fit,audit)
    @test fit.params≈[2.,.5] atol=1e-6
    @test fit.rss≈4.5 atol=1e-10
    @test fit.evaluations==length(fit.history)<=opts.slsqp_maxeval
    @test 0<fit.gradient_evaluations<=fit.evaluations
    @test fit.model_evaluations>fit.evaluations
    @test all(diff([r.best_rss for r in fit.history]).<=0)
    @test a.passed && a.agreement
    @test p.q0==before && fit.diagnostic.initial==p.q0
    @test all(p.lo.<=fit.params.<=p.hi)
    @test fit.status in ("FTOL_REACHED","XTOL_REACHED","SUCCESS")
    for (index,q) in fit.checkpoints
        index==0 && continue
        @test sum(abs2,p.predict(q).-p.target)≈fit.history[index].rss atol=1e-12
    end
    stopped=S.solve(p,merge(opts,(slsqp_maxeval=1,evaluation_checkpoints=[1])),audit)
    @test stopped.status=="MAXEVAL_REACHED" && stopped.evaluations==1
    @test !C.stationarity(stopped,audit).passed
end

@testset "Same molecular objectives, bounds and validity for both modes and profiles" begin
    opts=merge(S.settings(CONFIG),(slsqp_maxeval=1000,evaluation_checkpoints=[1,100,1000]))
    audit=C.settings(AUDIT); model=P.settings(joinpath(ROOT,"config","paired_acquisition.toml"))
    for profile in (:gaussian,:split), circular in (false,true)
        cfg=G.ChainSweepConfig(peak_profile=profile,chain_circular_sigmas=circular,chain_tilted_baseline=true)
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
        xs=collect(range(-2.,2.;length=21)); ys=collect(range(-.5,.5;length=7))
        x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
        p=zeros(G._chain_nparams(3,cfg)); p[1:3]=[.01,.001,-.002]
        pred=P.predictions(vcat(p,[.1,.003,.002,-.001]),3,x,y,axis,cfg,.03,.07)
        z=pred.mean.+0.0001 .* sin.(eachindex(x)); f=pred.fwd.+(z.-pred.mean); b=pred.bwd.+(z.-pred.mean)
        data=(;x,y,z,zfull=z,noise=.01,axisctx=axis,xs,ys,zimg=reshape(z,length(ys),length(xs)))
        initial=G.ChainModelResult(n=3,params=p,success=true,amp_min=.03,amp_range=.07)
        for mode in ("fused","paired")
            pr=S.problem(initial,data,cfg,f,b,mode,model)
            lm=P.fit(initial,data,cfg,f,b,mode,merge(model,(local_maxiter=3,));diagnostic=true)
            @test pr.q0==lm.diagnostic.initial && pr.lo==lm.lower && pr.hi==lm.upper
            @test pr.target==lm.diagnostic.target
            @test pr.predict(pr.q0)==lm.diagnostic.predict(pr.q0)
            sl=S.finish(S.solve(pr,opts,audit),initial,data,cfg,f,b,mode)
            @test sl.rss<=sl.initial_rss+1e-12
            @test sl.result.n==3 && sl.np==length(pr.q0) && sl.nd==length(pr.target)
            @test sl.gcv==sl.nd/(sl.nd-sl.np)^2*sl.rss
            @test sl.valid && all(sl.lower.<=sl.params.<=sl.upper)
            @test sl.predictions.mean≈(sl.predictions.fwd.+sl.predictions.bwd)./2 atol=1e-15
            @test initial.params==p
            strict=deepcopy(cfg); strict.residual_peak_snr_threshold=0.
            rejected=S.finish(sl,initial,data,strict,f,b,mode)
            @test !rejected.valid && occursin("residual high",rejected.reason)
            @test C.stationarity(sl,audit).agreement
        end
        @test_throws ErrorException S.problem(initial,data,cfg,fill(NaN,length(f)),b,"paired",model)
    end
end
