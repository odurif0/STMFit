using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_paired_shift.jl"))
const D=PairedShiftDiagnostic
const H=D.H
const G=H.G
const ROOT=dirname(@__DIR__)
const SETTINGS=joinpath(ROOT,"config","paired_shift.toml")
const OPT=H.settings(SETTINGS)
const AUDIT=H.C.settings(joinpath(ROOT,"config","paired_convergence.toml"))
const SOLVER=merge(H.S.settings(joinpath(ROOT,"config","paired_solver.toml")),(slsqp_maxeval=1000,max_time_s=30.,evaluation_checkpoints=[1,1000]))
const MODEL=H.P.settings(joinpath(ROOT,"config","paired_acquisition.toml"))
BLAS.set_num_threads(1)

@testset "Explicit label-free settings, deterministic bounded starts, spatial splits" begin
    @test OPT.start_seeds==[0,11,29,47]
    mktempdir() do dir
        for (section,key,v) in (("selection","expected_N",6),("model","labels","NKNNKN"),("model","shift_max_px",0.),
            ("model","start_seeds",[0,0]),("preprocessing","holdout_buffer_px",8),("selection","fold_shift_agreement_px",2.))
            t=TOML.parsefile(SETTINGS); t[section][key]=v; path=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,t),path,"w")
            @test_throws ErrorException H.settings(path)
        end
    end
    @test_throws ErrorException D.parse_cli(["--labels","forbidden"])
    q0=[.01,.002,-.001,.4,-.7,0.,0.,0.,0.,0.,0.]
    p=(;q0,lo=fill(-1.,11),hi=fill(1.,11),molecular_count=5)
    @test H.start(p,0,.02)==q0
    for seed in OPT.start_seeds[2:end]
        q=H.start(p,seed,.02)
        @test q==H.start(p,seed,.02) && q!=q0
        @test q[1:3]==q0[1:3] && q[6:end]==q0[6:end]
        @test maximum(abs.(q.-q0))<=.04
        @test all(p.lo.<=q.<=p.hi)
    end
    rows=repeat(collect(1:40);inner=40); cols=repeat(collect(1:40);outer=40)
    fs=H.folds(rows,cols,16,2)
    @test sort(vcat(fs[1].test,fs[2].test))==collect(eachindex(rows))
    for f in fs
        @test isempty(intersect(f.train,f.test)) && !isempty(f.train)
        for i in f.train
            @test minimum(max(abs(rows[i]-rows[j]),abs(cols[i]-cols[j])) for j in f.test)>2
        end
    end
    pixels=[Dict("x_nm"=>string(.04c),"y_nm"=>string(.05r),"row"=>string(r),"column"=>string(c)) for (r,c) in zip(rows,cols)]
    @test all(isapprox.(H.grid_step(pixels),(.04,.05)))
    pixels[10]["x_nm"]="2.123"
    @test_throws ErrorException H.grid_step(pixels)
    # A shared positive normalization leaves the solution, not the recorded initial RSS, invariant.
    simple=(q0=[0.],lo=[-1.],hi=[1.],predict=q->[q[1]],target=[.3])
    fit=H.S.solve(simple,SOLVER,AUDIT;objective_scale=2.)
    @test fit.params≈[.3] atol=1e-7
    @test fit.initial_rss≈.09
    @test_throws ErrorException H.S.solve(simple,SOLVER,AUDIT;objective_scale=0.)
end

function fixture(profile,circular;shift=(.37,-.28))
    cfg=G.ChainSweepConfig(peak_profile=profile,chain_circular_sigmas=circular,chain_tilted_baseline=true)
    n=2; axis=(origin=(0.,0.),axis=[.8,.6],perp=[-.6,.8],tmin=-.8,tmax=.8)
    xs=collect(range(-1.2,1.2;length=33)); ys=collect(range(-.8,.8;length=25))
    x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
    step=(xs[2]-xs[1],ys[2]-ys[1]); p=zeros(G._chain_nparams(n,cfg)); p[1:3]=[.01,.001,-.002]
    q=vcat(p,[.07,.0003,.0001,-.0002,shift...])
    pred=H.predictions(q,n,x,y,axis,cfg,.03,.07,step)
    z=pred.mean; data=(;x,y,z,zfull=z,noise=.01,axisctx=axis,xs,ys,zimg=reshape(z,length(ys),length(xs)))
    initial=G.ChainModelResult(n=n,params=p,success=true,amp_min=.03,amp_range=.07)
    return ((;initial,data,cfg,fwd=pred.fwd,bwd=pred.bwd),step,q)
end

@testset "Exact zero-shift nesting and joint recovery on synthetic data" begin
    for profile in (:gaussian,:split),circular in (true,false)
        ctx,step,q=fixture(profile,circular); k=length(ctx.initial.params)
        z=copy(q); z[end-1:end].=0
        native=H.P.predictions(z[1:k+4],ctx.initial.n,ctx.data.x,ctx.data.y,ctx.data.axisctx,ctx.cfg,.03,.07)
        @test H.predictions(z,ctx.initial.n,ctx.data.x,ctx.data.y,ctx.data.axisctx,ctx.cfg,.03,.07,step)==native
        free=H.problem(ctx,"shifted",MODEL,OPT,step); control=H.problem(ctx,"paired",MODEL,OPT,step)
        @test free.predict(free.q0)==control.predict(control.q0)
        @test free.target==control.target
        @test free.lo[1:k+4]==control.lo && free.hi[1:k+4]==control.hi
        @test free.lo[end-1:end]==[-1.,-1.] && free.hi[end-1:end]==[1.,1.]
        J=H.C.finite_jacobian(free.predict,q,free.lo,free.hi,AUDIT.finite_difference_relative_step)
        @test all(iszero,J[1:length(ctx.fwd),end-1:end])
        # Fixed finite differences independently check translation derivatives.
        qa=copy(q); qb=copy(q); qa[end]+=.0001; qb[end]-=.0001
        @test J[:,end]≈(free.predict(qa).-free.predict(qb))./.0002 rtol=1e-5 atol=1e-9
        fit=H.finish(H.S.solve(free,SOLVER,AUDIT),free,ctx,"shifted")
        @test fit.params[end-1:end]≈q[end-1:end] atol=0.002
        @test fit.rss<1e-8 && fit.valid
        @test fit.np==k+6 && fit.nd==2length(ctx.data.z)
        @test fit.gcv==fit.nd/(fit.nd-fit.np)^2*fit.rss
        @test fit.result.rss≈sum(abs2,ctx.data.z.-fit.predictions.mean)
        @test ctx.initial.params==q[1:k]
        strict=merge(ctx,(cfg=deepcopy(ctx.cfg),)); strict.cfg.residual_peak_snr_threshold=0.
        @test !H.finish(fit,free,strict,"shifted").valid
        f=H.problem(ctx,"shifted",MODEL,OPT,step;indices=1:2:length(ctx.data.z))
        altered=merge(ctx,(fwd=copy(ctx.fwd),bwd=copy(ctx.bwd)))
        altered.fwd[2:2:end].+=10; altered.bwd[2:2:end].-=10
        other=H.problem(altered,"shifted",MODEL,OPT,step;indices=1:2:length(ctx.data.z))
        @test f.target==other.target && f.predict(q)==other.predict(q)
    end
end

@testset "Selection ignores heldout outcomes until the fixed pilot check" begin
    rows=Dict{String,String}[]
    for c in D.CASES,m in ("fused","paired","shifted"),f in (m=="fused" ? (0,) : (0,1,2)),s in OPT.start_seeds
        push!(rows,Dict("id"=>"$(c.profile).$(c.family).$m.$f.$s","profile"=>c.profile,"family"=>c.family,"mode"=>m,"fold"=>string(f),"seed"=>string(s),
            "rss"=>string(1+s),"valid"=>"true","gcv"=>m=="shifted" ? "0.5" : "1.0",
            "heldout_rss"=>m=="shifted" ? "0.5" : "1.0","shift_x_px"=>"0.3","shift_y_px"=>"-0.2"))
    end
    summary=D.summaries(rows,OPT)
    @test length(summary.winners)==28 && all(r["seed"]=="0" for r in summary.winners)
    @test all(r["eligible"]=="true" for r in summary.pilots)
    altered=deepcopy(rows)
    for r in altered
        r["mode"]=="shifted" && r["fold"]=="1" && (r["heldout_rss"]="1.01")
    end
    other=D.summaries(altered,OPT)
    @test getindex.(other.winners,"id")==getindex.(summary.winners,"id")
    @test all(r["eligible"]=="false" for r in other.pilots)
    @test_throws ErrorException D.summaries(rows[2:end],OPT)
    launcher=joinpath(ROOT,"hpc","diagnose_paired_shift.sbatch")
    @test success(`bash -n $launcher`)
    @test occursin("SLURM_JOB_ID:?",read(launcher,String))
end
