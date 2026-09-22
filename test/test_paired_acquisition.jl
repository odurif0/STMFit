using Test, TOML, LinearAlgebra
include(joinpath(@__DIR__,"lib","paired_acquisition_fit.jl"))
const P=PairedAcquisitionFit
const G=P.G
const SETTINGS=joinpath(dirname(@__DIR__),"config","paired_acquisition.toml")
BLAS.set_num_threads(1)

@testset "Identifiable paired acquisition and unchanged mean geometry" begin
    options=P.settings(SETTINGS)
    for profile in (:gaussian,:split), circular in (false,true)
        cfg=G.ChainSweepConfig(peak_profile=profile,chain_circular_sigmas=circular,chain_tilted_baseline=true)
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
        xs=collect(range(-2.,2.;length=31)); ys=collect(range(-.5,.5;length=11))
        x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
        p=zeros(G._chain_nparams(3,cfg)); p[1:3]=[.01,.001,-.002]
        q=vcat(p,[.2,.003,.002,-.001])
        pred=P.predictions(q,3,x,y,axis,cfg,.03,.07)
        mean=G._chain_model_values(x,y,p,3,axis,cfg;amp_min=.03,amp_range=.07)
        @test pred.mean==mean
        @test (pred.fwd.+pred.bwd)./2≈mean rtol=1e-14
        @test pred.fwd.-pred.bwd≈2 .* (0.2 .* pred.lobes .+ 0.003 .+ 0.002 .* x .- 0.001 .* y) rtol=1e-14 atol=1e-16
        swapped=P.predictions(vcat(p,-q[end-3:end]),3,x,y,axis,cfg,.03,.07)
        @test swapped.fwd==pred.bwd && swapped.bwd==pred.fwd
        z=mean; zimg=reshape(z,length(ys),length(xs))
        d=(;xs,ys,x,y,z,zimg,zfull=z,noise=.01,axisctx=axis)
        initial=G.ChainModelResult(n=3,params=p,success=true,amp_min=.03,amp_range=.07)
        before=copy(p)
        fitted=P.fit(initial,d,cfg,pred.fwd,pred.bwd,"paired",options)
        @test fitted.converged && fitted.valid
        @test fitted.rss<1e-12
        @test fitted.params[end-3:end]≈q[end-3:end] atol=1e-7
        @test fitted.np==length(p)+4 && fitted.nd==2length(z)
        @test fitted.gcv==fitted.nd/(fitted.nd-fitted.np)^2*fitted.rss
        @test fitted.rss≈fitted.rss_fwd+fitted.rss_bwd atol=1e-25
        @test initial.params==before
        @test all(fitted.lower.<=fitted.params.<=fitted.upper)
        @test fitted.iterations<=options.local_maxiter
        fused=P.fit(initial,d,cfg,pred.fwd,pred.bwd,"fused",options)
        @test fused.valid && fused.result.n==initial.n
        @test fused.np==length(p) && fused.nd==length(z)
        @test fused.rss<1e-20
        @test_throws ErrorException P.fit(initial,d,cfg,pred.fwd .+ 0.1,pred.bwd,"paired",options)
        @test_throws ErrorException P.fit(initial,d,cfg,fill(NaN,length(z)),pred.bwd,"paired",options)
        @test_throws ErrorException P.fit(initial,d,cfg,pred.fwd,pred.bwd,"fallback",options)
        # Native mean validity is never rescued by acquisition nuisance terms.
        strict=deepcopy(cfg); strict.residual_peak_snr_threshold=0.
        bad=P.fit(initial,d,strict,pred.fwd,pred.bwd,"paired",options)
        @test !bad.valid && occursin("residual high",bad.reason)
    end
end

@testset "Strict paired settings reject label and count controls" begin
    raw=TOML.parsefile(SETTINGS)
    mktempdir() do dir
        path=joinpath(dir,"bad.toml")
        for (section,key,value) in (("selection","expected_N",6),("selection","count_policy","reselect"),
                ("model","gain_delta_max",1.),("model","local_maxiter",2.5),("model","labels","NKNNKN"))
            bad=deepcopy(raw); bad[section][key]=value
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException P.settings(path)
        end
    end
end
