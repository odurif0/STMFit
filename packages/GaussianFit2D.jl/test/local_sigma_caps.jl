using Random, LinearAlgebra

@testset "opt-in Gaussian width-conditioned gaps" begin
    G = GaussianFit2D
    rng = MersenneTwister(824)
    base = G.ChainSweepConfig(sigma_parallel_min_nm=0.191, sigma_parallel_max_nm=0.509,
        sigma_perp_min_nm=0.191, sigma_perp_max_nm=0.509,
        spacing_min_nm=0.35, spacing_max_nm=0.75, max_overlap=0.6)
    @test base.overlap_constraint == "global_sigma_max"
    @test G._effective_spacing_min_nm(base) == sqrt(-2log(0.6))*0.509
    localcfg = deepcopy(base); localcfg.overlap_constraint = "local_sigma_cap"
    @test G._effective_spacing_min_nm(localcfg) == 0.35
    axis = (origin=(0.,0.), axis=[1.,0.], perp=[0.,1.], tmin=0., tmax=1.3)
    @test !G._chain_can_fit_support(4, axis, base)
    @test G._chain_can_fit_support(4, axis, localcfg)
    for fieldvalue in ((:peak_profile,:split), (:overlap_constraint,"typo"),
                       (:max_overlap,0.), (:max_overlap,1.), (:spacing_max_nm,0.1),
                       (:sigma_parallel_min_nm,-0.1))
        bad=deepcopy(localcfg); setproperty!(bad,fieldvalue...)
        @test_throws ErrorException G._effective_spacing_min_nm(bad)
    end
    for n in (1,2,3,5,9), circular in (true,false), shared in (0,1,2), spacing in ("free","uniform","alternating")
        cfg=deepcopy(localcfg); cfg.chain_circular_sigmas=circular
        cfg.shared_sigma_types=shared; cfg.chain_spacing_model=spacing
        theta=0.71; ax=[cos(theta),sin(theta)]; perp=[-sin(theta),cos(theta)]
        ac=(origin=(0.1,-0.2),axis=ax,perp=perp,tmin=-0.25,tmax=max(n-1,1)*0.51)
        for trial in 1:5
            p=trial<=3 ? fill((-60.,0.,60.)[trial],G._chain_nparams(n,cfg)) :
                10randn(rng,G._chain_nparams(n,cfg))
            _,feats,ts,us,sp,sq=G._decode_chain(p,n,ac,cfg)
            @test all(0.35-1e-12 .<= diff(ts) .<= 0.75+1e-12)
            @test first(ts)>=ac.tmin-1e-12 && last(ts)<=ac.tmax+1e-12
            @test all(0.191 .<= sp .<= 0.509) && all(0.191 .<= sq .<= 0.509)
            independent=0.0
            for i in 1:n-1, j in i+1:n
                distance=hypot(ts[i]-ts[j],us[i]-us[j])
                value=exp(-distance^2/(2maximum((sp[i],sp[j],sq[i],sq[j]))^2))
                independent=max(independent,value)
                @test value<=0.6+1e-12
                # Actual directional Gaussian at the other center is no larger
                # than this radial envelope, in either direction.
                for k in (i,j)
                    g=exp(-0.5*((ts[i]-ts[j])^2/sp[k]^2+(us[i]-us[j])^2/sq[k]^2))
                    @test g<=value+1e-12
                end
            end
            @test G._chain_pair_overlap(feats,sp,sq)≈independent atol=1e-12
            caps=G._chain_local_sigma_caps(n,diff(ts),cfg)
            for i in 1:n
                k=G._chain_lobe_type(i,length(caps))
                @test max(sp[i],sq[i])<=caps[k]+1e-12
            end
            r=G.ChainModelResult(n=n,success=true,params=p)
            G._chain_metrics!(r,ac,cfg)
            @test r.overlap≈independent atol=1e-12
        end
    end
    # Minimum-width floor can dominate the literal spacing floor. The historical
    # max-spacing clipping must not hide an impossible local-width domain.
    c=deepcopy(localcfg); c.sigma_parallel_min_nm=c.sigma_perp_min_nm=0.45
    @test G._effective_spacing_min_nm(c)≈sqrt(-2log(0.6))*0.45
    c.spacing_max_nm=0.4
    @test_throws ErrorException G._effective_spacing_min_nm(c)
    # Initialization and circular nesting retain bounds and the intended shape.
    xs=collect(range(-1.,2.5;length=33)); ys=collect(range(-0.5,0.5;length=15))
    z=[0.1+exp(-(x-0.2)^2/0.08-y^2/0.08) for y in ys,x in xs]
    ac=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-0.2,tmax=1.3)
    c=deepcopy(localcfg); c.chain_circular_sigmas=true
    p,amin,arange=G._pack_chain_initial(xs,ys,z,4,ac,c)
    _,_,ts,_,sp,sq=G._decode_chain(p,4,ac,c;amp_min=amin,amp_range=arange)
    @test G._chain_can_fit_support(4,ac,c) && all(isfinite,p)
    @test maximum(sp)<=minimum(diff(ts))/sqrt(-2log(0.6))+1e-12 ||
        all(sp[i]<=min(i>1 ? ts[i]-ts[i-1] : Inf,i<4 ? ts[i+1]-ts[i] : Inf)/sqrt(-2log(0.6))+1e-12 for i in 1:4)
    e=deepcopy(c); e.chain_circular_sigmas=false
    prefix=1+4+G._chain_spacing_param_count(4,c)+4
    pe=vcat(p[1:prefix],p[prefix+1:end],p[prefix+1:end])
    xx=repeat(xs;inner=length(ys)); yy=repeat(ys;outer=length(xs))
    @test G._chain_model_values(xx,yy,p,4,ac,c;amp_min=amin,amp_range=arange)==
          G._chain_model_values(xx,yy,pe,4,ac,e;amp_min=amin,amp_range=arange)
    # Legacy mode's widths are still independent of the neighboring gaps.
    for circular in (true,false)
        c=deepcopy(base); c.chain_circular_sigmas=circular
        p=randn(rng,G._chain_nparams(2,c))
        _,feats,ts,us,sp,sq=G._decode_chain(p,2,ac,c)
        prefix=1+2+G._chain_spacing_param_count(2,c)+2
        @test sp==[0.191+(0.509-0.191)*G._rsigmoid(q) for q in p[prefix+1:prefix+2]]
        r=G.ChainModelResult(n=2,success=true,params=p); G._chain_metrics!(r,ac,c)
        @test r.overlap==G._chain_overlap(feats,sum(sp)/2,sum(sq)/2)
    end
end
