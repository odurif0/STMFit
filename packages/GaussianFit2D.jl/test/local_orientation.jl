using LinearAlgebra, Statistics
@testset "Deterministic local model axes, unchanged global and circular paths" begin
    G = GaussianFit2D
    for theta in (0.0,0.8,2.7), n in (1,2,3,7), degree in (1,2)
        ctx=(axis=[cos(theta),sin(theta)],)
        ts=collect(range(-1.0,1.0;length=max(n,2)))[1:n]
        us=.12 .+ .18ts .- .07ts.^2
        cfg=G.ChainSweepConfig(chain_peak_orientation="local_tangent",chain_tangent_degree=degree)
        axes=G._chain_peak_axes(ts,us,ctx,cfg)
        for (i,(a,b)) in enumerate(axes)
            slope=n==1 ? 0.0 : degree==1 || n==2 ? .18-.07sum(ts)/length(ts)*2 : .18-.14ts[i]
            @test a ≈ cos(theta+atan(slope)) atol=2e-14
            @test b ≈ sin(theta+atan(slope)) atol=2e-14
            @test hypot(a,b) ≈ 1.0 atol=4e-16
        end
        cfg.chain_circular_sigmas=true
        @test G._chain_peak_axes(ts,us,ctx,cfg)==fill(Tuple(ctx.axis),n)
        cfg.chain_peak_orientation="global"
        @test G._chain_peak_axes(ts,us,ctx,cfg)==fill(Tuple(ctx.axis),n)
    end
    ctx=(origin=(0.,0.),axis=[.8,.6],perp=[-.6,.8],tmin=-1.7,tmax=1.7)
    for circular in (false,true), profile in (:gaussian,:split)
        cfg=G.ChainSweepConfig(chain_circular_sigmas=circular,peak_profile=profile,chain_tilted_baseline=true)
        n=3; p=zeros(G._chain_nparams(n,cfg)); p[1:3].=[.01,.003,-.007]
        start=3+n+G._chain_spacing_param_count(n,cfg)
        p[start+1:start+n].=[.3,-.2,.45]
        if !circular; p[start+n+1:start+2n].=-.8; p[start+2n+1:start+3n].=.7; end
        profile==:split && (p[end-n+1:end].=[-.4,.2,.6])
        x=collect(range(-1.5,1.5;length=43)); y=.2sin.(3x)
        b0,fs,ts,us=G._decode_chain(p,n,ctx,cfg;amp_min=.03,amp_range=.08)
        old=@. b0+p[2]*x+p[3]*y
        for f in fs
            ax,ay=ctx.axis
            if profile==:gaussian
                @. old += f.amplitude*exp(-.5*((((x-f.x_nm)*ax+(y-f.y_nm)*ay)/f.sigma_x_nm)^2+(((x-f.x_nm)*(-ay)+(y-f.y_nm)*ax)/f.sigma_y_nm)^2))
            else
                for i in eachindex(x)
                    dt=(x[i]-f.x_nm)*ax+(y[i]-f.y_nm)*ay
                    du=-(x[i]-f.x_nm)*ay+(y[i]-f.y_nm)*ax
                    st=dt<0 ? f.sigma_x_nm/sqrt(f.skew_ratio) : f.sigma_x_nm*sqrt(f.skew_ratio)
                    old[i]+=f.amplitude*exp(-.5*((dt/st)^2+(du/f.sigma_y_nm)^2))
                end
            end
        end
        global_values=G._chain_model_values(x,y,p,n,ctx,cfg;amp_min=.03,amp_range=.08)
        @test global_values==old
        np=G._chain_nparams(n,cfg)
        cfg.chain_peak_orientation="local_tangent"
        local_values=G._chain_model_values(x,y,p,n,ctx,cfg;amp_min=.03,amp_range=.08)
        @test G._chain_nparams(n,cfg)==np
        if circular && profile==:gaussian
            @test local_values==global_values
        else
            @test maximum(abs,local_values-global_values)>1e-5
            # Independent unscaled Vandermonde fit and rotated coordinates.
            c=hcat(ones(n),ts,ts.^2)\us
            expected=@. b0+p[2]*x+p[3]*y
            for (j,f) in enumerate(fs)
                angle=atan(c[2]+2c[3]*ts[j]); ax=.8cos(angle)-.6sin(angle); ay=.6cos(angle)+.8sin(angle)
                for i in eachindex(x)
                    dt=(x[i]-f.x_nm)*ax+(y[i]-f.y_nm)*ay
                    du=-(x[i]-f.x_nm)*ay+(y[i]-f.y_nm)*ax
                    st=dt<0 ? f.sigma_x_nm/sqrt(f.skew_ratio) : f.sigma_x_nm*sqrt(f.skew_ratio)
                    expected[i]+=f.amplitude*exp(-.5*((dt/st)^2+(du/f.sigma_y_nm)^2))
                end
            end
            @test local_values≈expected atol=2e-16
            changed=copy(p); changed[start+1]+=.1
            _,_,tt,uu=G._decode_chain(changed,n,ctx,cfg;amp_min=.03,amp_range=.08)
            @test G._chain_peak_axes(ts,us,ctx,cfg)!=G._chain_peak_axes(tt,uu,ctx,cfg)
        end
    end
    cfg=G.ChainSweepConfig(chain_peak_orientation="wrong")
    @test_throws ErrorException G._chain_peak_axes([0.],[0.],ctx,cfg)
    cfg.chain_peak_orientation="local_tangent"; cfg.chain_tangent_degree=3
    @test_throws ErrorException G._chain_peak_axes([0.],[0.],ctx,cfg)
    cfg.chain_tangent_degree=2
    @test_throws ErrorException G._chain_peak_axes([0.,0.],[0.,.1],ctx,cfg)
end
