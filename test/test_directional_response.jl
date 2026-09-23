using Test, TOML, Random
include(joinpath(@__DIR__,"lib","directional_response.jl"))
const D=DirectionalResponse
const SETTINGS=joinpath(@__DIR__,"..","config","directional_response.toml")

"Independent analytic steady-state acquisition, no reuse of transfer kernel."
function analytic_pair(ny,nx,mu;dx=0.,dy=0.,gain=1.)
    a=mu/(1+mu); f=zeros(ny,nx); b=similar(f)
    waves=[(.21,.09,1.0),(.73,-.13,.7),(1.17,.23,.4)]
    for x in 1:nx,y in 1:ny
        f[y,x]=sum(v*real((1-a)/(1-a*exp(-im*w))*exp(im*w*mu)*exp(im*(w*x+k*y))) for (w,k,v) in waves)
        b[y,x]=gain*sum(v*real((1-a)/(1-a*exp(im*w))*exp(-im*w*mu)*exp(im*(w*(x-dx)+k*(y-dy)))) for (w,k,v) in waves)
    end
    return f,b
end

@testset "Response settings and acquisition-only panel" begin
    s=D.settings(SETTINGS)
    files=["scan$(lpad(i,3,'0')).sxm" for i in 1:146]
    @test getproperty.(D.panel(reverse(files),s),:index)==[1,37,74,110,146]
    @test_throws ErrorException D.panel(vcat(files,first(files)),s)
    @test_throws ErrorException D.panel(["../scan.sxm" for _ in 1:5],s)
    mktempdir() do dir
        for mutate in (t->t["model"]["mean_lag_px"]= [0.,.2,1.],
                       t->t["preprocessing"]["holdout_buffer_rows"]=0,
                       t->t["model"]["gain_min"]=-1,
                       t->t["selection"]["expected_N"]=6,
                       t->t["model"]["kernel_tail_l1"]=1e-2)
            t=TOML.parsefile(SETTINGS); mutate(t); p=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,t),p,"w")
            @test_throws ErrorException D.settings(p)
        end
    end
end

@testset "Relative kernel is not a smoothing shortcut" begin
    s=D.settings(SETTINGS)
    for mu in s.mean_lag_px,sign in (-1,1)
        k=D.kernel(mu,sign,s.kernel_tail_l1)
        @test k.tail<=s.kernel_tail_l1
        @test isapprox(sum(k.weights),1;atol=2s.kernel_tail_l1)
        @test isapprox(sum(k.offsets.*k.weights),0;atol=1e-9)
        @test isapprox(sum(abs2,k.weights),1;atol=2s.kernel_tail_l1)
        for w in range(-pi,pi;length=33)
            v=sum(k.weights.*exp.(im*w.*k.offsets))
            @test isapprox(abs(v),1;atol=2s.kernel_tail_l1)
        end
    end
    f,b=analytic_pair(20,420,2.)
    p=D.transfer(f,D.kernel(2.,1,s.kernel_tail_l1))
    q=D.transfer(b,D.kernel(2.,-1,s.kernel_tail_l1))
    ok=isfinite.(p).&isfinite.(q)
    @test maximum(abs.(p[ok].-b[ok]))<1e-10
    @test maximum(abs.(q[ok].-f[ok]))<1e-10
    @test D.transfer(f,D.kernel(0.,1,s.kernel_tail_l1))==f
    f[10,200]=NaN
    p=D.transfer(f,D.kernel(1.,1,s.kernel_tail_l1))
    for off in D.kernel(1.,1,s.kernel_tail_l1).offsets
        x=200-off
        1<=x<=size(p,2) && @test isnan(p[10,x])
    end
    @test D.sample([1. 3.;5. 7.],1.5,1.5)==4.
    @test D.sample([2. NaN;NaN NaN],1.,1.)==2.
end

@testset "Observed common support and source-safe folds" begin
    s=merge(D.settings(SETTINGS),(;min_rows=4,min_pixels=64,min_pixels_per_row=4,holdout_band_rows=32))
    a=ones(128,360); b=copy(a); a[23,180]=NaN; b[65,182]=NaN
    r=D.common_mask(a,b,s)
    @test count(r.mask)>0
    for y in axes(a,1),x in axes(a,2)
        r.mask[y,x] || continue
        @test all(isfinite,a[y-r.ry:y+r.ry,x-r.rx:x+r.rx])
        @test all(isfinite,b[y-r.ry:y+r.ry,x-r.rx:x+r.rx])
        for yy in (y-r.ry,y+r.ry)
            @test 1+mod(fld(yy-1,s.holdout_band_rows),2)==r.fold[y]
        end
    end
end

@testset "Recovery, reverse control, absence and heldout training isolation" begin
    s=merge(D.settings(SETTINGS),(;mean_lag_px=[0.,1.,2.,3.],residual_shift_max_px=2.,residual_shift_step_px=1.,
        min_rows=4,min_pixels=64,min_pixels_per_row=4,holdout_band_rows=24,holdout_buffer_rows=3))
    f,b=analytic_pair(96,360,2.;dx=1.,dy=-1.,gain=1.15)
    r=D.scan(f,b,s)
    @test r.summary.supported
    for p in filter(x->x.mode=="physical",r.chosen)
        @test (p.lag,p.dx,p.dy)==(2.,1.,-1.)
        @test p.test_rss<1e-9
        @test p.gain_bwd≈1.15
        @test p.gain_fwd≈1/1.15
    end
    zero=D.scan(f,copy(f),s)
    @test !zero.summary.supported
    @test all(p->p.lag==0 && p.dx==0 && p.dy==0,zero.chosen)
    # Reversing physical acquisition directions is detected, not accepted.
    f,b=analytic_pair(96,360,2.)
    reverse=D.scan(b,f,s)
    @test !reverse.summary.supported
    @test all(p->p.lag==2 && p.test_rss<1e-9,filter(p->p.mode=="reversed",reverse.chosen))
    # Only the training moments determine coefficients and candidate ranking.
    k=D.kernel(2.,1,s.kernel_tail_l1); kinv=D.kernel(2.,-1,s.kernel_tail_l1)
    support=D.common_mask(f,b,s)
    m=D.moments(f,b,D.transfer(b,kinv),D.transfer(f,k),0.,0.,support)
    v=D.scores(m,1,s); m.sums[:,2].*=100
    w=D.scores(m,1,s)
    @test (v.gain_fwd,v.gain_bwd,v.train_rss)==(w.gain_fwd,w.gain_bwd,w.train_rss)
    @test_throws ErrorException D.scan(ones(96,360),fill(NaN,96,360),s)
    @test_throws ErrorException D.scan(ones(96,360),ones(96,360),s)
    f2,b2=copy(f),copy(b)
    for y in axes(f,1)
        mod(fld(y-1,s.holdout_band_rows),2)==1 || continue
        f2[y,:].*=1.3; b2[y,:].*=0.7
    end
    original=D.scan(f,b,s); perturbed=D.scan(f2,b2,s)
    tr1=filter(p->p.fold==1,original.records); tr2=filter(p->p.fold==1,perturbed.records)
    @test getproperty.(tr1,:train_rss)==getproperty.(tr2,:train_rss)
    @test getproperty.(tr1,:test_rss)!=getproperty.(tr2,:test_rss)
end

include(joinpath(@__DIR__,"verify_directional_response.jl"))
const V=VerifyDirectionalResponse
const CLI=V.D
@testset "Mocked complete panel, saved-only independent verification and CLI isolation" begin
    mktempdir() do dir
        t=TOML.parsefile(SETTINGS)
        t["model"]["mean_lag_px"]=[0.,1.,2.,3.]
        t["model"]["residual_shift_max_px"]=2.
        t["model"]["residual_shift_step_px"]=1.
        t["selection"]["panel_quantiles"]=[0.,0.5,1.]
        t["selection"]["min_supporting_scans"]=2
        t["selection"]["min_rows"]=4; t["selection"]["min_pixels"]=64
        t["selection"]["min_pixels_per_row"]=4
        t["preprocessing"]["holdout_band_rows"]=24; t["preprocessing"]["holdout_buffer_rows"]=3
        cfg=joinpath(dir,"settings.toml"); open(io->TOML.print(io,t),cfg,"w")
        files=["scan$(i).sxm" for i in 1:3]
        for file in files; open(io->write(io,"synthetic fixture, not SXM"),joinpath(dir,file),"w"); end
        shifts=joinpath(dir,"shifts.tsv")
        # Metadata format is intentionally the production's exact two-column order.
        physical=joinpath(@__DIR__,"..","config","chitosan.toml")
        args=["--data-dir",dir,"--base-shifts",shifts,"--config",physical,"--settings",cfg,"--outdir",joinpath(dir,"out")]
        open(shifts,"w") do io
            println(io,"file\tbwd_sample_dx_px")
            for file in files; println(io,file,"\t0"); end
        end
        opts=CLI.parse_cli(vcat(args,["--dry-run"]))
        CLI.execute(opts;loader=(_...)->error("Dry-run loaded images"))
        @test !ispath(opts["--outdir"])
        @test_throws ErrorException CLI.parse_cli(vcat(args,["--expected-N","6"]))
        f,b=analytic_pair(96,360,2.;dx=1.,dy=-1.)
        loader=(_,pre)->(;a=f,b,xs=collect(1.:360.),ys=collect(1.:96.),step=(1.,1.))
        header=_ -> Dict("SCAN_TIME"=>"0.05 0.05","SCAN_ANGLE"=>"0.0","SCAN_DIR"=>"down","BIAS"=>"-0.3")
        actual=CLI.parse_cli(args); CLI.execute(actual;loader,header_reader=header)
        @test only(CLI.table(joinpath(actual["--outdir"],"panel_result.tsv")))["directional_hypothesis_supported"]=="true"
        V.verify(actual["--outdir"])
        @test_throws ErrorException CLI.parse_cli(args)
        # Unequal line times are reported, never replaced by a same-speed model.
        args[end]=joinpath(dir,"unequal_time")
        unequal=_ -> Dict("SCAN_TIME"=>"0.05 0.10")
        actual=CLI.parse_cli(args); CLI.execute(actual;loader=(_...)->error("Unequal-time images should not load"),header_reader=unequal)
        @test all(r->r["status"]=="failed",CLI.table(joinpath(actual["--outdir"],"summary.tsv")))
        V.verify(actual["--outdir"])
    end
end
