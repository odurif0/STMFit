using Test, Random, Statistics, LinearAlgebra, TOML
include(joinpath(@__DIR__,"diagnose_vacuum_surface_transfer.jl"))
const R=DiagnoseVacuumSurfaceTransfer
const F=R.F
const SETTINGS=joinpath(@__DIR__,"../config/vacuum_surface_transfer.toml")
const C=F.settings(SETTINGS)

@testset "Exact common-offset finite dictionary" begin
    rng=MersenneTwister(614)
    for trial in 1:100
        costs=[(;n=rand(rng,3:20),mu=randn(rng,4),v=rand(rng,4)) for _ in 1:3]
        got=F.common_offset(costs)
        expected=Inf
        for choices in Iterators.product((1:4 for _ in costs)...)
            b=sum(c.n*c.mu[j] for (c,j) in zip(costs,choices))/sum(c.n for c in costs)
            loss=sum(c.v[j]+c.n*(b-c.mu[j])^2 for (c,j) in zip(costs,choices))
            expected=min(expected,loss)
        end
        @test isapprox(got.loss,expected;rtol=1e-12,atol=1e-12)
        @test all(got.choices[j]==argmin(c.v.+c.n.*(got.offset.-c.mu).^2) for (j,c) in enumerate(costs))
        shift=11.3
        shifted=F.common_offset([(;c.n,mu=c.mu.+shift,c.v) for c in costs])
        @test isapprox(shifted.offset,got.offset+shift;atol=1e-12)
        @test isapprox(shifted.loss,got.loss;atol=1e-11)
    end
    same=F.common_offset([(;n=5,mu=[1.,1.,2.],v=[2.,0.,3.])])
    @test same.choices==[2] && same.offset==1 && same.loss==0
    single=F.common_offset([(;n=5,mu=[1.],v=[0.]),(;n=10,mu=[2.],v=[1.])])
    @test single.offset≈5/3
    @test single.segments==1
    @test_throws ErrorException F.common_offset(NamedTuple[])
    @test_throws ErrorException F.envelope([NaN],[1.],1)
end

@testset "Interpolation, physical height contrast, dictionary recovery" begin
    grid=collect(-.32:.04:.32)
    plane=[1+2t-3u for u in grid,t in grid]
    for t in -.25:.025:.25,u in -.25:.025:.25
        @test isapprox(F.surface_value(plane,t,u,.32,.04),1+2t-3u;atol=1e-14)
    end
    @test_throws ErrorException F.surface_value(plane,.4,0.,.32,.04)
    g0=[1+.025exp(-((t-.045)/.09)^2-(u/.12)^2)+.007sin(13t+9u) for u in grid,t in grid]
    g1=[1.03+.020exp(-(t/.14)^2-((u+.07)/.06)^2)-.008sin(13t-9u) for u in grid,t in grid]
    maps=Dict(0=>g0,1=>g1,-1=>(g0+g1)/2)
    coords=[(x,y) for x in -.15:.03:.15 for y in -.15:.03:.15 if hypot(x,y)<=.16]
    p=(;dx=first.(coords),dy=last.(coords),values=zeros(length(coords)))
    ss=F.states(C,"chemical")
    examples=[ss[7],ss[end-2],ss[133]]; offset=-.94
    pp=[merge(p,(;values=offset.+F.prediction(p,maps[s.type],s,C))) for s in examples]
    r=F.fit_dictionary(pp,maps,C,"chemical")
    @test r.fit.loss<1e-24
    @test r.fit.offset≈offset
    @test [r.ss[k].type for k in r.fit.choices]==getproperty.(examples,:type)
    @test all(isapprox(r.pred[j],p.values;atol=1e-12) for (j,p) in enumerate(pp))
    common=F.fit_dictionary(pp,maps,C,"common")
    @test common.fit.loss>r.fit.loss+1e-5
    @test mean(g1)-mean(g0)>0.02
    @test all(s.type==-1 for s in F.states(C,"common"))
    @test length(ss)==432
    @test all(length(c.mu)==432 && length(c.v)==432 for c in r.costs)
end

@testset "Calibration split is target-holdout invariant" begin
    c=deepcopy(C); p=c["preprocessing"]; q=c["selection"]
    p["calibration_block_rows"]=16; p["max_lag_nm"]=.15
    q["peak_neighborhood_nm"]=.04; q["min_band_rows"]=2
    q["min_band_pixels"]=30; q["min_row_pixels"]=8
    xs=collect(range(0,2.;length=96)); ys=collect(range(0,2.;length=128))
    rng=MersenneTwister(99); a=randn(rng,length(ys),length(xs)); b=fill(NaN,size(a)); lag=3
    co=[.17,.03,-.02]
    for x in 1:length(xs)-lag,y in eachindex(ys)
        b[y,x+lag]=a[y,x]+co[1]+co[2]*xs[x]+co[3]*ys[y]
    end
    got=F.calibrate_pair(a,b,xs,ys,c)
    @test got.reg.accepted
    @test got.reg.applied_dx_px==lag
    @test got.coeff≈co
    @test !any(got.masks.cal.&got.masks.test)
    @test count(got.masks.cal)==count(got.masks.test)>0
    altered=copy(b)
    altered[.!got.masks.cal,:].=9999randn(rng,count(.!got.masks.cal),length(xs))
    again=F.calibrate_pair(a,altered,xs,ys,c)
    @test isequal(got.reg,again.reg)
    @test got.coeff==again.coeff
    @test all(!got.reg.support[y,x] for y in findall(got.masks.test),x in eachindex(xs))
    @test_throws ErrorException F.calibrate_pair(a,fill(NaN,size(a)),xs,ys,c)
end

function toy_image(;missing_target=false)
    xs=collect(range(0,3.;length=96)); ys=collect(range(0,3.;length=96))
    z=[10+.004x-.002y+.1exp(-((x-.9)/.16)^2-((y-.9)/.16)^2)+
        .08exp(-((x-2.1)/.16)^2-((y-2.1)/.16)^2) for y in ys,x in xs]
    f=R.STMSXMIO.SXMChannel("Z","m","fwd",z*1e-9)
    b=R.STMSXMIO.SXMChannel("Z","m","bwd",(z.+.03)*1e-9)
    R.STMSXMIO.SXMImage("synthetic.sxm",Dict("BIAS"=>"-0.300"),96,96,(3.,3.),(0.,0.),missing_target ? [f] : [f,b])
end

@testset "Source-only extraction, observations and complete synthetic output" begin
    xx=collect(0:.1:3); zz=zeros(length(xx),length(xx)); zz[8:9,8:9].=1; zz[22,22]=1.2
    peaks=F.anchors(xx,xx,zz,zz,C)
    @test length(peaks)==2
    @test sort(getproperty.(peaks,:plateau_size))==[1,4]
    @test isempty(F.anchors(xx,xx,zero(zz),zero(zz),C))
    img=toy_image(); s=R.source_data(img,"fwd",C)
    @test length(s.aa)>=2
    @test length(s.pp)>=2
    @test length(vcat(getproperty.(s.pp,:inds)...))==length(unique(vcat(getproperty.(s.pp,:inds)...)))
    @test all(hypot(a.x-b.x,a.y-b.y)>2C["selection"]["patch_radius_nm"] for (i,a) in enumerate(s.aa) for b in s.aa[i+1:end])
    obs=R.source_data(toy_image(missing_target=true),"fwd",C)
    @test isequal(s,obs)
    grid=collect(-.32:.04:.32)
    a=[1+.03exp(-(t/.12)^2-(u/.12)^2) for u in grid,t in grid]
    b=[1.025+.02exp(-(t/.08)^2-((u+.05)/.1)^2) for u in grid,t in grid]
    maps=(heights=Dict(j=>Dict(-1=>(a+b)/2,0=>a,1=>b) for j in 1:3),isos=Dict(j=>Float64(j) for j in 1:3))
    mktempdir() do dir
        out=joinpath(dir,"test"); result=R.one_direction(img,"fwd",maps,C,out)
        @test result.patches==length(s.pp)
        models=R.V.tsv(joinpath(out,"models.tsv")); patches=R.V.tsv(joinpath(out,"patches.tsv"))
        observations=R.V.tsv(joinpath(out,"observations.tsv"))
        @test length(models)==6
        @test length(patches)==6length(s.pp)
        @test Set(r["arm"] for r in models)==Set(["common","chemical"])
        for model in models
            @test R.num(model,"training_sse_nm2")≈R.num(model,"optimizer_sse_nm2")
            @test R.integer(model,"test_pixels")>0
            @test R.num(model,"source_copy_sse_nm2")<1e-20
        end
        masks=F.row_masks(96,C)
        @test all(masks.test[R.integer(r,"row")] for r in observations if isfinite(R.num(r,"target_nm")))
        changed=copy(img.channels[2].data)
        changed[.!masks.cal,:].+=123e-9
        bwd=R.STMSXMIO.SXMChannel("Z","m","bwd",changed)
        altered=R.STMSXMIO.SXMImage(img.filepath,img.header,img.width,img.height,img.range_nm,img.offset_nm,[img.channels[1],bwd])
        other=joinpath(dir,"heldout_changed"); R.one_direction(altered,"fwd",maps,C,other)
        for file in ("calibration.toml","anchors.tsv","registration.tsv","registration_scores.tsv",
                     ("costs_$(j)_$(arm).tsv" for j in 1:3 for arm in ("common","chemical"))...)
            @test read(joinpath(out,file))==read(joinpath(other,file))
        end
        again=R.V.tsv(joinpath(other,"models.tsv"))
        for (a,b) in zip(models,again), key in ("offset_nm","training_sse_nm2","optimizer_sse_nm2","flat_height_nm")
            @test a[key]==b[key]
        end
        @test all(R.num(b,"test_sse_nm2")>R.num(a,"test_sse_nm2") for (a,b) in zip(models,again))
        @test_throws ErrorException R.one_direction(img,"fwd",maps,C,out)
        absent=joinpath(dir,"absent"); result=R.one_direction(toy_image(missing_target=true),"fwd",maps,C,absent)
        @test occursin("Missing target",result.status)
        @test all(R.integer(r,"test_pixels")==0 for r in R.V.tsv(joinpath(absent,"models.tsv")))
    end
    @test R.matching_bias(Dict("BIAS"=>"-0.3"),C)
    @test !R.matching_bias(Dict("BIAS"=>"0.3"),C)
    @test !R.matching_bias(Dict{String,String}(),C)
    @test_throws ErrorException R.main(["--truth","forbidden"])
    mktempdir() do dir
        for mutation in (c->(c["model"]["truth"]=0), c->(c["model"]["height_gain"]=2),
                         c->(c["selection"]["patch_radius_nm"]=.31),
                         c->(c["preprocessing"]["holdout_buffer_rows"]=32))
            changed=deepcopy(C); mutation(changed); path=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,changed),path,"w")
            @test_throws ErrorException F.settings(path)
        end
    end
end
