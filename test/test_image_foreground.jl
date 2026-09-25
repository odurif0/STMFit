# Julia 1.13 synthetic-only tests; optional argument is a NEW synthetic-metrics TSV.
using Test, Statistics, Random, LinearAlgebra, TOML
include(joinpath(@__DIR__,"diagnose_image_foreground.jl"))
const D=DiagnoseImageForeground
const G=D.F
const CFGPATH=joinpath(@__DIR__,"../config/image_foreground.toml")
const C=G.settings(CFGPATH)

function image_fixture(z;target=z.+0.03,span=3.)
    ny,nx=size(z)
    channels=[D.STMSXMIO.SXMChannel("Z","m","fwd",z.*1e-9),
        D.STMSXMIO.SXMChannel("Z","m","bwd",target.*1e-9)]
    D.STMSXMIO.SXMImage("synthetic.sxm",Dict("BIAS"=>"-0.300"),nx,ny,(span,span),(0.,0.),channels)
end

@testset "Exact Otsu criterion, constant and finite-only inputs" begin
    rng=MersenneTwister(71)
    for trial in 1:100
        values=randn(rng,rand(rng,3:70)); o=G.otsu(vcat(values,NaN,Inf))
        v=sort(values); scores=[(j/length(v))*((length(v)-j)/length(v))*(mean(v[1:j])-mean(v[j+1:end]))^2 for j in 1:length(v)-1]
        @test isapprox(o.score,maximum(scores);rtol=1e-12,atol=1e-14)
        @test count(<(o.threshold),values)==o.split
        @test o.n==length(values)
        moved=G.otsu(3.7values.+8.2)
        @test isapprox(moved.threshold,3.7o.threshold+8.2;atol=1e-12)
        @test isapprox(moved.score,3.7^2*o.score;rtol=1e-12)
    end
    o=G.otsu([-1.,-1.,1.,1.]); @test o.threshold==0 && o.split==2 && o.score==1
    @test G.otsu(ones(6)).split==0
    @test G.otsu([1.]).threshold==1
    @test isnan(G.otsu([NaN]).threshold)
    @test !any(G.segment(ones(4,6),0.,C).mask)
    @test G.segment(fill(NaN,4,6),0.,C).status=="unavailable"
end

@testset "All components and anisotropic box guard are retained" begin
    rng=MersenneTwister(14)
    for rx in (0,1,3,20),ry in (0,1,2,20)
        mask=rand(rng,9,11).<0.04
        expected=[any(mask[max(1,y-ry):min(end,y+ry),max(1,x-rx):min(end,x+rx)]) for y in 1:9,x in 1:11]
        @test G.box_guard(mask,rx,ry)==expected
    end
    m=falses(7,9); m[1,1]=m[2,2]=true; m[5:6,7:8].=true
    z=zeros(7,9); z[6,6]=NaN
    r=G.components(m,z,collect(1.:9.),collect(1.:7.))
    @test length(r.rows)==2 && sum(q.pixels for q in r.rows)==count(m)
    @test r.labels[1,1]==r.labels[2,2]==1
    @test r.rows[1].touches_edge && r.rows[2].touches_missing
    @test !r.rows[1].touches_missing && !r.rows[2].touches_edge
end

@testset "Joint background solve under unequal row coverage" begin
    xs=collect(range(0,3.;length=60)); ys=collect(range(0,2.;length=40))
    xn=(xs.-1.5)./3; levels=[.013y+.04sin(y) for y in ys]
    signal=zeros(40,60); signal[14:26,26:35].=0.17
    z=[levels[y]+.31xn[x]+signal[y,x] for y in eachindex(ys),x in eachindex(xs)]
    for y in 1:40; z[y,1:mod(y,18)].=NaN; end
    guard=falses(size(z)); guard[12:28,24:37].=true
    r=G.joint_background(z,xs,guard,C)
    @test r.status=="ok" && r.slope≈.31
    @test all(isfinite(r.corrected[i])==isfinite(z[i]) for i in eachindex(z))
    @test maximum(abs.(r.corrected[isfinite.(z)]-signal[isfinite.(z)]))<1e-12
    inds=findall(r.used); X=zeros(length(inds),41); response=[z[i] for i in inds]
    for (k,i) in enumerate(inds); X[k,1]=xn[i[2]]; X[k,1+i[1]]=1.; end
    beta=X\response
    @test beta[1]≈r.slope
    @test beta[2:end]≈getproperty.(r.rows,:row_level_centered_nm)
    @test abs(sum(r.corrected[r.used]))<1e-10
    shifted=G.joint_background(z.+0.7,xs,guard,C)
    @test shifted.slope≈r.slope
    @test maximum(abs.(shifted.corrected[isfinite.(z)]-r.corrected[isfinite.(z)]))<1e-12
    blocked=copy(guard); blocked[20,:].=true
    unavailable=G.joint_background(z,xs,blocked,C)
    @test unavailable.status=="ok" && !unavailable.rows[20].usable
    @test all(isnan,unavailable.corrected[20,:])
    no_background=G.joint_background(z,xs,trues(size(z)),C)
    @test no_background.status=="insufficient_background" && all(isnan,no_background.corrected)
end

function synthetic_cases()
    xs=collect(range(0,3.;length=96)); ys=copy(xs)
    baseline=[-64+.004x-.007y+.003sin(8y) for y in ys,x in xs]
    bump=[.1exp(-((x-.9)/.16)^2-((y-.9)/.16)^2)+.08exp(-((x-2.1)/.16)^2-((y-2.1)/.16)^2) for y in ys,x in xs]
    # Gaussian support has no sharp boundary; this contour is synthetic reporting only.
    truth=bump.>0.01
    rng=MersenneTwister(18); noise=.001randn(rng,size(baseline))
    partial=baseline+bump+noise; partial[1:8,:].=NaN; partial[44:49,10:27].=NaN
    row_signal=zeros(size(baseline)); row_signal[30:45,:].=0.04
    curved=[-64+.03cos(2pi*x/3)+.001sin(4y) for y in ys,x in xs]
    [(;name="two_bumps",image=image_fixture(baseline+bump+noise),truth),
     (;name="missing_observations",image=image_fixture(partial),truth),
     (;name="noisy_plane_rows",image=image_fixture(baseline+noise),truth=falses(size(baseline))),
     (;name="noiseless_plane_rows",image=image_fixture(baseline),truth=falses(size(baseline))),
     (;name="full_row_signal_unidentifiable",image=image_fixture(baseline+row_signal),truth=row_signal.>0),
     (;name="curved_background_false_foreground",image=image_fixture(curved),truth=falses(size(baseline)))]
end

const SYNTHETIC_METRICS=NamedTuple[]
@testset "Source-only masks, missingness and retained scientific failure cases" begin
    cases=synthetic_cases()
    for case in cases
        source=D.direction_data(case.image,"fwd",C); obs,r=source.obs,source.r
        @test r.background.status=="ok"
        @test all(!r.initial.mask[i] for i in eachindex(obs.raw) if !isfinite(obs.raw[i]))
        @test all(!r.final.mask[i] for i in eachindex(obs.raw) if !isfinite(obs.raw[i]))
        for (arm,mask,z) in (("initial",r.initial.mask,obs.z_smooth),("joint_background",r.final.mask,r.final_smooth))
            available=isfinite.(z); tp=count(mask.&case.truth); fp=count(mask.&.!case.truth); fn=count(available.&case.truth.&.!mask)
            push!(SYNTHETIC_METRICS,(;case=case.name,arm,raw_pixels=count(isfinite,obs.raw),available_pixels=count(available),
                foreground_pixels=count(mask),true_positive=tp,false_positive=fp,false_negative=fn,
                iou=tp+fp+fn>0 ? tp/(tp+fp+fn) : NaN))
        end
        if case.name=="two_bumps"
            @test length(r.initial_components.rows)==2 && length(r.final_components.rows)==2
            @test all(count(seg.mask.&case.truth)>0 for seg in (r.initial,r.final))
        elseif case.name=="missing_observations"
            @test all(isnan,r.final_smooth[1:8,:])
            @test all(!row.usable for row in r.background.rows[1:8])
        elseif case.name in ("noisy_plane_rows","noiseless_plane_rows","full_row_signal_unidentifiable")
            @test !any(r.initial.mask) && !any(r.final.mask)
        elseif case.name=="curved_background_false_foreground"
            # Deliberately retained model failure: positive curvature is not a molecule.
            @test any(r.initial.mask) && any(r.final.mask)
        end
    end
    source=D.direction_data(first(cases).image,"fwd",C)
    original=first(cases).image
    # Preserve the source bytes; an nm -> m round trip is itself a source mutation.
    channels=[original.channels[1],D.STMSXMIO.SXMChannel("Z","m","bwd",fill(99e-9,size(original.channels[1].data)))]
    changed=D.STMSXMIO.SXMImage(original.filepath,original.header,original.width,original.height,
        original.range_nm,original.offset_nm,channels)
    other=D.direction_data(changed,"fwd",C)
    @test isequal(source.r,other.r)
    mktempdir() do dir
        D.save_direction(dir,source,C)
        for (name,a,T) in (("raw_centered.f64",source.r.centered,Float64),("initial_core.u8",UInt8.(source.r.initial.mask),UInt8))
            @test isequal(D.get_array(dir,name,T,size(a)...),a)
        end
        @test length(D.t(joinpath(dir,"summary.tsv")))==2
        @test sum(D.i(r,"pixels") for r in D.t(joinpath(dir,"components.tsv")))==count(source.r.initial.mask)+count(source.r.final.mask)
    end
    @test_throws ErrorException G.analyse(fill(NaN,5,5),zeros(5,5),collect(1.:5),collect(1.:5),C)
    @test_throws ErrorException D.main(["--expected-N","6"])
    mktempdir() do dir
        for mutation in (c->(c["model"]["truth"]=0),c->(c["selection"]["connectivity"]=4),c->(c["preprocessing"]["stride"]=2))
            bad=deepcopy(C); mutation(bad); path=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException G.settings(path)
        end
    end
end

include(joinpath(@__DIR__,"verify_image_foreground.jl"))
const V=VerifyImageForeground
@testset "Saved synthetic verification and independent guard" begin
    rng=MersenneTwister(73)
    for rx in (0,1,7),ry in (0,2,7)
        mask=rand(rng,10,15).<.1
        @test V.guard_reference(mask,rx,ry)==G.box_guard(mask,rx,ry)
    end
    for case in synthetic_cases()
        mktempdir() do dir
            D.save_direction(dir,D.direction_data(case.image,"fwd",C),C)
            V.verify_direction(dir,C)
        end
    end
    a=(;file="fixture",view="fwd",identified=true,interval=1,mask="initial",group="retained",arm="local_plane",patches=2,test_pixels=10,test_sse_nm2=.01,rms_pm=1000sqrt(.001))
    rows=[a,merge(a,(view="bwd",identified=false,test_pixels=20,test_sse_nm2=.02))]
    summary=V.aggregate_predictions(rows)
    @test only(filter(r->r.subset=="all_scored" && r.interval==1 && r.mask=="initial" && r.group=="retained",summary)).test_pixels==30
    @test only(filter(r->r.subset=="identified" && r.interval==1 && r.mask=="initial" && r.group=="retained",summary)).test_pixels==10
end

if !isempty(ARGS)
    length(ARGS)==1 || error("Expected at most a new synthetic report path")
    ispath(only(ARGS)) && error("Synthetic report exists")
    D.R.table(only(ARGS),SYNTHETIC_METRICS)
end
