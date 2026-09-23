using Test, LinearAlgebra, Random, Statistics, TOML
include(joinpath(@__DIR__,"score_tangent_mold_templates.jl"))
const S=TangentMoldScorer
const T=S.T
const ROOT=dirname(@__DIR__)
const SETTINGS=joinpath(ROOT,"config/tangent_mold_projection.toml")
const OPT=T.settings(SETTINGS)
using STMSXMIO: _box_smooth
BLAS.set_num_threads(1)

function geometry(;source="ell",x=0.5,y=0.5,angle=0.0)
    Dict("file"=>"synthetic.sxm","N"=>"1","lobe"=>"1","x_nm"=>string(x),"y_nm"=>string(y),
        "axis_x"=>string(cos(angle)),"axis_y"=>string(sin(angle)),"source"=>source,
        "sigma_parallel_nm"=>"0.18","sigma_perp_nm"=>source=="circ" ? "0.18" : "0.14","skew_ratio"=>"1.0")
end

# Independent full-image construction and bilinear formula (not the optimized
# sparse native sampler). Also test the clipped box window at the image edge.
function full_reference(xs,ys,r,coords,radius; orientation=false)
    cx=T.number(r,"x_nm"); cy=T.number(r,"y_nm"); ax=T.number(r,"axis_x"); ay=T.number(r,"axis_y")
    sp=T.number(r,"sigma_parallel_nm"); su=T.number(r,"sigma_perp_nm")
    x=reshape(xs,1,:); y=reshape(ys,:,1)
    t=(x.-cx).*ax.+(y.-cy).*ay; u=-(x.-cx).*ay.+(y.-cy).*ax
    g=exp.(-0.5 .* ((t./sp).^2 .+ (u./su).^2))
    maps=[ones(size(g)),t,u,g,t./sp.*g,u./su.*g]
    r["source"]=="circ" ? push!(maps,((t./sp).^2 .+ (u./su).^2).*g) : append!(maps,[(t./sp).^2 .*g,(u./su).^2 .*g])
    if orientation && r["source"]=="ell" && sp!=su
        push!(maps,t.*u.*g.*(1/su^2-1/sp^2))
    end
    maps=_box_smooth.(maps,radius)
    out=fill(NaN,length(coords)^2,length(maps))
    for (k,(vv,tt)) in enumerate(Iterators.product(coords,coords))
        # product first coordinate is fastest: physical t, then u.
        tp=vv; up=tt; xx=cx+tp*ax-up*ay; yy=cy+tp*ay+up*ax
        i=searchsortedlast(xs,xx); j=searchsortedlast(ys,yy)
        1<=i<length(xs) && 1<=j<length(ys) || continue
        wx=(xx-xs[i])/(xs[i+1]-xs[i]); wy=(yy-ys[j])/(ys[j+1]-ys[j])
        for q in eachindex(maps)
            m=maps[q]; out[k,q]=(1-wx)*(1-wy)*m[j,i]+wx*(1-wy)*m[j,i+1]+(1-wx)*wy*m[j+1,i]+wx*wy*m[j+1,i+1]
        end
    end
    out
end

@testset "Orientation derivative: fixed pixels/frame, exact isotropy and unchanged base" begin
    opt=T.settings(joinpath(ROOT,"config/tangent_mold_orientation.toml"))
    @test opt.orientation && !OPT.orientation
    @test Base.structdiff(opt,(orientation=true,))==Base.structdiff(OPT,(orientation=false,))
    xs=collect(range(0,1,length=43)); ys=collect(range(0,1,length=39)); coords=collect(-.32:.04:.32)
    rng=MersenneTwister(924)
    for source in ("ell","circ"), angle in (0.0,0.63,-1.2), radius in (0,1,2), center in ((.5,.5),(.07,.12))
        r=geometry(;source,angle,x=center[1],y=center[2])
        B=T.native_design(xs,ys,r,coords,radius;orientation=true)
        R=full_reference(xs,ys,r,coords,radius;orientation=true)
        old=T.native_design(xs,ys,r,coords,radius)
        @test isequal(B[:,1:size(old,2)],old)
        @test size(B)==(289,source=="ell" ? 9 : 7)
        @test isnan.(B)==isnan.(R)
        good=isfinite.(R)
        @test B[good]≈R[good] rtol=3e-14 atol=3e-15
        rows=findall(isfinite,B[:,1]); D=B[rows,:]
        p=T.projector(D,opt)
        @test p.rank==size(D,2)
        y=randn(rng,length(rows)); c=randn(rng,length(rows)); M=hcat(-c,c)
        a=T.projected_score(y,M,D,opt)
        Q=Matrix(qr(D).Q)[:,1:size(D,2)]
        yp=y-Q*(Q'*y); mp=M-Q*(Q'*M)
        @test a.costs≈[-dot(yp,mp[:,j])/(norm(yp)*norm(mp[:,j])) for j in 1:2] atol=2e-11
        @test T.projected_score(y+D*randn(rng,size(D,2)),M+D*randn(rng,size(D,2),2),D,opt).costs≈a.costs atol=2e-11
        @test norm(T.project(D,p))<2e-12
    end
    # Vary only the model angle, never the fixed observation point.
    function gaussian(x,y,r,angle)
        dx=x-T.number(r,"x_nm"); dy=y-T.number(r,"y_nm")
        t=dx*cos(angle)+dy*sin(angle); u=-dx*sin(angle)+dy*cos(angle)
        exp(-.5*((t/T.number(r,"sigma_parallel_nm"))^2+(u/T.number(r,"sigma_perp_nm"))^2))
    end
    h=1e-5
    for angle in (0.0,.63,-1.2,pi), su in (.14,.22), (x,y) in ((.56,.61),(.3,.4),(.5,.5))
        r=merge(geometry(;angle),Dict("sigma_perp_nm"=>string(su)))
        b=T.basis_at(x,y,r;orientation=true)
        @test (gaussian(x,y,r,angle+h)-gaussian(x,y,r,angle-h))/(2h)≈b[9] rtol=2e-8 atol=1e-10
    end
    for source in ("ell","circ")
        r=merge(geometry(;source,angle=.63),Dict("sigma_perp_nm"=>"0.18"))
        @test !T.has_orientation(r)
        @test isequal(T.native_design(xs,ys,r,coords,1;orientation=true),T.native_design(xs,ys,r,coords,1))
    end
    near=merge(geometry(),Dict("sigma_perp_nm"=>string(nextfloat(.18))))
    B=T.native_design(xs,ys,near,coords,1;orientation=true)
    @test T.has_orientation(near) && size(B,2)==9
    @test norm(B[:,9])>0 && T.projector(B,opt).rank==9
    r=geometry(); B=T.native_design(xs,ys,r,coords,1;orientation=true)
    c=randn(rng,289); M=hcat(-c,c); y=randn(rng,289); y[1:40].=NaN
    a=T.projected_score(y,M,B,opt); old=T.projected_score(y,M,B[:,1:8],OPT)
    @test a.contrast_fraction<=old.contrast_fraction+1e-12
    @test a.patch_fraction<=old.patch_fraction+1e-12
    changed=copy(B); changed[1:40,:].=Inf
    @test T.projected_score(y,M,changed,opt).costs==a.costs
    @test T.projected_score(B*ones(9),M,B,opt).reason=="zero_projected_patch"
    @test_throws ErrorException T.projected_score(randn(rng,289),hcat(-B[:,9],B[:,9]),B,opt)
end

@testset "Exact native smoothed/sampled tangent span" begin
    xs=collect(range(0,1,length=43)); ys=collect(range(0,1,length=39)); coords=collect(-.32:.04:.32)
    for source in ("ell","circ"), angle in (0.0,0.63,-1.2), radius in (0,1,2), center in ((.5,.5),(.07,.12))
        r=geometry(;source,angle,x=center[1],y=center[2])
        B=T.native_design(xs,ys,r,coords,radius)
        R=full_reference(xs,ys,r,coords,radius)
        @test size(B)==size(R)==(289,source=="circ" ? 7 : 8)
        @test isnan.(B)==isnan.(R)
        good=isfinite.(R)
        @test B[good]≈R[good] rtol=2e-14 atol=2e-15
        @test all(B[isfinite.(B[:,1]),1].≈1)
    end
    r=geometry(); x=.56; y=.61; b=T.basis_at(x,y,r); h=1e-5
    g(r)=exp(-.5*((x-T.number(r,"x_nm"))^2/T.number(r,"sigma_parallel_nm")^2+(y-T.number(r,"y_nm"))^2/T.number(r,"sigma_perp_nm")^2))
    for (field,col,scale) in (("x_nm",5,.18),("y_nm",6,.14))
        a=copy(r); z=copy(r); a[field]=string(T.number(r,field)+h); z[field]=string(T.number(r,field)-h)
        @test (g(a)-g(z))/(2h)*scale≈b[col] rtol=1e-8
    end
    for (field,col) in (("sigma_parallel_nm",7),("sigma_perp_nm",8))
        a=copy(r); z=copy(r); a[field]=string(T.number(r,field)*exp(h)); z[field]=string(T.number(r,field)*exp(-h))
        @test (g(a)-g(z))/(2h)≈b[col] rtol=1e-8
    end
    @test_throws ErrorException T.native_design(xs,ys,merge(r,Dict("skew_ratio"=>"1.1")),coords,1)
    @test_throws ErrorException T.native_design(xs,ys,merge(r,Dict("source"=>"bad")),coords,1)
    @test_throws ErrorException T.native_design(xs,ys,merge(r,Dict("axis_x"=>"1.1")),coords,1)
    @test_throws ErrorException T.native_design(xs,ys,merge(r,Dict("sigma_perp_nm"=>"0")),coords,1)
end

@testset "Common support, numerical rank, invariance and contrast retention" begin
    rng=MersenneTwister(923); coords=collect(-.32:.04:.32)
    B=T.native_design(collect(0:.02:1),collect(0:.025:1),geometry(),coords,1)
    y=randn(rng,289); c=randn(rng,289); M=hcat(-c,c)
    p=T.projector(B,OPT); P=I-p.Q*p.Q'
    @test p.rank==8
    @test P*P≈P atol=1e-14
    @test P≈P' atol=1e-14
    @test norm(P*B)<1e-12
    @test T.project(y,p)≈P*y atol=1e-14
    a=T.projected_score(y,M,B,OPT)
    @test a.reason=="ok" && a.rank==8 && a.n==289
    @test 0<a.patch_fraction<=1 && 0<a.contrast_fraction<=1
    @test a.costs≈[-dot(P*y,P*M[:,j])/(norm(P*y)*norm(P*M[:,j])) for j in 1:2] atol=1e-14
    @test a.costs[1]≈-a.costs[2]
    for _ in 1:20
        b=T.projected_score(3y+B*randn(rng,8),M+B*randn(rng,8,2),B,OPT)
        @test b.costs≈a.costs atol=2e-13
        perm=randperm(rng,289)
        @test T.projected_score(y[perm],M[perm,:],B[perm,:],OPT).costs≈a.costs atol=2e-13
    end
    # QR is an independent construction when the finite-support design is full rank.
    partial=copy(y); partial[1:40].=NaN
    good=findall(isfinite,partial); Q=Matrix(qr(B[good,:]).Q)[:,1:8]
    yp=partial[good]-Q*(Q'*partial[good]); mp=M[good,:]-Q*(Q'*M[good,:])
    a=T.projected_score(partial,M,B,OPT)
    @test a.costs≈[-dot(yp,mp[:,j])/(norm(yp)*norm(mp[:,j])) for j in 1:2] atol=2e-13
    changed=copy(B); changed[1:40,:].=Inf
    @test T.projected_score(partial,M,changed,OPT).costs==a.costs
    partial[1:145].=NaN
    @test T.projected_score(partial,M,B,OPT).reason=="insufficient_native_support"
    partial[145]=y[145]
    @test T.projected_score(partial,M,B,OPT).n==145
    @test T.projected_score(B*ones(8),M,B,OPT).reason=="zero_projected_patch"
    @test_throws ErrorException T.projected_score(y,hcat(-B[:,2],B[:,2]),B,OPT)
    @test_throws ErrorException T.projected_score(y,fill(NaN,289,2),B,OPT)
    @test T.projector(hcat(B,B[:,4]),OPT).rank==8
    @test_throws ErrorException T.projector(zeros(20,8),OPT)
end

@testset "Strict configuration and compatible connected decoder" begin
    mktempdir() do dir
        for (section,key,value) in (("model","rank_rtol",0.0),("model","zero_norm_rtol",NaN),
            ("model","basis","target_gaussian_affine_neighbors"),
            ("selection","minimum_observed_fraction",0.4),("preprocessing","sampling","analytic_only"),("model","truth","forbidden"))
            cfg=TOML.parsefile(SETTINGS); cfg[section][key]=value
            p=joinpath(dir,"bad.toml"); open(io->TOML.print(io,cfg),p,"w")
            @test_throws ErrorException T.settings(p)
        end
    end
    @test_throws ErrorException S.parse_cli(["--truth","x"])
    @test_throws ErrorException S.parse_cli(["--expected-N","6"])
    @test S.parse_cli(["--help"])===nothing
    rng=MersenneTwister(15)
    records=[(file="s.sxm",lobe=j,amplitude=.1j,patch=randn(rng,25)) for j in 1:5]
    templates=Dict((t,p,m)=>randn(rng,25) for t in 0:1 for p in 0:1 for m in 0:1)
    opt=S.Connected.Options("","",nothing,"","res","ncc","contrast",0.,1.)
    a=S.Connected._decode_file(records,templates,nothing,opt)
    b=S.Connected._decode_file(records,templates,nothing,opt;scorer=(r,t)->S.Connected._score_patch(r.patch,t,"ncc"))
    @test a==b
    io=IOBuffer(); S.Connected.write_decoded(io,"s.sxm",records,a)
    @test length(split(chomp(String(take!(io))),'\n'))==5
end
