# Includes the unchanged target/orientation regression suite.
include(joinpath(@__DIR__,"test_tangent_mold_projection.jl"))
const FINITE=T.settings(joinpath(ROOT,"config/tangent_mold_finite.toml"))
const NEIGHBORS=T.settings(joinpath(ROOT,"config/tangent_mold_neighbors.toml"))

# Independent full native map for a neighboring Gaussian, sampled in target frame.
function neighbor_reference(xs,ys,target,neighbor,coords,radius)
    nx=T.number(neighbor,"x_nm"); ny=T.number(neighbor,"y_nm")
    ax=T.number(neighbor,"axis_x"); ay=T.number(neighbor,"axis_y")
    sp=T.number(neighbor,"sigma_parallel_nm"); su=T.number(neighbor,"sigma_perp_nm")
    G=[exp(-.5*(((x-nx)*ax+(y-ny)*ay)^2/sp^2+(-(x-nx)*ay+(y-ny)*ax)^2/su^2)) for y in ys,x in xs]
    smoothed=_box_smooth(G,radius)
    cx=T.number(target,"x_nm"); cy=T.number(target,"y_nm")
    tx=T.number(target,"axis_x"); ty=T.number(target,"axis_y")
    out=fill(NaN,length(coords)^2)
    k=0
    for u in coords,t in coords
        k+=1; x=cx+t*tx-u*ty; y=cy+t*ty+u*tx
        i=searchsortedlast(xs,x); j=searchsortedlast(ys,y)
        1<=i<length(xs) && 1<=j<length(ys) || continue
        wx=(x-xs[i])/(xs[i+1]-xs[i]); wy=(y-ys[j])/(ys[j+1]-ys[j])
        out[k]=(1-wx)*(1-wy)*smoothed[j,i]+wx*(1-wy)*smoothed[j,i+1]+
            (1-wx)*wy*smoothed[j+1,i]+wx*wy*smoothed[j+1,i+1]
    end
    out
end

@testset "Only adjacent amplitude columns; fixed old span and native target sampling" begin
    @test NEIGHBORS.adjacent_amplitudes && !NEIGHBORS.orientation && NEIGHBORS.omit_unavailable
    @test FINITE.omit_unavailable && !FINITE.adjacent_amplitudes && !FINITE.orientation
    @test Base.structdiff(NEIGHBORS,(adjacent_amplitudes=true,))==Base.structdiff(FINITE,(adjacent_amplitudes=false,))
    @test Base.structdiff(FINITE,(omit_unavailable=true,))==Base.structdiff(OPT,(omit_unavailable=false,))
    xs=collect(range(0,1,length=43)); ys=collect(range(0,1,length=39)); coords=collect(-.32:.04:.32)
    rng=MersenneTwister(92321)
    for source in ("ell","circ"), angle in (0.,.63,-1.2), radius in (0,1,2), center in ((.5,.5),(.07,.12))
        r=geometry(;source,angle,x=center[1],y=center[2])
        p=source=="circ" ? 7 : 8
        neighbors=[geometry(;source="ell",angle=angle+.2,x=center[1]-.25,y=center[2]+.05),
            geometry(;source="circ",angle=angle-.3,x=center[1]+.29,y=center[2]-.07)]
        old=T.native_design(xs,ys,r,coords,radius)
        for nn in 0:2
            ns=neighbors[1:nn]
            B=T.native_design(xs,ys,r,coords,radius;neighbors=ns)
            @test size(B)==(289,p+nn)
            @test isequal(B[:,1:p],old)
            for j in 1:nn
                ref=neighbor_reference(xs,ys,r,ns[j],coords,radius)
                @test isnan.(B[:,p+j])==isnan.(ref)
                good=findall(isfinite,ref)
                @test B[good,p+j]≈ref[good] rtol=4e-14 atol=4e-15
            end
            rows=findall(isfinite,B[:,1]); D=B[rows,:]
            y=randn(rng,length(rows)); c=randn(rng,length(rows)); molds=hcat(-c,c)
            result=T.projected_score(y,molds,D,NEIGHBORS)
            @test result.reason=="ok" && result.rank==p+nn
            Q=Matrix(qr(D).Q)[:,1:p+nn]; yp=y-Q*(Q'*y); Mp=molds-Q*(Q'*molds)
            @test result.costs≈[-dot(yp,Mp[:,j])/(norm(yp)*norm(Mp[:,j])) for j in 1:2] atol=2e-10
            @test T.projected_score(y+D*randn(rng,p+nn),molds+D*randn(rng,p+nn,2),D,NEIGHBORS).costs≈result.costs atol=2e-10
            prior=T.projected_score(y,molds,D[:,1:p],FINITE)
            @test result.patch_fraction<=prior.patch_fraction+2e-12
            @test result.contrast_fraction<=prior.contrast_fraction+2e-12
        end
    end
    r=geometry(); q=geometry(;x=.72,y=.56,angle=.4)
    # A duplicate amplitude is rank-redundant, not a new numerical direction.
    B=T.native_design(xs,ys,r,coords,1;neighbors=[r])
    @test T.projector(B,NEIGHBORS).rank==8
    @test_throws ErrorException T.native_design(xs,ys,r,coords,1;neighbors=[q,q,q])
    @test_throws ErrorException T.native_design(xs,ys,r,coords,1;orientation=true,neighbors=[q])
    @test_throws ErrorException T.native_design(xs,ys,r,coords,1;neighbors=[merge(q,Dict("sigma_parallel_nm"=>"0"))])
    for n in 1:5
        base=Dict(("a.sxm",i)=>merge(geometry(),Dict("lobe"=>string(i),"N"=>string(n))) for i in 1:n)
        for i in 1:n
            expected=filter(j->1<=j<=n,[i-1,i+1])
            @test parse.(Int,getindex.(T.adjacent_rows(base,"a.sxm",i,n),"lobe"))==expected
        end
    end
    @test_throws ErrorException T.adjacent_rows(Dict(),"a.sxm",0,3)
    mktempdir() do dir
        c=TOML.parsefile(joinpath(ROOT,"config/tangent_mold_neighbors.toml"))
        c["selection"]["missing_cost"]="fill_zero_scores"
        path=joinpath(dir,"bad.toml"); open(io->TOML.print(io,c),path,"w")
        @test_throws ErrorException T.settings(path)
    end
end
