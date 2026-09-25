#!/usr/bin/env julia
using Test, LinearAlgebra, SparseArrays, Random, TOML
include(joinpath(@__DIR__,"diagnose_shared_image_envelope.jl"))
const D=DiagnoseSharedImageEnvelope
const F=D.F
const C=F.settings(joinpath(@__DIR__,"../config/shared_image_envelope.toml"))
BLAS.set_num_threads(1)

# Independent scalar recursion, including a right endpoint approached from the left.
function scalar_basis(x,t,j,d)
    d==0 && return t[j]<=x<t[j+1] ? 1. : 0.
    left=t[j+d]-t[j]; right=t[j+d+1]-t[j+1]
    (left>0 ? (x-t[j])/left*scalar_basis(x,t,j,d-1) : 0.)+
        (right>0 ? (t[j+d+1]-x)/right*scalar_basis(x,t,j+1,d-1) : 0.)
end

@testset "Cubic basis, boundaries and tensor order" begin
    for span in (.12,1.,2.7), width in C["model"]["minimum_knot_intervals_nm"]
        x=collect(range(-.23,span-.23;length=29)); a=F.basis(x,width); B=Matrix(a.matrix)
        @test all(B.>=0) && all(isapprox.(sum(B;dims=2),1.;atol=1e-14,rtol=0))
        @test a.segments==max(1,floor(Int,span/width))
        @test first(B[1,:])==1 && last(B[end,:])==1
        for r in eachindex(x),j in axes(B,2)
            xx=r==length(x) ? prevfloat(x[r]) : x[r]
            @test isapprox(B[r,j],scalar_basis(xx,a.knots,j,3);rtol=1e-11,atol=1e-13)
        end
        greville=[sum(a.knots[j+1:j+3])/3 for j in axes(B,2)]
        @test B*greville≈x
    end
    bx=F.basis(collect(range(0.,2.;length=15)),.64).matrix
    by=F.basis(collect(range(0.,1.;length=12)),.64).matrix
    coeff=reshape(collect(1.:size(bx,2)*size(by,2)),size(by,2),size(bx,2))
    @test reshape(kron(bx,by)*vec(coeff),12,15)≈by*coeff*transpose(bx)
    @test_throws ErrorException F.basis([0.,0.],1.)
end

@testset "Masked background and exact conditional trace" begin
    rng=MersenneTwister(431); xs=collect(range(0.,2.;length=19)); ys=collect(range(0.,1.;length=15))
    valid=trues(length(ys),length(xs)); valid[5,:].=false; valid[8,4]=false
    guard=falses(size(valid)); guard[4:12,6:15].=true
    b=F.background_operator(valid,guard,(xs.-1.)./2)
    Z=Matrix(b.Z); Zb=Z[b.bg,:]; B=Z*((transpose(Zb)*Zb)\transpose(Zb))
    Bfull=zeros(length(b.inds),length(b.inds)); Bfull[:,b.bg]=B
    v=randn(rng,length(b.inds))
    @test F.background(b,v)≈Bfull*v
    @test tr(Bfull)≈b.rank
    @test Bfull*Bfull≈Bfull
    z=v-F.background(b,v)
    for width in C["model"]["minimum_knot_intervals_nm"]
        f=F.fit_candidate(z,valid,xs,ys,b,width,C); X=Matrix(f.X)
        P=X*((transpose(X)*X)\transpose(X)); H=Bfull+P-P*Bfull
        @test f.beta≈X\z
        @test f.pred≈P*z
        @test f.trace_pb≈tr(P*Bfull)
        @test f.df≈tr(H)
        @test isapprox(f.rss,sum(abs2,v-H*v);rtol=1e-12,atol=1e-12)
        @test f.gcv≈(sum(abs2,v-H*v)/length(v))/(1-tr(H)/length(v))^2
        # Probe zero, arbitrary contrast and a representable broad field.
        for delta in (zeros(length(v)),randn(rng,length(v)),X*ones(size(X,2)))
            p=F.probe(delta,b,f); residual=(I-H)*delta
            @test isapprox(p.residual_sse,sum(abs2,residual);rtol=1e-10,atol=1e-20)
            @test p.beta≈X\((I-Bfull)*delta) atol=1e-12
        end
        @test F.probe(X*ones(size(X,2)),b,f).residual_sse<1e-20
    end
    r=F.fit_all(z,valid,xs,ys,b,C)
    @test all(==("ok"),r.statuses)
    @test r.chosen==argmin([f.gcv for f in r.fits])
    @test F.fit_all(z,valid,xs,ys,b,C).fits[r.chosen].beta==r.fits[r.chosen].beta
    @test_throws ErrorException F.background_operator(valid,trues(size(valid)),xs)
    # Missing edge band can contain entire basis supports. Zero columns change
    # no prediction on observed pixels and must not disqualify the finer grid.
    edgevalid=copy(valid); edgevalid[end-5:end,:].=false
    edge=F.background_operator(edgevalid,guard,(xs.-1.)./2)
    values=randn(rng,count(edgevalid)); corrected=values-F.background(edge,values)
    f=F.fit_candidate(corrected,edgevalid,xs,ys,edge,.32,C)
    @test f.parameters<(f.segments_x+3)*(f.segments_y+3)
    @test length(f.active)==f.parameters && all(vec(sum(abs2,f.X;dims=1)).>0)
    @test f.pred≈Matrix(f.X)*(Matrix(f.X)\corrected)
end

@testset "Synthetic envelope, missingness and score-only target" begin
    xs=collect(range(0.,2.;length=35)); ys=collect(range(0.,1.4;length=29)); ny,nx=length(ys),length(xs)
    valid=trues(ny,nx); valid[3,4]=false; guard=falses(ny,nx); guard[8:22,10:25].=true
    b=F.background_operator(valid,guard,(xs.-1.)./2)
    raw=[.04sin(2pi*y)+.02x+.1exp(-((x-1)^2+(y-.7)^2)/.14) for y in ys,x in xs]
    z=raw[valid]-F.background(b,raw[valid]); result=F.fit_all(z,valid,xs,ys,b,C)
    selected=result.fits[result.chosen]; @test selected.gcv<=result.fits[1].gcv
    @test selected.rss<.02result.fits[1].rss
    a=(;meta=Dict("xs_nm"=>xs,"ys_nm"=>ys,"source_reference_nm"=>3.),raw,b,valid,bg=zeros(count(valid)))
    other=(;raw=copy(raw),meta=Dict("xs_nm"=>xs,"ys_nm"=>ys,"source_reference_nm"=>3.))
    base=D.R.F.settings(joinpath(@__DIR__,"../config/vacuum_surface_transfer.toml"))
    # Use a small synthetic row block without changing real-image settings.
    base["preprocessing"]["calibration_block_rows"]=8
    cal=Dict("applied_dx_px"=>1,"difference_plane"=>[.1,.02,.03])
    target=D.target_vector(a,other,cal,base); held=D.R.F.row_masks(ny,base).test
    for (k,q) in enumerate(b.inds)
        y,x=Tuple(q)
        if held[y] && x<nx
            @test target[k]≈raw[y,x+1]-.1-.02xs[x]-.03ys[y]
        else; @test isnan(target[k]); end
    end
    before=copy(selected.beta); other.raw[held,:].+=50
    changed=D.target_vector(a,other,cal,base)
    @test any(abs.(changed[isfinite.(changed)]-target[isfinite.(target)]).>1)
    @test selected.beta==before
    @test F.fit_all(z,valid,xs,ys,b,C).fits[result.chosen].beta==before
    groups=C["selection"]["groups"]; rows=D.score(selected.pred,changed,guard[valid],groups)
    @test rows[1].test_pixels==rows[2].test_pixels+rows[3].test_pixels
    @test rows[1].test_sse_nm2≈rows[2].test_sse_nm2+rows[3].test_sse_nm2
    delta=zeros(length(z)); delta[findall(guard[valid])].=.01
    p=F.probe(delta,b,selected)
    @test 0<=p.residual_sse<p.input_sse
end

@testset "Config and submission boundary" begin
    mktempdir() do tmp
        for (s,k,value) in (("selection","policy","choose_target_winner"),("model","extra",6),
            ("preprocessing","missing","fill_missing"),("model","minimum_knot_intervals_nm",[.32,.64,1.28]))
            c=deepcopy(C); c[s][k]=value; path=joinpath(tmp,"bad.toml")
            open(io->TOML.print(io,c),path,"w"); @test_throws ErrorException F.settings(path)
        end
    end
    batch=read(joinpath(@__DIR__,"../hpc/shared_image_envelope.sbatch"),String)
    @test occursin("#SBATCH --no-requeue",batch) && occursin("--time=01:00:00",batch)
    @test occursin("-t 4",batch) && occursin("SLURM_JOB_ID:?",batch)
end
