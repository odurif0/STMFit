using Test, LinearAlgebra, Random, Statistics, Printf, TOML
include(joinpath(@__DIR__,"lib","empirical_fisher_native.jl"))
using .EmpiricalFisherNative
const EF=EmpiricalFisherNative
const ROOT=dirname(@__DIR__)
const CONTROL_CONFIG=joinpath(ROOT,"config","unit_assignment_affine_residual.toml")
const AFFINE_CONFIG=joinpath(ROOT,"config","unit_assignment_affine_fisher.toml")
const CONTROL=load_fisher_config(CONTROL_CONFIG)
const AFFINE=load_fisher_config(AFFINE_CONFIG)
const GRID=fisher_grid(AFFINE)
const BASIS=EF._affine_disk_basis(GRID)
const U=[GRID.coords[div(i-1,GRID.side)+1] for i in GRID.disk_indices]
const T=[GRID.coords[mod(i-1,GRID.side)+1] for i in GRID.disk_indices]
const DESIGN=hcat(ones(length(U)),T,U)

@testset "Explicit projection config and fixed geometric basis" begin
    @test CONTROL.patch_projection == "none"
    @test AFFINE.patch_projection == "affine_disk"
    @test AFFINE.projection_zero_l1 == 1e-12
    @test size(BASIS) == (197,3)
    @test norm(BASIS'*BASIS-I) <= size(BASIS,1)*eps(Float64) # dot products over 197 pixels
    for key in ("fisher_patch_projection","fisher_projection_zero_l1")
        section=key=="fisher_patch_projection" ? "preprocessing" : "model"
        for value in (nothing,"guess",true,-1.0,NaN,Inf)
            cfg=TOML.parsefile(AFFINE_CONFIG)
            value===nothing ? delete!(cfg[section],key) : (cfg[section][key]=value)
            @test_throws ArgumentError load_fisher_config(cfg)
        end
    end
end

@testset "Affine projection algebra, signal loss and fixed mirror" begin
    rng=MersenneTwister(29)
    p=randn(rng,length(U)) .+ 3U .- 2T .+ 0.7
    original=copy(p)
    q=EF._project_disk(p,BASIS)
    normal=p-DESIGN*((DESIGN'*DESIGN)\(DESIGN'*p))
    @test q ≈ normal atol=2e-13
    @test EF._project_disk(q,BASIS) ≈ q atol=2e-13
    @test maximum(abs,DESIGN'*q) < 5e-13
    @test norm(p)^2 ≈ norm(q)^2+norm(p-q)^2 atol=5e-12
    @test p == original
    for coefficients in ((0.3,1.2,-0.7),(-2.,-3.,5.),(0.,0.,3.))
        plane=DESIGN*collect(coefficients)
        @test sum(abs,EF._project_disk(plane,BASIS)) < AFFINE.projection_zero_l1
        @test EF._project_disk(p+plane,BASIS) ≈ q atol=2e-13
    end
    shoulder=exp.(-0.5 .* (((U .- 0.12)./0.055).^2 .+ ((T .- 0.04)./0.07).^2))
    @test 0 < norm(EF._project_disk(shoulder,BASIS)) < norm(shoulder)
    @test EF._project_disk(flip_u_disk(p,GRID),BASIS) ≈ flip_u_disk(q,GRID) atol=2e-13
end

function fixture()
    rng=MersenneTwister(20260920)
    n=64
    full=0.04randn(rng,n,GRID.side^2)
    coords=GRID.coords
    bump=[exp(-30(t*t+u*u))+0.12sin(13t)*cos(9u) for u in coords for t in coords]
    plane=[0.4+0.7u-0.2t for u in coords for t in coords]
    for i in 1:n
        full[i,:] .+= (i<=44 ? -2.0 : 2.0).*bump .+ (i/n).*plane
    end
    table=PatchTable([("synthetic.sxm",i) for i in 1:n],full[:,GRID.disk_indices],
                     full[:,GRID.center_index],fill("",n),GRID)
    table,full
end

@testset "Same transform at fit and score, original amplitude anchor, no class count" begin
    patches,full=fixture()
    prepared=patches.X-(patches.X*BASIS)*BASIS'
    original=copy(patches.X)
    model=fit_fisher(patches.X,patches.amplitudes,AFFINE)
    reference=fit_fisher(prepared,patches.amplitudes,CONTROL)
    @test model.projection_basis == BASIS
    @test reference.projection_basis === nothing
    @test model.w_p == reference.w_p
    @test model.mid == reference.mid
    @test model.amplitude_means == reference.amplitude_means
    @test model.gmm.weights == reference.gmm.weights
    @test model.gmm.weights[1] != model.gmm.weights[2] # no imposed 50/50 composition
    @test patches.X == original
    for i in (1,3,17,46,63)
        x=patches.X[i,:]
        q=EF._project_disk(x,BASIS)
        expected=dot(q-model.mid,model.w_p)
        @test score(x,model) == expected
        @test score(x,model) ≈ score(prepared[i,:],reference) atol=2e-13
        @test score(x+DESIGN*[0.3,-1.2,0.8],model) ≈ expected atol=2e-13
        @test maxmirror_score(x,model,GRID) ≈ max(score(q,reference),score(flip_u_disk(q,GRID),reference)) atol=2e-13
    end
    scores=cv_scores(patches,AFFINE)
    fixed=cv_scores(PatchTable(patches.keys,prepared,patches.amplitudes,patches.invalid_reasons,GRID),CONTROL)
    @test [s.invalid_reason for s in scores] == [s.invalid_reason for s in fixed]
    @test [s.score for s in scores] ≈ [s.score for s in fixed] atol=2e-12
    # A held-out row cannot affect other held-out scores in its fold.
    changed=copy(patches.X); changed[1,:] .+= 0.1sin.(13T)
    altered=cv_scores(PatchTable(patches.keys,changed,patches.amplitudes,patches.invalid_reasons,GRID),AFFINE)
    @test [s.score for s in scores[3:2:end]] == [s.score for s in altered[3:2:end]]
    renamed=PatchTable([("renamed.sxm",i) for i in 1:64],patches.X,patches.amplitudes,patches.invalid_reasons,GRID)
    @test [s.score for s in cv_scores(renamed,AFFINE)] == [s.score for s in scores]
    planes=copy(patches.X); planes[2,:] = DESIGN*[0.7,1.2,-0.3]; planes[3,1]=NaN
    bad=cv_scores(PatchTable(patches.keys,planes,patches.amplitudes,patches.invalid_reasons,GRID),AFFINE)
    @test [(s.file,s.lobe) for s in bad] == patches.keys
    @test bad[2].invalid_reason == "zero_affine_patch_mass"
    @test bad[3].invalid_reason == "nonfinite_patch_input"
    @test isnan(bad[2].score) && isnan(bad[3].score)
    pure=repeat(reshape(DESIGN*[0.7,1.2,-0.3],1,:),64,1)
    allbad=cv_scores(PatchTable(patches.keys,pure,patches.amplitudes,patches.invalid_reasons,GRID),AFFINE)
    @test all(s->s.invalid_reason=="zero_affine_patch_mass" && isnan(s.score),allbad)
    @test_throws EF.FisherFitError fit_fisher(pure,patches.amplitudes,AFFINE)
    @test_throws DimensionMismatch fit_fisher(patches.X[:,1:20],patches.amplitudes,AFFINE)
    mktempdir() do dir
        path=joinpath(dir,"patches.tsv")
        open(path,"w") do io
            println(io,join(vcat(["file","lobe","amplitude"],[@sprintf("res_p%03d",i) for i in 1:289]),'\t'))
            for i in 1:64
                # Irrelevant metadata amplitude must not replace the original center pixel.
                println(io,join(vcat(["synthetic.sxm",string(i),"999"],string.(full[i,:])),'\t'))
            end
        end
        loaded=load_patches(path,"res",AFFINE)
        @test loaded.X == patches.X
        @test loaded.amplitudes == patches.amplitudes
        expected=joinpath(dir,"expected.tsv"); actual=joinpath(dir,"actual.tsv")
        write_scores(expected,scores)
        cmd=`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_empirical_fisher_native.jl")) --patches $path --prefix res --config $AFFINE_CONFIG --out $actual`
        @test success(cmd)
        @test read(expected) == read(actual)
        @test_throws ArgumentError write_scores(actual,scores)
    end
end
