# Includes the upstream regression and its independent synthetic SXM fixture.
include(joinpath(@__DIR__, "test_vacuum_surface_transfer.jl"))
include(joinpath(@__DIR__, "diagnose_vacuum_surface_controls.jl"))
include(joinpath(@__DIR__, "verify_vacuum_surface_controls.jl"))
const L = DiagnoseVacuumSurfaceControls
const U = L.C
const Q = VerifyVacuumSurfaceControls
const CONTROLS_PATH = joinpath(@__DIR__, "../config/vacuum_surface_controls.toml")
const CONTROL = U.settings(CONTROLS_PATH)

function control_fixture()
    grid=collect(-.32:.04:.32)
    common=[1+.03exp(-(t/.12)^2-(u/.13)^2) for u in grid,t in grid]
    contrast=[.025+.012exp(-((t-.07)/.06)^2-((u+.04)/.09)^2) for u in grid,t in grid]
    maps=Dict(-1=>common,0=>common-contrast,1=>common+contrast)
    coords=[(x,y) for x in -.15:.03:.15 for y in -.15:.03:.15 if hypot(x,y)<=.16]
    p=(;anchor=1,dx=first.(coords),dy=last.(coords),values=zeros(length(coords)))
    small=deepcopy(C); small["model"]["rotation_step_deg"]=360.; small["model"]["translations_nm"]=[0.]
    state=only(F.states(small,"common"))
    pp=[merge(p,(;anchor=k+1,values=F.prediction(p,maps[k],state,small) .- 0.96 .+ 0.03 .* p.dx .- 0.07 .* p.dy)) for k in 0:1]
    (;maps,p,small,pp)
end

@testset "Matched height contrast and source-only local plane" begin
    f=control_fixture(); r=U.fit_patches(f.pp,f.maps,f.small,CONTROL)
    @test r.delta≈mean((f.maps[1]-f.maps[0])/2)
    @test [p.state.type for p in r.fits["plane_chemical"]]==[0,1]
    @test all(norm(p.values-fit.pred)<1e-12 for (p,fit) in zip(f.pp,r.fits["plane_chemical"]))
    @test all(sum(abs2,p.values-fit.pred)>1e-5 for (p,fit) in zip(f.pp,r.fits["plane_common"]))
    for (j,p) in enumerate(f.pp)
        X=hcat(ones(length(p.values)),p.dx,p.dy)
        basis=U.plane_basis(p,1e-12)
        beta=X\p.values
        @test basis.X*U.remove_plane(p.values,basis).beta≈X*beta
        @test r.fits["local_plane"][j].pred≈X*beta
        @test r.fits["local_constant"][j].pred==fill(mean(p.values),length(p.values))
        for arm in ("plane_common","plane_chemical")
            fit=r.fits[arm][j]
            @test norm(X'*(p.values-fit.pred))<1e-11
            @test (fit.state.angle,fit.state.tx,fit.state.ty)==(0.,0.,0.)
        end
        for (k,state) in enumerate(r.hs)
            v=F.prediction(p,f.maps[-1],state,f.small).+(2state.type-1)*r.delta
            resid=p.values-v; center=mean(resid)
            @test isapprox(center,r.height_costs[j].mu[k];atol=1e-14)
            @test isapprox(sum(abs2,resid.-center),r.height_costs[j].v[k];atol=1e-18,rtol=1e-11)
        end
    end
    expected=Inf
    for types in Iterators.product(1:2,1:2)
        cost=r.height_costs; b=sum(c.n*c.mu[k] for (c,k) in zip(cost,types))/sum(c.n for c in cost)
        expected=min(expected,sum(c.v[k]+c.n*(b-c.mu[k])^2 for (c,k) in zip(cost,types)))
    end
    @test r.height.loss≈expected
    swapped=Dict(-1=>f.maps[-1],0=>f.maps[1],1=>f.maps[0])
    reverse=U.fit_patches(f.pp,swapped,f.small,CONTROL)
    for arm in ("height_contrast","plane_chemical")
        @test all(isapprox(a.pred,b.pred;atol=1e-13) for (a,b) in zip(r.fits[arm],reverse.fits[arm]))
        @test getproperty.(getproperty.(reverse.fits[arm],:state),:type)==1 .-getproperty.(getproperty.(r.fits[arm],:state),:type)
    end
    # A nuisance plane cannot create a new chemical contrast.
    shifted=[merge(p,(;values=p.values .+ 0.17 .+ 0.4 .* p.dx .- 0.3 .* p.dy)) for p in f.pp]
    again=U.fit_patches(shifted,f.maps,f.small,CONTROL)
    for arm in ("plane_common","plane_chemical"), j in eachindex(f.pp)
        @test r.fits[arm][j].state==again.fits[arm][j].state
        @test again.fits[arm][j].pred-r.fits[arm][j].pred ≈ 0.17 .+ 0.4 .* f.pp[j].dx .- 0.3 .* f.pp[j].dy
    end
    @test_throws ErrorException U.plane_basis((;values=ones(5),dx=ones(5),dy=zeros(5)),1e-12)
    @test_throws ErrorException U.fit_patches(NamedTuple[],f.maps,f.small,CONTROL)
end

@testset "Pure-height contrast disappears under the matched plane nuisance" begin
    f=control_fixture(); maps=Dict(-1=>f.maps[-1],0=>f.maps[-1].-.04,1=>f.maps[-1].+.04)
    result=U.fit_patches(f.pp,maps,f.small,CONTROL)
    for j in eachindex(f.pp)
        @test isapprox(result.fits["plane_common"][j].pred,result.fits["plane_chemical"][j].pred;atol=1e-13)
        @test isapprox(result.chemical_costs[j].losses[1],result.chemical_costs[j].losses[2];atol=1e-16)
    end
end

@testset "Frozen observation control workflow and target mutation" begin
    f=control_fixture(); maps=(heights=Dict(j=>f.maps for j in 1:3),isos=Dict(j=>Float64(j) for j in 1:3))
    mktempdir() do dir
        upstream=joinpath(dir,"upstream"); R.one_direction(toy_image(),"fwd",maps,C,upstream)
        out=joinpath(dir,"controls"); result=L.one_direction(upstream,maps,C,CONTROL,out)
        verified=Q.verify_direction(out,upstream,maps,C,CONTROL)
        @test verified.patches==result.patches && length(verified.models)==15
        rows=L.t(joinpath(out,"models.tsv")); patches=L.t(joinpath(out,"patches.tsv"))
        @test result.status=="ok" && result.patches>=2
        @test length(rows)==15 && length(patches)==15result.patches
        @test Set(r["arm"] for r in rows)==Set(CONTROL["model"]["arms"])
        # Full discrete common geometry is replayed, not silently replaced by a small test grid.
        @test all(isfinite(L.n(r,"training_sse_nm2")) for r in rows)
        observed=R.V.tsv(joinpath(upstream,"observations.tsv"))
        changed=[(;patch=L.i(r,"patch"),anchor=L.i(r,"anchor"),pixel=L.i(r,"pixel"),row=L.i(r,"row"),column=L.i(r,"column"),
            dx_nm=L.n(r,"dx_nm"),dy_nm=L.n(r,"dy_nm"),source_nm=L.n(r,"source_nm"),target_nm=L.n(r,"target_nm")+10.) for r in observed]
        R.table(joinpath(upstream,"observations.tsv"),changed)
        other=joinpath(dir,"target_changed"); L.one_direction(upstream,maps,C,CONTROL,other)
        for file in ("common_replay.tsv",("costs_$(j)_$(a).tsv" for j in 1:3 for a in ("height_contrast","plane_common","plane_chemical"))...)
            @test read(joinpath(out,file))==read(joinpath(other,file))
        end
        changedrows=L.t(joinpath(other,"models.tsv")); changedpatches=L.t(joinpath(other,"patches.tsv"))
        for (a,b) in zip(patches,changedpatches), key in keys(a)
            key=="test_sse_nm2" || @test a[key]==b[key]
        end
        @test all(L.n(b,"test_sse_nm2")>L.n(a,"test_sse_nm2") for (a,b) in zip(rows,changedrows))
        @test_throws ErrorException L.one_direction(upstream,maps,C,CONTROL,out)
    end
    @test_throws ErrorException L.main(["--truth","forbidden"])
    mktempdir() do dir
        for mutation in (c->(c["model"]["expected_N"]=6),c->(c["model"]["height_gain"]=2),
                         c->(c["preprocessing"]["target_use"]="fit"))
            bad=deepcopy(CONTROL); mutation(bad); file=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,bad),file,"w")
            @test_throws ErrorException U.settings(file)
        end
    end
end
