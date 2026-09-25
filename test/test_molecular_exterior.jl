#!/usr/bin/env julia
# No arguments: synthetic geometry. ROOT GEOMETRY [REPORT]: saved-output checks.
using Test, Random, Statistics
include(joinpath(@__DIR__,"lib/molecular_exterior.jl"))
const E=MolecularExterior
const Q=Rational{BigInt}
cross(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2])
norm2(a,b)=sum((a[k]-b[k])^2 for k in 1:2)

# Gift wrapping, not Andrew chains; ties choose the farthest collinear point.
function oracle_hull(points)
    ps=unique([(Q(p[1]),Q(p[2])) for p in points]); firstpoint=minimum(ps)
    result=typeof(firstpoint)[]; p=firstpoint
    while true
        push!(result,p); q=first(x for x in ps if x!=p)
        for r in ps
            turn=cross(p,q,r)
            (turn<0 || (turn==0 && norm2(p,r)>norm2(p,q))) && (q=r)
        end
        p=q; p==firstpoint && break
        length(result)>length(ps) && error("Oracle failed to close")
    end
    result
end

# Projection distance uses the expanded norm identity, not projected coordinates.
function oracle_distance(p,poly)
    inside=true; nearest=nothing
    for j in eachindex(poly)
        a=poly[j]; b=poly[mod1(j+1,length(poly))]; inside &= cross(a,b,p)>=0
        dot=sum((p[k]-a[k])*(b[k]-a[k]) for k in 1:2); length2=norm2(a,b)
        d=dot<=0 ? norm2(p,a) : dot>=length2 ? norm2(p,b) : norm2(p,a)-dot^2/length2
        nearest=nearest===nothing ? d : min(nearest,d)
    end
    inside ? Q(0) : nearest
end

function oracle_clearance(xy,poly,radius,cell)
    p=(mod(Q(xy[1]),Q(cell[1])),mod(Q(xy[2]),Q(cell[2])))
    # An intentionally larger explicit image enumeration than the implementation.
    minimum(oracle_distance((p[1]+ix*Q(cell[1]),p[2]+iy*Q(cell[2])),poly)
        for ix in -2:2 for iy in -2:2)-Q(radius)^2
end

@testset "Exact whole-molecule footprint geometry" begin
    square=[(2.,2.),(4.,2.),(4.,4.),(2.,4.)]
    expected=E.point.(square)
    @test E.hull(vcat(square,[(3.,3.),(3.,2.)],square))==expected
    @test E.hull(reverse(square))==expected==oracle_hull(square)
    @test_throws ErrorException E.hull([(0.,0.),(1.,0.),(2.,0.)])
    @test_throws ErrorException E.hull([(NaN,0.),(1.,0.),(0.,1.)])
    @test_throws ErrorException E.hull([(0.,0.),(1.,0.)])
    tiny=[(0.,0.),(1.,0.),(1.,nextfloat(0.))]
    @test E.hull(tiny)==oracle_hull(tiny) # no tolerance discards tiny orientation
    atoms=[(;z=6,position_nm=[p...,0.]) for p in square]
    push!(atoms,(;z=29,position_nm=[7.,7.,0.]))
    geo=(;cell_nm=[8.,8.,10.],atoms); f=E.footprint(geo,Dict(6=>.25,29=>3.))
    @test f.radius==1//4 && f.poly==expected && length(f.atoms)==4
    @test E.clearance((1.75,3.),f).boundary
    @test !E.clearance((1.75,3.),f).exterior
    @test E.clearance((prevfloat(1.75),3.),f).exterior
    @test !E.clearance((nextfloat(1.75),3.),f).exterior
    @test E.distance2(E.point((1.,1.)),f.poly)==2
    @test E.distance2(E.point((3.,3.)),f.poly)==0
    @test_throws ErrorException E.footprint((;cell_nm=[4.1,8.,10.],atoms),Dict(6=>.25,29=>3.))
    @test_throws ErrorException E.footprint((;cell_nm=[8.,8.,10.],atoms=atoms[end:end]),Dict(29=>3.))
    rng=Xoshiro(831974)
    for _ in 1:50
        points=[Tuple(rand(rng,0:128,2)./16) for _ in 1:20]
        poly=E.hull(points)
        @test poly==oracle_hull(points)==E.hull(shuffle(rng,points))
        for _ in 1:20
            p=E.point(Tuple(rand(rng,-64:192,2)./16))
            @test E.distance2(p,poly)==oracle_distance(p,poly)
        end
    end
    for _ in 1:200
        xy=E.point(Tuple(rand(rng,-128:256,2)./16)); c=E.clearance(xy,f)
        @test c.exterior==(oracle_clearance(xy,f.poly,f.radius,f.cell)>0)
        @test E.clearance((xy[1]+3f.cell[1],xy[2]-2f.cell[2]),f)==c
    end
end

if !isempty(ARGS)
    length(ARGS) in (2,3) || error("test_molecular_exterior.jl [ROOT SAVED_GEOMETRY [REPORT]]")
    include(joinpath(@__DIR__,"diagnose_molecular_exterior.jl"))
    const M=DiagnoseMolecularExterior
    root,prepared=ARGS[1:2]; g=M.saved_geometry(root,prepared)
    @testset "Independent metadata geometry and all saved memberships" begin
        oracle=Dict{String,Any}()
        hullrows=M.tsv(joinpath(prepared,"hulls.tsv"))
        for m in M.R.MOLECULES
            a=M.V.independent_geometry(g.src.run,m); atoms=filter(x->x.element!="Cu",a.atoms)
            poly=oracle_hull([x.position[1:2] for x in atoms]); radius=maximum(a.radii[x.element] for x in atoms)
            oracle[m]=(;poly,radius,cell=a.cell)
            hr=filter(r->r["molecule"]==m,hullrows)
            @test [(Q(M.n(r,"x_nm")),Q(M.n(r,"y_nm"))) for r in hr]==poly
            @test all(M.n(r,"radius_nm")==radius for r in hr)
            @test all(E.distance2(E.point(x.position[1:2]),poly)==0 for x in atoms)
        end
        old=M.tsv(joinpath(g.src.run,"grid.tsv"))
        @test length(g.rows)==length(old)==g.meta["sites"]
        @test M.i.(g.rows,"pixel")==M.i.(old,"pixel")
        for (r,o) in zip(g.rows,old)
            xy=(M.n(r,"x_nm"),M.n(r,"y_nm"))
            @test xy==(M.n(o,"x_nm"),M.n(o,"y_nm"))
            for m in M.R.MOLECULES
                f=oracle[m]; sign2=oracle_clearance(xy,f.poly,f.radius,f.cell)
                @test M.b(r,m*"_exterior")== (sign2>0)
                @test M.b(r,m*"_boundary")== (sign2==0)
                @test M.n(r,m*"_margin_nm")≈sqrt(Float64(sign2+Q(f.radius)^2))-f.radius atol=1e-15
            end
            @test M.b(r,"prior_eligible")==M.b(o,"eligible")
            @test M.b(r,"common_exterior")== (M.b(r,"glcn_exterior") && M.b(r,"glcnac_exterior"))
            @test M.b(r,"removed_by_envelope")== (M.b(o,"eligible") && !M.b(r,"common_exterior"))
            @test !M.b(r,"common_exterior") || M.b(o,"eligible")
        end
        for (key,column) in (("prior","prior_eligible"),("exterior","common_exterior"),("removed","removed_by_envelope"))
            @test count(r->M.b(r,column),g.rows)==g.meta[key]
        end
        @test g.meta["prior"]==g.meta["exterior"]+g.meta["removed"]
    end
    if length(ARGS)==3
        report=ARGS[3]; states=M.tsv(joinpath(report,"states.tsv")); paired=M.tsv(joinpath(report,"paired.tsv"))
        summary=M.tsv(joinpath(report,"paired_summary.tsv")); groups=g.s["selection"]["groups"]
        domains=g.s["selection"]["domains"]
        @testset "Saved-height grouping, complements and paired arithmetic" begin
            @test length(states)==2length(groups)*length(domains)*3
            @test length(summary)==length(groups)*length(domains)*3
            roots=Dict(m=>M.tsv(joinpath(g.src.run,m,"diagnostic/roots.tsv")) for m in M.R.MOLECULES)
            central=Dict(m=>M.tsv(joinpath(g.src.run,m,"molecular_surfaces.tsv")) for m in M.R.MOLECULES)
            ids=Dict(group=>Set(M.i(r,"pixel") for r in g.rows if M.b(r,group=="all_prior" ? "prior_eligible" : group)) for group in groups)
            lookup=Dict((m,r["domain"],M.i(r,"interval"),M.i(r,"pixel"))=>r for m in M.R.MOLECULES for r in roots[m])
            for r in paired
                d=r["domain"]; j=M.i(r,"interval"); p=M.i(r,"pixel"); group=r["group"]
                @test p in ids[group]
                a=lookup[("glcn",d,j,p)]; c=lookup[("glcnac",d,j,p)]
                @test M.b(a,"valid") && M.b(c,"valid")
                original=mean(M.n(x,"z_nm") for x in central["glcnac"] if M.i(x,"interval")==j)-
                    mean(M.n(x,"z_nm") for x in central["glcn"] if M.i(x,"interval")==j)
                @test M.n(r,"original_molecular_contrast_nm")≈original atol=1e-15
                @test M.n(r,"background_contrast_nm")==M.n(c,"z_nm")-M.n(a,"z_nm")
                @test M.n(r,"referenced_contrast_nm")==M.n(r,"original_molecular_contrast_nm")-M.n(r,"background_contrast_nm")
            end
            for s in states
                rr=filter(r->r["domain"]==s["domain"] && r["interval"]==s["interval"] && M.i(r,"pixel") in ids[s["group"]],roots[s["molecule"]])
                @test M.i(s,"eligible")==length(rr)
                @test M.i(s,"count")==count(r->M.b(r,"valid"),rr)
                @test M.i(s,"absent")==count(r->M.i(r,"root_count")==0 && M.i(r,"unresolved")==0,rr)
                @test M.i(s,"multiple")==count(r->M.i(r,"root_count")>1,rr)
                @test M.i(s,"unresolved")==count(r->M.i(r,"unresolved")>0,rr)
            end
            for s in summary
                rr=filter(r->all(r[k]==s[k] for k in ("group","domain","interval")),paired)
                @test M.i(s,"count")==length(rr)
                @test M.i(s,"eligible")==length(ids[s["group"]])
                expected=Set(p for p in ids[s["group"]] if all(M.b(lookup[(m,s["domain"],M.i(s,"interval"),p)],"valid") for m in M.R.MOLECULES))
                @test Set(M.i.(rr,"pixel"))==expected && allunique(M.i.(rr,"pixel"))
                values=BigFloat.(M.n.(rr,"referenced_contrast_nm"))
                for (key,fun) in (("mean_nm",mean),("median_nm",median),("min_nm",minimum),("max_nm",maximum),("std_nm",x->std(x;corrected=false)))
                    @test isempty(values) ? isnan(M.n(s,key)) : isapprox(M.n(s,key),Float64(fun(values));atol=1e-15,rtol=1e-14)
                end
                @test M.i(s,"negative")==count(<(0),values)
                @test M.i(s,"positive")==count(>(0),values)
                @test M.i(s,"zero")==count(==(0),values)
            end
            for d in domains, j in 1:3
                ss=filter(r->r["domain"]==d && M.i(r,"interval")==j,summary)
                @test length(ss)==3 && Set(r["group"] for r in ss)==Set(groups)
                for key in ("eligible","count","negative","positive","zero")
                    @test M.i(only(filter(r->r["group"]=="all_prior",ss)),key)==
                        sum(M.i(r,key) for r in ss if r["group"]!="all_prior")
                end
            end
        end
    end
end
