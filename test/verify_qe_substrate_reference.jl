#!/usr/bin/env julia
# Read-only independent geometry, partition, query and reference arithmetic.
# Does not reconstruct wavefunctions or claim a third physical calculation.
module VerifyQESubstrateReference
using Test, LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__,"qe_substrate_reference.jl"))
const R=QESubstrateReference
const V=R.V
const F=R.F
const IA=F.IA
n(r,k)=parse(Float64,r[k])
i(r,k)=parse(Int,r[k])
b(r,k)=parse(Bool,r[k])

function independent_geometry(run,m)
    xml=read(joinpath(run,m,"data-file-schema.xml"),String)
    output=match(r"<output>(.*?)</output>"s,xml).captures[1]
    atomtext=match(r"<atomic_positions>(.*?)</atomic_positions>"s,output).captures[1]
    atoms=[(;element=match(r"name=\"([^\"]+)\"",r.captures[1]).captures[1],
        position=parse.(Float64,split(r.captures[2]))*R.D.G.C.BOHR_NM)
        for r in eachmatch(r"<atom\b([^>]*)>(.*?)</atom>"s,atomtext)]
    cell=[parse.(Float64,split(match(Regex("<a$k>(.*?)</a$k>","s"),output).captures[1]))[k]*R.D.G.C.BOHR_NM for k in 1:3]
    fft=match(r"<fft_grid\b([^>]*)>",output).captures[1]
    nz=parse(Int,match(r"nr3=\"(\d+)\"",fft).captures[1])
    radii=Dict{String,Float64}()
    input=match(r"<input>(.*?)</input>"s,xml).captures[1]
    files=[strip(r.captures[1]) for r in eachmatch(r"<pseudo_file>(.*?)</pseudo_file>"s,input)]
    allunique(files) && all(f->basename(f)==f,files) || error("Invalid active pseudo list")
    for name in files
        file=joinpath(run,"pseudo",name)
        text=read(file,String); head=match(r"<PP_HEADER\b([^>]*)>"s,text).captures[1]
        element=strip(match(r"element=\"([^\"]+)\"",head).captures[1])
        haskey(radii,element) && error("Duplicate active element")
        radial=R.S.number.(split(match(r"<PP_R\b[^>]*>(.*?)</PP_R>"s,text).captures[1]))
        inds=[parse(Int,r.captures[1]) for r in eachmatch(r"cutoff_radius_index=\"(\d+)\"",text)]
        push!(inds,parse(Int,match(r"cutoff_r_index=\"(\d+)\"",text).captures[1]))
        radii[element]=radial[maximum(inds)]*R.D.G.C.BOHR_NM
    end
    (;atoms,cell,nz,radii)
end

# Explicit neighboring-cell enumeration, independent of minimum-image code.
function projected(g,xy)
    minimum(hypot(xy[1]-a.position[1]-ix*g.cell[1],xy[2]-a.position[2]-iy*g.cell[2])-g.radii[a.element]
        for a in g.atoms if a.element!="Cu" for ix in -1:1 for iy in -1:1)
end
function segment(g,xy,lo,hi)
    minimum(sqrt((xy[1]-a.position[1]-ix*g.cell[1])^2+(xy[2]-a.position[2]-iy*g.cell[2])^2+
        max(lo-a.position[3]-iz*g.cell[3],a.position[3]+iz*g.cell[3]-hi,0.)^2)-g.radii[a.element]
        for a in g.atoms for ix in -1:1 for iy in -1:1 for iz in -1:1)
end

function stats(v)
    isempty(v) && return (;count=0,mean_nm=NaN,median_nm=NaN,min_nm=NaN,max_nm=NaN,std_nm=NaN)
    (;count=length(v),mean_nm=mean(v),median_nm=median(v),min_nm=minimum(v),max_nm=maximum(v),std_nm=std(v;corrected=false))
end

function verify(run,out;witness=false)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    R.S.newdir(out); s=R.settings(joinpath(run,"settings.toml")); q=s["selection"]
    common=filter(r->r["source"]=="joint" && i(r,"interval")>0,V.tsv(joinpath(run,"common_intervals.tsv")))
    isos=Dict(i(r,"interval")=>n(r,"isovalue") for r in common)
    grid=V.tsv(joinpath(run,"grid.tsv")); ds=V.tsv(joinpath(run,"domains.tsv"))
    margins=V.tsv(joinpath(run,"segment_clearance.tsv")); geometries=Dict(m=>independent_geometry(run,m) for m in R.MOLECULES)
    rows=Dict{String,Any}(); summaries=NamedTuple[]; paired=NamedTuple[]
    @testset "Off-projection geometry and saved root evidence" begin
        @test length(grid)==prod(s["preprocessing"]["lateral_grid"])
        @test i.(grid,"pixel")==collect(1:length(grid))
        @test length(ds)==6 && Set(r["domain"] for r in ds)==Set(q["domains"])
        nx,ny=s["preprocessing"]["lateral_grid"]
        for row in grid
            ix,iy=i(row,"ix"),i(row,"iy")
            @test i(row,"pixel")==1+ix+iy*nx && 0<=ix<nx && 0<=iy<ny
            @test n(row,"fx")==(ix+.5)/nx && n(row,"fy")==(iy+.5)/ny
            @test n(row,"x_nm")==n(row,"fx")*geometries["glcn"].cell[1]
            @test n(row,"y_nm")==n(row,"fy")*geometries["glcn"].cell[2]
            p=[n(row,"x_nm"),n(row,"y_nm")]
            for m in R.MOLECULES; @test n(row,m*"_margin_nm")≈projected(geometries[m],p) atol=1e-14; end
            @test b(row,"eligible")==(n(row,"glcn_margin_nm")>0 && n(row,"glcnac_margin_nm")>0)
        end
        eligible=i.(filter(r->b(r,"eligible"),grid),"pixel")
        targets=witness ? eligible[1:1] : eligible
        @test length(margins)==6length(eligible)
        @test Set((r["molecule"],i(r,"pixel"),r["domain"]) for r in margins)==
            Set((m,p,d) for m in R.MOLECULES for p in eligible for d in q["domains"])
        for m in R.MOLECULES
            g=geometries[m]; domains=filter(r->r["molecule"]==m,ds)
            @test Set(r["domain"] for r in domains)==Set(q["domains"])
            for d in domains
                atoms=d["domain"]=="saved_molecular_half" ? g.atoms : filter(a->a.element=="Cu",g.atoms)
                lo=maximum(a.position[3]+g.radii[a.element] for a in atoms)
                hi=minimum(a.position[3]+g.cell[3]-g.radii[a.element] for a in atoms)
                @test n(d,"lower_paw_nm")==lo && n(d,"upper_paw_nm")==hi
                @test n(d,"midpoint_nm")== (lo+hi)/2
                # Old molecular domain used k*(L/nz), Cu domain uses k*L/nz.
                knots=d["domain"]=="saved_molecular_half" ? [k*(g.cell[3]/g.nz) for k in 0:g.nz-1] :
                    [k*g.cell[3]/g.nz for k in 0:g.nz-1]
                keep=filter(z->lo<z<hi && (d["domain"]=="full_cu_gap" || z<(lo+hi)/2),knots)
                @test n(d,"lo")==first(keep) && n(d,"hi")==last(keep)
            end
            for r in filter(r->r["molecule"]==m,margins)
                p=grid[i(r,"pixel")]; d=only(filter(x->x["domain"]==r["domain"],domains))
                actual=segment(g,[n(p,"x_nm"),n(p,"y_nm")],n(d,"lo"),n(d,"hi"))
                @test n(r,"clearance_nm")≈actual atol=1e-14
                @test actual>0
            end
            dir=joinpath(run,m,witness ? "witness" : "diagnostic")
            summary=TOML.parsefile(joinpath(dir,"summary.toml")); meta=TOML.parsefile(joinpath(run,m,"metadata.toml"))
            @test summary["witness"]==witness && all(values(summary["checks"]))
            @test summary["config_sha256"]==meta["config_sha256"]==V.digest(joinpath(run,"settings.toml"))
            @test summary["wfc_sha256"]==meta["wfc_sha256"]==R.D.WFC_SHA[findfirst(==(m),R.MOLECULES)]
            for (name,key) in (("data-file-schema.xml","xml_sha256"),("frame.tsv","frame_sha256"),("molecular_surfaces.tsv","molecular_surface_sha256"))
                @test V.digest(joinpath(run,m,name))==meta[key]
            end
            for (name,sha) in meta["geometry_sha256"]; @test V.digest(joinpath(run,name))==sha; end
            for (name,sha) in meta["pseudo_sha256"]; @test V.digest(joinpath(run,"pseudo",name))==sha; end
            roots=V.tsv(joinpath(dir,"roots.tsv")); leaves=V.tsv(joinpath(dir,"leaves.tsv"))
            queries=V.tsv(joinpath(dir,"direct_queries.tsv")); rows[m]=roots
            @test length(roots)==9length(targets)==summary["targets"]
            @test summary["sites"]==length(targets) && summary["eligible_sites"]==length(eligible)
            @test Set((i(r,"pixel"),r["domain"],i(r,"interval")) for r in roots)==
                Set((p,d,j) for p in targets for d in q["domains"] for j in keys(isos))
            partitions=Dict{Tuple{Int,String,Int},Vector{typeof(first(leaves))}}()
            for r in leaves; push!(get!(partitions,(i(r,"pixel"),r["domain"],i(r,"interval")),typeof(first(leaves))[]),r); end
            @test Set(keys(partitions))==Set((i(r,"pixel"),r["domain"],i(r,"interval")) for r in roots)
            expected_queries=Set((p,d,0,kind,0) for p in targets for d in q["domains"] for kind in ("lower","midpoint","upper"))
            for r in roots
                key=(i(r,"pixel"),r["domain"],i(r,"interval")); ls=partitions[key]
                d=only(filter(x->x["domain"]==r["domain"],domains)); iso=n(r,"isovalue")
                @test iso==isos[i(r,"interval")]
                @test n(first(ls),"lo")==n(d,"lo") && n(last(ls),"hi")==n(d,"hi")
                @test all(n(ls[k],"hi")==n(ls[k+1],"lo") for k in 1:length(ls)-1)
                @test i.(ls,"leaf")==collect(1:length(ls))
                @test i(r,"nodes")==2length(ls)-1
                for leaf in ls
                    @test n(leaf,"lo")<n(leaf,"hi") && n(leaf,"isovalue")==iso
                    status=leaf["status"]
                    if status=="excluded"
                        @test iso<n(leaf,"rho_lo") || iso>n(leaf,"rho_hi")
                    elseif status in ("unique","excluded_monotone","precision_limit")
                        @test (i(leaf,"direction")==-1 && n(leaf,"derivative_hi")<0) ||
                            (i(leaf,"direction")==1 && n(leaf,"derivative_lo")>0)
                        if status!="excluded_monotone"
                            @test n(leaf,"lo")<=n(leaf,"root_lo")<=n(leaf,"root_hi")<=n(leaf,"hi")
                            status=="unique" && @test n(leaf,"root_hi")-n(leaf,"root_lo")<=q["root_width_nm"]
                        end
                    else
                        @test status in ("unresolved","node_limit")
                    end
                    status=="unique" && push!(expected_queries,(key...,"root",i(leaf,"leaf")))
                end
                unique_rows=filter(x->x["status"]=="unique",ls)
                unresolved=count(x->x["status"] in ("unresolved","node_limit","precision_limit"),ls)
                @test i(r,"root_count")==length(unique_rows) && i(r,"unresolved")==unresolved
                @test i(r,"descending")==count(x->i(x,"direction")==-1,unique_rows)
                @test i(r,"ascending")==count(x->i(x,"direction")==1,unique_rows)
                valid=unresolved==0 && length(unique_rows)==1 && i(only(unique_rows),"direction")==-1
                @test b(r,"valid")==valid
                if valid
                    @test n(r,"z_lower_nm")==n(only(unique_rows),"root_lo") && n(r,"z_upper_nm")==n(only(unique_rows),"root_hi")
                    @test n(r,"z_nm")==IA.mid(F.I(n(r,"z_lower_nm"),n(r,"z_upper_nm")))
                else
                    @test all(isnan(n(r,k)) for k in ("z_nm","z_lower_nm","z_upper_nm"))
                end
            end
            @test length(queries)==length(expected_queries)
            @test Set((i(r,"pixel"),r["domain"],i(r,"interval"),r["kind"],i(r,"leaf")) for r in queries)==expected_queries
            # Reconstruct saved enclosed coefficients, not the source WFC.
            coeffs=V.tsv(joinpath(dir,"coefficients.tsv")); bands=summary["selected_bands"]
            bypixel=Dict{Int,Vector{typeof(first(coeffs))}}()
            for c in coeffs; push!(get!(bypixel,i(c,"pixel"),typeof(first(coeffs))[]),c); end
            @test Set(keys(bypixel))==Set(targets)
            series=Dict{Int,Any}()
            for p in targets
                cr=bypixel[p]; nk=1+maximum(i(c,"harmonic") for c in cr)
                @test length(cr)==nk*length(bands) && Set(i(c,"band") for c in cr)==Set(bands)
                @test length(Set((i(c,"band"),i(c,"harmonic")) for c in cr))==length(cr)
                a=fill(F.I(0),nk,length(bands)); bb=copy(a); omega=fill(F.I(0),nk)
                for c in cr
                    h=i(c,"harmonic")+1; j=findfirst(==(i(c,"band")),bands)
                    a[h,j]=F.I(n(c,"a_lo"),n(c,"a_hi")); bb[h,j]=F.I(n(c,"b_lo"),n(c,"b_hi"))
                    omega[h]=F.I(n(c,"omega_lo"),n(c,"omega_hi"))
                end
                series[p]=F.make_series(a,bb,omega,F.I.(summary["weights"]),F.I(summary["volume_bohr3"]),q["taylor_order"])
            end
            for r in queries
                p=grid[i(r,"pixel")]; d=only(filter(x->x["domain"]==r["domain"],domains))
                @test n(r,"x_nm")==n(p,"x_nm") && n(r,"y_nm")==n(p,"y_nm")
                if r["kind"]=="root"
                    ls=partitions[(i(r,"pixel"),r["domain"],i(r,"interval"))]; l=ls[i(r,"leaf")]
                    @test l["status"]=="unique"
                    @test n(r,"z_nm")==IA.mid(F.I(n(l,"root_lo"),n(l,"root_hi")))
                else
                    z=r["kind"]=="lower" ? n(d,"lo") : r["kind"]=="upper" ? n(d,"hi") : (n(d,"lo")+n(d,"hi"))/2
                    @test n(r,"z_nm")==z && i(r,"interval")==i(r,"leaf")==0
                end
                val=F.density(series[i(r,"pixel")],F.I(n(r,"z_nm")))
                @test n(r,"lower")==IA.inf(val) && n(r,"upper")==IA.sup(val)
                @test n(r,"direct")==n(r,"repeat") && isfinite(n(r,"direct")) && n(r,"direct")>0
                @test n(r,"distance")==max(0.,n(r,"lower")-n(r,"direct"),n(r,"direct")-n(r,"upper"))<=q["numeric_roundoff_atol"]
            end
            @test summary["valid_roots"]==count(r->b(r,"valid"),roots)
            @test summary["unresolved_leaves"]==sum(i(r,"unresolved") for r in roots)
            @test summary["max_direct_distance"]==maximum(n(r,"distance") for r in queries)
            mold=V.tsv(joinpath(run,m,"molecular_surfaces.tsv"))
            for domain in q["domains"], j in sort(collect(keys(isos)))
                group=filter(r->r["domain"]==domain && i(r,"interval")==j,roots)
                heights=n.(filter(r->b(r,"valid"),group),"z_nm")
                ms=filter(r->i(r,"interval")==j,mold)
                @test length(ms)==289 && all(r->b(r,"valid"),ms)
                molecular_mean=mean(n.(ms,"z_nm")); molecular_range=extrema(n.(ms,"z_nm"))
                st=stats(heights)
                push!(summaries,merge((;molecule=m,domain,interval=j,isovalue=isos[j],eligible=length(group),
                    complete=length(heights)==length(group),absent=count(r->i(r,"root_count")==0 && i(r,"unresolved")==0,group),
                    multiple=count(r->i(r,"root_count")>1,group),unresolved=sum(i(r,"unresolved") for r in group),
                    one_descending=count(r->i(r,"descending")==1 && i(r,"unresolved")==0,group),
                    molecular_mean_nm=molecular_mean,molecular_min_nm=molecular_range[1],molecular_max_nm=molecular_range[2],
                    relative_mean_nm=molecular_mean-st.mean_nm),st))
            end
            println(m,": ",length(targets)," sites; ",length(roots)," targets; ",length(queries)," independent direct queries checked")
            flush(stdout)
        end
        for domain in q["domains"], j in sort(collect(keys(isos)))
            nrows=Dict(i(r,"pixel")=>r for r in rows["glcn"] if r["domain"]==domain && i(r,"interval")==j)
            arows=Dict(i(r,"pixel")=>r for r in rows["glcnac"] if r["domain"]==domain && i(r,"interval")==j)
            both=filter(p->b(nrows[p],"valid") && b(arows[p],"valid"),targets)
            sm=[only(filter(r->r.molecule==m && r.domain==domain && r.interval==j,summaries)) for m in R.MOLECULES]
            raw_contrast=sm[2].molecular_mean_nm-sm[1].molecular_mean_nm
            for p in both
                bg=n(arows[p],"z_nm")-n(nrows[p],"z_nm")
                push!(paired,(;domain,interval=j,pixel=p,background_contrast_nm=bg,
                    original_molecular_contrast_nm=raw_contrast,referenced_contrast_nm=raw_contrast-bg))
            end
            println((;domain,interval=j,paired=length(both),eligible=length(targets),
                original_molecular_contrast_nm=raw_contrast,
                background=stats([n(arows[p],"z_nm")-n(nrows[p],"z_nm") for p in both])))
        end
    end
    R.table(joinpath(out,"summary.tsv"),summaries)
    R.C.Q.table(joinpath(out,"paired.tsv"),(:domain,:interval,:pixel,:background_contrast_nm,
        :original_molecular_contrast_nm,:referenced_contrast_nm),paired)
    nothing
end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    ARGS==["--help"] ? println("verify_qe_substrate_reference.jl RUN NEW_REPORT [--witness]") :
        (length(ARGS) in (2,3) && (length(ARGS)==2 || ARGS[3]=="--witness") ?
            VerifyQESubstrateReference.verify(ARGS[1:2]...;witness=length(ARGS)==3) : error("Expected saved run and new report"))
end
