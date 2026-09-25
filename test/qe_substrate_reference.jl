#!/usr/bin/env julia
# Diagnostic off-projection references. No image, label or assignment input.
module QESubstrateReference
using LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__,"qe_vacuum_continuous.jl"))
include(joinpath(@__DIR__,"lib/substrate_reference_geometry.jl"))
const C=QEVacuumContinuous
const D=C.D
const S=C.S
const F=C.F
const IA=F.IA
const G=SubstrateReferenceGeometry
const V=D.V
const MOLECULES=C.Q.MOLECULES
n(r,k)=parse(Float64,r[k])
i(r,k)=parse(Int,r[k])
table(file,rows)=C.Q.table(file,keys(first(rows)),rows)

function settings(file)
    s=TOML.parsefile(file)
    old=C.settings(joinpath(@__DIR__,"../config/qe_vacuum_continuous.toml"))
    model=merge(old["model"],Dict("reference_run"=>"qe/vacuum_continuous_20260925"))
    expected=Dict("model"=>model,"selection"=>Dict(
        "policy"=>"report_all_saved_isovalues_and_domains_no_reference_adoption",
        "support"=>"outside_union_of_both_periodic_molecular_paw_projections",
        "domains"=>["full_cu_gap","cu_half_gap","saved_molecular_half"],
        "require_reference_sha256"=>true,"require_exact_repeat"=>true,
        "interval_arithmetic_version"=>"1.0.12","taylor_order"=>12,"root_width_nm"=>1e-7,
        "max_depth"=>24,"max_nodes"=>20000,"newton_steps"=>32,
        "numeric_roundoff_atol"=>1e-13,"geometry_rtol"=>1e-12),
        "preprocessing"=>Dict("coordinates"=>"accepted_xml_cell_substrate_normal",
            "lateral_grid"=>[32,24],"lateral_sampling"=>"fractional_cell_centers",
            "vertical_knots"=>"native_grid_strictly_inside_cu_paw_planes",
            "clip_negative_values"=>false,"normalize_components"=>false,"fourier_block_points"=>32))
    s==expected || error("Not the frozen substrate-reference diagnostic")
    s
end

function geometry(run,s)
    old=C.settings(joinpath(@__DIR__,"../config/qe_vacuum_continuous.toml"))
    cases=[C.Q.geometry(run,m,old) for m in MOLECULES]
    grid=G.grid([c.geo.geometry for c in cases],[c.radii for c in cases],s["preprocessing"]["lateral_grid"])
    domains=[G.domains(c.geo.geometry,c.radii,c.geo.dims[3],c.domain) for c in cases]
    margins=NamedTuple[]
    for (k,m) in enumerate(MOLECULES), p in grid, d in domains[k]
        p.eligible || continue
        margin=G.segment_margin([p.x_nm,p.y_nm],d.lo,d.hi,cases[k].geo.geometry,cases[k].radii)
        margin>0 || error("A retained vertical segment intersects a PAW sphere")
        push!(margins,(;molecule=m,pixel=p.pixel,domain=d.domain,clearance_nm=margin))
    end
    (;cases,grid,domains,margins)
end

function prepare(root,config,outdir)
    s=settings(config); previous=joinpath(root,s["model"]["reference_run"])
    C.settings(joinpath(previous,"settings.toml"))
    sha=C.REFERENCE_SHA["common_intervals.tsv"]
    S.sha(joinpath(previous,"common_intervals.tsv"))==sha || error("Changed isovalues")
    pinned=TOML.parsefile(joinpath(root,"config/vacuum_surface_transfer.toml"))["model"]
    for (k,m) in enumerate(MOLECULES)
        dir=joinpath(previous,m)
        S.sha(joinpath(dir,"data-file-schema.xml"))==D.G.XML_SHA[k] || error("Changed state")
        S.sha(joinpath(dir,"frame.tsv"))==pinned[m*"_frame_sha256"] || error("Changed frame")
        S.sha(joinpath(dir,"continuous/surfaces.tsv"))==pinned[m*"_surface_sha256"] || error("Changed molecular surface")
        result=TOML.parsefile(joinpath(dir,"continuous/summary.toml"))
        result["certified_complete_maps"] && all(values(result["checks"])) || error("Incomplete molecular maps")
    end
    g=geometry(previous,s)
    S.newdir(outdir); cp(config,joinpath(outdir,"settings.toml"))
    cp(joinpath(previous,"pseudo"),joinpath(outdir,"pseudo"))
    cp(joinpath(previous,"common_intervals.tsv"),joinpath(outdir,"common_intervals.tsv"))
    table(joinpath(outdir,"grid.tsv"),g.grid)
    table(joinpath(outdir,"domains.tsv"),[merge((;molecule=m),d) for (k,m) in enumerate(MOLECULES) for d in g.domains[k]])
    table(joinpath(outdir,"segment_clearance.tsv"),g.margins)
    cp(joinpath(root,"hpc/qe_substrate_reference.sbatch"),joinpath(outdir,"run_substrate.sbatch"))
    for (k,m) in enumerate(MOLECULES)
        dir=joinpath(outdir,m); mkdir(dir); src=joinpath(previous,m)
        for name in ("data-file-schema.xml","frame.tsv"); cp(joinpath(src,name),joinpath(dir,name)); end
        cp(joinpath(src,"continuous/surfaces.tsv"),joinpath(dir,"molecular_surfaces.tsv"))
        meta=Dict("config_sha256"=>S.sha(config),"xml_sha256"=>D.G.XML_SHA[k],
            "wfc_sha256"=>D.WFC_SHA[k],"frame_sha256"=>S.sha(joinpath(dir,"frame.tsv")),
            "molecular_surface_sha256"=>S.sha(joinpath(dir,"molecular_surfaces.tsv")),
            "pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in g.cases[k].rs),
            "geometry_sha256"=>Dict(name=>S.sha(joinpath(outdir,name)) for name in
                ("grid.tsv","domains.tsv","segment_clearance.tsv","common_intervals.tsv")))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    end
    println("Geometry only: ",count(p->p.eligible,g.grid)," / ",length(g.grid),
        " paired off-projection sites; minimum segment clearance=",minimum(p.clearance_nm for p in g.margins)," nm")
    foreach(println,[merge((;molecule=m),d) for (k,m) in enumerate(MOLECULES) for d in g.domains[k]])
    nothing
end

function analyze(run,m,wfc_file;witness=false)
    k=findfirst(==(m),MOLECULES); k===nothing && error("Unknown molecule")
    s=settings(joinpath(run,"settings.toml")); q=s["selection"]; dir=joinpath(run,m)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed settings")
    for (name,sha) in meta["geometry_sha256"]; S.sha(joinpath(run,name))==sha || error("Changed geometry/isovalues"); end
    for (name,key) in (("data-file-schema.xml","xml_sha256"),("frame.tsv","frame_sha256"),
            ("molecular_surfaces.tsv","molecular_surface_sha256"))
        S.sha(joinpath(dir,name))==meta[key] || error("Changed $name")
    end
    g=geometry(run,s); case=g.cases[k]; domains=g.domains[k]
    Dict(r.pseudo=>r.sha256 for r in case.rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
    sites=filter(p->p.eligible,g.grid)
    witness && (sites=sites[1:1]) # fixed first GEOMETRY-eligible site; never density-selected
    isempty(sites) && error("No paired off-projection support")
    representatives=filter(r->r["source"]=="joint" && i(r,"interval")>0,V.tsv(joinpath(run,"common_intervals.tsv")))
    all(r->r["domain"]=="molecular_half_paw_gap",representatives) || error("Changed representatives")
    isos=n.(representatives,"isovalue"); ids=i.(representatives,"interval")
    out=joinpath(dir,witness ? "witness" : "diagnostic"); S.newdir(out)
    S.sha(wfc_file)==meta["wfc_sha256"] || error("Changed WFC")
    sp=S.spectrum(joinpath(dir,"data-file-schema.xml"),s)
    selected=[r.band for r in sp.rows if r.ildos_weight>0]
    weights=[r.ildos_weight for r in sp.rows if r.ildos_weight>0]
    w=D.W.read_wfc(wfc_file,selected)
    length(sp.rows)==w.nbnd && isapprox(w.cell,case.geo.cell;rtol=q["geometry_rtol"]) || error("WFC geometry")
    L=F.I(case.geo.cell[3,3])*F.I(D.G.C.BOHR_NM)
    series=Vector{Any}(undef,length(sites)); results=Vector{Any}(undef,length(sites))
    BLAS.set_num_threads(1)
    Threads.@threads for p in eachindex(sites)
        xy=[F.I(sites[p].x_nm)/F.I(D.G.C.BOHR_NM)/F.I(case.geo.cell[1,1]),
            F.I(sites[p].y_nm)/F.I(D.G.C.BOHR_NM)/F.I(case.geo.cell[2,2])]
        c=F.column(w,xy,L,weights,q["taylor_order"]); series[p]=c
        results[p]=[C.proof(c,d.lo,d.hi,iso,q) for d in domains, iso in isos]
    end
    println(m,": ",length(sites)," columns, all three domains/isovalues searched"); flush(stdout)
    repeated=all(isequal(results[p][di,j],C.proof(series[p],d.lo,d.hi,iso,q))
        for p in eachindex(sites) for (di,d) in enumerate(domains) for (j,iso) in enumerate(isos))
    println(m,": search repeat exact=",repeated); flush(stdout)
    coefficients=NamedTuple[]; leaves=NamedTuple[]; roots=NamedTuple[]; probes=NamedTuple[]
    for (p,site) in enumerate(sites)
        c=series[p]
        for harmonic in axes(c.a,1), band in axes(c.a,2)
            push!(coefficients,(;pixel=site.pixel,band=selected[band],harmonic=harmonic-1,
                a_lo=IA.inf(c.a[harmonic,band]),a_hi=IA.sup(c.a[harmonic,band]),
                b_lo=IA.inf(c.b[harmonic,band]),b_hi=IA.sup(c.b[harmonic,band]),
                omega_lo=IA.inf(c.omega[harmonic]),omega_hi=IA.sup(c.omega[harmonic])))
        end
        for (di,d) in enumerate(domains)
            for (kind,z) in (("lower",d.lo),("midpoint",(d.lo+d.hi)/2),("upper",d.hi))
                value=F.density(c,F.I(z))
                push!(probes,(;pixel=site.pixel,domain=d.domain,interval=0,kind,leaf=0,
                    x_nm=site.x_nm,y_nm=site.y_nm,z_nm=z,lower=IA.inf(value),upper=IA.sup(value)))
            end
            for (j,iso) in enumerate(isos)
                r=results[p][di,j]; interval=ids[j]
                for (leafno,leaf) in enumerate(r.leaves)
                    push!(leaves,merge((;pixel=site.pixel,domain=d.domain,interval,isovalue=iso,leaf=leafno),leaf))
                    leaf.status==:unique || continue
                    z=IA.mid(F.I(leaf.root_lo,leaf.root_hi)); value=F.density(c,F.I(z))
                    push!(probes,(;pixel=site.pixel,domain=d.domain,interval,kind="root",leaf=leafno,
                        x_nm=site.x_nm,y_nm=site.y_nm,z_nm=z,lower=IA.inf(value),upper=IA.sup(value)))
                end
                good=r.valid_root
                push!(roots,(;pixel=site.pixel,domain=d.domain,interval,isovalue=iso,
                    valid=good,root_count=length(r.roots),descending=count(x->x.direction==-1,r.roots),
                    ascending=count(x->x.direction==1,r.roots),unresolved=r.unresolved,nodes=r.nodes,
                    z_lower_nm=good ? only(r.roots).root_lo : NaN,
                    z_upper_nm=good ? only(r.roots).root_hi : NaN,
                    z_nm=good ? IA.mid(F.I(only(r.roots).root_lo,only(r.roots).root_hi)) : NaN))
            end
        end
    end
    table(joinpath(out,"coefficients.tsv"),coefficients)
    table(joinpath(out,"leaves.tsv"),leaves); table(joinpath(out,"roots.tsv"),roots)
    points=hcat([[r.x_nm,r.y_nm,r.z_nm] for r in probes]...)
    fractions=case.geo.cell\(points/D.G.C.BOHR_NM)
    block=s["preprocessing"]["fourier_block_points"]
    direct=D.W.density(w,fractions,weights;block_points=block)
    rep=D.W.density(w,fractions,weights;block_points=block,parallel=false)
    checked=[merge(r,(;direct=direct[j],repeat=rep[j],
        distance=max(0.,r.lower-direct[j],direct[j]-r.upper))) for (j,r) in enumerate(probes)]
    table(joinpath(out,"direct_queries.tsv"),checked)
    checks=Dict("root_repeat_exact"=>repeated,"direct_repeat_exact"=>direct==rep,
        "direct_enclosure_agreement"=>all(r->r.distance<=q["numeric_roundoff_atol"],checked),
        "finite_positive_direct"=>all(x->isfinite(x) && x>0,direct))
    summary=Dict("scope"=>"off_projection_diagnostic_not_clean_copper_or_current_calibration",
        "witness"=>witness,"checks"=>checks,"config_sha256"=>meta["config_sha256"],
        "wfc_sha256"=>meta["wfc_sha256"],"selected_bands"=>selected,"weights"=>weights,
        "volume_bohr3"=>w.volume,"sites"=>length(sites),"grid_sites"=>length(g.grid),
        "eligible_sites"=>count(p->p.eligible,g.grid),"targets"=>length(roots),
        "valid_roots"=>count(r->r.valid,roots),"unresolved_leaves"=>sum(r.unresolved for r in roots),
        "max_direct_distance"=>maximum(r.distance for r in checked))
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(m,": ",summary); flush(stdout)
    all(values(checks)) || error("Numerical verification failed; no adjustment or retry")
    nothing
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_substrate_reference.jl prepare ROOT CONFIG NEW_DIR | analyze RUN MOLECULE WFC | witness RUN MOLECULE WFC")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1] in ("analyze","witness") && return analyze(args[2:end]...;witness=args[1]=="witness")
    error("Expected prepare, analyze or fixed witness; no density, domain or point override")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QESubstrateReference.main()
