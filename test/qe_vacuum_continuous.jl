#!/usr/bin/env julia
# Enclosed continuous crossings at all three saved common representatives.
module QEVacuumContinuous
using LinearAlgebra, TOML
include(joinpath(@__DIR__,"qe_vacuum_crossings.jl"))
include(joinpath(@__DIR__,"lib/vacuum_fourier_intervals.jl"))
const Q=QEVacuumCrossings
const D=Q.D
const S=Q.S
const F=VacuumFourierIntervals
const IA=F.IA
const REFERENCE_SHA=Dict(
    "common_intervals.tsv"=>"2bfb5897b956c40c1455be21b4aba88677e9d1d44c23194de127921f6ea7f3c1",
    "surfaces.tsv"=>"c41b92a06a65b8ee4a5864fa5582be99515f72b759cc8c415084b643ce62142e",
    "glcn_profiles.tsv"=>"6dd23636f53c8d5b3477e6a144a93d373bbb6f1d4806e63e62b6ded0ec185e12",
    "glcnac_profiles.tsv"=>"a6105253cd43e9181f096312ad6a2d8d81ef6ee37fd504f9f38748f33e99112e")
number(r,key)=parse(Float64,r[key])

function settings(file)
    s=TOML.parsefile(file); old=Q.settings(joinpath(@__DIR__,"../config/qe_vacuum_crossings.toml"))
    expected=Dict("model"=>old["model"],"selection"=>Dict(
        "policy"=>"continuous_unique_descending_saved_representatives",
        "domain"=>"unchanged_molecular_half_native_endpoints",
        "require_reference_sha256"=>true,"require_exact_repeat"=>true,
        "repeat_scope"=>"interval_roots_and_direct_surface_values","require_outside_paw"=>true,
        "interval_arithmetic_version"=>"1.0.12","taylor_order"=>12,
        "root_width_nm"=>1e-7,"max_depth"=>24,"max_nodes"=>20000,"newton_steps"=>32,
        "numeric_roundoff_atol"=>1e-13,"geometry_rtol"=>1e-12),
        "preprocessing"=>Dict("coordinates"=>"accepted_xml_cell_substrate_normal",
            "harmonics"=>"all_stored_gamma_pairs","clip_negative_values"=>false,
            "normalize_components"=>false,"half_nm"=>.32,"step_nm"=>.04,"fourier_block_points"=>32))
    s==expected || error("Not the frozen continuous-crossing experiment")
    string(pkgversion(IA))==s["selection"]["interval_arithmetic_version"] || error("Changed interval library")
    s
end

function prepare(root,config,outdir)
    s=settings(config); previous=joinpath(root,"qe/vacuum_crossings_20260925")
    Q.settings(joinpath(previous,"settings.toml"))
    all(values(TOML.parsefile(joinpath(previous,"analysis/summary.toml"))["checks"])) || error("Previous comparison failed")
    for (name,sha) in REFERENCE_SHA
        S.sha(joinpath(previous,"analysis",name))==sha || error("Changed saved reference")
    end
    for (k,m) in enumerate(Q.MOLECULES)
        S.sha(joinpath(previous,m,"data-file-schema.xml"))==D.G.XML_SHA[k] || error("Changed accepted state")
        Q.geometry(previous,m,s)
    end
    S.newdir(outdir); cp(config,joinpath(outdir,"settings.toml"))
    cp(joinpath(previous,"pseudo"),joinpath(outdir,"pseudo"))
    cp(joinpath(previous,"geometry_before_density.tsv"),joinpath(outdir,"geometry.tsv"))
    for (name,_) in REFERENCE_SHA; cp(joinpath(previous,"analysis",name),joinpath(outdir,name)); end
    cp(joinpath(root,"hpc/qe_vacuum_continuous.sbatch"),joinpath(outdir,"run_continuous.sbatch"))
    for (k,m) in enumerate(Q.MOLECULES)
        dir=joinpath(outdir,m); mkdir(dir)
        for name in ("data-file-schema.xml","frame.tsv"); cp(joinpath(previous,m,name),joinpath(dir,name)); end
        _,rs=D.radii(outdir,dir)
        meta=Dict("config_sha256"=>S.sha(config),"xml_sha256"=>D.G.XML_SHA[k],
            "wfc_sha256"=>D.WFC_SHA[k],"frame_sha256"=>S.sha(joinpath(dir,"frame.tsv")),
            "pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in rs),"references"=>REFERENCE_SHA)
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    end
    println("Prepared all three saved representatives and unchanged molecular-half endpoints; no new density evaluation.")
end

function proof(c,lo,hi,iso,q)
    F.crossings(c,lo,hi,iso;root_width=q["root_width_nm"],max_depth=q["max_depth"],
        max_nodes=q["max_nodes"],newton_steps=q["newton_steps"])
end

function analyze(run,m,wfc_file)
    m in Q.MOLECULES || error("Unknown molecule")
    s=settings(joinpath(run,"settings.toml")); q=s["selection"]; p=s["preprocessing"]
    dir=joinpath(run,m); out=joinpath(dir,"continuous"); S.newdir(out)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed config")
    for (name,sha) in meta["references"]; S.sha(joinpath(run,name))==sha || error("Changed reference $name"); end
    xml=joinpath(dir,"data-file-schema.xml")
    S.sha(xml)==meta["xml_sha256"] && S.sha(wfc_file)==meta["wfc_sha256"] || error("Changed electronic state")
    S.sha(joinpath(dir,"frame.tsv"))==meta["frame_sha256"] || error("Changed frame")
    g=Q.geometry(run,m,s)
    Dict(r.pseudo=>r.sha256 for r in g.rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
    sp=S.spectrum(xml,s); selected=[r.band for r in sp.rows if r.ildos_weight>0]
    weights=[r.ildos_weight for r in sp.rows if r.ildos_weight>0]
    w=D.W.read_wfc(wfc_file,selected)
    length(sp.rows)==w.nbnd && isapprox(w.cell,g.geo.cell;rtol=q["geometry_rtol"]) || error("Changed WFC geometry")
    refs=D.V.tsv(joinpath(run,m*"_profiles.tsv"))
    saved=D.V.tsv(joinpath(run,"surfaces.tsv"))
    representatives=filter(r->r["source"]=="joint" && parse(Int,r["interval"])>0,D.V.tsv(joinpath(run,"common_intervals.tsv")))
    all(r->r["domain"]=="molecular_half_paw_gap",representatives) || error("Changed saved domains")
    isos=number.(representatives,"isovalue"); pixels=size(g.xy,1); nz=length(g.domain.knots)
    length(refs)==nz*pixels || error("Changed profile size")
    lo=first(g.domain.knots); hi=g.domain.knots[last(g.domain.half)]
    series=Vector{Any}(undef,pixels); results=Vector{Any}(undef,pixels)
    L=F.I(g.geo.cell[3,3])*F.I(D.G.C.BOHR_NM)
    BLAS.set_num_threads(1)
    Threads.@threads for pixel in 1:pixels
        xy=[F.I(g.xy[pixel,k])/F.I(D.G.C.BOHR_NM)/F.I(g.geo.cell[k,k]) for k in 1:2]
        c=F.column(w,xy,L,weights,q["taylor_order"]); series[pixel]=c
        results[pixel]=[proof(c,lo,hi,iso,q) for iso in isos]
    end
    println(m,": all ",pixels," columns regrouped and continuous searches completed"); flush(stdout)
    repeat_exact=true
    for pixel in 1:pixels, j in eachindex(isos)
        repeat_exact &= isequal(results[pixel][j],proof(series[pixel],lo,hi,isos[j],q))
    end
    println(m,": serial continuous-search repeat exact=",repeat_exact); flush(stdout)
    coeff_rows=NamedTuple[]; profile_rows=NamedTuple[]; leaves=NamedTuple[]; roots=NamedTuple[]
    for pixel in 1:pixels
        c=series[pixel]
        for k in axes(c.a,1), band in axes(c.a,2)
            push!(coeff_rows,(;pixel,band=selected[band],harmonic=k-1,
                a_lo=IA.inf(c.a[k,band]),a_hi=IA.sup(c.a[k,band]),
                b_lo=IA.inf(c.b[k,band]),b_hi=IA.sup(c.b[k,band]),
                omega_lo=IA.inf(c.omega[k]),omega_hi=IA.sup(c.omega[k])))
        end
        for k in 1:nz
            ref=refs[(pixel-1)*nz+k]; z=g.domain.knots[k]
            parse(Int,ref["pixel"])==pixel && number(ref,"z_nm")==z || error("Changed saved coordinates")
            all(number(ref,axis*"_nm")==g.xy[pixel,j] for (j,axis) in enumerate(("x","y"))) || error("Changed lateral frame")
            value=F.density(c,F.I(z)); direct=number(ref,"direct")
            distance=max(0.,IA.inf(value)-direct,direct-IA.sup(value))
            push!(profile_rows,(;pixel,z_nm=z,direct,lower=IA.inf(value),upper=IA.sup(value),distance))
        end
        for (j,iso) in enumerate(isos)
            r=results[pixel][j]; interval=parse(Int,representatives[j]["interval"])
            for leaf in r.leaves; push!(leaves,merge((;pixel,interval,isovalue=iso),leaf)); end
            old=only(filter(x->x["molecule"]==m && parse(Int,x["pixel"])==pixel && parse(Int,x["interval"])==interval,saved))
            good=r.valid_root
            zl=good ? only(r.roots).root_lo : NaN; zh=good ? only(r.roots).root_hi : NaN
            z=good ? IA.mid(F.I(zl,zh)) : NaN
            value=good ? F.density(c,F.I(z)) : F.I(0)
            push!(roots,(;pixel,interval,isovalue=iso,valid=good,unresolved=r.unresolved,
                root_count=length(r.roots),nodes=r.nodes,x_nm=g.xy[pixel,1],y_nm=g.xy[pixel,2],
                z_lower_nm=zl,z_upper_nm=zh,z_nm=z,old_z_nm=number(old,"z_nm"),
                density_lower=good ? IA.inf(value) : NaN,density_upper=good ? IA.sup(value) : NaN))
        end
    end
    Q.table(joinpath(out,"coefficients.tsv"),keys(first(coeff_rows)),coeff_rows)
    Q.table(joinpath(out,"profiles.tsv"),keys(first(profile_rows)),profile_rows)
    Q.table(joinpath(out,"leaves.tsv"),keys(first(leaves)),leaves)
    Q.table(joinpath(out,"roots.tsv"),keys(first(roots)),roots)
    good=findall(r->r.valid,roots); surface_rows=NamedTuple[]; surface_repeat=true
    if !isempty(good)
        points=hcat([[roots[i].x_nm,roots[i].y_nm,roots[i].z_nm] for i in good]...)
        fractions=g.geo.cell\(points/D.G.C.BOHR_NM)
        dv=D.W.density(w,fractions,weights;block_points=p["fourier_block_points"])
        rep=D.W.density(w,fractions,weights;block_points=p["fourier_block_points"],parallel=false)
        surface_repeat=dv==rep
        for (j,i) in enumerate(good)
            r=roots[i]; distance=max(0.,r.density_lower-dv[j],dv[j]-r.density_upper)
            push!(surface_rows,merge(r,(;direct=dv[j],repeat=rep[j],distance,
                relative_density_error=(dv[j]-r.isovalue)/r.isovalue)))
        end
    end
    surface_keys=(keys(first(roots))...,:direct,:repeat,:distance,:relative_density_error)
    Q.table(joinpath(out,"surfaces.tsv"),surface_keys,surface_rows)
    checks=Dict("root_repeat_exact"=>repeat_exact,"surface_repeat_exact"=>surface_repeat,
        "grouped_profiles_agree"=>all(r->r.distance<=q["numeric_roundoff_atol"],profile_rows),
        "direct_surfaces_agree"=>all(r->r.distance<=q["numeric_roundoff_atol"],surface_rows))
    summary=Dict("scope"=>"continuous_finite_basis_roots_not_experimental_calibration",
        "checks"=>checks,"certified_complete_maps"=>all(values(checks)) && all(r->r.valid,roots),
        "valid_roots"=>length(good),"targets"=>length(roots),"unresolved_leaves"=>sum(r.unresolved for r in roots),
        "max_profile_distance"=>maximum(r.distance for r in profile_rows),
        "max_abs_relative_density_error"=>isempty(surface_rows) ? NaN : maximum(abs(r.relative_density_error) for r in surface_rows),
        "config_sha256"=>S.sha(joinpath(run,"settings.toml")),"wfc_sha256"=>meta["wfc_sha256"],
        "volume_bohr3"=>w.volume,"weights"=>weights,"selected_bands"=>selected,
        "interval_arithmetic_version"=>string(pkgversion(IA)),"domain_lower_nm"=>lo,"domain_upper_nm"=>hi)
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(m,": ",summary); flush(stdout)
    all(values(checks)) || error("Saved numerical comparison failed; no retry or bound adjustment")
    nothing
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_vacuum_continuous.jl prepare ROOT CONFIG NEW_DIR | analyze RUN MOLECULE WFC")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1]=="analyze" && return analyze(args[2:end]...)
    error("Expected prepare or analyze; no label, isovalue or domain override")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEVacuumContinuous.main()
