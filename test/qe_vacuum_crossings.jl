#!/usr/bin/env julia
# Substrate-normal surfaces from geometry-defined vacuum domains; no fit/labels.
module QEVacuumCrossings
using LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__,"qe_vacuum_wavefunctions.jl"))
include(joinpath(@__DIR__,"lib/vacuum_crossings.jl"))
const D=QEVacuumWavefunctions
const S=D.S
const X=VacuumCrossings
const MOLECULES=("glcn","glcnac")

function settings(file)
    s=TOML.parsefile(file)
    base=D.settings(joinpath(@__DIR__,"../config/qe_vacuum_wavefunctions.toml"))
    expected=Dict("model"=>base["model"],
        "selection"=>Dict("policy"=>"all_common_unique_descending_intervals",
            "domains"=>["full_paw_gap","molecular_half_paw_gap"],
            "require_reference_sha256"=>true,"require_exact_repeat"=>true,
            "require_outside_paw"=>true,"geometry_rtol"=>1e-12,
            "surface_representative"=>"geometric_midpoint_each_open_interval"),
        "preprocessing"=>Dict("cube_order"=>"qe_last_axis_fast",
            "coordinates"=>"accepted_xml_cell_substrate_normal",
            "vertical_knots"=>"native_grid_strictly_inside_paw_planes",
            "clip_negative_values"=>false,"normalize_components"=>false,
            "half_nm"=>.32,"step_nm"=>.04,"fourier_block_points"=>32))
    s==expected || error("Not the frozen vacuum-crossing diagnostic")
    s
end

function table(file,keys,rows)
    open(file,"w") do io
        println(io,join(keys,'\t'))
        for r in rows; println(io,join([getproperty(r,k) for k in keys],'\t')); end
    end
end

function geometry(run,m,s)
    dir=joinpath(run,m); geo=D.xml_geometry(joinpath(dir,"data-file-schema.xml"))
    radii,rs=D.radii(run,dir); frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
    domain=X.domains(geo.geometry,radii,geo.dims[3]); p=s["preprocessing"]
    xy=X.lateral_points(frame,p["half_nm"],p["step_nm"])
    # The global PAW planes enclose every sphere and every interpolation vertex.
    # Check the lateral cube domain now, before any density payload is read.
    q=xy[:,1:2]./transpose(geo.geometry.cell_nm[1:2]).*transpose(geo.dims[1:2])
    all(0 .<= floor.(Int,q)) && all(floor.(Int,q).+1 .< transpose(geo.dims[1:2])) || error("Lateral domain")
    (;geo,radii,rs,frame,domain,xy)
end

function prepare(root,config,outdir)
    s=settings(config); previous=joinpath(root,"qe/vacuum_wavefunctions_20260925")
    rows=NamedTuple[]
    for (k,m) in enumerate(MOLECULES)
        dir=joinpath(previous,m)
        S.sha(joinpath(dir,"data-file-schema.xml"))==D.G.XML_SHA[k] || error("Changed XML")
        S.sha(joinpath(dir,"reference_smooth.cube"))==D.SMOOTH_SHA[k] || error("Changed cube")
        all(values(TOML.parsefile(joinpath(dir,"comparison/summary.toml"))["checks"])) || error("WFC verification incomplete")
        g=geometry(previous,m,s); d=g.domain
        push!(rows,(;molecule=m,lower_paw_nm=d.lower,midpoint_nm=d.midpoint,upper_paw_nm=d.upper,
            first_knot_nm=first(d.knots),half_last_knot_nm=d.knots[last(d.half)],
            last_knot_nm=last(d.knots),full_knots=length(d.knots),half_knots=length(d.half),
            ring_z_nm=g.frame.origin_nm[3]))
    end
    S.newdir(outdir); cp(config,joinpath(outdir,"settings.toml"))
    cp(joinpath(previous,"pseudo"),joinpath(outdir,"pseudo"))
    cp(joinpath(root,"hpc/qe_vacuum_crossings.sbatch"),joinpath(outdir,"run_crossings.sbatch"))
    table(joinpath(outdir,"geometry_before_density.tsv"),keys(first(rows)),rows)
    for (k,m) in enumerate(MOLECULES)
        dir=joinpath(outdir,m); mkdir(dir)
        for name in ("data-file-schema.xml","frame.tsv","reference_smooth.cube")
            cp(joinpath(previous,m,name),joinpath(dir,name))
        end
        _,rs=D.radii(outdir,dir)
        meta=Dict("config_sha256"=>S.sha(config),"xml_sha256"=>D.G.XML_SHA[k],
            "smooth_sha256"=>D.SMOOTH_SHA[k],"wfc_sha256"=>D.WFC_SHA[k],
            "frame_sha256"=>S.sha(joinpath(dir,"frame.tsv")),
            "pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in rs))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    end
    foreach(println,rows)
end

function load_case(run,m,s,wfc)
    dir=joinpath(run,m); meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
    for (name,key) in (("data-file-schema.xml","xml_sha256"),("reference_smooth.cube","smooth_sha256"),
            ("frame.tsv","frame_sha256"))
        S.sha(joinpath(dir,name))==meta[key] || error("Changed $name")
    end
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed config")
    S.sha(wfc)==meta["wfc_sha256"] || error("Changed WFC")
    g=geometry(run,m,s)
    Dict(r.pseudo=>r.sha256 for r in g.rs)==meta["pseudo_sha256"] || error("Changed UPF")
    sp=S.spectrum(joinpath(dir,"data-file-schema.xml"),s)
    selected=[r.band for r in sp.rows if r.ildos_weight>0]
    weights=[r.ildos_weight for r in sp.rows if r.ildos_weight>0]
    w=D.W.read_wfc(wfc,selected)
    length(sp.rows)==w.nbnd && isapprox(w.cell,g.geo.cell;rtol=s["selection"]["geometry_rtol"]) || error("WFC geometry")
    c=D.V.cube(joinpath(dir,"reference_smooth.cube"))
    c.dims==g.geo.dims && all(iszero,c.origin) || error("Changed cube grid")
    z=g.domain.knots; points=hcat([[p[1],p[2],h] for p in eachrow(g.xy) for h in z]...)
    fractions=g.geo.cell\(points/D.G.C.BOHR_NM)
    block=s["preprocessing"]["fourier_block_points"]
    dv=D.W.density(w,fractions,weights;block_points=block)
    repeat=D.W.density(w,fractions,weights;block_points=block,parallel=false)
    direct=reshape(dv,length(z),size(g.xy,1))
    interpolated=X.profiles(c,g.geo.geometry.cell_nm,g.xy,z)
    all(isfinite,dv) && all(>=(0),dv) && all(>=(0),interpolated) || error("Invalid vacuum density")
    (;g,w,weights,direct,interpolated,repeat_exact=dv==repeat,
        relative_l2=norm(direct-interpolated)/norm(direct),wfc_sha256=meta["wfc_sha256"])
end

function analyze(run,wfc_n,wfc_a)
    s=settings(joinpath(run,"settings.toml")); out=joinpath(run,"analysis"); S.newdir(out)
    BLAS.set_num_threads(1)
    cases=Dict{String,Any}(); common=Dict{Tuple{String,String},Vector{Tuple{Float64,Float64}}}()
    interval_rows=NamedTuple[]; diagnostics=NamedTuple[]; checks=Dict{String,Bool}()
    for (m,wfc) in zip(MOLECULES,(wfc_n,wfc_a))
        c=load_case(run,m,s,wfc); cases[m]=c; d=c.g.domain; z=d.knots
        checks[m*"_repeat_exact"]=c.repeat_exact
        rows=((;pixel=p,z_nm=z[k],x_nm=c.g.xy[p,1],y_nm=c.g.xy[p,2],
            cube=c.interpolated[k,p],direct=c.direct[k,p]) for p in axes(c.direct,2) for k in eachindex(z))
        table(joinpath(out,m*"_profiles.tsv"),(:pixel,:z_nm,:x_nm,:y_nm,:cube,:direct),rows)
        for domain in s["selection"]["domains"], source in ("cube","direct")
            inds=domain=="full_paw_gap" ? eachindex(z) : d.half
            vals=source=="cube" ? c.interpolated : c.direct
            key=(domain,source); joint=get(common,key,[(0.,Inf)])
            for p in axes(vals,2)
                v=vals[inds,p]; ranges=X.intervals(z[inds],v)
                joint=X.intersect_intervals(joint,ranges)
                push!(diagnostics,(;molecule=m,domain,source,pixel=p,
                    rising_segments=count(>(0),diff(v)),flat_segments=count(iszero,diff(v)),
                    intervals=length(ranges),minimum=minimum(v),maximum=maximum(v)))
                for (i,(lo,hi)) in enumerate(ranges)
                    push!(interval_rows,(;molecule=m,domain,source,pixel=p,interval=i,lo,hi))
                end
            end
            common[key]=joint
        end
        println(m,": ",size(c.direct)," profiles; direct repeat=",c.repeat_exact,
            "; full-gap interpolation relative L2=",c.relative_l2)
        flush(stdout)
    end
    table(joinpath(out,"column_intervals.tsv"),(:molecule,:domain,:source,:pixel,:interval,:lo,:hi),interval_rows)
    table(joinpath(out,"column_diagnostics.tsv"),keys(first(diagnostics)),diagnostics)
    common_rows=NamedTuple[]; surface_rows=NamedTuple[]
    for domain in s["selection"]["domains"]
        common[(domain,"joint")]=X.intersect_intervals(common[(domain,"cube")],common[(domain,"direct")])
        for source in ("cube","direct","joint")
            ranges=common[(domain,source)]
            isempty(ranges) && push!(common_rows,(;domain,source,interval=0,lo=NaN,hi=NaN,isovalue=NaN))
            for (i,(lo,hi)) in enumerate(ranges)
                iso=lo>0 ? exp((log(lo)+log(hi))/2) : hi/2
                lo<iso<hi || error("No representable interior isovalue")
                push!(common_rows,(;domain,source,interval=i,lo,hi,isovalue=iso))
                source=="joint" || continue
                for m in MOLECULES
                    c=cases[m]; z=c.g.domain.knots
                    inds=domain=="full_paw_gap" ? eachindex(z) : c.g.domain.half
                    roots=[X.crossing(z[inds],c.direct[inds,p],iso) for p in axes(c.direct,2)]
                    all(r->r.status==:found,roots) || error("Interval/crossing inconsistency")
                    points=hcat([[c.g.xy[p,1],c.g.xy[p,2],r.z] for (p,r) in enumerate(roots)]...)
                    fractions=c.g.geo.cell\(points/D.G.C.BOHR_NM)
                    block=s["preprocessing"]["fourier_block_points"]
                    dv=D.W.density(c.w,fractions,c.weights;block_points=block)
                    rep=D.W.density(c.w,fractions,c.weights;block_points=block,parallel=false)
                    checks[m*"_"*domain*"_surface_"*string(i)*"_repeat"]=dv==rep
                    for (p,r) in enumerate(roots)
                        rc=X.crossing(z[inds],c.interpolated[inds,p],iso)
                        rc.status==:found || error("Cube intersection inconsistency")
                        push!(surface_rows,(;domain,interval=i,molecule=m,pixel=p,isovalue=iso,
                            x_nm=points[1,p],y_nm=points[2,p],z_nm=r.z,
                            ring_relative_vertical_nm=r.z-c.g.frame.origin_nm[3],cube_z_nm=rc.z,
                            direct_at_surface=dv[p],relative_density_error=(dv[p]-iso)/iso))
                    end
                end
            end
        end
    end
    table(joinpath(out,"common_intervals.tsv"),(:domain,:source,:interval,:lo,:hi,:isovalue),common_rows)
    table(joinpath(out,"surfaces.tsv"),(:domain,:interval,:molecule,:pixel,:isovalue,:x_nm,:y_nm,
        :z_nm,:ring_relative_vertical_nm,:cube_z_nm,:direct_at_surface,:relative_density_error),surface_rows)
    result=Dict("scope"=>"geometry_defined_sampled_crossings_not_current_calibration",
        "checks"=>checks,"config_sha256"=>S.sha(joinpath(run,"settings.toml")),
        "surface_rows"=>length(surface_rows),"direct_between_knots_uniqueness_proven"=>false,
        "molecules"=>Dict(m=>Dict("profile_points"=>length(cases[m].direct),
            "interpolation_relative_l2"=>cases[m].relative_l2,"wfc_sha256"=>cases[m].wfc_sha256)
            for m in MOLECULES))
    open(io->TOML.print(io,result),joinpath(out,"summary.toml"),"w")
    foreach(println,common_rows)
    all(values(checks)) || error("Saved diagnostics fail exact repetition; no retry or adjustment")
    println("VACUUM_CROSSINGS_COMPLETE; no experimental current, mask, calibration or grade")
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_vacuum_crossings.jl prepare ROOT CONFIG NEW_DIR | analyze RUN WFC_GLCN WFC_GLCNAC")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1]=="analyze" && return analyze(args[2:end]...)
    error("Expected prepare or analyze; no labels, fit, current or cutoff options")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEVacuumCrossings.main()
