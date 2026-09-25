#!/usr/bin/env julia
# PP-only accepted-state potential diagnostic, never a chemical predictor.
module QEVacuumPotential
using TOML, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"qe_vacuum_wavefunctions.jl"))
include(joinpath(@__DIR__,"lib/qe_native_plot.jl"))
include(joinpath(@__DIR__,"lib/vacuum_potential.jl"))
const D=QEVacuumWavefunctions
const S=D.S
const V=D.V
const N=QENativePlot
const P=VacuumPotential
const MOLECULES=("glcn","glcnac")

function settings(file)
    s=TOML.parsefile(file)
    expected=Dict("model"=>Dict("qe_version"=>"7.4.1","source_run"=>"qe/substrate_reference_20260925",
        "sample_bias_ry"=>-0.0220495933,"spectral_window"=>"sharp_zero_temperature",
        "scf_acceptance_ry"=>5e-5,"total_plot_num"=>1,"electrostatic_plot_num"=>11,
        "glcn_density_sha256"=>"0636e18139cc6809be398eec5455153f3d32850437aa547f342c57a59ce299d1",
        "glcnac_density_sha256"=>"448de1eaea9760a222dcc994349684dcad38ff75202369285e12a1b15a0103ba"),
        "selection"=>Dict("policy"=>"all_native_planes_and_selected_bands_no_matching_plane_adoption",
            "whole_plane_domain"=>"strict_between_all_periodic_paw_bounding_planes",
            "report_lower_geometric_half"=>true,"require_exact_total_repeat"=>true,
            "adopt_reference"=>false,"propagate_wavefunctions"=>false),
        "preprocessing"=>Dict("primary_values"=>"native_filplot_first_axis_fast_with_padding",
            "independent_values"=>"cube_last_axis_fast","coordinates"=>"accepted_xml_cell_substrate_normal",
            "potential_units"=>"rydberg_energy","xml_energy_units"=>"hartree",
            "native_significant_digits"=>10,"cube_significant_digits"=>5,"native_alat_decimal_places"=>8,
            "native_tau_decimal_places"=>9,"cube_header_decimal_places"=>6,"geometry_rtol"=>1e-12,
            "numeric_roundoff_atol_ry"=>1e-13,"clip_potential"=>false,"normalize_potential"=>false))
    s==expected || error("Changed vacuum-potential experiment")
    s
end

function case_geometry(run,m)
    dir=joinpath(run,m); geo=D.xml_geometry(joinpath(dir,"data-file-schema.xml"))
    radii,rs=D.radii(run,dir); g=geo.geometry
    lower=maximum(a.position_nm[3]+radii[a.z] for a in g.atoms)
    upper=minimum(a.position_nm[3]+g.cell_nm[3]-radii[a.z] for a in g.atoms)
    lower<upper || error("No full-plane PAW-free gap")
    (;geo,radii,rs,lower,upper,midpoint=(lower+upper)/2)
end

function pp_input(m,stem,num)
    m in MOLECULES && stem in ("total","electrostatic","repeat_total") && num in (1,11) || error("Unsupported export")
    "&INPUTPP\n prefix='$(m)_central'\n outdir='../work/$m'\n filplot='$stem.dat'\n plot_num=$num\n/\n"*
        "&PLOT\n iflag=3\n output_format=6\n fileout='$stem.cube'\n/\n"
end

function prepare(root,config,out)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    s=settings(config); previous=joinpath(root,s["model"]["source_run"])
    cases=[case_geometry(previous,m) for m in MOLECULES]
    for (j,m) in enumerate(MOLECULES)
        meta=TOML.parsefile(joinpath(previous,m,"metadata.toml"))
        S.sha(joinpath(previous,m,"data-file-schema.xml"))==D.G.XML_SHA[j]==meta["xml_sha256"] || error("Changed accepted state")
        Dict(r.pseudo=>r.sha256 for r in cases[j].rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
        S.spectrum(joinpath(previous,m,"data-file-schema.xml"),s)
    end
    S.newdir(out); cp(config,joinpath(out,"settings.toml")); cp(joinpath(previous,"pseudo"),joinpath(out,"pseudo"))
    cp(joinpath(root,"hpc/qe_vacuum_potential.sbatch"),joinpath(out,"run_potential.sbatch"))
    open(joinpath(out,"accepted_states.sha256"),"w") do io
        for (j,m) in enumerate(MOLECULES)
            println(io,D.WFC_SHA[j],"  work/",m,"/",m,"_central.save/wfc1.dat")
            println(io,s["model"][m*"_density_sha256"],"  work/",m,"/",m,"_central.save/charge-density.dat")
        end
    end
    geometry=NamedTuple[]
    for (j,m) in enumerate(MOLECULES)
        dir=joinpath(out,m); mkdir(dir); cp(joinpath(previous,m,"data-file-schema.xml"),joinpath(dir,"data-file-schema.xml"))
        sp=S.spectrum(joinpath(dir,"data-file-schema.xml"),s); S.table(joinpath(dir,"bands.tsv"),sp.rows)
        for (stem,num) in (("total",s["model"]["total_plot_num"]),("electrostatic",s["model"]["electrostatic_plot_num"]),
                ("repeat_total",s["model"]["total_plot_num"]))
            write(joinpath(dir,stem*".in"),pp_input(m,stem,num))
        end
        c=cases[j]; nz=c.geo.dims[3]; L=c.geo.geometry.cell_nm[3]
        append!(geometry,[(;molecule=m,k,z_nm=k*L/nz,lower_paw_nm=c.lower,upper_paw_nm=c.upper,
            midpoint_nm=c.midpoint,paw_gap=c.lower<k*L/nz<c.upper,
            lower_half=c.lower<k*L/nz<c.midpoint) for k in 0:nz-1])
        meta=Dict("xml_sha256"=>D.G.XML_SHA[j],"wfc_sha256"=>D.WFC_SHA[j],"config_sha256"=>S.sha(config),
            "pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in c.rs),"fermi_ev"=>sp.fermi_ev,
            "scf_error_ry"=>sp.scf_error_ry,"bands_sha256"=>S.sha(joinpath(dir,"bands.tsv")))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    end
    S.table(joinpath(out,"geometry.tsv"),geometry)
    println("Prepared metadata only: ",count(r->r.paw_gap,geometry)," full-plane PAW-free native heights; no potential read")
end

function validate_grid(native,cube,g,p)
    native.dims==Tuple(g.geo.dims)==Tuple(cube.dims) || error("FFT dimensions differ")
    all(iszero,cube.origin) || error("Unexpected cube origin")
    cell=g.geo.cell; rt=p["geometry_rtol"]
    cell_bound=.5*10.0^(-p["native_alat_decimal_places"]).*abs.(native.at).+rt.*abs.(cell)
    all(abs.(native.cell_bohr.-cell).<=cell_bound) || error("Native cell outside printing precision")
    cube_bound=.5*10.0^(-p["cube_header_decimal_places"])
    all(abs.(cube.axes.-cell*Diagonal(1 ./ g.geo.dims)).<=cube_bound .+rt.*abs.(cube.axes)) || error("Cube axes")
    length(native.atoms)==length(g.geo.geometry.atoms) || error("Different atoms")
    for (a,b) in zip(native.atoms,g.geo.geometry.atoms)
        D.P.Z[a.element]==b.z || error("Different species")
        expected=b.position_nm/D.G.C.BOHR_NM
        delta=a.position_bohr-expected
        delta .-= diag(cell).*round.(delta./diag(cell))
        bound=.5*10.0^(-p["native_tau_decimal_places"])*native.alat .+
            .5*10.0^(-p["native_alat_decimal_places"]).*abs.(a.tau) .+rt.*diag(cell)
        all(abs.(delta).<=bound) || error("Native atomic coordinates")
    end
    maxerror=0.; maxratio=0.; count=0
    for k in axes(native.grid,3), j in axes(native.grid,2), i in axes(native.grid,1)
        a=native.grid[i,j,k]; b=cube.grid[k,j,i]
        bound=N.halfquantum(a,p["native_significant_digits"])+N.halfquantum(b,p["cube_significant_digits"])+p["numeric_roundoff_atol_ry"]
        err=abs(a-b); err<=bound || error("Native/cube mismatch at $((i,j,k)): $err > $bound")
        maxerror=max(maxerror,err); maxratio=max(maxratio,err/bound); count+=1
    end
    (;voxels=count,max_abs_difference_ry=maxerror,max_error_to_print_bound=maxratio)
end

function analyze(run,m)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    m in MOLECULES || error("Unknown accepted state")
    s=settings(joinpath(run,"settings.toml")); p=s["preprocessing"]; dir=joinpath(run,m)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml")); g=case_geometry(run,m)
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed configuration")
    S.sha(joinpath(dir,"data-file-schema.xml"))==meta["xml_sha256"] || error("Changed XML")
    Dict(r.pseudo=>r.sha256 for r in g.rs)==meta["pseudo_sha256"] || error("Changed active PAW")
    for ext in ("dat","cube")
        S.sha(joinpath(dir,"total.$ext"))==S.sha(joinpath(dir,"repeat_total.$ext")) || error("Inexact total repeat")
    end
    total=N.readplot(joinpath(dir,"total.dat")); electro=N.readplot(joinpath(dir,"electrostatic.dat"))
    total.plot_num==s["model"]["total_plot_num"] && electro.plot_num==s["model"]["electrostatic_plot_num"] || error("Wrong observable")
    checks=NamedTuple[]
    for (stem,native) in (("total",total),("electrostatic",electro))
        cube=V.cube(joinpath(dir,stem*".cube"))
        push!(checks,merge((;observable=stem),validate_grid(native,cube,g,p)))
    end
    sp=S.spectrum(joinpath(dir,"data-file-schema.xml"),s)
    bands=filter(r->r.ildos_weight>0,sp.rows); planes=NamedTuple[]; barriers=NamedTuple[]
    for k in axes(total.grid,3)
        z=(k-1)*g.geo.geometry.cell_nm[3]/g.geo.dims[3]
        stats=P.plane_stats(view(total.grid,:,:,k),view(electro.grid,:,:,k))
        location=(;molecule=m,k=k-1,z_nm=z,paw_gap=g.lower<z<g.upper,lower_half=g.lower<z<g.midpoint)
        push!(planes,merge(location,stats,(;mean_total_minus_fermi_ev=stats.mean_total_ry*S.RY_EV-sp.fermi_ev)))
        for b in bands
            push!(barriers,merge(location,(;band=b.band,energy_ev=b.energy_ev),P.barrier(stats,b.energy_ev,S.RY_EV,D.G.C.BOHR_NM)))
        end
    end
    out=joinpath(dir,"analysis"); S.newdir(out)
    S.table(joinpath(out,"planes.tsv"),planes); S.table(joinpath(out,"barriers.tsv"),barriers)
    S.table(joinpath(out,"export_comparison.tsv"),checks)
    hashes=Dict(stem*"."*ext=>S.sha(joinpath(dir,stem*"."*ext)) for stem in ("total","repeat_total","electrostatic") for ext in ("dat","cube"))
    open(io->TOML.print(io,Dict("config_sha256"=>meta["config_sha256"],"xml_sha256"=>meta["xml_sha256"],
        "export_sha256"=>hashes,"planes"=>length(planes),"selected_bands"=>length(bands),
        "paw_free_planes"=>count(r->r.paw_gap,planes),"exact_total_repeat"=>true,
        "matching_plane_adopted"=>false,"wavefunctions_propagated"=>false)),joinpath(out,"summary.toml"),"w")
    println(m,": verified ",sum(r.voxels for r in checks)," native/cube comparisons; ",length(planes)," planes, ",length(bands)," bands")
end

function main(args)
    args==["--help"] && return println("qe_vacuum_potential.jl prepare ROOT CONFIG NEW_RUN | analyze RUN MOLECULE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:4]...)
    length(args)==3 && args[1]=="analyze" && return analyze(args[2:3]...)
    error("Use --help")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEVacuumPotential.main(ARGS)
