#!/usr/bin/env julia
# Fixed-geometry 360/720-Ry sensitivity, never a predictor or mold selector.
module QEDensityCutoff
using TOML, Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"qe_vacuum_potential.jl"))
const R=QEVacuumPotential
const S=R.S; const N=R.N; const V=R.V
const STEMS=("density","total","electrostatic")
species_block(text)=S.onlymatch(r"<atomic_species\b[^>]*>(.*?)</atomic_species>"s,text,"species").captures[1]

function settings(file)
    s=TOML.parsefile(file)
    Set(keys(s))==Set(["model","selection","preprocessing"]) || error("Unexpected sections")
    m=s["model"]
    Set(keys(m))==Set(["qe_version","baseline_run","stock_pp","baseline_ecutrho_ry",
        "candidate_ecutrho_ry","ecutwfc_ry","functional","scf_acceptance_ry","sample_bias_ry",
        "startingpot","startingwfc","scf_max_seconds"]) || error("Unexpected model settings")
    (m["qe_version"],m["baseline_ecutrho_ry"],m["candidate_ecutrho_ry"],m["ecutwfc_ry"],
        m["functional"],m["scf_acceptance_ry"],m["sample_bias_ry"],m["startingpot"],m["startingwfc"]) ==
        ("7.4.1",360.,720.,50.,"PBE",5e-5,-0.0220495933,"atomic","file") || error("Outside approved scope")
    0<m["scf_max_seconds"]<=10800 || error("Insufficient reserved analysis time")
    s["selection"]==Dict("policy"=>"two_point_sensitivity_no_cutoff_or_mold_selection",
        "require_collected_converged_state"=>true,"require_exact_total_repeat"=>true,
        "adopt_reference"=>false,"propagate_wavefunctions"=>false) || error("Selection is forbidden")
    s["preprocessing"]==Dict("domain"=>"all_native_planes_with_unchanged_geometry_only_gap_tags",
        "cross_grid_interpolation"=>false,"clip_density"=>false,"normalize_density"=>false) || error("Changed domain")
    s
end

function state(xml,s,ecutrho)
    text=read(xml,String); inp=S.tag(text,"input"); out=S.tag(text,"output")
    sp=S.spectrum(xml,s)
    geo=R.D.xml_geometry(xml)
    basis=S.tag(out,"basis_set")
    geo.ecutwfc_ry==s["model"]["ecutwfc_ry"] && 2S.number(S.tag(basis,"ecutrho"))==ecutrho || error("Changed cutoffs")
    dft=S.tag(inp,"dft")
    S.tag(dft,"functional")==s["model"]["functional"] || error("Changed functional")
    S.tag(dft,"vdw_corr")=="grimme-d3" && S.tag(dft,"dftd3_version")=="3" &&
        S.tag(dft,"dftd3_threebody")=="true" || error("Changed dispersion")
    bands=S.tag(inp,"bands")
    S.number(S.tag(bands,"tot_charge"))==0 && S.tag(bands,"occupations")=="smearing" || error("Changed charge/occupations")
    sp.smearing_ev==.02S.RY_EV || error("Changed smearing")
    all(S.tag(S.tag(inp,"spin"),k)=="false" for k in ("lsda","noncolin","spinorbit")) || error("Changed spin")
    (;text,inp,out,sp,geo)
end

function input_text(xml,m,s)
    m in R.MOLECULES || error("Unknown molecule")
    a=state(xml,s,s["model"]["baseline_ecutrho_ry"]); cfg=s["model"]
    ec=S.tag(a.inp,"electron_control")
    S.tag(ec,"mixing_mode")=="plain" && S.number(S.tag(ec,"mixing_beta"))==.3 || error("Unexpected mixing")
    target=2S.number(S.tag(ec,"conv_thr")) # XML Hartree -> PW input Ry
    0<target<=cfg["scf_acceptance_ry"] || error("Original SCF target exceeds the common acceptance")
    species=collect(eachmatch(r"<species\s+name=\"([^\"]+)\">(.*?)</species>"s,species_block(a.inp)))
    structure=S.onlymatch(r"<atomic_structure\b[^>]*>(.*?)</atomic_structure>"s,a.out,"structure").captures[1]
    atoms=collect(eachmatch(r"<atom\s+name=\"([^\"]+)\"[^>]*>(.*?)</atom>"s,S.tag(structure,"atomic_positions")))
    io=IOBuffer()
    println(io,"&CONTROL\n calculation='scf'\n restart_mode='from_scratch'\n prefix='$(m)_central'\n outdir='./work'\n pseudo_dir='../pseudo'\n disk_io='low'\n tprnfor=.true.\n tstress=.true.\n max_seconds=$(cfg["scf_max_seconds"])\n/")
    println(io,"&SYSTEM\n ibrav=0\n nat=$(length(atoms))\n ntyp=$(length(species))\n nbnd=$(length(a.sp.rows))\n ecutwfc=$(cfg["ecutwfc_ry"])\n ecutrho=$(cfg["candidate_ecutrho_ry"])\n occupations='smearing'\n smearing='mv'\n degauss=0.02\n input_dft='PBE'\n vdw_corr='grimme-d3'\n dftd3_version=3\n dftd3_threebody=.true.\n/")
    println(io,"&ELECTRONS\n conv_thr=$target\n electron_maxstep=$(S.tag(ec,"max_nstep"))\n mixing_mode='plain'\n mixing_beta=0.3\n mixing_ndim=$(S.tag(ec,"mixing_ndim"))\n diagonalization='david'\n startingpot='$(cfg["startingpot"])'\n startingwfc='$(cfg["startingwfc"])'\n/\nATOMIC_SPECIES")
    for row in species
        body=row.captures[2]; pseudo=S.tag(body,"pseudo_file")
        basename(pseudo)==pseudo || error("Unsafe pseudo path")
        println(io,row.captures[1]," ",S.tag(body,"mass")," ",pseudo)
    end
    println(io,"CELL_PARAMETERS bohr")
    for i in 1:3; println(io,S.tag(S.tag(structure,"cell"),"a$i")); end
    println(io,"ATOMIC_POSITIONS bohr")
    for row in atoms; println(io,row.captures[1]," ",strip(row.captures[2])); end
    println(io,"K_POINTS gamma")
    String(take!(io))
end

function pp_input(m,stem)
    m in R.MOLECULES || error("Unknown molecule")
    num=stem=="density" ? 0 : stem in ("total","repeat_total") ? 1 : stem=="electrostatic" ? 11 : error("Unknown observable")
    "&INPUTPP\n prefix='$(m)_central'\n outdir='./work'\n filplot='$stem.dat'\n plot_num=$num\n/\n&PLOT\n iflag=3\n output_format=6\n fileout='$stem.cube'\n/\n"
end

function prepare(root,config,run)
    s=settings(config); old=joinpath(root,s["model"]["baseline_run"])
    pp=joinpath(root,s["model"]["stock_pp"]); isfile(pp) || error("Missing validated stock PP")
    inputs=String[]
    for (j,m) in enumerate(R.MOLECULES)
        xml=joinpath(old,m,"data-file-schema.xml")
        S.sha(xml)==R.D.G.XML_SHA[j] || error("Changed accepted geometry/state")
        # Checking the existing report does not re-read any scientific volume locally.
        report=TOML.parsefile(joinpath(old,m,"analysis/summary.toml"))
        report["exact_same_build_controls_and_repeats"]===true || error("Unqualified baseline")
        meta=TOML.parsefile(joinpath(old,m,"metadata.toml"))
        _,radii=R.D.radii(old,joinpath(old,m))
        Dict(r.pseudo=>r.sha256 for r in radii)==meta["pseudo_sha256"] || error("Changed PAW files")
        push!(inputs,input_text(xml,m,s))
    end
    S.newdir(run); cp(config,joinpath(run,"settings.toml"))
    cp(realpath(joinpath(old,"pseudo")),joinpath(run,"pseudo"))
    for (j,m) in enumerate(R.MOLECULES)
        dir=joinpath(run,m); mkdir(dir)
        cp(realpath(joinpath(old,m,"data-file-schema.xml")),joinpath(dir,"accepted.xml"))
        cp(joinpath(old,m,"controls/summary.toml"),joinpath(dir,"baseline_controls.toml"))
        write(joinpath(dir,"warmstart.sha256"),R.D.WFC_SHA[j]*"  work/$(m)_central.save/wfc1.dat\n")
        write(joinpath(dir,"pw_scf.in"),inputs[j])
        for stem in (STEMS...,"repeat_total"); write(joinpath(dir,stem*".in"),pp_input(m,stem)); end
        cp(joinpath(root,"hpc/qe_density_cutoff.sbatch"),joinpath(dir,"run_cutoff.sbatch"))
        open(joinpath(dir,"input_hashes.sha256"),"w") do io
            for file in ("accepted.xml","pw_scf.in","density.in","total.in","electrostatic.in","repeat_total.in",
                    "run_cutoff.sbatch","warmstart.sha256","baseline_controls.toml","../settings.toml")
                println(io,S.sha(joinpath(dir,file)),"  ",file)
            end
            for file in sort(readdir(joinpath(run,"pseudo")))
                println(io,S.sha(joinpath(run,"pseudo",file)),"  ../pseudo/",file)
            end
        end
    end
    open(io->println(io,S.sha(pp),"  ",s["model"]["stock_pp"]),joinpath(run,"stock_pp.sha256"),"w")
    println("Prepared two fixed-geometry inputs; no volumes loaded, no jobs submitted")
end

function check_scf(run,m)
    s=settings(joinpath(run,"settings.toml")); dir=joinpath(run,m)
    m in R.MOLECULES || error("Unknown molecule")
    src=joinpath(dir,"work",m*"_central.save/data-file-schema.xml")
    old=state(joinpath(dir,"accepted.xml"),s,s["model"]["baseline_ecutrho_ry"])
    a=state(src,s,s["model"]["candidate_ecutrho_ry"])
    length(a.geo.geometry.atoms)==length(old.geo.geometry.atoms) || error("Changed atom count")
    isapprox(a.geo.cell,old.geo.cell;rtol=1e-12,atol=0) || error("Changed cell")
    for (x,y) in zip(a.geo.geometry.atoms,old.geo.geometry.atoms)
        x.z==y.z && isapprox(x.position_nm,y.position_nm;rtol=0,atol=1e-12) || error("Changed geometry")
    end
    length(a.sp.rows)==length(old.sp.rows) || error("Changed band count")
    for key in ("conv_thr","mixing_mode","mixing_beta","mixing_ndim","max_nstep","diagonalization")
        S.tag(S.tag(a.inp,"electron_control"),key)==S.tag(S.tag(old.inp,"electron_control"),key) || error("Changed electronic setting: $key")
    end
    S.number(S.tag(S.tag(a.out,"band_structure"),"nelec"))==
        S.number(S.tag(S.tag(old.out,"band_structure"),"nelec")) || error("Changed electron count")
    split(species_block(a.inp))==split(species_block(old.inp)) || error("Changed species/pseudopotentials")
    a.geo.dims!=old.geo.dims || error("Density grid did not change")
    log=read(joinpath(dir,"pw_scf.out"),String)
    occursin("convergence has been achieved",log) && occursin("JOB DONE.",log) || error("Unfinished SCF")
    occursin("Starting wfcs from file",log) || error("Wavefunction warm start was not used")
    occursin("Initial potential from superposition of free atoms",log) || error("Unexpected initial potential")
    cp(realpath(src),joinpath(dir,"data-file-schema.xml"))
    println(m,": accepted 720-Ry SCF error=",a.sp.scf_error_ry," Ry; FFT=",a.geo.dims)
end

function plane_rows(density,total,electro,g,m,cutoff)
    size(density)==size(total)==size(electro) || error("Mismatched fields")
    all(isfinite,density) || error("Nonfinite density")
    rows=NamedTuple[]
    for k in axes(total,3)
        z=(k-1)*g.geo.geometry.cell_nm[3]/g.geo.dims[3]
        rho=vec(collect(view(density,:,:,k)))
        push!(rows,merge((;molecule=m,ecutrho_ry=cutoff,k=k-1,z_nm=z,
            paw_gap=g.lower<z<g.upper,lower_half=g.lower<z<g.midpoint),
            R.P.plane_stats(view(total,:,:,k),view(electro,:,:,k)),
            (;negative_density=count(<(0),rho),negative_density_fraction=count(<(0),rho)/length(rho),
                min_density=minimum(rho),mean_density=mean(rho),max_density=maximum(rho),
                std_density=std(rho;corrected=false),mean_negative_density=mean(x->max(-x,0.),rho))))
    end
    rows
end

function analyze(root,run,m)
    s=settings(joinpath(run,"settings.toml")); dir=joinpath(run,m); old=joinpath(root,s["model"]["baseline_run"])
    p=R.settings(joinpath(root,"config/qe_vacuum_potential.toml"))["preprocessing"]
    for ext in ("dat","cube")
        S.sha(joinpath(dir,"total.$ext"))==S.sha(joinpath(dir,"repeat_total.$ext")) || error("Total repeat differs")
    end
    checks=NamedTuple[]; rows=NamedTuple[]; summaries=NamedTuple[]; hashes=Dict{String,String}()
    for (name,cutoff,source,geometry_run) in (("baseline",s["model"]["baseline_ecutrho_ry"],joinpath(old,m,"stock"),old),
            ("candidate",s["model"]["candidate_ecutrho_ry"],dir,run))
        g=R.case_geometry(geometry_run,m); arrays=Dict{String,Any}()
        for stem in STEMS
            file=joinpath(source,stem*".dat"); cube=joinpath(source,stem*".cube")
            if name=="baseline"
                expected=TOML.parsefile(joinpath(dir,"baseline_controls.toml"))["control_exports_sha256"]
                S.sha(file)==expected["stock/$stem.dat"] && S.sha(cube)==expected["stock/$stem.cube"] || error("Changed baseline fields")
            end
            a=N.readplot(file); b=V.cube(cube)
            a.plot_num==(stem=="density" ? 0 : stem=="total" ? 1 : 11) || error("Wrong export")
            check=R.validate_grid(a,b,g,p)
            push!(checks,(;source=name,field=stem,unit=stem=="density" ? "electron_Bohr-3" : "Ry",
                voxels=check.voxels,max_abs_difference=check.max_abs_difference_ry,
                max_error_to_print_bound=check.max_error_to_print_bound))
            arrays[stem]=a.grid
            hashes[name*"/"*stem*".dat"]=S.sha(file); hashes[name*"/"*stem*".cube"]=S.sha(cube)
        end
        append!(rows,plane_rows(arrays["density"],arrays["total"],arrays["electrostatic"],g,m,cutoff))
        rho=arrays["density"]; voxel=prod(g.geo.geometry.cell_nm)/R.D.G.C.BOHR_NM^3/length(rho)
        push!(summaries,(;source=name,ecutrho_ry=cutoff,nx=g.geo.dims[1],ny=g.geo.dims[2],nz=g.geo.dims[3],
            lower_paw_nm=g.lower,upper_paw_nm=g.upper,midpoint_nm=g.midpoint,
            negative_charge_electrons=sum(x->max(-x,0.),rho)*voxel,charge_electrons=sum(rho)*voxel))
        empty!(arrays); GC.gc()
    end
    out=joinpath(dir,"analysis"); S.newdir(out)
    S.table(joinpath(out,"planes.tsv"),rows); S.table(joinpath(out,"grids.tsv"),summaries)
    S.table(joinpath(out,"export_comparison.tsv"),checks)
    open(io->TOML.print(io,Dict("config_sha256"=>S.sha(joinpath(run,"settings.toml")),
        "xml_sha256"=>S.sha(joinpath(dir,"data-file-schema.xml")),"export_sha256"=>hashes,
        "exact_total_repeat"=>true,"cross_grid_interpolation"=>false,"potential_adopted"=>false,
        "scf_acceptance_ry"=>s["model"]["scf_acceptance_ry"])),joinpath(out,"summary.toml"),"w")
    println(m,": ",sum(r.voxels for r in checks)," native/cube checks; all native planes retained")
end

function main(args)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_density_cutoff.jl prepare ROOT CONFIG NEW_RUN | check-scf RUN MOLECULE | analyze ROOT RUN MOLECULE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==3 && args[1]=="check-scf" && return check_scf(args[2:end]...)
    length(args)==4 && args[1]=="analyze" && return analyze(args[2:end]...)
    error("Use --help")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEDensityCutoff.main(ARGS)
