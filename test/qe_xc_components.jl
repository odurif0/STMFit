#!/usr/bin/env julia
# Accepted-state component attribution; never a predictor or potential repair.
module QEXCComponents
using TOML, Statistics
include(joinpath(@__DIR__, "qe_vacuum_potential.jl"))
const R=QEVacuumPotential
const S=R.S
const N=R.N
const V=R.V
const FIELDS=("rho_valence","rho_xc_input","grad2_xc_input","gga_active","lda","gga_local","gga_divergence")

function settings(file)
    actual=TOML.parsefile(file)
    expected=Dict("model"=>Dict("qe_version"=>"7.4.1",
        "qe_source_commit"=>"500de340b820e1cb8c05f2d8bb8fced102f377c1",
        "source_run"=>"qe/vacuum_potential_20260926","functional"=>"internal_nonmagnetic_PBE_1_4_3_4",
        "components"=>collect(FIELDS),"scf_acceptance_ry"=>5e-5),
        "selection"=>Dict("policy"=>"all_native_planes_and_source_gga_branches_no_physical_selection",
            "require_total_reference_exact"=>true,"require_density_control_exact"=>true,"require_component_repeat_exact"=>true,
            "adopt_reference"=>false,"modify_density"=>false,"propagate_wavefunctions"=>false),
        "preprocessing"=>Dict("rho_threshold_lda"=>1e-10,"rho_threshold_gga"=>1e-6,
            "gradient_squared_threshold_gga"=>1e-10,"native_significant_digits"=>10,
            "cube_significant_digits"=>5,"addition_roundoff_ulps"=>64,
            "normalize_components"=>false,"clip_density"=>false))
    actual==expected || error("Changed XC attribution experiment")
    actual
end

function extraction_input(m,repeat=false)
    m in R.MOLECULES || error("Unknown accepted state")
    dir=repeat ? "../../work/$m" : "../work/$m"
    "&INPUTPP\n prefix='$(m)_central'\n outdir='$dir'\n plot_num=1\n filplot='xc_total.dat'\n/\n" *
        "&PLOT\n iflag=3\n output_format=6\n fileout='xc_total.cube'\n/\n"
end

function cube_input(name)
    name in FIELDS || error("Unknown component")
    "&INPUTPP\n/\n&PLOT\n nfile=1\n filepp(1)='$name.dat'\n weight(1)=1\n iflag=3\n" *
        " output_format=6\n fileout='$name.cube'\n/\n"
end

function prepare(root,config,out)
    s=settings(config); old=joinpath(root,s["model"]["source_run"])
    R.settings(joinpath(old,"settings.toml"))
    S.newdir(out)
    cp(config,joinpath(out,"settings.toml"))
    for (source,dest) in (("settings.toml","potential_settings.toml"),("geometry.tsv","geometry.tsv"),
            ("accepted_states.sha256","accepted_states.sha256"),("pseudo","pseudo"))
        cp(joinpath(old,source),joinpath(out,dest))
    end
    cp(joinpath(root,"hpc/qe_xc_components.sbatch"),joinpath(out,"run_xc.sbatch"))
    for (j,m) in enumerate(R.MOLECULES)
        before=joinpath(old,m); dir=joinpath(out,m); mkdir(dir); mkdir(joinpath(dir,"repeat"))
        meta=TOML.parsefile(joinpath(before,"metadata.toml"))
        meta["xml_sha256"]==R.D.G.XML_SHA[j]==S.sha(joinpath(before,"data-file-schema.xml")) || error("Changed accepted XML")
        oldresult=TOML.parsefile(joinpath(before,"analysis/summary.toml"))
        meta["reference_hashes"]=oldresult["export_sha256"]
        meta["config_sha256"]=S.sha(config)
        meta["geometry_sha256"]=S.sha(joinpath(out,"geometry.tsv"))
        cp(joinpath(before,"data-file-schema.xml"),joinpath(dir,"data-file-schema.xml"))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
        write(joinpath(dir,"extract.in"),extraction_input(m))
        write(joinpath(dir,"repeat/extract.in"),extraction_input(m,true))
        write(joinpath(dir,"density.in"),replace(extraction_input(m),"plot_num=1"=>"plot_num=0", "xc_total"=>"density_control"))
        for field in FIELDS; write(joinpath(dir,"cube_"*field*".in"),cube_input(field)); end
    end
    println("Prepared frozen-state XC extraction; no density or component values read")
end

function additivity(total,electro,lda,local_gga,divergence,p)
    arrays=(total,electro,lda,local_gga,divergence)
    all(a->size(a)==size(total),arrays) || error("Different component dimensions")
    maxerror=0.; maxratio=0.
    for q in eachindex(total)
        values=map(a->a[q],arrays)
        all(isfinite,values) || error("Nonfinite component")
        err=abs(values[1]-values[2]-values[3]-values[4]-values[5])
        bound=sum(N.halfquantum(x,p["native_significant_digits"]) for x in values)+
            p["addition_roundoff_ulps"]*eps(Float64)*max(1.,sum(abs,values))
        err<=bound || error("XC additivity at $q: $err > $bound")
        maxerror=max(maxerror,err); maxratio=max(maxratio,err/bound)
    end
    (;voxels=length(total),max_abs_error_ry=maxerror,max_error_to_bound=maxratio)
end

function source_masks(rho,grad2,active,local_gga,p)
    all(a->size(a)==size(rho),(grad2,active,local_gga)) || error("Mask dimensions differ")
    definite=0; ambiguous=0
    for q in eachindex(rho)
        r=rho[q]; g=grad2[q]; a=active[q]; l=local_gga[q]
        all(isfinite,(r,g,a,l)) && g>=0 && a in (0.,1.) || error("Invalid source mask/value")
        a==0 && l!=0 && error("Inactive GGA has nonzero local term")
        rbound=N.halfquantum(r,p["native_significant_digits"])
        gbound=N.halfquantum(g,p["native_significant_digits"])
        rr=p["rho_threshold_gga"]; gg=p["gradient_squared_threshold_gga"]
        if abs(abs(r)-rr)<=rbound || abs(g-gg)<=gbound
            ambiguous+=1 # retain exact source mask; rounded values cannot adjudicate this case
        else
            (a==1)==(abs(r)>rr && g>gg) || error("GGA branch mismatch")
            definite+=1
        end
    end
    (;definite,ambiguous)
end

function description(v)
    isempty(v) && return (;samples=0,minimum=NaN,mean=NaN,maximum=NaN,spatial_sd=NaN)
    mu=mean(v)
    (;samples=length(v),minimum=minimum(v),mean=mu,maximum=maximum(v),spatial_sd=sqrt(mean(abs2,v.-mu)))
end

function analyze(root,run,m)
    m in R.MOLECULES || error("Unknown state")
    s=settings(joinpath(run,"settings.toml")); p=s["preprocessing"]
    old=joinpath(root,s["model"]["source_run"],m); dir=joinpath(run,m)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml")); g=R.case_geometry(run,m)
    ps=R.settings(joinpath(run,"potential_settings.toml"))["preprocessing"]
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed configuration")
    S.sha(joinpath(run,"geometry.tsv"))==meta["geometry_sha256"] || error("Changed geometry domains")
    S.sha(joinpath(dir,"data-file-schema.xml"))==meta["xml_sha256"] || error("Changed state")
    Dict(r.pseudo=>r.sha256 for r in g.rs)==meta["pseudo_sha256"] || error("Changed active PAW datasets")
    thresholds=TOML.parsefile(joinpath(dir,"xc_thresholds.toml"))
    thresholds==Dict(k=>p[k] for k in ("rho_threshold_lda","rho_threshold_gga","gradient_squared_threshold_gga")) || error("Changed QE thresholds")
    for ext in ("dat","cube")
        S.sha(joinpath(dir,"xc_total.$ext"))==meta["reference_hashes"]["total.$ext"]==S.sha(joinpath(old,"total.$ext")) || error("Total control differs")
        S.sha(joinpath(dir,"repeat/xc_total.$ext"))==S.sha(joinpath(dir,"xc_total.$ext")) || error("Total repeat differs")
    end
    S.sha(joinpath(old,"electrostatic.dat"))==meta["reference_hashes"]["electrostatic.dat"] || error("Changed electrostatic reference")
    S.sha(joinpath(dir,"rho_valence.dat"))==S.sha(joinpath(dir,"density_control.dat")) || error("Density control differs")
    analyze_fields(run,m,s,meta,g,joinpath(old,"electrostatic.dat"),"exact_total_reference_and_repeats")
end

# Shared arithmetic only. Each caller must first verify its own frozen controls;
# the original entrypoint above still requires the archived site bytes exactly.
function analyze_fields(run,m,s,meta,g,electrofile,control_key)
    p=s["preprocessing"]; dir=joinpath(run,m)
    ps=R.settings(joinpath(run,"potential_settings.toml"))["preprocessing"]
    values=Dict{String,Any}(); comparisons=NamedTuple[]; hashes=Dict{String,String}()
    for field in FIELDS
        file=joinpath(dir,field*".dat")
        S.sha(file)==S.sha(joinpath(dir,"repeat",field*".dat")) || error("Component repeat differs")
        native=N.readplot(file); cube=V.cube(joinpath(dir,field*".cube"))
        native.plot_num==(field in ("lda","gga_local","gga_divergence") ? 1 : 0) || error("Wrong component header")
        push!(comparisons,merge((;component=field),R.validate_grid(native,cube,g,ps)))
        values[field]=native.grid
        for ext in ("dat","cube"); hashes[field*"."*ext]=S.sha(joinpath(dir,field*"."*ext)); end
    end
    total=N.readplot(joinpath(dir,"xc_total.dat")); electro=N.readplot(electrofile)
    sums=additivity(total.grid,electro.grid,values["lda"],values["gga_local"],values["gga_divergence"],p)
    masks=source_masks(values["rho_xc_input"],values["grad2_xc_input"],values["gga_active"],values["gga_local"],p)
    values["xc_total"]=total.grid.-electro.grid
    rows=NamedTuple[]; counts=NamedTuple[]
    for k in axes(total.grid,3)
        z=(k-1)*g.geo.geometry.cell_nm[3]/g.geo.dims[3]
        location=(;molecule=m,k=k-1,z_nm=z,paw_gap=g.lower<z<g.upper,lower_half=g.lower<z<g.midpoint)
        active=view(values["gga_active"],:,:,k)
        rho=view(values["rho_xc_input"],:,:,k); valence=view(values["rho_valence"],:,:,k)
        push!(counts,merge(location,(;samples=length(active),gga_active=count(==(1),active),
            negative_valence=count(<(0),valence),negative_xc_input=count(<(0),rho),
            below_lda=count(x->abs(x)<=p["rho_threshold_lda"],rho),
            below_gga_density=count(x->abs(x)<=p["rho_threshold_gga"],rho))))
        for subset in ("all","gga_active","gga_inactive")
            mask=subset=="all" ? trues(size(active)) : active.==(subset=="gga_active" ? 1. : 0.)
            for field in (FIELDS...,"xc_total")
                vector=view(values[field],:,:,k)[mask]
                unit=field in ("lda","gga_local","gga_divergence","xc_total") ? "Ry" :
                    field=="gga_active" ? "binary" : field=="grad2_xc_input" ? "electron2_Bohr-8" : "electron_Bohr-3"
                push!(rows,merge(location,(;subset,component=field,unit),description(vector)))
            end
        end
    end
    out=joinpath(dir,"analysis"); S.newdir(out)
    S.table(joinpath(out,"components.tsv"),rows); S.table(joinpath(out,"counts.tsv"),counts)
    S.table(joinpath(out,"export_comparison.tsv"),comparisons)
    vol=prod(g.geo.geometry.cell_nm)/R.D.G.C.BOHR_NM^3/length(total.grid)
    summary=Dict("config_sha256"=>meta["config_sha256"],"component_sha256"=>hashes,
        "native_cube_voxels"=>sum(r.voxels for r in comparisons),"planes"=>length(counts),
        "additivity"=>Dict(string(k)=>v for (k,v) in pairs(sums)),
        "source_masks"=>Dict(string(k)=>v for (k,v) in pairs(masks)),
        "negative_valence_integral_electrons"=>sum(x->max(-x,0.),values["rho_valence"])*vol,
        control_key=>true,"density_modified"=>false,"potential_adopted"=>false)
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(m,": ",summary["native_cube_voxels"]," export comparisons; ",sums,"; source masks ",masks)
end

function main(args)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_xc_components.jl prepare ROOT CONFIG NEW_RUN | analyze ROOT RUN MOLECULE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1]=="analyze" && return analyze(args[2:end]...)
    error("Use --help")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEXCComponents.main(ARGS)
