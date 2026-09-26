#!/usr/bin/env julia
# Same-build PP controls; cross-build differences are measured, never fitted away.
module QEXCBuildControl
using TOML, Statistics
include(joinpath(@__DIR__,"qe_xc_components.jl"))
const C=QEXCComponents
const R=C.R; const S=C.S; const N=C.N; const V=C.V
const CONTROL_FILES=[stem*"."*ext for stem in ("stock/density","stock/total","stock/electrostatic","site/density") for ext in ("dat","cube")]

function settings(file)
    expected=C.settings(joinpath(@__DIR__,"../config/qe_xc_components.toml"))
    expected["model"]["paired_build"]="stock_then_extraction_same_configure_and_libraries"
    delete!(expected["selection"],"require_total_reference_exact")
    expected["selection"]["require_same_build_total_exact"]=true
    expected["selection"]["cross_build_policy"]="report_all_native_voxels_no_accuracy_threshold_or_adoption"
    actual=TOML.parsefile(file)
    actual==expected || error("Changed paired-build experiment")
    actual
end

function control_input(m,stem)
    m in R.MOLECULES || error("Unknown state")
    number=stem=="density" ? 0 : stem=="total" ? 1 : stem=="electrostatic" ? 11 : error("Unknown control")
    "&INPUTPP\n prefix='$(m)_central'\n outdir='../../work/$m'\n plot_num=$number\n filplot='$stem.dat'\n/\n" *
        "&PLOT\n iflag=3\n output_format=6\n fileout='$stem.cube'\n/\n"
end

function prepare(root,config,out)
    s=settings(config); old=joinpath(root,s["model"]["source_run"])
    R.settings(joinpath(old,"settings.toml"))
    S.newdir(out); cp(config,joinpath(out,"settings.toml"))
    for (src,dst) in (("settings.toml","potential_settings.toml"),("geometry.tsv","geometry.tsv"),
            ("accepted_states.sha256","accepted_states.sha256"),("pseudo","pseudo"))
        cp(joinpath(old,src),joinpath(out,dst))
    end
    cp(joinpath(root,"hpc/qe_xc_build_control.sbatch"),joinpath(out,"run_build_control.sbatch"))
    for (j,m) in enumerate(R.MOLECULES)
        before=joinpath(old,m); dir=joinpath(out,m)
        mkdir(dir)
        for sub in ("stock","site","repeat"); mkdir(joinpath(dir,sub)); end
        meta=TOML.parsefile(joinpath(before,"metadata.toml"))
        meta["xml_sha256"]==R.D.G.XML_SHA[j]==S.sha(joinpath(before,"data-file-schema.xml")) || error("Changed accepted XML")
        meta["reference_hashes"]=TOML.parsefile(joinpath(before,"analysis/summary.toml"))["export_sha256"]
        meta["config_sha256"]=S.sha(config); meta["geometry_sha256"]=S.sha(joinpath(out,"geometry.tsv"))
        meta["helper_sha256"]=S.sha(joinpath(root,"hpc/qe_xc_components.f90"))
        meta["patch_sha256"]=S.sha(joinpath(root,"hpc/patches/qe-7.4.1-xc-components.patch"))
        cp(joinpath(before,"data-file-schema.xml"),joinpath(dir,"data-file-schema.xml"))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
        write(joinpath(dir,"extract.in"),C.extraction_input(m))
        write(joinpath(dir,"repeat/extract.in"),C.extraction_input(m,true))
        for stem in ("density","total","electrostatic")
            write(joinpath(dir,"stock",stem*".in"),control_input(m,stem))
        end
        write(joinpath(dir,"site/density.in"),control_input(m,"density"))
        for field in C.FIELDS; write(joinpath(dir,"cube_"*field*".in"),C.cube_input(field)); end
    end
    println("Prepared paired stock/extraction controls; no density or potential read")
end

"Descriptive differences only: no cross-build acceptance tolerance."
function differences(a,b,p)
    size(a)==size(b) && !isempty(a) || error("Different or empty fields")
    all(isfinite,a) && all(isfinite,b) || error("Nonfinite field")
    delta=a.-b; d=abs.(delta); within=0
    # Loop avoids keeping another full grid of printing bounds.
    ratio=0.
    for q in eachindex(a)
        bound=N.halfquantum(a[q],p["native_significant_digits"])+N.halfquantum(b[q],p["native_significant_digits"])
        within+=d[q]<=bound
        ratio=max(ratio,bound>0 ? d[q]/bound : iszero(d[q]) ? 0. : Inf)
    end
    # Native grids are Cartesian SubArrays, even without physical padding.
    # Their generic reduction can accumulate serially; materialization gives
    # the same pairwise sum as a dense plane. Delta is already a dense array.
    normref=sum(abs2,collect(b)); normdiff=sum(abs2,delta)
    (;samples=length(a),different=count(!iszero,delta),maximum_absolute=maximum(d),
        mean_signed=mean(delta),mean_absolute=mean(d),rms=sqrt(normdiff/length(a)),
        relative_l2=normref>0 ? sqrt(normdiff/normref) : iszero(normdiff) ? 0. : Inf,
        within_combined_print_bound=within,max_error_to_combined_print_bound=ratio)
end

function native_pair(a,b)
    for key in (:dims,:padded_dims,:cell_bohr,:alat,:at,:plot_num,:cutoffs,:species,:atoms)
        getproperty(a,key)==getproperty(b,key) || error("Changed native header $key")
    end
end

function validated_metadata(root,run,m)
    m in R.MOLECULES || error("Unknown state")
    s=settings(joinpath(run,"settings.toml")); dir=joinpath(run,m)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml")); g=R.case_geometry(run,m)
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed configuration")
    S.sha(joinpath(run,"geometry.tsv"))==meta["geometry_sha256"] || error("Changed geometry")
    S.sha(joinpath(dir,"data-file-schema.xml"))==meta["xml_sha256"] || error("Changed accepted state")
    Dict(r.pseudo=>r.sha256 for r in g.rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
    S.sha(joinpath(root,"hpc/qe_xc_components.f90"))==meta["helper_sha256"] || error("Changed helper")
    S.sha(joinpath(root,"hpc/patches/qe-7.4.1-xc-components.patch"))==meta["patch_sha256"] || error("Changed patch")
    old=joinpath(root,s["model"]["source_run"],m)
    for stem in ("total","electrostatic"), ext in ("dat","cube")
        S.sha(joinpath(old,"$stem.$ext"))==meta["reference_hashes"]["$stem.$ext"] || error("Changed site archive")
    end
    (;s,meta,g,dir,old)
end

function identity_pairs()
    rows=[("stock/density.dat","rho_valence.dat"),("stock/total.dat","xc_total.dat"),
        ("stock/total.cube","xc_total.cube"),("xc_total.dat","repeat/xc_total.dat"),
        ("xc_total.cube","repeat/xc_total.cube"),("xc_thresholds.toml","repeat/xc_thresholds.toml")]
    append!(rows,[(field*".dat","repeat/"*field*".dat") for field in C.FIELDS])
    rows
end

function report(root,run,m)
    v=validated_metadata(root,run,m); (;s,meta,g,dir,old)=v; p=s["preprocessing"]
    thresholds=TOML.parsefile(joinpath(dir,"xc_thresholds.toml"))
    thresholds==Dict(k=>p[k] for k in ("rho_threshold_lda","rho_threshold_gga","gradient_squared_threshold_gga")) || error("Changed QE thresholds")
    identities=[(;a,b,sha_a=S.sha(joinpath(dir,a)),sha_b=S.sha(joinpath(dir,b)),exact=S.sha(joinpath(dir,a))==S.sha(joinpath(dir,b))) for (a,b) in identity_pairs()]
    rows=NamedTuple[]; totals=NamedTuple[]; exports=NamedTuple[]
    ps=R.settings(joinpath(run,"potential_settings.toml"))["preprocessing"]
    for stem in ("stock/density","stock/total","stock/electrostatic","site/density")
        a=N.readplot(joinpath(dir,stem*".dat")); cube=V.cube(joinpath(dir,stem*".cube"))
        push!(exports,merge((;observable=stem),R.validate_grid(a,cube,g,ps)))
    end
    comparisons=(("same_build_density","electron_Bohr-3",joinpath(dir,"stock/density.dat"),joinpath(dir,"rho_valence.dat")),
        ("same_build_total","Ry",joinpath(dir,"stock/total.dat"),joinpath(dir,"xc_total.dat")),
        ("cross_build_density","electron_Bohr-3",joinpath(dir,"stock/density.dat"),joinpath(dir,"site/density.dat")),
        ("cross_build_total","Ry",joinpath(dir,"stock/total.dat"),joinpath(old,"total.dat")),
        ("cross_build_electrostatic","Ry",joinpath(dir,"stock/electrostatic.dat"),joinpath(old,"electrostatic.dat")))
    for (comparison,unit,afile,bfile) in comparisons
        a=N.readplot(afile); b=N.readplot(bfile); native_pair(a,b)
        push!(totals,merge((;molecule=m,comparison,unit,sha_a=S.sha(afile),sha_b=S.sha(bfile)),differences(a.grid,b.grid,p)))
        for k in axes(a.grid,3)
            z=(k-1)*g.geo.geometry.cell_nm[3]/g.geo.dims[3]
            location=(;molecule=m,comparison,unit,k=k-1,z_nm=z,paw_gap=g.lower<z<g.upper,lower_half=g.lower<z<g.midpoint)
            push!(rows,merge(location,differences(view(a.grid,:,:,k),view(b.grid,:,:,k),p)))
        end
    end
    out=joinpath(dir,"controls"); S.newdir(out)
    S.table(joinpath(out,"identities.tsv"),identities); S.table(joinpath(out,"differences.tsv"),rows)
    S.table(joinpath(out,"totals.tsv"),totals); S.table(joinpath(out,"exports.tsv"),exports)
    summary=Dict("config_sha256"=>meta["config_sha256"],"same_build_controls_exact"=>all(r.exact for r in identities),
        "control_exports_sha256"=>Dict(file=>S.sha(joinpath(dir,file)) for file in CONTROL_FILES),
        "comparison_voxels"=>sum(r.samples for r in totals),"native_cube_voxels"=>sum(r.voxels for r in exports),
        "cross_build_acceptance_decision"=>false,"potential_adopted"=>false)
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(m,": ",summary)
end

function require_controls(run)
    for m in R.MOLECULES
        dir=joinpath(run,m)
        saved=TOML.parsefile(joinpath(dir,"controls/summary.toml"))
        saved["same_build_controls_exact"] || error("Same-build control failed for $m; keep full difference report")
        saved["config_sha256"]==S.sha(joinpath(run,"settings.toml")) || error("Changed control configuration")
        Set(keys(saved["control_exports_sha256"]))==Set(CONTROL_FILES) || error("Incomplete control exports")
        for file in CONTROL_FILES
            S.sha(joinpath(dir,file))==saved["control_exports_sha256"][file] || error("Changed control export $file")
        end
        records=V.tsv(joinpath(dir,"controls/identities.tsv"))
        Set((r["a"],r["b"]) for r in records)==Set(identity_pairs()) || error("Incomplete identity report")
        for r in records
            parse(Bool,r["exact"]) && r["sha_a"]==r["sha_b"]==S.sha(joinpath(dir,r["a"]))==S.sha(joinpath(dir,r["b"])) || error("Changed saved identity")
        end
        for (a,b) in identity_pairs()
            S.sha(joinpath(dir,a))==S.sha(joinpath(dir,b)) || error("Same-build identity failed: $m $a $b")
        end
    end
end

function analyze(root,run,m)
    require_controls(run)
    (;s,meta,g,dir)=validated_metadata(root,run,m)
    thresholds=TOML.parsefile(joinpath(dir,"xc_thresholds.toml"))
    thresholds==Dict(k=>s["preprocessing"][k] for k in ("rho_threshold_lda","rho_threshold_gga","gradient_squared_threshold_gga")) || error("Changed QE thresholds")
    C.analyze_fields(run,m,s,meta,g,joinpath(dir,"stock/electrostatic.dat"),"exact_same_build_controls_and_repeats")
end

function main(args)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_xc_build_control.jl prepare ROOT CONFIG NEW_RUN | report ROOT RUN MOLECULE | check RUN | analyze ROOT RUN MOLECULE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1]=="report" && return report(args[2:end]...)
    length(args)==2 && args[1]=="check" && return require_controls(args[2])
    length(args)==4 && args[1]=="analyze" && return analyze(args[2:end]...)
    error("Use --help")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEXCBuildControl.main(ARGS)
