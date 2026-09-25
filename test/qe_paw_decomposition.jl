#!/usr/bin/env julia
# PP-only attribution of signed density; never a mold exporter or predictor.
module QEPAWDecomposition
using TOML, Statistics
include(joinpath(@__DIR__,"qe_gamma_reconstruction.jl"))
const G=QEGammaReconstruction
const S=G.S
const FULL_SHA=["f9cb76e0fadb23a156af57f329f3e8d3985376021dbad270c8b535c997220911",
    "fd17389087a27fdf8d5ec2a6644ce0b2ace501b0b7d4c4de74389d573876ca35"]
const NAMES=("full_control","smooth_control","augmentation_control","full_ildos",
    "smooth_ildos","repeat_smooth","augmentation_ildos","repeat_augmentation")

function settings(file)
    expected=G.settings(joinpath(@__DIR__,"../config/qe_gamma_reconstruction.toml"))
    expected["model"]["decomposition"]="smooth_plus_reciprocal_augmentation"
    delete!(expected["selection"],"require_patched_repeat_sha256")
    expected["selection"]["require_corrected_total_sha256"]=true
    expected["selection"]["require_component_repeat_sha256"]=true
    expected["preprocessing"]["normalize_components"]=false
    actual=TOML.parsefile(file)
    actual==expected || error("Not the frozen PAW component experiment")
    actual
end

function additivity(total,smooth,augmentation,s)
    length(total)==length(smooth)==length(augmentation) || error("Grid length mismatch")
    digits=s["selection"]["cube_significant_digits"]; atol=s["selection"]["numeric_roundoff_atol"]
    violations=0; maxerror=0.; maxratio=0.
    for (t,p,a) in zip(total,smooth,augmentation)
        all(isfinite,(t,p,a)) || error("Nonfinite component")
        bound=(G.quantum(t,digits)+G.quantum(p,digits)+G.quantum(a,digits))/2+atol
        err=abs(t-p-a); violations+=err>bound
        maxerror=max(maxerror,err); maxratio=max(maxratio,err/bound)
    end
    (;violations,max_abs_error=maxerror,max_error_to_bound=maxratio)
end

function prepare(root,config,outdir)
    s=settings(config); old=joinpath(root,"qe/gamma_reconstruction_mkl_20260925")
    G.settings(joinpath(old,"settings.toml"))
    cases=[joinpath(old,m) for m in ("glcn","glcnac")]
    S.sha.([joinpath(d,"reference.cube") for d in cases])==G.CUBE_SHA || error("Wrong STM references")
    S.sha.([joinpath(d,"patched_ildos.cube") for d in cases])==FULL_SHA || error("Wrong corrected totals")
    S.sha.([joinpath(d,"data-file-schema.xml") for d in cases])==G.XML_SHA || error("Wrong accepted states")
    spectra=[S.spectrum(joinpath(d,"data-file-schema.xml"),s) for d in cases]
    S.require_plane_domains([joinpath(d,"reference.cube") for d in cases],
        [joinpath(d,"frame.tsv") for d in cases],s)
    S.newdir(outdir); cp(config,joinpath(outdir,"settings.toml"))
    for file in ("collected_wfc.sha256","pw_scf.in","pp_ldos.in","pseudo")
        cp(joinpath(old,file),joinpath(outdir,file))
    end
    cp(joinpath(root,"hpc/qe_paw_decomposition.sbatch"),joinpath(outdir,"run_pp_only.sbatch"))
    for (k,m) in enumerate(("glcn","glcnac"))
        case_dir=joinpath(outdir,m); mkdir(case_dir); previous=cases[k]; sp=spectra[k]
        for file in ("data-file-schema.xml","reference.cube","frame.tsv","bands.tsv","pp_ldos.in")
            cp(joinpath(previous,file),joinpath(case_dir,file))
        end
        cp(joinpath(previous,"patched_ildos.cube"),joinpath(case_dir,"reference_total.cube"))
        for name in NAMES
            plot=endswith(name,"control") ? 5 : 10
            write(joinpath(case_dir,name*".in"),S.pp_input(m*"_central","../work/$m",name,plot,sp,s))
        end
        meta=TOML.parsefile(joinpath(previous,"metadata.toml"))
        meta["scope"]="paw_component_attribution_only"
        meta["config_sha256"]=S.sha(config)
        meta["reference_total_sha256"]=FULL_SHA[k]
        meta["reference_stm_sha256"]=G.CUBE_SHA[k]
        open(io->TOML.print(io,meta),joinpath(case_dir,"metadata.toml"),"w")
    end
    println("Accepted states, fixed window and all six plane domains verified; no components selected.")
end

function compare(case_dir,config)
    s=settings(config); meta=TOML.parsefile(joinpath(case_dir,"metadata.toml"))
    meta["config_sha256"]==S.sha(config) || error("Changed settings")
    meta["xml_sha256"]==S.sha(joinpath(case_dir,"data-file-schema.xml")) || error("Changed state")
    sp=S.spectrum(joinpath(case_dir,"data-file-schema.xml"),s)
    expected=sum(r.ildos_weight for r in sp.rows)
    expected==meta["ildos_integral_expected"] || error("Changed state weight")
    meta["reference_total_sha256"]==S.sha(joinpath(case_dir,"reference_total.cube")) || error("Changed total reference")
    meta["reference_stm_sha256"]==S.sha(joinpath(case_dir,"reference.cube")) || error("Changed STM reference")
    files=Dict(n=>joinpath(case_dir,n*".cube") for n in NAMES)
    cubes=Dict(n=>S.QECubeMolds.read_qe_cube(files[n]) for n in NAMES)
    full=cubes["full_ildos"]
    all(c.dims==full.dims && c.origin_nm==full.origin_nm && c.axes_nm==full.axes_nm for c in values(cubes)) || error("Changed grid")
    sums=additivity(full.values,cubes["smooth_ildos"].values,cubes["augmentation_ildos"].values,s)
    stats=[merge((observable=n,),S.cube_stats(cubes[n])) for n in NAMES]
    hashes=Dict(n=>S.sha(files[n]) for n in NAMES)
    full_integral=only(r.integral for r in stats if r.observable=="full_ildos")
    checks=Dict("legacy_controls_exact"=>all(hashes[n]==meta["reference_stm_sha256"] for n in NAMES if endswith(n,"control")),
        "corrected_total_exact"=>hashes["full_ildos"]==meta["reference_total_sha256"],
        "smooth_repeat_exact"=>hashes["smooth_ildos"]==hashes["repeat_smooth"],
        "augmentation_repeat_exact"=>hashes["augmentation_ildos"]==hashes["repeat_augmentation"],
        "total_integral"=>abs(full_integral-expected)<=s["selection"]["ildos_integral_rtol"]*expected,
        "whole_volume_additivity"=>sums.violations==0)
    frame=G.C.read_frame(joinpath(case_dir,"frame.tsv"))
    planes=NamedTuple[]; pixels=NamedTuple[]; attribution=NamedTuple[]
    for height in s["preprocessing"]["diagnostic_heights_nm"]
        plane_values=Dict(n=>S.plane(cubes[n],frame,height,s) for n in NAMES)
        for n in NAMES
            v=plane_values[n]
            push!(planes,(observable=n,height_nm=height,values=length(v),negative=count(<(0),v),
                minimum=minimum(v),median=median(v),maximum=maximum(v),
                positive_sum=sum(max(0,x) for x in v),negative_abs_sum=-sum(min(0,x) for x in v)))
            append!(pixels,[(observable=n,height_nm=height,pixel=i,value=v[i]) for i in eachindex(v)])
        end
        t=plane_values["full_ildos"]; p=plane_values["smooth_ildos"]; a=plane_values["augmentation_ildos"]
        push!(attribution,(height_nm=height,total_negative=count(<(0),t),smooth_negative=count(<(0),p),
            augmentation_negative=count(<(0),a),negative_total_with_nonnegative_smooth=count(i->t[i]<0 && p[i]>=0,eachindex(t)),
            negative_total_with_negative_augmentation=count(i->t[i]<0 && a[i]<0,eachindex(t))))
    end
    out=joinpath(case_dir,"comparison"); S.newdir(out)
    S.table(joinpath(out,"cube_summary.tsv"),stats)
    S.table(joinpath(out,"plane_summary.tsv"),planes)
    S.table(joinpath(out,"planes.tsv"),pixels)
    S.table(joinpath(out,"attribution.tsv"),attribution)
    summary=Dict("scope"=>"components_only_no_mold_no_calibration_no_benchmark","checks"=>checks,
        "expected_total_integral"=>expected,"cube_sha256"=>hashes,
        "additivity"=>Dict(string(k)=>v for (k,v) in pairs(sums)))
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(basename(case_dir),": ",checks,"; additivity ",sums)
    all(Base.values(checks)) || error("Saved component comparison failed; no interpretation or promotion")
    nothing
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_paw_decomposition.jl prepare ROOT CONFIG NEW_DIR | compare CASE CONFIG")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==3 && args[1]=="compare" && return compare(args[2:end]...)
    error("Expected prepare or compare; no label options")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEPAWDecomposition.main()
