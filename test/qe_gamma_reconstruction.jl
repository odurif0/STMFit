#!/usr/bin/env julia
# Isolated PP-only scientific comparison. No STM images, fits, labels or molds.
module QEGammaReconstruction
using TOML, SHA, LinearAlgebra, Statistics, Printf
include(joinpath(@__DIR__,"qe_spectral_window.jl"))
const S=QESpectralWindow
const C=S.C
const QE_COMMIT="500de340b820e1cb8c05f2d8bb8fced102f377c1"
const XML_SHA=["0cbcb7f6b1f533c24b4e05f08f5bb5b279eb9255d4c0ac04e445c04d69b057c8",
    "3b25c35dad13754f505b0b923ec8845c3304e6a84b11f40cf0f90c2710ccc124"]
const CUBE_SHA=["80cd1d1fde94cf084cc7ea464d2bf065b36b2c015ea0bfaef8cebfee8ff88863",
    "40649ccd9eb6768444b8ff61bf4a639b3940cf926fe3eb42254eb024faf9b5bf"]
const WFC_SHA=["31f137f580bab932c453d8c3e79e9f19e993ce98ed3a1bf30c2c08f62ba85251",
    "5fbb2f8450f9b019d3d31aa69d7dc1f0220ff9b5cc990cd12f14ba0439ad004c"]
function settings(file)
    s=TOML.parsefile(file)
    expected=Dict("model"=>Dict("qe_version"=>"7.4.1","qe_commit"=>QE_COMMIT,
        "sample_bias_ry"=>-0.0220495933,"spectral_window"=>"sharp_zero_temperature",
        "scf_acceptance_ry"=>5e-5,"control_plot_num"=>5,"candidate_plot_num"=>10),
        "selection"=>Dict("policy"=>"no_mold_or_benchmark_selection",
            "require_archived_control_sha256"=>true,"require_patched_repeat_sha256"=>true,
            "ildos_integral_rtol"=>1e-4,"cube_significant_digits"=>5,"numeric_roundoff_atol"=>1e-13),
        "preprocessing"=>Dict("cube_order"=>"qe_last_axis_fast","clip_negative_values"=>false,
            "half_nm"=>.32,"step_nm"=>.04,"diagnostic_heights_nm"=>[.4,.5,.6]))
    s==expected || error("Not the frozen Gamma reconstruction experiment")
    s
end

"Spacing of a five-significant-digit E13.5 cube value; zero uses roundoff only."
quantum(x,digits) = iszero(x) ? 0.0 : 10.0^(floor(log10(abs(x)))-digits+1)
function quantized_difference(a,b,s)
    length(a)==length(b) || error("Grid length mismatch")
    digits=s["selection"]["cube_significant_digits"]; atol=s["selection"]["numeric_roundoff_atol"]
    violations=0; maxerror=0.0; maxratio=0.0
    for (x,y) in zip(a,b)
        bound=(quantum(x,digits)+quantum(y,digits))/2+atol
        err=abs(x-y); violations+=err>bound
        maxerror=max(maxerror,err); maxratio=max(maxratio,err/bound)
    end
    (;violations,max_abs_error=maxerror,max_error_to_bound=maxratio)
end
function half_pair_relation(stock,patched,s)
    length(stock)==length(patched) || error("Grid length mismatch")
    digits=s["selection"]["cube_significant_digits"]; atol=s["selection"]["numeric_roundoff_atol"]
    mean_patched=mean(patched)
    mean_quantum=sum(quantum(x,digits) for x in patched)/length(patched)
    violations=0; maxerror=0.0; maxratio=0.0
    for (x,y) in zip(stock,patched)
        # Error envelope from rounded stock, patched voxel, and patched mean.
        bound=quantum(x,digits)/2+quantum(y,digits)/4+mean_quantum/4+atol
        err=abs(x-(y+mean_patched)/2); violations+=err>bound
        maxerror=max(maxerror,err); maxratio=max(maxratio,err/bound)
    end
    (;violations,max_abs_error=maxerror,max_error_to_bound=maxratio,mean_patched)
end

function prepare(root,config,outdir)
    s=settings(config)
    grid_source=joinpath(root,"config/unit_assignment_patch_support.toml")
    historical=TOML.parsefile(grid_source)["model"]
    historical["mold_half_nm"]==s["preprocessing"]["half_nm"] || error("Not historical half-width")
    historical["mold_step_nm"]==s["preprocessing"]["step_nm"] || error("Not historical step")
    old=joinpath(root,"qe/spectral_window_20260925")
    xmls=[joinpath(root,"qe/glcn_restart5/qe_tmp/glcn_central.save/data-file-schema.xml"),
        joinpath(old,"work/glcnac/glcnac_central.save/data-file-schema.xml")]
    references=[joinpath(old,"$(m)_reference.cube") for m in ("glcn","glcnac")]
    frames=[joinpath(old,"$(m)_frame.tsv") for m in ("glcn","glcnac")]
    S.sha.(xmls)==XML_SHA || error("Wrong accepted electronic states")
    S.sha.(references)==CUBE_SHA || error("Wrong archived cubes")
    spectra=[S.spectrum(x,s) for x in xmls]
    S.require_plane_domains(references,frames,s) # Before any output or job.
    S.newdir(outdir)
    cp(config,joinpath(outdir,"settings.toml"))
    open(joinpath(outdir,"collected_wfc.sha256"),"w") do io
        for (m,hash) in zip(("glcn","glcnac"),WFC_SHA)
            println(io,hash,"  work/",m,"/",m,"_central.save/wfc1.dat")
        end
    end
    for (k,m) in enumerate(("glcn","glcnac"))
        case=joinpath(outdir,m); mkdir(case)
        sp=spectra[k]; prefix=m*"_central"
        cp(xmls[k],joinpath(case,"data-file-schema.xml"))
        cp(references[k],joinpath(case,"reference.cube"))
        cp(frames[k],joinpath(case,"frame.tsv"))
        S.table(joinpath(case,"bands.tsv"),sp.rows)
        for (name,num) in (("site_control",5),("stock_control",5),("patched_control",5),
                ("stock_ildos",10),("patched_ildos",10),("repeat_ildos",10))
            filename=name=="site_control" ? "pp_ldos.in" : name*".in"
            write(joinpath(case,filename),S.pp_input(prefix,"../work/$m",name,num,sp,s))
        end
        open(joinpath(case,"metadata.toml"),"w") do io
            TOML.print(io,Dict("scope"=>"gamma_reconstruction_only","accepted"=>sp.accepted,
                "xml_sha256"=>S.sha(xmls[k]),"config_sha256"=>S.sha(config),
                "scf_error_ry"=>sp.scf_error_ry,"fermi_ev"=>sp.fermi_ev,
                "emin_ev"=>sp.emin_ev,"emax_ev"=>sp.emax_ev,
                "ildos_bands"=>count(r->r.ildos_weight>0,sp.rows),
                "ildos_integral_expected"=>sum(r.ildos_weight for r in sp.rows)))
        end
    end
    # Provenance-only SCF input for the existing QE launcher's physical checks.
    scf=read(joinpath(root,"qe/glcnac/pw_scf_accept_plain.in"),String)
    write(joinpath(outdir,"pw_scf.in"),replace(scf,"./qe_tmp_accept_plain"=>"./work/glcnac"))
    write(joinpath(outdir,"pp_ldos.in"),S.pp_input("glcnac_central","./work/glcnac","site_control",5,spectra[2],s))
    cp(joinpath(root,"qe/glcnac/pseudo"),joinpath(outdir,"pseudo"))
    cp(joinpath(root,"hpc/qe_gamma_reconstruction.sbatch"),joinpath(outdir,"run_pp_only.sbatch"))
    open(io->TOML.print(io,Dict("scope"=>"pp_only_no_scf_no_benchmark",
        "grid_source_sha256"=>S.sha(grid_source),"xml_sha256"=>XML_SHA,
        "reference_sha256"=>CUBE_SHA,"wfc_sha256"=>WFC_SHA,
        "settings_sha256"=>S.sha(config))),joinpath(outdir,"inputs.toml"),"w")
    println("Both accepted states and all six 17x17 planes validated before preparation.")
end

function compare(case,config)
    s=settings(config); meta=TOML.parsefile(joinpath(case,"metadata.toml"))
    meta["config_sha256"]==S.sha(config) || error("Changed settings")
    meta["xml_sha256"]==S.sha(joinpath(case,"data-file-schema.xml")) || error("Changed XML")
    sp=S.spectrum(joinpath(case,"data-file-schema.xml"),s)
    expected=sum(r.ildos_weight for r in sp.rows)
    expected==meta["ildos_integral_expected"] || error("Wrong state weight")
    names=("site_control","stock_control","patched_control","stock_ildos","patched_ildos","repeat_ildos")
    files=Dict(n=>joinpath(case,n*".cube") for n in names)
    reference=joinpath(case,"reference.cube")
    cubes=Dict(n=>S.QECubeMolds.read_qe_cube(files[n]) for n in names)
    c=cubes["site_control"]
    all(a.dims==c.dims && a.origin_nm==c.origin_nm && a.axes_nm==c.axes_nm for a in values(cubes)) || error("Grid changed")
    stock_control=quantized_difference(c.values,cubes["stock_control"].values,s)
    relation=half_pair_relation(cubes["stock_ildos"].values,cubes["patched_ildos"].values,s)
    stats=[merge((observable=n,),S.cube_stats(cubes[n])) for n in names]
    checks=Dict("site_control_exact"=>S.sha(files["site_control"])==S.sha(reference),
        "rebuilt_control_within_rounding"=>stock_control.violations==0,
        "patch_leaves_stm_exact"=>S.sha(files["stock_control"])==S.sha(files["patched_control"]),
        "patched_repeat_exact"=>S.sha(files["patched_ildos"])==S.sha(files["repeat_ildos"]),
        "half_pair_relation"=>relation.violations==0,
        "ildos_integrals"=>all(abs(r.integral-expected)<=s["selection"]["ildos_integral_rtol"]*expected
            for r in stats if occursin("ildos",r.observable)))
    frame=C.read_frame(joinpath(case,"frame.tsv")); planes=NamedTuple[]; pixels=NamedTuple[]
    for h in s["preprocessing"]["diagnostic_heights_nm"], n in names
        v=S.plane(cubes[n],frame,h,s)
        push!(planes,(observable=n,height_nm=h,values=length(v),negative=count(<(0),v),
            minimum=minimum(v),median=median(v),maximum=maximum(v),
            positive_sum=sum(max(0,x) for x in v),negative_abs_sum=-sum(min(0,x) for x in v)))
        append!(pixels,[(observable=n,height_nm=h,pixel=i,value=v[i]) for i in eachindex(v)])
    end
    out=joinpath(case,"comparison"); S.newdir(out)
    S.table(joinpath(out,"cube_summary.tsv"),stats)
    S.table(joinpath(out,"plane_summary.tsv"),planes)
    S.table(joinpath(out,"planes.tsv"),pixels)
    hashes=Dict(n=>S.sha(files[n]) for n in names)
    summary=Dict("scope"=>"no_calibration_no_mold_no_benchmark","checks"=>checks,
        "expected_integral"=>expected,"cube_sha256"=>hashes,
        "stock_control_difference"=>Dict(string(k)=>v for (k,v) in pairs(stock_control)),
        "half_pair_relation"=>Dict(string(k)=>v for (k,v) in pairs(relation)))
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(basename(case),": ",checks,"; relation ",relation)
    all(values(checks)) || error("Saved Gamma comparison failed; do not promote")
    nothing
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_gamma_reconstruction.jl prepare ROOT CONFIG NEW_DIR | compare CASE CONFIG")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==3 && args[1]=="compare" && return compare(args[2:end]...)
    error("Expected prepare ROOT CONFIG NEW_DIR or compare CASE CONFIG; no label options")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEGammaReconstruction.main()
