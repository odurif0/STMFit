#!/usr/bin/env julia
# Two NSCF eigenvalue tolerances at one accepted density, not a new SCF or mold.
module QEDiagonalPrecision
using TOML, SHA, LinearAlgebra, Statistics
include(joinpath(@__DIR__, "qe_wavefunction_cutoff.jl"))
const B=QEWavefunctionCutoff
const X=B.X; const C=B.C; const D=B.D; const S=B.S; const W=B.W
const MOLECULES=B.MOLECULES
const TRIALS=("lo", "hi")
const STAGES=("reference", TRIALS...)

function settings(file)
    s=TOML.parsefile(file)
    expected=Dict(
        "model"=>Dict("qe_version"=>"7.4.1", "baseline_run"=>"qe/wavefunction_cutoff_20260927",
            "ecutwfc_ry"=>60., "ecutrho_ry"=>720., "functional"=>"PBE", "scf_acceptance_ry"=>5e-5,
            "sample_bias_ry"=>-0.0220495933, "spectral_window"=>"sharp_zero_temperature",
            "diago_thresholds_ry"=>[1e-10,1e-12], "diago_full_acc"=>true,
            "startingpot"=>"file", "startingwfc"=>"file", "nscf_max_seconds"=>3000),
        "selection"=>Dict("policy"=>"two_precision_fixed_density_no_mold_or_benchmark_selection",
            "require_exact_baseline_replay"=>true, "require_exact_parallel_repeat"=>true, "adopt_reference"=>false),
        "preprocessing"=>Dict("half_nm"=>.32, "step_nm"=>.04, "diagnostic_heights_nm"=>[.4,.5,.6],
            "fourier_block_points"=>32, "geometry_rtol"=>1e-12, "atom_position_atol_nm"=>1e-12,
            "clip_density"=>false, "normalize_density"=>false))
    # Preserve the old configuration exactly; the new control declares CG explicitly.
    if haskey(s["model"],"diagonalization") || haskey(s["model"],"diago_cg_maxiter")
        expected["model"]["diagonalization"]="cg"
        expected["model"]["diago_cg_maxiter"]=20
    end
    s==expected || error("Outside the bounded fixed-density precision comparison")
    s
end

iscg(s)=get(s["model"],"diagonalization",nothing)=="cg"

function threshold(s,trial)
    trial in TRIALS || error("Unknown NSCF precision")
    s["model"]["diago_thresholds_ry"][findfirst(==(trial),TRIALS)]
end

function reference_state(xml,s)
    a=C.state(xml,s,s["model"]["ecutrho_ry"])
    S.tag(S.tag(a.inp,"control_variables"),"calculation")=="scf" || error("Reference is not an SCF")
    if iscg(s)
        parse(Int,S.tag(S.tag(a.inp,"electron_control"),"diago_cg_maxiter"))==s["model"]["diago_cg_maxiter"] || error("Changed inherited CG iteration limit")
    end
    a
end

function input_text(xml,m,s,trial)
    reference_state(xml,s)
    q=deepcopy(s); q["model"]["baseline_ecutrho_ry"]=q["model"]["candidate_ecutrho_ry"]=720.
    q["model"]["scf_max_seconds"]=s["model"]["nscf_max_seconds"]
    # Same checked species/geometry writer; only the calculation/solver controls differ.
    text=C.input_text(xml,m,q)
    solver=iscg(s) ? " diagonalization='cg'\n diago_cg_maxiter=$(s["model"]["diago_cg_maxiter"])" : " diagonalization='david'"
    for (before,after) in ((" calculation='scf'", " calculation='nscf'"),
            (" pseudo_dir='../pseudo'", " pseudo_dir='../../pseudo'"),
            (" tprnfor=.true.", " tprnfor=.false."), (" tstress=.true.", " tstress=.false."),
            (" diagonalization='david'", solver*"\n diago_thr_init=$(threshold(s,trial))\n diago_full_acc=.true."))
        count(before,text)==1 || error("Ambiguous NSCF input replacement")
        text=replace(text,before=>after)
    end
    text
end

function verify_inputs(run,m)
    m in MOLECULES || error("Unknown molecule")
    s=settings(joinpath(run,"settings.toml")); dir=joinpath(run,m)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
    meta["molecule"]==m && meta["config_sha256"]==S.sha(joinpath(run,"settings.toml")) || error("Changed case/config")
    all(S.sha(joinpath(dir,f))==sha for (f,sha) in meta["input_sha256"]) || error("Changed prepared input")
    all(S.sha(joinpath(run,"pseudo",f))==sha for (f,sha) in meta["pseudo_sha256"]) || error("Changed PAW dataset")
    s,meta
end

function prepare(root,config,run)
    s=settings(config); old=joinpath(root,s["model"]["baseline_run"])
    prior=B.settings(joinpath(old,"settings.toml"))
    s["preprocessing"]==prior["preprocessing"] || error("Changed queries")
    sources=NamedTuple[]
    for m in MOLECULES
        _,meta=B.verify_inputs(old,m); dir=joinpath(old,m,"60")
        xml=joinpath(dir,"data-file-schema.xml"); reference_state(xml,s)
        report=TOML.parsefile(joinpath(dir,"analysis/summary.toml"))
        all(values(report["checks"])) && report["xml_sha256"]==S.sha(xml) || error("Unverified reference")
        manifest=joinpath(old,m,"converged_before_analysis.sha256")
        files=["data-file-schema.xml","wfc1.dat","charge-density.dat","paw.txt",sort(collect(keys(meta["pseudo_sha256"])))...]
        hashes=Dict(f=>X.manifest_hash(manifest,"work/$(m)_central.save/$f") for f in files)
        hashes["wfc1.dat"]==report["wfc_sha256"] && hashes["data-file-schema.xml"]==S.sha(xml) || error("Changed reference")
        all(hashes[f]==h for (f,h) in meta["pseudo_sha256"]) || error("Checkpoint PAW mismatch")
        push!(sources,(;m,dir,xml,hashes,meta,inputs=[input_text(xml,m,s,t) for t in TRIALS]))
    end
    S.newdir(run); cp(config,joinpath(run,"settings.toml"))
    cp(realpath(joinpath(old,"pseudo")),joinpath(run,"pseudo"))
    for src in sources
        dir=joinpath(run,src.m)
        for stage in STAGES
            mkpath(joinpath(dir,stage))
            cp(joinpath(src.dir,"frame.tsv"),joinpath(dir,stage,"frame.tsv"))
        end
        cp(src.xml,joinpath(dir,"reference/data-file-schema.xml"))
        for file in ("planes.tsv","queries.tsv")
            cp(joinpath(src.dir,"analysis",file),joinpath(dir,"reference_"*file))
        end
        for (trial,input) in zip(TRIALS,src.inputs); write(joinpath(dir,trial,"pw_nscf.in"),input); end
        cp(joinpath(root,"hpc/qe_diagonal_precision.sbatch"),joinpath(dir,"run_precision.sbatch"))
        files=["reference/data-file-schema.xml","reference/frame.tsv","reference_planes.tsv","reference_queries.tsv",
            "lo/frame.tsv","hi/frame.tsv","lo/pw_nscf.in","hi/pw_nscf.in","run_precision.sbatch"]
        meta=Dict("molecule"=>src.m,"config_sha256"=>S.sha(config),"source_sha256"=>src.hashes,
            "pseudo_sha256"=>src.meta["pseudo_sha256"],"input_sha256"=>Dict(f=>S.sha(joinpath(dir,f)) for f in files))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
        open(joinpath(dir,"input_hashes.sha256"),"w") do io
            for f in (files...,"metadata.toml","../settings.toml")
                println(io,S.sha(joinpath(dir,f)),"  ",f)
            end
            for f in sort(readdir(joinpath(run,"pseudo")))
                println(io,S.sha(joinpath(run,"pseudo",f)),"  ../pseudo/",f)
            end
        end
        verify_inputs(run,src.m)
    end
    println("Prepared two fixed-density precisions per geometry; no checkpoint copied, solve or submission")
end

workdir(run,m,trial)=joinpath(run,m,trial,"work",m*"_central.save")

function check_start(run,m,trial)
    s,meta=verify_inputs(run,m); threshold(s,trial)
    work=workdir(run,m,trial)
    all(S.sha(joinpath(work,f))==sha for (f,sha) in meta["source_sha256"]) || error("NSCF must start from its own identical source copy")
    reference_state(joinpath(work,"data-file-schema.xml"),s)
    println(m," ",trial,": identical original density, PAW state and orbitals")
end

function nscf_state(xml,s,old,trial)
    text=read(xml,String); inp=S.tag(text,"input"); out=S.tag(text,"output")
    # QE deliberately marks NSCF scf_conv=false; do not forge an accepted SCF.
    sp=S.spectrum(xml,s;require_accepted=false); geo=D.xml_geometry(xml)
    a=(;text,inp,out,sp,geo); X.same_geometry(old,a,s["preprocessing"])
    cv=S.tag(inp,"control_variables"); ec=S.tag(inp,"electron_control")
    S.tag(cv,"calculation")=="nscf" && S.tag(cv,"restart_mode")=="from_scratch" || error("Not the requested NSCF")
    S.number(S.tag(cv,"max_seconds"))==s["model"]["nscf_max_seconds"] || error("Changed time bound")
    all(S.tag(cv,k)=="false" for k in ("forces","stress")) || error("Unexpected force/stress run")
    # Unlike conv_thr, QE 7.4.1 writes diago_thr_init in input Ry, without /2.
    S.number(S.tag(ec,"diago_thr_init"))==threshold(s,trial) &&
        S.tag(ec,"diago_full_acc")=="true" || error("Changed diagonalization precision")
    expected_solver=iscg(s) ? s["model"]["diagonalization"] : S.tag(S.tag(old.inp,"electron_control"),"diagonalization")
    S.tag(ec,"diagonalization")==expected_solver || error("Changed diagonalization algorithm")
    if iscg(s)
        parse(Int,S.tag(ec,"diago_cg_maxiter"))==s["model"]["diago_cg_maxiter"] || error("Changed CG iteration limit")
    end
    for key in ("conv_thr","mixing_mode","mixing_beta","mixing_ndim","max_nstep")
        S.tag(ec,key)==S.tag(S.tag(old.inp,"electron_control"),key) || error("Changed inherited electronic setting: $key")
    end
    for key in ("basis","dft","bands","spin","k_points_IBZ")
        split(S.tag(inp,key))==split(S.tag(old.inp,key)) || error("Changed physical input: $key")
    end
    geo.dims==old.geo.dims && geo.ecutwfc_ry==old.geo.ecutwfc_ry || error("Changed basis/grid")
    S.tag(S.tag(out,"basis_set"),"ecutrho")==S.tag(S.tag(old.out,"basis_set"),"ecutrho") || error("Changed density cutoff")
    conv=S.tag(out,"convergence_info"); scf=S.tag(conv,"scf_conv")
    S.tag(conv,"wf_collected")=="true" && S.tag(scf,"convergence_achieved")=="false" &&
        parse(Int,S.tag(scf,"n_scf_steps"))==1 || error("Not a collected NSCF checkpoint")
    a
end

function check_nscf(run,m,trial)
    s,meta=verify_inputs(run,m); threshold(s,trial); dir=joinpath(run,m,trial)
    old=reference_state(joinpath(run,m,"reference/data-file-schema.xml"),s)
    work=workdir(run,m,trial); source=joinpath(work,"data-file-schema.xml")
    a=nscf_state(source,s,old,trial)
    for file in ("charge-density.dat","paw.txt",keys(meta["pseudo_sha256"])...)
        S.sha(joinpath(work,file))==meta["source_sha256"][file] || error("NSCF changed the frozen density/PAW state")
    end
    log=read(joinpath(dir,"pw_nscf.out"),String)
    if iscg(s)
        count(r"(?m)^\s*CG style diagonalization\s*$",log)==1 &&
            !occursin(r"Davidson diagonalization|PPCG style diagonalization",log) || error("Missing or unexpected CG solver marker")
    end
    for marker in ("The potential is recalculated from file", "Starting wfcs from file",
            "Band Structure Calculation", "End of band structure calculation", "JOB DONE.")
        count(marker,log)==1 || error("Missing or repeated NSCF marker: $marker")
    end
    occursin(r"(?i)eigenvalues not converged|convergence NOT achieved|Error in routine|Maximum CPU time|iteration #|recomputing them from scratch",log) && error("Unfinished or unexpected NSCF")
    ethr=S.onlymatch(r"ethr\s*=\s*([0-9.EeDd+-]+)",log,"NSCF threshold")
    S.number(ethr.captures[1])==threshold(s,trial) || error("Wrong achieved solver threshold")
    cp(source,joinpath(dir,"data-file-schema.xml");force=false)
    println(m," ",trial,": completed NSCF at ",threshold(s,trial)," Ry; original density/PAW unchanged; not a new SCF")
    nothing
end

function analyze(run,m,stage,wfc_file)
    stage in STAGES || error("Unknown stage")
    s,meta=verify_inputs(run,m); dir=joinpath(run,m,stage); out=joinpath(dir,"analysis")
    (ispath(out)||islink(out)) && error("Output already exists")
    p=s["preprocessing"]; old=reference_state(joinpath(run,m,"reference/data-file-schema.xml"),s)
    xml=joinpath(dir,"data-file-schema.xml")
    a=stage=="reference" ? old : nscf_state(xml,s,old,stage)
    wfc_hash=S.sha(wfc_file)
    expected=stage=="reference" ? meta["source_sha256"]["wfc1.dat"] :
        X.manifest_hash(joinpath(dir,"converged_before_analysis.sha256"),"work/$(m)_central.save/wfc1.dat")
    wfc_hash==expected || error("Changed reference or converged orbitals")
    if stage!="reference"
        S.sha(xml)==S.sha(joinpath(workdir(run,m,stage),"data-file-schema.xml")) || error("Changed NSCF metadata")
        for file in ("charge-density.dat","paw.txt")
            S.sha(joinpath(workdir(run,m,stage),file))==meta["source_sha256"][file] || error("Changed fixed density")
        end
    end
    selected=[r.band for r in a.sp.rows if r.ildos_weight>0]
    weights=[r.ildos_weight for r in a.sp.rows if r.ildos_weight>0]
    w=W.read_wfc(wfc_file,selected)
    w.nbnd==length(a.sp.rows) && isapprox(w.cell,a.geo.cell;rtol=p["geometry_rtol"],atol=0) || error("WFC/XML mismatch")
    maximum(sum(abs2,w.reciprocal*Float64.(w.miller);dims=1))<=s["model"]["ecutwfc_ry"]*(1+p["geometry_rtol"]) || error("WFC exceeds cutoff")
    all(4maximum(abs,view(w.miller,k,:))<a.geo.dims[k] for k in 1:3) || error("Product aliasing")
    miller_sha=bytes2hex(sha256(reinterpret(UInt8,vec(w.miller))))
    if stage!="reference"
        ref=TOML.parsefile(joinpath(run,m,"reference/analysis/summary.toml"))
        w.ng==ref["plane_waves"] && miller_sha==ref["miller_sha256"] || error("Changed ordered reciprocal basis")
    end
    radii,_=D.radii(run,dir); frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
    points=reduce(vcat,[S.plane_points(frame,h,s) for h in p["diagnostic_heights_nm"]])
    ns=D.native_vertices(a.geo.cell,a.geo.dims,points)
    pc=D.P.clearance(points,a.geo.geometry,radii)
    vc=D.P.clearance(transpose(a.geo.cell*ns.vertices)*D.G.C.BOHR_NM,a.geo.geometry,radii)
    min(pc.minimum_clearance_nm,vc.minimum_clearance_nm)>0 || error("Query enters PAW sphere")
    fractions=hcat(ns.vertices,ns.fractions); BLAS.set_num_threads(1)
    direct=W.density(w,fractions,weights;block_points=p["fourier_block_points"])
    repeated=W.density(w,fractions,weights;block_points=p["fourier_block_points"],parallel=false)
    nv=length(ns.indices); products=X.plane_products(direct[nv+1:end],points,p["diagnostic_heights_nm"])
    checks=Dict("outside_paw"=>true,"parallel_serial_exact"=>direct==repeated,
        "finite_nonnegative"=>all(x->isfinite(x)&&x>=0,direct),"unchanged_density_and_basis"=>true)
    if stage=="reference"
        ref=D.V.tsv(joinpath(run,m,"reference_queries.tsv"))
        checks["baseline_queries_exact"]=length(ref)==length(direct) && all(
            parse(Float64,r["density"])==direct[i] &&
            [parse(Float64,r[k]) for k in ("fractional_x","fractional_y","fractional_z")]==collect(fractions[:,i])
            for (i,r) in enumerate(ref))
        refplanes=D.V.tsv(joinpath(run,m,"reference_planes.tsv"))
        checks["baseline_planes_exact"]=length(refplanes)==length(products.rows) && all(
            parse(Float64,r["density"])==v.density && parse(Float64,r["height_nm"])==v.height_nm && parse(Int,r["pixel"])==v.pixel
            for (r,v) in zip(refplanes,products.rows))
    end
    normvalue=W.parseval(w,weights); checks["positive_finite_smooth_norm"]=isfinite(normvalue)&&normvalue>0
    S.newdir(out)
    S.table(joinpath(out,"planes.tsv"),products.rows); S.table(joinpath(out,"plane_summary.tsv"),products.stats)
    S.table(joinpath(out,"decay.tsv"),products.tails); S.table(joinpath(out,"bands.tsv"),a.sp.rows)
    S.table(joinpath(out,"queries.tsv"),[(;query=i,kind=i<=nv ? "native_vertex" : "plane",
        fractional_x=fractions[1,i],fractional_y=fractions[2,i],fractional_z=fractions[3,i],density=direct[i],repeat=repeated[i]) for i in eachindex(direct)])
    report=Dict("checks"=>checks,"molecule"=>m,"stage"=>stage,
        "state_kind"=>stage=="reference" ? "scf_reference" : "nscf",
        "diagonalization"=>S.tag(S.tag(a.inp,"electron_control"),"diagonalization"),
        "diago_thr_init_ry"=>stage=="reference" ? 0. : threshold(s,stage),
        "config_sha256"=>meta["config_sha256"],"xml_sha256"=>S.sha(xml),"wfc_sha256"=>wfc_hash,
        "source_density_sha256"=>meta["source_sha256"]["charge-density.dat"],
        "source_paw_sha256"=>meta["source_sha256"]["paw.txt"],"miller_sha256"=>miller_sha,
        "selected_bands"=>selected,"selected_weights"=>weights,"parseval_smooth_norm"=>normvalue,
        "plane_waves"=>w.ng,"native_grid"=>a.geo.dims,"native_vertices"=>nv,"plane_points"=>size(points,1),
        "source_scf_error_ry"=>old.sp.scf_error_ry,"fermi_ev"=>a.sp.fermi_ev,"emin_ev"=>a.sp.emin_ev,"emax_ev"=>a.sp.emax_ev,
        "point_clearance_nm"=>pc.minimum_clearance_nm,"vertex_clearance_nm"=>vc.minimum_clearance_nm,
        "new_cube_validation"=>false,"new_scf"=>false,"mold_adopted"=>false)
    open(io->TOML.print(io,report),joinpath(out,"summary.toml"),"w")
    println(m," ",stage,": ",checks,"; bands=",selected,"; smooth norm=",normvalue)
    all(values(checks)) || error("Analysis failed; no adjustment or adoption")
    nothing
end

function summarize(run,m)
    s,_=verify_inputs(run,m); dir=joinpath(run,m); rows=NamedTuple[]; stats=NamedTuple[]
    reports=Dict(stage=>TOML.parsefile(joinpath(dir,stage,"analysis/summary.toml")) for stage in STAGES)
    all(all(values(r["checks"])) for r in values(reports)) || error("Unverified stage")
    for r in values(reports), key in ("plane_waves","miller_sha256","source_density_sha256","source_paw_sha256","native_grid")
        r[key]==reports["reference"][key] || error("Unmatched basis/density/grid")
    end
    for (from,to) in (("reference","lo"),("reference","hi"),("lo","hi"))
        x,y=[D.V.tsv(joinpath(dir,stage,"analysis/planes.tsv")) for stage in (from,to)]
        length(x)==length(y) || error("Unpaired queries")
        for (u,v) in zip(x,y)
            all(u[k]==v[k] for k in ("height_nm","pixel","x_nm","y_nm","z_nm")) || error("Unmatched physical point")
            va=parse(Float64,u["density"]); vb=parse(Float64,v["density"])
            push!(rows,(;molecule=m,from,to,height_nm=parse(Float64,u["height_nm"]),pixel=parse(Int,u["pixel"]),
                density_from=va,density_to=vb,difference=vb-va,ratio=X.ratio(va,vb)))
        end
        for h in s["preprocessing"]["diagnostic_heights_nm"]
            pts=filter(r->r.from==from && r.to==to && r.height_nm==h,rows)
            va=[r.density_from for r in pts]; vb=[r.density_to for r in pts]
            push!(stats,(;molecule=m,from,to,height_nm=h,points=length(pts),relative_l2=norm(va)>0 ? norm(vb-va)/norm(va) : NaN,
                minimum_ratio=minimum(r.ratio for r in pts),median_ratio=median([r.ratio for r in pts]),
                maximum_ratio=maximum(r.ratio for r in pts),sum_ratio=X.ratio(sum(va),sum(vb))))
        end
    end
    out=joinpath(dir,"comparison"); S.newdir(out)
    S.table(joinpath(out,"paired_planes.tsv"),rows); S.table(joinpath(out,"summary.tsv"),stats)
    println(m,": ",length(rows)," paired points across three comparisons; primary comparison lo -> hi; no adoption")
    nothing
end

function main(args)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_diagonal_precision.jl prepare ROOT CONFIG NEW_RUN | check-start RUN MOLECULE lo|hi | check-nscf RUN MOLECULE lo|hi | analyze RUN MOLECULE reference|lo|hi WFC | summarize RUN MOLECULE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1]=="check-start" && return check_start(args[2:end]...)
    length(args)==4 && args[1]=="check-nscf" && return check_nscf(args[2:end]...)
    length(args)==5 && args[1]=="analyze" && return analyze(args[2:end]...)
    length(args)==3 && args[1]=="summarize" && return summarize(args[2:end]...)
    error("Use --help")
end

end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEDiagonalPrecision.main(ARGS)
