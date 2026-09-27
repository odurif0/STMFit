#!/usr/bin/env julia
# 50/60-Ry orbital-basis sensitivity at fixed 720-Ry density cutoff.
module QEWavefunctionCutoff
using TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__,"qe_cutoff_wavefunctions.jl"))
const X=QECutoffWavefunctions
const C=X.C; const D=X.D; const S=X.S; const W=X.W
const MOLECULES=X.MOLECULES

function settings(file)
    s=TOML.parsefile(file)
    expected=Dict(
        "model"=>Dict("qe_version"=>"7.4.1", "baseline_run"=>"qe/cutoff_wavefunctions_20260927",
            "baseline_scf_run"=>"qe/density_cutoff_20260926", "ecutwfc_ry"=>[50.,60.],
            "ecutrho_ry"=>720., "functional"=>"PBE", "scf_acceptance_ry"=>5e-5,
            "sample_bias_ry"=>-0.0220495933, "spectral_window"=>"sharp_zero_temperature",
            "startingpot"=>"atomic", "startingwfc"=>"atomic+random", "scf_max_seconds"=>12600),
        "selection"=>Dict("policy"=>"two_point_wavefunction_sensitivity_no_mold_or_benchmark_selection",
            "require_exact_baseline_replay"=>true, "require_exact_parallel_repeat"=>true, "adopt_reference"=>false),
        "preprocessing"=>Dict("half_nm"=>.32, "step_nm"=>.04, "diagnostic_heights_nm"=>[.4,.5,.6],
            "fourier_block_points"=>32, "geometry_rtol"=>1e-12, "atom_position_atol_nm"=>1e-12,
            "clip_density"=>false, "normalize_density"=>false))
    s==expected || error("Outside approved wavefunction-cutoff scope")
    s
end

function state(xml,s,cutoff)
    cutoff in s["model"]["ecutwfc_ry"] || error("Unexpected wavefunction cutoff")
    q=deepcopy(s); q["model"]["ecutwfc_ry"]=Float64(cutoff)
    a=C.state(xml,q,s["model"]["ecutrho_ry"])
    basis=S.tag(a.inp,"basis")
    2S.number(S.tag(basis,"ecutwfc"))==cutoff &&
        2S.number(S.tag(basis,"ecutrho"))==s["model"]["ecutrho_ry"] || error("Inconsistent input cutoffs")
    a
end

function input_text(xml,m,s)
    state(xml,s,50)
    # Reuse the checked geometry/species/electronic-setting writer, unchanged.
    q=deepcopy(s); q["model"]["ecutwfc_ry"]=50.
    q["model"]["baseline_ecutrho_ry"]=q["model"]["candidate_ecutrho_ry"]=720.
    text=C.input_text(xml,m,q)
    count(" ecutwfc=50.0\n",text)==1 || error("Ambiguous input cutoff")
    replace(text," ecutwfc=50.0\n"=>" ecutwfc=$(s["model"]["ecutwfc_ry"][2])\n")
end

function verify_inputs(run,m)
    m in MOLECULES || error("Unknown molecule")
    s=settings(joinpath(run,"settings.toml")); dir=joinpath(run,m)
    meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
    meta["molecule"]==m && meta["config_sha256"]==S.sha(joinpath(run,"settings.toml")) || error("Changed settings/case")
    all(S.sha(joinpath(dir,f))==sha for (f,sha) in meta["input_sha256"]) || error("Changed case input")
    all(S.sha(joinpath(run,"pseudo",f))==sha for (f,sha) in meta["pseudo_sha256"]) || error("Changed PAW dataset")
    s,meta
end

function prepare(root,config,run)
    s=settings(config); old=joinpath(root,s["model"]["baseline_run"])
    prior=X.settings(joinpath(old,"settings.toml"))
    s["preprocessing"]==prior["preprocessing"] || error("Changed physical queries")
    sources=NamedTuple[]
    for m in MOLECULES
        dir=joinpath(old,m,"720"); xml=joinpath(dir,"data-file-schema.xml")
        report=TOML.parsefile(joinpath(dir,"analysis/summary.toml"))
        meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
        all(values(report["checks"])) && S.sha(xml)==report["xml_sha256"] || error("Unverified baseline")
        report["wfc_sha256"]==meta["wfc_sha256"] || error("Unverified baseline orbitals")
        state(xml,s,50)
        _,rs=D.radii(old,dir)
        Dict(r.pseudo=>r.sha256 for r in rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
        push!(sources,(;m,dir,xml,meta,input=input_text(xml,m,s)))
    end
    S.newdir(run); cp(config,joinpath(run,"settings.toml"))
    cp(realpath(joinpath(old,"pseudo")),joinpath(run,"pseudo"))
    for src in sources
        dir=joinpath(run,src.m); mkpath(joinpath(dir,"50")); mkpath(joinpath(dir,"60"))
        cp(src.xml,joinpath(dir,"50/data-file-schema.xml"))
        for cutoff in (50,60); cp(joinpath(src.dir,"frame.tsv"),joinpath(dir,string(cutoff),"frame.tsv")); end
        for file in ("planes.tsv","queries.tsv")
            cp(joinpath(src.dir,"analysis",file),joinpath(dir,"reference_"*file))
        end
        write(joinpath(dir,"pw_scf.in"),src.input)
        cp(joinpath(root,"hpc/qe_wavefunction_cutoff.sbatch"),joinpath(dir,"run_basis.sbatch"))
        files=("50/data-file-schema.xml","50/frame.tsv","60/frame.tsv","reference_planes.tsv",
            "reference_queries.tsv","pw_scf.in","run_basis.sbatch")
        meta=Dict("molecule"=>src.m,"config_sha256"=>S.sha(config),
            "baseline_wfc_sha256"=>src.meta["wfc_sha256"],"pseudo_sha256"=>src.meta["pseudo_sha256"],
            "input_sha256"=>Dict(f=>S.sha(joinpath(dir,f)) for f in files))
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
    println("Prepared two new SCFs; no state copied into work and no job submitted")
end

function check_scf(run,m)
    s,_=verify_inputs(run,m); dir=joinpath(run,m)
    old=state(joinpath(dir,"50/data-file-schema.xml"),s,50)
    source=joinpath(dir,"work",m*"_central.save/data-file-schema.xml")
    a=state(source,s,60); X.same_geometry(old,a,s["preprocessing"])
    for key in ("conv_thr","mixing_mode","mixing_beta","mixing_ndim","max_nstep","diagonalization")
        S.tag(S.tag(old.inp,"electron_control"),key)==S.tag(S.tag(a.inp,"electron_control"),key) || error("Changed electronic setting: $key")
    end
    log=read(joinpath(dir,"pw_scf.out"),String)
    occursin("convergence has been achieved",log) && occursin("JOB DONE.",log) || error("Unfinished SCF")
    occursin("Initial potential from superposition of free atoms",log) || error("Wrong starting potential")
    occursin(r"Starting wfcs are [^\n]*random",log) && !occursin("Starting wfcs from file",log) || error("Wrong wavefunction initialization")
    cp(source,joinpath(dir,"60/data-file-schema.xml"))
    println(m,": accepted 60/720-Ry state; SCF error=",a.sp.scf_error_ry," Ry; native grid=",a.geo.dims)
end

function analyze(run,m,cutoff,wfc_file)
    s,meta=verify_inputs(run,m); dir=joinpath(run,m,string(cutoff)); out=joinpath(dir,"analysis")
    (ispath(out)||islink(out)) && error("Output already exists")
    a=state(joinpath(dir,"data-file-schema.xml"),s,cutoff); p=s["preprocessing"]
    old=state(joinpath(run,m,"50/data-file-schema.xml"),s,50)
    X.same_geometry(old,a,p)
    wfc_hash=S.sha(wfc_file)
    if cutoff==50
        wfc_hash==meta["baseline_wfc_sha256"] || error("Changed baseline orbitals")
    else
        source=joinpath(run,m,"work",m*"_central.save/data-file-schema.xml")
        S.sha(source)==S.sha(joinpath(dir,"data-file-schema.xml")) || error("Changed accepted state")
        expected=X.manifest_hash(joinpath(run,m,"converged_before_analysis.sha256"),"work/$(m)_central.save/wfc1.dat")
        wfc_hash==expected || error("Changed converged orbitals")
    end
    selected=[r.band for r in a.sp.rows if r.ildos_weight>0]; weights=[r.ildos_weight for r in a.sp.rows if r.ildos_weight>0]
    w=W.read_wfc(wfc_file,selected)
    w.nbnd==length(a.sp.rows) && isapprox(w.cell,a.geo.cell;rtol=p["geometry_rtol"],atol=0) || error("WFC/XML mismatch")
    maximum(sum(abs2,w.reciprocal*Float64.(w.miller);dims=1))<=cutoff*(1+p["geometry_rtol"]) || error("WFC exceeds cutoff")
    all(4maximum(abs,view(w.miller,k,:))<a.geo.dims[k] for k in 1:3) || error("Product aliasing")
    radii,_=D.radii(run,dir); frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
    points=reduce(vcat,[S.plane_points(frame,h,s) for h in p["diagnostic_heights_nm"]])
    ns=D.native_vertices(a.geo.cell,a.geo.dims,points)
    pc=D.P.clearance(points,a.geo.geometry,radii)
    vc=D.P.clearance(transpose(a.geo.cell*ns.vertices)*D.G.C.BOHR_NM,a.geo.geometry,radii)
    min(pc.minimum_clearance_nm,vc.minimum_clearance_nm)>0 || error("Query inside PAW sphere")
    fractions=hcat(ns.vertices,ns.fractions); BLAS.set_num_threads(1)
    direct=W.density(w,fractions,weights;block_points=p["fourier_block_points"])
    repeat=W.density(w,fractions,weights;block_points=p["fourier_block_points"],parallel=false)
    nv=length(ns.indices); products=X.plane_products(direct[nv+1:end],points,p["diagnostic_heights_nm"])
    checks=Dict("parallel_serial_exact"=>direct==repeat,"finite_nonnegative"=>all(x->isfinite(x)&&x>=0,direct),"outside_paw"=>true)
    if cutoff==50
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
        fractional_x=fractions[1,i],fractional_y=fractions[2,i],fractional_z=fractions[3,i],density=direct[i],repeat=repeat[i]) for i in eachindex(direct)])
    report=Dict("checks"=>checks,"molecule"=>m,"ecutwfc_ry"=>cutoff,"ecutrho_ry"=>s["model"]["ecutrho_ry"],
        "config_sha256"=>meta["config_sha256"],"xml_sha256"=>S.sha(joinpath(dir,"data-file-schema.xml")),"wfc_sha256"=>wfc_hash,
        "selected_bands"=>selected,"selected_weights"=>weights,"parseval_smooth_norm"=>normvalue,
        "plane_waves"=>w.ng,"native_grid"=>a.geo.dims,"native_vertices"=>nv,"plane_points"=>size(points,1),
        "scf_error_ry"=>a.sp.scf_error_ry,"fermi_ev"=>a.sp.fermi_ev,"emin_ev"=>a.sp.emin_ev,"emax_ev"=>a.sp.emax_ev,
        "point_clearance_nm"=>pc.minimum_clearance_nm,"vertex_clearance_nm"=>vc.minimum_clearance_nm,
        "new_cube_validation"=>false,"mold_adopted"=>false)
    open(io->TOML.print(io,report),joinpath(out,"summary.toml"),"w")
    println(m," ",cutoff," Ry: ",checks,"; plane waves=",w.ng,"; bands=",selected,"; smooth norm=",normvalue)
    all(values(checks)) || error("Analysis failed; no tolerance adjustment")
    nothing
end

function summarize(run,m)
    s,_=verify_inputs(run,m); dir=joinpath(run,m); out=joinpath(dir,"comparison")
    a,b=[TOML.parsefile(joinpath(dir,string(c),"analysis/summary.toml")) for c in (50,60)]
    all(values(a["checks"])) && all(values(b["checks"])) || error("Unverified analysis")
    b["plane_waves"]>a["plane_waves"] || error("Wavefunction basis did not expand")
    x,y=[D.V.tsv(joinpath(dir,string(c),"analysis/planes.tsv")) for c in (50,60)]
    length(x)==length(y) || error("Unpaired physical queries")
    rows=NamedTuple[]; stats=NamedTuple[]
    for (u,v) in zip(x,y)
        all(u[k]==v[k] for k in ("height_nm","pixel","x_nm","y_nm","z_nm")) || error("Unmatched physical point")
        va=parse(Float64,u["density"]); vb=parse(Float64,v["density"])
        push!(rows,(;molecule=m,height_nm=parse(Float64,u["height_nm"]),pixel=parse(Int,u["pixel"]),
            density_50=va,density_60=vb,difference=vb-va,ratio=X.ratio(va,vb)))
    end
    for h in s["preprocessing"]["diagnostic_heights_nm"]
        pts=filter(r->r.height_nm==h,rows); va=[r.density_50 for r in pts]; vb=[r.density_60 for r in pts]
        push!(stats,(;molecule=m,height_nm=h,points=length(pts),relative_l2=norm(va)>0 ? norm(vb-va)/norm(va) : NaN,
            minimum_ratio=minimum(r.ratio for r in pts),median_ratio=median([r.ratio for r in pts]),
            maximum_ratio=maximum(r.ratio for r in pts),sum_ratio=X.ratio(sum(va),sum(vb))))
    end
    S.newdir(out); S.table(joinpath(out,"paired_planes.tsv"),rows); S.table(joinpath(out,"summary.tsv"),stats)
    println(m,": ",length(rows)," paired physical points; expanded basis; no normalization or adoption")
end

function main(args)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_wavefunction_cutoff.jl prepare ROOT CONFIG NEW_RUN | check-scf RUN MOLECULE | analyze RUN MOLECULE 50|60 WFC | summarize RUN MOLECULE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==3 && args[1]=="check-scf" && return check_scf(args[2:end]...)
    length(args)==5 && args[1]=="analyze" && return analyze(args[2],args[3],parse(Int,args[4]),args[5])
    length(args)==3 && args[1]=="summarize" && return summarize(args[2:end]...)
    error("Use --help")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEWavefunctionCutoff.main(ARGS)
