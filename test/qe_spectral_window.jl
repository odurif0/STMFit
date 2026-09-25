#!/usr/bin/env julia
# A QE observable diagnostic, not a molecular predictor or mold exporter.
module QESpectralWindow
using SHA, TOML, Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "lib", "qe_cube_molds.jl"))
using .QECubeMolds
const C = QECubeMolds.CCMoldNative
const RY_EV = 13.605693122994 # QE 7.4.1 Modules/constants.f90
const HARTREE_EV = 2RY_EV
sha(path) = open(io -> bytes2hex(sha256(io)), path)

function settings(path)
    s = TOML.parsefile(path)
    Set(keys(s)) == Set(["model", "selection", "preprocessing"]) || error("Unexpected sections")
    m, q, p = s["model"], s["selection"], s["preprocessing"]
    Set(keys(m)) == Set(["qe_version", "sample_bias_ry", "spectral_window", "scf_acceptance_ry",
        "control_plot_num", "candidate_plot_num"]) || error("Unexpected model keys")
    m["qe_version"] == "7.4.1" && m["spectral_window"] == "sharp_zero_temperature" || error("Unsupported QE convention")
    m["control_plot_num"] == 5 && m["candidate_plot_num"] == 10 || error("Unexpected observables")
    m["sample_bias_ry"] == -0.0220495933 && m["scf_acceptance_ry"] == 5e-5 || error("Changed bias/acceptance")
    q == Dict("policy"=>"no_mold_or_benchmark_selection", "require_archived_control_sha256"=>true,
        "ildos_integral_rtol"=>1e-4) || error("Unexpected selection policy")
    Set(keys(p)) == Set(["cube_order", "clip_negative_values", "half_nm", "step_nm", "diagnostic_heights_nm"]) || error("Unexpected preprocessing")
    p["cube_order"] == "qe_last_axis_fast" && p["clip_negative_values"] === false || error("Invalid cube policy")
    p["half_nm"] == 0.64 && p["step_nm"] == 0.08 && p["diagnostic_heights_nm"] == [0.4,0.5,0.6] || error("Changed diagnostic grid")
    s
end

function onlymatch(rx, text, description)
    matches = collect(eachmatch(rx, text))
    length(matches) == 1 || error("Expected exactly one $description, found $(length(matches))")
    only(matches)
end
function tag(text, name)
    strip(onlymatch(Regex("<"*name*">(.*?)</"*name*">", "s"), text, name).captures[1])
end
number(x) = parse(Float64, replace(x, 'D'=>'E', 'd'=>'e'))

"Read only the audited single-Gamma, scalar QE 7.4.1 XML schema; energies are Hartree."
function spectrum(xml, s; require_accepted=true)
    text = read(xml, String)
    occursin("Units=\"Hartree atomic units\"",text) || error("Unknown XML energy units")
    occursin(r"<creator\b[^>]*VERSION=\"7\.4\.1\"", text) || error("Not QE 7.4.1 schema")
    out = tag(text, "output")
    tag(tag(out,"basis_set"),"gamma_only") == "true" || error("Not the Gamma-only algorithm")
    bands = tag(out, "band_structure")
    parse(Int, tag(bands, "nks")) == 1 || error("Only one Gamma point is supported")
    all(tag(bands,k) == "false" for k in ("lsda","noncolin","spinorbit")) || error("Only scalar nonmagnetic spectra supported")
    conv = tag(out,"convergence_info")
    scf = tag(conv,"scf_conv")
    accepted = tag(scf,"convergence_achieved") == "true" && tag(conv,"wf_collected") == "true"
    scf_error_ry = 2number(tag(scf,"scf_error"))
    accepted &= isfinite(scf_error_ry) && 0 <= scf_error_ry <= s["model"]["scf_acceptance_ry"]
    require_accepted && !accepted && error("Checkpoint is not accepted with collected wavefunctions")
    eig = tag(bands,"ks_energies")
    kp = onlymatch(r"<k_point\s+weight=\"([^\"]+)\">(.*?)</k_point>"s, eig, "Gamma point")
    all(iszero, number.(split(kp.captures[2]))) || error("Not Gamma")
    weight = number(kp.captures[1]); weight == 2.0 || error("Unexpected spin/k weight")
    ev = onlymatch(r"<eigenvalues\s+size=\"(\d+)\">(.*?)</eigenvalues>"s, eig, "eigenvalues")
    energies_ha = number.(split(ev.captures[2]))
    length(energies_ha) == parse(Int,ev.captures[1]) == parse(Int,tag(bands,"nbnd")) || error("Band count mismatch")
    all(isfinite,energies_ha) && issorted(energies_ha) || error("Invalid eigenvalues")
    sm = onlymatch(r"<smearing\s+degauss=\"([^\"]+)\">(.*?)</smearing>"s, bands, "smearing")
    strip(sm.captures[2]) == "mv" || error("Only the recorded cold smearing is supported")
    width_ha = number(sm.captures[1]); isfinite(width_ha) && width_ha > 0 || error("Invalid smearing")
    ef_ha = number(tag(bands,"fermi_energy")); isfinite(ef_ha) || error("Invalid Fermi energy")
    bias_ha = s["model"]["sample_bias_ry"]/2
    lo, hi = minmax(ef_ha,ef_ha+bias_ha)
    # Exactly the bands included by stm.f90; do not invent a tail beyond its 3*degauss cutoff.
    first_band = findlast(<(lo-3width_ha),energies_ha)
    last_band = findlast(<(hi+3width_ha),energies_ha)
    first_band === nothing && error("Legacy first_band would be uninitialized")
    last_band === nothing && error("No legacy band window")
    first_band += 1
    first_band <= last_band || error("Empty legacy band interval")
    rows = NamedTuple[]
    for (i,e) in enumerate(energies_ha)
        e in (lo,hi) && error("Legacy exact-boundary weight is not assigned by stm.f90")
        w = lo < e < hi ? weight : weight*cold_derivative(((e<lo ? lo : hi)-e)/width_ha)
        used = first_band <= i <= last_band
        push!(rows,(band=i,energy_ev=e*HARTREE_EV,relative_to_fermi_ev=(e-ef_ha)*HARTREE_EV,
            legacy_included=used,legacy_weight=used ? w : 0.0,
            ildos_weight=lo <= e <= hi ? weight : 0.0))
    end
    sum(r.ildos_weight for r in rows) > 0 || error("No states in physical bias window")
    (;rows,accepted,scf_error_ry,fermi_ev=ef_ha*HARTREE_EV,emin_ev=lo*HARTREE_EV,
        emax_ev=hi*HARTREE_EV,smearing_ev=width_ha*HARTREE_EV,first_band,last_band,
        energy_ry=2number(tag(tag(out,"total_energy"),"etot")))
end

cold_derivative(x) = exp(-min(200.0,(x-1/sqrt(2.0))^2))*(2-sqrt(2.0)*x)/sqrt(pi)

function newdir(path)
    (ispath(path) || islink(path)) && error("Output already exists: $path")
    mkpath(path)
end
function table(path, rows)
    isempty(rows) && error("Empty table")
    open(path,"w") do io
        println(io,join(string.(keys(first(rows))),'\t'))
        for row in rows
            println(io,join([x isa AbstractFloat ? @sprintf("%.17g",x) : string(x) for x in values(row)],'\t'))
        end
    end
end
function pp_input(prefix, checkpoint, stem, plot_num, sp, s)
    occursin(r"^[A-Za-z0-9_]+$",prefix) && occursin(r"^[A-Za-z0-9_]+$",stem) || error("Unsafe QE identifier")
    occursin(r"['\r\n]",checkpoint) && error("Unsupported checkpoint path")
    spectral = plot_num == 5 ? @sprintf("  sample_bias = %.17g\n",s["model"]["sample_bias_ry"]) :
        @sprintf("  emin = %.17g\n  emax = %.17g\n",sp.emin_ev,sp.emax_ev)
    plot_num in (5,10) || error("Unsupported PP observable")
    "&INPUTPP\n  prefix = '$prefix'\n  outdir = '$checkpoint'\n  filplot = '$stem.dat'\n  plot_num = $plot_num\n"*
        spectral*"/\n&PLOT\n  iflag = 3\n  output_format = 6\n  fileout = '$stem.cube'\n/\n"
end

function prepare_spectrum(xml, config, prefix, checkpoint, outdir)
    s = settings(config); sp = spectrum(xml,s)
    input_sha = sha(xml)
    newdir(outdir)
    table(joinpath(outdir,"bands.tsv"),sp.rows)
    for (stem,num) in (("control",5),("ildos",10))
        write(joinpath(outdir,"$stem.in"),pp_input(prefix,abspath(checkpoint),stem,num,sp,s))
    end
    metadata=Dict{String,Any}("scope"=>"spectral_diagnostic_only", "xml_sha256"=>input_sha,
        "config_sha256"=>sha(config),"accepted"=>sp.accepted,"scf_error_ry"=>sp.scf_error_ry,
        "energy_ry"=>sp.energy_ry,"fermi_ev"=>sp.fermi_ev,"emin_ev"=>sp.emin_ev,"emax_ev"=>sp.emax_ev,
        "smearing_ev"=>sp.smearing_ev,"legacy_first_band"=>sp.first_band,"legacy_last_band"=>sp.last_band,
        "legacy_negative_bands"=>count(r->r.legacy_weight<0,sp.rows),
        "legacy_negative_weight"=>sum(min(0,r.legacy_weight) for r in sp.rows),
        "legacy_positive_weight"=>sum(max(0,r.legacy_weight) for r in sp.rows),
        "ildos_bands"=>count(r->r.ildos_weight>0,sp.rows),
        "ildos_integral_expected"=>sum(r.ildos_weight for r in sp.rows))
    open(io->TOML.print(io,metadata),joinpath(outdir,"metadata.toml"),"w")
    cp(xml,joinpath(outdir,"data-file-schema.xml"))
    cp(config,joinpath(outdir,"settings.toml"))
    println(prefix,": ",metadata)
end

function cube_stats(c)
    all(isfinite,c.values) || error("Nonfinite cube")
    voxel_bohr3 = abs(det(c.axes_nm))/C.BOHR_NM^3
    (;values=length(c.values),negative=count(<(0),c.values),minimum=minimum(c.values),
        maximum=maximum(c.values),integral=sum(c.values)*voxel_bohr3,
        negative_abs_integral=-sum(x for x in c.values if x<0;init=0.0)*voxel_bohr3)
end
function plane_points(frame,height,s)
    p=s["preprocessing"]; grid=collect(-p["half_nm"]:p["step_nm"]:p["half_nm"])
    normal=C._normal(frame)
    reduce(vcat,[permutedims(frame.origin_nm+t*frame.t_axis+u*frame.u_axis+height*normal)
        for u in grid for t in grid])
end
function outside_count(origin,axes,dims,points)
    inverse=inv(axes)
    count(eachrow(points)) do row
        q=inverse*(row-origin)
        !all(0 <= q[k] <= dims[k]-1 for k in 1:3)
    end
end
function plane(c,frame,height,s)
    points=plane_points(frame,height,s)
    outside_count(c.origin_nm,c.axes_nm,c.dims,points)==0 || error("Diagnostic plane outside cube")
    C.sample_volume(c,points)
end

"Check reference geometry using only its six header lines, before any QE run."
function plane_domain(cubefile,framefile,s)
    header=open(cubefile) do io
        [split(readline(io)) for _ in 1:6]
    end
    origin=number.(header[3][2:4])*C.BOHR_NM
    axes=hcat([number.(header[k][2:4]) for k in 4:6]...)*C.BOHR_NM
    dims=[parse(Int,header[k][1]) for k in 4:6]
    all(>(0),dims) || error("Invalid QE reference dimensions")
    frame=C.read_frame(framefile)
    [(;height_nm=h,outside=outside_count(origin,axes,dims,plane_points(frame,h,s)))
        for h in s["preprocessing"]["diagnostic_heights_nm"]]
end
function require_plane_domains(cubes,frames,s)
    length(cubes)==length(frames) || error("Missing reference frame")
    failures=String[]
    for (cube,frame) in zip(cubes,frames), row in plane_domain(cube,frame,s)
        row.outside==0 || push!(failures,"$(basename(cube)): h=$(row.height_nm), $(row.outside) points outside")
    end
    isempty(failures) || error("Plane-domain preflight failed before preparation: "*join(failures,"; "))
    nothing
end
function compare(control,candidate,reference,framefile,metadatafile,config,outdir)
    s=settings(config); meta=TOML.parsefile(metadatafile)
    sha(control)==sha(reference) || error("Legacy control does not reproduce archived cube")
    meta["config_sha256"]==sha(config) || error("Changed spectral settings")
    c=QECubeMolds.read_qe_cube(control); a=QECubeMolds.read_qe_cube(candidate)
    c.dims==a.dims && c.origin_nm==a.origin_nm && c.axes_nm==a.axes_nm || error("Changed cube grid")
    cs,as=cube_stats(c),cube_stats(a)
    expected=meta["ildos_integral_expected"]
    abs(as.integral-expected)<=s["selection"]["ildos_integral_rtol"]*expected || error("ILDOS integral does not match selected state weight")
    frame=C.read_frame(framefile)
    rows=NamedTuple[]; samples=NamedTuple[]
    for h in s["preprocessing"]["diagnostic_heights_nm"]
        x,y=plane(c,frame,h,s),plane(a,frame,h,s)
        all(isfinite,x) && all(isfinite,y) || error("Unavailable plane")
        for (name,v) in (("control",x),("ildos",y))
            push!(rows,(observable=name,height_nm=h,negative=count(<(0),v),minimum=minimum(v),
                median=median(v),maximum=maximum(v),positive_sum=sum(max(0,z) for z in v),
                negative_abs_sum=-sum(min(0,z) for z in v)))
        end
        append!(samples,[(height_nm=h,pixel=i,control=x[i],ildos=y[i]) for i in eachindex(x)])
    end
    newdir(outdir)
    table(joinpath(outdir,"cube_summary.tsv"),[merge((observable="control",),cs),merge((observable="ildos",),as)])
    table(joinpath(outdir,"plane_summary.tsv"),rows); table(joinpath(outdir,"planes.tsv"),samples)
    inputs=Dict(k=>Dict("path"=>abspath(v),"sha256"=>sha(v)) for (k,v) in
        ("control"=>control,"candidate"=>candidate,"reference"=>reference,"frame"=>framefile,"metadata"=>metadatafile,"config"=>config))
    open(io->TOML.print(io,Dict("inputs"=>inputs,"scope"=>"no_current_calibration_no_mold_no_benchmark",
        "expected_integral"=>expected,"measured_integral"=>as.integral)),joinpath(outdir,"comparison.toml"),"w")
    println("Archived control exact; ILDOS integral ",as.integral," vs ",expected,
        "; native negative values: ",cs.negative," -> ",as.negative)
end

function prepare_run(root,config,outdir)
    s=settings(config)
    glcn=joinpath(root,"qe","glcn_restart5"); glcnac=joinpath(root,"qe","glcnac")
    scf=joinpath(glcnac,"pw_scf_accept_plain.in")
    input=read(scf,String)
    count("outdir = './qe_tmp_accept_plain'",input)==1 || error("Unexpected SCF source")
    input=replace(input,"outdir = './qe_tmp_accept_plain'"=>"outdir = './work/glcnac'")
    occursin("conv_thr = 5.0d-5",input) || error("Changed SCF acceptance")
    references=[joinpath(glcn,"glcn_central_ldos.cube"),joinpath(glcnac,"glcnac_central_ldos_accept_plain.cube")]
    expected=["80cd1d1fde94cf084cc7ea464d2bf065b36b2c015ea0bfaef8cebfee8ff88863",
        "40649ccd9eb6768444b8ff61bf4a639b3940cf926fe3eb42254eb024faf9b5bf"]
    sha.(references)==expected || error("Wrong archived reference cubes")
    require_plane_domains(references,[joinpath(glcn,"frame.tsv"),joinpath(glcnac,"frame.tsv")],s)
    newdir(outdir); mkdir(joinpath(outdir,"pseudo"))
    write(joinpath(outdir,"pw_scf.in"),input)
    cp(joinpath(root,"hpc","qe_spectral_window.sbatch"),joinpath(outdir,"run_scf_pp.sbatch"))
    # Existing launcher preflight keeps checking the unchanged plot_num=5 control.
    pp=read(joinpath(glcnac,"pp_ldos_accept_plain.in"),String)
    write(joinpath(outdir,"pp_ldos.in"),replace(pp,"./qe_tmp_accept_plain"=>"./work/glcnac"))
    for file in readdir(joinpath(glcnac,"pseudo"))
        endswith(file,".UPF") && cp(joinpath(glcnac,"pseudo",file),joinpath(outdir,"pseudo",file))
    end
    for (name,dir,reference) in zip(("glcn","glcnac"),(glcn,glcnac),references)
        cp(reference,joinpath(outdir,"$(name)_reference.cube"))
        cp(joinpath(dir,"frame.tsv"),joinpath(outdir,"$(name)_frame.tsv"))
    end
    cp(config,joinpath(outdir,"settings.toml"))
    open(io->TOML.print(io,Dict("accepted_scf_input_sha256"=>sha(scf),"settings_sha256"=>sha(config),
        "reference_sha256"=>expected,"scope"=>"one_hour_restore_control_then_ildos")),joinpath(outdir,"inputs.toml"),"w")
    println("Prepared new run directory ",outdir)
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    if isempty(args) || args==["--help"]
        println("qe_spectral_window.jl prepare-run --root R --settings S --outdir NEW\n",
            "  prepare-spectrum --xml X --settings S --prefix P --checkpoint D --outdir NEW\n",
            "  compare --control C --candidate C --reference C --frame F --metadata M --settings S --outdir NEW")
        return
    end
    mode=first(args); opts=Dict{String,String}()
    isodd(length(args)) || error("Expected key/value options")
    for i in 2:2:length(args)
        key=args[i]; startswith(key,"--") && !haskey(opts,key) || error("Invalid/repeated option")
        opts[key]=args[i+1]
    end
    names=mode=="prepare-run" ? ["root","settings","outdir"] : mode=="prepare-spectrum" ?
        ["xml","settings","prefix","checkpoint","outdir"] : mode=="compare" ?
        ["control","candidate","reference","frame","metadata","settings","outdir"] : error("Unknown mode")
    Set(keys(opts))==Set("--".*names) || error("Missing or forbidden options")
    fn=mode=="prepare-run" ? prepare_run : mode=="prepare-spectrum" ? prepare_spectrum : compare
    fn((opts["--"*k] for k in names)...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QESpectralWindow.main()
