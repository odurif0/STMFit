#!/usr/bin/env julia
# Compare the same physical plane queries from existing 360/720-Ry SCFs.
module QECutoffWavefunctions
using TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__, "qe_density_cutoff.jl"))
const C=QEDensityCutoff
const D=C.R.D; const S=D.S; const W=D.W
const MOLECULES=C.R.MOLECULES

function settings(file)
    s=TOML.parsefile(file)
    expected=Dict(
        "model"=>Dict("qe_version"=>"7.4.1", "baseline_run"=>"qe/vacuum_wavefunctions_20260925",
            "candidate_run"=>"qe/density_cutoff_20260926", "reference_settings"=>"config/qe_vacuum_wavefunctions.toml",
            "ecutrho_ry"=>[360.,720.], "ecutwfc_ry"=>50., "functional"=>"PBE",
            "sample_bias_ry"=>-0.0220495933, "spectral_window"=>"sharp_zero_temperature", "scf_acceptance_ry"=>5e-5),
        "selection"=>Dict("policy"=>"paired_saved_orbitals_no_mold_or_benchmark_selection",
            "require_exact_baseline_replay"=>true, "require_exact_parallel_repeat"=>true, "adopt_reference"=>false),
        "preprocessing"=>Dict("half_nm"=>.32, "step_nm"=>.04, "diagnostic_heights_nm"=>[.4,.5,.6],
            "fourier_block_points"=>32, "geometry_rtol"=>1e-12, "atom_position_atol_nm"=>1e-12,
            "clip_density"=>false, "normalize_density"=>false))
    s==expected || error("Outside approved saved-orbital comparison")
    s
end

function same_geometry(a,b,p)
    isapprox(a.geo.cell,b.geo.cell;rtol=p["geometry_rtol"],atol=0) || error("Changed cell")
    aa,bb=a.geo.geometry.atoms,b.geo.geometry.atoms
    length(aa)==length(bb) || error("Changed atom count")
    all(x.z==y.z && isapprox(x.position_nm,y.position_nm;rtol=0,atol=p["atom_position_atol_nm"])
        for (x,y) in zip(aa,bb)) || error("Changed atomic geometry")
    length(a.sp.rows)==length(b.sp.rows) || error("Changed available band count")
    S.number(S.tag(S.tag(a.out,"band_structure"),"nelec"))==
        S.number(S.tag(S.tag(b.out,"band_structure"),"nelec")) || error("Changed electron count")
    C.species_block(a.inp)==C.species_block(b.inp) ||
        split(C.species_block(a.inp))==split(C.species_block(b.inp)) || error("Changed species")
    true
end

function manifest_hash(file,path)
    rows=[split(line;limit=2) for line in eachline(file)]
    only(r[1] for r in rows if strip(r[2])==path)
end

function prepare(root,config,run)
    s=settings(config); model=s["model"]; ref=D.settings(joinpath(root,model["reference_settings"]))
    for key in ("sample_bias_ry","spectral_window","scf_acceptance_ry")
        model[key]==ref["model"][key] || error("Changed observable convention")
    end
    for key in ("half_nm","step_nm","diagnostic_heights_nm","fourier_block_points")
        s["preprocessing"][key]==ref["preprocessing"][key] || error("Changed physical queries")
    end
    old=joinpath(root,model["baseline_run"]); candidate=joinpath(root,model["candidate_run"])
    sources=NamedTuple[]
    for (j,m) in enumerate(MOLECULES)
        olddir=joinpath(old,m); newdir=joinpath(candidate,m)
        oldxml=joinpath(olddir,"data-file-schema.xml"); newxml=joinpath(newdir,"data-file-schema.xml")
        S.sha(oldxml)==D.G.XML_SHA[j] || error("Changed accepted baseline")
        previous=TOML.parsefile(joinpath(olddir,"comparison/summary.toml"))
        all(values(previous["checks"])) && previous["wfc_sha256"]==D.WFC_SHA[j] || error("Unverified baseline")
        report=TOML.parsefile(joinpath(newdir,"analysis/summary.toml"))
        report["exact_total_repeat"] && S.sha(newxml)==report["xml_sha256"] || error("Unverified 720-Ry state")
        _,rs=D.radii(candidate,newdir)
        Dict(r.pseudo=>r.sha256 for r in rs)==TOML.parsefile(joinpath(olddir,"metadata.toml"))["pseudo_sha256"] || error("Changed PAW datasets")
        a=C.state(oldxml,s,360.); b=C.state(newxml,s,720.)
        same_geometry(a,b,s["preprocessing"])
        newsha=manifest_hash(joinpath(newdir,"converged_before_pp.sha256"),"work/$(m)_central.save/wfc1.dat")
        push!(sources,(;m,olddir,oldxml,newxml,wfc_hashes=(D.WFC_SHA[j],newsha)))
    end
    S.newdir(run); cp(config,joinpath(run,"settings.toml"))
    cp(realpath(joinpath(candidate,"pseudo")),joinpath(run,"pseudo"))
    cp(joinpath(root,"hpc/qe_cutoff_wavefunctions.sbatch"),joinpath(run,"run_wavefunctions.sbatch"))
    for src in sources, (i,cutoff) in enumerate((360,720))
        dir=joinpath(run,src.m,string(cutoff)); mkpath(dir)
        cp(realpath(i==1 ? src.oldxml : src.newxml),joinpath(dir,"data-file-schema.xml"))
        cp(joinpath(src.olddir,"frame.tsv"),joinpath(dir,"frame.tsv"))
        cp(joinpath(src.olddir,"comparison/planes.tsv"),joinpath(dir,"reference_planes.tsv"))
        cp(joinpath(src.olddir,"comparison/native_vertices.tsv"),joinpath(dir,"reference_vertices.tsv"))
        geo=D.xml_geometry(joinpath(dir,"data-file-schema.xml")); radii,rs=D.radii(run,dir)
        frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
        points=reduce(vcat,[S.plane_points(frame,h,s) for h in s["preprocessing"]["diagnostic_heights_nm"]])
        ns=D.native_vertices(geo.cell,geo.dims,points)
        D.P.clearance(points,geo.geometry,radii).inside_or_boundary==0 || error("Plane enters PAW sphere")
        D.P.clearance(transpose(geo.cell*ns.vertices)*D.G.C.BOHR_NM,geo.geometry,radii).inside_or_boundary==0 || error("Vertex enters PAW sphere")
        meta=Dict("molecule"=>src.m,"ecutrho_ry"=>cutoff,"wfc_sha256"=>src.wfc_hashes[i],
            "config_sha256"=>S.sha(config),"pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in rs),
            "input_sha256"=>Dict(f=>S.sha(joinpath(dir,f)) for f in
                ("data-file-schema.xml","frame.tsv","reference_planes.tsv","reference_vertices.tsv")))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    end
    open(joinpath(run,"input_hashes.sha256"),"w") do io
        for (dir,_,files) in walkdir(run), file in sort(files)
            file=="input_hashes.sha256" && continue
            path=joinpath(dir,file); println(io,S.sha(path),"  ",relpath(path,run))
        end
    end
    println("Prepared four metadata-only cases; no orbital evaluation or job submission")
end

ratio(a,b)=a>0 ? b/a : NaN
decay(a,b,dz)=a>0 && b>0 ? log(a/b)/dz : NaN

function plane_products(values,points,heights)
    length(values)==size(points,1) && rem(length(values),length(heights))==0 || error("Incomplete planes")
    per=div(length(values),length(heights)); rows=NamedTuple[]; stats=NamedTuple[]; tails=NamedTuple[]
    for (j,h) in enumerate(heights)
        ix=(j-1)*per+1:j*per; v=values[ix]
        for (pixel,k) in enumerate(ix)
            push!(rows,(;height_nm=h,pixel,x_nm=points[k,1],y_nm=points[k,2],z_nm=points[k,3],density=v[pixel]))
        end
        push!(stats,(;height_nm=h,points=per,zeros=count(iszero,v),minimum=minimum(v),median=median(v),maximum=maximum(v),density_sum=sum(v)))
        j==1 && continue
        for pixel in 1:per
            a=values[(j-2)*per+pixel]; b=v[pixel]
            push!(tails,(;height_from_nm=heights[j-1],height_to_nm=h,pixel,lower_density=a,upper_density=b,
                rising=b>a,ratio=ratio(a,b),log_decay_per_nm=decay(a,b,h-heights[j-1])))
        end
    end
    (;rows,stats,tails)
end

function analyze(run,m,cutoff,wfc_file)
    m in MOLECULES && cutoff in (360,720) || error("Unknown case")
    s=settings(joinpath(run,"settings.toml")); p=s["preprocessing"]
    dir=joinpath(run,m,string(cutoff)); out=joinpath(dir,"analysis")
    (ispath(out) || islink(out)) && error("Output already exists")
    meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
    meta["molecule"]==m && meta["ecutrho_ry"]==cutoff || error("Wrong case metadata")
    S.sha(joinpath(run,"settings.toml"))==meta["config_sha256"] || error("Changed settings")
    all(S.sha(joinpath(dir,f))==sha for (f,sha) in meta["input_sha256"]) || error("Changed case input")
    S.sha(wfc_file)==meta["wfc_sha256"] || error("Changed collected orbitals")
    a=C.state(joinpath(dir,"data-file-schema.xml"),s,Float64(cutoff)); sp=a.sp; geo=a.geo
    selected=[r.band for r in sp.rows if r.ildos_weight>0]; weights=[r.ildos_weight for r in sp.rows if r.ildos_weight>0]
    w=W.read_wfc(wfc_file,selected)
    w.nbnd==length(sp.rows) && isapprox(w.cell,geo.cell;rtol=p["geometry_rtol"],atol=0) || error("WFC/XML mismatch")
    maximum(sum(abs2,w.reciprocal*Float64.(w.miller);dims=1))<=geo.ecutwfc_ry*(1+p["geometry_rtol"]) || error("WFC exceeds cutoff")
    all(4maximum(abs,view(w.miller,k,:))<geo.dims[k] for k in 1:3) || error("Product aliasing")
    radii,rs=D.radii(run,dir)
    Dict(r.pseudo=>r.sha256 for r in rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
    frame=D.G.C.read_frame(joinpath(dir,"frame.tsv")); heights=p["diagnostic_heights_nm"]
    points=reduce(vcat,[S.plane_points(frame,h,s) for h in heights])
    ns=D.native_vertices(geo.cell,geo.dims,points)
    pc=D.P.clearance(points,geo.geometry,radii)
    vc=D.P.clearance(transpose(geo.cell*ns.vertices)*D.G.C.BOHR_NM,geo.geometry,radii)
    min(pc.minimum_clearance_nm,vc.minimum_clearance_nm)>0 || error("Query inside PAW sphere")
    # Preserve the validated query ordering, including the native-vertex prefix.
    # At 720 Ry these vertices have no independently exported cube reference.
    fractions=hcat(ns.vertices,ns.fractions); BLAS.set_num_threads(1)
    direct=W.density(w,fractions,weights;block_points=p["fourier_block_points"])
    repeated=W.density(w,fractions,weights;block_points=p["fourier_block_points"],parallel=false)
    nv=length(ns.indices); plane_values=direct[nv+1:end]
    products=plane_products(plane_values,points,heights)
    checks=Dict("finite_nonnegative"=>all(x->isfinite(x)&&x>=0,direct),
        "parallel_serial_exact"=>direct==repeated,"outside_paw"=>true)
    if cutoff==360
        old=D.V.tsv(joinpath(dir,"reference_planes.tsv")); vertices=D.V.tsv(joinpath(dir,"reference_vertices.tsv"))
        checks["baseline_planes_exact"]=length(old)==length(products.rows) && all(
            parse(Float64,r["height_nm"])==v.height_nm && parse(Int,r["pixel"])==v.pixel &&
            parse(Float64,r["direct"])==v.density for (r,v) in zip(old,products.rows))
        checks["baseline_vertices_exact"]=length(vertices)==nv && all(
            Tuple(parse(Int,r[k]) for k in ("ix","iy","iz"))==ns.indices[i] && parse(Float64,r["direct"])==direct[i]
            for (i,r) in enumerate(vertices))
    end
    norm=W.parseval(w,weights); checks["positive_finite_smooth_norm"]=isfinite(norm)&&norm>0
    S.newdir(out)
    S.table(joinpath(out,"planes.tsv"),products.rows); S.table(joinpath(out,"plane_summary.tsv"),products.stats)
    S.table(joinpath(out,"decay.tsv"),products.tails); S.table(joinpath(out,"bands.tsv"),sp.rows)
    S.table(joinpath(out,"queries.tsv"),[(;query=i,kind=i<=nv ? "native_vertex" : "plane",
        fractional_x=fractions[1,i],fractional_y=fractions[2,i],fractional_z=fractions[3,i],density=direct[i],repeat=repeated[i]) for i in eachindex(direct)])
    summary=Dict("checks"=>checks,"molecule"=>m,"ecutrho_ry"=>cutoff,"config_sha256"=>meta["config_sha256"],
        "xml_sha256"=>meta["input_sha256"]["data-file-schema.xml"],"wfc_sha256"=>meta["wfc_sha256"],
        "scf_error_ry"=>sp.scf_error_ry,"fermi_ev"=>sp.fermi_ev,"emin_ev"=>sp.emin_ev,"emax_ev"=>sp.emax_ev,
        "selected_bands"=>selected,"selected_weights"=>weights,"parseval_smooth_norm"=>norm,
        "native_vertices"=>nv,"plane_points"=>length(plane_values),"point_clearance_nm"=>pc.minimum_clearance_nm,
        "vertex_clearance_nm"=>vc.minimum_clearance_nm,"new_cube_validation"=>false,"mold_adopted"=>false)
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    println(m," ",cutoff," Ry: ",checks,"; selected bands=",selected,"; smooth norm=",norm)
    all(values(checks)) || error("Saved-orbital check failed; no tolerance or state substitution")
    nothing
end

function summarize(run)
    settings(joinpath(run,"settings.toml")); out=joinpath(run,"comparison"); S.newdir(out)
    rows=NamedTuple[]; stats=NamedTuple[]
    for m in MOLECULES
        dirs=[joinpath(run,m,string(c),"analysis") for c in (360,720)]
        for dir in dirs
            all(values(TOML.parsefile(joinpath(dir,"summary.toml"))["checks"])) || error("Unverified case")
        end
        a,b=[D.V.tsv(joinpath(dir,"planes.tsv")) for dir in dirs]
        length(a)==length(b) || error("Unpaired points")
        localrows=NamedTuple[]
        for (x,y) in zip(a,b)
            all(x[k]==y[k] for k in ("height_nm","pixel","x_nm","y_nm","z_nm")) || error("Unmatched physical queries")
            va=parse(Float64,x["density"]); vb=parse(Float64,y["density"])
            push!(localrows,(;molecule=m,height_nm=parse(Float64,x["height_nm"]),pixel=parse(Int,x["pixel"]),
                density_360=va,density_720=vb,difference=vb-va,ratio=ratio(va,vb)))
        end
        append!(rows,localrows)
        for h in sort(unique(r.height_nm for r in localrows))
            pts=filter(r->r.height_nm==h,localrows); va=[r.density_360 for r in pts]; vb=[r.density_720 for r in pts]
            push!(stats,(;molecule=m,height_nm=h,points=length(pts),relative_l2=norm(vb-va)/norm(va),
                minimum_ratio=minimum(r.ratio for r in pts),median_ratio=median([r.ratio for r in pts]),
                maximum_ratio=maximum(r.ratio for r in pts),sum_ratio=ratio(sum(va),sum(vb))))
        end
    end
    S.table(joinpath(out,"paired_planes.tsv"),rows); S.table(joinpath(out,"summary.tsv"),stats)
    println("Compared ",length(rows)," unchanged physical points; no normalization or selection")
end

function main(args)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_cutoff_wavefunctions.jl prepare ROOT CONFIG NEW_RUN | analyze RUN MOLECULE CUTOFF WFC | summarize RUN")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==5 && args[1]=="analyze" && return analyze(args[2],args[3],parse(Int,args[4]),args[5])
    length(args)==2 && args[1]=="summarize" && return summarize(args[2])
    error("Use --help")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QECutoffWavefunctions.main(ARGS)
