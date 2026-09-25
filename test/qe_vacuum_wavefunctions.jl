#!/usr/bin/env julia
# Independent WFC -> vacuum density comparison; no mold, fit or label inputs.
module QEVacuumWavefunctions
using LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__,"lib/qe_gamma_wavefunctions.jl"))
include(joinpath(@__DIR__,"qe_paw_plane_clearance.jl"))
include(joinpath(@__DIR__,"verify_qe_spectral_window.jl"))
const W=QEGammaWavefunctions
const P=QEPAWPlaneClearance
const G=P.G
const S=G.S
const V=VerifyQESpectralWindow
const SMOOTH_SHA=["6fbaafb370c65fee6fa3281513683e490b04c6e36de518b94d74c4d94a9df84c",
    "df0c38626142f2a68d1d713b51be93e5e3260b62fb5899aa6730b09ee6344c47"]
const WFC_SHA=["31f137f580bab932c453d8c3e79e9f19e993ce98ed3a1bf30c2c08f62ba85251",
    "5fbb2f8450f9b019d3d31aa69d7dc1f0220ff9b5cc990cd12f14ba0439ad004c"]

function settings(file)
    s=TOML.parsefile(file)
    expected=G.settings(joinpath(@__DIR__,"../config/qe_gamma_reconstruction.toml"))
    for key in ("control_plot_num","candidate_plot_num"); delete!(expected["model"],key); end
    expected["model"]["observable"]="gamma_wavefunction_sum_outside_paw"
    expected["selection"]=Dict("policy"=>"no_mold_or_benchmark_selection",
        "require_reference_sha256"=>true,"require_exact_repeat"=>true,"require_outside_paw"=>true,
        "cube_significant_digits"=>5,"cube_header_decimal_places"=>6,
        "numeric_roundoff_atol"=>1e-13,"geometry_rtol"=>1e-12)
    expected["preprocessing"]["normalize_components"]=false
    expected["preprocessing"]["plane_coordinates"]="accepted_xml_cell"
    expected["preprocessing"]["fourier_block_points"]=32
    s==expected || error("Not the frozen vacuum-wavefunction experiment")
    s
end

function xml_geometry(file)
    out=S.tag(read(file,String),"output")
    structure=S.onlymatch(r"<atomic_structure\b([^>]*)>(.*?)</atomic_structure>"s,out,"structure")
    cell_text=S.tag(structure.captures[2],"cell")
    cell=hcat([S.number.(split(S.tag(cell_text,"a$i"))) for i in 1:3]...)
    cell==diagm(diag(cell)) && all(>(0),diag(cell)) || error("Unsupported cell")
    atoms=[(;z=P.Z[P.attribute(m.captures[1],"name")],
        position_nm=S.number.(split(m.captures[2]))*G.C.BOHR_NM)
        for m in eachmatch(r"<atom\b([^>]*)>(.*?)</atom>"s,S.tag(structure.captures[2],"atomic_positions"))]
    length(atoms)==parse(Int,P.attribute(structure.captures[1],"nat")) || error("Atom count mismatch")
    fft=S.onlymatch(r"<fft_grid\b([^>]*)>",out,"FFT grid").captures[1]
    dims=[parse(Int,P.attribute(fft,"nr$i")) for i in 1:3]
    (;cell,geometry=(;cell_nm=diag(cell)*G.C.BOHR_NM,atoms),dims,
        ecutwfc_ry=2S.number(S.tag(S.tag(out,"basis_set"),"ecutwfc")))
end

function radii(run_dir,case_dir)
    input=S.tag(read(joinpath(case_dir,"data-file-schema.xml"),String),"input")
    rows=NamedTuple[]
    for species in eachmatch(r"<species\b[^>]*>(.*?)</species>"s,input)
        name=strip(S.tag(species.captures[1],"pseudo_file"))
        basename(name)==name || error("Unsafe pseudo path")
        file=joinpath(run_dir,"pseudo",name)
        push!(rows,merge((;pseudo=name,sha256=S.sha(file)),P.paw_radius(file)))
    end
    Dict(r.z=>r.radius_nm for r in rows),rows
end

"All eight interpolation vertices, at precise XML native-grid coordinates."
function native_vertices(cell,dims,points_nm)
    fractions=cell\(transpose(points_nm)/G.C.BOHR_NM)
    indices=NTuple{3,Int}[]
    for f in eachcol(fractions)
        lo=floor.(Int,f.*dims)
        all(0 .<= lo) && all(lo.+1 .< dims) || error("Plane outside native cube")
        append!(indices,[Tuple(lo+[i,j,k]) for i in 0:1 for j in 0:1 for k in 0:1])
    end
    sort!(unique!(indices))
    vertices=hcat([collect(i)./dims for i in indices]...)
    (;indices,vertices,fractions)
end

function prepare(root,config,outdir)
    s=settings(config); previous=joinpath(root,"qe/paw_decomposition_20260925")
    S.sha(joinpath(previous,"settings.toml"))=="438820c143636986e90af005ed0ef6e31975177af1cc91f9c36d0383268a8592" || error("Changed component experiment")
    for (k,m) in enumerate(("glcn","glcnac"))
        dir=joinpath(previous,m)
        S.sha(joinpath(dir,"data-file-schema.xml"))==G.XML_SHA[k] || error("Changed XML")
        S.sha(joinpath(dir,"smooth_ildos.cube"))==SMOOTH_SHA[k] || error("Changed smooth reference")
        S.spectrum(joinpath(dir,"data-file-schema.xml"),s)
        all(values(TOML.parsefile(joinpath(dir,"comparison/summary.toml"))["checks"])) || error("Components not verified")
        geo=xml_geometry(joinpath(dir,"data-file-schema.xml")); frame=G.C.read_frame(joinpath(dir,"frame.tsv"))
        r,_=radii(previous,dir)
        for h in s["preprocessing"]["diagnostic_heights_nm"]
            pts=S.plane_points(frame,h,s); ns=native_vertices(geo.cell,geo.dims,pts)
            P.clearance(pts,geo.geometry,r).inside_or_boundary==0 || error("Plane inside PAW")
            P.clearance(transpose(geo.cell*ns.vertices)*G.C.BOHR_NM,geo.geometry,r).inside_or_boundary==0 || error("Native vertices inside PAW")
        end
    end
    S.newdir(outdir); cp(config,joinpath(outdir,"settings.toml"))
    for name in ("pw_scf.in","pp_ldos.in","pseudo"); cp(joinpath(previous,name),joinpath(outdir,name)); end
    cp(joinpath(root,"hpc/qe_vacuum_wavefunctions.sbatch"),joinpath(outdir,"run_fourier.sbatch"))
    for (k,m) in enumerate(("glcn","glcnac"))
        dir=joinpath(outdir,m); mkdir(dir); src=joinpath(previous,m)
        for name in ("data-file-schema.xml","frame.tsv","bands.tsv")
            cp(joinpath(src,name),joinpath(dir,name))
        end
        cp(joinpath(src,"smooth_ildos.cube"),joinpath(dir,"reference_smooth.cube"))
        _,rs=radii(outdir,dir)
        meta=Dict("scope"=>"vacuum_fourier_diagnostic_no_mold","config_sha256"=>S.sha(config),
            "xml_sha256"=>G.XML_SHA[k],"smooth_sha256"=>SMOOTH_SHA[k],"wfc_sha256"=>WFC_SHA[k],
            "pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in rs))
        open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
    end
    println("Accepted states, reference hashes, all plane domains and native-vertex PAW clearances pass.")
end

function compare(case_dir,config,wfc_file)
    out=joinpath(case_dir,"comparison")
    (ispath(out) || islink(out)) && error("Comparison already exists: $out")
    s=settings(config); p=s["preprocessing"]; q=s["selection"]
    meta=TOML.parsefile(joinpath(case_dir,"metadata.toml"))
    S.sha(config)==meta["config_sha256"] || error("Changed settings")
    xml=joinpath(case_dir,"data-file-schema.xml"); reference=joinpath(case_dir,"reference_smooth.cube")
    S.sha(xml)==meta["xml_sha256"] && S.sha(reference)==meta["smooth_sha256"] &&
        S.sha(wfc_file)==meta["wfc_sha256"] || error("Changed scientific input")
    sp=S.spectrum(xml,s); selected=[r.band for r in sp.rows if r.ildos_weight>0]
    weights=[r.ildos_weight for r in sp.rows if r.ildos_weight>0]
    w=W.read_wfc(wfc_file,selected); geo=xml_geometry(xml)
    length(sp.rows)==w.nbnd || error("WFC/XML band count differs")
    isapprox(geo.cell,w.cell;rtol=q["geometry_rtol"]) || error("WFC/XML cell differs")
    maximum(sum(abs2,w.reciprocal*Float64.(w.miller);dims=1))<=geo.ecutwfc_ry*(1+q["geometry_rtol"]) || error("WFC exceeds cutoff")
    all(4maximum(abs,view(w.miller,k,:))<geo.dims[k] for k in 1:3) || error("Dense grid aliases wavefunction products")
    c=V.cube(reference)
    c.dims==geo.dims && all(iszero,c.origin) || error("Changed native grid")
    exact_axes=geo.cell*Diagonal(1 ./ geo.dims)
    maximum(abs,c.axes-exact_axes)<=.5*10.0^(-q["cube_header_decimal_places"])+q["numeric_roundoff_atol"] || error("Cube header not rounded XML grid")
    r,rs=radii(dirname(case_dir),case_dir)
    Dict(x.pseudo=>x.sha256 for x in rs)==meta["pseudo_sha256"] || error("Changed PAW datasets")
    frame=G.C.read_frame(joinpath(case_dir,"frame.tsv")); heights=p["diagnostic_heights_nm"]
    points=reduce(vcat,[S.plane_points(frame,h,s) for h in heights])
    ns=native_vertices(geo.cell,geo.dims,points)
    point_clearance=P.clearance(points,geo.geometry,r)
    vertex_clearance=P.clearance(transpose(geo.cell*ns.vertices)*G.C.BOHR_NM,geo.geometry,r)
    min(point_clearance.minimum_clearance_nm,vertex_clearance.minimum_clearance_nm)>0 || error("Not a vacuum-only domain")
    BLAS.set_num_threads(1)
    allfractions=hcat(ns.vertices,ns.fractions)
    direct=W.density(w,allfractions,weights;block_points=p["fourier_block_points"])
    repeated=W.density(w,allfractions,weights;block_points=p["fourier_block_points"],parallel=false)
    nv=length(ns.indices); node_rows=NamedTuple[]
    for (i,(x,y,z)) in enumerate(ns.indices)
        observed=c.grid[z+1,y+1,x+1]
        bound=.5G.quantum(observed,q["cube_significant_digits"])+q["numeric_roundoff_atol"]
        push!(node_rows,(;ix=x,iy=y,iz=z,cube=observed,direct=direct[i],repeat=repeated[i],
            abs_error=abs(direct[i]-observed),bound))
    end
    quantum_sum=sum(v->G.quantum(v,q["cube_significant_digits"]),c.raw)
    voxel=det(geo.cell)/prod(geo.dims)
    cube_integral=sum(c.raw)*voxel; norm=W.parseval(w,weights)
    integral_bound=.5voxel*quantum_sum+q["numeric_roundoff_atol"]
    plane_rows=NamedTuple[]; summaries=NamedTuple[]; count_per_height=div(size(points,1),length(heights))
    exact_cube=merge(c,(;axes=exact_axes))
    for (j,h) in enumerate(heights)
        indices=(j-1)*count_per_height+1:j*count_per_height
        dv=direct[nv .+ indices]
        ev=[V.sample(exact_cube,vec(points[i,:])/G.C.BOHR_NM) for i in indices]
        hv=[V.sample(c,vec(points[i,:])/G.C.BOHR_NM) for i in indices]
        for (pixel,i) in enumerate(indices)
            push!(plane_rows,(;height_nm=h,pixel,direct=dv[pixel],exact_grid_interpolated=ev[pixel],
                rounded_header_interpolated=hv[pixel],repeat=repeated[nv+i]))
        end
        push!(summaries,(;height_nm=h,values=length(dv),zeros=count(iszero,dv),
            minimum=minimum(dv),median=median(dv),maximum=maximum(dv),density_sum=sum(dv),
            interpolation_relative_l2=LinearAlgebra.norm(dv-ev)/LinearAlgebra.norm(dv),
            header_rounding_relative_l2=LinearAlgebra.norm(hv-ev)/LinearAlgebra.norm(ev)))
    end
    decay=NamedTuple[]
    for j in 1:length(heights)-1
        first_plane=direct[nv .+ ((j-1)*count_per_height+1:j*count_per_height)]
        next_plane=direct[nv .+ (j*count_per_height+1:(j+1)*count_per_height)]
        for i in eachindex(first_plane)
            a,b=first_plane[i],next_plane[i]
            push!(decay,(;height_from_nm=heights[j],height_to_nm=heights[j+1],pixel=i,
                ratio=b/a,log_decay_per_nm=(a>0 && b>0 ? log(a/b)/(heights[j+1]-heights[j]) : NaN)))
        end
    end
    checks=Dict("native_vertex_reconstruction"=>all(x->x.abs_error<=x.bound,node_rows),
        "parseval_integral"=>abs(norm-cube_integral)<=integral_bound,
        "repeat_exact"=>direct==repeated,"nonnegative"=>all(>=(0),direct),
        "finite"=>all(isfinite,direct),"outside_paw"=>true)
    S.newdir(out)
    S.table(joinpath(out,"native_vertices.tsv"),node_rows)
    S.table(joinpath(out,"planes.tsv"),plane_rows); S.table(joinpath(out,"plane_summary.tsv"),summaries)
    S.table(joinpath(out,"decay.tsv"),decay); S.table(joinpath(out,"radii.tsv"),rs)
    result=Dict("scope"=>"vacuum_wavefunction_verification_no_mold_or_grade","checks"=>checks,
        "selected_bands"=>selected,"selected_weights"=>weights,"plane_points"=>size(points,1),
        "native_vertices"=>nv,"max_native_error_to_bound"=>maximum(x.abs_error/x.bound for x in node_rows),
        "native_violations"=>count(x->x.abs_error>x.bound,node_rows),"parseval_integral"=>norm,
        "cube_integral_exact_cell"=>cube_integral,"cube_integral_rounded_header"=>sum(c.raw)*abs(det(c.axes)),
        "integral_bound"=>integral_bound,"point_clearance_nm"=>point_clearance.minimum_clearance_nm,
        "vertex_clearance_nm"=>vertex_clearance.minimum_clearance_nm,
        "wfc_sha256"=>S.sha(wfc_file),"smooth_sha256"=>S.sha(reference),"config_sha256"=>S.sha(config))
    open(io->TOML.print(io,result),joinpath(out,"summary.toml"),"w")
    println(basename(case_dir),": ",checks,"; native vertices=",nv,"; max error/bound=",result["max_native_error_to_bound"])
    foreach(println,summaries)
    all(values(checks)) || error("Saved independent WFC comparison failed; no density or tolerance adjustment")
    nothing
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_vacuum_wavefunctions.jl prepare ROOT CONFIG NEW_DIR | compare CASE CONFIG WFC_FILE")
    length(args)==4 && args[1]=="prepare" && return prepare(args[2:end]...)
    length(args)==4 && args[1]=="compare" && return compare(args[2:end]...)
    error("Expected prepare or compare; no fitting or label options")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEVacuumWavefunctions.main()
