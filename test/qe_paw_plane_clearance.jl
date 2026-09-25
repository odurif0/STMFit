#!/usr/bin/env julia
# Header-only geometry diagnostic. No density correction, fit, mask or selection.
module QEPAWPlaneClearance
using LinearAlgebra, TOML
include(joinpath(@__DIR__, "qe_gamma_reconstruction.jl"))
const G = QEGammaReconstruction
const S = G.S
const Z = Dict("H"=>1, "C"=>6, "N"=>7, "O"=>8, "Cu"=>29)
attribute(text,key) = S.onlymatch(Regex("\\b"*key*"=\"([^\"]+)\""),text,key).captures[1]

function paw_radius(file)
    text = read(file,String)
    header = S.onlymatch(r"<PP_HEADER\b[^>]*>",text,"UPF header").match
    strip(attribute(header,"is_paw"))=="true" || error("Expected PAW dataset")
    element = strip(attribute(header,"element"))
    radial = S.number.(split(S.tag(text,"PP_R")))
    all(isfinite,radial) && all(>(0),diff(radial)) || error("Invalid radial mesh")
    length(radial)==parse(Int,attribute(header,"mesh_size")) || error("Wrong mesh size")
    aug = S.onlymatch(r"<PP_AUGMENTATION\b[^>]*>",text,"augmentation").match
    indices = [parse(Int,attribute(aug,"cutoff_r_index"))]
    for beta in eachmatch(r"<PP_BETA\.\d+\b[^>]*>",text)
        push!(indices,parse(Int,attribute(beta.match,"cutoff_radius_index")))
    end
    length(indices)>1 && all(i->1<=i<=length(radial),indices) || error("Invalid PAW radii")
    # QE 7.4.1 upflib/read_upf_new.f90 defines the outer integration radius by
    # kkbeta=max(maxval(kbeta),paw%iraug). Do not use the cutoff_r=-1 sentinel.
    outer = maximum(indices)
    (;element,z=Z[element],outer_index=outer,radius_nm=radial[outer]*G.C.BOHR_NM)
end

function cube_atoms(file)
    open(file) do io
        readline(io); readline(io)
        head = split(readline(io)); nat = parse(Int,head[1])
        nat>0 && length(head)==4 || error("Expected scalar QE cube header")
        axes = zeros(3,3); dims = zeros(Int,3)
        for k in 1:3
            line = split(readline(io)); dims[k] = parse(Int,line[1])
            axes[:,k] = S.number.(line[2:4])*G.C.BOHR_NM
        end
        all(>(0),dims) && all(>(0),diag(axes)) || error("Invalid cube axes")
        axes==diagm(diag(axes)) || error("Only the audited orthorhombic cells are supported")
        atoms = map(1:nat) do _
            line = split(readline(io))
            (;z=parse(Int,line[1]),position_nm=S.number.(line[3:5])*G.C.BOHR_NM)
        end
        all(a->all(isfinite,a.position_nm),atoms) || error("Nonfinite atom position")
        (;cell_nm=diag(axes).*dims,atoms)
    end
end

"Minimum distance to any periodically repeated atomic augmentation sphere."
function clearance(points,geometry,radii)
    cell = geometry.cell_nm
    all(>(0),cell) || error("Invalid periodic cell")
    margins = map(eachrow(points)) do point
        minimum(geometry.atoms) do atom
            delta = point-atom.position_nm
            delta -= cell.*round.(delta./cell)
            norm(delta)-radii[atom.z]
        end
    end
    all(isfinite,margins) || error("Nonfinite clearance")
    (;points=length(margins),inside_or_boundary=count(<=(0),margins),minimum_clearance_nm=minimum(margins))
end

function analyze(run_dir,outdir)
    config = joinpath(run_dir,"settings.toml"); settings = G.settings(config)
    rows = NamedTuple[]; radius_rows = NamedTuple[]
    for molecule in ("glcn","glcnac")
        case_dir = joinpath(run_dir,molecule)
        input = S.tag(read(joinpath(case_dir,"data-file-schema.xml"),String),"input")
        radii = Dict{Int,Float64}()
        for species in eachmatch(r"<species\b[^>]*>(.*?)</species>"s,input)
            file = strip(S.tag(species.captures[1],"pseudo_file"))
            basename(file)==file || error("Expected local pseudopotential filename")
            file = joinpath(run_dir,"pseudo",file); r = paw_radius(file)
            haskey(radii,r.z) && error("Duplicate species")
            radii[r.z] = r.radius_nm
            push!(radius_rows,merge((molecule=molecule,pseudo=basename(file),sha256=S.sha(file)),r))
        end
        geometry = cube_atoms(joinpath(case_dir,"reference.cube"))
        frame = G.C.read_frame(joinpath(case_dir,"frame.tsv"))
        for height in settings["preprocessing"]["diagnostic_heights_nm"]
            points = S.plane_points(frame,height,settings)
            result = clearance(points,geometry,radii)
            push!(rows,merge((molecule=molecule,height_nm=height),result))
        end
    end
    S.newdir(outdir)
    S.table(joinpath(outdir,"radii.tsv"),radius_rows)
    S.table(joinpath(outdir,"plane_clearance.tsv"),rows)
    for row in rows; println(row); end
    nothing
end

function main(args=ARGS)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    args==["--help"] && return println("qe_paw_plane_clearance.jl RUN_DIR NEW_OUTPUT_DIR (geometry only)")
    length(args)==2 || error("Expected run directory and new output directory")
    analyze(args...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && QEPAWPlaneClearance.main()
