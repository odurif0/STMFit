#!/usr/bin/env julia

# Extract the local STMFit mold frame from a slab+trimer XYZ.
#
# Use atom indices of the central unit to define:
#   origin = mean(origin_indices)
#   t_axis = axis_to - axis_from
#   u_axis = plane_index - origin, orthogonalized to t_axis
# Output is in nm and can be passed directly to cube_to_stm_maps.jl with --frame.

using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _parse_ints

struct Atom
    element::String
    xyz::Vector{Float64}
end

struct Options
    xyz::String
    origin_indices::Vector{Int}
    axis_from::Int
    axis_to::Int
    plane_index::Int
    out_tsv::String
    height_nm::Float64
end

function _parse_cli(args)
    xyz = ""
    origin_indices = Int[]
    axis_from = 0
    axis_to = 0
    plane_index = 0
    out_tsv = ""
    height_nm = 0.0
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--xyz"; xyz = args[i+1]; i += 2
        elseif startswith(arg, "--xyz="); xyz = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--origin-indices"; origin_indices = _parse_ints(args[i+1]); i += 2
        elseif startswith(arg, "--origin-indices="); origin_indices = _parse_ints(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--axis-from"; axis_from = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--axis-from="); axis_from = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--axis-to"; axis_to = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--axis-to="); axis_to = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--plane-index"; plane_index = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--plane-index="); plane_index = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--height-nm"; height_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--height-nm="); height_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/extract_qe_mold_frame.jl [options]

            Required:
              --xyz PATH                 Slab+trimer XYZ in angstrom
              --origin-indices I,J,...   Central-lobe atom indices for origin
              --axis-from I              Atom index defining start of backbone axis
              --axis-to J                Atom index defining end of backbone axis
              --plane-index K            Atom index on positive-u side

            Options:
              --out PATH                 Optional TSV output
              --height-nm FLOAT          Include sampling height in printed args [0]

            Output:
              origin_nm, t_axis, u_axis, normal_axis and ready-to-copy
              cube_to_stm_maps.jl argument fragments. Indices are 1-based XYZ
              atom-line indices. Truth labels and sequence are not used.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(xyz) && error("--xyz is required")
    isfile(xyz) || error("XYZ not found: $xyz")
    isempty(origin_indices) && error("--origin-indices is required")
    axis_from > 0 || error("--axis-from is required")
    axis_to > 0 || error("--axis-to is required")
    plane_index > 0 || error("--plane-index is required")
    return Options(xyz, origin_indices, axis_from, axis_to, plane_index, out_tsv, height_nm)
end

function _read_xyz(path::String)
    lines = readlines(path)
    length(lines) >= 2 || error("XYZ too short: $path")
    nat = parse(Int, strip(lines[1]))
    length(lines) >= nat + 2 || error("XYZ has fewer atom lines than declared")
    atoms = Atom[]
    for line in lines[3:(nat+2)]
        parts = split(strip(line))
        length(parts) >= 4 || error("Bad XYZ atom line: $line")
        push!(atoms, Atom(parts[1], [parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])]))
    end
    return atoms
end

function _check_indices(atoms, inds)
    n = length(atoms)
    for i in inds
        1 <= i <= n || error("Atom index $i outside 1:$n")
    end
end

function _fmt_vec(v)
    return join([@sprintf("%.10g", x) for x in v], ',')
end

function main()
    opt = _parse_cli(ARGS)
    atoms = _read_xyz(opt.xyz)
    _check_indices(atoms, vcat(opt.origin_indices, [opt.axis_from, opt.axis_to, opt.plane_index]))
    origin_a = zeros(3)
    for i in opt.origin_indices
        origin_a .+= atoms[i].xyz
    end
    origin_a ./= length(opt.origin_indices)
    t = atoms[opt.axis_to].xyz .- atoms[opt.axis_from].xyz
    norm(t) > 0 || error("axis-from and axis-to are identical")
    that = t ./ norm(t)
    u = atoms[opt.plane_index].xyz .- origin_a
    u = u .- dot(u, that) .* that
    norm(u) > 0 || error("plane-index is collinear with backbone axis")
    uhat = u ./ norm(u)
    nhat = cross(that, uhat)
    nhat ./= norm(nhat)
    origin_nm = origin_a ./ 10.0
    isempty(opt.out_tsv) || _ensure_parent(opt.out_tsv)
    if !isempty(opt.out_tsv)
        open(opt.out_tsv, "w") do io
            println(io, "key\tvalue")
            println(io, "origin_nm\t", _fmt_vec(origin_nm))
            println(io, "t_axis\t", _fmt_vec(that))
            println(io, "u_axis\t", _fmt_vec(uhat))
            println(io, "normal_axis\t", _fmt_vec(nhat))
            println(io, "height_nm\t", @sprintf("%.10g", opt.height_nm))
        end
    end
    println("QE mold frame")
    println("  xyz:       ", opt.xyz)
    println("  origin_nm: ", _fmt_vec(origin_nm))
    println("  t_axis:    ", _fmt_vec(that))
    println("  u_axis:    ", _fmt_vec(uhat))
    println("  normal:    ", _fmt_vec(nhat))
    !isempty(opt.out_tsv) && println("  output:    ", opt.out_tsv)
    println()
    println("cube_to_stm_maps args:")
    !isempty(opt.out_tsv) && println("  --frame TYPE:", opt.out_tsv)
    println("  --origin ", _fmt_vec(origin_nm), " --t-axis ", _fmt_vec(that), " --u-axis ", _fmt_vec(uhat), " --height-nm ", @sprintf("%.10g", opt.height_nm))
end

main()
