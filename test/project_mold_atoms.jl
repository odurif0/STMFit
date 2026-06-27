#!/usr/bin/env julia

# Project 3D atom/proxy coordinates into the aligned (t,u) patch frame used by
# connected mold templates.
#
# Input coordinate TSV format:
#   type  atom  element  x_nm  y_nm  z_nm  [weight]  [sigma_t_nm]  [sigma_u_nm]
# where type is 0=GlcN or 1=GlcNAc. Atom names must include the anchors used to
# define the local backbone frame.

using Printf
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

const DEFAULT_COORDS = "templates/chitosan_mold_coords.tsv"
const DEFAULT_OUT = "templates/chitosan_mold_atoms.tsv"

struct Options
    coords::String
    out_tsv::String
    origin_atoms::Vector{String}
    axis_from::String
    axis_to::String
    plane_atom::String
    sigma_t_nm::Float64
    sigma_u_nm::Float64
    z_decay_nm::Float64
end

function _parse_cli(args)
    coords = DEFAULT_COORDS
    out_tsv = DEFAULT_OUT
    origin_atoms = ["C1", "C4"]
    axis_from = "C1"
    axis_to = "C4"
    plane_atom = "C2"
    sigma_t_nm = 0.08
    sigma_u_nm = 0.08
    z_decay_nm = Inf
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--coords"; coords = args[i+1]; i += 2
        elseif startswith(arg, "--coords="); coords = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--origin-atoms"; origin_atoms = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--origin-atoms="); origin_atoms = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--axis-from"; axis_from = args[i+1]; i += 2
        elseif startswith(arg, "--axis-from="); axis_from = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--axis-to"; axis_to = args[i+1]; i += 2
        elseif startswith(arg, "--axis-to="); axis_to = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--plane-atom"; plane_atom = args[i+1]; i += 2
        elseif startswith(arg, "--plane-atom="); plane_atom = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--sigma-t-nm"; sigma_t_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--sigma-t-nm="); sigma_t_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--sigma-u-nm"; sigma_u_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--sigma-u-nm="); sigma_u_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--z-decay-nm"; z_decay_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--z-decay-nm="); z_decay_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/project_mold_atoms.jl [options]

            Options:
              --coords PATH        Input 3D coordinate TSV [$(DEFAULT_COORDS)]
              --out PATH           Output aligned atom/proxy TSV [$(DEFAULT_OUT)]
              --origin-atoms LIST  Comma-separated atoms whose mean is origin [C1,C4]
              --axis-from NAME     Backbone axis start atom [C1]
              --axis-to NAME       Backbone axis end atom [C4]
              --plane-atom NAME    Atom defining positive u side after removing t projection [C2]
              --sigma-t-nm FLOAT   Default projected atom blur along t [0.08]
              --sigma-u-nm FLOAT   Default projected atom blur along u [0.08]
              --z-decay-nm FLOAT   Optional exp(-Δz/z_decay) height weight; Inf disables [Inf]

            Input TSV columns:
              type, atom, element, x_nm, y_nm, z_nm, optional weight/sigma_t_nm/sigma_u_nm

            Output TSV columns:
              type, atom, t_nm, u_nm, weight, sigma_t_nm, sigma_u_nm
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(coords) || error("Coordinate TSV not found: $coords")
    return Options(coords, out_tsv, String.(strip.(origin_atoms)), axis_from, axis_to,
                   plane_atom, sigma_t_nm, sigma_u_nm, z_decay_nm)
end

function _default_element_weight(element::String)
    e = uppercase(strip(element))
    e == "H" && return 0.15
    e == "C" && return 1.0
    e == "N" && return 1.15
    e == "O" && return 1.3
    return 1.0
end

function _load_coords(path)
    _, rows = _read_tsv(path)
    by_type = Dict{Int,Vector{NamedTuple}}(0 => NamedTuple[], 1 => NamedTuple[])
    for row in rows
        for c in ("type", "atom", "element", "x_nm", "y_nm", "z_nm")
            haskey(row, c) || error("Coordinate TSV missing column: $c")
        end
        typ = parse(Int, row["type"])
        typ in (0, 1) || error("type must be 0 or 1")
        element = row["element"]
        weight = haskey(row, "weight") && !isempty(strip(row["weight"])) ? _parse_f(row["weight"]) : _default_element_weight(element)
        push!(by_type[typ], (
            atom=strip(row["atom"]),
            element=strip(element),
            xyz=[_parse_f(row["x_nm"]), _parse_f(row["y_nm"]), _parse_f(row["z_nm"])],
            weight=weight,
            sigt=haskey(row, "sigma_t_nm") ? _parse_f(get(row, "sigma_t_nm", "")) : NaN,
            sigu=haskey(row, "sigma_u_nm") ? _parse_f(get(row, "sigma_u_nm", "")) : NaN,
        ))
    end
    isempty(by_type[0]) && error("No coordinates for type=0")
    isempty(by_type[1]) && error("No coordinates for type=1")
    return by_type
end

function _atom_map(atoms)
    d = Dict{String,Vector{Float64}}()
    for a in atoms
        haskey(d, a.atom) && error("Duplicate atom name $(a.atom) in one type; atom names must be unique")
        d[a.atom] = a.xyz
    end
    return d
end

function _frame(atoms, opt::Options)
    coords = _atom_map(atoms)
    for name in vcat(opt.origin_atoms, [opt.axis_from, opt.axis_to, opt.plane_atom])
        haskey(coords, name) || error("Missing anchor atom $name for type frame")
    end
    origin = mean(reduce(hcat, [coords[name] for name in opt.origin_atoms]); dims=2)[:, 1]
    tvec = coords[opt.axis_to] .- coords[opt.axis_from]
    norm(tvec) > 0 || error("axis-from and axis-to are identical")
    that = tvec ./ norm(tvec)
    uvec = coords[opt.plane_atom] .- origin
    uvec = uvec .- dot(uvec, that) .* that
    norm(uvec) > 0 || error("plane-atom is collinear with axis")
    uhat = uvec ./ norm(uvec)
    return origin, that, uhat
end

function _project_type(atoms, typ::Int, opt::Options)
    origin, that, uhat = _frame(atoms, opt)
    zmax = maximum(a.xyz[3] for a in atoms)
    projected = NamedTuple[]
    for a in atoms
        r = a.xyz .- origin
        w = a.weight
        if isfinite(opt.z_decay_nm) && opt.z_decay_nm > 0
            w *= exp(-(zmax - a.xyz[3]) / opt.z_decay_nm)
        end
        sigt = isfinite(a.sigt) ? a.sigt : opt.sigma_t_nm
        sigu = isfinite(a.sigu) ? a.sigu : opt.sigma_u_nm
        push!(projected, (
            typ=typ,
            atom=a.atom,
            t=dot(r, that),
            u=dot(r, uhat),
            weight=w,
            sigt=sigt,
            sigu=sigu,
        ))
    end
    return projected
end

function main()
    opt = _parse_cli(ARGS)
    by_type = _load_coords(opt.coords)
    mkpath(dirname(opt.out_tsv))
    open(opt.out_tsv, "w") do io
        println(io, join(["type", "atom", "t_nm", "u_nm", "weight", "sigma_t_nm", "sigma_u_nm"], '\t'))
        for typ in (0, 1)
            for a in _project_type(by_type[typ], typ, opt)
                println(io, join([a.typ, a.atom, @sprintf("%.8g", a.t), @sprintf("%.8g", a.u),
                                  @sprintf("%.8g", a.weight), @sprintf("%.8g", a.sigt),
                                  @sprintf("%.8g", a.sigu)], '\t'))
            end
        end
    end
    println("Projected mold atoms")
    println("  coords: ", opt.coords)
    println("  output: ", opt.out_tsv)
    println("  origin: ", join(opt.origin_atoms, ","), "  axis: ", opt.axis_from, "→", opt.axis_to,
            "  +u: ", opt.plane_atom)
end

main()
