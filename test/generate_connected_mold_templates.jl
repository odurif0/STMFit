#!/usr/bin/env julia

# Generate connected mold-template TSVs from aligned geometric/proxy sites.
#
# Input site TSV format:
#   type  atom  t_nm  u_nm  weight  sigma_t_nm  sigma_u_nm
# where type is 0=GlcN or 1=GlcNAc. Coordinates are in the same aligned patch
# frame used by extract_lobe_patches.jl: t along the chain axis and u transverse.
#
# The script expands each base type into all unary connectivity variants:
#   type ∈ {0,1}, parity ∈ {0,1}, mirror ∈ {0,1}
# by applying parity and mirror flips to the coordinates, then writes the wide
# template TSV expected by score_connected_mold_templates.jl.
# If --bond-out is supplied, it also writes sliding pair templates:
#   left_type/right_type ∈ {00,01,10,11}, parity ∈ {0,1}, mirror ∈ {0,1}
# using concatenated left/right unary mold vectors (l_pNNN/r_pNNN columns).

using Printf
using Statistics

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

const DEFAULT_ATOMS = "templates/chitosan_geometric_sites.tsv"
const DEFAULT_OUT = "templates/chitosan_connected_molds.tsv"

struct Options
    atoms::String
    out_tsv::String
    bond_out_tsv::String
    half_nm::Float64
    step_nm::Float64
    parity_flip::String
    mirror_flip::String
    normalize::String
end

function _parse_cli(args)
    atoms = DEFAULT_ATOMS
    out_tsv = DEFAULT_OUT
    bond_out_tsv = ""
    half_nm = 0.32
    step_nm = 0.08
    parity_flip = "t"
    mirror_flip = "u"
    normalize = "zscore"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--atoms"; atoms = args[i+1]; i += 2
        elseif startswith(arg, "--atoms="); atoms = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--bond-out"; bond_out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--bond-out="); bond_out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--half-nm"; half_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--half-nm="); half_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--step-nm"; step_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--step-nm="); step_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--parity-flip"; parity_flip = args[i+1]; i += 2
        elseif startswith(arg, "--parity-flip="); parity_flip = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--mirror-flip"; mirror_flip = args[i+1]; i += 2
        elseif startswith(arg, "--mirror-flip="); mirror_flip = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--normalize"; normalize = args[i+1]; i += 2
        elseif startswith(arg, "--normalize="); normalize = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/generate_connected_mold_templates.jl [options]

            Options:
              --atoms PATH        Geometric/proxy site TSV [$(DEFAULT_ATOMS)]
              --out PATH          Output connected mold TSV [$(DEFAULT_OUT)]
              --bond-out PATH     Optional output sliding bond template TSV
              --half-nm FLOAT     Template half-size, must match patch extraction [0.32]
              --step-nm FLOAT     Template grid spacing, must match patch extraction [0.08]
              --parity-flip STR   none | t | u | both [t]
              --mirror-flip STR   none | t | u | both [u]
              --normalize STR     none | sum | max | zscore [zscore]

            Input TSV columns:
              type, atom, t_nm, u_nm, weight, sigma_t_nm, sigma_u_nm

            Output TSV columns:
              name, type, parity, mirror, p001, p002, ...

            Bond output TSV columns:
              name, left_type, right_type, parity, mirror, l_p001, ..., r_p001, ...
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(atoms) || error("Atom TSV not found: $atoms")
    parity_flip in ("none", "t", "u", "both") || error("Invalid --parity-flip")
    mirror_flip in ("none", "t", "u", "both") || error("Invalid --mirror-flip")
    normalize in ("none", "sum", "max", "zscore") || error("Invalid --normalize")
    return Options(atoms, out_tsv, bond_out_tsv, half_nm, step_nm, parity_flip, mirror_flip, normalize)
end

function _apply_flip(t, u, spec::String)
    spec == "none" && return t, u
    spec == "t" && return -t, u
    spec == "u" && return t, -u
    return -t, -u
end

function _normalize(v::Vector{Float64}, method::String)
    method == "none" && return v
    if method == "sum"
        s = sum(abs, v)
        return s > 0 ? v ./ s : v
    elseif method == "max"
        m = maximum(abs.(v))
        return m > 0 ? v ./ m : v
    end
    σ = std(v)
    σ = σ > 0 ? σ : 1.0
    return (v .- mean(v)) ./ σ
end

function _load_atoms(path)
    _, rows = _read_tsv(path)
    atoms = Dict{Int,Vector{NamedTuple}}(0 => NamedTuple[], 1 => NamedTuple[])
    required = ("type", "t_nm", "u_nm", "weight", "sigma_t_nm", "sigma_u_nm")
    for row in rows
        all(haskey(row, k) for k in required) || error("Atom TSV missing one of: $(join(required, ", "))")
        typ = parse(Int, row["type"])
        typ in (0, 1) || error("type must be 0 or 1")
        atom = get(row, "atom", "atom")
        push!(atoms[typ], (
            atom=atom,
            t=_parse_f(row["t_nm"]),
            u=_parse_f(row["u_nm"]),
            weight=_parse_f(row["weight"]),
            sigt=max(_parse_f(row["sigma_t_nm"]), eps(Float64)),
            sigu=max(_parse_f(row["sigma_u_nm"]), eps(Float64)),
        ))
    end
    isempty(atoms[0]) && error("No atoms for type=0 (GlcN)")
    isempty(atoms[1]) && error("No atoms for type=1 (GlcNAc)")
    return atoms
end

function _template_values(atoms, coords, parity::Int, mirror::Int, opt::Options)
    vals = Float64[]
    transformed = []
    for a in atoms
        t, u = a.t, a.u
        parity == 1 && ((t, u) = _apply_flip(t, u, opt.parity_flip))
        mirror == 1 && ((t, u) = _apply_flip(t, u, opt.mirror_flip))
        push!(transformed, (t=t, u=u, weight=a.weight, sigt=a.sigt, sigu=a.sigu))
    end
    for u in coords, t in coords
        v = 0.0
        for a in transformed
            v += a.weight * exp(-0.5 * (((t - a.t) / a.sigt)^2 + ((u - a.u) / a.sigu)^2))
        end
        push!(vals, v)
    end
    return _normalize(vals, opt.normalize)
end

function main()
    opt = _parse_cli(ARGS)
    atoms = _load_atoms(opt.atoms)
    coords = collect(-opt.half_nm:opt.step_nm:opt.half_nm)
    pix = [@sprintf("p%03d", i) for i in 1:(length(coords)^2)]
    mkpath(dirname(opt.out_tsv))
    open(opt.out_tsv, "w") do io
        println(io, join(vcat(["name", "type", "parity", "mirror"], pix), '\t'))
        for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
            vals = _template_values(atoms[typ], coords, parity, mirror, opt)
            name = @sprintf("%s_p%d_m%d", typ == 0 ? "GlcN" : "GlcNAc", parity, mirror)
            println(io, join(vcat([name, string(typ), string(parity), string(mirror)],
                                  [@sprintf("%.8g", v) for v in vals]), '\t'))
        end
    end
    if !isempty(opt.bond_out_tsv)
        mkpath(dirname(opt.bond_out_tsv))
        open(opt.bond_out_tsv, "w") do io
            left_pix = ["l_$(p)" for p in pix]
            right_pix = ["r_$(p)" for p in pix]
            println(io, join(vcat(["name", "left_type", "right_type", "parity", "mirror"], left_pix, right_pix), '\t'))
            for lt in (0, 1), rt in (0, 1), parity in (0, 1), mirror in (0, 1)
                left_vals = _template_values(atoms[lt], coords, parity, mirror, opt)
                # The right lobe has the opposite ring parity along a β-(1→4) chain.
                right_vals = _template_values(atoms[rt], coords, 1 - parity, mirror, opt)
                name = @sprintf("%d%d_p%d_m%d", lt, rt, parity, mirror)
                println(io, join(vcat([name, string(lt), string(rt), string(parity), string(mirror)],
                                      [@sprintf("%.8g", v) for v in left_vals],
                                      [@sprintf("%.8g", v) for v in right_vals]), '\t'))
            end
        end
    end
    println("Generated connected mold templates")
    println("  sites:     ", opt.atoms)
    println("  output:    ", opt.out_tsv)
    !isempty(opt.bond_out_tsv) && println("  bond out:  ", opt.bond_out_tsv)
    println("  grid:      ", length(coords), "x", length(coords), " (", length(pix), " pixels)")
    println("  parity:    ", opt.parity_flip)
    println("  mirror:    ", opt.mirror_flip)
    println("  normalize: ", opt.normalize)
end

main()
