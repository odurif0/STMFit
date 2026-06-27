#!/usr/bin/env julia

# Build deterministic initial XYZ geometries for beta-(1->4)-linked chitosan
# trimers used as starting points for QE relaxation. The structures are initial
# coordinates, not optimized geometries.

using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent

const DEFAULT_OUTDIR = "hpc/qe_molds"
const DEFAULT_SPACING = 5.65

struct Atom
    element::String
    label::String
    xyz::Vector{Float64}
end

mutable struct Builder
    atoms::Vector{Atom}
    index::Dict{String,Int}
    bonds::Vector{Tuple{Int,Int}}
    hydrogens::Dict{String,Tuple{Int,Float64}}
end

Builder() = Builder(Atom[], Dict{String,Int}(), Tuple{Int,Int}[], Dict{String,Tuple{Int,Float64}}())

function _add_atom!(b::Builder, element::String, label::String, xyz)
    haskey(b.index, label) && error("Duplicate atom label: $label")
    push!(b.atoms, Atom(element, label, collect(Float64, xyz)))
    b.index[label] = length(b.atoms)
    return length(b.atoms)
end

function _bond!(b::Builder, a::String, c::String)
    push!(b.bonds, (b.index[a], b.index[c]))
end

function _hydrogen_spec!(b::Builder, label::String, count::Int, length::Float64)
    b.hydrogens[label] = (count, length)
end

_unit(v) = norm(v) > 0 ? v ./ norm(v) : [0.0, 0.0, 1.0]

function _basis(axis)
    a = _unit(axis)
    ref = abs(dot(a, [0.0, 0.0, 1.0])) < 0.85 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    e1 = _unit(cross(a, ref))
    e2 = _unit(cross(a, e1))
    return a, e1, e2
end

function _hydrogen_dirs(base, n::Int, phase::Float64=0.0)
    a, e1, e2 = _basis(base)
    if n == 1
        return [a]
    elseif n == 2
        θ = 0.95
        r = cos(phase) .* e1 .+ sin(phase) .* e2
        return [_unit(cos(θ) .* a .+ sin(θ) .* r),
                _unit(cos(θ) .* a .- sin(θ) .* r)]
    elseif n == 3
        θ = 1.23
        return [_unit(cos(θ) .* a .+ sin(θ) .* (cos(φ + phase) .* e1 .+ sin(φ + phase) .* e2)) for φ in (0.0, 2π/3, 4π/3)]
    end
    error("unsupported hydrogen count: $n")
end

function _candidate_score(xyz, dirs, bond_length::Float64, obstacle_xyz)
    hs = [xyz .+ bond_length .* d for d in dirs]
    score = isempty(obstacle_xyz) ? Inf : minimum(norm(h .- o) for h in hs for o in obstacle_xyz)
    if length(hs) > 1
        score = min(score, minimum(norm(hs[i] .- hs[j]) for i in 1:length(hs)-1 for j in i+1:length(hs)))
    end
    return score
end

function _best_hydrogen_dirs(xyz, base, n::Int, bond_length::Float64, obstacle_xyz)
    n == 1 && return _hydrogen_dirs(base, n)
    best_dirs = _hydrogen_dirs(base, n)
    best_score = -Inf
    for phase in range(0, 2π; length=25)[1:end-1]
        dirs = _hydrogen_dirs(base, n, phase)
        score = _candidate_score(xyz, dirs, bond_length, obstacle_xyz)
        if score > best_score
            best_score = score
            best_dirs = dirs
        end
    end
    return best_dirs
end

function _add_residue!(b::Builder; prefix::String, cx::Float64, kind::Symbol, incoming::Bool, outgoing::Bool)
    z = 0.0
    coords = Dict(
        "C1" => [cx + 1.45,  0.00, z + 0.35],
        "C2" => [cx + 0.72,  1.24, z - 0.35],
        "C3" => [cx - 0.72,  1.24, z + 0.35],
        "C4" => [cx - 1.45,  0.00, z - 0.35],
        "C5" => [cx - 0.72, -1.24, z + 0.35],
        "O5" => [cx + 0.72, -1.24, z - 0.35],
        "C6" => [cx - 1.08, -2.52, z + 0.72],
        "O3" => [cx - 0.92,  2.55, z + 0.70],
        "O6" => [cx - 1.18, -3.75, z + 1.05],
        "N2" => [cx + 0.88,  2.55, z - 0.75],
    )
    for name in ("C1", "C2", "C3", "C4", "C5", "O5", "C6", "O3", "O6", "N2")
        element = startswith(name, "C") ? "C" : startswith(name, "O") ? "O" : "N"
        _add_atom!(b, element, "$(prefix)_$(name)", coords[name])
    end
    for (a, c) in (("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C4", "C5"),
                   ("C5", "O5"), ("O5", "C1"), ("C5", "C6"), ("C6", "O6"),
                   ("C3", "O3"), ("C2", "N2"))
        _bond!(b, "$(prefix)_$(a)", "$(prefix)_$(c)")
    end

    if !incoming
        _add_atom!(b, "O", "$(prefix)_O4H", [cx - 2.62, -0.22, z - 0.62])
        _bond!(b, "$(prefix)_C4", "$(prefix)_O4H")
        _hydrogen_spec!(b, "$(prefix)_O4H", 1, 0.98)
    end
    if !outgoing
        _add_atom!(b, "O", "$(prefix)_O1H", [cx + 2.62, -0.22, z + 0.62])
        _bond!(b, "$(prefix)_C1", "$(prefix)_O1H")
        _hydrogen_spec!(b, "$(prefix)_O1H", 1, 0.98)
    end

    if kind == :glcnac
        _add_atom!(b, "C", "$(prefix)_AcC", [cx + 1.02, 3.77, z - 0.58])
        _add_atom!(b, "O", "$(prefix)_AcO", [cx + 0.22, 4.55, z - 0.28])
        _add_atom!(b, "C", "$(prefix)_AcMe", [cx + 2.35, 4.25, z - 0.98])
        _bond!(b, "$(prefix)_N2", "$(prefix)_AcC")
        _bond!(b, "$(prefix)_AcC", "$(prefix)_AcO")
        _bond!(b, "$(prefix)_AcC", "$(prefix)_AcMe")
        _hydrogen_spec!(b, "$(prefix)_N2", 1, 1.02)
        _hydrogen_spec!(b, "$(prefix)_AcMe", 3, 1.09)
    else
        _hydrogen_spec!(b, "$(prefix)_N2", 2, 1.02)
    end

    for name in ("C1", "C2", "C3", "C4", "C5")
        _hydrogen_spec!(b, "$(prefix)_$(name)", 1, 1.09)
    end
    _hydrogen_spec!(b, "$(prefix)_C6", 2, 1.09)
    _hydrogen_spec!(b, "$(prefix)_O3", 1, 0.98)
    _hydrogen_spec!(b, "$(prefix)_O6", 1, 0.98)
    return nothing
end

function _add_glycosidic_link!(b::Builder, left::String, right::String)
    c1 = b.atoms[b.index["$(left)_C1"]].xyz
    c4 = b.atoms[b.index["$(right)_C4"]].xyz
    mid = 0.5 .* (c1 .+ c4) .+ [0.0, -0.08, 0.0]
    label = "$(left)$(right)_Olink"
    _add_atom!(b, "O", label, mid)
    _bond!(b, "$(left)_C1", label)
    _bond!(b, label, "$(right)_C4")
end

function _add_hydrogens!(b::Builder)
    heavy_bonds = copy(b.bonds)
    heavy_xyz = [a.xyz for a in b.atoms]
    for label in sort(collect(keys(b.hydrogens)))
        count, length = b.hydrogens[label]
        idx = b.index[label]
        xyz = b.atoms[idx].xyz
        neigh = Int[]
        for (i, j) in heavy_bonds
            i == idx && push!(neigh, j)
            j == idx && push!(neigh, i)
        end
        base = isempty(neigh) ? [0.0, 0.0, 1.0] : sum((xyz .- b.atoms[n].xyz for n in neigh); init=zeros(3))
        obstacle_xyz = [heavy_xyz[i] for i in eachindex(heavy_xyz) if i != idx]
        dirs = _best_hydrogen_dirs(xyz, base, count, length, obstacle_xyz)
        for k in 1:count
            hlabel = count == 1 ? "H_$(label)" : "H_$(label)_$(k)"
            hidx = _add_atom!(b, "H", hlabel, xyz .+ length .* dirs[k])
            push!(b.bonds, (idx, hidx))
        end
    end
end

function _build_trimer(central_kind::Symbol)
    b = Builder()
    s = DEFAULT_SPACING
    _add_residue!(b; prefix="L", cx=-s, kind=:glcn, incoming=false, outgoing=true)
    _add_residue!(b; prefix="C", cx=0.0, kind=central_kind, incoming=true, outgoing=true)
    _add_residue!(b; prefix="R", cx=s, kind=:glcn, incoming=true, outgoing=false)
    _add_glycosidic_link!(b, "L", "C")
    _add_glycosidic_link!(b, "C", "R")
    _add_hydrogens!(b)
    return b
end

function _bond_set(b::Builder)
    return Set((min(i, j), max(i, j)) for (i, j) in b.bonds)
end

function _distance_report(b::Builder)
    bonded = _bond_set(b)
    min_all = Inf
    min_nonbond = Inf
    min_pair = (0, 0)
    for i in 1:length(b.atoms)-1, j in i+1:length(b.atoms)
        d = norm(b.atoms[i].xyz .- b.atoms[j].xyz)
        min_all = min(min_all, d)
        if !((i, j) in bonded)
            if d < min_nonbond
                min_nonbond = d
                min_pair = (i, j)
            end
        end
    end
    return min_all, min_nonbond, min_pair
end

function _write_xyz(path::String, b::Builder, central_kind::Symbol)
    _ensure_parent(path)
    open(path, "w") do io
        println(io, length(b.atoms))
        println(io, "initial beta-(1->4) chitosan trimer; neighbors=GlcN; central=$(central_kind); coordinates=angstrom; label column is metadata")
        for a in b.atoms
            println(io, @sprintf("%-2s  % .10f  % .10f  % .10f  %s", a.element, a.xyz[1], a.xyz[2], a.xyz[3], a.label))
        end
    end
end

function _write_atoms_tsv(path::String, b::Builder)
    _ensure_parent(path)
    open(path, "w") do io
        println(io, "index\tlabel\telement\tx_ang\ty_ang\tz_ang")
        for (i, a) in enumerate(b.atoms)
            println(io, join([i, a.label, a.element, @sprintf("%.10f", a.xyz[1]), @sprintf("%.10f", a.xyz[2]), @sprintf("%.10f", a.xyz[3])], '\t'))
        end
    end
end

function _join_indices(b::Builder, labels)
    return join([b.index[l] for l in labels], ',')
end

function _write_index_tsv(path::String, b::Builder, xyz_name::String, central_kind::Symbol)
    _ensure_parent(path)
    origin_labels = ["C_C1", "C_C2", "C_C3", "C_C4", "C_C5", "C_O5"]
    origin = _join_indices(b, origin_labels)
    axis_from = b.index["C_C4"]
    axis_to = b.index["C_C1"]
    plane = b.index["C_N2"]
    default_offset = 8 * 8 * 4
    origin_default = join([parse(Int, x) + default_offset for x in split(origin, ',')], ',')
    open(path, "w") do io
        println(io, "key\tvalue")
        println(io, "xyz\t", xyz_name)
        println(io, "central_type\t", central_kind == :glcn ? "0" : "1")
        println(io, "neighbors\tGlcN,GlcN")
        println(io, "n_atoms\t", length(b.atoms))
        println(io, "origin_labels\t", join(origin_labels, ','))
        println(io, "origin_indices_bare\t", origin)
        println(io, "axis_from_label\tC_C4")
        println(io, "axis_from_index_bare\t", axis_from)
        println(io, "axis_to_label\tC_C1")
        println(io, "axis_to_index_bare\t", axis_to)
        println(io, "plane_index_label\tC_N2")
        println(io, "plane_index_bare\t", plane)
        println(io, "frame_command_bare\t--origin-indices $origin --axis-from $axis_from --axis-to $axis_to --plane-index $plane")
        println(io, "default_slab_offset_8x8x4\t", default_offset)
        println(io, "frame_command_after_default_slab\t--origin-indices $origin_default --axis-from $(axis_from + default_offset) --axis-to $(axis_to + default_offset) --plane-index $(plane + default_offset)")
        println(io, "note\tInitial unoptimized geometry; relax before extracting production STM/LDOS maps.")
    end
end

function _parse_cli(args)
    out_dir = DEFAULT_OUTDIR
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--out-dir"
            out_dir = args[i+1]; i += 2
        elseif startswith(arg, "--out-dir=")
            out_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/build_initial_chitosan_trimer_xyz.jl [options]

            Options:
              --out-dir PATH   Output directory [$(DEFAULT_OUTDIR)]

            Writes initial, unoptimized beta-(1->4) trimer XYZ files:
              glcn_central_trimer.xyz
              glcnac_central_trimer.xyz

            Neighbor units are GlcN in both structures. The central unit is GlcN
            or GlcNAc. Companion TSVs record atom labels and frame indices.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return out_dir
end

function main()
    out_dir = _parse_cli(ARGS)
    mkpath(out_dir)
    for (name, kind) in (("glcn", :glcn), ("glcnac", :glcnac))
        b = _build_trimer(kind)
        xyz = joinpath(out_dir, "$(name)_central_trimer.xyz")
        atoms_tsv = joinpath(out_dir, "$(name)_central_trimer_atoms.tsv")
        indices = joinpath(out_dir, "$(name)_central_trimer_indices.tsv")
        _write_xyz(xyz, b, kind)
        _write_atoms_tsv(atoms_tsv, b)
        _write_index_tsv(indices, b, basename(xyz), kind)
        min_all, min_nonbond, pair = _distance_report(b)
        println("Built initial $(name)-central trimer")
        println("  xyz:         ", xyz)
        println("  atoms:       ", length(b.atoms))
        println("  atom table:  ", atoms_tsv)
        println("  indices:     ", indices)
        println("  min distance: ", @sprintf("%.3f Å", min_all))
        println("  min nonbond:  ", @sprintf("%.3f Å", min_nonbond), " between ", b.atoms[pair[1]].label, " and ", b.atoms[pair[2]].label)
    end
end

main()
