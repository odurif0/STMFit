#!/usr/bin/env julia

# Validate the initial chitosan trimer and slab XYZ files used as QE starting
# geometries. This is a structural audit only; it does not use STM truth labels.

using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _parse_ints, _read_key_tsv

const DEFAULT_DIR = "hpc/qe_molds"
const DEFAULT_OUT = "hpc/qe_molds/structure_validation.tsv"

struct Atom
    element::String
    xyz::Vector{Float64}
    label::String
end

function _read_xyz(path::String)
    lines = readlines(path)
    length(lines) >= 2 || error("XYZ too short: $path")
    n = parse(Int, strip(lines[1]))
    length(lines) >= n + 2 || error("XYZ has fewer atom lines than declared: $path")
    atoms = Atom[]
    for (k, line) in enumerate(lines[3:(n+2)])
        parts = split(strip(line))
        length(parts) >= 4 || error("Bad XYZ atom line $k in $path: $line")
        label = length(parts) >= 5 ? parts[5] : "$(parts[1])$k"
        push!(atoms, Atom(parts[1], [parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])], label))
    end
    return atoms
end

function _formula(atoms)
    counts = Dict{String,Int}()
    for a in atoms
        counts[a.element] = get(counts, a.element, 0) + 1
    end
    elems = sort(collect(keys(counts)); by=e -> e == "C" ? "0" : e == "H" ? "1" : e)
    return join(["$e$(counts[e])" for e in elems], " ")
end

function _labels(atoms)
    return Set(a.label for a in atoms)
end

function _distance_stats(atoms)
    min_all = Inf
    min_heavy = Inf
    pair_all = (0, 0)
    pair_heavy = (0, 0)
    for i in 1:length(atoms)-1, j in i+1:length(atoms)
        d = norm(atoms[i].xyz .- atoms[j].xyz)
        if d < min_all
            min_all = d
            pair_all = (i, j)
        end
        if atoms[i].element != "H" && atoms[j].element != "H" && d < min_heavy
            min_heavy = d
            pair_heavy = (i, j)
        end
    end
    return min_all, pair_all, min_heavy, pair_heavy
end

function _center_of_labels(atoms, wanted)
    idx = Dict(a.label => i for (i, a) in enumerate(atoms))
    for label in wanted
        haskey(idx, label) || error("Missing label $label")
    end
    c = zeros(3)
    for label in wanted
        c .+= atoms[idx[label]].xyz
    end
    return c ./ length(wanted)
end

function _report_row(rows, structure, key, value)
    push!(rows, (structure=structure, key=key, value=string(value)))
end

function _validate_case!(rows, dir::String, name::String, expected_type::Int, expected_atoms::Int, expected_acetyl::Int)
    bare_path = joinpath(dir, "$(name)_central_trimer.xyz")
    slab_path = joinpath(dir, "$(name)_central_trimer_slab.xyz")
    index_path = joinpath(dir, "$(name)_central_trimer_indices.tsv")
    meta_path = joinpath(dir, "$(name)_central_trimer_slab_meta.tsv")
    for path in (bare_path, slab_path, index_path, meta_path)
        isfile(path) || error("Missing required file: $path")
    end

    bare = _read_xyz(bare_path)
    slab = _read_xyz(slab_path)
    idx = _read_key_tsv(index_path)
    meta = _read_key_tsv(meta_path)
    central_labels = split(idx["origin_labels"], ',')
    bare_indices = _parse_ints(idx["origin_indices_bare"])
    slab_offset = parse(Int, idx["default_slab_offset_8x8x4"])
    slab_indices = [i + slab_offset for i in bare_indices]

    labels = _labels(bare)
    acetyl_units = "C_AcC" in labels ? 1 : 0
    left_acetyl = "L_AcC" in labels ? 1 : 0
    right_acetyl = "R_AcC" in labels ? 1 : 0

    @assert length(bare) == expected_atoms
    @assert parse(Int, idx["central_type"]) == expected_type
    @assert acetyl_units == expected_acetyl
    @assert left_acetyl == 0 && right_acetyl == 0
    @assert length(slab) == parse(Int, meta["n_atoms"])
    @assert parse(Int, meta["n_atoms"]) == slab_offset + length(bare)

    for (bare_i, slab_i) in zip(bare_indices, slab_indices)
        @assert bare[bare_i].label == slab[slab_i].label
        @assert bare[bare_i].element == slab[slab_i].element
    end

    bare_min, bare_pair, bare_heavy, bare_heavy_pair = _distance_stats(bare)
    slab_min, slab_pair, slab_heavy, slab_heavy_pair = _distance_stats(slab)
    @assert bare_min > 0.75
    @assert bare_heavy > 1.10
    @assert slab_min > 0.75
    @assert slab_heavy > 1.10

    center = _center_of_labels(slab, central_labels)
    cell_a = [parse(Float64, x) for x in split(meta["cell_a"], ',')]
    cell_b = [parse(Float64, x) for x in split(meta["cell_b"], ',')]
    @assert abs(center[1] - cell_a[1] / 2) < 1e-6
    @assert abs(center[2] - cell_b[2] / 2) < 1e-6

    _report_row(rows, name, "bare_xyz", bare_path)
    _report_row(rows, name, "slab_xyz", slab_path)
    _report_row(rows, name, "bare_atoms", length(bare))
    _report_row(rows, name, "slab_atoms", length(slab))
    _report_row(rows, name, "bare_formula", _formula(bare))
    _report_row(rows, name, "central_type", expected_type)
    _report_row(rows, name, "central_acetyl_units", acetyl_units)
    _report_row(rows, name, "neighbor_acetyl_units", left_acetyl + right_acetyl)
    _report_row(rows, name, "bare_min_distance_ang", @sprintf("%.4f", bare_min))
    _report_row(rows, name, "bare_min_distance_pair", "$(bare[bare_pair[1]].label),$(bare[bare_pair[2]].label)")
    _report_row(rows, name, "bare_min_heavy_distance_ang", @sprintf("%.4f", bare_heavy))
    _report_row(rows, name, "bare_min_heavy_pair", "$(bare[bare_heavy_pair[1]].label),$(bare[bare_heavy_pair[2]].label)")
    _report_row(rows, name, "slab_min_distance_ang", @sprintf("%.4f", slab_min))
    _report_row(rows, name, "slab_min_distance_pair", "$(slab[slab_pair[1]].label),$(slab[slab_pair[2]].label)")
    _report_row(rows, name, "slab_min_heavy_distance_ang", @sprintf("%.4f", slab_heavy))
    _report_row(rows, name, "slab_min_heavy_pair", "$(slab[slab_heavy_pair[1]].label),$(slab[slab_heavy_pair[2]].label)")
    _report_row(rows, name, "central_ring_center_x_ang", @sprintf("%.8f", center[1]))
    _report_row(rows, name, "central_ring_center_y_ang", @sprintf("%.8f", center[2]))
    _report_row(rows, name, "frame_command_after_default_slab", idx["frame_command_after_default_slab"])
end

function _parse_cli(args)
    dir = DEFAULT_DIR
    out = DEFAULT_OUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--dir"
            dir = args[i+1]; i += 2
        elseif startswith(arg, "--dir=")
            dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            out = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/validate_chitosan_trimer_structures.jl [options]

            Options:
              --dir PATH   Directory containing generated trimer/slab files [$(DEFAULT_DIR)]
              --out PATH   Output validation TSV [$(DEFAULT_OUT)]

            Checks atom counts, central acetyl count, GlcN neighbors, minimum
            distances, slab label preservation, central-ring centering, and frame
            index consistency. Truth sequences are not used.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return dir, out
end

function main()
    dir, out = _parse_cli(ARGS)
    rows = NamedTuple[]
    _validate_case!(rows, dir, "glcn", 0, 69, 0)
    _validate_case!(rows, dir, "glcnac", 1, 74, 1)
    _ensure_parent(out)
    open(out, "w") do io
        println(io, "structure\tkey\tvalue")
        for row in rows
            println(io, join([row.structure, row.key, row.value], '\t'))
        end
    end
    println("Validated initial chitosan trimer structures")
    println("  dir:    ", dir)
    println("  report: ", out)
    println("  status: OK")
end

main()
