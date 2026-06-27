#!/usr/bin/env julia

# Extract the last ATOMIC_POSITIONS block from a Quantum ESPRESSO pw.x output and
# write it as XYZ. Intended for the relax -> SCF/frame handoff in DFT-STM mold
# calculations.

using Printf

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent

const BOHR_TO_ANG = 0.529177210903

struct Atom
    element::String
    xyz::NTuple{3,Float64}
end

struct Options
    qe_out::String
    xyz_out::String
    metadata::String
    alat_angstrom::Union{Nothing,Float64}
end

function _parse_cli(args)
    qe_out = ""
    xyz_out = ""
    metadata = ""
    alat_angstrom::Union{Nothing,Float64} = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--qe-out"; qe_out = args[i+1]; i += 2
        elseif startswith(arg, "--qe-out="); qe_out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; xyz_out = args[i+1]; i += 2
        elseif startswith(arg, "--out="); xyz_out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--metadata"; metadata = args[i+1]; i += 2
        elseif startswith(arg, "--metadata="); metadata = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--alat-angstrom"; alat_angstrom = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--alat-angstrom="); alat_angstrom = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/extract_qe_relaxed_xyz.jl [options]

            Required:
              --qe-out PATH          Quantum ESPRESSO pw.x output file
              --out PATH             Output XYZ path

            Options:
              --metadata PATH        Optional TSV with final cell if found
              --alat-angstrom FLOAT  Required only for ATOMIC_POSITIONS (alat)

            Notes:
              The script extracts the last ATOMIC_POSITIONS block and the last
              CELL_PARAMETERS block when present. Units angstrom and bohr are
              converted automatically; alat requires --alat-angstrom. It does
              not use truth labels or benchmark data.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(qe_out) && error("--qe-out is required")
    isempty(xyz_out) && error("--out is required")
    isfile(qe_out) || error("QE output not found: $qe_out")
    return Options(qe_out, xyz_out, metadata, alat_angstrom)
end

function _units_from_header(line::String)
    m = match(r"\(([^)]*)\)", lowercase(line))
    m === nothing && return "alat"
    u = strip(m.captures[1])
    occursin("angstrom", u) && return "angstrom"
    occursin("bohr", u) && return "bohr"
    occursin("crystal", u) && return "crystal"
    occursin("alat", u) && return "alat"
    return u
end

function _scale_for_units(units::String, opt::Options)
    units == "angstrom" && return 1.0
    units == "bohr" && return BOHR_TO_ANG
    if units == "alat"
        opt.alat_angstrom === nothing && error("ATOMIC_POSITIONS are in alat; pass --alat-angstrom")
        return opt.alat_angstrom
    end
    units == "crystal" && error("ATOMIC_POSITIONS (crystal) not supported by this extractor; rerun QE with angstrom positions")
    error("Unsupported ATOMIC_POSITIONS units: $units")
end

function _looks_like_atom_line(line::AbstractString)
    parts = split(strip(line))
    length(parts) >= 4 || return false
    occursin(r"^[A-Za-z][A-Za-z0-9_+-]*$", parts[1]) || return false
    return all(tryparse(Float64, parts[i]) !== nothing for i in 2:4)
end

function _parse_positions(lines, start_idx::Int, opt::Options)
    units = _units_from_header(lines[start_idx])
    scale = _scale_for_units(units, opt)
    atoms = Atom[]
    i = start_idx + 1
    while i <= length(lines)
        line = strip(lines[i])
        isempty(line) && break
        _looks_like_atom_line(line) || break
        parts = split(line)
        elem = replace(parts[1], r"[0-9_+-]+$" => "")
        push!(atoms, Atom(elem, (parse(Float64, parts[2]) * scale,
                                parse(Float64, parts[3]) * scale,
                                parse(Float64, parts[4]) * scale)))
        i += 1
    end
    isempty(atoms) && error("ATOMIC_POSITIONS block at line $start_idx has no atom rows")
    return atoms, units
end

function _parse_cell(lines, start_idx::Int, opt::Options)
    units = _units_from_header(lines[start_idx])
    scale = units == "angstrom" ? 1.0 : units == "bohr" ? BOHR_TO_ANG :
            units == "alat" ? _scale_for_units(units, opt) :
            error("CELL_PARAMETERS units $units unsupported")
    start_idx + 3 <= length(lines) || error("CELL_PARAMETERS block too short")
    cell = Vector{Float64}[]
    for line in lines[(start_idx+1):(start_idx+3)]
        parts = split(strip(line))
        length(parts) >= 3 || error("Bad CELL_PARAMETERS row: $line")
        push!(cell, [parse(Float64, parts[i]) * scale for i in 1:3])
    end
    return cell, units
end

function _extract(lines, opt::Options)
    last_atoms::Union{Nothing,Vector{Atom}} = nothing
    last_units = ""
    last_cell::Union{Nothing,Vector{Vector{Float64}}} = nothing
    last_cell_units = ""
    for (i, line) in enumerate(lines)
        s = strip(line)
        if startswith(s, "CELL_PARAMETERS")
            cell, units = _parse_cell(lines, i, opt)
            last_cell = cell
            last_cell_units = units
        elseif startswith(s, "ATOMIC_POSITIONS")
            atoms, units = _parse_positions(lines, i, opt)
            last_atoms = atoms
            last_units = units
        end
    end
    last_atoms === nothing && error("No ATOMIC_POSITIONS block found")
    return last_atoms, last_units, last_cell, last_cell_units
end

function _write_xyz(path::String, atoms, source::String, units::String)
    _ensure_parent(path)
    open(path, "w") do io
        println(io, length(atoms))
        println(io, "extracted_from=$(source) source_units=$(units)")
        for a in atoms
            x, y, z = a.xyz
            println(io, @sprintf("%-2s  %.10f  %.10f  %.10f", a.element, x, y, z))
        end
    end
end

function _write_metadata(path::String, cell, cell_units::String)
    isempty(path) && return
    _ensure_parent(path)
    open(path, "w") do io
        println(io, "key\tvalue")
        if cell === nothing
            println(io, "cell_found\tfalse")
        else
            println(io, "cell_found\ttrue")
            println(io, "source_cell_units\t", cell_units)
            println(io, "cell_a\t", join(cell[1], ','))
            println(io, "cell_b\t", join(cell[2], ','))
            println(io, "cell_c\t", join(cell[3], ','))
            println(io, "cell_a_arg\t--cell-a ", join(cell[1], ','))
            println(io, "cell_b_arg\t--cell-b ", join(cell[2], ','))
            println(io, "cell_c_arg\t--cell-c ", join(cell[3], ','))
        end
    end
end

function main()
    opt = _parse_cli(ARGS)
    lines = readlines(opt.qe_out)
    atoms, units, cell, cell_units = _extract(lines, opt)
    _write_xyz(opt.xyz_out, atoms, opt.qe_out, units)
    _write_metadata(opt.metadata, cell, cell_units)
    println("Extracted QE relaxed coordinates")
    println("  qe out:   ", opt.qe_out)
    println("  xyz:      ", opt.xyz_out)
    isempty(opt.metadata) || println("  metadata: ", opt.metadata)
    println("  atoms:    ", length(atoms))
    println("  units:    ", units, " -> angstrom")
    println("  cell:     ", cell === nothing ? "not found" : "found")
end

main()
