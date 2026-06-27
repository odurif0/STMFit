#!/usr/bin/env julia

# Replace the ATOMIC_POSITIONS block in a Quantum ESPRESSO input with coordinates
# from an XYZ file. Intended for relax -> SCF handoff.

using Printf

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent

struct Atom
    element::String
    xyz::NTuple{3,Float64}
end

struct Options
    input_qe::String
    xyz::String
    out_qe::String
end

function _parse_cli(args)
    input_qe = ""
    xyz = ""
    out_qe = ""
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input"; input_qe = args[i+1]; i += 2
        elseif startswith(arg, "--input="); input_qe = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--xyz"; xyz = args[i+1]; i += 2
        elseif startswith(arg, "--xyz="); xyz = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_qe = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_qe = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/update_qe_positions_from_xyz.jl [options]

            Required:
              --input PATH        Existing QE input containing ATOMIC_POSITIONS
              --xyz PATH          XYZ coordinates in angstrom
              --out PATH          Updated QE input path

            The output uses ATOMIC_POSITIONS angstrom and drops relaxation flags,
            which is appropriate for SCF/PP handoff after a relax calculation.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(input_qe) && error("--input is required")
    isempty(xyz) && error("--xyz is required")
    isempty(out_qe) && error("--out is required")
    isfile(input_qe) || error("QE input not found: $input_qe")
    isfile(xyz) || error("XYZ not found: $xyz")
    return Options(input_qe, xyz, out_qe)
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
        push!(atoms, Atom(parts[1], (parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4]))))
    end
    return atoms
end

function _section_starts(line::String)
    s = strip(line)
    return startswith(s, "K_POINTS") || startswith(s, "CELL_PARAMETERS") ||
           startswith(s, "ATOMIC_SPECIES") || startswith(s, "ATOMIC_FORCES") ||
           startswith(s, "CONSTRAINTS") || startswith(s, "OCCUPATIONS") ||
           startswith(s, "ATOMIC_VELOCITIES") || startswith(s, "&")
end

function _replace_positions(lines, atoms)
    idx = findfirst(line -> startswith(strip(line), "ATOMIC_POSITIONS"), lines)
    idx === nothing && error("No ATOMIC_POSITIONS block found in QE input")
    stop = idx + 1
    while stop <= length(lines)
        s = strip(lines[stop])
        if isempty(s) || _section_starts(lines[stop])
            break
        end
        stop += 1
    end
    out = String[]
    append!(out, lines[1:(idx-1)])
    push!(out, "ATOMIC_POSITIONS angstrom")
    for a in atoms
        x, y, z = a.xyz
        push!(out, @sprintf("  %-2s  %.10f  %.10f  %.10f", a.element, x, y, z))
    end
    append!(out, lines[stop:end])
    return out
end

function main()
    opt = _parse_cli(ARGS)
    atoms = _read_xyz(opt.xyz)
    lines = readlines(opt.input_qe)
    out_lines = _replace_positions(lines, atoms)
    _ensure_parent(opt.out_qe)
    open(opt.out_qe, "w") do io
        for line in out_lines
            println(io, line)
        end
    end
    println("Updated QE positions")
    println("  input: ", opt.input_qe)
    println("  xyz:   ", opt.xyz)
    println("  out:   ", opt.out_qe)
    println("  atoms: ", length(atoms))
end

main()
