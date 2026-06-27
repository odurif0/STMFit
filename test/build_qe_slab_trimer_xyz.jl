#!/usr/bin/env julia

# Build a simple Cu(100) slab + pre-built chitosan trimer XYZ for QE mold runs.
#
# This helper does not create the chitosan chemistry. It expects an already
# vetted, oriented trimer XYZ and only adds a reproducible Cu(100) slab below it.

using Printf
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _parse_ints, _parse_vec3

const DEFAULT_OUT = "hpc/qe_molds/trimer_slab.xyz"

struct Atom
    element::String
    xyz::NTuple{3,Float64}
    label::String
end

struct Options
    molecule::String
    out_xyz::String
    metadata::String
    nx::Int
    ny::Int
    layers::Int
    a_cu::Float64
    height_above_top::Float64
    vacuum::Float64
    center_xy::Bool
    center_indices::Vector{Int}
    shift::Vector{Float64}
end

function _parse_cli(args)
    molecule = ""
    out_xyz = DEFAULT_OUT
    metadata = ""
    nx = 8
    ny = 8
    layers = 4
    a_cu = 3.615
    height_above_top = 2.6
    vacuum = 18.0
    center_xy = true
    center_indices = Int[]
    shift = [0.0, 0.0, 0.0]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--molecule"; molecule = args[i+1]; i += 2
        elseif startswith(arg, "--molecule="); molecule = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_xyz = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_xyz = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--metadata"; metadata = args[i+1]; i += 2
        elseif startswith(arg, "--metadata="); metadata = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--nx"; nx = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--nx="); nx = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--ny"; ny = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--ny="); ny = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--layers"; layers = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--layers="); layers = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--a-cu"; a_cu = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--a-cu="); a_cu = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--height-above-top"; height_above_top = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--height-above-top="); height_above_top = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--vacuum"; vacuum = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--vacuum="); vacuum = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--shift"; shift = _parse_vec3(args[i+1]); i += 2
        elseif startswith(arg, "--shift="); shift = _parse_vec3(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--center-indices"; center_indices = _parse_ints(args[i+1]); i += 2
        elseif startswith(arg, "--center-indices="); center_indices = _parse_ints(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--no-center-xy"; center_xy = false; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/build_qe_slab_trimer_xyz.jl [options]

            Required:
              --molecule PATH          Oriented chitosan trimer XYZ in angstrom

            Options:
              --out PATH               Combined slab+trimer XYZ [$(DEFAULT_OUT)]
              --metadata PATH          Optional TSV with cell vectors and slab info
              --nx INT                 Cu surface repeats along x [8]
              --ny INT                 Cu surface repeats along y [8]
              --layers INT             Cu(100) layers [4]
              --a-cu FLOAT             fcc Cu lattice constant, angstrom [3.615]
              --height-above-top FLOAT Place molecule min-z this far above top Cu [2.6]
              --vacuum FLOAT           Vacuum above highest atom for suggested cell-c [18]
              --shift X,Y,Z            Additional molecule shift after placement [0,0,0]
              --center-indices I,J,...  Center these molecule atom indices in xy
              --no-center-xy           Do not center molecule over slab xy box

            Output:
              XYZ with slab atoms first, then molecule atoms. The comment line
              includes suggested QE cell vectors. The script does not build or
              rotate the beta-(1->4) trimer; supply that geometry explicitly.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(molecule) && error("--molecule is required")
    isfile(molecule) || error("Molecule XYZ not found: $molecule")
    nx > 1 || error("--nx must be > 1")
    ny > 1 || error("--ny must be > 1")
    layers > 1 || error("--layers must be > 1")
    a_cu > 0 || error("--a-cu must be positive")
    height_above_top > 0 || error("--height-above-top must be positive")
    vacuum > 0 || error("--vacuum must be positive")
    return Options(molecule, out_xyz, metadata, nx, ny, layers, a_cu,
                   height_above_top, vacuum, center_xy, center_indices, shift)
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
        label = length(parts) >= 5 ? parts[5] : "$(parts[1])$(length(atoms) + 1)"
        push!(atoms, Atom(parts[1], (parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])), label))
    end
    return atoms
end

function _cu100_slab(nx::Int, ny::Int, layers::Int, a_cu::Float64)
    s = a_cu / sqrt(2.0)
    dz = a_cu / 2.0
    atoms = Atom[]
    for layer in 0:(layers-1)
        offset = isodd(layer) ? s / 2.0 : 0.0
        z = layer * dz
        for ix in 0:(nx-1), iy in 0:(ny-1)
            label = @sprintf("Cu_L%d_X%d_Y%d", layer + 1, ix + 1, iy + 1)
            push!(atoms, Atom("Cu", (ix * s + offset, iy * s + offset, z), label))
        end
    end
    cell_a = [nx * s, 0.0, 0.0]
    cell_b = [0.0, ny * s, 0.0]
    top_z = (layers - 1) * dz
    return atoms, cell_a, cell_b, top_z
end

function _bounds(atoms)
    xs = [a.xyz[1] for a in atoms]
    ys = [a.xyz[2] for a in atoms]
    zs = [a.xyz[3] for a in atoms]
    return (xmin=minimum(xs), xmax=maximum(xs), ymin=minimum(ys), ymax=maximum(ys), zmin=minimum(zs), zmax=maximum(zs))
end

function _place_molecule(mol, cell_a, cell_b, top_z, opt::Options)
    b = _bounds(mol)
    if opt.center_xy && !isempty(opt.center_indices)
        for i in opt.center_indices
            1 <= i <= length(mol) || error("--center-indices atom $i outside 1:$(length(mol))")
        end
        cx = mean([mol[i].xyz[1] for i in opt.center_indices])
        cy = mean([mol[i].xyz[2] for i in opt.center_indices])
        dx = cell_a[1] / 2.0 - cx
        dy = cell_b[2] / 2.0 - cy
    else
        dx = opt.center_xy ? cell_a[1] / 2.0 - 0.5 * (b.xmin + b.xmax) : 0.0
        dy = opt.center_xy ? cell_b[2] / 2.0 - 0.5 * (b.ymin + b.ymax) : 0.0
    end
    dz = top_z + opt.height_above_top - b.zmin
    dx += opt.shift[1]
    dy += opt.shift[2]
    dz += opt.shift[3]
    return [Atom(a.element, (a.xyz[1] + dx, a.xyz[2] + dy, a.xyz[3] + dz), a.label) for a in mol]
end

function _write_xyz(path::String, atoms, cell_a, cell_b, cell_c)
    _ensure_parent(path)
    open(path, "w") do io
        println(io, length(atoms))
        println(io, @sprintf("cell_a=%.8f,%.8f,%.8f cell_b=%.8f,%.8f,%.8f cell_c=%.8f,%.8f,%.8f",
                            cell_a[1], cell_a[2], cell_a[3], cell_b[1], cell_b[2], cell_b[3], cell_c[1], cell_c[2], cell_c[3]))
        for a in atoms
            x, y, z = a.xyz
            println(io, @sprintf("%-2s  %.10f  %.10f  %.10f  %s", a.element, x, y, z, a.label))
        end
    end
end

function _write_metadata(path::String, opt::Options, cell_a, cell_b, cell_c, top_z, atoms)
    isempty(path) && return
    _ensure_parent(path)
    open(path, "w") do io
        println(io, "key\tvalue")
        println(io, "cell_a\t", join(cell_a, ','))
        println(io, "cell_b\t", join(cell_b, ','))
        println(io, "cell_c\t", join(cell_c, ','))
        println(io, "cell_a_arg\t--cell-a ", join(cell_a, ','))
        println(io, "cell_b_arg\t--cell-b ", join(cell_b, ','))
        println(io, "cell_c_arg\t--cell-c ", join(cell_c, ','))
        println(io, "fix_below_z_suggestion\t", top_z / max(opt.layers - 1, 1) + 1e-6)
        println(io, "top_cu_z\t", top_z)
        println(io, "n_atoms\t", length(atoms))
        println(io, "center_indices\t", isempty(opt.center_indices) ? "bbox" : join(opt.center_indices, ','))
    end
end

function main()
    opt = _parse_cli(ARGS)
    molecule = _read_xyz(opt.molecule)
    slab, cell_a, cell_b, top_z = _cu100_slab(opt.nx, opt.ny, opt.layers, opt.a_cu)
    placed = _place_molecule(molecule, cell_a, cell_b, top_z, opt)
    atoms = vcat(slab, placed)
    b = _bounds(atoms)
    cell_c = [0.0, 0.0, b.zmax + opt.vacuum]
    _write_xyz(opt.out_xyz, atoms, cell_a, cell_b, cell_c)
    _write_metadata(opt.metadata, opt, cell_a, cell_b, cell_c, top_z, atoms)
    println("Built Cu(100)+trimer XYZ")
    println("  molecule: ", opt.molecule)
    println("  output:   ", opt.out_xyz)
    isempty(opt.metadata) || println("  metadata: ", opt.metadata)
    println("  atoms:    slab=", length(slab), " molecule=", length(placed), " total=", length(atoms))
    println("  cell-a:   ", join(cell_a, ','))
    println("  cell-b:   ", join(cell_b, ','))
    println("  cell-c:   ", join(cell_c, ','))
    @printf("  top Cu z: %.6f Å\n", top_z)
end

main()
