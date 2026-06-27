#!/usr/bin/env julia

# Prepare Quantum ESPRESSO inputs for DFT-STM mold calculations from an XYZ file.
#
# The XYZ is expected to contain the complete model: Cu(100) slab + adsorbed
# beta-(1->4)-linked chitosan trimer. This script does not build or alter the
# chemistry; it only turns a vetted structure into QE input files.

using Printf

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_vec3, _read_key_tsv

const DEFAULT_OUTDIR = "qe"

const MASSES = Dict(
    "H" => 1.008,
    "C" => 12.011,
    "N" => 14.007,
    "O" => 15.999,
    "Cu" => 63.546,
)

const DEFAULT_PSEUDOS = Dict(
    "Cu" => "Cu.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "C" => "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "N" => "N.pbe-n-kjpaw_psl.1.0.0.UPF",
    "O" => "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "H" => "H.pbe-kjpaw_psl.1.0.0.UPF",
)

struct Atom
    element::String
    xyz::NTuple{3,Float64}
end

struct Options
    xyz::String
    out_dir::String
    prefix::String
    cell_a::Vector{Float64}
    cell_b::Vector{Float64}
    cell_c::Vector{Float64}
    pseudo_dir::String
    pseudo_map::Dict{String,String}
    ecutwfc::Float64
    ecutrho::Float64
    kpoints::NTuple{3,Int}
    emin_ev::Float64
    emax_ev::Float64
    fix_below_z::Union{Nothing,Float64}
    ntasks::Int
    mem_per_task_mb::Int
    walltime::String
end

function _parse_kpoints(s::AbstractString)
    vals = [parse(Int, strip(x)) for x in split(s, ',') if !isempty(strip(x))]
    length(vals) == 3 || error("expected kx,ky,kz, got: $s")
    return (vals[1], vals[2], vals[3])
end

function _parse_pseudo!(pseudo_map, spec::AbstractString)
    parts = split(spec, '=', limit=2)
    length(parts) == 2 || error("--pseudo expects ELEMENT=FILE.UPF")
    element = strip(parts[1])
    file = strip(parts[2])
    isempty(element) && error("empty pseudo element")
    isempty(file) && error("empty pseudo file")
    pseudo_map[element] = file
end

function _read_cell_metadata(path::String)
    isfile(path) || error("Cell metadata not found: $path")
    d = _read_key_tsv(path)
    for key in ("cell_a", "cell_b", "cell_c")
        haskey(d, key) || error("Cell metadata missing $key: $path")
    end
    return _parse_vec3(d["cell_a"]), _parse_vec3(d["cell_b"]), _parse_vec3(d["cell_c"])
end

function _parse_cli(args)
    xyz = ""
    out_dir = DEFAULT_OUTDIR
    prefix = "chitosan_mold"
    cell_a::Union{Nothing,Vector{Float64}} = nothing
    cell_b::Union{Nothing,Vector{Float64}} = nothing
    cell_c::Union{Nothing,Vector{Float64}} = nothing
    cell_metadata = ""
    pseudo_dir = "./pseudo"
    pseudo_map = copy(DEFAULT_PSEUDOS)
    ecutwfc = 50.0
    ecutrho = 360.0
    kpoints = (1, 1, 1)
    emin_ev = -1.0
    emax_ev = 0.0
    fix_below_z::Union{Nothing,Float64} = nothing
    ntasks = 8
    mem_per_task_mb = 12000
    walltime = "24:00:00"

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--xyz"; xyz = args[i+1]; i += 2
        elseif startswith(arg, "--xyz="); xyz = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out-dir"; out_dir = args[i+1]; i += 2
        elseif startswith(arg, "--out-dir="); out_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--prefix"; prefix = args[i+1]; i += 2
        elseif startswith(arg, "--prefix="); prefix = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--cell-a"; cell_a = _parse_vec3(args[i+1]); i += 2
        elseif startswith(arg, "--cell-a="); cell_a = _parse_vec3(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cell-b"; cell_b = _parse_vec3(args[i+1]); i += 2
        elseif startswith(arg, "--cell-b="); cell_b = _parse_vec3(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cell-c"; cell_c = _parse_vec3(args[i+1]); i += 2
        elseif startswith(arg, "--cell-c="); cell_c = _parse_vec3(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cell-metadata"; cell_metadata = args[i+1]; i += 2
        elseif startswith(arg, "--cell-metadata="); cell_metadata = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--pseudo-dir"; pseudo_dir = args[i+1]; i += 2
        elseif startswith(arg, "--pseudo-dir="); pseudo_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--pseudo"; _parse_pseudo!(pseudo_map, args[i+1]); i += 2
        elseif startswith(arg, "--pseudo="); _parse_pseudo!(pseudo_map, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--ecutwfc"; ecutwfc = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--ecutwfc="); ecutwfc = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--ecutrho"; ecutrho = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--ecutrho="); ecutrho = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--kpoints"; kpoints = _parse_kpoints(args[i+1]); i += 2
        elseif startswith(arg, "--kpoints="); kpoints = _parse_kpoints(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--emin-ev"; emin_ev = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--emin-ev="); emin_ev = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--emax-ev"; emax_ev = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--emax-ev="); emax_ev = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--fix-below-z"; fix_below_z = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--fix-below-z="); fix_below_z = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--ntasks"; ntasks = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--ntasks="); ntasks = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--mem-per-task-mb"; mem_per_task_mb = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--mem-per-task-mb="); mem_per_task_mb = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--walltime"; walltime = args[i+1]; i += 2
        elseif startswith(arg, "--walltime="); walltime = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/prepare_qe_mold_inputs.jl [options]

            Required:
              --xyz PATH             Complete slab+trimer XYZ in angstrom
              --cell-metadata PATH    Metadata from build_qe_slab_trimer_xyz.jl
                  or explicitly pass all of --cell-a/--cell-b/--cell-c

            Options:
              --out-dir PATH         Output directory [$(DEFAULT_OUTDIR)]
              --prefix STR           QE prefix [chitosan_mold]
              --cell-metadata PATH    Read cell_a/b/c from build_qe_slab_trimer_xyz metadata
              --cell-a X,Y,Z         Override QE cell vector A in angstrom
              --cell-b X,Y,Z         Override QE cell vector B in angstrom
              --cell-c X,Y,Z         Override QE cell vector C in angstrom
              --pseudo-dir PATH      QE pseudo_dir [./pseudo]
              --pseudo ELEM=FILE     Override pseudo file; repeatable
              --ecutwfc FLOAT        Wavefunction cutoff Ry [50]
              --ecutrho FLOAT        Charge-density cutoff Ry [360]
              --kpoints KX,KY,KZ     K-point grid; 1,1,1 writes K_POINTS gamma [1,1,1]
              --emin-ev FLOAT        LDOS lower bound for pp.x [-1]
              --emax-ev FLOAT        LDOS upper bound for pp.x [0]
              --fix-below-z FLOAT    Freeze atoms with z <= cutoff angstrom
              --ntasks INT           Slurm MPI tasks [8]
              --mem-per-task-mb INT  Slurm memory per MPI task MB [12000]
              --walltime HH:MM:SS    Slurm walltime [24:00:00]

            Output files:
              pw_relax.in, pw_scf.in, pp_ldos.in, run_qe_mold.sbatch

            The script does not generate molecular geometry. It only serializes
            a user-vetted Cu(100)+trimer structure into QE inputs. Truth labels
            and benchmark composition are not used.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isempty(xyz) && error("--xyz is required")
    isfile(xyz) || error("XYZ not found: $xyz")
    if !isempty(cell_metadata)
        md_a, md_b, md_c = _read_cell_metadata(cell_metadata)
        cell_a === nothing && (cell_a = md_a)
        cell_b === nothing && (cell_b = md_b)
        cell_c === nothing && (cell_c = md_c)
    end
    cell_a === nothing && error("--cell-a is required unless --cell-metadata supplies cell_a")
    cell_b === nothing && error("--cell-b is required unless --cell-metadata supplies cell_b")
    cell_c === nothing && error("--cell-c is required unless --cell-metadata supplies cell_c")
    ecutwfc > 0 || error("--ecutwfc must be positive")
    ecutrho > 0 || error("--ecutrho must be positive")
    ntasks > 0 || error("--ntasks must be positive")
    mem_per_task_mb > 0 || error("--mem-per-task-mb must be positive")
    return Options(xyz, out_dir, prefix, cell_a, cell_b, cell_c, pseudo_dir,
                   pseudo_map, ecutwfc, ecutrho, kpoints, emin_ev, emax_ev,
                   fix_below_z, ntasks, mem_per_task_mb, walltime)
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
        elem = parts[1]
        xyz = (parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4]))
        push!(atoms, Atom(elem, xyz))
    end
    return atoms
end

function _species(atoms)
    preferred = ["Cu", "C", "O", "H", "N"]
    present = Set(a.element for a in atoms)
    ordered = [e for e in preferred if e in present]
    append!(ordered, sort([e for e in present if !(e in preferred)]))
    return ordered
end

function _write_common_system(io, opt::Options, atoms, species)
    println(io, "&SYSTEM")
    println(io, "  ibrav = 0")
    println(io, "  nat = ", length(atoms))
    println(io, "  ntyp = ", length(species))
    println(io, @sprintf("  ecutwfc = %.8g", opt.ecutwfc))
    println(io, @sprintf("  ecutrho = %.8g", opt.ecutrho))
    println(io, "  occupations = 'smearing'")
    println(io, "  smearing = 'mv'")
    println(io, "  degauss = 0.02")
    println(io, "  input_dft = 'PBE'")
    println(io, "  vdw_corr = 'grimme-d3'")
    println(io, "/")
end

function _write_species(io, species, pseudo_map)
    println(io, "ATOMIC_SPECIES")
    for e in species
        mass = get(MASSES, e, 1.0)
        pseudo = get(pseudo_map, e, "$(e).UPF")
        println(io, @sprintf("  %-2s  %.6f  %s", e, mass, pseudo))
    end
end

function _write_cell(io, opt::Options)
    println(io, "CELL_PARAMETERS angstrom")
    for v in (opt.cell_a, opt.cell_b, opt.cell_c)
        println(io, @sprintf("  %.10f  %.10f  %.10f", v[1], v[2], v[3]))
    end
end

function _write_positions(io, atoms, fix_below_z)
    println(io, "ATOMIC_POSITIONS angstrom")
    for a in atoms
        x, y, z = a.xyz
        if fix_below_z !== nothing && z <= fix_below_z
            println(io, @sprintf("  %-2s  %.10f  %.10f  %.10f  0 0 0", a.element, x, y, z))
        else
            println(io, @sprintf("  %-2s  %.10f  %.10f  %.10f", a.element, x, y, z))
        end
    end
end

function _write_kpoints(io, opt::Options)
    kx, ky, kz = opt.kpoints
    if (kx, ky, kz) == (1, 1, 1)
        println(io, "K_POINTS gamma")
    else
        println(io, "K_POINTS automatic")
        println(io, "  $kx $ky $kz  0 0 0")
    end
end

function _write_pw(path, opt::Options, atoms, species; calculation::String, conv_thr::String)
    open(path, "w") do io
        println(io, "&CONTROL")
        println(io, "  calculation = '$calculation'")
        println(io, "  prefix = '$(opt.prefix)'")
        println(io, "  outdir = './qe_tmp'")
        println(io, "  pseudo_dir = '$(opt.pseudo_dir)'")
        println(io, "  tprnfor = .true.")
        println(io, "  tstress = .true.")
        println(io, "/")
        _write_common_system(io, opt, atoms, species)
        println(io, "&ELECTRONS")
        println(io, "  conv_thr = $conv_thr")
        println(io, "  mixing_beta = 0.3")
        println(io, "/")
        if calculation == "relax"
            println(io, "&IONS")
            println(io, "  ion_dynamics = 'bfgs'")
            println(io, "/")
        end
        _write_species(io, species, opt.pseudo_map)
        _write_cell(io, opt)
        _write_positions(io, atoms, calculation == "relax" ? opt.fix_below_z : nothing)
        _write_kpoints(io, opt)
    end
end

function _write_pp(path, opt::Options)
    open(path, "w") do io
        println(io, "&INPUTPP")
        println(io, "  prefix = '$(opt.prefix)'")
        println(io, "  outdir = './qe_tmp'")
        println(io, "  plot_num = 5")
        println(io, @sprintf("  emin = %.8g", opt.emin_ev))
        println(io, @sprintf("  emax = %.8g", opt.emax_ev))
        println(io, "  degauss_ldos = 0.02")
        println(io, "/")
        println(io, "&PLOT")
        println(io, "  iflag = 3")
        println(io, "  output_format = 6")
        println(io, "  fileout = '$(opt.prefix)_ldos.cube'")
        println(io, "/")
    end
end

function _write_sbatch(path, opt::Options)
    open(path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "#SBATCH --job-name=qe-$(opt.prefix)")
        println(io, "#SBATCH --nodes=1")
        println(io, "#SBATCH --ntasks-per-node=$(opt.ntasks)")
        println(io, "#SBATCH --cpus-per-task=1")
        println(io, "#SBATCH --time=$(opt.walltime)")
        println(io, "#SBATCH --mem=$(opt.ntasks * opt.mem_per_task_mb)MB")
        println(io, "#SBATCH --output=qe_$(opt.prefix)_%j.out")
        println(io, "#SBATCH --error=qe_$(opt.prefix)_%j.err")
        println(io)
        println(io, "set -euo pipefail")
        println(io, "module purge")
        println(io, "QE_COMPILER_MODULE=\${QE_COMPILER_MODULE:-intel/2024.0}")
        println(io, "QE_MPI_MODULE=\${QE_MPI_MODULE:-impi/2021.11}")
        println(io, "QE_MODULE=\${QE_MODULE:-qe/7.4.1}")
        println(io, "[[ -n \"\$QE_COMPILER_MODULE\" ]] && module load \"\$QE_COMPILER_MODULE\"")
        println(io, "[[ -n \"\$QE_MPI_MODULE\" ]] && module load \"\$QE_MPI_MODULE\"")
        println(io, "module load \"\$QE_MODULE\"")
        println(io, "JULIA_MODULE_VERSION=\${JULIA_MODULE_VERSION:-}")
        println(io, "module load \"julia\${JULIA_MODULE_VERSION:+/\$JULIA_MODULE_VERSION}\"")
        println(io, "command -v pw.x >/dev/null 2>&1 || { echo 'ERROR: pw.x not found after loading QE module stack' >&2; exit 1; }")
        println(io, "command -v julia >/dev/null 2>&1 || { echo 'ERROR: julia not found after loading Julia module' >&2; exit 1; }")
        println(io, "export OMP_NUM_THREADS=1")
        println(io, "QE_NTASKS=\${SLURM_NTASKS:-$(opt.ntasks)}")
        println(io, "echo \"QE_NTASKS=\$QE_NTASKS OMP_NUM_THREADS=\$OMP_NUM_THREADS SLURM_CPUS_ON_NODE=\${SLURM_CPUS_ON_NODE:-NA}\"")
        println(io, "STMFIT_ROOT=\${STMFIT_ROOT:-../..}")
        println(io, "rm -rf qe_tmp")
        println(io, "mkdir -p qe_tmp")
        println(io, "srun -n \"\$QE_NTASKS\" --cpu-bind=cores pw.x -in pw_relax.in > $(opt.prefix)_relax.out")
        println(io, "julia --project=\"\${STMFIT_ROOT}\" \"\${STMFIT_ROOT}/test/extract_qe_relaxed_xyz.jl\" --qe-out $(opt.prefix)_relax.out --out $(opt.prefix)_relaxed.xyz --metadata $(opt.prefix)_relaxed_meta.tsv")
        println(io, "julia --project=\"\${STMFIT_ROOT}\" \"\${STMFIT_ROOT}/test/update_qe_positions_from_xyz.jl\" --input pw_scf.in --xyz $(opt.prefix)_relaxed.xyz --out pw_scf_relaxed.in")
        println(io, "srun -n \"\$QE_NTASKS\" --cpu-bind=cores pw.x -in pw_scf_relaxed.in > $(opt.prefix)_scf.out")
        println(io, "srun -n \"\$QE_NTASKS\" --cpu-bind=cores pp.x -in pp_ldos.in > $(opt.prefix)_pp.out")
    end
end

function main()
    opt = _parse_cli(ARGS)
    atoms = _read_xyz(opt.xyz)
    species = _species(atoms)
    mkpath(opt.out_dir)
    _write_pw(joinpath(opt.out_dir, "pw_relax.in"), opt, atoms, species; calculation="relax", conv_thr="1.0d-6")
    _write_pw(joinpath(opt.out_dir, "pw_scf.in"), opt, atoms, species; calculation="scf", conv_thr="1.0d-7")
    _write_pp(joinpath(opt.out_dir, "pp_ldos.in"), opt)
    _write_sbatch(joinpath(opt.out_dir, "run_qe_mold.sbatch"), opt)
    println("Prepared QE mold inputs")
    println("  xyz:      ", opt.xyz)
    println("  out dir:  ", opt.out_dir)
    println("  prefix:   ", opt.prefix)
    println("  atoms:    ", length(atoms))
    println("  species:  ", join(species, ", "))
    println("  files:    pw_relax.in, pw_scf.in, pp_ldos.in, run_qe_mold.sbatch")
end

main()
