#!/usr/bin/env julia

# Post-process completed QE mold runs into STMFit connected mold templates.
# This assumes run_qe_mold.sbatch has produced relaxed XYZ files and LDOS cubes.

using Printf

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_ints, _read_key_tsv

const ROOT = normpath(joinpath(@__DIR__, ".."))

struct Options
    glcn_dir::String
    glcnac_dir::String
    height_nm::Union{Nothing,Float64}
    origin_indices::String
    axis_from::Int
    axis_to::Int
    plane_index::Int
    index_dir::String
    glcn_index_tsv::String
    glcnac_index_tsv::String
    cube_units::String
    half_nm::Float64
    step_nm::Float64
    maps_out::String
    templates_out::String
    bond_out::String
    dry_run::Bool
    check_only::Bool
end

function _parse_cli(args)
    glcn_dir = "qe/glcn"
    glcnac_dir = "qe/glcnac"
    height_nm::Union{Nothing,Float64} = nothing
    origin_indices = ""
    axis_from = 0
    axis_to = 0
    plane_index = 0
    index_dir = "hpc/qe_molds"
    glcn_index_tsv = ""
    glcnac_index_tsv = ""
    cube_units = "bohr"
    half_nm = 0.48
    step_nm = 0.08
    maps_out = "templates/chitosan_stm_maps.tsv"
    templates_out = "templates/chitosan_connected_molds_stm.tsv"
    bond_out = "templates/chitosan_connected_bond_molds_stm.tsv"
    dry_run = false
    check_only = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--glcn-dir"; glcn_dir = args[i+1]; i += 2
        elseif startswith(arg, "--glcn-dir="); glcn_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--glcnac-dir"; glcnac_dir = args[i+1]; i += 2
        elseif startswith(arg, "--glcnac-dir="); glcnac_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--height-nm"; height_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--height-nm="); height_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--origin-indices"; origin_indices = args[i+1]; i += 2
        elseif startswith(arg, "--origin-indices="); origin_indices = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--axis-from"; axis_from = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--axis-from="); axis_from = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--axis-to"; axis_to = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--axis-to="); axis_to = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--plane-index"; plane_index = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--plane-index="); plane_index = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--index-dir"; index_dir = args[i+1]; i += 2
        elseif startswith(arg, "--index-dir="); index_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--glcn-index-tsv"; glcn_index_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--glcn-index-tsv="); glcn_index_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--glcnac-index-tsv"; glcnac_index_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--glcnac-index-tsv="); glcnac_index_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--cube-units"; cube_units = lowercase(strip(args[i+1])); i += 2
        elseif startswith(arg, "--cube-units="); cube_units = lowercase(strip(split(arg, "=", limit=2)[2])); i += 1
        elseif arg == "--half-nm"; half_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--half-nm="); half_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--step-nm"; step_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--step-nm="); step_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--maps-out"; maps_out = args[i+1]; i += 2
        elseif startswith(arg, "--maps-out="); maps_out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--templates-out"; templates_out = args[i+1]; i += 2
        elseif startswith(arg, "--templates-out="); templates_out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--bond-out"; bond_out = args[i+1]; i += 2
        elseif startswith(arg, "--bond-out="); bond_out = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--dry-run"; dry_run = true; i += 1
        elseif arg == "--check-only"; check_only = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/finalize_qe_mold_workflow.jl [options]

            Required after QE has completed:
              --height-nm FLOAT       STM/LDOS sampling height in nm

            Options:
              --glcn-dir PATH         GlcN QE run directory [qe/glcn]
              --glcnac-dir PATH       GlcNAc QE run directory [qe/glcnac]
              --origin-indices LIST   Override relaxed slab central-ring indices [auto]
              --axis-from INT         Override frame axis start index [auto]
              --axis-to INT           Override frame axis end index [auto]
              --plane-index INT       Override positive-u side atom index [auto]
              --index-dir PATH        Directory containing *_trimer_indices.tsv [hpc/qe_molds]
              --glcn-index-tsv PATH   Override GlcN bare-index TSV [auto]
              --glcnac-index-tsv PATH Override GlcNAc bare-index TSV [auto]
              --cube-units STR        bohr | angstrom | nm [bohr]
              --half-nm FLOAT         Mold map half-size [0.48]
              --step-nm FLOAT         Mold map spacing [0.08]
              --maps-out PATH         Output STM map TSV [templates/chitosan_stm_maps.tsv]
              --templates-out PATH    Output connected unary mold TSV
              --bond-out PATH         Output connected bond mold TSV
              --check-only            Only check required QE outputs exist
              --dry-run               Print commands without running them

            This script uses no unit truth labels and imposes no composition prior.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    height_nm === nothing && !check_only && error("--height-nm is required unless --check-only")
    cube_units in ("bohr", "angstrom", "a", "nm") || error("invalid --cube-units")
    explicit = !isempty(origin_indices) || axis_from > 0 || axis_to > 0 || plane_index > 0
    if explicit
        !isempty(origin_indices) || error("--origin-indices is required when overriding frame indices")
        axis_from > 0 || error("--axis-from is required when overriding frame indices")
        axis_to > 0 || error("--axis-to is required when overriding frame indices")
        plane_index > 0 || error("--plane-index is required when overriding frame indices")
    end
    return Options(glcn_dir, glcnac_dir, height_nm, origin_indices, axis_from,
                   axis_to, plane_index, index_dir, glcn_index_tsv, glcnac_index_tsv,
                   cube_units, half_nm, step_nm, maps_out, templates_out, bond_out,
                   dry_run, check_only)
end

function _prefix(dir::String)
    pp = joinpath(dir, "pp_ldos.in")
    isfile(pp) || error("Missing pp_ldos.in: $pp")
    text = read(pp, String)
    m = match(r"(?m)^\s*prefix\s*=\s*'([^']+)'", text)
    m === nothing && error("Could not parse prefix from $pp")
    return m.captures[1]
end

function _paths(dir::String)
    p = _prefix(dir)
    return (
        prefix=p,
        relaxed=joinpath(dir, "$(p)_relaxed.xyz"),
        cube=joinpath(dir, "$(p)_ldos.cube"),
        frame=joinpath(dir, "frame.tsv"),
    )
end

function _xyz_atom_count(path::String)
    isfile(path) || error("Missing XYZ: $path")
    open(path) do io
        return parse(Int, strip(readline(io)))
    end
end

function _default_index_tsv(opt::Options, prefix::String)
    return joinpath(opt.index_dir, "$(prefix)_trimer_indices.tsv")
end

function _frame_args(paths, index_tsv::String, opt::Options)
    if !isempty(opt.origin_indices)
        return (origin_indices=opt.origin_indices, axis_from=opt.axis_from,
                axis_to=opt.axis_to, plane_index=opt.plane_index,
                source="explicit")
    end
    idx = _read_key_tsv(index_tsv)
    for key in ("n_atoms", "origin_indices_bare", "axis_from_index_bare",
                "axis_to_index_bare", "plane_index_bare")
        haskey(idx, key) || error("Index TSV missing $key: $index_tsv")
    end
    bare_n = parse(Int, idx["n_atoms"])
    total_n = _xyz_atom_count(paths.relaxed)
    offset = total_n - bare_n
    offset >= 0 || error("Relaxed XYZ $(paths.relaxed) has fewer atoms ($total_n) than bare trimer ($bare_n)")
    origins = _parse_ints(idx["origin_indices_bare"]) .+ offset
    axis_from = parse(Int, idx["axis_from_index_bare"]) + offset
    axis_to = parse(Int, idx["axis_to_index_bare"]) + offset
    plane_index = parse(Int, idx["plane_index_bare"]) + offset
    return (origin_indices=join(origins, ','), axis_from=axis_from,
            axis_to=axis_to, plane_index=plane_index,
            source="$(index_tsv) offset=$offset")
end

function _check(path::String)
    isfile(path) || error("Missing required QE output: $path")
end

function _run(cmd::Cmd, dry_run::Bool)
    if dry_run
        println("[dry-run] ", cmd)
    else
        run(cmd)
    end
end

function _script(path::String)
    return joinpath(ROOT, path)
end

function main()
    opt = _parse_cli(ARGS)
    glcn = _paths(opt.glcn_dir)
    glcnac = _paths(opt.glcnac_dir)
    glcn_index = isempty(opt.glcn_index_tsv) ? _default_index_tsv(opt, glcn.prefix) : opt.glcn_index_tsv
    glcnac_index = isempty(opt.glcnac_index_tsv) ? _default_index_tsv(opt, glcnac.prefix) : opt.glcnac_index_tsv
    for p in (glcn.relaxed, glcn.cube, glcnac.relaxed, glcnac.cube)
        _check(p)
    end
    println("QE outputs present")
    println("  GlcN relaxed:   ", glcn.relaxed)
    println("  GlcN cube:      ", glcn.cube)
    println("  GlcNAc relaxed: ", glcnac.relaxed)
    println("  GlcNAc cube:    ", glcnac.cube)
    opt.check_only && return

    height = opt.height_nm === nothing ? error("missing height") : opt.height_nm
    for (paths, typ, index_tsv) in ((glcn, 0, glcn_index), (glcnac, 1, glcnac_index))
        frame = _frame_args(paths, index_tsv, opt)
        _run(`$(Base.julia_cmd()) --project=$(ROOT) $(_script("test/extract_qe_mold_frame.jl")) --xyz $(paths.relaxed) --origin-indices $(frame.origin_indices) --axis-from $(frame.axis_from) --axis-to $(frame.axis_to) --plane-index $(frame.plane_index) --height-nm $height --out $(paths.frame)`, opt.dry_run)
        println("Frame ready for type $typ: ", paths.frame, "  (", frame.source, ")")
    end

    _run(`$(Base.julia_cmd()) --project=$(ROOT) $(_script("test/cube_to_stm_maps.jl")) --cube 0:$(glcn.cube) --frame 0:$(glcn.frame) --cube 1:$(glcnac.cube) --frame 1:$(glcnac.frame) --cube-units $(opt.cube_units) --half-nm $(opt.half_nm) --step-nm $(opt.step_nm) --out $(opt.maps_out)`, opt.dry_run)
    _run(`$(Base.julia_cmd()) --project=$(ROOT) $(_script("test/import_stm_mold_maps.jl")) --maps $(opt.maps_out) --out $(opt.templates_out) --bond-out $(opt.bond_out) --half-nm $(opt.half_nm) --step-nm $(opt.step_nm)`, opt.dry_run)

    println("Finalized QE STM molds")
    println("  maps:      ", opt.maps_out)
    println("  templates: ", opt.templates_out)
    println("  bonds:     ", opt.bond_out)
end

main()
