#!/usr/bin/env julia

# Diagnostic constant-current STM mold builder.
#
# Thin CLI over the T1 constant-current library (test/lib/constant_current_cube.jl)
# and the existing frame/map helpers (test/cube_to_stm_maps.jl). For each
# predeclared height in the fixed 0.40:0.05:0.60 nm sensitivity bracket, it
# extracts the first vacuum-side isosurface height map and per-column validity
# mask from the two LDOS cubes, then writes a provenance sidecar binding the
# observable contract (observable, cube hashes, bias, typed Cu reference frames,
# nominal/bracket heights, typed nominal isovalues, z spacing, crossing policy,
# invalid-mask hash, map hash, and the actual nominal isovalue used for each cube
# type).
#
# DIAGNOSTIC ONLY. The output is NOT the production stm_dft_v1 provider. It
# must not enter N_selected, fitting, calibration, thresholds, abstention, or
# the production selector. QE plot_num=5 is a discrete |psi_n(r)|^2 sum, never
# amperes; no nA-to-cube conversion is claimed. See
# .omo/plans/improve-unit-assignment-benchmark.md (T2) and
# docs/src/qe_stm_molds.md.

using Printf

# Reuse the T1 constant-current library AND the existing frame/map helpers
# (_read_cube, _read_frame, _frame, _sample_column_along_normal, _mask_path,
# _fmt, _auto_cc_resolution, _unit_factor, _parse_periodic_axes, FrameSpec,
# CubeSpec, DEFAULT_CC_HEIGHT_MIN/MAX_NM). Including the script is safe: its
# main() is guarded by abspath(PROGRAM_FILE) == @__FILE__. It also brings in
# ScriptUtils (_ensure_parent, _parse_vec3) and the ConstantCurrentCube module.
include(joinpath(@__DIR__, "cube_to_stm_maps.jl"))
using .ConstantCurrentCube: FrameRef, CrossingResult, first_vacuum_crossing,
                            isovalue_for_mean_height, DEFAULT_ISOVALUE_SCAN_INTERVALS

include(joinpath(@__DIR__, "lib", "qe_mold_provenance.jl"))
using .QEMoldProvenance: write_qe_mold_provenance

const DEFAULT_HEIGHTS = [0.40, 0.45, 0.50, 0.55, 0.60]
const DEFAULT_NOMINAL_HEIGHT = 0.50
const DEFAULT_HALF_NM = 0.32      # match the 9x9 production mold patch grid
const DEFAULT_STEP_NM = 0.08
const DEFAULT_BIAS_EV = -0.300
const DEFAULT_CUBE_UNITS = "bohr"
const CC_PROVIDER = "stm_dft_cc_diag"

struct CCCube
    typ::Int
    path::String
    frame::FrameSpec
end

struct CCOptions
    cubes::Vector{CCCube}
    out_prefix::String
    heights::Vector{Float64}
    nominal_height::Float64
    half_nm::Float64
    step_nm::Float64
    cube_units::String
    bias_ev::Float64
    isovalue::Union{Nothing,Float64}     # explicit nominal isovalue (skips calibration)
    z_spacing_nm::Union{Nothing,Float64} # override the auto z-step
    cc_height_min_nm::Float64
    cc_height_max_nm::Float64
    periodic_axes::NTuple{3,Bool}
    provenance_path::Union{Nothing,String}
    isovalue_scan_intervals::Int
end

function _parse_heights(raw::AbstractString)
    tokens = split(strip(raw), ',')
    vals = sort(unique(Float64(parse(Float64, strip(t))) for t in tokens))
    isempty(vals) && error("--heights must list at least one height")
    all(isfinite, vals) || error("--heights must be finite")
    all(>(0), vals) || error("--heights must be positive")
    return vals
end

function _cc_cube(typ::Int, path::AbstractString, frame_args, default_frame::FrameSpec,
                  frame_by_type::Dict{Int,FrameSpec})
    frame = _fill_missing(get(frame_by_type, typ, FrameSpec(nothing, nothing, nothing, nothing)), default_frame)
    frame = _fill_missing(frame, frame_args)
    frame.origin === nothing && error("No origin for cube type=$typ. Pass --origin, --origin$typ, --frame$typ, or --frame $typ:PATH")
    frame.t_axis === nothing && error("No t-axis for cube type=$typ. Pass --t-axis or --frame$typ")
    frame.u_axis === nothing && error("No u-axis for cube type=$typ. Pass --u-axis or --frame$typ")
    isfile(path) || error("Cube file not found ($typ): $path")
    return CCCube(typ, String(path), frame)
end

function _parse_cc_cli(args)
    # cubes are stored as (typ, path) during parsing and bound to frames at the
    # end so --origin/--t-axis/--u-axis/--frameN may appear in any order.
    cube_specs = Tuple{Int,String}[]
    out_prefix = ""
    heights = copy(DEFAULT_HEIGHTS)
    nominal_height = DEFAULT_NOMINAL_HEIGHT
    half_nm = DEFAULT_HALF_NM
    step_nm = DEFAULT_STEP_NM
    cube_units = DEFAULT_CUBE_UNITS
    bias_ev = DEFAULT_BIAS_EV
    isovalue::Union{Nothing,Float64} = nothing
    z_spacing_nm::Union{Nothing,Float64} = nothing
    cc_height_min_nm = DEFAULT_CC_HEIGHT_MIN_NM
    cc_height_max_nm = DEFAULT_CC_HEIGHT_MAX_NM
    periodic_axes = (false, false, false)
    provenance_path::Union{Nothing,String} = nothing
    isovalue_scan_intervals = DEFAULT_ISOVALUE_SCAN_INTERVALS
    default_frame = FrameSpec(nothing, nothing, nothing, nothing)
    frame_by_type = Dict{Int,FrameSpec}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--cube0"; push!(cube_specs, (0, args[i+1])); i += 2
        elseif startswith(arg, "--cube0="); push!(cube_specs, (0, split(arg, '=', limit=2)[2])); i += 1
        elseif arg == "--cube1"; push!(cube_specs, (1, args[i+1])); i += 2
        elseif startswith(arg, "--cube1="); push!(cube_specs, (1, split(arg, '=', limit=2)[2])); i += 1
        elseif arg == "--frame0"; f = _read_frame(args[i+1]); frame_by_type[0] = _fill_missing(f, get(frame_by_type, 0, default_frame)); i += 2
        elseif startswith(arg, "--frame0="); f = _read_frame(split(arg, '=', limit=2)[2]); frame_by_type[0] = _fill_missing(f, get(frame_by_type, 0, default_frame)); i += 1
        elseif arg == "--frame1"; f = _read_frame(args[i+1]); frame_by_type[1] = _fill_missing(f, get(frame_by_type, 1, default_frame)); i += 2
        elseif startswith(arg, "--frame1="); f = _read_frame(split(arg, '=', limit=2)[2]); frame_by_type[1] = _fill_missing(f, get(frame_by_type, 1, default_frame)); i += 1
        elseif arg == "--origin"; default_frame = _fill_missing(FrameSpec(_parse_vec3(args[i+1]), nothing, nothing, nothing), default_frame); i += 2
        elseif startswith(arg, "--origin="); default_frame = _fill_missing(FrameSpec(_parse_vec3(split(arg, '=', limit=2)[2]), nothing, nothing, nothing), default_frame); i += 1
        elseif arg == "--t-axis"; default_frame = _fill_missing(FrameSpec(nothing, _parse_vec3(args[i+1]), nothing, nothing), default_frame); i += 2
        elseif startswith(arg, "--t-axis="); default_frame = _fill_missing(FrameSpec(nothing, _parse_vec3(split(arg, '=', limit=2)[2]), nothing, nothing), default_frame); i += 1
        elseif arg == "--u-axis"; default_frame = _fill_missing(FrameSpec(nothing, nothing, _parse_vec3(args[i+1]), nothing), default_frame); i += 2
        elseif startswith(arg, "--u-axis="); default_frame = _fill_missing(FrameSpec(nothing, nothing, _parse_vec3(split(arg, '=', limit=2)[2]), nothing), default_frame); i += 1
        elseif arg == "--out-prefix"; out_prefix = args[i+1]; i += 2
        elseif startswith(arg, "--out-prefix="); out_prefix = split(arg, '=', limit=2)[2]; i += 1
        elseif arg == "--heights"; heights = _parse_heights(args[i+1]); i += 2
        elseif startswith(arg, "--heights="); heights = _parse_heights(split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--nominal-height"; nominal_height = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--nominal-height="); nominal_height = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--half-nm"; half_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--half-nm="); half_nm = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--step-nm"; step_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--step-nm="); step_nm = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--cube-units"; cube_units = lowercase(strip(args[i+1])); i += 2
        elseif startswith(arg, "--cube-units="); cube_units = lowercase(strip(split(arg, '=', limit=2)[2])); i += 1
        elseif arg == "--bias-ev"; bias_ev = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--bias-ev="); bias_ev = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--isovalue"; isovalue = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--isovalue="); isovalue = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--z-spacing-nm"; z_spacing_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--z-spacing-nm="); z_spacing_nm = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--cc-height-min-nm"; cc_height_min_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--cc-height-min-nm="); cc_height_min_nm = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--cc-height-max-nm"; cc_height_max_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--cc-height-max-nm="); cc_height_max_nm = parse(Float64, split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--periodic-axes"; periodic_axes = _parse_periodic_axes(args[i+1]); i += 2
        elseif startswith(arg, "--periodic-axes="); periodic_axes = _parse_periodic_axes(split(arg, '=', limit=2)[2]); i += 1
        elseif arg == "--provenance"; provenance_path = args[i+1]; i += 2
        elseif startswith(arg, "--provenance="); provenance_path = split(arg, '=', limit=2)[2]; i += 1
        elseif arg == "--isovalue-scan-intervals"; isovalue_scan_intervals = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--isovalue-scan-intervals="); isovalue_scan_intervals = parse(Int, split(arg, '=', limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/build_constant_current_stm_maps.jl [options]

            Diagnostic constant-current STM mold builder (T2). Generates the
            first vacuum-side isosurface height map + validity mask for each
            predeclared height and writes a bound provenance sidecar. NOT the
            production stm_dft_v1 provider.

            Required:
              --cube0 PATH               GlcN (type 0) LDOS cube
              --cube1 PATH               GlcNAc (type 1) LDOS cube
              --out-prefix PREFIX        Output prefix for maps/masks/provenance
              --origin X,Y,Z             Cu reference origin (nm) if no --frame
              --t-axis X,Y,Z             Chain-axis direction if no --frame
              --u-axis X,Y,Z             Transverse direction if no --frame

            Frame (per type, repeatable):
              --frame0 PATH              Frame TSV for type 0 (origin_nm/t_axis/u_axis)
              --frame1 PATH              Frame TSV for type 1

            Options:
              --heights LIST             Comma-separated heights nm [0.40,0.45,0.50,0.55,0.60]
              --nominal-height FLOAT     Nominal height nm [0.50]
              --half-nm FLOAT            Output half-size nm [0.32]
              --step-nm FLOAT            Output grid spacing nm [0.08]
              --cube-units STR           bohr | angstrom | nm [bohr]
              --bias-ev FLOAT            Sample bias eV recorded in provenance [-0.300]
              --isovalue FLOAT           Explicit nominal isovalue (skips calibration)
              --z-spacing-nm FLOAT       Override the auto z-step (nm)
              --cc-height-min-nm FLOAT   Crossing-search lower height [$(DEFAULT_CC_HEIGHT_MIN_NM)]
              --cc-height-max-nm FLOAT   Crossing-search upper height [$(DEFAULT_CC_HEIGHT_MAX_NM)]
              --periodic-axes STR        Must be none for provenance-bound output [none]
              --provenance PATH          Provenance sidecar path [<prefix>.provenance.toml]
              --isovalue-scan-intervals INT  Calibration ambiguity-scan intervals [$(DEFAULT_ISOVALUE_SCAN_INTERVALS)]

            Outputs (per height h, formatted hNNN where NNN=round(h*100)):
              {prefix}_h050.tsv          Map: type/t_nm/u_nm/value (value = height nm | NA)
              {prefix}_h050.mask.tsv     Mask: type/t_nm/u_nm/status
              {prefix}.provenance.toml   Bound provenance (nominal + bracket artifacts)

            The nominal height (0.50 nm) fixes the isovalue by the mean-height
            policy above the Cu reference unless --isovalue is given. Each bracket
            height calibrates its own isovalue. Missing first crossings stay NA
            in the map and the mask records the reason; they are never zero-filled.
            No nA-to-cube conversion is claimed.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    # Bind parsed frames to cubes now that all --origin/--frameN are known.
    cubes = CCCube[_cc_cube(typ, path, FrameSpec(nothing, nothing, nothing, nothing),
                            default_frame, frame_by_type) for (typ, path) in cube_specs]
    isempty(out_prefix) && error("--out-prefix is required")
    isempty(cubes) && error("Pass at least --cube0 and --cube1")
    length(cubes) == 2 || error("Exactly two cubes (type 0 and type 1) are required")
    Set(c.typ for c in cubes) == Set([0, 1]) || error("Cubes must be type 0 and type 1")
    half_nm > 0 || error("--half-nm must be positive")
    step_nm > 0 || error("--step-nm must be positive")
    cube_units in ("bohr", "angstrom", "a", "nm") || error("--cube-units must be bohr, angstrom, or nm")
    nominal_height in heights || error("--nominal-height ($nominal_height) must be in --heights ($heights)")
    isovalue === nothing || isfinite(isovalue) || error("--isovalue must be finite")
    z_spacing_nm === nothing || (z_spacing_nm > 0 && isfinite(z_spacing_nm)) || error("--z-spacing-nm must be positive finite")
    (isfinite(cc_height_min_nm) && isfinite(cc_height_max_nm) &&
     cc_height_max_nm > cc_height_min_nm) ||
        error("--cc-height-max-nm must exceed finite --cc-height-min-nm")
    any(periodic_axes) &&
        error("--periodic-axes must be none for provenance-bound constant-current maps")
    isovalue_scan_intervals > 0 ||
        error("--isovalue-scan-intervals must be a positive integer")
    return CCOptions(cubes, out_prefix, heights, nominal_height, half_nm, step_nm,
                     cube_units, bias_ev, isovalue, z_spacing_nm,
                     cc_height_min_nm, cc_height_max_nm,
                     periodic_axes, provenance_path, isovalue_scan_intervals)
end

_height_tag(h::Real) = @sprintf("%03d", round(Int, h * 100))
_map_path(prefix::AbstractString, h::Real) = string(prefix, "_h", _height_tag(h), ".tsv")

function _type_frame_entries(cubes::Vector{CCCube})
    return [(type=c.typ, origin_nm=c.frame.origin, t_axis=c.frame.t_axis,
             u_axis=c.frame.u_axis) for c in sort(cubes, by=c -> c.typ)]
end

function _atomic_write_rows(path::AbstractString, header::Vector{String}, rows::Vector{String})
    parent = dirname(abspath(path))
    isdir(parent) || error("output directory does not exist: $parent")
    tmp, io = mktemp(parent)
    committed = false
    try
        println(io, join(header, '\t'))
        isempty(rows) || println(io, join(rows, '\n'))
        flush(io)
        close(io)
        Base.Filesystem.rename(tmp, path)
        committed = true
    finally
        isopen(io) && close(io)
        !committed && ispath(tmp) && rm(tmp; force=true)
    end
    return path
end

# Build the constant-current map+mask for one cube at one height. Returns the
# grid rows (map, mask) and the isovalue actually used with its source label.
function _build_one_height(typ::Int, cube_path::AbstractString, frame::FrameSpec,
                           height::Float64, coords::Vector{Float64}, opt::CCOptions)
    cube = _read_cube(cube_path, opt.cube_units)
    that, uhat, nhat = _frame(frame.t_axis, frame.u_axis)
    s_step = opt.z_spacing_nm === nothing ? _auto_cc_resolution(cube) : opt.z_spacing_nm
    s_values = collect(opt.cc_height_min_nm:s_step:opt.cc_height_max_nm)
    length(s_values) >= 2 ||
        error("constant-current search range produces fewer than 2 samples; widen the cube z extent")
    if opt.isovalue !== nothing && height == opt.nominal_height
        iso = Float64(opt.isovalue); iso_source = "explicit"
    else
        frame_ref = FrameRef(frame.origin, frame.t_axis, frame.u_axis)
        iso = Float64(isovalue_for_mean_height(cube, frame_ref, height;
            scan_intervals=opt.isovalue_scan_intervals)); iso_source = "calibrated"
    end
    map_rows = String[]
    mask_rows = String[]
    counts = Dict{Symbol,Int}(:found => 0, :absent => 0, :ambiguous => 0,
                              :nonfinite => 0, :insufficient_z => 0)
    for u in coords, t in coords
        z_col, v_col = _sample_column_along_normal(cube, frame.origin, that, uhat, nhat,
                                                   t, u, s_values, opt.periodic_axes)
        result = isempty(z_col) ? CrossingResult(NaN, 0, :nonfinite, 0) :
                                first_vacuum_crossing(z_col, v_col, iso)
        counts[result.status] = get(counts, result.status, 0) + 1
        push!(map_rows, join([typ, @sprintf("%.8g", t), @sprintf("%.8g", u), _fmt(result.z)], '\t'))
        push!(mask_rows, join([typ, @sprintf("%.8g", t), @sprintf("%.8g", u), result.status], '\t'))
    end
    return map_rows, mask_rows, iso, iso_source, s_step, counts
end

# Wrapper exposing the per-type constant-current extraction for the synthetic
# gate tests (determinism, abstention, complement symmetry) without subprocess.
function _extract_surface_grid(cube_path::AbstractString, frame_origin, frame_t_axis,
                               frame_u_axis, height_nm::Real, half_nm::Real,
                               step_nm::Real, cube_units::AbstractString,
                               cc_height_min_nm::Real, cc_height_max_nm::Real,
                               z_spacing_nm::Union{Nothing,Real},
                               periodic_axes::NTuple{3,Bool};
                               isovalue_scan_intervals=DEFAULT_ISOVALUE_SCAN_INTERVALS)
    cube = _read_cube(cube_path, cube_units)
    that, uhat, nhat = _frame(frame_t_axis, frame_u_axis)
    s_step = z_spacing_nm === nothing ? _auto_cc_resolution(cube) : Float64(z_spacing_nm)
    s_values = collect(Float64(cc_height_min_nm):s_step:Float64(cc_height_max_nm))
    length(s_values) >= 2 || error("constant-current z range has fewer than 2 samples")
    frame_ref = FrameRef(frame_origin, frame_t_axis, frame_u_axis)
    iso = Float64(isovalue_for_mean_height(cube, frame_ref, height_nm;
        scan_intervals=isovalue_scan_intervals))
    coords = collect(-half_nm:step_nm:half_nm)
    heights = Float64[]; statuses = Symbol[]
    for u in coords, t in coords
        z_col, v_col = _sample_column_along_normal(cube, frame_origin, that, uhat, nhat,
                                                   t, u, s_values, periodic_axes)
        result = isempty(z_col) ? CrossingResult(NaN, 0, :nonfinite, 0) :
                                first_vacuum_crossing(z_col, v_col, iso)
        push!(heights, result.z)
        push!(statuses, result.status)
    end
    return coords, heights, statuses, iso, s_step
end

function _builder_output_set(opt::CCOptions)
    artifacts = [(height, _map_path(opt.out_prefix, height),
        _mask_path(_map_path(opt.out_prefix, height))) for height in opt.heights]
    provenance = opt.provenance_path === nothing ?
        string(opt.out_prefix, ".provenance.toml") : opt.provenance_path
    destinations = String[]
    for (_, map_path, mask_path) in artifacts
        push!(destinations, map_path, mask_path)
    end
    push!(destinations, provenance)
    absolute = abspath.(destinations)
    length(unique(absolute)) == length(absolute) ||
        error("constant-current output transaction destinations must be distinct")
    length(unique(dirname.(absolute))) == 1 ||
        error("constant-current map, mask, and provenance outputs must share one directory")
    _ensure_parent(first(destinations))
    return (; artifacts, provenance, destinations, absolute)
end

function _recover_builder_output_set(output_set)
    ordered = output_set.absolute
    gate = abspath(output_set.provenance)
    marker = _cc_txn_marker(gate)
    if _cc_txn_present(marker)
        id, _ = _cc_txn_read_marker(marker, _cc_txn_signature(ordered))
        for destination in ordered
            stage = _cc_txn_stage(destination, id)
            absent = _cc_txn_absent(destination, id)
            islink(stage) && error("refusing symlinked transaction stage sidecar: $stage")
            islink(absent) && error("refusing symlinked transaction absence sidecar: $absent")
            _cc_txn_present(stage) && !isfile(stage) &&
                error("invalid transaction stage sidecar: $stage")
            _cc_txn_present(absent) && !isfile(absent) &&
                error("invalid transaction absence sidecar: $absent")
            backup = _cc_txn_backup(destination, id)
            _cc_txn_present(backup) && !(isfile(backup) || islink(backup)) &&
                error("invalid transaction backup sidecar: $backup")
        end
    end
    _cc_txn_recover(ordered, gate)
    return nothing
end

function _run_build(opt::CCOptions; failpoint=nothing, interruptpoint=nothing)
    output_set = _builder_output_set(opt)
    _recover_builder_output_set(output_set)
    coords = collect(-opt.half_nm:opt.step_nm:opt.half_nm)
    built = NamedTuple[]
    nominal_map = nominal_mask = ""
    nominal_zsp = 0.0
    nominal_isovalues = Dict{Int,Tuple{Float64,String}}()
    for (height, map_path, mask_path) in output_set.artifacts
        map_rows = String[]
        mask_rows = String[]
        height_zsp::Union{Nothing,Float64} = nothing
        for spec in sort(opt.cubes, by=c -> c.typ)
            rows, mrows, iso, src, s_step, counts =
                _build_one_height(spec.typ, spec.path, spec.frame, height, coords, opt)
            append!(map_rows, rows); append!(mask_rows, mrows)
            height_zsp === nothing || isapprox(s_step, height_zsp; rtol=1e-12, atol=0.0) ||
                error("constant-current cubes must use the same z spacing; provenance records one shared value")
            height_zsp = s_step
            println("  type $(spec.typ) h=$(height) isovalue: $iso ($src); z_step=$s_step")
            println("    " * join(["$(String(k))=$(get(counts,k,0))" for k in keys(counts)], ' '))
            if height == opt.nominal_height
                nominal_isovalues[spec.typ] = (iso, src)
                nominal_zsp = s_step
            end
        end
        push!(built, (; height, map_path, mask_path, map_rows, mask_rows))
        height == opt.nominal_height && (nominal_map = map_path; nominal_mask = mask_path)
    end
    bracket_artifacts = [(item.height, item.map_path, item.mask_path)
        for item in built if item.height != opt.nominal_height]
    Set(keys(nominal_isovalues)) == Set([0, 1]) ||
        error("internal error: nominal isovalues were not recorded for both cube types")
    type_isovalues = [(type=typ, isovalue=nominal_isovalues[typ][1],
                        source=nominal_isovalues[typ][2]) for typ in (0, 1)]
    cube0 = [c for c in opt.cubes if c.typ == 0][1].path
    cube1 = [c for c in opt.cubes if c.typ == 1][1].path
    provenance = output_set.provenance
    _with_output_transaction(output_set.destinations; gate=provenance,
            failpoint=failpoint, interruptpoint=interruptpoint) do staged
        artifact_sources = Dict{String,String}()
        for item in built
            staged_map = staged[abspath(item.map_path)]
            staged_mask = staged[abspath(item.mask_path)]
            _atomic_write_rows(staged_map,
                ["type", "t_nm", "u_nm", "value"], item.map_rows)
            _atomic_write_rows(staged_mask,
                ["type", "t_nm", "u_nm", "status"], item.mask_rows)
            artifact_sources[item.map_path] = staged_map
            artifact_sources[item.mask_path] = staged_mask
        end
        write_qe_mold_provenance(staged[abspath(provenance)];
            provider=CC_PROVIDER, glcn_cube=cube0, glcnac_cube=cube1,
            maps=nominal_map, sample_bias_ev=opt.bias_ev, height_nm=opt.nominal_height,
            half_nm=opt.half_nm, step_nm=opt.step_nm, cube_units=opt.cube_units,
            observable="constant-current", nominal_height_nm=opt.nominal_height,
            bracket_heights_nm=opt.heights, type_isovalues=type_isovalues,
            z_spacing_nm=nominal_zsp, crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=nominal_mask, type_frames=_type_frame_entries(opt.cubes),
            bracket_artifacts=bracket_artifacts,
            isovalue_scan_intervals=opt.isovalue_scan_intervals,
            artifact_sources=artifact_sources)
    end
    println("Built diagnostic constant-current maps")
    println("  cubes:      $cube0, $cube1")
    println("  prefix:     $(opt.out_prefix)")
    println("  nominal:    h=$(opt.nominal_height) nm  -> $nominal_map / $nominal_mask")
    println("  heights:    " * join(["h=$h" for h in opt.heights], ", "))
    println("  provenance: $provenance")
    println("  isovalue scan intervals: $(opt.isovalue_scan_intervals)")
    println("  observable: constant-current (diagnostic; NOT stm_dft_v1)")
    return nothing
end

function _run_build(; failpoint=nothing, interruptpoint=nothing)
    return _run_build(_parse_cc_cli(ARGS);
        failpoint=failpoint, interruptpoint=interruptpoint)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run_build()
end
