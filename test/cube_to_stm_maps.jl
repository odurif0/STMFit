#!/usr/bin/env julia

# Convert Gaussian cube maps (for example Quantum ESPRESSO pp.x output) into the
# long STM/LDOS map TSV consumed by import_stm_mold_maps.jl.
#
# The script samples one or two cube files in an aligned local frame:
#   r(t,u) = origin + t * t_hat + u * u_hat + height * normal_hat
# and writes rows:
#   type  t_nm  u_nm  value
# where type=0 is GlcN and type=1 is GlcNAc.

using Printf
using LinearAlgebra
using SHA

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _parse_vec3

include(joinpath(@__DIR__, "lib", "constant_current_cube.jl"))
using .ConstantCurrentCube: CubeGrid, FrameRef, CrossingResult, SurfaceResult,
                            first_vacuum_crossing, constant_current_surface,
                            isovalue_for_mean_height,
                            DEFAULT_ISOVALUE_SCAN_INTERVALS

const DEFAULT_OUT = "templates/chitosan_stm_maps.tsv"
const DEFAULT_OBSERVABLE = "constant-height"
const SUPPORTED_OBSERVABLES = ("constant-height", "constant-current")
const DEFAULT_MEAN_HEIGHT_NM = 0.50
const DEFAULT_CC_HEIGHT_MIN_NM = -0.50
const DEFAULT_CC_HEIGHT_MAX_NM = 2.00

struct FrameSpec
    origin::Union{Nothing,Vector{Float64}}
    t_axis::Union{Nothing,Vector{Float64}}
    u_axis::Union{Nothing,Vector{Float64}}
    height_nm::Union{Nothing,Float64}
end

struct CubeSpec
    typ::Int
    path::String
    frame::FrameSpec
end

struct Options
    cubes::Vector{CubeSpec}
    out_tsv::String
    half_nm::Float64
    step_nm::Float64
    cube_units::String
    periodic_axes::NTuple{3,Bool}
    observable::String
    isovalue::Union{Nothing,Float64}
    mean_height_nm::Float64
    cc_height_min_nm::Float64
    cc_height_max_nm::Float64
    cc_height_resolution_nm::Union{Nothing,Float64}
    isovalue_scan_intervals::Int
end

function _parse_periodic_axes(raw::AbstractString)
    normalized = replace(lowercase(strip(raw)), "," => "")
    normalized == "none" && return (false, false, false)
    isempty(normalized) && throw(ArgumentError("--periodic-axes cannot be empty"))
    all(character -> character in ('x', 'y'), normalized) ||
        throw(ArgumentError("--periodic-axes supports diagnostic lateral axes x and y only"))
    length(unique(normalized)) == length(normalized) ||
        throw(ArgumentError("--periodic-axes contains a duplicate axis"))
    return ('x' in normalized, 'y' in normalized, false)
end

function _read_frame(path::String)
    isfile(path) || error("Frame TSV not found: $path")
    d = Dict{String,String}()
    for line in readlines(path)
        t = strip(line)
        isempty(t) && continue
        startswith(t, '#') && continue
        parts = split(t, '\t'; limit=2)
        length(parts) == 2 || continue
        d[strip(parts[1])] = strip(parts[2])
    end
    haskey(d, "origin_nm") || error("Frame TSV missing origin_nm: $path")
    haskey(d, "t_axis") || error("Frame TSV missing t_axis: $path")
    haskey(d, "u_axis") || error("Frame TSV missing u_axis: $path")
    return FrameSpec(
        _parse_vec3(d["origin_nm"]),
        _parse_vec3(d["t_axis"]),
        _parse_vec3(d["u_axis"]),
        haskey(d, "height_nm") ? parse(Float64, d["height_nm"]) : nothing,
    )
end

function _fill_missing(primary::FrameSpec, fallback::FrameSpec)
    return FrameSpec(
        primary.origin === nothing ? fallback.origin : primary.origin,
        primary.t_axis === nothing ? fallback.t_axis : primary.t_axis,
        primary.u_axis === nothing ? fallback.u_axis : primary.u_axis,
        primary.height_nm === nothing ? fallback.height_nm : primary.height_nm,
    )
end

function _apply_frame!(cubes::Vector{CubeSpec}, typ::Union{Nothing,Int}, frame::FrameSpec)
    for i in eachindex(cubes)
        if typ === nothing || cubes[i].typ == typ
            cubes[i] = CubeSpec(cubes[i].typ, cubes[i].path, _fill_missing(cubes[i].frame, frame))
        end
    end
end

function _parse_frame_arg(s::AbstractString)
    if startswith(s, "0:") || startswith(s, "1:")
        return parse(Int, s[1:1]), String(s[3:end])
    end
    return nothing, String(s)
end

function _typed_frame!(cubes::Vector{CubeSpec}, frame_by_type::Dict{Int,FrameSpec}, typ::Int, frame::FrameSpec)
    frame_by_type[typ] = _fill_missing(frame, get(frame_by_type, typ, FrameSpec(nothing, nothing, nothing, nothing)))
    _apply_frame!(cubes, typ, frame)
end

function _parse_cube_spec(s::AbstractString, default_frame::FrameSpec, frame_by_type::Dict{Int,FrameSpec})
    sep = findfirst(==(':'), s)
    sep === nothing && error("--cube expects TYPE:PATH or TYPE:PATH:ORIGIN")
    typ = round(Int, parse(Float64, strip(s[1:prevind(s, sep)])))
    typ in (0, 1) || error("cube type must be 0 or 1")
    rest = String(s[nextind(s, sep):end])
    path = strip(rest)
    frame = _fill_missing(get(frame_by_type, typ, FrameSpec(nothing, nothing, nothing, nothing)), default_frame)
    last_sep = findlast(==(':'), rest)
    if last_sep !== nothing && last_sep < lastindex(rest)
        maybe_origin = strip(rest[nextind(rest, last_sep):end])
        parsed_origin = try
            _parse_vec3(maybe_origin)
        catch
            nothing
        end
        if parsed_origin !== nothing
            path = strip(rest[1:prevind(rest, last_sep)])
            frame = _fill_missing(FrameSpec(parsed_origin, nothing, nothing, nothing), frame)
        end
    end
    isempty(path) && error("empty cube path in --cube spec: $s")
    return CubeSpec(typ, path, frame)
end

function _parse_cli(args)
    cubes = CubeSpec[]
    frame_by_type = Dict{Int,FrameSpec}()
    out_tsv = DEFAULT_OUT
    default_frame = FrameSpec(nothing, nothing, nothing, nothing)
    half_nm = 0.48
    step_nm = 0.08
    cube_units = "bohr"
    periodic_axes = (false, false, false)
    observable = DEFAULT_OBSERVABLE
    isovalue::Union{Nothing,Float64} = nothing
    mean_height_nm = DEFAULT_MEAN_HEIGHT_NM
    cc_height_min_nm = DEFAULT_CC_HEIGHT_MIN_NM
    cc_height_max_nm = DEFAULT_CC_HEIGHT_MAX_NM
    cc_height_resolution_nm::Union{Nothing,Float64} = nothing
    isovalue_scan_intervals = DEFAULT_ISOVALUE_SCAN_INTERVALS
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--cube"
            push!(cubes, _parse_cube_spec(args[i+1], default_frame, frame_by_type)); i += 2
        elseif startswith(arg, "--cube=")
            push!(cubes, _parse_cube_spec(split(arg, "=", limit=2)[2], default_frame, frame_by_type)); i += 1
        elseif arg == "--cube0"
            frame = _fill_missing(get(frame_by_type, 0, FrameSpec(nothing, nothing, nothing, nothing)), default_frame)
            push!(cubes, CubeSpec(0, args[i+1], frame)); i += 2
        elseif startswith(arg, "--cube0=")
            frame = _fill_missing(get(frame_by_type, 0, FrameSpec(nothing, nothing, nothing, nothing)), default_frame)
            push!(cubes, CubeSpec(0, split(arg, "=", limit=2)[2], frame)); i += 1
        elseif arg == "--cube1"
            frame = _fill_missing(get(frame_by_type, 1, FrameSpec(nothing, nothing, nothing, nothing)), default_frame)
            push!(cubes, CubeSpec(1, args[i+1], frame)); i += 2
        elseif startswith(arg, "--cube1=")
            frame = _fill_missing(get(frame_by_type, 1, FrameSpec(nothing, nothing, nothing, nothing)), default_frame)
            push!(cubes, CubeSpec(1, split(arg, "=", limit=2)[2], frame)); i += 1
        elseif arg == "--origin"
            frame = FrameSpec(_parse_vec3(args[i+1]), nothing, nothing, nothing)
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 2
        elseif startswith(arg, "--origin=")
            frame = FrameSpec(_parse_vec3(split(arg, "=", limit=2)[2]), nothing, nothing, nothing)
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 1
        elseif arg == "--origin0"
            _typed_frame!(cubes, frame_by_type, 0, FrameSpec(_parse_vec3(args[i+1]), nothing, nothing, nothing)); i += 2
        elseif startswith(arg, "--origin0=")
            _typed_frame!(cubes, frame_by_type, 0, FrameSpec(_parse_vec3(split(arg, "=", limit=2)[2]), nothing, nothing, nothing)); i += 1
        elseif arg == "--origin1"
            _typed_frame!(cubes, frame_by_type, 1, FrameSpec(_parse_vec3(args[i+1]), nothing, nothing, nothing)); i += 2
        elseif startswith(arg, "--origin1=")
            _typed_frame!(cubes, frame_by_type, 1, FrameSpec(_parse_vec3(split(arg, "=", limit=2)[2]), nothing, nothing, nothing)); i += 1
        elseif arg == "--frame"
            typ, path = _parse_frame_arg(args[i+1])
            frame = _read_frame(path)
            if typ === nothing
                default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame)
            else
                _typed_frame!(cubes, frame_by_type, typ, frame)
            end
            i += 2
        elseif startswith(arg, "--frame=")
            typ, path = _parse_frame_arg(split(arg, "=", limit=2)[2])
            frame = _read_frame(path)
            if typ === nothing
                default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame)
            else
                _typed_frame!(cubes, frame_by_type, typ, frame)
            end
            i += 1
        elseif arg == "--frame0"
            frame = _read_frame(args[i+1])
            _typed_frame!(cubes, frame_by_type, 0, frame); i += 2
        elseif startswith(arg, "--frame0=")
            frame = _read_frame(split(arg, "=", limit=2)[2])
            _typed_frame!(cubes, frame_by_type, 0, frame); i += 1
        elseif arg == "--frame1"
            frame = _read_frame(args[i+1])
            _typed_frame!(cubes, frame_by_type, 1, frame); i += 2
        elseif startswith(arg, "--frame1=")
            frame = _read_frame(split(arg, "=", limit=2)[2])
            _typed_frame!(cubes, frame_by_type, 1, frame); i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--t-axis"
            frame = FrameSpec(nothing, _parse_vec3(args[i+1]), nothing, nothing)
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 2
        elseif startswith(arg, "--t-axis=")
            frame = FrameSpec(nothing, _parse_vec3(split(arg, "=", limit=2)[2]), nothing, nothing)
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 1
        elseif arg == "--u-axis"
            frame = FrameSpec(nothing, nothing, _parse_vec3(args[i+1]), nothing)
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 2
        elseif startswith(arg, "--u-axis=")
            frame = FrameSpec(nothing, nothing, _parse_vec3(split(arg, "=", limit=2)[2]), nothing)
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 1
        elseif arg == "--height-nm"
            frame = FrameSpec(nothing, nothing, nothing, parse(Float64, args[i+1]))
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 2
        elseif startswith(arg, "--height-nm=")
            frame = FrameSpec(nothing, nothing, nothing, parse(Float64, split(arg, "=", limit=2)[2]))
            default_frame = _fill_missing(frame, default_frame); _apply_frame!(cubes, nothing, frame); i += 1
        elseif arg == "--half-nm"; half_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--half-nm="); half_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--step-nm"; step_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--step-nm="); step_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cube-units"; cube_units = lowercase(strip(args[i+1])); i += 2
        elseif startswith(arg, "--cube-units="); cube_units = lowercase(strip(split(arg, "=", limit=2)[2])); i += 1
        elseif arg == "--periodic-axes"; periodic_axes = _parse_periodic_axes(args[i+1]); i += 2
        elseif startswith(arg, "--periodic-axes="); periodic_axes = _parse_periodic_axes(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--observable"; observable = lowercase(strip(args[i+1])); i += 2
        elseif startswith(arg, "--observable="); observable = lowercase(strip(split(arg, "=", limit=2)[2])); i += 1
        elseif arg == "--isovalue"; isovalue = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--isovalue="); isovalue = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--mean-height-nm"; mean_height_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--mean-height-nm="); mean_height_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cc-height-min-nm"; cc_height_min_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--cc-height-min-nm="); cc_height_min_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cc-height-max-nm"; cc_height_max_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--cc-height-max-nm="); cc_height_max_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cc-height-resolution-nm"; cc_height_resolution_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--cc-height-resolution-nm="); cc_height_resolution_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--isovalue-scan-intervals"; isovalue_scan_intervals = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--isovalue-scan-intervals="); isovalue_scan_intervals = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/cube_to_stm_maps.jl [options]

            Required:
              --cube TYPE:PATH[:ORIGIN]  Cube map for type 0 or 1. Repeatable.
              --t-axis X,Y,Z             Chain-axis direction, unless frames provide it
              --u-axis X,Y,Z             Transverse direction, unless frames provide it

            Options:
              --origin X,Y,Z             Default central-lobe origin in nm
              --frame PATH               Default frame TSV from extract_qe_mold_frame.jl
              --frame TYPE:PATH          Per-type frame TSV, repeatable
              --frame0 PATH              Alias for --frame 0:PATH
              --frame1 PATH              Alias for --frame 1:PATH
              --out PATH                 Output map TSV [$(DEFAULT_OUT)]
              --height-nm FLOAT          Offset along normal t x u [0]
              --half-nm FLOAT            Output half-size [0.48]
              --step-nm FLOAT            Output grid spacing [0.08]
              --cube-units STR           bohr | angstrom | nm [bohr]
              --periodic-axes STR        Diagnostic wrapping: none | x | y | xy [none]
              --observable STR           constant-height | constant-current [$(DEFAULT_OBSERVABLE)]
              --isovalue FLOAT           Explicit LDOS isovalue for constant-current
              --mean-height-nm FLOAT     Target mean height for isovalue calibration [$(DEFAULT_MEAN_HEIGHT_NM)]
              --cc-height-min-nm FLOAT   Constant-current search minimum height [$(DEFAULT_CC_HEIGHT_MIN_NM)]
              --cc-height-max-nm FLOAT   Constant-current search maximum height [$(DEFAULT_CC_HEIGHT_MAX_NM)]
              --cc-height-resolution-nm FLOAT  Constant-current z-step [auto from cube]
              --isovalue-scan-intervals INT    Calibration ambiguity-scan intervals [$(DEFAULT_ISOVALUE_SCAN_INTERVALS)]

            Examples:
              julia --project=. test/cube_to_stm_maps.jl \
                --cube 0:glcn_ldos.cube:1.2,1.1,2.0 \
                --cube 1:glcnac_ldos.cube:1.2,1.1,2.0 \
                --t-axis 1,0,0 --u-axis 0,1,0 --height-nm 0.35

              julia --project=. test/cube_to_stm_maps.jl \
                --cube 0:glcn_ldos.cube --frame 0:qe/glcn/frame.tsv \
                --cube 1:glcnac_ldos.cube --frame 1:qe/glcnac/frame.tsv

            The cube values are sampled in an aligned local frame and written as
            type/t_nm/u_nm/value rows for import_stm_mold_maps.jl. Truth labels
            and benchmark composition are not used. Periodic sampling is
            diagnostic-only, wraps FFT-grid coordinates, and never wraps z.

            --observable constant-height (default) writes the trilinearly
            interpolated cube value at the height-nm plane above each frame
            origin. Byte-identical to the historical behavior.

            --observable constant-current writes the first vacuum-side
            isosurface height (nm above the frame origin along the normal) at
            each (t,u). The isovalue is calibrated from the 0.50 nm mean-height
            policy unless --isovalue is given; a sidecar .mask.tsv records
            per-column status. QE plot_num=5 is a discrete |psi|^2 sum, never
            amperes; no nA-to-cube conversion is claimed.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(cubes) && error("Pass at least one --cube TYPE:PATH[:ORIGIN]")
    for i in eachindex(cubes)
        frame = _fill_missing(
            cubes[i].frame,
            _fill_missing(
                get(frame_by_type, cubes[i].typ, FrameSpec(nothing, nothing, nothing, nothing)),
                _fill_missing(default_frame, FrameSpec(nothing, nothing, nothing, 0.0)),
            ),
        )
        frame.origin === nothing && error("No origin for cube type=$(cubes[i].typ). Pass --origin, --origin$(cubes[i].typ), --frame, --frame $(cubes[i].typ):PATH, or TYPE:PATH:x,y,z")
        frame.t_axis === nothing && error("No t-axis for cube type=$(cubes[i].typ). Pass --t-axis or --frame $(cubes[i].typ):PATH")
        frame.u_axis === nothing && error("No u-axis for cube type=$(cubes[i].typ). Pass --u-axis or --frame $(cubes[i].typ):PATH")
        frame.height_nm === nothing && error("No height for cube type=$(cubes[i].typ). Pass --height-nm or include height_nm in --frame")
        cubes[i] = CubeSpec(cubes[i].typ, cubes[i].path, frame)
    end
    for c in cubes
        isfile(c.path) || error("Cube file not found: $(c.path)")
    end
    cube_units in ("bohr", "angstrom", "a", "nm") || error("--cube-units must be bohr, angstrom, or nm")
    half_nm > 0 || error("--half-nm must be positive")
    step_nm > 0 || error("--step-nm must be positive")
    isovalue_scan_intervals > 0 || error("--isovalue-scan-intervals must be a positive integer")
    observable in SUPPORTED_OBSERVABLES ||
        error("--observable must be one of: $(join(SUPPORTED_OBSERVABLES, ", ")); got $observable")
    if observable == "constant-current"
        cc_height_max_nm > cc_height_min_nm ||
            error("--cc-height-max-nm must exceed --cc-height-min-nm")
        cc_height_resolution_nm === nothing || cc_height_resolution_nm > 0 ||
            error("--cc-height-resolution-nm must be positive")
        isovalue === nothing || isfinite(isovalue) ||
            error("--isovalue must be finite")
        mean_height_nm > 0 || error("--mean-height-nm must be positive")
    end
    return Options(cubes, out_tsv, half_nm, step_nm, cube_units, periodic_axes,
                   observable, isovalue, mean_height_nm,
                   cc_height_min_nm, cc_height_max_nm, cc_height_resolution_nm,
                   isovalue_scan_intervals)
end

function _unit_factor(units::String)
    units == "bohr" && return 0.0529177210903
    units in ("angstrom", "a") && return 0.1
    return 1.0
end

function _read_cube(path::String, units::String)
    lines = readlines(path)
    length(lines) >= 6 || error("Cube file too short: $path")
    scale = _unit_factor(units)
    fields = split(strip(lines[3]))
    length(fields) >= 4 || error("Bad cube atom-count/origin line")
    nat = abs(parse(Int, fields[1]))
    origin = [parse(Float64, fields[i]) * scale for i in 2:4]
    n = Int[]
    axes = zeros(3, 3)
    for ax in 1:3
        parts = split(strip(lines[3 + ax]))
        length(parts) >= 4 || error("Bad cube axis line $ax")
        push!(n, abs(parse(Int, parts[1])))
        axes[:, ax] .= [parse(Float64, parts[i]) * scale for i in 2:4]
    end
    value_start = 6 + nat + 1
    vals = Float64[]
    for line in lines[value_start:end]
        for tok in split(strip(line))
            isempty(tok) || push!(vals, parse(Float64, replace(tok, 'D' => 'E')))
        end
    end
    expected = prod(n)
    length(vals) < expected && error("Cube has $(length(vals)) values, expected $expected")
    length(vals) > expected && (vals = vals[1:expected])
    return CubeGrid(origin, axes, (n[1], n[2], n[3]), vals)
end

function _cube_value(c::CubeGrid, ix::Int, iy::Int, iz::Int)
    nx, ny, nz = c.n
    return c.values[((ix - 1) * ny + (iy - 1)) * nz + iz]
end

function _sample_cube(c::CubeGrid, r::Vector{Float64})
    q = c.axes \ (r .- c.origin)  # zero-based fractional grid coordinates
    nx, ny, nz = c.n
    if any(q .< 0) || q[1] > nx - 1 || q[2] > ny - 1 || q[3] > nz - 1
        return NaN
    end
    ix = floor(Int, q[1]) + 1
    iy = floor(Int, q[2]) + 1
    iz = floor(Int, q[3]) + 1
    ix == nx && (ix -= 1)
    iy == ny && (iy -= 1)
    iz == nz && (iz -= 1)
    (ix < 1 || iy < 1 || iz < 1 || ix >= nx || iy >= ny || iz >= nz) && return NaN
    tx = q[1] - (ix - 1)
    ty = q[2] - (iy - 1)
    tz = q[3] - (iz - 1)
    v = 0.0
    for dx in 0:1, dy in 0:1, dz in 0:1
        w = (dx == 1 ? tx : 1 - tx) * (dy == 1 ? ty : 1 - ty) * (dz == 1 ? tz : 1 - tz)
        v += w * _cube_value(c, ix + dx, iy + dy, iz + dz)
    end
    return v
end

function _sample_cube_periodic(c::CubeGrid, r::Vector{Float64},
                               periodic_axes::NTuple{3,Bool})
    periodic_axes[3] && throw(ArgumentError("periodic z sampling is not allowed"))
    q = c.axes \ (r .- c.origin)
    sizes = c.n
    for axis in 1:3
        if periodic_axes[axis]
            q[axis] = mod(q[axis], sizes[axis])
        elseif q[axis] < 0 || q[axis] > sizes[axis] - 1
            return NaN
        end
    end

    lower = [floor(Int, q[axis]) + 1 for axis in 1:3]
    for axis in 1:3
        if !periodic_axes[axis] && lower[axis] == sizes[axis]
            lower[axis] -= 1
        end
        if !periodic_axes[axis] && (lower[axis] < 1 || lower[axis] >= sizes[axis])
            return NaN
        end
    end
    fractions = [q[axis] - (lower[axis] - 1) for axis in 1:3]
    value = 0.0
    for dx in 0:1, dy in 0:1, dz in 0:1
        offsets = (dx, dy, dz)
        indices = ntuple(3) do axis
            candidate = lower[axis] + offsets[axis]
            periodic_axes[axis] ? mod1(candidate, sizes[axis]) : candidate
        end
        weight = prod(offsets[axis] == 1 ? fractions[axis] : 1 - fractions[axis]
            for axis in 1:3)
        value += weight * _cube_value(c, indices...)
    end
    return value
end

function _frame(t_axis::Vector{Float64}, u_axis::Vector{Float64})
    norm(t_axis) > 0 || error("--t-axis must be nonzero")
    that = t_axis ./ norm(t_axis)
    u = u_axis .- dot(u_axis, that) .* that
    norm(u) > 0 || error("--u-axis is collinear with --t-axis")
    uhat = u ./ norm(u)
    nhat = cross(that, uhat)
    norm(nhat) > 0 || error("invalid frame")
    nhat ./= norm(nhat)
    return that, uhat, nhat
end

_fmt(v) = isfinite(v) ? @sprintf("%.10g", v) : "NA"

function _mask_path(out_tsv::AbstractString)
    dir, base = splitdir(String(out_tsv))
    stem, ext = occursin('.', base) ? splitext(base) : (base, "")
    return joinpath(dir, "$(stem).mask$(ext)")
end

struct CubeTransactionInterruption <: Exception
    point::Symbol
end

Base.showerror(io::IO, err::CubeTransactionInterruption) =
    print(io, "simulated transaction interruption at ", err.point)

_cc_txn_present(path::AbstractString) = ispath(path) || islink(path)
_cc_txn_marker(gate::AbstractString) = string(gate, ".stmfit-txn")
_cc_txn_signature(destinations) = bytes2hex(sha256(join(abspath.(destinations), '\0')))
_cc_txn_stage(path, id) = string(path, ".", id, ".stmfit-txn.new", splitext(path)[2])
_cc_txn_backup(path, id) = string(path, ".", id, ".stmfit-txn.old")
_cc_txn_absent(path, id) = string(path, ".", id, ".stmfit-txn.absent")

function _cc_txn_rename(source::AbstractString, destination::AbstractString)
    rc = ccall(:rename, Cint, (Cstring, Cstring), source, destination)
    rc == 0 || Base.systemerror("rename $source to $destination", true)
    return nothing
end

function _cc_txn_fsync_file(path::AbstractString)
    open(path, "r") do io
        rc = ccall(:fsync, Cint, (Cint,), fd(io))
        rc == 0 || Base.systemerror("fsync $path", true)
    end
    return nothing
end

function _cc_txn_fsync_parent(path::AbstractString)
    directory = dirname(abspath(path))
    directory_fd = ccall(:open, Cint, (Cstring, Cint), directory, 0)
    directory_fd >= 0 || Base.systemerror("open directory $directory", true)
    try
        rc = ccall(:fsync, Cint, (Cint,), directory_fd)
        rc == 0 || Base.systemerror("fsync directory $directory", true)
    finally
        ccall(:close, Cint, (Cint,), directory_fd)
    end
    return nothing
end


function _cc_txn_write_marker(marker, id, phase, signature)
    temp_path, io = mktemp(dirname(marker); cleanup=false)
    try
        println(io, "stmfit-output-transaction-v1")
        println(io, "id\t", id)
        println(io, "phase\t", phase)
        println(io, "signature\t", signature)
        flush(io)
        rc = ccall(:fsync, Cint, (Cint,), fd(io))
        rc == 0 || Base.systemerror("fsync transaction marker", true)
        close(io)
        _cc_txn_rename(temp_path, marker)
        _cc_txn_fsync_parent(marker)
    finally
        isopen(io) && close(io)
        rm(temp_path; force=true)
    end
    return nothing
end

function _cc_txn_read_marker(marker, signature)
    islink(marker) && error("refusing symlinked transaction marker: $marker")
    lines = readlines(marker)
    length(lines) == 4 && lines[1] == "stmfit-output-transaction-v1" ||
        error("malformed output transaction marker: $marker")
    fields = Dict{String,String}()
    for line in lines[2:end]
        parts = split(line, '\t'; limit=2)
        length(parts) == 2 || error("malformed output transaction marker: $marker")
        fields[parts[1]] = parts[2]
    end
    get(fields, "signature", "") == signature ||
        error("output transaction destination set does not match marker: $marker")
    phase = get(fields, "phase", "")
    phase in ("prepared", "committed") ||
        error("invalid output transaction phase in $marker")
    id = get(fields, "id", "")
    occursin(r"^[0-9a-f]{16}$", id) || error("invalid output transaction id in $marker")
    return id, phase
end

function _cc_txn_cleanup(destinations, gate, id)
    for destination in destinations
        rm(_cc_txn_stage(destination, id); force=true)
        rm(_cc_txn_backup(destination, id); force=true)
        rm(_cc_txn_absent(destination, id); force=true)
    end
    rm(_cc_txn_marker(gate); force=true)
    _cc_txn_fsync_parent(gate)
    return nothing
end

function _cc_txn_restore_one(destination, id)
    backup = _cc_txn_backup(destination, id)
    absent = _cc_txn_absent(destination, id)
    if _cc_txn_present(backup)
        _cc_txn_present(destination) && rm(destination; force=true)
        _cc_txn_rename(backup, destination)
    elseif isfile(absent)
        _cc_txn_present(destination) && rm(destination; force=true)
        rm(absent; force=true)
    end
    return nothing
end

function _cc_txn_rollback(destinations, gate, id)
    gate_processed = _cc_txn_present(_cc_txn_backup(gate, id)) ||
        isfile(_cc_txn_absent(gate, id))
    gate_processed && _cc_txn_present(gate) && rm(gate; force=true)
    for destination in destinations
        destination == gate || _cc_txn_restore_one(destination, id)
    end
    _cc_txn_restore_one(gate, id)
    _cc_txn_fsync_parent(gate)
    _cc_txn_cleanup(destinations, gate, id)
    return nothing
end

function _cc_txn_recover(destinations, gate)
    marker = _cc_txn_marker(gate)
    _cc_txn_present(marker) || return nothing
    id, phase = _cc_txn_read_marker(marker, _cc_txn_signature(destinations))
    phase == "committed" ? _cc_txn_cleanup(destinations, gate, id) :
        _cc_txn_rollback(destinations, gate, id)
    return nothing
end

function _cc_txn_inject(failpoint, interruptpoint, point)
    interruptpoint == point && throw(CubeTransactionInterruption(point))
    failpoint == point && error("injected transaction failure at $point")
    return nothing
end

function _with_output_transaction(writer::Function, destinations;
        gate=last(destinations), failpoint=nothing, interruptpoint=nothing)
    ordered = unique(abspath.(String.(destinations)))
    isempty(ordered) && throw(ArgumentError("output transaction requires destinations"))
    gate = abspath(String(gate))
    gate in ordered || throw(ArgumentError("transaction gate must be a destination"))
    length(unique(dirname.(ordered))) == 1 ||
        throw(ArgumentError("transaction destinations must share one directory"))
    for destination in ordered
        isdir(destination) && !islink(destination) &&
            throw(ArgumentError("transaction destination is a directory: $destination"))
    end
    _cc_txn_recover(ordered, gate)

    id = bytes2hex(sha256(string(time_ns(), ':', getpid(), ':', objectid(writer))))[1:16]
    staged = Dict(destination => _cc_txn_stage(destination, id) for destination in ordered)
    marker = _cc_txn_marker(gate)
    signature = _cc_txn_signature(ordered)
    for destination in ordered
        for side_path in (staged[destination], _cc_txn_backup(destination, id),
                          _cc_txn_absent(destination, id))
            _cc_txn_present(side_path) && error("transaction side path already exists: $side_path")
        end
    end
    try
        writer(staged)
        all(destination -> isfile(staged[destination]) && !islink(staged[destination]), ordered) ||
            error("one or more staged constant-current outputs are missing or invalid")
        foreach(_cc_txn_fsync_file, values(staged))
        _cc_txn_write_marker(marker, id, "prepared", signature)
        _cc_txn_inject(failpoint, interruptpoint, :after_marker)

        backup_order = vcat([gate], filter(!=(gate), ordered))
        for (index, destination) in enumerate(backup_order)
            if _cc_txn_present(destination)
                _cc_txn_rename(destination, _cc_txn_backup(destination, id))
            else
                write(_cc_txn_absent(destination, id), "absent\n")
                _cc_txn_fsync_file(_cc_txn_absent(destination, id))
            end
            _cc_txn_fsync_parent(destination)
            _cc_txn_inject(failpoint, interruptpoint, Symbol("after_backup_$index"))
        end

        install_order = vcat(filter(!=(gate), ordered), [gate])
        for (index, destination) in enumerate(install_order)
            _cc_txn_rename(staged[destination], destination)
            _cc_txn_fsync_parent(destination)
            point = destination == gate ? :after_gate : Symbol("after_install_$index")
            _cc_txn_inject(failpoint, interruptpoint, point)
        end
        _cc_txn_write_marker(marker, id, "committed", signature)
        _cc_txn_inject(failpoint, interruptpoint, :after_commit)
        _cc_txn_cleanup(ordered, gate, id)
    catch err
        err isa CubeTransactionInterruption && rethrow()
        if _cc_txn_present(marker)
            try
                marker_id, phase = _cc_txn_read_marker(marker, signature)
                marker_id == id || error("output transaction marker changed during commit")
                phase == "committed" ? _cc_txn_cleanup(ordered, gate, id) :
                    _cc_txn_rollback(ordered, gate, id)
            catch rollback_error
                error("output transaction failed and recovery failed: $(sprint(showerror, err)); " *
                      "recovery: $(sprint(showerror, rollback_error))")
            end
        else
            foreach(path -> rm(path; force=true), values(staged))
        end
        rethrow()
    end
    return nothing
end

function _auto_cc_resolution(cube::CubeGrid)
    # Native z-step projected onto the cube z-axis (constant-current searches
    # along the frame normal, which for surface cubes is the z-axis).
    return abs(cube.axes[3, 3])
end

function _sample_column_along_normal(cube::CubeGrid, origin::Vector{Float64},
                                     that::Vector{Float64}, uhat::Vector{Float64},
                                     nhat::Vector{Float64}, t::Float64, u::Float64,
                                     s_values::Vector{Float64},
                                     periodic_axes::NTuple{3,Bool})
    # Sample the cube along the frame normal, then truncate to the contiguous
    # finite segment.  A line intersecting a convex cube region always produces
    # a single contiguous finite run; any internal NaN (hole) marks the column
    # as unusable and it is rejected as nonfinite downstream.
    v_all = Vector{Float64}(undef, length(s_values))
    for k in eachindex(s_values)
        s = s_values[k]
        r = origin .+ t .* that .+ u .* uhat .+ s .* nhat
        v_all[k] = any(periodic_axes) ?
            _sample_cube_periodic(cube, r, periodic_axes) : _sample_cube(cube, r)
    end
    first_finite = findfirst(isfinite, v_all)
    first_finite === nothing && return Float64[], Float64[]
    last_finite = findlast(isfinite, v_all)
    for k in first_finite:last_finite
        isfinite(v_all[k]) || return Float64[], Float64[]
    end
    return s_values[first_finite:last_finite], v_all[first_finite:last_finite]
end

function main()
    opt = _parse_cli(ARGS)
    coords = collect(-opt.half_nm:opt.step_nm:opt.half_nm)
    any(opt.periodic_axes) && println(stderr,
        "WARNING: periodic cube sampling is diagnostic-only; do not use this output as a production mold")
    _ensure_parent(opt.out_tsv)
    if opt.observable == "constant-height"
        _run_constant_height(opt, coords)
    else
        _run_constant_current(opt, coords)
    end
end

function _run_constant_height(opt::Options, coords::Vector{Float64})
    open(opt.out_tsv, "w") do io
        println(io, join(["type", "t_nm", "u_nm", "value"], '\t'))
        for spec in sort(opt.cubes, by=c -> c.typ)
            cube = _read_cube(spec.path, opt.cube_units)
            frame = spec.frame
            origin = frame.origin === nothing ? error("missing origin for cube type=$(spec.typ)") : frame.origin
            t_axis = frame.t_axis === nothing ? error("missing t-axis for cube type=$(spec.typ)") : frame.t_axis
            u_axis = frame.u_axis === nothing ? error("missing u-axis for cube type=$(spec.typ)") : frame.u_axis
            height_nm = frame.height_nm === nothing ? error("missing height for cube type=$(spec.typ)") : frame.height_nm
            that, uhat, nhat = _frame(t_axis, u_axis)
            for u in coords, t in coords
                r = origin .+ t .* that .+ u .* uhat .+ height_nm .* nhat
                val = any(opt.periodic_axes) ?
                    _sample_cube_periodic(cube, r, opt.periodic_axes) : _sample_cube(cube, r)
                println(io, join([spec.typ, @sprintf("%.8g", t), @sprintf("%.8g", u), _fmt(val)], '\t'))
            end
        end
    end
    println("Converted cube maps")
    println("  cubes:     ", join(["$(c.typ):$(c.path)" for c in opt.cubes], ", "))
    println("  out:       ", opt.out_tsv)
    println("  grid:      ", length(coords), "x", length(coords))
    println("  heights:   ", join(["$(c.typ):$(c.frame.height_nm)" for c in opt.cubes], ", "), " nm")
    println("  units:     ", opt.cube_units)
    periodic_label = join([axis for (axis, enabled) in zip(("x", "y", "z"), opt.periodic_axes) if enabled], "")
    println("  periodic:  ", isempty(periodic_label) ? "none" : periodic_label)
    println("  observable: ", opt.observable)
end

function _run_constant_current(opt::Options, coords::Vector{Float64})
    out_path = abspath(opt.out_tsv)
    mask_path = _mask_path(out_path)
    _ensure_parent(mask_path)
    _with_output_transaction([mask_path, out_path]; gate=out_path) do staged
        open(staged[out_path], "w") do io
            open(staged[mask_path], "w") do io_mask
            println(io, join(["type", "t_nm", "u_nm", "value"], '\t'))
            println(io_mask, join(["type", "t_nm", "u_nm", "status"], '\t'))
            for spec in sort(opt.cubes, by=c -> c.typ)
                cube = _read_cube(spec.path, opt.cube_units)
                frame = spec.frame
                origin = frame.origin === nothing ? error("missing origin for cube type=$(spec.typ)") : frame.origin
                t_axis = frame.t_axis === nothing ? error("missing t-axis for cube type=$(spec.typ)") : frame.t_axis
                u_axis = frame.u_axis === nothing ? error("missing u-axis for cube type=$(spec.typ)") : frame.u_axis
                that, uhat, nhat = _frame(t_axis, u_axis)
                s_step = opt.cc_height_resolution_nm === nothing ?
                    _auto_cc_resolution(cube) : opt.cc_height_resolution_nm
                s_values = collect(opt.cc_height_min_nm:s_step:opt.cc_height_max_nm)
                length(s_values) >= 2 ||
                    error("constant-current search range produces fewer than 2 samples; check --cc-height-* options")
                if opt.isovalue !== nothing
                    iso = opt.isovalue
                    iso_source = "explicit"
                else
                    frame_ref = FrameRef(origin, t_axis, u_axis)
                    iso = isovalue_for_mean_height(cube, frame_ref, opt.mean_height_nm;
                        scan_intervals=opt.isovalue_scan_intervals)
                    iso_source = "calibrated"
                end
                counts = Dict{Symbol,Int}(:found => 0, :absent => 0, :ambiguous => 0,
                                          :nonfinite => 0, :insufficient_z => 0)
                for u in coords, t in coords
                    z_col, v_col = _sample_column_along_normal(cube, origin, that, uhat, nhat,
                                                               t, u, s_values, opt.periodic_axes)
                    result = isempty(z_col) ?
                        CrossingResult(NaN, 0, :nonfinite, 0) :
                        first_vacuum_crossing(z_col, v_col, iso)
                    counts[result.status] = get(counts, result.status, 0) + 1
                    println(io, join([spec.typ, @sprintf("%.8g", t), @sprintf("%.8g", u),
                                      _fmt(result.z)], '\t'))
                    println(io_mask, join([spec.typ, @sprintf("%.8g", t), @sprintf("%.8g", u),
                                           result.status], '\t'))
                end
                println("  type $(spec.typ) isovalue:  $iso  ($iso_source)")
                println("  type $(spec.typ) diagnostics:")
                for status in (:found, :absent, :ambiguous, :nonfinite, :insufficient_z)
                    println("    $(String(status)):  ", get(counts, status, 0))
                end
            end
            end
        end
    end
    println("Converted cube maps (constant-current)")
    println("  cubes:     ", join(["$(c.typ):$(c.path)" for c in opt.cubes], ", "))
    println("  out:       ", opt.out_tsv)
    println("  mask:      ", mask_path)
    println("  grid:      ", length(coords), "x", length(coords))
    println("  mean-height policy:  ", opt.mean_height_nm, " nm")
    println("  isovalue scan intervals:  ", opt.isovalue_scan_intervals)
    println("  units:     ", opt.cube_units)
    periodic_label = join([axis for (axis, enabled) in zip(("x", "y", "z"), opt.periodic_axes) if enabled], "")
    println("  periodic:  ", isempty(periodic_label) ? "none" : periodic_label)
    println("  observable: ", opt.observable)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
