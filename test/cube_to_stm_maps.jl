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

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _parse_vec3

const DEFAULT_OUT = "templates/chitosan_stm_maps.tsv"

struct CubeGrid
    origin::Vector{Float64}
    axes::Matrix{Float64}
    n::NTuple{3,Int}
    values::Vector{Float64}
end

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
            and benchmark composition are not used.
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
    return Options(cubes, out_tsv, half_nm, step_nm, cube_units)
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

function main()
    opt = _parse_cli(ARGS)
    coords = collect(-opt.half_nm:opt.step_nm:opt.half_nm)
    _ensure_parent(opt.out_tsv)
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
                val = _sample_cube(cube, r)
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
end

main()
