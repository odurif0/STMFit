#!/usr/bin/env julia

using Printf

const SUPPORT_FIELDS = [:name, :n_points, :t_min_nm, :t_max_nm, :u_min_nm,
    :u_max_nm, :peak_t_nm, :peak_u_nm, :peak_abs_value, :peak_on_edge,
    :edge_l2_fraction, :edge_absmax_fraction]

function read_typed_maps(path::String)
    lines = readlines(path)
    isempty(lines) && error("empty map TSV: $path")
    header = split(first(lines), '\t')
    indices = Dict(name => findfirst(==(name), header) for name in ("type", "t_nm", "u_nm", "value"))
    all(index -> index !== nothing, values(indices)) ||
        error("map TSV requires type, t_nm, u_nm, value columns")
    maps = Dict{Int,Dict{Tuple{Float64,Float64},Float64}}()
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        typ = parse(Int, fields[something(indices["type"])])
        typ in (0, 1) || error("invalid type at $path:$line_number")
        coordinate = (parse(Float64, fields[something(indices["t_nm"])]),
            parse(Float64, fields[something(indices["u_nm"])]))
        value = parse(Float64, fields[something(indices["value"])])
        isfinite(value) || error("non-finite map value at $path:$line_number")
        typed = get!(maps, typ, Dict{Tuple{Float64,Float64},Float64}())
        haskey(typed, coordinate) && error("duplicate coordinate for type $typ at $path:$line_number")
        typed[coordinate] = value
    end
    return maps
end

function support_metrics(name::String, map_values::Dict{Tuple{Float64,Float64},Float64})
    isempty(map_values) && error("empty map for $name")
    ts = first.(keys(map_values))
    us = last.(keys(map_values))
    t_min, t_max = extrema(ts)
    u_min, u_max = extrema(us)
    is_edge(coordinate) = coordinate[1] in (t_min, t_max) || coordinate[2] in (u_min, u_max)
    peak = argmax(coordinate -> abs(map_values[coordinate]), collect(keys(map_values)))
    peak_abs = abs(map_values[peak])
    total_l2 = sum(abs2, values(map_values))
    edge_l2 = sum(abs2(value) for (coordinate, value) in map_values if is_edge(coordinate))
    edge_absmax = maximum(abs(value) for (coordinate, value) in map_values if is_edge(coordinate))
    return (; name, n_points=length(map_values), t_min_nm=t_min, t_max_nm=t_max,
        u_min_nm=u_min, u_max_nm=u_max, peak_t_nm=peak[1], peak_u_nm=peak[2],
        peak_abs_value=peak_abs, peak_on_edge=is_edge(peak),
        edge_l2_fraction=total_l2 > 0 ? edge_l2 / total_l2 : NaN,
        edge_absmax_fraction=peak_abs > 0 ? edge_absmax / peak_abs : NaN)
end

function audit_map_support(maps)
    all(haskey(maps, typ) for typ in (0, 1)) || error("typed maps must contain types 0 and 1")
    keys(maps[0]) == keys(maps[1]) || error("type maps use different coordinate grids")
    difference = Dict(coordinate => maps[1][coordinate] - maps[0][coordinate]
        for coordinate in keys(maps[0]))
    return [support_metrics("type_0", maps[0]), support_metrics("type_1", maps[1]),
        support_metrics("type_1_minus_type_0", difference)]
end

format_support(value::Bool) = value ? "true" : "false"
format_support(value::Real) = isfinite(Float64(value)) ? @sprintf("%.15g", Float64(value)) : "NaN"
format_support(value) = string(value)

function write_support_audit(path::String, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(SUPPORT_FIELDS, '\t'))
        for row in rows
            println(io, join((format_support(getproperty(row, field)) for field in SUPPORT_FIELDS), '\t'))
        end
    end
end

function parse_options(args)
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        flag = String(args[i])
        flag in ("--maps", "--out") || error("unknown argument: $flag")
        i < length(args) || error("$flag requires a value")
        values[flag] = String(args[i + 1])
        i += 2
    end
    haskey(values, "--maps") && haskey(values, "--out") || error("--maps and --out are required")
    return values["--maps"], values["--out"]
end

function main(args=ARGS)
    maps_path, out_path = parse_options(args)
    rows = audit_map_support(read_typed_maps(maps_path))
    write_support_audit(out_path, rows)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
