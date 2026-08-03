#!/usr/bin/env julia

# Diagnostic-only: measure edge-energy convergence across nested DFT mold frames.
# Requires all frames to share the same central crop (verified by max abs diff).
# Does not enter fitting, N_selected, or production unit assignment.

using Printf
using SHA

function _read_map_values(path::String)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty map TSV: $path"))
    header = split(first(lines), '\t')
    columns = Dict(name => findfirst(==(name), header)
        for name in ("type", "t_nm", "u_nm", "value"))
    all(!isnothing, Base.values(columns)) ||
        throw(ArgumentError("map TSV requires type, t_nm, u_nm, value columns"))
    map_values = Dict{Tuple{Int,Float64,Float64},Float64}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        typ = parse(Int, fields[something(columns["type"])])
        t = parse(Float64, fields[something(columns["t_nm"])])
        u = parse(Float64, fields[something(columns["u_nm"])])
        raw = strip(fields[something(columns["value"])])
        value = something(tryparse(Float64, raw), NaN)
        map_values[(typ, t, u)] = value
    end
    return map_values
end

function _map_info(map_values)
    ts = sort(unique(key[2] for key in keys(map_values)))
    us = sort(unique(key[3] for key in keys(map_values)))
    finite_count = count(isfinite, Base.values(map_values))
    return (; grid_n=length(ts),
        half_nm=max(maximum(abs, ts), maximum(abs, us)),
        step_nm=isempty(ts) ? NaN : minimum(diff(ts)),
        n_points=length(map_values),
        finite_fraction=finite_count / length(map_values),
        finite_count)
end

function _edge_fraction(map_values, half_nm, step_nm)
    coords = collect(-half_nm:step_nm:half_nm)
    is_edge(t, u) = abs(t) == half_nm || abs(u) == half_nm
    total_l2 = 0.0
    edge_l2 = 0.0
    edge_absmax = 0.0
    peak_value = 0.0
    peak_t, peak_u = NaN, NaN
    peak_on_edge = false
    for ((typ, t, u), value) in map_values
        isfinite(value) || continue
        total_l2 += value^2
        abs_val = abs(value)
        if abs_val > peak_value
            peak_value = abs_val
            peak_t, peak_u = t, u
            peak_on_edge = is_edge(t, u)
        end
        if is_edge(t, u)
            edge_l2 += value^2
            edge_absmax = max(edge_absmax, abs_val)
        end
    end
    return (; edge_l2_fraction=total_l2 > 0 ? edge_l2 / total_l2 : NaN,
        edge_absmax_fraction=peak_value > 0 ? edge_absmax / peak_value : NaN,
        peak_t_nm=peak_t, peak_u_nm=peak_u, peak_on_edge)
end

function _central_difference(reference, candidate)
    differences = Float64[]
    for (key, reference_value) in reference
        haskey(candidate, key) || return Inf
        candidate_value = candidate[key]
        (isfinite(reference_value) && isfinite(candidate_value)) || return Inf
        push!(differences, abs(reference_value - candidate_value))
    end
    return isempty(differences) ? Inf : maximum(differences)
end

const FIELDS = (:grid_n, :half_nm, :edge_l2_fraction_type0,
    :edge_l2_fraction_type1, :edge_l2_fraction_diff,
    :peak_on_edge_type0, :peak_on_edge_type1, :central_max_abs_difference,
    :finite_fraction)

function audit_frame_convergence(reference_path::String, candidate_paths::Vector{String})
    reference = _read_map_values(reference_path)
    ref_info = _map_info(reference)
    rows = []
    for path in candidate_paths
        candidate = _read_map_values(path)
        info = _map_info(candidate)
        # Split by type for per-type edge fractions
        ref0 = Dict(k => v for (k, v) in reference if k[1] == 0)
        ref1 = Dict(k => v for (k, v) in reference if k[1] == 1)
        can0 = Dict(k => v for (k, v) in candidate if k[1] == 0)
        can1 = Dict(k => v for (k, v) in candidate if k[1] == 1)
        ref_diff = Dict((0, k[2], k[3]) => ref1[(1, k[2], k[3])] - v
            for (k, v) in ref0 if haskey(ref1, (1, k[2], k[3])))
        can_diff = Dict((0, k[2], k[3]) => can1[(1, k[2], k[3])] - v
            for (k, v) in can0 if haskey(can1, (1, k[2], k[3])))

        edge0 = _edge_fraction(can0, info.half_nm, info.step_nm)
        edge1 = _edge_fraction(can1, info.half_nm, info.step_nm)
        edge_diff = _edge_fraction(can_diff, info.half_nm, info.step_nm)
        central_diff = _central_difference(reference, candidate)

        push!(rows, (; grid_n=info.grid_n, half_nm=info.half_nm,
            edge_l2_fraction_type0=edge0.edge_l2_fraction,
            edge_l2_fraction_type1=edge1.edge_l2_fraction,
            edge_l2_fraction_diff=edge_diff.edge_l2_fraction,
            peak_on_edge_type0=edge0.peak_on_edge,
            peak_on_edge_type1=edge1.peak_on_edge,
            central_max_abs_difference=central_diff,
            finite_fraction=info.finite_fraction))
    end
    return rows
end

function write_convergence_audit(path::String, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(FIELDS), '\t'))
        for row in rows
            values = []
            for field in FIELDS
                value = getproperty(row, field)
                push!(values, value isa Bool ? (value ? "true" : "false") :
                              value isa AbstractFloat ? @sprintf("%.15g", value) :
                              string(value))
            end
            println(io, join(values, '\t'))
        end
    end
end

function _frame_conv_options(args)
    options = Dict{String,String}()
    positional = String[]
    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg in ("--reference", "--out")
            i < length(args) || error("$arg requires a value")
            options[arg] = String(args[i + 1])
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end
    haskey(options, "--reference") || error("--reference is required")
    haskey(options, "--out") || error("--out is required")
    isempty(positional) && error("at least one candidate map path is required")
    return options["--reference"], positional, options["--out"]
end

function main(args=ARGS)
    reference, candidates, out = _frame_conv_options(args)
    rows = audit_frame_convergence(reference, candidates)
    write_convergence_audit(out, rows)
    println("Frame convergence audit")
    println("  reference: ", reference)
    println("  candidates: ", length(candidates))
    println("  output:    ", out)
    for row in rows
        println("  $(row.grid_n)x$(row.grid_n) ±$(row.half_nm) nm: ",
            "edge_L2 t0=$(@sprintf("%.4f", row.edge_l2_fraction_type0)) ",
            "t1=$(@sprintf("%.4f", row.edge_l2_fraction_type1)) ",
            "diff=$(@sprintf("%.4f", row.edge_l2_fraction_diff)) ",
            "cropΔ=$(@sprintf("%.2e", row.central_max_abs_difference))")
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
