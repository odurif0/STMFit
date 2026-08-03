#!/usr/bin/env julia

using Printf
using SHA
using TOML

function _frame_sha256(path::String)
    isfile(path) || throw(ArgumentError("artifact not found: $path"))
    return bytes2hex(open(SHA.sha256, path))
end

function read_typed_map(path::String)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty map TSV: $path"))
    header = split(first(lines), '\t')
    columns = Dict(name => findfirst(==(name), header)
        for name in ("type", "t_nm", "u_nm", "value"))
    all(!isnothing, values(columns)) ||
        throw(ArgumentError("map TSV requires type, t_nm, u_nm, value columns"))
    map_values = Dict{Tuple{Int,Float64,Float64},Float64}()
    for (offset, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        line_number = offset + 1
        typ = parse(Int, fields[something(columns["type"])])
        typ in (0, 1) || throw(ArgumentError("invalid type at $path:$line_number"))
        t = parse(Float64, fields[something(columns["t_nm"])])
        u = parse(Float64, fields[something(columns["u_nm"])])
        raw = strip(fields[something(columns["value"])])
        value = something(tryparse(Float64, raw), NaN)
        key = (typ, t, u)
        haskey(map_values, key) && throw(ArgumentError("duplicate map coordinate at $path:$line_number"))
        map_values[key] = value
    end
    isempty(map_values) && throw(ArgumentError("map TSV has no rows: $path"))
    ts = sort(unique(key[2] for key in keys(map_values)))
    us = sort(unique(key[3] for key in keys(map_values)))
    length(ts) == length(us) || throw(ArgumentError("map grid must be square"))
    expected = 2 * length(ts) * length(us)
    length(map_values) == expected || throw(ArgumentError("map grid is incomplete"))
    steps = vcat(diff(ts), diff(us))
    step_nm = isempty(steps) ? NaN : minimum(steps)
    finite_count = count(isfinite, values(map_values))
    return (; path, values=map_values, grid_n=length(ts),
        half_nm=max(maximum(abs, ts), maximum(abs, us)), step_nm,
        finite_count, total_count=length(map_values),
        finite_fraction=finite_count / length(map_values))
end

function _resolve_artifact_path(raw::String, provenance_path::String)
    isempty(raw) && return ""
    isabspath(raw) && return raw
    isfile(raw) && return normpath(raw)
    return normpath(joinpath(dirname(provenance_path), raw))
end

function read_frame_provenance(path::String)
    data = TOML.parsefile(path)
    required = ("maps_sha256", "glcn_cube_sha256", "glcnac_cube_sha256",
        "half_nm", "step_nm", "height_nm", "sample_bias_ev", "cube_units")
    all(key -> haskey(data, key), required) ||
        throw(ArgumentError("incomplete frame provenance: $path"))
    cube_paths = [_resolve_artifact_path(String(get(data, key, "")), path)
        for key in ("glcn_cube_path", "glcnac_cube_path")]
    return (; path, provider=String(get(data, "provider", "unknown")),
        maps_sha256=String(data["maps_sha256"]),
        glcn_cube_sha256=String(data["glcn_cube_sha256"]),
        glcnac_cube_sha256=String(data["glcnac_cube_sha256"]),
        half_nm=Float64(data["half_nm"]), step_nm=Float64(data["step_nm"]),
        height_nm=Float64(data["height_nm"]),
        sample_bias_ev=Float64(data["sample_bias_ev"]),
        cube_units=String(data["cube_units"]), cube_paths,
        cubes_present=length(cube_paths) == 2 && all(isfile, cube_paths))
end

function _artifact_summary(map, provenance)
    hash_ok = provenance !== nothing && _frame_sha256(map.path) == provenance.maps_sha256
    return (; grid_n=map.grid_n, half_nm=map.half_nm, step_nm=map.step_nm,
        finite_fraction=map.finite_fraction, maps_sha256=_frame_sha256(map.path),
        provenance_present=provenance !== nothing, provenance_hash_ok=hash_ok,
        cubes_present=provenance !== nothing && provenance.cubes_present)
end

function _central_difference(reference, candidate)
    differences = Float64[]
    for (key, reference_value) in reference.values
        haskey(candidate.values, key) || return Inf
        candidate_value = candidate.values[key]
        (isfinite(reference_value) && isfinite(candidate_value)) || return Inf
        push!(differences, abs(reference_value - candidate_value))
    end
    return isempty(differences) ? Inf : maximum(differences)
end

function audit_frame_compatibility(reference_maps::String, reference_provenance::Union{Nothing,String},
                                   candidate_maps::String, candidate_provenance::Union{Nothing,String})
    reference_map = read_typed_map(reference_maps)
    candidate_map = read_typed_map(candidate_maps)
    reference_prov = reference_provenance === nothing ? nothing : read_frame_provenance(reference_provenance)
    candidate_prov = candidate_provenance === nothing ? nothing : read_frame_provenance(candidate_provenance)
    reference = _artifact_summary(reference_map, reference_prov)
    candidate = _artifact_summary(candidate_map, candidate_prov)
    central_difference = _central_difference(reference_map, candidate_map)

    reason = if reference_prov === nothing
        "reference_provenance_missing"
    elseif candidate_prov === nothing
        "candidate_provenance_missing"
    elseif !reference.provenance_hash_ok || !candidate.provenance_hash_ok
        "maps_hash_mismatch"
    elseif reference_prov.glcn_cube_sha256 != candidate_prov.glcn_cube_sha256 ||
           reference_prov.glcnac_cube_sha256 != candidate_prov.glcnac_cube_sha256
        "cube_hash_mismatch"
    elseif reference_prov.cube_units != candidate_prov.cube_units ||
           !isapprox(reference_prov.height_nm, candidate_prov.height_nm; atol=1e-12) ||
           !isapprox(reference_prov.sample_bias_ev, candidate_prov.sample_bias_ev; atol=1e-12) ||
           !isapprox(reference_prov.step_nm, candidate_prov.step_nm; atol=1e-12)
        "extraction_conditions_mismatch"
    elseif candidate.half_nm <= reference.half_nm
        "candidate_not_wider"
    elseif reference.finite_fraction < 1.0 || candidate.finite_fraction < 1.0
        "nonfinite_map_values"
    elseif !isfinite(central_difference) || central_difference > 1e-12
        "central_crop_mismatch"
    else
        "compatible_nested_frame"
    end
    return (; compatible=reason == "compatible_nested_frame", reason,
        central_max_abs_difference=central_difference, reference, candidate)
end

function write_frame_compatibility(path::String, audit)
    mkpath(dirname(path))
    fields = ("compatible", "reason", "central_max_abs_difference",
        "reference_grid_n", "reference_half_nm", "reference_finite_fraction",
        "reference_provenance_present", "reference_maps_hash_ok", "reference_cubes_present",
        "candidate_grid_n", "candidate_half_nm", "candidate_finite_fraction",
        "candidate_provenance_present", "candidate_maps_hash_ok", "candidate_cubes_present")
    values = (audit.compatible, audit.reason, audit.central_max_abs_difference,
        audit.reference.grid_n, audit.reference.half_nm, audit.reference.finite_fraction,
        audit.reference.provenance_present, audit.reference.provenance_hash_ok,
        audit.reference.cubes_present, audit.candidate.grid_n, audit.candidate.half_nm,
        audit.candidate.finite_fraction, audit.candidate.provenance_present,
        audit.candidate.provenance_hash_ok, audit.candidate.cubes_present)
    open(path, "w") do io
        println(io, join(fields, '\t'))
        println(io, join(values, '\t'))
    end
end

function _frame_options(args)
    options = Dict{String,String}()
    i = 1
    while i <= length(args)
        flag = String(args[i])
        flag in ("--reference-maps", "--reference-provenance", "--candidate-maps",
            "--candidate-provenance", "--out") || error("unknown argument: $flag")
        i < length(args) || error("$flag requires a value")
        options[flag] = String(args[i + 1])
        i += 2
    end
    for required in ("--reference-maps", "--reference-provenance", "--candidate-maps", "--out")
        haskey(options, required) || error("$required is required")
    end
    return options
end

function main(args=ARGS)
    options = _frame_options(args)
    candidate_provenance = get(options, "--candidate-provenance", nothing)
    audit = audit_frame_compatibility(options["--reference-maps"],
        options["--reference-provenance"], options["--candidate-maps"], candidate_provenance)
    write_frame_compatibility(options["--out"], audit)
    return audit
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
