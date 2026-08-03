#!/usr/bin/env julia

# Diagnostic-only real-data controls layered on the frozen common/contrast
# pipeline. Registration runs exactly once; every perturbation and negative
# control reuses its frozen common image, mask, state, and identity.

if !isdefined(Main, :WholeRoiDiagnosticResult)
    include(joinpath(@__DIR__, "whole_roi_end_to_end.jl"))
end

using SHA
using TOML

struct FrozenControlScore
    label::String
    scoring::FrozenContrastResult
    registration_tsv::String
end

struct FrozenRealDiagnosticBundle
    nominal::WholeRoiDiagnosticResult
    perturbations::Vector{FrozenControlScore}
    controls::Vector{FrozenControlScore}
    registration_tsv::String
    perturbation_stable::Bool
    negative_control_beaten::Bool
    final_abstain::Bool
    final_abstention_reasons::Vector{String}
end

struct RealObservableContract
    mode::String
    provider::String
    provenance_path::String
    provenance_sha256::String
    maps_path::String
    maps_sha256::String
    invalid_mask_path::String
    invalid_mask_sha256::String
    molds_path::String
    molds_sha256::String
    mold_binding_path::String
    mold_binding_sha256::String
    valid_mask_pixels::Int
end

const REAL_CC_MAP_ENV = "STMFIT_CONSTANT_CURRENT_MAPS"
const REAL_CC_INVALID_MASK_ENV = "STMFIT_CONSTANT_CURRENT_INVALID_MASK"
const REAL_CC_MOLD_PROVENANCE_ENV = "STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE"
const REAL_CC_EXPECTED_PARITY_FLIP = "t"
const REAL_CC_EXPECTED_MIRROR_FLIP = "u"
const REAL_CC_EXPECTED_NORMALIZE = "zscore"

const _REAL_FORBIDDEN_FLAGS = Set([
    "--truth", "--grade", "--expected-n", "--control-sequence",
    "--manifest", "--benchmark-manifest", "--full145", "--control",
])

const _REAL_FORBIDDEN_ENV_TOKENS = (
    "TRUTH", "GRADE", "EXPECTED_N", "EXPECTEDN", "CONTROL_SEQUENCE",
    "CONTROLSEQUENCE", "MANIFEST", "FULL145",
)

function _sha256_file(path::AbstractString)
    isfile(path) || throw(ArgumentError("cannot hash missing file: $path"))
    return open(path, "r") do io
        bytes2hex(sha256(io))
    end
end

function _env_value(env, key::AbstractString)
    return haskey(env, key) ? String(env[key]) : ""
end

function _nonempty(value::AbstractString)
    stripped = strip(String(value))
    return isempty(stripped) ? nothing : stripped
end

function reject_forbidden_real_runtime_inputs(args;
        env=ENV, boundary::AbstractString="real diagnostic")
    for arg in args
        token = lowercase(replace(first(split(String(arg), '='; limit=2)), "_" => "-"))
        if token in _REAL_FORBIDDEN_FLAGS
            throw(ArgumentError(
                "$boundary rejects forbidden benchmark/truth input flag: $arg"))
        end
    end
    for (key, value) in env
        isempty(strip(String(value))) && continue
        upper = uppercase(String(key))
        startswith(upper, "STMFIT") || continue
        if any(token -> occursin(token, upper), _REAL_FORBIDDEN_ENV_TOKENS)
            throw(ArgumentError(
                "$boundary rejects forbidden benchmark/truth input environment variable: $key"))
        end
    end
    return nothing
end

function read_real_file_list(path::AbstractString)
    isfile(path) || throw(ArgumentError("real diagnostic file list not found: $path"))
    files = String[]
    for (line_no, raw) in enumerate(readlines(path))
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        basename(line) == line || throw(ArgumentError(
            "real diagnostic file list row $line_no must be an SXM basename, got: $line"))
        endswith(lowercase(line), ".sxm") || throw(ArgumentError(
            "real diagnostic file list row $line_no must name a .sxm file, got: $line"))
        push!(files, line)
    end
    isempty(files) && throw(ArgumentError("real diagnostic file list is empty: $path"))
    return files
end

function validate_real_diagnostic_paths(summary_path::AbstractString,
                                        mold_path::AbstractString,
                                        config_path::AbstractString,
                                        data_dir::AbstractString)
    isfile(summary_path) || throw(ArgumentError(
        "selection summary not found: $summary_path"))
    isfile(mold_path) || throw(ArgumentError("converged mold TSV not found: $mold_path"))
    isfile(config_path) || throw(ArgumentError("diagnostic config not found: $config_path"))
    isdir(data_dir) || throw(ArgumentError("STM data directory not found: $data_dir"))
    return nothing
end

function _required_key(table, key::AbstractString, label::AbstractString)
    haskey(table, key) || throw(ArgumentError("$label missing required key '$key'"))
    return table[key]
end

function _required_string(table, key::AbstractString, label::AbstractString)
    value = _required_key(table, key, label)
    value isa AbstractString ||
        throw(ArgumentError("$label.$key must be a string, got $(typeof(value))"))
    isempty(strip(value)) && throw(ArgumentError("$label.$key must not be empty"))
    return String(value)
end

function _required_float(table, key::AbstractString, label::AbstractString)
    value = _required_key(table, key, label)
    (value isa Real && !(value isa Bool)) ||
        throw(ArgumentError("$label.$key must be a non-Bool number, got $(typeof(value))"))
    out = Float64(value)
    isfinite(out) || throw(ArgumentError("$label.$key must be finite"))
    return out
end

function _required_hash(table, key::AbstractString, label::AbstractString)
    value = _required_string(table, key, label)
    occursin(r"^[0-9a-f]{64}$", value) ||
        throw(ArgumentError("$label.$key must be a lowercase SHA-256"))
    return value
end

function _same_float(actual::Real, expected::Real, label::AbstractString; atol=1e-9)
    isapprox(Float64(actual), Float64(expected); rtol=0.0, atol=atol) ||
        throw(ArgumentError("$label mismatch: provenance=$(Float64(actual)) config=$(Float64(expected))"))
    return nothing
end

function _resolve_artifact_path(raw::AbstractString, provenance_dir::AbstractString)
    path = String(raw)
    if isabspath(path)
        isfile(path) && return path
        staged = joinpath(provenance_dir, basename(path))
        isfile(staged) && return abspath(staged)
        return path
    end
    isfile(path) && return abspath(path)
    return abspath(joinpath(provenance_dir, path))
end

function _mask_valid_pixel_count(mask_path::AbstractString)
    lines = readlines(mask_path)
    isempty(lines) && throw(ArgumentError("constant-current validity mask is empty: $mask_path"))
    header = split(first(lines), '\t')
    required = ("type", "t_nm", "u_nm", "status")
    for col in required
        col in header || throw(ArgumentError(
            "constant-current validity mask missing column '$col': $mask_path"))
    end
    status_col = findfirst(==("status"), header)
    count_found = 0
    for (line_no, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        length(fields) >= status_col || throw(ArgumentError(
            "constant-current validity mask row $(line_no + 1) is truncated: $mask_path"))
        status = strip(fields[status_col])
        status == "found" || throw(ArgumentError(
            "constant-current validity mask has invalid isosurface at row $(line_no + 1): status=$status"))
        count_found += 1
    end
    count_found > 0 || throw(ArgumentError(
        "constant-current validity mask contains no valid pixels: $mask_path"))
    return count_found
end

function _require_same_grid_keys(map_path::AbstractString, mask_path::AbstractString)
    function keys_for(path, expected_value_col)
        lines = readlines(path)
        isempty(lines) && throw(ArgumentError("empty constant-current artifact: $path"))
        header = split(first(lines), '\t')
        for col in ("type", "t_nm", "u_nm", expected_value_col)
            col in header || throw(ArgumentError(
                "constant-current artifact missing column '$col': $path"))
        end
        type_col = findfirst(==("type"), header)
        t_col = findfirst(==("t_nm"), header)
        u_col = findfirst(==("u_nm"), header)
        return [(fields[type_col], fields[t_col], fields[u_col])
                for fields in (split(line, '\t') for line in lines[2:end])
                if !isempty(strip(join(fields, '\t')))]
    end
    keys_for(map_path, "value") == keys_for(mask_path, "status") ||
        throw(ArgumentError("constant-current map and validity-mask grid keys differ"))
    return nothing
end

function _validate_nominal_map_grid(map_path::AbstractString, sidecar, label)
    half_nm = _required_float(sidecar, "half_nm", label)
    step_nm = _required_float(sidecar, "step_nm", label)
    half_nm > 0 || throw(ArgumentError("$label.half_nm must be positive"))
    step_nm > 0 || throw(ArgumentError("$label.step_nm must be positive"))
    intervals_f = 2half_nm / step_nm
    intervals = round(Int, intervals_f)
    isapprox(intervals_f, intervals; rtol=0.0, atol=1e-8) || throw(ArgumentError(
        "constant-current map grid half_nm/step_nm do not form an integer grid"))
    intervals >= 1 || throw(ArgumentError(
        "constant-current map grid must contain at least two intervals"))

    lines = readlines(map_path)
    isempty(lines) && throw(ArgumentError("empty constant-current map: $map_path"))
    header = split(first(lines), '\t')
    for col in ("type", "t_nm", "u_nm", "value")
        col in header || throw(ArgumentError(
            "constant-current map missing column '$col': $map_path"))
    end
    type_col = findfirst(==("type"), header)
    t_col = findfirst(==("t_nm"), header)
    u_col = findfirst(==("u_nm"), header)
    seen = Dict(0 => Set{Tuple{Int,Int}}(), 1 => Set{Tuple{Int,Int}}())
    for (line_no, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        length(fields) >= maximum((type_col, t_col, u_col)) || throw(ArgumentError(
            "constant-current map row $(line_no + 1) is truncated: $map_path"))
        typ = tryparse(Int, strip(fields[type_col]))
        typ in (0, 1) || throw(ArgumentError(
            "constant-current map row $(line_no + 1) has invalid type: $(fields[type_col])"))
        t = tryparse(Float64, strip(fields[t_col]))
        u = tryparse(Float64, strip(fields[u_col]))
        (t !== nothing && u !== nothing && isfinite(t) && isfinite(u)) ||
            throw(ArgumentError("constant-current map row $(line_no + 1) has invalid grid coordinates"))
        ti = round(Int, (t + half_nm) / step_nm)
        ui = round(Int, (u + half_nm) / step_nm)
        0 <= ti <= intervals && 0 <= ui <= intervals || throw(ArgumentError(
            "constant-current map row $(line_no + 1) lies outside provenance grid"))
        isapprox(t, -half_nm + ti * step_nm; rtol=0.0, atol=1e-7) ||
            throw(ArgumentError("constant-current map t grid step mismatch at row $(line_no + 1)"))
        isapprox(u, -half_nm + ui * step_nm; rtol=0.0, atol=1e-7) ||
            throw(ArgumentError("constant-current map u grid step mismatch at row $(line_no + 1)"))
        push!(seen[typ], (ti, ui))
    end
    expected = (intervals + 1)^2
    for typ in (0, 1)
        length(seen[typ]) == expected || throw(ArgumentError(
            "constant-current map type $typ grid has $(length(seen[typ])) unique pixels, expected $expected from provenance half_nm/step_nm"))
    end
    return nothing
end

function _finite_array(table, key::AbstractString, label::AbstractString)
    vals = _required_key(table, key, label)
    vals isa AbstractVector || throw(ArgumentError("$label.$key must be an array"))
    return [_required_float(Dict("value" => value), "value", "$label.$key[$i]")
            for (i, value) in enumerate(vals)]
end

function _validate_bracket_heights(sidecar, observable, label, provenance_dir::AbstractString)
    expected = _finite_array(observable, "bracket_heights_nm", "[observable]")
    actual = _finite_array(sidecar, "bracket_heights_nm", label)
    sort(actual) == sort(expected) || throw(ArgumentError(
        "constant-current bracket_heights_nm mismatch between config and provenance"))
    nominal = _required_float(observable, "nominal_height_nm", "[observable]")
    _same_float(_required_float(sidecar, "nominal_height_nm", label), nominal,
        "constant-current nominal_height_nm"; atol=1e-12)
    non_nominal = sort(filter(h -> !isapprox(h, nominal; atol=1e-12, rtol=0.0),
        unique(expected)))
    artifacts = _required_key(sidecar, "bracket_artifacts", label)
    artifacts isa AbstractVector || throw(ArgumentError(
        "$label.bracket_artifacts must be an array"))
    artifact_heights = Float64[]
    for item in artifacts
        item isa AbstractDict || throw(ArgumentError(
            "$label.bracket_artifacts entries must be tables"))
        h = _required_float(item, "height_nm", "$label.bracket_artifacts")
        push!(artifact_heights, h)
        map_raw = _required_string(item, "map_path", "$label.bracket_artifacts[$h]")
        mask_raw = _required_string(item, "invalid_mask_path", "$label.bracket_artifacts[$h]")
        map_path = _resolve_artifact_path(map_raw, provenance_dir)
        mask_path = _resolve_artifact_path(mask_raw, provenance_dir)
        _sha256_file(map_path) == _required_hash(item, "map_sha256", "$label.bracket_artifacts[$h]") ||
            throw(ArgumentError("constant-current bracket map SHA-256 mismatch at height $h"))
        _sha256_file(mask_path) == _required_hash(item, "invalid_mask_sha256", "$label.bracket_artifacts[$h]") ||
            throw(ArgumentError("constant-current bracket validity-mask SHA-256 mismatch at height $h"))
    end
    sort!(artifact_heights)
    artifact_heights == non_nominal || throw(ArgumentError(
        "constant-current bracket_artifacts heights do not match non-nominal bracket heights"))
    return nothing
end

function _validate_type_frames(sidecar, label)
    frames = _required_key(sidecar, "type_frames", label)
    frames isa AbstractVector || throw(ArgumentError(
        "$label.type_frames must be an array"))
    seen = Set{Int}()
    for item in frames
        item isa AbstractDict || throw(ArgumentError(
            "$label.type_frames entries must be tables"))
        raw_type = _required_key(item, "type", "$label.type_frames")
        (raw_type isa Real && !(raw_type isa Bool)) || throw(ArgumentError(
            "$label.type_frames.type must be non-Bool numeric"))
        typ = Int(raw_type)
        Float64(raw_type) == Float64(typ) || throw(ArgumentError(
            "$label.type_frames.type must be an integer"))
        typ in (0, 1) || throw(ArgumentError(
            "$label.type_frames.type must be 0 or 1"))
        typ in seen && throw(ArgumentError(
            "$label.type_frames has duplicate type $typ"))
        push!(seen, typ)
        for key in ("origin_nm", "t_axis", "u_axis")
            vec = _required_key(item, key, "$label.type_frames[$typ]")
            vec isa AbstractVector || throw(ArgumentError(
                "$label.type_frames[$typ].$key must be a vector"))
            length(vec) == 3 || throw(ArgumentError(
                "$label.type_frames[$typ].$key must have length 3"))
            all(value -> value isa Real && !(value isa Bool) && isfinite(Float64(value)), vec) ||
                throw(ArgumentError(
                    "$label.type_frames[$typ].$key must be finite numbers"))
        end
    end
    seen == Set([0, 1]) || throw(ArgumentError(
        "$label.type_frames must contain exactly one frame for type 0 and type 1"))
    return nothing
end

function _validate_type_isovalues(sidecar, label)
    values = _required_key(sidecar, "type_isovalues", label)
    values isa AbstractVector || throw(ArgumentError(
        "$label.type_isovalues must be an array"))
    seen = Set{Int}()
    for item in values
        item isa AbstractDict || throw(ArgumentError(
            "$label.type_isovalues entries must be tables"))
        raw_type = _required_key(item, "type", "$label.type_isovalues")
        (raw_type isa Real && !(raw_type isa Bool)) || throw(ArgumentError(
            "$label.type_isovalues.type must be non-Bool numeric"))
        typ = Int(raw_type)
        Float64(raw_type) == Float64(typ) || throw(ArgumentError(
            "$label.type_isovalues.type must be an integer"))
        typ in (0, 1) || throw(ArgumentError(
            "$label.type_isovalues.type must be 0 or 1"))
        typ in seen && throw(ArgumentError(
            "$label.type_isovalues has duplicate type $typ"))
        push!(seen, typ)
        _required_float(item, "isovalue", "$label.type_isovalues[$typ]")
        if haskey(item, "source")
            _required_string(item, "source", "$label.type_isovalues[$typ]")
        end
    end
    seen == Set([0, 1]) || throw(ArgumentError(
        "$label.type_isovalues must contain exactly one isovalue for type 0 and type 1"))
    return nothing
end

function _validate_extraction_metadata(sidecar, observable, label)
    nominal = _required_float(observable, "nominal_height_nm", "[observable]")
    _same_float(_required_float(sidecar, "height_nm", label), nominal,
        "constant-current height_nm"; atol=1e-12)
    for key in ("half_nm", "step_nm", "z_spacing_nm")
        value = _required_float(sidecar, key, label)
        value > 0 || throw(ArgumentError("$label.$key must be positive"))
    end
    units = _required_string(sidecar, "cube_units", label)
    units in ("bohr", "angstrom", "a", "nm") ||
        throw(ArgumentError("$label.cube_units is unsupported: $units"))
    scan_intervals = _required_key(sidecar, "isovalue_scan_intervals", label)
    scan_intervals isa Integer && !(scan_intervals isa Bool) || throw(ArgumentError(
        "$label.isovalue_scan_intervals must be a positive integer"))
    scan_intervals > 0 || throw(ArgumentError(
        "$label.isovalue_scan_intervals must be a positive integer"))
    return nothing
end

function _binding_path(raw::Union{Nothing,AbstractString}, molds_path::AbstractString, env)
    raw !== nothing && return String(raw)
    env_raw = _nonempty(_env_value(env, REAL_CC_MOLD_PROVENANCE_ENV))
    env_raw !== nothing && return env_raw
    return string(molds_path, ".provenance.toml")
end

function _validate_mold_binding(sidecar, label, prov::AbstractString,
        actual_prov_sha::AbstractString, map_path::AbstractString,
        map_sha::AbstractString, molds_path::AbstractString,
        mold_binding_path::Union{Nothing,AbstractString}, env,
        half_nm::Union{Nothing,Real}, step_nm::Union{Nothing,Real})
    actual_molds_path = abspath(molds_path)
    isfile(actual_molds_path) || throw(ArgumentError(
        "constant-current mold TSV not found: $actual_molds_path"))
    binding_raw = _binding_path(mold_binding_path, actual_molds_path, env)
    binding = abspath(binding_raw)
    isfile(binding) || throw(ArgumentError(
        "constant-current mold binding sidecar not found: $binding"))
    binding_sha = _sha256_file(binding)
    data = TOML.parsefile(binding)
    bind_label = "constant-current mold binding"
    _required_string(data, "schema", bind_label) == "stmfit-constant-current-mold-binding-v1" ||
        throw(ArgumentError("constant-current mold binding schema mismatch"))
    _required_hash(data, "source_provenance_sha256", bind_label) == actual_prov_sha ||
        throw(ArgumentError("constant-current mold binding source provenance SHA-256 mismatch"))
    _required_hash(data, "source_maps_sha256", bind_label) == map_sha ||
        throw(ArgumentError("constant-current mold binding source map SHA-256 mismatch"))
    _required_hash(data, "molds_sha256", bind_label) == _sha256_file(actual_molds_path) ||
        throw(ArgumentError("constant-current mold TSV SHA-256 mismatch"))
    recorded_molds = _required_string(data, "molds_path", bind_label)
    basename(recorded_molds) == basename(actual_molds_path) || throw(ArgumentError(
        "constant-current mold binding path does not name the actual mold TSV"))
    recorded_map = _required_string(data, "source_maps_path", bind_label)
    basename(recorded_map) == basename(map_path) || throw(ArgumentError(
        "constant-current mold binding source map path mismatch"))
    recorded_prov = _required_string(data, "source_provenance_path", bind_label)
    basename(recorded_prov) == basename(prov) || throw(ArgumentError(
        "constant-current mold binding source provenance path mismatch"))
    if half_nm !== nothing
        _same_float(_required_float(data, "half_nm", bind_label), half_nm,
            "constant-current mold half_nm"; atol=1e-12)
    else
        _required_float(data, "half_nm", bind_label) > 0 || throw(ArgumentError(
            "$bind_label.half_nm must be positive"))
    end
    if step_nm !== nothing
        _same_float(_required_float(data, "step_nm", bind_label), step_nm,
            "constant-current mold step_nm"; atol=1e-12)
    else
        _required_float(data, "step_nm", bind_label) > 0 || throw(ArgumentError(
            "$bind_label.step_nm must be positive"))
    end
    _same_float(_required_float(data, "half_nm", bind_label),
        _required_float(sidecar, "half_nm", label),
        "constant-current mold/source half_nm"; atol=1e-12)
    _same_float(_required_float(data, "step_nm", bind_label),
        _required_float(sidecar, "step_nm", label),
        "constant-current mold/source step_nm"; atol=1e-12)
    _required_string(data, "parity_flip", bind_label) == REAL_CC_EXPECTED_PARITY_FLIP ||
        throw(ArgumentError("constant-current mold binding parity_flip mismatch"))
    _required_string(data, "mirror_flip", bind_label) == REAL_CC_EXPECTED_MIRROR_FLIP ||
        throw(ArgumentError("constant-current mold binding mirror_flip mismatch"))
    _required_string(data, "normalize", bind_label) == REAL_CC_EXPECTED_NORMALIZE ||
        throw(ArgumentError("constant-current mold binding normalize mismatch"))
    return actual_molds_path, _sha256_file(actual_molds_path), binding, binding_sha
end

function validate_real_observable_contract(config_path::AbstractString;
        env=ENV,
        provenance_path::Union{Nothing,AbstractString}=nothing,
        provenance_sha256::Union{Nothing,AbstractString}=nothing,
        maps_path::Union{Nothing,AbstractString}=nothing,
        invalid_mask_path::Union{Nothing,AbstractString}=nothing,
        molds_path::Union{Nothing,AbstractString}=nothing,
        mold_provenance_path::Union{Nothing,AbstractString}=nothing,
        half_nm::Union{Nothing,Real}=nothing,
        step_nm::Union{Nothing,Real}=nothing)
    isfile(config_path) || throw(ArgumentError(
        "diagnostic config not found: $config_path"))
    cfg = TOML.parsefile(config_path)
    haskey(cfg, "observable") || return nothing
    observable = cfg["observable"]
    observable isa AbstractDict ||
        throw(ArgumentError("[observable] must be a TOML table"))

    mode = _required_string(observable, "mode", "[observable]")
    mode == "constant-current" || throw(ArgumentError(
        "unsupported diagnostic observable mode: $mode"))
    provider = _required_string(observable, "provider", "[observable]")
    expected_schema = _required_string(observable, "provenance_schema", "[observable]")
    provenance_env = _required_string(observable, "provenance_path_env", "[observable]")
    provenance_sha_env = _required_string(observable, "provenance_sha256_env", "[observable]")
    map_hash_key = _required_string(observable, "maps_sha256_key", "[observable]")
    mask_hash_key = _required_string(observable, "invalid_mask_sha256_key", "[observable]")
    cube0_key = _required_string(observable, "glcn_cube_sha256_key", "[observable]")
    cube1_key = _required_string(observable, "glcnac_cube_sha256_key", "[observable]")

    prov_raw = provenance_path === nothing ?
        _nonempty(_env_value(env, provenance_env)) : String(provenance_path)
    prov_raw === nothing && throw(ArgumentError(
        "constant-current observable requires provenance path via $provenance_env"))
    prov = abspath(prov_raw)
    isfile(prov) || throw(ArgumentError(
        "constant-current provenance not found: $prov"))
    actual_prov_sha = _sha256_file(prov)
    expected_prov_sha = provenance_sha256 === nothing ?
        _nonempty(_env_value(env, provenance_sha_env)) : String(provenance_sha256)
    expected_prov_sha === nothing && throw(ArgumentError(
        "constant-current observable requires provenance SHA-256 via $provenance_sha_env"))
    actual_prov_sha == expected_prov_sha || throw(ArgumentError(
        "constant-current provenance SHA-256 mismatch"))

    sidecar = TOML.parsefile(prov)
    label = "constant-current provenance"
    _required_string(sidecar, "schema", label) == expected_schema ||
        throw(ArgumentError("constant-current provenance schema mismatch"))
    _required_string(sidecar, "provider", label) == provider ||
        throw(ArgumentError("constant-current provenance provider mismatch"))
    _required_string(sidecar, "observable", label) == mode ||
        throw(ArgumentError("constant-current provenance observable mismatch"))
    _same_float(_required_float(sidecar, "sample_bias_ev", label),
        _required_float(observable, "sample_bias_ev", "[observable]"),
        "constant-current sample_bias_ev"; atol=1e-9)
    _required_string(sidecar, "crossing_policy", label) ==
        _required_string(observable, "crossing_policy", "[observable]") ||
        throw(ArgumentError("constant-current crossing_policy mismatch"))
    provenance_dir = dirname(prov)
    _validate_extraction_metadata(sidecar, observable, label)
    _validate_bracket_heights(sidecar, observable, label, provenance_dir)
    _required_hash(sidecar, cube0_key, label) ==
        _required_hash(observable, "glcn_cube_sha256", "[observable]") ||
        throw(ArgumentError("constant-current GlcN cube hash mismatch"))
    _required_hash(sidecar, cube1_key, label) ==
        _required_hash(observable, "glcnac_cube_sha256", "[observable]") ||
        throw(ArgumentError("constant-current GlcNAc cube hash mismatch"))
    _validate_type_frames(sidecar, label)
    _validate_type_isovalues(sidecar, label)

    map_raw = maps_path === nothing ?
        _nonempty(_env_value(env, REAL_CC_MAP_ENV)) : String(maps_path)
    map_raw === nothing && haskey(sidecar, "maps_path") &&
        (map_raw = String(sidecar["maps_path"]))
    map_raw === nothing && throw(ArgumentError(
        "constant-current observable requires map path via $REAL_CC_MAP_ENV or provenance maps_path"))
    mask_raw = invalid_mask_path === nothing ?
        _nonempty(_env_value(env, REAL_CC_INVALID_MASK_ENV)) : String(invalid_mask_path)
    mask_raw === nothing && haskey(sidecar, "invalid_mask_path") &&
        (mask_raw = String(sidecar["invalid_mask_path"]))
    mask_raw === nothing && throw(ArgumentError(
        "constant-current observable requires validity mask path via $REAL_CC_INVALID_MASK_ENV or provenance invalid_mask_path"))

    map_path = _resolve_artifact_path(map_raw, provenance_dir)
    mask_path = _resolve_artifact_path(mask_raw, provenance_dir)
    map_sha = _sha256_file(map_path)
    mask_sha = _sha256_file(mask_path)
    map_sha == _required_hash(sidecar, map_hash_key, label) ||
        throw(ArgumentError("constant-current map SHA-256 mismatch"))
    mask_sha == _required_hash(sidecar, mask_hash_key, label) ||
        throw(ArgumentError("constant-current validity-mask SHA-256 mismatch"))
    _require_same_grid_keys(map_path, mask_path)
    _validate_nominal_map_grid(map_path, sidecar, label)
    valid_pixels = _mask_valid_pixel_count(mask_path)
    molds_path === nothing && throw(ArgumentError(
        "constant-current observable validation requires the actual mold TSV path"))
    actual_molds_path, molds_sha, binding, binding_sha = _validate_mold_binding(
        sidecar, label, prov, actual_prov_sha, map_path, map_sha,
        String(molds_path), mold_provenance_path, env, half_nm, step_nm)

    return RealObservableContract(mode, provider, prov, actual_prov_sha,
        map_path, map_sha, mask_path, mask_sha, actual_molds_path, molds_sha,
        binding, binding_sha, valid_pixels)
end

function _replace_frozen_input(input::FrozenScorerInput;
        direction=input.direction, phase=input.phase, mirror=input.mirror,
        theta_total=input.theta_total, shift_t_nm=input.shift_t_nm,
        shift_u_nm=input.shift_u_nm, blur_sigma_nm=input.blur_sigma_nm)
    return FrozenScorerInput(direction, phase, mirror, theta_total,
        shift_t_nm, shift_u_nm, blur_sigma_nm,
        copy(input.common_image), input.common_image_hash,
        copy(input.mask), input.mask_hash, input.grid_n, input.common_only_sse)
end

function _bounded_perturb(value::Float64, step::Float64, coarse)
    plus = value + step
    plus <= coarse.max + eps(Float64) && return min(plus, coarse.max)
    minus = value - step
    minus >= coarse.min - eps(Float64) && return max(minus, coarse.min)
    return value
end

function _perturbed_inputs(input::FrozenScorerInput, registration::FrozenRegistration, config)
    st = _bounded_perturb(input.shift_t_nm, config.shift_t_nm.fine.step,
                          config.shift_t_nm.coarse)
    su = _bounded_perturb(input.shift_u_nm, config.shift_u_nm.fine.step,
                          config.shift_u_nm.coarse)
    rd = _bounded_perturb(registration.rotation_deg, config.rotation_deg.fine.step,
                          config.rotation_deg.coarse)
    bl = _bounded_perturb(input.blur_sigma_nm, config.blur_sigma_nm.fine.step,
                          config.blur_sigma_nm.coarse)
    return [
        (label="shift_t_nm", input=_replace_frozen_input(input; shift_t_nm=st)),
        (label="shift_u_nm", input=_replace_frozen_input(input; shift_u_nm=su)),
        (label="rotation_deg", input=_replace_frozen_input(input;
            theta_total=input.theta_total + deg2rad(rd - registration.rotation_deg))),
        (label="blur_sigma_nm", input=_replace_frozen_input(input; blur_sigma_nm=bl)),
    ]
end

function _score_with_input(input, set, case, config, half_nm, step_nm)
    return score_frozen_contrast(input, set, case.observed, case.backbone,
        case.centers, case.xs, case.ys; half_nm=half_nm, step_nm=step_nm,
        config=config)
end

function _stable_perturbation(nominal, perturbed)
    perturbed.best.sequence == nominal.best.sequence || return false
    unstable_reasons = ("nonfinite_evidence", "best_runner_tie", "complement_tie")
    return !any(reason -> reason in perturbed.abstention_reasons, unstable_reasons)
end

function _strictly_better(value::Float64, control::Float64, config)
    isfinite(value) && isfinite(control) || return false
    return value > control &&
        !isapprox(value, control; rtol=config.tie_rtol, atol=config.tie_atol)
end

"""Run nominal scoring, four one-fine-step perturbations, and one shifted
contrast control without re-running registration."""
function run_frozen_real_controls(case::WholeRoiObservedCase,
                                  set::CommonContrastSet;
                                  config, half_nm::Real, step_nm::Real,
                                  control_shift_nm::Real)
    nominal = run_whole_roi_common_contrast(case, set; config=config,
        half_nm=half_nm, step_nm=step_nm)
    registration_tsv = freeze_registration_tsv(nominal.registration)
    frozen = frozen_scorer_input(nominal.registration; theta_base=case.theta_base)

    perturbations = FrozenControlScore[]
    for item in _perturbed_inputs(frozen, nominal.registration, config)
        scoring = _score_with_input(item.input, set, case, config, half_nm, step_nm)
        push!(perturbations, FrozenControlScore(item.label, scoring, registration_tsv))
    end
    perturbation_stable = all(item ->
        _stable_perturbation(nominal.scoring, item.scoring), perturbations)

    axis = (cos(case.theta_base), sin(case.theta_base))
    shifted_centers = [(x + Float64(control_shift_nm) * axis[1],
                        y + Float64(control_shift_nm) * axis[2])
                       for (x, y) in case.centers]
    shifted_case = WholeRoiObservedCase(case.observed, case.backbone,
        shifted_centers, case.xs, case.ys, case.theta_base)
    shifted_score = _score_with_input(frozen, set, shifted_case, config, half_nm, step_nm)
    controls = [FrozenControlScore("shifted_contrast", shifted_score, registration_tsv)]
    negative_control_beaten = _strictly_better(
        nominal.scoring.incremental_contrast_gain,
        shifted_score.incremental_contrast_gain, config)

    reasons = copy(nominal.scoring.abstention_reasons)
    isempty(nominal.boundary_parameters) || push!(reasons, "registration_boundary")
    perturbation_stable || push!(reasons, "transform_instability")
    negative_control_beaten || push!(reasons, "negative_control_not_beaten")
    unique!(reasons)
    return FrozenRealDiagnosticBundle(nominal, perturbations, controls,
        registration_tsv, perturbation_stable, negative_control_beaten,
        !isempty(reasons), reasons)
end
