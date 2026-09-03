module StructuredUniverse

using Dates
using Printf
using Random
using SHA
using Statistics
using TOML

module InputBoundary
include(joinpath(@__DIR__, "firewall.jl"))
end

export UniverseError,
       LobeKey,
       ScanRow,
       FrozenUniverse,
       NoiseScale,
       CANDIDATE_CONFIG_SHA256,
       MODEL_CONFIG_SHA256,
       FORWARD_PRODUCER_SHA256,
       BACKWARD_PRODUCER_SHA256,
       GRID_SHA256,
       GRID_COORDINATES,
       CONTRAST_FACTORS,
       NOISE_MAD_FRACTION,
       sha256_file,
       parse_scan_date,
       freeze_universe,
       validate_binding,
       bundle_bytes,
       publish_bundle,
       validate_bundle,
       validate_patch_keys,
       publish_patch_pair_receipt,
       drop_channel,
       shift_patch,
       contrast_patch,
       median_absolute_deviation,
       correlated_noise,
       freeze_noise_scale,
       apply_noise_patch

const CANDIDATE_CONFIG_SHA256 = "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const MODEL_CONFIG_SHA256 = "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const FORWARD_PRODUCER_SHA256 = "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad"
const BACKWARD_PRODUCER_SHA256 = "4f866a1e27289a0cf010a6ed4b5e0b92454dc55f97807d33cc95edb9e87182ba"
const GRID_SHA256 = "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"
const GRID_COORDINATES = (
    -0.32, -0.28, -0.24, -0.20, -0.16, -0.12, -0.08, -0.04, 0.00,
    0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32,
)
const CONTRAST_FACTORS = (0.95, 1.05)
const NOISE_MAD_FRACTION = 0.25
const GRID_SIDE = 17
const GRID_STEP_NM = 0.04
const SHIFT_NM = GRID_STEP_NM / 2
const OUTPUT_NAMES = (
    "folds.tsv",
    "patch_keys.tsv",
    "perturbations.tsv",
    "receipt.toml",
    "scan_universe.tsv",
    "seeds.tsv",
    "shards.tsv",
)
const FEATURE_COLUMNS = Set([
    "file", "N", "lobe", "amplitude", "x_nm", "y_nm", "t_nm", "u_nm",
    "sigma_parallel_nm", "sigma_perp_nm", "spacing_prev_nm", "spacing_next_nm",
    "spacing_asym", "amp_rel", "skew_ratio", "axis_x", "axis_y", "origin_x_nm",
    "origin_y_nm", "baseline", "tilt_x", "tilt_y", "gcv", "source", "integrated",
    "amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "bwd_neg_com_t",
    "bwd_neg_diag45", "split_log_skew", "chain_curv_deg", "elongation",
    "bwd_res_en", "bwd_res_com_u", "bwd_res_com_t", "bwd_res_sku", "bwd_res_skt",
    "bwd_res_asym_u", "bwd_res_asym_t", "bwd_res_kurt_u", "bwd_res_ring1",
    "bwd_res_ring2", "bwd_res_ring3", "diff_res_en", "diff_res_com_u",
    "diff_res_com_t", "diff_res_sku", "diff_res_skt", "diff_res_asym_u",
    "diff_res_asym_t", "diff_res_kurt_u", "diff_res_ring1", "diff_res_ring2",
    "diff_res_ring3", "diff_raw_en", "diff_raw_com_u", "diff_raw_com_t",
    "diff_raw_sku", "diff_raw_skt", "diff_raw_asym_u", "diff_raw_asym_t",
    "bwd_raw_en", "bwd_raw_com_u", "bwd_raw_com_t", "bwd_raw_sku",
    "bwd_raw_skt", "bwd_raw_asym_u", "bwd_raw_asym_t",
])

struct UniverseError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, error::UniverseError) =
    print(io, "BLOCKED ", error.code, ": ", error.message)

_fail(code::Symbol, message::AbstractString) = throw(UniverseError(code, String(message)))

struct LobeKey
    file::String
    lobe::Int
end

Base.isless(left::LobeKey, right::LobeKey) =
    (left.file, left.lobe) < (right.file, right.lobe)

struct ScanRow
    file::String
    date::String
    lobe_count::Int
end

struct BindingSnapshot
    role::String
    path::String
    sha256::String
    byte_count::Int
    device::UInt64
    inode::UInt64
    mtime::Float64
    ctime::Float64
end

struct FrozenUniverse
    root::String
    features_path::String
    candidate_config_path::String
    model_config_path::String
    feature_sha256::String
    candidate_config_sha256::String
    model_config_sha256::String
    source_sha256::String
    forward_producer_sha256::String
    backward_producer_sha256::String
    grid_sha256::String
    scans::Tuple{Vararg{ScanRow}}
    keys::Tuple{Vararg{LobeKey}}
    keys_sha256::String
    views::Tuple{Vararg{String}}
    bootstrap_seeds::Tuple{Vararg{Int}}
    shuffle_seeds::Tuple{Vararg{Int}}
    restart_seeds::Tuple{Vararg{Int}}
    snapshots::Tuple{Vararg{BindingSnapshot}}
end

struct PatchKeyValidation
    row_count::Int
    sha256::String
    keys_sha256::String
end

struct NoiseScale
    fold_date::String
    channel::String
    training_files::Tuple{Vararg{String}}
    training_scan_mads::Tuple{Vararg{Float64}}
    median_scan_mad::Float64
    target_mad::Float64
    binding_sha256::String
end

_sha256(bytes) = bytes2hex(sha256(bytes))
sha256_file(path::AbstractString) = _sha256(read(String(path)))

function _validate_hash(value::AbstractString, context::AbstractString)::String
    text = String(value)
    occursin(r"^[0-9a-f]{64}$", text) ||
        _fail(:invalid_hash, "$context must be a lowercase SHA-256")
    return text
end

function _canonical_root(root::AbstractString)::String
    supplied = String(root)
    isempty(supplied) && _fail(:invalid_root, "root is empty")
    isdir(supplied) || _fail(:invalid_root, "root is not a directory")
    islink(supplied) && _fail(:symlink_rejected, "root is a symlink")
    absolute = normpath(abspath(supplied))
    canonical = realpath(supplied)
    absolute == canonical || _fail(:noncanonical_path, "root is not canonical")
    return canonical
end

function _validate_relative(path::String, context::String)
    isempty(path) && _fail(:noncanonical_path, "$context is empty")
    occursin('\0', path) && _fail(:noncanonical_path, "$context contains NUL")
    occursin('\\', path) && _fail(:noncanonical_path, "$context uses a backslash")
    isabspath(path) && _fail(:path_escape, "$context must be root-relative")
    components = split(path, '/'; keepempty=true)
    any(component -> component in ("", ".", ".."), components) &&
        _fail(:path_escape, "$context contains a noncanonical component")
    normpath(path) == path || _fail(:noncanonical_path, "$context is not canonical")
    return nothing
end

function _relative_inside(root::String, absolute::String, context::String)::String
    relative = relpath(absolute, root)
    parts = splitpath(relative)
    (relative == "." || (!isempty(parts) && first(parts) == "..")) &&
        _fail(:path_escape, "$context escapes root")
    normalized = join(parts, "/")
    _validate_relative(normalized, context)
    return normalized
end

function _reject_symlink_components(root::String, relative::String, context::String)
    cursor = root
    for component in split(relative, '/')
        cursor = joinpath(cursor, component)
        islink(cursor) && _fail(:symlink_rejected, "$context contains a symlink component")
    end
    return nothing
end

function _resolve_existing(root::String, supplied::AbstractString, context::String)::Tuple{String,String}
    text = String(supplied)
    _validate_relative(text, context)
    absolute = normpath(joinpath(root, split(text, '/')...))
    relative = _relative_inside(root, absolute, context)
    _reject_symlink_components(root, relative, context)
    isfile(absolute) || _fail(:missing_input, "$context is not a regular file")
    return absolute, relative
end

function _reject_absolute_symlinks(path::String, context::String)
    absolute = normpath(abspath(path))
    cursor = first(splitdrive(absolute))
    isempty(cursor) && (cursor = "/")
    for component in splitpath(absolute)
        component == "/" && continue
        cursor = cursor == "/" ? joinpath(cursor, component) : joinpath(cursor, component)
        islink(cursor) && _fail(:symlink_rejected, "$context contains a symlink component")
    end
    return nothing
end

function _resolve_source(path::AbstractString)::String
    supplied = String(path)
    isabspath(supplied) || _fail(:noncanonical_path, "source path must be absolute")
    absolute = normpath(abspath(supplied))
    _reject_absolute_symlinks(absolute, "source path")
    isfile(absolute) || _fail(:missing_input, "source path is not a regular file")
    realpath(absolute) == absolute || _fail(:noncanonical_path, "source path is not canonical")
    return absolute
end

function _resolve_output(root::String, supplied::AbstractString, context::String)::String
    text = String(supplied)
    absolute = if isabspath(text)
        normpath(text)
    else
        _validate_relative(text, context)
        normpath(joinpath(root, split(text, '/')...))
    end
    relative = _relative_inside(root, absolute, context)
    parent_relative = dirname(relative)
    parent_relative == "." || _reject_symlink_components(root, parent_relative, context)
    parent = dirname(absolute)
    isdir(parent) || _fail(:missing_output_parent, "$context parent does not exist")
    islink(parent) && _fail(:symlink_rejected, "$context parent is a symlink")
    islink(absolute) && _fail(:symlink_rejected, "$context is a symlink")
    return absolute
end

function _snapshot(role::String, path::String)::Tuple{Vector{UInt8},BindingSnapshot}
    islink(path) && _fail(:symlink_rejected, "$role is a symlink")
    isfile(path) || _fail(:missing_input, "$role is not a regular file")
    before = stat(path)
    bytes = read(path)
    islink(path) && _fail(:symlink_rejected, "$role became a symlink")
    isfile(path) || _fail(:dependency_changed, "$role stopped being a regular file")
    after = stat(path)
    stable = before.device == after.device && before.inode == after.inode &&
             before.size == after.size && before.mtime == after.mtime &&
             before.ctime == after.ctime && length(bytes) == after.size
    stable || _fail(:dependency_changed, "$role changed while read")
    snapshot = BindingSnapshot(
        role,
        path,
        _sha256(bytes),
        length(bytes),
        UInt64(after.device),
        UInt64(after.inode),
        Float64(after.mtime),
        Float64(after.ctime),
    )
    return bytes, snapshot
end

function _verify_snapshot(snapshot::BindingSnapshot)
    bytes, current = _snapshot(snapshot.role, snapshot.path)
    same = current.sha256 == snapshot.sha256 &&
           current.byte_count == snapshot.byte_count &&
           current.device == snapshot.device && current.inode == snapshot.inode &&
           current.mtime == snapshot.mtime && current.ctime == snapshot.ctime
    same || _fail(:dependency_changed, "$(snapshot.role) changed after binding")
    return bytes
end

function validate_binding(frozen::FrozenUniverse)
    for snapshot in frozen.snapshots
        _verify_snapshot(snapshot)
    end
    return nothing
end

function parse_scan_date(file::AbstractString)::String
    text = String(file)
    strip(text) == text || _fail(:noncanonical_scan_name, "scan name has surrounding space")
    isempty(text) && _fail(:invalid_scan_name, "scan name is empty")
    all(isascii, text) || _fail(:noncanonical_scan_name, "scan name is not ASCII")
    (basename(text) == text && !occursin('/', text) && !occursin('\\', text)) ||
        _fail(:noncanonical_scan_name, "scan name must be a basename")
    tokens = [match.match for match in eachmatch(r"(?<![0-9])[0-9]{8}(?![0-9])", text)]
    isempty(tokens) && _fail(:invalid_scan_name, "scan name lacks a leading date token")
    length(tokens) == 1 || _fail(:ambiguous_scan_date, "scan name has multiple date tokens")
    startswith(text, only(tokens)) || _fail(:invalid_scan_name, "date token is not leading")
    occursin(r"^[0-9]{8}(?:_[A-Za-z0-9][A-Za-z0-9_.-]*)?\.sxm$", text) ||
        _fail(:noncanonical_scan_name, "scan name does not use the canonical form")
    token = only(tokens)
    try
        Date(token, dateformat"yyyymmdd")
    catch
        _fail(:invalid_scan_date, "scan date is invalid")
    end
    return token
end

function _parse_text_table(bytes::Vector{UInt8})
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "table is not valid UTF-8")
    startswith(text, '\ufeff') && _fail(:invalid_encoding, "table has a byte-order mark")
    occursin('\r', text) && _fail(:carriage_return, "table must use LF line endings")
    endswith(text, '\n') || _fail(:missing_final_lf, "table must end with LF")
    lines = split(text, '\n'; keepempty=true)
    pop!(lines)
    isempty(lines) && _fail(:empty_input, "table is empty")
    any(isempty, lines) && _fail(:empty_row, "table contains an empty row")
    header = split(first(lines), '\t'; keepempty=true)
    any(isempty, header) && _fail(:empty_column, "table has an empty column")
    length(header) == length(unique(header)) ||
        _fail(:duplicate_column, "table has duplicate columns")
    for column in header
        strip(column) == column || _fail(:noncanonical_column, "column has surrounding space")
        all(isascii, column) || _fail(:noncanonical_column, "column is not ASCII")
        occursin(r"^[A-Za-z][A-Za-z0-9_]*$", column) ||
            _fail(:noncanonical_column, "column name is not canonical")
        column == "N" || InputBoundary._has_forbidden_concept(column) &&
            _fail(:forbidden_input, "column is outside the construction boundary")
    end
    rows = Vector{Vector{String}}()
    for (row_index, line) in enumerate(lines[2:end])
        fields = split(line, '\t'; keepempty=true)
        length(fields) == length(header) ||
            _fail(:row_width_mismatch, "row $row_index has the wrong width")
        for field in fields
            strip(field) == field || _fail(:noncanonical_value, "row value has surrounding space")
            InputBoundary._has_forbidden_concept(field) &&
                _fail(:forbidden_input, "row value is outside the construction boundary")
        end
        push!(rows, fields)
    end
    isempty(rows) && _fail(:empty_input, "table has no data rows")
    return header, rows
end

function _positive_integer(text::String, code::Symbol, context::String)::Int
    if !occursin(r"^[1-9][0-9]*$", text)
        parsed = tryparse(Int, text)
        parsed !== nothing && parsed > 0 &&
            _fail(:noncanonical_integer, "$context is not canonically serialized")
        _fail(code, "$context is not a canonical positive integer")
    end
    value = tryparse(Int, text)
    value === nothing && _fail(code, "$context exceeds the integer range")
    return value
end

function _parse_features(bytes::Vector{UInt8})
    header, rows = _parse_text_table(bytes)
    for required in ("file", "N", "lobe")
        required in header || _fail(:missing_column, "feature table lacks $required")
    end
    all(column -> column in FEATURE_COLUMNS, header) ||
        _fail(:forbidden_input, "feature table contains a column outside the declared schema")
    file_index = findfirst(==("file"), header)::Int
    count_index = findfirst(==("N"), header)::Int
    lobe_index = findfirst(==("lobe"), header)::Int
    row_seen = Set{String}()
    key_seen = Dict{LobeKey,String}()
    per_file = Dict{String,Vector{Tuple{Int,Int}}}()
    for fields in rows
        encoded = join(fields, '\t')
        encoded in row_seen && _fail(:duplicate_row, "feature table repeats a row")
        push!(row_seen, encoded)
        file = fields[file_index]
        parse_scan_date(file)
        count = _positive_integer(fields[count_index], :invalid_integer, "lobe count")
        lobe = _positive_integer(fields[lobe_index], :invalid_integer, "lobe key")
        key = LobeKey(file, lobe)
        haskey(key_seen, key) && _fail(:duplicate_key, "feature table repeats a file/lobe key")
        key_seen[key] = encoded
        push!(get!(per_file, file, Tuple{Int,Int}[]), (lobe, count))
    end

    scans = ScanRow[]
    key_rows = LobeKey[]
    for file in sort!(collect(Base.keys(per_file)))
        entries = per_file[file]
        counts = unique(last(entry) for entry in entries)
        length(counts) == 1 || _fail(:lobe_count_mismatch, "scan has inconsistent lobe counts")
        count = only(counts)
        lobes = sort!(collect(first(entry) for entry in entries))
        lobes == collect(1:count) || _fail(:lobe_key_mismatch, "scan keys do not equal 1:lobe_count")
        length(entries) == count || _fail(:lobe_key_mismatch, "scan row count differs from lobe_count")
        push!(scans, ScanRow(file, parse_scan_date(file), count))
        append!(key_rows, [LobeKey(file, lobe) for lobe in lobes])
    end
    length(unique(row.date for row in scans)) >= 2 ||
        _fail(:insufficient_dates, "at least two scan dates are required")
    return Tuple(scans), Tuple(key_rows)
end

function _table(dictionary, key::String, context::String)
    value = get(dictionary, key, nothing)
    value isa AbstractDict || _fail(:malformed_config, "$context.$key is not a table")
    return value
end

function _config_string(dictionary, key::String, context::String)::String
    value = get(dictionary, key, nothing)
    value isa String || _fail(:malformed_config, "$context.$key is not a string")
    return value
end

function _config_int(dictionary, key::String, context::String)::Int
    value = get(dictionary, key, nothing)
    value isa Integer || _fail(:malformed_config, "$context.$key is not an integer")
    return Int(value)
end

function _config_float(dictionary, key::String, context::String)::Float64
    value = get(dictionary, key, nothing)
    value isa Real || _fail(:malformed_config, "$context.$key is not numerical")
    number = Float64(value)
    isfinite(number) || _fail(:malformed_config, "$context.$key is nonfinite")
    return number
end

function _config_strings(dictionary, key::String, context::String)::Vector{String}
    value = get(dictionary, key, nothing)
    value isa AbstractVector || _fail(:malformed_config, "$context.$key is not an array")
    all(item -> item isa String, value) ||
        _fail(:malformed_config, "$context.$key must contain strings")
    result = String.(value)
    length(result) == length(unique(result)) ||
        _fail(:malformed_config, "$context.$key contains duplicates")
    return result
end

function _grid_contract_bytes()::Vector{UInt8}
    text = "version=1\n" *
           "half_t_nm=0.32\n" *
           "half_u_nm=0.32\n" *
           "step_nm=0.04\n" *
           "order=u_outer,t_inner\n" *
           "coordinates_nm=-0.32,-0.28,-0.24,-0.20,-0.16,-0.12,-0.08,-0.04,0.00,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32\n" *
           "forward_prefix=res_p\n" *
           "backward_prefix=bwd_res_p\n" *
           "pixels=289\n"
    return Vector{UInt8}(codeunits(text))
end

function _parse_config(bytes::Vector{UInt8}, context::String)
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "$context is not valid UTF-8")
    try
        return TOML.parse(text)
    catch error
        _fail(:malformed_config, "$context cannot be parsed: $(sprint(showerror, error))")
    end
end

function _seed_range(table, prefix::String, context::String)::Tuple{Vararg{Int}}
    count = _config_int(table, prefix == "seed" ? "count" : "edge_shuffle_count", context)
    start_key = prefix == "seed" ? "seed_start" : "edge_shuffle_seed_start"
    last_key = prefix == "seed" ? "seed_stop" : "edge_shuffle_seed_stop"
    first_seed = _config_int(table, start_key, context)
    last_seed = _config_int(table, last_key, context)
    first_seed >= 0 || _fail(:malformed_config, "$context seed range starts below zero")
    last_seed >= first_seed || _fail(:malformed_config, "$context seed range is reversed")
    count == last_seed - first_seed + 1 ||
        _fail(:malformed_config, "$context seed count differs from its range")
    return Tuple(first_seed:last_seed)
end

function _source_bundle_sha(source_snapshots::Vector{BindingSnapshot})::String
    names = basename.(getfield.(source_snapshots, :path))
    length(names) == length(unique(names)) ||
        _fail(:source_collision, "source basenames are not unique")
    io = IOBuffer()
    for snapshot in sort(source_snapshots; by=item -> basename(item.path))
        println(io, basename(snapshot.path), '\t', snapshot.sha256)
    end
    return _sha256(take!(io))
end

function _key_identity_bytes(keys)::Vector{UInt8}
    io = IOBuffer()
    println(io, "schema_version\tfile\tlobe")
    for key in keys
        println(io, "1\t", key.file, '\t', key.lobe)
    end
    return take!(io)
end

function freeze_universe(
    root::AbstractString;
    features::AbstractString,
    feature_sha256::AbstractString,
    candidate_config::AbstractString,
    model_config::AbstractString,
    source_paths,
    candidate_config_sha256::AbstractString=CANDIDATE_CONFIG_SHA256,
    model_config_sha256::AbstractString=MODEL_CONFIG_SHA256,
)::FrozenUniverse
    canonical_root = _canonical_root(root)
    required_feature_hash = _validate_hash(feature_sha256, "feature hash")
    required_candidate_hash = _validate_hash(candidate_config_sha256, "candidate config hash")
    required_model_hash = _validate_hash(model_config_sha256, "model config hash")
    features_path, _ = _resolve_existing(canonical_root, features, "feature input")
    candidate_path, _ = _resolve_existing(canonical_root, candidate_config, "candidate config")
    model_path, _ = _resolve_existing(canonical_root, model_config, "model config")

    feature_bytes, feature_snapshot = _snapshot("feature input", features_path)
    candidate_bytes, candidate_snapshot = _snapshot("candidate config", candidate_path)
    model_bytes, model_snapshot = _snapshot("model config", model_path)
    feature_snapshot.sha256 == required_feature_hash ||
        _fail(:feature_hash_mismatch, "feature bytes do not match the declared hash")
    candidate_snapshot.sha256 == required_candidate_hash ||
        _fail(:config_hash_mismatch, "candidate config bytes do not match the required hash")
    model_snapshot.sha256 == required_model_hash ||
        _fail(:config_hash_mismatch, "model config bytes do not match the required hash")

    candidate = _parse_config(candidate_bytes, "candidate config")
    model = _parse_config(model_bytes, "model config")
    provenance = _table(candidate, "provenance", "candidate")
    declared_model_hash = _config_string(provenance, "model_config_sha256", "candidate.provenance")
    declared_model_hash == required_model_hash ||
        _fail(:config_binding_mismatch, "candidate does not bind the supplied model config")

    forward_relative = _config_string(provenance, "forward_patch_producer", "candidate.provenance")
    backward_relative = _config_string(provenance, "backward_patch_producer", "candidate.provenance")
    declared_forward_hash = _config_string(provenance, "forward_patch_producer_sha256", "candidate.provenance")
    declared_backward_hash = _config_string(provenance, "backward_patch_producer_sha256", "candidate.provenance")
    declared_forward_hash == FORWARD_PRODUCER_SHA256 ||
        _fail(:declaration_mismatch, "forward producer declaration is stale")
    declared_backward_hash == BACKWARD_PRODUCER_SHA256 ||
        _fail(:declaration_mismatch, "backward producer declaration is stale")
    forward_path, _ = _resolve_existing(canonical_root, forward_relative, "forward producer")
    backward_path, _ = _resolve_existing(canonical_root, backward_relative, "backward producer")
    _, forward_snapshot = _snapshot("forward producer", forward_path)
    _, backward_snapshot = _snapshot("backward producer", backward_path)
    forward_snapshot.sha256 == FORWARD_PRODUCER_SHA256 ||
        _fail(:source_hash_mismatch, "forward producer bytes are stale")
    backward_snapshot.sha256 == BACKWARD_PRODUCER_SHA256 ||
        _fail(:source_hash_mismatch, "backward producer bytes are stale")

    spacing_relative = _config_string(provenance, "spacing_source", "candidate.provenance")
    spacing_hash = _validate_hash(
        _config_string(provenance, "spacing_source_sha256", "candidate.provenance"),
        "spacing source hash",
    )
    spacing_path, _ = _resolve_existing(canonical_root, spacing_relative, "spacing source")
    _, spacing_snapshot = _snapshot("spacing source", spacing_path)
    spacing_snapshot.sha256 == spacing_hash ||
        _fail(:source_hash_mismatch, "spacing source bytes are stale")

    preprocessing = _table(model, "preprocessing", "model config")
    patch = _table(preprocessing, "patch", "model.preprocessing")
    _config_string(patch, "grid_sha256", "model.preprocessing.patch") == GRID_SHA256 ||
        _fail(:declaration_mismatch, "grid declaration is stale")
    _sha256(_grid_contract_bytes()) == GRID_SHA256 ||
        _fail(:internal_contract, "grid bytes do not match the frozen digest")
    _config_float(patch, "half_t_nm", "model.preprocessing.patch") == 0.32 ||
        _fail(:declaration_mismatch, "grid half_t_nm differs")
    _config_float(patch, "half_u_nm", "model.preprocessing.patch") == 0.32 ||
        _fail(:declaration_mismatch, "grid half_u_nm differs")
    _config_float(patch, "step_nm", "model.preprocessing.patch") == GRID_STEP_NM ||
        _fail(:declaration_mismatch, "grid step differs")
    _config_int(patch, "pixels_per_axis", "model.preprocessing.patch") == GRID_SIDE ||
        _fail(:declaration_mismatch, "grid side differs")
    _config_int(patch, "pixels", "model.preprocessing.patch") == GRID_SIDE^2 ||
        _fail(:declaration_mismatch, "grid pixel total differs")
    _config_string(patch, "order", "model.preprocessing.patch") == "u_outer,t_inner" ||
        _fail(:declaration_mismatch, "grid order differs")
    coordinates_any = get(patch, "coordinates_nm", nothing)
    coordinates_any isa AbstractVector ||
        _fail(:malformed_config, "model.preprocessing.patch.coordinates_nm is not an array")
    Tuple(Float64.(coordinates_any)) == GRID_COORDINATES ||
        _fail(:declaration_mismatch, "grid coordinates differ")
    _config_string(patch, "forward_prefix", "model.preprocessing.patch") == "res_p" ||
        _fail(:declaration_mismatch, "forward prefix differs")
    _config_string(patch, "backward_prefix", "model.preprocessing.patch") == "bwd_res_p" ||
        _fail(:declaration_mismatch, "backward prefix differs")

    producers = _table(preprocessing, "patch_producers", "model.preprocessing")
    _config_string(producers, "forward_source_sha256", "model.preprocessing.patch_producers") ==
        FORWARD_PRODUCER_SHA256 || _fail(:declaration_mismatch, "model forward digest differs")
    _config_string(producers, "backward_source_sha256", "model.preprocessing.patch_producers") ==
        BACKWARD_PRODUCER_SHA256 || _fail(:declaration_mismatch, "model backward digest differs")

    model_section = _table(model, "model", "model config")
    views_table = _table(model_section, "views", "model")
    views = Tuple(_config_strings(views_table, "names", "model.views"))
    isempty(views) && _fail(:malformed_config, "model.views.names is empty")
    restart_count = _config_int(model_section, "n_starts", "model")
    restart_count > 0 || _fail(:malformed_config, "model.n_starts is not positive")
    restart_seeds = Tuple(0:(restart_count - 1))
    selection = _table(model, "selection", "model config")
    bootstrap = _table(selection, "bootstrap", "model.selection")
    bootstrap_seeds = _seed_range(bootstrap, "seed", "model.selection.bootstrap")
    shuffle_seeds = _seed_range(bootstrap, "edge", "model.selection.bootstrap")

    scans, keys_tuple = _parse_features(feature_bytes)
    source_snapshots = BindingSnapshot[]
    for supplied in source_paths
        source_path = _resolve_source(String(supplied))
        _, snapshot = _snapshot("source $(basename(source_path))", source_path)
        push!(source_snapshots, snapshot)
    end
    isempty(source_snapshots) && _fail(:missing_input, "no source paths were supplied")
    source_hash = _source_bundle_sha(source_snapshots)
    snapshots = BindingSnapshot[
        feature_snapshot,
        candidate_snapshot,
        model_snapshot,
        forward_snapshot,
        backward_snapshot,
        spacing_snapshot,
    ]
    append!(snapshots, source_snapshots)
    unique_paths = Set{String}()
    for snapshot in snapshots
        snapshot.path in unique_paths && _fail(:source_collision, "a dependency path was bound twice")
        push!(unique_paths, snapshot.path)
    end
    key_hash = _sha256(_key_identity_bytes(keys_tuple))
    frozen = FrozenUniverse(
        canonical_root,
        features_path,
        candidate_path,
        model_path,
        required_feature_hash,
        required_candidate_hash,
        required_model_hash,
        source_hash,
        FORWARD_PRODUCER_SHA256,
        BACKWARD_PRODUCER_SHA256,
        GRID_SHA256,
        scans,
        keys_tuple,
        key_hash,
        views,
        bootstrap_seeds,
        shuffle_seeds,
        restart_seeds,
        Tuple(snapshots),
    )
    validate_binding(frozen)
    return frozen
end

_format_number(value::Float64) = value == 0.0 ? "0" : @sprintf("%.17g", value)

function _scan_universe_bytes(frozen::FrozenUniverse)::Vector{UInt8}
    io = IOBuffer()
    println(io, join((
        "schema_version", "file", "date", "lobe_count", "feature_sha256",
        "candidate_config_sha256", "model_config_sha256", "source_sha256",
        "forward_producer_sha256", "backward_producer_sha256", "grid_sha256",
    ), '\t'))
    for row in frozen.scans
        println(io, join((
            "1", row.file, row.date, string(row.lobe_count), frozen.feature_sha256,
            frozen.candidate_config_sha256, frozen.model_config_sha256,
            frozen.source_sha256, frozen.forward_producer_sha256,
            frozen.backward_producer_sha256, frozen.grid_sha256,
        ), '\t'))
    end
    return take!(io)
end

function _patch_key_bytes(frozen::FrozenUniverse, universe_sha256::String)::Vector{UInt8}
    io = IOBuffer()
    println(io, "schema_version\tfile\tlobe\tfeature_sha256\tuniverse_sha256\tkeys_sha256")
    for key in frozen.keys
        println(io, join((
            "1", key.file, string(key.lobe), frozen.feature_sha256,
            universe_sha256, frozen.keys_sha256,
        ), '\t'))
    end
    return take!(io)
end

function _fold_bytes(frozen::FrozenUniverse, universe_sha256::String)::Vector{UInt8}
    io = IOBuffer()
    println(io, join((
        "schema_version", "fold_date", "scan_file", "scan_date", "role", "lobe_count",
        "universe_sha256", "keys_sha256", "model_config_sha256", "source_sha256",
    ), '\t'))
    dates = sort!(unique(row.date for row in frozen.scans))
    for fold_date in dates, row in frozen.scans
        role = row.date == fold_date ? "heldout" : "training"
        println(io, join((
            "1", fold_date, row.file, row.date, role, string(row.lobe_count),
            universe_sha256, frozen.keys_sha256, frozen.model_config_sha256,
            frozen.source_sha256,
        ), '\t'))
    end
    return take!(io)
end

function _seed_groups(frozen::FrozenUniverse)
    return (
        "bootstrap" => frozen.bootstrap_seeds,
        "edge_shuffle" => frozen.shuffle_seeds,
        "noise" => frozen.bootstrap_seeds,
        "restart" => frozen.restart_seeds,
    )
end

function _seed_bytes(frozen::FrozenUniverse, universe_sha256::String)::Vector{UInt8}
    io = IOBuffer()
    println(io, join((
        "schema_version", "fold_date", "seed_family", "seed", "universe_sha256",
        "keys_sha256", "model_config_sha256", "source_sha256",
    ), '\t'))
    dates = sort!(unique(row.date for row in frozen.scans))
    for fold_date in dates, (family, seeds) in _seed_groups(frozen), seed in seeds
        println(io, join((
            "1", fold_date, family, string(seed), universe_sha256,
            frozen.keys_sha256, frozen.model_config_sha256, frozen.source_sha256,
        ), '\t'))
    end
    return take!(io)
end

function _perturbation_entries(frozen::FrozenUniverse)
    entries = NamedTuple[]
    push!(entries, (
        perturbation="unperturbed", view="NA", seed_family="NA", seed="NA",
        delta_t_nm="0", delta_u_nm="0", factor="1", noise_mad_fraction="0",
    ))
    for view in frozen.views
        push!(entries, (
            perturbation="channel_drop", view=view, seed_family="NA", seed="NA",
            delta_t_nm="0", delta_u_nm="0", factor="1", noise_mad_fraction="0",
        ))
    end
    for (name, delta_t, delta_u) in (
        ("shift_t_minus", -SHIFT_NM, 0.0),
        ("shift_t_plus", SHIFT_NM, 0.0),
        ("shift_u_minus", 0.0, -SHIFT_NM),
        ("shift_u_plus", 0.0, SHIFT_NM),
    )
        push!(entries, (
            perturbation=name, view="NA", seed_family="NA", seed="NA",
            delta_t_nm=_format_number(delta_t), delta_u_nm=_format_number(delta_u),
            factor="1", noise_mad_fraction="0",
        ))
    end
    for (name, factor) in zip(("contrast_low", "contrast_high"), CONTRAST_FACTORS)
        push!(entries, (
            perturbation=name, view="NA", seed_family="NA", seed="NA",
            delta_t_nm="0", delta_u_nm="0", factor=factor == 0.95 ? "0.95" : "1.05",
            noise_mad_fraction="0",
        ))
    end
    for seed in frozen.bootstrap_seeds
        push!(entries, (
            perturbation="noise", view="per_channel", seed_family="noise", seed=string(seed),
            delta_t_nm="0", delta_u_nm="0", factor="1",
            noise_mad_fraction="0.25",
        ))
    end
    for seed in frozen.restart_seeds
        push!(entries, (
            perturbation="seed_restart", view="NA", seed_family="restart", seed=string(seed),
            delta_t_nm="0", delta_u_nm="0", factor="1", noise_mad_fraction="0",
        ))
    end
    return entries
end

function _perturbation_bytes(frozen::FrozenUniverse, universe_sha256::String)::Vector{UInt8}
    io = IOBuffer()
    println(io, join((
        "schema_version", "fold_date", "scan_file", "perturbation", "view",
        "seed_family", "seed", "delta_t_nm", "delta_u_nm", "factor",
        "noise_mad_fraction", "universe_sha256", "keys_sha256",
        "model_config_sha256", "source_sha256",
    ), '\t'))
    entries = _perturbation_entries(frozen)
    for row in frozen.scans, entry in entries
        println(io, join((
            "1", row.date, row.file, entry.perturbation, entry.view,
            entry.seed_family, entry.seed, entry.delta_t_nm, entry.delta_u_nm,
            entry.factor, entry.noise_mad_fraction, universe_sha256,
            frozen.keys_sha256, frozen.model_config_sha256, frozen.source_sha256,
        ), '\t'))
    end
    return take!(io)
end

function _data_lines(bytes::Vector{UInt8})::Vector{String}
    lines = split(chomp(String(copy(bytes))), '\n')
    return lines[2:end]
end

function _shard_bytes(
    frozen::FrozenUniverse,
    universe_sha256::String,
    folds::Vector{UInt8},
    seeds::Vector{UInt8},
    perturbations::Vector{UInt8},
)::Vector{UInt8}
    records = NamedTuple[]
    for (kind, bytes) in (("fold", folds), ("seed", seeds), ("perturbation", perturbations))
        for line in _data_lines(bytes)
            digest = _sha256(codeunits(line * "\n"))
            identity = _sha256(codeunits(
                "shard-v1\nkind=$kind\nrow_sha256=$digest\nuniverse_sha256=$universe_sha256\n"))
            push!(records, (
                shard_id=identity[1:16], row_kind=kind, row_sha256=digest,
            ))
        end
    end
    sort!(records; by=record -> (record.row_kind, record.row_sha256, record.shard_id))
    length(records) == length(unique(record.shard_id for record in records)) ||
        _fail(:internal_contract, "shard identifiers collide")
    io = IOBuffer()
    println(io, join((
        "schema_version", "shard_id", "row_kind", "row_sha256", "universe_sha256",
        "keys_sha256", "model_config_sha256", "source_sha256",
    ), '\t'))
    for record in records
        println(io, join((
            "1", record.shard_id, record.row_kind, record.row_sha256,
            universe_sha256, frozen.keys_sha256, frozen.model_config_sha256,
            frozen.source_sha256,
        ), '\t'))
    end
    return take!(io)
end

_toml_quote(value::String) = repr(value)

function _receipt_bytes(
    frozen::FrozenUniverse,
    universe_sha256::String,
    artifacts::Dict{String,Vector{UInt8}},
)::Vector{UInt8}
    io = IOBuffer()
    fields = (
        "schema" => "stmfit-structured-scan-universe-receipt-v1",
        "schema_version" => 1,
        "feature_sha256" => frozen.feature_sha256,
        "candidate_config_sha256" => frozen.candidate_config_sha256,
        "model_config_sha256" => frozen.model_config_sha256,
        "source_sha256" => frozen.source_sha256,
        "forward_producer_sha256" => frozen.forward_producer_sha256,
        "backward_producer_sha256" => frozen.backward_producer_sha256,
        "grid_sha256" => frozen.grid_sha256,
        "universe_sha256" => universe_sha256,
        "keys_sha256" => frozen.keys_sha256,
    )
    for (key, value) in fields
        if value isa Integer
            println(io, key, " = ", value)
        else
            println(io, key, " = ", _toml_quote(String(value)))
        end
    end
    println(io, "\n[artifacts]")
    for name in sort!(collect(keys(artifacts)))
        println(io, _toml_quote(name), " = ", _toml_quote(_sha256(artifacts[name])))
    end
    return take!(io)
end

function bundle_bytes(frozen::FrozenUniverse)::Dict{String,Vector{UInt8}}
    validate_binding(frozen)
    universe = _scan_universe_bytes(frozen)
    universe_sha256 = _sha256(universe)
    keys = _patch_key_bytes(frozen, universe_sha256)
    folds = _fold_bytes(frozen, universe_sha256)
    seeds = _seed_bytes(frozen, universe_sha256)
    perturbations = _perturbation_bytes(frozen, universe_sha256)
    shards = _shard_bytes(frozen, universe_sha256, folds, seeds, perturbations)
    artifacts = Dict(
        "scan_universe.tsv" => universe,
        "patch_keys.tsv" => keys,
        "folds.tsv" => folds,
        "seeds.tsv" => seeds,
        "perturbations.tsv" => perturbations,
        "shards.tsv" => shards,
    )
    receipt = _receipt_bytes(frozen, universe_sha256, artifacts)
    result = copy(artifacts)
    result["receipt.toml"] = receipt
    return result
end

function _fsync_directory(path::String)
    descriptor = ccall(:open, Cint, (Cstring, Cint), path, 0)
    descriptor >= 0 || _fail(:fsync_failed, "could not open publication directory")
    try
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            _fail(:fsync_failed, "directory fsync failed")
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
    return nothing
end

function _write_synced(path::String, bytes::Vector{UInt8})
    ispath(path) && _fail(:publication_collision, "publication stage already exists")
    open(path, "w") do io
        write(io, bytes)
        flush(io)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 ||
            _fail(:fsync_failed, "file fsync failed")
    end
    isfile(path) && !islink(path) ||
        _fail(:publication_collision, "publication stage is not a regular file")
    return nothing
end

function _hashes(bytes_by_name::Dict{String,Vector{UInt8}})::Dict{String,String}
    return Dict(name => _sha256(bytes) for (name, bytes) in bytes_by_name)
end

function validate_bundle(frozen::FrozenUniverse, output_directory::AbstractString)::Dict{String,String}
    validate_binding(frozen)
    absolute = _resolve_output(frozen.root, output_directory, "output directory")
    isdir(absolute) || _fail(:missing_output, "output directory does not exist")
    islink(absolute) && _fail(:symlink_rejected, "output directory is a symlink")
    names = sort(readdir(absolute))
    names == collect(OUTPUT_NAMES) || _fail(:bundle_mismatch, "output file set is not exact")
    required = bundle_bytes(frozen)
    for name in names
        path = joinpath(absolute, name)
        isfile(path) && !islink(path) ||
            _fail(:bundle_mismatch, "output member is not a regular file")
        read(path) == required[name] || _fail(:bundle_mismatch, "output member bytes differ")
    end
    return _hashes(required)
end

function publish_bundle(
    frozen::FrozenUniverse,
    output_directory::AbstractString;
    failpoint::Union{Nothing,Symbol}=nothing,
)::Dict{String,String}
    failpoint in (nothing, :after_first_stage, :before_install, :after_install) ||
        _fail(:invalid_failpoint, "unsupported publication failpoint")
    validate_binding(frozen)
    destination = _resolve_output(frozen.root, output_directory, "output directory")
    if ispath(destination)
        return validate_bundle(frozen, destination)
    end
    bytes_by_name = bundle_bytes(frozen)
    parent = dirname(destination)
    stage = mktempdir(parent; prefix=".structured-universe-stage-", cleanup=false)
    installed = false
    try
        for (index, name) in enumerate(sort!(collect(keys(bytes_by_name))))
            _write_synced(joinpath(stage, name), bytes_by_name[name])
            index == 1 && failpoint == :after_first_stage &&
                _fail(:interrupted_publication, "simulated interruption after first stage")
        end
        _fsync_directory(stage)
        validate_binding(frozen)
        failpoint == :before_install &&
            _fail(:interrupted_publication, "simulated interruption before install")
        ispath(destination) && _fail(:publication_collision, "output appeared concurrently")
        Base.Filesystem.rename(stage, destination)
        installed = true
        _fsync_directory(parent)
        failpoint == :after_install &&
            _fail(:interrupted_publication, "simulated interruption after install")
        return validate_bundle(frozen, destination)
    finally
        !installed && ispath(stage) && rm(stage; recursive=true, force=true)
    end
end

function validate_patch_keys(frozen::FrozenUniverse, path::AbstractString)::PatchKeyValidation
    validate_binding(frozen)
    absolute = if isabspath(String(path))
        candidate = normpath(String(path))
        relative = _relative_inside(frozen.root, candidate, "patch table")
        _reject_symlink_components(frozen.root, relative, "patch table")
        isfile(candidate) || _fail(:missing_input, "patch table is missing")
        candidate
    else
        first(_resolve_existing(frozen.root, path, "patch table"))
    end
    bytes, _ = _snapshot("patch table", absolute)
    header, rows = _parse_text_table(bytes)
    for required in ("file", "lobe")
        required in header || _fail(:missing_column, "patch table lacks $required")
    end
    file_index = findfirst(==("file"), header)::Int
    lobe_index = findfirst(==("lobe"), header)::Int
    value_indices = setdiff(eachindex(header), (file_index, lobe_index))
    isempty(value_indices) && _fail(:missing_column, "patch table has no numerical payload")
    observed = LobeKey[]
    seen = Set{LobeKey}()
    for fields in rows
        file = fields[file_index]
        parse_scan_date(file)
        lobe = _positive_integer(fields[lobe_index], :invalid_integer, "patch lobe key")
        key = LobeKey(file, lobe)
        key in seen && _fail(:duplicate_patch_key, "patch table repeats a key")
        push!(seen, key)
        push!(observed, key)
        for index in value_indices
            value = fields[index]
            value == "NA" && continue
            parsed = tryparse(Float64, value)
            (parsed !== nothing && isfinite(parsed)) ||
                _fail(:invalid_patch_value, "patch payload is not finite or NA")
        end
    end
    observed_tuple = Tuple(observed)
    if Set(observed) == Set(frozen.keys) && length(observed) == length(frozen.keys)
        observed_tuple == frozen.keys || _fail(:patch_key_order, "patch key order is not canonical")
    else
        _fail(:patch_key_mismatch, "patch key set differs from the frozen set")
    end
    validate_binding(frozen)
    return PatchKeyValidation(length(observed), _sha256(bytes), frozen.keys_sha256)
end

function _patch_receipt_bytes(
    frozen::FrozenUniverse,
    forward::PatchKeyValidation,
    backward::PatchKeyValidation,
)::Vector{UInt8}
    io = IOBuffer()
    for (key, value) in (
        "schema" => "stmfit-structured-patch-key-receipt-v1",
        "schema_version" => 1,
        "feature_sha256" => frozen.feature_sha256,
        "keys_sha256" => frozen.keys_sha256,
        "forward_patch_sha256" => forward.sha256,
        "backward_patch_sha256" => backward.sha256,
        "forward_producer_sha256" => frozen.forward_producer_sha256,
        "backward_producer_sha256" => frozen.backward_producer_sha256,
        "grid_sha256" => frozen.grid_sha256,
        "model_config_sha256" => frozen.model_config_sha256,
        "source_sha256" => frozen.source_sha256,
    )
        value isa Integer ? println(io, key, " = ", value) :
            println(io, key, " = ", _toml_quote(String(value)))
    end
    return take!(io)
end

function _reconcile_patch_validation(
    channel::String,
    initial::PatchKeyValidation,
    current::PatchKeyValidation,
)
    (current.row_count == initial.row_count &&
     current.keys_sha256 == initial.keys_sha256) ||
        _fail(:patch_identity_changed,
              "$channel patch key identity changed during receipt publication")
    current.sha256 == initial.sha256 ||
        _fail(:patch_content_changed,
              "$channel patch bytes changed during receipt publication")
    return nothing
end

function _reconcile_patch_pair(
    frozen::FrozenUniverse,
    forward_path::AbstractString,
    backward_path::AbstractString,
    initial_forward::PatchKeyValidation,
    initial_backward::PatchKeyValidation,
)
    current_forward = validate_patch_keys(frozen, forward_path)
    current_backward = validate_patch_keys(frozen, backward_path)
    _reconcile_patch_validation("forward", initial_forward, current_forward)
    _reconcile_patch_validation("backward", initial_backward, current_backward)
    return nothing
end

function publish_patch_pair_receipt(
    frozen::FrozenUniverse,
    forward_path::AbstractString,
    backward_path::AbstractString,
    receipt_path::AbstractString;
    failpoint::Union{Nothing,Symbol}=nothing,
    _before_patch_reconciliation::Union{Nothing,Function}=nothing,
)::String
    failpoint in (nothing, :before_install) ||
        _fail(:invalid_failpoint, "unsupported patch receipt failpoint")
    forward = validate_patch_keys(frozen, forward_path)
    backward = validate_patch_keys(frozen, backward_path)
    bytes = _patch_receipt_bytes(frozen, forward, backward)
    destination = _resolve_output(frozen.root, receipt_path, "patch receipt")
    if ispath(destination)
        isfile(destination) && !islink(destination) ||
            _fail(:publication_collision, "patch receipt is not a regular file")
        read(destination) == bytes ||
            _fail(:publication_collision, "patch receipt already exists with different bytes")
        _before_patch_reconciliation === nothing || _before_patch_reconciliation()
        validate_binding(frozen)
        _reconcile_patch_pair(
            frozen,
            forward_path,
            backward_path,
            forward,
            backward,
        )
        return _sha256(bytes)
    end
    parent = dirname(destination)
    stage = joinpath(parent, ".$(basename(destination)).stage-$(randstring(20))")
    while ispath(stage)
        stage = joinpath(parent, ".$(basename(destination)).stage-$(randstring(20))")
    end
    installed = false
    try
        _write_synced(stage, bytes)
        failpoint == :before_install &&
            _fail(:interrupted_publication, "simulated interruption before patch receipt install")
        _before_patch_reconciliation === nothing || _before_patch_reconciliation()
        validate_binding(frozen)
        _reconcile_patch_pair(
            frozen,
            forward_path,
            backward_path,
            forward,
            backward,
        )
        ispath(destination) && _fail(:publication_collision, "patch receipt appeared concurrently")
        Base.Filesystem.rename(stage, destination)
        installed = true
        _fsync_directory(parent)
        return _sha256(bytes)
    finally
        !installed && ispath(stage) && rm(stage; force=true)
    end
end

function drop_channel(patch::AbstractMatrix)
    return fill(missing, size(patch))
end

function _bracket(value::Float64)
    clamped = clamp(value, first(GRID_COORDINATES), last(GRID_COORDINATES))
    if clamped == first(GRID_COORDINATES)
        return 1, 1, 0.0
    elseif clamped == last(GRID_COORDINATES)
        return GRID_SIDE, GRID_SIDE, 0.0
    end
    lower = findlast(coordinate -> coordinate <= clamped, GRID_COORDINATES)::Int
    upper = lower + 1
    weight = (clamped - GRID_COORDINATES[lower]) /
             (GRID_COORDINATES[upper] - GRID_COORDINATES[lower])
    return lower, upper, weight
end

function _interpolated_value(patch, u::Float64, t::Float64)
    t_lo, t_hi, t_weight = _bracket(t)
    u_lo, u_hi, u_weight = _bracket(u)
    terms = (
        ((1 - t_weight) * (1 - u_weight), patch[u_lo, t_lo]),
        (t_weight * (1 - u_weight), patch[u_lo, t_hi]),
        ((1 - t_weight) * u_weight, patch[u_hi, t_lo]),
        (t_weight * u_weight, patch[u_hi, t_hi]),
    )
    active = [(weight, value) for (weight, value) in terms if weight > 0.0]
    any(item -> item[2] === missing, active) && return missing
    values = Float64[item[2] for item in active]
    if any(!isfinite, values)
        length(values) == 1 && return only(values)
        all(isequal(first(values)), values) && return first(values)
        return NaN
    end
    return sum(weight * Float64(value) for (weight, value) in active)
end

function shift_patch(patch::AbstractMatrix, delta_t_nm::Real, delta_u_nm::Real)
    size(patch) == (GRID_SIDE, GRID_SIDE) ||
        _fail(:invalid_patch_shape, "shift requires a 17 by 17 patch")
    delta_t = Float64(delta_t_nm)
    delta_u = Float64(delta_u_nm)
    isfinite(delta_t) && isfinite(delta_u) ||
        _fail(:invalid_shift, "shift is nonfinite")
    (delta_t, delta_u) in (
        (-SHIFT_NM, 0.0), (SHIFT_NM, 0.0),
        (0.0, -SHIFT_NM), (0.0, SHIFT_NM),
    ) || _fail(:invalid_shift, "shift is outside the frozen set")
    output = Matrix{Union{Missing,Float64}}(undef, GRID_SIDE, GRID_SIDE)
    for (u_index, u) in enumerate(GRID_COORDINATES),
        (t_index, t) in enumerate(GRID_COORDINATES)
        output[u_index, t_index] = _interpolated_value(
            patch,
            u - delta_u,
            t - delta_t,
        )
    end
    return output
end

function contrast_patch(patch::AbstractMatrix, factor_value::Real)
    factor = Float64(factor_value)
    factor in CONTRAST_FACTORS || _fail(:invalid_contrast, "contrast factor is outside the frozen set")
    output = Matrix{Union{Missing,Float64}}(undef, size(patch))
    for index in eachindex(patch)
        value = patch[index]
        if value === missing
            output[index] = missing
        else
            number = Float64(value)
            output[index] = isfinite(number) ? factor * number : number
        end
    end
    return output
end

function _finite_values(values)
    result = Float64[]
    for value in values
        value === missing && continue
        number = Float64(value)
        isfinite(number) && push!(result, number)
    end
    return result
end

function median_absolute_deviation(values)::Float64
    finite = _finite_values(values)
    isempty(finite) && _fail(:missing_finite_values, "MAD input has no finite values")
    center = median(finite)
    return median(abs.(finite .- center))
end

function correlated_noise(seed::Integer; target_mad::Real)::Matrix{Float64}
    target = Float64(target_mad)
    isfinite(target) && target >= 0.0 ||
        _fail(:invalid_noise_scale, "noise MAD is not finite and nonnegative")
    rng = MersenneTwister(Int(seed))
    initial = randn(rng, GRID_SIDE, GRID_SIDE)
    along_t = similar(initial)
    for u in 1:GRID_SIDE, t in 1:GRID_SIDE
        along_t[u, t] = (
            initial[u, clamp(t - 1, 1, GRID_SIDE)] +
            2initial[u, t] +
            initial[u, clamp(t + 1, 1, GRID_SIDE)]
        ) / 4
    end
    filtered = similar(initial)
    for u in 1:GRID_SIDE, t in 1:GRID_SIDE
        filtered[u, t] = (
            along_t[clamp(u - 1, 1, GRID_SIDE), t] +
            2along_t[u, t] +
            along_t[clamp(u + 1, 1, GRID_SIDE), t]
        ) / 4
    end
    filtered .-= median(vec(filtered))
    current = median_absolute_deviation(vec(filtered))
    current > 0.0 || _fail(:invalid_noise_scale, "filtered noise has zero MAD")
    filtered .*= target / current
    return filtered
end

function _validate_patch_dictionary(frozen::FrozenUniverse, patches)
    supplied = Set{LobeKey}(keys(patches))
    supplied == Set(frozen.keys) ||
        _fail(:patch_key_mismatch, "noise patch keys differ from the frozen set")
    for key in frozen.keys
        size(patches[key]) == (GRID_SIDE, GRID_SIDE) ||
            _fail(:invalid_patch_shape, "noise scale requires 17 by 17 patches")
    end
    return nothing
end

function freeze_noise_scale(
    frozen::FrozenUniverse,
    fold_date::AbstractString,
    channel_value::AbstractString,
    patches,
)::NoiseScale
    validate_binding(frozen)
    fold = String(fold_date)
    fold in (row.date for row in frozen.scans) ||
        _fail(:unknown_fold, "noise scale fold is not frozen")
    channel = String(channel_value)
    isempty(channel) && _fail(:invalid_channel, "noise channel is empty")
    all(isascii, channel) || _fail(:invalid_channel, "noise channel is not ASCII")
    InputBoundary._has_forbidden_concept(channel) &&
        _fail(:invalid_channel, "noise channel is outside the construction boundary")
    _validate_patch_dictionary(frozen, patches)
    training_files = Tuple(row.file for row in frozen.scans if row.date != fold)
    isempty(training_files) && _fail(:unknown_fold, "noise fold has no training scans")
    scan_mads = Float64[]
    for file in training_files
        values = Float64[]
        for key in frozen.keys
            key.file == file || continue
            append!(values, _finite_values(patches[key]))
        end
        isempty(values) && _fail(:missing_finite_values, "training scan has no finite residual pixels")
        push!(scan_mads, median_absolute_deviation(values))
    end
    median_scan_mad = median(scan_mads)
    target = NOISE_MAD_FRACTION * median_scan_mad
    io = IOBuffer()
    println(io, "noise-scale-v1")
    println(io, "fold_date=", fold)
    println(io, "channel=", channel)
    println(io, "feature_sha256=", frozen.feature_sha256)
    println(io, "keys_sha256=", frozen.keys_sha256)
    for (file, value) in zip(training_files, scan_mads)
        println(io, "scan=", file, '\t', _format_number(value))
    end
    println(io, "target_mad=", _format_number(target))
    binding = _sha256(take!(io))
    return NoiseScale(
        fold,
        channel,
        training_files,
        Tuple(scan_mads),
        median_scan_mad,
        target,
        binding,
    )
end

function apply_noise_patch(patch::AbstractMatrix, scale::NoiseScale, seed::Integer)
    size(patch) == (GRID_SIDE, GRID_SIDE) ||
        _fail(:invalid_patch_shape, "noise requires a 17 by 17 patch")
    noise = correlated_noise(seed; target_mad=scale.target_mad)
    output = Matrix{Union{Missing,Float64}}(undef, size(patch))
    for index in eachindex(patch)
        value = patch[index]
        if value === missing
            output[index] = missing
        else
            number = Float64(value)
            output[index] = isfinite(number) ? number + noise[index] : number
        end
    end
    return output
end

end
