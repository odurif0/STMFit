module StructuredEdgeFeatures

using Printf
using SHA
using TOML

include(joinpath(@__DIR__, "universe.jl"))
using .StructuredUniverse

export EdgeFeatureError,
       EdgeBundle,
       CANDIDATE_CONFIG_SHA256,
       MODEL_CONFIG_SHA256,
       FORWARD_PRODUCER_SHA256,
       BACKWARD_PRODUCER_SHA256,
       GRID_SHA256,
       T2_REVIEW_SHA256,
       T3_REVIEW_SHA256,
       FORWARD_PATCH_HEADER,
       BACKWARD_PATCH_HEADER,
       EDGE_HEADER,
       NODE_HEADER,
       grid_contract_bytes,
       expected_producer_commands,
       build_edge_bundle,
       publish_edge_bundle,
       validate_edge_bundle

const TEST_DIRECTORY = normpath(joinpath(@__DIR__, "..", ".."))
const SOURCE_REPOSITORY_ROOT = normpath(joinpath(TEST_DIRECTORY, ".."))
const FIREWALL_SOURCE_PATH = joinpath(@__DIR__, "firewall.jl")
const UNIVERSE_SOURCE_PATH = joinpath(@__DIR__, "universe.jl")
const UNIVERSE_CLI_PATH = joinpath(TEST_DIRECTORY, "build_structured_scan_universe.jl")
const EDGE_SOURCE_PATH = abspath(@__FILE__)
const EDGE_CLI_PATH = joinpath(TEST_DIRECTORY, "build_label_free_edge_features.jl")

const T2_REVIEW_PATH = joinpath(
    SOURCE_REPOSITORY_ROOT,
    ".omo", "evidence", "structured-label-free-unit-assignment",
    "t2", "correction2", "review", "AdversarialVerify.json",
)
const T3_REVIEW_PATH = joinpath(
    SOURCE_REPOSITORY_ROOT,
    ".omo", "evidence", "structured-label-free-unit-assignment",
    "provenance-rebind", "phase3-t3-universe", "correction", "review", "AdversarialVerify.json",
)

const CANDIDATE_CONFIG_SHA256 =
    "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const MODEL_CONFIG_SHA256 =
    "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const FORWARD_PRODUCER_SHA256 =
    "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad"
const BACKWARD_PRODUCER_SHA256 =
    "4f866a1e27289a0cf010a6ed4b5e0b92454dc55f97807d33cc95edb9e87182ba"
const GRID_SHA256 =
    "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"
const T2_REVIEW_SHA256 =
    "5555459b99f2a848a53789b2e8e15eb683760ffc0d7056b9e38275f0d65b6aed"
const T3_REVIEW_SHA256 =
    "531d68e866d896809cdb0c1fde2ac8b934e06a2562ca1fecbb65d84bca5e9954"

const GRID_SIDE = 17
const PIXELS = 289
const NORM_FLOOR = 1.0e-12
const DOT_TOLERANCE = 1.0e-12
const SPACING_MIN_NM = 0.35
const SPACING_MAX_NM = 0.75

const PATCH_METADATA_HEADER = ["file", "lobe", "t_nm", "u_nm", "amplitude"]
const FORWARD_RAW_COLUMNS = [@sprintf("raw_p%03d", index) for index in 1:PIXELS]
const FORWARD_RESIDUAL_COLUMNS = [@sprintf("res_p%03d", index) for index in 1:PIXELS]
const BACKWARD_RAW_COLUMNS = [@sprintf("bwd_raw_p%03d", index) for index in 1:PIXELS]
const BACKWARD_RESIDUAL_COLUMNS = [@sprintf("bwd_res_p%03d", index) for index in 1:PIXELS]
const DIFFERENCE_RAW_COLUMNS = [@sprintf("diff_raw_p%03d", index) for index in 1:PIXELS]
const DIFFERENCE_RESIDUAL_COLUMNS = [@sprintf("diff_res_p%03d", index) for index in 1:PIXELS]

# The wide headers are validated as provenance from the frozen producers. Only
# FORWARD_RESIDUAL_COLUMNS and BACKWARD_RESIDUAL_COLUMNS are ever parsed into
# observations; the other payload blocks cannot be selected by any API.
const FORWARD_PATCH_HEADER = vcat(
    PATCH_METADATA_HEADER,
    FORWARD_RAW_COLUMNS,
    FORWARD_RESIDUAL_COLUMNS,
)
const BACKWARD_PATCH_HEADER = vcat(
    PATCH_METADATA_HEADER,
    BACKWARD_RAW_COLUMNS,
    BACKWARD_RESIDUAL_COLUMNS,
    DIFFERENCE_RAW_COLUMNS,
    DIFFERENCE_RESIDUAL_COLUMNS,
)

const EDGE_HEADER = [
    "file", "left_lobe", "right_lobe", "left_t_nm", "right_t_nm", "gap_nm",
    "left_segment_id", "right_segment_id", "edge_status", "split_reason",
    "corr_fwd", "corr_bwd", "grid_hash", "forward_patch_hash",
    "backward_patch_hash", "feature_hash", "config_hash", "source_hash",
]
const NODE_HEADER = [
    "file", "lobe", "t_nm", "segment_id", "node_status",
    "left_boundary_reason", "right_boundary_reason", "feature_hash",
    "config_hash", "source_hash",
]
const OUTPUT_NAMES = (
    "edge_observations.tsv",
    "node_segments.tsv",
    "receipt.toml",
)
const _PUBLICATION_INSTALL_HOOK = Ref{Union{Nothing,Function}}(nothing)
const _AT_FDCWD = Cint(-100)
const _RENAME_NOREPLACE = Cuint(1)
const _LINUX_EEXIST = 17
const _LINUX_ENOTEMPTY = 39
const _LINUX_EINVAL = 22
const _LINUX_ENOSYS = 38
const _LINUX_EOPNOTSUPP = 95
const PATCH_RECEIPT_KEYS = Set([
    "schema", "schema_version", "status", "channel", "execution_site",
    "julia_threads", "root", "cwd", "features_path", "data_path",
    "output_path", "source_path", "source_sha256", "feature_sha256",
    "output_sha256", "candidate_config_sha256", "model_config_sha256",
    "grid_sha256", "keys_sha256", "table_header_sha256", "row_count",
    "command",
])

struct EdgeFeatureError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, error::EdgeFeatureError) =
    print(io, "BLOCKED [", error.code, "]: ", error.message)

_fail(code::Symbol, message::AbstractString) =
    throw(EdgeFeatureError(code, String(message)))

struct FileSnapshot
    path::String
    sha256::String
    byte_count::Int
    device::UInt64
    inode::UInt64
    mtime::Float64
    ctime::Float64
end

struct NodeInput
    file::String
    lobe::Int
    t_nm::Float64
    t_text::String
    u_text::String
    amplitude_text::String
end

struct PatchVector
    status::Symbol
    normalized::Vector{Float64}
end

struct PatchTable
    snapshot::FileSnapshot
    header::Vector{String}
    vectors::Dict{Tuple{String,Int},Vector{Float64}}
end

struct EdgeDraft
    file::String
    left::NodeInput
    right::NodeInput
    gap_nm::Float64
    status::Symbol
    reason::Symbol
    corr_forward::Union{Nothing,Float64}
    corr_backward::Union{Nothing,Float64}
end

struct EdgeBundle
    files::Dict{String,Vector{UInt8}}
    frozen::StructuredUniverse.FrozenUniverse
    snapshots::Tuple{Vararg{FileSnapshot}}
end

_sha256(bytes) = bytes2hex(sha256(bytes))
_toml_quote(value::String) = repr(value)

function grid_contract_bytes()::Vector{UInt8}
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

function _canonical_existing_file(path::AbstractString, context::String)::String
    supplied = String(path)
    isempty(supplied) && _fail(:noncanonical_path, "$context path is empty")
    isabspath(supplied) || _fail(:noncanonical_path, "$context path is not absolute")
    absolute = normpath(abspath(supplied))
    isfile(absolute) || _fail(:missing_input, "$context is not a regular file")
    islink(absolute) && _fail(:symlink_rejected, "$context is a symlink")
    realpath(absolute) == absolute ||
        _fail(:noncanonical_path, "$context path is not canonical")
    return absolute
end

function _canonical_existing_directory(path::AbstractString, context::String)::String
    supplied = String(path)
    isempty(supplied) && _fail(:noncanonical_path, "$context path is empty")
    isabspath(supplied) || _fail(:noncanonical_path, "$context path is not absolute")
    absolute = normpath(abspath(supplied))
    isdir(absolute) || _fail(:missing_input, "$context is not a directory")
    islink(absolute) && _fail(:symlink_rejected, "$context is a symlink")
    realpath(absolute) == absolute ||
        _fail(:noncanonical_path, "$context path is not canonical")
    return absolute
end

function _resolve_existing(root::String, supplied::AbstractString, context::String)::String
    text = String(supplied)
    absolute = isabspath(text) ? normpath(text) : normpath(joinpath(root, text))
    relative = relpath(absolute, root)
    parts = splitpath(relative)
    (relative == "." || (!isempty(parts) && first(parts) == "..")) &&
        _fail(:path_escape, "$context escapes root")
    return _canonical_existing_file(absolute, context)
end

function _resolve_existing_directory(
    root::String,
    supplied::AbstractString,
    context::String,
)::String
    text = String(supplied)
    absolute = isabspath(text) ? normpath(text) : normpath(joinpath(root, text))
    relative = relpath(absolute, root)
    parts = splitpath(relative)
    (relative == "." || (!isempty(parts) && first(parts) == "..")) &&
        _fail(:path_escape, "$context escapes root")
    return _canonical_existing_directory(absolute, context)
end

function _resolve_output(root::String, supplied::AbstractString, context::String)::String
    text = String(supplied)
    isempty(text) && _fail(:noncanonical_path, "$context path is empty")
    absolute = isabspath(text) ? normpath(text) : normpath(joinpath(root, text))
    relative = relpath(absolute, root)
    parts = splitpath(relative)
    (relative == "." || (!isempty(parts) && first(parts) == "..")) &&
        _fail(:path_escape, "$context escapes root")
    parent = dirname(absolute)
    isdir(parent) || _fail(:missing_output_parent, "$context parent does not exist")
    islink(parent) && _fail(:symlink_rejected, "$context parent is a symlink")
    realpath(parent) == normpath(abspath(parent)) ||
        _fail(:noncanonical_path, "$context parent is not canonical")
    islink(absolute) && _fail(:symlink_rejected, "$context is a symlink")
    return absolute
end

function _snapshot(path::String, context::String)::Tuple{Vector{UInt8},FileSnapshot}
    absolute = _canonical_existing_file(path, context)
    before = stat(absolute)
    bytes = read(absolute)
    after = stat(absolute)
    stable = before.device == after.device && before.inode == after.inode &&
             before.size == after.size && before.mtime == after.mtime &&
             before.ctime == after.ctime && length(bytes) == after.size
    stable || _fail(:dependency_changed, "$context changed while read")
    snapshot = FileSnapshot(
        absolute,
        _sha256(bytes),
        length(bytes),
        UInt64(after.device),
        UInt64(after.inode),
        Float64(after.mtime),
        Float64(after.ctime),
    )
    return bytes, snapshot
end

function _verify_snapshot(snapshot::FileSnapshot)
    _, current = _snapshot(snapshot.path, basename(snapshot.path))
    current.sha256 == snapshot.sha256 &&
    current.byte_count == snapshot.byte_count &&
    current.device == snapshot.device &&
    current.inode == snapshot.inode &&
    current.mtime == snapshot.mtime &&
    current.ctime == snapshot.ctime ||
        _fail(:dependency_changed, "$(snapshot.path) changed after binding")
    return nothing
end

function _validate_hash(value, context::String)::String
    value isa String || _fail(:invalid_hash, "$context is not a string")
    occursin(r"^[0-9a-f]{64}$", value) ||
        _fail(:invalid_hash, "$context is not a lowercase SHA-256")
    return value
end

function _parse_toml(bytes::Vector{UInt8}, context::String)::Dict{String,Any}
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "$context is not valid UTF-8")
    startswith(text, '\ufeff') && _fail(:invalid_encoding, "$context has a byte-order mark")
    occursin('\r', text) && _fail(:carriage_return, "$context must use LF line endings")
    endswith(text, '\n') || _fail(:missing_final_lf, "$context must end with LF")
    parsed = try
        TOML.parse(text)
    catch error
        _fail(:malformed_toml, "$context cannot be parsed: $(sprint(showerror, error))")
    end
    return Dict{String,Any}(String(key) => value for (key, value) in parsed)
end

function _exact_keys(table::AbstractDict, expected::Set{String}, context::String)
    observed = Set(String(key) for key in keys(table))
    observed == expected || _fail(:patch_receipt_schema,
        "$context keys differ from the frozen schema")
    return nothing
end

function _required_string(table::AbstractDict, key::String, context::String)::String
    value = get(table, key, nothing)
    value isa String || _fail(:patch_receipt_schema, "$context.$key is not a string")
    return String(value)
end

function _required_integer(table::AbstractDict, key::String, context::String)::Int
    value = get(table, key, nothing)
    (value isa Integer && !(value isa Bool)) ||
        _fail(:patch_receipt_schema, "$context.$key is not an integer")
    return Int(value)
end

function _parse_table(
    bytes::Vector{UInt8},
    context::String;
    validate_width::Bool=true,
)
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "$context is not valid UTF-8")
    startswith(text, '\ufeff') && _fail(:invalid_encoding, "$context has a byte-order mark")
    occursin('\r', text) && _fail(:carriage_return, "$context must use LF line endings")
    endswith(text, '\n') || _fail(:missing_final_lf, "$context must end with LF")
    lines = split(text, '\n'; keepempty=true)
    pop!(lines)
    isempty(lines) && _fail(:empty_input, "$context is empty")
    any(isempty, lines) && _fail(:empty_row, "$context contains an empty row")
    header = String.(split(first(lines), '\t'; keepempty=true))
    any(isempty, header) && _fail(:empty_column, "$context contains an empty column")
    length(header) == length(unique(header)) ||
        _fail(:duplicate_column, "$context contains a duplicate column")
    rows = Vector{Vector{String}}()
    for (row_index, line) in enumerate(lines[2:end])
        fields = String.(split(line, '\t'; keepempty=true))
        !validate_width || length(fields) == length(header) ||
            _fail(:row_width_mismatch, "$context row $row_index has the wrong width")
        push!(rows, fields)
    end
    isempty(rows) && _fail(:empty_input, "$context has no data rows")
    return header, rows
end

function _canonical_positive_integer(text::String, context::String)::Int
    occursin(r"^[1-9][0-9]*$", text) ||
        _fail(:invalid_key, "$context is not a canonical positive integer")
    value = tryparse(Int, text)
    value === nothing && _fail(:invalid_key, "$context is outside the integer range")
    return value
end

function _finite_number(text::String, context::String)::Float64
    number = tryparse(Float64, text)
    number === nothing && _fail(:invalid_number, "$context is not numerical")
    isfinite(number) || _fail(:invalid_number, "$context is nonfinite")
    return number
end

function _patch_number(text::String, context::String)::Float64
    text == "NA" && return NaN
    number = tryparse(Float64, text)
    number === nothing && _fail(:invalid_patch_value, "$context is not numerical or NA")
    return number
end

function _validate_fixed_review(path::String, expected::String, context::String)::FileSnapshot
    bytes, snapshot = _snapshot(path, context)
    snapshot.sha256 == expected ||
        _fail(:dependency_hash_mismatch, "$context SHA-256 differs from the confirmed review")
    text = String(copy(bytes))
    isvalid(text) || _fail(:invalid_encoding, "$context bytes are not valid UTF-8")
    occursin("\"verdict\": \"confirmed\"", text) ||
        _fail(:dependency_verdict_mismatch, "$context is not a confirmed review")
    return snapshot
end

function _table_at(document::AbstractDict, keys::Vector{String}, context::String)
    current = document
    for key in keys
        value = get(current, key, nothing)
        value isa AbstractDict || _fail(:config_contract_mismatch,
            "$context.$key is not a table")
        current = value
    end
    return current
end

function _validate_model_contract(model_bytes::Vector{UInt8})
    model = _parse_toml(model_bytes, "model config")
    preprocessing = _table_at(model, ["preprocessing"], "model")
    get(preprocessing, "spacing_min_nm", nothing) == SPACING_MIN_NM ||
        _fail(:config_contract_mismatch, "spacing_min_nm differs from the frozen value")
    get(preprocessing, "spacing_max_nm", nothing) == SPACING_MAX_NM ||
        _fail(:config_contract_mismatch, "spacing_max_nm differs from the frozen value")
    get(preprocessing, "gap_policy", nothing) == "consecutive_observed_t_nm_only" ||
        _fail(:config_contract_mismatch, "gap policy differs from the frozen value")

    patch = _table_at(model, ["preprocessing", "patch"], "model")
    required_patch = Dict{String,Any}(
        "half_t_nm" => 0.32,
        "half_u_nm" => 0.32,
        "step_nm" => 0.04,
        "pixels_per_axis" => 17,
        "pixels" => 289,
        "order" => "u_outer,t_inner",
        "forward_prefix" => "res_p",
        "backward_prefix" => "bwd_res_p",
        "grid_sha256" => GRID_SHA256,
        "centering" => "per_lobe_channel_arithmetic_mean",
        "scaling" => "per_lobe_channel_unit_l2",
        "l2_norm_floor" => NORM_FLOOR,
        "correlation_roundoff_tolerance" => DOT_TOLERANCE,
        "serialization_format" => "%.17g",
        "reversal_index_formula" => "R(k)=17*(u_idx-1)+(18-t_idx)",
    )
    for (key, expected) in required_patch
        get(patch, key, nothing) == expected ||
            _fail(:config_contract_mismatch, "model.preprocessing.patch.$key differs")
    end
    get(patch, "coordinates_nm", nothing) == collect(StructuredUniverse.GRID_COORDINATES) ||
        _fail(:config_contract_mismatch, "patch coordinates differ")

    producers = _table_at(model, ["preprocessing", "patch_producers"], "model")
    required_producers = Dict{String,Any}(
        "forward_source" => "test/extract_lobe_patches.jl",
        "forward_source_sha256" => FORWARD_PRODUCER_SHA256,
        "forward_command" => "julia --project=\$ROOT \$ROOT/test/extract_lobe_patches.jl --features \$FEATURES --data-dir \$DATA --half-nm 0.32 --half-u-nm 0.32 --step-nm 0.04 --out \$FWD",
        "backward_source" => "test/extract_lobe_patches_bwd.jl",
        "backward_source_sha256" => BACKWARD_PRODUCER_SHA256,
        "backward_command" => "julia --project=\$ROOT \$ROOT/test/extract_lobe_patches_bwd.jl --features \$FEATURES --data-dir \$DATA --half-nm 0.32 --step-nm 0.04 --out \$BWD",
        "execution_site" => "Viper_Slurm_compute_node",
        "julia_threads" => 1,
    )
    for (key, expected) in required_producers
        get(producers, key, nothing) == expected ||
            _fail(:config_contract_mismatch, "model.preprocessing.patch_producers.$key differs")
    end
    return nothing
end

function _validate_candidate_contract(candidate_bytes::Vector{UInt8})
    candidate = _parse_toml(candidate_bytes, "candidate config")
    provenance = _table_at(candidate, ["provenance"], "candidate")
    required = Dict(
        "model_config_sha256" => MODEL_CONFIG_SHA256,
        "spacing_source" => "config/chitosan.toml",
        "forward_patch_producer" => "test/extract_lobe_patches.jl",
        "forward_patch_producer_sha256" => FORWARD_PRODUCER_SHA256,
        "backward_patch_producer" => "test/extract_lobe_patches_bwd.jl",
        "backward_patch_producer_sha256" => BACKWARD_PRODUCER_SHA256,
    )
    for (key, expected) in required
        get(provenance, key, nothing) == expected ||
            _fail(:config_contract_mismatch, "candidate.provenance.$key differs")
    end
    return nothing
end

function _source_snapshots()::Tuple{String,Vector{FileSnapshot}}
    entries = [
        "test/build_label_free_edge_features.jl" => EDGE_CLI_PATH,
        "test/lib/structured_assignment/edge_features.jl" => EDGE_SOURCE_PATH,
        "test/lib/structured_assignment/firewall.jl" => FIREWALL_SOURCE_PATH,
        "test/lib/structured_assignment/universe.jl" => UNIVERSE_SOURCE_PATH,
    ]
    snapshots = FileSnapshot[]
    io = IOBuffer()
    println(io, "structured-edge-source-v1")
    for (relative, path) in sort(entries; by=first)
        _, snapshot = _snapshot(path, relative)
        push!(snapshots, snapshot)
        println(io, relative, '=', snapshot.sha256)
    end
    return _sha256(take!(io)), snapshots
end

function _command_path(path::String, context::String)::String
    any(isspace, path) && _fail(:ambiguous_command_path,
        "$context contains whitespace and cannot match the frozen command")
    occursin('\0', path) && _fail(:ambiguous_command_path, "$context contains NUL")
    return path
end

function expected_producer_commands(
    root::AbstractString;
    features::AbstractString,
    data_dir::AbstractString,
    forward_patches::AbstractString,
    backward_patches::AbstractString,
)
    canonical_root = _command_path(_canonical_root(root), "root")
    feature_path = _command_path(
        _canonical_existing_file(abspath(String(features)), "feature table"),
        "feature table",
    )
    data_path = _command_path(
        _canonical_existing_directory(abspath(String(data_dir)), "data directory"),
        "data directory",
    )
    forward_path = _command_path(
        _canonical_existing_file(abspath(String(forward_patches)), "forward patch table"),
        "forward patch table",
    )
    backward_path = _command_path(
        _canonical_existing_file(abspath(String(backward_patches)), "backward patch table"),
        "backward patch table",
    )
    forward_source = _command_path(
        joinpath(canonical_root, "test", "extract_lobe_patches.jl"),
        "forward producer",
    )
    backward_source = _command_path(
        joinpath(canonical_root, "test", "extract_lobe_patches_bwd.jl"),
        "backward producer",
    )
    forward = "julia --project=$canonical_root $forward_source --features $feature_path " *
              "--data-dir $data_path --half-nm 0.32 --half-u-nm 0.32 " *
              "--step-nm 0.04 --out $forward_path"
    backward = "julia --project=$canonical_root $backward_source --features $feature_path " *
               "--data-dir $data_path --half-nm 0.32 --step-nm 0.04 " *
               "--out $backward_path"
    return (forward=forward, backward=backward)
end

function _parse_feature_nodes(
    bytes::Vector{UInt8},
    frozen::StructuredUniverse.FrozenUniverse,
)::Dict{Tuple{String,Int},NodeInput}
    header, rows = _parse_table(bytes, "feature table")
    required = ["file", "N", "lobe", "amplitude", "t_nm", "u_nm"]
    all(column -> column in header, required) ||
        _fail(:feature_schema_mismatch, "feature table lacks a required selected-lobe column")
    indices = Dict(column => findfirst(==(column), header)::Int for column in required)
    nodes = Dict{Tuple{String,Int},NodeInput}()
    for (row_index, fields) in enumerate(rows)
        file = fields[indices["file"]]
        StructuredUniverse.parse_scan_date(file)
        lobe = _canonical_positive_integer(fields[indices["lobe"]],
                                             "feature row $row_index lobe")
        key = (file, lobe)
        haskey(nodes, key) && _fail(:duplicate_feature_key, "feature table repeats a key")
        t_text = fields[indices["t_nm"]]
        u_text = fields[indices["u_nm"]]
        amplitude_text = fields[indices["amplitude"]]
        node = NodeInput(
            file,
            lobe,
            _finite_number(t_text, "feature row $row_index t_nm"),
            t_text,
            u_text,
            amplitude_text,
        )
        _finite_number(u_text, "feature row $row_index u_nm")
        _finite_number(amplitude_text, "feature row $row_index amplitude")
        nodes[key] = node
    end
    expected = Set((key.file, key.lobe) for key in frozen.keys)
    Set(keys(nodes)) == expected ||
        _fail(:feature_key_mismatch, "feature keys differ from the frozen universe")
    return nodes
end

function _residual_columns(header::Vector{String}, prefix::String)::Vector{String}
    return [column for column in header if startswith(column, prefix)]
end

function _validate_patch_header(channel::Symbol, header::Vector{String})
    expected_header = channel == :forward ? FORWARD_PATCH_HEADER : BACKWARD_PATCH_HEADER
    expected_residuals = channel == :forward ?
        FORWARD_RESIDUAL_COLUMNS : BACKWARD_RESIDUAL_COLUMNS
    prefix = channel == :forward ? "res_p" : "bwd_res_p"
    observed = _residual_columns(header, prefix)
    missing = setdiff(expected_residuals, observed)
    isempty(missing) || _fail(:missing_residual_column,
        "$(String(channel)) table is missing $(first(missing))")
    extra = setdiff(observed, expected_residuals)
    isempty(extra) || _fail(:extra_residual_column,
        "$(String(channel)) table has an extra residual column")
    observed == expected_residuals || _fail(:reordered_residual_columns,
        "$(String(channel)) residual columns are not in 001 through 289 order")
    header == expected_header || _fail(:patch_schema_mismatch,
        "$(String(channel)) patch table does not match the frozen producer schema")
    return nothing
end

function _parse_patch_table(
    path::String,
    channel::Symbol,
    nodes::Dict{Tuple{String,Int},NodeInput},
    frozen::StructuredUniverse.FrozenUniverse,
)::PatchTable
    bytes, snapshot = _snapshot(path, "$(String(channel)) patch table")
    header, rows = _parse_table(
        bytes,
        "$(String(channel)) patch table";
        validate_width=false,
    )
    _validate_patch_header(channel, header)
    residual_names = channel == :forward ?
        FORWARD_RESIDUAL_COLUMNS : BACKWARD_RESIDUAL_COLUMNS
    residual_indices = [findfirst(==(name), header)::Int for name in residual_names]
    vectors = Dict{Tuple{String,Int},Vector{Float64}}()
    observed_keys = Tuple{String,Int}[]
    for (row_index, fields) in enumerate(rows)
        length(fields) == length(header) || _fail(:row_width_mismatch,
            "$(String(channel)) patch row $row_index has the wrong width")
        file = fields[1]
        StructuredUniverse.parse_scan_date(file)
        lobe = _canonical_positive_integer(fields[2],
                                             "$(String(channel)) row $row_index lobe")
        key = (file, lobe)
        haskey(vectors, key) && _fail(:duplicate_patch_key,
            "$(String(channel)) patch table repeats a key")
        push!(observed_keys, key)
        node = get(nodes, key, nothing)
        node === nothing && continue
        fields[3] == node.t_text && fields[4] == node.u_text &&
        fields[5] == node.amplitude_text ||
            _fail(:patch_metadata_mismatch,
                "$(String(channel)) patch metadata differs from the feature row")
        for column_index in 6:length(fields)
            value = fields[column_index]
            value == "NA" && continue
            tryparse(Float64, value) === nothing &&
                _fail(:invalid_patch_value,
                    "$(String(channel)) row $row_index column $column_index is malformed")
        end
        vectors[key] = [
            _patch_number(fields[index],
                "$(String(channel)) row $row_index $(header[index])")
            for index in residual_indices
        ]
    end
    expected_keys = [(key.file, key.lobe) for key in frozen.keys]
    if length(observed_keys) != length(unique(observed_keys))
        _fail(:duplicate_patch_key, "$(String(channel)) patch table repeats a key")
    elseif Set(observed_keys) != Set(expected_keys)
        _fail(:patch_key_mismatch,
            "$(String(channel)) patch keys differ from the frozen universe")
    elseif observed_keys != expected_keys
        _fail(:patch_key_order,
            "$(String(channel)) patch keys are not in frozen canonical order")
    end
    length(vectors) == length(expected_keys) ||
        _fail(:patch_key_mismatch,
            "$(String(channel)) patch table lacks a frozen feature key")
    return PatchTable(snapshot, header, vectors)
end

function _validate_patch_receipt(
    path::String,
    channel::Symbol,
    root::String,
    features_path::String,
    patch::PatchTable,
    peer_patch_path::String,
    frozen::StructuredUniverse.FrozenUniverse,
)
    bytes, snapshot = _snapshot(path, "$(String(channel)) producer receipt")
    receipt = _parse_toml(bytes, "$(String(channel)) producer receipt")
    _exact_keys(receipt, PATCH_RECEIPT_KEYS, "$(String(channel)) producer receipt")
    _required_string(receipt, "schema", "producer receipt") ==
        "stmfit-structured-patch-producer-receipt-v1" ||
        _fail(:patch_receipt_schema, "producer receipt schema differs")
    _required_integer(receipt, "schema_version", "producer receipt") == 1 ||
        _fail(:patch_receipt_schema, "producer receipt version differs")
    _required_string(receipt, "status", "producer receipt") == "PASS" ||
        _fail(:patch_receipt_status, "producer receipt is not terminal PASS")
    _required_string(receipt, "channel", "producer receipt") == String(channel) ||
        _fail(:patch_receipt_schema, "producer receipt channel differs")
    _required_string(receipt, "execution_site", "producer receipt") ==
        "Viper_Slurm_compute_node" ||
        _fail(:producer_site_mismatch, "producer receipt execution site differs")
    _required_integer(receipt, "julia_threads", "producer receipt") == 1 ||
        _fail(:producer_threads_mismatch, "producer receipt thread count differs")

    receipt_root = _canonical_existing_directory(
        _required_string(receipt, "root", "producer receipt"), "receipt root")
    receipt_root == root || _fail(:producer_path_mismatch, "receipt root differs")
    receipt_cwd = _canonical_existing_directory(
        _required_string(receipt, "cwd", "producer receipt"), "receipt cwd")
    receipt_cwd == root || _fail(:producer_path_mismatch, "receipt cwd differs")
    receipt_features = _canonical_existing_file(
        _required_string(receipt, "features_path", "producer receipt"),
        "receipt feature table",
    )
    receipt_features == features_path ||
        _fail(:producer_path_mismatch, "receipt feature path differs")
    data_path = _canonical_existing_directory(
        _required_string(receipt, "data_path", "producer receipt"),
        "receipt data directory",
    )
    receipt_output = _canonical_existing_file(
        _required_string(receipt, "output_path", "producer receipt"),
        "receipt patch output",
    )
    receipt_output == patch.snapshot.path ||
        _fail(:producer_path_mismatch, "receipt output path differs")

    producer_relative = channel == :forward ?
        joinpath("test", "extract_lobe_patches.jl") :
        joinpath("test", "extract_lobe_patches_bwd.jl")
    producer_path = _canonical_existing_file(joinpath(root, producer_relative),
                                               "$(String(channel)) producer")
    receipt_source = _canonical_existing_file(
        _required_string(receipt, "source_path", "producer receipt"),
        "receipt producer source",
    )
    receipt_source == producer_path ||
        _fail(:producer_path_mismatch, "receipt source path differs")
    expected_source_hash = channel == :forward ?
        FORWARD_PRODUCER_SHA256 : BACKWARD_PRODUCER_SHA256
    _validate_hash(_required_string(receipt, "source_sha256", "producer receipt"),
                   "producer source hash") == expected_source_hash ||
        _fail(:producer_source_hash_mismatch, "receipt producer SHA-256 differs")
    patch.snapshot.sha256 == _validate_hash(
        _required_string(receipt, "output_sha256", "producer receipt"),
        "patch output hash",
    ) || _fail(:patch_hash_mismatch, "patch bytes differ from the producer receipt")
    frozen.feature_sha256 == _validate_hash(
        _required_string(receipt, "feature_sha256", "producer receipt"),
        "feature hash",
    ) || _fail(:feature_hash_mismatch, "producer receipt feature hash differs")
    _required_string(receipt, "candidate_config_sha256", "producer receipt") ==
        CANDIDATE_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "producer receipt candidate hash differs")
    _required_string(receipt, "model_config_sha256", "producer receipt") ==
        MODEL_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "producer receipt model hash differs")
    _required_string(receipt, "grid_sha256", "producer receipt") == GRID_SHA256 ||
        _fail(:grid_hash_mismatch, "producer receipt grid hash differs")
    _required_string(receipt, "keys_sha256", "producer receipt") == frozen.keys_sha256 ||
        _fail(:patch_key_mismatch, "producer receipt key hash differs")
    header_hash = _sha256(codeunits(join(patch.header, '\t') * "\n"))
    _required_string(receipt, "table_header_sha256", "producer receipt") == header_hash ||
        _fail(:patch_schema_hash_mismatch, "producer receipt table-header hash differs")
    _required_integer(receipt, "row_count", "producer receipt") == length(frozen.keys) ||
        _fail(:patch_key_mismatch, "producer receipt row count differs")

    commands = expected_producer_commands(
        root;
        features=features_path,
        data_dir=data_path,
        forward_patches=channel == :forward ? patch.snapshot.path : peer_patch_path,
        backward_patches=channel == :backward ? patch.snapshot.path : peer_patch_path,
    )
    expected_command = channel == :forward ? commands.forward : commands.backward
    _required_string(receipt, "command", "producer receipt") == expected_command ||
        _fail(:producer_command_mismatch, "producer command differs from the frozen invocation")
    return (snapshot=snapshot, command=expected_command, data_path=data_path)
end

function _normalize_patch(values::Vector{Float64})::PatchVector
    length(values) == PIXELS || _fail(:internal_contract, "patch vector does not have 289 values")
    all(isfinite, values) || return PatchVector(:nonfinite, Float64[])
    total = 0.0
    for value in values
        total += value
    end
    mean_value = total / PIXELS
    isfinite(mean_value) || return PatchVector(:nonfinite, Float64[])
    centered = Vector{Float64}(undef, PIXELS)
    squared_norm = 0.0
    for index in eachindex(values)
        value = values[index] - mean_value
        centered[index] = value
        squared_norm += value * value
    end
    norm_value = sqrt(squared_norm)
    isfinite(norm_value) || return PatchVector(:nonfinite, Float64[])
    norm_value <= NORM_FLOOR && return PatchVector(:zero_variance, Float64[])
    normalized = Vector{Float64}(undef, PIXELS)
    for index in eachindex(centered)
        normalized[index] = centered[index] / norm_value
    end
    all(isfinite, normalized) || return PatchVector(:nonfinite, Float64[])
    return PatchVector(:valid, normalized)
end

function _edge_patch_status(left::Union{Nothing,PatchVector}, right::Union{Nothing,PatchVector})
    (left === nothing || right === nothing) && return :missing
    (left.status == :nonfinite || right.status == :nonfinite) && return :nonfinite
    (left.status == :zero_variance || right.status == :zero_variance) &&
        return :zero_variance
    left.status == :valid && right.status == :valid ||
        _fail(:internal_contract, "unknown patch status")
    return :valid
end

function _guard_dot(value::Float64)
    (!isfinite(value) || value < -1.0 - DOT_TOLERANCE ||
     value > 1.0 + DOT_TOLERANCE) &&
        return (status=:out_of_range, value=value)
    return (status=:valid, value=clamp(value, -1.0, 1.0))
end

function _centered_dot(left::PatchVector, right::PatchVector)
    left.status == :valid && right.status == :valid ||
        return (status=:unavailable, value=NaN)
    value = 0.0
    for index in 1:PIXELS
        value += left.normalized[index] * right.normalized[index]
    end
    return _guard_dot(value)
end

function _select_split_reason(
    gap_ok::Bool,
    forward_status::Symbol,
    backward_status::Symbol,
    forward_correlation_status::Symbol,
    backward_correlation_status::Symbol,
)::Symbol
    !gap_ok && return :gap_out_of_range
    forward_status == :missing && return :missing_forward_key
    backward_status == :missing && return :missing_backward_key
    forward_status == :nonfinite && return :nonfinite_fwd
    backward_status == :nonfinite && return :nonfinite_bwd
    forward_status == :zero_variance && return :zero_variance_fwd
    backward_status == :zero_variance && return :zero_variance_bwd
    forward_correlation_status == :out_of_range && return :corr_out_of_range_fwd
    backward_correlation_status == :out_of_range && return :corr_out_of_range_bwd
    return :none
end

function _select_split_reason(
    gap_ok::Bool,
    forward_status::Symbol,
    backward_status::Symbol,
    forward_correlation_status::Symbol,
    backward_correlation_status::Symbol,
    ::Symbol,
    ::Symbol,
)::Symbol
    return _select_split_reason(
        gap_ok,
        forward_status,
        backward_status,
        forward_correlation_status,
        backward_correlation_status,
    )
end

function _reverse_patch_values(values::AbstractVector{<:Real})::Vector{Float64}
    length(values) == PIXELS || _fail(:invalid_patch_shape,
        "reversal requires exactly 289 values")
    output = Vector{Float64}(undef, PIXELS)
    for u_index in 1:GRID_SIDE, t_index in 1:GRID_SIDE
        index = GRID_SIDE * (u_index - 1) + t_index
        reversed_index = GRID_SIDE * (u_index - 1) + (GRID_SIDE + 1 - t_index)
        output[index] = Float64(values[reversed_index])
    end
    return output
end

function _format_number(value::Float64)::String
    isfinite(value) || _fail(:internal_contract, "attempted to serialize a nonfinite number")
    return @sprintf("%.17g", value)
end

function _segment_id(file::String, lobes)::String
    sequence = Int.(collect(lobes))
    isempty(sequence) && _fail(:internal_contract, "cannot identify an empty segment")
    first_lobe = min(first(sequence), last(sequence))
    last_lobe = max(first(sequence), last(sequence))
    bytes = codeunits(
        "segment-v1\n" *
        "file=$(basename(file))\n" *
        "first_lobe=$first_lobe\n" *
        "last_lobe=$last_lobe\n" *
        "grid_hash=$GRID_SHA256\n",
    )
    return _sha256(bytes)[1:16]
end

function _build_topology(
    nodes::Dict{Tuple{String,Int},NodeInput},
    forward_vectors::Dict{Tuple{String,Int},Vector{Float64}},
    backward_vectors::Dict{Tuple{String,Int},Vector{Float64}},
)
    normalized_forward = Dict(key => _normalize_patch(values)
                              for (key, values) in forward_vectors)
    normalized_backward = Dict(key => _normalize_patch(values)
                               for (key, values) in backward_vectors)
    by_file = Dict{String,Vector{NodeInput}}()
    for node in values(nodes)
        push!(get!(by_file, node.file, NodeInput[]), node)
    end
    drafts = EdgeDraft[]
    ordered_nodes = NodeInput[]
    node_segment = Dict{Tuple{String,Int},String}()
    node_connected = Dict{Tuple{String,Int},Bool}()
    left_boundary = Dict{Tuple{String,Int},String}()
    right_boundary = Dict{Tuple{String,Int},String}()

    for file in sort!(collect(keys(by_file)))
        chain = sort!(by_file[file]; by=node -> (node.t_nm, node.lobe))
        append!(ordered_nodes, chain)
        local_drafts = EdgeDraft[]
        for index in 1:max(0, length(chain) - 1)
            left = chain[index]
            right = chain[index + 1]
            gap = right.t_nm - left.t_nm
            gap_ok = SPACING_MIN_NM <= gap <= SPACING_MAX_NM
            left_key = (left.file, left.lobe)
            right_key = (right.file, right.lobe)
            left_forward = get(normalized_forward, left_key, nothing)
            right_forward = get(normalized_forward, right_key, nothing)
            left_backward = get(normalized_backward, left_key, nothing)
            right_backward = get(normalized_backward, right_key, nothing)
            forward_status = _edge_patch_status(left_forward, right_forward)
            backward_status = _edge_patch_status(left_backward, right_backward)
            forward_dot = forward_status == :valid ?
                _centered_dot(left_forward::PatchVector, right_forward::PatchVector) :
                (status=:unavailable, value=NaN)
            backward_dot = backward_status == :valid ?
                _centered_dot(left_backward::PatchVector, right_backward::PatchVector) :
                (status=:unavailable, value=NaN)
            reason = _select_split_reason(
                gap_ok,
                forward_status,
                backward_status,
                forward_dot.status,
                backward_dot.status,
            )
            eligible = reason == :none
            push!(local_drafts, EdgeDraft(
                file,
                left,
                right,
                gap,
                eligible ? :eligible : :split,
                reason,
                eligible ? forward_dot.value : nothing,
                eligible ? backward_dot.value : nothing,
            ))
        end

        segment_indices = ones(Int, length(chain))
        for index in eachindex(local_drafts)
            segment_indices[index + 1] = segment_indices[index] +
                                         (local_drafts[index].status == :split ? 1 : 0)
        end
        segment_lobes = Dict{Int,Vector{Int}}()
        for (node, segment_index) in zip(chain, segment_indices)
            push!(get!(segment_lobes, segment_index, Int[]), node.lobe)
        end
        segment_ids = Dict(index => _segment_id(file, lobes)
                           for (index, lobes) in segment_lobes)
        for (index, node) in enumerate(chain)
            key = (node.file, node.lobe)
            node_segment[key] = segment_ids[segment_indices[index]]
            left_reason = index == 1 ? "start" :
                (local_drafts[index - 1].status == :eligible ? "eligible" :
                 String(local_drafts[index - 1].reason))
            right_reason = index == length(chain) ? "end" :
                (local_drafts[index].status == :eligible ? "eligible" :
                 String(local_drafts[index].reason))
            left_boundary[key] = left_reason
            right_boundary[key] = right_reason
            connected_left = index > 1 && local_drafts[index - 1].status == :eligible
            connected_right = index < length(chain) && local_drafts[index].status == :eligible
            node_connected[key] = connected_left || connected_right
        end
        append!(drafts, local_drafts)
    end
    return drafts, ordered_nodes, node_segment, node_connected, left_boundary, right_boundary
end

function _edge_table_bytes(
    drafts::Vector{EdgeDraft},
    node_segment::Dict{Tuple{String,Int},String},
    forward_hash::String,
    backward_hash::String,
    feature_hash::String,
    source_hash::String,
)::Vector{UInt8}
    ordered = sort(drafts; by=edge -> (
        edge.file,
        edge.left.t_nm,
        edge.right.t_nm,
        edge.left.lobe,
        edge.right.lobe,
    ))
    io = IOBuffer()
    println(io, join(EDGE_HEADER, '\t'))
    for edge in ordered
        left_key = (edge.left.file, edge.left.lobe)
        right_key = (edge.right.file, edge.right.lobe)
        eligible = edge.status == :eligible
        values = [
            edge.file,
            string(edge.left.lobe),
            string(edge.right.lobe),
            _format_number(edge.left.t_nm),
            _format_number(edge.right.t_nm),
            _format_number(edge.gap_nm),
            node_segment[left_key],
            node_segment[right_key],
            String(edge.status),
            eligible ? "none" : String(edge.reason),
            eligible ? _format_number(edge.corr_forward::Float64) : "NA",
            eligible ? _format_number(edge.corr_backward::Float64) : "NA",
            GRID_SHA256,
            forward_hash,
            backward_hash,
            feature_hash,
            MODEL_CONFIG_SHA256,
            source_hash,
        ]
        println(io, join(values, '\t'))
    end
    return take!(io)
end

function _node_table_bytes(
    nodes::Vector{NodeInput},
    node_segment::Dict{Tuple{String,Int},String},
    node_connected::Dict{Tuple{String,Int},Bool},
    left_boundary::Dict{Tuple{String,Int},String},
    right_boundary::Dict{Tuple{String,Int},String},
    feature_hash::String,
    source_hash::String,
)::Vector{UInt8}
    ordered = sort(nodes; by=node -> (node.file, node.t_nm, node.lobe))
    io = IOBuffer()
    println(io, join(NODE_HEADER, '\t'))
    for node in ordered
        key = (node.file, node.lobe)
        println(io, join([
            node.file,
            string(node.lobe),
            _format_number(node.t_nm),
            node_segment[key],
            node_connected[key] ? "connected" : "isolated",
            left_boundary[key],
            right_boundary[key],
            feature_hash,
            MODEL_CONFIG_SHA256,
            source_hash,
        ], '\t'))
    end
    return take!(io)
end

function _receipt_bytes(
    edge_bytes::Vector{UInt8},
    node_bytes::Vector{UInt8},
    frozen::StructuredUniverse.FrozenUniverse,
    universe_receipt_hash::String,
    forward::PatchTable,
    backward::PatchTable,
    forward_receipt::FileSnapshot,
    backward_receipt::FileSnapshot,
    forward_command::String,
    backward_command::String,
    source_hash::String,
)::Vector{UInt8}
    edge_rows = count(==(UInt8('\n')), edge_bytes) - 1
    node_rows = count(==(UInt8('\n')), node_bytes) - 1
    fields = [
        "schema" => "stmfit-structured-edge-feature-receipt-v1",
        "schema_version" => 1,
        "status" => "PASS",
        "edge_observations_file" => "edge_observations.tsv",
        "node_segments_file" => "node_segments.tsv",
        "edge_header" => join(EDGE_HEADER, '\t'),
        "node_header" => join(NODE_HEADER, '\t'),
        "edge_header_sha256" => _sha256(codeunits(join(EDGE_HEADER, '\t') * "\n")),
        "node_header_sha256" => _sha256(codeunits(join(NODE_HEADER, '\t') * "\n")),
        "edge_observations_sha256" => _sha256(edge_bytes),
        "node_segments_sha256" => _sha256(node_bytes),
        "edge_row_count" => edge_rows,
        "node_row_count" => node_rows,
        "feature_sha256" => frozen.feature_sha256,
        "candidate_config_sha256" => CANDIDATE_CONFIG_SHA256,
        "model_config_sha256" => MODEL_CONFIG_SHA256,
        "universe_receipt_sha256" => universe_receipt_hash,
        "keys_sha256" => frozen.keys_sha256,
        "grid_sha256" => GRID_SHA256,
        "forward_patch_sha256" => forward.snapshot.sha256,
        "backward_patch_sha256" => backward.snapshot.sha256,
        "forward_patch_receipt_sha256" => forward_receipt.sha256,
        "backward_patch_receipt_sha256" => backward_receipt.sha256,
        "forward_producer_sha256" => FORWARD_PRODUCER_SHA256,
        "backward_producer_sha256" => BACKWARD_PRODUCER_SHA256,
        "forward_command" => forward_command,
        "backward_command" => backward_command,
        "source_sha256" => source_hash,
        "t2_review_sha256" => T2_REVIEW_SHA256,
        "t3_review_sha256" => T3_REVIEW_SHA256,
    ]
    io = IOBuffer()
    for (key, value) in fields
        if value isa Integer
            println(io, key, " = ", value)
        else
            println(io, key, " = ", _toml_quote(String(value)))
        end
    end
    return take!(io)
end

function _bundle_snapshots(
    universe_path::String,
    direct::Vector{FileSnapshot},
)::Vector{FileSnapshot}
    snapshots = copy(direct)
    for name in sort(readdir(universe_path))
        path = joinpath(universe_path, name)
        isfile(path) || _fail(:bundle_mismatch, "universe bundle member is not a file")
        _, snapshot = _snapshot(path, "universe $name")
        push!(snapshots, snapshot)
    end
    paths = [snapshot.path for snapshot in snapshots]
    length(paths) == length(unique(paths)) ||
        _fail(:source_collision, "a bound dependency path appears more than once")
    return snapshots
end

function build_edge_bundle(
    root::AbstractString;
    features::AbstractString,
    feature_sha256::AbstractString,
    candidate_config::AbstractString,
    model_config::AbstractString,
    universe_dir::AbstractString,
    forward_patches::AbstractString,
    backward_patches::AbstractString,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
)::EdgeBundle
    canonical_root = _canonical_root(root)
    _sha256(grid_contract_bytes()) == GRID_SHA256 ||
        _fail(:internal_contract, "grid contract bytes differ from the frozen digest")
    t2_snapshot = _validate_fixed_review(T2_REVIEW_PATH, T2_REVIEW_SHA256,
                                         "Todo 2 correction review")
    t3_snapshot = _validate_fixed_review(T3_REVIEW_PATH, T3_REVIEW_SHA256,
                                         "Todo 3 correction review")
    supplied_feature_hash = _validate_hash(String(feature_sha256), "feature hash")
    frozen = StructuredUniverse.freeze_universe(
        canonical_root;
        features=String(features),
        feature_sha256=supplied_feature_hash,
        candidate_config=String(candidate_config),
        model_config=String(model_config),
        source_paths=[UNIVERSE_SOURCE_PATH, UNIVERSE_CLI_PATH, FIREWALL_SOURCE_PATH],
        candidate_config_sha256=CANDIDATE_CONFIG_SHA256,
        model_config_sha256=MODEL_CONFIG_SHA256,
    )
    universe_path = _resolve_existing_directory(canonical_root, universe_dir,
                                                "universe bundle")
    StructuredUniverse.validate_bundle(frozen, universe_path)

    candidate_bytes, candidate_snapshot = _snapshot(frozen.candidate_config_path,
                                                     "candidate config")
    model_bytes, model_snapshot = _snapshot(frozen.model_config_path, "model config")
    candidate_snapshot.sha256 == CANDIDATE_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "candidate config SHA-256 differs")
    model_snapshot.sha256 == MODEL_CONFIG_SHA256 ||
        _fail(:config_hash_mismatch, "model config SHA-256 differs")
    _validate_candidate_contract(candidate_bytes)
    _validate_model_contract(model_bytes)

    feature_bytes, feature_snapshot = _snapshot(frozen.features_path, "feature table")
    feature_snapshot.sha256 == supplied_feature_hash ||
        _fail(:feature_hash_mismatch, "feature bytes changed after universe binding")
    nodes = _parse_feature_nodes(feature_bytes, frozen)
    forward_path = _resolve_existing(canonical_root, forward_patches,
                                     "forward patch table")
    backward_path = _resolve_existing(canonical_root, backward_patches,
                                      "backward patch table")
    forward = _parse_patch_table(forward_path, :forward, nodes, frozen)
    backward = _parse_patch_table(backward_path, :backward, nodes, frozen)
    Set(keys(forward.vectors)) == Set(keys(backward.vectors)) ||
        _fail(:patch_key_mismatch, "forward and backward patch keys differ")

    forward_receipt_path = _resolve_existing(canonical_root, forward_receipt,
                                             "forward producer receipt")
    backward_receipt_path = _resolve_existing(canonical_root, backward_receipt,
                                              "backward producer receipt")
    forward_provenance = _validate_patch_receipt(
        forward_receipt_path, :forward, canonical_root, frozen.features_path,
        forward, backward.snapshot.path, frozen,
    )
    backward_provenance = _validate_patch_receipt(
        backward_receipt_path, :backward, canonical_root, frozen.features_path,
        backward, forward.snapshot.path, frozen,
    )
    forward_provenance.data_path == backward_provenance.data_path ||
        _fail(:producer_path_mismatch,
              "forward and backward producer receipts bind different data directories")
    source_hash, source_snapshots = _source_snapshots()
    drafts, ordered_nodes, node_segment, node_connected, left_boundary, right_boundary =
        _build_topology(nodes, forward.vectors, backward.vectors)
    edge_bytes = _edge_table_bytes(
        drafts,
        node_segment,
        forward.snapshot.sha256,
        backward.snapshot.sha256,
        frozen.feature_sha256,
        source_hash,
    )
    node_bytes = _node_table_bytes(
        ordered_nodes,
        node_segment,
        node_connected,
        left_boundary,
        right_boundary,
        frozen.feature_sha256,
        source_hash,
    )
    universe_receipt_path = joinpath(universe_path, "receipt.toml")
    _, universe_receipt_snapshot = _snapshot(universe_receipt_path,
                                             "universe receipt")
    receipt_bytes = _receipt_bytes(
        edge_bytes,
        node_bytes,
        frozen,
        universe_receipt_snapshot.sha256,
        forward,
        backward,
        forward_provenance.snapshot,
        backward_provenance.snapshot,
        forward_provenance.command,
        backward_provenance.command,
        source_hash,
    )
    files = Dict(
        "edge_observations.tsv" => edge_bytes,
        "node_segments.tsv" => node_bytes,
        "receipt.toml" => receipt_bytes,
    )
    direct_snapshots = FileSnapshot[
        t2_snapshot,
        t3_snapshot,
        feature_snapshot,
        candidate_snapshot,
        model_snapshot,
        forward.snapshot,
        backward.snapshot,
        forward_provenance.snapshot,
        backward_provenance.snapshot,
    ]
    append!(direct_snapshots, source_snapshots)
    snapshots = _bundle_snapshots(universe_path, direct_snapshots)
    StructuredUniverse.validate_binding(frozen)
    for snapshot in snapshots
        _verify_snapshot(snapshot)
    end
    return EdgeBundle(files, frozen, Tuple(snapshots))
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

function _rename_noreplace(source::String, destination::String)
    Sys.islinux() || _fail(
        :atomic_publication_unsupported,
        "atomic no-replace publication requires Linux renameat2",
    )
    result = try
        ccall(
            :renameat2,
            Cint,
            (Cint, Cstring, Cint, Cstring, Cuint),
            _AT_FDCWD,
            source,
            _AT_FDCWD,
            destination,
            _RENAME_NOREPLACE,
        )
    catch error
        _fail(
            :atomic_publication_unsupported,
            "renameat2(RENAME_NOREPLACE) is unavailable: $(sprint(showerror, error))",
        )
    end
    result == 0 && return nothing

    error_number = Int(Libc.errno())
    error_number in (_LINUX_EEXIST, _LINUX_ENOTEMPTY) && _fail(
        :publication_collision,
        "edge bundle destination appeared concurrently",
    )
    error_number in (_LINUX_EINVAL, _LINUX_ENOSYS, _LINUX_EOPNOTSUPP) && _fail(
        :atomic_publication_unsupported,
        "renameat2(RENAME_NOREPLACE) is unsupported: $(Libc.strerror(error_number))",
    )
    _fail(
        :publication_install_failed,
        "atomic no-replace publication failed: $(Libc.strerror(error_number))",
    )
end

function _write_synced(path::String, bytes::Vector{UInt8})
    ispath(path) && _fail(:publication_collision, "publication stage member exists")
    open(path, "w") do io
        write(io, bytes)
        flush(io)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 ||
            _fail(:fsync_failed, "file fsync failed")
    end
    isfile(path) && !islink(path) ||
        _fail(:publication_collision, "publication stage member is not a regular file")
    return nothing
end

function _validate_bound_bundle(bundle::EdgeBundle)
    StructuredUniverse.validate_binding(bundle.frozen)
    for snapshot in bundle.snapshots
        _verify_snapshot(snapshot)
    end
    return nothing
end

function _run_publication_install_hook(stage::String, destination::String)
    hook = _PUBLICATION_INSTALL_HOOK[]
    hook === nothing || hook(stage, destination)
    return nothing
end

function _compare_output(destination::String, required::Dict{String,Vector{UInt8}})
    isdir(destination) && !islink(destination) ||
        _fail(:bundle_mismatch, "edge bundle destination is not a regular directory")
    sort(readdir(destination)) == collect(OUTPUT_NAMES) ||
        _fail(:bundle_mismatch, "edge bundle file set is not exact")
    for name in OUTPUT_NAMES
        path = joinpath(destination, name)
        isfile(path) && !islink(path) ||
            _fail(:bundle_mismatch, "edge bundle member is not a regular file")
        read(path) == required[name] ||
            _fail(:bundle_mismatch, "edge bundle member bytes differ")
    end
    return Dict(name => _sha256(required[name]) for name in OUTPUT_NAMES)
end

function _publish_edge_bundle(
    root::AbstractString,
    output_directory::AbstractString;
    features::AbstractString,
    feature_sha256::AbstractString,
    candidate_config::AbstractString,
    model_config::AbstractString,
    universe_dir::AbstractString,
    forward_patches::AbstractString,
    backward_patches::AbstractString,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
    failpoint::Union{Nothing,Symbol}=nothing,
)
    failpoint in (nothing, :before_install) ||
        _fail(:invalid_failpoint, "unsupported publication failpoint")
    canonical_root = _canonical_root(root)
    destination = _resolve_output(canonical_root, output_directory, "edge bundle")
    keywords = (
        features=features,
        feature_sha256=feature_sha256,
        candidate_config=candidate_config,
        model_config=model_config,
        universe_dir=universe_dir,
        forward_patches=forward_patches,
        backward_patches=backward_patches,
        forward_receipt=forward_receipt,
        backward_receipt=backward_receipt,
    )
    bundle = build_edge_bundle(canonical_root; keywords...)
    if ispath(destination)
        _validate_bound_bundle(bundle)
        return _compare_output(destination, bundle.files)
    end
    parent = dirname(destination)
    stage = mktempdir(parent; prefix=".structured-edge-stage-", cleanup=false)
    installed = false
    try
        for name in OUTPUT_NAMES
            _write_synced(joinpath(stage, name), bundle.files[name])
        end
        _fsync_directory(stage)
        failpoint == :before_install &&
            _fail(:interrupted_publication, "simulated interruption before install")
        _validate_bound_bundle(bundle)
        _run_publication_install_hook(stage, destination)
        _rename_noreplace(stage, destination)
        installed = true
        _fsync_directory(parent)
        return _compare_output(destination, bundle.files)
    finally
        !installed && ispath(stage) && rm(stage; recursive=true, force=true)
    end
end

function publish_edge_bundle(
    root::AbstractString,
    output_directory::AbstractString;
    features::AbstractString,
    feature_sha256::AbstractString,
    candidate_config::AbstractString,
    model_config::AbstractString,
    universe_dir::AbstractString,
    forward_patches::AbstractString,
    backward_patches::AbstractString,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
)
    return _publish_edge_bundle(
        root,
        output_directory;
        features=features,
        feature_sha256=feature_sha256,
        candidate_config=candidate_config,
        model_config=model_config,
        universe_dir=universe_dir,
        forward_patches=forward_patches,
        backward_patches=backward_patches,
        forward_receipt=forward_receipt,
        backward_receipt=backward_receipt,
    )
end

function validate_edge_bundle(
    root::AbstractString,
    output_directory::AbstractString;
    features::AbstractString,
    feature_sha256::AbstractString,
    candidate_config::AbstractString,
    model_config::AbstractString,
    universe_dir::AbstractString,
    forward_patches::AbstractString,
    backward_patches::AbstractString,
    forward_receipt::AbstractString,
    backward_receipt::AbstractString,
)
    canonical_root = _canonical_root(root)
    destination = _resolve_output(canonical_root, output_directory, "edge bundle")
    bundle = build_edge_bundle(
        canonical_root;
        features=features,
        feature_sha256=feature_sha256,
        candidate_config=candidate_config,
        model_config=model_config,
        universe_dir=universe_dir,
        forward_patches=forward_patches,
        backward_patches=backward_patches,
        forward_receipt=forward_receipt,
        backward_receipt=backward_receipt,
    )
    _validate_bound_bundle(bundle)
    return _compare_output(destination, bundle.files)
end

end
