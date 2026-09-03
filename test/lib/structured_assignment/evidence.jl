module StructuredEvidence

using SHA
using TOML

export ARTIFACT_INDEX_SCHEMA,
       GATE_RECEIPT_SCHEMA,
       TERMINAL_STATUSES,
       ArtifactSpec,
       ArtifactRecord,
       ArtifactIndex,
       GateReceipt,
       EvidenceResult,
       build_artifact_index,
       artifact_index_bytes,
       canonical_artifact_index_bytes,
       artifact_index_canonical_sha256,
       role_bundle_sha256,
       parse_artifact_index,
       parse_gate_receipt,
       gate_receipt_bytes,
       canonical_gate_receipt_bytes,
       gate_receipt_canonical_sha256,
       with_noncanonical_metadata,
       validate_artifact_index,
       validate_evidence_bundle,
       publish_evidence_bundle,
       publication_lock_path,
       command_identity_sha256,
       is_v2_selection_eligible

const ARTIFACT_INDEX_SCHEMA = "stmfit-structured-artifact-index-v1"
const GATE_RECEIPT_SCHEMA = "stmfit-structured-gate-receipt-v1"
const TERMINAL_STATUSES = (
    "PASS",
    "FAIL",
    "BLOCKED",
    "SKIPPED",
    "FOLLOWUP_REQUIRED",
)
const REQUIRED_ARTIFACT_ROLES = ("config", "source", "input", "output", "expected_rows")
const NONCANONICAL_METADATA_KEYS = ("repository_root", "generated_at", "exact_command")

const ARTIFACT_INDEX_HEADER = [
    "schema",
    "schema_version",
    "campaign_id",
    "model_id",
    "config_sha256",
    "source_sha256",
    "input_sha256",
    "output_sha256",
    "expected_rows_sha256",
    "terminal_status",
    "selection_eligible",
    "diagnostic_only",
    "reason_codes",
    "reasons",
    "artifact_role",
    "artifact_path",
    "byte_count",
    "artifact_sha256",
    "repository_root",
    "generated_at",
    "exact_command",
]
const CANONICAL_ARTIFACT_INDEX_HEADER = ARTIFACT_INDEX_HEADER[1:18]

struct EvidenceFailure <: Exception
    code::String
    message::String
end

Base.showerror(io::IO, failure::EvidenceFailure) =
    print(io, failure.code, ": ", failure.message)

struct ArtifactSpec
    role::String
    path::String
end

struct ArtifactRecord
    role::String
    path::String
    byte_count::Int
    sha256::String
end

struct ArtifactIndex
    schema::String
    schema_version::Int
    campaign_id::String
    model_id::String
    config_sha256::String
    source_sha256::String
    input_sha256::String
    output_sha256::String
    expected_rows_sha256::String
    terminal_status::String
    selection_eligible::Bool
    diagnostic_only::Bool
    reason_codes::Vector{String}
    reasons::Vector{String}
    artifacts::Vector{ArtifactRecord}
    metadata::Dict{String,String}
end

struct GateReceipt
    schema::String
    schema_version::Int
    campaign_id::String
    model_id::String
    config_sha256::String
    source_sha256::String
    input_sha256::String
    output_sha256::String
    expected_rows_sha256::String
    terminal_status::String
    selection_eligible::Bool
    diagnostic_only::Bool
    reason_codes::Vector{String}
    reasons::Vector{String}
    artifact_index_path::String
    artifact_index_sha256::String
    artifact_index_canonical_sha256::String
    canonical_sha256::String
    metadata::Dict{String,String}
end

struct EvidenceResult
    valid::Bool
    published::Bool
    terminal_status::String
    campaign_id::String
    model_id::String
    identity_verified::Bool
    command_identity_sha256::String
    durability_status::String
    warnings::Vector{String}
    reason_codes::Vector{String}
    reasons::Vector{String}
    artifact_index_path::String
    receipt_path::String
    artifact_index_sha256::String
    receipt_sha256::String
    canonical_sha256::String
end

_sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))
_sha256_hex(bytes::Base.CodeUnits) = bytes2hex(sha256(bytes))

function command_identity_sha256(root::String, exact_command::String)::String
    canonical_root = _canonical_root(root)
    _validate_plain_text(exact_command, "exact command")
    normalized = replace(exact_command, canonical_root => "\$ROOT")
    material = "stmfit-structured-command-identity-v1\ncommand=$normalized\n"
    return _sha256_hex(codeunits(material))
end

function _failure(code::String, message::String)
    throw(EvidenceFailure(code, message))
end

function _as_failure(error)::EvidenceFailure
    error isa EvidenceFailure && return error
    return EvidenceFailure("internal_error", sprint(showerror, error))
end

function _blocked(
    error;
    campaign_id::String="",
    model_id::String="",
    command_identity_sha256::String="",
    artifact_index_path::String="",
    receipt_path::String="",
)::EvidenceResult
    failure = _as_failure(error)
    return EvidenceResult(
        false,
        false,
        "BLOCKED",
        campaign_id,
        model_id,
        false,
        command_identity_sha256,
        "NOT_APPLICABLE",
        String[],
        [failure.code],
        [failure.message],
        artifact_index_path,
        receipt_path,
        "",
        "",
        "",
    )
end

function _valid_result(
    terminal_status::String;
    published::Bool,
    campaign_id::String,
    model_id::String,
    identity_verified::Bool,
    command_identity_sha256::String,
    durability_status::String="NOT_ASSERTED",
    warnings::Vector{String}=String[],
    reason_codes::Vector{String},
    reasons::Vector{String},
    artifact_index_path::String,
    receipt_path::String="",
    artifact_index_sha256::String,
    receipt_sha256::String="",
    canonical_sha256::String,
)::EvidenceResult
    durability_status in ("CONFIRMED", "NOT_ASSERTED", "WARNING") ||
        _failure("internal_error", "invalid durability status: $durability_status")
    durability_status == "WARNING" && isempty(warnings) &&
        _failure("internal_error", "WARNING durability status requires a warning")
    durability_status != "WARNING" && !isempty(warnings) &&
        _failure("internal_error", "durability warnings require WARNING status")
    return EvidenceResult(
        true,
        published,
        terminal_status,
        campaign_id,
        model_id,
        identity_verified,
        command_identity_sha256,
        durability_status,
        copy(warnings),
        copy(reason_codes),
        copy(reasons),
        artifact_index_path,
        receipt_path,
        artifact_index_sha256,
        receipt_sha256,
        canonical_sha256,
    )
end

function _toml_quote(value::String)::String
    io = IOBuffer()
    write(io, UInt8('"'))
    for character in value
        codepoint = Int(character)
        if character == '"'
            write(io, "\\\"")
        elseif character == '\\'
            write(io, "\\\\")
        elseif character == '\b'
            write(io, "\\b")
        elseif character == '\t'
            write(io, "\\t")
        elseif character == '\n'
            write(io, "\\n")
        elseif character == '\f'
            write(io, "\\f")
        elseif character == '\r'
            write(io, "\\r")
        elseif codepoint < 0x20 || codepoint == 0x7f
            write(io, "\\u", uppercase(string(codepoint; base=16, pad=4)))
        else
            write(io, string(character))
        end
    end
    write(io, UInt8('"'))
    return String(take!(io))
end

function _toml_string_array(values::Vector{String})::String
    return "[" * join((_toml_quote(value) for value in values), ", ") * "]"
end

function _parse_toml_string_array(text::AbstractString, context::String)::Vector{String}
    document = try
        TOML.parse("value = $(String(text))\n")
    catch error
        _failure("schema_invalid", "$context is not a TOML string array: $(sprint(showerror, error))")
    end
    Set(String(key) for key in keys(document)) == Set(["value"]) ||
        _failure("schema_invalid", "$context contains unexpected TOML keys")
    value = document["value"]
    value isa Vector || _failure("schema_invalid", "$context must be an array")
    all(item -> item isa String, value) ||
        _failure("schema_invalid", "$context must contain only strings")
    return String[String(item) for item in value]
end

function _validate_plain_text(
    value::String,
    context::String;
    allow_empty::Bool=false,
)::String
    !allow_empty && isempty(value) && _failure("invalid_text", "$context must not be empty")
    occursin('\0', value) && _failure("invalid_text", "$context contains NUL")
    for character in ('\t', '\n', '\r')
        occursin(character, value) &&
            _failure("invalid_text", "$context must be one line without TSV controls")
    end
    return value
end

function _validate_identifier(value::String, context::String)::String
    _validate_plain_text(value, context)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", value) ||
        _failure("invalid_identifier", "$context must use only A-Z, a-z, 0-9, dot, underscore, and hyphen")
    return value
end

function _validate_reason_code(value::String)::String
    _validate_plain_text(value, "reason code")
    occursin(r"^[a-z][a-z0-9_]*$", value) ||
        _failure("invalid_reasons", "reason code must be lowercase snake_case: $value")
    return value
end

function _validate_sha256(value::String, context::String)::String
    occursin(r"^[0-9a-f]{64}$", value) ||
        _failure("schema_invalid", "$context must be a lowercase SHA-256")
    return value
end

function _terminal_semantics(terminal_status::String)::Tuple{Bool,Bool}
    terminal_status in TERMINAL_STATUSES ||
        _failure("invalid_terminal_status", "terminal status must be exactly one of $(join(TERMINAL_STATUSES, '|'))")
    return terminal_status == "PASS", terminal_status == "FOLLOWUP_REQUIRED"
end

function _validate_reasons(reason_codes::Vector{String}, reasons::Vector{String})
    isempty(reason_codes) && _failure("invalid_reasons", "at least one reason code is required")
    length(reason_codes) == length(reasons) ||
        _failure("invalid_reasons", "reason_codes and reasons must have equal length")
    length(unique(reason_codes)) == length(reason_codes) ||
        _failure("invalid_reasons", "reason codes must be unique")
    foreach(_validate_reason_code, reason_codes)
    for reason in reasons
        _validate_plain_text(reason, "reason")
    end
    return nothing
end

function _validate_metadata(metadata::Dict{String,String})
    Set(keys(metadata)) == Set(NONCANONICAL_METADATA_KEYS) ||
        _failure("schema_invalid", "noncanonical metadata must contain exactly $(join(NONCANONICAL_METADATA_KEYS, ", "))")
    root = _validate_plain_text(metadata["repository_root"], "metadata.repository_root")
    isabspath(root) ||
        _failure("schema_invalid", "metadata.repository_root must be absolute")
    _validate_plain_text(metadata["generated_at"], "metadata.generated_at")
    _validate_plain_text(metadata["exact_command"], "metadata.exact_command")
    return nothing
end

function _canonical_root(root::String)::String
    _validate_plain_text(root, "repository root")
    isdir(root) || _failure("invalid_root", "repository root is not a directory: $root")
    islink(root) && _failure("symlink_rejected", "repository root is a symlink: $root")
    absolute = normpath(abspath(root))
    canonical = realpath(root)
    absolute == canonical ||
        _failure("noncanonical_path", "repository root must be its canonical real path")
    return canonical
end

function _validate_lexical_path(path::String, context::String)
    _validate_plain_text(path, context)
    occursin('\\', path) &&
        _failure("noncanonical_path", "$context must use '/' separators")
    components = split(path, '/'; keepempty=true)
    for (index, component) in enumerate(components)
        if component == ".."
            _failure("path_escape", "$context contains parent traversal: $path")
        elseif component == "."
            _failure("noncanonical_path", "$context contains a dot component: $path")
        elseif isempty(component) && !(isabspath(path) && index == 1)
            _failure("noncanonical_path", "$context contains an empty component: $path")
        end
    end
    normpath(path) == path ||
        _failure("noncanonical_path", "$context is not lexically canonical: $path")
    return nothing
end

function _relative_path(root::String, candidate::String, context::String)::String
    relative_native = relpath(candidate, root)
    parts = splitpath(relative_native)
    (relative_native == "." || (!isempty(parts) && first(parts) == "..")) &&
        _failure("path_escape", "$context escapes the repository root: $candidate")
    relative = join(parts, "/")
    _validate_relative_record_path(relative, context)
    return relative
end

function _validate_relative_record_path(path::String, context::String)::String
    isabspath(path) && _failure("schema_invalid", "$context must be repository-relative")
    _validate_lexical_path(path, context)
    isempty(path) && _failure("schema_invalid", "$context must not be empty")
    return path
end

function _candidate_path(root::String, supplied::String, context::String)::Tuple{String,String}
    _validate_lexical_path(supplied, context)
    candidate = if isabspath(supplied)
        supplied
    else
        joinpath(root, split(supplied, '/')...)
    end
    candidate = normpath(abspath(candidate))
    relative = _relative_path(root, candidate, context)
    return candidate, relative
end

function _reject_symlink_components(root::String, relative::String, context::String)
    cursor = root
    for component in split(relative, '/')
        cursor = joinpath(cursor, component)
        islink(cursor) &&
            _failure("symlink_rejected", "$context contains a symlink component: $relative")
    end
    return nothing
end

function _resolve_existing_file(
    root::String,
    supplied::String,
    context::String;
    missing_code::String="missing_artifact",
    not_regular_code::String="artifact_not_regular",
)::Tuple{String,String}
    candidate, relative = _candidate_path(root, supplied, context)
    _reject_symlink_components(root, relative, context)
    ispath(candidate) || _failure(missing_code, "$context is missing: $relative")
    isfile(candidate) || _failure(not_regular_code, "$context is not a regular file: $relative")
    return candidate, relative
end

function _resolve_directory(root::String, supplied::String, context::String)::Tuple{String,String}
    candidate, relative = _candidate_path(root, supplied, context)
    _reject_symlink_components(root, relative, context)
    isdir(candidate) || _failure("invalid_output_directory", "$context is not a directory: $relative")
    return candidate, relative
end

function _snapshot_bytes(path::String, context::String)::Vector{UInt8}
    islink(path) && _failure("symlink_rejected", "$context became a symlink")
    isfile(path) || _failure("artifact_not_regular", "$context is not a regular file")
    before = stat(path)
    bytes = read(path)
    islink(path) && _failure("symlink_rejected", "$context became a symlink while read")
    isfile(path) || _failure("artifact_not_regular", "$context stopped being a regular file while read")
    after = stat(path)
    stable = before.device == after.device &&
             before.inode == after.inode &&
             before.size == after.size &&
             before.mtime == after.mtime &&
             before.ctime == after.ctime
    stable || _failure("artifact_changed_during_read", "$context changed while its bytes were hashed")
    length(bytes) == after.size ||
        _failure("artifact_changed_during_read", "$context byte count changed while read")
    return bytes
end

function _record_from_spec(root::String, spec::ArtifactSpec)::ArtifactRecord
    _validate_reason_code(spec.role)
    absolute, relative = _resolve_existing_file(root, spec.path, "artifact $(spec.role)")
    bytes = _snapshot_bytes(absolute, "artifact $relative")
    return ArtifactRecord(spec.role, relative, length(bytes), _sha256_hex(bytes))
end

function _role_records(records::Vector{ArtifactRecord}, role::String)::Vector{ArtifactRecord}
    return sort!(filter(record -> record.role == role, records); by=record -> record.path)
end

function _role_bundle_bytes(records::Vector{ArtifactRecord}, role::String)::Vector{UInt8}
    selected = _role_records(records, role)
    isempty(selected) && _failure("schema_invalid", "artifact role $role has no records")
    io = IOBuffer()
    write(io, "stmfit-structured-artifact-role-bundle-v1\n")
    write(io, "role=$role\n")
    for record in selected
        write(io, "path=$(record.path)\n")
        write(io, "byte_count=$(record.byte_count)\n")
        write(io, "sha256=$(record.sha256)\n")
    end
    return take!(io)
end

role_bundle_sha256(index::ArtifactIndex, role::String) =
    _sha256_hex(_role_bundle_bytes(index.artifacts, role))

function _validate_artifact_structure(records::Vector{ArtifactRecord})
    isempty(records) && _failure("schema_invalid", "artifact index must not be empty")
    paths = String[]
    for record in records
        _validate_reason_code(record.role)
        _validate_relative_record_path(record.path, "artifact path")
        record.byte_count >= 0 ||
            _failure("schema_invalid", "artifact byte_count must be nonnegative")
        _validate_sha256(record.sha256, "artifact SHA-256")
        push!(paths, record.path)
    end
    length(unique(paths)) == length(paths) ||
        _failure("duplicate_artifact", "artifact index contains duplicate canonical paths")
    records == sort(records; by=record -> (record.role, record.path)) ||
        _failure("schema_invalid", "artifact records must be sorted by role and path")

    for role in REQUIRED_ARTIFACT_ROLES
        role_count = Base.count(record -> record.role == role, records)
        if role in ("config", "expected_rows")
            role_count == 1 ||
                _failure("schema_invalid", "artifact index requires exactly one $role artifact")
        else
            role_count >= 1 ||
                _failure("schema_invalid", "artifact index requires at least one $role artifact")
        end
    end
    return nothing
end

function _validate_index_structure(index::ArtifactIndex)
    index.schema == ARTIFACT_INDEX_SCHEMA ||
        _failure("schema_invalid", "artifact index schema is invalid")
    index.schema_version == 1 ||
        _failure("schema_invalid", "artifact index schema_version must equal 1")
    _validate_identifier(index.campaign_id, "campaign_id")
    _validate_identifier(index.model_id, "model_id")
    selection_eligible, diagnostic_only = _terminal_semantics(index.terminal_status)
    index.selection_eligible == selection_eligible ||
        _failure("schema_invalid", "selection_eligible contradicts terminal status")
    index.diagnostic_only == diagnostic_only ||
        _failure("schema_invalid", "diagnostic_only contradicts terminal status")
    _validate_reasons(index.reason_codes, index.reasons)
    _validate_metadata(index.metadata)
    _validate_artifact_structure(index.artifacts)

    config = only(_role_records(index.artifacts, "config"))
    expected_rows = only(_role_records(index.artifacts, "expected_rows"))
    expected = Dict(
        "config_sha256" => config.sha256,
        "source_sha256" => role_bundle_sha256(index, "source"),
        "input_sha256" => role_bundle_sha256(index, "input"),
        "output_sha256" => role_bundle_sha256(index, "output"),
        "expected_rows_sha256" => expected_rows.sha256,
    )
    for (field, digest) in expected
        actual = getfield(index, Symbol(field))
        _validate_sha256(actual, field)
        actual == digest ||
            _failure("schema_invalid", "$field is inconsistent with artifact records")
    end
    return nothing
end

function build_artifact_index(
    root::String,
    specs::Vector{ArtifactSpec};
    campaign_id::String,
    model_id::String,
    exact_command::String,
    generated_at::String,
    terminal_status::String,
    reason_codes::Vector{String},
    reasons::Vector{String},
)::ArtifactIndex
    canonical_root = _canonical_root(root)
    _validate_identifier(campaign_id, "campaign_id")
    _validate_identifier(model_id, "model_id")
    selection_eligible, diagnostic_only = _terminal_semantics(terminal_status)
    _validate_reasons(reason_codes, reasons)
    metadata = Dict(
        "repository_root" => canonical_root,
        "generated_at" => generated_at,
        "exact_command" => exact_command,
    )
    _validate_metadata(metadata)

    records = ArtifactRecord[_record_from_spec(canonical_root, spec) for spec in specs]
    sort!(records; by=record -> (record.role, record.path))
    _validate_artifact_structure(records)
    config = only(_role_records(records, "config"))
    expected_rows = only(_role_records(records, "expected_rows"))
    provisional = ArtifactIndex(
        ARTIFACT_INDEX_SCHEMA,
        1,
        campaign_id,
        model_id,
        config.sha256,
        repeat('0', 64),
        repeat('0', 64),
        repeat('0', 64),
        expected_rows.sha256,
        terminal_status,
        selection_eligible,
        diagnostic_only,
        copy(reason_codes),
        copy(reasons),
        records,
        metadata,
    )
    index = ArtifactIndex(
        provisional.schema,
        provisional.schema_version,
        provisional.campaign_id,
        provisional.model_id,
        provisional.config_sha256,
        role_bundle_sha256(provisional, "source"),
        role_bundle_sha256(provisional, "input"),
        role_bundle_sha256(provisional, "output"),
        provisional.expected_rows_sha256,
        provisional.terminal_status,
        provisional.selection_eligible,
        provisional.diagnostic_only,
        provisional.reason_codes,
        provisional.reasons,
        provisional.artifacts,
        provisional.metadata,
    )
    _validate_index_structure(index)
    return index
end

function _artifact_index_fields(
    index::ArtifactIndex,
    record::ArtifactRecord;
    canonical::Bool,
)::Vector{String}
    fields = [
        index.schema,
        string(index.schema_version),
        index.campaign_id,
        index.model_id,
        index.config_sha256,
        index.source_sha256,
        index.input_sha256,
        index.output_sha256,
        index.expected_rows_sha256,
        index.terminal_status,
        index.selection_eligible ? "true" : "false",
        index.diagnostic_only ? "true" : "false",
        _toml_string_array(index.reason_codes),
        _toml_string_array(index.reasons),
        record.role,
        record.path,
        string(record.byte_count),
        record.sha256,
    ]
    if !canonical
        append!(fields, [
            index.metadata["repository_root"],
            index.metadata["generated_at"],
            index.metadata["exact_command"],
        ])
    end
    return fields
end

function _serialize_artifact_index(index::ArtifactIndex; canonical::Bool)::Vector{UInt8}
    _validate_index_structure(index)
    header = canonical ? CANONICAL_ARTIFACT_INDEX_HEADER : ARTIFACT_INDEX_HEADER
    io = IOBuffer()
    write(io, join(header, '\t'), '\n')
    for record in index.artifacts
        write(io, join(_artifact_index_fields(index, record; canonical), '\t'), '\n')
    end
    return take!(io)
end

artifact_index_bytes(index::ArtifactIndex) =
    _serialize_artifact_index(index; canonical=false)

canonical_artifact_index_bytes(index::ArtifactIndex) =
    _serialize_artifact_index(index; canonical=true)

artifact_index_canonical_sha256(index::ArtifactIndex) =
    _sha256_hex(canonical_artifact_index_bytes(index))

function _parse_integer(text::AbstractString, context::String)::Int
    value = tryparse(Int, text)
    value === nothing && _failure("schema_invalid", "$context must be an integer")
    return value
end

function _parse_bool(text::AbstractString, context::String)::Bool
    text == "true" && return true
    text == "false" && return false
    _failure("schema_invalid", "$context must be true or false")
end

function parse_artifact_index(bytes::Vector{UInt8})::ArtifactIndex
    text = try
        String(copy(bytes))
    catch error
        _failure("schema_invalid", "artifact index is not valid UTF-8: $(sprint(showerror, error))")
    end
    startswith(text, "\ufeff") && _failure("schema_invalid", "artifact index must not contain a BOM")
    occursin('\r', text) && _failure("schema_invalid", "artifact index must use LF line endings")
    endswith(text, "\n") || _failure("schema_invalid", "artifact index must be LF-terminated")
    lines = split(text, '\n'; keepempty=true)
    pop!(lines) == "" || _failure("schema_invalid", "artifact index termination is invalid")
    length(lines) >= 2 || _failure("schema_invalid", "artifact index has no artifact rows")
    split(first(lines), '\t'; keepempty=true) == ARTIFACT_INDEX_HEADER ||
        _failure("schema_invalid", "artifact index header is invalid")

    first_fields = split(lines[2], '\t'; keepempty=true)
    length(first_fields) == length(ARTIFACT_INDEX_HEADER) ||
        _failure("schema_invalid", "artifact index row has the wrong field count")
    shared_indices = vcat(collect(1:14), collect(19:21))
    records = ArtifactRecord[]
    for (row_number, line) in enumerate(lines[2:end])
        fields = split(line, '\t'; keepempty=true)
        length(fields) == length(ARTIFACT_INDEX_HEADER) ||
            _failure("schema_invalid", "artifact index row $row_number has the wrong field count")
        all(fields[index] == first_fields[index] for index in shared_indices) ||
            _failure("schema_invalid", "artifact index row $row_number changes shared metadata")
        byte_count = _parse_integer(fields[17], "artifact byte_count")
        push!(records, ArtifactRecord(
            String(fields[15]),
            String(fields[16]),
            byte_count,
            String(fields[18]),
        ))
    end

    metadata = Dict(
        "repository_root" => String(first_fields[19]),
        "generated_at" => String(first_fields[20]),
        "exact_command" => String(first_fields[21]),
    )
    index = ArtifactIndex(
        String(first_fields[1]),
        _parse_integer(first_fields[2], "artifact index schema_version"),
        String(first_fields[3]),
        String(first_fields[4]),
        String(first_fields[5]),
        String(first_fields[6]),
        String(first_fields[7]),
        String(first_fields[8]),
        String(first_fields[9]),
        String(first_fields[10]),
        _parse_bool(first_fields[11], "selection_eligible"),
        _parse_bool(first_fields[12], "diagnostic_only"),
        _parse_toml_string_array(first_fields[13], "reason_codes"),
        _parse_toml_string_array(first_fields[14], "reasons"),
        records,
        metadata,
    )
    _validate_index_structure(index)
    artifact_index_bytes(index) == bytes ||
        _failure("schema_invalid", "artifact index is not in canonical full serialization")
    return index
end

function with_noncanonical_metadata(
    index::ArtifactIndex;
    repository_root::String,
    generated_at::String,
    exact_command::String,
)::ArtifactIndex
    updated = ArtifactIndex(
        index.schema,
        index.schema_version,
        index.campaign_id,
        index.model_id,
        index.config_sha256,
        index.source_sha256,
        index.input_sha256,
        index.output_sha256,
        index.expected_rows_sha256,
        index.terminal_status,
        index.selection_eligible,
        index.diagnostic_only,
        copy(index.reason_codes),
        copy(index.reasons),
        copy(index.artifacts),
        Dict(
            "repository_root" => repository_root,
            "generated_at" => generated_at,
            "exact_command" => exact_command,
        ),
    )
    _validate_index_structure(updated)
    return updated
end

function _validate_current_artifacts(root::String, index::ArtifactIndex)
    for record in index.artifacts
        absolute, relative = _resolve_existing_file(
            root,
            record.path,
            "artifact $(record.role)";
            missing_code="missing_artifact",
            not_regular_code="artifact_not_regular",
        )
        relative == record.path ||
            _failure("noncanonical_path", "artifact path changed under canonical resolution")
        bytes = _snapshot_bytes(absolute, "artifact $(record.path)")
        length(bytes) == record.byte_count ||
            _failure("artifact_size_mismatch", "artifact $(record.path) byte count is stale")
        _sha256_hex(bytes) == record.sha256 ||
            _failure("artifact_hash_mismatch", "artifact $(record.path) SHA-256 is stale")
    end
    return nothing
end

function validate_artifact_index(root::String, path::String)::EvidenceResult
    absolute = ""
    try
        canonical_root = _canonical_root(root)
        absolute, _ = _resolve_existing_file(
            canonical_root,
            path,
            "artifact index";
            missing_code="missing_artifact_index",
            not_regular_code="artifact_index_not_regular",
        )
        bytes = _snapshot_bytes(absolute, "artifact index")
        index = parse_artifact_index(bytes)
        _validate_current_artifacts(canonical_root, index)
        return _valid_result(
            index.terminal_status;
            published=false,
            campaign_id=index.campaign_id,
            model_id=index.model_id,
            identity_verified=false,
            command_identity_sha256="",
            reason_codes=index.reason_codes,
            reasons=index.reasons,
            artifact_index_path=absolute,
            artifact_index_sha256=_sha256_hex(bytes),
            canonical_sha256=artifact_index_canonical_sha256(index),
        )
    catch error
        return _blocked(error; artifact_index_path=absolute)
    end
end

function _receipt_from_index(
    root::String,
    index_path::String,
    index::ArtifactIndex,
    index_bytes::Vector{UInt8},
)::GateReceipt
    _, relative_index = _candidate_path(root, index_path, "artifact index path")
    raw_sha256 = _sha256_hex(index_bytes)
    basename(relative_index) == "artifact-index-$raw_sha256.tsv" ||
        _failure("schema_invalid", "artifact index filename must be content-addressed by current bytes")
    canonical_index_sha256 = artifact_index_canonical_sha256(index)
    provisional = GateReceipt(
        GATE_RECEIPT_SCHEMA,
        1,
        index.campaign_id,
        index.model_id,
        index.config_sha256,
        index.source_sha256,
        index.input_sha256,
        index.output_sha256,
        index.expected_rows_sha256,
        index.terminal_status,
        index.selection_eligible,
        index.diagnostic_only,
        copy(index.reason_codes),
        copy(index.reasons),
        relative_index,
        raw_sha256,
        canonical_index_sha256,
        repeat('0', 64),
        copy(index.metadata),
    )
    digest = _sha256_hex(canonical_gate_receipt_bytes(provisional))
    receipt = GateReceipt(
        provisional.schema,
        provisional.schema_version,
        provisional.campaign_id,
        provisional.model_id,
        provisional.config_sha256,
        provisional.source_sha256,
        provisional.input_sha256,
        provisional.output_sha256,
        provisional.expected_rows_sha256,
        provisional.terminal_status,
        provisional.selection_eligible,
        provisional.diagnostic_only,
        provisional.reason_codes,
        provisional.reasons,
        provisional.artifact_index_path,
        provisional.artifact_index_sha256,
        provisional.artifact_index_canonical_sha256,
        digest,
        provisional.metadata,
    )
    _validate_receipt_structure(receipt)
    return receipt
end

function _validate_receipt_structure(receipt::GateReceipt)
    receipt.schema == GATE_RECEIPT_SCHEMA ||
        _failure("schema_invalid", "gate receipt schema is invalid")
    receipt.schema_version == 1 ||
        _failure("schema_invalid", "gate receipt schema_version must equal 1")
    _validate_identifier(receipt.campaign_id, "campaign_id")
    _validate_identifier(receipt.model_id, "model_id")
    selection_eligible, diagnostic_only = _terminal_semantics(receipt.terminal_status)
    receipt.selection_eligible == selection_eligible ||
        _failure("schema_invalid", "receipt selection_eligible contradicts terminal status")
    receipt.diagnostic_only == diagnostic_only ||
        _failure("schema_invalid", "receipt diagnostic_only contradicts terminal status")
    _validate_reasons(receipt.reason_codes, receipt.reasons)
    _validate_metadata(receipt.metadata)
    for (name, digest) in (
        "config_sha256" => receipt.config_sha256,
        "source_sha256" => receipt.source_sha256,
        "input_sha256" => receipt.input_sha256,
        "output_sha256" => receipt.output_sha256,
        "expected_rows_sha256" => receipt.expected_rows_sha256,
        "artifact_index_sha256" => receipt.artifact_index_sha256,
        "artifact_index_canonical_sha256" => receipt.artifact_index_canonical_sha256,
        "canonical_sha256" => receipt.canonical_sha256,
    )
        _validate_sha256(digest, name)
    end
    _validate_relative_record_path(receipt.artifact_index_path, "artifact_index_path")
    basename(receipt.artifact_index_path) ==
        "artifact-index-$(receipt.artifact_index_sha256).tsv" ||
        _failure("schema_invalid", "receipt artifact_index_path is not content-addressed")
    expected_canonical = _sha256_hex(canonical_gate_receipt_bytes(receipt))
    receipt.canonical_sha256 == expected_canonical ||
        _failure("canonical_hash_mismatch", "gate receipt canonical SHA-256 is inconsistent")
    return nothing
end

function _write_receipt_common(io::IO, receipt::GateReceipt; canonical::Bool)
    write(io, "schema = ", _toml_quote(receipt.schema), '\n')
    write(io, "schema_version = $(receipt.schema_version)\n")
    write(io, "campaign_id = ", _toml_quote(receipt.campaign_id), '\n')
    write(io, "model_id = ", _toml_quote(receipt.model_id), '\n')
    write(io, "config_sha256 = ", _toml_quote(receipt.config_sha256), '\n')
    write(io, "source_sha256 = ", _toml_quote(receipt.source_sha256), '\n')
    write(io, "input_sha256 = ", _toml_quote(receipt.input_sha256), '\n')
    write(io, "output_sha256 = ", _toml_quote(receipt.output_sha256), '\n')
    write(io, "expected_rows_sha256 = ", _toml_quote(receipt.expected_rows_sha256), '\n')
    write(io, "terminal_status = ", _toml_quote(receipt.terminal_status), '\n')
    write(io, "selection_eligible = ", receipt.selection_eligible ? "true\n" : "false\n")
    write(io, "diagnostic_only = ", receipt.diagnostic_only ? "true\n" : "false\n")
    write(io, "reason_codes = ", _toml_string_array(receipt.reason_codes), '\n')
    write(io, "reasons = ", _toml_string_array(receipt.reasons), '\n')
    if canonical
        write(io, "artifact_index_canonical_sha256 = ",
            _toml_quote(receipt.artifact_index_canonical_sha256), '\n')
    else
        write(io, "artifact_index_path = ", _toml_quote(receipt.artifact_index_path), '\n')
        write(io, "artifact_index_sha256 = ", _toml_quote(receipt.artifact_index_sha256), '\n')
        write(io, "artifact_index_canonical_sha256 = ",
            _toml_quote(receipt.artifact_index_canonical_sha256), '\n')
        write(io, "canonical_sha256 = ", _toml_quote(receipt.canonical_sha256), '\n')
    end
    return nothing
end

function canonical_gate_receipt_bytes(receipt::GateReceipt)::Vector{UInt8}
    io = IOBuffer()
    _write_receipt_common(io, receipt; canonical=true)
    return take!(io)
end

function gate_receipt_bytes(receipt::GateReceipt)::Vector{UInt8}
    _validate_receipt_structure(receipt)
    io = IOBuffer()
    _write_receipt_common(io, receipt; canonical=false)
    write(io, "\n[metadata]\n")
    for key in NONCANONICAL_METADATA_KEYS
        write(io, key, " = ", _toml_quote(receipt.metadata[key]), '\n')
    end
    return take!(io)
end

gate_receipt_canonical_sha256(receipt::GateReceipt) =
    _sha256_hex(canonical_gate_receipt_bytes(receipt))

function _expect_keys(table::AbstractDict, expected::Vector{String}, context::String)
    actual = Set(String(key) for key in keys(table))
    wanted = Set(expected)
    actual == wanted ||
        _failure("schema_invalid", "$context keys differ; expected $(sort(expected)), got $(sort!(collect(actual)))")
    return nothing
end

function _expect_string(table::AbstractDict, key::String, context::String)::String
    value = get(table, key, nothing)
    value isa String || _failure("schema_invalid", "$context.$key must be a string")
    return value
end

function _expect_integer(table::AbstractDict, key::String, context::String)::Int
    value = get(table, key, nothing)
    value isa Integer && !(value isa Bool) ||
        _failure("schema_invalid", "$context.$key must be an integer")
    try
        return Int(value)
    catch
        _failure("schema_invalid", "$context.$key is outside Int range")
    end
end

function _expect_bool(table::AbstractDict, key::String, context::String)::Bool
    value = get(table, key, nothing)
    value isa Bool || _failure("schema_invalid", "$context.$key must be a boolean")
    return value
end

function _expect_string_array(table::AbstractDict, key::String, context::String)::Vector{String}
    value = get(table, key, nothing)
    value isa Vector || _failure("schema_invalid", "$context.$key must be an array")
    all(item -> item isa String, value) ||
        _failure("schema_invalid", "$context.$key must contain only strings")
    return String[String(item) for item in value]
end

function parse_gate_receipt(bytes::Vector{UInt8})::GateReceipt
    text = try
        String(copy(bytes))
    catch error
        _failure("schema_invalid", "gate receipt is not valid UTF-8: $(sprint(showerror, error))")
    end
    startswith(text, "\ufeff") && _failure("schema_invalid", "gate receipt must not contain a BOM")
    occursin('\r', text) && _failure("schema_invalid", "gate receipt must use LF line endings")
    endswith(text, "\n") || _failure("schema_invalid", "gate receipt must be LF-terminated")
    document = try
        TOML.parse(text)
    catch error
        _failure("schema_invalid", "gate receipt TOML parse failed: $(sprint(showerror, error))")
    end
    top_level = [
        "schema", "schema_version", "campaign_id", "model_id", "config_sha256",
        "source_sha256", "input_sha256", "output_sha256", "expected_rows_sha256",
        "terminal_status", "selection_eligible", "diagnostic_only", "reason_codes",
        "reasons", "artifact_index_path", "artifact_index_sha256",
        "artifact_index_canonical_sha256", "canonical_sha256", "metadata",
    ]
    _expect_keys(document, top_level, "gate receipt")
    metadata_value = get(document, "metadata", nothing)
    metadata_value isa AbstractDict ||
        _failure("schema_invalid", "gate receipt.metadata must be a TOML table")
    _expect_keys(metadata_value, collect(NONCANONICAL_METADATA_KEYS), "gate receipt.metadata")
    metadata = Dict(key => _expect_string(metadata_value, key, "metadata")
                    for key in NONCANONICAL_METADATA_KEYS)
    receipt = GateReceipt(
        _expect_string(document, "schema", "gate receipt"),
        _expect_integer(document, "schema_version", "gate receipt"),
        _expect_string(document, "campaign_id", "gate receipt"),
        _expect_string(document, "model_id", "gate receipt"),
        _expect_string(document, "config_sha256", "gate receipt"),
        _expect_string(document, "source_sha256", "gate receipt"),
        _expect_string(document, "input_sha256", "gate receipt"),
        _expect_string(document, "output_sha256", "gate receipt"),
        _expect_string(document, "expected_rows_sha256", "gate receipt"),
        _expect_string(document, "terminal_status", "gate receipt"),
        _expect_bool(document, "selection_eligible", "gate receipt"),
        _expect_bool(document, "diagnostic_only", "gate receipt"),
        _expect_string_array(document, "reason_codes", "gate receipt"),
        _expect_string_array(document, "reasons", "gate receipt"),
        _expect_string(document, "artifact_index_path", "gate receipt"),
        _expect_string(document, "artifact_index_sha256", "gate receipt"),
        _expect_string(document, "artifact_index_canonical_sha256", "gate receipt"),
        _expect_string(document, "canonical_sha256", "gate receipt"),
        metadata,
    )
    _validate_receipt_structure(receipt)
    gate_receipt_bytes(receipt) == bytes ||
        _failure("schema_invalid", "gate receipt is not in canonical full serialization")
    return receipt
end

function with_noncanonical_metadata(
    receipt::GateReceipt;
    repository_root::String,
    generated_at::String,
    exact_command::String,
)::GateReceipt
    updated = GateReceipt(
        receipt.schema,
        receipt.schema_version,
        receipt.campaign_id,
        receipt.model_id,
        receipt.config_sha256,
        receipt.source_sha256,
        receipt.input_sha256,
        receipt.output_sha256,
        receipt.expected_rows_sha256,
        receipt.terminal_status,
        receipt.selection_eligible,
        receipt.diagnostic_only,
        copy(receipt.reason_codes),
        copy(receipt.reasons),
        receipt.artifact_index_path,
        receipt.artifact_index_sha256,
        receipt.artifact_index_canonical_sha256,
        receipt.canonical_sha256,
        Dict(
            "repository_root" => repository_root,
            "generated_at" => generated_at,
            "exact_command" => exact_command,
        ),
    )
    _validate_receipt_structure(updated)
    return updated
end

function _receipt_index_consistency(receipt::GateReceipt, index::ArtifactIndex)
    fields = (
        :campaign_id,
        :model_id,
        :config_sha256,
        :source_sha256,
        :input_sha256,
        :output_sha256,
        :expected_rows_sha256,
        :terminal_status,
        :selection_eligible,
        :diagnostic_only,
        :reason_codes,
        :reasons,
        :metadata,
    )
    for field in fields
        getfield(receipt, field) == getfield(index, field) ||
            _failure("receipt_binding_mismatch", "gate receipt field $field differs from artifact index")
    end
    receipt.artifact_index_canonical_sha256 == artifact_index_canonical_sha256(index) ||
        _failure("receipt_binding_mismatch", "gate receipt canonical index hash differs from artifact index")
    return nothing
end

_receipt_filename(campaign_id::String, model_id::String) =
    "gate-receipt-$campaign_id-$model_id.toml"

function validate_evidence_bundle(
    root::String,
    receipt_path::String;
    expected_campaign_id::Union{Nothing,String}=nothing,
    expected_model_id::Union{Nothing,String}=nothing,
    expected_receipt_sha256::Union{Nothing,String}=nothing,
    expected_command_identity_sha256::Union{Nothing,String}=nothing,
    expected_command_sha256::Union{Nothing,String}=nothing,
)::EvidenceResult
    absolute_receipt = ""
    absolute_index = ""
    campaign_id = ""
    model_id = ""
    observed_command_identity = ""
    try
        canonical_root = _canonical_root(root)
        if expected_command_identity_sha256 !== nothing && expected_command_sha256 !== nothing &&
           expected_command_identity_sha256 != expected_command_sha256
            _failure("invalid_expected_identity",
                "expected command identity aliases disagree")
        end
        expected_command = expected_command_identity_sha256 === nothing ?
            expected_command_sha256 : expected_command_identity_sha256
        missing = String[]
        expected_campaign_id === nothing && push!(missing, "expected_campaign_id")
        expected_model_id === nothing && push!(missing, "expected_model_id")
        expected_receipt_sha256 === nothing && push!(missing, "expected_receipt_sha256")
        expected_command === nothing && push!(missing, "expected_command_identity_sha256")
        isempty(missing) || _failure("missing_expected_identity",
            "concrete receipt validation requires $(join(missing, ", "))")

        expected_campaign = something(expected_campaign_id)
        expected_model = something(expected_model_id)
        expected_receipt = something(expected_receipt_sha256)
        expected_command_digest = something(expected_command)
        _validate_identifier(expected_campaign, "expected_campaign_id")
        _validate_identifier(expected_model, "expected_model_id")
        occursin(r"^[0-9a-f]{64}$", expected_receipt) ||
            _failure("invalid_expected_identity",
                "expected_receipt_sha256 must be a lowercase SHA-256")
        occursin(r"^[0-9a-f]{64}$", expected_command_digest) ||
            _failure("invalid_expected_identity",
                "expected command identity must be a lowercase SHA-256")

        absolute_receipt, _ = _resolve_existing_file(
            canonical_root,
            receipt_path,
            "gate receipt";
            missing_code="missing_receipt",
            not_regular_code="receipt_not_regular",
        )
        receipt_bytes = _snapshot_bytes(absolute_receipt, "gate receipt")
        receipt = parse_gate_receipt(receipt_bytes)
        campaign_id = receipt.campaign_id
        model_id = receipt.model_id
        basename(absolute_receipt) ==
            _receipt_filename(receipt.campaign_id, receipt.model_id) ||
            _failure("receipt_path_mismatch",
                "gate receipt filename does not match its bound campaign and model")
        receipt.metadata["repository_root"] == canonical_root ||
            _failure("repository_root_mismatch",
                "gate receipt repository_root differs from the validator canonical root")
        receipt.campaign_id == expected_campaign ||
            _failure("campaign_identity_mismatch",
                "gate receipt campaign $(receipt.campaign_id) does not match expected $expected_campaign")
        receipt.model_id == expected_model ||
            _failure("model_identity_mismatch",
                "gate receipt model $(receipt.model_id) does not match expected $expected_model")
        receipt_sha256 = _sha256_hex(receipt_bytes)
        receipt_sha256 == expected_receipt ||
            _failure("receipt_sha256_mismatch",
                "gate receipt current full bytes do not match the trusted expected SHA-256")
        observed_command_identity = command_identity_sha256(
            canonical_root,
            receipt.metadata["exact_command"],
        )
        observed_command_identity == expected_command_digest ||
            _failure("command_identity_mismatch",
                "gate receipt exact command does not match the trusted command identity")
        absolute_index, relative_index = _resolve_existing_file(
            canonical_root,
            receipt.artifact_index_path,
            "artifact index";
            missing_code="missing_artifact_index",
            not_regular_code="artifact_index_not_regular",
        )
        dirname(absolute_index) == dirname(absolute_receipt) ||
            _failure("receipt_path_mismatch", "receipt and artifact index must share one directory")
        relative_index == receipt.artifact_index_path ||
            _failure("noncanonical_path", "receipt artifact index path changed under canonical resolution")
        index_bytes = _snapshot_bytes(absolute_index, "artifact index")
        _sha256_hex(index_bytes) == receipt.artifact_index_sha256 ||
            _failure("artifact_index_hash_mismatch", "artifact index current bytes do not match receipt")
        index = parse_artifact_index(index_bytes)
        _receipt_index_consistency(receipt, index)
        index.metadata["repository_root"] == canonical_root ||
            _failure("repository_root_mismatch",
                "artifact index repository_root differs from the validator canonical root")
        _validate_current_artifacts(canonical_root, index)
        return _valid_result(
            receipt.terminal_status;
            published=false,
            campaign_id=receipt.campaign_id,
            model_id=receipt.model_id,
            identity_verified=true,
            command_identity_sha256=observed_command_identity,
            reason_codes=receipt.reason_codes,
            reasons=receipt.reasons,
            artifact_index_path=absolute_index,
            receipt_path=absolute_receipt,
            artifact_index_sha256=receipt.artifact_index_sha256,
            receipt_sha256=receipt_sha256,
            canonical_sha256=receipt.canonical_sha256,
        )
    catch error
        return _blocked(error;
            campaign_id,
            model_id,
            command_identity_sha256=observed_command_identity,
            artifact_index_path=absolute_index,
            receipt_path=absolute_receipt)
    end
end

publication_lock_path(receipt_path::String) =
    joinpath(dirname(receipt_path), ".$(basename(receipt_path)).publish-lock")

function _fsync_file(path::String)
    open(path, "r") do io
        result = ccall(:fsync, Cint, (Cint,), fd(io))
        result == 0 || _failure("fsync_failed", "fsync failed for $path")
    end
    return nothing
end

function _fsync_directory(path::String)
    descriptor = ccall(:open, Cint, (Cstring, Cint), path, 0)
    descriptor >= 0 || _failure("fsync_failed", "could not open directory for fsync: $path")
    try
        result = ccall(:fsync, Cint, (Cint,), descriptor)
        result == 0 || _failure("fsync_failed", "directory fsync failed: $path")
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
    return nothing
end

function _write_stage(path::String, bytes::Vector{UInt8})
    ispath(path) && _failure("publication_collision", "publication stage already exists: $path")
    open(path, "w") do io
        write(io, bytes)
        flush(io)
        result = ccall(:fsync, Cint, (Cint,), fd(io))
        result == 0 || _failure("fsync_failed", "fsync failed for publication stage")
    end
    isfile(path) && !islink(path) ||
        _failure("publication_collision", "publication stage is not a regular file")
    return nothing
end

function _inject_failpoint(actual::Union{Nothing,Symbol}, expected::Symbol)
    actual == expected &&
        _failure("interrupted_publication", "simulated interruption at $expected")
    return nothing
end

function _inject_fsync_failpoint(actual::Union{Nothing,Symbol}, expected::Symbol)
    actual == expected &&
        _failure("fsync_failed", "simulated directory fsync failure at $expected")
    return nothing
end

function _receipt_references_index(
    root::String,
    receipt_path::String,
    index_path::String,
)::Bool
    isempty(root) && return false
    (isfile(receipt_path) && !islink(receipt_path)) || return false
    receipt = try
        parse_gate_receipt(read(receipt_path))
    catch
        return false
    end
    return normpath(joinpath(root, split(receipt.artifact_index_path, '/')...)) ==
           normpath(index_path)
end

function _acquire_publication_lock(lock_path::String)
    islink(lock_path) &&
        _failure("symlink_rejected", "publication lock path is a symlink")
    ispath(lock_path) &&
        _failure("concurrent_publication", "publication lock already exists")
    try
        mkdir(lock_path; mode=0o700)
    catch error
        ispath(lock_path) &&
            _failure("concurrent_publication", "publication lock was acquired concurrently")
        rethrow(error)
    end
    isdir(lock_path) && !islink(lock_path) ||
        _failure("symlink_rejected", "publication lock is not a real directory")
    return nothing
end

function publish_evidence_bundle(
    root::String,
    evidence_directory::String;
    campaign_id::String,
    model_id::String,
    artifacts::Vector{ArtifactSpec},
    exact_command::String,
    generated_at::String,
    terminal_status::String,
    reason_codes::Vector{String},
    reasons::Vector{String},
    failpoint::Union{Nothing,Symbol}=nothing,
)::EvidenceResult
    canonical_root = ""
    absolute_receipt = ""
    absolute_index = ""
    lock_path = ""
    lock_acquired = false
    created_index = false
    receipt_installed = false
    intended_index_sha256 = ""
    intended_receipt_sha256 = ""
    intended_command_identity = ""
    try
        failpoint in (nothing, :before_index_install, :after_index_install,
                      :before_receipt_install, :index_directory_fsync,
                      :receipt_directory_fsync) ||
            _failure("invalid_failpoint", "unsupported publication failpoint")
        canonical_root = _canonical_root(root)
        absolute_directory, _ = _resolve_directory(
            canonical_root,
            evidence_directory,
            "evidence directory",
        )
        index = build_artifact_index(
            canonical_root,
            artifacts;
            campaign_id,
            model_id,
            exact_command,
            generated_at,
            terminal_status,
            reason_codes,
            reasons,
        )
        index_bytes = artifact_index_bytes(index)
        index_sha256 = _sha256_hex(index_bytes)
        intended_index_sha256 = index_sha256
        absolute_index = joinpath(absolute_directory, "artifact-index-$index_sha256.tsv")
        absolute_receipt = joinpath(
            absolute_directory,
            _receipt_filename(index.campaign_id, index.model_id),
        )
        index_relative = _relative_path(canonical_root, absolute_index, "artifact index destination")
        receipt_relative = _relative_path(canonical_root, absolute_receipt, "gate receipt destination")
        any(record -> record.path in (index_relative, receipt_relative), index.artifacts) &&
            _failure("publication_collision", "evidence destinations overlap a bound artifact")

        for (path, context) in ((absolute_index, "artifact index destination"),
                                (absolute_receipt, "gate receipt destination"))
            islink(path) && _failure("symlink_rejected", "$context is a symlink")
            ispath(path) && !isfile(path) &&
                _failure("publication_collision", "$context is not a regular file")
        end

        lock_path = publication_lock_path(absolute_receipt)
        _acquire_publication_lock(lock_path)
        lock_acquired = true
        stat(lock_path).device == stat(absolute_directory).device ||
            _failure("cross_filesystem_publication", "publication stages are not on the destination filesystem")
        _validate_current_artifacts(canonical_root, index)

        receipt = _receipt_from_index(canonical_root, absolute_index, index, index_bytes)
        receipt_bytes = gate_receipt_bytes(receipt)
        intended_receipt_sha256 = _sha256_hex(receipt_bytes)
        intended_command_identity = command_identity_sha256(
            canonical_root,
            receipt.metadata["exact_command"],
        )
        parse_artifact_index(index_bytes)
        parse_gate_receipt(receipt_bytes)

        index_stage = joinpath(lock_path, "artifact-index.stage")
        receipt_stage = joinpath(lock_path, "gate-receipt.stage")
        _write_stage(index_stage, index_bytes)
        _write_stage(receipt_stage, receipt_bytes)
        _inject_failpoint(failpoint, :before_index_install)

        if ispath(absolute_index)
            islink(absolute_index) &&
                _failure("symlink_rejected", "content-addressed artifact index is a symlink")
            isfile(absolute_index) ||
                _failure("publication_collision", "content-addressed artifact index is not a regular file")
            existing = _snapshot_bytes(absolute_index, "content-addressed artifact index")
            existing == index_bytes ||
                _failure("content_address_collision", "artifact index path exists with different bytes")
            rm(index_stage; force=true)
        else
            islink(absolute_index) &&
                _failure("symlink_rejected", "artifact index destination became a symlink")
            Base.Filesystem.rename(index_stage, absolute_index)
            created_index = true
            _inject_fsync_failpoint(failpoint, :index_directory_fsync)
            _fsync_directory(absolute_directory)
        end
        _inject_failpoint(failpoint, :after_index_install)

        _validate_current_artifacts(canonical_root, index)
        islink(absolute_receipt) &&
            _failure("symlink_rejected", "gate receipt destination became a symlink")
        ispath(absolute_receipt) && !isfile(absolute_receipt) &&
            _failure("publication_collision", "gate receipt destination is not a regular file")
        _inject_failpoint(failpoint, :before_receipt_install)
        Base.Filesystem.rename(receipt_stage, absolute_receipt)
        receipt_installed = true
        _inject_fsync_failpoint(failpoint, :receipt_directory_fsync)
        _fsync_directory(absolute_directory)

        validation = validate_evidence_bundle(
            canonical_root,
            absolute_receipt;
            expected_campaign_id=receipt.campaign_id,
            expected_model_id=receipt.model_id,
            expected_receipt_sha256=intended_receipt_sha256,
            expected_command_identity_sha256=intended_command_identity,
        )
        validation.valid || _failure("publication_validation_failed",
            "published evidence did not independently validate: $(join(validation.reasons, "; "))")
        return _valid_result(
            receipt.terminal_status;
            published=true,
            campaign_id=receipt.campaign_id,
            model_id=receipt.model_id,
            identity_verified=true,
            command_identity_sha256=intended_command_identity,
            durability_status="CONFIRMED",
            reason_codes=receipt.reason_codes,
            reasons=receipt.reasons,
            artifact_index_path=absolute_index,
            receipt_path=absolute_receipt,
            artifact_index_sha256=index_sha256,
            receipt_sha256=intended_receipt_sha256,
            canonical_sha256=receipt.canonical_sha256,
        )
    catch error
        failure = _as_failure(error)
        if receipt_installed && !isempty(intended_receipt_sha256) &&
           !isempty(intended_command_identity)
            authoritative = validate_evidence_bundle(
                canonical_root,
                absolute_receipt;
                expected_campaign_id=campaign_id,
                expected_model_id=model_id,
                expected_receipt_sha256=intended_receipt_sha256,
                expected_command_identity_sha256=intended_command_identity,
            )
            if authoritative.valid && authoritative.identity_verified &&
               authoritative.artifact_index_sha256 == intended_index_sha256 &&
               authoritative.receipt_sha256 == intended_receipt_sha256
                warning = "$(failure.code) after authoritative receipt commit: $(failure.message)"
                return _valid_result(
                    authoritative.terminal_status;
                    published=true,
                    campaign_id=authoritative.campaign_id,
                    model_id=authoritative.model_id,
                    identity_verified=true,
                    command_identity_sha256=authoritative.command_identity_sha256,
                    durability_status="WARNING",
                    warnings=[warning],
                    reason_codes=authoritative.reason_codes,
                    reasons=authoritative.reasons,
                    artifact_index_path=authoritative.artifact_index_path,
                    receipt_path=authoritative.receipt_path,
                    artifact_index_sha256=authoritative.artifact_index_sha256,
                    receipt_sha256=authoritative.receipt_sha256,
                    canonical_sha256=authoritative.canonical_sha256,
                )
            end
        end
        if created_index && !isempty(absolute_index) && isfile(absolute_index) &&
           !_receipt_references_index(canonical_root, absolute_receipt, absolute_index)
            rm(absolute_index; force=true)
        end
        return _blocked(failure;
            campaign_id,
            model_id,
            command_identity_sha256=intended_command_identity,
            artifact_index_path=absolute_index,
            receipt_path=absolute_receipt)
    finally
        if lock_acquired && !isempty(lock_path) && ispath(lock_path)
            rm(lock_path; recursive=true, force=true)
        end
    end
end

is_v2_selection_eligible(::ArtifactIndex) = false

is_v2_selection_eligible(::GateReceipt) = false

is_v2_selection_eligible(result::EvidenceResult) =
    result.valid && result.terminal_status == "PASS" && result.identity_verified &&
    occursin(r"^[0-9a-f]{64}$", result.receipt_sha256) &&
    occursin(r"^[0-9a-f]{64}$", result.command_identity_sha256)

end # module StructuredEvidence
