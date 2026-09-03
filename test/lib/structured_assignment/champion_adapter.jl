# Immutable unary adapter for one frozen, hash-bound source.

using Printf
using SHA

const FROZEN_CHAMPION_SHA256 =
    "86c7b0985413cf4ceb61ee903ed9bd27e375b5cc2ed314fa1b15b9648dd4ecf5"
const FROZEN_CHAMPION_BYTE_COUNT = 27_451
const FROZEN_CHAMPION_ROW_COUNT = 892
const FROZEN_CHAMPION_KEY_ORDER_SHA256 =
    "ca21d9004b3ba95711d143359999b223fcbf6e69601454f5f9d1ccc04df6e1ba"
const _REPOSITORY_ROOT = realpath(joinpath(@__DIR__, "..", "..", ".."))
const FROZEN_CHAMPION_PATH = joinpath(
    _REPOSITORY_ROOT, "results", "unit_assignment",
    "best_labelfree_cc_soft_20260802.tsv")
const _FROZEN_COLUMNS = ("file", "lobe", "predicted", "confidence")
const _DENIED_DECLARATIONS = (
    "bench" * "mark",
    "tr" * "uth",
    "gr" * "ade",
    "gr" * "ader",
)

struct ChampionAdapterError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, err::ChampionAdapterError) =
    print(io, err.code, ": ", err.message)

struct FrozenChampionSource
    path::String
    sha256::String
    byte_count::Int
    row_count::Int
    key_order_sha256::String
end

struct BaselineUnaryKey
    file::String
    lobe::Int
end

struct BaselineUnaryEvidence
    predicted::Int
    confidence::Float64
    confidence_text::String
    probability_1::Float64
end

struct BaselineUnaryRecord
    key::BaselineUnaryKey
    unary::BaselineUnaryEvidence
end

struct BaselineUnaryStream{N} <: AbstractVector{BaselineUnaryRecord}
    source::FrozenChampionSource
    records::NTuple{N,BaselineUnaryRecord}
end

struct _SourceSnapshot
    path::String
    required_sha256::String
    bytes::Vector{UInt8}
    identity::NamedTuple
end

Base.IndexStyle(::Type{<:BaselineUnaryStream}) = IndexLinear()
Base.size(stream::BaselineUnaryStream) = (length(stream.records),)
Base.length(stream::BaselineUnaryStream) = length(stream.records)
Base.getindex(stream::BaselineUnaryStream, index::Int) = stream.records[index]

Base.:(==)(left::FrozenChampionSource, right::FrozenChampionSource) =
    left.path == right.path && left.sha256 == right.sha256 &&
    left.byte_count == right.byte_count && left.row_count == right.row_count &&
    left.key_order_sha256 == right.key_order_sha256
Base.:(==)(left::BaselineUnaryKey, right::BaselineUnaryKey) =
    left.file == right.file && left.lobe == right.lobe
Base.:(==)(left::BaselineUnaryEvidence, right::BaselineUnaryEvidence) =
    left.predicted == right.predicted && isequal(left.confidence, right.confidence) &&
    left.confidence_text == right.confidence_text &&
    isequal(left.probability_1, right.probability_1)
Base.:(==)(left::BaselineUnaryRecord, right::BaselineUnaryRecord) =
    left.key == right.key && left.unary == right.unary

_fail(code::Symbol, message::String) = throw(ChampionAdapterError(code, message))
_sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))

function _has_parent_component(path::String)
    normalized = replace(path, '\\' => '/')
    return any(==(".."), split(normalized, '/'; keepempty=true))
end

function _absolute_path(path::AbstractString)
    supplied = String(path)
    isempty(supplied) && _fail(:invalid_source_path, "source path is empty")
    occursin('\0', supplied) && _fail(:invalid_source_path, "source path contains NUL")
    _has_parent_component(supplied) &&
        _fail(:source_path_substitution, "parent traversal is not the frozen path")
    return normpath(isabspath(supplied) ? supplied : joinpath(_REPOSITORY_ROOT, supplied))
end

function _reject_symlink_components(path::String)
    parts = splitpath(abspath(path))
    cursor = first(parts)
    islink(cursor) && _fail(:symlink_rejected, "source path contains a symbolic link")
    for component in parts[2:end]
        cursor = joinpath(cursor, component)
        islink(cursor) && _fail(:symlink_rejected, "source path contains a symbolic link")
    end
    return nothing
end


function _resolve_source_path(path::AbstractString, required_path::AbstractString)
    candidate = _absolute_path(path)
    required = _absolute_path(required_path)
    _reject_symlink_components(candidate)
    _reject_symlink_components(required)
    candidate == required ||
        _fail(:source_path_substitution, "source path is not the frozen source")
    isfile(candidate) || _fail(:missing_source, "frozen source is not a regular file")
    return candidate
end

function _identity(info)
    return (
        device=info.device,
        inode=info.inode,
        mode=info.mode,
        size=info.size,
        mtime=info.mtime,
        ctime=info.ctime,
    )
end

function _read_source_snapshot(path::AbstractString, required_sha256::String)
    occursin(r"^[0-9a-f]{64}$", required_sha256) ||
        _fail(:invalid_source_hash, "required source digest is malformed")
    absolute = _absolute_path(path)
    _reject_symlink_components(absolute)
    isfile(absolute) || _fail(:missing_source, "source is not a regular file")
    before = _identity(stat(absolute))
    bytes = open(absolute, "r") do io
        opened = _identity(stat(io))
        opened == before || _fail(:source_changed, "source identity changed while opening")
        content = read(io)
        _identity(stat(io)) == opened ||
            _fail(:source_changed, "source changed during the first read")
        content
    end
    _identity(stat(absolute)) == before ||
        _fail(:source_changed, "source path changed after the first read")
    _sha256_hex(bytes) == required_sha256 ||
        _fail(:source_hash_mismatch, "source digest differs before parsing")
    return _SourceSnapshot(absolute, required_sha256, bytes, before)
end

function _verify_source_snapshot(snapshot::_SourceSnapshot)
    _reject_symlink_components(snapshot.path)
    isfile(snapshot.path) || _fail(:source_changed, "source disappeared before return")
    _identity(stat(snapshot.path)) == snapshot.identity ||
        _fail(:source_changed, "source identity changed before return")
    bytes = open(snapshot.path, "r") do io
        _identity(stat(io)) == snapshot.identity ||
            _fail(:source_changed, "source path was substituted before return")
        content = read(io)
        _identity(stat(io)) == snapshot.identity ||
            _fail(:source_changed, "source changed during the final read")
        content
    end
    bytes == snapshot.bytes || _fail(:source_changed, "source bytes changed before return")
    _sha256_hex(bytes) == snapshot.required_sha256 ||
        _fail(:source_changed, "source digest changed before return")
    return nothing
end

function _contains_denied_declaration(columns)
    for column in columns
        compact = replace(lowercase(String(column)), r"[^a-z0-9]+" => "")
        any(fragment -> occursin(fragment, compact), _DENIED_DECLARATIONS) && return true
    end
    return false
end

function _parse_confidence(value_text::AbstractString)
    text = String(value_text)
    text == strip(text) || _fail(:malformed_confidence, "confidence has whitespace")
    value = tryparse(Float64, text)
    value === nothing && _fail(:malformed_confidence, "confidence is not numerical")
    isfinite(value) || _fail(:nonfinite_confidence, "confidence is nonfinite")
    (value < 0.0 || value > 1.0 || signbit(value)) &&
        _fail(:confidence_out_of_range, "confidence is outside [0,1]")
    @sprintf("%.8f", value) == text ||
        _fail(:malformed_confidence, "confidence does not use frozen serialization")
    return value
end

function _unary_evidence(predicted::Int, confidence::Float64, confidence_text::String)
    predicted in (0, 1) || _fail(:invalid_prediction, "prediction must be binary")
    isfinite(confidence) || _fail(:nonfinite_confidence, "confidence is nonfinite")
    (confidence < 0.0 || confidence > 1.0 || signbit(confidence)) &&
        _fail(:confidence_out_of_range, "confidence is outside [0,1]")
    probability = predicted == 1 ? 0.5 + confidence / 2 : 0.5 - confidence / 2
    (isfinite(probability) && 0.0 <= probability <= 1.0) ||
        _fail(:invalid_probability, "reconstructed unary probability is invalid")
    return BaselineUnaryEvidence(predicted, confidence, confidence_text, probability)
end

function _parse_frozen_rows(bytes::Vector{UInt8})
    isvalid(String, bytes) || _fail(:invalid_encoding, "source is not valid UTF-8")
    isempty(bytes) && _fail(:malformed_source, "source is empty")
    last(bytes) == UInt8('\n') || _fail(:malformed_source, "source lacks final LF")
    UInt8('\r') in bytes && _fail(:malformed_source, "source contains CR bytes")
    text = String(copy(bytes))
    occursin('\0', text) && _fail(:malformed_source, "source contains NUL")
    lines = split(text, '\n'; keepempty=true)
    pop!(lines) == "" || _fail(:malformed_source, "source termination is malformed")
    isempty(lines) && _fail(:malformed_source, "source has no header")
    columns = split(first(lines), '\t'; keepempty=true)
    _contains_denied_declaration(columns) &&
        _fail(:forbidden_declaration, "source declares a forbidden field")
    Tuple(columns) == _FROZEN_COLUMNS ||
        _fail(:invalid_columns, "source columns differ from the frozen schema")

    records = BaselineUnaryRecord[]
    seen = Set{Tuple{String,Int}}()
    key_material = IOBuffer()
    for (row_index, line) in enumerate(lines[2:end])
        fields = split(line, '\t'; keepempty=true)
        length(fields) == 4 ||
            _fail(:invalid_row_columns, "row $row_index has the wrong column count")
        file = String(fields[1])
        (!isempty(file) && all(isascii, file) && basename(file) == file) ||
            _fail(:invalid_key, "row $row_index has an invalid file key")
        occursin(r"^[1-9][0-9]*$", fields[2]) ||
            _fail(:invalid_key, "row $row_index has an invalid lobe key")
        lobe = tryparse(Int, fields[2])
        lobe === nothing && _fail(:invalid_key, "row $row_index lobe is outside Int range")
        fields[3] in ("0", "1") ||
            _fail(:invalid_prediction, "row $row_index prediction is invalid")
        predicted = parse(Int, fields[3])
        confidence_text = String(fields[4])
        confidence = _parse_confidence(confidence_text)
        key_tuple = (file, lobe)
        key_tuple in seen && _fail(:duplicate_key, "row $row_index repeats a key")
        push!(seen, key_tuple)
        print(key_material, file, '\t', lobe, '\n')
        push!(records, BaselineUnaryRecord(
            BaselineUnaryKey(file, lobe),
            _unary_evidence(predicted, confidence, confidence_text),
        ))
    end
    length(records) == FROZEN_CHAMPION_ROW_COUNT ||
        _fail(:row_count_mismatch, "source row count differs from the frozen identity")
    _sha256_hex(take!(key_material)) == FROZEN_CHAMPION_KEY_ORDER_SHA256 ||
        _fail(:key_order_mismatch, "source key order differs from the frozen identity")
    return records
end

function _load_source_unaries(
    path::AbstractString,
    required_path::AbstractString,
    required_sha256::String,
    required_byte_count::Union{Nothing,Int}=nothing,
)
    source_path = _resolve_source_path(path, required_path)
    snapshot = _read_source_snapshot(source_path, required_sha256)
    records = _parse_frozen_rows(snapshot.bytes)
    source = FrozenChampionSource(
        source_path,
        required_sha256,
        length(snapshot.bytes),
        length(records),
        FROZEN_CHAMPION_KEY_ORDER_SHA256,
    )
    stream = BaselineUnaryStream(source, Tuple(records))
    required_byte_count === nothing || length(snapshot.bytes) == required_byte_count ||
        _fail(:source_size_mismatch, "source byte count differs from the frozen identity")
    _verify_source_snapshot(snapshot)
    return stream
end

function load_frozen_champion_unaries(
    path::AbstractString=FROZEN_CHAMPION_PATH,
)
    return _load_source_unaries(
        path,
        FROZEN_CHAMPION_PATH,
        FROZEN_CHAMPION_SHA256,
        FROZEN_CHAMPION_BYTE_COUNT,
    )
end
