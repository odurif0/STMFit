# Parse-boundary firewall for structured label-free unit assignment.

using Unicode

const JOIN_TOPOLOGY_METADATA = Set(["file", "lobe", "left_lobe", "right_lobe"])
const NUMERICAL_MODEL_FEATURES = Set([
    "amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated",
    "amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
    "bwd_neg_com_t", "bwd_neg_diag45", "split_log_skew", "corr_fwd", "corr_bwd",
])
const PROVENANCE_ONLY_FIELDS = Set([
    "source_path", "config_path", "input_path", "feature_path", "features_path",
    "forward_patch_path", "backward_patch_path", "output_path",
])
const _FORBIDDEN_COMPACT = [
    "benchmark", "truth", "grade", "grader", "labels", "labelpath",
    "knownlabel", "control", "controlsequence", "nknnkn", "010010", "101101",
    "expectedn", "targetn", "expectedcount", "targetcount",
    "expectedunitcount", "targetunitcount", "nselected",
    "lobecount", "chaincount", "chainlength", "classcount", "statecount",
    "composition", "topk", "occupancy", "position", "normalizedpos",
    "ordinal", "ordinalpos", "terminalpos", "centeredpos", "edgedistance",
    "lobeorder", "lobeindex", "runlength", "statefrequency", "transition", "switch",
    "stateprior", "classprior", "mixtureweight", "sequence",
]

struct StructuredFirewallError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, err::StructuredFirewallError) =
    print(io, err.code, ": ", err.message)

struct StructuredFeatureBatch
    feature_names::Vector{String}
    metadata::Vector{Dict{String,String}}
    matrix::Matrix{Float64}
end

_fail(code::Symbol, message::String) = throw(StructuredFirewallError(code, message))

function _hex_value(byte::UInt8)
    UInt8('0') <= byte <= UInt8('9') && return Int(byte - UInt8('0'))
    UInt8('a') <= byte <= UInt8('f') && return Int(byte - UInt8('a')) + 10
    UInt8('A') <= byte <= UInt8('F') && return Int(byte - UInt8('A')) + 10
    return -1
end

function _percent_decode_once(text::String)
    source = collect(codeunits(text))
    output = UInt8[]
    consumed_escape = false
    index = 1
    while index <= length(source)
        if source[index] == UInt8('%') && index + 2 <= length(source)
            high, low = _hex_value(source[index + 1]), _hex_value(source[index + 2])
            if high >= 0 && low >= 0
                push!(output, UInt8(16high + low))
                consumed_escape = true
                index += 3
                continue
            end
        end
        push!(output, source[index])
        index += 1
    end
    decoded = String(output)
    if !isvalid(decoded)
        consumed_escape &&
            _fail(:invalid_encoding, "percent escapes decode to invalid UTF-8")
        return text
    end
    return decoded
end

function _backslash_decode_once(text::String)
    source = collect(codeunits(text))
    output = UInt8[]
    index = 1
    while index <= length(source)
        if source[index] == UInt8('\\') && index < length(source)
            marker = source[index + 1]
            width = marker == UInt8('x') ? 2 : marker == UInt8('u') ? 4 :
                    marker == UInt8('U') ? 8 : 0
            if width > 0 && index + 1 + width <= length(source)
                codepoint = 0
                valid = true
                for offset in 1:width
                    digit = _hex_value(source[index + 1 + offset])
                    if digit < 0
                        valid = false
                        break
                    end
                    codepoint = 16codepoint + digit
                end
                if valid && codepoint <= 0x10ffff && !(0xd800 <= codepoint <= 0xdfff)
                    append!(output, codeunits(string(Char(codepoint))))
                    index += width + 2
                    continue
                end
            end
        end
        push!(output, source[index])
        index += 1
    end
    return String(output)
end

function _html_decode_once(text::String)
    source = collect(codeunits(text))
    output = UInt8[]
    index = 1
    while index <= length(source)
        if index + 3 <= length(source) && source[index] == UInt8('&') &&
           source[index + 1] == UInt8('#')
            cursor = index + 2
            base = 10
            if cursor <= length(source) && source[cursor] in (UInt8('x'), UInt8('X'))
                base = 16
                cursor += 1
            end
            digits_start = cursor
            codepoint = 0
            while cursor <= length(source) && source[cursor] != UInt8(';')
                digit = _hex_value(source[cursor])
                (digit < 0 || digit >= base) && break
                codepoint <= (0x10ffff - digit) ÷ base || break
                codepoint = base * codepoint + digit
                cursor += 1
            end
            if cursor > digits_start && cursor <= length(source) &&
               source[cursor] == UInt8(';') && codepoint <= 0x10ffff &&
               !(0xd800 <= codepoint <= 0xdfff)
                append!(output, codeunits(string(Char(codepoint))))
                index = cursor + 1
                continue
            end
        end
        push!(output, source[index])
        index += 1
    end
    return String(output)
end

_normalize_text(text::String) = Unicode.normalize(text; compat=true, casefold=true,
                                                  stripmark=true, stripignore=true)

function _decode_fixed_point(text::String, max_rounds::Int)
    max_rounds > 0 || _fail(:decode_not_converged, "decode round bound is exhausted")
    current = _normalize_text(text)
    for _ in 1:max_rounds
        decoded = _html_decode_once(
            _backslash_decode_once(_percent_decode_once(current)))
        candidate = _normalize_text(decoded)
        candidate == current && return current
        ncodeunits(candidate) < ncodeunits(current) ||
            _fail(:decode_not_converged, "canonical decoding made no bounded progress")
        current = candidate
    end
    _fail(:decode_not_converged, "canonical decoding did not reach a fixed point")
end

function _canonical_text(value::AbstractString)
    original = String(value)
    return _decode_fixed_point(original, max(1, ncodeunits(original)))
end

_compact(value::AbstractString) = replace(_canonical_text(value), r"[^a-z0-9]+" => "")

function _has_forbidden_concept(value::AbstractString)
    compact = _compact(value)
    compact == "n" && return true
    return any(token -> occursin(token, compact), _FORBIDDEN_COMPACT)
end

is_join_metadata(name::AbstractString) = String(name) in JOIN_TOPOLOGY_METADATA

function is_forbidden_feature(name::AbstractString)
    raw = strip(String(name))
    isempty(raw) && return true
    canonical = _canonical_text(raw)
    all(isascii, canonical) || return true
    raw in JOIN_TOPOLOGY_METADATA && return true
    return _has_forbidden_concept(raw)
end

function validate_declared_features(features)
    result = String[String(feature) for feature in features]
    length(result) == length(unique(result)) ||
        _fail(:duplicate_declared_feature, "declared numerical features contain duplicates")
    for feature in result
        feature in JOIN_TOPOLOGY_METADATA &&
            _fail(:metadata_as_feature, "join/topology metadata cannot be declared as a feature: $feature")
        is_forbidden_feature(feature) &&
            _fail(:forbidden_feature, "forbidden feature declaration: $feature")
        feature in NUMERICAL_MODEL_FEATURES ||
            _fail(:unknown_declared_feature, "feature is not in the physical-feature allowlist: $feature")
    end
    return result
end

function validate_feature_selection(selected, declared)
    allowed = Set(validate_declared_features(declared))
    result = String[String(feature) for feature in selected]
    length(result) == length(unique(result)) ||
        _fail(:duplicate_feature, "model feature selection contains duplicates")
    for feature in result
        feature in JOIN_TOPOLOGY_METADATA &&
            _fail(:metadata_as_feature, "join/topology metadata cannot enter a feature vector: $feature")
        is_forbidden_feature(feature) &&
            _fail(:forbidden_feature, "forbidden model feature: $feature")
        feature in allowed || _fail(:unknown_feature, "feature was not declared: $feature")
    end
    return result
end

function _validate_provenance_fields(fields)
    result = String[String(field) for field in fields]
    length(result) == length(unique(result)) ||
        _fail(:duplicate_provenance_field, "provenance field list contains duplicates")
    for field in result
        _has_forbidden_concept(field) &&
            _fail(:forbidden_provenance_field, "forbidden provenance field: $field")
        field in PROVENANCE_ONLY_FIELDS ||
            _fail(:unknown_provenance_field, "field is not provenance-only: $field")
    end
    return result
end

function validate_table_header(header, declared; provenance_fields=String[])
    columns = String[String(column) for column in header]
    length(columns) == length(unique(columns)) ||
        _fail(:duplicate_column, "input table contains duplicate columns")
    allowed_features = Set(validate_declared_features(declared))
    provenance = Set(_validate_provenance_fields(provenance_fields))
    for column in columns
        (column in JOIN_TOPOLOGY_METADATA || column in allowed_features ||
         column in provenance) && continue
        is_forbidden_feature(column) &&
            _fail(:forbidden_column, "input table contains a forbidden column: $column")
        _fail(:unknown_column, "input table contains an undeclared column: $column")
    end
    return columns
end

function _looks_like_path(value::AbstractString; basename_is_path::Bool=true)
    text = _canonical_text(strip(String(value)))
    isempty(text) && return false
    occursin(r"(^|[/\\])\.\.($|[/\\])", text) && return true
    occursin('/', text) && return true
    occursin('\\', text) && return true
    occursin(r"^[a-z]:", text) && return true
    occursin(r"^[a-z][a-z0-9+.-]*://", text) && return true
    startswith(text, '~') && return true
    return basename_is_path &&
           occursin(r"\.(tsv|csv|toml|json|jl|py|npy|sxm|txt|log|h5|hdf5)$", text)
end

function _reject_forbidden_value(column::String, value)
    value isa AbstractString || return
    canonical = _canonical_text(value)
    (!all(isascii, canonical) || _has_forbidden_concept(value)) &&
        _fail(:forbidden_value, "column $column contains a forbidden benchmark/control value")
end

function _metadata_value(column::String, value)
    text = strip(string(value))
    isempty(text) && _fail(:invalid_metadata, "metadata $column is empty")
    _reject_forbidden_value(column, text)
    if column == "file"
        _looks_like_path(text; basename_is_path=false) &&
            _fail(:path_value, "file metadata must be a basename, not a path")
        return text
    end
    occursin(r"^[1-9][0-9]*$", text) ||
        _fail(:invalid_metadata, "metadata $column must be a positive integer")
    return string(parse(Int, text))
end

function _numeric_value(column::String, value)
    _reject_forbidden_value(column, value)
    text = strip(string(value))
    text in ("", "NA") && return NaN
    _looks_like_path(text) &&
        _fail(:path_value, "numerical feature $column contains a path-valued cell")
    number = tryparse(Float64, text)
    number === nothing && _fail(:invalid_feature, "feature $column is not numerical")
    isfinite(number) || _fail(:nonfinite_feature, "feature $column is nonfinite")
    return number
end

function build_feature_batch(rows, selected, declared;
                             required_metadata=["file", "lobe"],
                             provenance_fields=String[])
    isempty(rows) && _fail(:empty_input, "feature table has no rows")
    feature_names = validate_feature_selection(selected, declared)
    required = String[String(field) for field in required_metadata]
    length(required) == length(unique(required)) ||
        _fail(:duplicate_metadata, "required metadata contains duplicates")
    all(field -> field in JOIN_TOPOLOGY_METADATA, required) ||
        _fail(:unknown_metadata, "required metadata must use the topology allowlist")
    header = sort!(String[string(key) for key in keys(first(rows))])
    validate_table_header(header, declared; provenance_fields=provenance_fields)
    all(field -> field in header, required) ||
        _fail(:missing_metadata, "feature table is missing required join/topology metadata")

    metadata_fields = [field for field in ("file", "lobe", "left_lobe", "right_lobe")
                       if field in header]
    metadata = Vector{Dict{String,String}}(undef, length(rows))
    matrix = Matrix{Float64}(undef, length(rows), length(feature_names))
    expected_columns = Set(header)
    provenance = Set(_validate_provenance_fields(provenance_fields))
    for (row_index, row) in enumerate(rows)
        Set(string(key) for key in keys(row)) == expected_columns ||
            _fail(:inconsistent_columns, "row $row_index has inconsistent columns")
        row_metadata = Dict{String,String}()
        for field in metadata_fields
            row_metadata[field] = _metadata_value(field, row[field])
        end
        for field in provenance
            _reject_forbidden_value(field, row[field])
        end
        for field in setdiff(expected_columns, union(JOIN_TOPOLOGY_METADATA, provenance))
            _numeric_value(field, row[field])
        end
        for (feature_index, feature) in enumerate(feature_names)
            matrix[row_index, feature_index] = _numeric_value(feature, row[feature])
        end
        metadata[row_index] = row_metadata
    end
    return StructuredFeatureBatch(feature_names, metadata, matrix)
end

function validate_cli_arguments(arguments;
                                allowed_flags,
                                value_flags=allowed_flags,
                                repeatable_flags=Set{String}())
    allowed = Set(String(flag) for flag in allowed_flags)
    valued = Set(String(flag) for flag in value_flags)
    repeatable = Set(String(flag) for flag in repeatable_flags)
    seen = Set{String}()
    index = 1
    while index <= length(arguments)
        argument = String(arguments[index])
        (!all(isascii, _canonical_text(argument)) || _has_forbidden_concept(argument)) &&
            _fail(:forbidden_cli, "forbidden benchmark/control CLI argument")
        startswith(argument, "--") ||
            _fail(:unexpected_positional, "unexpected positional CLI argument")
        pair = split(argument, '='; limit=2)
        flag = pair[1]
        flag in allowed || _fail(:unknown_flag, "unknown CLI flag: $flag")
        flag in seen && !(flag in repeatable) &&
            _fail(:duplicate_flag, "CLI flag may be supplied only once: $flag")
        push!(seen, flag)
        if length(pair) == 2
            flag in valued || _fail(:unexpected_value, "CLI flag takes no value: $flag")
            isempty(pair[2]) && _fail(:missing_value, "CLI flag requires a value: $flag")
        elseif flag in valued
            index < length(arguments) || _fail(:missing_value, "CLI flag requires a value: $flag")
            value = String(arguments[index + 1])
            (!all(isascii, _canonical_text(value)) || _has_forbidden_concept(value)) &&
                _fail(:forbidden_cli, "forbidden benchmark/control CLI value")
            startswith(value, "--") && _fail(:missing_value, "CLI flag requires a value: $flag")
            index += 1
        end
        index += 1
    end
    return nothing
end
