# Private implementation fragment for StructuredChainInference.
# This file intentionally has no module declaration and no product includes.

using Random
using Statistics
using LinearAlgebra

const _PLAN_SHA256 =
    "a9b386d613829e8f7e20b6e33f8e80898fa9a55f0b344dbb92ea64ac6804f3d0"
const _CANDIDATE_SHA256 =
    "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const _MODEL_SHA256 =
    "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const _GRID_SHA256 =
    "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"
const _FORWARD_PRODUCER_SHA256 =
    "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad"
const _BACKWARD_PRODUCER_SHA256 =
    "4f866a1e27289a0cf010a6ed4b5e0b92454dc55f97807d33cc95edb9e87182ba"
const _T2_REVIEW_SHA256 =
    "5555459b99f2a848a53789b2e8e15eb683760ffc0d7056b9e38275f0d65b6aed"
const _T3_REVIEW_SHA256 =
    "531d68e866d896809cdb0c1fde2ac8b934e06a2562ca1fecbb65d84bca5e9954"
const _T10_SOURCE_SHA256 =
    "285809e3711706d7a059ba2a0dd4139b44c1874dec79308e2762a63799aa954e"
const _T10_CLI_SHA256 =
    "c2133b0bc45e5ae68ab49f8ae39470c7bd2ee4fe85a52da548fbc3f84096022b"
const _T10_TEST_SHA256 =
    "21f4896cc1db508ef2dd9d5acc21196a51c7ef066e475468725db73ece98bb22"
const _T3_CLAIM_SHA256 =
    "625738b6cd79ad7b9083263c1bde4d3492d0b9720f58c791b9e16810ba692dfb"
const _T3_REVIEW_FILE_SHA256 = _T3_REVIEW_SHA256
const _T10_CLAIM_SHA256 =
    "895db471d9d607fb79084751d3d9cc37967ff67103750ff69d0c687194309d07"
const _T10_REVIEW_FILE_SHA256 =
    "02d421bb51dcdb4d658da791f14cc7a6cf5f6ced26f3f25216f0b6e6cb88cb77"
const _BASE_T11_CLAIM_SHA256 =
    "c7a9c05ee6bb1896b335170facb5aed657b1834b5448afe6c4b67ef0456a1e35"
const _BASE_T11_REVIEW_SHA256 =
    "e4457cb6678ecb018c8ccd038a1548565a4b50a5cbfd7e7833fa31007b0a50c0"
const _REJECTED_T11_CLAIM_SHA256 =
    "5eaaf639190caed62152cf0105cfd53b222806804432ae84e51a4a5cef5957c6"
const _REJECTED_T11_REVIEW_SHA256 =
    "f7c69427f782bb57870ad8c0834ef922b3d6b40d6c2075d5ba418347d7cf8648"
const _CORRECTION1_CLAIM_SHA256 =
    "83f623cfa35c92cc299e2f3165bd2cd6f97edf828c1c58815679ef6dfbe56cce"
const _CORRECTION1_CLAIM_CORRECTION_SHA256 =
    "ff2a0f9ff02ce07e37633e3c7bf6802d30e966b5837bc565bc1b1d4a378c6f67"
const _CORRECTION1_REVIEW_SHA256 =
    "ecb833d264ba0f9999865b569cf1be0d950a18891813a606765c3b12b695691a"
const _HISTORICAL_T11_CLAIM_SHA256 =
    "008f27ca2bdb56db7f470ffe4249279287556dcd4e3f227598d8cc212bc380c0"
const _HISTORICAL_T11_REVIEW_SHA256 =
    "ba89021d794181c2d10726f5385d74684a1596f875c6b757c1f2b064e3b18819"
const _HISTORICAL_T11_SOURCE_BUNDLE_SHA256 =
    "fd9cec0f28d824d9e1a2c7bafb0d0ec0b9cfb2eed7d634add7ab44d089421cad"
const _PRIOR_SCALE_T11_CLAIM_SHA256 =
    "5dc185994ecd9e117ba9f8dc782e70360daf397a2348e57fd65b2b14a4ffdfe2"
const _PRIOR_SCALE_T11_REVIEW_SHA256 =
    "93f986d6baa3dd8cf850328055c750b27b9873ae67da54f3360ae05c31b22e30"
const _PRIOR_SCALE_T11_SOURCE_BUNDLE_SHA256 =
    "7a024ceba895e670a6cf666b568b50e9bb32da110d1639c93e268d96ed5ecdaa"
const _EFFECTIVE_T11_CLAIM_SHA256 =
    "61b009a55a4ab403b4a2391f3a0a4d5b5bf72b3b21bc684703ac1b52dacd33cd"
const _EFFECTIVE_T11_REVIEW_SHA256 =
    "c02a8fdd553533d844b4d49f13b3cda2658bc02a38589092172d1a4b394cdb83"
const _EFFECTIVE_T11_SOURCE_BUNDLE_SHA256 =
    "b5e54d763523e937734d03742b7323a0f1da191591cf7683a17ff142f01c233b"
const _T11_SOURCE_SHA256 =
    "ca7ea9b05b55a5185a5f08834b4fec4e748830f3eb24cac1a45f387e4d768575"
const _T11_CLI_SHA256 =
    "77c7b2d5f1b9f73ae5c3bdfd310cbbad75771f5a38c5243488f4c2cc6db7ce27"
const _T11_SOURCE_BUNDLE_SHA256 =
    _EFFECTIVE_T11_SOURCE_BUNDLE_SHA256
const _T7_REVIEW_SHA256 =
    "8b0fb2c1a324a345de19d4e718278a52a6d568c0e474939cfb6bfe643e061635"
const _T11_HASH_KEYS = Set([
    "plan_sha256", "t7_review_sha256", "t10_review_sha256",
    "candidate_config_sha256", "model_config_sha256",
    "universe_receipt_sha256", "universe_sha256", "universe_keys_sha256",
    "edge_receipt_sha256", "edge_observations_sha256", "node_segments_sha256",
    "edge_source_sha256", "edge_feature_sha256", "t3_review_sha256",
    "input_sha256", "source_sha256", "normalization_sha256",
])

const _T10_EDGE_HEADER = [
    "file", "left_lobe", "right_lobe", "left_t_nm", "right_t_nm", "gap_nm",
    "left_segment_id", "right_segment_id", "edge_status", "split_reason",
    "corr_fwd", "corr_bwd", "grid_hash", "forward_patch_hash",
    "backward_patch_hash", "feature_hash", "config_hash", "source_hash",
]
const _T10_NODE_HEADER = [
    "file", "lobe", "t_nm", "segment_id", "node_status",
    "left_boundary_reason", "right_boundary_reason", "feature_hash",
    "config_hash", "source_hash",
]
const _T10_RECEIPT_KEYS = Set([
    "schema", "schema_version", "status", "edge_observations_file",
    "node_segments_file", "edge_header", "node_header", "edge_header_sha256",
    "node_header_sha256", "edge_observations_sha256", "node_segments_sha256",
    "edge_row_count", "node_row_count", "feature_sha256", "candidate_config_sha256",
    "model_config_sha256", "universe_receipt_sha256", "keys_sha256", "grid_sha256",
    "forward_patch_sha256", "backward_patch_sha256", "forward_patch_receipt_sha256",
    "backward_patch_receipt_sha256", "forward_producer_sha256",
    "backward_producer_sha256", "forward_command", "backward_command",
    "source_sha256", "t2_review_sha256", "t3_review_sha256",
])

const _GRAPH_HANDOFF_SCHEMA = "stmfit-structured-edge-model-handoff-v1"
const _GRAPH_HANDOFF_FILES = [
    "fitted_edge_models.tsv", "fitted_unary_nodes.tsv", "fitted_edge_transforms.tsv",
]
const _REPORT_FILES = [
    "models.tsv", "partitions.tsv", "ess.tsv", "starts.tsv", "scores.tsv",
    "shuffle.tsv", "bootstrap.tsv",
]
const _MODEL_RESULT_HEADER = [
    "schema_version", "model", "status", "reason", "outer_fold_count",
    "final_gate_sha256", "full_fit_sha256", "result_sha256",
]
const _PARTITION_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "status", "reason", "training_dates", "training_scan_count",
    "heldout_scan_count", "training_sha256", "fit_sha256", "score_sha256",
    "raw_heldout_sha256", "raw_at_score_sha256", "reversal_pass", "shuffle_pass",
]
const _ESS_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date",
    "category", "dates", "scans", "edges", "kish_ess", "minimum_support_scans",
    "free_parameters", "effective_observation_ratio", "sufficient", "reason",
    "fit_sha256",
]
const _START_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date", "fit_role",
    "start_index", "alpha", "status", "reason", "iterations", "objective_initial",
    "objective_final", "rejected_iteration", "rejected_objective", "trace_sha256",
    "fit_sha256",
]
const _SCORE_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date", "level", "file",
    "left_lobe", "right_lobe", "segment_id", "gain", "edge_count", "score_sha256",
]
const _SHUFFLE_HEADER = [
    "schema_version", "model", "scope", "outer_date", "heldout_date", "seed", "date_gain",
    "conditional_fit_sha256", "observed_control_gain", "upper_quantile", "passed", "score_sha256",
]
const _BOOTSTRAP_HEADER = [
    "schema_version", "model", "scope", "outer_date", "seed", "mean_gain", "lower_quantile",
    "every_date_positive", "every_shuffle_pass", "reversal_pass", "status", "reason", "gate_sha256",
]
const _GRAPH_CATEGORY_ORDER = ("00", "01", "10", "11")
const _GRAPH_STATE_ORDER = ("0", "1")
const _GRAPH_OUTPUT_NAMES = ("corr_fwd", "corr_bwd")
const _PRODUCER_MINIMUM_DATES = 5
const _PRODUCER_COMPLETE_DATE_COUNT = _PRODUCER_MINIMUM_DATES + 2
const _DENSITY_FAILURE_HOOK = Ref{Union{Nothing,Function}}(nothing)
const _DP_FAILURE_HOOK = Ref{Union{Nothing,Function}}(nothing)
const _SCALE_EIGENVALUE_FLOOR = 1.0e-4
const _SCALE_CONDITION_CAP = 1.0e4
const _GRAPH_COEFFICIENT_NAMES = (
    "intercept",
    "left_amp_prominence", "left_amp_neighbor_ratio", "left_integrated_prominence",
    "left_amp_rel", "left_bwd_neg_com_t", "left_bwd_neg_diag45", "left_split_log_skew",
    "right_amp_prominence", "right_amp_neighbor_ratio", "right_integrated_prominence",
    "right_amp_rel", "right_bwd_neg_com_t", "right_bwd_neg_diag45", "right_split_log_skew",
)
const _GRAPH_MODEL_HEADER = vcat([
    "schema_version", "fit_role", "fit_sha256", "model", "scope",
    "outer_date", "target_date", "training_dates", "reversed", "status",
    "reason", "training_sha256", "unary_fit_sha256", "reference_count",
    "references_sha256", "references", "unary_probability_sha256", "node_order_sha256",
    "residualizer_sha256", "residualizer_training_sha256", "residualizer_ridge",
    "student_t_nu", "category_order", "state_order", "null_mean_fwd", "null_mean_bwd",
    "null_scale_ff", "null_scale_fb", "null_scale_bb",
], ["conditional_$(c)_mean_$(o)" for c in _GRAPH_CATEGORY_ORDER for o in _GRAPH_OUTPUT_NAMES], [
    "conditional_scale_ff", "conditional_scale_fb", "conditional_scale_bb",
    "selected_start", "selected_alpha", "selected_objective", "selected_converged",
], ["$(n)_$(o)" for n in _GRAPH_COEFFICIENT_NAMES for o in _GRAPH_OUTPUT_NAMES])
const _GRAPH_UNARY_HEADER = [
    "schema_version", "fit_sha256", "model", "file", "date", "lobe", "t_nm",
    "node_identity", "state_order", "p0", "p1", "unary_fit_sha256",
    "unary_probability_sha256", "node_order_sha256",
]
const _GRAPH_EDGE_HEADER = [
    "schema_version", "fit_sha256", "model", "scope", "outer_date", "target_date",
    "file", "date", "left_lobe", "right_lobe", "left_t_nm", "right_t_nm",
    "segment_id", "ordinal", "edge_identity", "edge_row_raw_sha256",
    "residualizer_sha256", "unary_fit_sha256", "pred_fwd", "pred_bwd", "status", "reason",
]

struct _ChainBlocked <: Exception
    code::Symbol
    message::String
end
_blocked(code::Symbol, message) = throw(_ChainBlocked(code, String(message)))

struct _Snapshot
    path::String
    relative::String
    bytes::Vector{UInt8}
    sha256::String
    device::UInt64
    inode::UInt64
    size::Int64
    mode::UInt32
end

struct _Todo10Node
    file::String
    date::String
    lobe::Int
    t_nm::Float64
    t_text::String
    segment_id::String
    status::String
end

struct _Todo10Edge
    file::String
    date::String
    left_lobe::Int
    right_lobe::Int
    left_t_nm::Float64
    right_t_nm::Float64
    segment_id::String
    right_segment_id::String
    status::String
    reason::String
    corr_fwd::Union{Nothing,Float64}
    corr_bwd::Union{Nothing,Float64}
    raw_sha256::String
end

struct _Todo10Data
    nodes::Tuple
    edges::Tuple
    receipt_sha256::String
    provenance_sha256::String
    receipt_hashes::Dict{String,String}
end

struct _ReferenceData
    fit_role::String
    model::String
    scope::String
    outer_date::String
    target_date::String
    fit_sha256::String
end

struct _ModelData
    fit_sha256::String
    model::String
    fit_role::String
    scope::String
    outer_date::String
    target_date::String
    reversed::Bool
    status::String
    reason::String
    training_sha256::String
    unary_fit_sha256::String
    unary_probability_sha256::String
    node_order_sha256::String
    residualizer_sha256::String
    residualizer_training_sha256::String
    residualizer_ridge::Union{Nothing,Float64}
    null_mean::Vector{Float64}
    null_scale::Matrix{Float64}
    conditional_means::Matrix{Float64}
    conditional_scale::Matrix{Float64}
    residual_coefficients::Union{Nothing,Matrix{Float64}}
    training_dates::Tuple
    references::Tuple{Vararg{_ReferenceData}}
end

struct _UnaryData
    fit_sha256::String
    model::String
    file::String
    date::String
    lobe::Int
    t_nm::Float64
    p0::Float64
    p1::Float64
    unary_fit_sha256::String
    probability_sha256::String
    node_order_sha256::String
end

struct _TransformData
    fit_sha256::String
    model::String
    scope::String
    outer_date::String
    target_date::String
    file::String
    date::String
    left_lobe::Int
    right_lobe::Int
    left_t_nm::Float64
    right_t_nm::Float64
    segment_id::String
    ordinal::Int
    edge_identity::String
    raw_sha256::String
    residualizer_sha256::String
    unary_fit_sha256::String
    pred_fwd::Union{Nothing,Float64}
    pred_bwd::Union{Nothing,Float64}
    status::String
    reason::String
end

struct _Todo11Data
    report_dir::String
    receipt_sha256::String
    status::String
    reason::String
    models::Tuple
    unaries::Dict
    transforms::Dict
    provenance_sha256::String
    result_sha256::String
    normalization_sha256::String
end

function _sha(bytes)
    return bytes2hex(SHA.sha256(bytes))
end

function _hash_lines(lines)
    return _sha(codeunits(join(String.(collect(lines)), "\n") * "\n"))
end

function _fmt(value::Float64)
    isfinite(value) || _blocked(:invalid_number, "nonfinite value cannot be canonicalized")
    return @sprintf("%.17g", value)
end

function _sha_text(value, context)
    value isa String || _blocked(:schema_mismatch, "$context is not a string")
    occursin(r"^[0-9a-f]{64}$", value) ||
        _blocked(:invalid_hash, "$context is not a lowercase SHA-256: $value")
    return String(value)
end

function _required(table, key::String, context::String)
    haskey(table, key) || _blocked(:schema_mismatch, "$context.$key is absent")
    return table[key]
end

function _required_string(table, key::String, context::String)
    value = _required(table, key, context)
    value isa String || _blocked(:schema_mismatch, "$context.$key is not a string")
    return String(value)
end

function _required_int(table, key::String, context::String)
    value = _required(table, key, context)
    value isa Integer && !(value isa Bool) ||
        _blocked(:schema_mismatch, "$context.$key is not an integer")
    return Int(value)
end

function _exact_keys(table, expected::Set{String}, context::String)
    Set(String.(collect(keys(table)))) == expected ||
        _blocked(:schema_mismatch, "$context keys differ")
end

function _canonical_root(root::AbstractString)
    supplied = String(root)
    isempty(supplied) && _blocked(:invalid_root, "root is empty")
    isabspath(supplied) || _blocked(:noncanonical_path, "root must be absolute")
    isdir(supplied) || _blocked(:invalid_root, "root is not a directory")
    islink(supplied) && _blocked(:symlink_rejected, "root is a symlink")
    absolute = normpath(abspath(supplied))
    realpath(absolute) == absolute || _blocked(:noncanonical_path, "root is not canonical")
    _check_ancestors(absolute)
    return absolute
end

function _check_ancestors(path::String)
    parts = splitpath(normpath(abspath(path)))
    current = first(parts)
    for part in Iterators.drop(parts, 1)
        current = joinpath(current, part)
        islink(current) && _blocked(:symlink_rejected, "symlinked path ancestor: $current")
    end
end

function _resolve_under_root(root::String, supplied::AbstractString, context::String)
    text = String(supplied)
    isempty(text) && _blocked(:noncanonical_path, "$context path is empty")
    absolute = normpath(isabspath(text) ? text : joinpath(root, text))
    prefix = root == "/" ? "/" : root * "/"
    (absolute == root || startswith(absolute, prefix)) ||
        _blocked(:path_escape, "$context escapes canonical root")
    _check_ancestors(absolute)
    isfile(absolute) || _blocked(:missing_input, "$context is not a regular file")
    islink(absolute) && _blocked(:symlink_rejected, "$context is a symlink")
    realpath(absolute) == absolute || _blocked(:noncanonical_path, "$context is not canonical")
    return absolute
end

function _snapshot(root::String, path::AbstractString, context::String)
    absolute = _resolve_under_root(root, path, context)
    bytes = Vector{UInt8}(read(absolute))
    text = try
        String(copy(bytes))
    catch
        _blocked(:invalid_encoding, "$context is not UTF-8")
    end
    isvalid(text) || _blocked(:invalid_encoding, "$context is not UTF-8")
    occursin('\r', text) && _blocked(:carriage_return, "$context is not LF text")
    endswith(text, '\n') || _blocked(:missing_final_lf, "$context lacks final LF")
    file_stat = stat(absolute)
    return _Snapshot(absolute, relpath(absolute, root), bytes, _sha(bytes),
                     UInt64(file_stat.device), UInt64(file_stat.inode),
                     Int64(file_stat.size), UInt32(file_stat.mode))
end

function _resnapshot(snapshot::_Snapshot)
    isfile(snapshot.path) && !islink(snapshot.path) ||
        _blocked(:input_changed, "snapshot path is no longer a regular file: $(snapshot.relative)")
    _check_ancestors(snapshot.path)
    file_stat = stat(snapshot.path)
    UInt64(file_stat.device) == snapshot.device &&
    UInt64(file_stat.inode) == snapshot.inode &&
    Int64(file_stat.size) == snapshot.size &&
    UInt32(file_stat.mode) == snapshot.mode ||
        _blocked(:input_changed, "snapshot file identity changed: $(snapshot.relative)")
    bytes = Vector{UInt8}(read(snapshot.path))
    bytes == snapshot.bytes || _blocked(:input_changed, "snapshot changed before return: $(snapshot.relative)")
    _sha(bytes) == snapshot.sha256 || _blocked(:input_changed, "snapshot hash changed")
end

function _registry_snapshot(root::String, supplied::AbstractString, context::String,
                            registry::Dict{String,_Snapshot})
    absolute = _resolve_under_root(root, supplied, context)
    relative = relpath(absolute, root)
    haskey(registry, relative) || _blocked(:internal_snapshot, "snapshot was not registered: $context")
    snapshot = registry[relative]
    snapshot.path == absolute || _blocked(:source_collision, "snapshot path alias: $context")
    return snapshot
end

function _register_snapshot!(root::String, supplied::AbstractString, context::String,
                             registry::Dict{String,_Snapshot})
    absolute = _resolve_under_root(root, supplied, context)
    relative = relpath(absolute, root)
    haskey(registry, relative) && _blocked(:source_collision, "duplicate snapshot: $context")
    snapshot = _snapshot(root, supplied, context)
    registry[relative] = snapshot
    return snapshot
end

function _authority_snapshots(root::String)
    fixed = [
        (".omo/plans/structured-label-free-unit-assignment.md", _PLAN_SHA256),
        ("config/unit_assignment_structured_candidate.toml", _CANDIDATE_SHA256),
        ("config/unit_assignment_structured_model.toml", _MODEL_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/DoneClaim.json", _T3_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase3-t3-universe/correction/review/AdversarialVerify.json", _T3_REVIEW_FILE_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/DoneClaim.json", _T10_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/review/AdversarialVerify.json", _T10_REVIEW_FILE_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/correction/DoneClaim.json", _BASE_T11_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/correction/review/AdversarialVerify.json", _BASE_T11_REVIEW_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/DoneClaim.json", _REJECTED_T11_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/review/AdversarialVerify.json", _REJECTED_T11_REVIEW_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/DoneClaim.json", _CORRECTION1_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/claim-correction/ClaimCorrection.json", _CORRECTION1_CLAIM_CORRECTION_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction/review/AdversarialVerify.json", _CORRECTION1_REVIEW_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/DoneClaim.json", _HISTORICAL_T11_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/review/AdversarialVerify.json", _HISTORICAL_T11_REVIEW_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/correction2/canonical-t11-source-bundle.bytes", _HISTORICAL_T11_SOURCE_BUNDLE_SHA256),
        ("test/lib/structured_assignment/edge_admission.jl", _T11_SOURCE_SHA256),
        ("test/evaluate_structured_edge_admission.jl", _T11_CLI_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/DoneClaim.json", _PRIOR_SCALE_T11_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/review/AdversarialVerify.json", _PRIOR_SCALE_T11_REVIEW_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/canonical-t11-source-bundle.bytes", _PRIOR_SCALE_T11_SOURCE_BUNDLE_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/DoneClaim.json", _EFFECTIVE_T11_CLAIM_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/review/AdversarialVerify.json", _EFFECTIVE_T11_REVIEW_SHA256),
        (".omo/evidence/structured-label-free-unit-assignment/provenance-rebind/phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/canonical-t11-source-bundle.bytes", _EFFECTIVE_T11_SOURCE_BUNDLE_SHA256),
    ]
    snapshots = _Snapshot[]
    seen = Set{String}()
    for (relative, expected) in fixed
        snapshot = _snapshot(root, relative, "fixed authority $relative")
        snapshot.sha256 == expected ||
            _blocked(:dependency_hash_mismatch, "fixed authority differs: $relative")
        snapshot.relative in seen && _blocked(:source_collision, "duplicate authority snapshot")
        push!(seen, snapshot.relative)
        push!(snapshots, snapshot)
    end
    return snapshots
end

function _parse_toml(snapshot::_Snapshot, context::String)
    try
        return TOML.parse(String(copy(snapshot.bytes)))
    catch error
        _blocked(:schema_mismatch, "$context TOML parse failed: $(sprint(showerror, error))")
    end
end

function _parse_tsv(snapshot::_Snapshot, expected_header::Vector{String}, context::String)
    text = String(copy(snapshot.bytes))
    lines = split(chomp(text), '\n'; keepempty=true)
    isempty(lines) && _blocked(:empty_input, "$context is empty")
    header = split(first(lines), '\t'; keepempty=true)
    header == expected_header || _blocked(:header_mismatch, "$context header differs")
    rows = Vector{Vector{String}}()
    row_bytes = Vector{Vector{UInt8}}()
    for line in Iterators.drop(lines, 1)
        isempty(line) && _blocked(:schema_mismatch, "$context contains an empty row")
        fields = split(line, '\t'; keepempty=true)
        length(fields) == length(expected_header) ||
            _blocked(:schema_mismatch, "$context row width differs")
        push!(rows, fields)
        push!(row_bytes, Vector{UInt8}(codeunits(line * "\n")))
    end
    return header, rows, row_bytes
end

function _canonical_float(text::String, context::String; allow_na=false)
    allow_na && text == "NA" && return nothing
    value = tryparse(Float64, text)
    value === nothing && _blocked(:invalid_number, "$context is not numeric")
    isfinite(value) || _blocked(:invalid_number, "$context is nonfinite")
    _fmt(value) == text || _blocked(:noncanonical_number, "$context is not canonical")
    return value
end

function _positive_int(text::String, context::String)
    occursin(r"^[1-9][0-9]*$", text) || _blocked(:invalid_key, "$context is not canonical")
    value = tryparse(Int, text)
    value === nothing && _blocked(:invalid_key, "$context is outside Int")
    return value
end

function _nonnegative_int(text::String, context::String)
    occursin(r"^(0|[1-9][0-9]*)$", text) || _blocked(:invalid_key, "$context is not canonical")
    value = tryparse(Int, text)
    value === nothing && _blocked(:invalid_key, "$context is outside Int")
    return value
end

function _date(file::String)
    date_match = match(r"^([0-9]{8})", basename(file))
    date_match === nothing && _blocked(:date_mismatch, "scan filename has no canonical date: $file")
    return date_match.captures[1]
end

function _segment_id(file::String, first_lobe::Int, last_lobe::Int)
    lo, hi = minmax(first_lobe, last_lobe)
    bytes = codeunits("segment-v1\nfile=$(basename(file))\nfirst_lobe=$lo\nlast_lobe=$hi\ngrid_hash=$(_GRID_SHA256)\n")
    return _sha(bytes)[1:16]
end

function _key_hash(nodes)
    lines = ["schema_version\tfile\tlobe"]
    append!(lines, "1\t$(node.file)\t$(node.lobe)" for node in nodes)
    return _hash_lines(lines)
end

function _validate_todo10(root::String, edge_path, node_path, receipt_path,
                          registry::Dict{String,_Snapshot})
    edge = _registry_snapshot(root, edge_path, "Todo 10 edge observations", registry)
    node = _registry_snapshot(root, node_path, "Todo 10 node segments", registry)
    receipt = _registry_snapshot(root, receipt_path, "Todo 10 receipt", registry)
    dirs = Set(dirname(s.path) for s in (edge, node, receipt))
    length(dirs) == 1 || _blocked(:path_ambiguity, "Todo 10 files do not share one directory")
    basename(edge.path) == "edge_observations.tsv" || _blocked(:path_schema, "Todo 10 edge filename differs")
    basename(node.path) == "node_segments.tsv" || _blocked(:path_schema, "Todo 10 node filename differs")
    basename(receipt.path) == "receipt.toml" || _blocked(:path_schema, "Todo 10 receipt filename differs")
    parsed_receipt = _parse_toml(receipt, "Todo 10 receipt")
    _exact_keys(parsed_receipt, _T10_RECEIPT_KEYS, "Todo 10 receipt")
    _required_string(parsed_receipt, "schema", "Todo 10 receipt") ==
        "stmfit-structured-edge-feature-receipt-v1" || _blocked(:schema_mismatch, "Todo 10 schema differs")
    _required_int(parsed_receipt, "schema_version", "Todo 10 receipt") == 1 ||
        _blocked(:schema_mismatch, "Todo 10 schema version differs")
    _required_string(parsed_receipt, "status", "Todo 10 receipt") == "PASS" ||
        _blocked(:upstream_status, "Todo 10 receipt is not PASS")
    _required_string(parsed_receipt, "edge_observations_file", "Todo 10 receipt") == basename(edge.path) ||
        _blocked(:receipt_binding, "Todo 10 edge filename binding differs")
    _required_string(parsed_receipt, "node_segments_file", "Todo 10 receipt") == basename(node.path) ||
        _blocked(:receipt_binding, "Todo 10 node filename binding differs")
    _, edge_rows, edge_row_bytes = _parse_tsv(edge, _T10_EDGE_HEADER, "Todo 10 edge observations")
    _, node_rows, _ = _parse_tsv(node, _T10_NODE_HEADER, "Todo 10 node segments")
    _required_string(parsed_receipt, "edge_header", "Todo 10 receipt") == join(_T10_EDGE_HEADER, '\t') ||
        _blocked(:receipt_binding, "Todo 10 edge header differs")
    _required_string(parsed_receipt, "node_header", "Todo 10 receipt") == join(_T10_NODE_HEADER, '\t') ||
        _blocked(:receipt_binding, "Todo 10 node header differs")
    _required_string(parsed_receipt, "edge_header_sha256", "Todo 10 receipt") ==
        _sha(codeunits(join(_T10_EDGE_HEADER, '\t') * "\n")) || _blocked(:hash_mismatch, "Todo 10 edge header hash differs")
    _required_string(parsed_receipt, "node_header_sha256", "Todo 10 receipt") ==
        _sha(codeunits(join(_T10_NODE_HEADER, '\t') * "\n")) || _blocked(:hash_mismatch, "Todo 10 node header hash differs")
    _required_string(parsed_receipt, "edge_observations_sha256", "Todo 10 receipt") == edge.sha256 || _blocked(:hash_mismatch, "Todo 10 edge table hash differs")
    _required_string(parsed_receipt, "node_segments_sha256", "Todo 10 receipt") == node.sha256 || _blocked(:hash_mismatch, "Todo 10 node table hash differs")
    _required_int(parsed_receipt, "edge_row_count", "Todo 10 receipt") == length(edge_rows) || _blocked(:row_count, "Todo 10 edge row count differs")
    _required_int(parsed_receipt, "node_row_count", "Todo 10 receipt") == length(node_rows) || _blocked(:row_count, "Todo 10 node row count differs")
    _sha_text(_required_string(parsed_receipt, "feature_sha256", "Todo 10 receipt"), "Todo 10 feature hash")
    _required_string(parsed_receipt, "candidate_config_sha256", "Todo 10 receipt") == _CANDIDATE_SHA256 || _blocked(:hash_mismatch, "candidate config binding differs")
    _required_string(parsed_receipt, "model_config_sha256", "Todo 10 receipt") == _MODEL_SHA256 || _blocked(:hash_mismatch, "model config binding differs")
    for field in ("universe_receipt_sha256", "keys_sha256", "forward_patch_sha256", "backward_patch_sha256", "forward_patch_receipt_sha256", "backward_patch_receipt_sha256")
        _sha_text(_required_string(parsed_receipt, field, "Todo 10 receipt"), "Todo 10 $field")
    end
    _required_string(parsed_receipt, "grid_sha256", "Todo 10 receipt") == _GRID_SHA256 || _blocked(:hash_mismatch, "grid binding differs")
    _required_string(parsed_receipt, "forward_producer_sha256", "Todo 10 receipt") == _FORWARD_PRODUCER_SHA256 || _blocked(:hash_mismatch, "forward producer binding differs")
    _required_string(parsed_receipt, "backward_producer_sha256", "Todo 10 receipt") == _BACKWARD_PRODUCER_SHA256 || _blocked(:hash_mismatch, "backward producer binding differs")
    _required_string(parsed_receipt, "source_sha256", "Todo 10 receipt") == _T10_SOURCE_SHA256 || _blocked(:hash_mismatch, "Todo 10 source binding differs")
    _required_string(parsed_receipt, "t2_review_sha256", "Todo 10 receipt") == _T2_REVIEW_SHA256 || _blocked(:hash_mismatch, "Todo 2 review binding differs")
    _required_string(parsed_receipt, "t3_review_sha256", "Todo 10 receipt") == _T3_REVIEW_SHA256 || _blocked(:hash_mismatch, "Todo 3 review binding differs")

    nodes = _Todo10Node[]
    node_keys = Set{Tuple{String,Int}}()
    for (row_index, row) in enumerate(node_rows)
        file, lobe_text, t_text, segment, status = row[1:5]
        occursin(r"^[^/\\]+$", file) || _blocked(:invalid_key, "Todo 10 node file is not a basename")
        lobe = _positive_int(lobe_text, "Todo 10 node $row_index lobe")
        t = _canonical_float(t_text, "Todo 10 node $row_index t")
        occursin(r"^[0-9a-f]{16}$", segment) || _blocked(:segment_mismatch, "Todo 10 node segment is invalid")
        status in ("connected", "isolated") || _blocked(:schema_mismatch, "Todo 10 node status is invalid")
        key = (file, lobe)
        key in node_keys && _blocked(:duplicate_key, "Todo 10 node key repeats")
        push!(node_keys, key)
        _date(file)
        for reason in row[6:7]
            reason in ("start", "end", "eligible", "gap_out_of_range", "missing_forward_key", "missing_backward_key", "nonfinite_fwd", "nonfinite_bwd", "zero_variance_fwd", "zero_variance_bwd", "corr_out_of_range_fwd", "corr_out_of_range_bwd") ||
                _blocked(:schema_mismatch, "Todo 10 node boundary reason is invalid")
        end
        for (field, expected) in ((8, _required_string(parsed_receipt, "feature_sha256", "Todo 10 receipt")), (9, _MODEL_SHA256), (10, _T10_SOURCE_SHA256))
            row[field] == expected || _blocked(:hash_mismatch, "Todo 10 node row binding differs")
        end
        push!(nodes, _Todo10Node(file, _date(file), lobe, t, t_text, segment, status))
    end
    nodes == sort(nodes; by=n -> (n.file, n.t_nm, n.lobe)) || _blocked(:order_mismatch, "Todo 10 node order differs")
    _required_string(parsed_receipt, "keys_sha256", "Todo 10 receipt") == _key_hash(nodes) || _blocked(:hash_mismatch, "Todo 10 key order hash differs")

    by_key = Dict((n.file, n.lobe) => n for n in nodes)
    edge_data = _Todo10Edge[]
    expected_pairs = Tuple{String,Int,Int}[]
    for file in sort!(unique(n.file for n in nodes))
        chain = sort!([n for n in nodes if n.file == file]; by=n -> (n.t_nm, n.lobe))
        seen_segments = Set{String}()
        previous_segment = nothing
        for node in chain
            if previous_segment !== nothing && node.segment_id != previous_segment
                node.segment_id in seen_segments && _blocked(:segment_mismatch, "Todo 10 segment run is noncontiguous")
            end
            push!(seen_segments, node.segment_id)
            previous_segment = node.segment_id
        end
        for i in 1:max(0, length(chain) - 1)
            push!(expected_pairs, (file, chain[i].lobe, chain[i + 1].lobe))
        end
        groups = Dict{String,Vector{_Todo10Node}}()
        for n in chain
            push!(get!(groups, n.segment_id, _Todo10Node[]), n)
        end
        for group in values(groups)
            expected_segment = _segment_id(file, first(group).lobe, last(group).lobe)
            expected_segment == first(group).segment_id || _blocked(:segment_mismatch, "Todo 10 segment digest differs")
        end
    end
    observed_pairs = Tuple{String,Int,Int}[]
    precedence = ("gap_out_of_range", "missing_forward_key", "missing_backward_key", "nonfinite_fwd", "nonfinite_bwd", "zero_variance_fwd", "zero_variance_bwd", "corr_out_of_range_fwd", "corr_out_of_range_bwd")
    for (row_index, row) in enumerate(edge_rows)
        file, left_text, right_text, left_t_text, right_t_text, gap_text = row[1:6]
        left = _positive_int(left_text, "Todo 10 edge $row_index left")
        right = _positive_int(right_text, "Todo 10 edge $row_index right")
        left_t = _canonical_float(left_t_text, "Todo 10 edge $row_index left t")
        right_t = _canonical_float(right_t_text, "Todo 10 edge $row_index right t")
        gap = _canonical_float(gap_text, "Todo 10 edge $row_index gap")
        abs(gap - (right_t - left_t)) <= 0.0 || _blocked(:topology_mismatch, "Todo 10 edge gap differs")
        haskey(by_key, (file, left)) && haskey(by_key, (file, right)) || _blocked(:topology_mismatch, "Todo 10 edge endpoint absent")
        by_key[(file, left)].t_nm == left_t && by_key[(file, right)].t_nm == right_t || _blocked(:topology_mismatch, "Todo 10 edge t differs")
        push!(observed_pairs, (file, left, right))
        left_segment, right_segment, status, reason = row[7:10]
        left_segment == by_key[(file, left)].segment_id && right_segment == by_key[(file, right)].segment_id || _blocked(:segment_mismatch, "Todo 10 edge segment differs")
        status in ("eligible", "split") || _blocked(:schema_mismatch, "Todo 10 edge status invalid")
        if status == "eligible"
            0.35 <= gap <= 0.75 || _blocked(:topology_mismatch, "eligible Todo 10 edge gap is outside frozen range")
            left_segment == right_segment && reason == "none" || _blocked(:topology_mismatch, "eligible Todo 10 edge crosses a segment")
            corr_fwd = _canonical_float(row[11], "Todo 10 edge corr_fwd")
            corr_bwd = _canonical_float(row[12], "Todo 10 edge corr_bwd")
            -1.0 <= corr_fwd <= 1.0 && -1.0 <= corr_bwd <= 1.0 || _blocked(:topology_mismatch, "eligible Todo 10 correlation is outside [-1,1]")
        else
            expected_reason = gap < 0.35 || gap > 0.75 ? "gap_out_of_range" : "other"
            if expected_reason == "gap_out_of_range"
                reason == "gap_out_of_range" || _blocked(:split_reason, "out-of-range gap does not take precedence")
            else
                left_segment != right_segment && reason in precedence && reason != "gap_out_of_range" || _blocked(:split_reason, "Todo 10 split reason or boundary differs")
            end
            left_segment != right_segment || _blocked(:topology_mismatch, "split Todo 10 edge does not cross a segment")
            row[11] == "NA" && row[12] == "NA" || _blocked(:schema_mismatch, "split Todo 10 edge exposes correlation")
            corr_fwd = nothing
            corr_bwd = nothing
        end
        for (index, expected) in ((13, _GRID_SHA256), (16, _required_string(parsed_receipt, "feature_sha256", "Todo 10 receipt")), (17, _MODEL_SHA256), (18, _T10_SOURCE_SHA256))
            row[index] == expected || _blocked(:hash_mismatch, "Todo 10 edge row binding differs")
        end
        row[14] == _required_string(parsed_receipt, "forward_patch_sha256", "Todo 10 receipt") || _blocked(:hash_mismatch, "forward patch row binding differs")
        row[15] == _required_string(parsed_receipt, "backward_patch_sha256", "Todo 10 receipt") || _blocked(:hash_mismatch, "backward patch row binding differs")
        push!(edge_data, _Todo10Edge(file, _date(file), left, right, left_t, right_t, left_segment, right_segment, status, reason, corr_fwd, corr_bwd, _sha(edge_row_bytes[row_index])))
    end
    observed_pairs == expected_pairs || _blocked(:topology_mismatch, "Todo 10 edge pair order differs")
    edge_data == sort(edge_data; by=e -> (e.file, e.left_t_nm, e.right_t_nm, e.left_lobe, e.right_lobe)) || _blocked(:order_mismatch, "Todo 10 edge order differs")
    edge_by_key = Dict((e.file, e.left_lobe, e.right_lobe) => e for e in edge_data)
    for file in unique(n.file for n in nodes)
        chain = sort!([n for n in nodes if n.file == file]; by=n -> (n.t_nm, n.lobe))
        for (i, n) in enumerate(chain)
            left_reason = i == 1 ? "start" : (edge_by_key[(file, chain[i - 1].lobe, n.lobe)].status == "eligible" ? "eligible" : edge_by_key[(file, chain[i - 1].lobe, n.lobe)].reason)
            right_reason = i == length(chain) ? "end" : (edge_by_key[(file, n.lobe, chain[i + 1].lobe)].status == "eligible" ? "eligible" : edge_by_key[(file, n.lobe, chain[i + 1].lobe)].reason)
            row = only(node_rows[j] for j in eachindex(node_rows) if node_rows[j][1] == n.file && node_rows[j][2] == string(n.lobe))
            row[6] == left_reason && row[7] == right_reason || _blocked(:boundary_mismatch, "Todo 10 boundary consequence differs")
            has_eligible_neighbor = (i > 1 && edge_by_key[(file, chain[i - 1].lobe, n.lobe)].status == "eligible") ||
                                     (i < length(chain) && edge_by_key[(file, n.lobe, chain[i + 1].lobe)].status == "eligible")
            row[5] == (has_eligible_neighbor ? "connected" : "isolated") || _blocked(:boundary_mismatch, "Todo 10 node status differs")
        end
    end
    receipt_hashes = Dict{String,String}()
    for field in ("universe_receipt_sha256", "keys_sha256", "forward_patch_sha256",
                  "backward_patch_sha256", "forward_patch_receipt_sha256",
                  "backward_patch_receipt_sha256", "source_sha256",
                  "feature_sha256", "grid_sha256")
        receipt_hashes[field] = _required_string(parsed_receipt, field, "Todo 10 receipt")
    end
    receipt_hashes["edge_observations_sha256"] = edge.sha256
    receipt_hashes["node_segments_sha256"] = node.sha256
    key_hash = receipt_hashes["keys_sha256"]
    provenance = _hash_lines(["todo10-edge=$(edge.sha256)", "todo10-node=$(node.sha256)", "todo10-receipt=$(receipt.sha256)", "keys=$key_hash"])
    return _Todo10Data(Tuple(nodes), Tuple(edge_data), receipt.sha256, provenance, receipt_hashes)
end

_reference_sort_key(reference::_ReferenceData) =
    (reference.model, reference.fit_role == "partition" ? 0 : 1,
     reference.scope, reference.outer_date, reference.target_date)

function _validate_scale_feasibility(scale::Matrix{Float64}, context::String)
    size(scale) == (2, 2) || _blocked(:model_mismatch, "$context scale is not 2 by 2")
    scale[1, 2] == scale[2, 1] || _blocked(:model_mismatch, "$context scale is not symmetric")
    all(isfinite, scale) || _blocked(:model_mismatch, "$context scale is nonfinite")
    decomposition = try
        eigen(Symmetric(scale))
    catch
        _blocked(:model_mismatch, "$context scale eigendecomposition failed")
    end
    all(isfinite, decomposition.values) ||
        _blocked(:model_mismatch, "$context scale eigenvalues are nonfinite")
    eigenvalues = Float64.(decomposition.values)
    all(isfinite, eigenvalues) || _blocked(:model_mismatch, "$context scale eigenvalues are nonfinite")
    all(value -> value > 0.0, eigenvalues) || _blocked(:model_mismatch, "$context scale is not positive definite")
    minimum(eigenvalues) >= _SCALE_EIGENVALUE_FLOOR ||
        _blocked(:model_mismatch, "$context scale minimum eigenvalue is below floor")
    condition = maximum(eigenvalues) / minimum(eigenvalues)
    isfinite(condition) || _blocked(:model_mismatch, "$context scale condition number is nonfinite")
    condition <= _SCALE_CONDITION_CAP ||
        _blocked(:model_mismatch, "$context scale condition number exceeds cap")
    return nothing
end

function _parse_model(row::Vector{String}, columns::Dict{String,Int})
    getv(name) = row[columns[name]]
    _canonical_float(getv("schema_version"), "model schema"; allow_na=false) === nothing && nothing
    getv("schema_version") == "1" || _blocked(:schema_mismatch, "model schema version differs")
    fit_sha = _sha_text(getv("fit_sha256"), "model fit hash")
    role = getv("fit_role")
    role in ("partition", "full_refit") || _blocked(:reference_mismatch, "model fit role invalid")
    status = getv("status")
    status in ("PASS", "SKIPPED", "FAIL") || _blocked(:status_mismatch, "model status invalid")
    reversed = getv("reversed") == "true" ? true : getv("reversed") == "false" ? false : _blocked(:schema_mismatch, "model reversal flag invalid")
    getv("student_t_nu") == "8" || _blocked(:schema_mismatch, "model fixed nu differs")
    getv("category_order") == "00,01,10,11" || _blocked(:schema_mismatch, "model category order differs")
    getv("state_order") == "0,1" || _blocked(:schema_mismatch, "model state order differs")
    training_dates = getv("training_dates") == "NA" || isempty(getv("training_dates")) ? () : Tuple(String.(split(getv("training_dates"), ',')))
    all(occursin(r"^[0-9]{8}$", date) for date in training_dates) || _blocked(:schema_mismatch, "model training date is invalid")
    collect(training_dates) == sort!(unique(collect(training_dates))) || _blocked(:order_mismatch, "model training dates are not canonical")
    references_text = getv("references")
    references = isempty(references_text) ? String[] : String.(split(references_text, ';'))
    length(references) == _positive_int(getv("reference_count"), "model reference count") || _blocked(:reference_mismatch, "model reference count differs")
    length(references) == length(unique(references)) || _blocked(:reference_mismatch, "model reference list repeats an identity")
    _hash_lines([join(split(reference, '|'), '\t') for reference in references]) == getv("references_sha256") || _blocked(:reference_mismatch, "model reference hash differs")
    any(!isempty, references) || _blocked(:reference_mismatch, "model reference list is empty")
    reference_records = _ReferenceData[]
    for reference in references
        fields = split(reference, '|'; keepempty=true)
        length(fields) == 6 || _blocked(:reference_mismatch, "reference field count differs")
        fields[1] in ("partition", "full_refit") || _blocked(:reference_mismatch, "reference role is invalid")
        fields[2] == getv("model") && fields[6] == fit_sha || _blocked(:reference_mismatch, "reference fit/model identity differs")
        fields[1] == "partition" && fields[5] == "NA" && _blocked(:reference_mismatch, "partition reference has no target date")
        fields[1] == "full_refit" && (fields[3] != "full_refit" || fields[4] != "NA" || fields[5] != "NA") && _blocked(:reference_mismatch, "full-refit reference metadata differs")
        push!(reference_records, _ReferenceData(fields[1], fields[2], fields[3], fields[4], fields[5], fields[6]))
    end
    canonical_records = sort(copy(reference_records); by=_reference_sort_key)
    reference_records == canonical_records || _blocked(:reference_mismatch, "reference list is not canonically sorted")
    canonical = first(reference_records)
    canonical.fit_role == role && canonical.model == getv("model") && canonical.scope == getv("scope") && canonical.outer_date == getv("outer_date") && canonical.target_date == getv("target_date") || _blocked(:reference_mismatch, "canonical model metadata differs")
    null_mean = [_canonical_float(getv("null_mean_fwd"), "null mean fwd"; allow_na=true), _canonical_float(getv("null_mean_bwd"), "null mean bwd"; allow_na=true)]
    null_scale = [_canonical_float(getv("null_scale_ff"), "null scale ff"; allow_na=true), _canonical_float(getv("null_scale_fb"), "null scale fb"; allow_na=true), _canonical_float(getv("null_scale_bb"), "null scale bb"; allow_na=true)]
    conditional_means = fill(NaN, 4, 2)
    for (i, category) in enumerate(_GRAPH_CATEGORY_ORDER), (j, output) in enumerate(_GRAPH_OUTPUT_NAMES)
        value = _canonical_float(getv("conditional_$(category)_mean_$(output)"), "conditional mean"; allow_na=true)
        conditional_means[i, j] = value === nothing ? NaN : value
    end
    conditional_scale = [_canonical_float(getv("conditional_scale_ff"), "conditional scale ff"; allow_na=true), _canonical_float(getv("conditional_scale_fb"), "conditional scale fb"; allow_na=true), _canonical_float(getv("conditional_scale_bb"), "conditional scale bb"; allow_na=true)]
    if status == "PASS"
        all(isfinite, null_mean) && all(isfinite, conditional_means) || _blocked(:model_mismatch, "PASS model has nonfinite density parameters")
        conditional_means[2, :] == conditional_means[3, :] || _blocked(:model_mismatch, "mixed conditional means are not tied")
        _validate_scale_feasibility(
            [null_scale[1] null_scale[2]; null_scale[2] null_scale[3]], "null")
        _validate_scale_feasibility(
            [conditional_scale[1] conditional_scale[2];
             conditional_scale[2] conditional_scale[3]], "conditional")
    end
    residualizer_sha = getv("residualizer_sha256")
    residualizer_training = getv("residualizer_training_sha256")
    ridge = _canonical_float(getv("residualizer_ridge"), "residualizer ridge"; allow_na=true)
    coefficients = nothing
    if residualizer_sha != "NA"
        _sha_text(residualizer_sha, "residualizer hash")
        _sha_text(residualizer_training, "residualizer training hash")
        ridge === nothing && _blocked(:model_mismatch, "residualizer ridge is absent")
        coefficients = fill(NaN, length(_GRAPH_COEFFICIENT_NAMES), 2)
        for (i, name) in enumerate(_GRAPH_COEFFICIENT_NAMES), (j, output) in enumerate(_GRAPH_OUTPUT_NAMES)
            coefficients[i, j] = _canonical_float(getv("$(name)_$(output)"), "residualizer coefficient")
        end
    end
    finite_or_nan(value) = value === nothing ? NaN : Float64(value)
    training_sha = _sha_text(getv("training_sha256"), "model training hash")
    return _ModelData(fit_sha, getv("model"), role, getv("scope"), getv("outer_date"), getv("target_date"), reversed, status, getv("reason"), training_sha, _sha_text(getv("unary_fit_sha256"), "unary fit hash"), _sha_text(getv("unary_probability_sha256"), "unary probability hash"), _sha_text(getv("node_order_sha256"), "node order hash"), residualizer_sha == "NA" ? "NA" : _sha_text(residualizer_sha, "residualizer hash"), residualizer_training, ridge, [finite_or_nan(null_mean[1]), finite_or_nan(null_mean[2])], [finite_or_nan(null_scale[1]) finite_or_nan(null_scale[2]); finite_or_nan(null_scale[2]) finite_or_nan(null_scale[3])], conditional_means, [finite_or_nan(conditional_scale[1]) finite_or_nan(conditional_scale[2]); finite_or_nan(conditional_scale[2]) finite_or_nan(conditional_scale[3])], coefficients, training_dates, Tuple(reference_records))
end

function _validate_reference_semantics(model::_ModelData, todo10::_Todo10Data)
    all_dates = sort!(unique([node.date for node in todo10.nodes]))
    training_dates = collect(model.training_dates)
    training_dates == sort!(unique(training_dates)) || _blocked(:reference_mismatch, "model training dates are not unique")
    for reference in model.references
        reference.model == model.model && reference.fit_sha256 == model.fit_sha256 || _blocked(:reference_mismatch, "reference fit/model binding differs")
        if reference.fit_role == "full_refit"
            reference.scope == "full_refit" && reference.outer_date == "NA" && reference.target_date == "NA" || _blocked(:reference_mismatch, "full-refit reference scope differs")
            training_dates == all_dates || _blocked(:reference_mismatch, "full-refit training dates do not cover all Todo 10 dates")
        elseif reference.fit_role == "partition"
            reference.scope in ("outer_inner", "outer_score", "final_lodo") || _blocked(:reference_mismatch, "partition scope is not a T11 scope")
            reference.target_date in all_dates || _blocked(:reference_mismatch, "partition target date is absent")
            reference.target_date in training_dates && _blocked(:reference_mismatch, "partition target leaks into training dates")
            if reference.scope == "final_lodo"
                reference.outer_date == "NA" || _blocked(:reference_mismatch, "final-lodo outer date is not NA")
                training_dates == [date for date in all_dates if date != reference.target_date] ||
                    _blocked(:reference_mismatch, "final-lodo training-date exclusion differs")
            else
                reference.outer_date in all_dates || _blocked(:reference_mismatch, "outer partition date is absent")
                reference.outer_date in training_dates && _blocked(:reference_mismatch, "outer date leaks into training dates")
                reference.scope == "outer_inner" && reference.outer_date == reference.target_date && _blocked(:reference_mismatch, "inner target equals outer date")
                reference.scope == "outer_score" && reference.target_date != reference.outer_date && _blocked(:reference_mismatch, "outer score target differs from outer date")
                if reference.scope == "outer_inner"
                    expected = Tuple(date for date in all_dates if date != reference.outer_date && date != reference.target_date)
                    training_dates == collect(expected) || _blocked(:reference_mismatch, "inner training-date exclusion differs")
                else
                    expected = Tuple(date for date in all_dates if date != reference.outer_date)
                    training_dates == collect(expected) || _blocked(:reference_mismatch, "outer-score training-date exclusion differs")
                end
            end
        else
            _blocked(:reference_mismatch, "reference role is not a T11 role")
        end
    end
    return nothing
end

function _parse_unaries(snapshot::_Snapshot, models, todo10::_Todo10Data)
    _, rows, _ = _parse_tsv(snapshot, _GRAPH_UNARY_HEADER, "Todo 11 unary handoff")
    columns = Dict(name => i for (i, name) in enumerate(_GRAPH_UNARY_HEADER))
    result = Dict{Tuple{String,String,Int},_UnaryData}()
    for row in rows
        getv(name) = row[columns[name]]
        getv("schema_version") == "1" || _blocked(:schema_mismatch, "unary schema differs")
        fit = _sha_text(getv("fit_sha256"), "unary fit hash")
        haskey(models, fit) || _blocked(:reference_mismatch, "unary row references unknown fit")
        p0 = _canonical_float(getv("p0"), "unary p0")
        p1 = _canonical_float(getv("p1"), "unary p1")
        0.0 <= p0 <= 1.0 && 0.0 <= p1 <= 1.0 &&
            (isapprox(p0 + p1, 1.0; atol=1e-12, rtol=0.0) || (p0 == 0.0 && p1 == 0.0)) ||
            _blocked(:unary_mismatch, "unary probabilities are invalid")
        file = getv("file"); lobe = _positive_int(getv("lobe"), "unary lobe"); date = getv("date"); t = _canonical_float(getv("t_nm"), "unary t")
        date == _date(file) || _blocked(:unary_mismatch, "unary date differs from filename")
        getv("state_order") == "0,1" || _blocked(:unary_mismatch, "unary state order differs")
        expected_identity = join((file, date, string(lobe), _fmt(t)), '|')
        getv("node_identity") == expected_identity || _blocked(:unary_mismatch, "unary node identity differs")
        key = (fit, file, lobe)
        haskey(result, key) && _blocked(:duplicate_key, "unary key repeats")
        push!(result, key => _UnaryData(fit, getv("model"), file, date, lobe, t, p0, p1, _sha_text(getv("unary_fit_sha256"), "unary fit hash"), _sha_text(getv("unary_probability_sha256"), "unary probability hash"), _sha_text(getv("node_order_sha256"), "unary node order hash")))
    end
    for (fit, model) in models
        rows_for_fit = sort!([value for (key, value) in result if key[1] == fit]; by=u -> (u.file, u.t_nm, u.lobe))
        length(rows_for_fit) == length(todo10.nodes) || _blocked(:unary_coverage, "unary coverage differs from Todo 10 nodes")
        identities = [(u.file, u.date, u.lobe, u.t_nm) for u in rows_for_fit]
        node_hash = _hash_lines(vcat(["schema=structured-edge-todo10-node-order-v1"], [join((u.file, u.date, string(u.lobe), _fmt(u.t_nm)), '\t') for u in rows_for_fit]))
        all(u.node_order_sha256 == node_hash for u in rows_for_fit) && model.node_order_sha256 == node_hash || _blocked(:unary_mismatch, "unary node order binding differs")
        probability_hash = _matrix_hash([u.p1 for u in rows_for_fit])
        all(u.probability_sha256 == probability_hash for u in rows_for_fit) && model.unary_probability_sha256 == probability_hash || _blocked(:unary_mismatch, "unary probability binding differs")
        all(u.model == model.model && u.unary_fit_sha256 == model.unary_fit_sha256 for u in rows_for_fit) || _blocked(:unary_mismatch, "unary fit binding differs")
    end
    return result
end

function _matrix_hash(values)
    io = IOBuffer()
    println(io, length(values), 'x', 1)
    for value in values
        println(io, _fmt(Float64(value)))
    end
    return _sha(take!(io))
end

function _parse_transforms(snapshot::_Snapshot, models, todo10::_Todo10Data)
    _, rows, _ = _parse_tsv(snapshot, _GRAPH_EDGE_HEADER, "Todo 11 transform handoff")
    columns = Dict(name => i for (i, name) in enumerate(_GRAPH_EDGE_HEADER))
    result = Dict{Tuple{String,String,String,Int,Int},_TransformData}()
    row_keys = Tuple{String,String,String,Int,Int}[]
    edge_by_key = Dict((e.file, e.date, e.left_lobe, e.right_lobe) => e for e in todo10.edges if e.status == "eligible")
    for row in rows
        getv(name) = row[columns[name]]
        getv("schema_version") == "1" || _blocked(:schema_mismatch, "transform schema differs")
        fit = _sha_text(getv("fit_sha256"), "transform fit hash")
        haskey(models, fit) || _blocked(:reference_mismatch, "transform references unknown fit")
        model = models[fit]
        file = getv("file"); date = getv("date"); left = _positive_int(getv("left_lobe"), "transform left lobe"); right = _positive_int(getv("right_lobe"), "transform right lobe")
        left_t = _canonical_float(getv("left_t_nm"), "transform left t"); right_t = _canonical_float(getv("right_t_nm"), "transform right t")
        haskey(edge_by_key, (file, date, left, right)) || _blocked(:transform_coverage, "transform edge is not a Todo 10 eligible edge")
        edge = edge_by_key[(file, date, left, right)]
        edge.left_t_nm == left_t && edge.right_t_nm == right_t || _blocked(:transform_mismatch, "transform endpoint t differs")
        getv("segment_id") == edge.segment_id || _blocked(:transform_mismatch, "transform segment differs")
        getv("edge_identity") == join((file, string(left), string(right), edge.segment_id), '|') || _blocked(:transform_mismatch, "transform identity differs")
        _sha_text(getv("edge_row_raw_sha256"), "transform raw row hash") == edge.raw_sha256 || _blocked(:hash_mismatch, "transform raw row binding differs")
        getv("model") == model.model || _blocked(:reference_mismatch, "transform model binding differs")
        getv("scope") == model.scope && getv("outer_date") == model.outer_date && getv("target_date") == model.target_date || _blocked(:reference_mismatch, "transform must use canonical fit metadata")
        transform_reference = _ReferenceData(model.fit_role, model.model, model.scope, model.outer_date, model.target_date, fit)
        transform_reference in model.references || _blocked(:reference_mismatch, "transform reference binding differs")
        status = getv("status"); reason = getv("reason")
        status in ("eligible", "ineligible", "unavailable") || _blocked(:status_mismatch, "transform status invalid")
        expected_reason = status == "eligible" ? "none" : status == "ineligible" ? "nonfinite_endpoint_predictor" : "residualizer_unavailable"
        reason == expected_reason || _blocked(:status_mismatch, "transform status/reason differs")
        model.status == "PASS" && status == "unavailable" && _blocked(:status_mismatch, "PASS fit contains unavailable transform")
        pred_fwd = _canonical_float(getv("pred_fwd"), "transform pred_fwd"; allow_na=true)
        pred_bwd = _canonical_float(getv("pred_bwd"), "transform pred_bwd"; allow_na=true)
        if status == "eligible"
            pred_fwd !== nothing && pred_bwd !== nothing || _blocked(:transform_mismatch, "eligible transform prediction is absent")
            model.residualizer_sha256 != "NA" || _blocked(:model_mismatch, "eligible transform has no residualizer")
        else
            pred_fwd === nothing && pred_bwd === nothing || _blocked(:transform_mismatch, "noneligible transform exposes prediction")
        end
        residualizer = getv("residualizer_sha256")
        residualizer == model.residualizer_sha256 || _blocked(:hash_mismatch, "transform residualizer binding differs")
        unary_fit = _sha_text(getv("unary_fit_sha256"), "transform unary fit hash")
        unary_fit == model.unary_fit_sha256 || _blocked(:hash_mismatch, "transform unary binding differs")
        key = (fit, file, date, left, right)
        haskey(result, key) && _blocked(:duplicate_key, "transform key repeats")
        parsed = _TransformData(fit, getv("model"), getv("scope"), getv("outer_date"), getv("target_date"), file, date, left, right, left_t, right_t, getv("segment_id"), _positive_int(getv("ordinal"), "transform ordinal"), getv("edge_identity"), getv("edge_row_raw_sha256"), residualizer, unary_fit, pred_fwd, pred_bwd, status, reason)
        result[key] = parsed
        push!(row_keys, key)
    end
    for (fit, model) in models
        rows_for_fit = [value for (key, value) in result if key[1] == fit]
        length(rows_for_fit) == length(edge_by_key) || _blocked(:transform_coverage, "transform edge coverage differs")
        expected_keys = sort!([key for key in keys(result) if key[1] == fit]; by=key -> (key[2], key[3], key[4], key[5]))
        [key for key in row_keys if key[1] == fit] == expected_keys || _blocked(:order_mismatch, "transform row order differs")
    end
    return result
end

function _validate_todo11(root::String, model_path, unary_path, transform_path, receipt_path,
                          registry::Dict{String,_Snapshot}, todo10::_Todo10Data)
    model_snapshot = _registry_snapshot(root, model_path, "Todo 11 fitted edge models", registry)
    unary_snapshot = _registry_snapshot(root, unary_path, "Todo 11 fitted unary nodes", registry)
    transform_snapshot = _registry_snapshot(root, transform_path, "Todo 11 fitted edge transforms", registry)
    receipt_snapshot = _registry_snapshot(root, receipt_path, "Todo 11 admission receipt", registry)
    paths = (model_snapshot, unary_snapshot, transform_snapshot, receipt_snapshot)
    length(Set(dirname(s.path) for s in paths)) == 1 || _blocked(:path_ambiguity, "Todo 11 paths do not share one report directory")
    [basename(s.path) for s in paths] == ["fitted_edge_models.tsv", "fitted_unary_nodes.tsv", "fitted_edge_transforms.tsv", "receipt.toml"] || _blocked(:path_schema, "Todo 11 filename differs")
    report_dir = dirname(model_snapshot.path)
    expected_siblings = sort!(vcat(_REPORT_FILES, _GRAPH_HANDOFF_FILES, ["receipt.toml"]))
    sort!(readdir(report_dir)) == expected_siblings ||
        _blocked(:sibling_set, "Todo 11 sibling report file set differs")
    all(isfile(joinpath(report_dir, name)) && !islink(joinpath(report_dir, name)) for name in vcat(_REPORT_FILES, _GRAPH_HANDOFF_FILES)) || _blocked(:missing_input, "Todo 11 sibling artifact is absent")
    artifact_names = vcat(_REPORT_FILES, _GRAPH_HANDOFF_FILES)
    artifact_snapshots = Dict{String,_Snapshot}()
    for name in artifact_names
        artifact_path = joinpath(report_dir, name)
        relative = relpath(artifact_path, root)
        artifact_snapshots[name] = haskey(registry, relative) ? registry[relative] :
            _register_snapshot!(root, artifact_path, "Todo 11 $name", registry)
    end
    report_headers = Dict{String,Vector{String}}(
        "models.tsv" => _MODEL_RESULT_HEADER,
        "partitions.tsv" => _PARTITION_HEADER,
        "ess.tsv" => _ESS_HEADER,
        "starts.tsv" => _START_HEADER,
        "scores.tsv" => _SCORE_HEADER,
        "shuffle.tsv" => _SHUFFLE_HEADER,
        "bootstrap.tsv" => _BOOTSTRAP_HEADER,
        "fitted_edge_models.tsv" => _GRAPH_MODEL_HEADER,
        "fitted_unary_nodes.tsv" => _GRAPH_UNARY_HEADER,
        "fitted_edge_transforms.tsv" => _GRAPH_EDGE_HEADER,
    )
    report_rows = Dict{String,Vector{Vector{String}}}()
    for (name, header) in report_headers
        _, rows, _ = _parse_tsv(artifact_snapshots[name], header, "Todo 11 $name")
        report_rows[name] = rows
    end
    receipt = _parse_toml(receipt_snapshot, "Todo 11 receipt")
    expected_receipt_keys = Set(["schema", "schema_version", "status", "reason", "result_sha256", "model_count", "bootstrap_count", "shuffle_count", "fixed_nu", "node_priors", "pair_prior", "graph_handoff_schema", "graph_handoff_replay_sha256", "graph_handoff_files", "category_order", "state_order", "coefficient_order", "hashes", "models", "row_counts", "graph_handoff", "artifacts"])
    _exact_keys(receipt, expected_receipt_keys, "Todo 11 receipt")
    _required_string(receipt, "schema", "Todo 11 receipt") == "stmfit-structured-edge-admission-receipt-v2" || _blocked(:schema_mismatch, "Todo 11 receipt schema differs")
    _required_int(receipt, "schema_version", "Todo 11 receipt") == 2 || _blocked(:schema_mismatch, "Todo 11 receipt version differs")
    status = _required_string(receipt, "status", "Todo 11 receipt")
    status == "BLOCKED" && _blocked(:upstream_status, "Todo 11 receipt is BLOCKED")
    status in ("PASS", "SKIPPED", "FAIL") || _blocked(:status_mismatch, "Todo 11 status invalid")
    _sha_text(_required_string(receipt, "result_sha256", "Todo 11 receipt"), "Todo 11 overall result hash")
    _required_int(receipt, "fixed_nu", "Todo 11 receipt") == 8 || _blocked(:schema_mismatch, "Todo 11 nu differs")
    receipt["node_priors"] == [0.5, 0.5] || _blocked(:schema_mismatch, "Todo 11 node prior differs")
    receipt["pair_prior"] == [0.0, 0.0, 0.0, 0.0] || _blocked(:schema_mismatch, "Todo 11 pair prior differs")
    receipt["graph_handoff_schema"] == _GRAPH_HANDOFF_SCHEMA || _blocked(:schema_mismatch, "Todo 11 handoff schema differs")
    receipt["graph_handoff_files"] == _GRAPH_HANDOFF_FILES || _blocked(:schema_mismatch, "Todo 11 handoff file declaration differs")
    receipt["category_order"] == collect(_GRAPH_CATEGORY_ORDER) || _blocked(:schema_mismatch, "Todo 11 category order differs")
    receipt["state_order"] == collect(_GRAPH_STATE_ORDER) || _blocked(:schema_mismatch, "Todo 11 state order differs")
    expected_coefficient_order = ["$(name)_$(output)" for name in _GRAPH_COEFFICIENT_NAMES for output in _GRAPH_OUTPUT_NAMES]
    receipt["coefficient_order"] == expected_coefficient_order || _blocked(:schema_mismatch, "Todo 11 coefficient order differs")
    receipt["hashes"] isa AbstractDict || _blocked(:schema_mismatch, "Todo 11 dependency hashes are not a table")
    _exact_keys(receipt["hashes"], _T11_HASH_KEYS, "Todo 11 dependency hashes")
    all(begin _sha_text(String(value), "Todo 11 dependency hash"); true end for value in values(receipt["hashes"])) || _blocked(:hash_mismatch, "Todo 11 dependency hash is invalid")
    receipt["hashes"]["plan_sha256"] == _PLAN_SHA256 || _blocked(:hash_mismatch, "Todo 11 plan binding differs")
    receipt["hashes"]["t7_review_sha256"] == _T7_REVIEW_SHA256 || _blocked(:hash_mismatch, "Todo 7 review binding differs")
    receipt["hashes"]["t10_review_sha256"] == _T10_REVIEW_FILE_SHA256 || _blocked(:hash_mismatch, "Todo 10 review binding differs")
    receipt["hashes"]["candidate_config_sha256"] == _CANDIDATE_SHA256 || _blocked(:hash_mismatch, "candidate config binding differs")
    receipt["hashes"]["model_config_sha256"] == _MODEL_SHA256 || _blocked(:hash_mismatch, "model config binding differs")
    receipt["hashes"]["universe_receipt_sha256"] == todo10.receipt_hashes["universe_receipt_sha256"] || _blocked(:hash_mismatch, "universe receipt binding differs")
    receipt["hashes"]["universe_keys_sha256"] == todo10.receipt_hashes["keys_sha256"] || _blocked(:hash_mismatch, "universe key binding differs")
    receipt["hashes"]["edge_receipt_sha256"] == todo10.receipt_sha256 || _blocked(:hash_mismatch, "Todo 10 receipt binding differs")
    receipt["hashes"]["edge_observations_sha256"] == todo10.receipt_hashes["edge_observations_sha256"] || _blocked(:hash_mismatch, "Todo 10 edge binding differs")
    receipt["hashes"]["node_segments_sha256"] == todo10.receipt_hashes["node_segments_sha256"] || _blocked(:hash_mismatch, "Todo 10 node binding differs")
    receipt["hashes"]["edge_source_sha256"] == todo10.receipt_hashes["source_sha256"] || _blocked(:hash_mismatch, "Todo 10 source binding differs")
    receipt["hashes"]["edge_feature_sha256"] == todo10.receipt_hashes["feature_sha256"] || _blocked(:hash_mismatch, "Todo 10 feature binding differs")
    receipt["hashes"]["t3_review_sha256"] == _T3_REVIEW_SHA256 || _blocked(:hash_mismatch, "Todo 3 review binding differs")
    receipt["hashes"]["source_sha256"] == _T11_SOURCE_BUNDLE_SHA256 || _blocked(:hash_mismatch, "Todo 11 source bundle binding differs")
    _sha_text(receipt["hashes"]["normalization_sha256"], "Todo 11 normalization hash")
    _sha_text(receipt["graph_handoff_replay_sha256"], "Todo 11 replay hash")
    replay_hash = _graph_handoff_replay_hash(artifact_snapshots)
    receipt_replay = receipt["graph_handoff_replay_sha256"]
    replay_hash == receipt_replay || _blocked(:hash_mismatch, "Todo 11 graph replay bytes differ: $replay_hash != $receipt_replay")
    _exact_keys(receipt["artifacts"], Set(artifact_names), "Todo 11 artifact bindings")
    _exact_keys(receipt["row_counts"], Set(artifact_names), "Todo 11 row counts")
    for name in artifact_names
        _required(receipt["artifacts"], name, "Todo 11 artifacts") == artifact_snapshots[name].sha256 || _blocked(:hash_mismatch, "Todo 11 artifact hash differs: $name")
        _required(receipt["row_counts"], name, "Todo 11 row counts") == length(_parse_lines(artifact_snapshots[name].bytes)) - 1 || _blocked(:row_count, "Todo 11 row count differs: $name")
    end
    graph = receipt["graph_handoff"]
    _exact_keys(graph, Set(["fit_hashes", "full_refit_hashes", "model_row_count", "unary_node_row_count", "edge_transform_row_count", "schemas"]), "Todo 11 graph handoff")
    graph["schemas"] == Dict(name => _GRAPH_HANDOFF_SCHEMA for name in _GRAPH_HANDOFF_FILES) || _blocked(:schema_mismatch, "Todo 11 graph schema map differs")
    graph["model_row_count"] == receipt["row_counts"]["fitted_edge_models.tsv"] || _blocked(:row_count, "Todo 11 model row count differs")
    graph["unary_node_row_count"] == receipt["row_counts"]["fitted_unary_nodes.tsv"] || _blocked(:row_count, "Todo 11 unary row count differs")
    graph["edge_transform_row_count"] == receipt["row_counts"]["fitted_edge_transforms.tsv"] || _blocked(:row_count, "Todo 11 transform row count differs")
    model_header, model_rows, _ = _parse_tsv(artifact_snapshots["fitted_edge_models.tsv"], _GRAPH_MODEL_HEADER, "Todo 11 fitted edge models")
    model_columns = Dict(String(name) => i for (i, name) in enumerate(model_header))
    models = Dict{String,_ModelData}()
    for row in model_rows
        model = _parse_model(row, model_columns)
        haskey(models, model.fit_sha256) && _blocked(:duplicate_key, "Todo 11 fit row repeats")
        models[model.fit_sha256] = model
    end
    for model in values(models)
        _validate_reference_semantics(model, todo10)
    end
    all_dates = sort!(unique(node.date for node in todo10.nodes))
    if length(all_dates) >= _PRODUCER_COMPLETE_DATE_COUNT
        Set(model.model for model in values(models)) == Set(["C1", "C2"]) ||
            _blocked(:reference_mismatch, "Todo 11 model set is not exactly C1/C2")
        observed_logical = Set{Tuple{String,String,String,String,String}}()
        for model in values(models), reference in model.references
            push!(observed_logical, (reference.model, reference.fit_role,
                                     reference.scope, reference.outer_date,
                                     reference.target_date))
        end
        expected_logical = Set{Tuple{String,String,String,String,String}}()
        for model_id in ("C1", "C2")
            for outer in all_dates, target in all_dates
                target == outer && continue
                push!(expected_logical, (model_id, "partition", "outer_inner", outer, target))
            end
            for date in all_dates
                push!(expected_logical, (model_id, "partition", "outer_score", date, date))
                push!(expected_logical, (model_id, "partition", "final_lodo", "NA", date))
            end
            push!(expected_logical, (model_id, "full_refit", "full_refit", "NA", "NA"))
        end
        observed_logical == expected_logical ||
            _blocked(:row_count, "nested T11 logical graph is incomplete")
    end
    fit_models_by_id = Dict{String,Vector{_ModelData}}()
    for model in values(models)
        push!(get!(fit_models_by_id, model.model, _ModelData[]), model)
    end
    overall_rows = Dict{String,Vector{String}}()
    for row in report_rows["models.tsv"]
        length(row) == length(_MODEL_RESULT_HEADER) || _blocked(:schema_mismatch, "Todo 11 model result row width differs")
        haskey(overall_rows, row[2]) && _blocked(:duplicate_key, "Todo 11 overall model row repeats")
        row[1] == "1" || _blocked(:schema_mismatch, "Todo 11 overall model schema differs")
        row[3] in ("PASS", "SKIPPED", "FAIL") || _blocked(:status_mismatch, "Todo 11 overall model status invalid")
        _nonnegative_int(row[5], "Todo 11 outer fold count")
        _sha_text(row[6], "Todo 11 final gate hash")
        _sha_text(row[7], "Todo 11 full-fit hash")
        _sha_text(row[8], "Todo 11 overall result hash")
        overall_rows[row[2]] = row
    end
    receipt_models = receipt["models"]
    receipt_models isa AbstractDict || _blocked(:schema_mismatch, "Todo 11 overall model map is not a table")
    Set(String.(collect(keys(receipt_models)))) == Set(keys(overall_rows)) || _blocked(:reference_mismatch, "Todo 11 overall model map differs")
    Set(keys(overall_rows)) == Set(keys(fit_models_by_id)) || _blocked(:reference_mismatch, "Todo 11 fit model IDs differ")
    length(overall_rows) == _required_int(receipt, "model_count", "Todo 11 receipt") || _blocked(:row_count, "Todo 11 overall model count differs")
    overall_statuses = String[]
    for model_id in sort(collect(keys(overall_rows)))
        row = overall_rows[model_id]
        binding = split(String(receipt_models[model_id]), '|'; keepempty=true)
        length(binding) == 3 || _blocked(:reference_mismatch, "Todo 11 overall model receipt binding width differs")
        row[3] == binding[1] && row[4] == binding[2] && row[8] == binding[3] || _blocked(:reference_mismatch, "Todo 11 overall model receipt binding differs")
        push!(overall_statuses, row[3])
        fit_rows = fit_models_by_id[model_id]
        full_references = _ReferenceData[]
        for fit_model in fit_rows
            append!(full_references, [reference for reference in fit_model.references if reference.fit_role == "full_refit"])
        end
        length(full_references) == 1 || _blocked(:reference_mismatch, "Todo 11 overall model does not have exactly one full refit")
        row[7] == first(full_references).fit_sha256 || _blocked(:reference_mismatch, "Todo 11 overall full-fit binding differs")
    end
    serialized_partitions = _ReferenceData[]
    for model in values(models)
        append!(serialized_partitions, [reference for reference in model.references if reference.fit_role == "partition"])
    end
    partition_rows_seen = Set{_ReferenceData}()
    for row in report_rows["partitions.tsv"]
        length(row) == length(_PARTITION_HEADER) || _blocked(:schema_mismatch, "Todo 11 partition row width differs")
        fit = _sha_text(row[12], "Todo 11 partition fit hash")
        haskey(models, fit) || _blocked(:reference_mismatch, "Todo 11 partition references unknown fit")
        reference = _ReferenceData("partition", row[2], row[3], row[4], row[5], fit)
        reference in serialized_partitions || _blocked(:reference_mismatch, "Todo 11 partition row is not serialized")
        reference in partition_rows_seen && _blocked(:duplicate_key, "Todo 11 partition row repeats")
        push!(partition_rows_seen, reference)
    end
    Set(serialized_partitions) == partition_rows_seen || _blocked(:reference_mismatch, "Todo 11 partition reference coverage differs")
    ess_rows = report_rows["ess.tsv"]
    ess_seen = Set{Tuple{String,String,String,String,String,String}}()
    for row in ess_rows
        _nonnegative_int(row[9], "Todo 11 ESS edge count")
        fit = _sha_text(row[16], "Todo 11 ESS fit hash")
        haskey(models, fit) || _blocked(:reference_mismatch, "Todo 11 ESS references unknown fit")
        marker = (row[2], row[3], row[4], row[5], fit, row[6])
        marker in ess_seen && _blocked(:duplicate_key, "ESS category row repeats")
        push!(ess_seen, marker)
    end
    expected_fit_hashes = sort!(collect(keys(models)))
    graph["fit_hashes"] == expected_fit_hashes || _blocked(:reference_mismatch, "Todo 11 fit hash list differs")
    expected_full_refit_hashes = sort!(unique(collect(reference_fields[6] for model in values(models) for reference_fields in (_reference_fields(reference) for reference in model.references) if reference_fields[1] == "full_refit")))
    graph["full_refit_hashes"] == expected_full_refit_hashes || _blocked(:reference_mismatch, "Todo 11 full-refit hash list differs")
    unaries = _parse_unaries(artifact_snapshots["fitted_unary_nodes.tsv"], models, todo10)
    transforms = _parse_transforms(artifact_snapshots["fitted_edge_transforms.tsv"], models, todo10)
    result_hash, normalization_hash, _, _, _ = _validate_diagnostic_semantics!(report_rows, models, todo10, overall_rows, receipt, unaries, transforms)
    replay_value = receipt["graph_handoff_replay_sha256"]
    provenance = _hash_lines(vcat(["todo11-receipt=$(receipt_snapshot.sha256)", "graph=$replay_value"], ["$name=$(artifact_snapshots[name].sha256)" for name in sort(collect(keys(artifact_snapshots)))]))
    return _Todo11Data(report_dir, receipt_snapshot.sha256, status, _required_string(receipt, "reason", "Todo 11 receipt"), Tuple(sort!(collect(values(models)); by=model -> model.fit_sha256)), unaries, transforms, provenance, result_hash, normalization_hash)
end

function _parse_lines(bytes::Vector{UInt8})
    text = String(copy(bytes))
    endswith(text, '\n') || _blocked(:missing_final_lf, "artifact lacks final LF")
    return split(chomp(text), '\n'; keepempty=true)
end

function _graph_handoff_replay_hash(artifacts::Dict{String,_Snapshot})
    io = IOBuffer()
    write(io, codeunits("schema=structured-edge-graph-handoff-replay-v1\n"))
    for name in _GRAPH_HANDOFF_FILES
        snapshot = get(artifacts, name, nothing)
        snapshot isa _Snapshot || _blocked(:missing_input, "graph handoff artifact is absent: $name")
        write(io, codeunits("file=$(name)\nlength=$(length(snapshot.bytes))\n"))
        write(io, snapshot.bytes)
        write(io, UInt8('\n'))
    end
    return _sha(take!(io))
end

function _diag_bool(text::String, context::String)
    text in ("true", "false") || _blocked(:schema_mismatch, "$context is not a canonical boolean")
    return text == "true"
end

function _diag_dates(text::String, context::String)
    isempty(text) && _blocked(:schema_mismatch, "$context is empty")
    dates = String.(split(text, ','))
    all(occursin(r"^[0-9]{8}$", date) for date in dates) || _blocked(:schema_mismatch, "$context contains an invalid date")
    dates == sort!(unique(dates)) || _blocked(:order_mismatch, "$context is not canonical")
    return Tuple(dates)
end

function _diag_count(text::String, context::String)
    return _nonnegative_int(text, context)
end

function _diag_out_float(text::String, context::String; allow_na=true)
    return _canonical_float(text, context; allow_na=allow_na)
end

function _diag_type7(values::Vector{Float64}, probability::Float64)
    isempty(values) && _blocked(:row_count, "Type-7 quantile has no values")
    ordered = sort(copy(values))
    h = 1.0 + (length(ordered) - 1) * probability
    lower = floor(Int, h)
    upper = ceil(Int, h)
    lower == upper && return ordered[lower]
    return ordered[lower] + (h - lower) * (ordered[upper] - ordered[lower])
end

function _diag_bootstrap_values(scan_scores::Dict{String,Dict{String,Float64}})
    isempty(scan_scores) && _blocked(:row_count, "bootstrap scan scores are empty")
    values = Float64[]
    for seed in 0:499
        random = Random.MersenneTwister(seed)
        date_values = Float64[]
        for date in sort!(collect(keys(scan_scores)))
            scans = scan_scores[date]
            names = sort!(collect(keys(scans)))
            isempty(names) && _blocked(:row_count, "bootstrap date has no scans")
            sampled = [scans[names[rand(random, eachindex(names))]] for _ in eachindex(names)]
            all(isfinite, sampled) || _blocked(:invalid_number, "bootstrap sample is nonfinite")
            push!(date_values, sum(sampled) / length(sampled))
        end
        push!(values, sum(date_values) / length(date_values))
    end
    return values
end

_diag_number_text(value::Float64) = isfinite(value) ? _fmt(value) : string(value)

function _diag_invalid_score_hash(status::String, reason::String, fit::String, target::String)
    return _hash_lines([
        "schema=structured-edge-partition-score-v1",
        "status=$(status)",
        "reason=$(reason)",
        "fit=$(fit)",
        "target=$(target)",
    ])
end

function _diag_score_hash(fit::String, target::String, date_mean::Float64,
                          control_mean::Float64, shuffle_quantile::Float64,
                          shuffle_pass::Bool, scan_scores, gains,
                          shuffle_values)
    lines = [
        "schema=structured-edge-partition-score-v1",
        "fit=$(fit)",
        "target=$(target)",
        "date_mean=$(_diag_number_text(date_mean))",
        "control_mean=$(_diag_number_text(control_mean))",
        "shuffle_quantile=$(_diag_number_text(shuffle_quantile))",
        "shuffle_pass=$(shuffle_pass)",
    ]
    append!(lines, "scan=$(key):$(_diag_number_text(scan_scores[key]))" for key in sort!(collect(keys(scan_scores))))
    append!(lines, "gain=$(key):$(_diag_number_text(gains[key]))" for key in sort!(collect(keys(gains))))
    append!(lines, "shuffle=$(_diag_number_text(value))" for value in shuffle_values)
    return _hash_lines(lines)
end

function _diag_gate_hash(status::String, reason::String, dates::Vector{String},
                         date_means, bootstrap_values::Vector{Float64},
                         lower::Float64, every_positive::Bool,
                         every_shuffle::Bool, reversal::Bool)
    if status in ("BLOCKED", "FAIL", "SKIPPED") && isempty(bootstrap_values)
        return _hash_lines([
            "schema=structured-edge-gate-v1",
            "status=$(status)",
            "reason=$(reason)",
            "dates=$(join(dates, ','))",
        ])
    end
    lines = [
        "schema=structured-edge-gate-v1",
        "status=$(status)",
        "reason=$(reason)",
        "dates=$(join(dates, ','))",
        "bootstrap_lower=$(_diag_number_text(lower))",
        "every_positive=$(every_positive)",
        "every_shuffle=$(every_shuffle)",
        "reversal=$(reversal)",
    ]
    for date in dates
        push!(lines, "date=$(date):$(_diag_number_text(date_means[date]))")
    end
    append!(lines, "bootstrap=$(_diag_number_text(value))" for value in bootstrap_values)
    return _hash_lines(lines)
end

function _diag_edge_identity(edge::_Todo10Edge)
    return string(edge.file, '\t', edge.left_lobe, '\t', edge.right_lobe, '\t', edge.segment_id)
end

function _diag_raw_heldout_hash(todo10::_Todo10Data, target_date::String,
                                fit::String, transforms)
    lines = String[]
    for edge in todo10.edges
        edge.date == target_date && edge.status == "eligible" || continue
        transform = get(transforms, (fit, edge.file, edge.date, edge.left_lobe, edge.right_lobe), nothing)
        transform !== nothing && transform.status == "eligible" || continue
        push!(lines, string(_diag_edge_identity(edge), '=', edge.raw_sha256))
    end
    return isempty(lines) ? _sha(UInt8[]) : _hash_lines(lines)
end

function _diag_pair_weights(left::Float64, right::Float64)
    values = _diag_pair_weights_raw(left, right)
    abs(sum(values) - 1.0) <= 1e-12 || _blocked(:support_mismatch, "unary pair weights do not sum to one")
    return values
end

_diag_pair_weights_raw(left::Float64, right::Float64) =
    ((1.0 - left) * (1.0 - right), (1.0 - left) * right,
     left * (1.0 - right), left * right)

function _diag_eligible_edges(todo10::_Todo10Data, fit::String, transforms)
    result = _Todo10Edge[]
    for edge in todo10.edges
        edge.status == "eligible" || continue
        transform = get(transforms, (fit, edge.file, edge.date,
                                     edge.left_lobe, edge.right_lobe), nothing)
        transform !== nothing && transform.status == "eligible" || continue
        push!(result, edge)
    end
    return result
end

function _diag_edge_gain(model::_ModelData, edge::_Todo10Edge,
                         transform::_TransformData, unaries)
    transform.status == "eligible" || _blocked(:support_mismatch, "gain uses an ineligible transform")
    left = unaries[(model.fit_sha256, edge.file, edge.left_lobe)].p1
    right = unaries[(model.fit_sha256, edge.file, edge.right_lobe)].p1
    weights = model.reversed ? _diag_pair_weights(right, left) : _diag_pair_weights(left, right)
    raw = [edge.corr_fwd - transform.pred_fwd, edge.corr_bwd - transform.pred_bwd]
    null_density = _student_log_density(raw, model.null_mean, model.null_scale)
    values = 0.0
    values += weights[1] * (_student_log_density(raw, vec(model.conditional_means[1, :]), model.conditional_scale) - null_density)
    values += (weights[2] + weights[3]) * (_student_log_density(raw, vec(model.conditional_means[2, :]), model.conditional_scale) - null_density)
    values += weights[4] * (_student_log_density(raw, vec(model.conditional_means[4, :]), model.conditional_scale) - null_density)
    isfinite(values) || _infer_fail(:score_failure, "reconstructed edge gain is nonfinite")
    return values
end

function _diag_support_evidence(model::_ModelData, todo10::_Todo10Data,
                                unaries, transforms)
    training_edges = [edge for edge in _diag_eligible_edges(todo10, model.fit_sha256, transforms)
                      if edge.date in model.training_dates]
    weights = [_diag_pair_weights(
        unaries[(model.fit_sha256, edge.file, edge.left_lobe)].p1,
        unaries[(model.fit_sha256, edge.file, edge.right_lobe)].p1,
    ) for edge in training_edges]
    ess = [isempty(weights) ? 0.0 :
           sum(weight[category] for weight in weights)^2 /
           sum(weight[category]^2 for weight in weights) for category in 1:4]
    minimum_support = Int[]
    for category in 1:4
        counts = [length(Set(edge.file for edge in training_edges
                             if edge.date == date &&
                                _diag_pair_weights(
                                    unaries[(model.fit_sha256, edge.file, edge.left_lobe)].p1,
                                    unaries[(model.fit_sha256, edge.file, edge.right_lobe)].p1,
                                )[category] > 0.0))
                  for date in model.training_dates]
        push!(minimum_support, isempty(counts) ? 0 : Base.minimum(counts))
    end
    dates = isempty(training_edges) ? 0 : length(model.training_dates)
    scans = length(unique(edge.file for edge in training_edges))
    edges = length(training_edges)
    conditions = (
        dates >= _PRODUCER_MINIMUM_DATES,
        scans >= 20,
        edges >= 90,
        all(value -> value >= 10.0, ess),
        all(value -> value >= 3, minimum_support),
        edges / 9.0 > 10.0,
    )
    reason = !conditions[1] ? "insufficient_dates" :
             !conditions[2] ? "insufficient_scans" :
             !conditions[3] ? "insufficient_edges" :
             !conditions[4] ? "insufficient_category_ess" :
             !conditions[5] ? "insufficient_scan_support" :
             !conditions[6] ? "insufficient_parameter_ratio" : "ok"
    return (dates=dates, scans=scans, edges=edges, ess=ess,
            minimum=minimum_support, free=9, ratio=edges / 9.0,
            sufficient=all(conditions), reason=reason)
end

function _diag_float(text::String, context::String; allow_na=false)
    value = _canonical_float(text, context; allow_na=allow_na)
    return value
end

function _validate_diagnostic_semantics!(report_rows, models, todo10::_Todo10Data,
                                         overall_rows, receipt, unaries, transforms)
    references = Dict{Tuple{String,String,String,String,String},Tuple{_ReferenceData,_ModelData}}()
    for model in values(models)
        for reference in model.references
            key = (reference.model, reference.scope, reference.outer_date, reference.target_date, reference.fit_sha256)
            haskey(references, key) && _blocked(:duplicate_key, "serialized reference repeats")
            references[key] = (reference, model)
        end
    end
    all_dates = sort!(unique(node.date for node in todo10.nodes))
    legacy_contract = length(models) == 1 && all(length(model.references) == 1 &&
                          first(model.references).fit_role == "full_refit"
                          for model in values(models)) ||
                     (length(all_dates) < _PRODUCER_COMPLETE_DATE_COUNT &&
                      Set(model.model for model in values(models)) == Set(["C1"]) &&
                      any(reference.fit_role == "partition"
                          for model in values(models) for reference in model.references))
    if length(all_dates) < _PRODUCER_COMPLETE_DATE_COUNT && !legacy_contract
        Set(keys(overall_rows)) == Set(["C1", "C2"]) ||
            _blocked(:reference_mismatch, "early admission overall model IDs differ")
        all(row[3] == "SKIPPED" && row[4] == "insufficient_dates" &&
            row[5] == "0" for row in values(overall_rows)) ||
            _blocked(:status_mismatch, "four-date PASS or nested publication is invalid")
        isempty(report_rows["partitions.tsv"]) ||
            _blocked(:reference_mismatch, "early admission exposes partition rows")
    end
    logical_references = Set((reference.model, reference.fit_role, reference.scope,
                             reference.outer_date, reference.target_date)
                            for (reference, _) in values(references))
    model_ids = Set(model.model for model in values(models))
    model_ids == Set(["C1", "C2"]) ||
        (length(all_dates) < _PRODUCER_COMPLETE_DATE_COUNT &&
         all(model.model == "C1" for model in values(models)) ? nothing :
         _blocked(:reference_mismatch, "Todo 11 model set is not exactly C1/C2"))
    if length(all_dates) < _PRODUCER_COMPLETE_DATE_COUNT
        # A producer early-admission report has exactly one sparse full-refit
        # publication per model and no nested partitions.  Keep the historical
        # one-model/full-refit compatibility fixture isolated from this public
        # branch; it is not an authoritative admission report.
        if !legacy_contract
            length(models) == 2 || _blocked(:reference_mismatch, "early admission model count differs")
            all(model.model in ("C1", "C2") for model in values(models)) ||
                _blocked(:reference_mismatch, "early admission model IDs differ")
            all(model.status == "SKIPPED" && model.reason == "insufficient_dates" &&
                length(model.references) == 1 &&
                first(model.references).fit_role == "full_refit" for model in values(models)) ||
                _blocked(:status_mismatch, "early admission SKIPPED shape differs")
            isempty(logical_references) && _blocked(:reference_mismatch, "early admission has partition references")
        end
    elseif model_ids != Set(["C1", "C2"])
        _blocked(:reference_mismatch, "complete admission model IDs differ")
    end
    if !isempty(all_dates) && length(all_dates) >= _PRODUCER_COMPLETE_DATE_COUNT
        expected_logical = Set{Tuple{String,String,String,String,String}}()
        for model_data in values(models)
            model = model_data.model
            for outer in all_dates, target in all_dates
                target == outer && continue
                push!(expected_logical, (model, "partition", "outer_inner", outer, target))
            end
            for date in all_dates
                push!(expected_logical, (model, "partition", "outer_score", date, date))
                push!(expected_logical, (model, "partition", "final_lodo", "NA", date))
            end
            push!(expected_logical, (model, "full_refit", "full_refit", "NA", "NA"))
        end
        logical_references == expected_logical || _blocked(:row_count, "nested T11 logical graph is incomplete")
        for model_id in ("C1", "C2")
            model_references = [reference for (reference, _) in values(references)
                               if reference.model == model_id]
            count(reference -> reference.scope == "outer_inner", model_references) == 42 ||
                _blocked(:row_count, "outer-inner reference count differs")
            count(reference -> reference.scope == "outer_score", model_references) == 7 ||
                _blocked(:row_count, "outer-score reference count differs")
            count(reference -> reference.scope == "final_lodo", model_references) == 7 ||
                _blocked(:row_count, "final-lodo reference count differs")
            count(reference -> reference.fit_role == "full_refit", model_references) == 1 ||
                _blocked(:row_count, "full-refit reference count differs")
            length(unique(reference.fit_sha256 for reference in model_references)) == 29 ||
                _blocked(:row_count, "unique fit-row count differs")
            ab = only(reference.fit_sha256 for reference in model_references
                      if reference.scope == "outer_inner" &&
                         reference.outer_date == all_dates[1] &&
                         reference.target_date == all_dates[2])
            ba = only(reference.fit_sha256 for reference in model_references
                      if reference.scope == "outer_inner" &&
                         reference.outer_date == all_dates[2] &&
                         reference.target_date == all_dates[1])
            ab == ba || _blocked(:reference_mismatch, "A/B shared five-date fit differs")
            for date in all_dates
                outer = only(reference.fit_sha256 for reference in model_references
                             if reference.scope == "outer_score" && reference.outer_date == date)
                final = only(reference.fit_sha256 for reference in model_references
                             if reference.scope == "final_lodo" && reference.target_date == date)
                outer == final || _blocked(:reference_mismatch, "outer/final six-date fit alias differs")
            end
            full = only(reference.fit_sha256 for reference in model_references
                         if reference.fit_role == "full_refit")
            full in [reference.fit_sha256 for reference in model_references
                     if reference.scope == "outer_score"] &&
                _blocked(:reference_mismatch, "full refit is not distinct")
        end
    end
    partition_rows = Dict{Tuple{String,String,String,String,String},Vector{String}}()
    for row in report_rows["partitions.tsv"]
        row[1] == "1" || _blocked(:schema_mismatch, "partition schema version differs")
        fit = _sha_text(row[12], "partition fit hash")
        key = (row[2], row[3], row[4], row[5], fit)
        haskey(references, key) && references[key][1].fit_role == "partition" || _blocked(:reference_mismatch, "partition row reference differs")
        haskey(partition_rows, key) && _blocked(:duplicate_key, "partition row repeats")
        model = references[key][2]
        _diag_dates(row[8], "partition training dates") == model.training_dates || _blocked(:reference_mismatch, "partition training dates differ")
        row[11] == model.training_sha256 || _blocked(:hash_mismatch, "partition training hash differs")
        row[6] in ("PASS", "SKIPPED", "FAIL", "BLOCKED") || _blocked(:status_mismatch, "partition status invalid")
        _nonnegative_int(row[9], "partition training scan count")
        _nonnegative_int(row[10], "partition held-out scan count")
        for (index, field) in ((13, "partition score hash"), (14, "partition raw held-out hash"), (15, "partition raw-at-score hash"))
            _sha_text(row[index], field)
        end
        expected_training_scans = length(unique(node.file for node in todo10.nodes if node.date in model.training_dates))
        expected_heldout_scans = row[6] == "PASS" ?
            length(unique(edge.file for edge in _diag_eligible_edges(todo10, fit, transforms)
                       if edge.date == row[5])) : 0
        parse(Int, row[9]) == expected_training_scans || _blocked(:reference_mismatch, "partition training scan count differs")
        parse(Int, row[10]) == expected_heldout_scans || _blocked(:reference_mismatch, "partition held-out scan count differs")
        raw_hash = _diag_raw_heldout_hash(todo10, row[5], fit, transforms)
        row[14] == raw_hash && row[15] == raw_hash || _blocked(:hash_mismatch, "partition raw held-out hashes differ")
        _diag_bool(row[16], "partition reversal flag")
        _diag_bool(row[17], "partition shuffle flag")
        partition_rows[key] = row
    end
    expected_partition_keys = Set(key for key in keys(references) if references[key][1].fit_role == "partition")
    Set(keys(partition_rows)) == expected_partition_keys || _blocked(:row_count, "partition diagnostic coverage differs")

    ess_rows = Dict{Tuple{String,String,String,String,String},Vector{Vector{String}}}()
    ess_edge_counts = Dict{String,Int}()
    for row in report_rows["ess.tsv"]
        row[1] == "1" || _blocked(:schema_mismatch, "ESS schema version differs")
        fit = _sha_text(row[16], "ESS fit hash")
        key = (row[2], row[3], row[4], row[5], fit)
        haskey(references, key) || _blocked(:reference_mismatch, "ESS reference differs")
        dates_count = _diag_count(row[7], "ESS date count")
        scans_count = _diag_count(row[8], "ESS scan count")
        edge_count = _diag_count(row[9], "ESS edge count")
        category = _nonnegative_int(row[6], "ESS category")
        category in 0:3 || _blocked(:schema_mismatch, "ESS category outside 0:3")
        _diag_out_float(row[10], "ESS Kish ESS"; allow_na=true)
        _diag_count(row[11], "ESS minimum support scans")
        _diag_count(row[12], "ESS free parameters")
        _diag_out_float(row[13], "ESS observation ratio"; allow_na=true)
        _diag_bool(row[14], "ESS sufficient")
        isempty(row[15]) && _blocked(:schema_mismatch, "ESS reason is empty")
        model = references[key][2]
        if !legacy_contract
            support = _diag_support_evidence(model, todo10, unaries, transforms)
            dates_count == support.dates && scans_count == support.scans && edge_count == support.edges ||
                _blocked(:reference_mismatch, "ESS support counts differ from eligible samples")
            _nonnegative_int(row[12], "ESS free parameters") == support.free ||
                _blocked(:support_mismatch, "ESS free parameter count differs")
            ratio = _diag_out_float(row[13], "ESS observation ratio"; allow_na=true)
            ratio !== nothing && isapprox(ratio, support.ratio; atol=1e-12, rtol=1e-12) ||
                _blocked(:support_mismatch, "ESS observation ratio differs")
            value = _diag_out_float(row[10], "ESS Kish ESS"; allow_na=true)
            value !== nothing && isapprox(value, support.ess[category + 1]; atol=1e-10, rtol=1e-10) ||
                _blocked(:support_mismatch, "ESS Kish value differs from eligible weights")
            _diag_count(row[11], "ESS minimum support scans") == support.minimum[category + 1] ||
                _blocked(:support_mismatch, "ESS minimum scan support differs")
            row[14] == string(support.sufficient) ||
                _blocked(:status_mismatch, "ESS sufficiency differs from support evidence")
            row[15] == support.reason ||
                _blocked(:status_mismatch, "ESS first-failing reason differs")
        end
        if legacy_contract
            if model.status == "PASS"
                row[14] == "true" && row[15] == "ok" ||
                    _blocked(:status_mismatch, "legacy PASS ESS sufficiency/reason differs")
            end
        end
        push!(get!(ess_rows, key, Vector{Vector{String}}()), row)
        haskey(ess_edge_counts, fit) && ess_edge_counts[fit] != edge_count &&
            _blocked(:row_count, "ESS edge counts differ across references")
        ess_edge_counts[fit] = edge_count
    end
    expected_diag_keys = Set(keys(references))
    for key in expected_diag_keys
        rows = get(ess_rows, key, Vector{Vector{String}}())
        Set(_nonnegative_int(row[6], "ESS category") for row in rows) == Set(0:3) || _blocked(:row_count, "ESS category coverage differs")
        isempty(rows) && continue
        baseline = rows[1][7:9]
        all(row[7:9] == baseline && row[12:15] == rows[1][12:15] for row in rows) ||
            _blocked(:reference_mismatch, "ESS invariant support fields disagree")
    end
    for (fit, model) in models
        model.residualizer_sha256 == "NA" && continue
        haskey(ess_edge_counts, fit) || _blocked(:row_count, "residualizer lacks ESS evidence")
        _residualizer_sha(model, ess_edge_counts[fit]) == model.residualizer_sha256 ||
            _blocked(:hash_mismatch, "diagnostic residualizer digest differs")
    end

    start_seen = Set{Tuple{String,String,String,String,String,String,String}}()
    start_rows_by_key = Dict{Tuple{String,String,String,String,String},Vector{Vector{String}}}()
    graph_model_columns = Dict(name => index for (index, name) in enumerate(_GRAPH_MODEL_HEADER))
    graph_model_rows = Dict(row[3] => row for row in report_rows["fitted_edge_models.tsv"])
    for row in report_rows["starts.tsv"]
        row[1] == "1" || _blocked(:schema_mismatch, "start schema version differs")
        fit = _sha_text(row[17], "start fit hash")
        key = (row[2], row[3], row[4], row[5], "$(fit)")
        haskey(references, key) || _blocked(:reference_mismatch, "start reference differs")
        role = row[6]
        role in ("null", "conditional") || _blocked(:schema_mismatch, "start role differs")
        index = _nonnegative_int(row[7], "start index")
        if role == "null"
            index == 0 && row[8] == "NA" || _blocked(:schema_mismatch, "null start metadata differs")
        else
            index in 1:5 || _blocked(:schema_mismatch, "conditional start index differs")
            expected_alpha = ("0", "0.25", "0.5", "0.75", "1")[index]
            row[8] == expected_alpha || _blocked(:schema_mismatch, "conditional start alpha differs")
        end
        row[9] in ("PASS", "SKIPPED", "FAIL", "BLOCKED") || _blocked(:status_mismatch, "start status invalid")
        isempty(row[10]) && _blocked(:schema_mismatch, "start reason is empty")
        row[9] == "PASS" && row[10] != "ok" &&
            _blocked(:status_mismatch, "PASS start reason differs")
        _nonnegative_int(row[11], "start iterations")
        _diag_out_float(row[12], "start initial objective"; allow_na=true)
        _diag_out_float(row[13], "start final objective"; allow_na=true)
        _nonnegative_int(row[14], "start rejected iteration")
        _diag_out_float(row[15], "start rejected objective"; allow_na=true)
        _sha_text(row[16], "start trace hash")
        reference_key = (row[2], row[3], row[4], row[5], fit)
        push!(get!(start_rows_by_key, reference_key, Vector{Vector{String}}()), row)
        marker = (row[2], row[3], row[4], row[5], fit, role, string(index))
        marker in start_seen && _blocked(:duplicate_key, "start row repeats")
        push!(start_seen, marker)
    end
    expected_start_all = Set{Tuple{String,String,String,String,String,String,String}}()
    for key in expected_diag_keys
        model = references[key][2]
        fit = key[5]
        rows = get(start_rows_by_key, key, Vector{Vector{String}}())
        expected = Set{Tuple{String,String,String,String,String,String,String}}()
        if !isempty(rows)
            null_rows = [row for row in rows if row[6] == "null"]
            length(null_rows) == 1 || _blocked(:row_count, "null start coverage differs")
            conditional_indices = sort([_nonnegative_int(row[7], "conditional start index") for row in rows if row[6] == "conditional"])
            conditional_indices == collect(1:length(conditional_indices)) || _blocked(:order_mismatch, "conditional start coverage differs")
            for (role, indices) in (("null", [0]), ("conditional", conditional_indices))
                for index in indices
                    push!(expected, (key[1], key[2], key[3], key[4], fit, role, string(index)))
                end
            end
        end
        graph_row = get(graph_model_rows, fit, nothing)
        graph_row === nothing && _blocked(:reference_mismatch, "start fit is absent from graph model rows")
        selected_fields = [graph_row[graph_model_columns[name]] for name in ("selected_start", "selected_alpha", "selected_objective", "selected_converged")]
        null_published = graph_row[graph_model_columns["null_mean_fwd"]] != "NA"
        conditional_published = graph_row[graph_model_columns["conditional_00_mean_corr_fwd"]] != "NA"
        if null_published
            any(row[6] == "null" for row in rows) || _blocked(:row_count, "published null start is absent")
        else
            any(row[6] == "null" for row in rows) && _blocked(:row_count, "unpublished null start is present")
        end
        if conditional_published
            conditional_indices = sort([_nonnegative_int(row[7], "conditional start index")
                                        for row in rows if row[6] == "conditional"])
            conditional_indices == collect(1:5) || _blocked(:row_count, "published conditional start coverage differs")
        else
            any(row[6] == "conditional" for row in rows) &&
                _blocked(:row_count, "unpublished conditional start is present")
        end
        if !conditional_published
            all(field == "NA" for field in selected_fields) || _blocked(:reference_mismatch, "graph model exposes an un-emitted selected start")
        elseif model.status != "PASS" && all(field == "NA" for field in selected_fields)
            # A failed/skipped conditional stage may retain numeric means,
            # certified scale fields and per-start diagnostics without claiming
            # that any start was selected.
            nothing
        else
            selected_fields[1] == "NA" && _blocked(:reference_mismatch, "conditional starts lack selected start metadata")
            selected_index = _positive_int(selected_fields[1], "selected start index")
            selected_row = only(row for row in rows if row[6] == "conditional" && parse(Int, row[7]) == selected_index)
            selected_fields[2] == selected_row[8] || _blocked(:reference_mismatch, "selected start alpha differs")
            selected_fields[3] == selected_row[13] || _blocked(:reference_mismatch, "selected start objective differs")
            selected_fields[4] == (selected_row[9] == "PASS" ? "true" : "false") ||
                _blocked(:reference_mismatch, "selected start convergence differs")
        end
        if !legacy_contract && model.status == "PASS"
            null_published && conditional_published ||
                _blocked(:row_count, "PASS fit omits a required start stage")
        end
        intersect(expected, start_seen) == expected || _blocked(:row_count, "start diagnostic coverage differs")
        union!(expected_start_all, expected)
    end
    start_seen == expected_start_all || _blocked(:row_count, "start diagnostic coverage has extra rows")
    start_aliases = Dict{Tuple{String,String,Int},Vector{String}}()
    for rows in values(start_rows_by_key), row in rows
        alias_key = (_sha_text(row[17], "start fit alias"), row[6], parse(Int, row[7]))
        payload = row[7:16]
        if haskey(start_aliases, alias_key)
            start_aliases[alias_key] == payload ||
                _blocked(:reference_mismatch, "same-fit start aliases disagree")
        else
            start_aliases[alias_key] = payload
        end
    end

    partition_by_identity = Dict((row[2], row[3], row[4], row[5]) => row
                                 for row in values(partition_rows))
    score_seen = Set{Tuple{String,String,String,String,String,String,String,String,String,String}}()
    score_rows_by_key = Dict{Tuple{String,String,String,String,String},Vector{Vector{String}}}()
    for row in report_rows["scores.tsv"]
        score_sha = _sha_text(row[13], "score hash")
        identity = (row[2], row[3], row[4], row[5])
        haskey(partition_by_identity, identity) || _blocked(:reference_mismatch, "score reference differs")
        partition = partition_by_identity[identity]
        model = references[(identity..., partition[12])][2]
        key = (identity..., score_sha)
        score_sha == partition[13] || _blocked(:hash_mismatch, "score hash differs from partition")
        row[6] in ("edge", "scan", "date") || _blocked(:schema_mismatch, "score level differs")
        if partition[6] != "PASS"
            row[6] == "date" && row[7:10] == ["NA", "NA", "NA", "NA"] ||
                _blocked(:row_count, "sparse terminal score coverage differs")
            row[11] == "NA" || _blocked(:schema_mismatch, "sparse terminal date score is numeric")
            _nonnegative_int(row[12], "sparse terminal date edge count") == 0 ||
                _blocked(:row_count, "sparse terminal date edge count differs")
            marker = (key..., row[6], "NA", "NA", "NA", "NA")
            marker in score_seen && _blocked(:duplicate_key, "score row repeats")
            push!(score_seen, marker)
            push!(get!(score_rows_by_key, key, Vector{Vector{String}}()), row)
            continue
        end
        _diag_out_float(row[11], "score gain"; allow_na=false)
        edge_count = _positive_int(row[12], "score edge count")
        heldout_edges = [edge for edge in _diag_eligible_edges(todo10, partition[12], transforms)
                         if edge.date == row[5]]
        if row[6] == "edge"
            file = row[7]; left = _positive_int(row[8], "score edge left lobe"); right = _positive_int(row[9], "score edge right lobe")
            edge = get(Dict((item.file, item.left_lobe, item.right_lobe) => item for item in heldout_edges), (file, left, right), nothing)
            edge === nothing && _blocked(:reference_mismatch, "score edge key is not held-out")
            row[10] == edge.segment_id || _blocked(:reference_mismatch, "score edge segment differs")
            edge_count == 1 || _blocked(:row_count, "edge score count differs")
            transform = transforms[(partition[12], edge.file, edge.date,
                                    edge.left_lobe, edge.right_lobe)]
            expected_gain = _diag_edge_gain(model, edge, transform, unaries)
            value = _diag_out_float(row[11], "score edge gain"; allow_na=false)::Float64
            isapprox(value, expected_gain; atol=1e-12, rtol=1e-12) ||
                _blocked(:score_mismatch, "score edge gain differs: $(_fmt(value)) != $(_fmt(expected_gain))")
            marker = (key..., row[6], file, string(left), string(right), row[10])
        elseif row[6] == "scan"
            row[8] == "NA" && row[9] == "NA" && row[10] == "NA" || _blocked(:schema_mismatch, "scan score exposes edge key")
            file = row[7]
            scan_edges = [edge for edge in heldout_edges if edge.file == file]
            isempty(scan_edges) && _blocked(:reference_mismatch, "scan score key is not held-out")
            edge_count == length(scan_edges) || _blocked(:row_count, "scan score edge count differs")
            marker = (key..., row[6], file, "NA", "NA", "NA")
        else
            row[7:10] == ["NA", "NA", "NA", "NA"] || _blocked(:schema_mismatch, "date score exposes edge key")
            edge_count == length(heldout_edges) || _blocked(:row_count, "date score edge count differs")
            marker = (key..., row[6], "NA", "NA", "NA", "NA")
        end
        marker in score_seen && _blocked(:duplicate_key, "score row repeats")
        push!(score_seen, marker)
        push!(get!(score_rows_by_key, key, Vector{Vector{String}}()), row)
    end
    for (identity, partition) in partition_by_identity
        fit_key = (identity..., partition[13])
        rows = get(score_rows_by_key, fit_key, Vector{Vector{String}}())
        expected = Tuple{String,String,String,String,String,String,String,String,String,String}[]
        if partition[6] == "PASS"
            heldout_edges = [edge for edge in _diag_eligible_edges(todo10, partition[12], transforms)
                             if edge.date == identity[4]]
            append!(expected, (fit_key..., "edge", edge.file, string(edge.left_lobe), string(edge.right_lobe), edge.segment_id) for edge in heldout_edges)
            append!(expected, (fit_key..., "scan", file, "NA", "NA", "NA") for file in sort(unique(edge.file for edge in heldout_edges)))
            push!(expected, (fit_key..., "date", "NA", "NA", "NA", "NA"))
        else
            push!(expected, (fit_key..., "date", "NA", "NA", "NA", "NA"))
        end
        actual = [(fit_key..., row[6], row[7], row[8], row[9], row[10]) for row in rows]
        actual == expected || _blocked(:row_count, "score key coverage/order differs")
    end

    shuffle_seen = Dict{Tuple{String,String,String,String,String},Vector{Vector{String}}}()
    shuffle_summary = Dict{Tuple{String,String,String,String,String},NamedTuple}()
    shuffle_fit_bindings = Dict{Tuple{String,Int},String}()
    for row in report_rows["shuffle.tsv"]
        score_sha = _sha_text(row[12], "shuffle score hash")
        identity = (row[2], row[3], row[4], row[5])
        haskey(partition_by_identity, identity) || _blocked(:reference_mismatch, "shuffle reference differs")
        partition = partition_by_identity[identity]
        key = (identity..., score_sha)
        seed = _nonnegative_int(row[6], "shuffle seed"); seed in 0:499 || _blocked(:schema_mismatch, "shuffle seed outside 0:499")
        score_sha == partition[13] || _blocked(:hash_mismatch, "shuffle score hash differs")
        running = row[7] != "NA"
        if running
            _diag_out_float(row[7], "shuffle value"; allow_na=false)
            _sha_text(row[8], "shuffle conditional fit hash")
            _diag_out_float(row[9], "shuffle control mean"; allow_na=false)
            _diag_out_float(row[10], "shuffle quantile"; allow_na=false)
        else
            row[8] == "not_run" || _blocked(:schema_mismatch, "non-running shuffle exposes a fit hash")
            partition[6] == "PASS" ?
                (_diag_out_float(row[9], "non-running shuffle control mean"; allow_na=false); row[10] == "NA") :
                (row[9] == "NA" && row[10] == "NA") ||
                _blocked(:schema_mismatch, "non-running shuffle exposes invalid numerical values")
        end
        _diag_bool(row[11], "shuffle pass")
        push!(get!(shuffle_seen, key, Vector{Vector{String}}()), row)
    end
    for (identity, partition) in partition_by_identity
        key = (identity..., partition[13])
        rows = get(shuffle_seen, key, Vector{Vector{String}}())
        length(rows) == 500 || _blocked(:row_count, "shuffle seed coverage differs")
        [parse(Int, row[6]) for row in rows] == collect(0:499) || _blocked(:order_mismatch, "shuffle seed order differs")
        running = first(rows)[7] != "NA"
        all((row[7] != "NA") == running for row in rows) || _blocked(:status_mismatch, "shuffle run state differs")
        if running
            shuffle_values_local = [_diag_out_float(row[7], "shuffle value"; allow_na=false)::Float64 for row in rows]
            quantile = _diag_type7(shuffle_values_local, 0.975)
            controls = [_diag_out_float(row[9], "shuffle control mean"; allow_na=false)::Float64 for row in rows]
            quantiles = [_diag_out_float(row[10], "shuffle quantile"; allow_na=false)::Float64 for row in rows]
            all(isapprox(value, quantile; atol=1e-12, rtol=1e-12) for value in quantiles) ||
                _blocked(:hash_mismatch, "shuffle quantile differs")
            all(isapprox(value, first(controls); atol=1e-12, rtol=1e-12) for value in controls) ||
                _blocked(:reference_mismatch, "shuffle control means differ")
            pass_value = first(controls) > quantile
            all(_diag_bool(row[11], "shuffle pass") == pass_value for row in rows) || _blocked(:status_mismatch, "shuffle pass predicate differs")
            _diag_bool(partition[17], "partition shuffle pass") == pass_value || _blocked(:status_mismatch, "partition shuffle pass differs")
            fit_hashes = [_sha_text(row[8], "shuffle conditional fit hash") for row in rows]
            for row in rows
                seed_value = parse(Int, row[6])
                binding_key = (partition[12], seed_value)
                if haskey(shuffle_fit_bindings, binding_key)
                    shuffle_fit_bindings[binding_key] == row[8] ||
                        _blocked(:hash_mismatch, "shuffle conditional fit binding differs")
                else
                    shuffle_fit_bindings[binding_key] = row[8]
                end
            end
            published_quantile = _diag_out_float(first(rows)[10], "published shuffle quantile"; allow_na=false)::Float64
            shuffle_summary[key] = (running=true, quantile=published_quantile,
                                    control=first(controls), passed=pass_value,
                                    values=shuffle_values_local)
        else
            if partition[6] == "PASS"
                all(row[8] == "not_run" && row[10] == "NA" && row[11] == "true" &&
                    _diag_out_float(row[9], "non-running shuffle control mean"; allow_na=false) !== nothing
                    for row in rows) || _blocked(:status_mismatch, "non-running shuffle semantics differ")
                partition[17] == "true" || _blocked(:status_mismatch, "non-running partition shuffle pass differs")
            else
                all(row[8] == "not_run" && row[9] == "NA" && row[10] == "NA" && row[11] == "false" for row in rows) ||
                    _blocked(:status_mismatch, "sparse non-running shuffle semantics differ")
                partition[17] == "false" || _blocked(:status_mismatch, "sparse non-running partition shuffle pass differs")
            end
            if partition[6] == "PASS"
                control = _diag_out_float(first(rows)[9], "non-running shuffle control mean"; allow_na=false)::Float64
                shuffle_summary[key] = (running=false, quantile=-Inf,
                                        control=control, passed=true,
                                        values=fill(-Inf, 500))
            else
                shuffle_summary[key] = (running=false, quantile=-Inf,
                                        control=NaN, passed=false,
                                        values=fill(-Inf, 500))
            end
        end
    end

    for (identity, partition) in partition_by_identity
        fit_key = (identity..., partition[13])
        rows = get(score_rows_by_key, fit_key, Vector{Vector{String}}())
        if partition[6] != "PASS"
            length(rows) == 1 && rows[1][6] == "date" || _blocked(:row_count, "sparse terminal partition score coverage differs")
            fit_model = references[(identity..., partition[12])][2]
            expected_score = _diag_invalid_score_hash(fit_model.status, fit_model.reason,
                                                      partition[12], partition[5])
            partition[13] == expected_score || _blocked(:hash_mismatch, "sparse partition score hash differs")
            continue
        end
        summary = get(shuffle_summary, fit_key, nothing)
        summary === nothing && _blocked(:row_count, "partition shuffle summary is absent")
        scan_scores = Dict{String,Float64}()
        gains = Dict{String,Float64}()
        date_mean = nothing
        model = references[(identity..., partition[12])][2]
        eligible_edges = [edge for edge in _diag_eligible_edges(todo10, partition[12], transforms)
                          if edge.date == identity[4]]
        expected_gains = Dict(_diag_edge_identity(edge) =>
                              _diag_edge_gain(model, edge,
                                              transforms[(partition[12], edge.file, edge.date,
                                                          edge.left_lobe, edge.right_lobe)], unaries)
                              for edge in eligible_edges)
        expected_scan_values = Dict{String,Vector{Float64}}()
        for edge in eligible_edges
            push!(get!(expected_scan_values, edge.file, Float64[]),
                  expected_gains[_diag_edge_identity(edge)])
        end
        expected_scan_means = Dict(scan => mean(values)
                                   for (scan, values) in expected_scan_values)
        expected_date_mean = mean(expected_scan_means[scan]
                                  for scan in sort(collect(keys(expected_scan_means))))
        for row in rows
            value = _diag_out_float(row[11], "score gain"; allow_na=false)::Float64
            if row[6] == "edge"
                identity_text = join((row[7], row[8], row[9], row[10]), '\t')
                gains[identity_text] = value
            elseif row[6] == "scan"
                scan_scores[row[7]] = value
                haskey(expected_scan_means, row[7]) || _blocked(:reference_mismatch, "score scan is not reconstructed")
                isapprox(value, expected_scan_means[row[7]]; atol=1e-12, rtol=1e-12) ||
                    _blocked(:score_mismatch, string("score scan mean differs: ", _fmt(value), " != ", _fmt(expected_scan_means[row[7]])))
            else
                date_mean = value
            end
        end
        date_mean === nothing && _blocked(:row_count, "date score is absent")
        Set(keys(gains)) == Set(keys(expected_gains)) || _blocked(:row_count, "score edge coverage differs")
        for key in keys(expected_gains)
            isapprox(gains[key], expected_gains[key]; atol=1e-12, rtol=1e-12) ||
                _blocked(:score_mismatch, "score edge aggregate differs: $(_fmt(gains[key])) != $(_fmt(expected_gains[key]))")
        end
        Set(keys(scan_scores)) == Set(keys(expected_scan_means)) || _blocked(:row_count, "score scan coverage differs")
        isapprox(date_mean::Float64, expected_date_mean; atol=1e-12, rtol=1e-12) ||
            _blocked(:score_mismatch, "score date mean differs: $(_fmt(date_mean::Float64)) != $(_fmt(expected_date_mean))")
        expected_score = _diag_score_hash(partition[12], partition[5], date_mean::Float64,
            summary.control, summary.quantile, summary.passed, scan_scores, gains,
            summary.values)
        partition[13] == expected_score || _blocked(:hash_mismatch, "partition score hash differs: $(partition[13]) != $(expected_score)")
    end

    bootstrap_seen = Dict{Tuple{String,String,String},Vector{Vector{String}}}()
    for row in report_rows["bootstrap.tsv"]
        seed = _nonnegative_int(row[5], "bootstrap seed"); seed in 0:499 || _blocked(:schema_mismatch, "bootstrap seed outside 0:499")
        _sha_text(row[13], "bootstrap gate hash")
        _diag_out_float(row[6], "bootstrap mean"; allow_na=true)
        _diag_out_float(row[7], "bootstrap lower"; allow_na=true)
        for index in 8:10
            _diag_bool(row[index], "bootstrap boolean")
        end
        row[11] in ("PASS", "SKIPPED", "FAIL", "BLOCKED") || _blocked(:status_mismatch, "bootstrap status invalid")
        key = (row[2], row[3], row[4]); push!(get!(bootstrap_seen, key, Vector{Vector{String}}()), row)
    end
    gate_hashes = Dict{Tuple{String,String,String},String}()
    gate_statuses = Dict{Tuple{String,String,String},String}()
    gate_reasons = Dict{Tuple{String,String,String},String}()
    expected_gate_keys = Set{Tuple{String,String,String}}()
    for model in values(models)
        push!(expected_gate_keys, (model.model, "final_lodo", "NA"))
        for reference in model.references
            reference.fit_role == "partition" && reference.scope == "outer_inner" &&
                push!(expected_gate_keys, (reference.model, reference.scope, reference.outer_date))
        end
    end
    Set(keys(bootstrap_seen)) == expected_gate_keys || _blocked(:row_count, "bootstrap gate coverage differs")
    _required_int(receipt, "bootstrap_count", "Todo 11 receipt") == 500 || _blocked(:row_count, "bootstrap count differs")
    _required_int(receipt, "shuffle_count", "Todo 11 receipt") == 500 || _blocked(:row_count, "shuffle count differs")
    for key in sort!(collect(expected_gate_keys))
        rows = bootstrap_seen[key]
        [parse(Int, row[5]) for row in rows] == collect(0:499) || _blocked(:order_mismatch, "bootstrap seed order differs")
        partition_rows_for_gate = [row for (identity, row) in partition_by_identity
                                   if (identity[1], identity[2], identity[3]) == key]
        dates = sort(collect(row[5] for row in partition_rows_for_gate))
        if isempty(partition_rows_for_gate)
            key[2] == "final_lodo" && key[3] == "NA" || _blocked(:row_count, "gate has no serialized partition references")
            if legacy_contract
                degenerate_values = [_diag_out_float(row[6], "degenerate gate bootstrap value"; allow_na=false)::Float64 for row in rows]
                lower = _diag_out_float(first(rows)[7], "degenerate gate lower"; allow_na=false)::Float64
                all(row[8:10] == ["true", "true", "true"] && row[11] == "PASS" && row[12] == "ok" for row in rows) ||
                    _blocked(:status_mismatch, "degenerate gate status differs")
                expected_hash = _diag_gate_hash("PASS", "ok", String[], Dict{String,Float64}(), degenerate_values,
                                                lower, true, true, true)
                all(row[13] == expected_hash for row in rows) || _blocked(:hash_mismatch, "degenerate gate hash differs")
                gate_hashes[key] = expected_hash
                gate_statuses[key] = "PASS"
                gate_reasons[key] = "ok"
                continue
            end
            status = first(rows)[11]
            reason = first(rows)[12]
            status in ("BLOCKED", "FAIL", "SKIPPED") ||
                _blocked(:status_mismatch, "empty gate status differs")
            all(row[6:7] == ["NA", "NA"] && row[8:10] == ["false", "false", "false"] &&
                row[11] == status && row[12] == reason for row in rows) ||
                _blocked(:status_mismatch, "empty gate NA semantics differ")
            gate_model = only(model for model in values(models)
                              if model.model == key[1] && model.fit_role == "full_refit")
            gate_dates = key[2] == "final_lodo" ? collect(gate_model.training_dates) : String[]
            expected_hash = _diag_gate_hash(status, reason, gate_dates, Dict{String,Float64}(),
                                            Float64[], NaN, false, false, false)
            all(row[13] == expected_hash for row in rows) || _blocked(:hash_mismatch, "degenerate gate hash differs")
            gate_hashes[key] = expected_hash
            gate_statuses[key] = status
            gate_reasons[key] = reason
            continue
        end
        statuses = [row[6] for row in partition_rows_for_gate]
        status = any(==( "BLOCKED"), statuses) ? "BLOCKED" :
                 any(==( "FAIL"), statuses) ? "FAIL" :
                 any(==( "SKIPPED"), statuses) ? "SKIPPED" : "PASS"
        invalid = status != "PASS"
        if invalid
            reason = status == "BLOCKED" ? "partition_blocked" : status == "FAIL" ? "partition_failed" : "partition_skipped"
            expected_hash = _diag_gate_hash(status, reason, dates, Dict{String,Float64}(),
                                            Float64[], NaN, false, false, false)
            all(row[6] == "NA" && row[7] == "NA" && row[8:10] == ["false", "false", "false"] &&
                row[11] == status && row[12] == reason && row[13] == expected_hash for row in rows) ||
                _blocked(:status_mismatch, "invalid bootstrap gate semantics differ")
        else
            date_means = Dict{String,Float64}()
            scan_scores = Dict{String,Dict{String,Float64}}()
            for partition in partition_rows_for_gate
                fit_key = (partition[2], partition[3], partition[4], partition[5], partition[13])
                score_rows = get(score_rows_by_key, fit_key, Vector{Vector{String}}())
                scans = Dict(row[7] => (_diag_out_float(row[11], "bootstrap source scan score"; allow_na=false)::Float64)
                             for row in score_rows if row[6] == "scan")
                date_row = only(row for row in score_rows if row[6] == "date")
                date_means[partition[5]] = _diag_out_float(date_row[11], "bootstrap source date score"; allow_na=false)::Float64
                scan_scores[partition[5]] = scans
            end
            every_positive = all(value -> value > 0.0, values(date_means))
            every_shuffle = all(row[17] == "true" for row in partition_rows_for_gate)
            reversal = all(row[16] == "true" for row in partition_rows_for_gate)
            bootstrap_values = _diag_bootstrap_values(scan_scores)
            lower = _diag_type7(bootstrap_values, 0.025)
            status = every_positive && lower > 0.0 && every_shuffle && reversal ? "PASS" : "SKIPPED"
            reason = !every_positive ? "nonpositive_date" : lower <= 0.0 ? "bootstrap_not_positive" :
                     !every_shuffle ? "shuffle_control_failed" : !reversal ? "reversal_inequivalent" : "ok"
            expected_hash = _diag_gate_hash(status, reason, dates, date_means, bootstrap_values,
                                            lower, every_positive, every_shuffle, reversal)
            for (index, row) in enumerate(rows)
                row[6] == _fmt(bootstrap_values[index]) && row[7] == _fmt(lower) &&
                    row[8] == string(every_positive) && row[9] == string(every_shuffle) &&
                    row[10] == string(reversal) && row[11] == status && row[12] == reason &&
                    row[13] == expected_hash || _blocked(:hash_mismatch, "bootstrap gate reconstruction differs")
            end
        end
        gate_hashes[key] = expected_hash
        gate_statuses[key] = status
        gate_reasons[key] = reason
    end

    model_result_hashes = Dict{String,String}()
    component_statuses = Dict{String,Vector{String}}(model_id => String[] for model_id in keys(overall_rows))
    for (fit, model) in models
        push!(component_statuses[model.model], model.status)
    end
    for (key, row) in partition_rows
        key[2] == "outer_score" && push!(component_statuses[key[1]], row[6])
    end
    for (key, gate_status) in gate_statuses
        push!(component_statuses[key[1]], gate_status)
    end
    precedence(statuses) = any(==( "BLOCKED"), statuses) ? "BLOCKED" : any(==( "FAIL"), statuses) ? "FAIL" : any(==( "SKIPPED"), statuses) ? "SKIPPED" : "PASS"
    for model_id in keys(overall_rows)
        row = overall_rows[model_id]
        expected_status = precedence(component_statuses[model_id])
        expected_reason = expected_status == "PASS" ? "ok" : expected_status == "FAIL" ? "numerical_failure" :
                         expected_status == "BLOCKED" ? "integrity_failure" :
                         (!legacy_contract && length(all_dates) < _PRODUCER_COMPLETE_DATE_COUNT ?
                          "insufficient_dates" : "conditional_mechanism_not_admitted")
        row[3] == expected_status || _blocked(:status_mismatch, "overall/fitted status precedence differs")
        row[4] == expected_reason || _blocked(:status_mismatch, "overall model reason differs")
        full_fit = only(reference.fit_sha256 for model in values(models) for reference in model.references if reference.model == model_id && reference.fit_role == "full_refit")
        final_gate = get(gate_hashes, (model_id, "final_lodo", "NA"), nothing)
        final_gate === nothing && _blocked(:row_count, "final gate bootstrap coverage is absent")
        outer_dates = sort(unique(reference.outer_date for model in values(models) for reference in model.references if reference.model == model_id && reference.fit_role == "partition" && reference.scope == "outer_inner"))
        _nonnegative_int(row[5], "overall outer fold count") == length(outer_dates) || _blocked(:row_count, "overall outer fold count differs")
        outer_lines = String[]
        for outer_date in outer_dates
            inner_gate = get(gate_hashes, (model_id, "outer_inner", outer_date), nothing)
            inner_gate === nothing && _blocked(:row_count, "outer-inner gate coverage is absent")
            outer_gate_row = only(row for (key, row) in partition_rows if key[1] == model_id && key[2] == "outer_score" && key[3] == outer_date)
            push!(outer_lines, join(("outer", outer_date, inner_gate, outer_gate_row[13], outer_gate_row[12]), '\t'))
        end
        model_lines = if !legacy_contract && length(all_dates) < _PRODUCER_COMPLETE_DATE_COUNT
            ["model=$(model_id)", "status=$(row[3])", "reason=$(row[4])"]
        else
            lines = ["schema=structured-edge-model-admission-v1", "model=$(model_id)",
                     "status=$(row[3])", "reason=$(row[4])",
                     "final_gate=$(final_gate)", "full_fit=$(full_fit)"]
            append!(lines, outer_lines)
            lines
        end
        result_hash = _hash_lines(model_lines)
        row[6] == final_gate && row[7] == full_fit && row[8] == result_hash || _blocked(:hash_mismatch, "overall model result binding differs")
        model_result_hashes[model_id] = result_hash
    end
    receipt_status = precedence([row[3] for row in values(overall_rows)])
    receipt_reason = receipt_status == "PASS" ? "ok" : receipt_status == "FAIL" ? "numerical_failure" : receipt_status == "BLOCKED" ? "integrity_failure" : "conditional_mechanism_not_admitted"
    normalization = _sha_text(receipt["hashes"]["normalization_sha256"], "Todo 11 normalization hash")
    result_lines = [
        "schema=structured-edge-admission-report-v1",
        "status=$(receipt_status)",
        "reason=$(receipt_reason)",
        "normalization=$(normalization)",
    ]
    append!(result_lines, "model=$(model_id):$(model_result_hashes[model_id])" for model_id in sort!(collect(keys(model_result_hashes))))
    for key in sort!(collect(keys(receipt["hashes"])))
        key == "normalization_sha256" && continue
        push!(result_lines, "input=$(key):$(receipt["hashes"][key])")
    end
    result_hash = _hash_lines(result_lines)
    _required_string(receipt, "status", "Todo 11 receipt") == receipt_status || _blocked(:status_mismatch, "receipt status differs from diagnostic precedence")
    _required_string(receipt, "reason", "Todo 11 receipt") == receipt_reason || _blocked(:status_mismatch, "receipt reason differs from diagnostic precedence")
    _required_string(receipt, "result_sha256", "Todo 11 receipt") == result_hash || _blocked(:hash_mismatch, "Todo 11 result hash differs from diagnostics")
    return result_hash, normalization, model_result_hashes, receipt_status, receipt_reason
end

function _residualizer_sha(model::_ModelData, training_edge_count::Int)
    model.residualizer_sha256 == "NA" && return "NA"
    model.residualizer_training_sha256 != "NA" || _blocked(:model_mismatch, "residualizer training hash is absent")
    model.residualizer_ridge !== nothing || _blocked(:model_mismatch, "residualizer ridge is absent")
    training_edge_count > 0 || _blocked(:model_mismatch, "residualizer training edge count is absent")
    order = model.reversed ? vcat(1, collect(9:15), collect(2:8)) : collect(1:15)
    lines = String[
        "schema=structured-edge-residualizer-v1",
        "training_sha256=$(model.residualizer_training_sha256)",
        "training_row_count=$(training_edge_count)",
        "ridge=$(_fmt(model.residualizer_ridge::Float64))",
        "evaluation_order=$(join(order, ','))",
    ]
    coefficients = model.residual_coefficients
    coefficients isa Matrix{Float64} || _blocked(:model_mismatch, "residualizer coefficients are absent")
    for (name, row) in zip(_GRAPH_COEFFICIENT_NAMES, axes(coefficients, 1)),
        (output, column) in zip(_GRAPH_OUTPUT_NAMES, axes(coefficients, 2))
        push!(lines, "$(name)_$(output)=$(_fmt(coefficients[row, column]))")
    end
    return _hash_lines(lines)
end

function _reference_fields(reference::String)
    fields = split(reference, '|'; keepempty=true)
    length(fields) == 6 || _blocked(:reference_mismatch, "reference has wrong field count")
    fields
end

function _reference_fields(reference::_ReferenceData)
    return [reference.fit_role, reference.model, reference.scope, reference.outer_date,
            reference.target_date, reference.fit_sha256]
end

function _student_log_density(residual::Vector{Float64}, mean::Vector{Float64}, scale::Matrix{Float64})
    all(isfinite, residual) && all(isfinite, mean) && all(isfinite, scale) || _infer_fail(:density_failure, "density inputs are nonfinite")
    determinant = det(scale)
    determinant > 0.0 && isfinite(determinant) || _infer_fail(:density_failure, "density scale is not positive definite")
    delta = dot(residual - mean, scale \ (residual - mean))
    isfinite(delta) || _infer_fail(:density_failure, "density quadratic is nonfinite")
    return log(24.0) - log(6.0) - log(8.0 * π) - 0.5 * log(determinant) - 5.0 * log1p(delta / 8.0)
end

function _factor(model::_ModelData, edge::_Todo10Edge, transform::_TransformData)
    _FACTOR_BUILD_CALL_COUNT[] += 1
    if _DENSITY_FAILURE_HOOK[] !== nothing
        _DENSITY_FAILURE_HOOK[]()
        _infer_fail(:density_failure, "injected density-stage failure")
    end
    transform.status == "eligible" || return nothing
    raw = [edge.corr_fwd - transform.pred_fwd, edge.corr_bwd - transform.pred_bwd]
    null = _student_log_density(raw, model.null_mean, model.null_scale)
    ratios = ntuple(i -> _student_log_density(raw, vec(model.conditional_means[i, :]), model.conditional_scale) - null, 4)
    all(isfinite, ratios) || _infer_fail(:density_failure, "factor is nonfinite")
    return ratios
end

function _logsumexp(values)
    maximum_value = maximum(values)
    maximum_value == -Inf && return -Inf
    return maximum_value + log(sum(exp(value - maximum_value) for value in values))
end

@inline function _two_diff(left::Float64, right::Float64)::NTuple{2,Float64}
    difference = left - right
    b_virtual = left - difference
    a_virtual = difference + b_virtual
    b_roundoff = b_virtual - right
    a_roundoff = left - a_virtual
    return (difference, a_roundoff + b_roundoff)
end

@inline function _exp_loser_difference(loser::Float64, winner::Float64)::Float64
    delta_hi, delta_lo = _two_diff(loser, winner)
    delta_hi == -Inf && return 0.0

    return exp(delta_hi) * exp(delta_lo)
end

function _binary_softmax(l0::Float64, l1::Float64)::NTuple{2,Float64}
    valid_l0 = isfinite(l0) || l0 == -Inf
    valid_l1 = isfinite(l1) || l1 == -Inf
    valid_l0 && valid_l1 && !(l0 == -Inf && l1 == -Inf) ||
        _infer_fail(:numerical_failure, "binary softmax input is invalid")
    if l0 >= l1
        ratio = _exp_loser_difference(l1, l0)
        quotient = ratio / (1.0 + ratio)
        tail = fma(-ratio, quotient, ratio)
        return (1.0 - tail, tail)
    else
        ratio = _exp_loser_difference(l0, l1)
        quotient = ratio / (1.0 + ratio)
        tail = fma(-ratio, quotient, ratio)
        return (tail, 1.0 - tail)
    end
end

function _argmax_lower(values)
    best = 1
    for index in 2:length(values)
        values[index] > values[best] && (best = index)
    end
    return best
end

function _dp_block(unaries::Vector{NTuple{2,Float64}}, factors::Vector{NTuple{4,Float64}})
    _DP_CALL_COUNT[] += 1
    if _DP_FAILURE_HOOK[] !== nothing
        _DP_FAILURE_HOOK[]()
        _infer_fail(:dp_failure, "injected DP-stage failure")
    end
    n = length(unaries)
    n >= 1 || _infer_fail(:internal, "empty DP block")
    forward = fill(-Inf, n, 2)
    backward = fill(-Inf, n, 2)
    forward[1, :] .= unaries[1]
    for i in 2:n, state in 1:2
        factor = factors[i - 1]
        forward[i, state] = unaries[i][state] + _logsumexp([
            forward[i - 1, 1] + factor[(1 - 1) * 2 + state],
            forward[i - 1, 2] + factor[(2 - 1) * 2 + state],
        ])
    end
    backward[n, :] .= 0.0
    for i in (n - 1):-1:1, state in 1:2
        factor = factors[i]
        backward[i, state] = _logsumexp([
            factor[(state - 1) * 2 + 1] + unaries[i + 1][1] + backward[i + 1, 1],
            factor[(state - 1) * 2 + 2] + unaries[i + 1][2] + backward[i + 1, 2],
        ])
    end
    log_evidence = _logsumexp([forward[n, 1], forward[n, 2]])
    isfinite(log_evidence) || _infer_fail(:dp_failure, "no feasible chain state")
    log_marginals = forward .+ backward .- log_evidence
    marginals = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        marginals[i, :] .= _binary_softmax(log_marginals[i, 1], log_marginals[i, 2])
    end
    viterbi = fill(-Inf, n, 2)
    predecessors = fill(1, n, 2)
    viterbi[1, :] .= unaries[1]
    for i in 2:n, state in 1:2
        factor = factors[i - 1]
        values = [viterbi[i - 1, 1] + factor[state], viterbi[i - 1, 2] + factor[2 + state]]
        predecessor = _argmax_lower(values)
        predecessors[i, state] = predecessor
        viterbi[i, state] = unaries[i][state] + values[predecessor]
    end
    state = _argmax_lower([viterbi[n, 1], viterbi[n, 2]])
    viterbi_score = viterbi[n, state]
    labels = fill(0, n)
    labels[n] = state - 1
    for i in (n - 1):-1:1
        state = predecessors[i + 1, state]
        labels[i] = state - 1
    end
    return log_evidence, log_marginals, marginals, labels, viterbi_score
end

function _node_result(node::_Todo10Node, unary::_UnaryData, log_marginals, marginals, label)
    return ChainNodeResult(node.file, node.date, node.lobe, node.t_nm,
                           (log_marginals[1], log_marginals[2]),
                           (marginals[1], marginals[2]), label, unary.unary_fit_sha256)
end

function _model_for_reference(t11::_Todo11Data, reference::_ReferenceData)
    model = only(filter(item -> item.fit_sha256 == reference.fit_sha256, t11.models))
    model.model == reference.model && reference in model.references ||
        _blocked(:reference_mismatch, "reference does not match model row")
    return model
end

function _scan_unary_map!(cache, t11::_Todo11Data, fit_sha::String, file::String)
    key = (fit_sha, file)
    if !haskey(cache, key)
        _SCAN_DECODE_CALL_COUNT[] += 1
        cache[key] = Dict((item.file, item.lobe) => item for
                          (identity, item) in t11.unaries
                          if identity[1] == fit_sha && identity[2] == file)
    end
    return cache[key]
end

function _reference_result(reference::_ReferenceData, model::_ModelData,
                           todo10::_Todo10Data, t11::_Todo11Data,
                           factor_cache, scan_cache; unary_only=false)
    model.status == "FAIL" && return _reference_result_from_failure(reference)
    selected_nodes = reference.fit_role == "partition" ?
                     [node for node in todo10.nodes if node.date == reference.target_date] :
                     collect(todo10.nodes)
    isempty(selected_nodes) && _blocked(:reference_mismatch, "reference selects no scan nodes")
    for file in unique(node.file for node in selected_nodes)
        unary_map = _scan_unary_map!(scan_cache, t11, reference.fit_sha256, file)
        all(haskey(unary_map, (node.file, node.lobe)) for node in selected_nodes if node.file == file) ||
            _blocked(:unary_coverage, "reference unary coverage is incomplete")
    end
    unary_map(node) = _scan_unary_map!(scan_cache, t11, reference.fit_sha256, node.file)[(node.file, node.lobe)]
    edge_map = Dict((edge.file, edge.date, edge.left_lobe, edge.right_lobe) => edge for edge in todo10.edges)
    nodes_by_segment = Dict{String,Vector{_Todo10Node}}()
    for node in selected_nodes
        push!(get!(nodes_by_segment, node.segment_id, _Todo10Node[]), node)
    end
    blocks = ChainBlockResult[]
    usable_factor = false
    total_log_evidence = 0.0
    total_viterbi = 0.0
    for segment_id in sort!(collect(keys(nodes_by_segment)))
        chain = sort!(nodes_by_segment[segment_id]; by=node -> (node.file, node.t_nm, node.lobe))
        if unary_only || model.status == "SKIPPED"
            block = _unary_only_block(segment_id, chain,
                                      Dict((node.file, node.lobe) => unary_map(node) for node in chain))
            push!(blocks, block)
            total_log_evidence += block.log_evidence
            total_viterbi += block.viterbi_score
            continue
        end
        pieces = Vector{Tuple{Vector{_Todo10Node},Vector{ChainFactorResult},Vector{NTuple{4,Float64}}}}()
        current_nodes = [_Todo10Node(chain[1].file, chain[1].date, chain[1].lobe, chain[1].t_nm,
                                     chain[1].t_text, chain[1].segment_id, chain[1].status)]
        current_factors = ChainFactorResult[]
        current_values = NTuple{4,Float64}[]
        for index in 1:(length(chain) - 1)
            left = chain[index]
            right = chain[index + 1]
            edge = get(edge_map, (left.file, left.date, left.lobe, right.lobe), nothing)
            factor = nothing
            factor_result = nothing
            if edge !== nothing && edge.status == "eligible"
                transform = get(t11.transforms,
                                (reference.fit_sha256, edge.file, edge.date,
                                 edge.left_lobe, edge.right_lobe), nothing)
                transform === nothing && _blocked(:transform_coverage, "reference transform is absent")
                cache_key = (reference.fit_sha256, edge.file, edge.date,
                             edge.left_lobe, edge.right_lobe)
                if haskey(factor_cache, cache_key)
                    factor = factor_cache[cache_key]
                else
                    factor = _factor(model, edge, transform)
                    factor_cache[cache_key] = factor
                end
                factor === nothing && _infer_fail(:factor_failure, "eligible factor is unavailable")
                factor_result = ChainFactorResult(
                    edge.file, edge.date, edge.left_lobe, edge.right_lobe,
                    edge.left_t_nm, edge.right_t_nm, edge.segment_id, factor,
                    _factor_sha(edge, model, factor), edge.raw_sha256,
                    model.fit_sha256, model.model, model.unary_fit_sha256,
                    model.residualizer_sha256)
            end
            if factor === nothing
                push!(pieces, (current_nodes, current_factors, current_values))
                current_nodes = [_Todo10Node(right.file, right.date, right.lobe, right.t_nm,
                                             right.t_text, right.segment_id, right.status)]
                current_factors = ChainFactorResult[]
                current_values = NTuple{4,Float64}[]
            else
                push!(current_factors, factor_result)
                push!(current_values, factor)
                push!(current_nodes, right)
                usable_factor = true
            end
        end
        push!(pieces, (current_nodes, current_factors, current_values))
        for (natural_nodes, factor_results, natural_factors) in pieces
            natural_unaries = NTuple{2,Float64}[]
            for node in natural_nodes
                unary = unary_map(node)
                value = (unary.p0 == 0.0 ? -Inf : log(unary.p0),
                         unary.p1 == 0.0 ? -Inf : log(unary.p1))
                isfinite(value[1]) || isfinite(value[2]) || _infer_fail(:unary_failure, "unary has no feasible state")
                push!(natural_unaries, value)
            end
            natural_dp = _dp_block(natural_unaries, natural_factors)
            reversed_unaries = reverse(natural_unaries)
            reversed_factors = [_transpose_factor(factor) for factor in reverse(natural_factors)]
            reversed_dp = _dp_block(reversed_unaries, reversed_factors)
            _assert_reversal_equivalence(natural_dp, reversed_dp)
            primary = model.reversed ? reversed_dp : natural_dp
            log_marginals = model.reversed ? primary[2][end:-1:1, :] : primary[2]
            marginals = model.reversed ? primary[3][end:-1:1, :] : primary[3]
            labels = model.reversed ? primary[4][end:-1:1] : primary[4]
            node_results = ChainNodeResult[]
            for index in eachindex(natural_nodes)
                unary = unary_map(natural_nodes[index])
                push!(node_results, _node_result(natural_nodes[index], unary,
                                                 view(log_marginals, index, :),
                                                 view(marginals, index, :), labels[index]))
            end
            block = ChainBlockResult(segment_id, Tuple(node_results), Tuple(factor_results),
                                     primary[1], primary[5], Tuple(labels))
            push!(blocks, block)
            total_log_evidence += block.log_evidence
            total_viterbi += block.viterbi_score
        end
    end
    ref_status = usable_factor ? :PASS : :SKIPPED
    ref_reason = usable_factor ? :ok : :unary_only
    return ChainReferenceResult(reference.fit_sha256, reference.model, reference.fit_role,
                                reference.scope, reference.outer_date, reference.target_date,
                                ref_status, ref_reason, Tuple(blocks), total_log_evidence,
                                total_viterbi)
end
