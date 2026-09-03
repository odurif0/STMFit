#!/usr/bin/env julia

module StructuredStabilityAudit

using Base64
using Printf
using Random
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "lib", "structured_assignment", "universe.jl"))
include(joinpath(@__DIR__, "lib", "structured_assignment", "evidence.jl"))
include(joinpath(@__DIR__, "lib", "structured_assignment", "edge_features.jl"))
include(joinpath(@__DIR__, "lib", "structured_assignment", "chain_inference.jl"))

using .StructuredUniverse
using .StructuredEvidence
using .StructuredEdgeFeatures
using .StructuredChainInference

export audit_stability, report_files, publish_report, main

const CANDIDATE_CONFIG_SHA256 =
    "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const MODEL_CONFIG_SHA256 =
    "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const CHITOSAN_CONFIG_SHA256 =
    "a24460ab5114bd3907c063c1f6d7b678113037e084881a7fbbb211db9eeaf2d9"
const EVIDENCE_SOURCE_SHA256 =
    "3bdeeebbce06d86f43a77542185b499f858eb176deaf3461a524ef12111cd189"
const UNIVERSE_SOURCE_SHA256 =
    "a70a6e1382627fd584e5825dbbc6464d24ba5b7d96dabe7ecc36493c79084b18"
const EDGE_SOURCE_SHA256 =
    "c29fe06c7fc0ba33c4e168169fec5ced5a7907bbfcec626c183e3ced044da0f9"
const CHAIN_SOURCE_SHA256 =
    "cdc0c788d49298f50721b091535a238cd552454d9885404456c593e4625ea39f"
const EDGE_MODEL_SHA256 =
    "6d73aba926d5203e631cf904dc59804559f2b172ce2ff0f4f808fbf33b1690f3"

const T8_SOURCE_ENTRIES = (
    ("test/build_structured_unit_predictions.jl",
     "976459ce93d51d6384a63a7d5d28c0a72e4dcc65f9d61500e8a1d314b52b6585"),
    ("test/lib/hierarchical/emission_math.jl",
     "15792061ca23ae59bf6994319afc209bd8e95ca283754fc2844eddcaa1482762"),
    ("test/lib/hierarchical/emissions.jl",
     "ff98ef1e9ae77a4c0d25f4a53c7a974977e88582b70206a685d1ad79b6f8aaae"),
    ("test/lib/hierarchical/firewall.jl",
     "0502745ff8187b8dd926db3760e336404c068114e0da6bd493b284e482814bed"),
    ("test/lib/hierarchical/io_helpers.jl",
     "5e5e5bca2635c1da03250988d2e283cdcc1ef23a931570a11bbc54841d6918b5"),
    ("test/lib/hierarchical/loading.jl",
     "41171af123d489f8e3d4e20ba0ad2e55641a3cbb5768b3debc711c3f8262422a"),
    ("test/lib/hierarchical/nuisance.jl",
     "f256eb82966350544b55a23413417aa6eacda7cf6669c049f8dc90feb4fd2203"),
    ("test/lib/hierarchical/pipeline.jl",
     "31eeebc66ce9bb6072dc03cb5db95731909dce7d962beaaa528e2f8695bb95fa"),
    ("test/lib/hierarchical/views.jl",
     "40380fb586a7866e718225eac3fd610e4eba9a0e0c3bddb924f108655e5fdfcf"),
    ("test/lib/hierarchical_unit_assignment.jl",
     "4024d47b82077436a1c22a9511930effc5bca110f4c9eba5daef487b923513a8"),
    ("test/lib/structured_assignment/champion_adapter.jl",
     "470a4c6676eaa3b0a00a51a64f9fce96ac7282fd1adb84c921dc3b967f6aa83c"),
    ("test/lib/structured_assignment/firewall.jl",
     "a189898d31759352d9dac0de8ff281d6001e1dc2db94df2cb725ff597537a9bb"),
    ("test/lib/structured_assignment/robust_emissions.jl",
     "0ca4a3863a4b9ed13aee5cd3a0229df4f2764a8524b764eeb0fd27238f14f241"),
    ("test/lib/structured_unit_assignment.jl",
     "8c74f21924381146d58dd9611d5a75578e7b691db3c16c847164549fd48e30b9"),
)

const NODE_HEADER = [
    "schema_version", "evidence_scope", "fold_date", "scan_date", "file", "lobe", "model",
    "control_class", "perturbation", "view", "seed_family", "seed", "target", "posterior_1",
    "hard_output", "node_status", "invalid_reason", "views_used", "selected_unary_model",
    "graph_enabled", "fit_sha256", "decision_sha256", "edge_table_sha256", "graph_result_sha256",
    "transform_sha256",
]
const EDGE_HEADER = [
    "schema_version", "evidence_scope", "fold_date", "scan_date", "file", "control_class",
    "perturbation", "view", "seed_family", "seed", "target", "execution_direction",
    "left_lobe", "right_lobe", "left_t_nm", "right_t_nm", "gap_nm", "left_segment_id",
    "right_segment_id", "edge_status", "split_reason", "corr_fwd", "corr_bwd",
    "edge_table_sha256", "transform_sha256",
]
const TRANSFORM_HEADER = [
    "schema_version", "evidence_scope", "fold_date", "scan_date", "file", "lobe", "control_class",
    "perturbation", "view", "seed_family", "seed", "target", "channel", "baseline_table_sha256",
    "perturbed_table_sha256", "baseline_row_sha256", "perturbed_row_sha256", "noise_scale_sha256",
    "artifact_path", "transform_sha256",
]
const EXPECTED_HEADER = [
    "schema_version", "evidence_scope", "row_kind", "fold_date", "scan_date", "file", "lobe", "model",
    "control_class", "perturbation", "view", "seed_family", "seed", "target", "channel", "left_lobe",
    "right_lobe", "universe_sha256",
]
const LOBE_HEADER = [
    "schema_version", "evidence_scope", "fold_date", "scan_date", "file", "lobe", "model",
    "control_class", "perturbation", "view", "seed_family", "seed", "target", "baseline_posterior_1",
    "perturbed_posterior_1", "baseline_hard", "perturbed_hard", "baseline_status", "perturbed_status",
    "baseline_reason", "perturbed_reason", "baseline_classified", "perturbed_classified", "posterior_abs_delta",
    "hard_agree", "baseline_row_sha256", "perturbed_row_sha256",
]
const SCAN_HEADER = [
    "schema_version", "evidence_scope", "fold_date", "scan_date", "file", "model", "control_class",
    "perturbation", "view", "seed_family", "seed", "target", "node_count", "hard_agree_count", "loss",
    "baseline_classified", "perturbed_classified", "baseline_coverage", "perturbed_coverage",
    "posterior_pair_count", "posterior_correlation", "posterior_correlation_status", "abstention_count",
    "failure_count",
]
const TOPOLOGY_HEADER = [
    "schema_version", "evidence_scope", "fold_date", "scan_date", "file", "control_class", "perturbation",
    "view", "seed_family", "seed", "target", "left_lobe", "right_lobe", "baseline_edge_status",
    "perturbed_edge_status", "baseline_split_reason", "perturbed_split_reason", "baseline_corr_fwd",
    "perturbed_corr_fwd", "baseline_corr_bwd", "perturbed_corr_bwd", "baseline_left_segment_id",
    "baseline_right_segment_id", "perturbed_left_segment_id", "perturbed_right_segment_id", "equivalent",
]
const BOOTSTRAP_HEADER = [
    "schema_version", "evidence_scope", "control_class", "perturbation", "view", "seed_family", "seed",
    "target", "bootstrap_seed", "date_count", "scan_count", "candidate_mean_loss", "c1_mean_loss",
    "mean_loss_difference",
]
const REASON_HEADER = [
    "schema_version", "evidence_scope", "level", "model", "control_class", "perturbation", "view",
    "seed_family", "seed", "target", "status", "reason", "count",
]
const SUMMARY_HEADER = [
    "schema_version", "evidence_scope", "control_class", "perturbation", "view", "seed_family", "seed",
    "target", "gate_role", "scan_count", "node_count", "candidate_mean_loss", "c1_mean_loss",
    "mean_loss_difference", "bootstrap_upper", "candidate_coverage", "c1_coverage", "hard_agreement",
    "posterior_correlation_defined_scans", "reversal_max_abs_delta", "topology_equivalent", "computed_status",
    "publication_status", "reason",
]

const PATCH_HEADER = ["schema_version", "patch_sha256", "channel", "pixel_count", "values_base64"]
const NOISE_HEADER = [
    "schema_version", "fold_date", "channel", "training_file", "training_scan_mad", "median_scan_mad",
    "target_mad", "heldout_excluded", "scale_sha256", "provenance_sha256", "evidence_scope",
]
const EXECUTION_FILES = (
    "node_runs.tsv", "edge_runs.tsv", "transform_runs.tsv", "artifact-index.tsv", "execution-receipt.toml",
)
const UNIVERSE_FILES = (
    "folds.tsv", "patch_keys.tsv", "perturbations.tsv", "receipt.toml", "scan_universe.tsv", "seeds.tsv",
    "shards.tsv",
)
const OUTPUT_NAMES = (
    "expected_rows.tsv", "lobe_pairs.tsv", "scan_metrics.tsv", "topology_pairs.tsv", "bootstrap.tsv",
    "reason_counts.tsv", "summary.tsv",
)
const EXECUTION_RECEIPT_KEYS = Set([
    "schema", "schema_version", "status", "evidence_scope", "stage", "attempt_id", "execution_site",
    "julia_threads", "feature_sha256", "candidate_config_sha256", "model_config_sha256",
    "universe_receipt_sha256", "expected_rows_sha256", "node_runs_sha256", "edge_runs_sha256",
    "transform_runs_sha256", "artifact_index_sha256", "producer_source_sha256", "producer_command_sha256",
    "producer_registry_sha256",
])
const EXECUTION_COMMAND_TOKENS = (
    "--root=\$ROOT",
    "--features=\$FEATURES",
    "--feature-sha256=\$FEATURE_SHA256",
    "--candidate-config=config/unit_assignment_structured_candidate.toml",
    "--model-config=config/unit_assignment_structured_model.toml",
    "--universe-dir=\$UNIVERSE_DIR",
    "--execution-dir=\$EXECUTION_DIR",
    "--out-dir=\$OUT_DIR",
)
const EXECUTION_CAMPAIGN = "structured-stability"
const EXECUTION_MODEL = "C1-meta"
const EXECUTION_GENERATED_AT = "1970-01-01T00:00:00Z"
const PATCH_MATERIAL_SUFFIX = "patch-rows.tsv"
const NOISE_MATERIAL_SUFFIX = "noise-scales.tsv"
const PUBLIC_MATERIAL_SUFFIX = "public"
const PUBLIC_TODO3_NAMES = ("forward.tsv", "backward.tsv", "forward-receipt.toml", "backward-receipt.toml")
const PUBLIC_TODO10_NAMES = ("edge_observations.tsv", "node_segments.tsv", "receipt.toml")
const PUBLIC_TODO11_NAMES = (
    "models.tsv", "partitions.tsv", "ess.tsv", "starts.tsv", "scores.tsv", "shuffle.tsv", "bootstrap.tsv",
    "fitted_edge_models.tsv", "fitted_unary_nodes.tsv", "fitted_edge_transforms.tsv", "receipt.toml",
)
const PUBLIC_FEATURE_REQUIRED_COLUMNS = ("file", "N", "lobe", "amplitude", "t_nm", "u_nm")
const PUBLIC_FEATURE_ALLOWED_COLUMNS = Set(String.(collect(StructuredUniverse.FEATURE_COLUMNS)))
const PUBLIC_TODO12_REASON_MAP = Dict(
    "graph_inference" => "graph_inference",
    "unary_only" => "unary_only",
    "factor_failure" => "factor_failure",
    "unary_failure" => "unary_failure",
    "reversal_mismatch" => "reversal_mismatch",
    "internal_inference" => "internal_inference",
    "density_failure" => "density_failure",
    "dp_failure" => "dp_failure",
    "score_failure" => "score_failure",
    "internal" => "internal_inference",
    "todo11_fail" => "todo11_fail",
    "internal_validation" => "internal_validation",
    "reference_mismatch" => "todo12_binding",
    "schema_mismatch" => "schema",
    "hash_mismatch" => "hash",
    "status_mismatch" => "node_status",
    "row_count" => "completeness",
    "model_mismatch" => "model",
    "transform_mismatch" => "transform_mismatch",
    "transform_coverage" => "completeness",
    "unary_coverage" => "completeness",
    "topology_mismatch" => "topology",
    "order_mismatch" => "schema",
    "invalid_hash" => "hash",
    "invalid_root" => "path_schema",
    "noncanonical_path" => "path_schema",
    "symlink_rejected" => "symlink_path",
    "missing_input" => "missing_input",
    "invalid_encoding" => "invalid_encoding",
    "carriage_return" => "schema",
    "missing_final_lf" => "schema",
    "input_changed" => "stale_input",
    "internal_snapshot" => "internal_validation",
    "empty_input" => "missing_input",
    "invalid_key" => "schema",
    "date_mismatch" => "schema",
    "path_ambiguity" => "path_schema",
    "upstream_status" => "todo10_validation",
    "receipt_binding" => "transform_binding",
    "segment_mismatch" => "topology",
    "boundary_mismatch" => "topology",
    "split_reason" => "topology",
)

function _execution_command()
    tokens = collect(String.(EXECUTION_COMMAND_TOKENS))
    parsed = _parse_arguments(tokens)
    Set(keys(parsed)) == VALUE_FLAGS ||
        error("execution command token list does not cover the CLI contract")
    return "julia --project=. test/audit_structured_unit_assignment_stability.jl " * join(tokens, " ")
end
const ALLOWED_STATUSES = Set(("ok", "abstained", "failed", "eligible", "split"))
const ALLOWED_NODE_REASONS = Set(("none", "all_views_missing", "nonfinite_posterior", "numerical_failure"))
const ALLOWED_EDGE_REASONS = Set(("none", "gap_out_of_range", "missing_forward_key", "missing_backward_key",
                                 "nonfinite_fwd", "nonfinite_bwd", "zero_variance_fwd", "zero_variance_bwd",
                                 "corr_out_of_range_fwd", "corr_out_of_range_bwd"))
const ALLOWED_EXECUTION_DIRECTIONS = Set(("forward", "backward"))
const ALLOWED_METADATA_VALUES = Set(("synthetic_fixture", "real_viper", "stability_audit", "PASS", "C1-meta"))
const ALLOWED_REASON_CODES = Set((
    "path_invalid", "path_substitution", "path_escape", "symlink_path", "hardlink_path", "special_file",
    "changed_during_read", "stale_input", "missing_path", "missing_input", "schema", "header_mismatch",
    "invalid_encoding", "invalid_integer", "invalid_number", "noncanonical_number", "hash", "feature_hash",
    "dependency_hash_mismatch", "source_collision", "scope", "execution_status", "stage", "execution_site",
    "execution_threads", "attempt_id", "fileset", "artifact_path", "artifact_index", "artifact_role",
    "artifact_closure", "identity", "control_set", "universe_binding", "universe_receipt", "axis",
    "material_reference", "noise_material", "missing_material", "duplicate", "completeness", "expected_rows",
    "fold_binding", "model", "node_status", "node_reason", "node_schema", "edge_status", "edge_schema",
    "execution_direction", "todo10_projection", "todo12_binding", "todo12_projection", "hash_binding",
    "transform_binding", "transform_projection", "reversal_inference", "reversal_mismatch", "bootstrap",
    "internal_error", "unregistered_real_producer", "synthetic_fixture_only",
    "stability_gate_failure", "publication_blocked", "firewall", "forbidden_cli", "missing_bindings",
    "todo10_validation", "todo12_validation", "public_material", "transform_mismatch", "path_schema",
    "topology", "feature_key_mismatch", "duplicate_key", "empty_hidden_directory", "missing_output_parent",
    "atomic_publication_unsupported", "publication_collision", "atomic_install_failed", "fsync_failed",
    "cross_filesystem_publication", "interrupted_publication", "publication_identity",
    "publication_validation_failed", "report_schema", "invalid_failpoint", "invalid_expected_identity",
    "publication_blocked", "graph_inference", "todo11_fail", "unary_only", "internal_validation",
    "none", "all_views_missing", "nonfinite_posterior", "numerical_failure", "gap_out_of_range",
    "missing_forward_key", "missing_backward_key", "nonfinite_fwd", "nonfinite_bwd", "zero_variance_fwd",
    "zero_variance_bwd", "corr_out_of_range_fwd", "corr_out_of_range_bwd", "ok", "partition_blocked",
    "partition_failed", "partition_skipped", "nonpositive_date", "bootstrap_not_positive",
    "shuffle_control_failed", "reversal_inequivalent", "conditional_mechanism_not_admitted",
    "insufficient_dates", "integrity_failure", "nonfinite_endpoint_predictor", "residualizer_unavailable",
    "factor_failure", "unary_failure", "reversal_mismatch", "internal_inference", "density_failure",
    "dp_failure", "score_failure",
))

const _SOURCE_REASON_VOCABULARY = Set(ALLOWED_REASON_CODES)
isempty(setdiff(_SOURCE_REASON_VOCABULARY, ALLOWED_REASON_CODES)) ||
    error("structured stability reason vocabulary is not closed")
isempty(setdiff(Set(values(PUBLIC_TODO12_REASON_MAP)), ALLOWED_REASON_CODES)) ||
    error("public Todo 12 reason vocabulary is not closed")

struct AuditError <: Exception
    status::Symbol
    reason::String
    message::String
end
Base.showerror(io::IO, error::AuditError) =
    print(io, error.status, " [", error.reason, "]: ", error.message)

struct PublicationError <: Exception
    state::String
    reason::String
    message::String
end
Base.showerror(io::IO, error::PublicationError) =
    print(io, error.state, " [", error.reason, "]: ", error.message)

struct Binding
    role::String
    relative::String
    path::String
    bytes::Vector{UInt8}
    sha256::String
    device::UInt64
    inode::UInt64
    size::Int64
    mtime::Float64
    ctime::Float64
end

"""One immutable public predecessor surface for one declared control axis."""
struct PublicAxisMaterial
    axis_key::Tuple
    todo3_forward::String
    todo3_backward::String
    todo3_forward_receipt::String
    todo3_backward_receipt::String
    todo10_edge::String
    todo10_node::String
    todo10_receipt::String
    todo11_model::String
    todo11_unary::String
    todo11_transform::String
    todo11_receipt::String
    todo10_sha256::String
    todo11_model_sha256::String
    todo11_unary_sha256::String
    todo11_transform_sha256::String
    todo11_receipt_sha256::String
    report::Union{Nothing,ChainInferenceReport}
    result_sha256::String
    status::String
    reason::String
end

const _SNAPSHOT_CACHE = Ref{Union{Nothing,Dict{String,Binding}}}(nothing)
const _REPORT_PROGRESS = Ref{Union{Nothing,NamedTuple}}(nothing)

struct AuditReport
    status::String
    reason::String
    computed_status::String
    computed_reason::String
    publication_status::String
    publication_reason::String
    selection_eligible::Bool
    diagnostic_only::Bool
    files::Dict{String,Vector{UInt8}}
    source_hashes::Dict{String,String}
    metrics::Dict{String,String}
    bindings::Vector{Binding}
    specs::Vector{ArtifactSpec}
end

_sha(bytes) = bytes2hex(sha256(bytes))
_sha_text(text::String) = _sha(codeunits(text))
_hash_re(value) = occursin(r"^[0-9a-f]{64}$", value)
_fmt(value::Real) = @sprintf("%.17g", Float64(value))
_bool(value::Bool) = value ? "true" : "false"

function _fail(status::Symbol, reason::AbstractString, message::AbstractString)
    String(reason) in _SOURCE_REASON_VOCABULARY ||
        throw(ArgumentError("unregistered structured stability reason: $reason"))
    throw(AuditError(status, String(reason), String(message)))
end
_blocked(reason, message) = _fail(:BLOCKED, reason, message)
_scientific(reason, message) = _fail(:FAIL, reason, message)

function _canonical_root(root::AbstractString)::String
    supplied = String(root)
    isempty(supplied) && _blocked("path_invalid", "root is empty")
    isdir(supplied) || _blocked("path_invalid", "root is not a directory")
    islink(supplied) && _blocked("symlink_path", "root is a symlink")
    absolute = normpath(abspath(supplied))
    realpath(absolute) == absolute || _blocked("path_invalid", "root is not canonical")
    return absolute
end

function _relative_path(root::String, supplied::AbstractString, context::String; directory=false)
    text = String(supplied)
    isempty(text) && _blocked("path_invalid", "$context is empty")
    isabspath(text) && _blocked("path_substitution", "$context must be root-relative")
    occursin('\0', text) && _blocked("path_invalid", "$context contains NUL")
    occursin('\\', text) && _blocked("path_invalid", "$context uses a backslash")
    parts = split(text, '/'; keepempty=true)
    any(part -> part in ("", ".", ".."), parts) &&
        _blocked("path_escape", "$context is not canonical")
    absolute = normpath(joinpath(root, parts...))
    relpath(absolute, root) == text || _blocked("path_substitution", "$context is not root-relative")
    cursor = root
    for part in parts
        cursor = joinpath(cursor, part)
        islink(cursor) && _blocked("symlink_path", "$context contains a symlink")
    end
    directory ? (isdir(absolute) || _blocked("missing_path", "$context is not a directory")) :
        (isfile(absolute) || _blocked("missing_path", "$context is not a regular file"))
    return absolute, text
end

function _source_path(root::String, relative::String, context::String)
    absolute, _ = _relative_path(root, relative, context)
    realpath(absolute) == absolute || _blocked("path_invalid", "$context is not canonical")
    return absolute
end

function _snapshot(root::String, relative::String, role::String; expected=nothing, fresh=false)::Binding
    absolute, canonical = _relative_path(root, relative, role)
    cache = _SNAPSHOT_CACHE[]
    if !fresh && cache !== nothing && haskey(cache, canonical)
        binding = cache[canonical]
        expected !== nothing && binding.sha256 != expected &&
            _blocked("dependency_hash_mismatch", "$role SHA-256 differs")
        _progress_record_binding!(binding)
        return binding
    end
    islink(absolute) && _blocked("symlink_path", "$role is a symlink")
    st_before = stat(absolute)
    (st_before.mode & 0o170000) == 0o100000 ||
        _blocked("special_file", "$role is not a regular file")
    getproperty(st_before, :nlink) == 1 || _blocked("hardlink_path", "$role is a hardlink")
    bytes = Vector{UInt8}(read(absolute))
    islink(absolute) && _blocked("symlink_path", "$role became a symlink")
    st_after = stat(absolute)
    (st_after.mode & 0o170000) == 0o100000 || _blocked("special_file", "$role became a special file")
    getproperty(st_after, :nlink) == 1 || _blocked("hardlink_path", "$role became a hardlink")
    stable = st_before.device == st_after.device && st_before.inode == st_after.inode &&
             st_before.size == st_after.size && st_before.mtime == st_after.mtime &&
             st_before.ctime == st_after.ctime && length(bytes) == st_after.size
    stable || _blocked("changed_during_read", "$role changed during read")
    digest = _sha(bytes)
    expected !== nothing && digest != expected &&
        _blocked("dependency_hash_mismatch", "$role SHA-256 differs")
    binding = Binding(String(role), canonical, absolute, bytes, digest, UInt64(st_after.device),
                      UInt64(st_after.inode), Int64(st_after.size), Float64(st_after.mtime),
                      Float64(st_after.ctime))
    cache !== nothing && !fresh && (cache[canonical] = binding)
    _progress_record_binding!(binding)
    return binding
end

# `_verify` needs the root which is deliberately not part of the public report.
function _verify(root::String, binding::Binding)
    current = _snapshot(root, binding.relative, binding.role; fresh=true)
    (current.sha256 == binding.sha256 && current.device == binding.device &&
     current.inode == binding.inode && current.size == binding.size &&
     current.mtime == binding.mtime && current.ctime == binding.ctime) ||
        _blocked("stale_input", "$(binding.role) changed after binding")
    return nothing
end

function _progress_init!(scope::String="unknown")
    _REPORT_PROGRESS[] = (scope=scope, files=nothing, bindings=Binding[], specs=ArtifactSpec[],
                          metrics=Dict{String,String}())
    return nothing
end

function _progress_binding_role(binding::Binding)
    startswith(binding.role, "source ") && return "source"
    binding.role == "candidate config" && return "config"
    return "input"
end

function _progress_record_binding!(binding::Binding; role=nothing)
    progress = _REPORT_PROGRESS[]
    progress === nothing && return nothing
    bindings = copy(progress.bindings)
    existing = findfirst(item -> item.relative == binding.relative, bindings)
    if existing === nothing
        push!(bindings, binding)
    else
        bindings[existing].sha256 == binding.sha256 ||
            _blocked("source_collision", "one canonical path has conflicting snapshots")
    end
    specs = copy(progress.specs)
    spec_role = role === nothing ? _progress_binding_role(binding) : String(role)
    spec_index = findfirst(item -> item.path == binding.relative, specs)
    if spec_index === nothing
        push!(specs, ArtifactSpec(spec_role, binding.relative))
    elseif specs[spec_index].role != spec_role
        # A source/config role is authoritative over the generic input role.
        if specs[spec_index].role == "input" && spec_role in ("source", "config")
            specs[spec_index] = ArtifactSpec(spec_role, binding.relative)
        elseif spec_role != "input"
            _blocked("artifact_role", "acquired binding role differs")
        end
    end
    _REPORT_PROGRESS[] = (scope=progress.scope, files=progress.files, bindings=bindings,
                          specs=specs, metrics=progress.metrics)
    return nothing
end

function _progress_update!(; scope=nothing, files=nothing, bindings=nothing, specs=nothing, metrics=nothing)
    progress = _REPORT_PROGRESS[]
    progress === nothing && return nothing
    _REPORT_PROGRESS[] = (scope=scope === nothing ? progress.scope : String(scope),
                          files=files === nothing ? progress.files : files,
                          bindings=bindings === nothing ? progress.bindings : bindings,
                          specs=specs === nothing ? progress.specs : specs,
                          metrics=metrics === nothing ? progress.metrics : metrics)
    return nothing
end

function _parse_tsv(bytes::Vector{UInt8}, header::Vector{String}, context::String)
    text = try
        String(copy(bytes))
    catch error
        _blocked("invalid_encoding", "$context is not UTF-8: $(sprint(showerror, error))")
    end
    isvalid(text) || _blocked("invalid_encoding", "$context is not UTF-8")
    startswith(text, '\ufeff') && _blocked("schema", "$context has a BOM")
    occursin('\r', text) && _blocked("schema", "$context uses CR line endings")
    endswith(text, '\n') || _blocked("schema", "$context lacks a final LF")
    lines = split(text, '\n'; keepempty=true)
    pop!(lines) == "" || _blocked("schema", "$context termination is invalid")
    isempty(lines) && _blocked("schema", "$context is empty")
    split(first(lines), '\t'; keepempty=true) == header ||
        _blocked("header_mismatch", "$context header differs")
    rows = Vector{Vector{String}}()
    for (row_number, line) in enumerate(lines[2:end])
        isempty(line) && _blocked("schema", "$context has an empty row $row_number")
        fields = String.(split(line, '\t'; keepempty=true))
        length(fields) == length(header) ||
            _blocked("schema", "$context row $row_number has the wrong width")
        any(occursin('\0', field) for field in fields) &&
            _blocked("schema", "$context contains NUL")
        push!(rows, fields)
    end
    return rows
end

function _parse_toml(bytes::Vector{UInt8}, context::String)
    text = try
        String(copy(bytes))
    catch error
        _blocked("invalid_encoding", "$context is not UTF-8: $(sprint(showerror, error))")
    end
    isvalid(text) || _blocked("invalid_encoding", "$context is not UTF-8")
    occursin('\r', text) && _blocked("schema", "$context uses CR line endings")
    endswith(text, '\n') || _blocked("schema", "$context lacks a final LF")
    try
        return TOML.parse(text)
    catch error
        _blocked("schema", "$context is not valid TOML: $(sprint(showerror, error))")
    end
end

function _exact_keys(table::AbstractDict, expected::Set{String}, context::String)
    Set(String.(collect(keys(table)))) == expected ||
        _blocked("schema", "$context keys differ")
    return nothing
end

function _string(table, key::String, context::String)::String
    value = get(table, key, nothing)
    value isa String || _blocked("schema", "$context.$key is not a string")
    return String(value)
end

function _integer(table, key::String, context::String)::Int
    value = get(table, key, nothing)
    value isa Integer && !(value isa Bool) || _blocked("schema", "$context.$key is not an integer")
    try
        return Int(value)
    catch
        _blocked("schema", "$context.$key is outside Int range")
    end
end

function _canonical_int(text::String, context::String; positive=false)::Int
    expression = positive ? r"^[1-9][0-9]*$" : r"^(0|[1-9][0-9]*)$"
    occursin(expression, text) || _blocked("invalid_integer", "$context is not canonical")
    value = tryparse(Int, text)
    value === nothing && _blocked("invalid_integer", "$context is outside Int range")
    return value
end

function _canonical_float(text::String, context::String; allow_na=false)::Union{Nothing,Float64}
    allow_na && text == "NA" && return nothing
    value = tryparse(Float64, text)
    value === nothing && _blocked("invalid_number", "$context is not numeric")
    isfinite(value) || _blocked("invalid_number", "$context is not finite")
    _fmt(value) == text || _blocked("noncanonical_number", "$context is not canonical")
    return value
end

function _public_float(text::String, context::String; allow_na=false)::Union{Nothing,Float64}
    allow_na && text == "NA" && return nothing
    value = tryparse(Float64, text)
    value === nothing && _blocked("invalid_number", "$context is not numeric")
    isfinite(value) || _blocked("invalid_number", "$context is not finite")
    @sprintf("%.7g", value) == text || _blocked("noncanonical_number", "$context is not canonical public wire text")
    return value
end

function _tsv(header::Vector{String}, rows)::Vector{UInt8}
    io = IOBuffer()
    println(io, join(header, '\t'))
    for row in rows
        println(io, join(String.(row), '\t'))
    end
    return take!(io)
end

function _source_authorities(root::String)::Tuple{Vector{Binding},Vector{ArtifactSpec}}
    entries = Tuple{String,Union{Nothing,String}}[
        ("config/unit_assignment_structured_model.toml", MODEL_CONFIG_SHA256),
        ("config/chitosan.toml", CHITOSAN_CONFIG_SHA256),
        ("test/audit_structured_unit_assignment_stability.jl", nothing),
        ("test/lib/structured_assignment/evidence.jl", EVIDENCE_SOURCE_SHA256),
        ("test/lib/structured_assignment/universe.jl", UNIVERSE_SOURCE_SHA256),
        ("test/lib/structured_assignment/edge_features.jl", EDGE_SOURCE_SHA256),
        ("test/lib/structured_assignment/chain_inference.jl", CHAIN_SOURCE_SHA256),
        ("test/lib/structured_assignment/edge_model.jl", EDGE_MODEL_SHA256),
    ]
    append!(entries, T8_SOURCE_ENTRIES)
    bindings = Binding[]
    seen = Set{String}()
    for (relative, expected) in entries
        relative in seen && _blocked("source_collision", "authority appears twice: $relative")
        push!(seen, relative)
        push!(bindings, _snapshot(root, relative, "source $relative"; expected))
    end
    specs = ArtifactSpec[ArtifactSpec("source", binding.relative) for binding in bindings]
    return bindings, specs
end

function _config_binding(root::String)::Binding
    return _snapshot(root, "config/unit_assignment_structured_candidate.toml", "candidate config";
                     expected=CANDIDATE_CONFIG_SHA256)
end

function _artifact_relative(execution_relative::String, suffix::String)::String
    execution_relative == "." && _blocked("path_schema", "execution directory cannot be root")
    startswith(suffix, "/") && _blocked("path_schema", "execution artifact suffix is absolute")
    occursin(r"(^|/)\.\.(/|$)", suffix) && _blocked("path_escape", "execution artifact suffix escapes")
    return "$execution_relative/artifacts/$suffix"
end

function _safe_material_file(root::String, execution_relative::String, suffix::String, role::String)
    relative = _artifact_relative(execution_relative, suffix)
    return _snapshot(root, relative, role)
end

function _axis(control_class, perturbation, view, seed_family, seed, target)
    return (control_class=String(control_class), perturbation=String(perturbation), view=String(view),
            seed_family=String(seed_family), seed=String(seed), target=String(target))
end

function _controls(frozen::FrozenUniverse)
    controls = NamedTuple[]
    push!(controls, _axis("baseline", "unperturbed", "NA", "NA", "NA", "NA"))
    for view in frozen.views
        push!(controls, _axis("unary_view_drop", "channel_drop", view, "NA", "NA", "unary_view"))
    end
    for perturbation in ("shift_t_minus", "shift_t_plus", "shift_u_minus", "shift_u_plus",
                         "contrast_low", "contrast_high")
        push!(controls, _axis("image_perturbation", perturbation, "NA", "NA", "NA",
                              "both_patch_channels"))
    end
    for seed in 0:499
        push!(controls, _axis("image_perturbation", "noise", "per_channel", "noise", string(seed),
                              "both_patch_channels"))
    end
    for seed in 0:4
        push!(controls, _axis("restart", "seed_restart", "NA", "restart", string(seed),
                              "optimizer_restart"))
    end
    push!(controls, _axis("patch_channel_drop", "patch_channel_drop", "forward_patch", "NA", "NA",
                          "forward_patch"))
    push!(controls, _axis("patch_channel_drop", "patch_channel_drop", "backward_patch", "NA", "NA",
                          "backward_patch"))
    push!(controls, _axis("reversal", "reversal", "NA", "NA", "NA", "chain_t"))
    return controls
end

_axis_key(axis) = (axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed, axis.target)
_axis_sort(axis) = _axis_key(axis)

function _axis_dir(axis)
    return join((axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed, axis.target), "/")
end

function _dynamic_tsv(bytes::Vector{UInt8}, context::String)
    text = try
        String(copy(bytes))
    catch error
        _blocked("invalid_encoding", "$context is not UTF-8: $(sprint(showerror, error))")
    end
    isempty(text) && _blocked("schema", "$context is empty")
    first_line = first(split(text, '\n'; keepempty=true))
    header = String.(split(first_line, '\t'; keepempty=true))
    return header, _parse_tsv(bytes, header, context)
end

function _feature_finite(text::String, context::String)::Float64
    value = tryparse(Float64, text)
    value === nothing && _blocked("invalid_number", "$context is not numeric")
    isfinite(value) || _blocked("invalid_number", "$context is not finite")
    return value
end

function _parse_feature_t(root::String, relative::String)
    binding = _snapshot(root, relative, "feature input")
    header, rows = _dynamic_tsv(binding.bytes, "feature input")
    length(header) == length(unique(header)) || _blocked("schema", "feature header repeats a column")
    all(column -> column in PUBLIC_FEATURE_ALLOWED_COLUMNS, header) ||
        _blocked("schema", "feature header contains a forbidden column")
    all(column -> column in header for column in PUBLIC_FEATURE_REQUIRED_COLUMNS) ||
        _blocked("schema", "feature header lacks a required public column")
    indices = Dict(column => findfirst(==(column), header)::Int for column in PUBLIC_FEATURE_REQUIRED_COLUMNS)
    values = Dict{Tuple{String,Int},Float64}()
    u_values = Dict{Tuple{String,Int},Float64}()
    metadata = Dict{Tuple{String,Int},NamedTuple}()
    for (number, row) in enumerate(rows)
        file = row[indices["file"]]
        StructuredUniverse.parse_scan_date(file)
        n_text = row[indices["N"]]
        n_lobes = _canonical_int(n_text, "feature row $number N"; positive=true)
        lobe_text = row[indices["lobe"]]
        lobe = _canonical_int(lobe_text, "feature row $number lobe"; positive=true)
        lobe <= n_lobes || _blocked("feature_key_mismatch", "feature lobe exceeds declared N")
        key = (file, lobe)
        haskey(values, key) && _blocked("duplicate_key", "feature key repeats")
        t_text = row[indices["t_nm"]]
        u_text = row[indices["u_nm"]]
        amplitude_text = row[indices["amplitude"]]
        values[key] = _feature_finite(t_text, "feature row $number t_nm")
        u_values[key] = _feature_finite(u_text, "feature row $number u_nm")
        _feature_finite(amplitude_text, "feature row $number amplitude")
        metadata[key] = (file=file, n=n_text, lobe=lobe_text, t_nm=t_text,
                         u_nm=u_text, amplitude=amplitude_text)
    end
    return binding, values, u_values, metadata
end

function _observed_order(frozen::FrozenUniverse, t_values)
    result = Dict{String,Vector{Int}}()
    for scan in frozen.scans
        lobes = [key.lobe for key in frozen.keys if key.file == scan.file]
        sort!(lobes; by=lobe -> (t_values[(scan.file, lobe)], lobe))
        result[scan.file] = lobes
    end
    return result
end

function _involution(order::Vector{Int})
    result = Dict{Int,Int}()
    for index in eachindex(order)
        result[order[index]] = order[length(order) - index + 1]
    end
    return result
end

function _universe_sha(root::String, universe_relative::String)
    binding = _snapshot(root, joinpath(universe_relative, "scan_universe.tsv") |> x -> replace(x, '\\' => '/'),
                        "universe scan table")
    return binding.sha256, binding
end

function _parse_execution_receipt(root::String, execution_relative::String; real_only=false)
    receipt_binding = _snapshot(root, "$execution_relative/execution-receipt.toml", "execution receipt")
    receipt = _parse_toml(receipt_binding.bytes, "execution receipt")
    _exact_keys(receipt, EXECUTION_RECEIPT_KEYS, "execution receipt")
    _string(receipt, "schema", "execution receipt") == "stmfit-structured-stability-execution-v2" ||
        _blocked("schema", "execution receipt schema differs")
    _integer(receipt, "schema_version", "execution receipt") == 2 ||
        _blocked("schema", "execution receipt version differs")
    _string(receipt, "status", "execution receipt") == "PASS" ||
        _blocked("execution_status", "execution receipt is not complete")
    scope = _string(receipt, "evidence_scope", "execution receipt")
    scope in ("synthetic_fixture", "real_viper") || _blocked("scope", "execution scope differs")
    real_only && scope != "real_viper" && _blocked("scope", "real preflight scope differs")
    _string(receipt, "stage", "execution receipt") == "stability_audit" ||
        _blocked("stage", "execution stage differs")
    _string(receipt, "execution_site", "execution receipt") == scope ||
        _blocked("execution_site", "execution site differs")
    _integer(receipt, "julia_threads", "execution receipt") == 1 ||
        _blocked("execution_threads", "execution thread count differs")
    attempt = _string(receipt, "attempt_id", "execution receipt")
    occursin(r"^t14-[A-Za-z0-9_.-]{1,48}$", attempt) || _blocked("attempt_id", "attempt id is not canonical")
    scope in ("synthetic_fixture", "real_viper") || _blocked("scope", "execution scope differs")
    for key in ("feature_sha256", "candidate_config_sha256", "model_config_sha256", "universe_receipt_sha256",
                "expected_rows_sha256", "node_runs_sha256", "edge_runs_sha256", "transform_runs_sha256",
                "artifact_index_sha256", "producer_source_sha256", "producer_command_sha256",
                "producer_registry_sha256")
        _hash_re(_string(receipt, key, "execution receipt")) ||
            _blocked("hash", "execution receipt $key is not a lowercase SHA-256")
    end
    return receipt_binding, receipt, scope
end

function _execution_shape(root::String, execution_relative::String)
    execution_path, _ = _relative_path(root, execution_relative, "execution directory"; directory=true)
    names = sort(readdir(execution_path))
    sort(vcat(collect(EXECUTION_FILES), ["artifacts"])) == names ||
        _blocked("fileset", "execution directory file set is not exact")
    artifacts = joinpath(execution_path, "artifacts")
    isdir(artifacts) && !islink(artifacts) || _blocked("artifact_path", "execution artifacts is not a directory")
    for (directory, subdirectories, files) in walkdir(execution_path; follow_symlinks=false)
        for name in subdirectories
            path = joinpath(directory, name)
            islink(path) && _blocked("symlink_path", "execution directory contains a symlink")
            st = lstat(path)
            (st.mode & 0o170000) == 0o40000 || _blocked("special_file", "execution directory contains a special directory")
            isempty(readdir(path)) &&
                _blocked("empty_hidden_directory", "execution directory contains an empty directory")
        end
        for name in files
            path = joinpath(directory, name)
            islink(path) && _blocked("symlink_path", "execution directory contains a symlink")
            st = lstat(path)
            (st.mode & 0o170000) == 0o100000 || _blocked("special_file", "execution directory contains a special file")
            st.nlink == 1 || _blocked("hardlink_path", "execution directory contains a hardlink")
        end
    end
    return execution_path
end

function _public_relative_paths(execution_relative::String, controls)
    public = _artifact_relative(execution_relative, PUBLIC_MATERIAL_SUFFIX)
    paths = String[]
    for axis in controls
        axis_dir = _axis_dir(axis)
        append!(paths, ["$public/todo3/$axis_dir/$name" for name in PUBLIC_TODO3_NAMES])
        append!(paths, ["$public/todo10/$axis_dir/$name" for name in PUBLIC_TODO10_NAMES])
        append!(paths, ["$public/todo11/$axis_dir/$name" for name in PUBLIC_TODO11_NAMES])
    end
    return sort!(paths)
end

function _execution_material_inventory(root::String, execution_relative::String)
    material = _artifact_relative(execution_relative, "")
    material = chop(material; tail=1)
    material_path = joinpath(root, split(material, '/')...)
    isdir(material_path) && !islink(material_path) ||
        _blocked("artifact_path", "execution material root is not a directory")
    root_stat = lstat(material_path)
    (root_stat.mode & 0o170000) == 0o40000 ||
        _blocked("special_file", "execution material root is special")
    isempty(readdir(material_path)) &&
        _blocked("empty_hidden_directory", "execution material root is empty")
    files = Set{String}()
    directories = Set{String}([material])
    for (directory, subdirectories, names) in walkdir(material_path; follow_symlinks=false)
        for name in subdirectories
            path = joinpath(directory, name)
            islink(path) && _blocked("symlink_path", "execution material directory is a symlink")
            st = lstat(path)
            (st.mode & 0o170000) == 0o40000 || _blocked("special_file", "execution material directory is special")
            isempty(readdir(path)) &&
                _blocked("empty_hidden_directory", "execution material has an empty directory")
            push!(directories, replace(relpath(path, root), '\\' => '/'))
        end
        for name in names
            path = joinpath(directory, name)
            islink(path) && _blocked("symlink_path", "execution material is a symlink")
            st = lstat(path)
            (st.mode & 0o170000) == 0o100000 || _blocked("special_file", "execution material is special")
            st.nlink == 1 || _blocked("hardlink_path", "execution material is a hardlink")
            push!(files, replace(relpath(path, root), '\\' => '/'))
        end
    end
    return (files=files, directories=directories)
end

function _recursive_parent_closure(files::Set{String}, artifacts_root::String)
    closure = Set{String}([artifacts_root])
    for file in files
        parent = dirname(file)
        while true
            push!(closure, replace(parent, '\\' => '/'))
            parent == artifacts_root && break
            parent = dirname(parent)
            relpath(parent, artifacts_root) == ".." &&
                _blocked("artifact_closure", "expected material path escapes artifacts root")
        end
    end
    all(artifacts_root in closure for _ in files) ||
        error("structured stability recursive parent closure lost artifacts root")
    return closure
end

function _registry_hash(records)
    rows = sort([join((record.role, record.path, string(record.byte_count), record.sha256), '\t')
                 for record in records])
    return _sha_text(join(rows, "\n") * (isempty(rows) ? "" : "\n"))
end

function _execution_index(root::String, execution_relative::String, frozen::FrozenUniverse,
                          feature_relative::String, feature_sha::String, receipt, controls,
                          source_specs, universe_relative::String)
    index_relative = "$execution_relative/artifact-index.tsv"
    index_binding = _snapshot(root, index_relative, "execution artifact index")
    result = StructuredEvidence.validate_artifact_index(root, index_relative)
    result.valid || _blocked("artifact_index", join(result.reasons, "; "))
    index = try
        StructuredEvidence.parse_artifact_index(index_binding.bytes)
    catch error
        _blocked("artifact_index", sprint(showerror, error))
    end
    index.campaign_id == EXECUTION_CAMPAIGN || _blocked("identity", "execution campaign differs")
    index.model_id == EXECUTION_MODEL || _blocked("identity", "execution model differs")
    index.metadata["repository_root"] == root || _blocked("identity", "execution root differs")
    index.metadata["generated_at"] == EXECUTION_GENERATED_AT ||
        _blocked("identity", "execution generated_at differs")
    index.metadata["exact_command"] == _execution_command() ||
        _blocked("identity", "execution command differs")
    receipt["artifact_index_sha256"] == index_binding.sha256 ||
        _blocked("hash", "execution receipt does not bind artifact index")
    receipt["feature_sha256"] == feature_sha || _blocked("hash", "execution feature hash differs")
    receipt["candidate_config_sha256"] == CANDIDATE_CONFIG_SHA256 ||
        _blocked("hash", "execution candidate hash differs")
    receipt["model_config_sha256"] == MODEL_CONFIG_SHA256 ||
        _blocked("hash", "execution model hash differs")

    records = index.artifacts
    paths = Set(record.path for record in records)
    length(paths) == length(records) || _blocked("artifact_index", "execution index repeats a path")
    expected_relative = _artifact_relative(execution_relative, "expected_rows.tsv")
    material_paths = Set(vcat([expected_relative, _artifact_relative(execution_relative, PATCH_MATERIAL_SUFFIX),
                               _artifact_relative(execution_relative, NOISE_MATERIAL_SUFFIX)],
                              _public_relative_paths(execution_relative, controls)))
    expected_roles = Dict{String,String}()
    expected_roles["config/unit_assignment_structured_candidate.toml"] = "config"
    for spec in source_specs
        haskey(expected_roles, spec.path) && _blocked("artifact_role", "source closure repeats a path")
        expected_roles[spec.path] = "source"
    end
    expected_roles[expected_relative] = "expected_rows"
    for name in ("node_runs.tsv", "edge_runs.tsv", "transform_runs.tsv")
        expected_roles["$execution_relative/$name"] = "output"
    end
    input_paths = Set{String}([feature_relative, "$execution_relative/artifact-index.tsv",
                               "$execution_relative/execution-receipt.toml", _artifact_relative(execution_relative, PATCH_MATERIAL_SUFFIX),
                               _artifact_relative(execution_relative, NOISE_MATERIAL_SUFFIX)])
    append!(input_paths, ["$universe_relative/$name" for name in UNIVERSE_FILES])
    union!(input_paths, Set(_public_relative_paths(execution_relative, controls)))
    for path in input_paths
        haskey(expected_roles, path) && _blocked("artifact_role", "artifact path has two roles: $path")
        expected_roles[path] = "input"
    end
    observed_roles = Dict(record.path => record.role for record in records)
    observed_roles == expected_roles || _blocked("artifact_role", "execution index path-to-role map is not exact")
    material_inventory = _execution_material_inventory(root, execution_relative)
    material_inventory.files == material_paths || _blocked("artifact_closure", "execution material path closure differs")
    artifacts_root = _artifact_relative(execution_relative, "") |> x -> chop(x; tail=1)
    expected_directories = _recursive_parent_closure(material_paths, artifacts_root)
    material_inventory.directories == expected_directories ||
        _blocked("artifact_closure", "execution material directory closure differs")
    for record in records
        _relative_path(root, record.path, "indexed artifact")
        _hash_re(record.sha256) || _blocked("hash", "indexed artifact digest is malformed")
        binding = _snapshot(root, record.path, "indexed $(record.role) $(record.path)")
        binding.sha256 == record.sha256 && binding.size == record.byte_count ||
            _blocked("artifact_index", "indexed artifact bytes differ: $(record.path)")
    end

    StructuredEvidence.role_bundle_sha256(index, "source") == receipt["producer_source_sha256"] ||
        _blocked("hash", "producer source digest is not the exact indexed source bundle")
    StructuredEvidence.command_identity_sha256(root, _execution_command()) == receipt["producer_command_sha256"] ||
        _blocked("hash", "producer command digest is not the exact indexed command")
    _registry_hash(records) == receipt["producer_registry_sha256"] ||
        _blocked("hash", "producer registry digest is not the exact indexed registry")

    expected_binding = _snapshot(root, expected_relative, "execution expected rows")
    expected_binding.sha256 == receipt["expected_rows_sha256"] ||
        _blocked("hash", "execution expected-row hash differs")
    return (index=index, index_binding=index_binding, expected_binding=expected_binding,
            expected_relative=expected_relative, records=records)
end

function _run_bindings(root::String, execution_relative::String, receipt)
    result = Dict{String,Binding}()
    for name in ("node_runs.tsv", "edge_runs.tsv", "transform_runs.tsv")
        binding = _snapshot(root, "$execution_relative/$name", "execution $name")
        result[name] = binding
        receipt_key = Dict("node_runs.tsv" => "node_runs_sha256",
                           "edge_runs.tsv" => "edge_runs_sha256",
                           "transform_runs.tsv" => "transform_runs_sha256")[name]
        binding.sha256 == receipt[receipt_key] || _blocked("hash", "execution $name hash differs")
    end
    return result
end

function _patch_hash(channel::String, values)
    return _sha_text("structured-patch-row-v1\nchannel=$channel\npixels=289\nvalues=" *
                    join((ismissing(value) ? "NA" : _fmt(value) for value in values), ',') * "\n")
end

function _transform_hash(row)
    return _sha_text("structured-transform-v1\n$(row.fold_date)\n$(row.scan_date)\n$(row.file)\n$(row.lobe)\n" *
                    "$(row.control_class)\n$(row.perturbation)\n$(row.view)\n$(row.seed_family)\n$(row.seed)\n" *
                    "$(row.target)\n$(row.channel)\n$(row.baseline_row_sha256)\n$(row.perturbed_row_sha256)\n" *
                    "$(row.noise_scale_sha256)\n")
end

function _parse_patch_material(root::String, execution_relative::String, index, receipt)
    relative = _artifact_relative(execution_relative, PATCH_MATERIAL_SUFFIX)
    binding = _safe_material_file(root, execution_relative, PATCH_MATERIAL_SUFFIX, "patch material")
    any(record -> record.path == relative && record.role == "input", index.artifacts) ||
        _blocked("artifact_role", "patch material is not an indexed input")
    rows = _parse_tsv(binding.bytes, PATCH_HEADER, "patch material")
    patches = Dict{String,Tuple{String,Vector{Union{Missing,Float64}}}}()
    for (number, row) in enumerate(rows)
        row[1] == "1" || _blocked("schema", "patch material version differs")
        channel = row[3]
        channel in ("forward", "backward") || _blocked("schema", "patch channel differs")
        _canonical_int(row[4], "patch material pixel count"; positive=true) == 289 ||
            _blocked("schema", "patch material pixel count differs")
        digest = row[2]
        _hash_re(digest) || _blocked("hash", "patch material hash is malformed")
        encoded = try
            base64decode(row[5])
        catch error
            _blocked("schema", "patch material base64 is invalid: $(sprint(showerror, error))")
        end
        base64encode(encoded) == row[5] || _blocked("schema", "patch material base64 is not canonical")
        text = String(encoded)
        values = Union{Missing,Float64}[]
        for value in split(text, ','; keepempty=true)
            value == "NA" ? push!(values, missing) :
                push!(values, something(_canonical_float(value, "patch material value")))
        end
        length(values) == 289 || _blocked("schema", "patch material row $number has 289 values")
        _patch_hash(channel, values) == digest || _blocked("hash", "patch material row hash differs")
        haskey(patches, digest) && patches[digest] != (channel, values) &&
            _blocked("duplicate", "patch material hash has conflicting bytes")
        patches[digest] = (channel, values)
    end
    isempty(patches) && _blocked("missing_material", "patch material is empty")
    return binding, patches
end

function _noise_provenance(scale)
    return _sha_text("structured-noise-provenance-v1\nfold_date=$(scale.fold_date)\nchannel=$(scale.channel)\n" *
                    "training_files=$(join(scale.training_files, ','))\n" *
                    "scan_mads=$(join(_fmt.(scale.training_scan_mads), ','))\n" *
                    "median_scan_mad=$(_fmt(scale.median_scan_mad))\n" *
                    "target_mad=$(_fmt(scale.target_mad))\nheldout_excluded=true\n")
end

function _matrix(values)
    length(values) == 17^2 || _blocked("schema", "patch material does not contain 289 values")
    # Producer serialization is u_outer,t_inner.  The matrix consumed by the
    # public perturbation API is indexed [u,t], hence this transpose is part of
    # the frozen wire contract rather than a presentation choice.
    return permutedims(reshape(values, 17, 17))
end

function _encode_matrix(matrix)
    size(matrix) == (17, 17) || _blocked("schema", "physical patch is not 17 by 17")
    return collect(vec(permutedims(matrix)))
end

function _patch_geometry_invariant()
    wire = Float64[1000u + 10t + 0.001u * t for u in 1:17, t in 1:17]
    serialized = collect(vec(wire))
    physical = _matrix(serialized)
    _encode_matrix(physical) == serialized ||
        error("structured stability patch geometry encode/decode invariant failed")
    physical[1, 2] != physical[2, 1] ||
        error("structured stability patch geometry collapsed u and t axes")
    shifted_t = StructuredUniverse.shift_patch(physical, 0.02, 0.0)
    shifted_u = StructuredUniverse.shift_patch(physical, 0.0, 0.02)
    shifted_t != shifted_u ||
        error("structured stability patch geometry confused t and u shifts")
    reversed_t = physical[:, end:-1:1]
    _matrix(_encode_matrix(reversed_t)) == reversed_t ||
        error("structured stability patch geometry reversed the wrong axis")
    _matrix(_encode_matrix(reversed_t)) != physical[end:-1:1, :] ||
        error("structured stability patch geometry reversal is not t-local")
    return nothing
end

_patch_geometry_invariant()

function _material_rows(root::String, execution_relative::String, index, frozen, controls,
                        t_values, patches, patch_binding, noise_binding, noise_rows, execution_scope)
    transform_binding = _snapshot(root, "$execution_relative/transform_runs.tsv", "transform runs")
    transform_rows = _parse_tsv(transform_binding.bytes, TRANSFORM_HEADER, "transform runs")
    by_digest = Dict{String,Tuple{String,Vector{Union{Missing,Float64}}}}(patches)
    transform_map = Dict{Tuple{String,String,String,String,String,String,String,String,Int,String},NamedTuple}()
    observed_scope = nothing
    for (number, row) in enumerate(transform_rows)
        row[1] == "1" || _blocked("schema", "transform row $number version differs")
        row_scope = row[2]
        row_scope in ("synthetic_fixture", "real_viper") || _blocked("scope", "transform scope differs")
        observed_scope === nothing ? (observed_scope = row_scope) : observed_scope == row_scope ||
            _blocked("scope", "transform rows mix execution scopes")
        row_scope == execution_scope || _blocked("scope", "transform row scope differs")
        fold = row[3]; scan_date = row[4]; file = row[5]
        fold == scan_date || _blocked("fold_binding", "transform rows must use one scan fold, not fold-by-scan expansion")
        lobe = _canonical_int(row[6], "transform row lobe"; positive=true)
        axis = _axis(row[7], row[8], row[9], row[10], row[11], row[12])
        axis in controls || _blocked("axis", "transform row uses an undeclared control")
        channel = row[13]
        channel in ("forward", "backward") || _blocked("axis", "transform channel differs")
        row[14] == patch_binding.sha256 && row[15] == patch_binding.sha256 ||
            _blocked("hash", "transform table binding differs")
        row[19] == _artifact_relative(execution_relative, PATCH_MATERIAL_SUFFIX) ||
            _blocked("path_schema", "transform material path differs")
        base_digest, pert_digest = row[16], row[17]
        haskey(by_digest, base_digest) && haskey(by_digest, pert_digest) ||
            _blocked("material_reference", "transform references an absent patch row")
        by_digest[base_digest][1] == channel && by_digest[pert_digest][1] == channel ||
            _blocked("material_reference", "transform channel does not match patch bytes")
        _hash_re(row[18]) || (row[18] == "NA" || _blocked("hash", "noise scale hash is malformed"))
        key = (fold, scan_date, file, axis.control_class, axis.perturbation, axis.view,
               axis.seed_family, axis.seed, lobe, channel)
        haskey(transform_map, key) && _blocked("duplicate", "transform key repeats")
        named = (schema_version=row[1], evidence_scope=row_scope, fold_date=fold, scan_date=scan_date,
                 file=file, lobe=lobe, axis=axis, channel=channel, baseline_row_sha256=base_digest,
                 perturbed_row_sha256=pert_digest, noise_scale_sha256=row[18], artifact_path=row[19],
                 transform_sha256=row[20])
        _transform_hash((fold_date=fold, scan_date=scan_date, file=file, lobe=row[6],
                        control_class=axis.control_class, perturbation=axis.perturbation, view=axis.view,
                        seed_family=axis.seed_family, seed=axis.seed, target=axis.target, channel=channel,
                        baseline_row_sha256=base_digest, perturbed_row_sha256=pert_digest,
                        noise_scale_sha256=row[18])) == row[20] ||
            _blocked("hash", "transform row hash differs")
        transform_map[key] = named
    end
    observed_scope == execution_scope || _blocked("scope", "transform execution scope differs")

    # Baseline bytes are a global keyed map.  Fold rows are merely references to
    # that map; they are not a second fold-by-scan materialization of patches.
    baseline_by_key = Dict{Tuple{String,Int,String},String}()
    for item in values(transform_map)
        item.axis.control_class == "baseline" || continue
        key = (item.file, item.lobe, item.channel)
        previous = get(baseline_by_key, key, item.baseline_row_sha256)
        previous == item.baseline_row_sha256 ||
            _blocked("material_reference", "baseline patch bytes vary for one physical key")
        baseline_by_key[key] = previous
    end
    expected_baseline_keys = Set((key.file, key.lobe, channel)
                                 for key in frozen.keys for channel in ("forward", "backward"))
    length(expected_baseline_keys) == 2 * length(frozen.keys) ||
        error("structured stability baseline channel key invariant failed")
    Set(keys(baseline_by_key)) == expected_baseline_keys ||
        _blocked("material_reference", "global baseline patch map is incomplete or has extras")
    all(haskey(by_digest, digest) for digest in values(baseline_by_key)) ||
        _blocked("material_reference", "global baseline map references absent patch bytes")
    # The public transform functions are used to derive the expected bytes only
    # after the complete global baseline map has been sealed.  No fit,
    # posterior, graph matrix, or topology is reconstructed here.
    scales = Dict{Tuple{String,String},NoiseScale}()
    noise_groups = Dict{Tuple{String,String},Vector{Vector{String}}}()
    expected_noise_keys = Set((fold, channel, scan.file)
                              for fold in unique(scan.date for scan in frozen.scans),
                                  channel in ("forward", "backward"),
                                  scan in frozen.scans if scan.date != fold)
    observed_noise_keys = Set{Tuple{String,String,String}}()
    for (number, row) in enumerate(noise_rows)
        row[1] == "1" || _blocked("schema", "noise row $number version differs")
        row[11] == execution_scope || _blocked("scope", "noise row $number scope differs")
        fold, channel, training_file = row[2], row[3], row[4]
        StructuredUniverse.parse_scan_date(fold)
        channel in ("forward", "backward") || _blocked("noise_material", "noise channel is not allowlisted")
        StructuredUniverse.parse_scan_date(training_file)
        fold != StructuredUniverse.parse_scan_date(training_file) ||
            _blocked("noise_material", "noise training file is held out incorrectly")
        row[8] == "true" || _blocked("noise_material", "noise heldout_excluded is not true")
        noise_key = (fold, channel, training_file)
        noise_key in expected_noise_keys || _blocked("noise_material", "noise group is extra")
        noise_key in observed_noise_keys && _blocked("duplicate_key", "noise group repeats")
        push!(observed_noise_keys, noise_key)
        push!(get!(noise_groups, (fold, channel), Vector{String}[]), row)
    end
    observed_noise_keys == expected_noise_keys ||
        _blocked("completeness", "noise groups are incomplete or have extras")
    for ((fold, channel), rows) in noise_groups
        patch_dict = Dict{LobeKey,Any}()
        for key in frozen.keys
            digest = baseline_by_key[(key.file, key.lobe, channel)]
            patch_dict[key] = _matrix(by_digest[digest][2])
        end
        scale = StructuredUniverse.freeze_noise_scale(frozen, fold, channel, patch_dict)
        expected_files = Set(row[4] for row in rows)
        expected_files == Set(scan.file for scan in frozen.scans if scan.date != fold) ||
            _blocked("noise_material", "noise training-file set differs")
        length(rows) == length(expected_files) || _blocked("noise_material", "noise training file repeats")
        for row in rows
            _canonical_float(row[5], "noise training MAD")
            _canonical_float(row[6], "noise median MAD")
            _canonical_float(row[7], "noise target MAD")
            row[9] == scale.binding_sha256 || _blocked("noise_material", "noise scale hash differs")
            row[10] == _noise_provenance(scale) || _blocked("noise_material", "noise provenance differs")
            training_index = findfirst(==(row[4]), scale.training_files)
            training_index === nothing && _blocked("noise_material", "noise training file is absent")
            expected_mad = scale.training_scan_mads[training_index]
            row[5] == _fmt(expected_mad) && row[6] == _fmt(scale.median_scan_mad) &&
                row[7] == _fmt(scale.target_mad) || _blocked("noise_material", "noise scale values differ")
        end
        scales[(fold, channel)] = scale
    end

    expected = Dict{Tuple{String,String,String,String,String,String,String,String,Int,String},NamedTuple}()
    ordered = _observed_order(frozen, t_values)
    for scan in frozen.scans, axis in controls
        involution = _involution(ordered[scan.file])
        for output_lobe in ordered[scan.file]
            source_lobe = axis.control_class == "reversal" ? involution[output_lobe] : output_lobe
            for channel in ("forward", "backward")
                base_key = (scan.file, source_lobe, channel)
                haskey(baseline_by_key, base_key) || _blocked("material_reference", "baseline transform is incomplete")
                base_matrix = _matrix(by_digest[baseline_by_key[base_key]][2])
                perturbed = if axis.control_class in ("baseline", "unary_view_drop", "restart")
                    base_matrix
                elseif axis.control_class == "patch_channel_drop"
                    axis.view == "forward_patch" && channel == "forward" ? StructuredUniverse.drop_channel(base_matrix) :
                    axis.view == "backward_patch" && channel == "backward" ? StructuredUniverse.drop_channel(base_matrix) : base_matrix
                elseif axis.control_class == "reversal"
                    _matrix(by_digest[baseline_by_key[(scan.file, source_lobe, channel)]][2])[:, end:-1:1]
                elseif axis.perturbation == "noise"
                    haskey(scales, (scan.date, channel)) || _blocked("noise_material", "noise scale is absent")
                    StructuredUniverse.apply_noise_patch(base_matrix, scales[(scan.date, channel)], parse(Int, axis.seed))
                elseif axis.perturbation == "shift_t_minus"
                    StructuredUniverse.shift_patch(base_matrix, -0.02, 0.0)
                elseif axis.perturbation == "shift_t_plus"
                    StructuredUniverse.shift_patch(base_matrix, 0.02, 0.0)
                elseif axis.perturbation == "shift_u_minus"
                    StructuredUniverse.shift_patch(base_matrix, 0.0, -0.02)
                elseif axis.perturbation == "shift_u_plus"
                    StructuredUniverse.shift_patch(base_matrix, 0.0, 0.02)
                elseif axis.perturbation == "contrast_low"
                    StructuredUniverse.contrast_patch(base_matrix, 0.95)
                else
                    StructuredUniverse.contrast_patch(base_matrix, 1.05)
                end
                base_values = _encode_matrix(base_matrix); pert_values = _encode_matrix(perturbed)
                expected_key = (scan.date, scan.date, scan.file, axis.control_class, axis.perturbation,
                                axis.view, axis.seed_family, axis.seed, output_lobe, channel)
                expected[expected_key] = (base=_patch_hash(channel, base_values), perturbed=_patch_hash(channel, pert_values),
                                          noise=axis.perturbation == "noise" ? scales[(scan.date, channel)].binding_sha256 : "NA")
            end
        end
    end
    Set(keys(transform_map)) == Set(keys(expected)) || _blocked("completeness", "transform rows are incomplete or extra")
    for key in keys(expected)
        row = transform_map[key]; wanted = expected[key]
        row.baseline_row_sha256 == wanted.base && row.perturbed_row_sha256 == wanted.perturbed &&
            row.noise_scale_sha256 == wanted.noise && row.transform_sha256 == _transform_hash(
                (fold_date=key[1], scan_date=key[2], file=key[3], lobe=string(key[9]),
                 control_class=key[4], perturbation=key[5], view=key[6], seed_family=key[7], seed=key[8],
                 target=row.axis.target, channel=key[10], baseline_row_sha256=row.baseline_row_sha256,
                 perturbed_row_sha256=row.perturbed_row_sha256, noise_scale_sha256=row.noise_scale_sha256)) ||
            _blocked("transform_mismatch", "materialized transform differs from the public transform")
    end
    return (binding=transform_binding, rows=transform_rows, map=transform_map, expected=expected,
            scales=scales)
end

function _expected_rows_bytes(frozen::FrozenUniverse, controls, t_values, scope::String,
                              universe_sha::String)::Vector{UInt8}
    rows = Vector{Vector{String}}()
    ordered = _observed_order(frozen, t_values)
    for scan in frozen.scans, axis in controls, model in ("C1", "meta")
        for lobe in ordered[scan.file]
            push!(rows, ["1", scope, "node", scan.date, scan.date, scan.file, string(lobe), model,
                         axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed,
                         axis.target, "NA", "NA", "NA", universe_sha])
        end
    end
    for scan in frozen.scans, axis in controls
        order = ordered[scan.file]
        involution = _involution(order)
        for index in 1:max(0, length(order) - 1)
            left, right = order[index:index + 1]
            if axis.control_class == "reversal"
                left, right = involution[right], involution[left]
            end
            push!(rows, ["1", scope, "edge", scan.date, scan.date, scan.file, "NA", "NA",
                         axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed,
                         axis.target, "NA", string(left), string(right), universe_sha])
        end
    end
    for scan in frozen.scans, axis in controls, lobe in ordered[scan.file], channel in ("forward", "backward")
        push!(rows, ["1", scope, "transform", scan.date, scan.date, scan.file, string(lobe), "NA",
                     axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed,
                     axis.target, channel, "NA", "NA", universe_sha])
    end
    sort!(rows; by=Tuple)
    return _tsv(EXPECTED_HEADER, rows)
end

function _row_axis(row, offset::Int)
    return _axis(row[offset], row[offset + 1], row[offset + 2], row[offset + 3], row[offset + 4], row[offset + 5])
end

function _validate_execution_rows(root::String, execution_relative::String, run_bindings, receipt,
                                  frozen, controls, t_values, expected_binding, expected_bytes,
                                  scope, transform, universe_sha)
    expected_rows = _parse_tsv(expected_binding.bytes, EXPECTED_HEADER, "expected rows")
    expected_binding.bytes == expected_bytes ||
        _blocked("expected_rows", "producer expected_rows bytes differ from the independent derivation")
    expected_binding.sha256 == _sha(expected_bytes) ||
        _blocked("expected_rows", "independent expected_rows hash differs")
    expected_nodes = Set{Tuple}(); expected_edges = Set{Tuple}(); expected_transform = Dict{Tuple,Int}()
    for (number, row) in enumerate(expected_rows)
        row[1] == "1" && row[2] == scope || _blocked("expected_rows", "expected row $number scope differs")
        row[3] in ("node", "edge", "transform") || _blocked("expected_rows", "expected row kind differs")
        row[4] == row[5] || _blocked("expected_rows", "expected row $number fold differs")
        row[18] == universe_sha || _blocked("expected_rows", "expected row $number universe hash differs")
        axis = _row_axis(row, 9)
        axis in controls || _blocked("expected_rows", "expected row $number uses an undeclared axis")
        lobe = row[7] == "NA" ? 0 : _canonical_int(row[7], "expected row lobe"; positive=true)
        if row[3] == "node"
            row[8] in ("C1", "meta") || _blocked("expected_rows", "expected node model differs")
            row[15] == "NA" && row[16] == "NA" && row[17] == "NA" ||
                _blocked("expected_rows", "expected node adjacency metadata differs")
            key = (row[4], row[6], lobe, row[8], _axis_key(axis))
            haskey(expected_nodes, key) && _blocked("expected_rows", "expected node key repeats")
            push!(expected_nodes, key)
        elseif row[3] == "edge"
            row[7] == "NA" && row[8] == "NA" && row[15] == "NA" ||
                _blocked("expected_rows", "expected edge node metadata differs")
            left = _canonical_int(row[16], "expected edge left"; positive=true)
            right = _canonical_int(row[17], "expected edge right"; positive=true)
            key = (row[4], row[6], left, right, _axis_key(axis))
            haskey(expected_edges, key) && _blocked("expected_rows", "expected edge key repeats")
            push!(expected_edges, key)
        elseif row[3] == "transform"
            row[8] == "NA" && row[15] in ("forward", "backward") && row[16] == "NA" && row[17] == "NA" ||
                _blocked("expected_rows", "expected transform metadata differs")
            key = (row[4], row[6], lobe, _axis_key(axis), row[15])
            haskey(expected_transform, key) && _blocked("duplicate_key", "expected transform key repeats")
            expected_transform[key] = 1
        else
            _blocked("expected_rows", "expected row kind differs")
        end
    end
    _tsv(EXPECTED_HEADER, expected_rows) == expected_bytes ||
        _blocked("expected_rows", "expected row order is not canonical")
    required_node_count = sum(length([key for key in frozen.keys if key.file == scan.file]) for scan in frozen.scans) * length(controls) * 2
    length(expected_nodes) == required_node_count || _blocked("completeness", "expected node rows are incomplete")
    length(expected_transform) == required_node_count || _blocked("completeness", "expected transform keys are incomplete")
    node_rows = _parse_tsv(run_bindings["node_runs.tsv"].bytes, NODE_HEADER, "node runs")
    edge_rows = _parse_tsv(run_bindings["edge_runs.tsv"].bytes, EDGE_HEADER, "edge runs")
    transform_rows = transform.rows
    node_map = Dict{Tuple,NamedTuple}(); edge_map = Dict{Tuple,NamedTuple}()
    for (number, row) in enumerate(node_rows)
        row[1] == "1" && row[2] == scope || _blocked("schema", "node row scope differs")
        row[3] == row[4] || _blocked("fold_binding", "node fold and scan date differ")
        lobe = _canonical_int(row[6], "node row lobe"; positive=true)
        model = row[7]; model in ("C1", "meta") || _blocked("model", "node model differs")
        axis = _row_axis(row, 8)
        key = (row[3], row[5], lobe, model, _axis_key(axis))
        key in expected_nodes || _blocked("completeness", "node row is not expected")
        haskey(node_map, key) && _blocked("duplicate", "node key repeats")
        posterior = _canonical_float(row[14], "node posterior"; allow_na=true)
        status = row[16]
        status in ("ok", "abstained", "failed") || _blocked("node_status", "node status differs")
        if status == "ok"
            posterior === nothing && _blocked("node_status", "ok node has no posterior")
            row[15] == (posterior >= 0.5 ? "1" : "0") || _blocked("node_status", "hard output is not posterior-derived")
            row[17] == "none" || _blocked("node_status", "ok node has an invalid reason")
        else
            posterior !== nothing && _blocked("node_status", "abstained/failed node has a posterior")
            row[15] == "?" || _blocked("node_status", "abstained/failed node is not unknown")
            row[17] in ALLOWED_NODE_REASONS && row[17] != "none" ||
                _blocked("node_status", "invalid node reason is empty")
        end
        row[17] in ALLOWED_NODE_REASONS || _blocked("node_reason", "node reason is not allowlisted")
        0 <= _canonical_int(row[18], "node views_used") <= 4 || _blocked("node_schema", "views_used is outside the frozen view range")
        row[19] in ("C1", "C2", "NA") || _blocked("node_schema", "selected unary model is not allowlisted")
        if status == "failed"
            row[19] == "NA" || _blocked("node_schema", "failed node must have no selected unary model")
        elseif model == "C1"
            row[19] == "C1" || _blocked("node_schema", "C1 node selected a different unary model")
        else
            row[19] in ("C1", "C2") || _blocked("node_schema", "meta node selected an invalid unary model")
        end
        row[20] in ("true", "false") || _blocked("node_schema", "graph_enabled is not boolean")
        for index in (21, 22, 23, 24, 25)
            _hash_re(row[index]) || _blocked("hash", "node binding hash is malformed")
        end
        tr_key = (row[3], row[4], row[5], axis.control_class, axis.perturbation, axis.view,
                  axis.seed_family, axis.seed, lobe, "forward")
        haskey(transform.map, tr_key) || _blocked("transform_binding", "node transform is absent")
        haskey(transform.map, (row[3], row[4], row[5], axis.control_class, axis.perturbation, axis.view,
                               axis.seed_family, axis.seed, lobe, "backward")) ||
            _blocked("transform_binding", "backward transform is absent")
        row[25] == transform.map[tr_key].transform_sha256 || _blocked("transform_binding", "node transform hash differs")
        node_map[key] = (row=row, posterior=posterior, hard=row[15], status=status, reason=row[17], axis=axis)
    end
    Set(keys(node_map)) == expected_nodes || _blocked("completeness", "node rows are incomplete or extra")

    for (number, row) in enumerate(edge_rows)
        row[1] == "1" && row[2] == scope || _blocked("schema", "edge row scope differs")
        row[3] == row[4] || _blocked("fold_binding", "edge fold and scan date differ")
        axis = _row_axis(row, 6)
        row[12] in ALLOWED_EXECUTION_DIRECTIONS || _blocked("execution_direction", "edge direction is not allowlisted")
        left = _canonical_int(row[13], "edge left lobe"; positive=true)
        right = _canonical_int(row[14], "edge right lobe"; positive=true)
        key = (row[3], row[5], left, right, _axis_key(axis))
        key in expected_edges || _blocked("completeness", "edge row is not expected")
        haskey(edge_map, key) && _blocked("duplicate", "edge key repeats")
        lt = something(_canonical_float(row[15], "edge left t")); rt = something(_canonical_float(row[16], "edge right t"))
        gap = something(_canonical_float(row[17], "edge gap")); gap == rt - lt || _blocked("topology", "edge gap differs")
        status = row[20]; status in ("eligible", "split") || _blocked("edge_status", "edge status differs")
        if status == "eligible"
            row[21] == "none" || _blocked("edge_status", "eligible edge reason differs")
            _canonical_float(row[22], "forward correlation"); _canonical_float(row[23], "backward correlation")
            row[18] == row[19] || _blocked("topology", "eligible edge crosses a segment")
        else
            row[22] == "NA" && row[23] == "NA" || _blocked("edge_status", "split edge exposes a correlation")
            row[18] != row[19] || _blocked("topology", "split edge does not cross a segment")
            row[21] in ALLOWED_EDGE_REASONS && row[21] != "none" || _blocked("edge_status", "split reason differs")
        end
        for index in (24, 25)
            _hash_re(row[index]) || _blocked("hash", "edge binding hash is malformed")
        end
        edge_map[key] = (row=row, axis=axis)
    end
    Set(keys(edge_map)) == expected_edges || _blocked("completeness", "edge rows are incomplete or extra")
    return (nodes=node_map, edges=edge_map, transform=transform)
end

function _correlation(values_a, values_b)
    isempty(values_a) && return nothing
    length(values_a) < 2 && return nothing
    mean_a, mean_b = mean(values_a), mean(values_b)
    numerator = sum((a - mean_a) * (b - mean_b) for (a, b) in zip(values_a, values_b))
    denominator = sqrt(sum((a - mean_a)^2 for a in values_a) * sum((b - mean_b)^2 for b in values_b))
    denominator == 0.0 && return nothing
    return numerator / denominator
end

function _type7(values, probability::Float64)
    isempty(values) && _blocked("bootstrap", "bootstrap has no values")
    ordered = sort(Float64.(values))
    position = 1 + (length(ordered) - 1) * probability
    lower = floor(Int, position); upper = ceil(Int, position)
    lower == upper ? ordered[lower] : ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])
end

function _public_material(root::String, execution_relative::String, index, controls)
    paths = Set(record.path for record in index.records)
    registry = Dict{Tuple,PublicAxisMaterial}()
    for axis in controls
        base = _artifact_relative(execution_relative, PUBLIC_MATERIAL_SUFFIX)
        todo3 = "$base/todo3/$(_axis_dir(axis))"
        todo10 = "$base/todo10/$(_axis_dir(axis))"
        todo11 = "$base/todo11/$(_axis_dir(axis))"
        declared = vcat(["$todo3/$name" for name in PUBLIC_TODO3_NAMES],
                        ["$todo10/$name" for name in PUBLIC_TODO10_NAMES],
                        ["$todo11/$name" for name in PUBLIC_TODO11_NAMES])
        all(path in paths for path in declared) ||
            _blocked("public_material", "public material is incomplete for one control axis")
        for path in declared
            _snapshot(root, path, "public material $path")
        end
        registry[_axis_key(axis)] = PublicAxisMaterial(
            _axis_key(axis), "$todo3/forward.tsv", "$todo3/backward.tsv",
            "$todo3/forward-receipt.toml", "$todo3/backward-receipt.toml",
            "$todo10/edge_observations.tsv", "$todo10/node_segments.tsv", "$todo10/receipt.toml",
            "$todo11/fitted_edge_models.tsv", "$todo11/fitted_unary_nodes.tsv",
            "$todo11/fitted_edge_transforms.tsv", "$todo11/receipt.toml",
            "", "", "", "", "", nothing, "", "", "")
    end
    return registry
end

function _canonical_public_reason(reason)
    raw = String(reason)
    mapped = get(PUBLIC_TODO12_REASON_MAP, raw, nothing)
    mapped === nothing && _blocked("todo12_validation", "public Todo 12 reason is not in the closed vocabulary")
    return mapped
end

function _call_public_predecessors(root::String, registry, controls, feature_relative, feature_sha,
                                   candidate_relative, model_relative, universe_relative)
    # Todo 10's public validator is used only on the supplied producer material;
    # it never receives a caller-created observation vector.
    for axis in controls
        material = registry[_axis_key(axis)]
        # The producer paths are declared by the material registry.  A missing
        # producer input is a closed failure rather than a local reconstruction.
        forward = material.todo3_forward
        backward = material.todo3_backward
        forward_receipt = material.todo3_forward_receipt
        backward_receipt = material.todo3_backward_receipt
        try
            StructuredEdgeFeatures.validate_edge_bundle(
                root, dirname(material.todo10_edge);
                features=feature_relative,
                feature_sha256=feature_sha,
                candidate_config=candidate_relative,
                model_config=model_relative,
                universe_dir=universe_relative,
                forward_patches=forward,
                backward_patches=backward,
                forward_receipt=forward_receipt,
                backward_receipt=backward_receipt,
            )
        catch error
            _blocked("todo10_validation", sprint(showerror, error))
        end
        todo10_edge = _snapshot(root, material.todo10_edge, "Todo 10 edge observations")
        todo10_node = _snapshot(root, material.todo10_node, "Todo 10 node segments")
        todo10_receipt = _snapshot(root, material.todo10_receipt, "Todo 10 receipt")
        _snapshot(root, material.todo11_model, "Todo 11 fitted edge models")
        _snapshot(root, material.todo11_unary, "Todo 11 fitted unary nodes")
        _snapshot(root, material.todo11_transform, "Todo 11 fitted edge transforms")
        _snapshot(root, material.todo11_receipt, "Todo 11 receipt")
        try
            report = StructuredChainInference.infer_structured_chains(
                root;
                edge_observations=material.todo10_edge,
                node_segments=material.todo10_node,
                todo10_receipt=material.todo10_receipt,
                fitted_edge_models=material.todo11_model,
                fitted_unary_nodes=material.todo11_unary,
                fitted_edge_transforms=material.todo11_transform,
                todo11_receipt=material.todo11_receipt,
            )
            report.status in (:PASS, :FAIL, :BLOCKED, :SKIPPED) ||
                _blocked("todo12_validation", "Todo 12 status is not a public status")
            public_reason = _canonical_public_reason(report.reason)
            report.status == :FAIL &&
                any(reference.status != :FAIL || !isempty(reference.blocks) for reference in report.references) &&
                _blocked("todo12_validation", "Todo 12 FAIL evidence is inconsistent")
            model_binding = _snapshot(root, material.todo11_model, "Todo 11 fitted edge models")
            unary_binding = _snapshot(root, material.todo11_unary, "Todo 11 fitted unary nodes")
            transform_binding = _snapshot(root, material.todo11_transform, "Todo 11 fitted edge transforms")
            receipt_binding = _snapshot(root, material.todo11_receipt, "Todo 11 receipt")
            public_status = report.status == :SKIPPED ? "BLOCKED" : String(report.status)
            material = PublicAxisMaterial(
                material.axis_key, material.todo3_forward, material.todo3_backward,
                material.todo3_forward_receipt, material.todo3_backward_receipt,
                material.todo10_edge, material.todo10_node, material.todo10_receipt,
                material.todo11_model, material.todo11_unary, material.todo11_transform,
                material.todo11_receipt, todo10_edge.sha256, model_binding.sha256,
                unary_binding.sha256, transform_binding.sha256, receipt_binding.sha256,
                report, report.result_sha256, public_status, public_reason)
            registry[_axis_key(axis)] = material
        catch error
            _blocked("todo12_validation", sprint(showerror, error))
        end
    end
    return registry
end

function _bind_public_hashes(execution, registry, controls)
    for axis in controls
        material = registry[_axis_key(axis)]
        axis_key = _axis_key(axis)
        for item in values(execution.nodes)
            item.axis === nothing && continue
            _axis_key(item.axis) == axis_key || continue
            row = item.row
            (row[21] == material.todo11_model_sha256 ||
             any(reference.fit_sha256 == row[21] for reference in something(material.report).references)) ||
                _blocked("hash_binding", "node fit hash is not the indexed Todo 11 or reference fit")
            row[22] == material.todo11_receipt_sha256 || _blocked("hash_binding", "node decision hash is not the indexed Todo 11 receipt")
            row[23] == material.todo10_sha256 || _blocked("hash_binding", "node edge hash is not the indexed Todo 10 table")
            row[24] == material.result_sha256 || _blocked("hash_binding", "node graph hash is not the public Todo 12 result")
        end
        for item in values(execution.edges)
            _axis_key(item.axis) == axis_key || continue
            row = item.row
            row[24] == material.todo10_sha256 || _blocked("hash_binding", "edge hash is not the indexed Todo 10 table")
            row[25] == material.todo11_transform_sha256 || _blocked("hash_binding", "edge transform hash is not the indexed Todo 11 transform")
        end
    end
    return nothing
end

function _compare_public_edge_truth(root::String, execution, registry, controls, scope)
    for axis in controls
        material = registry[_axis_key(axis)]
        public_rows = _parse_tsv(_snapshot_bytes(root, material.todo10_edge, "Todo 10 edge observations"),
                                 StructuredEdgeFeatures.EDGE_HEADER, "Todo 10 edge observations")
        public_map = Dict{Tuple{String,Int,Int},Vector{String}}()
        for row in public_rows
            left = _canonical_int(row[2], "public edge left"; positive=true)
            right = _canonical_int(row[3], "public edge right"; positive=true)
            key = (row[1], left, right)
            haskey(public_map, key) && _blocked("todo10_projection", "public edge key repeats")
            public_map[key] = row
        end
        observed = Set{Tuple{String,Int,Int}}()
        for item in values(execution.edges)
            _axis_key(item.axis) == _axis_key(axis) || continue
            row = item.row
            row[2] == scope || _blocked("todo10_projection", "execution edge scope differs")
            row[12] == (axis.control_class == "reversal" ? "backward" : "forward") ||
                _blocked("execution_direction", "execution edge direction differs")
            left = _canonical_int(row[13], "execution edge left"; positive=true)
            right = _canonical_int(row[14], "execution edge right"; positive=true)
            key = (row[5], left, right)
            haskey(public_map, key) || _blocked("todo10_projection", "execution edge is absent from public Todo 10")
            public_row = public_map[key]
            observed_key = (public_row[1], _canonical_int(public_row[2], "public edge left"; positive=true),
                            _canonical_int(public_row[3], "public edge right"; positive=true))
            push!(observed, observed_key)
            execution_fields = [row[5], row[13], row[14], row[15], row[16], row[17], row[18], row[19]]
            public_fields = public_row[1:8]
            execution_fields == public_fields ||
                _blocked("todo10_projection", "execution edge endpoints or topology differ from public Todo 10")
            row[20:23] == public_row[9:12] ||
                _blocked("todo10_projection", "execution edge status, split, or correlations differ from public Todo 10")
        end
        Set(keys(public_map)) == observed || _blocked("todo10_projection", "public Todo 10 edge key coverage differs")
    end
    return nothing
end

function _snapshot_bytes(root::String, relative::String, role::String)
    return _snapshot(root, relative, role).bytes
end

function _public_patch_payload_indices(header::Vector{String}, channel::String)
    prefix = channel == "forward" ? "res_p" : "bwd_res_p"
    expected = [@sprintf("%s%03d", prefix, index) for index in 1:289]
    observed = [name for name in header if startswith(name, prefix)]
    observed == expected || _blocked("todo10_projection", "$channel residual payload header differs")
    indices = [findfirst(==(name), header)::Int for name in expected]
    indices == collect(first(indices):last(indices)) ||
        _blocked("todo10_projection", "$channel residual payload block is not contiguous")
    return indices
end

function _compare_public_patch_truth(root::String, execution, registry, controls, frozen, t_values,
                                    feature_metadata)
    baseline = Dict{Tuple{String,Int,String},String}()
    for item in values(execution.transform.map)
        item.axis.control_class == "baseline" || continue
        key = (item.file, item.lobe, item.channel)
        old = get(baseline, key, item.baseline_row_sha256)
        old == item.baseline_row_sha256 || _blocked("transform_projection", "baseline transform identity varies")
        baseline[key] = old
    end
    expected_baseline = Set((key.file, key.lobe, channel)
                            for key in frozen.keys for channel in ("forward", "backward"))
    Set(keys(baseline)) == expected_baseline ||
        _blocked("transform_projection", "public baseline binding coverage differs")
    ordered = _observed_order(frozen, t_values)
    for axis in controls
        material = registry[_axis_key(axis)]
        for (channel, path, header) in (
            ("forward", material.todo3_forward, StructuredEdgeFeatures.FORWARD_PATCH_HEADER),
            ("backward", material.todo3_backward, StructuredEdgeFeatures.BACKWARD_PATCH_HEADER))
            rows = _parse_tsv(_snapshot_bytes(root, path, "Todo 3 $channel patch table"), header,
                              "Todo 3 $channel patch table")
            payload_indices = _public_patch_payload_indices(header, channel)
            seen = Set{Tuple{String,Int}}()
            for row in rows
                file = row[1]; lobe = _canonical_int(row[2], "Todo 3 patch lobe"; positive=true)
                key = (file, lobe)
                key in seen && _blocked("duplicate_key", "Todo 3 patch key repeats")
                push!(seen, key)
                values = Union{Missing,Float64}[]
                for field in row[payload_indices]
                    public_value = _public_float(field, "Todo 3 patch value"; allow_na=true)
                    push!(values, public_value === nothing ? missing : public_value)
                end
                haskey(ordered, file) || _blocked("transform_projection", "Todo 3 patch file is not frozen")
                lobe in ordered[file] || _blocked("transform_projection", "Todo 3 patch lobe is not frozen")
                feature = get(feature_metadata, key, nothing)
                feature === nothing && _blocked("transform_projection", "Todo 3 patch metadata key is not frozen")
                (row[1] == feature.file && row[2] == feature.lobe && row[3] == feature.t_nm &&
                 row[4] == feature.u_nm && row[5] == feature.amplitude) ||
                    _blocked("transform_projection", "Todo 3 patch metadata differs from the feature row")
                source_lobe = axis.control_class == "reversal" ? _involution(ordered[file])[lobe] : lobe
                expected_baseline = get(baseline, (file, source_lobe, channel), nothing)
                expected_baseline === nothing && _blocked("transform_projection", "Todo 3 patch key is not frozen")
                date = StructuredUniverse.parse_scan_date(file)
                axis_key = material.axis_key
                transform_key = (date, date, file, axis_key[1], axis_key[2], axis_key[3],
                                 axis_key[4], axis_key[5], lobe, channel)
                expected_transform = get(execution.transform.map, transform_key, nothing)
                expected_transform === nothing && _blocked("transform_projection", "Todo 3 patch key is not declared")
                expected_transform.baseline_row_sha256 == expected_baseline ||
                    _blocked("transform_projection", "Todo 3 baseline is not the global baseline material")
                _patch_hash(channel, values) == expected_transform.perturbed_row_sha256 ||
                    _blocked("transform_projection", "Todo 3 patch bytes differ from transform material")
            end
            expected_keys = Set((key.file, key.lobe) for key in frozen.keys)
            seen == expected_keys || _blocked("transform_projection", "Todo 3 patch key coverage differs")
        end
    end
    return nothing
end

_public_reference_identity(reference) =
    (reference.fit_sha256, reference.model, reference.fit_role, reference.scope,
     reference.outer_date, reference.target_date)

function _public_failure_nodes(execution, registry, axis)
    material = registry[_axis_key(axis)]
    rows = [item.row for item in values(execution.nodes) if _axis_key(item.axis) == _axis_key(axis)]
    expected_keys = Set((row[4], row[5], _canonical_int(row[6], "failed execution lobe"; positive=true), row[7])
                        for row in rows)
    length(rows) == length(expected_keys) && !isempty(rows) ||
        _blocked("todo12_projection", "public failure node coverage is incomplete")
    for row in rows
        row[16] == "failed" && row[15] == "?" && row[14] == "NA" && row[19] == "NA" &&
            row[20] == "false" && row[17] in ALLOWED_NODE_REASONS && row[17] != "none" ||
            _blocked("todo12_projection", "public failure node is not the exact failed tuple")
        views = _canonical_int(row[18], "failed execution views_used")
        0 <= views <= 4 || _blocked("todo12_projection", "failed execution views_used is outside the node schema")
        all(_hash_re(row[index]) for index in 21:25) ||
            _blocked("todo12_projection", "failed execution node identity is incomplete")
        (row[21] == material.todo11_model_sha256 ||
         any(reference.fit_sha256 == row[21] for reference in something(material.report).references)) ||
            _blocked("todo12_projection", "failed execution fit identity has success residue")
        row[22] == material.todo11_receipt_sha256 && row[23] == material.todo10_sha256 &&
            row[24] == material.result_sha256 ||
            _blocked("todo12_projection", "failed execution decision or graph identity differs")
    end
    return rows
end

function _public_report_integrity(report)
    identities = Set{Tuple}()
    for reference in report.references
        identity = _public_reference_identity(reference)
        identity in identities && _blocked("duplicate_key", "public Todo 12 reference identity repeats")
        push!(identities, identity)
        block_ids = Set{String}()
        for block in reference.blocks
            block.segment_id in block_ids && _blocked("duplicate_key", "public Todo 12 block identity repeats")
            push!(block_ids, block.segment_id)
            node_keys = Set((node.file, node.date, node.lobe, node.t_nm) for node in block.nodes)
            length(node_keys) == length(block.nodes) || _blocked("duplicate_key", "public Todo 12 node identity repeats")
        end
    end
    return identities
end

function _public_node_truth(execution, registry, controls, frozen, t_values)
    ordered = _observed_order(frozen, t_values)
    for axis in controls
        material = registry[_axis_key(axis)]
        report = material.report
        report === nothing && _blocked("todo12_binding", "public Todo 12 report is absent")
        mapped_public_reason = _canonical_public_reason(report.reason)
        identities = _public_report_integrity(report)
        if material.status in ("FAIL", "BLOCKED")
            expected_status = material.status == "FAIL" ? (:FAIL,) : (:BLOCKED, :SKIPPED)
            report.status in expected_status && material.reason == mapped_public_reason ||
                _blocked("todo12_validation", "public Todo 12 status/reason differs")
            _public_failure_nodes(execution, registry, axis)
            isempty(report.blocks) &&
                all(reference.status != :PASS && isempty(reference.blocks) for reference in report.references) ||
                _blocked("todo12_validation", "public Todo 12 failure report contains success residue")
            material.status == "FAIL" &&
                all(reference.status == :FAIL &&
                    (reference.reason == Symbol(mapped_public_reason) || reference.reason == :todo11_fail)
                    for reference in report.references) || material.status == "BLOCKED" ||
                _blocked("todo12_validation", "public Todo 12 FAIL evidence is not complete")
            continue
        end
        report.status == :PASS || _blocked("todo12_validation", "public Todo 12 status was not retained")
        selected_reference_ids = Set{Tuple}()
        expected_selected_nodes = Set{Tuple}()
        observed_selected_nodes = Set{Tuple}()
        fallback_to_c1 = false
        for item in values(execution.nodes)
            _axis_key(item.axis) == _axis_key(axis) || continue
            row = item.row
            file, lobe = row[5], _canonical_int(row[6], "execution node lobe"; positive=true)
            date = row[4]
            haskey(ordered, file) || _blocked("todo12_projection", "execution node file is not frozen")
            public_t = t_values[(file, lobe)]
            selected_model = row[19]
            selected_model in ("C1", "C2") ||
                _blocked("todo12_binding", "successful execution selected unary model is not C1 or C2")
            row[7] == "C1" && selected_model == "C1" || row[7] == "meta" ||
                _blocked("todo12_projection", "C1 execution row does not select C1")
            row[7] == "meta" && selected_model == "C1" && (fallback_to_c1 = true)
            model_hints = Set([selected_model])
            candidates = NamedTuple[]
            for reference in report.references
                reference.model ∈ model_hints || continue
                reference.fit_role == "partition" || continue
                reference.target_date == date || continue
                for block in reference.blocks, node in block.nodes
                    node.file == file && node.date == date && node.lobe == lobe && node.t_nm == public_t || continue
                    push!(candidates, (reference=reference, block=block, node=node,
                                       graph_enabled=!isempty(block.factors),
                                       identity=_public_reference_identity(reference)))
                end
            end
            exact_fit_candidates = [candidate for candidate in candidates if row[21] == candidate.reference.fit_sha256]
            candidates = isempty(exact_fit_candidates) ?
                         [candidate for candidate in candidates if row[21] == material.todo11_model_sha256] :
                         exact_fit_candidates
            length(candidates) == 1 || _blocked("todo12_binding", "execution node does not resolve to exactly one public reference")
            resolved = only(candidates)
            resolved.identity in identities || _blocked("todo12_binding", "selected public reference is not declared")
            resolved.reference.status in (:PASS, :SKIPPED) ||
                _blocked("todo12_projection", "execution node resolved to a noncomplete public reference")
            push!(selected_reference_ids, resolved.identity)
            push!(expected_selected_nodes, (resolved.identity, file, date, lobe, public_t))
            public_node = resolved.node
            push!(observed_selected_nodes, (resolved.identity, public_node.file, public_node.date,
                                            public_node.lobe, public_node.t_nm))
            public_node.marginals[2] == something(_canonical_float(row[14], "execution posterior")) ||
                _blocked("todo12_projection", "execution posterior differs from public Todo 12")
            row[15] == string(public_node.label) || _blocked("todo12_projection", "execution hard output differs from public Todo 12")
            row[16] == "ok" && row[17] == "none" ||
                _blocked("todo12_projection", "execution node status differs from complete public Todo 12")
            row[20] == _bool(resolved.graph_enabled) || _blocked("todo12_projection", "graph identity differs from public Todo 12")
        end
        observed_selected_nodes == expected_selected_nodes ||
            _blocked("todo12_projection", "selected public Todo 12 node coverage differs")
        public_selected_nodes = Set{Tuple}()
        for reference in report.references
            reference.fit_role == "partition" || continue
            identity = _public_reference_identity(reference)
            if identity in selected_reference_ids
                for block in reference.blocks, node in block.nodes
                    push!(public_selected_nodes, (identity, node.file, node.date, node.lobe, node.t_nm))
                end
            elseif !(fallback_to_c1 && reference.model == "C2")
                _blocked("todo12_projection", "public Todo 12 report contains an orphan reference")
            end
        end
        public_selected_nodes == expected_selected_nodes ||
            _blocked("todo12_projection", "selected public Todo 12 reference coverage differs")
    end
    return nothing
end

function _public_reversal_check(registry, controls, frozen, t_values, execution)
    baseline_axis = only(filter(axis -> axis.control_class == "baseline", controls))
    reversal_axis = only(filter(axis -> axis.control_class == "reversal", controls))
    base = get(registry, _axis_key(baseline_axis), nothing)
    reversed = get(registry, _axis_key(reversal_axis), nothing)
    (base === nothing || reversed === nothing) && _blocked("reversal_inference", "public Todo 12 reversal material is absent")
    (base.status in ("FAIL", "BLOCKED") || reversed.status in ("FAIL", "BLOCKED")) && return Inf
    base_report = something(base.report); reverse_report = something(reversed.report)
    # The producer fit/provenance hashes are axis-local.  Reversal identity is
    # therefore physical and structural, not a byte identity between fits.
    reference_key(reference) = (reference.model, reference.fit_role, reference.scope,
                                reference.outer_date, reference.target_date)
    base_refs = Dict{Tuple,Any}(); reverse_refs = Dict{Tuple,Any}()
    for reference in base_report.references
        key = reference_key(reference)
        haskey(base_refs, key) && return Inf
        base_refs[key] = reference
    end
    for reference in reverse_report.references
        key = reference_key(reference)
        haskey(reverse_refs, key) && return Inf
        reverse_refs[key] = reference
    end
    Set(keys(base_refs)) == Set(keys(reverse_refs)) || return Inf
    maximum_delta = 0.0
    observed = _observed_order(frozen, t_values)
    for key in sort!(collect(keys(base_refs)))
        base_reference = base_refs[key]; reverse_reference = reverse_refs[key]
        base_reference.status == reverse_reference.status && base_reference.reason == reverse_reference.reason || return Inf
        abs(base_reference.log_evidence - reverse_reference.log_evidence) <= 1e-12 || return Inf
        abs(base_reference.viterbi_score - reverse_reference.viterbi_score) <= 1e-12 || return Inf
        paired_blocks = Tuple{Any,Any}[]
        unused_reverse = collect(reverse_reference.blocks)
        for natural in base_reference.blocks
            expected_keys = Set((node.file, node.date, _involution(observed[node.file])[node.lobe],
                                 t_values[(node.file, _involution(observed[node.file])[node.lobe])])
                                for node in natural.nodes)
            matches = [block for block in unused_reverse
                       if Set((node.file, node.date, node.lobe, node.t_nm) for node in block.nodes) == expected_keys]
            length(matches) == 1 || return Inf
            matched = only(matches)
            match_index = findfirst(block -> block === matched, unused_reverse)
            match_index === nothing && return Inf
            deleteat!(unused_reverse, match_index)
            push!(paired_blocks, (natural, matched))
        end
        isempty(unused_reverse) || return Inf
        for (natural, reversed_block) in paired_blocks
            abs(natural.log_evidence - reversed_block.log_evidence) <= 1e-12 || return Inf
            abs(natural.viterbi_score - reversed_block.viterbi_score) <= 1e-12 || return Inf
            base_nodes = Dict{Tuple,Any}(); reverse_nodes = Dict{Tuple,Any}()
            for node in natural.nodes
                node_key = (node.file, node.date, node.lobe, node.t_nm)
                haskey(base_nodes, node_key) && return Inf
                base_nodes[node_key] = node
            end
            for node in reversed_block.nodes
                node_key = (node.file, node.date, node.lobe, node.t_nm)
                haskey(reverse_nodes, node_key) && return Inf
                reverse_nodes[node_key] = node
            end
            expected_reverse_keys = Set{Tuple}()
            for (file, date, lobe, _) in keys(base_nodes)
                inv = _involution(observed[file])
                push!(expected_reverse_keys, (file, date, inv[lobe], t_values[(file, inv[lobe])]))
            end
            Set(keys(reverse_nodes)) == expected_reverse_keys || return Inf
            for (node_key, left) in base_nodes
                file, date, lobe, _ = node_key
                inv_lobe = _involution(observed[file])[lobe]
                right = reverse_nodes[(file, date, inv_lobe, t_values[(file, inv_lobe)])]
                for index in 1:2
                    maximum_delta = max(maximum_delta, abs(left.marginals[index] - right.marginals[index]))
                    maximum_delta = max(maximum_delta, abs(left.log_marginals[index] - right.log_marginals[index]))
                end
                left.label == right.label || return Inf
            end
            Tuple(reverse([node.label for node in natural.nodes])) == reversed_block.labels || return Inf
            base_factors = Dict{Tuple,Any}(); reverse_factors = Dict{Tuple,Any}()
            for factor in natural.factors
                factor_key = (factor.file, factor.date, factor.left_lobe, factor.right_lobe,
                              factor.left_t_nm, factor.right_t_nm)
                haskey(base_factors, factor_key) && return Inf
                base_factors[factor_key] = factor
            end
            for factor in reversed_block.factors
                factor_key = (factor.file, factor.date, factor.left_lobe, factor.right_lobe,
                              factor.left_t_nm, factor.right_t_nm)
                haskey(reverse_factors, factor_key) && return Inf
                reverse_factors[factor_key] = factor
            end
            expected_factor_keys = Set{Tuple}()
            for (file, date, left, right, _, _) in keys(base_factors)
                inv = _involution(observed[file])
                push!(expected_factor_keys, (file, date, inv[right], inv[left],
                                              t_values[(file, inv[right])], t_values[(file, inv[left])]))
            end
            Set(keys(reverse_factors)) == expected_factor_keys || return Inf
            factor_segment_map = Dict{String,String}()
            factor_segment_consistent = true
            for (factor_key, left_factor) in base_factors
                file, date, left, right, _, _ = factor_key
                inv = _involution(observed[file])
                right_factor = reverse_factors[(file, date, inv[right], inv[left],
                                                t_values[(file, inv[right])], t_values[(file, inv[left])])]
                factor_segment_consistent &= _bind_segment_mapping!(factor_segment_map,
                                                                      left_factor.segment_id,
                                                                      right_factor.segment_id)
                left_factor.model == right_factor.model || return Inf
                left_factor.left_t_nm == t_values[(file, left)] &&
                    left_factor.right_t_nm == t_values[(file, right)] || return Inf
                right_factor.left_t_nm == t_values[(file, inv[right])] &&
                    right_factor.right_t_nm == t_values[(file, inv[left])] || return Inf
                for index in 1:4
                    reverse_index = index == 1 ? 1 : index == 2 ? 3 : index == 3 ? 2 : 4
                    maximum_delta = max(maximum_delta, abs(left_factor.factor[index] -
                                                            right_factor.factor[reverse_index]))
                end
            end
            factor_segment_consistent || return Inf
        end
    end
    abs(base_report.log_evidence - reverse_report.log_evidence) <= 1e-12 || return Inf
    abs(base_report.viterbi_score - reverse_report.viterbi_score) <= 1e-12 || return Inf
    return maximum_delta
end

function _bind_segment_mapping!(mapping::Dict{String,String}, source::String, destination::String)
    old = get(mapping, source, nothing)
    if old === nothing
        mapping[source] = destination
        return true
    end
    return old == destination
end

function _topology_number_equal(left::String, right::String, context::String)
    left == "NA" && right == "NA" && return true
    (left == "NA") == (right == "NA") || return false
    a = something(_canonical_float(left, "$context baseline"))
    b = something(_canonical_float(right, "$context perturbed"))
    return isapprox(a, b; atol=1e-12, rtol=1e-12)
end

function _topology_equivalent(base, pert, left::Int, right::Int, output_left::Int,
                              output_right::Int, file::String, t_values, segment_map,
                              control_class::String)
    base[20] == pert[20] && base[21] == pert[21] || return false
    if control_class == "reversal"
        get(segment_map, base[18], nothing) == pert[19] || return false
        get(segment_map, base[19], nothing) == pert[18] || return false
    else
        get(segment_map, base[18], nothing) == pert[18] || return false
        get(segment_map, base[19], nothing) == pert[19] || return false
    end
    _topology_number_equal(base[22], pert[22], "topology forward correlation") || return false
    _topology_number_equal(base[23], pert[23], "topology backward correlation") || return false
    base_left_t = something(_canonical_float(base[15], "topology baseline left t"))
    base_right_t = something(_canonical_float(base[16], "topology baseline right t"))
    pert_left_t = something(_canonical_float(pert[15], "topology perturbed left t"))
    pert_right_t = something(_canonical_float(pert[16], "topology perturbed right t"))
    base_left_t == t_values[(file, left)] && base_right_t == t_values[(file, right)] || return false
    pert_left_t == t_values[(file, output_left)] && pert_right_t == t_values[(file, output_right)] || return false
    return true
end

function _metrics(execution, frozen, controls, t_values, registry)
    nodes, edges = execution.nodes, execution.edges
    lobe_rows = Vector{Vector{String}}(); scan_rows = Vector{Vector{String}}(); topology_rows = Vector{Vector{String}}()
    summaries = Dict{Tuple,NamedTuple}()
    scans_by_date = Dict{String,Vector{String}}()
    for scan in frozen.scans
        push!(get!(scans_by_date, scan.date, String[]), scan.file)
    end
    for axis in controls
        axis_key = _axis_key(axis)
        for model in ("C1", "meta")
            by_date = Dict{String,Vector{NTuple{6,Float64}}}()
            numerical_failures = 0
            for scan in frozen.scans
                keys_for_scan = [(scan.date, scan.file, lobe, model, axis_key) for lobe in
                                 _observed_order(frozen, t_values)[scan.file]]
                base_axis = only(filter(item -> item.control_class == "baseline", controls))
                order = _observed_order(frozen, t_values)[scan.file]
                involution = _involution(order)
                baseline_keys = [(scan.date, scan.file,
                                  axis.control_class == "reversal" ? involution[lobe] : lobe,
                                  model, _axis_key(base_axis)) for lobe in order]
                agrees = 0; base_class = 0; pert_class = 0; abstentions = 0; failures = 0
                posterior_a = Float64[]; posterior_b = Float64[]
                for (key, base_key) in zip(keys_for_scan, baseline_keys)
                    p = nodes[key]; b = nodes[base_key]
                    p_classified = p.status == "ok" && p.hard != "?"
                    b_classified = b.status == "ok" && b.hard != "?"
                    p_classified && (pert_class += 1); b_classified && (base_class += 1)
                    !p_classified && (abstentions += 1)
                    p.status == "failed" && (failures += 1)
                    p.status == "failed" && (numerical_failures += 1)
                    agreement = p_classified && b_classified && p.hard == b.hard
                    agreement && (agrees += 1)
                    p.posterior !== nothing && b.posterior !== nothing &&
                        (push!(posterior_a, something(b.posterior)); push!(posterior_b, something(p.posterior)))
                    push!(lobe_rows, ["1", execution.scope, scan.date, scan.date, scan.file, string(key[3]), model,
                                      axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed,
                                      axis.target, b.posterior === nothing ? "NA" : _fmt(something(b.posterior)),
                                      p.posterior === nothing ? "NA" : _fmt(something(p.posterior)), b.hard, p.hard,
                                      b.status, p.status, b.reason, p.reason, _bool(b_classified), _bool(p_classified),
                                      b.posterior === nothing || p.posterior === nothing ? "NA" : _fmt(abs(something(p.posterior) - something(b.posterior))),
                                      _bool(agreement), _sha_text(join(b.row, '\t') * "\n"), _sha_text(join(p.row, '\t') * "\n")])
                end
                corr = _correlation(posterior_a, posterior_b)
                loss = 1.0 - agrees / length(keys_for_scan)
                push!(scan_rows, ["1", execution.scope, scan.date, scan.date, scan.file, model, axis.control_class,
                                  axis.perturbation, axis.view, axis.seed_family, axis.seed, axis.target,
                                  string(length(keys_for_scan)), string(agrees), _fmt(loss), string(base_class),
                                  string(pert_class), _fmt(base_class / length(keys_for_scan)),
                                  _fmt(pert_class / length(keys_for_scan)), string(length(posterior_a)),
                                  corr === nothing ? "NA" : _fmt(corr), corr === nothing ? "undefined" : "defined",
                                  string(abstentions), string(failures)])
                push!(get!(by_date, scan.date, NTuple{6,Float64}[]),
                      (loss, base_class / length(keys_for_scan), pert_class / length(keys_for_scan),
                       agrees / length(keys_for_scan), corr === nothing ? 0.0 : 1.0, failures))
            end
            date_values = [by_date[date] for date in sort!(collect(keys(by_date)))]
            summaries[(axis_key, model)] = (
                losses=[mean(item[1] for item in values) for values in date_values],
                base_coverage=[mean(item[2] for item in values) for values in date_values],
                pert_coverage=[mean(item[3] for item in values) for values in date_values],
                agrees=[mean(item[4] for item in values) for values in date_values],
                corr_count=sum(round(Int, item[5]) for values in date_values for item in values),
                nodes=sum(length(_observed_order(frozen, t_values)[scan.file]) for scan in frozen.scans),
                numerical_failures=numerical_failures,
                abstentions=sum(item[6] for values in date_values for item in values))
        end
        for scan in frozen.scans
            order = _observed_order(frozen, t_values)[scan.file]
            involution = _involution(order)
            base_axis = only(filter(item -> item.control_class == "baseline", controls))
            edge_pairs = Tuple{Int,Int,Int,Int}[]
            segment_map = Dict{String,String}()
            segment_consistent = true
            for index in 1:max(0, length(order) - 1)
                left, right = order[index:index + 1]
                output_left, output_right = axis.control_class == "reversal" ? (involution[right], involution[left]) : (left, right)
                key = (scan.date, scan.file, output_left, output_right, axis_key)
                base = edges[(scan.date, scan.file, left, right, _axis_key(base_axis))].row
                pert = edges[key].row
                push!(edge_pairs, (left, right, output_left, output_right))
                if axis.control_class == "reversal"
                    segment_consistent &= _bind_segment_mapping!(segment_map, base[18], pert[19])
                    segment_consistent &= _bind_segment_mapping!(segment_map, base[19], pert[18])
                else
                    segment_consistent &= _bind_segment_mapping!(segment_map, base[18], pert[18])
                    segment_consistent &= _bind_segment_mapping!(segment_map, base[19], pert[19])
                end
            end
            for (left, right, output_left, output_right) in edge_pairs
                base = edges[(scan.date, scan.file, left, right, _axis_key(base_axis))].row
                pert = edges[(scan.date, scan.file, output_left, output_right, axis_key)].row
                equivalent = segment_consistent &&
                    _topology_equivalent(base, pert, left, right, output_left, output_right,
                                         scan.file, t_values, segment_map, axis.control_class)
                push!(topology_rows, ["1", execution.scope, scan.date, scan.date, scan.file, axis.control_class,
                                      axis.perturbation, axis.view, axis.seed_family, axis.seed, axis.target,
                                      string(left), string(right), base[20], pert[20], base[21], pert[21], base[22],
                                      pert[22], base[23], pert[23], base[18], base[19], pert[18], pert[19],
                                      _bool(equivalent)])
        end
    end
    end
    sort!(lobe_rows; by=Tuple); sort!(scan_rows; by=Tuple); sort!(topology_rows; by=Tuple)
    return (lobe=_tsv(LOBE_HEADER, lobe_rows), scan=_tsv(SCAN_HEADER, scan_rows), topology=_tsv(TOPOLOGY_HEADER, topology_rows), summaries=summaries, lobe_rows=lobe_rows, scan_rows=scan_rows, topology_rows=topology_rows)
end

function _bootstrap(metrics, frozen, controls, execution_scope)
    rows = Vector{Vector{String}}(); upper = Dict{Tuple,Float64}()
    for axis in controls
        axis.control_class in ("image_perturbation", "unary_view_drop", "patch_channel_drop", "restart") || continue
        key = _axis_key(axis)
        c1_rows = filter(row -> Tuple(row[7:12]) == key && row[6] == "C1", metrics.scan_rows)
        meta_rows = filter(row -> Tuple(row[7:12]) == key && row[6] == "meta", metrics.scan_rows)
        by_date = Dict{String,Vector{Tuple{String,Float64,Float64}}}()
        for row in c1_rows
            meta = only(filter(candidate -> candidate[5] == row[5], meta_rows))
            push!(get!(by_date, row[3], Tuple{String,Float64,Float64}[]), (row[5], parse(Float64, row[15]), parse(Float64, meta[15])))
        end
        date_names = sort(collect(keys(by_date)))
        for seed in 0:499
            rng = MersenneTwister(seed); c1_dates=Float64[]; meta_dates=Float64[]; scan_count=0
            for date in date_names
                scans = by_date[date]; samples = [scans[rand(rng, 1:length(scans))] for _ in scans]
                push!(c1_dates, mean(item[2] for item in samples)); push!(meta_dates, mean(item[3] for item in samples)); scan_count += length(samples)
            end
            c1 = mean(c1_dates); meta = mean(meta_dates); diff = meta - c1
            push!(rows, ["1", execution_scope, axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed, axis.target, string(seed), string(length(date_names)), string(scan_count), _fmt(meta), _fmt(c1), _fmt(diff)])
            upper[key] = get(upper, key, NaN)
        end
        differences = [parse(Float64, row[14]) for row in rows if Tuple(row[3:8]) == key]
        upper[key] = _type7(differences, 0.975)
    end
    sort!(rows; by=Tuple)
    return _tsv(BOOTSTRAP_HEADER, rows), upper
end

function _reason_bytes(metrics, execution_scope)
    counts = Dict{Tuple{String,String,String,String,String,String,String,String,String,String,String},Int}()
    for row in metrics.lobe_rows
        for (level, model, status, reason) in (("node", row[7], row[19], row[21]), ("node", row[7], row[20], row[22]))
            key = (level, model, row[8], row[9], row[10], row[11], row[12], row[13], status, reason, "")
            counts[key] = get(counts, key, 0) + 1
        end
    end
    rows = Vector{Vector{String}}()
    for (key, count) in sort(collect(counts); by=first)
        level, model, class, perturbation, view, family, seed, target, status, reason, _ = key
        push!(rows, ["1", execution_scope, level, model, class, perturbation, view, family, seed, target, status, reason, string(count)])
    end
    _tsv(REASON_HEADER, rows)
end

function _summary(metrics, bootstrap, upper, frozen, controls, execution_scope, computed_status, computed_reason,
                  reversal_delta, reversal_topology)
    rows = Vector{Vector{String}}()
    for axis in controls
        key = _axis_key(axis)
        summary = get(metrics.summaries, (key, "meta"), nothing)
        summary === nothing && continue
        mean_loss = mean(summary.losses)
        c1_summary = metrics.summaries[(key, "C1")]
        c1_loss = mean(c1_summary.losses)
        difference = mean_loss - c1_loss
        gate_role = axis.control_class == "restart" ? "report_only_restart" :
                    axis.control_class == "reversal" ? "reversal_equivalence" : "noninferiority"
        local_status = if computed_status == "PASS" || axis.control_class == "restart"
            "SKIPPED"
        elseif axis.control_class == "reversal"
            reversal_topology && reversal_delta <= 1e-12 ? "PASS" : "FAIL"
        elseif axis.control_class in ("image_perturbation", "unary_view_drop", "patch_channel_drop")
            get(upper, key, Inf) <= 0.0 ? "PASS" : "FAIL"
        else
            "PASS"
        end
        publication_status = axis.control_class == "restart" ? "SKIPPED" :
                             (computed_status == "PASS" ? "SKIPPED" : computed_status)
        push!(rows, ["1", execution_scope, axis.control_class, axis.perturbation, axis.view, axis.seed_family, axis.seed,
                     axis.target, gate_role, string(length(frozen.scans)), string(summary.nodes), _fmt(mean_loss),
                     _fmt(c1_loss), _fmt(difference), haskey(upper, key) ? _fmt(upper[key]) : "NA",
                     _fmt(mean(summary.pert_coverage)), _fmt(mean(c1_summary.pert_coverage)), _fmt(mean(summary.agrees)),
                     string(summary.corr_count), axis.control_class == "reversal" && isfinite(reversal_delta) ? _fmt(reversal_delta) : "NA",
                     axis.control_class == "reversal" ? _bool(reversal_topology && isfinite(reversal_delta) && reversal_delta <= 1e-12) : "NA",
                     local_status, publication_status, computed_reason])
    end
    isempty(rows) && push!(rows, ["1", execution_scope, "NA", "NA", "NA", "NA", "NA", "NA", "blocked", "0", "0", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "0", "NA", "NA", computed_status, computed_status, computed_reason])
    sort!(rows; by=Tuple)
    _tsv(SUMMARY_HEADER, rows)
end

function _assert_reason_correspondence(scientific_status::String, scientific_reason::String,
                                       summary_bytes::Vector{UInt8})
    expected_status = scientific_status == "PASS" ? "SKIPPED" : scientific_status
    expected_reason = scientific_status == "PASS" ? "synthetic_fixture_only" : scientific_reason
    rows = _parse_tsv(summary_bytes, SUMMARY_HEADER, "summary")
    !isempty(rows) && all((row[24] == expected_status || row[24] == "SKIPPED") &&
                          row[25] == expected_reason for row in rows) ||
        _blocked("integrity_failure", "structured stability summary/report/publication reason correspondence failed")
    scientific_reason == expected_reason ||
        _blocked("integrity_failure", "structured stability scientific reason correspondence failed")
    return (status=expected_status, reason=expected_reason)
end

function _empty_files(scope::String, status::String, reason::String)
    files = Dict{String,Vector{UInt8}}()
    files["expected_rows.tsv"] = _tsv(EXPECTED_HEADER, Vector{Vector{String}}())
    files["lobe_pairs.tsv"] = _tsv(LOBE_HEADER, Vector{Vector{String}}())
    files["scan_metrics.tsv"] = _tsv(SCAN_HEADER, Vector{Vector{String}}())
    files["topology_pairs.tsv"] = _tsv(TOPOLOGY_HEADER, Vector{Vector{String}}())
    files["bootstrap.tsv"] = _tsv(BOOTSTRAP_HEADER, Vector{Vector{String}}())
    files["reason_counts.tsv"] = _tsv(REASON_HEADER, Vector{Vector{String}}())
    files["summary.tsv"] = _tsv(SUMMARY_HEADER, [["1", scope, "NA", "NA", "NA", "NA", "NA", "NA", "blocked", "0", "0", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "0", "NA", "NA", status, status, reason]])
    return files
end

function _report_from_error(error; scope="unknown", bindings=Binding[], specs=ArtifactSpec[], files=nothing,
                            metrics=Dict{String,String}())
    failure = error isa AuditError ? error : AuditError(:BLOCKED, "internal_error", sprint(showerror, error))
    status = String(failure.status)
    reason = failure.reason in ALLOWED_REASON_CODES ? failure.reason : "internal_error"
    payloads = files === nothing ? _empty_files(scope, status, reason) : files
    return AuditReport(status, reason, status, reason, status, reason, false, false,
                       payloads, Dict(binding.relative => binding.sha256 for binding in bindings),
                       merge(Dict("scope" => scope), metrics), bindings, specs)
end

function _audit_core(root::String; features::String, feature_sha256::String, candidate_config::String,
                     model_config::String, universe_dir::String, execution_dir::String)
    canonical_root = _canonical_root(root)
    candidate_config == "config/unit_assignment_structured_candidate.toml" ||
        _blocked("path_substitution", "candidate config is not the frozen path")
    model_config == "config/unit_assignment_structured_model.toml" ||
        _blocked("path_substitution", "model config is not the frozen path")
    feature_relative = (_relative_path(canonical_root, features, "feature input"))[2]
    universe_relative = (_relative_path(canonical_root, universe_dir, "universe directory"; directory=true))[2]
    execution_relative = (_relative_path(canonical_root, execution_dir, "execution directory"; directory=true))[2]
    _hash_re(feature_sha256) || _blocked("feature_hash", "feature SHA-256 is malformed")
    authorities, source_specs = _source_authorities(canonical_root)
    candidate_binding = _config_binding(canonical_root)
    feature_binding = _snapshot(canonical_root, feature_relative, "feature input";
                                expected=feature_sha256)
    receipt_binding, receipt, scope = _parse_execution_receipt(canonical_root, execution_relative)
    _progress_update!(scope=scope)

    # A real producer is not registered in this lane.  Only the immutable source
    # closure and a small, exact execution envelope are inspected before blocking.
    if scope == "real_viper"
        _execution_shape(canonical_root, execution_relative)
        index_binding = _snapshot(canonical_root, "$execution_relative/artifact-index.tsv", "real execution artifact index")
        _hash(index_binding.bytes) == receipt["artifact_index_sha256"] || _blocked("hash", "real receipt index hash differs")
        index = StructuredEvidence.parse_artifact_index(index_binding.bytes)
        index.campaign_id == EXECUTION_CAMPAIGN && index.model_id == EXECUTION_MODEL ||
            _blocked("identity", "real execution registry identity differs")
        files = _empty_files(scope, "BLOCKED", "unregistered_real_producer")
        summary_text = String(copy(files["summary.tsv"]))
        files["summary.tsv"] = Vector{UInt8}(replace(summary_text, "\tunknown\t" => "\t$scope\t"))
        bindings = vcat(authorities, [candidate_binding, feature_binding, receipt_binding, index_binding])
        specs = vcat(ArtifactSpec[ArtifactSpec("config", candidate_binding.relative)], source_specs,
                     ArtifactSpec[ArtifactSpec("input", item.relative) for item in bindings if item.relative != candidate_binding.relative && item.relative ∉ [source.relative for source in authorities]])
        return AuditReport("BLOCKED", "unregistered_real_producer", "BLOCKED", "unregistered_real_producer",
                           "BLOCKED", "unregistered_real_producer", false, false, files,
                           Dict(item.relative => item.sha256 for item in bindings), Dict("scope" => scope), bindings, specs)
    end

    scope == "synthetic_fixture" || _blocked("scope", "unsupported execution scope")
    frozen = try
        StructuredUniverse.freeze_universe(
            canonical_root;
            features=feature_relative,
            feature_sha256=feature_sha256,
            candidate_config=candidate_config,
            model_config=model_config,
            source_paths=[joinpath(canonical_root, "test", "lib", "structured_assignment", "universe.jl"),
                         joinpath(canonical_root, "test", "build_structured_scan_universe.jl"),
                         joinpath(canonical_root, "test", "lib", "structured_assignment", "firewall.jl")],
        )
    catch error
        _blocked("universe_binding", sprint(showerror, error))
    end
    StructuredUniverse.validate_bundle(frozen, universe_relative)
    feature_binding, t_values, u_values, feature_metadata = _parse_feature_t(canonical_root, feature_relative)
    Set(keys(t_values)) == Set((key.file, key.lobe) for key in frozen.keys) ||
        _blocked("feature_key_mismatch", "feature t keys differ from frozen universe")
    Set(keys(u_values)) == Set(keys(t_values)) && Set(keys(feature_metadata)) == Set(keys(t_values)) ||
        _blocked("feature_key_mismatch", "feature physical metadata keys differ from frozen universe")
    receipt["universe_receipt_sha256"] ==
        _sha(_snapshot(canonical_root, "$universe_relative/receipt.toml", "universe receipt").bytes) ||
        _blocked("hash", "execution universe receipt hash differs")
    controls = _controls(frozen)
    length(controls) == 519 || _blocked("control_set", "stability control declaration count differs")
    execution_path = _execution_shape(canonical_root, execution_relative)
    execution = _execution_index(canonical_root, execution_relative, frozen, feature_relative, feature_sha256, receipt,
                                 controls, source_specs, universe_relative)
    run_bindings = _run_bindings(canonical_root, execution_relative, receipt)
    patch_binding, patches = _parse_patch_material(canonical_root, execution_relative, execution, receipt)
    noise_relative = _artifact_relative(execution_relative, NOISE_MATERIAL_SUFFIX)
    noise_binding = _safe_material_file(canonical_root, execution_relative, NOISE_MATERIAL_SUFFIX, "noise material")
    noise_rows = _parse_tsv(noise_binding.bytes, NOISE_HEADER, "noise material")
    transform = _material_rows(canonical_root, execution_relative, execution, frozen, controls, t_values,
                               patches, patch_binding, noise_binding, noise_rows, scope)
    universe_sha, universe_binding = _universe_sha(canonical_root, universe_relative)
    expected_bytes = _expected_rows_bytes(frozen, controls, t_values, scope, universe_sha)
    execution_rows = _validate_execution_rows(canonical_root, execution_relative, run_bindings, receipt, frozen,
                                               controls, t_values, execution.expected_binding, expected_bytes,
                                               scope, transform, universe_sha)
    registry = _public_material(canonical_root, execution_relative, execution, controls)
    registry = _call_public_predecessors(canonical_root, registry, controls, feature_relative, feature_sha256,
                                          candidate_config, model_config, universe_relative)
    _bind_public_hashes(execution_rows, registry, controls)
    _compare_public_patch_truth(canonical_root, execution_rows, registry, controls, frozen, t_values,
                                feature_metadata)
    _compare_public_edge_truth(canonical_root, execution_rows, registry, controls, scope)
    _public_node_truth(execution_rows, registry, controls, frozen, t_values)
    reversal_delta = _public_reversal_check(registry, controls, frozen, t_values, execution_rows)
    metrics = _metrics(execution_rows, frozen, controls, t_values, registry)
    bootstrap_bytes, upper = _bootstrap(metrics, frozen, controls, scope)
    gated_upper = [value for (axis_key, value) in upper if
                   only(filter(axis -> _axis_key(axis) == axis_key, controls)).control_class != "restart"]
    scientific_axis(axis_key) = only(filter(axis -> _axis_key(axis) == axis_key, controls)).control_class != "restart"
    public_failure = any(scientific_axis(axis_key) && material.status == "FAIL"
                         for (axis_key, material) in registry)
    public_failure_reasons = [material.reason for (axis_key, material) in registry
                              if scientific_axis(axis_key) && material.status == "FAIL"]
    public_blocked = [(axis_key, material.reason) for (axis_key, material) in registry
                      if scientific_axis(axis_key) && material.status == "BLOCKED"]
    numerical_failure = any(scientific_axis(key[1]) && summary.numerical_failures > 0
                            for (key, summary) in metrics.summaries)
    reversal_rows = [row for row in metrics.topology_rows if row[6] == "reversal"]
    reversal_topology = !isempty(reversal_rows) && all(row[26] == "true" for row in reversal_rows)
    scientific_status = !isempty(public_blocked) ? "BLOCKED" :
                        public_failure || numerical_failure || any(value > 0.0 for value in gated_upper) ||
                        !reversal_topology || reversal_delta > 1e-12 ? "FAIL" : "PASS"
    scientific_reason = !isempty(public_blocked) ? first(public_blocked)[2] :
                        !isempty(public_failure_reasons) ? first(public_failure_reasons) :
                        numerical_failure ? "numerical_failure" :
                        reversal_delta > 1e-12 ? "reversal_mismatch" :
                        !reversal_topology ? "topology" :
                        scientific_status == "PASS" ? "synthetic_fixture_only" : "stability_gate_failure"
    summary_bytes = _summary(metrics, bootstrap_bytes, upper, frozen, controls, scope, scientific_status,
                             scientific_reason, reversal_delta, reversal_topology)
    publication_identity = _assert_reason_correspondence(scientific_status, scientific_reason, summary_bytes)
    files = Dict(
        "expected_rows.tsv" => expected_bytes,
        "lobe_pairs.tsv" => metrics.lobe,
        "scan_metrics.tsv" => metrics.scan,
        "topology_pairs.tsv" => metrics.topology,
        "bootstrap.tsv" => bootstrap_bytes,
        "reason_counts.tsv" => _reason_bytes(metrics, scope),
        "summary.tsv" => summary_bytes,
    )
    bindings = Binding[]
    append!(bindings, authorities)
    append!(bindings, [candidate_binding, feature_binding, universe_binding, patch_binding, noise_binding,
                       execution.index_binding, execution.expected_binding, receipt_binding])
    append!(bindings, collect(values(run_bindings)))
    append!(bindings, [transform.binding])
    for record in execution.records
        record.path in (binding.relative for binding in bindings) && continue
        push!(bindings, _snapshot(canonical_root, record.path, "execution $(record.role) $(record.path)"; expected=record.sha256))
    end
    ordered_bindings = Dict{String,Binding}()
    for binding in bindings
        if haskey(ordered_bindings, binding.relative)
            ordered_bindings[binding.relative].sha256 == binding.sha256 ||
                _blocked("source_collision", "one canonical path has conflicting snapshots")
        else
            ordered_bindings[binding.relative] = binding
        end
    end
    bindings = [ordered_bindings[path] for path in sort(collect(keys(ordered_bindings)))]
    specs = ArtifactSpec[]
    push!(specs, ArtifactSpec("config", candidate_binding.relative))
    append!(specs, source_specs)
    append!(specs, [ArtifactSpec(record.role, record.path) for record in sort(execution.records; by=record -> (record.role, record.path))])
    unique_specs = Dict{String,ArtifactSpec}()
    for spec in specs
        if haskey(unique_specs, spec.path)
            unique_specs[spec.path].role == spec.role || _blocked("source_collision", "artifact spec role differs")
        else
            unique_specs[spec.path] = spec
        end
    end
    specs = [unique_specs[path] for path in sort(collect(keys(unique_specs)))]
    _progress_update!(scope=scope, files=files, bindings=bindings, specs=specs,
                      metrics=Dict("universe_sha256" => universe_sha,
                                   "reversal_max_abs_delta" => _fmt(reversal_delta)))
    for binding in bindings
        _verify(canonical_root, binding)
    end
    # The final evidence bundle owns the published expected rows and the other
    # six computed payloads; the execution expected rows remain an input.
    return AuditReport(publication_identity.status, publication_identity.reason, scientific_status, scientific_reason,
                       publication_identity.status, publication_identity.reason, false, false, files,
                       Dict(binding.relative => binding.sha256 for binding in bindings),
                       Dict("scope" => scope, "universe_sha256" => universe_sha, "reversal_max_abs_delta" => _fmt(reversal_delta)),
                       bindings, specs)
end

function _percent_decode_fixed(text::String, context::String)::String
    current = text
    for _ in 1:4
        occursin('%', current) || return current
        bytes = UInt8[]; index = firstindex(current)
        while index <= lastindex(current)
            character = current[index]
            if character == '%'
                index + 2 <= lastindex(current) || _blocked("firewall", "$context has malformed percent encoding")
                digits = current[index + 1:index + 2]
                all(character -> character in "0123456789abcdefABCDEF", digits) ||
                    _blocked("firewall", "$context has malformed percent encoding")
                push!(bytes, parse(UInt8, digits; base=16)); index += 3
            else
                encoded = codeunits(string(character))
                append!(bytes, encoded); index = nextind(current, index)
            end
        end
        decoded = try
            String(bytes)
        catch
            _blocked("firewall", "$context percent decoding is not UTF-8")
        end
        isascii(decoded) || _blocked("firewall", "$context percent decoding is not ASCII")
        decoded == current && return current
        current = decoded
    end
    occursin('%', current) && _blocked("firewall", "$context percent encoding does not reach a fixed point")
    return current
end

function _firewall_value(value::AbstractString, context::String; path_value=false)
    text = String(value)
    all(isascii, text) || _blocked("firewall", "$context is not ASCII")
    occursin('\0', text) && _blocked("firewall", "$context contains NUL")
    decoded = _percent_decode_fixed(text, context)
    all(character -> Int(character) >= 0x20 && Int(character) != 0x7f, decoded) ||
        _blocked("firewall", "$context contains a control character")
    compact = replace(lowercase(decoded), r"[^a-z0-9]+" => "")
    allowed_hash = _hash_re(text)
    forbidden = ("truth", "truths", "label", "labels", "benchmark", "benchmarks", "grade", "grader",
                 "graders", "count", "counts", "classcount", "classcounts", "expectedn", "nselected",
                 "nknnkn", "sequence", "sequences", "motif", "motifs", "composition", "compositions",
                 "classprior", "classpriors", "transition", "transitions", "position", "positions",
                 "order", "orders", "index", "indices", "rank", "ranks", "lobepos", "topk",
                 "runlength", "runlengths", "groundtruth", "gold")
    !allowed_hash && any(token -> occursin(token, compact), forbidden) &&
        _blocked("firewall", "$context contains a forbidden scientific concept")
    if path_value
        isabspath(decoded) && context != "root" &&
            _blocked("firewall", "$context must be root-relative")
        occursin('\\', decoded) && _blocked("firewall", "$context uses a backslash")
        parts = split(decoded, '/'; keepempty=true)
        any(part -> part in ("", ".", ".."), parts) && context != "root" &&
            _blocked("firewall", "$context has a noncanonical path value")
        occursin(r"(^|/)\.\.(/|$)", decoded) && _blocked("firewall", "$context escapes its root")
    end
    return nothing
end

function _audit_impl(root::String; features::String, feature_sha256::String, candidate_config::String,
                     model_config::String, universe_dir::String, execution_dir::String)
    canonical_root = _canonical_root(root)
    _progress_init!()
    authorities, source_specs = _source_authorities(canonical_root)
    candidate_binding = _config_binding(canonical_root)
    try
        return _audit_core(canonical_root; features=features, feature_sha256=feature_sha256,
                           candidate_config=candidate_config, model_config=model_config,
                           universe_dir=universe_dir, execution_dir=execution_dir)
    catch error
        progress = _REPORT_PROGRESS[]
        progress === nothing && return _report_from_error(error)
        stable = Binding[]
        for binding in progress.bindings
            try
                _verify(canonical_root, binding)
                push!(stable, binding)
            catch
                # A changed input is not safely acquired and must not be
                # claimed as part of a retained failure closure.
            end
        end
        stable_paths = Set(binding.relative for binding in stable)
        stable_specs = [spec for spec in progress.specs if spec.path in stable_paths]
        isempty(stable_specs) && !isempty(stable) &&
            push!(stable_specs, ArtifactSpec("input", first(stable).relative))
        safe_files = length(stable) == length(progress.bindings) ? progress.files : nothing
        return _report_from_error(error; scope=progress.scope, bindings=stable,
                                  specs=stable_specs, files=safe_files,
                                  metrics=progress.metrics)
    end
end

function audit_stability(root::AbstractString; features::AbstractString, feature_sha256::AbstractString,
                          candidate_config::AbstractString, model_config::AbstractString,
                          universe_dir::AbstractString, execution_dir::AbstractString)::AuditReport
    old_cache = _SNAPSHOT_CACHE[]
    try
        values = (String(root), String(features), String(feature_sha256), String(candidate_config),
                  String(model_config), String(universe_dir), String(execution_dir))
        for (value, context) in zip(values, ("root", "features", "feature_sha256", "candidate_config",
                                             "model_config", "universe_dir", "execution_dir"))
            _firewall_value(value, context; path_value=context != "feature_sha256")
        end
        _SNAPSHOT_CACHE[] = Dict{String,Binding}()
        _REPORT_PROGRESS[] = nothing
        return _audit_impl(values[1]; features=values[2], feature_sha256=values[3],
                           candidate_config=values[4], model_config=values[5],
                           universe_dir=values[6], execution_dir=values[7])
    catch error
        progress = _REPORT_PROGRESS[]
        if progress === nothing
            return _report_from_error(error)
        end
        stable = Binding[]
        canonical_root = try
            _canonical_root(String(root))
        catch
            ""
        end
        if !isempty(canonical_root)
            for binding in progress.bindings
                try
                    _verify(canonical_root, binding)
                    push!(stable, binding)
                catch
                end
            end
        end
        stable_paths = Set(binding.relative for binding in stable)
        stable_specs = [spec for spec in progress.specs if spec.path in stable_paths]
        isempty(stable_specs) && !isempty(stable) &&
            push!(stable_specs, ArtifactSpec("input", first(stable).relative))
        safe_files = length(stable) == length(progress.bindings) ? progress.files : nothing
        return _report_from_error(error; scope=progress.scope, bindings=stable,
                                  specs=stable_specs, files=safe_files,
                                  metrics=progress.metrics)
    finally
        _SNAPSHOT_CACHE[] = old_cache
        _REPORT_PROGRESS[] = nothing
    end
end

function report_files(report::AuditReport)
    Set(keys(report.files)) == Set(OUTPUT_NAMES) || throw(PublicationError("BLOCKED", "report_schema", "report payload set differs"))
    return Dict(name => copy(report.files[name]) for name in OUTPUT_NAMES)
end

function _output_directory(root::String, supplied::AbstractString)
    text = String(supplied)
    _firewall_value(text, "out_dir"; path_value=true)
    isempty(text) && throw(PublicationError("BLOCKED", "path_invalid", "output directory is empty"))
    isabspath(text) && throw(PublicationError("BLOCKED", "path_substitution", "output directory must be root-relative"))
    parts = split(text, '/'; keepempty=true)
    any(part -> part in ("", ".", ".."), parts) && throw(PublicationError("BLOCKED", "path_invalid", "output directory is not canonical"))
    absolute = normpath(joinpath(root, parts...))
    cursor = root
    for part in parts
        cursor = joinpath(cursor, part)
        islink(cursor) && throw(PublicationError("BLOCKED", "symlink_path", "output path contains a symlink component"))
    end
    parent = dirname(absolute)
    isdir(parent) || throw(PublicationError("BLOCKED", "missing_output_parent", "output parent does not exist"))
    islink(parent) && throw(PublicationError("BLOCKED", "symlink_path", "output parent is a symlink"))
    st = stat(parent)
    (st.mode & 0o170000) == 0o40000 || throw(PublicationError("BLOCKED", "special_file", "output parent is not a directory"))
    return absolute, text
end

function _safe_publication_member(path::String, context::String)
    isfile(path) && !islink(path) ||
        throw(PublicationError("COLLISION", "special_file", "$context is not a regular file"))
    lstat(path).nlink == 1 ||
        throw(PublicationError("COLLISION", "hardlink_path", "$context is a hardlink"))
    return nothing
end

function _publication_state(root::String, out_relative::String)
    absolute, _ = _output_directory(root, out_relative)
    ispath(absolute) || return (:absent, absolute, nothing, nothing)
    isdir(absolute) && !islink(absolute) ||
        throw(PublicationError("COLLISION", "publication_collision", "output is not a directory"))
    names = sort(readdir(absolute))
    output_names = sort(collect(OUTPUT_NAMES))
    if names == output_names
        for name in names
            _safe_publication_member(joinpath(absolute, name), "payload $name")
        end
        return (:payload_only, absolute, nothing, nothing)
    end
    index_names = filter(name -> occursin(r"^artifact-index-[0-9a-f]{64}\.tsv$", name), names)
    receipt_name = "gate-receipt-$(EXECUTION_CAMPAIGN)-$(EXECUTION_MODEL).toml"
    complete_names = sort(vcat(collect(OUTPUT_NAMES), index_names, [receipt_name]))
    length(index_names) == 1 && names == complete_names ||
        throw(PublicationError("COLLISION", "publication_collision", "output is neither absent, payload-only, nor complete"))
    for name in names
        _safe_publication_member(joinpath(absolute, name), "publication member $name")
    end
    return (:complete, absolute, only(index_names), receipt_name)
end

function _rename_noreplace(source::String, destination::String)
    Sys.islinux() || throw(PublicationError("PRECOMMIT", "atomic_publication_unsupported", "renameat2 is required"))
    result = try
        ccall(:renameat2, Cint, (Cint, Cstring, Cint, Cstring, Cuint), Cint(-100), source,
              Cint(-100), destination, Cuint(1))
    catch error
        throw(PublicationError("PRECOMMIT", "atomic_publication_unsupported", sprint(showerror, error)))
    end
    result == 0 && return nothing
    code = Int(Libc.errno())
    code in (17, 39) && throw(PublicationError("COLLISION", "publication_collision", "output appeared concurrently"))
    throw(PublicationError("PRECOMMIT", "atomic_install_failed", "renameat2 failed: $(Libc.strerror(code))"))
end

function _fsync_directory(path::String; committed=false)
    descriptor = ccall(:open, Cint, (Cstring, Cint), path, 0)
    descriptor >= 0 || throw(PublicationError(committed ? "COMMITTED_WARNING" : "PRECOMMIT", "fsync_failed", "could not open directory"))
    try
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            throw(PublicationError(committed ? "COMMITTED_WARNING" : "PRECOMMIT", "fsync_failed", "directory fsync failed"))
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
end

function _materialize_payloads(root::String, out_relative::String, payloads; failpoint=nothing)
    state, absolute, _, _ = _publication_state(root, out_relative)
    if state == :complete
        throw(PublicationError("COLLISION", "publication_collision", "complete publication must be independently replay-validated"))
    elseif state == :payload_only
        for name in OUTPUT_NAMES
            path = joinpath(absolute, name)
            _safe_publication_member(path, "published payload $name")
            read(path) == payloads[name] || throw(PublicationError("COLLISION", "publication_collision", "published payload differs"))
        end
        return absolute
    end
    stage = mktempdir(dirname(absolute); prefix=".structured-stability-stage-")
    stat(stage).device == stat(dirname(absolute)).device ||
        throw(PublicationError("PRECOMMIT", "cross_filesystem_publication", "stage is not on the destination filesystem"))
    installed = false
    try
        for name in OUTPUT_NAMES
            path = joinpath(stage, name)
            open(path, "w") do io
                write(io, payloads[name]); flush(io)
                ccall(:fsync, Cint, (Cint,), fd(io)) == 0 || throw(PublicationError("PRECOMMIT", "fsync_failed", "payload fsync failed"))
            end
            isfile(path) && !islink(path) && lstat(path).nlink == 1 ||
                throw(PublicationError("PRECOMMIT", "special_file", "staged payload is unsafe"))
        end
        _fsync_directory(stage)
        failpoint == :before_install && throw(PublicationError("PRECOMMIT", "interrupted_publication", "simulated interruption before install"))
        _rename_noreplace(stage, absolute)
        installed = true
        _fsync_directory(dirname(absolute); committed=true)
        return absolute
    finally
        !installed && ispath(stage) && rm(stage; recursive=true, force=true)
    end
end

const EVIDENCE_FAILPOINTS = Set{Union{Nothing,Symbol}}((
    nothing, :before_index_install, :after_index_install, :before_receipt_install,
    :index_directory_fsync, :receipt_directory_fsync,
))

function _validate_publication_request(report::AuditReport, failpoint)
    failpoint in union(EVIDENCE_FAILPOINTS, Set{Union{Nothing,Symbol}}((:before_install,))) ||
        throw(PublicationError("BLOCKED", "invalid_failpoint", "unsupported publication failpoint"))
    seen = Dict{String,String}()
    config_count = 0; source_count = 0; input_count = 0
    for spec in report.specs
        isempty(spec.path) && throw(PublicationError("BLOCKED", "artifact_role", "publication spec path is empty"))
        haskey(seen, spec.path) &&
            throw(PublicationError("BLOCKED", "artifact_role", "publication specs contain a duplicate path"))
        seen[spec.path] = spec.role
        spec.role == "config" && (config_count += 1)
        spec.role == "source" && (source_count += 1)
        spec.role == "input" && (input_count += 1)
    end
    config_count == 1 || throw(PublicationError("BLOCKED", "missing_bindings", "publication requires exactly one config"))
    source_count >= 1 || throw(PublicationError("BLOCKED", "missing_bindings", "publication requires a source closure"))
    input_count >= 1 || throw(PublicationError("BLOCKED", "missing_bindings", "publication requires an input closure"))
    return nothing
end

function _acquire_outer_publication_lock(root::String, out_relative::String)
    absolute, _ = _output_directory(root, out_relative)
    lock_path = joinpath(dirname(absolute), ".structured-stability-$(basename(absolute)).publication-lock")
    islink(lock_path) && throw(PublicationError("COLLISION", "symlink_path", "outer publication lock is a symlink"))
    ispath(lock_path) && throw(PublicationError("COLLISION", "publication_collision", "outer publication lock is held"))
    try
        mkdir(lock_path; mode=0o700)
    catch error
        ispath(lock_path) && throw(PublicationError("COLLISION", "publication_collision", "outer publication lock was acquired concurrently"))
        rethrow(error)
    end
    token = bytes2hex(rand(UInt8, 32))
    owner_path = joinpath(lock_path, "owner")
    try
        open(owner_path, "w") do io
            write(io, token); flush(io)
            ccall(:fsync, Cint, (Cint,), fd(io)) == 0 ||
                throw(PublicationError("PRECOMMIT", "fsync_failed", "outer publication lock fsync failed"))
        end
    catch error
        candidate = (path=lock_path, owner=owner_path, token=token)
        _release_outer_publication_lock(candidate)
        rethrow(error)
    end
    return (path=lock_path, owner=owner_path, token=token)
end

function _release_outer_publication_lock(lock)
    lock === nothing && return nothing
    owns = try
        isdir(lock.path) && !islink(lock.path) && isfile(lock.owner) && !islink(lock.owner) &&
            read(lock.owner, String) == lock.token
    catch
        false
    end
    owns && rm(lock.path; recursive=true, force=true)
    return nothing
end

function _publication_plan(out_relative::String, report::AuditReport)
    terminal = report.publication_status
    terminal in StructuredEvidence.TERMINAL_STATUSES ||
        throw(PublicationError("BLOCKED", "publication_blocked", "report publication status is not terminal"))
    reason_code = report.publication_reason
    reason_code in ALLOWED_REASON_CODES && occursin(r"^[a-z][a-z0-9_]*$", reason_code) ||
        throw(PublicationError("BLOCKED", "publication_blocked", "publication reason is not in the closed vocabulary"))
    specs = copy(report.specs)
    append!(specs, [ArtifactSpec("expected_rows", "$out_relative/expected_rows.tsv")])
    append!(specs, [ArtifactSpec("output", "$out_relative/$name") for name in OUTPUT_NAMES if name != "expected_rows.tsv"])
    unique_specs = Dict{String,ArtifactSpec}()
    for spec in specs
        if haskey(unique_specs, spec.path)
            unique_specs[spec.path].role == spec.role ||
                throw(PublicationError("BLOCKED", "artifact_role", "publication artifact role differs"))
        else
            unique_specs[spec.path] = spec
        end
    end
    exact_specs = [unique_specs[path] for path in sort(collect(keys(unique_specs)))]
    command = _execution_command()
    return (terminal=terminal, reason_code=reason_code, reasons=[reason_code], specs=exact_specs,
            command=command)
end

function _validate_complete_publication(root::String, out_relative::String, payloads, report::AuditReport,
                                        plan, index_name::String, receipt_name::String)
    receipt_relative = "$out_relative/$receipt_name"
    index_relative = "$out_relative/$index_name"
    absolute_output = joinpath(root, split(out_relative, '/')...)
    for name in OUTPUT_NAMES
        read(joinpath(absolute_output, name)) == payloads[name] ||
            throw(PublicationError("COLLISION", "publication_collision",
                                   "complete publication payload differs from the requested report"))
    end
    receipt_bytes = _snapshot(root, receipt_relative, "published gate receipt").bytes
    receipt_sha = _sha(receipt_bytes)
    command_identity = StructuredEvidence.command_identity_sha256(root, plan.command)
    validation = StructuredEvidence.validate_evidence_bundle(
        root, receipt_relative;
        expected_campaign_id=EXECUTION_CAMPAIGN,
        expected_model_id=EXECUTION_MODEL,
        expected_receipt_sha256=receipt_sha,
        expected_command_identity_sha256=command_identity,
    )
    validation.valid && validation.identity_verified ||
        throw(PublicationError("COLLISION", "publication_validation_failed",
                               "complete publication receipt is not independently valid"))
    absolute_index = joinpath(root, split(index_relative, '/')...)
    actual_index_bytes = read(absolute_index)
    expected_index = try
        StructuredEvidence.build_artifact_index(root, plan.specs;
            campaign_id=EXECUTION_CAMPAIGN, model_id=EXECUTION_MODEL,
            exact_command=plan.command, generated_at=EXECUTION_GENERATED_AT,
            terminal_status=plan.terminal, reason_codes=plan.reasons, reasons=plan.reasons)
    catch error
        throw(PublicationError("COLLISION", "publication_validation_failed",
                              "complete publication cannot be rebuilt: $(sprint(showerror, error))"))
    end
    expected_index_bytes = StructuredEvidence.artifact_index_bytes(expected_index)
    actual_index_bytes == expected_index_bytes ||
        throw(PublicationError("COLLISION", "publication_collision",
                               "complete publication index differs from the requested report"))
    _sha(actual_index_bytes) == replace(replace(index_name, r"^artifact-index-" => ""), r"\.tsv$" => "") ||
        throw(PublicationError("COLLISION", "publication_collision", "index filename is not content addressed"))
    receipt = StructuredEvidence.parse_gate_receipt(receipt_bytes)
    receipt.metadata["exact_command"] == plan.command && receipt.terminal_status == plan.terminal &&
        receipt.reason_codes == plan.reasons && receipt.reasons == plan.reasons ||
        throw(PublicationError("COLLISION", "publication_validation_failed",
                               "complete publication receipt does not match the requested report"))
    return validation
end

function publish_report(root::AbstractString, out_dir::AbstractString, report::AuditReport;
                        failpoint::Union{Nothing,Symbol}=nothing)
    _validate_publication_request(report, failpoint)
    _firewall_value(String(root), "root"; path_value=true)
    _firewall_value(String(out_dir), "out_dir"; path_value=true)
    canonical_root = _canonical_root(String(root))
    payloads = report_files(report)
    _, out_relative = _output_directory(canonical_root, String(out_dir))
    for binding in report.bindings
        _verify(canonical_root, binding)
    end
    plan = _publication_plan(out_relative, report)
    outer_lock = nothing
    try
        outer_lock = _acquire_outer_publication_lock(canonical_root, out_relative)
        state, _, index_name, receipt_name = _publication_state(canonical_root, out_relative)
        if state == :complete
            return _validate_complete_publication(canonical_root, out_relative, payloads, report,
                                                  plan, index_name, receipt_name)
        end
        _materialize_payloads(canonical_root, out_relative, payloads; failpoint=failpoint)
        state_after, _, index_after, receipt_after = _publication_state(canonical_root, out_relative)
        if state_after == :complete
            return _validate_complete_publication(canonical_root, out_relative, payloads, report,
                                                  plan, index_after, receipt_after)
        end
        state_after == :payload_only ||
            throw(PublicationError("PRECOMMIT", "publication_collision", "payload installation did not reach payload-only state"))
        for binding in report.bindings
            _verify(canonical_root, binding)
        end
        state_before_delegate, _, index_before_delegate, receipt_before_delegate =
            _publication_state(canonical_root, out_relative)
        if state_before_delegate == :complete
            return _validate_complete_publication(canonical_root, out_relative, payloads, report,
                                                  plan, index_before_delegate, receipt_before_delegate)
        end
        state_before_delegate == :payload_only ||
            throw(PublicationError("PRECOMMIT", "publication_collision", "publication state changed before delegation"))
        evidence_failpoint = failpoint == :before_install ? nothing : failpoint
        result = StructuredEvidence.publish_evidence_bundle(
            canonical_root, out_relative;
            campaign_id=EXECUTION_CAMPAIGN,
            model_id=EXECUTION_MODEL,
            artifacts=plan.specs,
            exact_command=plan.command,
            generated_at=EXECUTION_GENERATED_AT,
            terminal_status=plan.terminal,
            reason_codes=plan.reasons,
            reasons=plan.reasons,
            failpoint=evidence_failpoint,
        )
        evidence_failpoint !== nothing && !result.valid && return result
        result.valid && result.identity_verified && result.campaign_id == EXECUTION_CAMPAIGN &&
            result.model_id == EXECUTION_MODEL ||
            throw(PublicationError("PRECOMMIT", "publication_identity", "public evidence publisher returned an invalid identity"))
        expected_command_identity = StructuredEvidence.command_identity_sha256(canonical_root, plan.command)
        validation = StructuredEvidence.validate_evidence_bundle(
            canonical_root, result.receipt_path;
            expected_campaign_id=EXECUTION_CAMPAIGN,
            expected_model_id=EXECUTION_MODEL,
            expected_receipt_sha256=result.receipt_sha256,
            expected_command_identity_sha256=expected_command_identity,
        )
        validation.valid && validation.identity_verified && validation.campaign_id == EXECUTION_CAMPAIGN &&
            validation.model_id == EXECUTION_MODEL ||
            throw(PublicationError("COMMITTED_WARNING", "publication_validation_failed", "independent evidence validation failed"))
        return result
    finally
        _release_outer_publication_lock(outer_lock)
    end
end

const VALUE_FLAGS = Set([
    "--root", "--features", "--feature-sha256", "--candidate-config", "--model-config", "--universe-dir",
    "--execution-dir", "--out-dir",
])
const ALLOWED_FLAGS = union(VALUE_FLAGS, Set(["--help"]))

function _cli_text_safe(text::String)
    all(isascii, text) || throw(ArgumentError("CLI text is not ASCII"))
    occursin('\0', text) && throw(ArgumentError("CLI text contains NUL"))
    return nothing
end

function _parse_arguments(arguments::Vector{String})
    values = Dict{String,String}(); index = 1
    while index <= length(arguments)
        argument = arguments[index]; _cli_text_safe(argument)
        _firewall_value(argument, "CLI token")
        startswith(argument, "--") || throw(ArgumentError("unexpected positional argument"))
        pair = split(argument, '='; limit=2); flag = pair[1]
        flag in ALLOWED_FLAGS || throw(ArgumentError("unknown option: $flag"))
        haskey(values, flag) && throw(ArgumentError("duplicate option: $flag"))
        if flag == "--help"
            length(pair) == 1 || throw(ArgumentError("--help takes no value"))
            values[flag] = "true"; index += 1; continue
        end
        value = if length(pair) == 2
            isempty(pair[2]) && throw(ArgumentError("missing value for $flag")); pair[2]
        else
            index < length(arguments) || throw(ArgumentError("missing value for $flag"))
            index += 1; arguments[index]
        end
        _cli_text_safe(String(value)); _firewall_value(String(value), "CLI value")
        values[flag] = String(value); index += 1
    end
    haskey(values, "--help") && length(values) == 1 && return values
    for flag in VALUE_FLAGS
        haskey(values, flag) || throw(ArgumentError("required option is absent: $flag"))
    end
    return values
end

function _usage(io::IO=stdout)
    println(io, "Usage: julia --project=. test/audit_structured_unit_assignment_stability.jl \\")
    println(io, "  --root PATH --features PATH --feature-sha256 HEX \\")
    println(io, "  --candidate-config PATH --model-config PATH --universe-dir PATH \\")
    println(io, "  --execution-dir PATH --out-dir PATH")
end

function main(arguments::Vector{String}=copy(ARGS))::Int
    try
        values = _parse_arguments(arguments)
        if haskey(values, "--help")
            _usage(); return 0
        end
        report = audit_stability(values["--root"]; features=values["--features"],
                                 feature_sha256=values["--feature-sha256"],
                                 candidate_config=values["--candidate-config"],
                                 model_config=values["--model-config"],
                                 universe_dir=values["--universe-dir"],
                                 execution_dir=values["--execution-dir"])
        result = publish_report(values["--root"], values["--out-dir"], report)
        println("computed_status=", report.computed_status)
        println("terminal_status=", result.terminal_status)
        return result.valid ? 0 : 2
    catch error
        showerror(stderr, error); println(stderr)
        _usage(stderr)
        return 2
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

end # module StructuredStabilityAudit
