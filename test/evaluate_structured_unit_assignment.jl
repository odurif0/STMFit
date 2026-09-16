#!/usr/bin/env julia

module StructuredUnitAssignmentEvaluator
using Printf, Random, SHA, Statistics, TOML, LinearAlgebra
include(joinpath(@__DIR__, "lib", "structured_unit_assignment.jl"))
include(joinpath(@__DIR__, "lib", "structured_assignment", "edge_admission.jl"))
include(joinpath(@__DIR__, "lib", "structured_assignment", "chain_inference.jl"))
const T8 = StructuredUnitAssignment
const T11 = StructuredEdgeAdmission
const T12 = StructuredChainInference

export EvaluatorError, ViewJointEvidence, NodeEvidence, EdgeNullEvidence,
       GraphEvidence, InnerFoldEvidence, OuterFoldEvidence, EvaluationInput,
       EvaluatorReport, evaluate, report_files,
       UNARY_SELECTION_REASON_CODES, UnarySelectionReference,
       UnarySelectionFoldEvidence, UnarySelectionScanEvidence,
       UnarySelectionDateEvidence, UnarySelectionBootstrapEvidence,
       UnarySelectionDecision, UnarySelectionEvidence,
       produce_unary_selection_evidence

struct EvaluatorError <: Exception
    status::Symbol
    reason::Symbol
    message::String
    function EvaluatorError(status::Symbol, reason::Symbol, message::String)
        status in (:BLOCKED, :FAIL) || throw(ArgumentError("invalid evaluator error status"))
        new(status, reason, message)
    end
end

Base.showerror(io::IO, error::EvaluatorError) =
    print(io, error.status, " [", error.reason, "]: ", error.message)

struct ViewJointEvidence
    name::String
    c1_log_joint::NTuple{2,Float64}
    c2_log_joint::NTuple{2,Float64}
end

struct NodeEvidence
    file::String
    date::String
    lobe::Int
    t_nm::Float64
    views::Vector{ViewJointEvidence}
end

struct EdgeNullEvidence
    file::String
    date::String
    left_lobe::Int
    right_lobe::Int
    segment_id::String
    null_loglik::Float64
end

struct GraphEvidence
    model::String
    selected_unary_fit_sha256::String
    t11_fit_sha256::String
    common_null_sha256::String
    t11_gate_status::Symbol
    t11_gate_sha256::String
    t12_status::Symbol
    t12_reason::Symbol
    t12_result_sha256::String
    scan_logz::Dict{String,Float64}
    node_p1::Dict{Tuple{String,Int},Float64}
end

struct InnerFoldEvidence
    heldout_date::String
    training_dates::Vector{String}
    c1_status::Symbol
    c1_reason::Symbol
    c1_fit_sha256::String
    c2_status::Symbol
    c2_reason::Symbol
    c2_fit_sha256::String
    nodes::Vector{NodeEvidence}
    eligible_edge_counts::Dict{String,Int}
    c1_edge_fit_sha256::String
    c2_edge_fit_sha256::String
end

InnerFoldEvidence(heldout_date, training_dates, c1_status, c1_reason, c1_fit_sha256,
                  c2_status, c2_reason, c2_fit_sha256, nodes, eligible_edge_counts) =
    InnerFoldEvidence(heldout_date, training_dates, c1_status, c1_reason,
                      c1_fit_sha256, c2_status, c2_reason, c2_fit_sha256,
                      nodes, eligible_edge_counts, "", "")

struct OuterFoldEvidence
    outer_date::String
    training_dates::Vector{String}
    inner_folds::Vector{InnerFoldEvidence}
    c1_status::Symbol
    c1_reason::Symbol
    c1_fit_sha256::String
    c2_status::Symbol
    c2_reason::Symbol
    c2_fit_sha256::String
    common_null_sha256::String
    eligible_edge_count::Int
    nodes::Vector{NodeEvidence}
    edges::Vector{EdgeNullEvidence}
    graphs::Dict{String,GraphEvidence}
end

OuterFoldEvidence(outer_date, training_dates, inner_folds,
                  c1_status, c1_reason, c1_fit_sha256,
                  c2_status, c2_reason, c2_fit_sha256,
                  common_null_sha256, nodes, edges, graphs) =
    OuterFoldEvidence(outer_date, training_dates, inner_folds,
                      c1_status, c1_reason, c1_fit_sha256,
                      c2_status, c2_reason, c2_fit_sha256,
                      common_null_sha256, length(edges), nodes, edges, graphs)

struct EvaluationInput
    authority_sha256::String
    universe_sha256::String
    scan_dates::Vector{Tuple{String,String}}
    folds::Vector{OuterFoldEvidence}
end

struct EvaluatorReport
    status::Symbol
    reason::Symbol
    folds::Vector{NamedTuple}
    nodes::Vector{NamedTuple}
    scans::Vector{NamedTuple}
    dates::Vector{NamedTuple}
    signs::Vector{NamedTuple}
    events::Vector{NamedTuple}
    metrics::NamedTuple
    authority_sha256::String
    universe_sha256::String
    source_sha256::String
    bootstrap_values::Vector{Float64}
    result_sha256::String
    provenance_sha256::String
end

# T13 is deliberately a unary-only, typed evidence seam.  These records carry
# identities and observations, but no labels, expected counts, priors, or model
# selection inputs.  Keep the public records independent of the report TSV
# structs above: the evaluator report remains the frozen historical product.
const UNARY_SELECTION_REASON_CODES = (
    :ok, :c2_unavailable_fallback_c1, :nonpositive_date, :bootstrap_not_positive,
    :scope_unavailable,
    :unary_numerical_failure, :denominator_nonpositive, :bootstrap_replay_mismatch,
    :authority_or_factor_mismatch, :schema_mismatch,
    :incomplete_or_stale_evidence, :all_views_missing_incident_edge,
    :selection_reference_mismatch,
)

struct UnarySelectionReference
    model::String
    fit_role::String
    scope::String
    outer_date::String
    target_date::String
    fit_sha256::String
    unary_fit_sha256::String
end

struct UnarySelectionFoldEvidence
    decision_scope::String
    outer_date::String
    heldout_date::String
    training_dates::Vector{String}
    c1_status::Symbol
    c1_reason::Symbol
    c1_reference::Union{Nothing,UnarySelectionReference}
    c2_status::Symbol
    c2_reason::Symbol
    c2_reference::Union{Nothing,UnarySelectionReference}
end

struct UnarySelectionScanEvidence
    decision_scope::String
    outer_date::String
    date::String
    file::String
    node_count::Int
    eligible_edge_count::Int
    denominator::Int
    c1_log_evidence::Float64
    c2_log_evidence::Float64
    delta::Float64
end

struct UnarySelectionDateEvidence
    decision_scope::String
    outer_date::String
    date::String
    scan_count::Int
    date_mean::Float64
    strict_positivity::Bool
end

struct UnarySelectionBootstrapEvidence
    scope::String
    outer_date::String
    seed::Int
    mean_delta::Float64
end

struct UnarySelectionDecision
    scope::String
    outer_date::String
    status::Symbol
    reason::Symbol
    selected_model::String
    selected_reference::Union{Nothing,UnarySelectionReference}
    c1_reference::Union{Nothing,UnarySelectionReference}
    c2_reference::Union{Nothing,UnarySelectionReference}
    fold_count::Int
    scan_count::Int
    bootstrap_count::Int
    bootstrap_lower::Float64
    every_date_positive::Bool
    decision_hash::String
    evidence_hash::String
    upstream_reason::Union{Nothing,Symbol}
end

# The names used by receipts in adjacent seams are *_sha256.  Keep the short
# names as the storage contract and provide read-only aliases for callers that
# use the longer spelling.
function Base.getproperty(decision::UnarySelectionDecision, name::Symbol)
    name === :model && return getfield(decision, :selected_model)
    name === :hash && return getfield(decision, :decision_hash)
    name === :lower_bound && return getfield(decision, :bootstrap_lower)
    name === :candidate_references &&
        return (getfield(decision, :c1_reference), getfield(decision, :c2_reference))
    name === :decision_sha256 && return getfield(decision, :decision_hash)
    name === :evidence_sha256 && return getfield(decision, :evidence_hash)
    return getfield(decision, name)
end

function Base.getproperty(fold::UnarySelectionFoldEvidence, name::Symbol)
    name === :scope && return getfield(fold, :decision_scope)
    return getfield(fold, name)
end

function Base.getproperty(scan::UnarySelectionScanEvidence, name::Symbol)
    name === :scope && return getfield(scan, :decision_scope)
    return getfield(scan, name)
end

function Base.getproperty(date::UnarySelectionDateEvidence, name::Symbol)
    name === :strict_positive && return getfield(date, :strict_positivity)
    return getfield(date, name)
end

function Base.getproperty(row::UnarySelectionBootstrapEvidence, name::Symbol)
    name === :delta && return getfield(row, :mean_delta)
    return getfield(row, name)
end

struct UnarySelectionEvidence
    status::Symbol
    reason::Symbol
    upstream_reason::Union{Nothing,Symbol}
    outer_decisions::Vector{UnarySelectionDecision}
    final_decision::UnarySelectionDecision
    folds::Vector{UnarySelectionFoldEvidence}
    scans::Vector{UnarySelectionScanEvidence}
    dates::Vector{UnarySelectionDateEvidence}
    bootstrap::Vector{UnarySelectionBootstrapEvidence}
    bindings::NamedTuple
    authority_sha256::String
    evidence_hash::String
end

function Base.getproperty(evidence::UnarySelectionEvidence, name::Symbol)
    name === :outer && return getfield(evidence, :outer_decisions)
    name === :final && return getfield(evidence, :final_decision)
    name === :scan_evidence && return getfield(evidence, :scans)
    name === :date_evidence && return getfield(evidence, :dates)
    name === :bootstrap_evidence && return getfield(evidence, :bootstrap)
    bindings = getfield(evidence, :bindings)
    name in propertynames(bindings) && return getproperty(bindings, name)
    return getfield(evidence, name)
end

const _SCHEMA = "structured_" * "la" * "bel_free_unit_assignment_evaluator_v1"
const _VIEWS = ("base_local", "base_local+bwd_neg_com_t",
                "base_local+bwd_neg_diag45", "base_local+split_log_skew")
const _FILES = ("bootstrap.tsv", "dates.tsv", "events.tsv", "folds.tsv",
                "nodes.tsv", "receipt.toml", "scans.tsv", "signs.tsv",
                "summary.tsv")
const _BOOTSTRAP_SEEDS = 0:499
const _PRECEDENCE = Dict(:BLOCKED => 4, :FAIL => 3, :SKIPPED => 2, :PASS => 1)

const _FOLD_HEADER = [
    "schema_version", "outer_date", "training_dates", "c1_unary_fit_sha256",
    "c2_unary_fit_sha256", "inner_unary_status", "inner_unary_reason",
    "inner_unary_decision_sha256", "selected_model", "selected_unary_fit_sha256",
    "t11_gate_status", "t11_gate_sha256", "t11_outer_fit_sha256", "t12_status",
    "t12_reason", "t12_result_sha256", "graph_enabled", "fold_status", "fold_reason",
]
const _NODE_HEADER = [
    "schema_version", "outer_date", "date", "file", "lobe", "t_nm",
    "selected_model", "selected_unary_fit_sha256", "graph_enabled", "views_used",
    "c1_u", "selected_u", "p1", "hard_output", "entropy",
    "view_agreement_numerator", "view_agreement_denominator", "status", "reason",
]
const _SCAN_HEADER = [
    "schema_version", "outer_date", "date", "file", "selected_model",
    "graph_enabled", "node_count", "eligible_edge_count", "denominator",
    "c1_node_loglik", "selected_node_loglik", "common_null_loglik", "graph_logz",
    "c1_loglik", "meta_loglik", "d_s", "score_sha256", "status", "reason",
]
const _DATE_HEADER = ["schema_version", "date", "scan_count", "d_k", "status", "reason"]
const _BOOTSTRAP_HEADER = ["schema_version", "seed", "mean_d"]
const _SIGN_HEADER = ["schema_version", "mask", "signs", "t_epsilon", "in_upper_tail"]
const _EVENT_HEADER = [
    "schema_version", "ordinal", "outer_date", "event", "consumed",
    "terminal_status", "consequence", "reason", "formula_count",
    "bootstrap_count", "sign_count", "graph_count",
]
const _SUMMARY_HEADER = [
    "schema_version", "status", "reason", "date_count", "scan_count",
    "node_count", "eligible_edge_count", "t_obs", "sign_tail_count",
    "sign_denominator", "sign_p", "bootstrap_lower", "every_date_positive",
    "coverage", "entropy", "view_agreement_numerator", "view_agreement_denominator",
    "view_agreement", "unary_parameter_count", "residualizer_parameter_count",
    "null_parameter_count", "conditional_edge_parameter_count", "unary_meta_total",
    "graph_meta_total", "graph_enabled_fold_count",
]

_hash_bytes(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))
_hash_text(text::AbstractString) = _hash_bytes(Vector{UInt8}(codeunits(String(text))))
_f(value::Real) = isfinite(Float64(value)) ? @sprintf("%.17g", Float64(value)) : "NA"
_b(value::Bool) = value ? "true" : "false"
_status_text(value::Symbol) = String(value)

function _fail(status::Symbol, reason::Symbol, message::AbstractString)
    throw(EvaluatorError(status, reason, String(message)))
end

function _logsumexp(values::NTuple{2,Float64})::Float64
    maximum(values) == -Inf && return -Inf
    m = max(values[1], values[2])
    isfinite(m) || _fail(:FAIL, :unary_numerical_failure, "nonfinite unary joint")
    return m + log(exp(values[1] - m) + exp(values[2] - m))
end

function _logsumexp(values::AbstractVector{<:Real})::Float64
    isempty(values) && return -Inf
    m = maximum(Float64.(values))
    m == -Inf && return -Inf
    isfinite(m) || _fail(:FAIL, :unary_numerical_failure, "nonfinite log-sum-exp input")
    return m + log(sum(exp(Float64(value) - m) for value in values))
end

function _canonical_lines(lines)
    io = IOBuffer()
    for line in lines
        println(io, String(line))
    end
    return take!(io)
end

function _type7(values::Vector{Float64}, probability::Float64=0.025)::Float64
    isempty(values) && _fail(:FAIL, :bootstrap_replay_mismatch, "empty quantile population")
    ordered = sort(copy(values))
    h = 1.0 + (length(ordered) - 1) * probability
    lower = floor(Int, h)
    upper = ceil(Int, h)
    lower == upper && return ordered[lower]
    gamma = h - lower
    return (1.0 - gamma) * ordered[lower] + gamma * ordered[upper]
end

function _view_available(view::ViewJointEvidence, model::Symbol)::Bool
    joints = model == :C1 ? view.c1_log_joint : view.c2_log_joint
    any(isinf, joints) &&
        _fail(:FAIL, :unary_numerical_failure, "infinite unary joint")
    return isfinite(joints[1]) && isfinite(joints[2])
end

function _unary(node::NodeEvidence, model::Symbol)
    eta = Float64[]
    active = ViewJointEvidence[]
    c2_unavailable = false
    for view in node.views
        c1 = view.c1_log_joint
        c2 = view.c2_log_joint
        c1_missing = all(isnan, c1)
        c2_missing = all(isnan, c2)
        c1_partial = any(isnan, c1) && !c1_missing
        c2_partial = any(isnan, c2) && !c2_missing
        !(c1_partial || c2_partial) ||
            _fail(:FAIL, :unary_numerical_failure, "partial unary view evidence")
        any(isinf, c1) || any(isinf, c2) ?
            _fail(:FAIL, :unary_numerical_failure, "infinite unary view evidence") : nothing
        c1_missing && continue
        c2_missing && (c2_unavailable = true)
        model == :C2 && c2_missing && continue
        joints = model == :C1 ? c1 : c2
        all(isfinite, joints) || continue
        a = _logsumexp(joints)
        raw = exp(joints[2] - a)
        isfinite(raw) || _fail(:FAIL, :unary_numerical_failure, "nonfinite view posterior")
        clipped = clamp(raw, 1.0e-12, 1.0 - 1.0e-12)
        push!(eta, a + log1p(-clipped))
        push!(eta, a + log(clipped))
        push!(active, view)
    end
    isempty(active) && return (u=0.0, p1=0.5, views=0, active=active,
                               c2_unavailable=c2_unavailable)
    length(eta) == 2 * length(active) ||
        _fail(:FAIL, :unary_numerical_failure, "unary opinion pool width differs")
    eta0 = mean(eta[index] for index in 1:2:length(eta))
    eta1 = mean(eta[index] for index in 2:2:length(eta))
    u = _logsumexp((eta0, eta1))
    isfinite(u) || _fail(:FAIL, :unary_numerical_failure, "unary score is nonfinite")
    p1 = exp(eta1 - u)
    isfinite(p1) || _fail(:FAIL, :unary_numerical_failure, "unary posterior is nonfinite")
    return (u=u, p1=p1, views=length(active), active=active,
            c2_unavailable=c2_unavailable)
end

function _raw_view_posterior(view::ViewJointEvidence, model::Symbol)::Float64
    all(isnan, view.c1_log_joint) && return NaN
    _view_available(view, :C1) || return NaN
    model == :C2 && !_view_available(view, :C2) && return NaN
    joints = model == :C1 ? view.c1_log_joint : view.c2_log_joint
    _view_available(view, model) || return NaN
    return exp(joints[2] - _logsumexp(joints))
end

function _hard(p1::Float64)::String
    p1 > 0.5 && return "1"
    p1 < 0.5 && return "0"
    return "?"
end

function _entropy(p1::Float64)::Float64
    (p1 == 0.0 || p1 == 1.0) && return 0.0
    (0.0 <= p1 <= 1.0) || _fail(:FAIL, :unary_numerical_failure, "posterior outside [0,1]")
    return -(p1 * log2(p1) + (1.0 - p1) * log2(1.0 - p1))
end

function _agreement(node::NodeEvidence, model::Symbol)
    probabilities = Float64[]
    for view in node.views
        p = _raw_view_posterior(view, model)
        isfinite(p) && push!(probabilities, p)
    end
    numerator = 0
    denominator = 0
    for left in 1:length(probabilities), right in (left + 1):length(probabilities)
        denominator += 1
        ((probabilities[left] - 0.5) * (probabilities[right] - 0.5) > 0.0) &&
            (numerator += 1)
    end
    return numerator, denominator
end

function _validate_hash(text::AbstractString, context::AbstractString)
    occursin(r"^[0-9a-f]{64}$", String(text)) ||
        _fail(:BLOCKED, :authority_hash_mismatch, "$context is not a SHA-256")
end

function _validate_input(input::EvaluationInput)
    _validate_hash(input.authority_sha256, "authority hash")
    _validate_hash(input.universe_sha256, "universe hash")
    isempty(input.scan_dates) && _fail(:BLOCKED, :schema_mismatch, "scan universe is empty")
    scans = sort(copy(input.scan_dates); by=pair -> (pair[2], pair[1]))
    scans == input.scan_dates || _fail(:BLOCKED, :schema_mismatch, "scan dates are not canonical")
    length(unique(scans)) == length(scans) ||
        _fail(:BLOCKED, :schema_mismatch, "scan dates contain duplicates")
    all(!isempty(file) && occursin(r"^[0-9]{8}$", date) for (file, date) in scans) ||
        _fail(:BLOCKED, :schema_mismatch, "scan date identity is malformed")
    dates = Dict(file => date for (file, date) in scans)
    all_dates = sort(unique(last.(scans)))
    length(input.folds) == length(unique(fold.outer_date for fold in input.folds)) ||
        _fail(:BLOCKED, :schema_mismatch, "outer fold dates are duplicated")
    sort([fold.outer_date for fold in input.folds]) == all_dates ||
        _fail(:BLOCKED, :schema_mismatch, "outer fold coverage differs from universe")
    for fold in input.folds
        fold.outer_date in values(dates) ||
            _fail(:BLOCKED, :schema_mismatch, "outer fold date is absent from universe")
        expected_training = [date for date in all_dates if date != fold.outer_date]
        fold.training_dates == expected_training ||
            _fail(:BLOCKED, :schema_mismatch, "training dates are not canonical")
        fold.eligible_edge_count >= 0 ||
            _fail(:BLOCKED, :schema_mismatch, "eligible edge count is invalid")
        length(fold.edges) <= fold.eligible_edge_count ||
            _fail(:BLOCKED, :schema_mismatch, "edge evidence count exceeds eligible count")
        expected_outer_files = Set(file for (file, date) in scans if date == fold.outer_date)
        haskey(fold.graphs, "C1") && haskey(fold.graphs, "C2") &&
            length(fold.graphs) == 2 ||
            _fail(:BLOCKED, :schema_mismatch, "graph keys are not exactly C1 and C2")
        for node in fold.nodes
            haskey(dates, node.file) ||
                _fail(:BLOCKED, :schema_mismatch, "node file is absent from universe")
            dates[node.file] == node.date ||
                _fail(:BLOCKED, :schema_mismatch, "node date differs from universe")
            node.date == fold.outer_date ||
                _fail(:BLOCKED, :schema_mismatch, "outer node is not held out")
            node.lobe > 0 && isfinite(node.t_nm) ||
                _fail(:BLOCKED, :schema_mismatch, "node identity is invalid")
            names = [view.name for view in node.views]
            length(names) == length(unique(names)) ||
                _fail(:BLOCKED, :schema_mismatch, "node view names are duplicated")
            all(name in _VIEWS for name in names) ||
                _fail(:BLOCKED, :schema_mismatch, "node view name is not configured")
        end
        outer_keys = [(node.file, node.lobe) for node in fold.nodes]
        length(outer_keys) == length(unique(outer_keys)) ||
            _fail(:BLOCKED, :schema_mismatch, "outer nodes are duplicated")
        Set(node.file for node in fold.nodes) == expected_outer_files ||
            _fail(:BLOCKED, :schema_mismatch, "outer node scan coverage differs from universe")
        edges = Set{Tuple{String,Int,Int,String}}()
        nodekeys = Set(outer_keys)
        for edge in fold.edges
            haskey(dates, edge.file) ||
                _fail(:BLOCKED, :schema_mismatch, "edge file is absent from universe")
            dates[edge.file] == edge.date == fold.outer_date ||
                _fail(:BLOCKED, :schema_mismatch, "edge date differs from universe")
            (edge.file, edge.left_lobe) in nodekeys && (edge.file, edge.right_lobe) in nodekeys ||
                _fail(:BLOCKED, :schema_mismatch, "edge endpoint is absent")
            edge.left_lobe != edge.right_lobe && isfinite(edge.null_loglik) ||
                _fail(:FAIL, :common_null_numerical_failure, "edge null score is invalid")
            key = (edge.file, edge.left_lobe, edge.right_lobe, edge.segment_id)
            key in edges && _fail(:BLOCKED, :schema_mismatch, "edge is duplicated")
            push!(edges, key)
        end
        expected_inner_dates = fold.training_dates
        [inner.heldout_date for inner in fold.inner_folds] == expected_inner_dates ||
            _fail(:BLOCKED, :schema_mismatch, "inner heldout topology differs from outer training")
        for inner in fold.inner_folds
            inner.heldout_date in values(dates) ||
                _fail(:BLOCKED, :schema_mismatch, "inner date is absent from universe")
            expected_inner_training = [date for date in fold.training_dates
                                       if date != inner.heldout_date]
            inner.training_dates == expected_inner_training ||
                _fail(:BLOCKED, :schema_mismatch, "inner training topology differs")
            expected_inner_files = Set(file for (file, date) in scans
                                       if date == inner.heldout_date)
            for (file, count) in inner.eligible_edge_counts
                haskey(dates, file) ||
                    _fail(:BLOCKED, :schema_mismatch, "inner edge file is absent from universe")
                dates[file] == inner.heldout_date && count >= 0 ||
                    _fail(:BLOCKED, :schema_mismatch, "inner edge count is invalid")
            end
            for node in inner.nodes
                haskey(dates, node.file) ||
                    _fail(:BLOCKED, :schema_mismatch, "inner node file is absent from universe")
                dates[node.file] == node.date == inner.heldout_date ||
                    _fail(:BLOCKED, :schema_mismatch, "inner node date is invalid")
            end
            inner_keys = [(node.file, node.lobe) for node in inner.nodes]
            length(inner_keys) == length(unique(inner_keys)) ||
                _fail(:BLOCKED, :schema_mismatch, "inner nodes are duplicated")
            Set(node.file for node in inner.nodes) == expected_inner_files ||
                _fail(:BLOCKED, :schema_mismatch, "inner node scan coverage differs")
            Set(keys(inner.eligible_edge_counts)) == expected_inner_files ||
                _fail(:BLOCKED, :schema_mismatch, "inner edge-count scan coverage differs")
        end
    end
    return nothing
end

function _scan_groups(nodes::Vector{NodeEvidence}, edges::Vector{EdgeNullEvidence})
    grouped_nodes = Dict{String,Vector{NodeEvidence}}()
    for node in nodes
        push!(get!(grouped_nodes, node.file, NodeEvidence[]), node)
    end
    grouped_edges = Dict{String,Vector{EdgeNullEvidence}}()
    for edge in edges
        push!(get!(grouped_edges, edge.file, EdgeNullEvidence[]), edge)
    end
    return grouped_nodes, grouped_edges
end

function _selection_reference(model::String, fit_role::String, scope::String,
                              outer_date::String, target_date::String,
                              fit_sha256::String, unary_fit_sha256::String)
    return UnarySelectionReference(model, fit_role, scope, outer_date, target_date,
                                   fit_sha256, unary_fit_sha256)
end

function _inner_reference(fold::InnerFoldEvidence, model::String,
                          decision_scope::String, outer_date::String)
    if model == "C1"
        return _selection_reference("C1", "partition",
                                    decision_scope == "final" ? "final_lodo" : "outer_inner",
                                    decision_scope == "final" ? "NA" : outer_date,
                                    fold.heldout_date, fold.c1_edge_fit_sha256,
                                    fold.c1_fit_sha256)
    end
    return _selection_reference("C2", "partition",
                                decision_scope == "final" ? "final_lodo" : "outer_inner",
                                decision_scope == "final" ? "NA" : outer_date,
                                fold.heldout_date, fold.c2_edge_fit_sha256,
                                fold.c2_fit_sha256)
end

function _selection_fold(fold::InnerFoldEvidence, decision_scope::String,
                         outer_date::String)
    return UnarySelectionFoldEvidence(
        decision_scope, outer_date, fold.heldout_date, copy(fold.training_dates),
        fold.c1_status, fold.c1_reason,
        _inner_reference(fold, "C1", decision_scope, outer_date),
        fold.c2_status, fold.c2_reason,
        _inner_reference(fold, "C2", decision_scope, outer_date),
    )
end

function _selection_reason(status::Symbol, reason::Symbol)
    status == :PASS &&
        return reason in (:ok, :c2_unavailable_fallback_c1, :nonpositive_date,
                          :bootstrap_not_positive) ? reason : :ok
    status == :SKIPPED && return :scope_unavailable
    status == :FAIL &&
        return reason in (:unary_numerical_failure, :denominator_nonpositive,
                          :bootstrap_replay_mismatch) ? reason : :unary_numerical_failure
    status == :BLOCKED &&
        return reason in (:authority_or_factor_mismatch, :schema_mismatch,
                          :incomplete_or_stale_evidence,
                          :all_views_missing_incident_edge,
                          :selection_reference_mismatch) ? reason :
               :authority_or_factor_mismatch
    return :schema_mismatch
end

function _selection_decision(scope::String, outer_date::String, status::Symbol,
                             reason::Symbol, model::String,
                             selected_reference, c1_reference, c2_reference,
                             fold_count::Int, scan_count::Int,
                             bootstrap_count::Int, lower::Float64,
                             every_date_positive::Bool, decision_hash::String,
                             evidence_hash::String, upstream_reason=nothing)
    closed_status = status in (:PASS, :SKIPPED, :FAIL, :BLOCKED) ? status : :BLOCKED
    closed_reason = _selection_reason(closed_status, reason)
    if closed_status != :PASS
        model = ""
        selected_reference = nothing
    end
    deterministic_hash = _selection_decision_hash(
        scope, outer_date, closed_status, closed_reason, model,
        selected_reference, c1_reference, c2_reference, fold_count,
        scan_count, bootstrap_count, lower, every_date_positive,
        upstream_reason, decision_hash)
    deterministic_evidence_hash = _selection_decision_evidence_hash(
        evidence_hash, selected_reference, c1_reference, c2_reference)
    return UnarySelectionDecision(
        scope, outer_date, closed_status, closed_reason, model,
        selected_reference, c1_reference, c2_reference, fold_count, scan_count,
        bootstrap_count, lower, every_date_positive, deterministic_hash,
        deterministic_evidence_hash,
        upstream_reason === nothing ? nothing : Symbol(upstream_reason),
    )
end

function _selection_scan_hash(scan::UnarySelectionScanEvidence)
    return join((scan.decision_scope, scan.outer_date, scan.date, scan.file, scan.node_count,
                 scan.eligible_edge_count, scan.denominator,
                 _f(scan.c1_log_evidence), _f(scan.c2_log_evidence),
                 _f(scan.delta)), '\t')
end

function _selection_evidence_hash(folds, scans, dates, bootstrap)
    lines = String["schema=structured-unary-selection-evidence-v1"]
    append!(lines, join((fold.decision_scope, fold.outer_date, fold.heldout_date,
                         join(fold.training_dates, ','), fold.c1_status,
                         fold.c1_reason, fold.c1_reference === nothing ? "" :
                         _selection_reference_identity(fold.c1_reference), fold.c2_status,
                         fold.c2_reason, fold.c2_reference === nothing ? "" :
                         _selection_reference_identity(fold.c2_reference)), '\t') for fold in folds)
    append!(lines, _selection_scan_hash(scan) for scan in scans)
    append!(lines, join((date.decision_scope, date.outer_date, date.date, date.scan_count, _f(date.date_mean),
                         date.strict_positivity), '\t') for date in dates)
    append!(lines, join((row.scope, row.outer_date, row.seed, _f(row.mean_delta)), '\t')
            for row in bootstrap)
    return _length_digest(lines)
end

function _selection_reference_identity(reference)
    reference === nothing && return ""
    return join((reference.model, reference.fit_role, reference.scope,
                 reference.outer_date, reference.target_date,
                 reference.fit_sha256, reference.unary_fit_sha256), '\t')
end

function _selection_decision_hash(scope::String, outer_date::String,
                                  status::Symbol, reason::Symbol,
                                  selected_model::String, selected_reference,
                                  c1_reference, c2_reference, fold_count::Int,
                                  scan_count::Int, bootstrap_count::Int,
                                  lower::Float64, every_date_positive::Bool,
                                  upstream_reason, upstream_hash::String="")
    return _length_digest([
        "schema=structured-unary-selection-decision-v2",
        "scope=$scope",
        "outer_date=$outer_date",
        "status=$(String(status))",
        "reason=$(String(reason))",
        "selected_model=$selected_model",
        "selected_reference=$(_selection_reference_identity(selected_reference))",
        "c1_reference=$(_selection_reference_identity(c1_reference))",
        "c2_reference=$(_selection_reference_identity(c2_reference))",
        "fold_count=$(fold_count)",
        "scan_count=$(scan_count)",
        "bootstrap_count=$(bootstrap_count)",
        "bootstrap_lower=$(_f(lower))",
        "every_date_positive=$(_b(every_date_positive))",
        "upstream_reason=$(upstream_reason === nothing ? "" : String(upstream_reason))",
        "upstream_hash=$upstream_hash",
    ])
end

function _selection_decision_evidence_hash(upstream_hash::String,
                                           selected_reference,
                                           c1_reference, c2_reference)
    return _length_digest([
        "schema=structured-unary-selection-decision-evidence-v3",
        "upstream_evidence_hash=$upstream_hash",
        "selected_reference=$(_selection_reference_identity(selected_reference))",
        "c1_reference=$(_selection_reference_identity(c1_reference))",
        "c2_reference=$(_selection_reference_identity(c2_reference))",
    ])
end

function _inner_gate_evidence_impl(folds::Vector{InnerFoldEvidence};
                                   decision_scope::String="outer",
                                   outer_date::String="",
                                   c1_reference=nothing, c2_reference=nothing)
    decision_scope in ("outer", "final") ||
        _fail(:BLOCKED, :schema_mismatch, "invalid unary selection decision scope")
    if !isempty(folds)
        representative = first(sort(copy(folds); by=fold -> fold.heldout_date))
        c1_reference === nothing &&
            (c1_reference = _inner_reference(representative, "C1",
                                             decision_scope, outer_date))
        c2_reference === nothing &&
            (c2_reference = _inner_reference(representative, "C2",
                                             decision_scope, outer_date))
    end
    selection_folds = UnarySelectionFoldEvidence[
        _selection_fold(fold, decision_scope, outer_date) for fold in folds]
    sort!(selection_folds; by=fold -> (fold.heldout_date, fold.outer_date))
    empty_decision = () -> _selection_decision(
        decision_scope, outer_date, :SKIPPED, :scope_unavailable, "", nothing,
        c1_reference, c2_reference, length(folds), 0, 0, NaN, false, "",
        _selection_evidence_hash(selection_folds,
                                 UnarySelectionScanEvidence[],
                                 UnarySelectionDateEvidence[],
                                 UnarySelectionBootstrapEvidence[]),
    )
    isempty(folds) && return (
        decision=empty_decision(), folds=selection_folds,
        scans=UnarySelectionScanEvidence[], dates=UnarySelectionDateEvidence[],
        bootstrap=UnarySelectionBootstrapEvidence[])

    all_c1 = true
    c2_available = true
    for fold in folds
        fold.c1_status == :PASS && continue
        all_c1 = false
        fold.c1_status == :SKIPPED && continue
        fold.c1_status == :FAIL &&
            _fail(:FAIL, :unary_numerical_failure, "inner C1 failed")
        fold.c1_status == :BLOCKED &&
            _fail(:BLOCKED, :authority_or_factor_mismatch, "inner C1 blocked")
        _fail(:BLOCKED, :schema_mismatch, "unknown inner C1 status")
    end
    if !all_c1
        decision = _selection_decision(
            decision_scope, outer_date, :SKIPPED, :scope_unavailable, "", nothing,
            c1_reference, c2_reference, length(folds), 0, 0, NaN, false, "",
            _selection_evidence_hash(selection_folds,
                                     UnarySelectionScanEvidence[],
                                     UnarySelectionDateEvidence[],
                                     UnarySelectionBootstrapEvidence[]),
        )
        return (decision=decision, folds=selection_folds,
                scans=UnarySelectionScanEvidence[], dates=UnarySelectionDateEvidence[],
                bootstrap=UnarySelectionBootstrapEvidence[])
    end
    for fold in folds
        if fold.c2_status == :SKIPPED
            c2_available = false
        elseif fold.c2_status == :FAIL
            _fail(:FAIL, :unary_numerical_failure, "inner C2 failed")
        elseif fold.c2_status == :BLOCKED
            _fail(:BLOCKED, :authority_or_factor_mismatch, "inner C2 blocked")
        elseif fold.c2_status != :PASS
            _fail(:BLOCKED, :schema_mismatch, "unknown inner C2 status")
        end
    end
    any(_unary(node, :C2).c2_unavailable for fold in folds for node in fold.nodes) &&
        (c2_available = false)
    if !c2_available
        selected = c1_reference
        decision = _selection_decision(
            decision_scope, outer_date, :PASS, :c2_unavailable_fallback_c1, "C1",
            selected, c1_reference, c2_reference, length(folds), 0, 0, NaN, false,
            "", _selection_evidence_hash(selection_folds,
                                          UnarySelectionScanEvidence[],
                                          UnarySelectionDateEvidence[],
                                          UnarySelectionBootstrapEvidence[]),
        )
        return (decision=decision, folds=selection_folds,
                scans=UnarySelectionScanEvidence[], dates=UnarySelectionDateEvidence[],
                bootstrap=UnarySelectionBootstrapEvidence[])
    end

    scan_values = Dict{String,Dict{String,Float64}}()
    scans = UnarySelectionScanEvidence[]
    decision_records = String["schema=structured-evaluator-inner-decision-v2"]
    for fold in folds
        nodes, _ = _scan_groups(fold.nodes, EdgeNullEvidence[])
        haskey(scan_values, fold.heldout_date) &&
            _fail(:BLOCKED, :schema_mismatch, "inner heldout date is duplicated")
        date_values = Dict{String,Float64}()
        for file in sort(collect(keys(nodes)))
            any(_unary(node, :C1).views == 0 for node in nodes[file]) &&
                get(fold.eligible_edge_counts, file, 0) > 0 &&
                _fail(:BLOCKED, :all_views_missing_incident_edge,
                      "all-view-missing inner node has an eligible edge")
            n = length(nodes[file])
            e = get(fold.eligible_edge_counts, file, 0)
            denominator = n + e
            denominator > 0 ||
                _fail(:FAIL, :denominator_nonpositive, "inner denominator is not positive")
            ordered_nodes = sort(nodes[file]; by=node -> (node.t_nm, node.lobe))
            # Do not rewrite this as (sum C2 - sum C1): the historical
            # decision arithmetic is a sum of per-node differences.
            delta = sum(_unary(node, :C2).u - _unary(node, :C1).u
                        for node in ordered_nodes; init=0.0) / denominator
            c1_log = sum(_unary(node, :C1).u for node in ordered_nodes; init=0.0)
            c2_log = sum(_unary(node, :C2).u for node in ordered_nodes; init=0.0)
            date_values[file] = delta
            push!(scans, UnarySelectionScanEvidence(
                decision_scope, outer_date, fold.heldout_date, file, n, e, denominator,
                c1_log, c2_log, delta))
            push!(decision_records,
                  "scan\t$(fold.heldout_date)\t$file\t$(_f(delta))\t$denominator")
        end
        Set(keys(nodes)) == Set(keys(fold.eligible_edge_counts)) ||
            _fail(:BLOCKED, :schema_mismatch, "inner edge scan set differs from nodes")
        scan_values[fold.heldout_date] = date_values
    end
    dates = UnarySelectionDateEvidence[]
    date_means = Float64[]
    for date in sort(collect(keys(scan_values)))
        values = scan_values[date]
        date_mean = mean(values[file] for file in sort(collect(keys(values))))
        push!(date_means, date_mean)
        push!(dates, UnarySelectionDateEvidence(decision_scope, outer_date, date,
                                                length(values), date_mean,
                                                date_mean > 0.0))
    end
    bootstrap_groups = Dict(date => [scan_values[date][file]
                                     for file in sort(collect(keys(scan_values[date])))]
                            for date in sort(collect(keys(scan_values))))
    bootstrap_values = _bootstrap_scan_values(bootstrap_groups)
    lower = _type7(bootstrap_values)
    bootstrap = UnarySelectionBootstrapEvidence[
        UnarySelectionBootstrapEvidence(decision_scope, outer_date, Int(seed), value)
        for (seed, value) in zip(_BOOTSTRAP_SEEDS, bootstrap_values)]
    sort!(scans; by=scan -> (scan.date, scan.file))
    append!(decision_records, "bootstrap\t$seed\t$(_f(value))"
            for (seed, value) in zip(_BOOTSTRAP_SEEDS, bootstrap_values))
    digest = _length_digest(decision_records)
    every_positive = all(>(0.0), date_means)
    reason = !every_positive ? :nonpositive_date : lower > 0.0 ? :ok : :bootstrap_not_positive
    model = reason == :ok ? "C2" : "C1"
    selected = model == "C2" ? c2_reference : c1_reference
    evidence_hash = _selection_evidence_hash(selection_folds, scans, dates, bootstrap)
    decision = _selection_decision(
        decision_scope, outer_date, :PASS, reason, model, selected,
        c1_reference, c2_reference, length(folds), length(scans),
        length(bootstrap), lower, every_positive, digest, evidence_hash,
    )
    return (decision=decision, folds=selection_folds, scans=scans, dates=dates,
            bootstrap=bootstrap)
end

function _inner_gate_evidence(folds::Vector{InnerFoldEvidence};
                              _propagate_errors::Bool=false, kwargs...)
    try
        return _inner_gate_evidence_impl(folds; kwargs...)
    catch error
        error isa EvaluatorError || rethrow()
        _propagate_errors && rethrow()
        scope = get(kwargs, :decision_scope, "outer")
        outer_date = get(kwargs, :outer_date, "")
        c1_reference = get(kwargs, :c1_reference, nothing)
        c2_reference = get(kwargs, :c2_reference, nothing)
        if !isempty(folds)
            representative = first(sort(copy(folds); by=fold -> fold.heldout_date))
            c1_reference === nothing &&
                (c1_reference = _inner_reference(representative, "C1",
                                                 scope, outer_date))
            c2_reference === nothing &&
                (c2_reference = _inner_reference(representative, "C2",
                                                 scope, outer_date))
        end
        folds_evidence = UnarySelectionFoldEvidence[
            _selection_fold(fold, scope, outer_date) for fold in folds]
        sort!(folds_evidence; by=fold -> (fold.heldout_date, fold.outer_date))
        decision = _selection_decision(
            scope, outer_date, error.status, error.reason, "", nothing,
            c1_reference, c2_reference, length(folds), 0, 0, NaN, false, "",
            _selection_evidence_hash(folds_evidence,
                                     UnarySelectionScanEvidence[],
                                     UnarySelectionDateEvidence[],
                                     UnarySelectionBootstrapEvidence[]),
            error.reason,
        )
        return (decision=decision, folds=folds_evidence,
                scans=UnarySelectionScanEvidence[], dates=UnarySelectionDateEvidence[],
                bootstrap=UnarySelectionBootstrapEvidence[])
    end
end

function _inner_gate(folds::Vector{InnerFoldEvidence})
    evidence = _inner_gate_evidence(folds; _propagate_errors=true)
    decision = evidence.decision
    return (status=decision.status, reason=decision.reason,
            model=decision.selected_model, hash=
            (decision.reason == :c2_unavailable_fallback_c1 ||
             decision.status == :SKIPPED ? "" : decision.decision_hash))
end

function _bootstrap_scan_values(scan_values::Dict{String,Vector{Float64}})::Vector{Float64}
    isempty(scan_values) && _fail(:FAIL, :bootstrap_replay_mismatch, "no scan values for bootstrap")
    dates = sort(collect(keys(scan_values)))
    all(!isempty(scan_values[file]) for file in dates) ||
        _fail(:FAIL, :bootstrap_replay_mismatch, "empty scan group for bootstrap")
    result = Float64[]
    for seed in _BOOTSTRAP_SEEDS
        rng = MersenneTwister(Int(seed))
        date_means = Float64[]
        for date in dates
            values = scan_values[date]
            sampled = [values[rand(rng, 1:length(values))] for _ in eachindex(values)]
            push!(date_means, mean(sampled))
        end
        push!(result, mean(date_means))
    end
    return result
end

function _inner_decision_hash(folds::Vector{InnerFoldEvidence}, decision)
    lines = String["schema=structured-evaluator-inner-decision-v1",
                   "model=$(decision.model)", "reason=$(decision.reason)"]
    for fold in sort(copy(folds); by=fold -> fold.heldout_date)
        push!(lines, join((fold.heldout_date, fold.c1_status, fold.c1_reason,
                           fold.c1_fit_sha256, fold.c2_status, fold.c2_reason,
                           fold.c2_fit_sha256), '\t'))
    end
    push!(lines, "decision=$(decision.hash)")
    return _hash_bytes(_canonical_lines(lines))
end

function _graph_state(fold::OuterFoldEvidence, model::String)
    haskey(fold.graphs, model) ||
        _fail(:BLOCKED, :schema_mismatch, "graph evidence is absent for $model")
    graph = fold.graphs[model]
    graph.model == model ||
        _fail(:BLOCKED, :authority_or_factor_mismatch, "graph model binding differs")
    expected_unary = model == "C1" ? fold.c1_fit_sha256 : fold.c2_fit_sha256
    graph.selected_unary_fit_sha256 == expected_unary ||
        _fail(:BLOCKED, :authority_or_factor_mismatch, "graph unary fit binding differs")
    graph.common_null_sha256 == fold.common_null_sha256 ||
        _fail(:BLOCKED, :authority_or_factor_mismatch, "graph common-null binding differs")
    graph.t11_gate_status in (:PASS, :SKIPPED, :FAIL, :BLOCKED) ||
        _fail(:BLOCKED, :schema_mismatch, "graph gate status is invalid")
    graph.t12_status in (:PASS, :SKIPPED, :FAIL, :BLOCKED) ||
        _fail(:BLOCKED, :schema_mismatch, "graph inference status is invalid")
    graph.t11_gate_status == :BLOCKED &&
        _fail(:BLOCKED, :authority_or_factor_mismatch, "graph admission is blocked")
    graph.t12_status == :BLOCKED &&
        _fail(:BLOCKED, :authority_or_factor_mismatch, "graph inference is blocked")
    return graph
end

function _selected_graph(fold::OuterFoldEvidence, selected_model::String)
    graph = _graph_state(fold, selected_model)
    graph.t11_gate_status == :FAIL &&
        _fail(:FAIL, :t11_numerical_failure, "selected edge admission failed")
    graph.t12_status == :FAIL &&
        _fail(:FAIL, :t12_failure, "selected chain inference failed")
    enabled = graph.model == selected_model &&
              graph.t11_gate_status == :PASS && graph.t12_status == :PASS
    return graph, enabled
end

function _score_fold(fold::OuterFoldEvidence, selected_model::String)
    selected_model in ("C1", "C2") ||
        _fail(:BLOCKED, :schema_mismatch, "selected unary model is invalid")
    fold.c1_status == :PASS || begin
        fold.c1_status == :FAIL && _fail(:FAIL, :unary_numerical_failure, "outer C1 failed")
        fold.c1_status == :BLOCKED && _fail(:BLOCKED, :authority_or_factor_mismatch, "outer C1 blocked")
        return (status=:SKIPPED, reason=:scope_unavailable, folds=NamedTuple[],
                nodes=NamedTuple[], scans=NamedTuple[], graph_enabled=false,
                edge_count=0)
    end
    effective_model = selected_model
    fallback = false
    if selected_model == "C2"
        fold.c2_status == :FAIL && _fail(:FAIL, :unary_numerical_failure, "outer C2 failed")
        fold.c2_status == :BLOCKED && _fail(:BLOCKED, :authority_or_factor_mismatch, "outer C2 blocked")
        fold.c2_status == :SKIPPED && (effective_model = "C1"; fallback = true)
        !fallback && any(_unary(node, :C2).c2_unavailable for node in fold.nodes) &&
            (effective_model = "C1"; fallback = true)
    end
    fold.eligible_edge_count > length(fold.edges) &&
        return (status=:SKIPPED, reason=:common_null_unavailable,
                effective_model=effective_model, folds=NamedTuple[],
                nodes=NamedTuple[], scans=NamedTuple[], graph_enabled=false,
                edge_count=0)
    graph, graph_enabled = _selected_graph(fold, effective_model)
    grouped_nodes, grouped_edges = _scan_groups(fold.nodes, fold.edges)
    Set(keys(grouped_nodes)) == Set(keys(grouped_edges)) ∪ Set(keys(grouped_nodes)) ||
        _fail(:BLOCKED, :schema_mismatch, "scan grouping failed")
    node_keys = Set((node.file, node.lobe) for node in fold.nodes)
    for node in fold.nodes
        _unary(node, :C1).views == 0 && _unary(node, :C2).views == 0 &&
            any(edge -> edge.file == node.file &&
                        (edge.left_lobe == node.lobe || edge.right_lobe == node.lobe),
                fold.edges) &&
            _fail(:BLOCKED, :all_views_missing_incident_edge,
                  "all-view-missing node has an eligible edge")
    end
    if graph_enabled
        Set(keys(graph.node_p1)) == node_keys ||
            _fail(:BLOCKED, :t12_reference_mismatch, "graph node cardinality differs")
        all(isfinite(value) && 0.0 <= value <= 1.0 for value in values(graph.node_p1)) ||
            _fail(:FAIL, :t12_failure, "graph marginal is invalid")
        Set(keys(graph.scan_logz)) == Set(keys(grouped_nodes)) ||
            _fail(:BLOCKED, :t12_reference_mismatch, "graph scan cardinality differs")
        all(isfinite(value) for value in values(graph.scan_logz)) ||
            _fail(:FAIL, :t12_failure, "graph normalizer is invalid")
    end
    score_rows = NamedTuple[]
    node_rows = NamedTuple[]
    scan_rows = NamedTuple[]
    edge_count = 0
    for file in sort(collect(keys(grouped_nodes)))
        nodes = sort(copy(grouped_nodes[file]); by=node -> (node.t_nm, node.lobe))
        edges = sort(copy(get(grouped_edges, file, EdgeNullEvidence[]));
                     by=edge -> (edge.left_lobe, edge.right_lobe, edge.segment_id))
        edge_count += length(edges)
        for node in nodes
            c1 = _unary(node, :C1)
            selected = _unary(node, Symbol(effective_model))
            p1 = graph_enabled ? graph.node_p1[(node.file, node.lobe)] : selected.p1
            numerator, denominator = _agreement(node, Symbol(effective_model))
            reason = selected.views == 0 ? :all_views_missing : :ok
            push!(node_rows, (
                outer_date=fold.outer_date, date=node.date, file=node.file,
                lobe=node.lobe, t_nm=node.t_nm, selected_model=effective_model,
                selected_fit=effective_model == "C1" ? fold.c1_fit_sha256 : fold.c2_fit_sha256,
                graph_enabled=graph_enabled, views_used=selected.views, c1_u=c1.u,
                selected_u=selected.u, p1=p1, hard=_hard(p1), entropy=_entropy(p1),
                agreement_numerator=numerator, agreement_denominator=denominator,
                status=:PASS, reason=reason,
            ))
        end
        node_log_c1 = sum(_unary(node, :C1).u for node in nodes)
        node_log_selected = sum(_unary(node, Symbol(effective_model)).u for node in nodes)
        null_log = sum((edge.null_loglik for edge in edges); init=0.0)
        denominator = length(nodes) + length(edges)
        denominator > 0 || _fail(:FAIL, :denominator_nonpositive, "outer scan denominator is not positive")
        graph_logz = graph_enabled ? graph.scan_logz[file] : 0.0
        c1_loglik = node_log_c1 + null_log
        meta_loglik = node_log_selected + null_log + (graph_enabled ? graph_logz : 0.0)
        d_s = (meta_loglik - c1_loglik) / denominator
        isfinite(d_s) || _fail(:FAIL, :numerical_failure, "scan difference is nonfinite")
        score_records = String[
            "schema=structured-evaluator-scan-score-v2", fold.outer_date, file,
            effective_model, effective_model == "C1" ? fold.c1_fit_sha256 : fold.c2_fit_sha256,
            fold.common_null_sha256, graph.t11_fit_sha256, graph.t11_gate_sha256,
            String(graph.t12_status), graph.t12_result_sha256, _f(d_s), string(denominator),
        ]
        append!(score_records, begin
            selected_unary = _unary(node, Symbol(effective_model))
            p = graph_enabled ? graph.node_p1[(node.file, node.lobe)] : selected_unary.p1
            "node\t$(node.file)\t$(node.lobe)\t$(_f(node.t_nm))\t$(_f(_unary(node, :C1).u))\t$(_f(selected_unary.u))\t$(_f(p))"
        end for node in nodes)
        append!(score_records, "edge\t$(edge.file)\t$(edge.left_lobe)\t$(edge.right_lobe)\t$(edge.segment_id)\t$(_f(edge.null_loglik))"
                for edge in edges)
        score_sha = _hash_bytes(_canonical_lines(score_records))
        push!(scan_rows, (
            outer_date=fold.outer_date, date=fold.outer_date, file=file,
            selected_model=effective_model, graph_enabled=graph_enabled,
            node_count=length(nodes), eligible_edge_count=length(edges),
            denominator=denominator, c1_node_loglik=node_log_c1,
            selected_node_loglik=node_log_selected, common_null_loglik=null_log,
            graph_logz=graph_logz, c1_loglik=c1_loglik, meta_loglik=meta_loglik,
            d_s=d_s, score_sha256=score_sha, status=:PASS, reason=:ok,
        ))
    end
    return (status=:PASS, reason=fallback ? :c2_unavailable_fallback_c1 : :ok,
            effective_model=effective_model, folds=NamedTuple[], nodes=node_rows,
            scans=scan_rows, graph_enabled=graph_enabled, edge_count=edge_count)
end

function _empty_metrics()
    return (coverage=NaN, entropy=NaN, agreement_numerator=0,
            agreement_denominator=0, agreement=NaN, unary_parameter_count=0,
            residualizer_parameter_count=0, null_parameter_count=0,
            conditional_edge_parameter_count=0, unary_meta_total=0,
            graph_meta_total=0, graph_enabled_fold_count=0,
            eligible_edge_count=0, t_obs=NaN, sign_tail_count=0,
            sign_denominator=0, sign_p=NaN, bootstrap_lower=NaN,
            every_date_positive=false)
end

function _terminal_event(status::Symbol, reason::Symbol, outer_date::String="")
    status in (:PASS, :SKIPPED, :FAIL, :BLOCKED) ||
        throw(ArgumentError("invalid terminal event status"))
    return (ordinal=1, outer_date=outer_date, event="terminal", consumed=false,
            terminal_status=status, consequence=reason, reason=reason,
            formula_count=0, bootstrap_count=0, sign_count=0, graph_count=0)
end

function _empty_report(status::Symbol, reason::Symbol;
                       authority::String="", universe::String="", source::String="",
                       input_records::Vector{String}=String[],
                       input::Union{Nothing,EvaluationInput}=nothing,
                       outer_date::String="")::EvaluatorReport
    metrics = _empty_metrics()
    if input !== nothing
        source, provenance = _report_bindings(input)
    else
        provenance = _length_digest(vcat([
            "schema=structured-evaluator-empty-report-provenance-v2",
            authority, universe, String(status), String(reason), source,
        ], input_records))
    end
    provisional = EvaluatorReport(status, reason, NamedTuple[], NamedTuple[], NamedTuple[],
                                  NamedTuple[], NamedTuple[],
                                  NamedTuple[_terminal_event(status, reason, outer_date)], metrics,
                                  authority, universe, source, Float64[], "", provenance)
    result = _result_binding_hash(provenance, _tsv_files(provisional))
    return EvaluatorReport(status, reason, provisional.folds, provisional.nodes,
                           provisional.scans, provisional.dates, provisional.signs,
                           provisional.events, metrics, authority, universe, source,
                           Float64[], result, provenance)
end

function _exhaustive_sign_test(means::Vector{Float64})
    k = length(means)
    1 <= k <= 13 ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence,
              "sign-test date count is outside the supported range")
    all(isfinite, means) ||
        _fail(:FAIL, :numerical_failure, "sign-test means are nonfinite")
    t_obs = sum(means) / k
    isfinite(t_obs) ||
        _fail(:FAIL, :numerical_failure, "sign-test observed statistic is nonfinite")
    denominator = 2^k
    rows = NamedTuple[]
    tail = 0
    for mask in 0:(denominator - 1)
        signs = [((mask >> (index - 1)) & 1) == 1 ? 1 : -1
                 for index in eachindex(means)]
        t_epsilon = sum(signs[index] * means[index]
                        for index in eachindex(means)) / k
        isfinite(t_epsilon) ||
            _fail(:FAIL, :numerical_failure, "sign-test statistic is nonfinite")
        in_upper_tail = t_epsilon >= t_obs
        in_upper_tail && (tail += 1)
        push!(rows, (mask=mask, signs=join(signs, ','), t_epsilon=t_epsilon,
                     in_upper_tail=in_upper_tail))
    end
    denominator > 0 || _fail(:FAIL, :numerical_failure, "sign-test denominator is invalid")
    p = tail / denominator
    isfinite(p) || _fail(:FAIL, :numerical_failure, "sign-test p-value is nonfinite")
    return (t_obs=t_obs, rows=rows, tail=tail, denominator=denominator, p=p)
end

function _evaluate_valid(input::EvaluationInput)::EvaluatorReport
    context = get(_PRODUCTION_CONTEXT, input, nothing)
    context === nothing || _verify_context(context)
    _validate_input(input)
    sorted_folds = sort(copy(input.folds); by=fold -> fold.outer_date)
    fold_rows = NamedTuple[]
    node_rows = NamedTuple[]
    scan_rows = NamedTuple[]
    all_inner_events = NamedTuple[]
    fold_outcomes = NamedTuple[]
    for fold in sorted_folds
        try
            decision = _inner_gate(fold.inner_folds)
            decision_hash = _inner_decision_hash(fold.inner_folds, decision)
            if decision.status != :PASS
                push!(fold_outcomes, (outer_date=fold.outer_date,
                                      status=decision.status, reason=decision.reason))
                continue
            end
            selected = decision.model
            result = _score_fold(fold, selected)
            if result.status != :PASS
                push!(fold_outcomes, (outer_date=fold.outer_date,
                                      status=result.status, reason=result.reason))
                continue
            end
            effective_selected = result.effective_model
            graph = _graph_state(fold, effective_selected)
            graph_enabled = result.graph_enabled
            t12_status = graph.t12_status
            t12_reason = graph.t12_reason
            selected_fit = effective_selected == "C1" ? fold.c1_fit_sha256 : fold.c2_fit_sha256
            t11_status = graph.t11_gate_status
            t11_hash = graph.t11_gate_sha256
            t11_fit = graph.t11_fit_sha256
            fold_status = result.status
            fold_reason = result.reason
            push!(fold_rows, (
                outer_date=fold.outer_date, training_dates=join(fold.training_dates, ','),
                c1_fit=fold.c1_fit_sha256, c2_fit=fold.c2_fit_sha256,
                inner_status=decision.status, inner_reason=decision.reason,
                inner_hash=decision_hash, selected_model=effective_selected,
                selected_fit=selected_fit, t11_gate_status=t11_status,
                t11_gate_sha256=t11_hash, t11_outer_fit_sha256=t11_fit,
                t12_status=t12_status, t12_reason=t12_reason,
                t12_result_sha256=graph.t12_result_sha256, graph_enabled=graph_enabled,
                fold_status=fold_status, fold_reason=fold_reason,
            ))
            append!(node_rows, result.nodes)
            append!(scan_rows, result.scans)
            push!(all_inner_events, (outer_date=fold.outer_date, decision=decision))
            push!(fold_outcomes, (outer_date=fold.outer_date,
                                  status=result.status, reason=result.reason))
        catch error
            error isa EvaluatorError || rethrow()
            push!(fold_outcomes, (outer_date=fold.outer_date,
                                  status=error.status, reason=error.reason))
        end
    end
    nonpass = [outcome for outcome in fold_outcomes if outcome.status != :PASS]
    if !isempty(nonpass)
        winner = first(sort(nonpass; by=outcome ->
            (-_PRECEDENCE[outcome.status], outcome.outer_date, String(outcome.reason))))
        return _empty_report(winner.status, winner.reason;
                             authority=input.authority_sha256,
                             universe=input.universe_sha256,
                             input_records=_input_manifest_records(input), input=input,
                             outer_date=winner.outer_date)
    end
    expected_files = Set(first.(input.scan_dates))
    Set(row.file for row in scan_rows) == expected_files ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence, "scan score coverage differs from universe")
    length(scan_rows) == length(expected_files) ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence, "scan score rows are duplicated or absent")
    by_date = Dict{String,Vector{Float64}}()
    for row in scan_rows
        push!(get!(by_date, row.date, Float64[]), row.d_s)
    end
    dates = sort(collect(keys(by_date)))
    date_rows = NamedTuple[]
    date_means = Float64[]
    for date in dates
        values = by_date[date]
        d = mean(values)
        push!(date_means, d)
        push!(date_rows, (date=date, scan_count=length(values), d_k=d,
                          status=d > 0.0 ? :PASS : :SKIPPED,
                          reason=d > 0.0 ? :ok : :nonpositive_date))
    end
    final_bootstrap = _bootstrap_scan_values(by_date)
    lower = _type7(final_bootstrap)
    every_positive = all(>(0.0), date_means)
    sign_test = _exhaustive_sign_test(date_means)
    t_obs = sign_test.t_obs
    sign_rows = sign_test.rows
    tail = sign_test.tail
    sign_denominator = sign_test.denominator
    sign_p = sign_test.p
    terminal_status = if !every_positive || !(lower > 0.0) || !(sign_p < 0.05)
        :SKIPPED
    else
        :PASS
    end
    terminal_reason = !every_positive ? :nonpositive_date :
                      !(lower > 0.0) ? :bootstrap_not_positive :
                      !(sign_p < 0.05) ? :sign_not_significant : :ok
    selected_nodes = length(node_rows)
    classified = count(row -> row.views_used > 0 && row.p1 != 0.5, node_rows)
    coverage = selected_nodes == 0 ? NaN : classified / selected_nodes
    entropy = selected_nodes == 0 ? NaN : mean(row.entropy for row in node_rows)
    agreement_numerator = sum(row.agreement_numerator for row in node_rows)
    agreement_denominator = sum(row.agreement_denominator for row in node_rows)
    agreement = agreement_denominator == 0 ? NaN : agreement_numerator / agreement_denominator
    graph_folds = count(row -> row.graph_enabled, fold_rows)
    edge_count = sum(row.eligible_edge_count for row in scan_rows)
    metrics = (coverage=coverage, entropy=entropy,
               agreement_numerator=agreement_numerator,
               agreement_denominator=agreement_denominator, agreement=agreement,
               unary_parameter_count=76, residualizer_parameter_count=30,
               null_parameter_count=5, conditional_edge_parameter_count=9,
               unary_meta_total=111, graph_meta_total=120,
               graph_enabled_fold_count=graph_folds, eligible_edge_count=edge_count,
               t_obs=t_obs, sign_tail_count=tail, sign_denominator=sign_denominator,
               sign_p=sign_p, bootstrap_lower=lower, every_date_positive=every_positive)
    events = NamedTuple[]
    ordinal = 0
    for item in all_inner_events
        ordinal += 1
        push!(events, (ordinal=ordinal, outer_date=item.outer_date,
                       event="inner_selection", consumed=true, terminal_status=:NONE,
                       consequence="continue", reason=item.decision.reason,
                       formula_count=11, bootstrap_count=500, sign_count=0,
                       graph_count=0))
    end
    ordinal += 1
    push!(events, (ordinal=ordinal, outer_date="", event="outer_scoring",
                   consumed=true, terminal_status=:NONE, consequence="continue",
                   reason=:fixed_observations, formula_count=4,
                   bootstrap_count=0, sign_count=0, graph_count=graph_folds))
    ordinal += 1
    push!(events, (ordinal=ordinal, outer_date="", event="final_gate",
                   consumed=false, terminal_status=terminal_status,
                   consequence=terminal_status == :PASS ? :complete_all_guards : terminal_reason,
                   reason=terminal_reason, formula_count=5, bootstrap_count=500,
                   sign_count=sign_denominator, graph_count=graph_folds))
    source_sha, provenance = _report_bindings(input)
    provisional = EvaluatorReport(terminal_status, terminal_reason, fold_rows, node_rows,
                                  scan_rows, date_rows, sign_rows, events, metrics,
                                  input.authority_sha256, input.universe_sha256,
                                  source_sha, final_bootstrap, "", provenance)
    result_hash = _result_binding_hash(provenance, _tsv_files(provisional))
    return EvaluatorReport(terminal_status, terminal_reason, fold_rows, node_rows,
                           scan_rows, date_rows, sign_rows, events, metrics,
                           input.authority_sha256, input.universe_sha256,
                           source_sha, final_bootstrap, result_hash, provenance)
end

function evaluate(input::EvaluationInput)::EvaluatorReport
    try
        return _evaluate_valid(input)
    catch error
        error isa EvaluatorError || rethrow()
        return _empty_report(error.status, error.reason;
                             authority=input.authority_sha256,
                             universe=input.universe_sha256,
                             input_records=_input_manifest_records(input), input=input)
    end
end

function _field(value)
    value isa Symbol && return String(value)
    value isa Bool && return _b(value)
    value isa AbstractFloat && return _f(value)
    value isa Integer && return string(value)
    value === nothing && return "NA"
    return String(value)
end

function _table_bytes(header::Vector{String}, rows::Vector{<:AbstractVector})::Vector{UInt8}
    io = IOBuffer()
    println(io, join(header, '\t'))
    for row in rows
        length(row) == length(header) || throw(ArgumentError("report row width differs"))
        fields = [_field(value) for value in row]
        any(occursin(r"[\t\r\n]", value) for value in fields) &&
            throw(ArgumentError("report field contains a delimiter"))
        println(io, join(fields, '\t'))
    end
    return take!(io)
end

function _fold_bytes(report::EvaluatorReport)
    rows = Vector{Vector{Any}}()
    for row in report.folds
        push!(rows, Any[
            1, row.outer_date, row.training_dates, row.c1_fit, row.c2_fit,
            row.inner_status, row.inner_reason, row.inner_hash, row.selected_model,
            row.selected_fit, row.t11_gate_status, row.t11_gate_sha256,
            row.t11_outer_fit_sha256, row.t12_status, row.t12_reason,
            row.t12_result_sha256, row.graph_enabled, row.fold_status, row.fold_reason,
        ])
    end
    return _table_bytes(_FOLD_HEADER, rows)
end

function _node_bytes(report::EvaluatorReport)
    rows = Vector{Vector{Any}}()
    for row in report.nodes
        push!(rows, Any[
            1, row.outer_date, row.date, row.file, row.lobe, row.t_nm,
            row.selected_model, row.selected_fit, row.graph_enabled, row.views_used,
            row.c1_u, row.selected_u, row.p1, row.hard, row.entropy,
            row.agreement_numerator, row.agreement_denominator, row.status, row.reason,
        ])
    end
    return _table_bytes(_NODE_HEADER, rows)
end

function _scan_bytes(report::EvaluatorReport)
    rows = Vector{Vector{Any}}()
    for row in report.scans
        push!(rows, Any[
            1, row.outer_date, row.date, row.file, row.selected_model,
            row.graph_enabled, row.node_count, row.eligible_edge_count, row.denominator,
            row.c1_node_loglik, row.selected_node_loglik, row.common_null_loglik,
            row.graph_logz, row.c1_loglik, row.meta_loglik, row.d_s,
            row.score_sha256, row.status, row.reason,
        ])
    end
    return _table_bytes(_SCAN_HEADER, rows)
end

function _date_bytes(report::EvaluatorReport)
    rows = [Any[1, row.date, row.scan_count, row.d_k, row.status, row.reason]
            for row in report.dates]
    return _table_bytes(_DATE_HEADER, rows)
end

function _bootstrap_bytes(report::EvaluatorReport)
    values = length(report.bootstrap_values) == length(_BOOTSTRAP_SEEDS) ?
             report.bootstrap_values : fill(NaN, length(_BOOTSTRAP_SEEDS))
    rows = [Any[1, seed, values[seed + 1]] for seed in _BOOTSTRAP_SEEDS]
    return _table_bytes(_BOOTSTRAP_HEADER, rows)
end

function _sign_bytes(report::EvaluatorReport)
    rows = [Any[1, row.mask, row.signs, row.t_epsilon, row.in_upper_tail]
            for row in report.signs]
    return _table_bytes(_SIGN_HEADER, rows)
end

function _event_bytes(report::EvaluatorReport)
    rows = [Any[
        1, row.ordinal, row.outer_date, row.event, row.consumed,
        row.terminal_status, row.consequence, row.reason, row.formula_count,
        row.bootstrap_count, row.sign_count, row.graph_count,
    ] for row in report.events]
    return _table_bytes(_EVENT_HEADER, rows)
end

function _summary_bytes(report::EvaluatorReport)
    metrics = report.metrics
    status = report.status
    reason = report.reason
    row = Any[
        1, status, reason, length(report.dates), length(report.scans), length(report.nodes),
        metrics.eligible_edge_count, metrics.t_obs, metrics.sign_tail_count,
        metrics.sign_denominator, metrics.sign_p, metrics.bootstrap_lower,
        metrics.every_date_positive, metrics.coverage, metrics.entropy,
        metrics.agreement_numerator, metrics.agreement_denominator, metrics.agreement,
        metrics.unary_parameter_count, metrics.residualizer_parameter_count,
        metrics.null_parameter_count, metrics.conditional_edge_parameter_count,
        metrics.unary_meta_total, metrics.graph_meta_total,
        metrics.graph_enabled_fold_count,
    ]
    return _table_bytes(_SUMMARY_HEADER, [row])
end

function _tsv_files(report::EvaluatorReport)::Dict{String,Vector{UInt8}}
    files = Dict{String,Vector{UInt8}}(
        "bootstrap.tsv" => _bootstrap_bytes(report),
        "dates.tsv" => _date_bytes(report),
        "events.tsv" => _event_bytes(report),
        "folds.tsv" => _fold_bytes(report),
        "nodes.tsv" => _node_bytes(report),
        "scans.tsv" => _scan_bytes(report),
        "signs.tsv" => _sign_bytes(report),
        "summary.tsv" => _summary_bytes(report),
    )
    Set(keys(files)) == setdiff(Set(_FILES), Set(["receipt.toml"])) ||
        throw(ArgumentError("incomplete report file set"))
    return files
end

function _canonical_tsv_hash(tsv::Dict{String,Vector{UInt8}})::String
    io = IOBuffer()
    for name in sort(collect(keys(tsv)))
        name_bytes = Vector{UInt8}(codeunits(name))
        write(io, UInt64(length(name_bytes))); write(io, name_bytes)
        write(io, UInt64(length(tsv[name]))); write(io, tsv[name])
    end
    return _hash_bytes(take!(io))
end

function _result_binding_hash(provenance::String,
                              tsv::Dict{String,Vector{UInt8}})::String
    tsv_hash = _canonical_tsv_hash(tsv)
    return _length_digest([
        "schema=structured-evaluator-result-binding-v1",
        "provenance_sha256=$provenance",
        "result_tsv_sha256=$tsv_hash",
    ])
end

function _legacy_receipt_bytes(report::EvaluatorReport,
                               tsv::Dict{String,Vector{UInt8}})
    io = IOBuffer()
    println(io, "schema = ", repr(_SCHEMA * "_receipt_v2"))
    println(io, "schema_version = 2")
    println(io, "status = ", repr(String(report.status)))
    println(io, "reason = ", repr(String(report.reason)))
    println(io, "authority_sha256 = ", repr(report.authority_sha256))
    println(io, "authority_manifest_sha256 = ", repr(report.authority_sha256))
    println(io, "universe_sha256 = ", repr(report.universe_sha256))
    println(io, "input_universe_sha256 = ", repr(report.universe_sha256))
    println(io, "evaluator_source_sha256 = ", repr(report.source_sha256))
    println(io, "evaluator_config_sha256 = ", repr(_report_config_sha256(report)))
    println(io, "preclosure_plan_sha256 = ", repr(_PLAN_SHA256))
    println(io, "gateclosure_sha256 = ", repr(_GATECLOSURE_SHA256))
    println(io, "preclosure_boulder_sha256 = ", repr(_BOULDER_SHA256))
    println(io, "live_plan_boulder_runtime_authority = false")
    println(io, "review_sha256 = ", repr(_REVIEW_SHA256))
    println(io, "closure_v1_claim_sha256 = ", repr(_CLOSURE_V1_CLAIM_SHA256))
    println(io, "closure_v1_review_sha256 = ", repr(_CLOSURE_V1_REVIEW_SHA256))
    println(io, "closure_v1_root_manifest_sha256 = ", repr(_CLOSURE_V1_ROOT_MANIFEST_SHA256))
    println(io, "closure_v1_publication_receipt_sha256 = ",
            repr(_CLOSURE_V1_PUBLICATION_RECEIPT_SHA256))
    println(io, "closure_v1_publication_receipt_sidecar_sha256 = ",
            repr(_CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_SHA256))
    println(io, "result_sha256 = ", repr(report.result_sha256))
    println(io, "provenance_sha256 = ", repr(report.provenance_sha256))
    println(io, "provenance_manifest_sha256 = ", repr(report.provenance_sha256))
    println(io, "result_tsv_sha256 = ", repr(_canonical_tsv_hash(tsv)))
    println(io, "file_count = 9")
    println(io, "serialization_encoding = \"UTF-8\"")
    println(io, "serialization_line_ending = \"LF\"")
    println(io, "serialization_final_lf = true")
    println(io, "serialization_float_format = \"%.17g\"")
    println(io, "serialization_missing_token = \"NA\"")
    println(io, "serialization_boolean_format = \"lowercase\"")
    println(io, "files = ", repr(collect(_FILES)))
    for name in sort(collect(keys(tsv)))
        println(io)
        println(io, "[artifacts.", repr(name), "]")
        println(io, "sha256 = ", repr(_hash_bytes(tsv[name])))
        println(io, "bytes = ", length(tsv[name]))
        println(io, "data_rows = ", max(0, count(==(UInt8('\n')), tsv[name]) - 1))
    end
    return take!(io)
end

# True successor activated report receipt.  It names the activation, spec,
# review, source and Boulder/publication bindings and deliberately omits the
# historical GateClosure/closure-v1/plan authority fields, which do not certify
# an activated run.  The earlier transition `publication` record is a distinct
# upstream artifact; this receipt cannot certify its own publication.
function _activated_receipt_bytes(report::EvaluatorReport,
                                  tsv::Dict{String,Vector{UInt8}},
                                  activation)
    io = IOBuffer()
    println(io, "schema = ", repr(_T13_ACTIVATED_REPORT_SCHEMA))
    println(io, "schema_version = 2")
    println(io, "status = ", repr(String(report.status)))
    println(io, "reason = ", repr(String(report.reason)))
    println(io, "authority_sha256 = ", repr(report.authority_sha256))
    println(io, "universe_sha256 = ", repr(report.universe_sha256))
    println(io, "evaluator_source_sha256 = ", repr(report.source_sha256))
    println(io, "evaluator_config_sha256 = ", repr(_report_config_sha256(report)))
    _t13_activation_receipt_lines(io, activation)
    println(io, "result_sha256 = ", repr(report.result_sha256))
    println(io, "provenance_sha256 = ", repr(report.provenance_sha256))
    println(io, "result_tsv_sha256 = ", repr(_canonical_tsv_hash(tsv)))
    println(io, "file_count = 9")
    println(io, "serialization_encoding = \"UTF-8\"")
    println(io, "serialization_line_ending = \"LF\"")
    println(io, "serialization_final_lf = true")
    println(io, "serialization_float_format = \"%.17g\"")
    println(io, "serialization_missing_token = \"NA\"")
    println(io, "serialization_boolean_format = \"lowercase\"")
    println(io, "files = ", repr(collect(_FILES)))
    for name in sort(collect(keys(tsv)))
        println(io)
        println(io, "[artifacts.", repr(name), "]")
        println(io, "sha256 = ", repr(_hash_bytes(tsv[name])))
        println(io, "bytes = ", length(tsv[name]))
        println(io, "data_rows = ", max(0, count(==(UInt8('\n')), tsv[name]) - 1))
    end
    return take!(io)
end

function _receipt_bytes(report::EvaluatorReport, tsv::Dict{String,Vector{UInt8}};
                        activation=nothing)
    activation === nothing && return _legacy_receipt_bytes(report, tsv)
    return _activated_receipt_bytes(report, tsv, activation)
end

function report_files(report::EvaluatorReport;
                      activation=nothing)::Dict{String,Vector{UInt8}}
    tsv = _tsv_files(report)
    files = copy(tsv)
    files["receipt.toml"] = _receipt_bytes(report, tsv; activation=activation)
    sort(collect(keys(files))) == sort(collect(_FILES)) ||
        throw(ArgumentError("report file set differs from the frozen contract"))
    return Dict(name => files[name] for name in _FILES)
end

struct _Snapshot
    path::String
    bytes::Vector{UInt8}
    sha256::String
    device::UInt64
    inode::UInt64
    nlink::UInt64
end

# Dedicated Main-control object for the owner-selected "Main Boulder only"
# model.  It is deliberately separate from the ordinary R-relative snapshot
# vector: the execution root R owns code/inputs/outputs, while the single live
# `.omo/boulder.json` lives under the canonical Main root M.  The object holds
# the one fixed control-file snapshot plus the approved modes and the captured
# M / M/.omo directory identities.  It is an integrity/scope object, not a
# signer or process attestation; the externally reviewed launch pin remains the
# trust anchor for which M was reviewed.
struct _MainControl
    root::String
    file::_Snapshot
    root_mode::UInt
    parent_mode::UInt
    file_mode::UInt
    root_identity::Tuple{UInt64,UInt64}
    parent_identity::Tuple{UInt64,UInt64}
    preimage_sha256::String
    postimage_sha256::String
end

# Activation administrative bindings.  These are integrity/scope identities for
# an externally reviewed activation launch; they are not a signer identity and
# they never carry a digest of the receipt itself (that would be a cycle).
struct _ActivationBindings
    scope::String
    receipt_path::String
    receipt_sha256::String
    spec_path::String
    spec_sha256::String
    source_review_path::String
    source_review_sha256::String
    transition_review_path::String
    transition_review_sha256::String
    publication_path::String
    publication_sha256::String
    boulder_preimage_path::String
    boulder_preimage_sha256::String
    boulder_postimage_path::String
    boulder_postimage_sha256::String
    entrypoint_path::String
    entrypoint_sha256::String
    destination::String
    modes::Dict{String,UInt}
    root::String
    command_body::Vector{String}
    directory_identities::Dict{String,Tuple{UInt64,UInt64}}
    control::_MainControl
end

struct _ProductionContext
    root::String
    snapshots::Vector{_Snapshot}
    inventories::Vector{Tuple{String,Vector{String}}}
    manifest_sha256::String
    activation::Union{Nothing,_ActivationBindings}
end

_ProductionContext(root::String, snapshots::Vector{_Snapshot},
                   inventories::Vector{Tuple{String,Vector{String}}},
                   manifest_sha256::String) =
    _ProductionContext(root, snapshots, inventories, manifest_sha256, nothing)

struct _ActivationAuthorization
    root::String
    bindings::_ActivationBindings
    snapshots::Vector{_Snapshot}
    inventories::Vector{Tuple{String,Vector{String}}}
    spec_files::Vector{NamedTuple}
    spec_directories::Vector{NamedTuple}
end

const _PRODUCTION_CONTEXT = IdDict{Any,_ProductionContext}()
const _PRODUCTION_SELECTION_METADATA = IdDict{Any,Any}()

function _report_config_sha256(report::EvaluatorReport)::String
    for (input, context) in _PRODUCTION_CONTEXT
        input.authority_sha256 == report.authority_sha256 || continue
        runtime_v4 = joinpath(context.root, _RUNTIME_V4_CONFIG_PATH)
        any(snapshot.path == runtime_v4 for snapshot in context.snapshots) &&
            return _RUNTIME_V4_CONFIG_SHA256
        runtime_v3 = joinpath(context.root, _RUNTIME_V3_CONFIG_PATH)
        any(snapshot.path == runtime_v3 for snapshot in context.snapshots) &&
            return _RUNTIME_V3_CONFIG_SHA256
        runtime_v2 = joinpath(context.root, _RUNTIME_V2_CONFIG_PATH)
        any(snapshot.path == runtime_v2 for snapshot in context.snapshots) &&
            return _RUNTIME_V2_CONFIG_SHA256
        runtime_v1 = joinpath(context.root, _RUNTIME_V1_CONFIG_PATH)
        any(snapshot.path == runtime_v1 for snapshot in context.snapshots) &&
            return _RUNTIME_V1_CONFIG_SHA256
    end
    return _CONFIG_SHA256
end

function _length_digest(records::Vector{String})::String
    io = IOBuffer()
    for record in records
        bytes = Vector{UInt8}(codeunits(record))
        write(io, UInt64(length(bytes)))
        write(io, bytes)
    end
    return _hash_bytes(take!(io))
end

const _CONFIG_SHA256 =
    "ec0546096b3c4742cd86d8c3d40788a5894d20fc5a4702b0318581f10f8b0b90"
const _HISTORICAL_CONFIG_PATH = "config/unit_assignment_structured_evaluator.toml"
const _RUNTIME_V1_CONFIG_PATH = "config/unit_assignment_structured_evaluator_runtime.toml"
const _RUNTIME_V1_CONFIG_SHA256 =
    "4a355c7ec5572976f5f54eaf14f2a01e8c91a9304d789e60b844388ac20939a4"
const _RUNTIME_V2_CONFIG_PATH = "config/unit_assignment_structured_evaluator_runtime_v2.toml"
const _RUNTIME_V2_CONFIG_SHA256 =
    "69da829faf71a00ffec07e976c5c95245bbb57fd1ac7d0d98e2d54d8064175ae"
const _RUNTIME_V3_CONFIG_PATH =
    "config/unit_assignment_structured_evaluator_runtime_v3.toml"
const _RUNTIME_V3_CONFIG_SHA256 =
    "a0a04794b346f351c61a281384869bcde3f378aa0836c9271197d4004488586e"
const _RUNTIME_V4_CONFIG_PATH =
    "config/unit_assignment_structured_evaluator_runtime_v4.toml"
const _RUNTIME_V4_CONFIG_SHA256 =
    "abe22ed5047c898f594067d4a54cbe8fcc99c15fef17b429fa366f64c255547b"
const _RUNTIME_CONFIG_PATH = _RUNTIME_V4_CONFIG_PATH
const _RUNTIME_CONFIG_SHA256 = _RUNTIME_V4_CONFIG_SHA256
const _RUNTIME_V4_JULIA = v"1.12.7"
const _RUNTIME_V4_EXECUTABLE_SHA256 =
    "582718c20c563d824b88a08e2a9d0afcf5015b3053f425697cecc79ec4e9ea36"
const _RUNTIME_V4_SYSIMAGE_SHA256 =
    "03041fa63edac21123a7d4ced06364243c5b9653f4ceb6ba618f69876eefa66f"
const _RUNTIME_V4_BOULDER_SHA256 =
    "e509c919734c8110ab7ad4a7e2d475a552a0f2fca5f3533c9863683be1ac9133"
const _RUNTIME_V4_AUTHORITY_SCHEMA = "schema=structured-evaluator-authority-v6"

# Runtime-v4 is an administrative Gate5 candidate.  Keep its complete
# authority binding in this source so a TOML table cannot silently downgrade a
# path/hash/mode check to a type-only check.  The config byte digest is checked
# before TOML parsing; these values are the second, semantic fail-closed
# boundary used by the audit and synthetic contract seams.
const _RUNTIME_V4_AUTHORITY_BINDINGS = (
    ("t11_source_tar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/v7-run/input/source.tar",
     "t11_source_tar_sha256", "aaf515357901c3276edf867c04ecfd0229b494f836ddd72ef3f0effed9f367cc",
     "t11_source_tar_mode", "0444"),
    ("t11_source_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/v7-run/input/source-manifest.tsv",
     "t11_source_manifest_sha256", "2e7696824e208ba0c3ad960a08e452f81d79f892c776a1ebd0b28cbf1a196714",
     "t11_source_manifest_mode", "0444"),
    ("t11_source_symlink_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/v7-run/input/source-symlink-manifest.tsv",
     "t11_source_symlink_manifest_sha256", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
     "t11_source_symlink_manifest_mode", "0444"),
    ("t11_test_path", "test/test_structured_edge_admission.jl",
     "t11_test_sha256", "bff69a629ef99056bdf888478b15633abc921c29cd268d5625c7de9f8ed97e3c",
     "t11_test_mode", "0644"),
    ("t11_cli_path", "test/evaluate_structured_edge_admission.jl",
     "t11_cli_sha256", "77c7b2d5f1b9f73ae5c3bdfd310cbbad75771f5a38c5243488f4c2cc6db7ce27",
     "t11_cli_mode", "0644"),
    ("t11_claim_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/T11RuntimeAuthority.json",
     "t11_claim_sha256", "36ae27e177025c39cf7c41af0bf6ac4f27e43b263c04628de53941a07b6e937a",
     "t11_claim_mode", "0444"),
    ("t11_claim_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/T11RuntimeAuthority.sha256",
     "t11_claim_sidecar_sha256", "7e30c20c92dddcf621565fab73585befcc59856479a6e6d0230c6403106c9a99",
     "t11_claim_sidecar_mode", "0444"),
    ("t11_review_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/review/AdversarialVerify.json",
     "t11_review_sha256", "dbf3f11801bc6ec4e06afb73d516523c1c88c577493891e81c0fbee9f660e537",
     "t11_review_mode", "0644"),
    ("t11_review_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/review/AdversarialVerify.sha256",
     "t11_review_sidecar_sha256", "eed9d31c7585701ca99d6ab4416368e97f6128f7420b0023e5945df67f074c0c",
     "t11_review_sidecar_mode", "0644"),
    ("t11_publication_receipt_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/PublicationReceipt.json",
     "t11_publication_receipt_sha256", "38b4546324843d478b6d14e0fdeaa49ced5aa0d9c81011f9bf5ac864a34e1f26",
     "t11_publication_receipt_mode", "0444"),
    ("t11_publication_receipt_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/PublicationReceipt.sha256",
     "t11_publication_receipt_sidecar_sha256", "b3bbcdacb0dcb929807a48e8ca0c079005761845db7a38fdcc93ce9871384d9f",
     "t11_publication_receipt_sidecar_mode", "0444"),
    ("t11_payload_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/payload-files.sha256",
     "t11_payload_manifest_sha256", "6aa2aad9c19b0fded88d0717c9321e13f2d5d3b10cde871246105a817744a36f",
     "t11_payload_manifest_mode", "0444"),
    ("t11_payload_manifest_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/payload-files.sha256.sha256",
     "t11_payload_manifest_sidecar_sha256", "232a33d0bf4ee8ad5fda05f738e85de316e05933f5c3a82a834ca271b038130a",
     "t11_payload_manifest_sidecar_mode", "0444"),
    ("t12_source_tar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/v9-input/source.tar",
     "t12_source_tar_sha256", "e3538df88145753fc4eb65b1b5a0032e01391ef35f2ce54f16a1aede837b52b6",
     "t12_source_tar_mode", "0400"),
    ("t12_source_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/v9-input/source-manifest.tsv",
     "t12_source_manifest_sha256", "ab9deb0c68f349a1908d9dba81a5cc5187d85332846ba9f0e9855b88bff9a3b2",
     "t12_source_manifest_mode", "0400"),
    ("t12_source_symlink_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/v9-input/source-symlink-manifest.tsv",
     "t12_source_symlink_manifest_sha256", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
     "t12_source_symlink_manifest_mode", "0400"),
    ("t12_edge_model_path", "test/lib/structured_assignment/edge_model.jl",
     "t12_edge_model_sha256", "216571a5c74ba45293ab05b24c3363b05cb77fba6bcf3c690adcb57b42e3a146",
     "t12_edge_model_mode", "0644"),
    ("t12_chain_path", "test/lib/structured_assignment/chain_inference.jl",
     "t12_chain_sha256", "4c5b84339e61a1bc6dadba12373b8b98745f30a3c4359a5c8e92a70e24d96598",
     "t12_chain_mode", "0644"),
    ("t12_test_path", "test/test_structured_chain_inference.jl",
     "t12_test_sha256", "8e0a2d6a4e965bb0e00f290299264f0a75ead553c6d04a95f87ad7815205daeb",
     "t12_test_mode", "0644"),
    ("t12_claim_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/claim/DoneClaim.json",
     "t12_claim_sha256", "ab0ca8cfe5651998d3e002d3a158c58a113b5f00f5eb3feb88d6091d8d1b984f",
     "t12_claim_mode", "0400"),
    ("t12_claim_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/claim/DoneClaim.json.sha256",
     "t12_claim_sidecar_sha256", "84a07bcc25ae0106209652e1f1dafa5f7e057981f6d98b277d821a9d3be7820f",
     "t12_claim_sidecar_mode", "0400"),
    ("t12_review_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/review/AdversarialVerify.json",
     "t12_review_sha256", "28242978aa521a02c0783441605b09e733019b64ef6f5a5c55334c4a63df7946",
     "t12_review_mode", "0400"),
    ("t12_review_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/review/AdversarialVerify.sha256",
     "t12_review_sidecar_sha256", "3b5dd69d17bb930ea5a3d8247fd1b8a68bf7ceff045ef56be8138ca09e5b13a1",
     "t12_review_sidecar_mode", "0400"),
    ("t12_publication_receipt_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/PublicationReceipt.json",
     "t12_publication_receipt_sha256", "47d011da68b8029d38b3a55a11d6888dcd9fd05d82f99a7b721e8f42b51b6980",
     "t12_publication_receipt_mode", "0400"),
    ("t12_publication_receipt_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/PublicationReceipt.sha256",
     "t12_publication_receipt_sidecar_sha256", "536c00588a44381c426cf464ceb69fad9cbf0f187055c6260566f530787e8cd1",
     "t12_publication_receipt_sidecar_mode", "0400"),
    ("t12_payload_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/payload-manifest.tsv",
     "t12_payload_manifest_sha256", "4775cde16d87bb71958016e98947f9e72fc2ee9ac8926f56f3d226314f331b96",
     "t12_payload_manifest_mode", "0400"),
    ("t12_payload_manifest_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/payload-manifest.tsv.sha256",
     "t12_payload_manifest_sidecar_sha256", "450038587e0a058cd55a11948e7e5938059b797243306e222d7588f3fbbdd7a0",
     "t12_payload_manifest_sidecar_mode", "0400"),
    ("t12_evidence_manifest_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence-manifest.tsv",
     "t12_evidence_manifest_sha256", "7f48993887f70733d0330962052383513ef7a01c41a3559808570b993bfc42b8",
     "t12_evidence_manifest_mode", "0400"),
    ("t12_evidence_manifest_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence-manifest.tsv.sha256",
     "t12_evidence_manifest_sidecar_sha256", "3a4f8da410fd1daf3af8f6118c88fbdf2f040d5b4a0f832f98dbe67fdc27359f",
     "t12_evidence_manifest_sidecar_mode", "0400"),
    ("t12_runtime_waiver_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/runtime-waiver/RuntimeComparatorWaiver.json",
     "t12_runtime_waiver_sha256", "16cd577ac099ee2561c279e980bcd33ad8fd07ae1d5820bab84862efdcc5eec2",
     "t12_runtime_waiver_mode", "0400"),
    ("t12_runtime_waiver_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/runtime-waiver/RuntimeComparatorWaiver.sha256",
     "t12_runtime_waiver_sidecar_sha256", "08a67f2340d073f485d5b69068e336be07d247299f44cd573ef9ec034150e143",
     "t12_runtime_waiver_sidecar_mode", "0400"),
    ("gate4_expected_boulder_postimage_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-authority-rebind-v3/transition/BoulderExpectedPostimage.json",
     "gate4_expected_boulder_postimage_sha256", _RUNTIME_V4_BOULDER_SHA256,
     "gate4_expected_boulder_postimage_mode", "0400"),
    ("gate4_expected_boulder_postimage_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-authority-rebind-v3/transition/BoulderExpectedPostimage.json.sha256",
     "gate4_expected_boulder_postimage_sidecar_sha256", "d038fbdd7adad660c3a9089ce9b13de39a2173788d0db4ad01c2d467c672a4d0",
     "gate4_expected_boulder_postimage_sidecar_mode", "0400"),
    ("gate4_transition_receipt_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-boulder-transition-v1/TransitionReceipt.json",
     "gate4_transition_receipt_sha256", "87818b63b132d8856dec22cb389cc936283bb0827cfa4884cbfd6e8c169de092",
     "gate4_transition_receipt_mode", "0400"),
    ("gate4_transition_receipt_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-boulder-transition-v1/TransitionReceipt.json.sha256",
     "gate4_transition_receipt_sidecar_sha256", "5cc4e9737ff71e21d95222fc5a60ad4ef4ace43a24e808c9b3c1c36c7b44de33",
     "gate4_transition_receipt_sidecar_mode", "0400"),
    ("gate4_post_review_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-boulder-transition-v1-receipt-recovery-v2/post-review/IndependentPostExecutionReview.json",
     "gate4_post_review_sha256", "80a133ba8b13d8ac29b1f8d9dcbb07d5c10ee19d47dc60e5105d0897b8032e13",
     "gate4_post_review_mode", "0400"),
    ("gate4_post_review_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-boulder-transition-v1-receipt-recovery-v2/post-review/IndependentPostExecutionReview.json.sha256",
     "gate4_post_review_sidecar_sha256", "bc2acaee30213ed2eb8ab4b34600d2b8beb3287443b81f80d223452d60b9d381",
     "gate4_post_review_sidecar_mode", "0400"),
    ("gate4_oracle_final_review_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-boulder-transition-v1-receipt-recovery-v2/post-review/OracleFinalVerify.json",
     "gate4_oracle_final_review_sha256", "cab6c866c256ce160c13b4694531f661b2db9cbfac5dd47b34cf69b073635fbd",
     "gate4_oracle_final_review_mode", "0400"),
    ("gate4_oracle_final_review_sidecar_path", ".omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate4-boulder-transition-v1-receipt-recovery-v2/post-review/OracleFinalVerify.json.sha256",
     "gate4_oracle_final_review_sidecar_sha256", "06381b69be68e4abec1aad5f324632ef6b2345c3c8f5f69cca5226e85d2b6751",
     "gate4_oracle_final_review_sidecar_mode", "0400"),
)
const _RUNTIME_V4_AUTHORITY_NONBINDING_KEYS = (
    "proposal_path", "proposal_sha256", "plan_path", "plan_sha256",
    "candidate_path", "candidate_sha256", "model_path", "model_sha256",
    "t8_facade_path", "t8_facade_sha256", "t8_source_bundle_path",
    "t8_source_bundle_sha256", "t8_claim_path", "t8_claim_sha256",
    "t8_review_path", "t8_review_sha256",
)
const _RUNTIME_T12_EDGE_MODEL_SHA256 =
    "13c0ea81ccc44ce1154290bd9d7997bf8ddedba8fed262a4e7b2f949835435c5"
const _RUNTIME_T12_ROOT_MANIFEST_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t10-source-bundle-authority-correction-v1/root-manifest.sha256"
const _RUNTIME_T12_ROOT_MANIFEST_SHA256 =
    "bc2bdd2c9e8b618eb90b8f3b86af96a1d009dcc56eca43cae4cfed9b381b20f8"
const _RUNTIME_T12_CLAIM_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t10-source-bundle-authority-correction-v1/DoneClaim.json"
const _RUNTIME_T12_CLAIM_SHA256 =
    "35adeb66c9a5ac0202bfa8c3c49fe6c89cb6e6915a0bebaeae33f6c64faba35c"
const _RUNTIME_T12_REVIEW_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t10-source-bundle-authority-correction-v1/review/AdversarialVerify.json"
const _RUNTIME_T12_REVIEW_SHA256 =
    "f1babbacf371f835c3d766d5c517a13a9c546a07ae37209a35191bdf9348a504"
const _RUNTIME_T12_PUBLICATION_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t10-source-bundle-authority-correction-v1-publication-receipt.json"
const _RUNTIME_T12_PUBLICATION_SHA256 =
    "a58b7240dd56c277ff786261ad5abf73815dc8b46cca89716a263e89a6756c7c"
const _RUNTIME_T12_CANONICAL_TODO10_BUNDLE_SHA256 =
    "285809e3711706d7a059ba2a0dd4139b44c1874dec79308e2762a63799aa954e"
const _RUNTIME_T12_CORRECTED_TEST_SHA256 =
    "fb909c2f969910856664ac247bdab1b92ddad5c785856bd27d3ca305bf7f4bbf"
const _RUNTIME_V2_T12_EDGE_MODEL_SHA256 =
    "5e3b1c0371bb0387024714d420d25a96a9d076aa3f0ba6637a10899837a6c8b1"
const _RUNTIME_V2_T12_ROOT_MANIFEST_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t11-input-role-authority-correction-v2/root-manifest.sha256"
const _RUNTIME_V2_T12_ROOT_MANIFEST_SHA256 =
    "d079f2cee7e0558a93b1bfda156c979e86deb988beb3ef0cb553759fe4bb9626"
const _RUNTIME_V2_T12_CLAIM_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t11-input-role-authority-correction-v2/DoneClaim.json"
const _RUNTIME_V2_T12_CLAIM_SHA256 =
    "0b1a26903851f574f00991ef6e8e8b89e100712a07b59d8346e637161f4c3583"
const _RUNTIME_V2_T12_REVIEW_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t11-input-role-authority-correction-v2/review/AdversarialVerify.json"
const _RUNTIME_V2_T12_REVIEW_SHA256 =
    "d4f984cb070757da30ac713928c41473164911f4b1d95aa48ff8fae9534084fb"
const _RUNTIME_V2_T12_PUBLICATION_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "t11-input-role-authority-correction-v2-publication-receipt.json"
const _RUNTIME_V2_T12_PUBLICATION_SHA256 =
    "7e16b95824687731f160d0e3f13aa1b7c5b460e60989900df0f2093582286c9e"
const _RUNTIME_V2_T12_CORRECTED_TEST_SHA256 =
    "b8513fbbdd4c217dd2715527b40b40a736083e6a60d118cad2a2a0d8e669b9c9"
const _RUNTIME_V3_T12_EDGE_MODEL_SHA256 =
    "216571a5c74ba45293ab05b24c3363b05cb77fba6bcf3c690adcb57b42e3a146"
const _RUNTIME_V3_T12_CHAIN_SHA256 =
    "4c5b84339e61a1bc6dadba12373b8b98745f30a3c4359a5c8e92a70e24d96598"
const _RUNTIME_V3_T12_CORRECTED_TEST_SHA256 =
    "8e0a2d6a4e965bb0e00f290299264f0a75ead553c6d04a95f87ad7815205daeb"
const _RUNTIME_V3_T12_ROOT =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "stable-binary-marginal-normalization-v3"
const _RUNTIME_V3_T12_ROOT_MANIFEST_PATH = _RUNTIME_V3_T12_ROOT * "/root-manifest.sha256"
const _RUNTIME_V3_T12_ROOT_MANIFEST_SHA256 =
    "45f846e35f2965b50d13f4242add20a7d3e4e25e9ada94e619fa12003ed3b03f"
const _RUNTIME_V3_T12_CLAIM_PATH = _RUNTIME_V3_T12_ROOT * "/DoneClaim.json"
const _RUNTIME_V3_T12_CLAIM_SHA256 =
    "af8f0cecc26613856a479cc0db661934de34329158428e92bebcf024579148f4"
const _RUNTIME_V3_T12_REVIEW_PATH = _RUNTIME_V3_T12_ROOT * "/review/AdversarialVerify.json"
const _RUNTIME_V3_T12_REVIEW_SHA256 =
    "f14d2099dc62bc9fa19dcc88b977e53e99bf1090a2b346fabd8c6ed138149efc"
const _RUNTIME_V3_T12_PUBLICATION_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t12/" *
    "stable-binary-marginal-normalization-v3-publication-receipt.json"
const _RUNTIME_V3_T12_PUBLICATION_SHA256 =
    "bca7050c0985fa6ce870a3d660ec3244333d584d22d28a0e450a0fe613a3a118"
const _RUNTIME_V3_T12_EDGE_PRODUCT_PATH = _RUNTIME_V3_T12_ROOT * "/products/v3/edge_model.jl"
const _RUNTIME_V3_T12_CHAIN_PRODUCT_PATH = _RUNTIME_V3_T12_ROOT * "/products/v3/chain_inference.jl"
const _RUNTIME_V3_T12_TEST_PRODUCT_PATH =
    _RUNTIME_V3_T12_ROOT * "/products/v3/test_structured_chain_inference.jl"
const _RUNTIME_V3_T12_ROOT_MANIFEST_ENTRIES = 92
const _RUNTIME_V3_T12_ROOT_FILE_COUNT = 93
const _RUNTIME_V3_T12_ROOT_DIRECTORY_COUNT = 32
const _RUNTIME_V3_T12_ROOT_BYTES = 62981603
const _RUNTIME_V3_AUTHORITY_SCHEMA = "schema=structured-evaluator-authority-v5"
const _PLAN_SHA256 =
    "a9b386d613829e8f7e20b6e33f8e80898fa9a55f0b344dbb92ea64ac6804f3d0"
const _BOULDER_SHA256 =
    "10d685cac291c63cb12bc70b61b01ab5ed51873f2d7130abf1d377ef1fbb40b0"
const _GATECLOSURE_SHA256 =
    "e9c613c6b6609df900e61c02339742e6bb75053d8f977767392adb42a9795895"
const _REVIEW_SHA256 =
    "f42412c8f84ee3f1c4003d19962a6442bde2fb170dc0cc53a7fb1336588197b8"

# These are the immutable implementation-closure authorities.  The Plan and
# Boulder hashes above are used only as explicit preclosure provenance;
# neither live file is a runtime authority snapshot after Todo13 closure.
const _CLOSURE_V1_ROOT =
    ".omo/evidence/structured-label-free-unit-assignment/t13/implementation/closure-v1"
const _CLOSURE_V1_CLAIM_PATH = _CLOSURE_V1_ROOT * "/DoneClaim.json"
const _CLOSURE_V1_CLAIM_SIDECAR_PATH = _CLOSURE_V1_ROOT * "/DoneClaim.sha256"
const _CLOSURE_V1_REVIEW_PATH = _CLOSURE_V1_ROOT * "/review/AdversarialVerify.json"
const _CLOSURE_V1_REVIEW_SIDECAR_PATH = _CLOSURE_V1_ROOT * "/review/AdversarialVerify.sha256"
const _CLOSURE_V1_ROOT_MANIFEST_PATH = _CLOSURE_V1_ROOT * "/root-manifest.sha256"
const _CLOSURE_V1_PUBLICATION_RECEIPT_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t13/implementation/" *
    "closure-v1-publication-receipt.json"
const _CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_PATH =
    ".omo/evidence/structured-label-free-unit-assignment/t13/implementation/" *
    "closure-v1-publication-receipt.sha256"
const _CLOSURE_V1_CLAIM_SHA256 =
    "30cdfbf9be939104f76b93aa0ed610b17c5db20079af0ae32dfe781cecccc183"
const _CLOSURE_V1_CLAIM_SIDECAR_SHA256 =
    "7d642697b5940bfae7ef6255abf827a562fc483effb59d4be218642e658c75f1"
const _CLOSURE_V1_REVIEW_SHA256 =
    "1b383bda991f113b31923700a57731af56834c4ef0c2e7b3e8743e6ee24981c0"
const _CLOSURE_V1_REVIEW_SIDECAR_SHA256 =
    "e1067e8aa61efbd4311c5ece1f924a1ece575b13dd859b3863ee9e12a9f78ad4"
const _CLOSURE_V1_ROOT_MANIFEST_SHA256 =
    "a8f40846213cd976e4ffdd1c8d40e078718f30a917d9ad5ffd4fd38ef3a76aef"
const _CLOSURE_V1_PUBLICATION_RECEIPT_SHA256 =
    "59cf1c60006c549c6df6c8c07e9e9915494d5ed08d392095fa9a1cbfa5896249"
const _CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_SHA256 =
    "f2896081e273b9a9a9d267a7bbc162702f62ba9f81a64d3aae00fac8bfbac64f"

function _v4_table(document, name::String)
    table = get(document, name, nothing)
    table isa AbstractDict ||
        _fail(:BLOCKED, :policy_schema_mismatch, "runtime-v4 $name table is absent")
    return table
end

function _v4_value(table, key::String, expected, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :policy_schema_mismatch, "$context.$key is absent")
    value = table[key]
    if expected isa AbstractVector
        value isa AbstractVector ||
            _fail(:BLOCKED, :policy_type_mismatch, "$context.$key has the wrong type")
        collect(value) == expected ||
            _fail(:BLOCKED, :policy_value_mismatch, "$context.$key differs")
    else
        value isa typeof(expected) ||
            _fail(:BLOCKED, :policy_type_mismatch, "$context.$key has the wrong type")
        value == expected ||
            _fail(:BLOCKED, :policy_value_mismatch, "$context.$key differs")
    end
    return value
end

function _v4_path_value(path::String, context::String)
    (!isempty(path) && !isabspath(path) &&
     !any(==(".."), split(replace(path, '\\' => '/'), '/'; keepempty=true))) ||
        _fail(:BLOCKED, :authority_path_invalid, "invalid path for $context")
    return path
end

function _v4_bool(table, key::String, expected::Bool, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :policy_schema_mismatch, "$context.$key is absent")
    table[key] isa Bool ||
        _fail(:BLOCKED, :policy_type_mismatch, "$context.$key is not boolean")
    table[key] == expected ||
        _fail(:BLOCKED, :policy_value_mismatch, "$context.$key is not fail-closed")
    return nothing
end

function _v4_mode(mode, context::String)::String
    mode isa String && occursin(r"^0[0-7]{3}$", mode) ||
        _fail(:BLOCKED, :policy_type_mismatch, "$context mode is invalid")
    return mode
end

function _validate_runtime_v4_authority_config(authority)
    authority isa AbstractDict ||
        _fail(:BLOCKED, :policy_schema_mismatch, "runtime-v4 authority table is absent")
    expected_keys = Set(String.(collect(_RUNTIME_V4_AUTHORITY_NONBINDING_KEYS)))
    for (path_key, _, hash_key, _, mode_key, _) in _RUNTIME_V4_AUTHORITY_BINDINGS
        push!(expected_keys, path_key)
        push!(expected_keys, hash_key)
        push!(expected_keys, mode_key)
    end
    actual_keys = Set(String(key) for key in keys(authority))
    actual_keys == expected_keys ||
        _fail(:BLOCKED, :authority_binding_mismatch,
              "runtime-v4 authority key set differs")
    path_keys = Set{String}()
    for key0 in keys(authority)
        key = String(key0)
        endswith(key, "_path") || continue
        push!(path_keys, key)
        hash_key = key[1:end-5] * "_sha256"
        haskey(authority, hash_key) ||
            _fail(:BLOCKED, :authority_binding_mismatch, "missing hash for $key")
        path = authority[key]
        digest = authority[hash_key]
        path isa String && digest isa String ||
            _fail(:BLOCKED, :authority_binding_mismatch, "binding has the wrong type for $key")
        _v4_path_value(path, "authority.$key")
        _validate_hash(digest, "authority.$hash_key")
        mode_key = key[1:end-5] * "_mode"
        haskey(authority, mode_key) && _v4_mode(authority[mode_key], "authority.$key")
    end
    for key0 in keys(authority)
        key = String(key0)
        endswith(key, "_mode") || continue
        path_key = key[1:end-5] * "_path"
        haskey(authority, path_key) ||
            _fail(:BLOCKED, :authority_binding_mismatch, "mode has no path for $key")
    end
    for (path_key, expected_path, hash_key, expected_hash, mode_key, expected_mode) in
        _RUNTIME_V4_AUTHORITY_BINDINGS
        _v4_value(authority, path_key, expected_path, "authority")
        _v4_value(authority, hash_key, expected_hash, "authority")
        _v4_value(authority, mode_key, expected_mode, "authority")
    end
    return nothing
end

function _validate_runtime_v4_document(document)
    evaluator = _v4_table(document, "evaluator")
    _v4_value(evaluator, "schema", "structured_label_free_unit_assignment_evaluator_v1_correction2", "evaluator")
    _v4_value(evaluator, "version", 2, "evaluator")
    _v4_value(evaluator, "runtime_generation", 4, "evaluator")
    _v4_bool(evaluator, "allow_extra_keys", false, "evaluator")
    _v4_value(evaluator, "status_precedence", ["BLOCKED", "FAIL", "SKIPPED", "PASS"], "evaluator")
    _v4_value(evaluator, "authority_state", "gate5_implementation_candidate_pending_gate2_review", "evaluator")
    _v4_bool(evaluator, "authorizes_todo13", false, "evaluator")

    prerequisite = _v4_table(document, "prerequisite")
    _v4_value(prerequisite, "scope", "policy_and_evidence_only", "prerequisite")
    _v4_value(prerequisite, "julia", "1.12.7", "prerequisite")
    for key in ("execution_used", "real_data_used", "labels_used", "expected_N_used",
                "NKNNKN_used", "class_counts_used", "composition_prior_used",
                "equivalence_used", "benchmark_or_grader_used")
        _v4_bool(prerequisite, key, false, "prerequisite")
    end

    history = _v4_table(document, "history")
    for (key, value) in (
        ("original_config_path", _HISTORICAL_CONFIG_PATH),
        ("original_config_sha256", _CONFIG_SHA256),
        ("runtime_v1_config_path", _RUNTIME_V1_CONFIG_PATH),
        ("runtime_v1_config_sha256", _RUNTIME_V1_CONFIG_SHA256),
        ("runtime_v2_config_path", _RUNTIME_V2_CONFIG_PATH),
        ("runtime_v2_config_sha256", _RUNTIME_V2_CONFIG_SHA256),
        ("runtime_v3_config_path", _RUNTIME_V3_CONFIG_PATH),
        ("runtime_v3_config_sha256", _RUNTIME_V3_CONFIG_SHA256),
        ("scope", "history_only"),
    )
        _v4_value(history, key, value, "history")
    end
    _v4_bool(history, "runtime_equivalence", false, "history")

    runtime = _v4_table(document, "runtime")
    _v4_value(runtime, "schema", "structured_evaluator_runtime_authority_v4", "runtime")
    _v4_value(runtime, "version", 4, "runtime")
    _v4_value(runtime, "julia", "1.12.7", "runtime")
    _v4_value(runtime, "executable_sha256", _RUNTIME_V4_EXECUTABLE_SHA256, "runtime")
    _v4_value(runtime, "sysimage_sha256", _RUNTIME_V4_SYSIMAGE_SHA256, "runtime")
    for key in ("cross_runtime_compare", "julia_1_12_6_equivalence", "equivalence")
        _v4_bool(runtime, key, false, "runtime")
    end
    _v4_value(runtime, "executable", "/raven/u/system/soft/SLE_15/packages/x86_64/julia/1.12.7/bin/julia", "runtime")
    _v4_value(runtime, "sysimage", "/raven/u/system/soft/SLE_15/packages/x86_64/julia/1.12.7/lib/julia/sys.so", "runtime")

    state = _v4_table(document, "state")
    for (key, value) in (("gate4", "final"), ("next", "Gate 2 Oracle integration review"))
        _v4_value(state, key, value, "state")
    end
    _v4_bool(state, "gate5_effective", true, "state")
    _v4_bool(state, "gate5_work_started", true, "state")
    for key in ("evaluator_execution_ready", "evaluator_runtime_ready", "runtime_v4_certified",
                "todo13_authorized", "t14_authorized", "downstream_authorized")
        _v4_bool(state, key, false, "state")
    end

    outputs = _v4_table(document, "outputs")
    _v4_value(outputs, "state", "ABSENT_AT_GATE5_CONFIG_PUBLICATION", "outputs")
    _v4_bool(outputs, "execution", false, "outputs")
    _v4_value(outputs, "expected_count", 9, "outputs")
    _v4_value(outputs, "present_count", 0, "outputs")
    _v4_value(outputs, "expected_names", collect(_FILES), "outputs")
    _v4_value(outputs, "absent_names", collect(_FILES), "outputs")
    _v4_value(outputs, "present_files", String[], "outputs")
    _v4_bool(outputs, "output_hashes_recorded", false, "outputs")

    _validate_runtime_v4_authority_config(_v4_table(document, "authority"))
    return nothing
end

function _runtime_v4_config_bytes_guard(bytes::AbstractVector{UInt8})
    digest = _hash_bytes(bytes)
    digest == _RUNTIME_V4_CONFIG_SHA256 ||
        _fail(:BLOCKED, :authority_hash_mismatch, "runtime-v4 config bytes differ")
    return nothing
end

function _runtime_v4_value(runtime, key::String)
    value = runtime isa AbstractDict ? get(runtime, key, nothing) : getproperty(runtime, Symbol(key))
    value === nothing && _fail(:BLOCKED, :policy_schema_mismatch, "runtime.$key is absent")
    return value
end

function _runtime_v4_resolve_file(path::AbstractString, context::String,
                                   reason::Symbol;
                                   allow_final_symlink::Bool=false)::String
    text = String(path)
    isempty(text) && _fail(:BLOCKED, reason, "$context is absent")
    ispath(text) || _fail(:BLOCKED, reason, "$context is absent")
    isfile(text) || _fail(:BLOCKED, reason, "$context is not regular")
    !allow_final_symlink && islink(text) &&
        _fail(:BLOCKED, reason, "$context has a symlinked final component")
    canonical = try
        realpath(text)
    catch error
        _fail(:BLOCKED, reason,
              "$context cannot be resolved: $(sprint(showerror, error))")
    end
    isfile(canonical) || _fail(:BLOCKED, reason, "$context is not regular")
    try
        open(canonical, "r") do io
            read(io, 1)
        end
    catch error
        _fail(:BLOCKED, reason,
              "$context is unreadable: $(sprint(showerror, error))")
    end
    return canonical
end

function _runtime_v4_file_hash(path::AbstractString, context::String;
                               reason::Symbol=:runtime_executable_mismatch,
                               allow_final_symlink::Bool=false)::String
    canonical = _runtime_v4_resolve_file(path, context, reason;
                                         allow_final_symlink=allow_final_symlink)
    before = _identity(canonical)
    bytes = try
        Vector{UInt8}(read(canonical))
    catch error
        _fail(:BLOCKED, reason,
              "$context is unreadable: $(sprint(showerror, error))")
    end
    _identity(canonical) == before ||
        _fail(:BLOCKED, reason, "$context changed during read")
    return _hash_bytes(bytes)
end

function _runtime_v4_sysimage_path(image_pointer)
    image_pointer == C_NULL &&
        _fail(:BLOCKED, :runtime_sysimage_mismatch,
              "active Julia sysimage is unavailable")
    image_pointer isa Ptr{UInt8} ||
        _fail(:BLOCKED, :runtime_sysimage_mismatch,
              "active Julia sysimage pointer has the wrong type")
    path = try
        unsafe_string(image_pointer)
    catch error
        _fail(:BLOCKED, :runtime_sysimage_mismatch,
              "active Julia sysimage path is unavailable: $(sprint(showerror, error))")
    end
    isempty(path) &&
        _fail(:BLOCKED, :runtime_sysimage_mismatch,
              "active Julia sysimage path is empty")
    return path
end

function _runtime_v4_identity_paths(runtime;
                                    version=VERSION,
                                    configured_executable_path::AbstractString,
                                    configured_sysimage_path::AbstractString,
                                    active_executable_path::AbstractString,
                                    active_sysimage_path::AbstractString,
                                    proc_executable_path::AbstractString,
                                    allow_final_symlinks::Bool=true)
    version == _RUNTIME_V4_JULIA ||
        _fail(:BLOCKED, :runtime_version_mismatch, "Julia 1.12.7 is required")
    _runtime_v4_value(runtime, "julia") == "1.12.7" ||
        _fail(:BLOCKED, :runtime_version_mismatch, "runtime Julia version differs")
    expected_executable = _runtime_v4_value(runtime, "executable_sha256")
    expected_sysimage = _runtime_v4_value(runtime, "sysimage_sha256")
    _validate_hash(expected_executable, "runtime executable hash")
    _validate_hash(expected_sysimage, "runtime sysimage hash")

    configured_executable = _runtime_v4_resolve_file(
        configured_executable_path, "configured Julia executable",
        :runtime_executable_mismatch; allow_final_symlink=allow_final_symlinks)
    active_executable = _runtime_v4_resolve_file(
        active_executable_path, "active Julia executable", :runtime_executable_mismatch;
        allow_final_symlink=allow_final_symlinks)
    proc_executable = _runtime_v4_resolve_file(
        proc_executable_path, "running Julia executable", :runtime_executable_mismatch;
        allow_final_symlink=true)
    configured_executable == active_executable == proc_executable ||
        _fail(:BLOCKED, :runtime_executable_mismatch,
              "configured, active, and running Julia executable identities differ")

    configured_sysimage = _runtime_v4_resolve_file(
        configured_sysimage_path, "configured Julia sysimage",
        :runtime_sysimage_mismatch; allow_final_symlink=allow_final_symlinks)
    active_sysimage = _runtime_v4_resolve_file(
        active_sysimage_path, "active Julia sysimage", :runtime_sysimage_mismatch;
        allow_final_symlink=allow_final_symlinks)
    configured_sysimage == active_sysimage ||
        _fail(:BLOCKED, :runtime_sysimage_mismatch,
              "configured and active Julia sysimage identities differ")

    configured_executable_hash = _runtime_v4_file_hash(
        configured_executable, "configured Julia executable";
        reason=:runtime_executable_mismatch)
    active_executable_hash = _runtime_v4_file_hash(
        active_executable, "active Julia executable";
        reason=:runtime_executable_mismatch)
    proc_executable_hash = _runtime_v4_file_hash(
        proc_executable, "running Julia executable";
        reason=:runtime_executable_mismatch)
    configured_executable_hash == active_executable_hash == proc_executable_hash ==
        expected_executable ||
        _fail(:BLOCKED, :runtime_executable_mismatch,
              "active Julia executable bytes differ")

    configured_sysimage_hash = _runtime_v4_file_hash(
        configured_sysimage, "configured Julia sysimage";
        reason=:runtime_sysimage_mismatch)
    active_sysimage_hash = _runtime_v4_file_hash(
        active_sysimage, "active Julia sysimage";
        reason=:runtime_sysimage_mismatch)
    configured_sysimage_hash == active_sysimage_hash == expected_sysimage ||
        _fail(:BLOCKED, :runtime_sysimage_mismatch,
              "active Julia sysimage bytes differ")
    return nothing
end

function _validate_runtime_v4_identity(runtime; version=VERSION)
    version == _RUNTIME_V4_JULIA ||
        _fail(:BLOCKED, :runtime_version_mismatch, "Julia 1.12.7 is required")
    configured_executable = _runtime_v4_value(runtime, "executable")
    configured_sysimage = _runtime_v4_value(runtime, "sysimage")
    configured_executable isa AbstractString && configured_sysimage isa AbstractString ||
        _fail(:BLOCKED, :policy_type_mismatch, "runtime executable/sysimage paths are invalid")
    active_executable = try
        command = Base.julia_cmd().exec[1]
        command isa AbstractString || (command = String(command))
        isabspath(command) ? String(command) : joinpath(Sys.BINDIR, String(command))
    catch error
        _fail(:BLOCKED, :runtime_executable_mismatch,
              "active Julia executable is unavailable: $(sprint(showerror, error))")
    end
    Sys.isunix() && ispath("/proc/self/exe") ||
        _fail(:BLOCKED, :runtime_executable_mismatch,
              "running Julia executable is unavailable")
    active_sysimage = _runtime_v4_sysimage_path(Base.JLOptions().image_file)
    return _runtime_v4_identity_paths(runtime;
        version=version,
        configured_executable_path=configured_executable,
        configured_sysimage_path=configured_sysimage,
        active_executable_path=active_executable,
        active_sysimage_path=active_sysimage,
        proc_executable_path="/proc/self/exe",
        allow_final_symlinks=false)
end

function _runtime_v4_expected_mode(path::String)
    for (path_key, expected_path, _, _, _, expected_mode) in _RUNTIME_V4_AUTHORITY_BINDINGS
        endswith(path, "/" * expected_path) && return parse(UInt, expected_mode; base=8)
    end
    return nothing
end

function _validate_runtime_v4_live_hashes(; current_boulder_sha256,
                                          expected_boulder_sha256=_RUNTIME_V4_BOULDER_SHA256,
                                          gate4_transition_receipt_sha256,
                                          gate4_post_review_sha256,
                                          gate4_oracle_final_review_sha256)
    current_boulder_sha256 == expected_boulder_sha256 ||
        _fail(:BLOCKED, :authority_hash_mismatch, "current Boulder bytes differ")
    gate4_transition_receipt_sha256 ==
        "87818b63b132d8856dec22cb389cc936283bb0827cfa4884cbfd6e8c169de092" ||
        _fail(:BLOCKED, :authority_hash_mismatch, "Gate4 transition receipt binding differs")
    gate4_post_review_sha256 ==
        "80a133ba8b13d8ac29b1f8d9dcbb07d5c10ee19d47dc60e5105d0897b8032e13" ||
        _fail(:BLOCKED, :authority_hash_mismatch, "Gate4 post-review binding differs")
    gate4_oracle_final_review_sha256 ==
        "cab6c866c256ce160c13b4694531f661b2db9cbfac5dd47b34cf69b073635fbd" ||
        _fail(:BLOCKED, :authority_hash_mismatch, "Gate4 Oracle-review binding differs")
    return nothing
end

function _validate_runtime_v4_execution_gate(document)
    evaluator = _v4_table(document, "evaluator")
    state = _v4_table(document, "state")
    outputs = _v4_table(document, "outputs")
    candidate = get(evaluator, "authorizes_todo13", true) === false &&
                get(state, "evaluator_execution_ready", true) === false &&
                get(state, "runtime_v4_certified", true) === false &&
                get(outputs, "execution", true) === false &&
                get(outputs, "present_count", -1) == 0
    candidate &&
        _fail(:BLOCKED, :todo_state_mismatch,
              "runtime-v4 Gate5 candidate is not authorized for evaluator execution")
    _fail(:BLOCKED, :todo_state_mismatch,
          "runtime-v4 execution state is not an authorized certified state")
    return nothing
end

function _identity(path::String)
    info = stat(path)
    return (UInt64(info.device), UInt64(info.inode), UInt64(info.nlink))
end

function _reject_link_components(path::String)
    parts = splitpath(abspath(path))
    cursor = first(parts)
    islink(cursor) && _fail(:BLOCKED, :authority_symlink, "symlinked authority ancestor")
    for part in parts[2:end]
        cursor = joinpath(cursor, part)
        islink(cursor) && _fail(:BLOCKED, :authority_symlink, "symlinked authority path")
    end
    return nothing
end

function _repo_path(root::String, supplied::AbstractString, context::String)
    text = String(supplied)
    isempty(text) && _fail(:BLOCKED, :authority_path_invalid, "$context is empty")
    occursin('\0', text) && _fail(:BLOCKED, :authority_path_invalid, "$context contains NUL")
    any(==(".."), split(replace(text, '\\' => '/'), '/'; keepempty=true)) &&
        _fail(:BLOCKED, :authority_path_invalid, "$context contains parent traversal")
    absolute = normpath(isabspath(text) ? text : joinpath(root, text))
    prefix = root == "/" ? "/" : root * "/"
    (absolute == root || startswith(absolute, prefix)) ||
        _fail(:BLOCKED, :authority_path_invalid, "$context escapes root")
    _reject_link_components(absolute)
    ispath(absolute) ||
        _fail(:BLOCKED, :authority_path_invalid, "$context is absent")
    realpath(absolute) == absolute ||
        _fail(:BLOCKED, :authority_path_mismatch, "$context is not canonical")
    return absolute
end

function _snapshot(root::String, supplied::AbstractString, expected::String,
                   context::String, seen::Dict{Tuple{UInt64,UInt64},String})::_Snapshot
    path = _repo_path(root, supplied, context)
    isfile(path) || _fail(:BLOCKED, :authority_path_invalid, "$context is not a regular file")
    identity = _identity(path)
    identity[3] == 1 || _fail(:BLOCKED, :authority_hardlink, "$context is hard-linked")
    key = (identity[1], identity[2])
    haskey(seen, key) && _fail(:BLOCKED, :authority_identity_collision,
                               "$context shares a file identity with $(seen[key])")
    seen[key] = path
    before = identity
    bytes = open(path, "r") do io
        _identity(path) == before || _fail(:BLOCKED, :authority_snapshot_changed, "$context changed on open")
        read(io)
    end
    _identity(path) == before ||
        _fail(:BLOCKED, :authority_snapshot_changed, "$context changed after read")
    digest = _hash_bytes(bytes)
    digest == expected ||
        _fail(:BLOCKED, :authority_hash_mismatch, "$context bytes differ")
    return _Snapshot(path, Vector{UInt8}(bytes), digest, identity[1], identity[2], identity[3])
end

function _verify_snapshot(snapshot::_Snapshot, context::String)
    _reject_link_components(snapshot.path)
    isfile(snapshot.path) || _fail(:BLOCKED, :authority_snapshot_changed, "$context disappeared")
    _is_immutable_closure_authority_path(snapshot.path) &&
        _readonly_mode(snapshot.path, context)
    identity = _identity(snapshot.path)
    (identity[1], identity[2], identity[3]) ==
        (snapshot.device, snapshot.inode, snapshot.nlink) ||
        _fail(:BLOCKED, :authority_snapshot_changed, "$context identity changed")
    _runtime_v3_immutable_path(snapshot.path) &&
        _runtime_v3_mode(snapshot.path, UInt(0o444), context)
    runtime_v4_mode = _runtime_v4_expected_mode(snapshot.path)
    runtime_v4_mode === nothing || _runtime_v3_mode(snapshot.path, runtime_v4_mode, context)
    read(snapshot.path) == snapshot.bytes ||
        _fail(:BLOCKED, :authority_snapshot_changed, "$context bytes changed")
    return nothing
end

function _authority_entry(table, path_key::String, hash_key::String)
    haskey(table, path_key) && haskey(table, hash_key) ||
        _fail(:BLOCKED, :authority_binding_mismatch, "authority binding is incomplete")
    path = table[path_key]
    digest = table[hash_key]
    path isa String && digest isa String ||
        _fail(:BLOCKED, :authority_binding_mismatch, "authority binding has the wrong type")
    _validate_hash(digest, "authority.$hash_key")
    return String(path), String(digest)
end

function _check_plan(root::String, path::String, seen)
    absolute = _repo_path(root, path, "plan")
    isfile(absolute) || _fail(:BLOCKED, :authority_path_invalid, "plan is absent")
    actual = _hash_bytes(read(absolute))
    actual == _PLAN_SHA256 && return _snapshot(root, path, _PLAN_SHA256, "plan", seen)
    bytes = String(read(absolute))
    marker = "- [x] 13."
    count(marker, bytes) == 1 ||
        _fail(:BLOCKED, :authority_binding_mismatch, "checked Todo13 marker is absent")
    reversed = replace(bytes, marker => "- [ ] 13."; count=1)
    _hash_text(reversed) == _PLAN_SHA256 ||
        _fail(:BLOCKED, :authority_binding_mismatch, "plan is not marker-reversible")
    return _snapshot(root, path, actual, "plan", seen)
end

function _json_true(text::AbstractString, key::AbstractString)::Bool
    return occursin(Regex("\\\"" * String(key) * "\\\"\\s*:\\s*true"), String(text))
end

function _json_false(text::AbstractString, key::AbstractString)::Bool
    return occursin(Regex("\\\"" * String(key) * "\\\"\\s*:\\s*false"), String(text))
end

function _json_string(text::AbstractString, key::AbstractString,
                      value::AbstractString)::Bool
    # The accepted authorities are hash-bound snapshots.  This small parser is
    # deliberately only a field-binding check; it does not treat arbitrary
    # JSON text as an authority.
    return occursin(Regex("\\\"" * String(key) * "\\\"\\s*:\\s*\\\"" *
                         String(value) * "\\\""), String(text))
end

function _json_number(text::AbstractString, key::AbstractString, value::Int)::Bool
    return occursin(Regex("\\\"" * String(key) * "\\\"\\s*:\\s*" *
                         string(value) * "(?:\\s*[,}])"), String(text))
end

function _json_number_text(text::AbstractString, key::AbstractString,
                           value::AbstractString)::Bool
    return occursin(Regex("\\\"" * String(key) * "\\\"\\s*:\\s*" *
                          String(value) * "(?:\\s*[,}])"), String(text))
end

function _require_authority(condition::Bool, message::AbstractString)
    condition || _fail(:BLOCKED, :authority_binding_mismatch, message)
    return nothing
end

function _validate_sidecar(snapshot::_Snapshot, expected::String, target::String,
                           context::String)
    String(copy(snapshot.bytes)) == "$expected  $target\n" ||
        _fail(:BLOCKED, :authority_binding_mismatch, "$context sidecar binding differs")
    return nothing
end

function _readonly_mode(path::String, context::String)
    mode = UInt(stat(path).mode)
    mode & UInt(0o222) == 0 ||
        _fail(:BLOCKED, :authority_path_mismatch, "$context is writable")
    return nothing
end

function _authority_config_kind(root::String, supplied::AbstractString)::Symbol
    text = String(supplied)
    absolute = normpath(isabspath(text) ? text : joinpath(root, text))
    historical = normpath(joinpath(root, _HISTORICAL_CONFIG_PATH))
    runtime_v1 = normpath(joinpath(root, _RUNTIME_V1_CONFIG_PATH))
    runtime_v2 = normpath(joinpath(root, _RUNTIME_V2_CONFIG_PATH))
    runtime_v3 = normpath(joinpath(root, _RUNTIME_V3_CONFIG_PATH))
    runtime_v4 = normpath(joinpath(root, _RUNTIME_V4_CONFIG_PATH))
    absolute == historical && return :historical
    absolute == runtime_v1 && return :runtime
    absolute == runtime_v2 && return :runtime_v2
    absolute == runtime_v3 && return :runtime_v3
    absolute == runtime_v4 && return :runtime_v4
    _fail(:BLOCKED, :authority_path_mismatch,
          "evaluator config path is not an accepted authority path")
end

function _preflight_config_context(options::Dict{String,String};
                                   activation::Union{Nothing,_ActivationAuthorization}=nothing)::_ProductionContext
    root = options["--root"]
    isabspath(root) ||
        _fail(:BLOCKED, :authority_path_invalid, "root must be absolute")
    isdir(root) ||
        _fail(:BLOCKED, :authority_path_invalid, "root is not a directory")
    root = normpath(root)
    realpath(root) == root ||
        _fail(:BLOCKED, :authority_path_invalid, "root is not canonical")
    VERSION == _RUNTIME_V4_JULIA ||
        _fail(:BLOCKED, :runtime_version_mismatch, "Julia 1.12.7 is required")

    supplied = options["--evaluator-config"]
    kind = _authority_config_kind(root, supplied)
    kind == :runtime_v4 ||
        _fail(:BLOCKED, :authority_path_mismatch,
              "only the runtime-v4 evaluator config is executable")
    relative = _RUNTIME_V4_CONFIG_PATH
    config_path = _repo_path(root, supplied, "evaluator config")
    config_path == normpath(joinpath(root, relative)) ||
        _fail(:BLOCKED, :authority_path_mismatch, "evaluator config path differs")
    activation === nothing ||
        _validate_activation_receiving_options(activation, root, options)
    authority = _authority_snapshots(root, supplied; activation=activation)
    return _assemble_production_context(root, options, authority;
                                        activation=activation)
end

function _historical_blocker_context(options::Dict{String,String})
    root = options["--root"]
    isabspath(root) || return nothing
    isdir(root) || return nothing
    root = normpath(root)
    realpath(root) == root || return nothing
    supplied = options["--evaluator-config"]
    kind = try
        _authority_config_kind(root, supplied)
    catch
        return nothing
    end
    kind == :runtime_v4 && return nothing
    relative, expected = if kind == :historical
        (_HISTORICAL_CONFIG_PATH, _CONFIG_SHA256)
    elseif kind == :runtime
        (_RUNTIME_V1_CONFIG_PATH, _RUNTIME_V1_CONFIG_SHA256)
    elseif kind == :runtime_v2
        (_RUNTIME_V2_CONFIG_PATH, _RUNTIME_V2_CONFIG_SHA256)
    else
        (_RUNTIME_V3_CONFIG_PATH, _RUNTIME_V3_CONFIG_SHA256)
    end
    config_path = try
        _repo_path(root, supplied, "evaluator config")
    catch
        return nothing
    end
    config_path == normpath(joinpath(root, relative)) || return nothing
    snapshot = try
        _snapshot(root, supplied, expected, "evaluator config",
                  Dict{Tuple{UInt64,UInt64},String}())
    catch
        return nothing
    end
    inventories = Tuple{String,Vector{String}}[]
    return _ProductionContext(root, [snapshot], inventories,
                              _manifest_digest(root, [snapshot], inventories))
end

# Successor Main-control-only schemas.  The previous v1 activation route is
# rejected by schema equality: there is no fallback from v2 to v1, and the v1
# bytes certify only the earlier execution-root-owned Boulder design.
const _T13_ACTIVATION_SCHEMA = "stmfit_t13_activation_receipt_v2"
const _T13_EXECUTION_SPEC_SCHEMA = "stmfit_t13_execution_spec_v2"
# Fixed Main-control role and filename.  No CLI override, arbitrary control
# filename, multi-root whitelist or fallback discovery is accepted.
const _T13_CONTROL_ROLE = "main_boulder"
const _T13_CONTROL_FILENAME = ".omo/boulder.json"
const _T13_CONTROL_PARENT = ".omo"
const _T13_CONTROL_KEYS = ("role", "root", "root_mode", "parent_mode", "file_mode")
const _T13_SYNTHETIC_SCOPE = "T13_SYNTHETIC_EXECUTION_ONLY"
const _T13_SYNTHETIC_ENVIRONMENT = "STMFIT_T13_SYNTHETIC_EXECUTION_ONLY"
const _T13_GATE5_PUBLICATION_PATH =
    ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/runtime_v4/" *
    "t13-gate5-implementation-evidence-v1/PublicationReceipt.json"
const _T13_GATE5_PUBLICATION_SHA256 =
    "ae1af8edb963c632988a0719a1e7b0def2c2ad689628f91435aa238b355aad83"
const _T13_GATE2_ORACLE_REVIEW_PATH =
    ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/runtime_v4/" *
    "t13-gate5-gate2-oracle-review-v1/Gate2OracleReview.json"
const _T13_GATE2_ORACLE_REVIEW_SHA256 =
    "fe0b60720af63655f5415def50ef79a3096858840ed9d4f386e24569544e6fbe"
const _T13_GATE2_ORACLE_REVIEW_MODE = UInt(0o444)
const _T13_GATE5_PUBLICATION_MODE = UInt(0o444)
# New activation control records must be reviewed read-only.  This is a
# captured-mode comparison, not a sealing operation: the code never chmods.
const _T13_CONTROL_RECORD_MODE = UInt(0o444)
# Successor authority-manifest version for an activated T13 context.  The
# non-activated runtime-v4 manifest keeps its historical v6 schema.
const _T13_ACTIVATION_AUTHORITY_SCHEMA =
    "schema=structured-evaluator-authority-v8-t13-main-control"
const _T13_ACTIVATED_REPORT_SCHEMA = _SCHEMA * "_activation_receipt_v2"
const _T13_ACTIVATED_BLOCKER_SCHEMA = _SCHEMA * "_activation_blocker_receipt_v2"
const _T13_LAUNCH_KEYS = ("active_project", "cwd", "startup_file", "history_file",
                          "threads", "blas_threads", "depot_path", "load_path",
                          "offline", "forbidden_environment")
const _T13_FORBIDDEN_ENVIRONMENT = (
    "STMFIT_DATA_DIR", "STMFIT_T13_PUBLIC_REVALIDATION", "JULIA_DEBUG",
)
const _T13_REBIND_PREFIX =
    ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/provenance-rebind/"
# Fixed provenance-only reads derived from the unchanged T8/T11/T12 readers.
const _T13_FIXED_PROVENANCE_FILES = (
    ".omo/plans/structured-label-free-unit-assignment.md",
    ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/t7/correction/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase4-t10-edge-features/review/AdversarialVerify.json",
    ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/t2/correction2/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase3-t3-universe/correction/review/AdversarialVerify.json",
    "test/build_structured_unit_predictions.jl",
    "test/test_structured_robust_emissions.jl",
    "test/test_structured_edge_features.jl",
    "test/build_label_free_edge_features.jl",
    _T13_REBIND_PREFIX * "phase3-t3-universe/correction/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase4-t10-edge-features/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/correction/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/correction/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/correction/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/correction/claim-correction/ClaimCorrection.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/correction/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/correction2/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/correction2/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/correction2/canonical-t11-source-bundle.bytes",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/scale-certification-correction/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/scale-certification-correction/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/scale-certification-correction/canonical-t11-source-bundle.bytes",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/DoneClaim.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/review/AdversarialVerify.json",
    _T13_REBIND_PREFIX * "phase5-t11-admission/graph-handoff-extension/scale-certification-correction/checked-plan-rebind/canonical-t11-source-bundle.bytes",
)
const _T13_HIERARCHICAL_DIRECTORY = "test/lib/hierarchical"
const _T13_PRODUCER_SCRIPTS = (
    (path="test/extract_lobe_patches.jl",
     sha256=T11.StructuredEdgeFeatures.FORWARD_PRODUCER_SHA256),
    (path="test/extract_lobe_patches_bwd.jl",
     sha256=T11.StructuredEdgeFeatures.BACKWARD_PRODUCER_SHA256),
)

# The permitted activation command body.  The externally supplied
# --activation-sha256 digest is deliberately absent: it is the launch trust
# anchor, and recording it here would make the spec/receipt chain cyclic.
const _T13_ACTIVATION_COMMAND_FLAGS = (
    "--root", "--evaluator-config", "--features", "--candidate-config",
    "--model-config", "--universe-dir", "--edge-dir", "--forward-receipt",
    "--backward-receipt", "--admission-dir", "--out-dir",
    "--activation-receipt",
)

# Consumed include closure of the evaluator entrypoint at this revision.
const _T13_SOURCE_CLOSURE = (
    "test/evaluate_structured_unit_assignment.jl",
    "test/lib/structured_unit_assignment.jl",
    "test/lib/hierarchical_unit_assignment.jl",
    "test/lib/hierarchical/io_helpers.jl",
    "test/lib/hierarchical/firewall.jl",
    "test/lib/hierarchical/loading.jl",
    "test/lib/hierarchical/nuisance.jl",
    "test/lib/hierarchical/emission_math.jl",
    "test/lib/hierarchical/emissions.jl",
    "test/lib/hierarchical/views.jl",
    "test/lib/hierarchical/pipeline.jl",
    "test/lib/structured_assignment/firewall.jl",
    "test/lib/structured_assignment/champion_adapter.jl",
    "test/lib/structured_assignment/robust_emissions.jl",
    "test/lib/structured_assignment/edge_admission.jl",
    "test/lib/structured_assignment/edge_features.jl",
    "test/lib/structured_assignment/universe.jl",
    "test/lib/structured_assignment/chain_inference.jl",
    "test/lib/structured_assignment/edge_model.jl",
)

function _t13_relative_path_value(text::AbstractString, context::String)
    value = String(text)
    isempty(value) && _fail(:BLOCKED, :activation_path_invalid, "$context is empty")
    occursin('\0', value) &&
        _fail(:BLOCKED, :activation_path_invalid, "$context contains NUL")
    isabspath(value) &&
        _fail(:BLOCKED, :activation_path_invalid, "$context must be repository-relative")
    occursin('\\', value) &&
        _fail(:BLOCKED, :activation_path_invalid, "$context uses a backslash")
    parts = split(value, '/'; keepempty=true)
    any(part -> isempty(part) || part == "." || part == "..", parts) &&
        _fail(:BLOCKED, :activation_path_invalid, "$context is not normalized")
    normpath(value) == value ||
        _fail(:BLOCKED, :activation_path_invalid, "$context is not canonical")
    return value
end

function _t13_relative_path(table, key::String, context::String)
    return _t13_relative_path_value(_t13_string(table, key, context),
                                    "$context.$key")
end

function _t13_exact_keys(table, expected, context::String)
    table isa AbstractDict ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context is absent or not a table")
    actual = Set(String(key) for key in keys(table))
    expected_set = Set(String(key) for key in expected)
    actual == expected_set ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context key set differs")
    return table
end

function _t13_table(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa AbstractDict ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not a table")
    return value
end

function _t13_table_array(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa AbstractVector ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not an array")
    for entry in value
        entry isa AbstractDict ||
            _fail(:BLOCKED, :activation_schema_mismatch,
                  "$context.$key has a non-table entry")
    end
    return value
end

function _t13_string(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa AbstractString ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not a string")
    text = String(value)
    isempty(text) &&
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is empty")
    return text
end

function _t13_string_array(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa AbstractVector ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not an array")
    result = String[]
    for item in value
        item isa AbstractString ||
            _fail(:BLOCKED, :activation_schema_mismatch,
                  "$context.$key has a non-string entry")
        push!(result, String(item))
    end
    length(result) == length(unique(result)) ||
        _fail(:BLOCKED, :activation_schema_mismatch,
              "$context.$key has duplicate entries")
    return result
end

function _t13_hash(table, key::String, context::String)
    text = _t13_string(table, key, context)
    occursin(r"^[0-9a-f]{64}$", text) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not a SHA-256")
    return text
end

function _t13_mode(table, key::String, context::String)
    text = _t13_string(table, key, context)
    occursin(r"^0[0-7]{3}$", text) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not a mode")
    return parse(UInt, text; base=8)
end

_t13_mode_text(mode::UInt) = lpad(string(mode, base=8), 4, '0')

# Dedicated Main-control reader.  It exposes only the fixed
# `.omo/boulder.json` filename under the canonical Main root and may reuse the
# existing `_snapshot(M, ...)` helper, but never an arbitrary caller-selected
# filename.  Directory identities and modes are acquired before the file read,
# the file is checked regular/direct/single-link/digest, and after the read the
# file identity/link count/mode plus the directories are rechecked.  No chmod
# repair and no identity refresh.
function _read_main_control(control, receipt_control_root::String,
                            preimage_sha256::String, postimage_sha256::String)
    # `control` is the parsed `[control]` NamedTuple produced by
    # `_validate_t13_execution_spec_document`, which already enforced the exact
    # key set, canonical absolute root and mode syntax.
    role = String(control.role)
    root = String(control.root)
    root_mode = UInt(control.root_mode)
    parent_mode = UInt(control.parent_mode)
    file_mode = UInt(control.file_mode)
    role == _T13_CONTROL_ROLE ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation Main control role differs")
    (isabspath(root) && normpath(root) == root) ||
        _fail(:BLOCKED, :activation_path_invalid,
              "activation Main control root is not canonical absolute")
    root == receipt_control_root ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation Main control root differs from the receipt binding")
    parent = joinpath(root, _T13_CONTROL_PARENT)
    isdir(root) && !islink(root) ||
        _fail(:BLOCKED, :activation_control_mismatch,
              "activation Main control root is absent or not a real directory")
    _runtime_v3_mode(root, root_mode, "activation Main control root")
    root_identity = (UInt64(stat(root).device), UInt64(stat(root).inode))
    isdir(parent) && !islink(parent) ||
        _fail(:BLOCKED, :activation_control_mismatch,
              "activation Main control parent is absent or not a real directory")
    _runtime_v3_mode(parent, parent_mode, "activation Main control parent")
    parent_identity = (UInt64(stat(parent).device), UInt64(stat(parent).inode))
    # Acquire the file identity/mode/link state before the read.  A mode or
    # link failure is classified before an unreadable file can raise a raw
    # system error.
    file_path = joinpath(root, _T13_CONTROL_FILENAME)
    isfile(file_path) && !islink(file_path) ||
        _fail(:BLOCKED, :activation_control_mismatch,
              "activation Main control file is absent or not regular")
    file_info = stat(file_path)
    UInt64(file_info.nlink) == 1 ||
        _fail(:BLOCKED, :authority_hardlink,
              "activation Main control file is hard-linked")
    _runtime_v3_mode(file_path, file_mode, "activation Main control file")
    file_identity = (UInt64(file_info.device), UInt64(file_info.inode))
    seen = Dict{Tuple{UInt64,UInt64},String}()
    file = _snapshot(root, _T13_CONTROL_FILENAME, postimage_sha256,
                     "activation Main control file", seen)
    # Post-read recheck of the file identity/link count/approved mode and of the
    # captured root/parent identities and modes.  The reader returns only after
    # all of these still match the pre-read acquisition; identities are never
    # refreshed or rebased.
    _main_control_file_recheck(file.path, file_mode, file_identity,
                               UInt64(file_info.nlink),
                               "activation Main control file")
    _main_control_directory_recheck(root, root_mode, root_identity,
                                    "activation Main control root")
    _main_control_directory_recheck(parent, parent_mode, parent_identity,
                                    "activation Main control parent")
    return _MainControl(root, file, root_mode, parent_mode, file_mode,
                        root_identity, parent_identity, preimage_sha256,
                        postimage_sha256)
end

function _main_control_directory_recheck(path::String, mode::UInt,
                                         identity::Tuple{UInt64,UInt64},
                                         context::String)
    isdir(path) && !islink(path) ||
        _fail(:BLOCKED, :activation_snapshot_changed, "$context disappeared")
    _runtime_v3_mode(path, mode, context)
    (UInt64(stat(path).device), UInt64(stat(path).inode)) == identity ||
        _fail(:BLOCKED, :activation_snapshot_changed,
              "$context identity changed")
    return nothing
end

# Shared post-read file recheck for the dedicated Main-control object.  It
# reacquires the captured file identity, link count and approved mode after a
# read, so a same-bytes replacement or a mode change between the pre-read
# acquisition and the post-read revalidation is rejected.  It never refreshes
# the captured identity and never repairs a mode.
function _main_control_file_recheck(path::String, mode::UInt,
                                    identity::Tuple{UInt64,UInt64},
                                    nlink::UInt64, context::String)
    isfile(path) && !islink(path) ||
        _fail(:BLOCKED, :activation_snapshot_changed, "$context disappeared")
    info = stat(path)
    UInt64(info.nlink) == 1 ||
        _fail(:BLOCKED, :authority_hardlink, "$context is hard-linked")
    (UInt64(info.device), UInt64(info.inode), UInt64(info.nlink)) ==
        (identity[1], identity[2], nlink) ||
        _fail(:BLOCKED, :authority_snapshot_changed,
              "$context identity changed")
    _runtime_v3_mode(path, mode, context)
    return nothing
end

function _verify_main_control(control::_MainControl)
    parent = joinpath(control.root, _T13_CONTROL_PARENT)
    # Bracket the file read with the captured root/parent identity and mode
    # checks: the same dedicated object is revalidated before and after the
    # read, and the file identity/link count/mode is reacquired after the read.
    _main_control_directory_recheck(control.root, control.root_mode,
                                    control.root_identity,
                                    "activation Main control root")
    _main_control_directory_recheck(parent, control.parent_mode,
                                    control.parent_identity,
                                    "activation Main control parent")
    _verify_snapshot(control.file, "activation Main control file")
    _main_control_file_recheck(control.file.path, control.file_mode,
                               (control.file.device, control.file.inode),
                               control.file.nlink,
                               "activation Main control file")
    _main_control_directory_recheck(control.root, control.root_mode,
                                    control.root_identity,
                                    "activation Main control root")
    _main_control_directory_recheck(parent, control.parent_mode,
                                    control.parent_identity,
                                    "activation Main control parent")
    return nothing
end

function _t13_bool(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa Bool ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not boolean")
    return value
end

function _t13_int(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa Integer ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not an integer")
    return Int(value)
end

function _t13_false(table, key::String, context::String)
    haskey(table, key) ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is absent")
    value = table[key]
    value isa Bool ||
        _fail(:BLOCKED, :activation_schema_mismatch, "$context.$key is not boolean")
    value === false ||
        _fail(:BLOCKED, :activation_permission_mismatch, "$context.$key is not denied")
    return nothing
end

function _t13_parse_toml(snapshot::_Snapshot, context::String)
    return try
        TOML.parse(String(copy(snapshot.bytes)))
    catch error
        _fail(:BLOCKED, :activation_schema_mismatch,
              "$context is not valid TOML: $(sprint(showerror, error))")
    end
end

function _t13_canonical_command(options::Dict{String,String})
    body = String[]
    for flag in _T13_ACTIVATION_COMMAND_FLAGS
        haskey(options, flag) ||
            _fail(:BLOCKED, :activation_command_mismatch,
                  "activation command is missing $flag")
        push!(body, flag)
        push!(body, options[flag])
    end
    return body
end

function _t13_relative_from_option(root::String, value::AbstractString)
    text = String(value)
    isempty(text) && _fail(:BLOCKED, :activation_path_invalid, "activation input is empty")
    absolute = normpath(isabspath(text) ? text : joinpath(root, text))
    prefix = root == "/" ? "/" : root * "/"
    (absolute == root || startswith(absolute, prefix)) ||
        _fail(:BLOCKED, :activation_path_invalid, "activation input escapes root")
    relative = relpath(absolute, root)
    return _t13_relative_path_value(relative, "activation input")
end

function _t13_input_relative_paths(root::String, options::Dict{String,String})
    return [_t13_relative_from_option(root, options[flag])
            for flag in ("--features", "--candidate-config", "--model-config",
                         "--forward-receipt", "--backward-receipt")]
end

function _t13_directory_relative_paths(root::String, options::Dict{String,String})
    return [_t13_relative_from_option(root, options[flag])
            for flag in ("--universe-dir", "--edge-dir", "--admission-dir")]
end

function _t13_authority_binding(stem::String)
    target = stem * "_path"
    for (path_key, path, _, hash, _, _) in _RUNTIME_V4_AUTHORITY_BINDINGS
        path_key == target && return (path=path, sha256=hash)
    end
    _fail(:BLOCKED, :activation_predecessor_mismatch,
          "unknown activation predecessor $stem")
end

function _t13_check_predecessor(table, stem::String, expected_path::String,
                                expected_sha256::String, context::String)
    path = _t13_relative_path(table, stem * "_path", context)
    digest = _t13_hash(table, stem * "_sha256", context)
    path == expected_path && digest == expected_sha256 ||
        _fail(:BLOCKED, :activation_predecessor_mismatch,
              "$context.$stem differs from the fixed predecessor binding")
    return nothing
end

function _validate_t13_activation_receipt_document(document, receipt_rel::String)
    context = "activation receipt"
    _t13_exact_keys(document,
                    ("schema", "scope", "spec", "source_review",
                     "transition_review", "boulder", "publication"), context)
    _t13_string(document, "schema", context) == _T13_ACTIVATION_SCHEMA ||
        _fail(:BLOCKED, :activation_schema_mismatch, "activation receipt schema differs")
    scope = _t13_string(document, "scope", context)
    scope == _T13_SYNTHETIC_SCOPE ||
        _fail(:BLOCKED, :activation_scope_mismatch, "activation receipt scope differs")
    spec = _t13_exact_keys(_t13_table(document, "spec", context),
                           ("path", "sha256"), "$context.spec")
    source_review = _t13_exact_keys(_t13_table(document, "source_review", context),
                                    ("path", "sha256"), "$context.source_review")
    transition_review = _t13_exact_keys(
        _t13_table(document, "transition_review", context),
        ("path", "sha256"), "$context.transition_review")
    boulder = _t13_exact_keys(_t13_table(document, "boulder", context),
                              ("control_root",
                               "preimage_path", "preimage_sha256",
                               "postimage_path", "postimage_sha256"),
                              "$context.boulder")
    publication = _t13_exact_keys(_t13_table(document, "publication", context),
                                  ("path", "sha256"), "$context.publication")
    referenced = String[
        _t13_relative_path(spec, "path", "$context.spec"),
        _t13_relative_path(source_review, "path", "$context.source_review"),
        _t13_relative_path(transition_review, "path", "$context.transition_review"),
        _t13_relative_path(publication, "path", "$context.publication"),
    ]
    all(path -> path != receipt_rel, referenced) ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation receipt references itself")
    control_root = _t13_string(boulder, "control_root", "$context.boulder")
    (isabspath(control_root) && normpath(control_root) == control_root) ||
        _fail(:BLOCKED, :activation_path_invalid,
              "activation Boulder control_root is not canonical absolute")
    preimage_path = _t13_relative_path(boulder, "preimage_path", "$context.boulder")
    preimage_sha256 = _t13_hash(boulder, "preimage_sha256", "$context.boulder")
    postimage_path = _t13_relative_path(boulder, "postimage_path", "$context.boulder")
    postimage_sha256 = _t13_hash(boulder, "postimage_sha256", "$context.boulder")
    preimage_path == _T13_CONTROL_FILENAME ||
        _fail(:BLOCKED, :activation_state_mismatch,
              "activation Boulder preimage path is not the fixed Gate4 state path")
    postimage_path == _T13_CONTROL_FILENAME ||
        _fail(:BLOCKED, :activation_state_mismatch,
              "activation Boulder postimage path is not the fixed Main-control state path")
    preimage_sha256 == _RUNTIME_V4_BOULDER_SHA256 ||
        _fail(:BLOCKED, :activation_state_mismatch,
              "activation Boulder preimage is not the fixed Gate4 state")
    postimage_sha256 != preimage_sha256 ||
        _fail(:BLOCKED, :activation_state_mismatch,
              "activation Boulder postimage equals the Gate4 preimage")
    return (scope=scope,
            control_root=control_root,
            spec_path=_t13_relative_path(spec, "path", "$context.spec"),
            spec_sha256=_t13_hash(spec, "sha256", "$context.spec"),
            source_review_path=_t13_relative_path(source_review, "path",
                                                  "$context.source_review"),
            source_review_sha256=_t13_hash(source_review, "sha256",
                                           "$context.source_review"),
            transition_review_path=_t13_relative_path(transition_review, "path",
                                                      "$context.transition_review"),
            transition_review_sha256=_t13_hash(transition_review, "sha256",
                                               "$context.transition_review"),
            boulder_preimage_path=preimage_path,
            boulder_preimage_sha256=preimage_sha256,
            boulder_postimage_path=postimage_path,
            boulder_postimage_sha256=postimage_sha256,
            publication_path=_t13_relative_path(publication, "path",
                                                "$context.publication"),
            publication_sha256=_t13_hash(publication, "sha256",
                                         "$context.publication"))
end

function _validate_t13_execution_spec_document(document, root::String,
                                               options::Dict{String,String})
    context = "activation execution spec"
    _t13_exact_keys(document,
                    ("schema", "scope", "root", "control", "entrypoint", "config",
                     "project", "manifest", "source", "command", "launch",
                     "destination", "inputs", "dependencies", "predecessors",
                     "permissions"),
                    context)
    _t13_string(document, "schema", context) == _T13_EXECUTION_SPEC_SCHEMA ||
        _fail(:BLOCKED, :activation_schema_mismatch,
              "activation execution spec schema differs")
    scope = _t13_string(document, "scope", context)
    scope == _T13_SYNTHETIC_SCOPE ||
        _fail(:BLOCKED, :activation_scope_mismatch,
              "activation execution spec scope differs")

    root_table = _t13_exact_keys(_t13_table(document, "root", context),
                                 ("path",), "$context.root")
    root_path = _t13_string(root_table, "path", "$context.root")
    (isabspath(root_path) && normpath(root_path) == root_path) ||
        _fail(:BLOCKED, :activation_root_mismatch,
              "activation execution spec root is not canonical")
    root_path == root ||
        _fail(:BLOCKED, :activation_root_mismatch,
              "activation execution spec root differs")

    # Main-control table: exactly the fixed role, one canonical absolute Main
    # root and the three approved modes.  The spec does not contain a future
    # Boulder postimage digest; the receipt carries it, so the spec/receipt
    # chain stays acyclic.
    control = _t13_exact_keys(_t13_table(document, "control", context),
                              _T13_CONTROL_KEYS, "$context.control")
    control_role = _t13_string(control, "role", "$context.control")
    control_root = _t13_string(control, "root", "$context.control")
    (isabspath(control_root) && normpath(control_root) == control_root) ||
        _fail(:BLOCKED, :activation_path_invalid,
              "activation Main control root is not canonical absolute")
    control_root_mode = _t13_mode(control, "root_mode", "$context.control")
    control_parent_mode = _t13_mode(control, "parent_mode", "$context.control")
    control_file_mode = _t13_mode(control, "file_mode", "$context.control")

    entrypoint = _t13_exact_keys(_t13_table(document, "entrypoint", context),
                                 ("path", "sha256"), "$context.entrypoint")
    entrypoint_path = _t13_relative_path(entrypoint, "path", "$context.entrypoint")
    entrypoint_path == "test/evaluate_structured_unit_assignment.jl" ||
        _fail(:BLOCKED, :activation_entrypoint_mismatch,
              "activation execution spec entrypoint differs")
    entrypoint_sha256 = _t13_hash(entrypoint, "sha256", "$context.entrypoint")

    config = _t13_exact_keys(_t13_table(document, "config", context),
                             ("path", "sha256"), "$context.config")
    config_path = _t13_relative_path(config, "path", "$context.config")
    config_sha256 = _t13_hash(config, "sha256", "$context.config")
    config_path == _RUNTIME_V4_CONFIG_PATH &&
        config_sha256 == _RUNTIME_V4_CONFIG_SHA256 ||
        _fail(:BLOCKED, :activation_config_mismatch,
              "activation execution spec config differs from the unchanged runtime-v4 config")

    project = _t13_exact_keys(_t13_table(document, "project", context),
                              ("path", "sha256"), "$context.project")
    project_path = _t13_relative_path(project, "path", "$context.project")
    project_sha256 = _t13_hash(project, "sha256", "$context.project")
    project_path == "Project.toml" ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation Project.toml path differs")
    manifest = _t13_exact_keys(_t13_table(document, "manifest", context),
                               ("path", "sha256"), "$context.manifest")
    manifest_path = _t13_relative_path(manifest, "path", "$context.manifest")
    manifest_sha256 = _t13_hash(manifest, "sha256", "$context.manifest")
    manifest_path == "Manifest.toml" ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation Manifest.toml path differs")

    source = _t13_exact_keys(_t13_table(document, "source", context),
                             ("files",), "$context.source")
    source_entries = _t13_table_array(source, "files", "$context.source")
    length(source_entries) == length(_T13_SOURCE_CLOSURE) ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "activation source closure cardinality differs")
    source_files = NamedTuple[]
    source_paths = Set{String}()
    for entry in source_entries
        entry_table = _t13_exact_keys(entry, ("path", "sha256", "mode"),
                                      "$context.source.files[]")
        path = _t13_relative_path(entry_table, "path", "$context.source.files[]")
        path in source_paths &&
            _fail(:BLOCKED, :activation_source_mismatch,
                  "activation source closure contains a duplicate")
        push!(source_paths, path)
        push!(source_files, (path=path,
                             sha256=_t13_hash(entry_table, "sha256",
                                              "$context.source.files[]"),
                             mode=_t13_mode(entry_table, "mode",
                                            "$context.source.files[]")))
    end
    source_paths == Set(String.(collect(_T13_SOURCE_CLOSURE))) ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "activation source closure differs from the consumed include closure")

    command = _t13_exact_keys(_t13_table(document, "command", context),
                              ("body", "environment"), "$context.command")
    body = _t13_string_array(command, "body", "$context.command")
    body == _t13_canonical_command(options) ||
        _fail(:BLOCKED, :activation_command_mismatch,
              "activation command body differs from the approved invocation")
    environment = _t13_table_array(command, "environment", "$context.command")
    length(environment) == 1 ||
        _fail(:BLOCKED, :activation_command_mismatch,
              "activation environment must record exactly one restricted entry")
    environment_entry = _t13_exact_keys(environment[1], ("name", "value"),
                                        "$context.command.environment[]")
    environment_name = _t13_string(environment_entry, "name",
                                   "$context.command.environment[]")
    environment_value = _t13_string(environment_entry, "value",
                                    "$context.command.environment[]")
    environment_name == _T13_SYNTHETIC_ENVIRONMENT &&
        environment_value == "1" &&
        get(ENV, environment_name, "") == environment_value ||
        _fail(:BLOCKED, :activation_command_mismatch,
              "activation environment is not the approved synthetic scope")

    launch = _t13_exact_keys(_t13_table(document, "launch", context),
                             _T13_LAUNCH_KEYS, "$context.launch")
    launch_active_project = _t13_string(launch, "active_project", "$context.launch")
    (isabspath(launch_active_project) &&
     normpath(launch_active_project) == launch_active_project) ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "activation launch active_project is not canonical")
    launch_cwd = _t13_string(launch, "cwd", "$context.launch")
    (isabspath(launch_cwd) && normpath(launch_cwd) == launch_cwd) ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "activation launch cwd is not canonical")
    launch_startup_file = _t13_bool(launch, "startup_file", "$context.launch")
    launch_history_file = _t13_bool(launch, "history_file", "$context.launch")
    launch_threads = _t13_int(launch, "threads", "$context.launch")
    launch_blas_threads = _t13_int(launch, "blas_threads", "$context.launch")
    launch_depot_path = _t13_string(launch, "depot_path", "$context.launch")
    launch_load_path = _t13_string_array(launch, "load_path", "$context.launch")
    launch_offline = _t13_bool(launch, "offline", "$context.launch")
    launch_forbidden = _t13_string_array(launch, "forbidden_environment",
                                         "$context.launch")
    Set(launch_forbidden) == Set(String.(collect(_T13_FORBIDDEN_ENVIRONMENT))) &&
        length(launch_forbidden) == length(_T13_FORBIDDEN_ENVIRONMENT) ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "activation launch forbidden_environment differs from the fixed set")

    destination = _t13_exact_keys(_t13_table(document, "destination", context),
                                  ("path", "reports"), "$context.destination")
    destination_path = _t13_relative_path(destination, "path", "$context.destination")
    reports = _t13_string_array(destination, "reports", "$context.destination")
    Set(reports) == Set(String.(collect(_FILES))) &&
        length(reports) == length(_FILES) ||
        _fail(:BLOCKED, :activation_destination_mismatch,
              "activation report name set differs from the nine frozen report names")
    destination_absolute = _publication_destination(
        root, options["--out-dir"], "activation destination")
    normpath(joinpath(root, destination_path)) == destination_absolute ||
        _fail(:BLOCKED, :activation_destination_mismatch,
              "activation destination differs from the approved command")

    inputs = _t13_exact_keys(_t13_table(document, "inputs", context),
                             ("files", "directories"), "$context.inputs")
    input_entries = _t13_table_array(inputs, "files", "$context.inputs")
    length(input_entries) == 5 ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation input file scope cardinality differs")
    spec_files = NamedTuple[]
    input_paths = Set{String}()
    for entry in input_entries
        entry_table = _t13_exact_keys(entry, ("path", "sha256", "mode"),
                                      "$context.inputs.files[]")
        path = _t13_relative_path(entry_table, "path", "$context.inputs.files[]")
        path in input_paths &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation input file scope contains a duplicate")
        push!(input_paths, path)
        push!(spec_files, (path=path,
                           sha256=_t13_hash(entry_table, "sha256",
                                            "$context.inputs.files[]"),
                           mode=_t13_mode(entry_table, "mode",
                                          "$context.inputs.files[]")))
    end
    input_paths == Set(_t13_input_relative_paths(root, options)) ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation input file scope differs from the approved command")
    directory_entries = _t13_table_array(inputs, "directories", "$context.inputs")
    length(directory_entries) == 3 ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation input directory scope cardinality differs")
    spec_directories = NamedTuple[]
    directory_paths = Set{String}()
    for entry in directory_entries
        entry_table = _t13_exact_keys(entry, ("path", "mode", "members"),
                                      "$context.inputs.directories[]")
        path = _t13_relative_path(entry_table, "path",
                                  "$context.inputs.directories[]")
        path in directory_paths &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation input directory scope contains a duplicate")
        push!(directory_paths, path)
        directory_mode = _t13_mode(entry_table, "mode",
                                   "$context.inputs.directories[]")
        member_entries = _t13_table_array(entry_table, "members",
                                          "$context.inputs.directories[]")
        members = NamedTuple[]
        member_names = Set{String}()
        for member_entry in member_entries
            member_table = _t13_exact_keys(
                member_entry, ("name", "sha256", "mode"),
                "$context.inputs.directories[].members[]")
            name = _t13_string(member_table, "name",
                               "$context.inputs.directories[].members[]")
            (!isempty(name) && name == basename(name) &&
             !occursin('\\', name) && name ∉ (".", "..")) ||
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "activation directory inventory contains a non-member name")
            name in member_names &&
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "activation directory inventory contains a duplicate member")
            push!(member_names, name)
            push!(members, (name=name,
                            sha256=_t13_hash(member_table, "sha256",
                                             "$context.inputs.directories[].members[]"),
                            mode=_t13_mode(member_table, "mode",
                                           "$context.inputs.directories[].members[]")))
        end
        push!(spec_directories, (path=path, mode=directory_mode, members=members))
    end
    directory_paths == Set(_t13_directory_relative_paths(root, options)) ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation input directory scope differs from the approved command")

    dependencies = _t13_exact_keys(_t13_table(document, "dependencies", context),
                                   ("files", "directories", "identity_paths"),
                                   "$context.dependencies")
    dependency_files = NamedTuple[]
    dependency_file_paths = Set{String}()
    for entry in _t13_table_array(dependencies, "files", "$context.dependencies")
        entry_table = _t13_exact_keys(entry, ("path", "sha256", "mode"),
                                      "$context.dependencies.files[]")
        path = _t13_relative_path(entry_table, "path",
                                  "$context.dependencies.files[]")
        path in dependency_file_paths &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation dependency file set contains a duplicate")
        push!(dependency_file_paths, path)
        push!(dependency_files, (path=path,
                                 sha256=_t13_hash(entry_table, "sha256",
                                                  "$context.dependencies.files[]"),
                                 mode=_t13_mode(entry_table, "mode",
                                                "$context.dependencies.files[]")))
    end
    dependency_directories = NamedTuple[]
    dependency_directory_paths = Set{String}()
    for entry in _t13_table_array(dependencies, "directories", "$context.dependencies")
        entry_table = _t13_exact_keys(entry, ("path", "mode", "members"),
                                      "$context.dependencies.directories[]")
        path = _t13_relative_path(entry_table, "path",
                                  "$context.dependencies.directories[]")
        path in dependency_directory_paths &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation dependency directory set contains a duplicate")
        push!(dependency_directory_paths, path)
        directory_mode = _t13_mode(entry_table, "mode",
                                   "$context.dependencies.directories[]")
        members = NamedTuple[]
        member_names = Set{String}()
        for member_entry in _t13_table_array(
                entry_table, "members", "$context.dependencies.directories[]")
            member_table = _t13_exact_keys(
                member_entry, ("name", "sha256", "mode"),
                "$context.dependencies.directories[].members[]")
            name = _t13_string(member_table, "name",
                               "$context.dependencies.directories[].members[]")
            (!isempty(name) && name == basename(name) &&
             !occursin('\\', name) && name ∉ (".", "..")) ||
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "activation dependency directory has a non-member name")
            name in member_names &&
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "activation dependency directory has a duplicate member")
            push!(member_names, name)
            push!(members, (name=name,
                            sha256=_t13_hash(member_table, "sha256",
                                             "$context.dependencies.directories[].members[]"),
                            mode=_t13_mode(member_table, "mode",
                                           "$context.dependencies.directories[].members[]")))
        end
        push!(dependency_directories, (path=path, mode=directory_mode,
                                       members=members))
    end
    identity_paths = NamedTuple[]
    identity_seen = Set{String}()
    for entry in _t13_table_array(dependencies, "identity_paths",
                                  "$context.dependencies")
        entry_table = _t13_exact_keys(entry, ("path", "kind"),
                                      "$context.dependencies.identity_paths[]")
        path = _t13_relative_path(entry_table, "path",
                                  "$context.dependencies.identity_paths[]")
        kind = _t13_string(entry_table, "kind",
                           "$context.dependencies.identity_paths[]")
        kind in ("file", "directory") ||
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation identity path kind is invalid")
        path in identity_seen &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation identity path set contains a duplicate")
        push!(identity_seen, path)
        push!(identity_paths, (path=path, kind=kind))
    end

    predecessors = _t13_exact_keys(
        _t13_table(document, "predecessors", context),
        ("gate4_expected_boulder_postimage_path",
         "gate4_expected_boulder_postimage_sha256",
         "gate4_transition_receipt_path", "gate4_transition_receipt_sha256",
         "gate4_post_review_path", "gate4_post_review_sha256",
         "gate4_oracle_final_review_path", "gate4_oracle_final_review_sha256",
         "gate5_publication_receipt_path", "gate5_publication_receipt_sha256",
         "gate2_oracle_review_path", "gate2_oracle_review_sha256"),
        "$context.predecessors")
    predecessor_bindings = Dict{String,NamedTuple}()
    for key in ("gate4_expected_boulder_postimage", "gate4_transition_receipt",
                "gate4_post_review", "gate4_oracle_final_review")
        predecessor_bindings[key] = _t13_authority_binding(key)
    end
    predecessor_bindings["gate5_publication_receipt"] =
        (path=_T13_GATE5_PUBLICATION_PATH, sha256=_T13_GATE5_PUBLICATION_SHA256)
    predecessor_bindings["gate2_oracle_review"] =
        (path=_T13_GATE2_ORACLE_REVIEW_PATH, sha256=_T13_GATE2_ORACLE_REVIEW_SHA256)
    for (key, binding) in predecessor_bindings
        _t13_check_predecessor(predecessors, key, binding.path, binding.sha256,
                               "$context.predecessors")
    end
    gate4_expected = predecessor_bindings["gate4_expected_boulder_postimage"]
    gate4_transition = predecessor_bindings["gate4_transition_receipt"]
    gate4_post_review = predecessor_bindings["gate4_post_review"]
    gate4_oracle = predecessor_bindings["gate4_oracle_final_review"]
    gate5_publication = predecessor_bindings["gate5_publication_receipt"]
    gate2_oracle = predecessor_bindings["gate2_oracle_review"]

    permissions = _t13_exact_keys(
        _t13_table(document, "permissions", context),
        ("real_campaign", "public_producer", "t14", "downstream", "grading",
         "runtime_equivalence"), "$context.permissions")
    for key in ("real_campaign", "public_producer", "t14", "downstream", "grading",
                "runtime_equivalence")
        _t13_false(permissions, key, "$context.permissions")
    end

    return (scope=scope,
            control=(role=control_role, root=control_root,
                     root_mode=control_root_mode,
                     parent_mode=control_parent_mode,
                     file_mode=control_file_mode),
            entrypoint_path=entrypoint_path,
            entrypoint_sha256=entrypoint_sha256,
            project_path=project_path,
            project_sha256=project_sha256,
            manifest_path=manifest_path,
            manifest_sha256=manifest_sha256,
            destination_path=destination_path,
            command_body=body,
            source_files=source_files,
            spec_files=spec_files,
            spec_directories=spec_directories,
            dependency_files=dependency_files,
            dependency_directories=dependency_directories,
            identity_paths=identity_paths,
            launch=(active_project=launch_active_project,
                    cwd=launch_cwd,
                    startup_file=launch_startup_file,
                    history_file=launch_history_file,
                    threads=launch_threads,
                    blas_threads=launch_blas_threads,
                    depot_path=launch_depot_path,
                    load_path=launch_load_path,
                    offline=launch_offline,
                    forbidden_environment=launch_forbidden),
            gate4_expected_boulder_postimage_path=gate4_expected.path,
            gate4_expected_boulder_postimage_sha256=gate4_expected.sha256,
            gate4_transition_receipt_path=gate4_transition.path,
            gate4_transition_receipt_sha256=gate4_transition.sha256,
            gate4_post_review_path=gate4_post_review.path,
            gate4_post_review_sha256=gate4_post_review.sha256,
            gate4_oracle_final_review_path=gate4_oracle.path,
            gate4_oracle_final_review_sha256=gate4_oracle.sha256,
            gate5_publication_receipt_path=gate5_publication.path,
            gate5_publication_receipt_sha256=gate5_publication.sha256,
            gate2_oracle_review_path=gate2_oracle.path,
            gate2_oracle_review_sha256=gate2_oracle.sha256)
end

function _validate_t13_launch_environment(launch, root::String)
    expected_project = normpath(joinpath(root, "Project.toml"))
    active = Base.active_project()
    active === nothing &&
        _fail(:BLOCKED, :activation_launch_mismatch, "no active project is set")
    realpath(active) == expected_project ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "active project is not the bound-root Project.toml")
    launch.active_project == expected_project ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "recorded active project differs from the bound root")
    launch.cwd == pwd() ||
        _fail(:BLOCKED, :activation_launch_mismatch, "recorded cwd differs")
    launch.startup_file == (Base.JLOptions().startupfile == 2) ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "startup-file setting differs from the recorded binding")
    launch.history_file == (Base.JLOptions().historyfile == 0) ||
        _fail(:BLOCKED, :activation_launch_mismatch,
              "history-file setting differs from the recorded binding")
    launch.threads == Base.Threads.nthreads() ||
        _fail(:BLOCKED, :activation_launch_mismatch, "Julia thread count differs")
    launch.blas_threads == LinearAlgebra.BLAS.get_num_threads() ||
        _fail(:BLOCKED, :activation_launch_mismatch, "BLAS thread count differs")
    launch.depot_path == first(DEPOT_PATH) ||
        _fail(:BLOCKED, :activation_launch_mismatch, "depot path differs")
    launch.load_path == String.(LOAD_PATH) ||
        _fail(:BLOCKED, :activation_launch_mismatch, "load path differs")
    launch.offline == (get(ENV, "JULIA_PKG_OFFLINE", "") in ("true", "1")) ||
        _fail(:BLOCKED, :activation_launch_mismatch, "offline setting differs")
    for name in _T13_FORBIDDEN_ENVIRONMENT
        haskey(ENV, name) && !isempty(ENV[name]) &&
            _fail(:BLOCKED, :activation_launch_mismatch,
                  "forbidden environment variable is set: $name")
    end
    return nothing
end

function _t13_snapshot_at(snapshots::Vector{_Snapshot}, root::String,
                          relative::String, context::String)
    absolute = normpath(joinpath(root, relative))
    matches = [snapshot for snapshot in snapshots if snapshot.path == absolute]
    length(matches) == 1 ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "$context snapshot is absent or duplicated")
    return only(matches)
end

function _t13_receipt_string(document, key::String, context::String)
    haskey(document, key) ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "$context.$key is absent")
    value = document[key]
    value isa AbstractString ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "$context.$key is not a string")
    return String(value)
end

function _t13_dependency_hash(entries, path::String)
    matches = [entry.sha256 for entry in entries if entry.path == path]
    length(matches) == 1 ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "declared dependency is absent or duplicated: $path")
    return only(matches)
end

function _t13_snapshot_approved_directory(root::String, entry, seen,
                                          modes::Dict{String,UInt};
                                          existing::Union{Nothing,Vector{_Snapshot}}=nothing)
    directory = _repo_path(root, entry.path, "activation input directory")
    isdir(directory) && !islink(directory) ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation input directory is absent or not a directory")
    _runtime_v3_mode(directory, entry.mode,
                     "activation input directory $(entry.path)")
    actual_names = sort(readdir(directory))
    approved_names = sort(String[member.name for member in entry.members])
    actual_names == approved_names ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation input directory membership differs: $(entry.path)")
    snapshots = _Snapshot[]
    for member in sort(entry.members; by=member -> member.name)
        member_path = joinpath(directory, member.name)
        islink(member_path) &&
            _fail(:BLOCKED, :activation_symlink,
                  "activation input member is a symlink: $(entry.path)/$(member.name)")
        absolute = normpath(joinpath(root, entry.path, member.name))
        overlap = existing === nothing ? _Snapshot[] :
                  [snapshot for snapshot in existing if snapshot.path == absolute]
        if !isempty(overlap)
            length(overlap) == 1 ||
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "overlapping directory member is duplicated: $(entry.path)/$(member.name)")
            only(overlap).sha256 == member.sha256 ||
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "overlapping directory member hash differs: $(entry.path)/$(member.name)")
            _runtime_v3_mode(absolute, member.mode,
                             "activation input member $(entry.path)/$(member.name)")
            modes[absolute] = member.mode
            continue
        end
        snapshot = _snapshot(root, joinpath(entry.path, member.name), member.sha256,
                             "activation input member $(entry.path)/$(member.name)",
                             seen)
        _runtime_v3_mode(snapshot.path, member.mode,
                         "activation input member $(entry.path)/$(member.name)")
        modes[snapshot.path] = member.mode
        push!(snapshots, snapshot)
    end
    modes[directory] = entry.mode
    return directory, actual_names, snapshots
end

function _authorize_t13_activation(options::Dict{String,String})
    root = options["--root"]
    isabspath(root) || _fail(:BLOCKED, :authority_path_invalid, "root must be absolute")
    isdir(root) || _fail(:BLOCKED, :authority_path_invalid, "root is not a directory")
    root = normpath(root)
    realpath(root) == root ||
        _fail(:BLOCKED, :authority_path_invalid, "root is not canonical")
    VERSION == _RUNTIME_V4_JULIA ||
        _fail(:BLOCKED, :runtime_version_mismatch, "Julia 1.12.7 is required")

    expected_receipt = String(options["--activation-sha256"])
    occursin(r"^[0-9a-f]{64}$", expected_receipt) ||
        _fail(:BLOCKED, :activation_receipt_mismatch,
              "activation receipt digest is not a SHA-256")
    receipt_rel = _t13_relative_path_value(options["--activation-receipt"],
                                           "activation receipt")
    seen = Dict{Tuple{UInt64,UInt64},String}()
    receipt = _snapshot(root, receipt_rel, expected_receipt, "activation receipt", seen)
    receipt_document = _t13_parse_toml(receipt, "activation receipt")
    receipt_bindings = _validate_t13_activation_receipt_document(receipt_document,
                                                                 receipt_rel)
    spec = _snapshot(root, receipt_bindings.spec_path, receipt_bindings.spec_sha256,
                     "activation execution spec", seen)
    spec_document = _t13_parse_toml(spec, "activation execution spec")
    spec_bindings = _validate_t13_execution_spec_document(spec_document, root, options)
    _validate_t13_launch_environment(spec_bindings.launch, root)

    entrypoint_absolute = _repo_path(root, spec_bindings.entrypoint_path,
                                     "activation entrypoint")
    running_entrypoint = realpath(@__FILE__)
    running_entrypoint == entrypoint_absolute ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "running evaluator entrypoint does not belong to the bound root")
    isfile(running_entrypoint) && !islink(running_entrypoint) ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "running evaluator entrypoint is not a regular file")
    _hash_bytes(read(running_entrypoint)) == spec_bindings.entrypoint_sha256 ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "running evaluator entrypoint bytes differ from the approved source")

    # The receipt's control root must agree with the spec's Main-control root.
    # The receipt preimage must be the fixed Gate4 state; the receipt postimage
    # is the reviewed digest that the live Main file must match.  The dedicated
    # control reader is the only place the fixed Main filename is read.
    receipt_bindings.control_root == spec_bindings.control.root ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation receipt control_root differs from the spec control root")
    receipt_bindings.boulder_preimage_sha256 == _RUNTIME_V4_BOULDER_SHA256 ||
        _fail(:BLOCKED, :activation_state_mismatch,
              "activation Boulder preimage is not the fixed Gate4 state")
    control = _read_main_control(spec_bindings.control,
                                 receipt_bindings.control_root,
                                 receipt_bindings.boulder_preimage_sha256,
                                 receipt_bindings.boulder_postimage_sha256)

    _runtime_v3_mode(receipt.path, _T13_CONTROL_RECORD_MODE,
                     "activation receipt")
    _runtime_v3_mode(spec.path, _T13_CONTROL_RECORD_MODE,
                     "activation execution spec")
    source_review = _snapshot(root, receipt_bindings.source_review_path,
                              receipt_bindings.source_review_sha256,
                              "activation source review", seen)
    _runtime_v3_mode(source_review.path, _T13_CONTROL_RECORD_MODE,
                     "activation source review")
    transition_review = _snapshot(root, receipt_bindings.transition_review_path,
                                  receipt_bindings.transition_review_sha256,
                                  "activation transition review", seen)
    _runtime_v3_mode(transition_review.path, _T13_CONTROL_RECORD_MODE,
                     "activation transition review")
    publication = _snapshot(root, receipt_bindings.publication_path,
                            receipt_bindings.publication_sha256,
                            "activation transition publication", seen)
    _runtime_v3_mode(publication.path, _T13_CONTROL_RECORD_MODE,
                     "activation transition publication")
    # The Main-control object above already bound the one fixed
    # `.omo/boulder.json` under the reviewed Main root.  The execution root R is
    # never consulted for control bytes; an R-local file is a non-consumed
    # decoy, and no R Boulder snapshot is created or aliased here.
    predecessors = _Snapshot[]
    for (path, digest, mode, label) in (
        (spec_bindings.gate4_expected_boulder_postimage_path,
         spec_bindings.gate4_expected_boulder_postimage_sha256,
         _runtime_v4_expected_mode(
             joinpath(root, spec_bindings.gate4_expected_boulder_postimage_path)),
         "Gate4 expected Boulder postimage"),
        (spec_bindings.gate4_transition_receipt_path,
         spec_bindings.gate4_transition_receipt_sha256,
         _runtime_v4_expected_mode(
             joinpath(root, spec_bindings.gate4_transition_receipt_path)),
         "Gate4 transition receipt"),
        (spec_bindings.gate4_post_review_path,
         spec_bindings.gate4_post_review_sha256,
         _runtime_v4_expected_mode(
             joinpath(root, spec_bindings.gate4_post_review_path)),
         "Gate4 post-execution review"),
        (spec_bindings.gate4_oracle_final_review_path,
         spec_bindings.gate4_oracle_final_review_sha256,
         _runtime_v4_expected_mode(
             joinpath(root, spec_bindings.gate4_oracle_final_review_path)),
         "Gate4 Oracle final review"),
        (spec_bindings.gate5_publication_receipt_path,
         spec_bindings.gate5_publication_receipt_sha256,
         _T13_GATE5_PUBLICATION_MODE,
         "Gate5 publication receipt"),
        (spec_bindings.gate2_oracle_review_path,
         spec_bindings.gate2_oracle_review_sha256,
         _T13_GATE2_ORACLE_REVIEW_MODE,
         "Gate2 Oracle review"),
    )
        mode === nothing &&
            _fail(:BLOCKED, :activation_predecessor_mismatch,
                  "$label has no fixed predecessor mode")
        snapshot = _snapshot(root, path, digest, label, seen)
        _runtime_v3_mode(snapshot.path, mode, label)
        push!(predecessors, snapshot)
    end
    project = _snapshot(root, spec_bindings.project_path, spec_bindings.project_sha256,
                        "activation Project.toml", seen)
    manifest = _snapshot(root, spec_bindings.manifest_path,
                         spec_bindings.manifest_sha256,
                         "activation Manifest.toml", seen)
    source_snapshots = _Snapshot[]
    for entry in spec_bindings.source_files
        snapshot = _snapshot(root, entry.path, entry.sha256,
                             "activation source $(entry.path)", seen)
        _runtime_v3_mode(snapshot.path, entry.mode, "activation source $(entry.path)")
        push!(source_snapshots, snapshot)
    end
    entrypoint_matches = [snapshot for snapshot in source_snapshots
                          if snapshot.path == entrypoint_absolute]
    length(entrypoint_matches) == 1 ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "activation source closure does not bind the entrypoint exactly once")
    only(entrypoint_matches).sha256 == spec_bindings.entrypoint_sha256 ||
        _fail(:BLOCKED, :activation_source_mismatch,
              "activation entrypoint binding differs from the source closure")

    snapshots = vcat(_Snapshot[receipt, spec, source_review, transition_review,
                               publication],
                     predecessors, _Snapshot[project, manifest], source_snapshots)
    modes = Dict{String,UInt}()
    for snapshot in snapshots
        modes[snapshot.path] = UInt(stat(snapshot.path).mode) & UInt(0o777)
    end
    for entry in spec_bindings.source_files
        modes[normpath(joinpath(root, entry.path))] = entry.mode
    end
    # Bind every approved input byte and mode before any publication-capable
    # context exists.  The loader must not re-hash these inputs with observed
    # values; it revalidates these snapshots instead.
    input_snapshots = _Snapshot[]
    input_inventories = Tuple{String,Vector{String}}[]
    input_directories = String[]
    for entry in spec_bindings.spec_files
        snapshot = _snapshot(root, entry.path, entry.sha256,
                             "activation input file $(entry.path)", seen)
        _runtime_v3_mode(snapshot.path, entry.mode,
                         "activation input file $(entry.path)")
        modes[snapshot.path] = entry.mode
        push!(input_snapshots, snapshot)
    end
    for entry in spec_bindings.spec_directories
        directory, names, members = _t13_snapshot_approved_directory(
            root, entry, seen, modes)
        push!(input_directories, directory)
        push!(input_inventories, (relpath(directory, root), names))
        append!(input_snapshots, members)
    end
    directory_identities = Dict{String,Tuple{UInt64,UInt64}}()
    # F3: finite receipt-derived and fixed provenance dependencies.  The
    # expected closed set is derived from the bound producer receipts and the
    # unchanged predecessor readers; the spec must declare exactly that set.
    forward_rel = _t13_relative_from_option(root, options["--forward-receipt"])
    backward_rel = _t13_relative_from_option(root, options["--backward-receipt"])
    forward_snapshot = _t13_snapshot_at(input_snapshots, root, forward_rel,
                                        "forward producer receipt")
    backward_snapshot = _t13_snapshot_at(input_snapshots, root, backward_rel,
                                         "backward producer receipt")
    forward_receipt = _t13_parse_toml(forward_snapshot, "forward producer receipt")
    backward_receipt = _t13_parse_toml(backward_snapshot, "backward producer receipt")
    for (document, label) in ((forward_receipt, "forward"),
                              (backward_receipt, "backward"))
        _t13_receipt_string(document, "root", "$label producer receipt") == root ||
            _fail(:BLOCKED, :activation_dependency_mismatch,
                  "$label producer receipt root differs")
        _t13_receipt_string(document, "cwd", "$label producer receipt") == root ||
            _fail(:BLOCKED, :activation_dependency_mismatch,
                  "$label producer receipt cwd differs")
    end
    features_rel = _t13_relative_from_option(
        root, _t13_receipt_string(forward_receipt, "features_path",
                                  "forward producer receipt"))
    features_rel == _t13_relative_from_option(
        root, _t13_receipt_string(backward_receipt, "features_path",
                                  "backward producer receipt")) ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "producer feature table paths differ")
    forward_output_rel = _t13_relative_from_option(
        root, _t13_receipt_string(forward_receipt, "output_path",
                                  "forward producer receipt"))
    backward_output_rel = _t13_relative_from_option(
        root, _t13_receipt_string(backward_receipt, "output_path",
                                  "backward producer receipt"))
    data_path_rel = _t13_relative_from_option(
        root, _t13_receipt_string(forward_receipt, "data_path",
                                  "forward producer receipt"))
    data_path_rel == _t13_relative_from_option(
        root, _t13_receipt_string(backward_receipt, "data_path",
                                  "backward producer receipt")) ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "producer data paths differ")
    for (document, script, digest, label) in (
            (forward_receipt, _T13_PRODUCER_SCRIPTS[1].path,
             _T13_PRODUCER_SCRIPTS[1].sha256, "forward"),
            (backward_receipt, _T13_PRODUCER_SCRIPTS[2].path,
             _T13_PRODUCER_SCRIPTS[2].sha256, "backward"))
        _t13_relative_from_option(
            root, _t13_receipt_string(document, "source_path",
                                      "$label producer receipt")) == script ||
            _fail(:BLOCKED, :activation_dependency_mismatch,
                  "$label producer source_path differs")
        _t13_receipt_string(document, "source_sha256",
                            "$label producer receipt") == digest ||
            _fail(:BLOCKED, :activation_dependency_mismatch,
                  "$label producer source hash differs")
    end
    expected_files = Set{String}(String.(collect(_T13_FIXED_PROVENANCE_FILES)))
    push!(expected_files, features_rel)
    push!(expected_files, forward_output_rel)
    push!(expected_files, backward_output_rel)
    declared_files = Set(String(entry.path)
                         for entry in spec_bindings.dependency_files)
    declared_files == expected_files ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "declared dependency file set differs from the derived closed set")
    declared_directories = Set(String(entry.path)
                               for entry in spec_bindings.dependency_directories)
    declared_directories == Set([_T13_HIERARCHICAL_DIRECTORY]) ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "declared dependency directory set differs from the derived set")
    expected_identity = Set([(_T13_PRODUCER_SCRIPTS[1].path, "file"),
                             (_T13_PRODUCER_SCRIPTS[2].path, "file"),
                             (data_path_rel, "directory")])
    declared_identity = Set([(String(entry.path), String(entry.kind))
                             for entry in spec_bindings.identity_paths])
    declared_identity == expected_identity ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "declared identity path set differs from the derived closed set")
    universe_receipt = _t13_snapshot_at(
        input_snapshots, root,
        joinpath(_t13_relative_from_option(root, options["--universe-dir"]),
                 "receipt.toml"), "universe receipt")
    edge_receipt = _t13_snapshot_at(
        input_snapshots, root,
        joinpath(_t13_relative_from_option(root, options["--edge-dir"]),
                 "receipt.toml"), "edge receipt")
    universe_document = _t13_parse_toml(universe_receipt, "universe receipt")
    edge_document = _t13_parse_toml(edge_receipt, "edge receipt")
    features_hash = _t13_dependency_hash(spec_bindings.dependency_files,
                                         features_rel)
    features_hash == _t13_receipt_string(forward_receipt, "feature_sha256",
                                         "forward producer receipt") ==
        _t13_receipt_string(backward_receipt, "feature_sha256",
                            "backward producer receipt") ==
        _t13_receipt_string(universe_document, "feature_sha256",
                            "universe receipt") ==
        _t13_receipt_string(edge_document, "feature_sha256", "edge receipt") ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "feature table hash cross-binding differs")
    _t13_dependency_hash(spec_bindings.dependency_files, forward_output_rel) ==
        _t13_receipt_string(edge_document, "forward_patch_sha256",
                            "edge receipt") ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "forward patch hash cross-binding differs")
    _t13_dependency_hash(spec_bindings.dependency_files, backward_output_rel) ==
        _t13_receipt_string(edge_document, "backward_patch_sha256",
                            "edge receipt") ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "backward patch hash cross-binding differs")
    forward_snapshot.sha256 ==
        _t13_receipt_string(edge_document, "forward_patch_receipt_sha256",
                            "edge receipt") ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "forward producer receipt cross-binding differs")
    backward_snapshot.sha256 ==
        _t13_receipt_string(edge_document, "backward_patch_receipt_sha256",
                            "edge receipt") ||
        _fail(:BLOCKED, :activation_dependency_mismatch,
              "backward producer receipt cross-binding differs")
    for entry in spec_bindings.dependency_files
        absolute = normpath(joinpath(root, entry.path))
        overlap = [snapshot for snapshot in input_snapshots
                   if snapshot.path == absolute]
        if !isempty(overlap)
            length(overlap) == 1 ||
                _fail(:BLOCKED, :activation_dependency_mismatch,
                      "overlapping dependency is duplicated: $(entry.path)")
            only(overlap).sha256 == entry.sha256 ||
                _fail(:BLOCKED, :activation_dependency_mismatch,
                      "overlapping dependency hash differs: $(entry.path)")
            _runtime_v3_mode(absolute, entry.mode,
                             "activation dependency $(entry.path)")
            modes[absolute] = entry.mode
            continue
        end
        snapshot = _snapshot(root, entry.path, entry.sha256,
                             "activation dependency $(entry.path)", seen)
        _runtime_v3_mode(snapshot.path, entry.mode,
                         "activation dependency $(entry.path)")
        modes[snapshot.path] = entry.mode
        push!(input_snapshots, snapshot)
    end
    for entry in spec_bindings.dependency_directories
        directory, names, members = _t13_snapshot_approved_directory(
            root, entry, seen, modes;
            existing=vcat(snapshots, input_snapshots))
        directory in input_directories || push!(input_directories, directory)
        push!(input_inventories, (relpath(directory, root), names))
        append!(input_snapshots, members)
    end
    for entry in spec_bindings.identity_paths
        absolute = _repo_path(root, entry.path, "activation identity path")
        if entry.kind == "file"
            isfile(absolute) && !islink(absolute) ||
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "activation identity file is not regular: $(entry.path)")
        else
            isdir(absolute) && !islink(absolute) ||
                _fail(:BLOCKED, :activation_binding_mismatch,
                      "activation identity directory is absent: $(entry.path)")
        end
        info = stat(absolute)
        modes[absolute] = UInt(info.mode) & UInt(0o777)
        directory_identities[absolute] = (UInt64(info.device), UInt64(info.inode))
    end
    snapshots = vcat(snapshots, input_snapshots)
    # The dedicated Main-control file must not double as an approved ordinary
    # input (including the M == R case where the R-local control file is also a
    # declared input).  Reject the overlap before any authorization is usable.
    _reject_main_control_overlap(snapshots, control, "activation authorization")
    # Retain captured directory identities (root, destination parent, approved
    # input directories) so a replacement directory is rejected even when its
    # members are hardlinks to the same files.
    destination_absolute = _publication_destination(
        root, options["--out-dir"], "activation destination")
    for directory in vcat(String[root, dirname(destination_absolute)],
                          input_directories)
        isdir(directory) && !islink(directory) ||
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "activation directory is absent or not a directory: $directory")
        info = stat(directory)
        directory_identities[directory] = (UInt64(info.device), UInt64(info.inode))
    end
    bindings = _ActivationBindings(
        spec_bindings.scope,
        receipt.path, receipt.sha256,
        spec.path, spec.sha256,
        source_review.path, source_review.sha256,
        transition_review.path, transition_review.sha256,
        publication.path, publication.sha256,
        receipt_bindings.boulder_preimage_path, receipt_bindings.boulder_preimage_sha256,
        control.file.path, control.postimage_sha256,
        entrypoint_absolute, spec_bindings.entrypoint_sha256,
        spec_bindings.destination_path, modes,
        root, spec_bindings.command_body, directory_identities,
        control)
    # Aggregate revalidation before any authorization is returned: snapshots,
    # modes and captured directory identities must still hold.
    for snapshot in snapshots
        _verify_snapshot(snapshot, "activation authorization snapshot")
    end
    _verify_activation_bindings(bindings)
    return _ActivationAuthorization(root, bindings, snapshots,
                                    input_inventories,
                                    spec_bindings.spec_files,
                                    spec_bindings.spec_directories)
end

function _validate_activation_receiving_options(activation::_ActivationAuthorization,
                                                root::String,
                                                options::Dict{String,String})
    activation.root == root ||
        _fail(:BLOCKED, :activation_binding_mismatch,
              "activation root differs from the receiving root")
    _t13_canonical_command(options) == activation.bindings.command_body ||
        _fail(:BLOCKED, :activation_command_mismatch,
              "activation command differs from the receiving options")
    destination = _publication_destination(root, options["--out-dir"],
                                           "activation destination")
    normpath(joinpath(root, activation.bindings.destination)) == destination ||
        _fail(:BLOCKED, :activation_destination_mismatch,
              "activation destination differs from the receiving options")
    return nothing
end

function _assemble_production_context(root::String, options::Dict{String,String},
                                      authority;
                                      activation::Union{Nothing,_ActivationAuthorization}=nothing)::_ProductionContext
    if activation === nothing
        input_snapshots, inventories = _snapshot_context(root, options, authority)
        context = _ProductionContext(root, vcat(authority.snapshots, input_snapshots),
                                     vcat(authority.inventories, inventories), "")
    else
        context = _ProductionContext(
            root,
            _merge_snapshots(authority.snapshots, activation.snapshots,
                             "activation binding"),
            vcat(authority.inventories, activation.inventories), "",
            activation.bindings)
    end
    # The context is only returned after full revalidation; callers must never
    # treat an assembled-but-unverified context as publication-capable.
    _verify_context(context)
    return context
end

function _t13_activation_receipt_lines(io::IO, activation::_ActivationBindings)
    println(io, "activation_schema = \"", _T13_ACTIVATION_SCHEMA, "\"")
    println(io, "activation_authority_schema = ",
            repr(_T13_ACTIVATION_AUTHORITY_SCHEMA))
    println(io, "activation_scope = ", repr(activation.scope))
    println(io, "activation_receipt_sha256 = ", repr(activation.receipt_sha256))
    println(io, "activation_spec_sha256 = ", repr(activation.spec_sha256))
    println(io, "activation_source_review_sha256 = ",
            repr(activation.source_review_sha256))
    println(io, "activation_transition_review_sha256 = ",
            repr(activation.transition_review_sha256))
    println(io, "activation_boulder_preimage_sha256 = ",
            repr(activation.boulder_preimage_sha256))
    println(io, "activation_boulder_postimage_sha256 = ",
            repr(activation.boulder_postimage_sha256))
    # Main-control provenance: actual role/root/fixed file/digest/modes.  This
    # records which Main root and file the run consumed; it is integrity/scope
    # provenance, not signer or process attestation.
    println(io, "activation_control_role = ", repr(_T13_CONTROL_ROLE))
    println(io, "activation_control_root = ", repr(activation.control.root))
    println(io, "activation_control_file = ", repr(_T13_CONTROL_FILENAME))
    println(io, "activation_control_sha256 = ",
            repr(activation.control.file.sha256))
    println(io, "activation_control_root_mode = ",
            repr(_t13_mode_text(activation.control.root_mode)))
    println(io, "activation_control_parent_mode = ",
            repr(_t13_mode_text(activation.control.parent_mode)))
    println(io, "activation_control_file_mode = ",
            repr(_t13_mode_text(activation.control.file_mode)))
    println(io, "activation_publication_sha256 = ",
            repr(activation.publication_sha256))
    println(io, "activation_entrypoint_sha256 = ",
            repr(activation.entrypoint_sha256))
    println(io, "activation_destination = ", repr(activation.destination))
    return nothing
end

function _runtime_contract(kind::Symbol)
    kind == :runtime_v4 && return (
        generation=4,
        julia="1.12.7",
        executable_sha256=_RUNTIME_V4_EXECUTABLE_SHA256,
        sysimage_sha256=_RUNTIME_V4_SYSIMAGE_SHA256,
        boulder_sha256=_RUNTIME_V4_BOULDER_SHA256,
        transition_receipt_sha256=
            "87818b63b132d8856dec22cb389cc936283bb0827cfa4884cbfd6e8c169de092",
        post_review_sha256=
            "80a133ba8b13d8ac29b1f8d9dcbb07d5c10ee19d47dc60e5105d0897b8032e13",
        oracle_final_review_sha256=
            "cab6c866c256ce160c13b4694531f661b2db9cbfac5dd47b34cf69b073635fbd",
    )
    kind == :runtime && return (
        edge_model_sha256=_RUNTIME_T12_EDGE_MODEL_SHA256,
        chain_sha256="cdc0c788d49298f50721b091535a238cd552454d9885404456c593e4625ea39f",
        edge_model_product_path="products/edge_model.corrected.jl",
        corrected_test_sha256=_RUNTIME_T12_CORRECTED_TEST_SHA256,
        root_manifest_path=_RUNTIME_T12_ROOT_MANIFEST_PATH,
        root_manifest_sha256=_RUNTIME_T12_ROOT_MANIFEST_SHA256,
        claim_path=_RUNTIME_T12_CLAIM_PATH,
        claim_sha256=_RUNTIME_T12_CLAIM_SHA256,
        review_path=_RUNTIME_T12_REVIEW_PATH,
        review_sha256=_RUNTIME_T12_REVIEW_SHA256,
        publication_path=_RUNTIME_T12_PUBLICATION_PATH,
        publication_sha256=_RUNTIME_T12_PUBLICATION_SHA256,
        manifest_entries=33,
        claim_schema="stmfit-t12-t10-source-bundle-authority-correction-done-claim-v1",
        review_schema="stmfit-t12-t10-source-bundle-authority-correction-review-v1",
        publication_schema="stmfit-t12-authority-correction-publication-receipt-v1",
    )
    kind == :runtime_v2 && return (
        edge_model_sha256=_RUNTIME_V2_T12_EDGE_MODEL_SHA256,
        chain_sha256="cdc0c788d49298f50721b091535a238cd552454d9885404456c593e4625ea39f",
        edge_model_product_path="products/edge_model.corrected-v2.jl",
        corrected_test_sha256=_RUNTIME_V2_T12_CORRECTED_TEST_SHA256,
        root_manifest_path=_RUNTIME_V2_T12_ROOT_MANIFEST_PATH,
        root_manifest_sha256=_RUNTIME_V2_T12_ROOT_MANIFEST_SHA256,
        claim_path=_RUNTIME_V2_T12_CLAIM_PATH,
        claim_sha256=_RUNTIME_V2_T12_CLAIM_SHA256,
        review_path=_RUNTIME_V2_T12_REVIEW_PATH,
        review_sha256=_RUNTIME_V2_T12_REVIEW_SHA256,
        publication_path=_RUNTIME_V2_T12_PUBLICATION_PATH,
        publication_sha256=_RUNTIME_V2_T12_PUBLICATION_SHA256,
        manifest_entries=37,
        claim_schema="stmfit-t12-t11-input-role-authority-correction-done-claim-v2",
        review_schema="stmfit-t12-t11-input-role-authority-adversarial-review-v2",
        publication_schema="stmfit-t12-authority-correction-publication-receipt-v2",
    )
    kind == :runtime_v3 && return (
        edge_model_sha256=_RUNTIME_V3_T12_EDGE_MODEL_SHA256,
        chain_sha256=_RUNTIME_V3_T12_CHAIN_SHA256,
        edge_model_product_path="products/v3/edge_model.jl",
        chain_product_path="products/v3/chain_inference.jl",
        test_product_path="products/v3/test_structured_chain_inference.jl",
        corrected_test_sha256=_RUNTIME_V3_T12_CORRECTED_TEST_SHA256,
        root_manifest_path=_RUNTIME_V3_T12_ROOT_MANIFEST_PATH,
        root_manifest_sha256=_RUNTIME_V3_T12_ROOT_MANIFEST_SHA256,
        claim_path=_RUNTIME_V3_T12_CLAIM_PATH,
        claim_sha256=_RUNTIME_V3_T12_CLAIM_SHA256,
        review_path=_RUNTIME_V3_T12_REVIEW_PATH,
        review_sha256=_RUNTIME_V3_T12_REVIEW_SHA256,
        publication_path=_RUNTIME_V3_T12_PUBLICATION_PATH,
        publication_sha256=_RUNTIME_V3_T12_PUBLICATION_SHA256,
        manifest_entries=_RUNTIME_V3_T12_ROOT_MANIFEST_ENTRIES,
        claim_schema="stmfit-t12-stable-binary-marginal-normalization-done-claim-v3",
        review_schema="stmfit-t12-stable-binary-marginal-normalization-review-v3",
        publication_schema="stmfit-t12-stable-binary-marginal-normalization-publication-receipt-v3",
    )
    _fail(:BLOCKED, :authority_binding_mismatch, "runtime authority generation differs")
end

function _runtime_correction_bindings(contract)
    return (
        ("t12_authority_root_manifest_path", contract.root_manifest_path,
         "t12_authority_root_manifest_sha256", contract.root_manifest_sha256),
        ("t12_authority_claim_path", contract.claim_path,
         "t12_authority_claim_sha256", contract.claim_sha256),
        ("t12_authority_review_path", contract.review_path,
         "t12_authority_review_sha256", contract.review_sha256),
        ("t12_authority_publication_receipt_path", contract.publication_path,
         "t12_authority_publication_receipt_sha256", contract.publication_sha256),
    )
end

function _validate_runtime_authority_config(authority, kind::Symbol)
    kind == :runtime_v4 && return _validate_runtime_v4_authority_config(authority)
    contract = _runtime_contract(kind)
    edge_path, edge_hash = _authority_entry(
        authority, "t12_edge_model_path", "t12_edge_model_sha256")
    edge_path == "test/lib/structured_assignment/edge_model.jl" &&
        edge_hash == contract.edge_model_sha256 ||
        _fail(:BLOCKED, :authority_binding_mismatch,
              "runtime corrected edge-model binding differs")
    if kind == :runtime_v3
        chain_path, chain_hash = _authority_entry(
            authority, "t12_chain_path", "t12_chain_sha256")
        chain_path == "test/lib/structured_assignment/chain_inference.jl" &&
            chain_hash == contract.chain_sha256 ||
            _fail(:BLOCKED, :authority_binding_mismatch,
                  "runtime corrected chain-inference binding differs")
        test_path, test_hash = _authority_entry(
            authority, "t12_test_path", "t12_test_sha256")
        test_path == "test/test_structured_chain_inference.jl" &&
            test_hash == contract.corrected_test_sha256 ||
            _fail(:BLOCKED, :authority_binding_mismatch,
                  "runtime corrected T12 test binding differs")
    end
    expected_keys = Set{String}()
    for (path_key, expected_path, hash_key, expected_hash) in _runtime_correction_bindings(contract)
        push!(expected_keys, path_key)
        push!(expected_keys, hash_key)
        path, digest = _authority_entry(authority, path_key, hash_key)
        path == expected_path && digest == expected_hash ||
            _fail(:BLOCKED, :authority_binding_mismatch,
                  "runtime correction binding differs for $path_key")
    end
    actual_keys = Set(String(key) for key in keys(authority)
                      if startswith(String(key), "t12_authority_"))
    actual_keys == expected_keys ||
        _fail(:BLOCKED, :authority_binding_mismatch,
              "runtime correction authority key set differs")
    return nothing
end

function _runtime_snapshot(snapshots::Vector{_Snapshot}, root::String,
                           relative::String, context::String)::_Snapshot
    absolute = normpath(joinpath(root, relative))
    matches = [snapshot for snapshot in snapshots if snapshot.path == absolute]
    length(matches) == 1 ||
        _fail(:BLOCKED, :authority_binding_mismatch,
              "$context snapshot is absent or duplicated")
    return only(matches)
end

function _runtime_v3_mode(path::String, expected::UInt, context::String)
    actual = UInt(stat(path).mode) & UInt(0o777)
    actual == expected ||
        _fail(:BLOCKED, :authority_path_mismatch,
              "$context mode differs: expected $(string(expected, base=8))")
    return nothing
end

function _runtime_v3_immutable_path(path::String)::Bool
    root_suffix = "/" * _RUNTIME_V3_T12_ROOT
    return occursin(root_suffix * "/", path) ||
           endswith(path, "/" * _RUNTIME_V3_T12_PUBLICATION_PATH)
end

function _runtime_v3_correction_directory(path::String)::Bool
    root_suffix = "/" * _RUNTIME_V3_T12_ROOT
    return endswith(path, root_suffix) || occursin(root_suffix * "/", path)
end

function _runtime_v3_manifest_path(path::String)::String
    isempty(path) && _fail(:BLOCKED, :authority_path_invalid,
                           "runtime-v3 manifest path is empty")
    occursin('\0', path) &&
        _fail(:BLOCKED, :authority_path_invalid,
              "runtime-v3 manifest path contains NUL")
    isabspath(path) &&
        _fail(:BLOCKED, :authority_path_invalid,
              "runtime-v3 manifest path is absolute")
    occursin('\\', path) &&
        _fail(:BLOCKED, :authority_path_invalid,
              "runtime-v3 manifest path uses a backslash")
    parts = split(path, '/'; keepempty=true)
    any(part -> isempty(part) || part == "." || part == "..", parts) &&
        _fail(:BLOCKED, :authority_path_invalid,
              "runtime-v3 manifest path is not normalized")
    normpath(path) == path ||
        _fail(:BLOCKED, :authority_path_invalid,
              "runtime-v3 manifest path is not canonical")
    return path
end

function _runtime_v3_walk(path::String, directories::Vector{String},
                          files::Vector{String})
    islink(path) && _fail(:BLOCKED, :authority_symlink,
                          "runtime-v3 correction path is a symlink")
    isdir(path) || _fail(:BLOCKED, :authority_bundle_mismatch,
                         "runtime-v3 correction directory is absent")
    _runtime_v3_mode(path, UInt(0o555), "runtime-v3 correction directory")
    push!(directories, path)
    for name in sort(readdir(path))
        child = joinpath(path, name)
        islink(child) && _fail(:BLOCKED, :authority_symlink,
                               "runtime-v3 correction member is a symlink")
        if isdir(child)
            _runtime_v3_walk(child, directories, files)
        elseif isfile(child)
            push!(files, child)
        else
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "runtime-v3 correction member is not regular")
        end
    end
    return nothing
end

function _runtime_v3_bind_snapshot!(snapshots::Vector{_Snapshot}, root::String,
                                    relative::String, expected::String,
                                    context::String,
                                    seen::Dict{Tuple{UInt64,UInt64},String})::_Snapshot
    path = _repo_path(root, relative, context)
    matches = [snapshot for snapshot in snapshots if snapshot.path == path]
    length(matches) <= 1 ||
        _fail(:BLOCKED, :authority_identity_collision,
              "$context snapshot is duplicated")
    if isempty(matches)
        snapshot = _snapshot(root, relative, expected, context, seen)
        push!(snapshots, snapshot)
    else
        snapshot = only(matches)
        snapshot.sha256 == expected ||
            _fail(:BLOCKED, :authority_hash_mismatch, "$context bytes differ")
        _identity(path) == (snapshot.device, snapshot.inode, snapshot.nlink) ||
            _fail(:BLOCKED, :authority_snapshot_changed, "$context identity differs")
        read(path) == snapshot.bytes ||
            _fail(:BLOCKED, :authority_snapshot_changed, "$context bytes changed")
    end
    _runtime_v3_immutable_path(path) &&
        _runtime_v3_mode(path, UInt(0o444), context)
    return snapshot
end

function _validate_runtime_v3_bundle(root::String,
                                     snapshots::Vector{_Snapshot},
                                     seen::Dict{Tuple{UInt64,UInt64},String})
    bundle_root = _repo_path(root, _RUNTIME_V3_T12_ROOT,
                             "runtime-v3 correction root")
    directories = String[]
    files = String[]
    _runtime_v3_walk(bundle_root, directories, files)
    length(directories) == _RUNTIME_V3_T12_ROOT_DIRECTORY_COUNT ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime-v3 correction directory count differs")
    length(files) == _RUNTIME_V3_T12_ROOT_FILE_COUNT ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime-v3 correction file count differs")
    sum(Int(stat(path).size) for path in files) == _RUNTIME_V3_T12_ROOT_BYTES ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime-v3 correction byte count differs")

    manifest = _runtime_v3_bind_snapshot!(
        snapshots, root, _RUNTIME_V3_T12_ROOT_MANIFEST_PATH,
        _RUNTIME_V3_T12_ROOT_MANIFEST_SHA256, "runtime-v3 root manifest", seen)
    lines = split(chomp(String(copy(manifest.bytes))), '\n'; keepempty=false)
    length(lines) == _RUNTIME_V3_T12_ROOT_MANIFEST_ENTRIES ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime-v3 root manifest entry count differs")
    expected = Dict{String,String}()
    for line in lines
        match_result = match(r"^([0-9a-f]{64})  (.+)$", line)
        match_result === nothing &&
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "runtime-v3 root manifest row is malformed")
        digest, relative = match_result.captures
        relative = _runtime_v3_manifest_path(String(relative))
        haskey(expected, relative) &&
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "runtime-v3 root manifest contains a duplicate path")
        expected[relative] = String(digest)
    end
    actual = Set(relpath(path, bundle_root) for path in files
                 if path != joinpath(bundle_root, "root-manifest.sha256"))
    Set(keys(expected)) == actual ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime-v3 root manifest fileset differs")
    for relative in sort(collect(keys(expected)))
        _runtime_v3_bind_snapshot!(snapshots, root,
                                   joinpath(_RUNTIME_V3_T12_ROOT, relative),
                                   expected[relative],
                                   "runtime-v3 member $relative", seen)
    end
    _runtime_v3_bind_snapshot!(
        snapshots, root, _RUNTIME_V3_T12_PUBLICATION_PATH,
        _RUNTIME_V3_T12_PUBLICATION_SHA256,
        "runtime-v3 external publication receipt", seen)
    for path in files
        _runtime_v3_mode(path, UInt(0o444), "runtime-v3 correction member")
    end
    inventories = Tuple{String,Vector{String}}[]
    for directory in sort(directories)
        push!(inventories, (relpath(directory, root), sort(readdir(directory))))
    end
    return inventories
end

function _validate_runtime_correction_semantics(root::String,
                                                snapshots::Vector{_Snapshot},
                                                kind::Symbol)
    contract = _runtime_contract(kind)
    root_manifest = _runtime_snapshot(snapshots, root, contract.root_manifest_path,
                                      "runtime root manifest")
    claim = _runtime_snapshot(snapshots, root, contract.claim_path,
                              "runtime correction claim")
    review = _runtime_snapshot(snapshots, root, contract.review_path,
                               "runtime correction review")
    publication = _runtime_snapshot(snapshots, root, contract.publication_path,
                                    "runtime publication receipt")

    manifest_lines = split(chomp(String(copy(root_manifest.bytes))), '\n'; keepempty=false)
    length(manifest_lines) == contract.manifest_entries ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime correction root manifest entry count differs")
    required_manifest_lines = kind == :runtime_v3 ? (
        "$(contract.claim_sha256)  DoneClaim.json",
        "$(contract.review_sha256)  review/AdversarialVerify.json",
        "$(contract.edge_model_sha256)  $(contract.edge_model_product_path)",
        "$(contract.chain_sha256)  $(contract.chain_product_path)",
        "$(contract.corrected_test_sha256)  $(contract.test_product_path)",
    ) : (
        "$(contract.claim_sha256)  DoneClaim.json",
        "$(contract.review_sha256)  review/AdversarialVerify.json",
        "$(contract.edge_model_sha256)  $(contract.edge_model_product_path)",
    )
    all(line in manifest_lines for line in required_manifest_lines) ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "runtime correction root manifest members differ")

    claim_text = String(copy(claim.bytes))
    claim_specific = if kind == :runtime_v3
        _json_string(claim_text, "schema", contract.claim_schema) &&
        _json_string(claim_text, "edge_model_sha256", contract.edge_model_sha256) &&
        _json_string(claim_text, "chain_inference_sha256", contract.chain_sha256) &&
        _json_string(claim_text, "test_sha256", contract.corrected_test_sha256) &&
        _json_true(claim_text, "user_authorized") &&
        _json_true(claim_text, "canonical_and_provenance_unchanged") &&
        _json_false(claim_text, "factors_changed") &&
        _json_false(claim_text, "labels_changed") &&
        _json_false(claim_text, "log_evidence_changed") &&
        _json_false(claim_text, "model_changed") &&
        _json_false(claim_text, "recurrence_changed") &&
        _json_false(claim_text, "threshold_changed") &&
        _json_false(claim_text, "topology_changed") &&
        _json_false(claim_text, "viterbi_changed") &&
        _json_false(claim_text, "benchmark_labels_used") &&
        _json_false(claim_text, "production_data_used") &&
        _json_false(claim_text, "t13_runtime_authorized") &&
        _json_false(claim_text, "t14_authorized") &&
        _json_false(claim_text, "downstream_authorized")
    elseif kind == :runtime
        _json_string(claim_text, "canonical_todo10_bundle_sha256",
                     _RUNTIME_T12_CANONICAL_TODO10_BUNDLE_SHA256) &&
        _json_false(claim_text, "inference_logic_changed") &&
        _json_false(claim_text, "t13_authorized") &&
        _json_false(claim_text, "t14_authorized") &&
        _json_false(claim_text, "downstream_authorized")
    else
        _json_string(claim_text, "removed_invalid_equality",
                     "Todo11 input_sha256 == Todo10 feature_sha256") &&
        _json_false(claim_text, "t13_runtime_authorized") &&
        _json_false(claim_text, "t14_authorized") &&
        _json_false(claim_text, "downstream_todos_authorized")
    end
    _require_authority(
        (kind == :runtime_v3 || _json_string(claim_text, "schema", contract.claim_schema)) &&
        _json_string(claim_text, "status", "PASS") &&
        _json_true(claim_text, "terminal") &&
        (kind == :runtime_v3 ||
         (_json_string(claim_text, "class", "authority_identity_only") &&
          _json_string(claim_text, "corrected_edge_model_sha256",
                       contract.edge_model_sha256) &&
          _json_string(claim_text, "corrected_test_sha256",
                       contract.corrected_test_sha256) &&
          _json_false(claim_text, "numerics_changed") &&
          _json_false(claim_text, "scientific_behavior_changed") &&
          _json_false(claim_text, "thresholds_changed"))) &&
        claim_specific,
        "runtime correction claim semantics differ")

    review_text = String(copy(review.bytes))
    _require_authority(
        _json_string(review_text, "schema", contract.review_schema) &&
        _json_string(review_text, "status", "PASS") &&
        _json_string(review_text, "verdict", "PASS") &&
        _json_true(review_text, "terminal") &&
        _json_number_text(review_text, "confidence", "0.999") &&
        (kind == :runtime_v3 ?
         (_json_true(review_text, "t12_v3_publication_authorized") &&
          _json_false(review_text, "t13_runtime_authorized") &&
          _json_false(review_text, "t14_authorized") &&
          _json_false(review_text, "downstream_authorized")) :
         (_json_false(review_text, kind == :runtime ? "t13_authorized" :
                                                     "t13_runtime_authorized") &&
          _json_false(review_text, "t14_authorized") &&
          _json_false(review_text, kind == :runtime ? "downstream_authorized" :
                                                      "downstream_todos_authorized"))),
        "runtime correction review semantics differ")

    publication_text = String(copy(publication.bytes))
    _require_authority(
        _json_string(publication_text, "schema", contract.publication_schema) &&
        _json_string(publication_text, "status", "PASS") &&
        _json_string(publication_text, "publication_state", "committed_verified") &&
        _json_string(publication_text, "root_manifest_sha256",
                     contract.root_manifest_sha256) &&
        _json_string(publication_text, "claim_sha256", contract.claim_sha256) &&
        _json_string(publication_text, "review_sha256", contract.review_sha256) &&
        (kind == :runtime_v3 ?
         (_json_number(publication_text, "root_manifest_entries", contract.manifest_entries) &&
          _json_string(publication_text, "final_path", _RUNTIME_V3_T12_ROOT) &&
          _json_true(publication_text, "identity_matches_checkpoint") &&
          _json_true(publication_text, "no_post_rename_root_writes") &&
          _json_false(publication_text, "t13_runtime_authorized") &&
          _json_false(publication_text, "t14_authorized") &&
          _json_false(publication_text, "downstream_authorized")) :
         (_json_false(publication_text, "t13_authorized") &&
          _json_false(publication_text, "t14_authorized") &&
          _json_false(publication_text, "downstream_authorized"))),
        "runtime publication receipt semantics differ")
    if kind == :runtime_v3
        for (path, context) in (
            (contract.root_manifest_path, "runtime v3 root manifest"),
            (contract.claim_path, "runtime v3 claim"),
            (contract.review_path, "runtime v3 review"),
            (contract.publication_path, "runtime v3 publication receipt"),
        )
            artifact = _runtime_snapshot(snapshots, root, path, context)
            _readonly_mode(artifact.path, context)
        end
        for (path, context) in (
            (_RUNTIME_V3_T12_EDGE_PRODUCT_PATH, "runtime v3 edge product"),
            (_RUNTIME_V3_T12_CHAIN_PRODUCT_PATH, "runtime v3 chain product"),
            (_RUNTIME_V3_T12_TEST_PRODUCT_PATH, "runtime v3 test product"),
        )
            product = _runtime_snapshot(snapshots, root, path, context)
            _readonly_mode(product.path, context)
        end
    end
    return nothing
end

function _is_closure_path(path::String)::Bool
    suffix = "/" * _CLOSURE_V1_ROOT
    # Snapshot paths are absolute, so the closure marker is not at the start
    # of paths in an isolated fixture.  Keep the boundary check exact: a
    # similarly named sibling must not inherit the closure mode policy.
    return endswith(path, suffix) || occursin(suffix * "/", path)
end

function _is_immutable_closure_authority_path(path::String)::Bool
    return _is_closure_path(path) ||
           endswith(path, "/" * _CLOSURE_V1_PUBLICATION_RECEIPT_PATH) ||
           endswith(path, "/" * _CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_PATH)
end

function _closure_directory(root::String, relative::String, context::String)
    supplied = isempty(relative) ? _CLOSURE_V1_ROOT : _CLOSURE_V1_ROOT * "/" * relative
    path = _repo_path(root, supplied, context)
    isdir(path) && !islink(path) ||
        _fail(:BLOCKED, :authority_path_invalid, "$context is not a real directory")
    _readonly_mode(path, context)
    return path
end

function _validate_closure_manifest(root::String, manifest::_Snapshot,
                                   seen::Dict{Tuple{UInt64,UInt64},String})
    lines = split(chomp(String(copy(manifest.bytes))), '\n'; keepempty=false)
    length(lines) == 26 ||
        _fail(:BLOCKED, :authority_bundle_mismatch, "closure root manifest entry count differs")
    seen_paths = Set{String}()
    prefix = _CLOSURE_V1_ROOT * "/"
    expected_files = Dict{String,Vector{String}}(
        "" => ["root-manifest.sha256"], "guards" => String[], "review" => String[])
    members = Dict{String,_Snapshot}(_CLOSURE_V1_ROOT_MANIFEST_PATH => manifest)
    _readonly_mode(manifest.path, "closure-v1 root manifest")
    for line in lines
        match_result = match(r"^([0-9a-f]{64})  (.+)$", line)
        match_result === nothing &&
            _fail(:BLOCKED, :authority_bundle_mismatch, "closure root manifest row is malformed")
        expected, supplied = match_result.captures
        startswith(supplied, prefix) ||
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "closure root manifest escapes its root")
        supplied in seen_paths &&
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "closure root manifest contains a duplicate path")
        push!(seen_paths, supplied)
        relative = supplied[length(prefix) + 1:end]
        parts = split(relative, '/'; keepempty=false)
        length(parts) in (1, 2) ||
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "closure root manifest member has an invalid depth")
        directory = length(parts) == 1 ? "" : parts[1]
        directory in ("", "guards", "review") ||
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "closure root manifest member uses an unexpected directory")
        push!(expected_files[directory], parts[end])
        members[String(supplied)] = _snapshot(root, String(supplied), String(expected),
                                      "closure root manifest member", seen)
        _readonly_mode(members[String(supplied)].path, "closure root manifest member")
    end
    length(members) == 27 ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "closure bundle does not contain exactly 27 files")
    root_path = _closure_directory(root, "", "closure-v1 root")
    expected_root = sort(vcat(expected_files[""], ["guards", "review"]))
    sort(readdir(root_path)) == expected_root ||
        _fail(:BLOCKED, :authority_bundle_mismatch,
              "closure-v1 root inventory differs")
    inventories = Tuple{String,Vector{String}}[]
    push!(inventories, (relpath(root_path, root), copy(expected_root)))
    for directory in ("guards", "review")
        path = _closure_directory(root, directory, "closure-v1/$directory")
        expected = sort(expected_files[directory])
        sort(readdir(path)) == expected ||
            _fail(:BLOCKED, :authority_bundle_mismatch,
                  "closure-v1/$directory inventory differs")
        push!(inventories, (relpath(path, root), copy(expected)))
    end
    return (members=members,
            snapshots=sort(collect(values(members)); by=snapshot -> snapshot.path),
            inventories=inventories)
end

function _validate_closure_claim(claim::_Snapshot)
    text = String(copy(claim.bytes))
    _require_authority(_json_string(text, "schema", "stmfit-t13-structured-evaluator-done-claim-v1"),
                       "closure DoneClaim schema differs")
    _require_authority(_json_string(text, "status", "PASS") &&
                       _json_true(text, "terminal"),
                       "closure DoneClaim is not a terminal implementation pass")
    _require_authority(_json_true(text, "closure_pending_plan_and_boulder") &&
                       _json_true(text, "post_publication_receipt_required"),
                       "closure DoneClaim does not record the prerequisite history")
    _require_authority(_json_string(text, "evaluator_config_sha256", _CONFIG_SHA256) &&
                       _json_string(text, "evaluator_sha256", "7b4717077c24f122298e61cd944d55ab439ef278b32a76ec8a6fa8049735dd2e") &&
                       _json_string(text, "test_sha256", "86d758174968851e50d64fced8369749ac128f23384b879b856e72b27a0dbe78") &&
                       _json_string(text, "preclosure_plan_sha256", _PLAN_SHA256) &&
                       _json_string(text, "preclosure_boulder_sha256", _BOULDER_SHA256),
                       "closure DoneClaim product provenance differs")
    return nothing
end

function _validate_closure_review(review::_Snapshot)
    text = String(copy(review.bytes))
    _require_authority(_json_string(text, "schema", "stmfit-t13-structured-evaluator-adversarial-review-v1") &&
                       _json_string(text, "status", "PASS") &&
                       _json_string(text, "verdict", "PASS") &&
                       _json_true(text, "terminal") &&
                       _json_true(text, "read_only") &&
                       _json_string(text, "materialization_owner", "parent"),
                       "closure independent review is not an accepted pass")
    _require_authority(_json_true(text, "todo13_implementation_confirmed") &&
                       _json_true(text, "todo13_may_be_checked") &&
                       _json_true(text, "boulder_closure_authorized"),
                       "closure independent review does not authorize the closed implementation")
    _require_authority(_json_string(text, "sha256", _CLOSURE_V1_CLAIM_SHA256),
                       "closure independent review bindings differ")
    return nothing
end

function _validate_publication_receipt(receipt::_Snapshot)
    text = String(copy(receipt.bytes))
    _require_authority(_json_string(text, "schema", "stmfit-t13-closure-v1-publication-receipt-v1") &&
                       _json_string(text, "status", "PASS") &&
                       _json_true(text, "review_valid") &&
                       _json_true(text, "todo13_implementation_confirmed") &&
                       _json_true(text, "todo13_may_be_checked") &&
                       _json_false(text, "final_root_files_modified_after_rename"),
                       "closure publication receipt is not a verified pass")
    _require_authority(_json_string(text, "final_path", _CLOSURE_V1_ROOT) &&
                       _json_string(text, "publication_state",
                                    "committed_verified_after_interrupted_parent_check"),
                       "closure publication target differs")
    _require_authority(_json_string(text, "claim_sha256", _CLOSURE_V1_CLAIM_SHA256) &&
                       _json_string(text, "review_sha256", _CLOSURE_V1_REVIEW_SHA256) &&
                       _json_string(text, "review_sidecar_sha256", _CLOSURE_V1_REVIEW_SIDECAR_SHA256) &&
                       _json_string(text, "config_sha256", _CONFIG_SHA256) &&
                       _json_string(text, "evaluator_sha256", "7b4717077c24f122298e61cd944d55ab439ef278b32a76ec8a6fa8049735dd2e") &&
                       _json_string(text, "test_sha256", "86d758174968851e50d64fced8369749ac128f23384b879b856e72b27a0dbe78") &&
                       _json_string(text, "plan_preclosure_sha256", _PLAN_SHA256) &&
                       _json_string(text, "boulder_preclosure_sha256", _BOULDER_SHA256),
                       "closure publication product bindings differ")
    _require_authority(_json_string(text, "root_manifest_sha256", _CLOSURE_V1_ROOT_MANIFEST_SHA256) &&
                       _json_true(text, "exact_fileset") &&
                       _json_number(text, "file_count", 27),
                       "closure publication root binding differs")
    return nothing
end

function _check_closure(root::String, seen)
    root_manifest = _snapshot(root, _CLOSURE_V1_ROOT_MANIFEST_PATH,
                              _CLOSURE_V1_ROOT_MANIFEST_SHA256,
                              "closure-v1 root manifest", seen)
    bundle = _validate_closure_manifest(root, root_manifest, seen)
    claim = bundle.members[_CLOSURE_V1_CLAIM_PATH]
    claim_sidecar = bundle.members[_CLOSURE_V1_CLAIM_SIDECAR_PATH]
    review = bundle.members[_CLOSURE_V1_REVIEW_PATH]
    review_sidecar = bundle.members[_CLOSURE_V1_REVIEW_SIDECAR_PATH]
    publication = _snapshot(root, _CLOSURE_V1_PUBLICATION_RECEIPT_PATH,
                            _CLOSURE_V1_PUBLICATION_RECEIPT_SHA256,
                            "closure-v1 publication receipt", seen)
    publication_sidecar = _snapshot(root, _CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_PATH,
                                    _CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_SHA256,
                                    "closure-v1 publication receipt sidecar", seen)
    _readonly_mode(publication.path, "closure-v1 publication receipt")
    _readonly_mode(publication_sidecar.path, "closure-v1 publication receipt sidecar")
    _validate_sidecar(claim_sidecar, _CLOSURE_V1_CLAIM_SHA256,
                      _CLOSURE_V1_CLAIM_PATH, "closure-v1 DoneClaim")
    _validate_sidecar(review_sidecar, _CLOSURE_V1_REVIEW_SHA256,
                      _CLOSURE_V1_REVIEW_PATH, "closure-v1 review")
    _validate_sidecar(publication_sidecar, _CLOSURE_V1_PUBLICATION_RECEIPT_SHA256,
                      _CLOSURE_V1_PUBLICATION_RECEIPT_PATH,
                      "closure-v1 publication receipt")
    _validate_closure_claim(claim)
    _validate_closure_review(review)
    _validate_publication_receipt(publication)
    return (claim=claim, claim_sidecar=claim_sidecar, review=review,
            review_sidecar=review_sidecar, root_manifest=root_manifest,
            publication=publication, publication_sidecar=publication_sidecar,
            snapshots=bundle.snapshots, inventories=bundle.inventories)
end

function _check_gate(root::String, seen)
    gate_path = ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/t13/" *
                "evaluator-policy-v1/correction6-GateClosure.json"
    review_path = ".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/t13/" *
                  "evaluator-policy-v1/correction6-review/AdversarialVerify.json"
    gate = _snapshot(root, gate_path, _GATECLOSURE_SHA256, "GateClosure", seen)
    review = _snapshot(root, review_path, _REVIEW_SHA256, "Todo13 review", seen)
    gate_text = String(copy(gate.bytes))
    _require_authority(_json_string(gate_text, "schema", "stmfit-t13-correction6-gateclosure-v1") &&
                       _json_true(gate_text, "authorizes_todo13") &&
                       _json_true(gate_text, "authorization_effective"),
                       "GateClosure does not authorize Todo13 effectively")
    review_text = String(copy(review.bytes))
    _require_authority(_json_string(review_text, "schema",
                                    "stmfit-t13-correction6-independent-review-v1") &&
                       _json_string(review_text, "verdict", "PASS") &&
                       _json_true(review_text, "terminal") &&
                       _json_true(review_text, "todo13_products_absent") &&
                       _json_true(review_text, "todo13_unchecked"),
                       "Todo13 prerequisite review is not an accepted history")
    closure = _check_closure(root, seen)
    return (gate=gate, review=review, closure=closure)
end

function _authority_snapshots(root::String, config_path::String;
                              activation::Union{Nothing,_ActivationAuthorization}=nothing)
    seen = Dict{Tuple{UInt64,UInt64},String}()
    config_kind = _authority_config_kind(root, config_path)
    config_hash = config_kind == :historical ? _CONFIG_SHA256 :
                  config_kind == :runtime ? _RUNTIME_V1_CONFIG_SHA256 :
                  config_kind == :runtime_v2 ? _RUNTIME_V2_CONFIG_SHA256 :
                  config_kind == :runtime_v3 ? _RUNTIME_V3_CONFIG_SHA256 :
                  _RUNTIME_V4_CONFIG_SHA256
    config = _snapshot(root, config_path, config_hash, "evaluator config", seen)
    return _collect_authority_snapshots(root, config_path, config_kind, config, seen;
                                        activation=activation)
end

function _collect_authority_snapshots(root::String, config_path::String,
                                      config_kind::Symbol, config::_Snapshot,
                                      seen::Dict{Tuple{UInt64,UInt64},String};
                                      activation::Union{Nothing,_ActivationAuthorization}=nothing)
    config_kind == :runtime_v4 && _runtime_v4_config_bytes_guard(config.bytes)
    document = try
        TOML.parse(String(copy(config.bytes)))
    catch error
        _fail(:BLOCKED, :policy_schema_mismatch, sprint(showerror, error))
    end
    if config_kind == :runtime_v4
        _validate_runtime_v4_document(document)
        _validate_runtime_v4_identity(document["runtime"])
        # The default denial stays in force.  Only an authorization that has
        # already passed the activation receipt/spec/scope checks may continue.
        activation === nothing && _validate_runtime_v4_execution_gate(document)
    else
        get(document, "evaluator", nothing) isa AbstractDict ||
            _fail(:BLOCKED, :policy_schema_mismatch, "evaluator table is absent")
        prerequisite = get(document, "prerequisite", nothing)
        prerequisite isa AbstractDict ||
            _fail(:BLOCKED, :policy_schema_mismatch, "prerequisite table is absent")
        # These fields describe the pre-implementation gate.  They are retained
        # as typed historical configuration, but are not runtime product-state
        # invariants after the accepted implementation closure.
        get(prerequisite, "todo13_products_absent", nothing) isa Bool ||
            _fail(:BLOCKED, :policy_schema_mismatch,
                  "historical Todo13 product state is not boolean")
    end
    authority = document["authority"]
    authority isa AbstractDict ||
        _fail(:BLOCKED, :authority_binding_mismatch, "authority table is absent")
    config_kind != :historical && _validate_runtime_authority_config(authority, config_kind)
    snapshots = _Snapshot[config]
    authority_keys = sort(String.(collect(keys(authority))))
    claim_key = "t12_claim_path"
    claim_key in authority_keys &&
        (authority_keys = [claim_key; filter(key -> key != claim_key, authority_keys)])
    for key in authority_keys
        endswith(key, "_path") || continue
        hash_key = key[1:end-5] * "_sha256"
        haskey(authority, hash_key) ||
            _fail(:BLOCKED, :authority_binding_mismatch, "missing hash for $key")
        path, digest = _authority_entry(authority, key, hash_key)
        if key == "plan_path"
            digest == _PLAN_SHA256 ||
                _fail(:BLOCKED, :authority_binding_mismatch, "plan binding differs")
            # The Plan binding is pre-implementation history.  Do not read or
            # snapshot the live Plan here: a checked marker or an unrelated
            # later administrative edit must not revoke runtime authority.
        else
            # The published T12 identity-only correction includes the current
            # test snapshot.  The additive runtime config remains a byte-copy
            # of the historical policy for that field, so its stale historical
            # binding is checked while the corrected product is bound from the
            # immutable correction authority below.
            if config_kind == :historical &&
               key in ("t12_edge_model_path", "t12_chain_path", "t12_test_path")
                # These are stale pre-rebind historical declarations.  The
                # historical config remains non-authoritative because its
                # claim is checked first and fails closed; test-only synthetic
                # rebinds may exercise the otherwise unchanged historical body.
                continue
            elseif config_kind in (:runtime, :runtime_v2) &&
               key in ("t12_edge_model_path", "t12_chain_path", "t12_test_path")
                expected_live_digest = key == "t12_edge_model_path" ?
                    _runtime_contract(config_kind).edge_model_sha256 :
                    key == "t12_chain_path" ?
                    "cdc0c788d49298f50721b091535a238cd552454d9885404456c593e4625ea39f" :
                    "376847c800fec317acdc33444798869c65c87b4bb6c159c594594082cf6f6b71"
                digest == expected_live_digest ||
                    _fail(:BLOCKED, :authority_binding_mismatch,
                          "historical T12 live binding differs")
                continue
            end
            snapshot = _snapshot(root, path, digest, "authority.$key", seen)
            if config_kind == :runtime_v4
                mode_key = key[1:end-5] * "_mode"
                haskey(authority, mode_key) &&
                    _runtime_v3_mode(snapshot.path,
                                     parse(UInt, authority[mode_key]; base=8),
                                     "authority.$key")
            end
            push!(snapshots, snapshot)
        end
    end
    if config_kind == :runtime_v3
        product_bindings = (
            (_RUNTIME_V3_T12_EDGE_PRODUCT_PATH, _RUNTIME_V3_T12_EDGE_MODEL_SHA256),
            (_RUNTIME_V3_T12_CHAIN_PRODUCT_PATH, _RUNTIME_V3_T12_CHAIN_SHA256),
            (_RUNTIME_V3_T12_TEST_PRODUCT_PATH, _RUNTIME_V3_T12_CORRECTED_TEST_SHA256),
        )
        for (path, digest) in product_bindings
            push!(snapshots, _snapshot(root, path, digest, "runtime-v3.$path", seen))
        end
    end
    # The evaluator configuration binds the downstream bundles; these producer
    # members are the explicit T3/T7/T10 authority boundary checked by T11.
    extras = [
        ("test/lib/structured_assignment/universe.jl",
         "a70a6e1382627fd584e5825dbbc6464d24ba5b7d96dabe7ecc36493c79084b18"),
        ("test/lib/structured_assignment/robust_emissions.jl",
         "0ca4a3863a4b9ed13aee5cd3a0229df4f2764a8524b764eeb0fd27238f14f241"),
        ("test/test_structured_robust_emissions.jl",
         "c3deb8e82c2ba2a8720474d91f2294898b92372cf0bcbf53e1470e92a4260824"),
        ("test/lib/structured_assignment/edge_features.jl",
         "c29fe06c7fc0ba33c4e168169fec5ced5a7907bbfcec626c183e3ced044da0f9"),
        ("test/build_" * "la" * "bel_free_edge_features.jl",
         "c2133b0bc45e5ae68ab49f8ae39470c7bd2ee4fe85a52da548fbc3f84096022b"),
        ("test/test_structured_edge_features.jl",
         "21f4896cc1db508ef2dd9d5acc21196a51c7ef066e475468725db73ece98bb22"),
        (".omo/evidence/structured-" * "la" * "bel-free-unit-assignment/provenance-rebind/phase4-t10-edge-features/review/AdversarialVerify.json",
         "02d421bb51dcdb4d658da791f14cc7a6cf5f6ced26f3f25216f0b6e6cb88cb77"),
    ]
    for (path, digest) in extras
        push!(snapshots, _snapshot(root, path, digest, "producer.$path", seen))
    end
    source_path = joinpath(root, "test", "evaluate_structured_unit_assignment.jl")
    source = _current_snapshot(source_path, "evaluator source", seen)
    source === nothing || push!(snapshots, source)
    correction_inventories = Tuple{String,Vector{String}}[]
    if config_kind == :runtime_v3
        correction_inventories = _validate_runtime_v3_bundle(root, snapshots, seen)
    end
    config_kind != :historical && config_kind != :runtime_v4 &&
        _validate_runtime_correction_semantics(root, snapshots, config_kind)
    _, plan_hash = _authority_entry(authority, "plan_path", "plan_sha256")
    plan_hash == _PLAN_SHA256 ||
        _fail(:BLOCKED, :authority_binding_mismatch, "plan binding differs")
    if config_kind == :runtime_v4
        if activation === nothing
            # The native live bound-root Boulder file is the only accepted path
            # for the non-activated runtime-v4 route.
            boulder_sha256 = _RUNTIME_V4_BOULDER_SHA256
            boulder = _snapshot(root, ".omo/boulder.json", boulder_sha256,
                                "current Gate4 Boulder", seen)
            # The live Boulder snapshot must remain part of the revalidated
            # context; previously it was verified here but never propagated to
            # the caller.
            push!(snapshots, boulder)
        else
            # Main-control route: the single control file lives under the
            # reviewed Main root M, not under the execution root R.  It is
            # carried as the dedicated tagged control object and is never added
            # to the ordinary R-relative snapshot vector; R's own
            # `.omo/boulder.json`, if present, is a non-consumed decoy.
            control = activation.bindings.control
            _verify_main_control(control)
            boulder_sha256 = control.postimage_sha256
            boulder = control.file
        end
        transition = _runtime_snapshot(snapshots, root,
                                       authority["gate4_transition_receipt_path"],
                                       "Gate4 transition receipt")
        post_review = _runtime_snapshot(snapshots, root,
                                        authority["gate4_post_review_path"],
                                        "Gate4 post-review")
        oracle = _runtime_snapshot(snapshots, root,
                                   authority["gate4_oracle_final_review_path"],
                                   "Gate4 Oracle final review")
        _validate_runtime_v4_live_hashes(
            current_boulder_sha256=boulder.sha256,
            expected_boulder_sha256=boulder_sha256,
            gate4_transition_receipt_sha256=transition.sha256,
            gate4_post_review_sha256=post_review.sha256,
            gate4_oracle_final_review_sha256=oracle.sha256)
        gate = (boulder=boulder, transition=transition, post_review=post_review,
                oracle=oracle)
        closure = nothing
        inventories = correction_inventories
        authority_schema = activation === nothing ? _RUNTIME_V4_AUTHORITY_SCHEMA :
                           _T13_ACTIVATION_AUTHORITY_SCHEMA
    else
        gate = _check_gate(root, seen)
        closure = gate.closure
        append!(snapshots, (gate.gate, gate.review))
        append!(snapshots, closure.snapshots)
        push!(snapshots, closure.publication)
        push!(snapshots, closure.publication_sidecar)
        inventories = vcat(correction_inventories, closure.inventories)
        authority_schema = config_kind == :runtime_v3 ? _RUNTIME_V3_AUTHORITY_SCHEMA :
            config_kind == :runtime_v2 ?
            "schema=structured-evaluator-authority-v4" :
            config_kind == :runtime ?
            "schema=structured-evaluator-authority-v3" :
            "schema=structured-evaluator-authority-v2"
    end
    control = activation === nothing ? nothing : activation.bindings.control
    authority_sha256 = _manifest_digest(root, snapshots, inventories,
                                        [authority_schema]; control=control)
    return (snapshots=snapshots, document=document, gate=gate,
            closure=closure, inventories=inventories,
            authority_sha256=authority_sha256, config_kind=config_kind,
            config_sha256=config.sha256)
end

function _status(value::AbstractString, context::String)::Symbol
    value in ("PASS", "SKIPPED", "FAIL", "BLOCKED") ||
        _fail(:BLOCKED, :schema_mismatch, "$context status is invalid")
    return Symbol(value)
end

function _map_predecessor_error(error, reason::Symbol, context::String)
    if error isa T8.StructuredUnaryError || error isa T11.EdgeAdmissionError
        _fail(:BLOCKED, reason, "$context: $(sprint(showerror, error))")
    end
    rethrow()
end

function _view_joint(fit_view, normalized, index::Int)::NTuple{2,Float64}
    fit_view.status == :failed &&
        _fail(:FAIL, :unary_numerical_failure, "predecessor unary view failed")
    fit_view.status == :ok || return (NaN, NaN)
    values = view(normalized.normalized_predictors, index, fit_view.feature_indices)
    all(isfinite, values) || return (NaN, NaN)
    component_log = if fit_view.family == :C1
        component -> log(0.5) +
            T8.HierarchicalUnitAssignment._log_gaussian_diag(
                values, view(fit_view.means, component, :),
                view(fit_view.scales, component, :))
    elseif fit_view.family == :C2
        component -> log(0.5) + T11.StructuredRobustEmissions.log_student_t_diag(
            values, view(fit_view.means, component, :),
            view(fit_view.scales, component, :))
    else
        return (NaN, NaN)
    end
    high = fit_view.high_component
    high in (1, 2) || return (NaN, NaN)
    other = 3 - high
    result = (Float64(component_log(other)), Float64(component_log(high)))
    all(isfinite, result) ||
        _fail(:FAIL, :unary_numerical_failure, "predecessor density is nonfinite")
    return result
end

function _node_evidence(normalized, fit1, fit2, date::String)
    views1 = Dict(view.name => view for view in fit1.unary.views)
    views2 = Dict(view.name => view for view in fit2.unary.views)
    nodes = NodeEvidence[]
    for (index, node) in enumerate(normalized.data.nodes)
        node.date == date || continue
        pairs = ViewJointEvidence[]
        for name in _VIEWS
            j1 = haskey(views1, name) ? _view_joint(views1[name], normalized, index) : (NaN, NaN)
            j2 = haskey(views2, name) ? _view_joint(views2[name], normalized, index) : (NaN, NaN)
            # A C1-inactive view is structural omission for both models.  A
            # C1-active/C2-missing view remains visible to C1 and marks C2
            # unavailable for the enclosing fold; it is never a mixed node score.
            if all(isnan, j1)
                j1 = (NaN, NaN)
                j2 = (NaN, NaN)
            elseif all(isfinite, j1) && all(isnan, j2)
                nothing
            elseif all(isfinite, j1) && all(isfinite, j2)
                nothing
            else
                _fail(:FAIL, :unary_numerical_failure, "invalid production view mask")
            end
            push!(pairs, ViewJointEvidence(name, j1, j2))
        end
        push!(nodes, NodeEvidence(node.file, node.date, node.lobe, node.t_nm, pairs))
    end
    return nodes
end

function _residualizer_equal(left, right)
    (left === nothing) == (right === nothing) || return false
    left === nothing && return true
    return left.ridge == right.ridge && left.coefficients == right.coefficients &&
           left.training_row_count == right.training_row_count &&
           left.training_sha256 == right.training_sha256 &&
           left.evaluation_order == right.evaluation_order
end

function _sample_identity(sample)
    edge = sample.edge
    return _edge_identity(edge)
end

function _edge_identity(edge)
    return (edge.file, edge.date, edge.left_lobe, edge.right_lobe, edge.segment_id)
end

function _authoritative_transform_states(fit)
    states = Any[]
    seen = Set{Tuple{String,String,Int,Int,String}}()
    for state in fit.transform_edges
        status = String(state.status)
        reason = state.reason isa Symbol ? state.reason : Symbol(String(state.reason))
        finite = all(isfinite, state.endpoint_predictors)
        include = if status == "eligible" && reason == :none
            finite || _fail(:BLOCKED, :authority_or_factor_mismatch,
                            "eligible transform has nonfinite predictors")
            true
        elseif status == "unavailable" && reason == :residualizer_unavailable
            finite || _fail(:BLOCKED, :authority_or_factor_mismatch,
                            "unavailable transform has nonfinite predictors")
            true
        elseif status == "ineligible" && reason == :nonfinite_endpoint_predictor
            !finite || _fail(:BLOCKED, :authority_or_factor_mismatch,
                             "ineligible transform has finite predictors")
            false
        else
            _fail(:BLOCKED, :authority_or_factor_mismatch,
                  "unknown transform status/reason combination")
        end
        include || continue
        identity = _edge_identity(state.edge)
        identity in seen &&
            _fail(:BLOCKED, :authority_or_factor_mismatch, "duplicate authoritative transform edge")
        push!(seen, identity)
        push!(states, state)
    end
    sort!(states; by=state -> (
        state.edge.file, state.edge.date, state.edge.left_lobe,
        state.edge.right_lobe, state.edge.segment_id,
    ))
    return states
end

_eligible_transform_states(fit) = _authoritative_transform_states(fit)

function _authoritative_samples(fit)
    samples = Dict{Tuple{String,String,Int,Int,String},Any}()
    for sample in fit.samples
        identity = _sample_identity(sample)
        haskey(samples, identity) &&
            _fail(:BLOCKED, :authority_or_factor_mismatch, "duplicate transformed edge sample")
        all(isfinite, sample.residual) ||
            _fail(:FAIL, :common_null_numerical_failure, "transformed edge residual is nonfinite")
        samples[identity] = sample
    end
    states = _eligible_transform_states(fit)
    state_ids = Set(_edge_identity(state.edge) for state in states)
    sample_ids = Set(keys(samples))
    issubset(sample_ids, state_ids) ||
        _fail(:BLOCKED, :authority_or_factor_mismatch,
              "transformed edge samples include a non-authoritative edge")
    return samples, states
end

function _append_residualizer_records!(records::Vector{String}, prefix::String, fit)
    residualizer = fit.residualizer
    residualizer === nothing && (push!(records, "$prefix\tresidualizer=none"); return records)
    push!(records, "$prefix\tresidualizer_sha256=$(residualizer.training_sha256)")
    push!(records, "$prefix\tresidualizer_ridge=$(_f(residualizer.ridge))")
    push!(records, "$prefix\tresidualizer_rows=$(residualizer.training_row_count)")
    push!(records, "$prefix\tresidualizer_order=$(join(residualizer.evaluation_order, ','))")
    append!(records, "$prefix\tresidualizer_coefficient=$(_f(value))"
            for value in residualizer.coefficients)
    return records
end

function _append_null_fit_records!(records::Vector{String}, prefix::String, fit)
    null_fit = fit.null_fit
    null_fit === nothing && (push!(records, "$prefix\tnull=none"); return records)
    push!(records, "$prefix\tnull_status=$(null_fit.status)")
    push!(records, "$prefix\tnull_reason=$(null_fit.reason)")
    append!(records, "$prefix\tnull_mean=$(_f(value))" for value in null_fit.mean)
    append!(records, "$prefix\tnull_scale=$(_f(value))" for value in null_fit.scale)
    push!(records, "$prefix\tnull_objective=$(_f(null_fit.objective))")
    return records
end

function _common_null_cache_key(c1fit, c2fit, c1states, c2states,
                                c1samples, c2samples)::String
    records = String["schema=structured-evaluator-common-null-reference-v4"]
    for (prefix, fit, states, samples) in (("c1", c1fit, c1states, c1samples),
                                           ("c2", c2fit, c2states, c2samples))
        dates = sort(unique(String.(fit.training_dates)))
        push!(records, "$prefix\tfit_sha256=$(fit.fit_sha256)")
        push!(records, "$prefix\ttraining_dates=$(join(dates, ','))")
        push!(records, "$prefix\tunary_status=$(fit.unary.status)")
        push!(records, "$prefix\tunary_reason=$(fit.unary.reason)")
        _append_residualizer_records!(records, prefix, fit)
        _append_null_fit_records!(records, prefix, fit)
        training_states = [state for state in states if state.edge.date in dates]
        for state in training_states
            identity = _edge_identity(state.edge)
            sample = get(samples, identity, nothing)
            push!(records, "$prefix\tedge\t$(join(string.(identity), '|'))\tstatus=$(state.status)\treason=$(state.reason)\tavailable=$(sample !== nothing)")
            sample === nothing && continue
            push!(records, "$prefix\traw_sha256=$(sample.raw_sha256)")
            append!(records, "$prefix\tendpoint=$(_f(value))" for value in sample.endpoint_predictors)
            append!(records, "$prefix\tresidual=$(_f(value))" for value in sample.residual)
        end
    end
    return _length_digest(records)
end

function _common_null(c1fit, c2fit, target_date::String;
                      cache::Union{Nothing,Dict}=nothing)
    c1fit.training_dates == c2fit.training_dates ||
        _fail(:BLOCKED, :authority_or_factor_mismatch, "common-null training dates differ")
    c1sample_map, c1states = _authoritative_samples(c1fit)
    c2sample_map, c2states = _authoritative_samples(c2fit)
    c1ids = [_edge_identity(state.edge) for state in c1states]
    c2ids = [_edge_identity(state.edge) for state in c2states]
    if c2fit.unary.status == "PASS"
        _residualizer_equal(c1fit.residualizer, c2fit.residualizer) ||
            _fail(:BLOCKED, :authority_or_factor_mismatch, "common-null residualizers differ")
        c1ids == c2ids ||
            _fail(:BLOCKED, :authority_or_factor_mismatch, "eligible edge identities differ")
        Set(keys(c1sample_map)) == Set(keys(c2sample_map)) ||
            _fail(:BLOCKED, :authority_or_factor_mismatch, "eligible sample availability differs")
        for identity in sort(collect(keys(c1sample_map)))
            c1sample_map[identity].residual == c2sample_map[identity].residual ||
                _fail(:BLOCKED, :authority_or_factor_mismatch, "eligible residual vectors differ")
        end
    end
    key = _common_null_cache_key(c1fit, c2fit, c1states, c2states,
                                 c1sample_map, c2sample_map)
    training_dates = Set(String.(c1fit.training_dates))
    training_states = [state for state in c1states if state.edge.date in training_dates]
    missing_training = [state for state in training_states
                        if !haskey(c1sample_map, _edge_identity(state.edge))]
    training = [c1sample_map[_edge_identity(state.edge)] for state in training_states
                if haskey(c1sample_map, _edge_identity(state.edge))]
    if cache !== nothing && haskey(cache, key)
        null, digest = cache[key]
    else
        null = if !isempty(missing_training) || isempty(training)
            nothing
        elseif c1fit.null_fit !== nothing && c1fit.null_fit.status == :ok
            c1fit.null_fit
        else
            residuals = reduce(vcat, (reshape(sample.residual, 1, :) for sample in training))
            fitted = T11.fit_state_independent(residuals)
            fitted.status == :ok ||
                _fail(:FAIL, :common_null_numerical_failure, "common null fit failed")
            fitted
        end
        digest = if null === nothing
            ""
        else
            records = String["schema=structured-evaluator-common-null-v4", "reference=$key"]
            append!(records, "mean=$(_f(value))" for value in null.mean)
            append!(records, "scale=$(_f(value))" for value in null.scale)
            _length_digest(records)
        end
        cache !== nothing && (cache[key] = (null, digest))
    end
    if null !== nothing && c2fit.unary.status == "PASS" &&
       c2fit.null_fit !== nothing && c2fit.null_fit.status == :ok
        c2fit.null_fit.mean == null.mean && c2fit.null_fit.scale == null.scale ||
            _fail(:BLOCKED, :authority_or_factor_mismatch, "published common null differs")
    end
    has_heldout_edges = any(state.edge.date == target_date for state in c1states)
    common_available = null !== nothing || !has_heldout_edges
    return null, digest, common_available
end

function _edge_evidence(fit, null, target_date::String)
    null === nothing && return EdgeNullEvidence[]
    samples, states = _authoritative_samples(fit)
    heldout_states = [state for state in states if state.edge.date == target_date]
    all(haskey(samples, _edge_identity(state.edge)) for state in heldout_states) ||
        return EdgeNullEvidence[]
    edges = EdgeNullEvidence[]
    for state in heldout_states
        edge = state.edge
        sample = samples[_edge_identity(edge)]
        value = T11.log_student_t_full(sample.residual, null.mean, null.scale)
        isfinite(value) || _fail(:FAIL, :common_null_numerical_failure, "common null score failed")
        push!(edges, EdgeNullEvidence(edge.file, edge.date, edge.left_lobe,
                                      edge.right_lobe, edge.segment_id, value))
    end
    return sort(edges; by=edge -> (edge.file, edge.left_lobe, edge.right_lobe, edge.segment_id))
end

function _edge_counts(fit, target_date::String)
    counts = Dict{String,Int}()
    for identity in fit.unary.node_identities
        file, date = identity[1], identity[2]
        date == target_date && (counts[file] = 0)
    end
    for state in _eligible_transform_states(fit)
        state.edge.date == target_date || continue
        counts[state.edge.file] = get(counts, state.edge.file, 0) + 1
    end
    return counts
end

function _partition_by_date(gate, date::String)
    matches = [partition for partition in gate.partitions if partition.target_date == date]
    length(matches) == 1 ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence, "partition cardinality differs")
    return only(matches)
end

function _placeholder_graph(model::String, unary_hash::String, fit_hash::String,
                            null_hash::String, gate)
    return GraphEvidence(model, unary_hash, fit_hash, null_hash,
                         _status(gate.status, "T11 inner gate"),
                         gate.gate_sha256, :SKIPPED, :unary_only, "",
                         Dict{String,Float64}(), Dict{Tuple{String,Int},Float64}())
end

function _directory_path(root::String, supplied::AbstractString, context::String)
    path = _repo_path(root, supplied, context)
    isdir(path) || _fail(:BLOCKED, :authority_path_invalid, "$context is not a directory")
    return path
end

function _current_snapshot(path::String, context::String, seen::Dict{Tuple{UInt64,UInt64},String})
    _reject_link_components(path)
    isfile(path) || _fail(:BLOCKED, :incomplete_or_stale_evidence, "$context is absent")
    identity = _identity(path)
    identity[3] == 1 || _fail(:BLOCKED, :authority_hardlink, "$context is hard-linked")
    key = (identity[1], identity[2])
    if haskey(seen, key)
        seen[key] == path && return nothing
        _fail(:BLOCKED, :authority_identity_collision,
              "$context shares a file identity with $(seen[key])")
    end
    seen[key] = path
    bytes = Vector{UInt8}(read(path))
    _identity(path) == identity ||
        _fail(:BLOCKED, :authority_snapshot_changed, "$context changed during read")
    return _Snapshot(path, bytes, _hash_bytes(bytes), identity[1], identity[2], identity[3])
end

function _snapshot_directory(directory::String, context::String,
                             seen::Dict{Tuple{UInt64,UInt64},String})
    _reject_link_components(directory)
    names = sort(readdir(directory))
    snapshots = _Snapshot[]
    for name in names
        path = joinpath(directory, name)
        isfile(path) || islink(path) ||
            _fail(:BLOCKED, :schema_mismatch, "$context contains a non-file member")
        snap = _current_snapshot(path, "$context/$name", seen)
        snap === nothing || push!(snapshots, snap)
    end
    return snapshots
end

function _manifest_records(root::String, snapshots::Vector{_Snapshot},
                           inventories::Vector{Tuple{String,Vector{String}}},
                           extra::Vector{String}=String[];
                           control::Union{Nothing,_MainControl}=nothing)
    records = String["schema=structured-evaluator-provenance-manifest-v2"]
    for snapshot in sort(copy(snapshots); by=snapshot -> relpath(snapshot.path, root))
        relative = relpath(snapshot.path, root)
        # Ordinary records are strictly R-relative.  A path outside the
        # execution root must never be smuggled into the ordinary snapshot
        # vector; the one Main-control record is appended separately below.
        (relative == ".." || startswith(relative, "../")) &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "ordinary manifest path escapes the execution root")
        push!(records, "file\t$relative\t$(snapshot.sha256)\t$(length(snapshot.bytes))")
    end
    for (directory, names) in sort(copy(inventories); by=first)
        push!(records, "directory\t$directory\t$(join(sort(names), ','))")
    end
    if control !== nothing
        # Exactly one fixed tagged Main-control record.  It records the
        # canonical Main root, the fixed relative control filename, the file
        # hash/bytes and the approved modes.  No inode numbers enter the
        # deterministic content digest.  Defensively, the control file must not
        # also appear among the ordinary snapshots that precede the tagged
        # record.
        _reject_main_control_overlap(snapshots, control, "provenance manifest")
        push!(records, join((
            "main_control",
            "role=$(_T13_CONTROL_ROLE)",
            "root=$(control.root)",
            "file=$(_T13_CONTROL_FILENAME)",
            "root_mode=$(_t13_mode_text(control.root_mode))",
            "parent_mode=$(_t13_mode_text(control.parent_mode))",
            "file_mode=$(_t13_mode_text(control.file_mode))",
            "sha256=$(control.file.sha256)",
            "bytes=$(length(control.file.bytes))",
        ), '\t'))
    end
    append!(records, extra)
    return records
end

function _manifest_digest(root::String, snapshots::Vector{_Snapshot},
                          inventories::Vector{Tuple{String,Vector{String}}},
                          extra::Vector{String}=String[];
                          control::Union{Nothing,_MainControl}=nothing)::String
    return _length_digest(_manifest_records(root, snapshots, inventories, extra;
                                            control=control))
end

function _snapshot_context(root::String, options::Dict{String,String},
                           authority::NamedTuple)::Tuple{Vector{_Snapshot},Vector{Tuple{String,Vector{String}}}}
    seen = Dict{Tuple{UInt64,UInt64},String}(
        (snapshot.device, snapshot.inode) => snapshot.path
        for snapshot in authority.snapshots)
    snapshots = _Snapshot[]
    inventories = Tuple{String,Vector{String}}[]
    for (path, context) in (
        (_repo_path(root, options["--features"], "features"), "features"),
        (_repo_path(root, options["--candidate-config"], "candidate config"), "candidate config"),
        (_repo_path(root, options["--model-config"], "model config"), "model config"),
        (_repo_path(root, options["--forward-receipt"], "forward receipt"), "forward receipt"),
        (_repo_path(root, options["--backward-receipt"], "backward receipt"), "backward receipt"),
    )
        snap = _current_snapshot(path, context, seen)
        snap === nothing || push!(snapshots, snap)
    end
    for (supplied, context) in ((options["--universe-dir"], "universe"),
                                (options["--edge-dir"], "edge"),
                                (options["--admission-dir"], "admission"))
        directory = _directory_path(root, supplied, context)
        names = sort(readdir(directory))
        push!(inventories, (relpath(directory, root), names))
        append!(snapshots, _snapshot_directory(directory, context, seen))
    end
    return snapshots, inventories
end

function _merge_snapshots(base::Vector{_Snapshot}, extra::Vector{_Snapshot},
                          context::String)
    by_path = Dict{String,_Snapshot}(snapshot.path => snapshot for snapshot in base)
    for snapshot in extra
        existing = get(by_path, snapshot.path, nothing)
        if existing === nothing
            by_path[snapshot.path] = snapshot
        elseif existing.sha256 != snapshot.sha256 ||
               (existing.device, existing.inode, existing.nlink) !=
               (snapshot.device, snapshot.inode, snapshot.nlink) ||
               existing.bytes != snapshot.bytes
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "$context bytes or identity changed during binding")
        end
    end
    merged = collect(values(by_path))
    sort!(merged; by=snapshot -> snapshot.path)
    return merged
end

function _verify_activation_bindings(bindings::_ActivationBindings)
    for (path, expected) in bindings.modes
        (isfile(path) && !islink(path)) || (isdir(path) && !islink(path)) ||
            _fail(:BLOCKED, :activation_snapshot_changed,
                  "activation binding disappeared: $path")
        actual = UInt(stat(path).mode) & UInt(0o777)
        actual == expected ||
            _fail(:BLOCKED, :activation_mode_mismatch,
                  "activation binding mode changed: $path")
    end
    for (path, expected) in bindings.directory_identities
        (isfile(path) || isdir(path)) && !islink(path) ||
            _fail(:BLOCKED, :activation_snapshot_changed,
                  "activation identity path disappeared: $path")
        info = stat(path)
        (UInt64(info.device), UInt64(info.inode)) == expected ||
            _fail(:BLOCKED, :activation_snapshot_changed,
                  "activation directory identity changed: $path")
    end
    # The same captured Main-control object is revalidated on every receiving
    # boundary: its file bytes/mode/link/inode and the M / M/.omo directory
    # identities and modes.  No identity is rebased or repaired.
    _verify_main_control(bindings.control)
    return nothing
end

# The single live Main-control file is a dedicated control-plane object and
# must never double as an ordinary R-relative snapshot.  Overlap is rejected by
# the captured fixed pathname and by the captured dev/inode (a different
# pathname that aliases the same file is equally rejected).  Nothing is
# silently dropped: an overlapping input is a hard rejection.
function _reject_main_control_overlap(snapshots::Vector{_Snapshot},
                                      control::_MainControl,
                                      context::String)
    control_key = (control.file.device, control.file.inode)
    for snapshot in snapshots
        snapshot.path == control.file.path &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "$context contains the Main control file as an ordinary snapshot")
        (snapshot.device, snapshot.inode) == control_key &&
            _fail(:BLOCKED, :activation_binding_mismatch,
                  "$context ordinary snapshot shares the Main control file identity")
    end
    return nothing
end

function _verify_context(context::_ProductionContext)
    for snapshot in context.snapshots
        _verify_snapshot(snapshot, "production input snapshot")
    end
    for (directory, expected) in context.inventories
        # Directory names are repository-relative and are resolved only for
        # revalidation; membership is compared before any predecessor call.
        path = _repo_path(context.root, directory, "production directory inventory")
        isdir(path) && !islink(path) ||
            _fail(:BLOCKED, :authority_snapshot_changed, "production directory disappeared")
        _is_closure_path(path) && _readonly_mode(path, "production directory inventory")
        _runtime_v3_correction_directory(path) &&
            _runtime_v3_mode(path, UInt(0o555), "runtime-v3 correction directory")
        sort(readdir(path)) == sort(expected) ||
            _fail(:BLOCKED, :authority_snapshot_changed, "production directory membership changed")
    end
    if context.activation !== nothing
        _verify_activation_bindings(context.activation)
        _reject_main_control_overlap(context.snapshots, context.activation.control,
                                     "production context")
    end
    return nothing
end

function _find_outer(result, date::String)
    matches = [fold for fold in result.outer_folds if fold.outer_date == date]
    length(matches) == 1 ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence, "outer fold cardinality differs")
    return only(matches)
end

function _make_outer_fold(normalized, c1outer, c2outer, outer_date::String,
                          common_null_cache::Dict)
    c1partition = c1outer.outer_score
    c2partition = c2outer.outer_score
    c1fit = c1partition.fit
    c2fit = c2partition.fit
    common_null, null_hash, common_available = _common_null(
        c1fit, c2fit, outer_date; cache=common_null_cache)
    eligible_edge_count = sum(values(_edge_counts(c1fit, outer_date)); init=0)
    common_available == (common_null !== nothing || eligible_edge_count == 0) ||
        _fail(:BLOCKED, :authority_or_factor_mismatch,
              "common-null availability is inconsistent with held-out edges")
    nodes = _node_evidence(normalized, c1fit, c2fit, outer_date)
    edges = _edge_evidence(c1fit, common_null, outer_date)
    inner = InnerFoldEvidence[]
    inner_dates = sort(unique(partition.target_date for partition in c1outer.inner_gate.partitions))
    for date in inner_dates
        c1p = _partition_by_date(c1outer.inner_gate, date)
        c2p = _partition_by_date(c2outer.inner_gate, date)
        _common_null(c1p.fit, c2p.fit, date; cache=common_null_cache)
        push!(inner, InnerFoldEvidence(
            date, copy(c1p.training_dates), _status(c1p.fit.unary.status, "inner C1 unary"),
            c1p.fit.unary.reason, c1p.fit.unary.fit_sha256,
            _status(c2p.fit.unary.status, "inner C2 unary"), c2p.fit.unary.reason,
            c2p.fit.unary.fit_sha256, _node_evidence(normalized, c1p.fit, c2p.fit, date),
            _edge_counts(c1p.fit, date), c1p.fit.fit_sha256, c2p.fit.fit_sha256,
        ))
    end
    graphs = Dict{String,GraphEvidence}(
        "C1" => _placeholder_graph("C1", c1fit.unary.fit_sha256,
                                    c1fit.fit_sha256, null_hash, c1outer.inner_gate),
        "C2" => _placeholder_graph("C2", c2fit.unary.fit_sha256,
                                    c2fit.fit_sha256, null_hash, c2outer.inner_gate),
    )
    return OuterFoldEvidence(
        outer_date, copy(c1fit.training_dates), inner,
        _status(c1fit.unary.status, "outer C1 unary"), c1fit.unary.reason,
        c1fit.unary.fit_sha256, _status(c2fit.unary.status, "outer C2 unary"),
        c2fit.unary.reason, c2fit.unary.fit_sha256, isempty(null_hash) ? "" : null_hash,
        eligible_edge_count, nodes, edges, graphs,
    )
end

_make_outer_fold(normalized, c1outer, c2outer, outer_date::String) =
    _make_outer_fold(normalized, c1outer, c2outer, outer_date, Dict{String,Any}())

function _reference_for(report, model::String, outer_date::String)
    matches = [reference for reference in report.references
               if reference.model == model && reference.fit_role == "partition" &&
                  reference.scope == "outer_score" && reference.outer_date == outer_date &&
                  reference.target_date == outer_date]
    length(matches) == 1 ||
        _fail(:BLOCKED, :t12_reference_mismatch,
              "T12 reference cardinality differs for $model/$outer_date")
    return only(matches)
end

function _graph_from_reference(graph::GraphEvidence, reference, report, nodes, edges)
    reference.fit_sha256 == graph.t11_fit_sha256 ||
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 fit does not match T11 fit")
    grouped_logz = Dict{String,Float64}()
    marginals = Dict{Tuple{String,Int},Float64}()
    factor_keys = Set{Tuple{String,Int,Int,String}}()
    block_ids = Set{String}()
    block_logz = Float64[]
    for block in reference.blocks
        block.segment_id in block_ids &&
            _fail(:BLOCKED, :t12_reference_mismatch, "T12 block is duplicated")
        push!(block_ids, block.segment_id)
        block_files = Set{String}()
        block_dates = Set{String}()
        for factor in block.factors
            push!(block_files, factor.file)
            push!(block_dates, factor.date)
            key = (factor.file, factor.left_lobe, factor.right_lobe, factor.segment_id)
            key in factor_keys &&
                _fail(:BLOCKED, :t12_reference_mismatch, "T12 factor is duplicated")
            push!(factor_keys, key)
            factor.unary_fit_sha256 == graph.selected_unary_fit_sha256 ||
                _fail(:BLOCKED, :t12_reference_mismatch, "T12 unary reference differs")
            factor.fit_sha256 == graph.t11_fit_sha256 ||
                _fail(:BLOCKED, :t12_reference_mismatch, "T12 factor fit differs")
        end
        for node in block.nodes
            push!(block_files, node.file)
            push!(block_dates, node.date)
            key = (node.file, node.lobe)
            haskey(marginals, key) &&
                _fail(:BLOCKED, :t12_reference_mismatch, "T12 node is duplicated")
            marginals[key] = node.marginals[2]
        end
        length(block_files) == 1 && block_dates == Set([reference.target_date]) ||
            _fail(:BLOCKED, :t12_reference_mismatch, "T12 block mixes scan identity")
        file = only(block_files)
        grouped_logz[file] = get(grouped_logz, file, 0.0) + block.log_evidence
        push!(block_logz, block.log_evidence)
    end
    reference.log_evidence == sum(block_logz; init=0.0) ||
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 reference log evidence differs")
    expected_nodes = Set((node.file, node.lobe) for node in nodes)
    Set(keys(marginals)) == expected_nodes ||
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 node set differs")
    expected_scans = Set(node.file for node in nodes)
    Set(keys(grouped_logz)) == expected_scans ||
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 scan set differs")
    expected_factors = Set((edge.file, edge.left_lobe, edge.right_lobe, edge.segment_id)
                           for edge in edges)
    factor_keys == expected_factors ||
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 factor set differs")
    all(isfinite(value) && 0.0 <= value <= 1.0 for value in values(marginals)) ||
        _fail(:FAIL, :t12_failure, "T12 marginal is invalid")
    all(isfinite(value) for value in values(grouped_logz)) ||
        _fail(:FAIL, :t12_failure, "T12 log normalizer is invalid")
    return GraphEvidence(graph.model, graph.selected_unary_fit_sha256,
                         graph.t11_fit_sha256, graph.common_null_sha256,
                         graph.t11_gate_status, graph.t11_gate_sha256,
                         reference.status, reference.reason, report.result_sha256,
                         grouped_logz, marginals)
end

function _manifest_node_records(prefix::String, node::NodeEvidence)
    records = String[
        "node\t$prefix\t$(node.file)\t$(node.date)\t$(node.lobe)\t$(_f(node.t_nm))\t$(length(node.views))",
    ]
    append!(records, begin
        "node_view\t$prefix\t$(node.file)\t$(node.date)\t$(node.lobe)\t$index\t$(view.name)\t$(_f(view.c1_log_joint[1]))\t$(_f(view.c1_log_joint[2]))\t$(_f(view.c2_log_joint[1]))\t$(_f(view.c2_log_joint[2]))"
    end for (index, view) in enumerate(node.views))
    return records
end

function _manifest_edge_record(prefix::String, edge::EdgeNullEvidence)
    return "edge\t$prefix\t$(edge.file)\t$(edge.date)\t$(edge.left_lobe)\t$(edge.right_lobe)\t$(edge.segment_id)\t$(_f(edge.null_loglik))"
end

function _input_manifest_records(input::EvaluationInput)
    records = String["schema=structured-evaluator-input-evidence-v3"]
    append!(records, "scan\t$(pair[1])\t$(pair[2])"
            for pair in sort(copy(input.scan_dates); by=pair -> (pair[2], pair[1])))
    for fold in sort(copy(input.folds); by=fold -> fold.outer_date)
        push!(records, join((
            "outer", fold.outer_date, join(sort(copy(fold.training_dates)), ','),
            String(fold.c1_status), String(fold.c1_reason), fold.c1_fit_sha256,
            String(fold.c2_status), String(fold.c2_reason), fold.c2_fit_sha256,
            fold.common_null_sha256, string(fold.eligible_edge_count),
        ), '\t'))
        for inner in sort(copy(fold.inner_folds); by=inner -> inner.heldout_date)
            push!(records, join((
                "inner", fold.outer_date, inner.heldout_date,
                join(sort(copy(inner.training_dates)), ','),
                String(inner.c1_status), String(inner.c1_reason), inner.c1_fit_sha256,
                String(inner.c2_status), String(inner.c2_reason), inner.c2_fit_sha256,
            ), '\t'))
            append!(records, "inner_edge\t$(fold.outer_date)\t$(inner.heldout_date)\t$file\t$count"
                    for (file, count) in sort(collect(inner.eligible_edge_counts); by=first))
            inner_prefix = "$(fold.outer_date)|inner|$(inner.heldout_date)"
            for node in sort(copy(inner.nodes); by=node -> (node.file, node.date, node.lobe, node.t_nm))
                append!(records, _manifest_node_records(inner_prefix, node))
            end
        end
        outer_prefix = "$(fold.outer_date)|outer"
        for node in sort(copy(fold.nodes); by=node -> (node.file, node.date, node.lobe, node.t_nm))
            append!(records, _manifest_node_records(outer_prefix, node))
        end
        for edge in sort(copy(fold.edges); by=edge ->
                         (edge.file, edge.date, edge.left_lobe, edge.right_lobe, edge.segment_id))
            push!(records, _manifest_edge_record(outer_prefix, edge))
        end
        for model in sort(collect(keys(fold.graphs)))
            graph = fold.graphs[model]
            push!(records, join((
                "graph", fold.outer_date, model, graph.selected_unary_fit_sha256,
                graph.t11_fit_sha256, graph.common_null_sha256,
                String(graph.t11_gate_status), graph.t11_gate_sha256,
                String(graph.t12_status), String(graph.t12_reason),
                graph.t12_result_sha256,
            ), '\t'))
            append!(records, "graph_scan\t$(fold.outer_date)\t$model\t$file\t$(_f(value))"
                    for (file, value) in sort(collect(graph.scan_logz); by=first))
            append!(records, "graph_node\t$(fold.outer_date)\t$model\t$(key[1])\t$(key[2])\t$(_f(value))"
                    for (key, value) in sort(collect(graph.node_p1); by=first))
        end
    end
    return records
end

function _report_bindings(input::EvaluationInput)
    source_path = abspath(joinpath(@__DIR__, "evaluate_structured_unit_assignment.jl"))
    context = get(_PRODUCTION_CONTEXT, input, nothing)
    if context !== nothing
        matches = [snapshot for snapshot in context.snapshots
                   if snapshot.path == source_path]
        length(matches) == 1 ||
            throw(ArgumentError("evaluator source snapshot is absent"))
        snapshot = only(matches)
        return snapshot.sha256, context.manifest_sha256
    end
    source_sha = _hash_bytes(read(source_path))
    provenance = _length_digest(vcat([
        "schema=structured-evaluator-provenance-v3",
        input.authority_sha256, input.universe_sha256, "source_sha256=$source_sha",
    ], _input_manifest_records(input)))
    return source_sha, provenance
end

function _convert_production_input(root::String, options::Dict{String,String},
                                   admission, data, normalized,
                                   chain_report)::EvaluationInput
    chain_report.status in (:PASS, :SKIPPED, :FAIL, :BLOCKED) ||
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 returned an unknown status")
    c1 = admission.models["C1"]
    c2 = admission.models["C2"]
    outer_dates = sort(unique(node.date for node in data.nodes))
    common_null_cache = Dict{String,Any}()
    folds = OuterFoldEvidence[]
    for date in outer_dates
        c1outer = _find_outer(c1, date)
        c2outer = _find_outer(c2, date)
        push!(folds, _make_outer_fold(normalized, c1outer, c2outer, date,
                                      common_null_cache))
    end
    if chain_report.status == :BLOCKED
        _fail(:BLOCKED, :t12_reference_mismatch, "T12 validation was blocked")
    end
    if chain_report.status == :PASS
        references = [reference for reference in chain_report.references
                      if reference.fit_role == "partition" &&
                         reference.scope == "outer_score"]
        length(references) == 2 * length(folds) ||
            _fail(:BLOCKED, :t12_reference_mismatch, "T12 outer reference count differs")
        reference_keys = [(reference.model, reference.outer_date, reference.target_date)
                          for reference in references]
        length(reference_keys) == length(unique(reference_keys)) ||
            _fail(:BLOCKED, :t12_reference_mismatch, "T12 outer reference is duplicated")
    end
    updated = OuterFoldEvidence[]
    for fold in folds
        graphs = Dict{String,GraphEvidence}()
        for model in ("C1", "C2")
            if chain_report.status == :FAIL
                base = fold.graphs[model]
                graphs[model] = GraphEvidence(
                    base.model, base.selected_unary_fit_sha256, base.t11_fit_sha256,
                    base.common_null_sha256, base.t11_gate_status, base.t11_gate_sha256,
                    :FAIL, chain_report.reason, chain_report.result_sha256,
                    Dict{String,Float64}(), Dict{Tuple{String,Int},Float64}())
            elseif chain_report.status == :SKIPPED
                base = fold.graphs[model]
                graphs[model] = GraphEvidence(
                    base.model, base.selected_unary_fit_sha256, base.t11_fit_sha256,
                    base.common_null_sha256, base.t11_gate_status, base.t11_gate_sha256,
                    :SKIPPED, chain_report.reason, chain_report.result_sha256,
                    Dict{String,Float64}(), Dict{Tuple{String,Int},Float64}())
            else
                reference = _reference_for(chain_report, model, fold.outer_date)
                if reference.status == :PASS && fold.graphs[model].t11_gate_status == :PASS
                    graphs[model] = _graph_from_reference(
                        fold.graphs[model], reference, chain_report, fold.nodes, fold.edges)
                else
                    base = fold.graphs[model]
                    graphs[model] = GraphEvidence(
                        base.model, base.selected_unary_fit_sha256, base.t11_fit_sha256,
                        base.common_null_sha256, base.t11_gate_status, base.t11_gate_sha256,
                        reference.status, reference.reason, chain_report.result_sha256,
                        Dict{String,Float64}(), Dict{Tuple{String,Int},Float64}())
                end
            end
        end
        push!(updated, OuterFoldEvidence(
            fold.outer_date, fold.training_dates, fold.inner_folds,
            fold.c1_status, fold.c1_reason, fold.c1_fit_sha256,
            fold.c2_status, fold.c2_reason, fold.c2_fit_sha256,
            fold.common_null_sha256, fold.eligible_edge_count, fold.nodes, fold.edges, graphs,
        ))
    end
    scans = sort([(node.file, node.date) for node in data.nodes]; by=pair -> (pair[2], pair[1]))
    scans = unique(scans)
    universe = get(data.hashes, "universe_sha256", "")
    _validate_hash(universe, "universe hash")
    authority = get(data.hashes, "source_sha256", "")
    _validate_hash(authority, "production authority hash")
    return EvaluationInput(authority, universe, scans, updated)
end

function _load_production_input(options::Dict{String,String};
                                activation::Union{Nothing,_ActivationAuthorization}=nothing)::EvaluationInput
    # Keep the Gate5 candidate guard at the input boundary as well as in main:
    # the producer seam must not become a bypass around the preflight block.
    # The preflight returns the single complete, already-revalidated context
    # (ordinary inputs, or activation plus approved inputs); reuse it instead of
    # assembling a second, weaker context here.
    context = _preflight_config_context(options; activation=activation)
    root = options["--root"]
    isabspath(root) || _fail(:BLOCKED, :authority_path_invalid, "root must be absolute")
    isdir(root) || _fail(:BLOCKED, :authority_path_invalid, "root is not a directory")
    root = normpath(root)
    realpath(root) == root || _fail(:BLOCKED, :authority_path_invalid, "root is not canonical")
    config_path = _repo_path(root, options["--evaluator-config"], "evaluator config")
    config_path == joinpath(root, _RUNTIME_V4_CONFIG_PATH) ||
        _fail(:BLOCKED, :authority_path_mismatch, "evaluator config path differs")
    authority = _authority_snapshots(root, options["--evaluator-config"];
                                     activation=activation)
    candidate = _repo_path(root, options["--candidate-config"], "candidate config")
    model = _repo_path(root, options["--model-config"], "model config")
    features = _repo_path(root, options["--features"], "features")
    universe_dir = _directory_path(root, options["--universe-dir"], "universe")
    edge_dir = _directory_path(root, options["--edge-dir"], "edge bundle")
    admission_dir = _directory_path(root, options["--admission-dir"], "admission output")
    _repo_path(root, options["--forward-receipt"], "forward receipt")
    _repo_path(root, options["--backward-receipt"], "backward receipt")
    _verify_context(context)
    contract = try
        result = T8.load_contract(model)
        T8.source_hash()
        result
    catch error
        _map_predecessor_error(error, :authority_binding_mismatch, "T8 contract")
    end
    _verify_context(context)
    data = try
        T11.load_admission_data(root;
            features=options["--features"], candidate_config=options["--candidate-config"],
            model_config=options["--model-config"], universe_dir=options["--universe-dir"],
            edge_dir=options["--edge-dir"], forward_receipt=options["--forward-receipt"],
            backward_receipt=options["--backward-receipt"])
    catch error
        _map_predecessor_error(error, :authority_binding_mismatch, "T11 input")
    end
    _verify_context(context)
    admission = try
        T11.evaluate_admission(data)
    catch error
        _map_predecessor_error(error, :authority_or_factor_mismatch, "T11 admission")
    end
    _verify_context(context)
    generated = try
        T11.report_files(admission)
    catch error
        _map_predecessor_error(error, :incomplete_or_stale_evidence, "T11 report")
    end
    sort(readdir(admission_dir)) == sort(collect(keys(generated))) ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence, "admission file set differs")
    for name in sort(collect(keys(generated)))
        path = joinpath(admission_dir, name)
        isfile(path) && !islink(path) && read(path) == generated[name] ||
            _fail(:BLOCKED, :incomplete_or_stale_evidence, "admission bytes differ for $name")
    end
    admission.status == "BLOCKED" &&
        _fail(:BLOCKED, :authority_or_factor_mismatch, "T11 admission is blocked")
    admission.status in ("PASS", "SKIPPED", "FAIL") ||
        _fail(:BLOCKED, :authority_or_factor_mismatch, "T11 returned an unknown status")
    _verify_context(context)
    normalized = try
        T11.normalize_admission_data(data)
    catch error
        _map_predecessor_error(error, :authority_or_factor_mismatch, "T11 normalization")
    end
    _verify_context(context)
    chain_report = admission.status == "FAIL" ?
        (status=:FAIL, reason=:t11_numerical_failure,
         result_sha256=_hash_text("structured-evaluator-t11-terminal-failure")) :
        T12.infer_structured_chains(root;
            edge_observations=joinpath(edge_dir, "edge_observations.tsv"),
            node_segments=joinpath(edge_dir, "node_segments.tsv"),
            todo10_receipt=joinpath(edge_dir, "receipt.toml"),
            fitted_edge_models=joinpath(admission_dir, "fitted_edge_models.tsv"),
            fitted_unary_nodes=joinpath(admission_dir, "fitted_unary_nodes.tsv"),
            fitted_edge_transforms=joinpath(admission_dir, "fitted_edge_transforms.tsv"),
            todo11_receipt=joinpath(admission_dir, "receipt.toml"))
    _verify_context(context)
    input = _convert_production_input(root, options, admission, data, normalized, chain_report)
    input = EvaluationInput(authority.authority_sha256, input.universe_sha256,
                            input.scan_dates, input.folds)
    for snapshot in authority.snapshots
        _verify_snapshot(snapshot, "authority snapshot")
    end
    _verify_context(context)
    manifest = _manifest_digest(root, context.snapshots, context.inventories,
                                vcat(["authority=$(input.authority_sha256)",
                                      "universe=$(input.universe_sha256)"],
                                     _input_manifest_records(input));
                                control=context.activation === nothing ? nothing :
                                        context.activation.control)
    context = _ProductionContext(context.root, context.snapshots, context.inventories,
                                 manifest, context.activation)
    _verify_context(context)
    for snapshot in authority.snapshots
        _verify_snapshot(snapshot, "authority snapshot before return")
    end
    _verify_context(context)
    _PRODUCTION_CONTEXT[input] = context
    _PRODUCTION_SELECTION_METADATA[input] = (
        admission=admission, chain_report=chain_report, data=data,
        normalized=normalized,
        t12_result_sha256=hasproperty(chain_report, :result_sha256) ?
                          String(chain_report.result_sha256) : "",
        t12_provenance_sha256=hasproperty(chain_report, :provenance_sha256) ?
                              String(chain_report.provenance_sha256) : "",
    )
    _ = contract
    return input
end

function _selection_outer_reference(fold::OuterFoldEvidence, model::String)
    graph = get(fold.graphs, model, nothing)
    graph === nothing &&
        _fail(:BLOCKED, :selection_reference_mismatch,
              "outer-score graph identity is absent for $model/$(fold.outer_date)")
    fit = graph.t11_fit_sha256
    unary = model == "C1" ? fold.c1_fit_sha256 : fold.c2_fit_sha256
    graph.selected_unary_fit_sha256 == unary ||
        _fail(:BLOCKED, :selection_reference_mismatch,
              "outer-score unary identity differs for $model/$(fold.outer_date)")
    return _selection_reference(model, "partition", "outer_score", fold.outer_date,
                                fold.outer_date, fit, unary)
end

function _selection_effective_outer_model(fold::OuterFoldEvidence, model::String)
    fold.c1_status == :FAIL &&
        _fail(:FAIL, :unary_numerical_failure, "outer C1 failed")
    fold.c1_status == :BLOCKED &&
        _fail(:BLOCKED, :authority_or_factor_mismatch, "outer C1 blocked")
    fold.c1_status == :PASS || return nothing
    model == "C2" || return "C1"
    fold.c2_status == :FAIL &&
        _fail(:FAIL, :unary_numerical_failure, "outer C2 failed")
    fold.c2_status == :BLOCKED &&
        _fail(:BLOCKED, :authority_or_factor_mismatch, "outer C2 blocked")
    fold.c2_status == :SKIPPED && return "C1"
    any(_unary(node, :C2).c2_unavailable for node in fold.nodes) && return "C1"
    return "C2"
end

function _selection_full_reference(result, model::String)
    fit = result.full_fit
    return _selection_reference(model, "full_refit", "full_refit", "NA", "NA",
                                fit.fit_sha256, fit.unary.fit_sha256)
end

function _selection_full_unary_outcome(result)
    status = _status(String(result.full_fit.unary.status), "full unary selection")
    reason = Symbol(result.full_fit.unary.reason)
    return status, reason
end

function _selection_rebind_decision(decision::UnarySelectionDecision,
                                    c1_reference, c2_reference, selected_reference;
                                    selected_model=decision.selected_model,
                                    reason=decision.reason)
    selected_model = decision.status == :PASS ? selected_model : ""
    selected = decision.status == :PASS ? selected_reference : nothing
    return _selection_decision(
        decision.scope, decision.outer_date, decision.status, reason,
        selected_model, selected, c1_reference, c2_reference,
        decision.fold_count, decision.scan_count, decision.bootstrap_count,
        decision.bootstrap_lower, decision.every_date_positive,
        decision.decision_hash, decision.evidence_hash,
        decision.upstream_reason,
    )
end

function _selection_final_partition(result, date::String)
    matches = [partition for partition in result.final_gate.partitions
               if partition.target_date == date]
    length(matches) == 1 ||
        _fail(:BLOCKED, :selection_reference_mismatch,
              "final_lodo partition cardinality differs for $(result.model_id)/$date")
    partition = only(matches)
    partition.fit.model_id == result.model_id ||
        _fail(:BLOCKED, :selection_reference_mismatch,
              "final_lodo partition model identity differs")
    return partition
end

function _selection_final_folds(input::EvaluationInput, metadata)
    models = metadata.admission.models
    haskey(models, "C1") && haskey(models, "C2") ||
        _fail(:BLOCKED, :selection_reference_mismatch,
              "final_lodo model identities are incomplete")
    all_dates = sort(unique(last.(input.scan_dates)))
    final_folds = InnerFoldEvidence[]
    final_c1 = UnarySelectionReference[]
    final_c2 = UnarySelectionReference[]
    final_full_c1 = _selection_full_reference(models["C1"], "C1")
    final_full_c2 = _selection_full_reference(models["C2"], "C2")
    logical_ids = Set{Tuple{String,String,String,String,String}}()
    for (model, result) in (("C1", models["C1"]), ("C2", models["C2"]))
        result.full_fit.model_id == model ||
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "full_refit model identity differs")
        sort(String.(result.full_fit.training_dates)) == all_dates &&
            result.full_fit.unary.training_dates == all_dates ||
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "full_refit training identity differs")
    end
    for date in all_dates
        c1 = _selection_final_partition(models["C1"], date)
        c2 = _selection_final_partition(models["C2"], date)
        expected_training = [value for value in all_dates if value != date]
        sort(String.(c1.training_dates)) == expected_training &&
            sort(String.(c2.training_dates)) == expected_training ||
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "final_lodo training dates are not the exact complement")
        sort(String.(c1.fit.training_dates)) == expected_training &&
            sort(String.(c2.fit.training_dates)) == expected_training ||
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "final unary training dates are not the exact complement")
        c1.fit.unary.training_dates == expected_training &&
            c2.fit.unary.training_dates == expected_training ||
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "final unary identities have mismatched training dates")
        c1.target_date == date && c2.target_date == date ||
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "final_lodo target identity differs")
        r1 = _selection_reference("C1", "partition", "final_lodo", "NA", date,
                                  c1.fit.fit_sha256, c1.fit.unary.fit_sha256)
        r2 = _selection_reference("C2", "partition", "final_lodo", "NA", date,
                                  c2.fit.fit_sha256, c2.fit.unary.fit_sha256)
        for reference in (r1, r2)
            identity = (reference.model, reference.fit_role, reference.scope,
                        reference.outer_date, reference.target_date)
            identity in logical_ids &&
                _fail(:BLOCKED, :selection_reference_mismatch,
                      "final reference identity is duplicated")
            push!(logical_ids, identity)
        end
        push!(final_c1, r1)
        push!(final_c2, r2)
        push!(final_folds, InnerFoldEvidence(
            date, expected_training,
            _status(c1.fit.unary.status, "final C1 unary"), c1.fit.unary.reason,
            c1.fit.unary.fit_sha256,
            _status(c2.fit.unary.status, "final C2 unary"), c2.fit.unary.reason,
            c2.fit.unary.fit_sha256,
            _node_evidence(metadata.normalized, c1.fit, c2.fit, date),
            _edge_counts(c1.fit, date), c1.fit.fit_sha256, c2.fit.fit_sha256,
        ))
    end
    full_ids = Set{Tuple{String,String,String,String,String}}()
    for reference in (final_full_c1, final_full_c2)
        identity = (reference.model, reference.fit_role, reference.scope,
                    reference.outer_date, reference.target_date)
        identity in full_ids &&
            _fail(:BLOCKED, :selection_reference_mismatch,
                  "full_refit reference identity is duplicated")
        push!(full_ids, identity)
    end
    length(final_c1) == length(all_dates) && length(final_c2) == length(all_dates) ||
        _fail(:BLOCKED, :selection_reference_mismatch,
              "final_lodo reference coverage differs")
    return final_folds, final_c1, final_c2, final_full_c1, final_full_c2
end

function _selection_binding_digest(context::_ProductionContext, path::String,
                                   label::String="binding")
    matches = [snapshot for snapshot in context.snapshots if snapshot.path == path]
    length(matches) == 1 ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence,
              "$label snapshot is absent or ambiguous: $path")
    return only(matches).sha256
end

function _selection_bindings(root::String, options::Dict{String,String}, input::EvaluationInput)
    context = get(_PRODUCTION_CONTEXT, input, nothing)
    context === nothing &&
        _fail(:BLOCKED, :incomplete_or_stale_evidence,
              "unary selection production context is absent")
    metadata = get(_PRODUCTION_SELECTION_METADATA, input, nothing)
    metadata === nothing &&
        _fail(:BLOCKED, :incomplete_or_stale_evidence,
              "unary selection provenance metadata is absent")
    file_digest(flag) = _selection_binding_digest(
        context, _repo_path(root, options[flag], flag), flag)
    directory_digest(flag, key) = begin
        directory = _directory_path(root, options[flag], flag)
        inventories = [(normpath(joinpath(context.root, relative)), expected)
                       for (relative, expected) in context.inventories
                       if normpath(joinpath(context.root, relative)) == directory]
        length(inventories) == 1 ||
            _fail(:BLOCKED, :incomplete_or_stale_evidence,
                  "$flag directory inventory is absent or ambiguous")
        expected = Set(joinpath(directory, String(name)) for name in only(inventories)[2])
        members = [snapshot for snapshot in context.snapshots
                   if startswith(snapshot.path, directory * "/")]
        Set(snapshot.path for snapshot in members) == expected ||
            _fail(:BLOCKED, :incomplete_or_stale_evidence,
                  "$flag directory snapshot membership is incomplete")
        records = String["schema=structured-unary-selection-directory-v1", key]
        append!(records, "$(relpath(snapshot.path, directory))\t$(snapshot.sha256)"
                for snapshot in sort(members; by=snapshot -> snapshot.path))
        _length_digest(records)
    end
    return (
        evaluator_config_sha256=file_digest("--evaluator-config"),
        evaluator_source_sha256=_selection_binding_digest(
            context, joinpath(root, "test", "evaluate_structured_unit_assignment.jl"),
            "evaluator source"),
        features_sha256=file_digest("--features"),
        candidate_config_sha256=file_digest("--candidate-config"),
        model_config_sha256=file_digest("--model-config"),
        universe_dir_sha256=directory_digest("--universe-dir", "universe_dir_sha256"),
        edge_dir_sha256=directory_digest("--edge-dir", "edge_dir_sha256"),
        forward_receipt_sha256=file_digest("--forward-receipt"),
        backward_receipt_sha256=file_digest("--backward-receipt"),
        admission_dir_sha256=directory_digest("--admission-dir", "admission_dir_sha256"),
        authority_sha256=input.authority_sha256,
        universe_sha256=input.universe_sha256,
        provenance_manifest_sha256=context.manifest_sha256,
        t12_result_sha256=metadata.t12_result_sha256,
        t12_provenance_sha256=metadata.t12_provenance_sha256,
        todo12_result_sha256=metadata.t12_result_sha256,
        todo12_provenance_sha256=metadata.t12_provenance_sha256,
    )
end

function _selection_semantic_path_guard(value::AbstractString, context::String)
    text = lowercase(replace(String(value), '\\' => '/'))
    forbidden = ("label", "truth", "expected", "benchmark", "groundtruth",
                 "ground-truth", "n_selected", "n-selected", "count")
    any(token -> occursin(token, text), forbidden) &&
        _fail(:BLOCKED, :authority_path_invalid,
              "$context is not an allowed label-free unary input")
    return nothing
end

const _UNARY_SELECTION_STATUS_PRECEDENCE = Dict(
    :PASS => 0, :SKIPPED => 1, :FAIL => 2, :BLOCKED => 3)

function _selection_terminal_decision(outer_decisions::Vector{UnarySelectionDecision},
                                      final_decision::UnarySelectionDecision)
    candidates = [(decision, false) for decision in outer_decisions]
    push!(candidates, (final_decision, true))
    sort!(candidates; by=pair -> (-get(_UNARY_SELECTION_STATUS_PRECEDENCE,
                                      pair[1].status, 3), pair[2] ? 0 : 1,
                                   pair[1].outer_date, String(pair[1].reason)))
    decision = first(candidates)[1]
    return (status=decision.status, reason=decision.reason,
            upstream_reason=decision.upstream_reason)
end

function _selection_top_hash(outer_decisions::Vector{UnarySelectionDecision},
                             final_decision::UnarySelectionDecision,
                             bindings)
    lines = String["schema=structured-unary-selection-top-v3",
                   "final_status=$(final_decision.status)",
                   "final_reason=$(final_decision.reason)",
                   "final_upstream=$(final_decision.upstream_reason)",
                   "final_model=$(final_decision.selected_model)",
                   "final_reference=$(_selection_reference_identity(final_decision.selected_reference))",
                   "final_c1_reference=$(_selection_reference_identity(final_decision.c1_reference))",
                   "final_c2_reference=$(_selection_reference_identity(final_decision.c2_reference))",
                   "final_decision_hash=$(final_decision.decision_hash)",
                   "final_evidence_hash=$(final_decision.evidence_hash)"]
    append!(lines, "outer\t$(decision.outer_date)\t$(decision.status)\t$(decision.reason)\t$(decision.selected_model)\t$(decision.decision_hash)\t$(decision.evidence_hash)\t$(_selection_reference_identity(decision.selected_reference))\t$(_selection_reference_identity(decision.c1_reference))\t$(_selection_reference_identity(decision.c2_reference))"
            for decision in sort(copy(outer_decisions); by=decision -> decision.outer_date))
    append!(lines, "binding\t$(name)\t$(getfield(bindings, name))"
            for name in propertynames(bindings))
    return _length_digest(lines)
end

function _produce_unary_selection_evidence_from_context(input::EvaluationInput,
                                                        metadata, bindings)
    outer_decisions = UnarySelectionDecision[]
    outer_folds = UnarySelectionFoldEvidence[]
    outer_scans = UnarySelectionScanEvidence[]
    outer_dates = UnarySelectionDateEvidence[]
    outer_bootstrap = UnarySelectionBootstrapEvidence[]
    for fold in sort(copy(input.folds); by=fold -> fold.outer_date)
        c1 = _selection_outer_reference(fold, "C1")
        c2 = _selection_outer_reference(fold, "C2")
        result = _inner_gate_evidence(fold.inner_folds;
            decision_scope="outer", outer_date=fold.outer_date,
            c1_reference=c1, c2_reference=c2)
        outer_decision = result.decision
        if outer_decision.status == :PASS
            effective_model = _selection_effective_outer_model(
                fold, outer_decision.selected_model)
            if effective_model === nothing
                outer_decision = _selection_decision(
                    "outer", fold.outer_date, :SKIPPED, :scope_unavailable,
                    "", nothing, c1, c2, outer_decision.fold_count,
                    outer_decision.scan_count, outer_decision.bootstrap_count,
                    outer_decision.bootstrap_lower,
                    outer_decision.every_date_positive, "",
                    outer_decision.evidence_hash, outer_decision.upstream_reason)
            else
                effective_reference = effective_model == "C1" ? c1 : c2
                outer_decision = _selection_rebind_decision(
                    outer_decision, c1, c2, effective_reference;
                    selected_model=effective_model,
                    reason=effective_model == "C1" && outer_decision.selected_model == "C2" ?
                           :c2_unavailable_fallback_c1 : outer_decision.reason)
            end
        end
        push!(outer_decisions, outer_decision)
        append!(outer_folds, result.folds)
        append!(outer_scans, result.scans)
        append!(outer_dates, result.dates)
        append!(outer_bootstrap, result.bootstrap)
    end

    final_folds = UnarySelectionFoldEvidence[]
    final_scans = UnarySelectionScanEvidence[]
    final_dates = UnarySelectionDateEvidence[]
    final_bootstrap = UnarySelectionBootstrapEvidence[]
    final_decision = try
        folds, c1_refs, c2_refs, full_c1, full_c2 =
            _selection_final_folds(input, metadata)
        result = _inner_gate_evidence(folds;
            decision_scope="final", outer_date="NA",
            c1_reference=full_c1, c2_reference=full_c2)
        final_folds = result.folds
        final_scans = result.scans
        final_dates = result.dates
        final_bootstrap = result.bootstrap
        if result.decision.status != :PASS
            _selection_decision("final", "NA", result.decision.status,
                               result.decision.reason, "", nothing,
                               full_c1, full_c2, result.decision.fold_count,
                               result.decision.scan_count,
                               result.decision.bootstrap_count,
                               result.decision.bootstrap_lower,
                               result.decision.every_date_positive, "",
                               result.decision.evidence_hash,
                               result.decision.upstream_reason)
        else
        selected_model = result.decision.selected_model
        selected_result = metadata.admission.models[selected_model]
        selected_status, selected_reason = _selection_full_unary_outcome(selected_result)
        if selected_status == :PASS
            selected_full = selected_model == "C2" ? full_c2 : full_c1
            _selection_rebind_decision(result.decision, full_c1, full_c2, selected_full;
                                       selected_model=selected_model)
        else
            executable_reason = selected_status == :SKIPPED ? :scope_unavailable :
                                selected_status == :FAIL ? :unary_numerical_failure :
                                :authority_or_factor_mismatch
            _selection_decision("final", "NA", selected_status, executable_reason,
                               "", nothing, full_c1, full_c2,
                               result.decision.fold_count, result.decision.scan_count,
                               result.decision.bootstrap_count,
                               result.decision.bootstrap_lower,
                               result.decision.every_date_positive, "",
                               result.decision.evidence_hash, selected_reason)
        end
        end
    catch error
        error isa EvaluatorError || rethrow()
        _selection_decision("final", "NA", error.status, error.reason, "", nothing,
                           nothing, nothing, 0, 0, 0, NaN, false, "", "",
                           error.reason)
    end
    all_folds = sort(vcat(outer_folds, final_folds);
                     by=fold -> (fold.decision_scope, fold.outer_date, fold.heldout_date))
    all_scans = sort(vcat(outer_scans, final_scans);
                     by=scan -> (scan.decision_scope, scan.outer_date, scan.date, scan.file))
    all_dates = sort(vcat(outer_dates, final_dates);
                     by=date -> (date.decision_scope, date.outer_date, date.date))
    all_bootstrap = sort(vcat(outer_bootstrap, final_bootstrap);
                         by=row -> (row.scope, row.outer_date, row.seed))
    canonical_outer = sort(outer_decisions; by=decision -> decision.outer_date)
    terminal = _selection_terminal_decision(canonical_outer, final_decision)
    evidence_hash = _selection_top_hash(canonical_outer, final_decision, bindings)
    return UnarySelectionEvidence(terminal.status, terminal.reason,
                                  terminal.upstream_reason, canonical_outer,
                                  final_decision, all_folds, all_scans, all_dates,
                                  all_bootstrap, bindings, input.authority_sha256,
                                  evidence_hash)
end

"""Produce immutable, label-free unary model-selection evidence.

This seam intentionally has no output path and no selector, threshold, score,
matrix, label, prior, or seed argument.  All source and input authority checks
are performed by the existing production loader before evidence is assembled.
"""
function produce_unary_selection_evidence(root::AbstractString;
        evaluator_config::AbstractString,
        features::AbstractString,
        candidate_config::AbstractString,
        model_config::AbstractString,
        universe_dir::AbstractString,
        edge_dir::AbstractString,
        forward_receipt::AbstractString,
        backward_receipt::AbstractString,
        admission_dir::AbstractString)::UnarySelectionEvidence
    root_text = String(root)
    cli_pairs = (
        ("--root", root_text),
        ("--evaluator-config", String(evaluator_config)),
        ("--features", String(features)),
        ("--candidate-config", String(candidate_config)),
        ("--model-config", String(model_config)),
        ("--universe-dir", String(universe_dir)),
        ("--edge-dir", String(edge_dir)),
        ("--forward-receipt", String(forward_receipt)),
        ("--backward-receipt", String(backward_receipt)),
        ("--admission-dir", String(admission_dir)),
    )
    isabspath(root_text) ||
        _fail(:BLOCKED, :authority_path_invalid, "root must be absolute")
    isdir(root_text) ||
        _fail(:BLOCKED, :authority_path_invalid, "root is not a directory")
    cli_arguments = String[value for pair in cli_pairs for value in pair]
    cli_options = _selection_cli_options(cli_arguments)
    root_text = cli_options["--root"]
    root_path = normpath(root_text)
    isdir(root_path) ||
        _fail(:BLOCKED, :authority_path_invalid, "root is not a directory")
    realpath(root_path) == root_path ||
        _fail(:BLOCKED, :authority_path_invalid, "root is not canonical")
    evaluator_config = cli_options["--evaluator-config"]
    features = cli_options["--features"]
    candidate_config = cli_options["--candidate-config"]
    model_config = cli_options["--model-config"]
    universe_dir = cli_options["--universe-dir"]
    edge_dir = cli_options["--edge-dir"]
    forward_receipt = cli_options["--forward-receipt"]
    backward_receipt = cli_options["--backward-receipt"]
    admission_dir = cli_options["--admission-dir"]
    values = (evaluator_config, features, candidate_config, model_config,
              universe_dir, edge_dir, forward_receipt, backward_receipt,
              admission_dir)
    labels = ("evaluator config", "features", "candidate config", "model config",
              "universe", "edge", "forward receipt", "backward receipt",
              "admission")
    for (value, label) in zip(values, labels)
        isempty(String(value)) &&
            _fail(:BLOCKED, :authority_path_invalid, "$label is empty")
        isabspath(String(value)) &&
            _fail(:BLOCKED, :authority_path_invalid, "$label must be repository-relative")
        _selection_semantic_path_guard(value, label)
    end
    options = Dict{String,String}(
        "--root" => root_path,
        "--evaluator-config" => String(evaluator_config),
        "--features" => String(features),
        "--candidate-config" => String(candidate_config),
        "--model-config" => String(model_config),
        "--universe-dir" => String(universe_dir),
        "--edge-dir" => String(edge_dir),
        "--forward-receipt" => String(forward_receipt),
        "--backward-receipt" => String(backward_receipt),
        "--admission-dir" => String(admission_dir),
    )
    input = _load_production_input(options)
    context = _PRODUCTION_CONTEXT[input]
    metadata = get(_PRODUCTION_SELECTION_METADATA, input, nothing)
    metadata === nothing &&
        _fail(:BLOCKED, :incomplete_or_stale_evidence,
              "unary selection loader metadata is absent")
    bindings = _selection_bindings(root_path, options, input)
    evidence = _produce_unary_selection_evidence_from_context(input, metadata, bindings)
    _verify_context(context)
    return evidence
end

const _VALUE_FLAGS = Set([
    "--root", "--evaluator-config", "--features", "--candidate-config",
    "--model-config", "--universe-dir", "--edge-dir", "--forward-receipt",
    "--backward-receipt", "--admission-dir", "--out-dir",
])
const _ACTIVATION_VALUE_FLAGS = Set(["--activation-receipt", "--activation-sha256"])
const _CLI_FLAGS = union(_VALUE_FLAGS, _ACTIVATION_VALUE_FLAGS, Set(["--help"]))

function _clean_cli_value(flag::AbstractString, value::String)::String
    all(isascii, value) && !occursin('%', value) && !occursin('\\', value) ||
        _fail(:BLOCKED, :cli_error, "noncanonical encoded path value for $flag")
    return value
end

function _parse_cli_options(arguments::Vector{String}; require_output::Bool=true)::Dict{String,String}
    try
        T8.StructuredFirewall.validate_cli_arguments(
            arguments; allowed_flags=_CLI_FLAGS,
            value_flags=union(_VALUE_FLAGS, _ACTIVATION_VALUE_FLAGS),
            repeatable_flags=Set{String}())
    catch error
        error isa T8.StructuredFirewall.StructuredFirewallError || rethrow()
        _fail(:BLOCKED, :cli_error, sprint(showerror, error))
    end
    options = Dict{String,String}()
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        startswith(argument, "--") ||
            _fail(:BLOCKED, :cli_error, "positional arguments are not accepted")
        pair = split(argument, '='; limit=2)
        flag = pair[1]
        flag in _CLI_FLAGS || _fail(:BLOCKED, :cli_error, "unknown option $flag")
        haskey(options, flag) && _fail(:BLOCKED, :cli_error, "duplicate option $flag")
        if flag == "--help"
            length(pair) == 1 || _fail(:BLOCKED, :cli_error, "--help is standalone")
            options[flag] = "true"
            index += 1
        elseif length(pair) == 2
            isempty(pair[2]) && _fail(:BLOCKED, :cli_error, "empty value for $flag")
            options[flag] = _clean_cli_value(flag, pair[2])
            index += 1
        else
            index < length(arguments) ||
                _fail(:BLOCKED, :cli_error, "missing value for $flag")
            next = arguments[index + 1]
            startswith(next, "--") &&
                _fail(:BLOCKED, :cli_error, "missing value for $flag")
            options[flag] = _clean_cli_value(flag, next)
            index += 2
        end
    end
    haskey(options, "--help") && length(options) == 1 ||
        (haskey(options, "--help") &&
         _fail(:BLOCKED, :cli_error, "--help accepts no evaluation inputs"))
    haskey(options, "--help") && return options
    required_flags = require_output ? _VALUE_FLAGS : setdiff(_VALUE_FLAGS, Set(["--out-dir"]))
    for flag in sort(collect(required_flags))
        haskey(options, flag) || _fail(:BLOCKED, :cli_error, "missing option $flag")
    end
    haskey(options, "--activation-receipt") == haskey(options, "--activation-sha256") ||
        _fail(:BLOCKED, :cli_error,
              "--activation-receipt and --activation-sha256 must be supplied together")
    return options
end

parse_cli(arguments::Vector{String})::Dict{String,String} =
    _parse_cli_options(arguments; require_output=true)

_selection_cli_options(arguments::Vector{String})::Dict{String,String} =
    _parse_cli_options(arguments; require_output=false)

function _usage(io::IO=stdout)
    println(io, "Usage: julia --project=. test/evaluate_structured_unit_assignment.jl \\")
    println(io, "  --root PATH --evaluator-config PATH --features PATH \\")
    println(io, "  --candidate-config PATH --model-config PATH --universe-dir PATH \\")
    println(io, "  --edge-dir PATH --forward-receipt PATH --backward-receipt PATH \\")
    println(io, "  --admission-dir PATH --out-dir PATH")
    println(io, "  [--activation-receipt REPO_RELATIVE_PATH --activation-sha256 EXPECTED_SHA256]")
end

struct _PublicationError <: Exception
    state::Symbol
    reason::Symbol
    message::String
    function _PublicationError(state::Symbol, reason::Symbol, message::AbstractString)
        state in (:not_committed, :ambiguous) ||
            throw(ArgumentError("invalid publication error state"))
        new(state, reason, String(message))
    end
end

Base.showerror(io::IO, error::_PublicationError) =
    print(io, "publication ", error.state, " [", error.reason, "]: ", error.message)

const _PUBLICATION_HOOK = Ref{Union{Nothing,Function}}(nothing)

function _publication_hook(event::Symbol, args...)
    hook = _PUBLICATION_HOOK[]
    hook === nothing && return nothing
    return hook(event, args...)
end

function _publication_identity(path::String)
    info = stat(path)
    return (UInt64(info.device), UInt64(info.inode))
end

function _publication_regular_file(path::String, context::String)
    isfile(path) && !islink(path) ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "$context is not a regular file"))
    info = stat(path)
    UInt64(info.nlink) == 1 ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "$context is hard-linked"))
    UInt(info.mode) & UInt(0o777) == UInt(0o644) ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "$context mode differs"))
    return nothing
end

function _verify_publication_bytes(path::String, expected::Vector{UInt8},
                                   context::String, state::Symbol)
    isfile(path) && !islink(path) ||
        throw(_publication_error(state, :publication_failure,
                                 "$context is not a regular file"))
    UInt(stat(path).mode) & UInt(0o777) == UInt(0o644) ||
        throw(_publication_error(state, :publication_failure,
                                 "$context mode differs"))
    before = _publication_identity(path)
    bytes = read(path)
    _publication_identity(path) == before ||
        throw(_publication_error(state, :publication_failure,
                                 "$context identity changed while reading"))
    bytes == expected ||
        throw(_publication_error(state, :publication_failure,
                                 "$context bytes differ"))
    return nothing
end

function _publication_directory(path::String, mode::UInt, context::String)
    isdir(path) && !islink(path) ||
        _fail(:BLOCKED, :publication_failure, "$context is not a real directory")
    UInt(stat(path).mode) & UInt(0o777) == mode ||
        _fail(:BLOCKED, :publication_failure, "$context mode differs")
    return nothing
end

function _publication_file_contract(files::Dict{String,Vector{UInt8}})
    names = String.(collect(keys(files)))
    all(name -> !isempty(name) && name == basename(name) &&
                 !occursin('\0', name) && !occursin('/', name) &&
                 !occursin('\\', name) && name ∉ (".", ".."), names) ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "publication contains an unsafe basename"))
    expected = Set(String.(collect(_FILES)))
    actual = Set(names)
    actual == expected || actual == Set(["receipt.toml"]) ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "publication file set is not an accepted product"))
    return sort(filter(!=("receipt.toml"), names))
end

function _publication_error(state::Symbol, reason::Symbol, message::AbstractString)
    state in (:not_committed, :ambiguous) ||
        throw(ArgumentError("invalid publication error state"))
    return _PublicationError(state, reason, String(message))
end

function _publication_destination(root::String, supplied::AbstractString, context::String)
    text = String(supplied)
    isempty(text) && _fail(:BLOCKED, :authority_path_invalid, "$context is empty")
    occursin('\0', text) && _fail(:BLOCKED, :authority_path_invalid, "$context contains NUL")
    any(==(".."), split(replace(text, '\\' => '/'), '/'; keepempty=true)) &&
        _fail(:BLOCKED, :authority_path_invalid, "$context contains parent traversal")
    destination = normpath(isabspath(text) ? text : joinpath(root, text))
    prefix = root == "/" ? "/" : root * "/"
    (destination == root || startswith(destination, prefix)) ||
        _fail(:BLOCKED, :authority_path_invalid, "$context escapes root")
    parent = dirname(destination)
    isdir(parent) || _fail(:BLOCKED, :authority_path_invalid, "$context parent is absent")
    _reject_link_components(parent)
    if ispath(destination)
        islink(destination) && _fail(:BLOCKED, :authority_symlink, "$context is a symlink")
        realpath(destination) == destination ||
            _fail(:BLOCKED, :authority_path_mismatch, "$context is not canonical")
    end
    return destination
end

function _fsync_directory(path::String; state::Symbol=:not_committed)
    _publication_hook(:before_fsync_directory, path, state)
    descriptor = ccall(:open, Cint, (Cstring, Cint), path, 0)
    descriptor >= 0 || throw(_publication_error(state, :publication_failure,
                                                 "cannot open directory for fsync"))
    try
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            throw(_publication_error(state, :publication_failure, "directory fsync failed"))
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
    _publication_hook(:after_fsync_directory, path, state)
end

function _write_exclusive(path::String, bytes::Vector{UInt8})
    flags = Cint(0x0001 | 0x0040 | 0x0080 | 0x80000) # WRONLY|CREAT|EXCL|CLOEXEC
    descriptor = ccall(:open, Cint, (Cstring, Cint, Cint), path, flags, Cint(0o644))
    descriptor >= 0 || throw(_publication_error(:not_committed, :publication_failure,
                                                 "exclusive stage creation failed"))
    committed = false
    try
        chmod(path, 0o644)
        offset = 0
        while offset < length(bytes)
            written = GC.@preserve bytes ccall(
                :write, Cssize_t, (Cint, Ptr{UInt8}, Csize_t), descriptor,
                pointer(bytes, offset + 1), length(bytes) - offset)
            written > 0 || throw(_publication_error(:not_committed, :publication_failure,
                                                     "staged write was partial"))
            offset += Int(written)
        end
        _publication_hook(:before_fsync_file, path)
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            throw(_publication_error(:not_committed, :publication_failure,
                                     "staged file fsync failed"))
        _publication_hook(:after_fsync_file, path)
        chmod(path, 0o644)
        committed = true
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
    committed || ispath(path) || nothing
    return nothing
end

function _rename_noreplace(source::String, destination::String)
    Sys.islinux() || throw(_publication_error(:not_committed, :publication_failure,
                                              "renameat2 is required"))
    forced_errno = _publication_hook(:rename_noreplace, source, destination)
    forced_errno === nothing || (forced_errno isa Integer) ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "rename hook returned an invalid errno"))
    result = try
        forced_errno === nothing ? ccall(:renameat2, Cint,
            (Cint, Cstring, Cint, Cstring, Cuint), Cint(-100), source,
            Cint(-100), destination, Cuint(1)) : Cint(-1)
    catch error
        throw(_publication_error(:not_committed, :publication_failure, sprint(showerror, error)))
    end
    result == 0 && return nothing
    errno = forced_errno === nothing ? Int(Base.Libc.errno()) : Int(forced_errno)
    source_intact = isdir(source) && !islink(source)
    destination_absent = !ispath(destination)
    state = source_intact && destination_absent ? :not_committed : :ambiguous
    errno == 17 && source_intact &&
        throw(_publication_error(:not_committed, :destination_exists,
                                 "destination appeared concurrently"))
    errno in (22, 38, 95) && source_intact && destination_absent &&
        throw(_publication_error(:not_committed, :rename_unsupported,
                                 "renameat2 failed with errno $errno"))
    throw(_publication_error(state, :publication_failure,
                             "renameat2 failed with errno $errno"))
end

function _same_publication(destination::String, files::Dict{String,Vector{UInt8}})
    try
        _publication_directory(destination, UInt(0o700), "published output")
        sort(readdir(destination)) == sort(collect(keys(files))) || return false
        identities = Set{Tuple{UInt64,UInt64}}()
        for name in sort(collect(keys(files)))
            path = joinpath(destination, name)
            _publication_regular_file(path, "published $name")
            identity = _publication_identity(path)
            identity in identities && return false
            push!(identities, identity)
            bytes = read(path)
            _publication_identity(path) == identity || return false
            bytes == files[name] || return false
        end
        return true
    catch
        return false
    end
end

function _verify_publication_identity(destination::String,
                                      identity::Tuple{UInt64,UInt64})
    isdir(destination) && !islink(destination) ||
        throw(_publication_error(:ambiguous, :publication_failure,
                                 "reserved destination is not a real directory"))
    UInt(stat(destination).mode) & UInt(0o777) == UInt(0o700) ||
        throw(_publication_error(:ambiguous, :publication_failure,
                                 "reserved destination mode differs"))
    _publication_identity(destination) == identity ||
        throw(_publication_error(:ambiguous, :publication_failure,
                                 "reserved destination identity changed"))
    return nothing
end

function _verify_stage_identity(stage::String, identity::Tuple{UInt64,UInt64})
    isdir(stage) && !islink(stage) ||
        throw(_publication_error(:ambiguous, :publication_failure,
                                 "private publication stage is not a real directory"))
    UInt(stat(stage).mode) & UInt(0o777) == UInt(0o700) ||
        throw(_publication_error(:not_committed, :publication_failure,
                                 "private publication stage mode differs"))
    _publication_identity(stage) == identity ||
        throw(_publication_error(:ambiguous, :publication_failure,
                                 "private publication stage identity changed"))
    return nothing
end

function _verify_publication_inventory(directory::String, names::Vector{String},
                                       context::String)
    sort(readdir(directory)) == sort(names) ||
        throw(_publication_error(:ambiguous, :publication_failure,
                                 "$context inventory differs"))
    for name in names
        _publication_regular_file(joinpath(directory, name), "$context/$name")
    end
    return nothing
end

function _cleanup_private_stage(stage::String,
                                identity::Tuple{UInt64,UInt64})
    ispath(stage) || return nothing
    _publication_hook(:before_stage_cleanup, stage, identity)
    isdir(stage) && !islink(stage) || return nothing
    _publication_identity(stage) == identity || return nothing
    rm(stage; recursive=true, force=true)
    return nothing
end

function _link_noreplace(source::String, destination::String, state::Symbol)
    _publication_hook(:before_hardlink, source, destination, state)
    try
        hardlink(source, destination)
    catch error
        ispath(destination) &&
            throw(_publication_error(state, :destination_exists,
                                     "publication destination member appeared"))
        throw(_publication_error(state, :publication_failure, sprint(showerror, error)))
    end
    return nothing
end

function _verify_context_for_publication(context::Union{Nothing,_ProductionContext},
                                         state::Symbol)
    context === nothing && return nothing
    try
        _verify_context(context)
    catch error
        error isa _PublicationError && rethrow()
        reason = error isa EvaluatorError ? error.reason : :publication_failure
        throw(_publication_error(state, reason, sprint(showerror, error)))
    end
    return nothing
end

function _publish_fallback(destination::String, parent::String, stage::String,
                           stage_identity::Tuple{UInt64,UInt64},
                           files::Dict{String,Vector{UInt8}},
                           names::Vector{String},
                           context::Union{Nothing,_ProductionContext}=nothing)
    state = :not_committed
    try
        _verify_stage_identity(stage, stage_identity)
        _verify_publication_inventory(stage, vcat(names, ["receipt.toml"]), "staged output")
        !ispath(destination) ||
            _fail(:BLOCKED, :publication_collision, "destination appeared before reservation")
        _publication_hook(:before_reserve, stage, destination, parent)
        try
            mkdir(destination; mode=0o700)
        catch error
            if ispath(destination)
                _same_publication(destination, files) &&
                    (_verify_context_for_publication(context, :not_committed);
                     return :committed_verified)
                _fail(:BLOCKED, :publication_collision, "concurrent destination differs")
            end
            throw(_publication_error(:not_committed, :publication_failure,
                                     sprint(showerror, error)))
        end
        _publication_hook(:after_reserve_create, stage, destination, parent)
        chmod(destination, 0o700)
        destination_identity = _publication_identity(destination)
        _verify_publication_identity(destination, destination_identity)
        _fsync_directory(parent)
        _verify_publication_identity(destination, destination_identity)
        _publication_hook(:after_reserve, stage, destination, destination_identity)
        _verify_stage_identity(stage, stage_identity)

        for name in names
            _verify_publication_identity(destination, destination_identity)
            source = joinpath(stage, name)
            target = joinpath(destination, name)
            _publication_regular_file(source, "staged $name")
            _verify_publication_bytes(source, files[name], "staged $name", :not_committed)
            _link_noreplace(source, target, :not_committed)
            _verify_publication_bytes(target, files[name], "published $name", :not_committed)
            _publication_hook(:after_payload_link, name, destination, destination_identity)
            _verify_publication_identity(destination, destination_identity)
            _verify_publication_bytes(target, files[name], "published $name", :not_committed)
            rm(source; force=true)
            _publication_regular_file(target, "published $name")
            _verify_publication_bytes(target, files[name], "published $name", :not_committed)
            _verify_publication_identity(destination, destination_identity)
        end
        _fsync_directory(destination)
        _fsync_directory(stage)
        _verify_publication_identity(destination, destination_identity)
        _verify_stage_identity(stage, stage_identity)
        _verify_publication_inventory(destination, names, "precommit destination")
        _verify_publication_inventory(stage, ["receipt.toml"], "precommit stage")
        _verify_publication_bytes(joinpath(stage, "receipt.toml"),
                                  files["receipt.toml"], "staged receipt", :not_committed)
        _verify_context_for_publication(context, :not_committed)

        _publication_hook(:before_receipt_link, stage, destination, destination_identity)
        _verify_publication_identity(destination, destination_identity)
        _verify_context_for_publication(context, :not_committed)
        _verify_stage_identity(stage, stage_identity)
        _verify_publication_bytes(joinpath(stage, "receipt.toml"),
                                  files["receipt.toml"], "staged receipt", :not_committed)
        receipt_destination = joinpath(destination, "receipt.toml")
        try
            _link_noreplace(joinpath(stage, "receipt.toml"), receipt_destination,
                            :not_committed)
        catch error
            if error isa _PublicationError && error.reason == :destination_exists
                state = :ambiguous
                _publication_hook(:receipt_collision, destination, state)
                throw(_publication_error(:ambiguous, :destination_exists,
                                         "receipt appeared concurrently"))
            end
            rethrow()
        end
        state = :ambiguous
        _publication_hook(:after_receipt_link, destination, destination_identity)
        _verify_publication_bytes(receipt_destination, files["receipt.toml"],
                                  "published receipt", :ambiguous)
        _verify_publication_identity(destination, destination_identity)
        _fsync_directory(destination; state=:ambiguous)
        rm(joinpath(stage, "receipt.toml"); force=true)
        _verify_publication_identity(destination, destination_identity)
        _fsync_directory(destination; state=:ambiguous)
        _verify_stage_identity(stage, stage_identity)
        isempty(readdir(stage)) ||
            throw(_publication_error(:ambiguous, :publication_failure,
                                     "private publication stage is not empty"))
        rm(stage; recursive=true, force=true)
        _fsync_directory(parent; state=:ambiguous)
        _publication_hook(:before_final_verify, destination, destination_identity)
        _verify_publication_identity(destination, destination_identity)
        _same_publication(destination, files) ||
            throw(_publication_error(:ambiguous, :publication_failure,
                                     "published bytes failed replay"))
        _verify_context_for_publication(context, :ambiguous)
        return :committed_verified
    catch error
        error isa _PublicationError && rethrow()
        error isa EvaluatorError && rethrow()
        throw(_publication_error(state, :publication_failure, sprint(showerror, error)))
    finally
        ispath(stage) && _cleanup_private_stage(stage, stage_identity)
    end
end

function _publish_atomic(root::String, supplied::AbstractString,
                         files::Dict{String,Vector{UInt8}};
                         context::Union{Nothing,_ProductionContext}=nothing)
    names = _publication_file_contract(files)
    destination = _publication_destination(root, supplied, "output directory")
    if ispath(destination)
        _same_publication(destination, files) ||
            _fail(:BLOCKED, :publication_collision, "existing output bytes differ")
        _verify_context_for_publication(context, :not_committed)
        return :committed_verified
    end
    parent = dirname(destination)
    stage = try
        mktempdir(parent; prefix=".stmfit-t13-evaluator-phase1-", cleanup=false)
    catch error
        throw(_publication_error(:not_committed, :publication_failure, sprint(showerror, error)))
    end
    chmod(stage, 0o700)
    stage_identity = _publication_identity(stage)
    state = :not_committed
    handed_to_fallback = false
    try
        for name in vcat(names, ["receipt.toml"])
            _write_exclusive(joinpath(stage, name), files[name])
        end
        _fsync_directory(stage)
        _verify_stage_identity(stage, stage_identity)
        _verify_context_for_publication(context, :not_committed)
        try
            _rename_noreplace(stage, destination)
        catch error
            if error isa _PublicationError && error.reason == :rename_unsupported
                _verify_stage_identity(stage, stage_identity)
                !ispath(destination) ||
                    _fail(:BLOCKED, :publication_collision,
                          "destination appeared during rename fallback")
                handed_to_fallback = true
                return _publish_fallback(destination, parent, stage, stage_identity,
                                         files, names, context)
            elseif error isa _PublicationError && error.reason == :destination_exists
                if _same_publication(destination, files)
                    # Preferred-rename destination-exists fast path: matching
                    # bytes alone do not prove the reviewed context still
                    # holds.  Revalidate before returning committed_verified so
                    # a precommit Main-control invalidation cannot be reported
                    # as verified publication.
                    _verify_context_for_publication(context, :not_committed)
                    return :committed_verified
                end
                _fail(:BLOCKED, :publication_collision, "concurrent destination differs")
            end
            error isa _PublicationError && (state = error.state; rethrow())
            throw(_publication_error(state, :publication_failure, sprint(showerror, error)))
        end
        state = :ambiguous
        _publication_hook(:after_preferred_rename, destination, state)
        destination_identity = _publication_identity(destination)
        destination_identity == stage_identity ||
            throw(_publication_error(:ambiguous, :publication_failure,
                                     "preferred publication identity differs"))
        _fsync_directory(parent; state=:ambiguous)
        _same_publication(destination, files) ||
            throw(_publication_error(:ambiguous, :publication_failure,
                                     "published bytes failed replay"))
        _verify_context_for_publication(context, :ambiguous)
        return :committed_verified
    catch error
        error isa _PublicationError && (state = error.state; rethrow())
        error isa EvaluatorError && rethrow()
        throw(_publication_error(state, :publication_failure, sprint(showerror, error)))
    finally
        !handed_to_fallback && state == :not_committed &&
            _cleanup_private_stage(stage, stage_identity)
    end
end

function _context_evaluator_config_sha256(context::_ProductionContext)
    historical_path = normpath(joinpath(context.root, _HISTORICAL_CONFIG_PATH))
    runtime_path = normpath(joinpath(context.root, _RUNTIME_V1_CONFIG_PATH))
    runtime_v2_path = normpath(joinpath(context.root, _RUNTIME_V2_CONFIG_PATH))
    runtime_v3_path = normpath(joinpath(context.root, _RUNTIME_V3_CONFIG_PATH))
    runtime_v4_path = normpath(joinpath(context.root, _RUNTIME_V4_CONFIG_PATH))
    expected = Dict(
        historical_path => _CONFIG_SHA256,
        runtime_path => _RUNTIME_V1_CONFIG_SHA256,
        runtime_v2_path => _RUNTIME_V2_CONFIG_SHA256,
        runtime_v3_path => _RUNTIME_V3_CONFIG_SHA256,
        runtime_v4_path => _RUNTIME_V4_CONFIG_SHA256,
    )
    matches = [snapshot for snapshot in context.snapshots
               if haskey(expected, snapshot.path)]
    length(matches) == 1 ||
        _fail(:BLOCKED, :incomplete_or_stale_evidence,
              "evaluator config snapshot is absent or ambiguous")
    snapshot = only(matches)
    _verify_snapshot(snapshot, "evaluator config snapshot")
    digest = _hash_bytes(snapshot.bytes)
    snapshot.sha256 == digest ||
        _fail(:BLOCKED, :authority_snapshot_changed,
              "evaluator config snapshot digest differs")
    digest == expected[snapshot.path] ||
        _fail(:BLOCKED, :authority_hash_mismatch,
              "evaluator config snapshot bytes differ")
    return digest
end

function _legacy_blocker_files(error::EvaluatorError,
                               evaluator_config_sha256::String)
    io = IOBuffer()
    println(io, "schema = ", repr(_SCHEMA * "_blocker_receipt_v2"))
    println(io, "schema_version = 2")
    println(io, "status = \"BLOCKED\"")
    println(io, "reason = ", repr(String(error.reason)))
    println(io, "message_sha256 = ", repr(_hash_text(error.message)))
    println(io, "evaluator_config_sha256 = ", repr(evaluator_config_sha256))
    println(io, "preclosure_plan_sha256 = ", repr(_PLAN_SHA256))
    println(io, "gateclosure_sha256 = ", repr(_GATECLOSURE_SHA256))
    println(io, "preclosure_boulder_sha256 = ", repr(_BOULDER_SHA256))
    println(io, "live_plan_boulder_runtime_authority = false")
    println(io, "review_sha256 = ", repr(_REVIEW_SHA256))
    println(io, "closure_v1_claim_sha256 = ", repr(_CLOSURE_V1_CLAIM_SHA256))
    println(io, "closure_v1_review_sha256 = ", repr(_CLOSURE_V1_REVIEW_SHA256))
    println(io, "closure_v1_root_manifest_sha256 = ", repr(_CLOSURE_V1_ROOT_MANIFEST_SHA256))
    println(io, "closure_v1_publication_receipt_sha256 = ",
            repr(_CLOSURE_V1_PUBLICATION_RECEIPT_SHA256))
    println(io, "closure_v1_publication_receipt_sidecar_sha256 = ",
            repr(_CLOSURE_V1_PUBLICATION_RECEIPT_SIDECAR_SHA256))
    println(io, "work_rows = 0")
    println(io, "bootstrap_rows = 0")
    println(io, "sign_rows = 0")
    return Dict("receipt.toml" => take!(io))
end

# Successor activated blocker receipt; the historical closure-v1 and GateClosure
# authority fields are intentionally absent because they do not certify an
# activated run.
function _activated_blocker_files(error::EvaluatorError,
                                  evaluator_config_sha256::String,
                                  activation::_ActivationBindings)
    io = IOBuffer()
    println(io, "schema = ", repr(_T13_ACTIVATED_BLOCKER_SCHEMA))
    println(io, "schema_version = 2")
    println(io, "status = \"BLOCKED\"")
    println(io, "reason = ", repr(String(error.reason)))
    println(io, "message_sha256 = ", repr(_hash_text(error.message)))
    println(io, "evaluator_config_sha256 = ", repr(evaluator_config_sha256))
    _t13_activation_receipt_lines(io, activation)
    println(io, "work_rows = 0")
    println(io, "bootstrap_rows = 0")
    println(io, "sign_rows = 0")
    return Dict("receipt.toml" => take!(io))
end

function _blocker_files(error::EvaluatorError;
                        context::Union{Nothing,_ProductionContext}=nothing)
    evaluator_config_sha256 = context === nothing ? _CONFIG_SHA256 :
                              _context_evaluator_config_sha256(context)
    activation = context === nothing ? nothing : context.activation
    activation === nothing &&
        return _legacy_blocker_files(error, evaluator_config_sha256)
    return _activated_blocker_files(error, evaluator_config_sha256, activation)
end

function main(arguments::Vector{String}=copy(ARGS))::Int
    options = try
        parse_cli(arguments)
    catch error
        error isa EvaluatorError || rethrow()
        showerror(stderr, error)
        println(stderr)
        return 2
    end
    haskey(options, "--help") && (_usage(); return 0)
    # The activation receipt/spec/scope checks happen before any production
    # context exists, so an invalid authorization cannot publish a blocker
    # receipt to the supplied destination.
    activation = nothing
    if haskey(options, "--activation-receipt")
        activation = try
            _authorize_t13_activation(options)
        catch error
            if error isa EvaluatorError
                showerror(stderr, error)
                println(stderr)
                return 2
            end
            rethrow()
        end
    end
    context = nothing
    try
        context = _preflight_config_context(options; activation=activation)
        input = _load_production_input(options; activation=activation)
        context = get(_PRODUCTION_CONTEXT, input, nothing)
        report = evaluate(input)
        if report.status == :BLOCKED
            blocked = EvaluatorError(:BLOCKED, report.reason, "evaluation blocked")
            _publish_atomic(options["--root"], options["--out-dir"],
                            _blocker_files(blocked; context=context); context=context)
            showerror(stderr, blocked)
            println(stderr)
            return 2
        end
        report.status in (:PASS, :SKIPPED, :FAIL) ||
            _fail(:BLOCKED, :publication_failure, "invalid terminal report status")
        context = get(_PRODUCTION_CONTEXT, input, nothing)
        context === nothing || _verify_context(context)
        files = report_files(report;
                             activation=context === nothing ? nothing :
                             context.activation)
        context === nothing || _verify_context(context)
        _publish_atomic(options["--root"], options["--out-dir"], files; context=context)
        println("status=", report.status)
        println("reason=", report.reason)
        println("result_sha256=", report.result_sha256)
        return 0
    catch error
        if error isa _PublicationError
            showerror(stderr, error)
            println(stderr)
            return error.state == :ambiguous ? 1 : 2
        end
        if error isa EvaluatorError
            if error.status == :BLOCKED
                context === nothing && (context = _historical_blocker_context(options))
                if context !== nothing
                    try
                        root = options["--root"]
                        destination = options["--out-dir"]
                        _publish_atomic(root, destination,
                                        _blocker_files(error; context=context);
                                        context=context)
                    catch publication_error
                        showerror(stderr, publication_error)
                        println(stderr)
                    end
                end
                showerror(stderr, error)
                println(stderr)
                return 2
            end
            showerror(stderr, error)
            println(stderr)
            return 1
        end
        rethrow()
    end
end

end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(StructuredUnitAssignmentEvaluator.main(copy(ARGS)))
end
