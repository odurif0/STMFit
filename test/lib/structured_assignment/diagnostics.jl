# allow: SIZE_OK — Todo 9's exclusive scope requires one self-contained diagnostics module.
module StructuredAssignmentDiagnostics

using LinearAlgebra
using Printf
using Random
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "..", "hierarchical_unit_assignment.jl"))
include(joinpath(@__DIR__, "robust_emissions.jl"))

const HUA = HierarchicalUnitAssignment
const SRE = StructuredRobustEmissions
const ROOT = realpath(joinpath(@__DIR__, "..", "..", ".."))
const MODEL_CONFIG_PATH = joinpath(ROOT, "config", "unit_assignment_structured_model.toml")
const MODEL_CONFIG_SHA256 =
    "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const UNARY_REVIEW_PATH = joinpath(
    ROOT, ".omo", "evidence", "structured-label-free-unit-assignment", "t8",
    "correction", "review", "AdversarialVerify.json")
const UNARY_REVIEW_SHA256 =
    "cd9218523df0c07fc4d7283b9c92d557d35233bcb2e46279d4fe03be4aad6144"
const CONFIRMED_DEPENDENCY_HASHES = Dict(
    "test/lib/structured_assignment/robust_emissions.jl" =>
        "0ca4a3863a4b9ed13aee5cd3a0229df4f2764a8524b764eeb0fd27238f14f241",
    "test/lib/structured_assignment/firewall.jl" =>
        "a189898d31759352d9dac0de8ff281d6001e1dc2db94df2cb725ff597537a9bb",
    "test/lib/structured_unit_assignment.jl" =>
        "8c74f21924381146d58dd9611d5a75578e7b691db3c16c847164549fd48e30b9",
    "test/build_structured_unit_predictions.jl" =>
        "976459ce93d51d6384a63a7d5d28c0a72e4dcc65f9d61500e8a1d314b52b6585",
    "test/lib/hierarchical_unit_assignment.jl" =>
        "4024d47b82077436a1c22a9511930effc5bca110f4c9eba5daef487b923513a8",
    "test/lib/hierarchical/loading.jl" =>
        "41171af123d489f8e3d4e20ba0ad2e55641a3cbb5768b3debc711c3f8262422a",
    "test/lib/hierarchical/nuisance.jl" =>
        "f256eb82966350544b55a23413417aa6eacda7cf6669c049f8dc90feb4fd2203",
    "test/lib/hierarchical/views.jl" =>
        "40380fb586a7866e718225eac3fd610e4eba9a0e0c3bddb924f108655e5fdfcf",
)
const DIAGNOSTIC_NAMES = ["mfa_q1", "scan_effects", "view_asymmetry"]
const CATEGORY_NAMES = ["00", "01", "10", "11"]
const TERMINAL_STATES = ["FOLLOWUP_REQUIRED", "SKIPPED", "BLOCKED"]
const SYNTHETIC_PAIR_PRIOR::NTuple{4,Float64} = (0.0, 0.0, 0.0, 0.0)
const REPORT_COLUMNS = [
    "record_type", "outer_heldout_date", "inner_date", "mechanism", "category",
    "status", "reason", "seed", "statistic", "threshold", "lower", "upper",
    "raw_p", "holm_p", "dimensions", "free_parameters", "inherited_parameters",
    "dates", "scans", "nodes", "edges", "effective_observations",
    "parameter_observation_ratio", "kish_ess", "support_scans",
    "scale_min_eigenvalue", "scale_condition_number", "config_hash", "feature_hash",
    "fold_hash", "input_receipt_hash", "unary_review_hash", "source_hash",
    "data_hash", "edge_hash", "detail_hash",
]

export DiagnosticError, FileSnapshot, DiagnosticPolicy, DiagnosticConfig, ResidualScan,
       PartitionEvidence, EdgeCategorySupport, FeasibilityInput, FeasibilityRow,
       DateFactorResult, MechanismResult, PartitionResult, HolmResult,
       NullCampaignRecord, NullCampaignSummary, FeatureDataset, FoldPlan,
       InputReceipt, DiagnosticRun, load_diagnostic_config, load_feature_dataset,
       load_fold_plan, load_input_receipt, run_diagnostics,
       first_component_statistic, scan_icc, holm_adjust, kish_ess,
       diagonal_two_component_parameter_count, edge_admission_parameter_count,
       scale_is_feasible, feasibility_table, clopper_pearson_upper,
       evaluate_partition, run_null_campaigns, null_campaign_hash, quantile_type7,
       model_source_hash, REPORT_COLUMNS, report_tsv_bytes, blocked_report_bytes,
       write_report

struct DiagnosticError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, error::DiagnosticError) =
    print(io, error.code, ": ", error.message)

struct FileSnapshot
    path::String
    bytes::Vector{UInt8}
    sha256::String
end

struct DiagnosticPolicy
    schema::String
    residual_model::String
    residual_formula::String
    residual_fallback::String
    residual_failure::String
    icc_estimator::String
    icc_date_centering::String
    icc_feature_pooling::String
    icc_unbalanced_group_size::String
    view_score::String
    view_contrast::String
    channel_drop_features::Vector{String}
    channel_drop_refit::String
    channel_drop_statistic::String
    channel_drop_same_sign::String
    quantile_method::String
    permutation_tail::String
    view_tail::String
    threshold_equality::String
    holm_order::Vector{String}
    holm_reject_equality::Bool
end

struct DiagnosticConfig
    path::String
    sha256::String
    view_names::Vector{String}
    view_features::Vector{Vector{String}}
    feature_names::Vector{String}
    permutation_count::Int
    permutation_quantile::Float64
    bootstrap_seeds::Vector{Int}
    bootstrap_lower_quantile::Float64
    bootstrap_upper_quantile::Float64
    interval_level::Float64
    alpha::Float64
    null_campaign_seeds::Vector{Int}
    null_dates::Int
    null_scans::Int
    null_chain_min::Int
    null_chain_max::Int
    clopper_confidence::Float64
    clopper_upper_max::Float64
    minimum_dates::Int
    minimum_scans::Int
    minimum_category_ess::Float64
    minimum_edges::Int
    minimum_support_scans::Int
    observations_per_parameter::Float64
    scale_eigenvalue_floor::Float64
    scale_condition_cap::Float64
    policy::DiagnosticPolicy
    snapshot::FileSnapshot
end

struct ResidualScan
    date::String
    scan::String
    values::Matrix{Float64}
end

struct PartitionEvidence
    outer_heldout_date::String
    residual_scans::Vector{ResidualScan}
    view_score_difference::Dict{String,Float64}
    dropout_difference::Dict{String,Float64}
end

struct EdgeCategorySupport
    name::String
    sum_weights::Float64
    sum_squared_weights::Float64
    support_scans_by_date::Dict{String,Int}
end

struct FeasibilityInput
    dates::Vector{String}
    scans::Int
    nodes::Int
    edges::Int
    categories::Vector{EdgeCategorySupport}
    scale_eigenvalues::Vector{Float64}
end

struct FeasibilityRow
    model::String
    dimensions::Vector{Int}
    free_parameters::Int
    inherited_parameters::Int
    dates::Int
    scans::Int
    nodes::Int
    edges::Int
    effective_observations::Int
    parameter_observation_ratio::Float64
    category_names::Vector{String}
    category_ess::Vector{Float64}
    category_support_scans::Vector{Int}
    minimum_category_ess::Float64
    minimum_support_scans::Int
    minimum_scale_eigenvalue::Float64
    scale_condition_number::Float64
    status::String
    reason::String
end

struct DateFactorResult
    date::String
    observed::Float64
    threshold::Float64
    excess::Float64
    raw_p::Float64
    permutation_values::Vector{Float64}
end

struct MechanismResult
    name::String
    status::String
    reason::String
    observed::Float64
    threshold::Float64
    lower::Float64
    upper::Float64
    raw_p::Float64
    holm_p::Float64
    date_results::Vector{DateFactorResult}
    permutation_values::Vector{Float64}
    bootstrap_values::Vector{Float64}
end

struct PartitionResult
    outer_heldout_date::String
    status::String
    mechanisms::Vector{MechanismResult}
    feasibility::Vector{FeasibilityRow}
    evidence_sha256::String
end

struct HolmResult
    names::Vector{String}
    raw_p::Vector{Float64}
    adjusted_p::Vector{Float64}
    rejected::Vector{Bool}
    sorted_indices::Vector{Int}
end

struct SyntheticNullScan
    date::String
    scan::String
    values::Matrix{Float64}
    amplitudes::Vector{Float64}
    edge_values::Matrix{Float64}
end

struct SyntheticNullCampaign
    seed::Int
    dates::Vector{String}
    scans::Vector{SyntheticNullScan}
    nodes::Int
    edges::Int
    latent_ones::Int
    data_sha256::String
    edge_sha256::String
end

struct SyntheticOuterDecision
    outer_heldout_date::String
    training_sha256::String
    inner_decision_sha256::String
    refit_sha256::String
    unary_model::String
    edge_admitted::Bool
    graph_enabled::Bool
    fit_invocations::Int
    scan_scores::Vector{Float64}
    mean_score::Float64
    bootstrap_lower::Float64
    status::String
    evidence_sha256::String
end

struct SyntheticMetaDecision
    status::String
    outer_folds::Vector{SyntheticOuterDecision}
    date_scores::Vector{Float64}
    sign_flip_p::Float64
    pair_prior::NTuple{4,Float64}
    sha256::String
end

struct NullCampaignRecord
    seed::Int
    nodes::Int
    edges::Int
    latent_ones::Int
    statuses::Vector{String}
    diagnostic_trigger::Bool
    meta_status::String
    meta_false_pass::Bool
    meta_sha256::String
    data_sha256::String
    edge_sha256::String
    result_sha256::String
end

struct NullCampaignSummary
    records::Vector{NullCampaignRecord}
    failures::Int
    campaigns::Int
    clopper_pearson_upper::Float64
    accepted::Bool
    sha256::String
end

struct FeatureDataset
    records::Vector{HUA.LobeRecord}
    date_by_scan::Dict{String,String}
    snapshot::FileSnapshot
end

struct FoldPlan
    dates::Vector{String}
    pairs::Vector{Tuple{String,String}}
    snapshot::FileSnapshot
end

struct ReceiptOuter
    outer_heldout_date::String
    edges::Int
    scale_eigenvalues::Vector{Float64}
    categories::Vector{EdgeCategorySupport}
end

struct InputReceipt
    by_outer::Dict{String,ReceiptOuter}
    snapshot::FileSnapshot
    unary_review::FileSnapshot
end

struct DiagnosticRun
    status::String
    reason::String
    aggregate_statuses::Vector{String}
    partitions::Vector{PartitionResult}
    null_campaigns::NullCampaignSummary
    config_hash::String
    feature_hash::String
    fold_hash::String
    input_receipt_hash::String
    unary_review_hash::String
    source_hash::String
    run_sha256::String
    snapshots::Vector{FileSnapshot}
end

_fail(code::Symbol, message::AbstractString) =
    throw(DiagnosticError(code, String(message)))
sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))

function _reject_symlink_components(path::String)
    cursor = first(splitpath(abspath(path)))
    islink(cursor) && _fail(:symlink_rejected, "path contains a symbolic link")
    for component in splitpath(abspath(path))[2:end]
        cursor = joinpath(cursor, component)
        islink(cursor) && _fail(:symlink_rejected, "path contains a symbolic link")
    end
    return nothing
end

function _snapshot(path::AbstractString; expected_path::Union{Nothing,String}=nothing,
                   expected_hash::Union{Nothing,String}=nothing)
    supplied = String(path)
    isempty(supplied) && _fail(:invalid_path, "path is empty")
    occursin('\0', supplied) && _fail(:invalid_path, "path contains NUL")
    any(==(".."), split(replace(supplied, '\\' => '/'), '/'; keepempty=true)) &&
        _fail(:path_escape, "parent traversal is forbidden")
    absolute = normpath(isabspath(supplied) ? supplied : joinpath(ROOT, supplied))
    _reject_symlink_components(absolute)
    isfile(absolute) || _fail(:missing_input, "required input is absent")
    expected_path === nothing || absolute == expected_path ||
        _fail(:path_substitution, "immutable input path differs")
    bytes = read(absolute)
    digest = sha256_hex(bytes)
    expected_hash === nothing || digest == expected_hash ||
        _fail(:hash_mismatch, "immutable input bytes differ")
    return FileSnapshot(absolute, bytes, digest)
end

function _verify(snapshot::FileSnapshot)
    current = _snapshot(snapshot.path)
    current.sha256 == snapshot.sha256 && current.bytes == snapshot.bytes ||
        _fail(:stale_input, "input changed after validation")
    return nothing
end

function _table(document::AbstractDict, name::String, context::String)
    value = get(document, name, nothing)
    value isa AbstractDict || _fail(:config_contract, "$context.$name is absent")
    return value
end

function _expect(table::AbstractDict, name::String, expected, context::String)
    get(table, name, nothing) == expected ||
        _fail(:config_contract, "$context.$name differs from the frozen declaration")
    return nothing
end

function _finite_float(value, context::String)
    value isa Real && !(value isa Bool) || _fail(:config_contract, "$context is not numeric")
    result = Float64(value)
    isfinite(result) || _fail(:config_contract, "$context is not finite")
    return result
end

function _positive_integer(value, context::String)
    value isa Integer && !(value isa Bool) ||
        _fail(:config_contract, "$context is not an integer")
    result = Int(value)
    result > 0 || _fail(:config_contract, "$context must be positive")
    return result
end

function load_diagnostic_config(path::AbstractString=MODEL_CONFIG_PATH)
    snapshot = _snapshot(path; expected_path=MODEL_CONFIG_PATH,
                         expected_hash=MODEL_CONFIG_SHA256)
    document = try
        TOML.parse(String(copy(snapshot.bytes)))
    catch error
        _fail(:config_parse, sprint(showerror, error))
    end
    model = _table(document, "model", "config")
    views = _table(model, "views", "model")
    edge = _table(model, "edge", "model")
    graph = _table(model, "graph", "model")
    selection = _table(document, "selection", "config")
    bootstrap = _table(selection, "bootstrap", "selection")
    scoring = _table(selection, "scoring", "selection")
    unary_gate = _table(selection, "unary_gate", "selection")
    edge_gate = _table(selection, "edge_gate", "selection")
    diagnostics = _table(selection, "diagnostics", "selection")
    policy_table = _table(diagnostics, "policy", "selection.diagnostics")
    primary = _table(selection, "primary_test", "selection")
    feasibility = _table(selection, "feasibility", "selection")
    preprocessing = _table(document, "preprocessing", "config")

    _expect(model, "schema", "structured_label_free_unit_assignment_model_v2", "model")
    _expect(model, "version", 2, "model")
    _expect(model, "model_ids", ["C1", "C2", "edge_admission", "admitted_graph"],
            "model")
    _expect(model, "node_priors", [0.5, 0.5], "model")
    _expect(model, "mixture_weights", [0.5, 0.5], "model")
    _expect(model, "student_t_nu", SRE.FIXED_NU, "model")
    _expect(model, "covariance_floor", SRE.SCALE_FLOOR, "model")
    _expect(model, "n_starts", length(SRE.START_QUANTILE_DELTAS), "model")
    _expect(model, "unary_start_quantile_deltas", collect(SRE.START_QUANTILE_DELTAS),
            "model")
    _expect(model, "max_iter", SRE.MAX_ITER, "model")
    _expect(model, "tol", SRE.CONVERGENCE_TOL, "model")
    _expect(model, "objective_decrease_tolerance", SRE.DECREASE_TOL, "model")
    _expect(views, "fusion", "equal_logit_mean", "model.views")
    _expect(views, "fixed_equal_weights", [0.25, 0.25, 0.25, 0.25], "model.views")
    _expect(edge, "basis", ["corr_fwd", "corr_bwd"], "model.edge")
    _expect(edge, "dimension", 2, "model.edge")
    _expect(edge, "category_order", CATEGORY_NAMES, "model.edge")
    _expect(edge, "student_t_nu", SRE.FIXED_NU, "model.edge")
    _expect(edge, "mixed_state_mean", "tied_01_10", "model.edge")
    _expect(edge, "scale", "shared_full_2x2", "model.edge")
    _expect(edge, "pair_prior", collect(SYNTHETIC_PAIR_PRIOR), "model.edge")
    _expect(edge, "pair_prior_policy", "zero_data_independent", "model.edge")
    _expect(graph, "inference", "exact_log_domain_two_state_chain", "model.graph")
    _expect(graph, "enable_condition", "verified_edge_admission_pass", "model.graph")
    _expect(graph, "reversal_equivariance_required", true, "model.graph")
    _expect(selection, "meta_algorithm", "nested_C1_C2_edge_admission_admitted_graph",
            "selection")
    _expect(selection, "common_reference", "C1_with_state_independent_edge_null",
            "selection")
    _expect(unary_gate, "comparison", "C2_minus_C1", "selection.unary_gate")
    _expect(unary_gate, "nested_inner_leave_date_out", true, "selection.unary_gate")
    _expect(unary_gate, "every_inner_date_mean_positive", true,
            "selection.unary_gate")
    _expect(unary_gate, "bootstrap_lower_bound_strictly_positive", true,
            "selection.unary_gate")
    _expect(edge_gate, "comparison",
            "conditional_edge_density_minus_state_independent_edge_null",
            "selection.edge_gate")
    _expect(edge_gate, "nested_inner_leave_date_out", true, "selection.edge_gate")
    _expect(edge_gate, "every_inner_date_mean_positive", true,
            "selection.edge_gate")
    _expect(edge_gate, "bootstrap_lower_bound_strictly_positive", true,
            "selection.edge_gate")
    _expect(edge_gate, "observed_above_shuffle_quantile", true,
            "selection.edge_gate")
    _expect(edge_gate, "reversal_equivariance_required", true,
            "selection.edge_gate")
    _expect(scoring, "numerical_failure", "FAIL", "selection.scoring")
    _expect(diagnostics, "follow_up_only", true, "selection.diagnostics")
    _expect(diagnostics, "mechanisms", DIAGNOSTIC_NAMES, "selection.diagnostics")
    _expect(diagnostics, "allowed_outcomes", TERMINAL_STATES, "selection.diagnostics")
    _expect(diagnostics, "outer_training_partitions_only", true,
            "selection.diagnostics")
    _expect(diagnostics, "multiple_testing", "Holm", "selection.diagnostics")
    _expect(policy_table, "schema", "structured_diagnostics_policy_v1",
            "selection.diagnostics.policy")
    _expect(policy_table, "residual_model", "fixed_nu8_diagonal_student_t_two",
            "selection.diagnostics.policy")
    _expect(policy_table, "residual_formula",
            "x_minus_posterior_weighted_component_mean", "selection.diagnostics.policy")
    _expect(policy_table, "residual_fallback", "matched_one_component_student_t_mean",
            "selection.diagnostics.policy")
    _expect(policy_table, "residual_failure", "BLOCKED", "selection.diagnostics.policy")
    _expect(policy_table, "icc_estimator", "ICC_1_1_oneway_random_unbalanced",
            "selection.diagnostics.policy")
    _expect(policy_table, "icc_date_centering", "row_weighted_per_date_per_feature",
            "selection.diagnostics.policy")
    _expect(policy_table, "icc_feature_pooling", "equal_feature_summed_squares",
            "selection.diagnostics.policy")
    _expect(policy_table, "icc_unbalanced_group_size",
            "n0=(N-sum(n_s^2)/N)/(K-1)", "selection.diagnostics.policy")
    _expect(policy_table, "view_score", "mean_lobes(log_predictive_density/dimension)",
            "selection.diagnostics.policy")
    _expect(policy_table, "view_contrast",
            "mean(bwd_com,bwd_diag)-mean(base,split)", "selection.diagnostics.policy")
    _expect(policy_table, "channel_drop_features",
            ["bwd_neg_com_t", "bwd_neg_diag45"], "selection.diagnostics.policy")
    _expect(policy_table, "channel_drop_refit", "fresh_inner_training_base_local",
            "selection.diagnostics.policy")
    _expect(policy_table, "channel_drop_statistic",
            "mean(full_bwd_scores)-refit_base_score", "selection.diagnostics.policy")
    _expect(policy_table, "channel_drop_same_sign", "strict_nonzero_each_inner_date",
            "selection.diagnostics.policy")
    _expect(policy_table, "quantile_method", "Hyndman_Fan_type_7_explicit",
            "selection.diagnostics.policy")
    _expect(policy_table, "permutation_tail", "upper_inclusive_plus_one",
            "selection.diagnostics.policy")
    _expect(policy_table, "view_tail", "two_sided_zero_inclusive_plus_one",
            "selection.diagnostics.policy")
    _expect(policy_table, "threshold_equality", "SKIPPED",
            "selection.diagnostics.policy")
    _expect(policy_table, "holm_order", DIAGNOSTIC_NAMES,
            "selection.diagnostics.policy")
    _expect(policy_table, "holm_reject_equality", true,
            "selection.diagnostics.policy")
    _expect(preprocessing, "normalization", "per_scan_median_mad", "preprocessing")
    _expect(preprocessing, "normalization_scope", "each_scan_uses_only_its_own_lobes",
            "preprocessing")

    view_names = String.(views["names"])
    base = String.(views["base_local"])
    backward = String.(views["backward_descriptors"])
    split = String.(views["split_descriptor"])
    view_features = [
        copy(base), vcat(base, backward[1:1]), vcat(base, backward[2:2]),
        vcat(base, split),
    ]
    length(view_names) == length(view_features) == 4 ||
        _fail(:config_contract, "exactly four declared views are required")
    feature_names = sort!(unique(reduce(vcat, view_features)))

    bootstrap_count = _positive_integer(bootstrap["count"], "selection.bootstrap.count")
    first_seed = Int(bootstrap["seed_start"])
    last_seed = Int(bootstrap["seed_stop"])
    seeds = collect(first_seed:last_seed)
    length(seeds) == bootstrap_count ||
        _fail(:config_contract, "bootstrap seed range is incomplete")
    _expect(bootstrap, "seed_family", "integer_range_0_499", "selection.bootstrap")
    _expect(bootstrap, "resample_unit", "whole_scan_within_date", "selection.bootstrap")

    permutation_count = _positive_integer(
        diagnostics["permutation_count"], "selection.diagnostics.permutation_count")
    permutation_count <= length(seeds) ||
        _fail(:config_contract, "permutation seeds exceed frozen seed family")
    null_count = _positive_integer(
        diagnostics["null_campaign_count"], "selection.diagnostics.null_campaign_count")
    null_count == length(seeds) ||
        _fail(:config_contract, "null campaign count must equal the frozen seed family")
    interval_level = _finite_float(
        diagnostics["view_interval_level"], "selection.diagnostics.view_interval_level")
    0.0 < interval_level < 1.0 ||
        _fail(:config_contract, "view interval level is outside (0,1)")
    permutation_quantile = _finite_float(
        diagnostics["permutation_upper_quantile"],
        "selection.diagnostics.permutation_upper_quantile")
    bootstrap_lower = _finite_float(
        _table(selection, "unary_gate", "selection")["bootstrap_lower_quantile"],
        "selection.unary_gate.bootstrap_lower_quantile")
    bootstrap_upper = _finite_float(
        _table(selection, "noninferiority", "selection")["bootstrap_upper_quantile"],
        "selection.noninferiority.bootstrap_upper_quantile")
    alpha = _finite_float(primary["alpha"], "selection.primary_test.alpha")

    contract = DiagnosticConfig(
        snapshot.path, snapshot.sha256, view_names, view_features, feature_names,
        permutation_count, permutation_quantile, seeds, bootstrap_lower,
        bootstrap_upper, interval_level, alpha, seeds,
        _positive_integer(diagnostics["null_campaign_min_dates"],
                          "selection.diagnostics.null_campaign_min_dates"),
        _positive_integer(diagnostics["null_campaign_min_scans"],
                          "selection.diagnostics.null_campaign_min_scans"),
        _positive_integer(diagnostics["null_campaign_chain_length_min"],
                          "selection.diagnostics.null_campaign_chain_length_min"),
        _positive_integer(diagnostics["null_campaign_chain_length_max"],
                          "selection.diagnostics.null_campaign_chain_length_max"),
        _finite_float(diagnostics["clopper_pearson_confidence"],
                      "selection.diagnostics.clopper_pearson_confidence"),
        _finite_float(diagnostics["clopper_pearson_upper_max"],
                      "selection.diagnostics.clopper_pearson_upper_max"),
        _positive_integer(feasibility["minimum_dates"],
                          "selection.feasibility.minimum_dates"),
        _positive_integer(feasibility["minimum_scans"],
                          "selection.feasibility.minimum_scans"),
        _finite_float(feasibility["minimum_category_kish_ess"],
                      "selection.feasibility.minimum_category_kish_ess"),
        _positive_integer(feasibility["minimum_distinct_edges"],
                          "selection.feasibility.minimum_distinct_edges"),
        _positive_integer(feasibility["minimum_support_scans_per_training_date"],
                          "selection.feasibility.minimum_support_scans_per_training_date"),
        _finite_float(feasibility["effective_observations_per_free_parameter"],
                      "selection.feasibility.effective_observations_per_free_parameter"),
        _finite_float(feasibility["edge_scale_eigenvalue_floor"],
                      "selection.feasibility.edge_scale_eigenvalue_floor"),
        _finite_float(feasibility["edge_scale_condition_number_cap"],
                      "selection.feasibility.edge_scale_condition_number_cap"),
        DiagnosticPolicy(
            String(policy_table["schema"]),
            String(policy_table["residual_model"]),
            String(policy_table["residual_formula"]),
            String(policy_table["residual_fallback"]),
            String(policy_table["residual_failure"]),
            String(policy_table["icc_estimator"]),
            String(policy_table["icc_date_centering"]),
            String(policy_table["icc_feature_pooling"]),
            String(policy_table["icc_unbalanced_group_size"]),
            String(policy_table["view_score"]),
            String(policy_table["view_contrast"]),
            String.(policy_table["channel_drop_features"]),
            String(policy_table["channel_drop_refit"]),
            String(policy_table["channel_drop_statistic"]),
            String(policy_table["channel_drop_same_sign"]),
            String(policy_table["quantile_method"]),
            String(policy_table["permutation_tail"]),
            String(policy_table["view_tail"]),
            String(policy_table["threshold_equality"]),
            String.(policy_table["holm_order"]),
            policy_table["holm_reject_equality"] === true,
        ),
        snapshot,
    )
    contract.null_chain_min <= contract.null_chain_max ||
        _fail(:config_contract, "null chain-length range is reversed")
    _verify(snapshot)
    return contract
end

function first_component_statistic(values::AbstractMatrix)
    matrix = Matrix{Float64}(values)
    n, p = size(matrix)
    n >= 2 || _fail(:underpowered, "factor statistic needs at least two observations")
    p >= 1 || _fail(:malformed_input, "factor statistic needs a positive dimension")
    all(isfinite, matrix) || _fail(:nonfinite_input, "factor values are nonfinite")
    centered = matrix .- mean(matrix; dims=1)
    covariance = Symmetric((centered' * centered) / (n - 1))
    eigenvalues = eigvals(covariance)
    all(isfinite, eigenvalues) || _fail(:numerical_failure, "factor eigenvalues are nonfinite")
    return maximum(eigenvalues)
end

function quantile_type7(values::AbstractVector, probability::Real)
    data = sort!(Float64.(collect(values)))
    isempty(data) && _fail(:underpowered, "Type-7 quantile needs at least one value")
    all(isfinite, data) || _fail(:nonfinite_input, "Type-7 quantile values are nonfinite")
    q = Float64(probability)
    isfinite(q) && 0.0 <= q <= 1.0 ||
        _fail(:malformed_input, "Type-7 quantile probability is outside [0,1]")
    q == 0.0 && return first(data)
    q == 1.0 && return last(data)
    h = 1.0 + (length(data) - 1) * q
    lower_index = floor(Int, h)
    gamma = h - lower_index
    return (1.0 - gamma) * data[lower_index] + gamma * data[lower_index + 1]
end

function _strict_diagnostic_positive(value::Real)
    return isfinite(value) && Float64(value) > 0.0
end

function _view_contrast(scores::AbstractVector{<:Real}, config::DiagnosticConfig)
    length(scores) == length(config.view_names) ||
        _fail(:config_contract, "view score count differs from configured views")
    backward = findall(name -> occursin("bwd", lowercase(name)), config.view_names)
    non_backward = findall(name -> !occursin("bwd", lowercase(name)), config.view_names)
    length(backward) == length(non_backward) == 2 ||
        _fail(:config_contract, "frozen channel view groups differ")
    values = Float64.(scores)
    all(isfinite, values) || _fail(:nonfinite_input, "view scores are nonfinite")
    return mean(values[index] for index in backward) -
           mean(values[index] for index in non_backward)
end

function _concatenate(scans::Vector{ResidualScan})
    isempty(scans) && _fail(:underpowered, "no residual scans were supplied")
    dimension = size(first(scans).values, 2)
    dimension > 0 || _fail(:malformed_input, "residual dimension is zero")
    all(size(scan.values, 2) == dimension for scan in scans) ||
        _fail(:malformed_input, "residual dimensions differ")
    values = reduce(vcat, (scan.values for scan in scans))
    labels = reduce(vcat, (fill(scan.scan, size(scan.values, 1)) for scan in scans))
    dates = reduce(vcat, (fill(scan.date, size(scan.values, 1)) for scan in scans))
    return values, labels, dates
end

function scan_icc(values::AbstractMatrix, scans::AbstractVector{<:AbstractString},
                  dates::AbstractVector{<:AbstractString})
    matrix = Matrix{Float64}(values)
    n, p = size(matrix)
    length(scans) == length(dates) == n ||
        _fail(:malformed_input, "ICC labels do not match observations")
    n >= 3 && p >= 1 || _fail(:underpowered, "ICC has insufficient observations")
    all(isfinite, matrix) || _fail(:nonfinite_input, "ICC values are nonfinite")

    centered = copy(matrix)
    for date in sort!(unique(String.(dates)))
        indices = findall(==(date), String.(dates))
        isempty(indices) && continue
        centered[indices, :] .-= mean(centered[indices, :]; dims=1)
    end
    groups = sort!(unique(String.(scans)))
    k = length(groups)
    k >= 2 && n > k || _fail(:underpowered, "ICC needs at least two non-singleton scans")
    grand = mean(centered; dims=1)
    ss_between = 0.0
    ss_within = 0.0
    sizes = Int[]
    for scan in groups
        indices = findall(==(scan), String.(scans))
        push!(sizes, length(indices))
        group_mean = mean(centered[indices, :]; dims=1)
        ss_between += length(indices) * sum(abs2, group_mean .- grand)
        ss_within += sum(abs2, centered[indices, :] .- group_mean)
    end
    mean_between = ss_between / (k - 1)
    mean_within = ss_within / (n - k)
    effective_group_size = (n - sum(size^2 for size in sizes) / n) / (k - 1)
    denominator = mean_between + (effective_group_size - 1.0) * mean_within
    denominator > 0.0 || return 0.0
    value = (mean_between - mean_within) / denominator
    return clamp(value, -1.0, 1.0)
end

function holm_adjust(names::Vector{String}, raw_p::Vector{Float64}, alpha::Real)
    length(names) == length(raw_p) || _fail(:malformed_input, "Holm inputs differ")
    length(names) == length(unique(names)) || _fail(:malformed_input, "Holm names repeat")
    all(value -> isfinite(value) && 0.0 <= value <= 1.0, raw_p) ||
        _fail(:malformed_input, "Holm p-values are invalid")
    level = Float64(alpha)
    0.0 < level < 1.0 || _fail(:malformed_input, "Holm alpha is invalid")
    count = length(raw_p)
    order = sortperm(1:count; by=index -> (raw_p[index], index))
    adjusted = zeros(Float64, count)
    running = 0.0
    for (rank, index) in enumerate(order)
        running = max(running, (count - rank + 1) * raw_p[index])
        adjusted[index] = min(1.0, running)
    end
    rejected = falses(count)
    open_sequence = true
    for (rank, index) in enumerate(order)
        local_reject = raw_p[index] <= level / (count - rank + 1)
        rejected[index] = open_sequence && local_reject
        open_sequence &= local_reject
    end
    return HolmResult(copy(names), copy(raw_p), adjusted, rejected, order)
end

function kish_ess(sum_weights::Real, sum_squared_weights::Real)
    first = Float64(sum_weights)
    second = Float64(sum_squared_weights)
    isfinite(first) && isfinite(second) && first >= 0.0 && second >= 0.0 ||
        _fail(:malformed_input, "Kish inputs are invalid")
    second > 0.0 || return 0.0
    result = first * first / second
    isfinite(result) || _fail(:numerical_failure, "Kish ESS is nonfinite")
    return result
end

diagonal_two_component_parameter_count(dimension::Integer) = begin
    dimension > 0 || _fail(:malformed_input, "model dimension must be positive")
    4 * Int(dimension)
end

function edge_admission_parameter_count(dimension::Integer=2)
    dimension == 2 || _fail(:config_contract, "edge dimension must remain two")
    means = 3 * Int(dimension) # 00, tied 01/10, 11
    shared_symmetric_scale = Int(dimension) * (Int(dimension) + 1) ÷ 2
    return means + shared_symmetric_scale
end

function scale_is_feasible(eigenvalues::AbstractVector, floor::Real, cap::Real)
    values = Float64.(eigenvalues)
    length(values) == 2 || return false
    all(isfinite, values) && all(>(0.0), values) || return false
    minimum(values) >= Float64(floor) || return false
    return maximum(values) / minimum(values) <= Float64(cap)
end

function _feasibility_reason(input::FeasibilityInput, config::DiagnosticConfig)
    reasons = String[]
    length(input.dates) >= config.minimum_dates || push!(reasons, "insufficient_dates")
    input.scans >= config.minimum_scans || push!(reasons, "insufficient_scans")
    input.edges >= config.minimum_edges || push!(reasons, "insufficient_distinct_edges")
    names = sort([category.name for category in input.categories])
    names == sort(CATEGORY_NAMES) || push!(reasons, "incomplete_categories")
    for category in input.categories
        kish_ess(category.sum_weights, category.sum_squared_weights) >=
            config.minimum_category_ess || push!(reasons, "insufficient_category_ess")
        for date in input.dates
            get(category.support_scans_by_date, date, 0) >= config.minimum_support_scans ||
                push!(reasons, "insufficient_category_scan_support")
        end
    end
    edge_parameters = edge_admission_parameter_count()
    input.edges >= ceil(Int, config.observations_per_parameter * edge_parameters) ||
        push!(reasons, "insufficient_effective_observations")
    scale_is_feasible(input.scale_eigenvalues, config.scale_eigenvalue_floor,
                      config.scale_condition_cap) || push!(reasons, "unsupported_edge_scale")
    return unique(reasons)
end

function feasibility_table(input::FeasibilityInput, config::DiagnosticConfig)
    length(input.dates) == length(unique(input.dates)) ||
        _fail(:malformed_input, "feasibility dates repeat")
    input.scans >= 0 && input.nodes >= 0 && input.edges >= 0 ||
        _fail(:malformed_input, "feasibility counts are negative")
    view_dimensions = length.(config.view_features)
    unary_parameters = sum(diagonal_two_component_parameter_count, view_dimensions)
    edge_parameters = edge_admission_parameter_count()
    category_values = [kish_ess(item.sum_weights, item.sum_squared_weights)
                       for item in input.categories]
    minimum_ess = isempty(category_values) ? 0.0 : minimum(category_values)
    support_values = [get(item.support_scans_by_date, date, 0)
                      for item in input.categories for date in input.dates]
    minimum_support = isempty(support_values) ? 0 : minimum(support_values)
    minimum_eigenvalue = length(input.scale_eigenvalues) == 2 &&
                         all(isfinite, input.scale_eigenvalues) ?
                         minimum(input.scale_eigenvalues) : NaN
    condition_number = length(input.scale_eigenvalues) == 2 &&
                       all(value -> isfinite(value) && value > 0.0,
                           input.scale_eigenvalues) ?
                       maximum(input.scale_eigenvalues) / minimum(input.scale_eigenvalues) : Inf
    edge_reasons = _feasibility_reason(input, config)
    edge_status = isempty(edge_reasons) ? "SUPPORTED" : "SKIPPED"
    edge_reason = isempty(edge_reasons) ? "ok" : join(edge_reasons, ",")

    ratio(parameters, observations) = observations > 0 ? parameters / observations : Inf
    common = (length(input.dates), input.scans, input.nodes, input.edges,
              minimum_ess, minimum_support, minimum_eigenvalue, condition_number)
    return FeasibilityRow[
        FeasibilityRow("C1", view_dimensions, unary_parameters, 0, common[1], common[2],
                       common[3], common[4], input.nodes,
                       ratio(unary_parameters, input.nodes), CATEGORY_NAMES,
                       category_values, support_values, common[5], common[6],
                       common[7], common[8], input.nodes > 0 ? "SUPPORTED" : "SKIPPED",
                       input.nodes > 0 ? "descriptive_unary_support" : "no_nodes"),
        FeasibilityRow("C2", view_dimensions, unary_parameters, 0, common[1], common[2],
                       common[3], common[4], input.nodes,
                       ratio(unary_parameters, input.nodes), CATEGORY_NAMES,
                       category_values, support_values, common[5], common[6],
                       common[7], common[8], input.nodes > 0 ? "SUPPORTED" : "SKIPPED",
                       input.nodes > 0 ? "descriptive_unary_support" : "no_nodes"),
        FeasibilityRow("edge_admission", [2], edge_parameters, 0, common[1], common[2],
                       common[3], common[4], input.edges,
                       ratio(edge_parameters, input.edges), CATEGORY_NAMES,
                       category_values, support_values, common[5], common[6],
                       common[7], common[8], edge_status, edge_reason),
        FeasibilityRow("admitted_graph", [2], 0, edge_parameters, common[1], common[2],
                       common[3], common[4], input.nodes + input.edges,
                       ratio(edge_parameters, input.nodes + input.edges), CATEGORY_NAMES,
                       category_values, support_values, common[5], common[6], common[7],
                       common[8], edge_status,
                       isempty(edge_reasons) ? "inherits_verified_edge_support" : edge_reason),
    ]
end

function _binomial_cdf(failures::Int, campaigns::Int, probability::Float64)
    failures < 0 && return 0.0
    failures >= campaigns && return 1.0
    probability <= 0.0 && return 1.0
    probability >= 1.0 && return 0.0
    term = (1.0 - probability)^campaigns
    total = term
    for count in 0:(failures - 1)
        term *= (campaigns - count) / (count + 1) * probability / (1.0 - probability)
        total += term
    end
    return clamp(total, 0.0, 1.0)
end

function clopper_pearson_upper(failures::Integer, campaigns::Integer,
                               confidence::Real)
    n = Int(campaigns)
    x = Int(failures)
    n > 0 && 0 <= x <= n || _fail(:malformed_input, "binomial counts are invalid")
    level = Float64(confidence)
    0.0 < level < 1.0 || _fail(:malformed_input, "confidence is invalid")
    x == n && return 1.0
    target = 1.0 - level
    lower = 0.0
    upper = 1.0
    for _ in 1:100
        midpoint = (lower + upper) / 2.0
        if _binomial_cdf(x, n, midpoint) > target
            lower = midpoint
        else
            upper = midpoint
        end
    end
    return (lower + upper) / 2.0
end

function _seeded_rng(seed::Int, parts::AbstractString...)
    material = "structured-diagnostics-seed-v1\nseed=$seed\n" *
               join(("part=$(String(part))\n" for part in parts))
    digest = sha256(collect(codeunits(material)))
    value = zero(UInt64)
    for byte in digest[1:8]
        value = (value << 8) | UInt64(byte)
    end
    return MersenneTwister(UInt32(value % UInt64(typemax(UInt32))))
end

function _permuted_factor(scans::Vector{ResidualScan}, seed::Int, context::String)
    rng = _seeded_rng(seed, "factor", context)
    permuted = Matrix{Float64}[]
    for scan in scans
        matrix = copy(scan.values)
        for column in axes(matrix, 2)
            order = randperm(rng, size(matrix, 1))
            matrix[:, column] = matrix[order, column]
        end
        push!(permuted, matrix)
    end
    return first_component_statistic(reduce(vcat, permuted))
end

function _resample_scans(scans::Vector{ResidualScan}, seed::Int, context::String)
    rng = _seeded_rng(seed, "bootstrap", context)
    output = ResidualScan[]
    for date in sort!(unique(scan.date for scan in scans))
        candidates = [scan for scan in scans if scan.date == date]
        isempty(candidates) && _fail(:incomplete_fold, "bootstrap date has no scans")
        for draw in eachindex(candidates)
            selected = candidates[rand(rng, eachindex(candidates))]
            push!(output, ResidualScan(date, "$(selected.scan)#draw$(draw)",
                                       copy(selected.values)))
        end
    end
    return output
end

function _permuted_scan_labels(scans::Vector{ResidualScan}, seed::Int, context::String)
    values, labels, dates = _concatenate(scans)
    rng = _seeded_rng(seed, "scan-label", context)
    permuted = copy(labels)
    for date in sort!(unique(dates))
        indices = findall(==(date), dates)
        permuted[indices] = labels[indices][randperm(rng, length(indices))]
    end
    return scan_icc(values, permuted, dates)
end

function _date_balanced_mean(values::Dict{String,Float64}, scans::Vector{ResidualScan})
    date_means = Float64[]
    for date in sort!(unique(scan.date for scan in scans))
        names = sort!(unique(scan.scan for scan in scans if scan.date == date))
        isempty(names) && _fail(:incomplete_fold, "view date has no scans")
        all(haskey(values, name) for name in names) ||
            _fail(:incomplete_fold, "view evidence omits a scan")
        push!(date_means, mean(values[name] for name in names))
    end
    return mean(date_means)
end

function _bootstrap_view(values::Dict{String,Float64}, scans::Vector{ResidualScan},
                         seed::Int, context::String)
    rng = _seeded_rng(seed, "view-bootstrap", context)
    date_means = Float64[]
    for date in sort!(unique(scan.date for scan in scans))
        names = sort!(unique(scan.scan for scan in scans if scan.date == date))
        sampled = [names[rand(rng, eachindex(names))] for _ in eachindex(names)]
        push!(date_means, mean(values[name] for name in sampled))
    end
    return mean(date_means)
end

function _skipped_mechanism(name::String, reason::String)
    return MechanismResult(name, "SKIPPED", reason, NaN, NaN, NaN, NaN, 1.0,
                           1.0, DateFactorResult[], Float64[], Float64[])
end

function _partition_hash(outer::String, mechanisms::Vector{MechanismResult},
                         feasibility::Vector{FeasibilityRow})
    io = IOBuffer()
    println(io, "structured-diagnostic-partition-v1")
    println(io, "outer=", outer)
    for result in mechanisms
        println(io, join((result.name, result.status, result.reason,
                          bitstring(result.observed), bitstring(result.threshold),
                          bitstring(result.lower), bitstring(result.upper),
                          bitstring(result.raw_p), bitstring(result.holm_p)), '\t'))
        for date in result.date_results
            println(io, join((date.date, bitstring(date.observed),
                              bitstring(date.threshold), bitstring(date.excess),
                              bitstring(date.raw_p)), '\t'))
        end
        for value in result.permutation_values
            println(io, "perm\t", bitstring(value))
        end
        for value in result.bootstrap_values
            println(io, "boot\t", bitstring(value))
        end
    end
    for row in feasibility
        println(io, join((row.model, join(row.dimensions, ','), row.free_parameters,
                          row.inherited_parameters, row.dates, row.scans, row.nodes,
                          row.edges, row.effective_observations,
                          bitstring(row.parameter_observation_ratio),
                          join(row.category_names, ','),
                          join(bitstring.(row.category_ess), ','),
                          join(row.category_support_scans, ','), row.status,
                          row.reason), '\t'))
    end
    return sha256_hex(take!(io))
end

function evaluate_partition(evidence::PartitionEvidence, support::FeasibilityInput,
                            config::DiagnosticConfig)
    scans = evidence.residual_scans
    isempty(scans) && _fail(:incomplete_fold, "partition has no residual scans")
    scan_names = [scan.scan for scan in scans]
    length(scan_names) == length(unique(scan_names)) ||
        _fail(:duplicate_scan, "partition repeats a scan")
    all(size(scan.values, 1) >= 2 && all(isfinite, scan.values) for scan in scans) ||
        _fail(:nonfinite_input, "residual scan is empty or nonfinite")
    dates = sort!(unique(scan.date for scan in scans))
    sort(support.dates) == dates || _fail(:incomplete_fold, "feasibility dates differ")
    support.scans == length(scans) || _fail(:incomplete_fold, "scan count differs")
    support.nodes == sum(size(scan.values, 1) for scan in scans) ||
        _fail(:incomplete_fold, "node count differs")
    Set(keys(evidence.view_score_difference)) == Set(scan_names) ||
        _fail(:incomplete_fold, "view scores differ from scan universe")
    Set(keys(evidence.dropout_difference)) == Set(scan_names) ||
        _fail(:incomplete_fold, "dropout scores differ from scan universe")
    all(isfinite, Base.values(evidence.view_score_difference)) ||
        _fail(:nonfinite_input, "view scores are nonfinite")
    all(isfinite, Base.values(evidence.dropout_difference)) ||
        _fail(:nonfinite_input, "dropout scores are nonfinite")
    table = feasibility_table(support, config)

    if length(dates) < config.minimum_dates || length(scans) < config.minimum_scans
        mechanisms = [_skipped_mechanism(name, "underpowered_partition")
                      for name in DIAGNOSTIC_NAMES]
        digest = _partition_hash(evidence.outer_heldout_date, mechanisms, table)
        return PartitionResult(evidence.outer_heldout_date, "SKIPPED", mechanisms,
                               table, digest)
    end

    date_results = DateFactorResult[]
    thresholds = Dict{String,Float64}()
    for date in dates
        date_scans = [scan for scan in scans if scan.date == date]
        observed = first_component_statistic(reduce(vcat, (scan.values for scan in date_scans)))
        permutations = [_permuted_factor(
            date_scans, seed, "$(evidence.outer_heldout_date):$date")
                        for seed in config.bootstrap_seeds[1:config.permutation_count]]
        threshold = quantile_type7(permutations, config.permutation_quantile)
        raw_p = (1 + count(value -> value >= observed, permutations)) /
                (length(permutations) + 1)
        thresholds[date] = threshold
        push!(date_results, DateFactorResult(
            date, observed, threshold, observed - threshold, raw_p, permutations))
    end
    factor_bootstrap = Float64[]
    for seed in config.bootstrap_seeds
        sampled = _resample_scans(scans, seed, "factor:$(evidence.outer_heldout_date)")
        excesses = Float64[]
        for date in dates
            matrices = [scan.values for scan in sampled if scan.date == date]
            push!(excesses, first_component_statistic(reduce(vcat, matrices)) -
                              thresholds[date])
        end
        push!(factor_bootstrap, minimum(excesses))
    end
    factor_lower = quantile_type7(factor_bootstrap, config.bootstrap_lower_quantile)
    factor_p = maximum(result.raw_p for result in date_results)

    residual_values, labels, row_dates = _concatenate(scans)
    observed_icc = scan_icc(residual_values, labels, row_dates)
    scan_permutations = [_permuted_scan_labels(
        scans, seed, "$(evidence.outer_heldout_date)")
                         for seed in config.bootstrap_seeds[1:config.permutation_count]]
    scan_threshold = quantile_type7(scan_permutations, config.permutation_quantile)
    scan_p = (1 + count(value -> value >= observed_icc, scan_permutations)) /
             (length(scan_permutations) + 1)
    scan_bootstrap = [begin
        sampled = _resample_scans(scans, seed, "scan:$(evidence.outer_heldout_date)")
        matrix, sampled_labels, sampled_dates = _concatenate(sampled)
        scan_icc(matrix, sampled_labels, sampled_dates) - scan_threshold
    end for seed in config.bootstrap_seeds]
    scan_lower = quantile_type7(scan_bootstrap, config.bootstrap_lower_quantile)

    observed_view = _date_balanced_mean(evidence.view_score_difference, scans)
    view_bootstrap = [_bootstrap_view(
        evidence.view_score_difference, scans, seed,
        "$(evidence.outer_heldout_date)") for seed in config.bootstrap_seeds]
    interval_tail = (1.0 - config.interval_level) / 2.0
    view_lower = quantile_type7(view_bootstrap, interval_tail)
    view_upper = quantile_type7(view_bootstrap, 1.0 - interval_tail)
    lower_tail = (1 + count(<=(0.0), view_bootstrap)) / (length(view_bootstrap) + 1)
    upper_tail = (1 + count(>=(0.0), view_bootstrap)) / (length(view_bootstrap) + 1)
    view_p = min(1.0, 2.0 * min(lower_tail, upper_tail))

    holm = holm_adjust(config.policy.holm_order, [factor_p, scan_p, view_p], config.alpha)
    factor_raw = all(result.observed > result.threshold for result in date_results) &&
                 factor_lower > 0.0
    scan_raw = observed_icc > scan_threshold && scan_lower > 0.0
    view_sign = sign(observed_view)
    date_view_signs = [sign(mean(evidence.view_score_difference[scan.scan]
                                 for scan in scans if scan.date == date)) for date in dates]
    date_dropout_signs = [sign(mean(evidence.dropout_difference[scan.scan]
                                    for scan in scans if scan.date == date)) for date in dates]
    interval_excludes_zero = view_lower > 0.0 || view_upper < 0.0
    view_raw = interval_excludes_zero && view_sign != 0.0 &&
               _strict_diagnostic_positive(abs(view_lower)) &&
               _strict_diagnostic_positive(abs(view_upper)) &&
               all(==(view_sign), date_view_signs) &&
               all(==(view_sign), date_dropout_signs)

    statuses = [factor_raw && holm.rejected[1], scan_raw && holm.rejected[2],
                view_raw && holm.rejected[3]]
    mechanisms = MechanismResult[
        MechanismResult("mfa_q1", statuses[1] ? "FOLLOWUP_REQUIRED" : "SKIPPED",
                        statuses[1] ? "holm_factor_trigger" : "factor_trigger_absent",
                        minimum(result.observed for result in date_results),
                        maximum(result.threshold for result in date_results),
                        factor_lower, NaN, factor_p, holm.adjusted_p[1], date_results,
                        reduce(vcat, (result.permutation_values for result in date_results)),
                        factor_bootstrap),
        MechanismResult("scan_effects", statuses[2] ? "FOLLOWUP_REQUIRED" : "SKIPPED",
                        statuses[2] ? "holm_scan_trigger" : "scan_trigger_absent",
                        observed_icc, scan_threshold, scan_lower, NaN, scan_p,
                        holm.adjusted_p[2], DateFactorResult[], scan_permutations,
                        scan_bootstrap),
        MechanismResult("view_asymmetry", statuses[3] ? "FOLLOWUP_REQUIRED" : "SKIPPED",
                        statuses[3] ? "holm_view_trigger" : "view_trigger_absent",
                        observed_view, 0.0, view_lower, view_upper, view_p,
                        holm.adjusted_p[3], DateFactorResult[], Float64[],
                        view_bootstrap),
    ]
    partition_status = any(statuses) ? "FOLLOWUP_REQUIRED" : "SKIPPED"
    digest = _partition_hash(evidence.outer_heldout_date, mechanisms, table)
    return PartitionResult(evidence.outer_heldout_date, partition_status, mechanisms,
                           table, digest)
end

function _diagonal_student_t8_row(rng::AbstractRNG, dimension::Int)
    dimension > 0 || _fail(:malformed_input, "t8 dimension must be positive")
    denominator = sqrt(sum(abs2, randn(rng, SRE.FIXED_NU)) / SRE.FIXED_NU)
    return randn(rng, dimension) ./ denominator
end

function _bivariate_student_t8_covariance(rng::AbstractRNG)
    first_normal = randn(rng)
    second_normal = 0.3 * first_normal + sqrt(1.0 - 0.3^2) * randn(rng)
    denominator = sqrt(sum(abs2, randn(rng, SRE.FIXED_NU)) / SRE.FIXED_NU)
    covariance_scale = sqrt((SRE.FIXED_NU - 2.0) / SRE.FIXED_NU)
    return covariance_scale * first_normal / denominator,
           covariance_scale * second_normal / denominator
end

function _normalize_scan(values::Matrix{Float64})
    output = similar(values)
    for column in axes(values, 2)
        median_value = median(view(values, :, column))
        mad = median(abs.(view(values, :, column) .- median_value))
        mad > 0.0 && isfinite(mad) || _fail(:numerical_failure, "null scan has zero MAD")
        output[:, column] = (values[:, column] .- median_value) ./ mad
    end
    return output
end

function _null_view_scores(values::Matrix{Float64}, config::DiagnosticConfig)
    scores = Float64[]
    index_by_name = Dict(name => index for (index, name) in enumerate(config.feature_names))
    for features in config.view_features
        indices = [index_by_name[name] for name in features]
        score = mean(begin
            row_score = 0.0
            for index in indices
                row_score += SRE.log_student_t_diag([values[row, index]], [0.0], [1.0])
            end
            row_score / length(indices)
        end for row in axes(values, 1))
        push!(scores, score)
    end
    return scores
end

function _null_view_difference(values::Matrix{Float64}, config::DiagnosticConfig)
    return _view_contrast(_null_view_scores(values, config), config)
end

function _null_dropout_difference(values::Matrix{Float64}, config::DiagnosticConfig)
    scores = _null_view_scores(values, config)
    backward = findall(name -> occursin("bwd", lowercase(name)), config.view_names)
    dropout_indices = _view_indices(config, config.view_features[1])
    dropout_values = values[:, dropout_indices]
    amplitudes = vec(values[:, 1])
    dropout_fit = SRE.fit_student_t_two_component(dropout_values, amplitudes)
    dropout_score = if dropout_fit.status == SRE.ROBUST_VALID
        mean(begin
            first = log(0.5) + SRE.log_student_t_diag(
                view(dropout_values, row, :), view(dropout_fit.means, 1, :),
                view(dropout_fit.scales, 1, :))
            second = log(0.5) + SRE.log_student_t_diag(
                view(dropout_values, row, :), view(dropout_fit.means, 2, :),
                view(dropout_fit.scales, 2, :))
            _logsumexp(first, second) / size(dropout_values, 2)
        end for row in axes(dropout_values, 1))
    elseif dropout_fit.one_component.status == SRE.ROBUST_VALID
        mean(SRE.log_student_t_diag(view(dropout_values, row, :),
                                    dropout_fit.one_component.mean,
                                    dropout_fit.one_component.scale) /
             size(dropout_values, 2) for row in axes(dropout_values, 1))
    else
        _fail(:numerical_failure, "synthetic dropout fit failed")
    end
    return mean(scores[index] for index in backward) - dropout_score
end

abstract type _SyntheticViewFit end

struct _SyntheticGaussianViewFit <: _SyntheticViewFit
    fit::HUA.TwoComponentFit
    high_component::Int
end

struct _SyntheticStudentViewFit <: _SyntheticViewFit
    fit::SRE.StudentTTwoFit
end

struct _SyntheticUnaryFit
    model::String
    views::Vector{_SyntheticViewFit}
    training_dates::Vector{String}
    training_sha256::String
    fit_sha256::String
end

function _synthetic_scans(campaign::SyntheticNullCampaign,
                          dates::AbstractVector{<:AbstractString})
    selected = Set(String.(dates))
    scans = [scan for scan in campaign.scans if scan.date in selected]
    sort!(scans; by=scan -> (scan.date, scan.scan))
    return scans
end

function _synthetic_training_hash(scans::Vector{SyntheticNullScan})
    io = IOBuffer()
    println(io, "structured-null-training-v1")
    for scan in scans
        println(io, scan.date, '\t', scan.scan, '\t', size(scan.values, 1))
        for value in scan.values
            println(io, "x\t", bitstring(value))
        end
        for value in scan.amplitudes
            println(io, "a\t", bitstring(value))
        end
        for value in scan.edge_values
            println(io, "e\t", bitstring(value))
        end
    end
    return sha256_hex(take!(io))
end

function _view_indices(config::DiagnosticConfig, features::Vector{String})
    by_name = Dict(name => index for (index, name) in enumerate(config.feature_names))
    all(haskey(by_name, feature) for feature in features) ||
        _fail(:config_contract, "synthetic view contains an unknown feature")
    return [by_name[feature] for feature in features]
end

function _synthetic_view_matrix(scans::Vector{SyntheticNullScan}, indices::Vector{Int})
    isempty(scans) && _fail(:underpowered, "synthetic unary training set is empty")
    return reduce(vcat, (scan.values[:, indices] for scan in scans))
end

function _synthetic_amplitudes(scans::Vector{SyntheticNullScan})
    isempty(scans) && _fail(:underpowered, "synthetic amplitude training set is empty")
    return reduce(vcat, (scan.amplitudes for scan in scans))
end

function _high_amplitude_component(responsibilities::Matrix{Float64},
                                   amplitudes::Vector{Float64})
    size(responsibilities, 1) == length(amplitudes) ||
        _fail(:incomplete_fold, "synthetic orientation rows differ")
    assignments = [responsibilities[row, 1] >= responsibilities[row, 2] ? 1 : 2
                   for row in axes(responsibilities, 1)]
    counts = [count(==(component), assignments) for component in 1:2]
    all(>(0), counts) || return 0
    means = [mean(amplitudes[index] for index in eachindex(amplitudes)
                  if assignments[index] == component) for component in 1:2]
    all(isfinite, means) && means[1] != means[2] || return 0
    return means[1] > means[2] ? 1 : 2
end

function _synthetic_unary_hash(model::String, training_sha256::String,
                               views::Vector{_SyntheticViewFit})
    io = IOBuffer()
    println(io, "structured-null-unary-fit-v1")
    println(io, "model=", model)
    println(io, "training=", training_sha256)
    for (index, view_fit) in enumerate(views)
        println(io, "view=", index)
        if view_fit isa _SyntheticGaussianViewFit
            fit = view_fit.fit
            println(io, "family=gaussian\thigh=", view_fit.high_component)
            for value in fit.means
                println(io, "mean\t", bitstring(value))
            end
            for value in fit.vars
                println(io, "scale\t", bitstring(value))
            end
        else
            fit = (view_fit::_SyntheticStudentViewFit).fit
            println(io, "family=student_t\thigh=", fit.high_amplitude_component)
            for value in fit.means
                println(io, "mean\t", bitstring(value))
            end
            for value in fit.scales
                println(io, "scale\t", bitstring(value))
            end
        end
    end
    return sha256_hex(take!(io))
end

function _fit_synthetic_unary(campaign::SyntheticNullCampaign,
                              training_dates::Vector{String}, model::String,
                              config::DiagnosticConfig, fit_counter::Base.RefValue{Int})
    model in ("C1", "C2") || _fail(:config_contract, "unknown synthetic unary model")
    dates = sort!(unique(copy(training_dates)))
    scans = _synthetic_scans(campaign, dates)
    training_sha256 = _synthetic_training_hash(scans)
    amplitudes = _synthetic_amplitudes(scans)
    views = _SyntheticViewFit[]
    for features in config.view_features
        fit_counter[] += 1
        indices = _view_indices(config, features)
        values = _synthetic_view_matrix(scans, indices)
        if model == "C1"
            fit = try
                HUA.fit_em_two_component(
                    values, 0; n_starts=5, cov_floor=SRE.SCALE_FLOOR,
                    max_iter=SRE.MAX_ITER, tol=SRE.CONVERGENCE_TOL)
            catch
                return nothing
            end
            fit.converged && fit.monotone && isfinite(fit.loglik) || return nothing
            responsibilities = HUA.responsibilities(fit, values)
            all(isfinite, responsibilities) || return nothing
            high_component = _high_amplitude_component(responsibilities, amplitudes)
            high_component in (1, 2) || return nothing
            push!(views, _SyntheticGaussianViewFit(fit, high_component))
        else
            fit = try
                SRE.fit_student_t_two_component(values, amplitudes)
            catch
                return nothing
            end
            fit.status == SRE.ROBUST_VALID || return nothing
            fit.high_amplitude_component in (1, 2) || return nothing
            push!(views, _SyntheticStudentViewFit(fit))
        end
    end
    fit_sha256 = _synthetic_unary_hash(model, training_sha256, views)
    return _SyntheticUnaryFit(model, views, dates, training_sha256, fit_sha256)
end

function _synthetic_state_logs(fit::_SyntheticUnaryFit, scan::SyntheticNullScan,
                               config::DiagnosticConfig)
    length(fit.views) == length(config.view_features) ||
        _fail(:incomplete_fold, "synthetic unary view count differs")
    logs = zeros(Float64, size(scan.values, 1), 2)
    view_weight = 1.0 / length(fit.views)
    for (view_fit, features) in zip(fit.views, config.view_features)
        indices = _view_indices(config, features)
        values = scan.values[:, indices]
        high_component = view_fit isa _SyntheticGaussianViewFit ?
            view_fit.high_component :
            (view_fit::_SyntheticStudentViewFit).fit.high_amplitude_component
        low_component = 3 - high_component
        for row in axes(values, 1)
            if view_fit isa _SyntheticGaussianViewFit
                gaussian = view_fit.fit
                logs[row, 1] += view_weight * HUA._log_gaussian_diag(
                    view(values, row, :), view(gaussian.means, low_component, :),
                    view(gaussian.vars, low_component, :))
                logs[row, 2] += view_weight * HUA._log_gaussian_diag(
                    view(values, row, :), view(gaussian.means, high_component, :),
                    view(gaussian.vars, high_component, :))
            else
                student = (view_fit::_SyntheticStudentViewFit).fit
                logs[row, 1] += view_weight * SRE.log_student_t_diag(
                    view(values, row, :), view(student.means, low_component, :),
                    view(student.scales, low_component, :))
                logs[row, 2] += view_weight * SRE.log_student_t_diag(
                    view(values, row, :), view(student.means, high_component, :),
                    view(student.scales, high_component, :))
            end
        end
    end
    all(isfinite, logs) || _fail(:numerical_failure, "synthetic unary scores are nonfinite")
    return logs
end

function _synthetic_posteriors(logs::Matrix{Float64})
    probabilities = similar(logs)
    prior = log(0.5)
    for row in axes(logs, 1)
        normalizer = _logsumexp(prior + logs[row, 1], prior + logs[row, 2])
        isfinite(normalizer) || _fail(:numerical_failure, "synthetic posterior is nonfinite")
        probabilities[row, 1] = exp(prior + logs[row, 1] - normalizer)
        probabilities[row, 2] = exp(prior + logs[row, 2] - normalizer)
    end
    return probabilities
end

function _synthetic_node_loglik(fit::_SyntheticUnaryFit, scan::SyntheticNullScan,
                                config::DiagnosticConfig)
    logs = _synthetic_state_logs(fit, scan, config)
    prior = log(0.5)
    return sum(_logsumexp(prior + logs[row, 1], prior + logs[row, 2])
               for row in axes(logs, 1))
end

_synthetic_observation_count(scan::SyntheticNullScan) =
    size(scan.values, 1) + size(scan.edge_values, 1)

function _whole_scan_bootstrap_lower(gains::Dict{String,Float64},
                                     scans::Vector{SyntheticNullScan},
                                     dates::Vector{String}, config::DiagnosticConfig,
                                     context::String)
    by_scan = Dict(scan.scan => scan for scan in scans)
    Set(keys(gains)) == Set(keys(by_scan)) ||
        _fail(:incomplete_fold, "synthetic bootstrap scan universe differs")
    values = Float64[]
    for seed in config.bootstrap_seeds
        rng = _seeded_rng(seed, "meta-bootstrap", context)
        date_means = Float64[]
        for date in dates
            names = sort!([scan.scan for scan in scans if scan.date == date])
            isempty(names) && _fail(:incomplete_fold, "synthetic bootstrap date is empty")
            sampled = [names[rand(rng, eachindex(names))] for _ in eachindex(names)]
            push!(date_means, mean(gains[name] for name in sampled))
        end
        push!(values, mean(date_means))
    end
    return quantile_type7(values, config.bootstrap_lower_quantile), values
end

function _select_synthetic_unary(campaign::SyntheticNullCampaign, outer::String,
                                 config::DiagnosticConfig,
                                 fit_counter::Base.RefValue{Int})
    inner_dates = [date for date in campaign.dates if date != outer]
    c1_fits = Dict{String,_SyntheticUnaryFit}()
    c2_fits = Dict{String,_SyntheticUnaryFit}()
    gains = Dict{String,Float64}()
    valid_c2 = true
    io = IOBuffer()
    println(io, "structured-null-inner-unary-v1")
    println(io, "outer=", outer)
    for inner in inner_dates
        training_dates = [date for date in inner_dates if date != inner]
        c1 = _fit_synthetic_unary(campaign, training_dates, "C1", config, fit_counter)
        c1 === nothing && return (
            status="FAIL", model="C1", fits=Dict{String,_SyntheticUnaryFit}(),
            hash=sha256_hex(collect(codeunits("c1-fit-failed:$outer:$inner\n"))),
            lower=-Inf)
        c1_fits[inner] = c1
        c2 = _fit_synthetic_unary(campaign, training_dates, "C2", config, fit_counter)
        if c2 === nothing
            valid_c2 = false
            println(io, inner, "\t", c1.fit_sha256, "\tC2_INVALID")
            continue
        end
        c2_fits[inner] = c2
        println(io, inner, '\t', c1.fit_sha256, '\t', c2.fit_sha256)
        for scan in campaign.scans
            scan.date == inner || continue
            denominator = _synthetic_observation_count(scan)
            gains[scan.scan] = (_synthetic_node_loglik(c2, scan, config) -
                                _synthetic_node_loglik(c1, scan, config)) / denominator
            println(io, scan.scan, '\t', bitstring(gains[scan.scan]))
        end
    end
    lower = -Inf
    selected = "C1"
    if valid_c2
        heldout_scans = [scan for scan in campaign.scans if scan.date in inner_dates]
        date_positive = all(mean(gains[scan.scan] for scan in heldout_scans
                                 if scan.date == date) > 0.0 for date in inner_dates)
        lower, bootstrap = _whole_scan_bootstrap_lower(
            gains, heldout_scans, inner_dates, config, "unary:$outer")
        println(io, "bootstrap_lower=", bitstring(lower))
        println(io, "bootstrap_hash=", sha256_hex(collect(codeunits(
            join(bitstring.(bootstrap), '\n') * "\n"))))
        selected = date_positive && lower > 0.0 ? "C2" : "C1"
    end
    println(io, "selected=", selected)
    selected_fits = selected == "C2" ? c2_fits : c1_fits
    return (status="OK", model=selected, fits=selected_fits,
            hash=sha256_hex(take!(io)), lower=lower)
end

struct _SyntheticEdgeNullFit
    mean::Vector{Float64}
    scale::Matrix{Float64}
end

struct _SyntheticEdgeFit
    means::Matrix{Float64}
    scale::Matrix{Float64}
    null:: _SyntheticEdgeNullFit
    sha256::String
end

function _project_edge_scale(scale::AbstractMatrix, config::DiagnosticConfig)
    size(scale) == (2, 2) || _fail(:numerical_failure, "edge scale is not 2 by 2")
    symmetric = Symmetric((Matrix{Float64}(scale) + Matrix{Float64}(scale)') / 2.0)
    decomposition = eigen(symmetric)
    all(isfinite, decomposition.values) ||
        _fail(:numerical_failure, "edge scale eigenvalues are nonfinite")
    values = max.(decomposition.values, config.scale_eigenvalue_floor)
    minimum_value = minimum(values)
    values = min.(values, minimum_value * config.scale_condition_cap)
    projected = decomposition.vectors * Diagonal(values) * decomposition.vectors'
    return Matrix{Float64}(Symmetric((projected + projected') / 2.0))
end

function _edge_quadratic(value::AbstractVector, mean_value::AbstractVector,
                         scale::AbstractMatrix)
    length(value) == length(mean_value) == 2 && size(scale) == (2, 2) ||
        _fail(:numerical_failure, "edge density dimensions differ")
    first = Float64(scale[1, 1])
    cross = 0.5 * (Float64(scale[1, 2]) + Float64(scale[2, 1]))
    second = Float64(scale[2, 2])
    determinant = first * second - cross * cross
    isfinite(determinant) && determinant > 0.0 ||
        _fail(:numerical_failure, "edge scale determinant is nonpositive")
    dx = Float64(value[1]) - Float64(mean_value[1])
    dy = Float64(value[2]) - Float64(mean_value[2])
    quadratic = (second * dx * dx - 2.0 * cross * dx * dy + first * dy * dy) /
                determinant
    return max(quadratic, 0.0), determinant
end

function _log_student_t_full2(value::AbstractVector, mean_value::AbstractVector,
                              scale::AbstractMatrix)
    quadratic, determinant = _edge_quadratic(value, mean_value, scale)
    nu = Float64(SRE.FIXED_NU)
    log_normalizer = SRE._loggamma_integer_or_half_twice(SRE.FIXED_NU + 2) -
                     SRE._loggamma_integer_or_half_twice(SRE.FIXED_NU) -
                     log(nu * π) - 0.5 * log(determinant)
    result = log_normalizer - 0.5 * (nu + 2.0) * log1p(quadratic / nu)
    isfinite(result) || _fail(:numerical_failure, "edge density is nonfinite")
    return result
end

function _fit_edge_null(values::Matrix{Float64}, config::DiagnosticConfig)
    size(values, 1) >= 2 && size(values, 2) == 2 && all(isfinite, values) ||
        return nothing
    count = size(values, 1)
    mean_value = [median(view(values, :, column)) for column in 1:2]
    centered = values .- reshape(mean_value, 1, :)
    scale = _project_edge_scale((centered' * centered) / count, config)
    for _ in 1:SRE.MAX_ITER
        latent = Vector{Float64}(undef, count)
        for row in axes(values, 1)
            quadratic, _ = _edge_quadratic(view(values, row, :), mean_value, scale)
            latent[row] = (SRE.FIXED_NU + 2.0) / (SRE.FIXED_NU + quadratic)
        end
        mass = sum(latent)
        isfinite(mass) && mass > 0.0 || return nothing
        proposed_mean = [sum(latent[row] * values[row, column]
                             for row in axes(values, 1)) / mass for column in 1:2]
        proposed_scale = zeros(Float64, 2, 2)
        for row in axes(values, 1)
            difference = collect(view(values, row, :)) .- proposed_mean
            proposed_scale .+= latent[row] .* (difference * difference')
        end
        proposed_scale = _project_edge_scale(proposed_scale / count, config)
        change = max(maximum(abs, proposed_mean .- mean_value),
                     maximum(abs, proposed_scale .- scale))
        mean_value = proposed_mean
        scale = proposed_scale
        change <= SRE.CONVERGENCE_TOL && break
    end
    all(isfinite, mean_value) && all(isfinite, scale) || return nothing
    return _SyntheticEdgeNullFit(mean_value, scale)
end

function _edge_pair_weights(probabilities::Matrix{Float64})
    size(probabilities, 2) == 2 ||
        _fail(:numerical_failure, "synthetic unary probabilities are not binary")
    edges = max(size(probabilities, 1) - 1, 0)
    weights = Matrix{Float64}(undef, edges, 4)
    for edge in 1:edges
        left_zero, left_one = probabilities[edge, 1], probabilities[edge, 2]
        right_zero, right_one = probabilities[edge + 1, 1], probabilities[edge + 1, 2]
        weights[edge, 1] = left_zero * right_zero
        weights[edge, 2] = left_zero * right_one
        weights[edge, 3] = left_one * right_zero
        weights[edge, 4] = left_one * right_one
    end
    all(isfinite, weights) || _fail(:numerical_failure, "edge weights are nonfinite")
    return weights
end

function _edge_training_material(campaign::SyntheticNullCampaign,
                                 training_dates::Vector{String},
                                 unary::_SyntheticUnaryFit,
                                 config::DiagnosticConfig)
    scans = _synthetic_scans(campaign, training_dates)
    values = reduce(vcat, (scan.edge_values for scan in scans))
    weight_blocks = Matrix{Float64}[]
    support = Dict(name => Dict(date => Set{String}() for date in training_dates)
                   for name in CATEGORY_NAMES)
    for scan in scans
        probabilities = _synthetic_posteriors(_synthetic_state_logs(unary, scan, config))
        weights = _edge_pair_weights(probabilities)
        size(weights, 1) == size(scan.edge_values, 1) ||
            _fail(:incomplete_fold, "edge and endpoint rows differ")
        push!(weight_blocks, weights)
        for category in eachindex(CATEGORY_NAMES)
            sum(view(weights, :, category)) > 0.0 &&
                push!(support[CATEGORY_NAMES[category]][scan.date], scan.scan)
        end
    end
    weights = reduce(vcat, weight_blocks)
    return scans, values, weights, support
end

function _edge_support_is_sufficient(campaign::SyntheticNullCampaign,
                                     training_dates::Vector{String},
                                     config::DiagnosticConfig)
    scans = _synthetic_scans(campaign, training_dates)
    edge_count = sum(size(scan.edge_values, 1) for scan in scans)
    return length(training_dates) >= config.minimum_dates &&
           length(scans) >= config.minimum_scans &&
           edge_count >= config.minimum_edges &&
           edge_count >= ceil(Int, config.observations_per_parameter *
                                    edge_admission_parameter_count())
end

function _edge_objective(values::Matrix{Float64}, weights::Matrix{Float64},
                         means::Matrix{Float64}, scale::Matrix{Float64})
    size(values, 1) == size(weights, 1) && size(weights, 2) == 4 &&
        size(means) == (4, 2) || _fail(:numerical_failure, "edge objective dimensions differ")
    total = 0.0
    for row in axes(values, 1), category in 1:4
        total += weights[row, category] * _log_student_t_full2(
            view(values, row, :), view(means, category, :), scale)
    end
    return total
end

function _expanded_edge_means(group_means::Matrix{Float64})
    size(group_means) == (3, 2) ||
        _fail(:numerical_failure, "edge group means differ")
    return vcat(group_means[1:1, :], group_means[2:2, :],
                group_means[2:2, :], group_means[3:3, :])
end

function _run_edge_start(values::Matrix{Float64}, weights::Matrix{Float64},
                         null::_SyntheticEdgeNullFit, alpha::Float64,
                         config::DiagnosticConfig)
    group_weights = hcat(weights[:, 1], weights[:, 2] + weights[:, 3], weights[:, 4])
    empirical = Matrix{Float64}(undef, 3, 2)
    for group in 1:3
        mass = sum(view(group_weights, :, group))
        isfinite(mass) && mass > 0.0 || return nothing
        for column in 1:2
            empirical[group, column] = sum(group_weights[row, group] * values[row, column]
                                            for row in axes(values, 1)) / mass
        end
    end
    group_means = (1.0 - alpha) .* repeat(reshape(null.mean, 1, :), 3, 1) .+
                  alpha .* empirical
    scale = copy(null.scale)
    means = _expanded_edge_means(group_means)
    objective = _edge_objective(values, weights, means, scale)
    isfinite(objective) || return nothing
    converged = false
    for _ in 1:SRE.MAX_ITER
        latent = Matrix{Float64}(undef, size(values, 1), 3)
        for row in axes(values, 1), group in 1:3
            quadratic, _ = _edge_quadratic(
                view(values, row, :), view(group_means, group, :), scale)
            latent[row, group] = (SRE.FIXED_NU + 2.0) /
                                 (SRE.FIXED_NU + quadratic)
        end
        proposed_means = similar(group_means)
        for group in 1:3
            effective = view(group_weights, :, group) .* view(latent, :, group)
            mass = sum(effective)
            isfinite(mass) && mass > 0.0 || return nothing
            for column in 1:2
                proposed_means[group, column] =
                    sum(effective[row] * values[row, column]
                        for row in axes(values, 1)) / mass
            end
        end
        proposed_scale = zeros(Float64, 2, 2)
        for row in axes(values, 1), group in 1:3
            difference = collect(view(values, row, :)) .-
                         collect(view(proposed_means, group, :))
            proposed_scale .+= group_weights[row, group] * latent[row, group] .*
                              (difference * difference')
        end
        proposed_scale ./= size(values, 1)
        shrinkage = min(0.5, 3.0 / size(values, 1))
        proposed_scale = (1.0 - shrinkage) .* proposed_scale .+
                         shrinkage .* null.scale
        proposed_scale = _project_edge_scale(proposed_scale, config)
        proposed_expanded = _expanded_edge_means(proposed_means)
        proposed_objective = _edge_objective(
            values, weights, proposed_expanded, proposed_scale)
        proposed_objective + SRE.DECREASE_TOL >= objective || return nothing
        change = abs(proposed_objective - objective)
        group_means = proposed_means
        means = proposed_expanded
        scale = proposed_scale
        objective = proposed_objective
        if change <= SRE.CONVERGENCE_TOL * (1.0 + abs(objective))
            converged = true
            break
        end
    end
    converged || return nothing
    return (means=means, scale=scale, objective=objective)
end

function _synthetic_edge_hash(means::Matrix{Float64}, scale::Matrix{Float64},
                              null::_SyntheticEdgeNullFit,
                              training_sha256::String)
    io = IOBuffer()
    println(io, "structured-null-edge-fit-v1")
    println(io, "training=", training_sha256)
    println(io, "pair_prior=", join(bitstring.(SYNTHETIC_PAIR_PRIOR), ','))
    for value in means
        println(io, "mean\t", bitstring(value))
    end
    for value in scale
        println(io, "scale\t", bitstring(value))
    end
    for value in null.mean
        println(io, "null_mean\t", bitstring(value))
    end
    for value in null.scale
        println(io, "null_scale\t", bitstring(value))
    end
    return sha256_hex(take!(io))
end

function _fit_synthetic_edge(campaign::SyntheticNullCampaign,
                             training_dates::Vector{String},
                             unary::_SyntheticUnaryFit,
                             config::DiagnosticConfig)
    _edge_support_is_sufficient(campaign, training_dates, config) || return nothing
    scans, values, weights, support = _edge_training_material(
        campaign, training_dates, unary, config)
    for category in eachindex(CATEGORY_NAMES)
        sum_weights = sum(view(weights, :, category))
        sum_squared = sum(abs2, view(weights, :, category))
        kish_ess(sum_weights, sum_squared) >= config.minimum_category_ess || return nothing
        all(length(support[CATEGORY_NAMES[category]][date]) >=
            config.minimum_support_scans for date in training_dates) || return nothing
    end
    null = _fit_edge_null(values, config)
    null === nothing && return nothing
    best = nothing
    for alpha in (0.0, 0.25, 0.5, 0.75, 1.0)
        candidate = _run_edge_start(values, weights, null, alpha, config)
        candidate === nothing && continue
        if best === nothing || candidate.objective > best.objective + 1.0e-12
            best = candidate
        end
    end
    best === nothing && return nothing
    eigenvalues = eigvals(Symmetric(best.scale))
    scale_is_feasible(eigenvalues, config.scale_eigenvalue_floor,
                      config.scale_condition_cap) || return nothing
    training_sha256 = _synthetic_training_hash(scans)
    digest = _synthetic_edge_hash(best.means, best.scale, null, training_sha256)
    return _SyntheticEdgeFit(best.means, best.scale, null, digest)
end

function _edge_gain_with_weights(fit::_SyntheticEdgeFit, scan::SyntheticNullScan,
                                 weights::Matrix{Float64})
    size(weights) == (size(scan.edge_values, 1), 4) ||
        _fail(:incomplete_fold, "held-out edge weights differ")
    gains = Vector{Float64}(undef, size(scan.edge_values, 1))
    for edge in axes(scan.edge_values, 1)
        null_score = _log_student_t_full2(
            view(scan.edge_values, edge, :), fit.null.mean, fit.null.scale)
        conditional = sum(weights[edge, category] * _log_student_t_full2(
            view(scan.edge_values, edge, :), view(fit.means, category, :), fit.scale)
                          for category in 1:4)
        gains[edge] = conditional - null_score
    end
    return gains
end

function _heldout_edge_weights(unary::_SyntheticUnaryFit, scan::SyntheticNullScan,
                               config::DiagnosticConfig)
    probabilities = _synthetic_posteriors(_synthetic_state_logs(unary, scan, config))
    return _edge_pair_weights(probabilities)
end

function _edge_reversal_equivariant(fit::_SyntheticEdgeFit,
                                    scan::SyntheticNullScan,
                                    weights::Matrix{Float64})
    forward = sum(_edge_gain_with_weights(fit, scan, weights))
    reversed_scan = SyntheticNullScan(
        scan.date, scan.scan, reverse(scan.values; dims=1), reverse(scan.amplitudes),
        reverse(scan.edge_values; dims=1))
    reversed_weights = reverse(weights[:, [1, 3, 2, 4]]; dims=1)
    reversed = sum(_edge_gain_with_weights(fit, reversed_scan, reversed_weights))
    return isapprox(forward, reversed; atol=1.0e-12, rtol=1.0e-12)
end

function _synthetic_edge_admission(campaign::SyntheticNullCampaign, outer::String,
                                   unary_decision, config::DiagnosticConfig)
    outer_training_dates = [date for date in campaign.dates if date != outer]
    io = IOBuffer()
    println(io, "structured-null-inner-edge-v1")
    println(io, "outer=", outer)
    println(io, "unary=", unary_decision.model)
    println(io, "pair_prior=", join(bitstring.(SYNTHETIC_PAIR_PRIOR), ','))
    if !_edge_support_is_sufficient(campaign, outer_training_dates, config)
        println(io, "decision=DISABLED\treason=outer_training_feasibility")
        return (admitted=false, hash=sha256_hex(take!(io)))
    end

    gains = Dict{String,Float64}()
    fits = Dict{String,_SyntheticEdgeFit}()
    weights_by_scan = Dict{String,Matrix{Float64}}()
    scans_by_name = Dict(scan.scan => scan for scan in campaign.scans)
    reversal_ok = true
    for inner in outer_training_dates
        training_dates = [date for date in outer_training_dates if date != inner]
        _edge_support_is_sufficient(campaign, training_dates, config) || begin
            println(io, "decision=DISABLED\treason=inner_training_feasibility\tinner=", inner)
            return (admitted=false, hash=sha256_hex(take!(io)))
        end
        unary_fit = get(unary_decision.fits, inner, nothing)
        unary_fit === nothing && begin
            println(io, "decision=DISABLED\treason=missing_inner_unary\tinner=", inner)
            return (admitted=false, hash=sha256_hex(take!(io)))
        end
        edge_fit = _fit_synthetic_edge(campaign, training_dates, unary_fit, config)
        edge_fit === nothing && begin
            println(io, "decision=DISABLED\treason=edge_fit_failed\tinner=", inner)
            return (admitted=false, hash=sha256_hex(take!(io)))
        end
        println(io, inner, '\t', unary_fit.fit_sha256, '\t', edge_fit.sha256)
        for scan in campaign.scans
            scan.date == inner || continue
            weights = _heldout_edge_weights(unary_fit, scan, config)
            scan_gain = sum(_edge_gain_with_weights(edge_fit, scan, weights)) /
                        _synthetic_observation_count(scan)
            gains[scan.scan] = scan_gain
            fits[scan.scan] = edge_fit
            weights_by_scan[scan.scan] = weights
            reversal_ok &= _edge_reversal_equivariant(edge_fit, scan, weights)
            println(io, scan.scan, '\t', bitstring(scan_gain))
        end
    end
    heldout_scans = [scan for scan in campaign.scans if scan.date in outer_training_dates]
    date_positive = all(mean(gains[scan.scan] for scan in heldout_scans
                             if scan.date == date) > 0.0 for date in outer_training_dates)
    lower, bootstrap = _whole_scan_bootstrap_lower(
        gains, heldout_scans, outer_training_dates, config, "edge:$outer")
    observed = mean(mean(gains[scan.scan] for scan in heldout_scans
                         if scan.date == date) for date in outer_training_dates)
    shuffle_values = Float64[]
    for seed in config.bootstrap_seeds
        rng = _seeded_rng(seed, "edge-shuffle", outer)
        shuffled = Dict{String,Float64}()
        for scan in heldout_scans
            weights = weights_by_scan[scan.scan]
            edge_count = size(weights, 1)
            shift = edge_count <= 1 ? 0 : rand(rng, 0:(edge_count - 1))
            shifted = circshift(weights, (shift, 0))
            shuffled[scan.scan] = sum(_edge_gain_with_weights(
                fits[scan.scan], scans_by_name[scan.scan], shifted)) /
                _synthetic_observation_count(scan)
        end
        push!(shuffle_values, mean(mean(shuffled[scan.scan] for scan in heldout_scans
                                         if scan.date == date)
                                   for date in outer_training_dates))
    end
    shuffle_threshold = quantile_type7(shuffle_values, config.permutation_quantile)
    admitted = date_positive && lower > 0.0 && observed > shuffle_threshold && reversal_ok
    println(io, "date_positive=", date_positive)
    println(io, "bootstrap_lower=", bitstring(lower))
    println(io, "bootstrap_hash=", sha256_hex(collect(codeunits(
        join(bitstring.(bootstrap), '\n') * "\n"))))
    println(io, "observed=", bitstring(observed))
    println(io, "shuffle_threshold=", bitstring(shuffle_threshold))
    println(io, "shuffle_hash=", sha256_hex(collect(codeunits(
        join(bitstring.(shuffle_values), '\n') * "\n"))))
    println(io, "reversal=", reversal_ok)
    println(io, "decision=", admitted ? "ENABLED" : "DISABLED")
    return (admitted=admitted, hash=sha256_hex(take!(io)))
end

function _edge_null_loglik(fit::_SyntheticEdgeFit, scan::SyntheticNullScan)
    return sum(_log_student_t_full2(view(scan.edge_values, edge, :),
                                    fit.null.mean, fit.null.scale)
               for edge in axes(scan.edge_values, 1))
end

function _synthetic_graph_loglik(unary::_SyntheticUnaryFit,
                                 edge_fit::_SyntheticEdgeFit,
                                 scan::SyntheticNullScan,
                                 config::DiagnosticConfig)
    SYNTHETIC_PAIR_PRIOR == (0.0, 0.0, 0.0, 0.0) ||
        _fail(:config_contract, "synthetic graph pair prior is not zero")
    logs = _synthetic_state_logs(unary, scan, config)
    prior = log(0.5)
    dynamic = [prior + logs[1, 1], prior + logs[1, 2]]
    for edge in axes(scan.edge_values, 1)
        null_score = _log_student_t_full2(
            view(scan.edge_values, edge, :), edge_fit.null.mean, edge_fit.null.scale)
        factors = [_log_student_t_full2(
            view(scan.edge_values, edge, :), view(edge_fit.means, category, :),
            edge_fit.scale) - null_score for category in 1:4]
        next_zero = prior + logs[edge + 1, 1] + _logsumexp(
            dynamic[1] + factors[1] + SYNTHETIC_PAIR_PRIOR[1],
            dynamic[2] + factors[3] + SYNTHETIC_PAIR_PRIOR[3])
        next_one = prior + logs[edge + 1, 2] + _logsumexp(
            dynamic[1] + factors[2] + SYNTHETIC_PAIR_PRIOR[2],
            dynamic[2] + factors[4] + SYNTHETIC_PAIR_PRIOR[4])
        dynamic = [next_zero, next_one]
    end
    return _logsumexp(dynamic[1], dynamic[2]) + _edge_null_loglik(edge_fit, scan)
end

function _exact_sign_flip_p(values::Vector{Float64})
    isempty(values) && return 1.0
    all(isfinite, values) || return 1.0
    observed = mean(values)
    tail = 0
    assignments = 1 << length(values)
    for mask in 0:(assignments - 1)
        statistic = mean(((mask >> (index - 1)) & 1 == 1 ? 1.0 : -1.0) *
                         values[index] for index in eachindex(values))
        statistic >= observed && (tail += 1)
    end
    return tail / assignments
end

function _failed_synthetic_outer(outer::String, training_sha256::String,
                                 inner_hash::String, refit_hash::String,
                                 unary_model::String, fit_invocations::Int,
                                 reason::String)
    material = "structured-null-outer-failure-v1\nouter=$outer\nreason=$reason\n" *
               "training=$training_sha256\ninner=$inner_hash\nrefit=$refit_hash\n"
    return SyntheticOuterDecision(
        outer, training_sha256, inner_hash, refit_hash, unary_model, false, false,
        fit_invocations, Float64[], -Inf, -Inf, "FAIL",
        sha256_hex(collect(codeunits(material))))
end

function _synthetic_outer_decision(campaign::SyntheticNullCampaign, outer::String,
                                   config::DiagnosticConfig)
    training_dates = [date for date in campaign.dates if date != outer]
    training_scans = _synthetic_scans(campaign, training_dates)
    training_sha256 = _synthetic_training_hash(training_scans)
    fit_counter = Ref(0)
    unary_decision = _select_synthetic_unary(campaign, outer, config, fit_counter)
    if unary_decision.status != "OK"
        return _failed_synthetic_outer(
            outer, training_sha256, unary_decision.hash, unary_decision.hash, "C1",
            fit_counter[], "inner_unary_failure")
    end
    c1 = _fit_synthetic_unary(campaign, training_dates, "C1", config, fit_counter)
    c1 === nothing && return _failed_synthetic_outer(
        outer, training_sha256, unary_decision.hash, unary_decision.hash,
        unary_decision.model, fit_counter[], "outer_c1_failure")
    selected = unary_decision.model == "C1" ? c1 :
        _fit_synthetic_unary(campaign, training_dates, "C2", config, fit_counter)
    selected === nothing && return _failed_synthetic_outer(
        outer, training_sha256, unary_decision.hash, c1.fit_sha256,
        unary_decision.model, fit_counter[], "outer_selected_unary_failure")

    edge_decision = _synthetic_edge_admission(campaign, outer, unary_decision, config)
    graph_enabled = edge_decision.admitted
    edge_fit = graph_enabled ?
        _fit_synthetic_edge(campaign, training_dates, selected, config) : nothing
    graph_enabled && edge_fit === nothing && return _failed_synthetic_outer(
        outer, training_sha256, unary_decision.hash, selected.fit_sha256,
        unary_decision.model, fit_counter[], "outer_edge_refit_failure")
    refit_material = "structured-null-outer-refit-v1\nouter=$outer\n" *
                     "training=$training_sha256\nc1=$(c1.fit_sha256)\n" *
                     "selected=$(selected.fit_sha256)\nedge=$(edge_decision.hash)\n"
    refit_hash = sha256_hex(collect(codeunits(refit_material)))
    scan_scores = Float64[]
    heldout_scans = [scan for scan in campaign.scans if scan.date == outer]
    gains = Dict{String,Float64}()
    for scan in heldout_scans
        difference = if unary_decision.model == "C1" && !graph_enabled
            0.0
        elseif graph_enabled
            reference = _synthetic_node_loglik(c1, scan, config) +
                        _edge_null_loglik(edge_fit, scan)
            _synthetic_graph_loglik(selected, edge_fit, scan, config) - reference
        else
            _synthetic_node_loglik(selected, scan, config) -
            _synthetic_node_loglik(c1, scan, config)
        end
        score = difference / _synthetic_observation_count(scan)
        isfinite(score) || return _failed_synthetic_outer(
            outer, training_sha256, unary_decision.hash, refit_hash,
            unary_decision.model, fit_counter[], "outer_score_nonfinite")
        gains[scan.scan] = score
        push!(scan_scores, score)
    end
    mean_score = mean(scan_scores)
    bootstrap_lower, bootstrap = _whole_scan_bootstrap_lower(
        gains, heldout_scans, [outer], config, "outer:$outer")
    status = mean_score > 0.0 && bootstrap_lower > 0.0 ? "PASS" : "FAIL"
    io = IOBuffer()
    println(io, "structured-null-outer-decision-v1")
    println(io, "outer=", outer)
    println(io, "training=", training_sha256)
    println(io, "inner=", unary_decision.hash)
    println(io, "refit=", refit_hash)
    println(io, "unary=", unary_decision.model)
    println(io, "edge_admitted=", edge_decision.admitted)
    println(io, "graph_enabled=", graph_enabled)
    println(io, "pair_prior=", join(bitstring.(SYNTHETIC_PAIR_PRIOR), ','))
    for score in scan_scores
        println(io, "scan_score=", bitstring(score))
    end
    println(io, "mean=", bitstring(mean_score))
    println(io, "bootstrap_lower=", bitstring(bootstrap_lower))
    println(io, "bootstrap_hash=", sha256_hex(collect(codeunits(
        join(bitstring.(bootstrap), '\n') * "\n"))))
    println(io, "status=", status)
    return SyntheticOuterDecision(
        outer, training_sha256, unary_decision.hash, refit_hash,
        unary_decision.model, edge_decision.admitted, graph_enabled, fit_counter[],
        scan_scores, mean_score, bootstrap_lower, status, sha256_hex(take!(io)))
end

function _complete_nested_meta(campaign::SyntheticNullCampaign,
                               config::DiagnosticConfig;
                               outer_dates::Vector{String}=copy(campaign.dates))
    requested = sort!(unique(copy(outer_dates)))
    all(date -> date in campaign.dates, requested) ||
        _fail(:malformed_fold, "synthetic outer date is unknown")
    length(requested) == length(outer_dates) ||
        _fail(:duplicate_fold, "synthetic outer dates repeat")
    folds = [_synthetic_outer_decision(campaign, outer, config) for outer in requested]
    date_scores = [fold.mean_score for fold in folds]
    sign_flip_p = _exact_sign_flip_p(date_scores)
    complete = requested == campaign.dates
    terminal_pass = complete && all(fold.status == "PASS" for fold in folds) &&
                    sign_flip_p < config.alpha
    status = terminal_pass ? "PASS" : "FAIL"
    io = IOBuffer()
    println(io, "structured-null-complete-meta-v1")
    println(io, "seed=", campaign.seed)
    println(io, "complete=", complete)
    println(io, "pair_prior=", join(bitstring.(SYNTHETIC_PAIR_PRIOR), ','))
    for fold in folds
        println(io, fold.outer_heldout_date, '\t', fold.evidence_sha256)
    end
    println(io, "sign_flip_p=", bitstring(sign_flip_p))
    println(io, "status=", status)
    return SyntheticMetaDecision(status, folds, date_scores, sign_flip_p,
                                 SYNTHETIC_PAIR_PRIOR, sha256_hex(take!(io)))
end

function _classify_null_campaign(statuses::Vector{String}, meta_status::String)
    all(status -> status in TERMINAL_STATES, statuses) ||
        _fail(:malformed_input, "diagnostic null status is invalid")
    meta_status in ("PASS", "FAIL") ||
        _fail(:malformed_input, "complete-meta null status is invalid")
    return (diagnostic_trigger=any(==("FOLLOWUP_REQUIRED"), statuses),
            meta_false_pass=meta_status == "PASS")
end

function _meta_failure_count(statuses::AbstractVector{<:AbstractString})
    all(status -> String(status) in ("PASS", "FAIL"), statuses) ||
        _fail(:malformed_input, "complete-meta status collection is invalid")
    return count(status -> String(status) == "PASS", statuses)
end

function _null_support(dates::Vector{String}, scans_per_date::Int, nodes::Int,
                       edges::Int, category_counts::Dict{String,Tuple{Float64,Float64}})
    categories = EdgeCategorySupport[]
    for name in CATEGORY_NAMES
        sums = category_counts[name]
        support = Dict(date => scans_per_date for date in dates)
        push!(categories, EdgeCategorySupport(name, sums[1], sums[2], support))
    end
    return FeasibilityInput(dates, length(dates) * scans_per_date, nodes, edges,
                            categories, [0.7, 1.3])
end

function _generate_null_campaign(seed::Int, config::DiagnosticConfig)
    rng = MersenneTwister(seed)
    dates = [@sprintf("%08d", 20300101 + index - 1) for index in 1:config.null_dates]
    config.null_scans % config.null_dates == 0 ||
        _fail(:config_contract, "null scans must divide evenly across dates")
    scans_per_date = config.null_scans ÷ config.null_dates
    scans = SyntheticNullScan[]
    data_material = IOBuffer()
    edge_material = IOBuffer()
    latent_ones = 0
    edges = 0
    for date in dates
        for scan_index in 1:scans_per_date
            scan = "$(date)_null_$(lpad(scan_index, 2, '0')).sxm"
            chain_length = rand(rng, config.null_chain_min:config.null_chain_max)
            p = length(config.feature_names)
            base = Matrix{Float64}(undef, chain_length, p)
            shifts = 0.5 .* randn(rng, p)
            scales = exp.(0.1 .* randn(rng, p))
            for row in 1:chain_length
                draw = _diagonal_student_t8_row(rng, p)
                for column in 1:p
                    base[row, column] = shifts[column] + scales[column] * draw[column]
                end
            end
            normalized = _normalize_scan(base)
            for value in base
                println(data_material, bitstring(value))
            end
            states = rand(rng, Bool, chain_length)
            latent_ones += count(identity, states)
            edge_values = Matrix{Float64}(undef, chain_length - 1, 2)
            for edge in 1:(chain_length - 1)
                first, second = _bivariate_student_t8_covariance(rng)
                edge_values[edge, 1] = first
                edge_values[edge, 2] = second
                println(edge_material, bitstring(first), '\t', bitstring(second))
                edges += 1
            end
            push!(scans, SyntheticNullScan(
                date, scan, normalized, copy(base[:, 1]), edge_values))
        end
    end
    nodes = sum(size(scan.values, 1) for scan in scans)
    return SyntheticNullCampaign(
        seed, dates, scans, nodes, edges, latent_ones,
        sha256_hex(take!(data_material)), sha256_hex(take!(edge_material)))
end

function _campaign(seed::Int, config::DiagnosticConfig)
    campaign = _generate_null_campaign(seed, config)
    scans_per_date = config.null_scans ÷ config.null_dates
    residual_scans = [ResidualScan(scan.date, scan.scan, scan.values)
                      for scan in campaign.scans]
    view_differences = Dict(scan.scan => _null_view_difference(scan.values, config)
                            for scan in campaign.scans)
    dropout_differences = Dict(
        scan.scan => _null_dropout_difference(scan.values, config)
        for scan in campaign.scans)
    category_sum = Dict(name => 0.25 * campaign.edges for name in CATEGORY_NAMES)
    category_sum_squared = Dict(name => 0.25^2 * campaign.edges
                                for name in CATEGORY_NAMES)
    category_counts = Dict(name => (category_sum[name], category_sum_squared[name])
                           for name in CATEGORY_NAMES)
    support = _null_support(campaign.dates, scans_per_date, campaign.nodes,
                            campaign.edges, category_counts)
    evidence = PartitionEvidence("synthetic-null-$seed", residual_scans,
                                 view_differences, dropout_differences)
    result = evaluate_partition(evidence, support, config)
    statuses = [mechanism.status for mechanism in result.mechanisms]
    meta = _complete_nested_meta(campaign, config)
    classification = _classify_null_campaign(statuses, meta.status)
    result_material = "structured-null-campaign-v2\nseed=$seed\n" *
                      "nodes=$(campaign.nodes)\nedges=$(campaign.edges)\n" *
                      "latent_ones=$(campaign.latent_ones)\n" *
                      "statuses=$(join(statuses, ','))\npartition=$(result.evidence_sha256)\n" *
                      "diagnostic_trigger=$(classification.diagnostic_trigger)\n" *
                      "meta_status=$(meta.status)\n" *
                      "meta_false_pass=$(classification.meta_false_pass)\n" *
                      "meta=$(meta.sha256)\n" *
                      "data=$(campaign.data_sha256)\nedge=$(campaign.edge_sha256)\n"
    result_hash = sha256_hex(collect(codeunits(result_material)))
    return NullCampaignRecord(
        seed, campaign.nodes, campaign.edges, campaign.latent_ones, statuses,
        classification.diagnostic_trigger, meta.status, classification.meta_false_pass,
        meta.sha256, campaign.data_sha256, campaign.edge_sha256, result_hash)
end

function null_campaign_hash(records::Vector{NullCampaignRecord})
    io = IOBuffer()
    println(io, "structured-null-campaigns-v2")
    for record in records
        println(io, join((record.seed, record.nodes, record.edges, record.latent_ones,
                          join(record.statuses, ','), record.diagnostic_trigger,
                          record.meta_status, record.meta_false_pass, record.meta_sha256,
                          record.data_sha256, record.edge_sha256,
                          record.result_sha256), '\t'))
    end
    return sha256_hex(take!(io))
end

function run_null_campaigns(config::DiagnosticConfig)
    records = [_campaign(seed, config) for seed in config.null_campaign_seeds]
    failures = _meta_failure_count([record.meta_status for record in records])
    upper = clopper_pearson_upper(failures, length(records),
                                  config.clopper_confidence)
    accepted = upper <= config.clopper_upper_max
    digest = null_campaign_hash(records)
    return NullCampaignSummary(records, failures, length(records), upper, accepted, digest)
end

function _strict_lines(snapshot::FileSnapshot, context::String)
    bytes = snapshot.bytes
    isvalid(String, bytes) || _fail(:invalid_encoding, "$context is not UTF-8")
    isempty(bytes) && _fail(:empty_input, "$context is empty")
    last(bytes) == UInt8('\n') || _fail(:malformed_input, "$context lacks final LF")
    UInt8('\r') in bytes && _fail(:malformed_input, "$context contains CR")
    UInt8('\0') in bytes && _fail(:malformed_input, "$context contains NUL")
    lines = split(String(copy(bytes)), '\n'; keepempty=true)
    pop!(lines) == "" || _fail(:malformed_input, "$context termination differs")
    any(isempty, lines) && _fail(:malformed_input, "$context contains a blank row")
    return lines
end

function _parse_table(snapshot::FileSnapshot, expected_columns::Vector{String},
                      context::String)
    lines = _strict_lines(snapshot, context)
    length(lines) >= 2 || _fail(:empty_input, "$context contains no rows")
    header = String.(split(first(lines), '\t'; keepempty=true))
    header == expected_columns || _fail(:invalid_columns, "$context header differs")
    length(header) == length(unique(header)) ||
        _fail(:duplicate_column, "$context repeats a column")
    rows = Dict{String,String}[]
    for (index, line) in enumerate(lines[2:end])
        fields = String.(split(line, '\t'; keepempty=true))
        length(fields) == length(header) ||
            _fail(:invalid_row_columns, "$context row $index has the wrong width")
        any(value -> occursin('\t', value) || occursin('\n', value), fields) &&
            _fail(:malformed_input, "$context row $index contains a delimiter")
        push!(rows, Dict(header[position] => fields[position]
                         for position in eachindex(header)))
    end
    return rows
end

function _parse_date(scan::String)
    basename(scan) == scan || _fail(:path_value, "scan key must be a basename")
    matched = match(r"^(\d{8})(?:\D|$)", scan)
    matched === nothing && _fail(:missing_date, "scan key lacks one leading date")
    return String(only(matched.captures))
end

function _parse_int(text::AbstractString, context::String; minimum::Int=typemin(Int))
    value = tryparse(Int, String(text))
    value === nothing && _fail(:malformed_input, "$context is not an integer")
    value >= minimum || _fail(:malformed_input, "$context is below its minimum")
    return value
end

function _parse_float(text::AbstractString, context::String; minimum::Float64=-Inf)
    value = tryparse(Float64, String(text))
    value === nothing && _fail(:malformed_input, "$context is not numeric")
    isfinite(value) || _fail(:nonfinite_input, "$context is nonfinite")
    value >= minimum || _fail(:malformed_input, "$context is below its minimum")
    return value
end

function load_feature_dataset(path::AbstractString, config::DiagnosticConfig)
    snapshot = _snapshot(path)
    columns = vcat(["file", "lobe", "amplitude"], config.feature_names)
    rows = _parse_table(snapshot, columns, "feature table")
    keys_seen = Set{Tuple{String,Int}}()
    date_by_scan = Dict{String,String}()
    for (index, row) in enumerate(rows)
        scan = row["file"]
        isempty(scan) && _fail(:malformed_input, "feature row $index has an empty scan")
        date = _parse_date(scan)
        lobe = _parse_int(row["lobe"], "feature row $index lobe"; minimum=1)
        key = (scan, lobe)
        key in keys_seen && _fail(:duplicate_key, "feature table repeats a key")
        push!(keys_seen, key)
        date_by_scan[scan] = date
        _parse_float(row["amplitude"], "feature row $index amplitude")
        for feature in config.feature_names
            _parse_float(row[feature], "feature row $index $feature")
        end
    end
    ordered_keys = [(row["file"], _parse_int(row["lobe"], "feature lobe"; minimum=1))
                    for row in rows]
    ordered_keys == sort(ordered_keys) ||
        _fail(:malformed_input, "feature rows are not sorted by file and lobe")
    parsed = try
        HUA.load_records(snapshot.path)
    catch error
        _fail(:feature_load, sprint(showerror, error))
    end
    length(parsed) == length(rows) || _fail(:key_mismatch, "feature loader changed rows")
    Set((record.file, record.lobe) for record in parsed) == keys_seen ||
        _fail(:key_mismatch, "feature loader changed keys")
    _verify(snapshot)
    return FeatureDataset(parsed, date_by_scan, snapshot)
end

function load_fold_plan(path::AbstractString, dataset::FeatureDataset)
    snapshot = _snapshot(path)
    rows = _parse_table(snapshot, ["outer_heldout_date", "inner_heldout_date"],
                        "fold table")
    dates = sort!(unique(Base.values(dataset.date_by_scan)))
    pairs = Tuple{String,String}[]
    for (index, row) in enumerate(rows)
        outer = row["outer_heldout_date"]
        inner = row["inner_heldout_date"]
        occursin(r"^\d{8}$", outer) && occursin(r"^\d{8}$", inner) ||
            _fail(:malformed_fold, "fold row $index has an invalid date")
        outer != inner || _fail(:malformed_fold, "outer and inner dates coincide")
        outer in dates && inner in dates || _fail(:malformed_fold, "fold date is unknown")
        push!(pairs, (outer, inner))
    end
    pairs == sort(pairs) || _fail(:malformed_fold, "fold rows are not sorted")
    length(pairs) == length(unique(pairs)) || _fail(:duplicate_fold, "fold rows repeat")
    expected = sort([(outer, inner) for outer in dates for inner in dates if inner != outer])
    pairs == expected || _fail(:incomplete_fold, "fold topology is incomplete")
    _verify(snapshot)
    return FoldPlan(dates, pairs, snapshot)
end

function _parse_support(text::String, expected_dates::Vector{String}, context::String)
    entries = isempty(text) ? String[] : String.(split(text, ','; keepempty=true))
    support = Dict{String,Int}()
    encountered = String[]
    for entry in entries
        pair = split(entry, ':'; keepempty=true)
        length(pair) == 2 || _fail(:malformed_receipt, "$context support is malformed")
        date = pair[1]
        haskey(support, date) && _fail(:duplicate_receipt, "$context support repeats a date")
        support[date] = _parse_int(pair[2], "$context support count"; minimum=0)
        push!(encountered, date)
    end
    encountered == sort(expected_dates) ||
        _fail(:incomplete_receipt, "$context support dates differ")
    return support
end

function load_input_receipt(path::AbstractString, config::DiagnosticConfig,
                            dataset::FeatureDataset, folds::FoldPlan)
    snapshot = _snapshot(path)
    columns = [
        "outer_heldout_date", "category", "sum_weights", "sum_squared_weights",
        "distinct_eligible_edges", "scale_eigenvalue_1", "scale_eigenvalue_2",
        "support_scans_by_training_date", "config_sha256", "feature_sha256",
        "fold_sha256", "unary_review_sha256",
    ]
    rows = _parse_table(snapshot, columns, "diagnostic input receipt")
    unary = _snapshot(UNARY_REVIEW_PATH; expected_path=UNARY_REVIEW_PATH,
                      expected_hash=UNARY_REVIEW_SHA256)
    grouped = Dict{String,Vector{Dict{String,String}}}()
    for row in rows
        row["config_sha256"] == config.sha256 ||
            _fail(:stale_hash, "receipt config hash differs")
        row["feature_sha256"] == dataset.snapshot.sha256 ||
            _fail(:stale_hash, "receipt feature hash differs")
        row["fold_sha256"] == folds.snapshot.sha256 ||
            _fail(:stale_hash, "receipt fold hash differs")
        row["unary_review_sha256"] == unary.sha256 ||
            _fail(:stale_hash, "receipt unary dependency hash differs")
        push!(get!(grouped, row["outer_heldout_date"], Dict{String,String}[]), row)
    end
    [(row["outer_heldout_date"], row["category"]) for row in rows] ==
        sort([(row["outer_heldout_date"], row["category"]) for row in rows]) ||
        _fail(:malformed_receipt, "receipt rows are not sorted")
    sort!(collect(keys(grouped))) == folds.dates ||
        _fail(:incomplete_receipt, "receipt outer folds differ")
    by_outer = Dict{String,ReceiptOuter}()
    for outer in folds.dates
        outer_rows = grouped[outer]
        length(outer_rows) == length(CATEGORY_NAMES) ||
            _fail(:incomplete_receipt, "receipt category count differs")
        [row["category"] for row in outer_rows] == CATEGORY_NAMES ||
            _fail(:malformed_receipt, "receipt category order differs")
        training_dates = [date for date in folds.dates if date != outer]
        edges_values = Set(_parse_int(row["distinct_eligible_edges"], "edge count";
                                      minimum=0) for row in outer_rows)
        first_eigenvalues = Set(_parse_float(row["scale_eigenvalue_1"], "scale eigenvalue";
                                             minimum=0.0) for row in outer_rows)
        second_eigenvalues = Set(_parse_float(row["scale_eigenvalue_2"], "scale eigenvalue";
                                              minimum=0.0) for row in outer_rows)
        length(edges_values) == length(first_eigenvalues) == length(second_eigenvalues) == 1 ||
            _fail(:incompatible_receipt, "outer receipt declarations disagree")
        categories = EdgeCategorySupport[]
        for row in outer_rows
            support = _parse_support(row["support_scans_by_training_date"],
                                     training_dates, "receipt $(row["category"])")
            push!(categories, EdgeCategorySupport(
                row["category"],
                _parse_float(row["sum_weights"], "sum_weights"; minimum=0.0),
                _parse_float(row["sum_squared_weights"], "sum_squared_weights";
                             minimum=0.0),
                support,
            ))
        end
        edges = only(edges_values)
        total_weight = sum(category.sum_weights for category in categories)
        tolerance = 8 * eps(max(1.0, Float64(edges)))
        abs(total_weight - edges) <= tolerance ||
            _fail(:incompatible_receipt, "soft category weights do not sum to edges")
        actual_scans_by_date = Dict(date => length(unique(
            record.file for record in dataset.records
            if dataset.date_by_scan[record.file] == date)) for date in training_dates)
        for category in categories
            category.sum_weights <= edges + tolerance ||
                _fail(:incompatible_receipt, "category weight exceeds edge count")
            category.sum_squared_weights <= category.sum_weights^2 + tolerance ||
                _fail(:incompatible_receipt, "category squared weights are impossible")
            kish_ess(category.sum_weights, category.sum_squared_weights) <=
                edges + tolerance ||
                _fail(:incompatible_receipt, "category Kish ESS exceeds edges")
            for date in training_dates
                category.support_scans_by_date[date] <= actual_scans_by_date[date] ||
                    _fail(:incompatible_receipt, "category support exceeds scan count")
            end
        end
        eigenvalues = [only(first_eigenvalues), only(second_eigenvalues)]
        issorted(eigenvalues) ||
            _fail(:incompatible_receipt, "scale eigenvalues are not sorted")
        training_records = [record for record in dataset.records
                            if dataset.date_by_scan[record.file] != outer]
        maximum_edges = sum(max(count(record -> record.file == scan, training_records) - 1, 0)
                            for scan in unique(record.file for record in training_records))
        edges <= maximum_edges ||
            _fail(:incompatible_receipt, "eligible edges exceed consecutive pairs")
        by_outer[outer] = ReceiptOuter(
            outer, edges, eigenvalues, categories)
    end
    _verify(snapshot)
    _verify(unary)
    return InputReceipt(by_outer, snapshot, unary)
end

function _fit_residual_model(train::Matrix{Float64}, amplitudes::Vector{Float64})
    size(train, 1) == length(amplitudes) || _fail(:incomplete_fold, "fit rows differ")
    size(train, 1) >= 4 || _fail(:underpowered, "inner fit has too few rows")
    fit = SRE.fit_student_t_two_component(train, amplitudes)
    if fit.status == SRE.ROBUST_VALID
        return (:two, fit)
    elseif fit.one_component.status == SRE.ROBUST_VALID
        return (:one, fit.one_component)
    end
    _fail(:numerical_failure, "inner fixed-nu residual fit failed")
end

function _residuals(model, heldout::Matrix{Float64})
    kind, fit = model
    if kind == :two
        responsibilities = SRE.student_t_responsibilities(heldout, fit.means, fit.scales)
        return heldout - responsibilities * fit.means
    end
    return heldout .- reshape(fit.mean, 1, :)
end

function _logsumexp(first::Float64, second::Float64)
    maximum_value = max(first, second)
    return maximum_value + log(exp(first - maximum_value) + exp(second - maximum_value))
end

function _score_rows(model, heldout::Matrix{Float64})
    kind, fit = model
    scores = Vector{Float64}(undef, size(heldout, 1))
    if kind == :two
        for row in axes(heldout, 1)
            first = log(0.5) + SRE.log_student_t_diag(
                view(heldout, row, :), view(fit.means, 1, :), view(fit.scales, 1, :))
            second = log(0.5) + SRE.log_student_t_diag(
                view(heldout, row, :), view(fit.means, 2, :), view(fit.scales, 2, :))
            scores[row] = _logsumexp(first, second) / size(heldout, 2)
        end
    else
        for row in axes(heldout, 1)
            scores[row] = SRE.log_student_t_diag(
                view(heldout, row, :), fit.mean, fit.scale) / size(heldout, 2)
        end
    end
    return scores
end

function _crossfit_partition(dataset::FeatureDataset, outer::String,
                             config::DiagnosticConfig)
    records = dataset.records
    record_dates = [dataset.date_by_scan[record.file] for record in records]
    residual_scans = ResidualScan[]
    view_difference = Dict{String,Float64}()
    dropout_difference = Dict{String,Float64}()
    union_matrix, _ = HUA.feature_matrix(records, config.feature_names)
    union_normalized = HUA.normalize_per_scan(records, union_matrix, config.feature_names)
    view_normalized = Matrix{Float64}[]
    for features in config.view_features
        matrix, _ = HUA.feature_matrix(records, features)
        push!(view_normalized, HUA.normalize_per_scan(records, matrix, features))
    end
    length(config.view_features) == 4 ||
        _fail(:config_contract, "dropout requires exactly the four configured views")
    dropout_features = copy(config.view_features[1])
    isempty(intersect(dropout_features, config.policy.channel_drop_features)) ||
        _fail(:config_contract, "dropout features overlap base_local")
    all(feature in config.feature_names for feature in config.policy.channel_drop_features) ||
        _fail(:config_contract, "dropout feature is absent from configured union")
    dropout_matrix, _ = HUA.feature_matrix(records, dropout_features)
    dropout_normalized = HUA.normalize_per_scan(records, dropout_matrix, dropout_features)
    inner_dates = sort!(unique(date for date in record_dates if date != outer))
    for inner in inner_dates
        training = findall(date -> date != outer && date != inner, record_dates)
        heldout = findall(==(inner), record_dates)
        isempty(training) && _fail(:underpowered, "inner training fold is empty")
        isempty(heldout) && _fail(:incomplete_fold, "inner held-out fold is empty")
        all(isfinite, union_normalized[training, :]) &&
            all(isfinite, union_normalized[heldout, :]) ||
            _fail(:nonfinite_input, "inner union view is invalid")
        model = _fit_residual_model(
            union_normalized[training, :], [records[index].amplitude for index in training])
        heldout_residuals = _residuals(model, union_normalized[heldout, :])
        heldout_scans = sort!(unique(records[index].file for index in heldout))
        for scan in heldout_scans
            local_positions = findall(position -> records[heldout[position]].file == scan,
                                      eachindex(heldout))
            push!(residual_scans, ResidualScan(
                inner, scan, Matrix{Float64}(heldout_residuals[local_positions, :])))
        end

        per_view_scan = [Dict{String,Float64}() for _ in config.view_names]
        for view_index in eachindex(config.view_names)
            normalized = view_normalized[view_index]
            all(isfinite, normalized[training, :]) && all(isfinite, normalized[heldout, :]) ||
                _fail(:nonfinite_input, "inner required view is invalid")
            view_model = _fit_residual_model(
                normalized[training, :], [records[index].amplitude for index in training])
            scores = _score_rows(view_model, normalized[heldout, :])
            for scan in heldout_scans
                positions = findall(position -> records[heldout[position]].file == scan,
                                    eachindex(heldout))
                per_view_scan[view_index][scan] = mean(scores[positions])
            end
        end
        all(isfinite, dropout_normalized[training, :]) &&
            all(isfinite, dropout_normalized[heldout, :]) ||
            _fail(:nonfinite_input, "dropout view is invalid")
        dropout_model = _fit_residual_model(
            dropout_normalized[training, :], [records[index].amplitude for index in training])
        dropout_scores = _score_rows(dropout_model, dropout_normalized[heldout, :])
        dropout_per_scan = Dict{String,Float64}()
        backward = findall(name -> occursin("bwd", lowercase(name)), config.view_names)
        for scan in heldout_scans
            view_scores = [per_view_scan[index][scan] for index in eachindex(per_view_scan)]
            difference = _view_contrast(view_scores, config)
            positions = findall(position -> records[heldout[position]].file == scan,
                                eachindex(heldout))
            dropout_per_scan[scan] = mean(dropout_scores[positions])
            view_difference[scan] = difference
            dropout_difference[scan] = mean(per_view_scan[index][scan] for index in backward) -
                                       dropout_per_scan[scan]
        end
    end
    return PartitionEvidence(outer, residual_scans, view_difference, dropout_difference)
end

function _run_hash(status::String, reason::String, aggregate::Vector{String},
                   partitions::Vector{PartitionResult}, nulls::NullCampaignSummary,
                   hashes::Vector{String})
    io = IOBuffer()
    println(io, "structured-diagnostic-run-v1")
    println(io, "status=", status)
    println(io, "reason=", reason)
    println(io, "aggregate=", join(aggregate, ','))
    println(io, "null=", nulls.sha256)
    for digest in hashes
        println(io, "input=", digest)
    end
    for partition in partitions
        println(io, "partition=", partition.outer_heldout_date, ':',
                partition.evidence_sha256)
    end
    return sha256_hex(take!(io))
end

function run_diagnostics(config_path::AbstractString, feature_path::AbstractString,
                         fold_path::AbstractString, receipt_path::AbstractString)
    config = load_diagnostic_config(config_path)
    dataset = load_feature_dataset(feature_path, config)
    folds = load_fold_plan(fold_path, dataset)
    receipt = load_input_receipt(receipt_path, config, dataset, folds)
    nulls = run_null_campaigns(config)
    source = model_source_hash()
    snapshots = [config.snapshot, dataset.snapshot, folds.snapshot, receipt.snapshot,
                 receipt.unary_review]
    foreach(_verify, snapshots)
    input_hashes = [config.sha256, dataset.snapshot.sha256, folds.snapshot.sha256,
                    receipt.snapshot.sha256, receipt.unary_review.sha256, source]
    if !nulls.accepted
        aggregate = fill("BLOCKED", length(DIAGNOSTIC_NAMES))
        digest = _run_hash("BLOCKED", "null_campaign_gate_failed", aggregate,
                           PartitionResult[], nulls, input_hashes)
        return DiagnosticRun(
            "BLOCKED", "null_campaign_gate_failed", aggregate, PartitionResult[], nulls,
            config.sha256, dataset.snapshot.sha256, folds.snapshot.sha256,
            receipt.snapshot.sha256, receipt.unary_review.sha256, source, digest, snapshots)
    end

    partitions = PartitionResult[]
    for outer in folds.dates
        evidence = _crossfit_partition(dataset, outer, config)
        training_records = [record for record in dataset.records
                            if dataset.date_by_scan[record.file] != outer]
        training_scans = unique(record.file for record in training_records)
        declaration = receipt.by_outer[outer]
        support = FeasibilityInput(
            [date for date in folds.dates if date != outer],
            length(training_scans),
            length(training_records),
            declaration.edges,
            declaration.categories,
            declaration.scale_eigenvalues,
        )
        push!(partitions, evaluate_partition(evidence, support, config))
    end
    aggregate = String[]
    for (index, _) in enumerate(DIAGNOSTIC_NAMES)
        statuses = [partition.mechanisms[index].status for partition in partitions]
        push!(aggregate, all(==("FOLLOWUP_REQUIRED"), statuses) ?
              "FOLLOWUP_REQUIRED" : "SKIPPED")
    end
    status = any(==("FOLLOWUP_REQUIRED"), aggregate) ? "FOLLOWUP_REQUIRED" : "SKIPPED"
    reason = status == "FOLLOWUP_REQUIRED" ? "outer_partitions_agree" :
             "followup_criterion_absent"
    foreach(_verify, snapshots)
    digest = _run_hash(status, reason, aggregate, partitions, nulls, input_hashes)
    return DiagnosticRun(
        status, reason, aggregate, partitions, nulls, config.sha256,
        dataset.snapshot.sha256, folds.snapshot.sha256, receipt.snapshot.sha256,
        receipt.unary_review.sha256, source, digest, snapshots)
end

_field(value::Integer) = string(value)
_field(value::Bool) = value ? "true" : "false"
_field(value::AbstractString) = String(value)
_field(value::Real) = isfinite(value) ? @sprintf("%.17g", Float64(value)) : "NA"

function _empty_report_row()
    return Dict(column => "NA" for column in REPORT_COLUMNS)
end

function _append_report_row!(rows::Vector{Dict{String,String}}, values::Pair...)
    row = _empty_report_row()
    for (key, value) in values
        key in REPORT_COLUMNS || _fail(:invalid_output, "unknown report field $key")
        text = _field(value)
        (occursin('\t', text) || occursin('\n', text)) &&
            _fail(:invalid_output, "report field contains a delimiter")
        row[key] = text
    end
    push!(rows, row)
    return nothing
end

function _report_rows(run::DiagnosticRun)
    rows = Dict{String,String}[]
    hashes = [
        "config_hash" => run.config_hash,
        "feature_hash" => run.feature_hash,
        "fold_hash" => run.fold_hash,
        "input_receipt_hash" => run.input_receipt_hash,
        "unary_review_hash" => run.unary_review_hash,
        "source_hash" => run.source_hash,
    ]
    _append_report_row!(rows,
        "record_type" => "receipt", "status" => run.status, "reason" => run.reason,
        "detail_hash" => run.run_sha256, hashes...)
    for (index, name) in enumerate(DIAGNOSTIC_NAMES)
        _append_report_row!(rows,
            "record_type" => "aggregate", "mechanism" => name,
            "status" => run.aggregate_statuses[index], "reason" => "all_outer_folds",
            "detail_hash" => run.run_sha256, hashes...)
    end
    _append_report_row!(rows,
        "record_type" => "null_summary",
        "status" => run.null_campaigns.accepted ? "SKIPPED" : "BLOCKED",
        "reason" => run.null_campaigns.accepted ? "null_gate_passed" : "null_gate_failed",
        "statistic" => run.null_campaigns.failures,
        "threshold" => run.null_campaigns.clopper_pearson_upper,
        "nodes" => run.null_campaigns.campaigns,
        "detail_hash" => run.null_campaigns.sha256, hashes...)
    for record in run.null_campaigns.records
        _append_report_row!(rows,
            "record_type" => "null_campaign", "seed" => record.seed,
            "status" => record.diagnostic_trigger ? "FOLLOWUP_REQUIRED" : "SKIPPED",
            "reason" => join(record.statuses, ','), "nodes" => record.nodes,
            "edges" => record.edges, "statistic" => record.latent_ones,
            "data_hash" => record.data_sha256, "edge_hash" => record.edge_sha256,
            "detail_hash" => record.result_sha256, hashes...)
        _append_report_row!(rows,
            "record_type" => "null_meta", "seed" => record.seed,
            "mechanism" => "complete_nested_meta", "status" => record.meta_status,
            "reason" => "complete_C1_C2_edge_graph_decision",
            "statistic" => record.meta_false_pass,
            "nodes" => record.nodes, "edges" => record.edges,
            "data_hash" => record.data_sha256, "edge_hash" => record.edge_sha256,
            "detail_hash" => record.meta_sha256, hashes...)
    end
    config = load_diagnostic_config()
    for partition in run.partitions
        _append_report_row!(rows,
            "record_type" => "partition", "outer_heldout_date" => partition.outer_heldout_date,
            "status" => partition.status, "reason" => "outer_training_only",
            "detail_hash" => partition.evidence_sha256, hashes...)
        for mechanism in partition.mechanisms
            _append_report_row!(rows,
                "record_type" => "mechanism", "outer_heldout_date" => partition.outer_heldout_date,
                "mechanism" => mechanism.name, "status" => mechanism.status,
                "reason" => mechanism.reason, "statistic" => mechanism.observed,
                "threshold" => mechanism.threshold, "lower" => mechanism.lower,
                "upper" => mechanism.upper, "raw_p" => mechanism.raw_p,
                "holm_p" => mechanism.holm_p, "detail_hash" => partition.evidence_sha256,
                hashes...)
            for date_result in mechanism.date_results
                _append_report_row!(rows,
                    "record_type" => "date_statistic",
                    "outer_heldout_date" => partition.outer_heldout_date,
                    "inner_date" => date_result.date, "mechanism" => mechanism.name,
                    "status" => mechanism.status, "reason" => "inner_date_only",
                    "statistic" => date_result.observed,
                    "threshold" => date_result.threshold, "lower" => date_result.excess,
                    "raw_p" => date_result.raw_p,
                    "detail_hash" => partition.evidence_sha256, hashes...)
                for (index, value) in enumerate(date_result.permutation_values)
                    _append_report_row!(rows,
                        "record_type" => "permutation",
                        "outer_heldout_date" => partition.outer_heldout_date,
                        "inner_date" => date_result.date,
                        "mechanism" => mechanism.name,
                        "status" => mechanism.status, "reason" => "within_scan_columns",
                        "seed" => config.bootstrap_seeds[index], "statistic" => value,
                        "detail_hash" => partition.evidence_sha256, hashes...)
                end
            end
            if mechanism.name == "scan_effects"
                for (index, value) in enumerate(mechanism.permutation_values)
                    _append_report_row!(rows,
                        "record_type" => "permutation",
                        "outer_heldout_date" => partition.outer_heldout_date,
                        "mechanism" => mechanism.name, "status" => mechanism.status,
                        "reason" => "within_date_scan_labels",
                        "seed" => config.bootstrap_seeds[index], "statistic" => value,
                        "detail_hash" => partition.evidence_sha256, hashes...)
                end
            end
            for (index, value) in enumerate(mechanism.bootstrap_values)
                _append_report_row!(rows,
                    "record_type" => "bootstrap",
                    "outer_heldout_date" => partition.outer_heldout_date,
                    "mechanism" => mechanism.name, "status" => mechanism.status,
                    "reason" => "whole_scan_within_date",
                    "seed" => config.bootstrap_seeds[index], "statistic" => value,
                    "detail_hash" => partition.evidence_sha256, hashes...)
            end
        end
        for feasibility in partition.feasibility
            _append_report_row!(rows,
                "record_type" => "feasibility",
                "outer_heldout_date" => partition.outer_heldout_date,
                "mechanism" => feasibility.model, "status" => feasibility.status,
                "reason" => feasibility.reason,
                "dimensions" => join(feasibility.dimensions, ','),
                "free_parameters" => feasibility.free_parameters,
                "inherited_parameters" => feasibility.inherited_parameters,
                "dates" => feasibility.dates, "scans" => feasibility.scans,
                "nodes" => feasibility.nodes, "edges" => feasibility.edges,
                "effective_observations" => feasibility.effective_observations,
                "parameter_observation_ratio" =>
                    feasibility.parameter_observation_ratio,
                "kish_ess" => feasibility.minimum_category_ess,
                "support_scans" => feasibility.minimum_support_scans,
                "scale_min_eigenvalue" => feasibility.minimum_scale_eigenvalue,
                "scale_condition_number" => feasibility.scale_condition_number,
                "detail_hash" => partition.evidence_sha256, hashes...)
            date_count = feasibility.dates
            for category_index in eachindex(feasibility.category_names)
                offset = (category_index - 1) * date_count
                support_slice = date_count == 0 ? Int[] :
                    feasibility.category_support_scans[(offset + 1):(offset + date_count)]
                _append_report_row!(rows,
                    "record_type" => "category_feasibility",
                    "outer_heldout_date" => partition.outer_heldout_date,
                    "mechanism" => feasibility.model,
                    "category" => feasibility.category_names[category_index],
                    "status" => feasibility.status, "reason" => feasibility.reason,
                    "kish_ess" => feasibility.category_ess[category_index],
                    "support_scans" => isempty(support_slice) ? 0 : minimum(support_slice),
                    "detail_hash" => partition.evidence_sha256, hashes...)
            end
        end
    end
    return rows
end

function report_tsv_bytes(run::DiagnosticRun)
    rows = _report_rows(run)
    io = IOBuffer()
    println(io, join(REPORT_COLUMNS, '\t'))
    for row in rows
        println(io, join((row[column] for column in REPORT_COLUMNS), '\t'))
    end
    return take!(io)
end

function blocked_report_bytes(error::Exception;
                              config_hash::String="NA",
                              feature_hash::String="NA",
                              fold_hash::String="NA",
                              input_receipt_hash::String="NA",
                              unary_review_hash::String="NA",
                              source_hash::String="NA")
    code = error isa DiagnosticError ? String(error.code) : "unexpected_error"
    error_hash = sha256_hex(collect(codeunits(sprint(showerror, error))))
    material = join((code, error_hash, config_hash, feature_hash, fold_hash, input_receipt_hash,
                     unary_review_hash, source_hash), '\n') * "\n"
    reason_hash = sha256_hex(collect(codeunits(material)))
    row = _empty_report_row()
    row["record_type"] = "receipt"
    row["status"] = "BLOCKED"
    row["reason"] = code
    row["config_hash"] = config_hash
    row["feature_hash"] = feature_hash
    row["fold_hash"] = fold_hash
    row["input_receipt_hash"] = input_receipt_hash
    row["unary_review_hash"] = unary_review_hash
    row["source_hash"] = source_hash
    row["detail_hash"] = reason_hash
    io = IOBuffer()
    println(io, join(REPORT_COLUMNS, '\t'))
    println(io, join((row[column] for column in REPORT_COLUMNS), '\t'))
    return take!(io)
end

function _rename_noreplace(source::String, destination::String)
    Sys.islinux() || _fail(:atomic_noreplace_unsupported,
                          "atomic no-replace publication requires Linux")
    result = ccall(:renameat2, Cint,
                   (Cint, Cstring, Cint, Cstring, Cuint),
                   Cint(-100), source, Cint(-100), destination, Cuint(1))
    result == 0 && return nothing
    error_number = Base.Libc.errno()
    error_number in (17, 21, 39) &&
        _fail(:output_collision, "output destination already has an owner")
    _fail(:atomic_publish_failure, "renameat2 failed with errno=$error_number")
end

function write_report(path::AbstractString, bytes::Vector{UInt8};
                      snapshots::Vector{FileSnapshot}=FileSnapshot[],
                      expected_source_hash::Union{Nothing,String}=nothing)
    supplied = String(path)
    isempty(supplied) && _fail(:invalid_path, "output path is empty")
    occursin('\0', supplied) && _fail(:invalid_path, "output path contains NUL")
    any(==(".."), split(replace(supplied, '\\' => '/'), '/'; keepempty=true)) &&
        _fail(:path_escape, "output parent traversal is forbidden")
    destination = normpath(isabspath(supplied) ? supplied : joinpath(ROOT, supplied))
    parent = dirname(destination)
    isdir(parent) || _fail(:missing_output_parent, "output parent is absent")
    _reject_symlink_components(parent)
    last(bytes) == UInt8('\n') || _fail(:invalid_output, "report lacks final LF")
    temporary, io = mktemp(parent; cleanup=false)
    committed = false
    try
        write(io, bytes)
        flush(io)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 ||
            _fail(:fsync_failed, "report fsync failed")
        close(io)
        foreach(_verify, snapshots)
        expected_source_hash === nothing || model_source_hash() == expected_source_hash ||
            _fail(:stale_source, "diagnostic source changed before publication")
        _rename_noreplace(temporary, destination)
        committed = true
    finally
        isopen(io) && close(io)
        !committed && ispath(temporary) && rm(temporary; force=true)
    end
    return destination
end

function model_source_hash()
    paths = [
        joinpath(@__DIR__, "diagnostics.jl"),
        joinpath(@__DIR__, "robust_emissions.jl"),
        joinpath(@__DIR__, "firewall.jl"),
        joinpath(@__DIR__, "..", "hierarchical_unit_assignment.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "emission_math.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "emissions.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "firewall.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "io_helpers.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "loading.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "nuisance.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "pipeline.jl"),
        joinpath(@__DIR__, "..", "hierarchical", "views.jl"),
        joinpath(ROOT, "test", "lib", "structured_unit_assignment.jl"),
        joinpath(ROOT, "test", "build_structured_unit_predictions.jl"),
        joinpath(ROOT, "test", "diagnose_structured_assignment_residuals.jl"),
    ]
    io = IOBuffer()
    println(io, "structured-diagnostics-source-v1")
    for path in sort(paths)
        normalized = normpath(path)
        snapshot = _snapshot(normalized)
        relative = replace(relpath(normalized, ROOT), '\\' => '/')
        expected = get(CONFIRMED_DEPENDENCY_HASHES, relative, nothing)
        expected === nothing || snapshot.sha256 == expected ||
            _fail(:dependency_hash, "confirmed dependency bytes changed")
        println(io, relative, '=', snapshot.sha256)
    end
    return sha256_hex(take!(io))
end

end
