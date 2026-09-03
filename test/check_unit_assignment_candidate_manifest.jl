#!/usr/bin/env julia

# Checker for config/unit_assignment_candidate.toml.
#
# Enforces the T0 contract from `.omo/plans/improve-unit-assignment-benchmark.md`:
#   1. FIREWALL — forbidden benchmark motif/encoding, truth/count column names,
#      and benchmark truth/grade paths may appear only inside [grader_only].
#   2. STRUCTURE — every T0 decision is present with the exact frozen value.
#   3. STATE — grade_status in {locked, frozen_once, graded}; provenance.status
#      and candidate.frozen_hash are consistent with grade_status.
#   4. HASH BINDING — when grade_status = "frozen_once" or "graded",
#      candidate.frozen_hash must equal SHA-256 of the exact non-grader source
#      bytes, excluding only [grader_only] plus the frozen_hash and grade_status
#      assignment text needed for self-reference and the frozen_once -> graded
#      lifecycle transition. Other comments and line endings remain bound.
#
# CLI:
#   julia --project=. test/check_unit_assignment_candidate_manifest.jl \
#       --config config/unit_assignment_candidate.toml [--expect-locked | --expect-frozen-once]
#
# Exit 0 on PASS with a parsed-status summary on stdout; exit nonzero on FAIL
# with diagnostics on stderr. PASS is decided from parsed values and exit code,
# not from log wording.

using SHA, TOML

const FORBIDDEN_TOKENS = ["NKNNKN", "010010", "sequence", "expected_N", "target_N"]
const FORBIDDEN_PATH_FRAGMENTS = [
    "benchmarks/",
    "results/benchmark_grades",
    "report_unit_assignment_benchmark",
    "grade_unit_assignment",
    "_unit_sequences.tsv",
]

struct Options
    config::String
    model::Union{Nothing,String}
    lifecycle_dir::Union{Nothing,String}
    expect_locked::Bool
    expect_frozen_once::Bool
    expect_grade_reserved::Bool
    expect_graded::Bool
    expect_unfrozen::Bool
end

function _fail(msg)
    println(stderr, "Unit-assignment candidate manifest: FAIL")
    println(stderr, msg)
    return 1
end

function _parse_cli(args)
    config = ""
    model = nothing
    lifecycle_dir = nothing
    expect_locked = false
    expect_frozen_once = false
    expect_grade_reserved = false
    expect_graded = false
    expect_unfrozen = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i < length(args) || error("--config requires a value")
            config = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config = split(arg, "="; limit=2)[2]; i += 1
        elseif arg in ("--model", "--model-config")
            i < length(args) || error("$arg requires a value")
            model === nothing || error("model config may be supplied only once")
            model = args[i+1]; i += 2
        elseif startswith(arg, "--model=") || startswith(arg, "--model-config=")
            model === nothing || error("model config may be supplied only once")
            model = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--lifecycle-dir"
            i < length(args) || error("--lifecycle-dir requires a value")
            lifecycle_dir === nothing || error("--lifecycle-dir may be supplied only once")
            lifecycle_dir = args[i+1]; i += 2
        elseif startswith(arg, "--lifecycle-dir=")
            lifecycle_dir === nothing || error("--lifecycle-dir may be supplied only once")
            lifecycle_dir = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--expect-locked"
            expect_locked = true; i += 1
        elseif arg == "--expect-frozen-once"
            expect_frozen_once = true; i += 1
        elseif arg == "--expect-grade-reserved"
            expect_grade_reserved = true; i += 1
        elseif arg == "--expect-graded"
            expect_graded = true; i += 1
        elseif arg == "--expect-unfrozen"
            expect_unfrozen = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/check_unit_assignment_candidate_manifest.jl --config PATH [--model PATH] [--lifecycle-dir PATH] [EXPECTATION]

            Validates the T0 firewall, structure, state, and hash-binding contract
            of config/unit_assignment_candidate.toml. PASS is decided from parsed
            values and exit code, not from log wording. Mutations after a frozen
            hash are rejected via the exact-source SHA-256 comparison.
            """)
            return nothing
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(config) && error("--config is required")
    return Options(config, model, lifecycle_dir, expect_locked, expect_frozen_once,
        expect_grade_reserved, expect_graded, expect_unfrozen)
end

# Walk a dotted path through nested Dicts; returns (value, ok).
function _get(parsed, path::AbstractString)
    parts = split(path, '.')
    cur = parsed
    for p in parts
        cur isa AbstractDict && haskey(cur, p) || return nothing, false
        cur = cur[p]
    end
    return cur, true
end

# Split a single-line TOML source line at its first comment marker outside a
# basic or literal string. Multiline strings are rejected separately, making
# table and assignment boundaries unambiguous without maintaining a TOML lexer.
function _split_comment(line::AbstractString)
    code = IOBuffer()
    comment = IOBuffer()
    in_basic = false
    in_literal = false
    escaped = false
    found_comment = false
    for c in line
        if found_comment
            write(comment, c)
        elseif in_basic
            write(code, c)
            if escaped
                escaped = false
            elseif c == '\\'
                escaped = true
            elseif c == '"'
                in_basic = false
            end
        elseif in_literal
            write(code, c)
            c == '\'' && (in_literal = false)
        elseif c == '#'
            found_comment = true
            write(comment, c)
        else
            write(code, c)
            c == '"' && (in_basic = true)
            c == '\'' && (in_literal = true)
        end
    end
    return String(take!(code)), String(take!(comment))
end

function _line_ending(line::String)
    endswith(line, "\r\n") && return "\r\n"
    endswith(line, "\n") && return "\n"
    endswith(line, "\r") && return "\r"
    return ""
end

function _assignment_key(code::String)
    m = match(r"^\s*(?:([A-Za-z0-9_-]+)|\"([^\"]+)\"|'([^']+)')\s*=", code)
    m === nothing && return nothing
    return something(m.captures[1], m.captures[2], m.captures[3])
end

# A deliberately conservative source contract. TOML multiline strings and
# non-simple table headers are rejected so text that resembles a table or
# assignment can never alter the firewall/hash boundary.
function _source_contract(text::String)
    section = ""
    grader_tables = 0
    frozen_hash_assignments = 0
    violations = String[]
    hash_buf = IOBuffer()

    for (lineno, raw) in enumerate(eachline(IOBuffer(text); keep=true))
        ending = _line_ending(raw)
        body = isempty(ending) ? raw : chop(raw; tail=length(ending))
        code, comment = _split_comment(body)

        if occursin("\"\"\"", code) || occursin("'''", code)
            push!(violations, "line $lineno: TOML multiline strings are not permitted in this manifest")
        end

        stripped = strip(code)
        header = match(r"^\[\s*([A-Za-z0-9_.-]+)\s*\]$", stripped)
        if startswith(stripped, '[') && endswith(stripped, ']') && header === nothing
            push!(violations, "line $lineno: unsupported table-header syntax; use simple [name] headers")
            section = "__invalid_header__"
        elseif header !== nothing
            section = header.captures[1]
            if section == "grader_only"
                grader_tables += 1
            elseif startswith(section, "grader_only.")
                push!(violations, "line $lineno: nested [grader_only.*] tables are forbidden")
            end
        end

        in_grader = section == "grader_only"
        if !in_grader
            lowered = lowercase(string(code, comment))
            for tok in FORBIDDEN_TOKENS
                occursin(lowercase(tok), lowered) && push!(violations, "line $lineno: forbidden token '$tok' outside [grader_only]: $(strip(body))")
            end
            for frag in FORBIDDEN_PATH_FRAGMENTS
                occursin(lowercase(frag), lowered) && push!(violations, "line $lineno: forbidden path fragment '$frag' outside [grader_only]: $(strip(body))")
            end
        end

        key = _assignment_key(code)
        if section == "candidate" && key == "frozen_hash"
            frozen_hash_assignments += 1
        end

        in_grader && continue
        if section == "candidate" && key in ("frozen_hash", "grade_status")
            # Exclude only lifecycle assignment text. Retain an inline comment
            # and the exact line ending so all other source bytes stay bound.
            write(hash_buf, comment)
            write(hash_buf, ending)
        else
            write(hash_buf, raw)
        end
    end

    grader_tables == 1 || push!(violations, "manifest must contain exactly one real [grader_only] table; found $grader_tables")
    return violations, frozen_hash_assignments, bytes2hex(sha256(take!(hash_buf)))
end

function _semantic_string_violations(value::AbstractString, path::String, kind::String)
    lowered = lowercase(value)
    violations = String[]
    for tok in FORBIDDEN_TOKENS
        occursin(lowercase(tok), lowered) &&
            push!(violations, "parsed $kind at $path contains forbidden token '$tok'")
    end
    for frag in FORBIDDEN_PATH_FRAGMENTS
        occursin(lowercase(frag), lowered) &&
            push!(violations, "parsed $kind at $path contains forbidden path fragment '$frag'")
    end
    return violations
end

# TOML.parse has already decoded \uXXXX and \UXXXXXXXX escapes. Walk every key
# and string value outside the sole top-level [grader_only] table so semantic
# tokens cannot evade the exact-source firewall by changing their spelling.
function _semantic_firewall(parsed::AbstractDict)
    violations = String[]

    function walk(value, path::String; root::Bool=false)
        if value isa AbstractDict
            for (key, child) in pairs(value)
                key_text = string(key)
                root && key_text == "grader_only" && continue
                key_path = isempty(path) ? key_text : "$path.$key_text"
                append!(violations, _semantic_string_violations(key_text, key_path, "key"))
                walk(child, key_path)
            end
        elseif value isa AbstractVector
            for (index, child) in pairs(value)
                walk(child, "$path[$index]")
            end
        elseif value isa AbstractString
            append!(violations, _semantic_string_violations(value, path, "string value"))
        end
    end

    walk(parsed, ""; root=true)
    return violations
end

const EXACT_CHECKS = [
    ("candidate.name", "hierarchical_equalprior_constant_current"),
    ("candidate.plan", ".omo/plans/improve-unit-assignment-benchmark.md"),
    ("candidate.version", 1),
    ("denominator.files", 145),
    ("denominator.control_positions_per_file", 6),
    ("denominator.control_positions_total", 870),
    ("features.policy", "equal_priors"),
    ("views.policy", "equal_weights"),
    ("bootstrap.count", 500),
    ("bootstrap.seeds_range_start", 0),
    ("bootstrap.seeds_range_stop", 499),
    ("bootstrap.resample_unit", "whole_scan_within_training_dates"),
    ("date_parser.rule", "leading_yyyymmdd_token"),
    ("date_parser.ambiguous", "fail"),
    ("date_parser.missing", "fail"),
    ("constant_current.mean_height_policy_nm", 0.50),
    ("constant_current.sensitivity_bracket_nm", [0.40, 0.45, 0.50, 0.55, 0.60]),
    ("constant_current.bracket_policy", "diagnostics_only"),
    ("real_gates.registration_interior_to_bounds", true),
    ("real_gates.beat_shifted_and_rotated_controls", true),
    ("real_gates.forward_backward_agreement", true),
    ("real_gates.fixed_perturbation_agreement", true),
    ("real_gates.abstain_on_invalid_isosurface", true),
    ("real_gates.abstain_on_registration_boundary", true),
    ("real_gates.do_not_widen_bounds_or_select_bracket_from_outputs", true),
    ("promotion.denominator_files", 145),
    ("promotion.denominator_control_positions", 870),
    ("promotion.honest_correct_min", 677),
    ("promotion.physical_accuracy_classified_min_fraction", 0.789),
    ("promotion.exact_chains_min", 18),
]

const CONTAIN_CHECKS = [
    ("features.hierarchical_emission_policy", "two-component"),
    ("features.hierarchical_emission_policy", "0.5"),
    ("features.hierarchical_emission_policy", "no occupancy regularizer"),
    ("one_vs_two_gate.pass_condition", "positive"),
    ("one_vs_two_gate.pass_condition", "95%"),
]

const REQUIRED_KEYS = String[
    "candidate", "denominator", "firewall", "provenance", "features", "views",
    "bootstrap", "date_parser", "constant_current", "one_vs_two_gate",
    "real_gates", "ranking", "promotion", "grader_only",
]
const MANDATORY_PROVENANCE_PATHS = String[
    "provenance.feature_tsves",
    "provenance.constant_current_cubes",
    "provenance.generated_maps",
    "provenance.molds",
    "provenance.configs",
]

_is_sha256_hex(s) = typeof(s) == String && length(s) == 64 && all(c -> c in "0123456789abcdef", s)

# Exact source for the frozen hash. Source-contract validation has already made
# section boundaries unambiguous. Only grader-only bytes and lifecycle
# assignment text are excluded; comments and line endings remain bound.
function _exact_source_hash(text::String)
    _, _, digest = _source_contract(text)
    return digest
end

function _require_ranking(parsed)
    order, ok = _get(parsed, "ranking.order")
    ok && order isa AbstractVector && !isempty(order) ||
        return ["ranking.order must be a non-empty list"]
    return String[]
end

function _check_view_weights(parsed)
    # Equal weights: views.policy must be equal_weights AND views.weight is a
    # finite number in (0,1]. The actual per-view weight is pinned at T0; T4 may
    # raise the count only by adding views that share this same weight.
    errs = String[]
    pol, ok = _get(parsed, "views.policy")
    (ok && pol == "equal_weights") || push!(errs, "views.policy must be \"equal_weights\", got $(ok ? pol : "<missing>")")
    w, ok2 = _get(parsed, "views.weight")
    (ok2 && w isa Number && isfinite(w) && w > 0.0 && w <= 1.0) || push!(errs, "views.weight must be a finite number in (0,1], got $(ok2 ? w : "<missing>")")
    return errs
end

function _check_locked(parsed, status, frozen_hash_assignments::Int)
    errs = String[]
    prov, ok = _get(parsed, "provenance.status")
    (ok && prov == "pending") || push!(errs, "locked candidate must have provenance.status = \"pending\", got $(ok ? prov : "<missing>")")
    fh, has_fh = _get(parsed, "candidate.frozen_hash")
    !has_fh || push!(errs, "locked candidate must not yet declare candidate.frozen_hash")
    frozen_hash_assignments == 0 || push!(errs, "locked candidate must contain no real candidate.frozen_hash assignment")
    return errs
end

function _check_frozen_provenance(parsed, status::String)
    errs = String[]
    for path in MANDATORY_PROVENANCE_PATHS
        value, ok = _get(parsed, path)
        if !ok
            push!(errs, "$status candidate is missing $path; every mandatory provenance field must be a lowercase 64-hex SHA-256")
        elseif !_is_sha256_hex(value)
            push!(errs, "$path must be a lowercase 64-hex SHA-256 for a $status candidate, got $(repr(value))")
        end
    end
    return errs
end

function _check_frozen_binding(parsed, text::String, status::String, frozen_hash_assignments::Int)
    errs = String[]
    prov, ok = _get(parsed, "provenance.status")
    (ok && prov == "frozen") || push!(errs, "$status candidate must have provenance.status = \"frozen\", got $(ok ? prov : "<missing>")")
    append!(errs, _check_frozen_provenance(parsed, status))
    fh, has_fh = _get(parsed, "candidate.frozen_hash")
    has_fh || (push!(errs, "$status candidate must declare candidate.frozen_hash"); return errs)
    frozen_hash_assignments == 1 || push!(errs, "$status candidate must contain exactly one real candidate.frozen_hash assignment; found $frozen_hash_assignments")
    _is_sha256_hex(fh) || push!(errs, "candidate.frozen_hash must be a 64-character lowercase-hex SHA-256, got: $fh")
    recomputed = _exact_source_hash(text)
    fh == recomputed || push!(errs, "Frozen hash mismatch: stored=$fh recomputed=$recomputed. Mutation after the freeze is rejected.")
    return errs
end

const V2_MODEL_SHA256 = "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const V2_CANDIDATE_SHA256 = "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const V2_FORBIDDEN_TERMS = [
    "benchmark", "truth", "motif", "labels", "label_path", "label_column",
    "expected_n", "expected n", "target_n", "target n", "class_count",
    "class count", "composition", "top_k", "top-k", "top k", "position",
    "parity", "run_length", "run-length", "run length", "transition", "switch",
    "learned_class_weight", "learned_view_weight", "class weights", "rank",
]

const V2_CANDIDATE_SCHEMA = Dict(
    "" => ["candidate", "provenance", "freeze", "grader_only"],
    "candidate" => ["name", "plan", "version", "model_config", "lifecycle_schema"],
    "provenance" => [
        "v1_candidate", "v1_candidate_sha256", "model_config_sha256",
        "checker_source", "spacing_source", "spacing_source_sha256", "forward_patch_producer",
        "forward_patch_producer_sha256", "backward_patch_producer",
        "backward_patch_producer_sha256",
    ],
    "freeze" => [
        "candidate_hash_scope", "model_hash_scope", "canonical_hash_exceptions",
        "receipt_chain", "receipt_files", "receipt_states",
        "prior_receipt_hash_required", "rewrite_prior_receipt", "mutable_toml_state",
    ],
    "grader_only" => [
        "denominator_manifest", "denominator_manifest_sha256", "grade_script",
        "grade_script_sha256", "profile", "denominator_files",
        "denominator_control_positions", "lobes_classified_min",
        "physical_accuracy_numerator", "physical_accuracy_denominator",
        "exact_chains_min", "honest_correct_min", "strict_improvement_required",
    ],
)

const V2_RECEIPT_SCHEMAS = Dict(
    "structured_freeze_receipt_v1" => [
        "schema", "schema_version", "state", "campaign_id", "attempt_id",
        "candidate_path", "candidate_sha256", "model_path", "model_sha256",
        "prediction_path", "prediction_sha256", "source_bundle_sha256",
        "universe_sha256", "evidence_index_sha256", "eligibility_receipt_sha256",
        "grader_source_sha256", "previous_receipt", "previous_receipt_sha256",
    ],
    "structured_grade_reservation_v1" => [
        "schema", "schema_version", "state", "campaign_id", "attempt_id",
        "candidate_sha256", "model_sha256", "prediction_path", "prediction_sha256",
        "repository_root", "project_path", "project_sha256", "grade_wrapper_path",
        "grade_wrapper_sha256", "report_script_path", "report_script_sha256",
        "exact_command", "working_directory", "report_destination", "environment_keys",
        "previous_receipt", "previous_receipt_sha256",
    ],
    "structured_grade_final_v1" => [
        "schema", "schema_version", "state", "outcome", "error_reason",
        "campaign_id", "attempt_id", "candidate_sha256", "model_sha256", "profile",
        "prediction_path", "prediction_sha256", "process_exit_code", "invocation_count",
        "summary_tsv_path", "summary_tsv_sha256", "report_md_path", "report_md_sha256",
        "lobe_position_errors_tsv_path", "lobe_position_errors_tsv_sha256",
        "grade_tsv_path", "grade_tsv_sha256", "stdout_path", "stdout_sha256",
        "stderr_path", "stderr_sha256", "previous_receipt", "previous_receipt_sha256",
    ],
)

const V2_RECEIPT_HASH_FIELDS = Dict(
    "structured_freeze_receipt_v1" => [
        "candidate_sha256", "model_sha256", "prediction_sha256",
        "source_bundle_sha256", "universe_sha256", "evidence_index_sha256",
        "eligibility_receipt_sha256", "grader_source_sha256",
    ],
    "structured_grade_reservation_v1" => [
        "candidate_sha256", "model_sha256", "prediction_sha256", "project_sha256",
        "grade_wrapper_sha256", "report_script_sha256", "previous_receipt_sha256",
    ],
    "structured_grade_final_v1" => [
        "candidate_sha256", "model_sha256", "prediction_sha256", "summary_tsv_sha256",
        "report_md_sha256", "lobe_position_errors_tsv_sha256", "grade_tsv_sha256",
        "stdout_sha256", "stderr_sha256", "previous_receipt_sha256",
    ],
)

const V2_RECEIPT_ABSOLUTE_PATH_FIELDS = Dict(
    "structured_freeze_receipt_v1" => [
        "candidate_path", "model_path", "prediction_path",
    ],
    "structured_grade_reservation_v1" => [
        "prediction_path", "repository_root", "project_path", "grade_wrapper_path",
        "report_script_path", "working_directory", "report_destination",
    ],
    "structured_grade_final_v1" => [
        "prediction_path", "summary_tsv_path", "report_md_path",
        "lobe_position_errors_tsv_path", "grade_tsv_path", "stdout_path", "stderr_path",
    ],
)

const V2_RECEIPT_IDENTITY_FIELDS = (
    "campaign_id", "attempt_id", "prediction_path", "prediction_sha256",
)
const V2_PRELAUNCH_ERROR_REASON = "crash_after_reservation_before_launch"

const V2_MODEL_SCHEMA = Dict(
    "" => ["model", "selection", "preprocessing"],
    "model" => [
        "schema", "version", "model_ids", "reference_model_id",
        "challenger_model_id", "edge_admission_model_id", "admitted_graph_model_id",
        "reference_family", "challenger_family", "node_priors", "mixture_weights",
        "mixture_weight_policy", "covariance", "covariance_floor", "student_t_nu",
        "n_starts", "unary_start_quantile_deltas", "max_iter", "tol",
        "objective_decrease_tolerance", "physical_orientation_feature", "views",
        "edge", "graph",
    ],
    "model.views" => [
        "names", "base_local", "backward_descriptors", "split_descriptor", "fusion",
        "fixed_equal_weights", "missing_view_policy", "all_views_missing_posterior",
        "all_views_missing_score_increment",
    ],
    "model.edge" => [
        "basis", "dimension", "category_order", "student_t_nu",
        "conditional_start_alphas", "mixed_state_mean", "scale", "scale_shrink_target",
        "scale_shrinkage_formula", "eigenvalue_floor", "condition_number_cap",
        "edge_weight", "pair_prior", "pair_prior_policy", "responsibility_policy",
        "pair_weight_formula", "kish_ess_formula", "fixed_weight_objective",
        "heldout_gain_formula", "max_iter", "tol", "objective_decrease_tolerance",
        "start_tie_tolerance", "start_tie_policy",
    ],
    "model.graph" => [
        "inference", "enable_condition", "eligible_edge_factor", "split_gap_policy",
        "isolated_node_policy", "reversal_equivariance_required",
    ],
    "selection" => [
        "schema", "version", "meta_algorithm", "common_reference", "unary_fallback",
        "edge_fallback", "terminal_outcomes", "heldout_partition", "scan_aggregation",
        "date_aggregation", "bootstrap", "scoring", "unary_gate", "edge_gate",
        "primary_test", "noninferiority", "feasibility", "diagnostics",
    ],
    "selection.bootstrap" => [
        "count", "seed_start", "seed_stop", "seed_family", "resample_unit",
        "edge_shuffle_count", "edge_shuffle_seed_start", "edge_shuffle_seed_stop",
        "edge_shuffle_seed_family", "sign_flip",
    ],
    "selection.scoring" => [
        "observation_count_formula", "scan_score_formula", "node_only_edge_score",
        "score_every_frozen_node", "score_every_eligible_edge",
        "all_view_abstention_score_increment", "all_view_abstention_posterior",
        "numerical_failure",
    ],
    "selection.unary_gate" => [
        "comparison", "nested_inner_leave_date_out", "every_inner_date_mean_positive",
        "bootstrap_lower_quantile", "bootstrap_lower_bound_strictly_positive",
        "failure_fallback",
    ],
    "selection.edge_gate" => [
        "comparison", "nested_inner_leave_date_out", "every_inner_date_mean_positive",
        "bootstrap_lower_quantile", "bootstrap_lower_bound_strictly_positive",
        "shuffle_upper_quantile", "observed_above_shuffle_quantile",
        "reversal_equivariance_required", "failure_fallback",
    ],
    "selection.primary_test" => [
        "single_primary_test", "claim", "date_statistic", "observed_statistic",
        "sign_statistic", "p_value_rule", "ties_in_tail", "monte_carlo_plus_one",
        "alpha", "per_date_positive_guard",
    ],
    "selection.noninferiority" => [
        "loss_formula", "difference_formula", "resample_unit",
        "bootstrap_upper_quantile", "maximum_allowed_loss_difference",
        "seed_restart_reported_separately",
    ],
    "selection.feasibility" => [
        "minimum_dates", "minimum_scans", "minimum_category_kish_ess",
        "minimum_distinct_edges", "minimum_support_scans_per_training_date",
        "effective_observations_per_free_parameter", "edge_scale_eigenvalue_floor",
        "edge_scale_condition_number_cap", "insufficient_outcome",
    ],
    "selection.diagnostics" => [
        "follow_up_only", "mechanisms", "allowed_outcomes", "policy",
        "outer_training_partitions_only", "multiple_testing", "permutation_count",
        "permutation_upper_quantile", "view_interval_level", "null_campaign_count",
        "null_campaign_min_dates", "null_campaign_min_scans",
        "null_campaign_chain_length_min", "null_campaign_chain_length_max",
        "clopper_pearson_confidence", "clopper_pearson_upper_max",
    ],
    "selection.diagnostics.policy" => [
        "schema", "residual_model", "residual_formula", "residual_fallback",
        "residual_failure", "icc_estimator", "icc_date_centering",
        "icc_feature_pooling", "icc_unbalanced_group_size", "view_score",
        "view_contrast", "channel_drop_features", "channel_drop_refit",
        "channel_drop_statistic", "channel_drop_same_sign", "quantile_method",
        "permutation_tail", "view_tail", "threshold_equality", "holm_order",
        "holm_reject_equality",
    ],
    "preprocessing" => [
        "schema", "version", "normalization", "normalization_scope", "spacing_source",
        "spacing_min_nm", "spacing_max_nm", "gap_policy", "channel", "stride",
        "flatten", "smooth_radius_px", "patch", "patch_producers", "residualization",
    ],
    "preprocessing.patch" => [
        "basis", "half_t_nm", "half_u_nm", "step_nm", "pixels_per_axis", "pixels",
        "order", "coordinates_nm", "forward_prefix", "backward_prefix", "grid_sha256",
        "centering", "scaling", "l2_norm_floor", "correlation_roundoff_tolerance",
        "serialization_format", "reversal_index_formula",
    ],
    "preprocessing.patch_producers" => [
        "forward_source", "forward_source_sha256", "forward_command", "backward_source",
        "backward_source_sha256", "backward_command", "execution_site", "julia_threads",
        "producer_value_normalization", "producer_nonpositive_std_fallback",
        "producer_serialization_format", "interpolation", "out_of_bounds",
        "forward_residual", "backward_residual",
    ],
    "preprocessing.residualization" => [
        "features", "endpoint_order", "design", "intercept_penalized",
        "predictor_count", "ridge_formula", "coefficient_formula",
        "heldout_residual_formula", "fit_scope",
    ],
)

const V2_CANDIDATE_EXACT_CHECKS = [
    ("candidate.name", "structured_label_free_unit_assignment"),
    ("candidate.plan", ".omo/plans/structured-label-free-unit-assignment.md"),
    ("candidate.version", 2),
    ("candidate.model_config", "config/unit_assignment_structured_model.toml"),
    ("candidate.lifecycle_schema", "sidecar_receipts_v1"),
    ("provenance.v1_candidate", "config/unit_assignment_candidate.toml"),
    ("provenance.v1_candidate_sha256", "9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa"),
    ("provenance.model_config_sha256", V2_MODEL_SHA256),
    ("provenance.checker_source", "test/check_unit_assignment_candidate_manifest.jl"),
    ("provenance.spacing_source", "config/chitosan.toml"),
    ("provenance.forward_patch_producer", "test/extract_lobe_patches.jl"),
    ("provenance.forward_patch_producer_sha256", "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad"),
    ("provenance.backward_patch_producer", "test/extract_lobe_patches_bwd.jl"),
    ("provenance.backward_patch_producer_sha256", "4f866a1e27289a0cf010a6ed4b5e0b92454dc55f97807d33cc95edb9e87182ba"),
    ("freeze.candidate_hash_scope", "exact_file_bytes"),
    ("freeze.model_hash_scope", "exact_file_bytes"),
    ("freeze.canonical_hash_exceptions", Any[]),
    ("freeze.receipt_chain", ["absence", "freeze_receipt", "grade_reservation", "grade_final"]),
    ("freeze.receipt_files", ["freeze_receipt.toml", "grade_reservation.toml", "grade_final.toml"]),
    ("freeze.receipt_states", ["unfrozen", "frozen_once", "grade_reserved", "graded"]),
    ("freeze.prior_receipt_hash_required", true),
    ("freeze.rewrite_prior_receipt", false),
    ("freeze.mutable_toml_state", false),
    ("grader_only.denominator_manifest", "benchmarks/chitosan_6mer_counting_confirmed.toml"),
    ("grader_only.grade_script", "test/report_unit_assignment_benchmark.jl"),
    ("grader_only.profile", "structured_v2"),
    ("grader_only.denominator_files", 145),
    ("grader_only.denominator_control_positions", 870),
    ("grader_only.lobes_classified_min", 854),
    ("grader_only.physical_accuracy_numerator", 677),
    ("grader_only.physical_accuracy_denominator", 854),
    ("grader_only.exact_chains_min", 36),
    ("grader_only.honest_correct_min", 677),
    ("grader_only.strict_improvement_required", true),
]

const V2_MODEL_CORE_CHECKS = [
    ("model.version", 2),
    ("model.model_ids", ["C1", "C2", "edge_admission", "admitted_graph"]),
    ("model.node_priors", [0.5, 0.5]),
    ("model.mixture_weights", [0.5, 0.5]),
    ("model.student_t_nu", 8),
    ("model.covariance_floor", 1.0e-4),
    ("model.unary_start_quantile_deltas", [0.25, 0.40, 0.45, 0.475, 0.49]),
    ("model.max_iter", 200),
    ("model.tol", 1.0e-6),
    ("model.objective_decrease_tolerance", 1.0e-10),
    ("model.edge.basis", ["corr_fwd", "corr_bwd"]),
    ("model.edge.conditional_start_alphas", [0.0, 0.25, 0.5, 0.75, 1.0]),
    ("model.edge.edge_weight", 1.0),
    ("model.edge.pair_prior", [0.0, 0.0, 0.0, 0.0]),
    ("selection.bootstrap.count", 500),
    ("selection.bootstrap.edge_shuffle_count", 500),
    ("selection.primary_test.single_primary_test", true),
    ("selection.primary_test.monte_carlo_plus_one", false),
    ("selection.feasibility.minimum_distinct_edges", 90),
    ("selection.diagnostics.follow_up_only", true),
    ("selection.diagnostics.policy.schema", "structured_diagnostics_policy_v1"),
    ("selection.diagnostics.policy.residual_model", "fixed_nu8_diagonal_student_t_two"),
    ("selection.diagnostics.policy.residual_formula",
        "x_minus_posterior_weighted_component_mean"),
    ("selection.diagnostics.policy.residual_fallback",
        "matched_one_component_student_t_mean"),
    ("selection.diagnostics.policy.residual_failure", "BLOCKED"),
    ("selection.diagnostics.policy.icc_estimator",
        "ICC_1_1_oneway_random_unbalanced"),
    ("selection.diagnostics.policy.icc_date_centering",
        "row_weighted_per_date_per_feature"),
    ("selection.diagnostics.policy.icc_feature_pooling",
        "equal_feature_summed_squares"),
    ("selection.diagnostics.policy.icc_unbalanced_group_size",
        "n0=(N-sum(n_s^2)/N)/(K-1)"),
    ("selection.diagnostics.policy.view_score",
        "mean_lobes(log_predictive_density/dimension)"),
    ("selection.diagnostics.policy.view_contrast",
        "mean(bwd_com,bwd_diag)-mean(base,split)"),
    ("selection.diagnostics.policy.channel_drop_features",
        ["bwd_neg_com_t", "bwd_neg_diag45"]),
    ("selection.diagnostics.policy.channel_drop_refit",
        "fresh_inner_training_base_local"),
    ("selection.diagnostics.policy.channel_drop_statistic",
        "mean(full_bwd_scores)-refit_base_score"),
    ("selection.diagnostics.policy.channel_drop_same_sign",
        "strict_nonzero_each_inner_date"),
    ("selection.diagnostics.policy.quantile_method",
        "Hyndman_Fan_type_7_explicit"),
    ("selection.diagnostics.policy.permutation_tail",
        "upper_inclusive_plus_one"),
    ("selection.diagnostics.policy.view_tail",
        "two_sided_zero_inclusive_plus_one"),
    ("selection.diagnostics.policy.threshold_equality", "SKIPPED"),
    ("selection.diagnostics.policy.holm_order",
        ["mfa_q1", "scan_effects", "view_asymmetry"]),
    ("selection.diagnostics.policy.holm_reject_equality", true),
    ("preprocessing.spacing_source", "config/chitosan.toml"),
    ("preprocessing.spacing_min_nm", 0.35),
    ("preprocessing.spacing_max_nm", 0.75),
    ("preprocessing.patch.grid_sha256", "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"),
    ("preprocessing.patch.l2_norm_floor", 1.0e-12),
    ("preprocessing.patch.correlation_roundoff_tolerance", 1.0e-12),
    ("preprocessing.patch.serialization_format", "%.17g"),
    ("preprocessing.patch_producers.forward_source_sha256", "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad"),
    ("preprocessing.patch_producers.backward_source_sha256", "4f866a1e27289a0cf010a6ed4b5e0b92454dc55f97807d33cc95edb9e87182ba"),
    ("preprocessing.residualization.ridge_formula", "max(1e-12,1e-6*trace(X[:,2:end]'*X[:,2:end])/14)"),
]

function _v2_schema_errors(parsed::AbstractDict, schema::Dict{String,Vector{String}}, label::String)
    errs = String[]
    for path in sort!(collect(keys(schema)))
        table = if isempty(path)
            parsed
        else
            value, ok = _get(parsed, path)
            if !ok || !(value isa AbstractDict)
                push!(errs, "$label table $path is missing or not a table")
                continue
            end
            value
        end
        actual = Set(String(key) for key in keys(table))
        expected = Set(schema[path])
        actual == expected || push!(errs,
            "$label keys at $(isempty(path) ? "<root>" : path) differ; expected $(sort!(collect(expected))), got $(sort!(collect(actual)))")
    end
    return errs
end

function _v2_exact_errors(parsed::AbstractDict, checks, label::String)
    errs = String[]
    for (path, expected) in checks
        value, ok = _get(parsed, path)
        ok && isequal(value, expected) || push!(errs,
            "$label $path must equal $(repr(expected)), got $(ok ? repr(value) : "<missing>")")
    end
    return errs
end

function _v2_source_violations(text::String; allow_grader::Bool)
    section = ""
    grader_tables = 0
    violations = String[]
    for (lineno, raw) in enumerate(eachline(IOBuffer(text); keep=true))
        ending = _line_ending(raw)
        body = isempty(ending) ? raw : chop(raw; tail=length(ending))
        code, comment = _split_comment(body)
        if occursin("\"\"\"", code) || occursin("'''", code)
            push!(violations, "line $lineno: TOML multiline strings are not permitted")
        end
        stripped = strip(code)
        header = match(r"^\[\s*([A-Za-z0-9_.-]+)\s*\]$", stripped)
        if startswith(stripped, '[') && endswith(stripped, ']') && header === nothing
            push!(violations, "line $lineno: unsupported table-header syntax")
            section = "__invalid_header__"
        elseif header !== nothing
            section = header.captures[1]
            section == "grader_only" && (grader_tables += 1)
            startswith(section, "grader_only.") &&
                push!(violations, "line $lineno: nested grader-only tables are forbidden")
        end
        in_grader = allow_grader && section == "grader_only"
        if !in_grader
            lowered = lowercase(string(code, comment))
            for term in V2_FORBIDDEN_TERMS
                occursin(term, lowered) && push!(violations,
                    "line $lineno: forbidden term '$term' outside the external-only namespace")
            end
        end
    end
    if allow_grader
        grader_tables == 1 || push!(violations,
            "candidate must contain exactly one real [grader_only] table; found $grader_tables")
    elseif grader_tables != 0
        push!(violations, "model config must not contain a grader-only table")
    end
    return violations
end

function _v2_semantic_violations(parsed::AbstractDict; allow_grader::Bool)
    violations = String[]
    function walk(value, path::String; root::Bool=false)
        if value isa AbstractDict
            for (key, child) in pairs(value)
                key_text = string(key)
                root && allow_grader && key_text == "grader_only" && continue
                child_path = isempty(path) ? key_text : "$path.$key_text"
                lowered_key = lowercase(key_text)
                for term in V2_FORBIDDEN_TERMS
                    occursin(term, lowered_key) && push!(violations,
                        "parsed key at $child_path contains forbidden term '$term'")
                end
                walk(child, child_path)
            end
        elseif value isa AbstractVector
            for (index, child) in pairs(value)
                walk(child, "$path[$index]")
            end
        elseif value isa AbstractString
            lowered_value = lowercase(value)
            for term in V2_FORBIDDEN_TERMS
                occursin(term, lowered_value) && push!(violations,
                    "parsed string at $path contains forbidden term '$term'")
            end
        end
    end
    walk(parsed, ""; root=true)
    return violations
end

_sha256_file(path::String) = bytes2hex(sha256(read(path)))

function _v2_symlink_component(path::String)
    lexical_path = isabspath(path) ? path : string(pwd(), '/', path)
    components = splitpath(lexical_path)
    isempty(components) && return nothing
    cursor = first(components)
    islink(cursor) && return cursor
    ispath(cursor) || return nothing
    for component in @view components[2:end]
        isempty(component) && continue
        component == "." && continue
        if component == ".."
            cursor = dirname(cursor)
            continue
        end
        cursor = joinpath(cursor, component)
        islink(cursor) && return cursor
        ispath(cursor) || return nothing
    end
    return nothing
end

function _v2_regular_file(path::String, label::String)
    symlink_component = _v2_symlink_component(path)
    symlink_component === nothing ||
        return nothing, "$label path contains a symlink component: $symlink_component"
    isfile(path) || return nothing, "$label not found: $path"
    return path, nothing
end

function _v2_registered_path(relative::String)
    isabspath(relative) && return nothing
    any(part -> part == "..", split(replace(relative, '\\' => '/'), '/')) && return nothing
    root = realpath(joinpath(@__DIR__, ".."))
    candidate = normpath(joinpath(root, relative))
    relpath(candidate, root) == ".." && return nothing
    cursor = root
    for part in split(replace(relative, '\\' => '/'), '/')
        isempty(part) && continue
        cursor = joinpath(cursor, part)
        islink(cursor) && return nothing
    end
    isfile(candidate) || return nothing
    return candidate
end

function _v2_registered_hash_errors(parsed::AbstractDict)
    bindings = [
        ("provenance.v1_candidate", "provenance.v1_candidate_sha256"),
        ("provenance.spacing_source", "provenance.spacing_source_sha256"),
        ("provenance.forward_patch_producer", "provenance.forward_patch_producer_sha256"),
        ("provenance.backward_patch_producer", "provenance.backward_patch_producer_sha256"),
        ("grader_only.denominator_manifest", "grader_only.denominator_manifest_sha256"),
        ("grader_only.grade_script", "grader_only.grade_script_sha256"),
    ]
    errs = String[]
    for (path_field, hash_field) in bindings
        relative, path_ok = _get(parsed, path_field)
        digest, hash_ok = _get(parsed, hash_field)
        if !path_ok || !(relative isa String)
            push!(errs, "$path_field must be a registered relative source path")
            continue
        end
        source = _v2_registered_path(relative)
        source === nothing && (push!(errs, "$path_field is absent, escaped, symlinked, or unregistered"); continue)
        if !hash_ok || !_is_sha256_hex(digest)
            push!(errs, "$hash_field must be a lowercase SHA-256")
            continue
        end
        actual = _sha256_file(source)
        actual == digest || push!(errs, "$hash_field mismatch for $relative: declared=$digest actual=$actual")
    end
    return errs
end

function _v2_receipt(path::String, expected_schema::String, expected_state::String,
                     candidate_hash::String, model_hash::String,
                     previous::String, previous_hash::String)
    symlink_component = _v2_symlink_component(path)
    symlink_component === nothing ||
        return nothing, ["lifecycle receipt path contains a symlink component: $symlink_component"]
    isfile(path) || return nothing, ["missing lifecycle receipt: $path"]
    parsed = try
        TOML.parsefile(path)
    catch e
        return nothing, ["lifecycle receipt TOML parse error in $path: $(sprint(showerror, e))"]
    end
    required = Dict(
        "schema" => expected_schema,
        "state" => expected_state,
        "candidate_sha256" => candidate_hash,
        "model_sha256" => model_hash,
        "previous_receipt" => previous,
        "previous_receipt_sha256" => previous_hash,
    )
    errs = String[]
    allowed_keys = get(V2_RECEIPT_SCHEMAS, expected_schema, nothing)
    if allowed_keys === nothing
        push!(errs, "unsupported lifecycle receipt schema: $expected_schema")
        return parsed, errs
    end
    actual_keys = Set(String(key) for key in keys(parsed))
    expected_keys = Set(allowed_keys)
    actual_keys == expected_keys || push!(errs,
        "receipt $(basename(path)) keys differ; expected $(sort!(collect(expected_keys))), got $(sort!(collect(actual_keys)))")
    for (key, expected) in required
        value = get(parsed, key, nothing)
        isequal(value, expected) || push!(errs,
            "receipt $(basename(path)) field $key must equal $(repr(expected)), got $(repr(value))")
    end
    get(parsed, "schema_version", nothing) === 1 || push!(errs,
        "receipt $(basename(path)) schema_version must equal 1")
    for key in V2_RECEIPT_HASH_FIELDS[expected_schema]
        value = get(parsed, key, nothing)
        _is_sha256_hex(value) || push!(errs,
            "receipt $(basename(path)) field $key must be a lowercase SHA-256")
    end
    for key in V2_RECEIPT_ABSOLUTE_PATH_FIELDS[expected_schema]
        value = get(parsed, key, nothing)
        value isa String && isabspath(value) || push!(errs,
            "receipt $(basename(path)) field $key must be an absolute path")
    end
    for key in ("campaign_id", "attempt_id")
        value = get(parsed, key, nothing)
        value isa String && !isempty(value) || push!(errs,
            "receipt $(basename(path)) field $key must be a nonempty string")
    end
    if expected_schema == "structured_grade_reservation_v1"
        command = get(parsed, "exact_command", nothing)
        command isa String && !isempty(command) || push!(errs,
            "receipt $(basename(path)) exact_command must be a nonempty string")
        environment_keys = get(parsed, "environment_keys", nothing)
        allowed_environment = Set([
            "PATH", "HOME", "USER", "JULIA_DEPOT_PATH", "JULIA_LOAD_PATH",
            "LANG", "LC_ALL", "GKSwstype",
        ])
        if !(environment_keys isa Vector) ||
           !all(value -> value isa String && value in allowed_environment, environment_keys) ||
           length(environment_keys) != length(unique(environment_keys))
            push!(errs, "receipt $(basename(path)) environment_keys must be a unique allowed subset")
        end
    elseif expected_schema == "structured_grade_final_v1"
        outcome = get(parsed, "outcome", nothing)
        outcome in ("success", "error") || push!(errs,
            "receipt $(basename(path)) outcome must be success or error")
        error_reason = get(parsed, "error_reason", nothing)
        error_reason isa String || push!(errs,
            "receipt $(basename(path)) error_reason must be a string")
        process_exit_code = get(parsed, "process_exit_code", nothing)
        process_exit_code isa Integer && !(process_exit_code isa Bool) || push!(errs,
            "receipt $(basename(path)) process_exit_code must be an integer")
        invocation_count = get(parsed, "invocation_count", nothing)
        invocation_count isa Integer && !(invocation_count isa Bool) || push!(errs,
            "receipt $(basename(path)) invocation_count must be an integer")
        get(parsed, "profile", nothing) == "structured_v2" || push!(errs,
            "receipt $(basename(path)) profile must equal structured_v2")
        if outcome == "success"
            error_reason == "" || push!(errs,
                "receipt $(basename(path)) successful outcome requires an empty error_reason")
            process_exit_code === 0 || push!(errs,
                "receipt $(basename(path)) successful outcome requires process_exit_code = 0")
            invocation_count === 1 || push!(errs,
                "receipt $(basename(path)) successful outcome requires invocation_count = 1")
        elseif outcome == "error"
            error_reason isa String && !isempty(error_reason) || push!(errs,
                "receipt $(basename(path)) error outcome requires a nonempty error_reason")
            if invocation_count === 0
                error_reason == V2_PRELAUNCH_ERROR_REASON || push!(errs,
                    "receipt $(basename(path)) zero-invocation error must record the post-reservation prelaunch failure")
                process_exit_code === -1 || push!(errs,
                    "receipt $(basename(path)) zero-invocation prelaunch error requires process_exit_code = -1")
            elseif invocation_count === 1
                error_reason != V2_PRELAUNCH_ERROR_REASON || push!(errs,
                    "receipt $(basename(path)) launched error must not use the prelaunch-only error_reason")
                process_exit_code !== -1 || push!(errs,
                    "receipt $(basename(path)) launched error must not use process_exit_code = -1")
            elseif invocation_count !== 1
                push!(errs,
                    "receipt $(basename(path)) error outcome requires invocation_count = 0 for prelaunch failure or 1 after launch")
            end
        end
    end
    return parsed, errs
end

function _v2_receipt_continuity_errors(
    receipt::AbstractDict, frozen::AbstractDict, receipt_name::String,
)
    errs = String[]
    for key in V2_RECEIPT_IDENTITY_FIELDS
        expected = get(frozen, key, nothing)
        actual = get(receipt, key, nothing)
        isequal(actual, expected) || push!(errs,
            "receipt $receipt_name field $key must preserve frozen value $(repr(expected)), got $(repr(actual))")
    end
    return errs
end

function _v2_reservation_binding_errors(
    reservation::AbstractDict, frozen::AbstractDict,
)
    prediction_path = get(frozen, "prediction_path", nothing)
    prediction_path isa String && isabspath(prediction_path) ||
        return ["freeze receipt prediction_path must be absolute before grade reservation"]

    root = realpath(joinpath(@__DIR__, ".."))
    project_path = joinpath(root, "Project.toml")
    grade_wrapper_path = joinpath(root, "test", "grade_frozen_structured_candidate.jl")
    report_script_path = joinpath(root, "test", "report_unit_assignment_benchmark.jl")
    report_destination = joinpath(dirname(prediction_path), "grade")
    exact_command = "julia --project=$root $report_script_path --full145-own-n " *
        "--profile structured_v2=$prediction_path --outdir $report_destination"
    required = [
        "repository_root" => root,
        "project_path" => project_path,
        "grade_wrapper_path" => grade_wrapper_path,
        "report_script_path" => report_script_path,
        "working_directory" => root,
        "report_destination" => report_destination,
        "exact_command" => exact_command,
    ]
    errs = String[]
    for (key, expected) in required
        actual = get(reservation, key, nothing)
        isequal(actual, expected) || push!(errs,
            "receipt grade_reservation.toml field $key must equal canonical $(repr(expected)), got $(repr(actual))")
    end
    return errs
end

function _v2_lifecycle(opt::Options, candidate_hash::String, model_hash::String)
    opt.lifecycle_dir === nothing && return "unfrozen", String[]
    directory = something(opt.lifecycle_dir)
    symlink_component = _v2_symlink_component(directory)
    symlink_component === nothing || return "unfrozen",
        ["lifecycle directory path contains a symlink component: $symlink_component"]
    isdir(directory) || return "unfrozen", ["lifecycle directory not found: $directory"]
    allowed = Set(["freeze_receipt.toml", "grade_reservation.toml", "grade_final.toml"])
    unknown = sort!([name for name in readdir(directory) if endswith(name, ".toml") && !(name in allowed)])
    isempty(unknown) || return "unfrozen", ["unregistered lifecycle TOML files: $(join(unknown, ", "))"]
    freeze_path = joinpath(directory, "freeze_receipt.toml")
    reservation_path = joinpath(directory, "grade_reservation.toml")
    final_path = joinpath(directory, "grade_final.toml")
    has_freeze = islink(freeze_path) || ispath(freeze_path)
    has_reservation = islink(reservation_path) || ispath(reservation_path)
    has_final = islink(final_path) || ispath(final_path)
    !has_freeze && (has_reservation || has_final) &&
        return "unfrozen", ["lifecycle receipts are out of order: freeze_receipt.toml is absent"]
    has_final && !has_reservation &&
        return "frozen_once", ["lifecycle receipts are out of order: grade_reservation.toml is absent"]
    !has_freeze && return "unfrozen", String[]
    freeze, errs = _v2_receipt(freeze_path, "structured_freeze_receipt_v1", "frozen_once",
        candidate_hash, model_hash, "absence", "absence")
    isempty(errs) || return "frozen_once", errs
    !has_reservation && return "frozen_once", String[]
    freeze_hash = _sha256_file(freeze_path)
    reservation, errs = _v2_receipt(reservation_path, "structured_grade_reservation_v1", "grade_reserved",
        candidate_hash, model_hash, "freeze_receipt.toml", freeze_hash)
    if isempty(errs)
        append!(errs, _v2_receipt_continuity_errors(
            something(reservation), something(freeze), "grade_reservation.toml"))
        append!(errs, _v2_reservation_binding_errors(
            something(reservation), something(freeze)))
    end
    isempty(errs) || return "grade_reserved", errs
    !has_final && return "grade_reserved", String[]
    reservation_hash = _sha256_file(reservation_path)
    final, errs = _v2_receipt(final_path, "structured_grade_final_v1", "graded",
        candidate_hash, model_hash, "grade_reservation.toml", reservation_hash)
    isempty(errs) && append!(errs, _v2_receipt_continuity_errors(
        something(final), something(freeze), "grade_final.toml"))
    isempty(errs) || return "graded", errs
    return "graded", String[]
end

function _run_v2(opt::Options, candidate_text::String, candidate::AbstractDict)
    opt.model === nothing && return _fail("Version 2 candidate requires --model PATH")
    _, config_error = _v2_regular_file(opt.config, "Candidate config")
    config_error === nothing || return _fail(config_error)
    model_path = something(opt.model)
    _, model_error = _v2_regular_file(model_path, "Model config")
    model_error === nothing || return _fail(model_error)
    model_text = read(model_path, String)
    model = try
        TOML.parse(model_text)
    catch e
        return _fail("TOML parse error in $model_path: $(sprint(showerror, e))")
    end

    violations = _v2_source_violations(candidate_text; allow_grader=true)
    append!(violations, _v2_semantic_violations(candidate; allow_grader=true))
    append!(violations, _v2_source_violations(model_text; allow_grader=false))
    append!(violations, _v2_semantic_violations(model; allow_grader=false))
    isempty(violations) || return _fail(join(vcat(["Source/firewall violations:"], violations), "\n  "))

    errs = String[]
    append!(errs, _v2_schema_errors(candidate, V2_CANDIDATE_SCHEMA, "candidate"))
    append!(errs, _v2_schema_errors(model, V2_MODEL_SCHEMA, "model"))
    append!(errs, _v2_exact_errors(candidate, V2_CANDIDATE_EXACT_CHECKS, "candidate"))
    append!(errs, _v2_exact_errors(model, V2_MODEL_CORE_CHECKS, "model"))
    candidate_hash = bytes2hex(sha256(codeunits(candidate_text)))
    candidate_hash == V2_CANDIDATE_SHA256 || push!(errs,
        "candidate config exact-byte SHA-256 must equal frozen $V2_CANDIDATE_SHA256, got $candidate_hash")
    model_hash = bytes2hex(sha256(codeunits(model_text)))
    model_hash == V2_MODEL_SHA256 || push!(errs,
        "model config exact-byte SHA-256 must equal frozen $V2_MODEL_SHA256, got $model_hash")
    declared_model_hash, has_model_hash = _get(candidate, "provenance.model_config_sha256")
    has_model_hash && declared_model_hash == model_hash || push!(errs,
        "candidate provenance.model_config_sha256 does not bind the supplied model bytes")
    append!(errs, _v2_registered_hash_errors(candidate))

    expectations = count(identity, [opt.expect_frozen_once, opt.expect_grade_reserved,
        opt.expect_graded, opt.expect_unfrozen])
    expectations <= 1 || push!(errs, "at most one v2 lifecycle expectation may be supplied")
    opt.expect_locked && push!(errs, "--expect-locked is valid only for version 1")
    lifecycle_state, lifecycle_errors = _v2_lifecycle(opt, candidate_hash, model_hash)
    append!(errs, lifecycle_errors)
    opt.expect_unfrozen && lifecycle_state != "unfrozen" &&
        push!(errs, "--expect-unfrozen requires lifecycle state unfrozen, got $lifecycle_state")
    opt.expect_frozen_once && lifecycle_state != "frozen_once" &&
        push!(errs, "--expect-frozen-once requires lifecycle state frozen_once, got $lifecycle_state")
    opt.expect_grade_reserved && lifecycle_state != "grade_reserved" &&
        push!(errs, "--expect-grade-reserved requires lifecycle state grade_reserved, got $lifecycle_state")
    opt.expect_graded && lifecycle_state != "graded" &&
        push!(errs, "--expect-graded requires lifecycle state graded, got $lifecycle_state")
    isempty(errs) || return _fail(join(vcat(["Contract violations:"], errs), "\n  "))

    println("Unit-assignment candidate manifest: OK")
    println("  config:            ", opt.config)
    println("  candidate_version: 2")
    println("  model:             ", model_path)
    println("  candidate_sha256:  ", candidate_hash)
    println("  model_sha256:      ", model_hash)
    println("  lifecycle_state: ", lifecycle_state)
    println("  firewall:          pass (external-only namespace isolated)")
    println("  status:            ok")
    return 0
end

function _run_v1(opt::Options)
    isfile(opt.config) || return _fail("Config not found: $(opt.config)")
    text = read(opt.config, String)
    parsed = try
        TOML.parse(text)
    catch e
        return _fail("TOML parse error in $(opt.config): $(sprint(showerror, e))")
    end

    # 1. Unambiguous source boundary and firewall
    violations, frozen_hash_assignments, _ = _source_contract(text)
    append!(violations, _semantic_firewall(parsed))
    isempty(violations) || return _fail(join(vcat(["Source/firewall violations:"], violations), "\n  "))

    # 2. Required top-level sections
    missing_sections = [k for k in REQUIRED_KEYS if !haskey(parsed, k)]
    isempty(missing_sections) || return _fail("Missing required sections: $(join(missing_sections, ", "))")

    # 3. Required exact and substring field checks
    errs = String[]
    for (path, expected) in EXACT_CHECKS
        v, ok = _get(parsed, path)
        ok && v == expected || push!(errs, "$path must equal $(repr(expected)), got $(ok ? repr(v) : "<missing>")")
    end
    for (path, needle) in CONTAIN_CHECKS
        v, ok = _get(parsed, path)
        ok && v isa AbstractString && occursin(needle, v) || push!(errs, "$path must contain $(repr(needle)), got $(ok ? repr(v) : "<missing>")")
    end
    append!(errs, _require_ranking(parsed))
    append!(errs, _check_view_weights(parsed))

    # 4. State consistency
    status, ok = _get(parsed, "candidate.grade_status")
    ok && status in ("locked", "frozen_once", "graded") || push!(errs, "candidate.grade_status must be one of locked|frozen_once|graded, got $(ok ? repr(status) : "<missing>")")
    if ok && status == "locked"
        append!(errs, _check_locked(parsed, status, frozen_hash_assignments))
    elseif ok && status in ("frozen_once", "graded")
        append!(errs, _check_frozen_binding(parsed, text, status, frozen_hash_assignments))
    end

    # 5. Expectation flags
    if opt.expect_locked
        (ok && status == "locked") || push!(errs, "--expect-locked requires grade_status = \"locked\", got $(ok ? repr(status) : "<missing>")")
    end
    if opt.expect_frozen_once
        (ok && status == "frozen_once") || push!(errs, "--expect-frozen-once requires grade_status = \"frozen_once\", got $(ok ? repr(status) : "<missing>")")
    end

    isempty(errs) || return _fail(join(vcat(["Contract violations:"], errs), "\n  "))

    # PASS
    prov, _ = _get(parsed, "provenance.status")
    println("Unit-assignment candidate manifest: OK")
    println("  config:        ", opt.config)
    println("  grade_status:  ", status)
    println("  provenance:    ", prov)
    println("  firewall:      pass (forbidden tokens confined to [grader_only])")
    if status in ("frozen_once", "graded")
        fh, _ = _get(parsed, "candidate.frozen_hash")
        println("  frozen_hash:   ", fh)
    end
    println("  status:        ok")
    return 0
end

function _dispatch_run(opt::Options)
    isfile(opt.config) || return _fail("Config not found: $(opt.config)")
    text = read(opt.config, String)
    parsed = try
        TOML.parse(text)
    catch e
        return _fail("TOML parse error in $(opt.config): $(sprint(showerror, e))")
    end
    version, ok = _get(parsed, "candidate.version")
    if !ok
        violations, _, _ = _source_contract(text)
        append!(violations, _semantic_firewall(parsed))
        isempty(violations) || return _fail(join(
            vcat(["Source/firewall violations:"], violations), "\n  "))
        return _fail("candidate.version is required for exact version dispatch")
    end
    if version === 1
        if opt.model !== nothing || opt.lifecycle_dir !== nothing || opt.expect_grade_reserved ||
           opt.expect_graded || opt.expect_unfrozen
            return _fail("Version 1 does not accept version-2 model, lifecycle, or expectation options")
        end
        return _run_v1(opt)
    elseif version === 2
        return _run_v2(opt, text, parsed)
    end
    return _fail("candidate.version must be exact integer 1 or 2, got $(repr(version))")
end

function main(args=ARGS)
    opt = try
        _parse_cli(args)
    catch e
        println(stderr, "Argument error: ", sprint(showerror, e))
        return 2
    end
    opt === nothing && return 0
    return _dispatch_run(opt)
end

exit(main())
