#!/usr/bin/env julia

using SHA, TOML, Test

const ROOT = dirname(@__DIR__)
const CHECKER = joinpath(@__DIR__, "check_unit_assignment_candidate_manifest.jl")
const V1_CONFIG = joinpath(ROOT, "config", "unit_assignment_candidate.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_structured_candidate.toml")
const MODEL = joinpath(ROOT, "config", "unit_assignment_structured_model.toml")
const V1_SHA256 = "9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa"
const V2_CANDIDATE_SHA256 = "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const CANONICAL_ROOT = realpath(ROOT)
const FIXTURE_PREDICTION = joinpath(
    CANONICAL_ROOT, "results", "structured_unit_assignment", "fixture", "predictions.tsv")
const FIXTURE_GRADE = joinpath(dirname(FIXTURE_PREDICTION), "grade")
const GRADE_WRAPPER = joinpath(CANONICAL_ROOT, "test", "grade_frozen_structured_candidate.jl")
const REPORT_SCRIPT = joinpath(CANONICAL_ROOT, "test", "report_unit_assignment_benchmark.jl")
const FIXTURE_GRADE_COMMAND =
    "julia --project=$CANONICAL_ROOT $REPORT_SCRIPT --full145-own-n " *
    "--profile structured_v2=$FIXTURE_PREDICTION --outdir $FIXTURE_GRADE"

sha256_hex(bytes) = bytes2hex(sha256(bytes))
file_sha256(path) = sha256_hex(read(path))

function replace_once_each(text::String, replacements...)
    updated = text
    for replacement in replacements
        source = first(replacement)
        length(findall(source, updated)) == 1 ||
            error("expected exactly one receipt field occurrence: $source")
        updated = replace(updated, replacement; count=1)
    end
    return updated
end

function run_capture(cmd::Cmd)
    out = IOBuffer()
    err = IOBuffer()
    process = run(pipeline(cmd; stdout=out, stderr=err), wait=false)
    wait(process)
    return process.exitcode, String(take!(out)), String(take!(err))
end

checker_cmd() = `$(Base.julia_cmd()) --project=$ROOT $CHECKER`

function with_model_hash(candidate_text::String, model_text::String)
    digest = sha256_hex(codeunits(model_text))
    return replace(
        candidate_text,
        r"(?m)^model_config_sha256\s*=\s*\"[0-9a-f]{64}\"$" =>
            "model_config_sha256 = \"$digest\"";
        count=1,
    )
end

function run_fixture(
    candidate_text::String,
    model_text::String;
    prepare::Function=(tmp, candidate_path, model_path) -> nothing,
    flags::Vector{String}=String[],
)
    mktempdir() do tmp
        candidate_path = joinpath(tmp, "candidate.toml")
        model_path = joinpath(tmp, "model.toml")
        write(candidate_path, candidate_text)
        write(model_path, model_text)
        prepare(tmp, candidate_path, model_path)
        command = Cmd(vcat(
            collect(checker_cmd()),
            ["--config", candidate_path, "--model", model_path],
            flags,
        ))
        code, out, err = run_capture(command)
        return code, out, err, sort(readdir(tmp))
    end
end

function run_model_mutation(mutation::Function)
    model_text = mutation(read(MODEL, String))
    candidate_text = with_model_hash(read(CANDIDATE, String), model_text)
    return run_fixture(candidate_text, model_text)
end

function receipt_text(schema, state, candidate_hash, model_hash, previous, previous_hash)
    if schema == "structured_freeze_receipt_v1"
        return """
        schema = "$schema"
        schema_version = 1
        state = "$state"
        campaign_id = "fixture-campaign"
        attempt_id = "fixture-attempt"
        candidate_path = "$(abspath(CANDIDATE))"
        candidate_sha256 = "$candidate_hash"
        model_path = "$(abspath(MODEL))"
        model_sha256 = "$model_hash"
        prediction_path = "$FIXTURE_PREDICTION"
        prediction_sha256 = "$(repeat('1', 64))"
        source_bundle_sha256 = "$(repeat('2', 64))"
        universe_sha256 = "$(repeat('3', 64))"
        evidence_index_sha256 = "$(repeat('4', 64))"
        eligibility_receipt_sha256 = "$(repeat('5', 64))"
        grader_source_sha256 = "$(repeat('6', 64))"
        previous_receipt = "$previous"
        previous_receipt_sha256 = "$previous_hash"
        """
    elseif schema == "structured_grade_reservation_v1"
        return """
        schema = "$schema"
        schema_version = 1
        state = "$state"
        campaign_id = "fixture-campaign"
        attempt_id = "fixture-attempt"
        candidate_sha256 = "$candidate_hash"
        model_sha256 = "$model_hash"
        prediction_path = "$FIXTURE_PREDICTION"
        prediction_sha256 = "$(repeat('1', 64))"
        repository_root = "$CANONICAL_ROOT"
        project_path = "$(joinpath(CANONICAL_ROOT, "Project.toml"))"
        project_sha256 = "$(repeat('7', 64))"
        grade_wrapper_path = "$GRADE_WRAPPER"
        grade_wrapper_sha256 = "$(repeat('8', 64))"
        report_script_path = "$REPORT_SCRIPT"
        report_script_sha256 = "$(repeat('9', 64))"
        exact_command = "$FIXTURE_GRADE_COMMAND"
        working_directory = "$CANONICAL_ROOT"
        report_destination = "$FIXTURE_GRADE"
        environment_keys = ["PATH", "HOME", "USER", "JULIA_DEPOT_PATH", "JULIA_LOAD_PATH", "LANG", "LC_ALL", "GKSwstype"]
        previous_receipt = "$previous"
        previous_receipt_sha256 = "$previous_hash"
        """
    elseif schema == "structured_grade_final_v1"
        return """
        schema = "$schema"
        schema_version = 1
        state = "$state"
        outcome = "success"
        error_reason = ""
        campaign_id = "fixture-campaign"
        attempt_id = "fixture-attempt"
        candidate_sha256 = "$candidate_hash"
        model_sha256 = "$model_hash"
        profile = "structured_v2"
        prediction_path = "$FIXTURE_PREDICTION"
        prediction_sha256 = "$(repeat('1', 64))"
        process_exit_code = 0
        invocation_count = 1
        summary_tsv_path = "$(joinpath(ROOT, "results", "structured_unit_assignment", "fixture", "grade", "summary.tsv"))"
        summary_tsv_sha256 = "$(repeat('a', 64))"
        report_md_path = "$(joinpath(ROOT, "results", "structured_unit_assignment", "fixture", "grade", "report.md"))"
        report_md_sha256 = "$(repeat('b', 64))"
        lobe_position_errors_tsv_path = "$(joinpath(ROOT, "results", "structured_unit_assignment", "fixture", "grade", "lobe_position_errors.tsv"))"
        lobe_position_errors_tsv_sha256 = "$(repeat('c', 64))"
        grade_tsv_path = "$(joinpath(ROOT, "results", "structured_unit_assignment", "fixture", "grade", "grade.tsv"))"
        grade_tsv_sha256 = "$(repeat('d', 64))"
        stdout_path = "$(joinpath(ROOT, "results", "structured_unit_assignment", "fixture", "grade", "stdout.log"))"
        stdout_sha256 = "$(repeat('e', 64))"
        stderr_path = "$(joinpath(ROOT, "results", "structured_unit_assignment", "fixture", "grade", "stderr.log"))"
        stderr_sha256 = "$(repeat('f', 64))"
        previous_receipt = "$previous"
        previous_receipt_sha256 = "$previous_hash"
        """
    end
    error("unsupported receipt schema: $schema")
end

function run_receipt_fixture(;
    reservation_mutation::Function=identity,
    final_mutation::Function=identity,
    expectation::String="--expect-graded",
    include_final::Bool=true,
)
    candidate_hash = file_sha256(CANDIDATE)
    model_hash = file_sha256(MODEL)
    return mktempdir() do tmp
        lifecycle = joinpath(tmp, "lifecycle")
        mkdir(lifecycle)
        freeze_path = joinpath(lifecycle, "freeze_receipt.toml")
        write(freeze_path, receipt_text(
            "structured_freeze_receipt_v1", "frozen_once", candidate_hash, model_hash,
            "absence", "absence",
        ))
        reservation_path = joinpath(lifecycle, "grade_reservation.toml")
        reservation = receipt_text(
            "structured_grade_reservation_v1", "grade_reserved", candidate_hash,
            model_hash, "freeze_receipt.toml", file_sha256(freeze_path),
        )
        write(reservation_path, reservation_mutation(reservation))
        if include_final
            final = receipt_text(
                "structured_grade_final_v1", "graded", candidate_hash, model_hash,
                "grade_reservation.toml", file_sha256(reservation_path),
            )
            write(joinpath(lifecycle, "grade_final.toml"), final_mutation(final))
        end
        return run_capture(`$(checker_cmd()) --config $CANDIDATE --model $MODEL --lifecycle-dir $lifecycle $expectation`)
    end
end

@testset "structured candidate/model artifacts" begin
    @test isfile(CHECKER)
    @test isfile(V1_CONFIG)
    @test file_sha256(V1_CONFIG) == V1_SHA256
    @test isfile(CANDIDATE)
    @test isfile(MODEL)

    candidate = TOML.parsefile(CANDIDATE)
    model = TOML.parsefile(MODEL)
    @test candidate["candidate"]["version"] === 2
    @test candidate["candidate"]["lifecycle_schema"] == "sidecar_receipts_v1"
    @test Set(keys(candidate)) == Set(["candidate", "provenance", "freeze", "grader_only"])
    @test Set(keys(model)) == Set(["model", "selection", "preprocessing"])
    @test file_sha256(CANDIDATE) == V2_CANDIDATE_SHA256
    @test candidate["provenance"]["model_config_sha256"] == file_sha256(MODEL)
    @test candidate["provenance"]["v1_candidate_sha256"] == V1_SHA256
    @test candidate["provenance"]["checker_source"] ==
          "test/check_unit_assignment_candidate_manifest.jl"
    @test candidate["provenance"]["spacing_source_sha256"] ==
          file_sha256(joinpath(ROOT, "config", "chitosan.toml"))
    @test candidate["freeze"]["canonical_hash_exceptions"] == []
    @test candidate["freeze"]["mutable_toml_state"] === false
    @test candidate["freeze"]["receipt_states"] ==
          ["unfrozen", "frozen_once", "grade_reserved", "graded"]
    grader = candidate["grader_only"]
    @test grader["lobes_classified_min"] === 854
    @test grader["physical_accuracy_numerator"] === 677
    @test grader["physical_accuracy_denominator"] === 854
    @test grader["exact_chains_min"] === 36
    @test grader["honest_correct_min"] === 677
    @test grader["strict_improvement_required"] === true
end

@testset "exact structured model contract" begin
    document = TOML.parsefile(MODEL)
    model = document["model"]
    selection = document["selection"]
    preprocessing = document["preprocessing"]

    @test model["model_ids"] == ["C1", "C2", "edge_admission", "admitted_graph"]
    @test model["reference_model_id"] == "C1"
    @test model["challenger_model_id"] == "C2"
    @test model["node_priors"] == [0.5, 0.5]
    @test model["mixture_weights"] == [0.5, 0.5]
    @test model["mixture_weight_policy"] == "fixed_not_estimated"
    @test model["student_t_nu"] === 8
    @test model["covariance_floor"] === 1.0e-4
    @test model["unary_start_quantile_deltas"] == [0.25, 0.40, 0.45, 0.475, 0.49]
    @test model["n_starts"] === 5
    @test model["max_iter"] === 200
    @test model["tol"] === 1.0e-6
    @test model["objective_decrease_tolerance"] === 1.0e-10

    views = model["views"]
    @test views["names"] == [
        "base_local",
        "base_local+bwd_neg_com_t",
        "base_local+bwd_neg_diag45",
        "base_local+split_log_skew",
    ]
    @test views["base_local"] == [
        "amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
    ]
    @test views["fixed_equal_weights"] == [0.25, 0.25, 0.25, 0.25]

    edge = model["edge"]
    @test edge["basis"] == ["corr_fwd", "corr_bwd"]
    @test edge["dimension"] === 2
    @test edge["student_t_nu"] === 8
    @test edge["conditional_start_alphas"] == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test edge["mixed_state_mean"] == "tied_01_10"
    @test edge["scale"] == "shared_full_2x2"
    @test edge["scale_shrink_target"] == "state_independent_null"
    @test edge["eigenvalue_floor"] === 1.0e-4
    @test edge["condition_number_cap"] === 1.0e4
    @test edge["edge_weight"] === 1.0
    @test edge["pair_prior"] == [0.0, 0.0, 0.0, 0.0]
    @test edge["pair_prior_policy"] == "zero_data_independent"
    @test occursin("w_ec", edge["fixed_weight_objective"])
    @test occursin("g_e", edge["heldout_gain_formula"])

    bootstrap = selection["bootstrap"]
    @test bootstrap["count"] === 500
    @test bootstrap["seed_start"] === 0
    @test bootstrap["seed_stop"] === 499
    @test bootstrap["edge_shuffle_count"] === 500
    @test bootstrap["sign_flip"] == "exhaustive_all_date_sign_assignments"
    @test selection["primary_test"]["single_primary_test"] === true
    @test selection["primary_test"]["p_value_rule"] ==
          "count(T_epsilon>=T_obs)/2^K"
    @test selection["primary_test"]["monte_carlo_plus_one"] === false
    @test selection["unary_gate"]["bootstrap_lower_quantile"] === 0.025
    @test selection["edge_gate"]["shuffle_upper_quantile"] === 0.975
    @test selection["noninferiority"]["bootstrap_upper_quantile"] === 0.975
    @test selection["noninferiority"]["maximum_allowed_loss_difference"] === 0.0

    feasibility = selection["feasibility"]
    @test feasibility["minimum_dates"] === 5
    @test feasibility["minimum_scans"] === 20
    @test feasibility["minimum_category_kish_ess"] === 10.0
    @test feasibility["minimum_distinct_edges"] === 90
    @test feasibility["minimum_support_scans_per_training_date"] === 3
    @test feasibility["effective_observations_per_free_parameter"] === 10.0
    @test selection["diagnostics"]["follow_up_only"] === true
    @test selection["diagnostics"]["mechanisms"] ==
          ["mfa_q1", "scan_effects", "view_asymmetry"]
    policy = selection["diagnostics"]["policy"]
    @test Set(keys(policy)) == Set([
        "schema", "residual_model", "residual_formula", "residual_fallback",
        "residual_failure", "icc_estimator", "icc_date_centering",
        "icc_feature_pooling", "icc_unbalanced_group_size", "view_score",
        "view_contrast", "channel_drop_features", "channel_drop_refit",
        "channel_drop_statistic", "channel_drop_same_sign", "quantile_method",
        "permutation_tail", "view_tail", "threshold_equality", "holm_order",
        "holm_reject_equality",
    ])
    @test policy["schema"] == "structured_diagnostics_policy_v1"
    @test policy["residual_model"] == "fixed_nu8_diagonal_student_t_two"
    @test policy["residual_formula"] ==
          "x_minus_posterior_weighted_component_mean"
    @test policy["residual_fallback"] == "matched_one_component_student_t_mean"
    @test policy["residual_failure"] == "BLOCKED"
    @test policy["icc_estimator"] == "ICC_1_1_oneway_random_unbalanced"
    @test policy["icc_date_centering"] == "row_weighted_per_date_per_feature"
    @test policy["icc_feature_pooling"] == "equal_feature_summed_squares"
    @test policy["icc_unbalanced_group_size"] ==
          "n0=(N-sum(n_s^2)/N)/(K-1)"
    @test policy["view_score"] == "mean_lobes(log_predictive_density/dimension)"
    @test policy["view_contrast"] == "mean(bwd_com,bwd_diag)-mean(base,split)"
    @test policy["channel_drop_features"] == ["bwd_neg_com_t", "bwd_neg_diag45"]
    @test policy["channel_drop_refit"] == "fresh_inner_training_base_local"
    @test policy["channel_drop_statistic"] ==
          "mean(full_bwd_scores)-refit_base_score"
    @test policy["channel_drop_same_sign"] == "strict_nonzero_each_inner_date"
    @test policy["quantile_method"] == "Hyndman_Fan_type_7_explicit"
    @test policy["permutation_tail"] == "upper_inclusive_plus_one"
    @test policy["view_tail"] == "two_sided_zero_inclusive_plus_one"
    @test policy["threshold_equality"] == "SKIPPED"
    @test policy["holm_order"] == ["mfa_q1", "scan_effects", "view_asymmetry"]
    @test policy["holm_reject_equality"] === true

    @test preprocessing["normalization"] == "per_scan_median_mad"
    @test preprocessing["spacing_source"] == "config/chitosan.toml"
    @test preprocessing["spacing_min_nm"] === 0.35
    @test preprocessing["spacing_max_nm"] === 0.75
    patch = preprocessing["patch"]
    @test patch["basis"] == ["corr_fwd", "corr_bwd"]
    @test patch["l2_norm_floor"] === 1.0e-12
    @test patch["correlation_roundoff_tolerance"] === 1.0e-12
    @test patch["serialization_format"] == "%.17g"
    @test patch["grid_sha256"] ==
          "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"
    @test preprocessing["residualization"]["ridge_formula"] ==
          "max(1e-12,1e-6*trace(X[:,2:end]'*X[:,2:end])/14)"
end

@testset "version dispatch preserves v1 and closes unknown versions" begin
    code1, out1, err1 = run_capture(`$(checker_cmd()) --config $V1_CONFIG --expect-locked`)
    @test code1 == 0
    @test isempty(err1)
    @test occursin("grade_status:  locked", out1)
    @test file_sha256(V1_CONFIG) == V1_SHA256

    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    for mutation in (
        text -> replace(text, r"(?m)^version\s*=\s*2\s*$" => ""; count=1),
        text -> replace(text, "version = 2" => "version = 3"; count=1),
        text -> replace(text, "version = 2" => "version = \"2\""; count=1),
        text -> replace(text, "version = 2" => "version = true"; count=1),
    )
        code, out, err, _ = run_fixture(mutation(candidate_text), model_text)
        @test code != 0
        @test isempty(out)
        @test occursin("version", lowercase(err))
    end
end

@testset "real v2 pair is read-only and accepted" begin
    before_candidate = read(CANDIDATE)
    before_model = read(MODEL)
    code, out, err = run_capture(`$(checker_cmd()) --config $CANDIDATE --model $MODEL`)
    @test code == 0
    @test isempty(err)
    @test occursin("candidate_version: 2", out)
    @test occursin("lifecycle_state: unfrozen", out)
    @test read(CANDIDATE) == before_candidate
    @test read(MODEL) == before_model
end

@testset "immutable sidecar lifecycle derives state through reservation" begin
    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    candidate_hash = sha256_hex(codeunits(candidate_text))
    model_hash = sha256_hex(codeunits(model_text))

    mktempdir() do tmp
        candidate_path = joinpath(tmp, "candidate.toml")
        model_path = joinpath(tmp, "model.toml")
        lifecycle = joinpath(tmp, "lifecycle")
        write(candidate_path, candidate_text)
        write(model_path, model_text)
        mkdir(lifecycle)
        frozen_candidate = read(candidate_path)
        frozen_model = read(model_path)

        cmd = `$(checker_cmd()) --config $candidate_path --model $model_path --lifecycle-dir $lifecycle`
        code0, out0, err0 = run_capture(cmd)
        @test code0 == 0
        @test isempty(err0)
        @test occursin("lifecycle_state: unfrozen", out0)

        freeze = receipt_text(
            "structured_freeze_receipt_v1", "frozen_once", candidate_hash, model_hash,
            "absence", "absence",
        )
        freeze_path = joinpath(lifecycle, "freeze_receipt.toml")
        write(freeze_path, freeze)
        code1, out1, err1 = run_capture(`$cmd --expect-frozen-once`)
        @test code1 == 0
        @test isempty(err1)
        @test occursin("lifecycle_state: frozen_once", out1)

        reservation = receipt_text(
            "structured_grade_reservation_v1", "grade_reserved", candidate_hash,
            model_hash, "freeze_receipt.toml", file_sha256(freeze_path),
        )
        reservation_path = joinpath(lifecycle, "grade_reservation.toml")
        write(reservation_path, reservation)
        code2, out2, err2 = run_capture(`$cmd --expect-grade-reserved`)
        @test code2 == 0
        @test isempty(err2)
        @test occursin("lifecycle_state: grade_reserved", out2)
        @test read(candidate_path) == frozen_candidate
        @test read(model_path) == frozen_model

        final = receipt_text(
            "structured_grade_final_v1", "graded", candidate_hash, model_hash,
            "grade_reservation.toml", file_sha256(reservation_path),
        )
        write(joinpath(lifecycle, "grade_final.toml"), final)
        code_final, out_final, err_final = run_capture(`$cmd --expect-graded`)
        @test code_final == 0
        @test isempty(err_final)
        @test occursin("lifecycle_state: graded", out_final)
        @test read(candidate_path) == frozen_candidate
        @test read(model_path) == frozen_model

        broken = replace(reservation, file_sha256(freeze_path) => repeat('0', 64); count=1)
        write(reservation_path, broken)
        code3, out3, err3 = run_capture(cmd)
        @test code3 != 0
        @test isempty(out3)
        @test occursin("previous", lowercase(err3)) || occursin("hash", lowercase(err3))
        @test read(candidate_path) == frozen_candidate
        @test read(model_path) == frozen_model
    end
end

@testset "exact model mutations fail after recalculated outer hash" begin
    mutations = Dict{String,Function}(
        "learned class weights" => text -> replace(
            text, "[model]" => "[model]\nlearned_class_weights = [0.5, 0.5]"; count=1),
        "nonzero pair prior" => text -> replace(
            text, "pair_prior = [0.0, 0.0, 0.0, 0.0]" =>
                "pair_prior = [0.0, 0.2, 0.2, 0.0]"; count=1),
        "changed nu" => text -> replace(text, "student_t_nu = 8" => "student_t_nu = 7"; count=1),
        "changed unary start" => text -> replace(
            text, "unary_start_quantile_deltas = [0.25, 0.40, 0.45, 0.475, 0.49]" =>
                "unary_start_quantile_deltas = [0.20, 0.40, 0.45, 0.475, 0.49]";
            count=1,
        ),
        "changed edge weight" => text -> replace(
            text, "edge_weight = 1.0" => "edge_weight = 0.5"; count=1),
        "changed spacing source" => text -> replace(
            text, "spacing_source = \"config/chitosan.toml\"" =>
                "spacing_source = \"config/template.toml\""; count=1),
        "changed grid hash" => text -> replace(
            text,
            "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e" =>
                repeat('0', 64);
            count=1,
        ),
        "changed producer hash" => text -> replace(
            text,
            "f79258197cd123f833e0541647c4e7107044149deb558ba765fbeee0545737ad" =>
                repeat('1', 64);
            count=1,
        ),
        "unknown model id" => text -> replace(
            text, "\"admitted_graph\"]" => "\"admitted_graph\", \"mfa_q1\"]"; count=1),
        "alternate rank" => text -> replace(text, "[model]" => "[model]\nrank = 1"; count=1),
        "unregistered source" => text -> replace(
            text, "forward_source = \"test/extract_lobe_patches.jl\"" =>
                "forward_source = \"test/unregistered_patch_source.jl\"";
            count=1,
        ),
    )
    for (case, mutation) in sort!(collect(mutations); by=first)
        code, out, err, entries = run_model_mutation(mutation)
        @test code != 0
        @test isempty(out)
        @test !isempty(err)
        @test entries == ["candidate.toml", "model.toml"]
    end
end

@testset "diagnostics policy schema and semantic mutations fail closed" begin
    mutations = [
        "missing policy table" => (text -> replace(
            text, "[selection.diagnostics.policy]" =>
                "[selection.diagnostics.missing_policy]"; count=1)),
        "extra policy key" => (text -> replace(
            text, "quantile_method = \"Hyndman_Fan_type_7_explicit\"" =>
                "quantile_method = \"Hyndman_Fan_type_7_explicit\"\npolicy_extra = true";
            count=1)),
        "missing policy key" => (text -> replace(
            text, "quantile_method = \"Hyndman_Fan_type_7_explicit\"\n" => "";
            count=1)),
        "wrong channel-drop refit semantic" => (text -> replace(
            text, "channel_drop_refit = \"fresh_inner_training_base_local\"" =>
                "channel_drop_refit = \"view_contrast_alias\""; count=1)),
        "wrong quantile/tie semantic" => (text -> replace(
            text, "view_tail = \"two_sided_zero_inclusive_plus_one\"" =>
                "view_tail = \"two_sided_strict\""; count=1)),
    ]
    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    for (case, mutation) in mutations
        @testset "$case" begin
            code, out, err, entries = run_model_mutation(mutation)
            @test code != 0
            @test isempty(out)
            @test !isempty(err)
            @test occursin("policy", lowercase(err))
            @test entries == ["candidate.toml", "model.toml"]
        end
    end
    @test candidate_text == read(CANDIDATE, String)
    @test model_text == read(MODEL, String)
end

@testset "candidate boundary, lifecycle state, and source hashes fail closed" begin
    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    mutations = (
        text -> replace(text, "[candidate]" => "[candidate]\ncurrent_state = \"eligible\""; count=1),
        text -> replace(text, "[provenance]" =>
            "[provenance]\nbenchmark_path = \"results/benchmark_grades/x.tsv\""; count=1),
        text -> replace(text, V1_SHA256 => repeat('0', 64); count=1),
        text -> replace(
            text,
            r"(?m)^spacing_source_sha256\s*=\s*\"[0-9a-f]{64}\"$" =>
                "spacing_source_sha256 = \"$(repeat('0', 64))\"";
            count=1,
        ),
        text -> replace(text, "model_config = \"config/unit_assignment_structured_model.toml\"" =>
            "model_config = \"../escaped.toml\""; count=1),
    )
    for mutation in mutations
        code, out, err, _ = run_fixture(mutation(candidate_text), model_text)
        @test code != 0
        @test isempty(out)
        @test !isempty(err)
    end

    # The real external-only namespace intentionally contains otherwise
    # forbidden vocabulary, yet the unmodified pair remains valid.
    @test occursin("benchmarks/", candidate_text)
    @test occursin("report_unit_assignment_benchmark", candidate_text)
    code, _, err, _ = run_fixture(candidate_text, model_text)
    @test code == 0
    @test isempty(err)
end

@testset "malformed, duplicate, partial, and symlink inputs fail closed" begin
    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    malformed = (
        candidate_text * "\nbroken = [\n",
        replace(candidate_text, "version = 2" => "version = 2\nversion = 2"; count=1),
        first(candidate_text, length(candidate_text) ÷ 2),
    )
    for candidate_fixture in malformed
        code, out, err, _ = run_fixture(candidate_fixture, model_text)
        @test code != 0
        @test isempty(out)
        @test !isempty(err)
    end

    partial_model = first(model_text, length(model_text) ÷ 2)
    code0, out0, err0, _ = run_fixture(
        with_model_hash(candidate_text, partial_model), partial_model)
    @test code0 != 0
    @test isempty(out0)
    @test !isempty(err0)

    mktempdir() do tmp
        candidate_path = joinpath(tmp, "candidate.toml")
        model_path = joinpath(tmp, "model.toml")
        candidate_link = joinpath(tmp, "candidate-link.toml")
        model_link = joinpath(tmp, "model-link.toml")
        write(candidate_path, candidate_text)
        write(model_path, model_text)
        symlink(candidate_path, candidate_link)
        symlink(model_path, model_link)
        code1, out1, err1 = run_capture(
            `$(checker_cmd()) --config $candidate_link --model $model_path`)
        @test code1 != 0
        @test isempty(out1)
        @test occursin("symlink", lowercase(err1))
        code2, out2, err2 = run_capture(
            `$(checker_cmd()) --config $candidate_path --model $model_link`)
        @test code2 != 0
        @test isempty(out2)
        @test occursin("symlink", lowercase(err2))
    end
end

@testset "correction regressions reject symlink ancestors and mutable sidecars" begin
    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    candidate_hash = sha256_hex(codeunits(candidate_text))
    model_hash = sha256_hex(codeunits(model_text))

    mktempdir() do tmp
        actual = joinpath(tmp, "actual")
        mkdir(actual)
        candidate_path = joinpath(actual, "candidate.toml")
        model_path = joinpath(actual, "model.toml")
        write(candidate_path, candidate_text)
        write(model_path, model_text)
        linked = joinpath(tmp, "linked")
        symlink(actual, linked)
        code, out, err = run_capture(`$(checker_cmd()) --config $(joinpath(linked, "candidate.toml")) --model $(joinpath(linked, "model.toml"))`)
        @test code != 0
        @test isempty(out)
        @test occursin("symlink", lowercase(err))
    end

    mktempdir() do tmp
        actual = joinpath(tmp, "actual")
        lifecycle = joinpath(actual, "lifecycle")
        mkpath(lifecycle)
        freeze_path = joinpath(lifecycle, "freeze_receipt.toml")
        write(freeze_path, receipt_text(
            "structured_freeze_receipt_v1", "frozen_once", candidate_hash, model_hash,
            "absence", "absence",
        ))
        linked = joinpath(tmp, "linked")
        symlink(actual, linked)
        code, out, err = run_capture(`$(checker_cmd()) --config $CANDIDATE --model $MODEL --lifecycle-dir $(joinpath(linked, "lifecycle")) --expect-frozen-once`)
        @test code != 0
        @test isempty(out)
        @test occursin("symlink", lowercase(err))
    end

    mktempdir() do tmp
        lifecycle = joinpath(tmp, "lifecycle")
        mkdir(lifecycle)
        symlink(joinpath(tmp, "absent-freeze.toml"), joinpath(lifecycle, "freeze_receipt.toml"))
        code, out, err = run_capture(`$(checker_cmd()) --config $CANDIDATE --model $MODEL --lifecycle-dir $lifecycle`)
        @test code != 0
        @test isempty(out)
        @test occursin("symlink", lowercase(err))
    end

    mktempdir() do tmp
        lifecycle = joinpath(tmp, "lifecycle")
        mkdir(lifecycle)
        freeze_path = joinpath(lifecycle, "freeze_receipt.toml")
        freeze = receipt_text(
            "structured_freeze_receipt_v1", "frozen_once", candidate_hash, model_hash,
            "absence", "absence",
        ) * "current_state = \"eligible\"\n"
        write(freeze_path, freeze)
        code, out, err = run_capture(`$(checker_cmd()) --config $CANDIDATE --model $MODEL --lifecycle-dir $lifecycle --expect-frozen-once`)
        @test code != 0
        @test isempty(out)
        @test occursin("current_state", err) || occursin("keys", lowercase(err))
    end

    mktempdir() do tmp
        lifecycle = joinpath(tmp, "lifecycle")
        mkdir(lifecycle)
        freeze_path = joinpath(lifecycle, "freeze_receipt.toml")
        write(freeze_path, receipt_text(
            "structured_freeze_receipt_v1", "frozen_once", candidate_hash, model_hash,
            "absence", "absence",
        ))
        reservation_path = joinpath(lifecycle, "grade_reservation.toml")
        write(reservation_path, receipt_text(
            "structured_grade_reservation_v1", "grade_reserved", candidate_hash,
            model_hash, "freeze_receipt.toml", file_sha256(freeze_path),
        ))
        final = receipt_text(
            "structured_grade_final_v1", "graded", candidate_hash, model_hash,
            "grade_reservation.toml", file_sha256(reservation_path),
        ) * "current_state = \"eligible\"\n"
        write(joinpath(lifecycle, "grade_final.toml"), final)
        code, out, err = run_capture(`$(checker_cmd()) --config $CANDIDATE --model $MODEL --lifecycle-dir $lifecycle --expect-graded`)
        @test code != 0
        @test isempty(out)
        @test occursin("current_state", err) || occursin("keys", lowercase(err))
    end
end

@testset "correction2 rejects lexical symlink components before dotdot" begin
    candidate_text = read(CANDIDATE, String)
    model_text = read(MODEL, String)
    mktempdir() do tmp
        visible = joinpath(tmp, "visible")
        target = joinpath(visible, "target")
        mkpath(target)
        write(joinpath(visible, "candidate.toml"), candidate_text)
        write(joinpath(visible, "model.toml"), model_text)
        link = joinpath(visible, "link")
        symlink(target, link)
        candidate_path = string(link, "/../candidate.toml")
        model_path = string(link, "/../model.toml")
        code, out, err = run_capture(
            `$(checker_cmd()) --config $candidate_path --model $model_path`)
        @test code != 0
        @test isempty(out)
        @test occursin("symlink", lowercase(err)) || occursin("lexical", lowercase(err))
    end
end

@testset "correction2 freezes receipt identity and prediction continuity" begin
    alternate_prediction = joinpath(dirname(FIXTURE_PREDICTION), "other-predictions.tsv")
    reservation_mutations = [
        "campaign_id" => (text -> replace(
            text, "campaign_id = \"fixture-campaign\"" =>
                "campaign_id = \"other-campaign\""; count=1)),
        "attempt_id" => (text -> replace(
            text, "attempt_id = \"fixture-attempt\"" =>
                "attempt_id = \"other-attempt\""; count=1)),
        "prediction_path" => (text -> replace(
            text, "prediction_path = \"$FIXTURE_PREDICTION\"" =>
                "prediction_path = \"$alternate_prediction\""; count=1)),
        "prediction_sha256" => (text -> replace(
            text, "prediction_sha256 = \"$(repeat('1', 64))\"" =>
                "prediction_sha256 = \"$(repeat('0', 64))\""; count=1)),
    ]
    for (field, mutation) in reservation_mutations
        code, out, err = run_receipt_fixture(
            reservation_mutation=mutation,
            expectation="--expect-grade-reserved",
            include_final=false,
        )
        @test code != 0
        @test isempty(out)
        @test occursin(field, err) || occursin("frozen", lowercase(err))
    end

    final_mutations = [
        "campaign_id" => (text -> replace(
            text, "campaign_id = \"fixture-campaign\"" =>
                "campaign_id = \"other-campaign\""; count=1)),
        "attempt_id" => (text -> replace(
            text, "attempt_id = \"fixture-attempt\"" =>
                "attempt_id = \"other-attempt\""; count=1)),
        "prediction_path" => (text -> replace(
            text, "prediction_path = \"$FIXTURE_PREDICTION\"" =>
                "prediction_path = \"$alternate_prediction\""; count=1)),
        "prediction_sha256" => (text -> replace(
            text, "prediction_sha256 = \"$(repeat('1', 64))\"" =>
                "prediction_sha256 = \"$(repeat('0', 64))\""; count=1)),
    ]
    for (field, mutation) in final_mutations
        code, out, err = run_receipt_fixture(final_mutation=mutation)
        @test code != 0
        @test isempty(out)
        @test occursin(field, err) || occursin("frozen", lowercase(err))
    end
end

@testset "correction2 binds the canonical grade reservation" begin
    alternate_root = joinpath(CANONICAL_ROOT, "other-root")
    mutations = [
        "repository_root" => (text -> replace(
            text, "repository_root = \"$CANONICAL_ROOT\"" =>
                "repository_root = \"$alternate_root\""; count=1)),
        "project_path" => (text -> replace(
            text, "project_path = \"$(joinpath(CANONICAL_ROOT, "Project.toml"))\"" =>
                "project_path = \"$(joinpath(CANONICAL_ROOT, "OtherProject.toml"))\"";
            count=1)),
        "grade_wrapper_path" => (text -> replace(
            text, "grade_wrapper_path = \"$GRADE_WRAPPER\"" =>
                "grade_wrapper_path = \"$(joinpath(CANONICAL_ROOT, "test", "other_wrapper.jl"))\"";
            count=1)),
        "report_script_path" => (text -> replace(
            text, "report_script_path = \"$REPORT_SCRIPT\"" =>
                "report_script_path = \"$(joinpath(CANONICAL_ROOT, "test", "other_report.jl"))\"";
            count=1)),
        "working_directory" => (text -> replace(
            text, "working_directory = \"$CANONICAL_ROOT\"" =>
                "working_directory = \"$(joinpath(CANONICAL_ROOT, "test"))\"";
            count=1)),
        "report_destination" => (text -> replace(
            text, "report_destination = \"$FIXTURE_GRADE\"" =>
                "report_destination = \"$(joinpath(dirname(FIXTURE_PREDICTION), "other-grade"))\"";
            count=1)),
        "exact_command" => (text -> replace(
            text, "exact_command = \"$FIXTURE_GRADE_COMMAND\"" =>
                "exact_command = \"true\""; count=1)),
    ]
    for (field, mutation) in mutations
        code, out, err = run_receipt_fixture(
            reservation_mutation=mutation,
            expectation="--expect-grade-reserved",
            include_final=false,
        )
        @test code != 0
        @test isempty(out)
        @test occursin(field, err) || occursin("canonical", lowercase(err))
    end
end

@testset "correction2 enforces phase-aware terminal outcomes" begin
    prelaunch_error = text -> replace_once_each(
        text,
        "outcome = \"success\"" => "outcome = \"error\"",
        "error_reason = \"\"" => "error_reason = \"crash_after_reservation_before_launch\"",
        "process_exit_code = 0" => "process_exit_code = -1",
        "invocation_count = 1" => "invocation_count = 0",
    )
    code0, out0, err0 = run_receipt_fixture(final_mutation=prelaunch_error)
    @test code0 == 0
    @test isempty(err0)
    @test occursin("lifecycle_state: graded", out0)

    launched_error = text -> replace_once_each(
        text,
        "outcome = \"success\"" => "outcome = \"error\"",
        "error_reason = \"\"" => "error_reason = \"grader_failed\"",
        "process_exit_code = 0" => "process_exit_code = 1",
    )
    code1, out1, err1 = run_receipt_fixture(final_mutation=launched_error)
    @test code1 == 0
    @test isempty(err1)
    @test occursin("lifecycle_state: graded", out1)

    parser_error = text -> replace_once_each(
        text,
        "outcome = \"success\"" => "outcome = \"error\"",
        "error_reason = \"\"" => "error_reason = \"invalid_grade_artifacts\"",
    )
    code2, out2, err2 = run_receipt_fixture(final_mutation=parser_error)
    @test code2 == 0
    @test isempty(err2)
    @test occursin("lifecycle_state: graded", out2)

    invalid_finals = [
        text -> replace(text, "process_exit_code = 0" => "process_exit_code = 1"; count=1),
        text -> replace(text, "error_reason = \"\"" => "error_reason = \"unexpected\""; count=1),
        text -> replace_once_each(
            text,
            "outcome = \"success\"" => "outcome = \"error\"",
            "error_reason = \"\"" => "error_reason = \"wrong_phase\"",
            "process_exit_code = 0" => "process_exit_code = -1",
            "invocation_count = 1" => "invocation_count = 0",
        ),
        text -> replace(text, "invocation_count = 1" => "invocation_count = 2"; count=1),
    ]
    for mutation in invalid_finals
        code, out, err = run_receipt_fixture(final_mutation=mutation)
        @test code != 0
        @test isempty(out)
        @test !isempty(err)
    end
end

@testset "correction3 rejects launched receipts with prelaunch-only fields" begin
    contradictions = [
        "prelaunch-only reason after launch" => (text -> replace_once_each(
            text,
            "outcome = \"success\"" => "outcome = \"error\"",
            "error_reason = \"\"" =>
                "error_reason = \"crash_after_reservation_before_launch\"",
            "process_exit_code = 0" => "process_exit_code = 1",
        )),
        "prelaunch-only exit sentinel after launch" => (text -> replace_once_each(
            text,
            "outcome = \"success\"" => "outcome = \"error\"",
            "error_reason = \"\"" => "error_reason = \"grader_failed\"",
            "process_exit_code = 0" => "process_exit_code = -1",
        )),
    ]
    for (case, mutation) in contradictions
        @testset "$case" begin
            code, out, err = run_receipt_fixture(final_mutation=mutation)
            @test code != 0
            @test isempty(out)
            @test !isempty(err)
        end
    end
end

exit(0)
