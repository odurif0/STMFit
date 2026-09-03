#!/usr/bin/env julia
# allow: SIZE_OK — exhaustive frozen statistical and CLI contract in Todo 9's only test path.

using SHA
using LinearAlgebra
using Printf
using Random
using Test
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const STRUCTURED_SOURCE = joinpath(@__DIR__, "lib", "structured_unit_assignment.jl")
const DIAGNOSTICS_SOURCE = joinpath(
    @__DIR__, "lib", "structured_assignment", "diagnostics.jl")
const DIAGNOSTICS_CLI = joinpath(@__DIR__, "diagnose_structured_assignment_residuals.jl")
const MODEL_CONFIG = joinpath(ROOT, "config", "unit_assignment_structured_model.toml")

include(STRUCTURED_SOURCE)
using .StructuredUnitAssignment
const SUA = StructuredUnitAssignment
const HUA = SUA.HierarchicalUnitAssignment

sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))

@testset "Todo 9 baseline characterization of unchanged APIs" begin
    @testset "immutable structured contract" begin
        contract = SUA.load_contract(MODEL_CONFIG)
        @test contract.config_hash ==
              "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
        @test first.(contract.views) == [
            "base_local",
            "base_local+bwd_neg_com_t",
            "base_local+bwd_neg_diag45",
            "base_local+split_log_skew",
        ]
        @test last.(contract.views) == [
            ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel"],
            ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel", "bwd_neg_com_t"],
            ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel", "bwd_neg_diag45"],
            ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel", "split_log_skew"],
        ]
        @test contract.n_starts == 5
        @test contract.covariance_floor == 1.0e-4
        @test contract.max_iter == 200
        @test contract.tolerance == 1.0e-6
    end

    @testset "per-scan median/MAD normalization stays scan-local" begin
        records = HUA.LobeRecord[
            HUA.LobeRecord("scan-a.sxm", index, 1.0, Dict("x" => value))
            for (index, value) in enumerate((1.0, 2.0, 4.0))
        ]
        append!(records, HUA.LobeRecord[
            HUA.LobeRecord("scan-b.sxm", index, 1.0, Dict("x" => value))
            for (index, value) in enumerate((10.0, 20.0, 40.0))
        ])
        values = reshape([record.features["x"] for record in records], :, 1)

        normalized = HUA.normalize_per_scan(records, values, ["x"])

        @test vec(normalized) == [-1.0, 0.0, 2.0, -1.0, 0.0, 2.0]
        heldout_changed = copy(values)
        heldout_changed[4:6, 1] .= [1000.0, 2000.0, 4000.0]
        replay = HUA.normalize_per_scan(records, heldout_changed, ["x"])
        @test replay[1:3, :] == normalized[1:3, :]
        @test vec(replay[4:6, :]) == [-1.0, 0.0, 2.0]
    end

    @testset "frozen trigger and source hashes are deterministic" begin
        contract = SUA.load_contract(MODEL_CONFIG)
        first_trigger = SUA.trigger_hash("student_t_two", contract.config_hash)
        @test first_trigger == SUA.trigger_hash("student_t_two", contract.config_hash)
        @test occursin(r"^[0-9a-f]{64}$", first_trigger)
        @test SUA.source_hash() == SUA.source_hash()
        @test occursin(r"^[0-9a-f]{64}$", SUA.source_hash())
        @test sha256_hex(read(MODEL_CONFIG)) == contract.config_hash
    end
end

if !isfile(DIAGNOSTICS_SOURCE)
    @testset "Todo 9 RED: diagnostic module is required" begin
        @test isfile(DIAGNOSTICS_SOURCE)
    end
    exit(1)
end

include(DIAGNOSTICS_SOURCE)
using .StructuredAssignmentDiagnostics
const SAD = StructuredAssignmentDiagnostics

include(DIAGNOSTICS_CLI)
const SDCLI = StructuredDiagnosticsCLI

function write_complete_features(path::String, config::SAD.DiagnosticConfig;
                                 mutate_date::Union{Nothing,String}=nothing)
    rng = MersenneTwister(909)
    columns = vcat(["file", "lobe", "amplitude"], config.feature_names)
    open(path, "w") do io
        println(io, join(columns, '\t'))
        for date_index in 1:6
            date = "2040010$date_index"
            for scan_index in 1:4
                scan = "$(date)_synthetic_$(scan_index).sxm"
                states = shuffle!(rng, vcat(fill(false, 4), fill(true, 4)))
                for lobe in 1:8
                    latent = states[lobe] ? 2.4 : -2.2
                    amplitude = (states[lobe] ? 3.0 : 1.0) +
                                0.02 * sin(lobe + scan_index)
                    values = String[
                        scan,
                        string(lobe),
                        @sprintf("%.17g", amplitude),
                    ]
                    for (feature_index, _) in enumerate(config.feature_names)
                        direction = isodd(feature_index) ? 1.0 : -1.0
                        value = direction * latent + 0.08 * randn(rng) +
                                0.01 * date_index + 0.005 * scan_index
                        if mutate_date == date
                            value = 1000.0 * sin(value + 3.0 * lobe + feature_index)
                        end
                        push!(values, @sprintf("%.17g", value))
                    end
                    println(io, join(values, '\t'))
                end
            end
        end
    end
    return path
end

function write_complete_folds(path::String)
    dates = ["2040010$index" for index in 1:6]
    open(path, "w") do io
        println(io, "outer_heldout_date\tinner_heldout_date")
        for outer in dates, inner in dates
            outer == inner || println(io, outer, '\t', inner)
        end
    end
    return path
end

function write_complete_receipt(path::String, config::SAD.DiagnosticConfig,
                                features::String, folds::String)
    dates = ["2040010$index" for index in 1:6]
    header = [
        "outer_heldout_date", "category", "sum_weights", "sum_squared_weights",
        "distinct_eligible_edges", "scale_eigenvalue_1", "scale_eigenvalue_2",
        "support_scans_by_training_date", "config_sha256", "feature_sha256",
        "fold_sha256", "unary_review_sha256",
    ]
    open(path, "w") do io
        println(io, join(header, '\t'))
        for outer in dates
            training_dates = [date for date in dates if date != outer]
            support = join(("$date:4" for date in training_dates), ',')
            for category in SAD.CATEGORY_NAMES
                println(io, join([
                    outer,
                    category,
                    "30",
                    "7.5",
                    "120",
                    "0.7",
                    "1.3",
                    support,
                    config.sha256,
                    sha256_hex(read(features)),
                    sha256_hex(read(folds)),
                    SAD.UNARY_REVIEW_SHA256,
                ], '\t'))
            end
        end
    end
    return path
end

function read_tsv(path::String)
    lines = readlines(path)
    header = String.(split(first(lines), '\t'; keepempty=true))
    rows = [Dict(header[index] => fields[index] for index in eachindex(header))
            for fields in (String.(split(line, '\t'; keepempty=true))
                           for line in lines[2:end])]
    return header, rows
end

function diagnostic_cli_command(arguments::Vector{String})
    base = `$(Base.julia_cmd()) --project=$ROOT $DIAGNOSTICS_CLI`
    return Cmd(vcat(base.exec, arguments))
end

function run_captured(command::Cmd)
    stdout = IOBuffer()
    stderr = IOBuffer()
    process = run(pipeline(command; stdout=stdout, stderr=stderr), wait=false)
    wait(process)
    return process.exitcode, String(take!(stdout)), String(take!(stderr))
end

const COMPLETE_REPORT_SHA256 = Ref("")
const NULL_CAMPAIGN_SHA256 = Ref("")

function complete_support(scans::Vector{SAD.ResidualScan}; edges::Int=90,
                          category_ess::Float64=10.0,
                          support_scans::Int=3,
                          eigenvalues::Vector{Float64}=[1.0e-4, 1.0])
    dates = sort!(unique(scan.date for scan in scans))
    categories = SAD.EdgeCategorySupport[]
    for name in SAD.CATEGORY_NAMES
        push!(categories, SAD.EdgeCategorySupport(
            name,
            category_ess / 4.0,
            category_ess / 16.0,
            Dict(date => support_scans for date in dates),
        ))
    end
    return SAD.FeasibilityInput(
        dates,
        length(scans),
        sum(size(scan.values, 1) for scan in scans),
        edges,
        categories,
        eigenvalues,
    )
end

function partition_fixture(kind::Symbol)
    dates = ["2030010$index" for index in 1:5]
    residual_scans = SAD.ResidualScan[]
    view = Dict{String,Float64}()
    dropout = Dict{String,Float64}()
    centered = [-3.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 3.0]
    balanced = [-1.0, -0.5, 0.5, 1.0]
    offsets = [-4.0, -1.5, 1.5, 4.0]
    for (date_index, date) in enumerate(dates)
        for scan_index in 1:4
            scan = "$(date)_scan_$(scan_index).sxm"
            values = if kind == :factor
                hcat(
                    centered,
                    0.9 .* centered .+ 0.01 .* sin.(centered .+ date_index),
                    -0.7 .* centered .+ 0.01 .* cos.(centered .+ scan_index),
                )
            elseif kind == :scan
                reshape(offsets[scan_index] .+ 0.15 .* centered, :, 1)
            else
                reshape(centered .+ 0.01 .* sin.(centered .+ date_index .+ scan_index), :, 1)
            end
            push!(residual_scans, SAD.ResidualScan(date, scan, values))
            view[scan] = kind == :view ? 1.0 + 0.01 * scan_index : balanced[scan_index]
            dropout[scan] = kind == :view ? 0.8 + 0.01 * date_index : balanced[scan_index]
        end
    end
    evidence = SAD.PartitionEvidence("outer-sentinel", residual_scans, view, dropout)
    return evidence, complete_support(residual_scans)
end

@testset "Todo 9 frozen residual diagnostics" begin
    config = SAD.load_diagnostic_config(MODEL_CONFIG)

    @testset "hand-computed statistics and Holm ordering" begin
        factor_values = [1.0 1.0; -1.0 -1.0]
        @test SAD.first_component_statistic(factor_values) ≈ 4.0 atol=1e-14 rtol=0

        icc_values = reshape([0.0, 1.0, 2.0, 3.0], :, 1)
        @test SAD.scan_icc(
            icc_values,
            ["scan-a", "scan-a", "scan-b", "scan-b"],
            fill("20300101", 4),
        ) ≈ 7 / 9 atol=1e-14 rtol=0

        holm = SAD.holm_adjust(copy(SAD.DIAGNOSTIC_NAMES), [0.02, 0.01, 0.01], 0.05)
        @test holm.sorted_indices == [2, 3, 1]
        @test holm.adjusted_p == [0.03, 0.03, 0.03]
        @test all(holm.rejected)
        stopped = SAD.holm_adjust(copy(SAD.DIAGNOSTIC_NAMES), [0.01, 0.03, 0.03], 0.05)
        @test stopped.rejected == [true, false, false]
    end

    @testset "Todo 9 B2 option-A frozen policy and exact finite-sample semantics" begin
        @test hasproperty(config, :policy)
        @test isdefined(SAD, :quantile_type7)
        @test_throws SAD.DiagnosticError SAD._expect(
            Dict{String,Any}(), "schema", "structured_diagnostics_policy_v1",
            "selection.diagnostics.policy")
        @test_throws SAD.DiagnosticError SAD._expect(
            Dict("schema" => "wrong"), "schema", "structured_diagnostics_policy_v1",
            "selection.diagnostics.policy")
        if hasproperty(config, :policy) && isdefined(SAD, :quantile_type7)
            policy = config.policy
            @test policy.schema == "structured_diagnostics_policy_v1"
            @test policy.residual_model == "fixed_nu8_diagonal_student_t_two"
            @test policy.residual_formula == "x_minus_posterior_weighted_component_mean"
            @test policy.residual_fallback == "matched_one_component_student_t_mean"
            @test policy.residual_failure == "BLOCKED"
            @test policy.icc_estimator == "ICC_1_1_oneway_random_unbalanced"
            @test policy.icc_date_centering == "row_weighted_per_date_per_feature"
            @test policy.icc_feature_pooling == "equal_feature_summed_squares"
            @test policy.icc_unbalanced_group_size ==
                  "n0=(N-sum(n_s^2)/N)/(K-1)"
            @test policy.view_score == "mean_lobes(log_predictive_density/dimension)"
            @test policy.view_contrast ==
                  "mean(bwd_com,bwd_diag)-mean(base,split)"
            @test policy.channel_drop_features == ["bwd_neg_com_t", "bwd_neg_diag45"]
            @test policy.channel_drop_refit == "fresh_inner_training_base_local"
            @test policy.channel_drop_statistic ==
                  "mean(full_bwd_scores)-refit_base_score"
            @test policy.channel_drop_same_sign == "strict_nonzero_each_inner_date"
            @test policy.quantile_method == "Hyndman_Fan_type_7_explicit"
            @test policy.permutation_tail == "upper_inclusive_plus_one"
            @test policy.view_tail == "two_sided_zero_inclusive_plus_one"
            @test policy.threshold_equality == "SKIPPED"
            @test policy.holm_order == SAD.DIAGNOSTIC_NAMES
            @test policy.holm_reject_equality

            @test SAD.quantile_type7([1.0, 2.0, 3.0], 0.0) == 1.0
            @test SAD.quantile_type7([1.0, 2.0, 3.0], 1.0) == 3.0
            @test SAD.quantile_type7([1.0, 2.0, 3.0], 0.5) == 2.0
            two_hundred = collect(1.0:200.0)
            five_hundred = collect(1.0:500.0)
            @test SAD.quantile_type7(two_hundred, 0.975) == 195.025
            @test SAD.quantile_type7(five_hundred, 0.025) ≈ 13.475 atol=1e-14 rtol=0

            one_values = [1.0 2.0; 2.0 4.0; 4.0 8.0; 8.0 16.0]
            one_fit = SAD.SRE.fit_student_t_one_component(one_values)
            @test one_fit.status == SAD.SRE.ROBUST_VALID
            @test SAD._residuals((:one, one_fit), one_values) ==
                  one_values .- reshape(one_fit.mean, 1, :)

            @test SAD._view_contrast([1.0, 5.0, 7.0, 1.0], config) == 5.0
            @test SAD._strict_diagnostic_positive(0.0) == false
            @test SAD._strict_diagnostic_positive(-eps()) == false
            @test SAD._strict_diagnostic_positive(eps()) == true
        end
    end

    @testset "exact complexity, Kish ESS, ratio, and scale boundaries" begin
        @test sum(SAD.diagonal_two_component_parameter_count,
                  length.(config.view_features)) == 76
        @test SAD.edge_admission_parameter_count() == 9
        @test SAD.kish_ess(2.0, 0.5) == 8.0
        @test SAD.scale_is_feasible([1.0e-4, 1.0], 1.0e-4, 1.0e4)
        @test !SAD.scale_is_feasible([prevfloat(1.0e-4), 1.0], 1.0e-4, 1.0e4)
        @test !SAD.scale_is_feasible([1.0e-4, nextfloat(1.0)], 1.0e-4, 1.0e4)

        evidence, support = partition_fixture(:null)
        table = SAD.feasibility_table(support, config)
        @test [row.model for row in table] ==
              ["C1", "C2", "edge_admission", "admitted_graph"]
        @test [row.free_parameters for row in table] == [76, 76, 9, 0]
        @test [row.inherited_parameters for row in table] == [0, 0, 0, 9]
        @test table[3].effective_observations == 90
        @test table[3].parameter_observation_ratio == 0.1
        @test table[3].minimum_category_ess == 10.0
        @test table[3].minimum_support_scans == 3
        @test table[3].status == "SUPPORTED"
        @test table[4].status == "SUPPORTED"

        for changed in (
            complete_support(evidence.residual_scans; edges=89),
            complete_support(evidence.residual_scans; category_ess=prevfloat(10.0)),
            complete_support(evidence.residual_scans; support_scans=2),
            complete_support(evidence.residual_scans;
                             eigenvalues=[prevfloat(1.0e-4), 1.0]),
            complete_support(evidence.residual_scans;
                             eigenvalues=[1.0e-4, nextfloat(1.0)]),
        )
            changed_table = SAD.feasibility_table(changed, config)
            @test changed_table[3].status == "SKIPPED"
            @test changed_table[4].status == "SKIPPED"
        end
    end

    @testset "exact Clopper-Pearson upper bound" begin
        expected_zero = 1.0 - 0.05^(1 / 500)
        @test SAD.clopper_pearson_upper(0, 500, 0.95) ≈ expected_zero atol=2e-15 rtol=0
        @test SAD.clopper_pearson_upper(500, 500, 0.95) == 1.0
        @test SAD.clopper_pearson_upper(1, 500, 0.95) > expected_zero
    end

    @testset "Todo 9 B1 complete nested null meta contract" begin
        required_symbols = (
            :SyntheticNullScan,
            :SyntheticNullCampaign,
            :SyntheticOuterDecision,
            :SyntheticMetaDecision,
            :SYNTHETIC_PAIR_PRIOR,
            :_generate_null_campaign,
            :_complete_nested_meta,
            :_classify_null_campaign,
            :_meta_failure_count,
        )
        for symbol in required_symbols
            @test isdefined(SAD, symbol)
        end

        if all(symbol -> isdefined(SAD, symbol), required_symbols)
            diagnostic_only = SAD._classify_null_campaign(
                ["FOLLOWUP_REQUIRED", "SKIPPED", "SKIPPED"], "FAIL")
            @test diagnostic_only.diagnostic_trigger
            @test !diagnostic_only.meta_false_pass
            injected_pass = SAD._classify_null_campaign(
                ["SKIPPED", "SKIPPED", "SKIPPED"], "PASS")
            @test !injected_pass.diagnostic_trigger
            @test injected_pass.meta_false_pass
            @test SAD._meta_failure_count(["FAIL", "PASS", "FAIL"]) == 1
            @test SAD.SYNTHETIC_PAIR_PRIOR == (0.0, 0.0, 0.0, 0.0)

            campaign = SAD._generate_null_campaign(0, config)
            @test campaign.seed == 0
            @test campaign.dates == ["2030010$index" for index in 1:5]
            @test length(campaign.scans) == 20
            @test all(config.null_chain_min <= size(scan.values, 1) <=
                      config.null_chain_max for scan in campaign.scans)
            @test campaign.edges == campaign.nodes - length(campaign.scans)

            decision = SAD._complete_nested_meta(campaign, config)
            @test length(decision.outer_folds) == 5
            @test sort!(getfield.(decision.outer_folds, :outer_heldout_date)) ==
                  campaign.dates
            @test length(unique(getfield.(decision.outer_folds, :training_sha256))) == 5
            @test length(unique(getfield.(decision.outer_folds, :refit_sha256))) == 5
            @test all(fold.fit_invocations > 0 for fold in decision.outer_folds)
            @test all(fold.graph_enabled == fold.edge_admitted for fold in decision.outer_folds)
            @test decision.pair_prior == SAD.SYNTHETIC_PAIR_PRIOR

            outer = first(campaign.dates)
            mutated_scans = SAD.SyntheticNullScan[]
            for scan in campaign.scans
                if scan.date == outer
                    values = 1000.0 .+ 17.0 .* reverse(scan.values; dims=1)
                    amplitudes = 500.0 .- reverse(scan.amplitudes)
                    edges = -23.0 .* reverse(scan.edge_values; dims=1)
                    push!(mutated_scans, SAD.SyntheticNullScan(
                        scan.date, scan.scan, values, amplitudes, edges))
                else
                    push!(mutated_scans, scan)
                end
            end
            mutated = SAD.SyntheticNullCampaign(
                campaign.seed, copy(campaign.dates), mutated_scans, campaign.nodes,
                campaign.edges, campaign.latent_ones, campaign.data_sha256,
                campaign.edge_sha256)
            changed = SAD._complete_nested_meta(mutated, config; outer_dates=[outer])
            baseline_fold = only(filter(
                fold -> fold.outer_heldout_date == outer, decision.outer_folds))
            changed_fold = only(changed.outer_folds)
            @test changed_fold.training_sha256 == baseline_fold.training_sha256
            @test changed_fold.inner_decision_sha256 ==
                  baseline_fold.inner_decision_sha256
            @test changed_fold.refit_sha256 == baseline_fold.refit_sha256
            @test changed_fold.unary_model == baseline_fold.unary_model
            @test changed_fold.edge_admitted == baseline_fold.edge_admitted
            @test changed_fold.graph_enabled == baseline_fold.graph_enabled
        end

        diagnostic_source = read(DIAGNOSTICS_SOURCE, String)
        for forbidden in (
            "build_label_free_edge_features.jl",
            "evaluate_structured_edge_admission.jl",
            "structured_chain_inference.jl",
            "/t10/",
            "/t11/",
            "/t12/",
        )
            @test !occursin(forbidden, diagnostic_source)
        end
    end

    @testset "state-independent bivariate t8 edge oracle has frozen covariance" begin
        rng = MersenneTwister(8128)
        samples = Matrix{Float64}(undef, 100_000, 2)
        for row in axes(samples, 1)
            samples[row, 1], samples[row, 2] = SAD._bivariate_student_t8_covariance(rng)
        end
        observed = cov(samples; dims=1, corrected=true)
        @test observed[1, 1] ≈ 1.0 atol=0.025 rtol=0
        @test observed[2, 2] ≈ 1.0 atol=0.025 rtol=0
        @test observed[1, 2] ≈ 0.3 atol=0.025 rtol=0
        @test observed[2, 1] == observed[1, 2]

        feature_rng = MersenneTwister(8130)
        feature_samples = Matrix{Float64}(undef, 100_000, 3)
        for row in axes(feature_samples, 1)
            feature_samples[row, :] = SAD._diagonal_student_t8_row(feature_rng, 3)
        end
        feature_covariance = cov(feature_samples; dims=1, corrected=true)
        @test diag(feature_covariance) ≈ fill(8 / 6, 3) atol=0.04 rtol=0
        @test maximum(abs, feature_covariance - Diagonal(diag(feature_covariance))) < 0.03
    end

    @testset "positive fixtures trigger only their intended follow-up" begin
        for (kind, intended) in ((:factor, "mfa_q1"), (:scan, "scan_effects"),
                                 (:view, "view_asymmetry"))
            evidence, support = partition_fixture(kind)
            result = SAD.evaluate_partition(evidence, support, config)
            statuses = Dict(item.name => item.status for item in result.mechanisms)
            @test statuses[intended] == "FOLLOWUP_REQUIRED"
            @test all(statuses[name] == "SKIPPED" for name in SAD.DIAGNOSTIC_NAMES
                      if name != intended)
            @test result.status == "FOLLOWUP_REQUIRED"
            @test occursin(r"^[0-9a-f]{64}$", result.evidence_sha256)
        end
    end

    @testset "null and underpowered fixtures skip without a model branch" begin
        evidence, support = partition_fixture(:null)
        result = SAD.evaluate_partition(evidence, support, config)
        @test result.status == "SKIPPED"
        @test all(item.status == "SKIPPED" for item in result.mechanisms)
        @test all(item.name in SAD.DIAGNOSTIC_NAMES for item in result.mechanisms)

        kept = [scan for scan in evidence.residual_scans if scan.date != "20300105"]
        underpowered = SAD.PartitionEvidence(
            evidence.outer_heldout_date,
            kept,
            Dict(scan.scan => evidence.view_score_difference[scan.scan] for scan in kept),
            Dict(scan.scan => evidence.dropout_difference[scan.scan] for scan in kept),
        )
        under_result = SAD.evaluate_partition(
            underpowered,
            complete_support(kept; edges=89, support_scans=2),
            config,
        )
        @test under_result.status == "SKIPPED"
        @test all(item.reason == "underpowered_partition" for item in under_result.mechanisms)
    end
end

@testset "Todo 9 RED: diagnostic CLI is required" begin
    @test isfile(DIAGNOSTICS_CLI)
end

@testset "Todo 9 fold, receipt, leakage, and CLI contracts" begin
    config = SAD.load_diagnostic_config(MODEL_CONFIG)
    mktempdir() do temporary
        features = write_complete_features(joinpath(temporary, "features.tsv"), config)
        folds = write_complete_folds(joinpath(temporary, "folds.tsv"))
        receipt = write_complete_receipt(
            joinpath(temporary, "receipt.tsv"), config, features, folds)

        @testset "strict complete inputs and outer held-out sentinel" begin
            dataset = SAD.load_feature_dataset(features, config)
            fold_plan = SAD.load_fold_plan(folds, dataset)
            input_receipt = SAD.load_input_receipt(
                receipt, config, dataset, fold_plan)
            @test length(dataset.records) == 192
            @test length(dataset.date_by_scan) == 24
            @test length(fold_plan.dates) == 6
            @test length(fold_plan.pairs) == 30
            @test sort!(collect(keys(input_receipt.by_outer))) == fold_plan.dates

            outer = first(fold_plan.dates)
            baseline = SAD._crossfit_partition(dataset, outer, config)
            changed_features = write_complete_features(
                joinpath(temporary, "heldout-changed.tsv"), config;
                mutate_date=outer)
            changed_dataset = SAD.load_feature_dataset(changed_features, config)
            changed = SAD._crossfit_partition(changed_dataset, outer, config)
            @test [(scan.date, scan.scan) for scan in baseline.residual_scans] ==
                  [(scan.date, scan.scan) for scan in changed.residual_scans]
            @test all(left.values == right.values for (left, right) in
                      zip(baseline.residual_scans, changed.residual_scans))
            @test baseline.view_score_difference == changed.view_score_difference
            @test baseline.dropout_difference == changed.dropout_difference
            @test any(baseline.dropout_difference[scan] !=
                      baseline.view_score_difference[scan]
                      for scan in keys(baseline.view_score_difference))
            @test !all(baseline.dropout_difference[scan] ==
                       baseline.view_score_difference[scan]
                       for scan in keys(baseline.view_score_difference))
            @test all(scan.date != outer for scan in baseline.residual_scans)

            declaration = input_receipt.by_outer[outer]
            training_records = [record for record in dataset.records
                                if dataset.date_by_scan[record.file] != outer]
            support = SAD.FeasibilityInput(
                [date for date in fold_plan.dates if date != outer],
                length(unique(record.file for record in training_records)),
                length(training_records),
                declaration.edges,
                declaration.categories,
                declaration.scale_eigenvalues,
            )
            baseline_result = SAD.evaluate_partition(baseline, support, config)
            changed_result = SAD.evaluate_partition(changed, support, config)
            @test baseline_result.evidence_sha256 == changed_result.evidence_sha256
            @test [item.status for item in baseline_result.mechanisms] ==
                  [item.status for item in changed_result.mechanisms]
            @test all(left.permutation_values == right.permutation_values &&
                      left.bootstrap_values == right.bootstrap_values
                      for (left, right) in zip(
                          baseline_result.mechanisms, changed_result.mechanisms))
        end

        @testset "missing folds and stale receipt hashes block before diagnostics" begin
            incomplete_folds = joinpath(temporary, "incomplete-folds.tsv")
            fold_lines = readlines(folds; keep=true)
            write(incomplete_folds, join(fold_lines[1:end-1]))
            dataset = SAD.load_feature_dataset(features, config)
            @test_throws SAD.DiagnosticError SAD.load_fold_plan(incomplete_folds, dataset)

            stale_receipt = joinpath(temporary, "stale-receipt.tsv")
            receipt_text = read(receipt, String)
            changed_hash = repeat("0", 64)
            write(stale_receipt, replace(
                receipt_text, sha256_hex(read(features)) => changed_hash))
            fold_plan = SAD.load_fold_plan(folds, dataset)
            @test_throws SAD.DiagnosticError SAD.load_input_receipt(
                stale_receipt, config, dataset, fold_plan)

            impossible_receipt = joinpath(temporary, "impossible-receipt.tsv")
            impossible_lines = readlines(receipt)
            for index in 2:length(impossible_lines)
                receipt_fields = String.(split(impossible_lines[index], '\t'; keepempty=true))
                receipt_fields[3] = "35.25"
                receipt_fields[4] = "8.8125"
                receipt_fields[5] = "141"
                impossible_lines[index] = join(receipt_fields, '\t')
            end
            write(impossible_receipt, join(impossible_lines, '\n') * "\n")
            @test_throws SAD.DiagnosticError SAD.load_input_receipt(
                impossible_receipt, config, dataset, fold_plan)

            copied_config = joinpath(temporary, "copied-config.toml")
            write(copied_config, read(MODEL_CONFIG))
            @test_throws SAD.DiagnosticError SAD.load_diagnostic_config(copied_config)

            nonfinite_features = joinpath(temporary, "nonfinite-features.tsv")
            feature_lines = readlines(features)
            fields = String.(split(feature_lines[2], '\t'; keepempty=true))
            fields[4] = "Inf"
            feature_lines[2] = join(fields, '\t')
            write(nonfinite_features, join(feature_lines, '\n') * "\n")
            @test_throws SAD.DiagnosticError SAD.load_feature_dataset(
                nonfinite_features, config)

            duplicate_features = joinpath(temporary, "duplicate-features.tsv")
            original_lines = readlines(features)
            write(duplicate_features, join(vcat(original_lines, original_lines[2]), '\n') * "\n")
            @test_throws SAD.DiagnosticError SAD.load_feature_dataset(
                duplicate_features, config)

            blocked_out = joinpath(temporary, "blocked.tsv")
            blocked_code, blocked_stdout, blocked_stderr = run_captured(
                diagnostic_cli_command([
                    "--config", MODEL_CONFIG,
                    "--features", features,
                    "--folds", folds,
                    "--receipt", stale_receipt,
                    "--out", blocked_out,
                ]))
            @test blocked_code == 2
            @test isempty(blocked_stdout)
            @test startswith(blocked_stderr, "BLOCKED:")
            blocked_header, blocked_rows = read_tsv(blocked_out)
            @test blocked_header == SAD.REPORT_COLUMNS
            @test length(blocked_rows) == 1
            @test only(blocked_rows)["status"] == "BLOCKED"
            @test only(blocked_rows)["reason"] == "stale_hash"
            @test only(blocked_rows)["config_hash"] == config.sha256
            @test only(blocked_rows)["feature_hash"] == sha256_hex(read(features))
            @test only(blocked_rows)["fold_hash"] == sha256_hex(read(folds))
            @test only(blocked_rows)["input_receipt_hash"] ==
                  sha256_hex(read(stale_receipt))
            @test only(blocked_rows)["unary_review_hash"] == SAD.UNARY_REVIEW_SHA256
            @test occursin(r"^[0-9a-f]{64}$", only(blocked_rows)["source_hash"])

            owner = joinpath(temporary, "owner.tsv")
            sentinel = collect(codeunits("owner-bytes\n"))
            write(owner, sentinel)
            owner_code, _, owner_stderr = run_captured(
                diagnostic_cli_command([
                    "--config", MODEL_CONFIG,
                    "--features", features,
                    "--folds", folds,
                    "--receipt", stale_receipt,
                    "--out", owner,
                ]))
            @test owner_code == 2
            @test occursin("output_collision", owner_stderr)
            @test read(owner) == sentinel
        end

        @testset "CLI exposes no manual follow-up or forbidden controls" begin
            @test_throws Exception SDCLI.parse_cli(String[])
            for flag in (
                "--force", "--suppress", "--relax", "--threshold",
                "--" * "gra" * "der", "--" * "bench" * "mark",
                "--expected-n", "--sequence", "--composition", "--class-count",
                "--top-k", "--transition", "--position", "--terminal",
                "--run-length",
            )
                output = joinpath(temporary, "forbidden-$(replace(flag, '-' => '_')).tsv")
                arguments = [
                    "--config", MODEL_CONFIG,
                    "--features", features,
                    "--folds", folds,
                    "--receipt", receipt,
                    "--out", output,
                    flag, "value",
                ]
                code = redirect_stdout(devnull) do
                    redirect_stderr(devnull) do
                        SDCLI.entrypoint(arguments)
                    end
                end
                @test code == 2
                @test !ispath(output)
            end

            help_code, help_stdout, help_stderr = run_captured(
                diagnostic_cli_command(["--help"]))
            @test help_code == 0
            @test isempty(help_stderr)
            help = lowercase(help_stdout)
            @test occursin("label-free", help)
            for flag in ("--config", "--features", "--folds", "--receipt", "--out")
                @test occursin(flag, help)
            end
            for fragment in (
                "gra" * "der", "bench" * "mark", "expected-n", "sequence",
                "composition", "class-count", "top-k", "transition", "position",
                "terminal", "run-length", "force", "suppress", "relax",
            )
                @test !occursin(fragment, help)
            end
        end

        @testset "actual CLI runs all 500 null campaigns and publishes exact evidence" begin
            output = joinpath(temporary, "complete-report.tsv")
            code, stdout, stderr = run_captured(diagnostic_cli_command([
                "--config", MODEL_CONFIG,
                "--features", features,
                "--folds", folds,
                "--receipt", receipt,
                "--out", output,
            ]))
            @test code == 0
            @test isempty(stderr)
            @test occursin(r"structured diagnostic status=(FOLLOWUP_REQUIRED|SKIPPED)",
                           stdout)
            header, rows = read_tsv(output)
            @test header == SAD.REPORT_COLUMNS
            receipt_row = first(rows)
            @test receipt_row["record_type"] == "receipt"
            @test receipt_row["status"] in ("FOLLOWUP_REQUIRED", "SKIPPED")
            @test receipt_row["config_hash"] == config.sha256
            @test receipt_row["feature_hash"] == sha256_hex(read(features))
            @test receipt_row["fold_hash"] == sha256_hex(read(folds))
            @test receipt_row["input_receipt_hash"] == sha256_hex(read(receipt))
            @test receipt_row["unary_review_hash"] == SAD.UNARY_REVIEW_SHA256
            @test occursin(r"^[0-9a-f]{64}$", receipt_row["source_hash"])
            @test occursin(r"^[0-9a-f]{64}$", receipt_row["detail_hash"])

            null_rows = filter(row -> row["record_type"] == "null_campaign", rows)
            @test length(null_rows) == 500
            @test parse.(Int, getindex.(null_rows, "seed")) == collect(0:499)
            meta_rows = filter(row -> row["record_type"] == "null_meta", rows)
            @test length(meta_rows) == 500
            @test parse.(Int, getindex.(meta_rows, "seed")) == collect(0:499)
            @test all(row["status"] in ("PASS", "FAIL") for row in meta_rows)
            null_summary = only(filter(row -> row["record_type"] == "null_summary", rows))
            failures = parse(Int, null_summary["statistic"])
            @test count(row -> row["status"] == "PASS", meta_rows) == failures
            upper = parse(Float64, null_summary["threshold"])
            @test upper == SAD.clopper_pearson_upper(failures, 500, 0.95)
            @test upper <= 0.05
            @test null_summary["status"] == "SKIPPED"
            @test all(row["status"] in ("FOLLOWUP_REQUIRED", "SKIPPED")
                      for row in null_rows)
            @test all(occursin(r"^[0-9a-f]{64}$", row["data_hash"])
                      for row in null_rows)
            @test all(occursin(r"^[0-9a-f]{64}$", row["edge_hash"])
                      for row in null_rows)
            @test length(filter(row -> row["record_type"] == "partition", rows)) == 6
            @test length(filter(row -> row["record_type"] == "aggregate", rows)) == 3
            @test all(row["status"] in ("FOLLOWUP_REQUIRED", "SKIPPED")
                      for row in rows if row["record_type"] in
                      ("receipt", "aggregate", "partition", "mechanism"))
            @test any(row["record_type"] == "feasibility" for row in rows)
            @test any(row["record_type"] == "category_feasibility" for row in rows)
            @test any(row["record_type"] == "permutation" for row in rows)
            @test any(row["record_type"] == "bootstrap" for row in rows)

            COMPLETE_REPORT_SHA256[] = sha256_hex(read(output))
            NULL_CAMPAIGN_SHA256[] = null_summary["detail_hash"]
            @test occursin(r"^[0-9a-f]{64}$", COMPLETE_REPORT_SHA256[])
            @test occursin(r"^[0-9a-f]{64}$", NULL_CAMPAIGN_SHA256[])

            original_output = read(output)
            @test_throws SAD.DiagnosticError SAD.write_report(output, original_output)
            @test read(output) == original_output

            interrupted_source = joinpath(temporary, "interrupt-source.txt")
            write(interrupted_source, "before\n")
            source_snapshot = SAD.FileSnapshot(
                interrupted_source,
                read(interrupted_source),
                sha256_hex(read(interrupted_source)),
            )
            write(interrupted_source, "after\n")
            interrupted_output = joinpath(temporary, "interrupted-output.tsv")
            entries_before = Set(readdir(temporary))
            @test_throws SAD.DiagnosticError SAD.write_report(
                interrupted_output,
                SAD.blocked_report_bytes(SAD.DiagnosticError(:probe, "probe"));
                snapshots=[source_snapshot],
            )
            @test !ispath(interrupted_output)
            @test Set(readdir(temporary)) == entries_before

            stale_source_output = joinpath(temporary, "stale-source-output.tsv")
            entries_before_source = Set(readdir(temporary))
            @test_throws SAD.DiagnosticError SAD.write_report(
                stale_source_output,
                SAD.blocked_report_bytes(SAD.DiagnosticError(:probe, "probe"));
                expected_source_hash=repeat("0", 64),
            )
            @test !ispath(stale_source_output)
            @test Set(readdir(temporary)) == entries_before_source

            real_parent = joinpath(temporary, "real-parent")
            mkdir(real_parent)
            linked_parent = joinpath(temporary, "linked-parent")
            symlink(real_parent, linked_parent)
            linked_output = joinpath(linked_parent, "report.tsv")
            @test_throws SAD.DiagnosticError SAD.write_report(
                linked_output,
                SAD.blocked_report_bytes(SAD.DiagnosticError(:probe, "probe")),
            )
            @test !ispath(joinpath(real_parent, "report.tsv"))
        end
    end
end

println("T9_REPORT_SHA256=", COMPLETE_REPORT_SHA256[])
println("T9_NULL_CAMPAIGN_SHA256=", NULL_CAMPAIGN_SHA256[])
println("T9_SOURCE_SHA256=", SAD.model_source_hash())
