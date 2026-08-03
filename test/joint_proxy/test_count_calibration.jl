#!/usr/bin/env julia

# TDD tests for the synthetic count-calibration facade (Todo 5).
#
# Scope guardrails:
# - synthetic CandidateViewReport fixtures only
# - no real SXM / summary / manifest / truth imports
# - count calibration only; no type/mold evidence is mixed into N

using Test

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "count_calibration.jl"))
using .JointProxyCountCalibration:
    CountCalibrationCase, CountCalibrationError, CountCalibrationModel,
    CountPosteriorRow, CountCalibrationDiagnostics, CountReliabilityBin,
    score_count_report, calibrate_count_temperature, predict_count

struct MiniCandidate
    n::Int
    joint_gcv::Float64
    valid::Bool
end

struct MiniReport
    candidates::Vector{MiniCandidate}
end

_cand(n::Int, gcv::Float64; valid::Bool=true) = MiniCandidate(n, gcv, valid)

function _report(ns::Vector{Int}, gcvs::Vector{Float64}; order::Vector{Int}=collect(1:length(ns)))
    cands = [_cand(ns[i], gcvs[i]) for i in order]
    return MiniReport(cands)
end

function _case(case_id::String, split::Symbol, true_n::Int, report;
              seed::Integer=0x17, config_hash::String="cfg-a")
    return CountCalibrationCase(case_id, split, true_n, report, UInt64(seed), config_hash)
end

function _good_cases()
    cases = CountCalibrationCase[]
    push!(cases, _case("cal-1", :calibration, 3, _report([2,3,4], [0.20, 0.02, 0.30])))
    push!(cases, _case("cal-2", :calibration, 4, _report([3,4,5], [0.40, 0.01, 0.28])))
    push!(cases, _case("held-1", :heldout, 3, _report([2,3,4], [0.30, 0.04, 0.33])))
    push!(cases, _case("held-2", :heldout, 4, _report([3,4,5], [0.42, 0.02, 0.25])))
    return cases
end

@testset "JointProxyCountCalibration — Todo 5" begin
    @testset "red-first guardrails before implementation" begin
        @test_throws Exception calibrate_count_temperature(CountCalibrationCase[];
            seed=0x17, config_hash="cfg-a")
    end

    @testset "happy synthetic fit, scoring, and replay" begin
        cases = _good_cases()
        model = calibrate_count_temperature(cases; seed=0x17, config_hash="cfg-a")
        @test model isa CountCalibrationModel
        @test isfinite(model.temperature) && model.temperature > 0
        @test isfinite(model.confidence_threshold)
        @test 0.0 <= model.confidence_threshold <= 1.0
        @test isfinite(model.heldout_nll)
        @test isfinite(model.baseline_temp1_nll)
        @test isfinite(model.calibration_nll)
        @test isfinite(model.heldout_nll)
        @test model.heldout_nll <= model.baseline_temp1_nll + 1e-9
        @test !isempty(model.reliability_bins)
        @test all(b -> b isa CountReliabilityBin, model.reliability_bins)

        scored = score_count_report(cases[1].report, model.temperature;
            threshold=model.confidence_threshold, true_n=cases[1].true_n)
        @test scored.rows isa Vector{CountPosteriorRow}
        @test isapprox(sum(r.probability for r in scored.rows), 1.0; atol=1e-12)
        @test first(scored.rows).rank == 1
        @test scored.true_rank == 1
        @test scored.predicted_n == cases[1].true_n
        @test !scored.abstained
        @test predict_count(cases[1].report, model; true_n=cases[1].true_n).predicted_n == cases[1].true_n
        
        replay = calibrate_count_temperature(cases; seed=0x17, config_hash="cfg-a")
        @test replay.temperature == model.temperature
        @test replay.confidence_threshold == model.confidence_threshold
        @test replay.heldout_nll == model.heldout_nll
    end

    @testset "true N ranks first under low-noise synthetic scores" begin
        report = _report([2,3,4,5], [0.35, 0.01, 0.12, 0.20])
        scored = score_count_report(report, 0.75; true_n=3)
        @test scored.rows[1].n == 3
        @test scored.rows[1].rank == 1
        @test scored.rows[1].probability == maximum(r.probability for r in scored.rows)
    end

    @testset "tied scores abstain rather than follow row order" begin
        left = _report([3,4], [0.10, 0.10], order=[1,2])
        right = _report([3,4], [0.10, 0.10], order=[2,1])
        scored_left = score_count_report(left, 1.0; threshold=0.50, true_n=3)
        scored_right = score_count_report(right, 1.0; threshold=0.50, true_n=3)
        @test scored_left.abstained
        @test scored_right.abstained
        @test scored_left.predicted_n === missing
        @test scored_right.predicted_n === missing
    end

    @testset "explicit failure modes" begin
        @test_throws CountCalibrationError calibrate_count_temperature([
            _case("bad-true", :calibration, 9, _report([2,3,4], [0.2, 0.1, 0.3]))
        ]; seed=0x17, config_hash="cfg-a")

        @test_throws CountCalibrationError calibrate_count_temperature([
            _case("all-inf", :calibration, 3, _report([3,4], [Inf, Inf]))
        ]; seed=0x17, config_hash="cfg-a")

        @test_throws CountCalibrationError calibrate_count_temperature([
            _case("one-cand", :calibration, 3, _report([3], [0.1]))
        ]; seed=0x17, config_hash="cfg-a")

        @test_throws CountCalibrationError calibrate_count_temperature(
            _good_cases(); seed=0x99, config_hash="cfg-a")
        @test_throws CountCalibrationError calibrate_count_temperature(
            _good_cases(); seed=0x17, config_hash="cfg-b")
    end
end
