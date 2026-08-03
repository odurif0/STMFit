#!/usr/bin/env julia

# TDD tests for the synthetic type-calibration facade (Todo 6).
#
# Scope guardrails:
# - synthetic TypePosteriorResult / ProxyEnsemble fixtures only
# - no real SXM / summary / manifest / truth imports
# - type calibration only; no real-label tuning of exact-chain, composition,
#   flip, family, height, or confidence rules

using Test

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "proxy_registry.jl"))
include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "type_posterior.jl"))
include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "type_calibration.jl"))

using .JointProxyRegistry: ProxySource, ProxyTemplate, ProxyEntry, ProxyEnsemble
using .JointProxyTypePosterior: TypePosteriorPriors, TypePosteriorResult,
                                 TypePosteriorLobeEvidence,
                                 infer_type_posterior
using .JointProxyTypeCalibration: TypeCalibrationCase, TypeCalibrationError,
                                  TypeCalibrationModel, TypeCasePosterior,
                                  TypePosteriorRow, TypeReliabilityBin,
                                  TypeCalibrationDiagnostics,
                                  score_type_case, calibrate_type_temperature,
                                  predict_type

const FAMILY_NAMES = ("geom", "stm")

_basis(i, n) = [j == i ? 1.0 : 0.0 for j in 1:n]

function _template_vec(type::Int, parity::Int, mirror::Int; source_bias::Float64=0.0)
    idx = 1 + type + 2 * parity + 4 * mirror
    v = _basis(idx, 8)
    if source_bias != 0.0
        v[8] += source_bias
    end
    return v
end

function _ensemble(; swap::Bool=false, identical::Bool=false, bias::Float64=0.0)
    entries = ProxyEntry[]
    tmpls = ProxyTemplate[]
    for typ in 0:1, parity in 0:1, mirror in 0:1
        t = identical ? typ : (swap ? 1 - typ : typ)
        push!(tmpls, ProxyTemplate(typ, parity, mirror,
                                   identical ? fill(1.0, 8) : _template_vec(t, parity, mirror; source_bias=bias)))
    end
    push!(entries, ProxyEntry(ProxySource("geom", "geom", "geom.tsv", "sha_geom", 0.5, 0.5, true), tmpls))
    push!(entries, ProxyEntry(ProxySource("stm", "stm_h050", "stm.tsv", "sha_stm", 0.5, 0.5, true), tmpls))
    return ProxyEnsemble(entries, 0.32, 0.08, 9, 81, "payload", String[])
end

function _evidence(n::Int; duplicate::Bool=false, reversed::Bool=false, swap_second::Bool=false,
                   identical::Bool=false, bias::Float64=0.15)
    true_types = isodd(n) ? [mod1(i, 2) - 1 for i in 1:n] : [0 for _ in 1:n]
    if reversed
        true_types = reverse(true_types)
    end
    lobes = TypePosteriorLobeEvidence[]
    for (i, typ) in enumerate(true_types)
        parity = mod(i - 1, 2)
        true_t = _template_vec(typ, parity, 0)
        false_t = _template_vec(1 - typ, parity, 0)
        fwd = identical ? fill(1.0, 8) : (0.75 .* true_t .+ 0.25 .* false_t)
        bwd = duplicate ? copy(fwd) : (swap_second ? (0.20 .* true_t .+ 0.80 .* false_t) : (0.90 .* true_t .+ 0.10 .* false_t))
        push!(lobes, TypePosteriorLobeEvidence(Dict("fwd" => fwd, "bwd" => bwd)))
    end
    rho = duplicate ? 1.0 : 0.0
    return lobes, true_types, rho
end

function _case(case_id::String, split::Symbol, true_types::Vector{Int};
               seed::Integer=0x17, config_hash::String="cfg-a",
               manifest_hash::String="man-a", source_hash::String="src-a",
               control::Symbol=:ok, family_reports::Vector{Pair{String,TypePosteriorResult}})
    return TypeCalibrationCase(case_id=case_id, split=split, true_types=true_types,
        family_reports=family_reports, seed=UInt64(seed), config_hash=config_hash,
        manifest_hash=manifest_hash, source_hash=source_hash, control=control)
end

function _family_reports(n::Int; swap_geom::Bool=false, swap_stm::Bool=false,
                         identical::Bool=false, duplicate::Bool=false,
                         reversed::Bool=false, bias::Float64=0.15)
    lobes, true_types, rho = _evidence(n; duplicate=duplicate, reversed=reversed,
        swap_second=swap_stm, identical=identical, bias=bias)
    geom = infer_type_posterior(lobes, _ensemble(; swap=swap_geom, identical=identical, bias=bias);
                                priors=TypePosteriorPriors(), rho=rho)
    stm = infer_type_posterior(lobes, _ensemble(; swap=swap_stm, identical=identical, bias=bias);
                               priors=TypePosteriorPriors(), rho=rho)
    return true_types, [FAMILY_NAMES[1] => geom, FAMILY_NAMES[2] => stm]
end

function _good_cases()
    cal_true, cal_reports = _family_reports(4)
    held_true, held_reports = _family_reports(4; bias=0.20)
    null_true, null_reports = _family_reports(3; identical=true)
    conflict_true, conflict_reports = _family_reports(4; swap_stm=true)
    return TypeCalibrationCase[
        _case("cal-1", :calibration, cal_true; family_reports=cal_reports),
        _case("held-1", :heldout, held_true; family_reports=held_reports),
        _case("held-null", :heldout, null_true; control=:null, family_reports=null_reports),
        _case("held-conflict", :heldout, conflict_true; control=:family_conflict, family_reports=conflict_reports),
    ]
end

@testset "JointProxyTypeCalibration — Todo 6" begin
    @testset "red-first guardrails before implementation" begin
        @test_throws Exception calibrate_type_temperature(TypeCalibrationCase[];
            seed=0x17, config_hash="cfg-a", manifest_hash="man-a", source_hash="src-a")
    end

    @testset "happy synthetic fit, scoring, and replay" begin
        cases = _good_cases()
        model = calibrate_type_temperature(cases; seed=0x17, config_hash="cfg-a",
            manifest_hash="man-a", source_hash="src-a")
        @test model isa TypeCalibrationModel
        @test isfinite(model.temperature) && model.temperature > 0
        @test isfinite(model.confidence_threshold)
        @test 0.0 <= model.confidence_threshold <= 1.0
        @test isfinite(model.heldout_nll)
        @test isfinite(model.baseline_temp1_nll)
        @test isfinite(model.calibration_nll)
        @test model.heldout_nll <= model.baseline_temp1_nll + 1e-9
        @test !isempty(model.reliability_bins)
        @test all(b -> b isa TypeReliabilityBin, model.reliability_bins)

        scored = score_type_case(cases[1], model)
        @test scored isa TypeCasePosterior
        @test scored.rows isa Vector{TypePosteriorRow}
        @test isapprox(sum(r.probability_0 + r.probability_1 for r in scored.rows), length(scored.rows); atol=1e-12)
        @test all(isapprox(r.probability_0 + r.probability_1, 1.0; atol=1e-12) for r in scored.rows)
        @test all(0.0 .<= [r.confidence for r in scored.rows] .<= 1.0)
        @test scored.predicted_types == cases[1].true_types
        @test !any(ismissing, scored.predicted_types)

        replay = calibrate_type_temperature(cases; seed=0x17, config_hash="cfg-a",
            manifest_hash="man-a", source_hash="src-a")
        @test replay.temperature == model.temperature
        @test replay.confidence_threshold == model.confidence_threshold
        @test replay.heldout_nll == model.heldout_nll
    end

    @testset "confidence is monotone in separation and low-noise types cross threshold" begin
        weak_true, weak_reports = _family_reports(3; identical=true)
        strong_true, strong_reports = _family_reports(3; bias=0.30)
        weak = score_type_case(_case("weak", :heldout, weak_true; control=:null, family_reports=weak_reports),
            calibrate_type_temperature(_good_cases(); seed=0x17, config_hash="cfg-a",
                manifest_hash="man-a", source_hash="src-a"))
        strong = score_type_case(_case("strong", :heldout, strong_true; family_reports=strong_reports),
            calibrate_type_temperature(_good_cases(); seed=0x17, config_hash="cfg-a",
                manifest_hash="man-a", source_hash="src-a"))
        @test strong.rows[1].confidence >= weak.rows[1].confidence
        @test all(!ismissing, predict_type(_case("strong-pred", :heldout, strong_true; family_reports=strong_reports),
            calibrate_type_temperature(_good_cases(); seed=0x17, config_hash="cfg-a",
                manifest_hash="man-a", source_hash="src-a")).predicted_types)
    end

    @testset "mixed confidence emits strong lobes and ablates only conflicted weak lobe" begin
        model = calibrate_type_temperature(_good_cases(); seed=0x17, config_hash="cfg-a",
            manifest_hash="man-a", source_hash="src-a")
        mixed_true, mixed_reports = _family_reports(4; bias=0.30)
        mixed_reports[1].second.lobe_marginals[2, 1] = 0.10
        mixed_reports[1].second.lobe_marginals[2, 2] = 0.90
        mixed_reports[2].second.lobe_marginals[2, 1] = 0.90
        mixed_reports[2].second.lobe_marginals[2, 2] = 0.10
        mixed = score_type_case(_case("mixed", :heldout, mixed_true; family_reports=mixed_reports), model)
        @test ismissing(mixed.predicted_types[2])
        @test all(!ismissing, mixed.predicted_types[[1, 3, 4]])
        @test occursin("lobe_2_family_conflict", mixed.flags)
    end

    @testset "explicit abstention modes" begin
        model = calibrate_type_temperature(_good_cases(); seed=0x17, config_hash="cfg-a",
            manifest_hash="man-a", source_hash="src-a")
        null_case = _case("null", :heldout, _good_cases()[3].true_types; control=:null,
            family_reports=_good_cases()[3].family_reports)
        conflict_case = _case("conflict", :heldout, _good_cases()[4].true_types; control=:family_conflict,
            family_reports=_good_cases()[4].family_reports)
        stale_case = _case("stale", :heldout, _good_cases()[1].true_types; source_hash="src-stale",
            family_reports=_good_cases()[1].family_reports)
        stale_config_case = _case("stale-config", :heldout, _good_cases()[1].true_types; config_hash="cfg-stale",
            family_reports=_good_cases()[1].family_reports)
        uncalibrated = predict_type(_good_cases()[1], nothing)
        @test all(ismissing, uncalibrated.predicted_types)
        @test all(ismissing, predict_type(null_case, model).predicted_types)
        @test all(ismissing, predict_type(conflict_case, model).predicted_types)
        @test all(ismissing, predict_type(stale_case, model).predicted_types)
        stale_cfg = predict_type(stale_config_case, model)
        @test all(ismissing, stale_cfg.predicted_types)
        @test occursin("stale_config", stale_cfg.flags)
    end

    @testset "explicit failure modes" begin
        @test_throws TypeCalibrationError calibrate_type_temperature(TypeCalibrationCase[
            _case("bad-seed", :calibration, [0, 1]; seed=0x99, family_reports=_good_cases()[1].family_reports)
        ]; seed=0x17, config_hash="cfg-a", manifest_hash="man-a", source_hash="src-a")

        @test_throws TypeCalibrationError calibrate_type_temperature(TypeCalibrationCase[
            _case("only-heldout", :heldout, [0, 1]; family_reports=_good_cases()[1].family_reports)
        ]; seed=0x17, config_hash="cfg-a", manifest_hash="man-a", source_hash="src-a")

        @test_throws TypeCalibrationError calibrate_type_temperature(TypeCalibrationCase[
            _case("only-cal", :calibration, [0, 1]; family_reports=_good_cases()[1].family_reports)
        ]; seed=0x17, config_hash="cfg-a", manifest_hash="man-a", source_hash="src-a")
    end
end
