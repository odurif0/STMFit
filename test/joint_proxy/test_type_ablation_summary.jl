#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "summarize_joint_proxy_type_ablation.jl"))

@testset "joint-proxy type-ablation summary" begin
    synthetic = [
        Dict("injected_family" => "stm_dft_v1", "ablation" => "stm_dft_v1",
            "noise_sigma" => "0", "correct" => "true", "true_margin" => "0.2"),
        Dict("injected_family" => "stm_dft_v1", "ablation" => "stm_dft_v1",
            "noise_sigma" => "0", "correct" => "true", "true_margin" => "0.1"),
    ]
    synth_metrics = synthetic_metrics(synthetic)
    recovery = only(filter(row -> row.metric == "recovery", synth_metrics))
    @test recovery.value == 1.0
    @test recovery.n == 2

    real = [
        Dict("patch_kind" => "residual", "ablation" => "combined", "p1" => "0.8",
            "argmax_type" => "1", "confidence" => "0.8"),
        Dict("patch_kind" => "residual", "ablation" => "combined", "p1" => "0.4",
            "argmax_type" => "0", "confidence" => "0.6"),
    ]
    summaries = [
        Dict("file" => "a.sxm", "status" => "ok", "patch_kind" => "residual",
            "ablation" => "combined", "stm_dft_v1_mass" => "0.75"),
        Dict("file" => "bad.sxm", "status" => "nonfinite_evidence", "patch_kind" => "residual",
            "ablation" => "combined", "stm_dft_v1_mass" => "NaN"),
    ]
    metrics = real_metrics(real, summaries)
    argmax1 = only(filter(row -> row.metric == "argmax1_fraction", metrics))
    family_mass_row = only(filter(row -> row.metric == "mean_stm_dft_v1_mass", metrics))
    @test argmax1.value == 0.5
    @test argmax1.n == 2
    @test family_mass_row.value == 0.75
    @test any(row -> row.metric == "status_ok" && row.value == 1, metrics)
    @test any(row -> row.metric == "status_nonfinite_evidence" && row.value == 1, metrics)

    geometric_only = [Dict("patch_kind" => "image", "ablation" => "geometric", "p1" => "0.5",
        "argmax_type" => "0", "confidence" => "0.5")]
    geometric_summary = [Dict("file" => "a.sxm", "status" => "ok", "patch_kind" => "image",
        "ablation" => "geometric", "stm_dft_v1_mass" => "0")]
    geometric_metrics = real_metrics(geometric_only, geometric_summary)
    @test !any(row -> row.metric == "mean_stm_dft_v1_mass", geometric_metrics)

    stability_rows = [Dict(
        "file" => "a.sxm", "n" => "2", "lobe" => "1", "patch_kind" => patch_kind,
        "ablation" => "combined_matched", "p1" => "0.5", "argmax_type" => argmax_type,
        "confidence" => "0.5",
    ) for (patch_kind, argmax_type) in zip(
        ("image", "raw", "residual", "neg_residual"), ("0", "1", "1", "0"))]
    stability_metrics = real_metrics(stability_rows, geometric_summary)
    image_residual = only(filter(row -> row.metric == "image_residual_agreement", stability_metrics))
    raw_residual = only(filter(row -> row.metric == "raw_residual_agreement", stability_metrics))
    @test image_residual.value == 0.0
    @test raw_residual.value == 1.0

    mktempdir() do dir
        path = joinpath(dir, "metrics.tsv")
        write_metrics(path, vcat(synth_metrics, metrics))
        @test length(readlines(path)) == length(synth_metrics) + length(metrics) + 1
    end
end
