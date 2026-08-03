#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "diagnose_joint_proxy_view_preprocessing.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

@testset "joint-proxy view preprocessing diagnostic" begin
    options = parse_preprocessing_options([
        "--config", CONFIG,
        "--data-dir", "/tmp/data",
        "--file", "scan.sxm",
        "--out", "/tmp/view.tsv",
        "--flatten-modes", "none,plane,rows,plane+rows",
        "--max-shift-steps", "1",
    ])
    @test options.flatten_modes == ["none", "plane", "rows", "plane+rows"]
    @test options.max_shift_steps == 1
    @test_throws ErrorException parse_preprocessing_options([
        "--config", CONFIG, "--data-dir", "/tmp", "--file", "../scan.sxm",
        "--out", "x.tsv", "--flatten-modes", "none", "--max-shift-steps", "1",
    ])
    @test_throws ErrorException parse_preprocessing_options([
        "--config", CONFIG, "--data-dir", "/tmp", "--file", "scan.sxm",
        "--out", "x.tsv", "--flatten-modes", "none,none", "--max-shift-steps", "1",
    ])
    @test_throws ErrorException parse_preprocessing_options([
        "--config", CONFIG, "--data-dir", "/tmp", "--file", "scan.sxm",
        "--out", "x.tsv", "--flatten-modes", "invalid", "--max-shift-steps", "1",
    ])

    registry = load_registry(CONFIG)
    fwd = collect(1.0:registry.npix)
    bwd = 1.5 .* fwd .+ 0.2
    candidate = (; n=1, lobes=[(; index=1)])
    evidence = [TypePosteriorLobeEvidence(Dict("fwd" => fwd, "bwd" => bwd))]
    rows = preprocessing_transfer_rows("scan.sxm", "none", candidate, evidence,
        registry, 1)
    @test length(rows) == 1
    @test rows[1].flatten_mode == "none"
    @test rows[1].n == 1
    @test rows[1].lobe == 1
    @test rows[1].corr ≈ 1.0 atol=1e-12
    @test rows[1].slope ≈ 1.5 atol=1e-12
    @test rows[1].affine_nrmse < 1e-12

    base = Inference._default_pattern_config("scan.sxm", "/tmp")
    plane = pattern_with_flatten(base, "plane")
    @test plane.flatten == "plane"
    @test base.flatten == "plane+rows"
end
