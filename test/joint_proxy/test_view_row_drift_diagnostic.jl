#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "diagnose_joint_proxy_view_row_drift.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

@testset "joint-proxy view row-drift diagnostic" begin
    options = parse_row_drift_options([
        "--config", CONFIG, "--data-dir", "/tmp/data", "--file", "scan.sxm",
        "--out", "/tmp/rows.tsv", "--flatten-mode", "plane+rows",
        "--max-lag-px", "3", "--min-pairs", "6",
    ])
    @test options.flatten_mode == "plane+rows"
    @test options.max_lag_px == 3
    @test options.min_pairs == 6
    @test_throws ErrorException parse_row_drift_options([
        "--config", CONFIG, "--data-dir", "/tmp", "--file", "../scan.sxm",
        "--out", "x.tsv", "--flatten-mode", "none", "--max-lag-px", "1",
        "--min-pairs", "3",
    ])
    @test_throws ErrorException parse_row_drift_options([
        "--config", CONFIG, "--data-dir", "/tmp", "--file", "scan.sxm",
        "--out", "x.tsv", "--flatten-mode", "invalid", "--max-lag-px", "1",
        "--min-pairs", "3",
    ])

    fwd = [0.0, 2.0, -1.0, 3.0, 0.5, -2.0, 4.0, 1.0, -0.5, 2.5,
        -3.0, 0.7, 1.8, -1.4, 3.6, 0.2]
    bwd = fill(NaN, length(fwd))
    bwd[3:end] .= 1.4 .* fwd[1:end-2] .- 0.3
    metrics = row_lag_metrics(fwd, bwd, 3, 6)
    @test metrics.status == "ok"
    @test metrics.best_lag_px == 2
    @test metrics.best_corr ≈ 1.0 atol=1e-12
    @test metrics.slope ≈ 1.4 atol=1e-12
    @test metrics.affine_nrmse < 1e-12

    insufficient = row_lag_metrics(fwd[1:4], bwd[1:4], 1, 6)
    @test insufficient.status == "insufficient_pairs"

    mask = falses(5, 7)
    mask[2:4, 3:6] .= true
    @test roi_bounds(mask) == (2:4, 3:6)
end
