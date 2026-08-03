#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "diagnose_joint_proxy_view_registration.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

@testset "joint-proxy view registration diagnostic" begin
    options = parse_registration_options([
        "--config", CONFIG, "--data-dir", "/tmp/data", "--file", "scan.sxm",
        "--out", "/tmp/registered.tsv", "--flatten-mode", "plane+rows",
        "--max-lag-px", "8", "--min-pairs", "20",
    ])
    @test options.max_lag_px == 8
    @test options.min_pairs == 20

    z = reshape(collect(1.0:18.0), 3, 6)
    view_data = Main.JointProxyCandidateViews.ViewData("bwd", collect(1.0:6.0),
        collect(1.0:3.0), z, 0.1, true)
    shifted = shift_view_x(view_data, -2)
    @test all(isnan, shifted.z[:, 1:2])
    @test shifted.z[:, 3:6] == z[:, 1:4]
    @test view_data.z == z

    lag_rows = [(; status="ok", best_lag_px=-4), (; status="ok", best_lag_px=-3),
        (; status="ok", best_lag_px=-4), (; status="insufficient_pairs", best_lag_px=0)]
    @test robust_global_lag(lag_rows) == -4
    @test_throws ErrorException robust_global_lag([(; status="insufficient_pairs", best_lag_px=0)])
end
