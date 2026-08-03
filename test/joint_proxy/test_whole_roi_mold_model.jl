#!/usr/bin/env julia

using Test

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "proxy_registry.jl"))
include(joinpath(dirname(@__DIR__), "diagnose_joint_proxy_whole_roi.jl"))

const WR_REG = JointProxyRegistry

function synthetic_entry()
    type_pixels = Dict(
        0 => [0.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 0.0],
        1 => [0.0, 0.0, 0.0,
              0.0, 0.0, 1.0,
              0.0, 0.0, 0.0],
    )
    templates = WR_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        push!(templates, WR_REG.ProxyTemplate(typ, parity, mirror, copy(type_pixels[typ])))
    end
    source = WR_REG.ProxySource("synthetic", "whole_roi_test", "", "", 0.0, 1.0, true)
    return WR_REG.ProxyEntry(source, templates)
end

@testset "whole-ROI assembled mold diagnostic" begin
    entry = synthetic_entry()
    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    center = [(0.0, 0.0)]

    single = assemble_whole_roi(entry, [0], center, xs, ys;
        half_nm=0.2, step_nm=0.2, theta=0.0,
        direction=0, phase=0, mirror=0)
    doubled = assemble_whole_roi(entry, [0, 0], [(0.0, 0.0), (0.0, 0.0)], xs, ys;
        half_nm=0.2, step_nm=0.2, theta=0.0,
        direction=0, phase=0, mirror=0)
    @test doubled == 2 .* single
    @test single[findfirst(==(0.0), ys), findfirst(==(0.0), xs)] == 1.0
    @test all(single[:, firstindex(xs)] .== 0.0)

    centers = [(-0.2, 0.0), (0.2, 0.0), (0.6, 0.0)]
    xs_search = collect(-0.6:0.05:1.0)
    true_sequence = [0, 1, 0]
    true_mold = assemble_whole_roi(entry, true_sequence, centers, xs_search, ys;
        half_nm=0.2, step_nm=0.2, theta=0.0,
        direction=0, phase=0, mirror=0)
    backbone = [exp(-((x + 0.1)^2 + y^2) / 0.35) for y in ys, x in xs_search]
    observed = 0.7 .+ 1.4 .* backbone .+ 0.55 .* true_mold

    fit = fit_whole_roi_nuisance(observed, backbone, true_mold)
    @test fit.sse ≤ 1e-20
    @test fit.background ≈ 0.7 atol=1e-10
    @test fit.backbone_coefficient ≈ 1.4 atol=1e-10
    @test fit.mold_coefficient ≈ 0.55 atol=1e-10

    result = search_whole_roi_sequences(observed, backbone, entry, centers, xs_search, ys;
        half_nm=0.2, step_nm=0.2, theta=0.0,
        directions=(0,), phases=(0,), mirrors=(0,), max_n=8)
    @test result.sequence == true_sequence
    @test result.fit.sse ≤ 1e-20

    @test_throws ArgumentError assemble_whole_roi(entry, [2], center, xs, ys;
        half_nm=0.2, step_nm=0.2)
    @test_throws ArgumentError search_whole_roi_sequences(observed, backbone, entry,
        centers, xs_search, ys; half_nm=0.2, step_nm=0.2, max_n=2)
    @test_throws ArgumentError search_whole_roi_sequences(observed, backbone, entry,
        centers, xs_search, ys; half_nm=0.2, step_nm=0.2, max_n=13)

    root = normpath(joinpath(@__DIR__, "..", ".."))
    registry = WR_REG.load_registry(joinpath(root, "config", "joint_proxy_molds.toml"))
    dft_entry = only(item for item in registry.entries if item.source.family == "stm_dft_v1")
    tracked = assemble_whole_roi(dft_entry, [0, 1], [(-0.2, 0.0), (0.2, 0.0)],
        xs_search, ys; half_nm=registry.grid_half_nm, step_nm=registry.grid_step_nm,
        direction=0, phase=0, mirror=0)
    @test all(isfinite, tracked)
    @test maximum(tracked) > minimum(tracked)
end
