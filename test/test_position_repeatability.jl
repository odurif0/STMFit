# Synthetic tests of the alignment used by test/report_position_repeatability.jl.
#   julia --project=. test/test_position_repeatability.jl
using Test, Random, Statistics
include(joinpath(@__DIR__, "report_position_repeatability.jl"))

@testset "rigid and Procrustes alignment" begin
    base = [(0.66 * (i - 1), 0.1 * (-1)^i) for i in 1:6]
    rotate(p, a, t) = (cos(a) * p[1] - sin(a) * p[2] + t[1], sin(a) * p[1] + cos(a) * p[2] + t[2])
    angles = Dict("a" => 0.0, "b" => 0.3, "c" => -0.45, "d" => 0.1)
    shifts = Dict("a" => (0.0, 0.0), "b" => (1.2, -0.4), "c" => (-0.7, 2.0), "d" => (0.05, 0.02))
    P = Dict(f => [rotate(p, angles[f], shifts[f]) for p in base] for f in keys(angles))
    scans = sort(collect(keys(angles)))
    Q, rot = align_procrustes(P, scans, 6)
    m = lobe_means(Q, scans, 6)
    @test maximum(hypot(Q[f][r][1] - m[r][1], Q[f][r][2] - m[r][2]) for f in scans, r in 1:6) < 1e-9
    for f in scans
        @test isapprox(rot[f] - rot["a"], -angles[f]; atol=1e-9)   # rotation that undoes each scan's
    end
    # translation-only shapes: rigid alignment removes them exactly
    T = Dict(f => [(p[1] + shifts[f][1], p[2] + shifts[f][2]) for p in base] for f in scans)
    R = align_translation(T, scans, 6); mt = lobe_means(R, scans, 6)
    @test maximum(hypot(R[f][r][1] - mt[r][1], R[f][r][2] - mt[r][2]) for f in scans, r in 1:6) < 1e-12
    # noise of known size is recovered by the pooled SD
    rng = MersenneTwister(1); σ = 0.07
    N = Dict(f => [(p[1] + σ * randn(rng), p[2] + σ * randn(rng)) for p in base] for f in ["s$i" for i in 1:400])
    sc = sort(collect(keys(N))); mn = lobe_means(N, sc, 6)
    ss = sum((N[f][r][1] - mn[r][1])^2 for f in sc, r in 1:6)
    @test isapprox(pooled_sd(ss, 6 * (length(sc) - 1)), σ; rtol=0.05)
    @test isnan(pooled_sd(1.0, 0))
end
