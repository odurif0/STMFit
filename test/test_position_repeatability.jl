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

@testset "per-lobe uncertainty source" begin
    proc(t, r, pos, ss) = Dict{String,Any}("track" => t, "rank" => r, "N" => 3, "scans" => 3, "mode" => "procrustes",
        "position" => pos, "ss_along" => ss, "ss_across" => 4ss, "dof" => 2.0)
    rows = [proc(1, 1, "end", 0.02), proc(1, 2, "interior", 0.008), proc(1, 3, "end", 0.02)]
    geo = Dict((f, l) => Dict("x_nm" => "1", "y_nm" => "2") for f in ("a", "b", "c", "lone"), l in 1:3)
    R = (rows=rows, used=Dict(1 => ["a", "b", "c"]), geo=geo,
         phys=Dict((f, l) => "track1_$l" for f in ("a", "b", "c"), l in 1:3),
         tracks=Dict("a" => 1, "b" => 1, "c" => 1, "lone" => 2), counts=Dict(f => 3 for f in ("a", "b", "c", "lone")),
         lobes_abs=Dict(f => [(0.0, 0.0), (0.66, 0.0), (1.32, 0.0)] for f in ("a", "b", "c", "lone")))
    out = position_rows(R)
    @test length(out) == 12
    a2 = only(r for r in out if r["file"] == "a" && r["lobe"] == "2")
    @test a2["sd_source"] == "track" && a2["sd_along_nm"] == @sprintf("%.4f", sqrt(0.008 / 2))
    lone1 = only(r for r in out if r["file"] == "lone" && r["lobe"] == "1")
    @test lone1["sd_source"] == "cohort" && lone1["track_scans"] == "1"
    @test lone1["sd_along_nm"] == @sprintf("%.4f", sqrt(0.04 / 4))   # pooled over the two tracked end lobes
    R0 = merge(R, (rows=Dict{String,Any}[],))
    @test all(r["sd_source"] == "none" && r["sd_along_nm"] == "NA" for r in position_rows(R0) if r["file"] == "lone")
end
