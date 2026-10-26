# Joint shape fit of repeated scans (campaign J, stage 2) on synthetic data.
#   julia --project=. test/test_joint_shape.jl
using Test, Random, Statistics
using STMFit.JointShape

@testset "rigid transforms and molecule frame" begin
    T = (0.4, 1.5, -2.0)
    src = [(0.0, 0.0), (1.0, 0.2), (2.1, -0.3), (0.5, 1.1)]
    dst = [apply_rigid(T, p) for p in src]
    θ, tx, ty = rigid_fit(src, dst)
    @test isapprox(θ, 0.4; atol=1e-10) && isapprox(tx, 1.5; atol=1e-10) && isapprox(ty, -2.0; atol=1e-10)
    back = [apply_rigid(invert_rigid(T), p) for p in dst]
    @test all(isapprox(collect(a), collect(b); atol=1e-10) for (a, b) in zip(back, src))
    pts = [(1.0 + 0.6k * cos(0.3), 2.0 + 0.6k * sin(0.3)) for k in 0:5]
    F = molecule_frame(pts)
    tu = [to_frame(F, p) for p in pts]
    @test issorted(first.(tu)) && all(abs(x[2]) < 1e-9 for x in tu)
    @test all(isapprox(collect(from_frame(F, x...)), collect(p); atol=1e-9) for (x, p) in zip(tu, pts))
    @test_throws ArgumentError rigid_fit([(0.0, 0.0)], [(1.0, 1.0)])
end

@testset "stack of normalised scans" begin
    ts = collect(-1.0:0.02:1.0); us = collect(-0.8:0.02:0.8)
    f(gain, off) = (t, u) -> off + gain * exp(-(t^2 + u^2) / (2 * 0.3^2))
    img, cnt, used = stack_scans([f(2.0, 5.0), f(0.5, -1.0), (t, u) -> t > 0.9 ? NaN : f(1.0, 0.0)(t, u)], ts, us; chain_t=(-0.5, 0.5))
    @test used == 3
    i0, j0 = argmin(abs.(us)), argmin(abs.(ts))
    @test img[i0, j0] > 0.9                                    # gains and offsets are removed before the mean
    @test maximum(cnt) == 3 && minimum(cnt) == 2               # the third scan misses t > 0.9
end

@testset "joint fit recovers known units" begin
    rng = MersenneTwister(7)
    ts = collect(-2.2:0.02:2.2); us = collect(-0.8:0.02:0.8)
    truth = [LobeFit(-1.5 + 0.62k, 0.03 * (-1)^k, 1.0 + 0.3 * (k == 2), k == 3 ? 0.36 : 0.26, k == 1 ? 0.22 : 0.30, k == 4 ? 0.3 : 0.0) for k in 0:5]
    T = [t for u in us for t in ts]; U = [u for u in us for t in ts]
    clean = lobe_model(T, U, truth, (0.1, 0.0, 0.02))
    img = permutedims(reshape(clean .+ 0.01 .* randn(rng, length(clean)), length(ts), length(us)))
    init = [(l.t + 0.08 * randn(rng), 0.0) for l in truth]
    lobes, bg, info = fit_lobes(ts, us, img, init)
    @test info.converged && info.rms_rel < 0.02
    @test all(abs(a.t - b.t) < 0.02 && abs(a.u - b.u) < 0.02 for (a, b) in zip(lobes, truth))
    @test all(abs(a.sigma_par - b.sigma_par) < 0.02 && abs(a.sigma_perp - b.sigma_perp) < 0.02 for (a, b) in zip(lobes, truth))
    @test abs(lobes[5].theta - 0.3) < 0.1
    @test isapprox(lobes[3].amplitude / lobes[1].amplitude, 1.3; rtol=0.05)
end
