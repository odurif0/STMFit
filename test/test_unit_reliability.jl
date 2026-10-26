# Synthetic tests of the benchmark-free reliability of per-unit signals.
#   julia --project=. test/test_unit_reliability.jl
using Test, Random, Statistics, Dates
using STMFit.UnitReliability

@testset "variance components recover known variances" begin
    rng = MersenneTwister(1)
    σb, σw = 0.3, 0.2
    # balanced: 400 lobes, 5 scans each
    groups = [(μ = σb * randn(rng); [μ + σw * randn(rng) for _ in 1:5]) for _ in 1:400]
    b, w, icc = variance_components(moments(groups))
    @test isapprox(sqrt(b), σb; rtol=0.06)
    @test isapprox(sqrt(w), σw; rtol=0.06)
    @test isapprox(icc, σb^2 / (σb^2 + σw^2); atol=0.03)
    # unbalanced: 3 to 12 scans per lobe
    groups = [(μ = σb * randn(rng); [μ + σw * randn(rng) for _ in 1:rand(rng, 3:12)]) for _ in 1:600]
    b, w, icc = variance_components(moments(groups))
    @test isapprox(sqrt(b), σb; rtol=0.06)
    @test isapprox(sqrt(w), σw; rtol=0.06)
    # no between-lobe variance: ICC near 0 (the estimate may be slightly negative)
    flat = [[σw * randn(rng) for _ in 1:5] for _ in 1:400]
    @test abs(variance_components(moments(flat))[3]) < 0.04
    # Moments add up like the pooled design
    a, c = groups[1:300], groups[301:end]
    @test all(isapprox.(variance_components(moments(a) + moments(c)), variance_components(moments(groups)); rtol=1e-9))
    # undefined designs
    @test all(isnan, variance_components(moments([[1.0, 2.0, 3.0]])))
    @test all(isnan, variance_components(Moments()))
end

@testset "bootstrap interval over tracks covers the truth and ratio detects a difference" begin
    rng = MersenneTwister(2)
    track(σb; lobes=6, scans=5) = moments([(μ = σb * randn(rng); [μ + 0.2 * randn(rng) for _ in 1:scans]) for _ in 1:lobes])
    a = [track(0.3) for _ in 1:60]
    B, W, I = bootstrap_components(a; draws=400, rng)
    @test isapprox(median(B), 0.3^2; rtol=0.25)                             # the bootstrap centres on the estimate
    @test quantile(B, 0.005) < 0.3^2 < quantile(B, 0.995)                   # a wide interval covers the truth
    @test quantile(I, 0.005) < 0.3^2 / (0.3^2 + 0.2^2) < quantile(I, 0.995)
    @test 0 < quantile(B, 0.975) - quantile(B, 0.025) < 0.1                 # and is not degenerate
    same = [track(0.3) for _ in 1:60]
    r, lo, hi = ratio_interval(a, same; draws=400, rng)
    @test lo < 1 < hi
    larger = [track(0.5) for _ in 1:60]
    r, lo, hi = ratio_interval(larger, a; draws=400, rng)
    @test lo > 1 && isapprox(r, (0.5 / 0.3)^2; rtol=0.35)
end

@testset "capped files are evenly spaced and keep the ends" begin
    @test capped_files(collect(1:5), 8) == collect(1:5)
    @test capped_files(collect(1:43), 8) == [1, 7, 13, 19, 25, 31, 37, 43]
    c = capped_files(collect(1:43), 8)
    @test length(c) == 8 && c[1] == 1 && c[end] == 43 && issorted(c)
end

# Synthetic run: tracks of scans of a 7-lobe chain with a known end profile, a scan gain,
# a persistent lobe effect and noise. Lobes of a scan whose count differs have no physical id.
function synthetic_run(; tracks=40, scans=6, σb=0.15, σw=0.1, seed=3)
    rng = MersenneTwister(seed)
    profile = Dict(1 => -0.2, 2 => 0.05, 3 => 0.08, 4 => 0.0)
    trk = Dict{String,Vector{String}}(); track_of = Dict{String,String}(); lobes = Dict{String,Vector{Lobe}}()
    n_scan = Dict{String,Int}(); n_final = Dict{String,Int}(); times = Dict{String,DateTime}()
    t0 = DateTime(2026, 1, 1)
    for t in 1:tracks
        effect = σb .* randn(rng, 7)
        for s in 1:scans
            f = "t$(t)_s$(s).sxm"; gain = 0.5 * randn(rng)
            push!(get!(trk, string(t), String[]), f); track_of[f] = string(t)
            odd = t == 1 && s == 1                       # a scan with another count, no physical lobes
            n = odd ? 5 : 7; n_scan[f] = n; n_final[f] = n
            times[f] = t0 + Minute(5 * s + 1000 * t)
            lobes[f] = [Lobe(k, min(k, n + 1 - k), odd ? "" : "track$(t)_$(k)",
                             exp(gain + profile[min(min(k, n + 1 - k), 4)] + (odd ? 0.0 : effect[k]) + σw * randn(rng)),
                             exp(-1.0), exp(-1.0)) for k in 1:n]
        end
    end
    return RunUnits("synthetic", trk, track_of, lobes, n_scan, n_final), times, profile
end

@testset "profile, residual groups and variance components on a synthetic run" begin
    ru, times, truth = synthetic_run()
    prof = distance_profile(ru, :amplitude)
    @test length(prof) == 4
    # The scan median removes a lobe-independent offset: compare differences to the d >= 4 class.
    @test isapprox(prof[1] - prof[4], truth[1] - truth[4]; atol=0.03)
    @test isapprox(prof[2] - prof[4], truth[2] - truth[4]; atol=0.03)
    tg = track_groups(ru, :amplitude; profile=prof)
    @test length(tg) == 40
    # d >= 2 keeps 5 lobes per scan; the unlabelled scan of track 1 is not used
    @test all(length(g) == 5 for (t, g) in tg if t != "1")
    @test all(length(x) == 6 for (t, g) in tg if t != "1" for x in g)
    @test all(length(x) == 5 for x in tg["1"])
    m = sum((moments(g) for g in values(tg)); init=Moments())
    b, w, icc = variance_components(m)
    @test isapprox(sqrt(b), 0.15; rtol=0.15)
    @test isapprox(sqrt(w), 0.1; rtol=0.1)
    # a cap below the number of scans lowers the observations per lobe
    tg3 = track_groups(ru, :amplitude; profile=prof, cap=3)
    @test all(length(x) == 3 for g in values(tg3) for x in g)
end

@testset "time-lag correlations separate a persistent pattern from noise" begin
    ru, times, _ = synthetic_run(σb=0.2, σw=0.05)
    ps = pair_correlations(ru, :amplitude, times)
    @test !isempty(ps) && all(p -> p[2] >= 0, ps)
    @test mean(p[3] for p in ps) > 0.8              # same persistent pattern in every scan
    ls = lag_summary(ps, [(0.0, 20.0), (20.0, Inf)]; draws=200, rng=MersenneTwister(4))
    @test length(ls) >= 1 && all(x -> x[4] <= x[3] <= x[5], ls)
    ru0, times0, _ = synthetic_run(σb=0.0, σw=0.1, seed=5)
    @test abs(mean(p[3] for p in pair_correlations(ru0, :amplitude, times0))) < 0.1
end

@testset "autocorrelation of profiles and its shuffle null" begin
    alt = [(-1.0)^k for k in 1:8]
    @test profile_autocorrelation(alt, 1) < -0.8
    @test profile_autocorrelation(alt, 2) > 0.5
    rng = MersenneTwister(6)
    P = [randn(rng, 7) for _ in 1:80]
    obs, null = autocorrelation_null(P, 1; draws=300, rng)
    @test quantile(null, 0.025) < obs < quantile(null, 0.975)          # random profiles sit inside the null
    obs2, null2 = autocorrelation_null([alt for _ in 1:20], 1; draws=300, rng)
    @test obs2 < quantile(null2, 0.025)                                 # an alternating profile does not
end

@testset "count dispersion" begin
    ru, _, _ = synthetic_run(tracks=4, scans=3)
    d = count_dispersion(ru)
    @test [x[1] for x in d] == ["1", "2", "3", "4"]
    @test [x[3] for x in d] == [2, 0, 0, 0]                              # track 1 has one scan with 5 lobes
end

@testset "normal distribution function" begin
    @test isapprox(normal_cdf(0), 0.5; atol=1e-9)
    @test isapprox(normal_cdf(1.0), 0.8413447; atol=2e-7)
    @test isapprox(normal_cdf(-1.96), 0.0249979; atol=2e-7)
    @test isapprox(normal_cdf(0.25) + normal_cdf(-0.25), 1; atol=1e-12)
    @test normal_cdf(8) > 0.999999 && normal_cdf(-8) < 1e-6
end

@testset "call robustness across models" begin
    @test call_robustness(0.0) == :dim && call_robustness(0.2) == :dim
    @test call_robustness(0.8) == :bright && call_robustness(1.0) == :bright
    @test call_robustness(0.5) == :unstable && call_robustness(0.21) == :unstable && call_robustness(0.79) == :unstable
    @test call_robustness(NaN) == :unknown
    @test call_robustness(0.7; threshold=0.7) == :bright && call_robustness(0.3; threshold=0.7) == :dim
    @test_throws ArgumentError call_robustness(0.5; threshold=0.5)
end
