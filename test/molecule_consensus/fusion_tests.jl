# Synthetic tests for latent-class molecule fusion.
#   julia --project=. test/molecule_consensus/fusion_tests.jl
using Test, Random, TOML, Statistics
include(joinpath(@__DIR__, "..", "lib", "molecule_fusion.jl"))
using .MoleculeFusion

const ST = MoleculeFusion.FusionSettings(3, 500, 1e-10, 0.5, 0.1, 0.9)

@testset "EM recovers detection and false-call rates" begin
    rng = MersenneTwister(7)
    z = rand(rng, 400) .< 0.35; n = rand(rng, 3:20, 400)
    k = [sum(rand(rng, n[i]) .< (z[i] ? 0.6 : 0.03)) for i in 1:400]
    fit = fit_latent_class(k, n, ST)
    @test abs(fit.theta1 - 0.6) < 0.05 && abs(fit.theta0 - 0.03) < 0.02 && abs(fit.pi - 0.35) < 0.07
    acc = mean((fit.posterior .>= 0.5) .== z)
    @test acc > 0.95
    # a minority of positive calls with rare false calls is evidence for type 1
    @test fit_latent_class(vcat(k, 3), vcat(n, 10), ST).posterior[end] > 0.9
    # start-independent within the declared neutral family
    fit2 = fit_latent_class(k, n, MoleculeFusion.FusionSettings(3, 500, 1e-10, 0.3, 0.2, 0.7))
    @test isapprox(fit.theta1, fit2.theta1; atol=1e-6) && isapprox(fit.pi, fit2.pi; atol=1e-6)
    @test_throws ErrorException fit_latent_class([1, 5], [2, 4], ST)
end


@testset "physical lobes are ranked along a fixed absolute axis" begin
    base = [(0.6i, 0.15 * (-1)^i) for i in 1:6]
    rot(p, a) = (cos(a) * p[1] - sin(a) * p[2], sin(a) * p[1] + cos(a) * p[2])
    order = ["a", "b", "c", "d"]
    tracks = Dict("a" => 1, "b" => 1, "c" => 1, "d" => 1)
    counts = Dict("a" => 6, "b" => 6, "c" => 6, "d" => 5)
    lobes = Dict("a" => base,
                 "b" => reverse([(p[1] + 1.3, p[2] - 0.7) for p in base]),     # drift, reversed lobe order
                 "c" => [rot((p[1] + 0.2, p[2]), 0.3) for p in base],           # small rotation
                 "d" => base[1:5])                                              # different count: not fused
    k = physical_lobes(order, tracks, counts, lobes)
    @test [k[("a", l)] for l in 1:6] == ["track1_$r" for r in 1:6]
    @test [k[("b", l)] for l in 1:6] == ["track1_$r" for r in 6:-1:1]
    @test [k[("c", l)] for l in 1:6] == ["track1_$r" for r in 1:6]
    @test !haskey(k, ("d", 1))
end

@testset "settings" begin
    repo = dirname(dirname(@__DIR__))
    st = load_fusion_settings(TOML.parsefile(joinpath(repo, "config", "molecule_consensus.toml")))
    @test st.min_scans == 3 && st.init_theta0 < st.init_theta1
    bad = TOML.parsefile(joinpath(repo, "config", "molecule_consensus.toml")); bad["selection"]["fusion"] = "majority"
    @test_throws ErrorException load_fusion_settings(bad)
end
