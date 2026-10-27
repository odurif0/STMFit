# Label-free typing of units in two families (campaign T) on synthetic data.
#   julia --project=. test/test_unit_typing.jl
using Test, Random, Statistics, LinearAlgebra
using STMFit.UnitTyping

@testset "mixtures and BIC separate one family from two" begin
    rng = MersenneTwister(1)
    two = vcat(randn(rng, 120, 2) .+ [3.0 0.0], randn(rng, 80, 2) .- [3.0 0.0])
    one = randn(rng, 200, 2)
    for (X, expect_two) in ((two, true), (one, false))
        M1 = fit_mixture(X, 1); M2 = fit_mixture(X, 2)
        d = mixture_bic(M1, size(X, 1)) - mixture_bic(M2, size(X, 1))
        @test (d >= 6) == expect_two
    end
    M2 = fit_mixture(two, 2)
    @test sort(M2.weights) ≈ [0.4, 0.6] atol=0.05
    P = posterior(M2, two)
    @test all(isapprox.(sum(P; dims=2), 1))
    lab = P[:, argmax(M2.means[:, 1])] .> 0.5
    @test mean(lab[1:120]) > 0.97 && mean(lab[121:end]) < 0.03
end

@testset "reproducible subspace keeps what repeats" begin
    rng = MersenneTwister(2)
    n = 300; trait = randn(rng, n)
    A = hcat(trait .+ 0.3randn(rng, n), randn(rng, n), randn(rng, n))     # only feature 1 carries a unit trait
    B = hcat(trait .+ 0.3randn(rng, n), randn(rng, n), randn(rng, n))
    V, vals = reproducible_subspace(A, B)
    @test size(V, 2) == 1 && abs(V[1, 1]) > 0.99
    @test vals[1] > 0.8 && all(abs.(vals[2:end]) .< 0.2)
    @test_throws ArgumentError reproducible_subspace(A, B[1:10, :])
end

@testset "kappa and robust z" begin
    @test cohen_kappa([true, true, false, false], [true, true, false, false]) == 1.0
    @test isapprox(cohen_kappa([true, false, true, false], [true, true, false, false]), 0.0; atol=1e-12)
    z = robust_z([1.0, 2.0, 3.0, 4.0, 100.0, NaN])
    @test z[3] == 0.0 && isnan(z[end]) && z[5] > 10
    @test robust_z([2.0, 2.0, 2.0]) == [0.0, 0.0, 0.0]
end

@testset "contrastive direction and winsorising" begin
    rng = MersenneTwister(3)
    # control: shared variation on feature 1; application: the same plus extra variation on feature 2
    Cc = [1.0 0.0 0.0; 0.0 0.1 0.0; 0.0 0.0 0.5]
    Ca = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.5]
    v, λ = contrastive_direction(Ca, Cc)
    @test abs(v[2]) > 0.99 && λ > 5
    @test minimum(eigvals(psd_shrink([1.0 0.0; 0.0 -1.0]))) > 0
    w = winsorize(vcat(collect(1.0:100.0), 1000.0))
    @test maximum(w) < 100 && minimum(w) > 1
end

@testset "two calibrated types against a continuum" begin
    rng = MersenneTwister(4)
    n = 300; τ = 0.5; noise = fill(0.05, n)
    types = rand(rng, n) .< 0.4
    s2 = ifelse.(types, 1.5, -1.0) .+ τ .* randn(rng, n) .+ sqrt.(noise) .* randn(rng, n)
    s1 = 1.2 .* randn(rng, n) .+ sqrt.(noise) .* randn(rng, n)              # same spread, no types
    for (s, expect) in ((s2, true), (s1, false))
        v = τ^2 .+ noise
        T = fit_two_types(s, v); C = fit_continuum(s, noise)
        dbic = (-2C.loglik + 2log(n)) - (-2T.loglik + 3log(n))
        @test (dbic >= 6 && min(T.π, 1 - T.π) >= 0.1) == expect
    end
    T = fit_two_types(s2, τ^2 .+ noise)
    hi = T.μ1 > T.μ2 ? T.π : 1 - T.π
    @test isapprox(hi, 0.4; atol=0.07)
    p = type_posterior(s2, τ^2 .+ noise, T.π, T.μ1, T.μ2)
    lab = T.μ1 > T.μ2 ? p .> 0.5 : p .< 0.5
    @test mean(lab .== types) > 0.95
    C = fit_continuum(1.2 .* randn(rng, 2000), fill(0.0, 2000))
    @test isapprox(C.τ2, 1.44; rtol=0.1)
end
