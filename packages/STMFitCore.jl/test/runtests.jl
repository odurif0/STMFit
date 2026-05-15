using Test
using STMFitCore
using Random

@testset "STMFitCore constraints" begin
    c = PhysicalChainConstraints(spacing_min_nm=0.35, spacing_max_nm=0.75,
                                 fwhm_min_nm=0.45, fwhm_max_nm=1.20,
                                 max_overlap=0.60)
    @test sigma_from_fwhm(1.20) ≈ 1.20 / FWHM_TO_SIGMA
    @test effective_spacing_min(c) ≈ 0.5150802014248211 atol=1e-12
    @test support_n_bounds(4.2, c; n_min_config=2, n_max_config=14) == (5, 9)
    @test chain_can_fit_support(9, 4.2, c)
    @test !chain_can_fit_support(10, 4.2, c)
    @test endpoint_overrun([0.0, 1.0, 2.0], 0.0, 2.0) == 0.0
    @test endpoint_overrun([-0.1, 1.0, 2.2], 0.0, 2.0) ≈ 0.3
end

@testset "Condition number κ" begin
    # Well-separated peaks: κ → 1
    @test overlap_condition_number(10.0, 0.2) ≈ 1.0 atol=0.01
    @test overlap_condition_number(10.0, 0.2) >= 1.0

    # Coincident: κ = Inf
    @test overlap_condition_number(0.0, 0.2) == Inf

    # Known value: d = σ → ρ = exp(-0.25), κ = (1+exp(-0.25))/(1-exp(-0.25))
    kappa_known = (1.0 + exp(-0.25)) / (1.0 - exp(-0.25))
    @test overlap_condition_number(0.5, 0.5) ≈ kappa_known atol=0.01

    # Penalty zero below threshold
    @test kappa_penalty(1.0; kappa_max=25.0) == 0.0
    @test kappa_penalty(25.0; kappa_max=25.0) == 0.0

    # Quadratic ramp: κ = 50, κ_max = 25 → (25/25)² = 1.0
    @test kappa_penalty(50.0; kappa_max=25.0, weight=1.0) ≈ 1.0 atol=1e-12

    # Disabled
    @test kappa_penalty(100.0; kappa_max=0.0) == 0.0

    # adjacent_kappa_max: empty
    @test adjacent_kappa_max(Float64[], Float64[]) == 1.0

    # Single pair: d=0.5, σ=(0.3+0.3)/2=0.3 → uses max(0.3,0.3)=0.3
    kappa_pair = overlap_condition_number(0.5, 0.3)
    @test adjacent_kappa_max([0.5], [0.3, 0.3]) ≈ kappa_pair atol=0.01
end

@testset "Residual diagnostics" begin
    # Durbin-Watson on pure noise ≈ 2
    rng = MersenneTwister(42)
    noise = randn(rng, 200)
    dw, dw_p = durbin_watson(noise)
    @test 1.5 < dw < 2.5
    @test dw_p > 0.01

    # DW on autocorrelated series << 2
    ar1 = cumsum(randn(rng, 200)) * 0.1
    dw_ar, _ = durbin_watson(ar1)
    @test dw_ar < 1.5

    # DW edge cases
    @test durbin_watson(Float64[]) === (NaN, NaN)
    @test durbin_watson([1.0, 2.0]) === (NaN, NaN)

    # Runs test on alternating signs → many runs
    alt = repeat([1.0, -1.0], 50)
    rn, re, rp = runs_test(alt)
    @test rn == 100
    @test rp < 0.05

    # Runs on uniform signs → 1 run
    @test runs_test(ones(10)) === (1, NaN, NaN)

    # Round-trip compute_residual_diagnostics
    rd = compute_residual_diagnostics(noise)
    @test isfinite(rd.durbin_watson)
    @test isfinite(rd.durbin_watson_p)
    @test isfinite(rd.residual_rms)
    @test isfinite(rd.residual_max)
    @test rd.runs_n > 0
    @test isfinite(rd.runs_expected)
    @test isfinite(rd.runs_p)
end
