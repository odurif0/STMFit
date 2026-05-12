using Test
using STMFitCore

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

    # Known value: d = σ → ρ = exp(-0.5), κ = (1+exp(-0.5))/(1-exp(-0.5))
    kappa_known = (1.0 + exp(-0.5)) / (1.0 - exp(-0.5))
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
