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
