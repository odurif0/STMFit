using Test
using GaussianFit2D

@testset "split-width chain profile" begin
    axisctx = (origin=(0.0, 0.0), axis=[1.0, 0.0], perp=[0.0, 1.0], tmin=-1.0, tmax=1.0)
    ccfg_gauss = GaussianFit2D.ChainSweepConfig(
        n_min=1, n_max=1,
        sigma_parallel_min_nm=0.10, sigma_parallel_max_nm=0.50,
        sigma_perp_min_nm=0.10, sigma_perp_max_nm=0.50,
        chain_tilted_baseline=false,
        peak_profile=:gaussian,
    )
    ccfg_split = deepcopy(ccfg_gauss)
    ccfg_split.peak_profile = :split
    ccfg_split.skew_ratio_max = 2.0

    sigma_raw = GaussianFit2D._rlogit((0.20 - 0.10) / (0.50 - 0.10))
    p_gauss = [0.0, log(2.0), 0.0, 0.0, sigma_raw, sigma_raw]
    p_split = vcat(p_gauss, 0.0) # _rsigmoid(0)=0.5 => skew_ratio=1

    x = collect(range(-0.6, 0.6; length=31))
    y = collect(range(-0.3, 0.3; length=31))
    z_gauss = GaussianFit2D._chain_model_values(x, y, p_gauss, 1, axisctx, ccfg_gauss)
    z_split = GaussianFit2D._chain_model_values(x, y, p_split, 1, axisctx, ccfg_split)
    @test z_split ≈ z_gauss atol=1e-12 rtol=1e-12

    _b, feats, _ts, _us, _spars, _sperps = GaussianFit2D._decode_chain(p_split, 1, axisctx, ccfg_split)
    @test length(feats) == 1
    @test feats[1].skew_ratio ≈ 1.0 atol=1e-12
end
