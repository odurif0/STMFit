# Axial support threshold rules on synthetic axial profiles.
@testset "support threshold: half-maximum cap" begin
    G = GaussianFit2D
    t = collect(range(0.0, 8.0; length=4000))
    # A chain filling most of the axis: plateau 0.15 with a 0.10 modulation,
    # short ramps to background 0 at both ends.
    z = [0.15 * clamp(min(ti - 0.3, 7.7 - ti) / 0.3, 0, 1) * (0.8 + 0.2cos(2pi * ti / 0.9)) for ti in t]
    legacy = G.ChainSweepConfig(support_noise_k=2.5, support_padding_nm=0.0, support_min_length_nm=1.0)
    capped = deepcopy(legacy); capped.support_threshold_rule = "half_maximum_cap"
    lo0, hi0, m0 = G._active_t_support(t, z, legacy)
    lo1, hi1, m1 = G._active_t_support(t, z, capped; background=0.0)
    @test m1.threshold <= m0.threshold
    @test m1.threshold <= 0.5 * m1.peak + 1e-12
    @test hi1 - lo1 >= hi0 - lo0
    @test hi1 - lo1 > 6.5                      # the whole plateau
    @test_throws ErrorException G._active_t_support(t, z, capped)   # background required
    bad = deepcopy(capped); bad.support_threshold_fraction = 1.0
    @test_throws ErrorException G._active_t_support(t, z, bad; background=0.0)
    # Background below the legacy baseline only lowers the threshold; a legacy
    # threshold already under the cap is kept unchanged.
    short = [0.15 * exp(-0.5 * ((ti - 4.0) / 0.9)^2) for ti in t]
    _, _, ms0 = G._active_t_support(t, short, legacy)
    _, _, ms1 = G._active_t_support(t, short, capped; background=0.0)
    @test ms1.threshold == min(ms0.threshold, 0.5 * ms1.peak)
    oddrule = deepcopy(legacy); oddrule.support_threshold_rule = "quantile"
    @test_throws ErrorException G._active_t_support(t, z, oddrule)
    img = [1.0 2.0; 3.0 4.0]; mask = BitMatrix([true false; false true])
    @test G._off_roi_background(img, mask) == 2.5
    @test isnan(G._off_roi_background(img, trues(2, 2)))
end
