using Statistics, STMSXMIO

function fusion_image(fwd, bwd; forward_unit="nm", backward_unit="nm")
    SXMImage("synthetic_fusion.sxm", Dict{String,String}(), size(fwd, 2), size(fwd, 1),
        (4.0, 3.0), (0.0, 0.0),
        [SXMChannel("Z", forward_unit, "fwd", fwd),
         SXMChannel("Z", backward_unit, "bwd", bwd)])
end

@testset "Fused fitting treats both scan directions symmetrically" begin
    f = [exp(-((y-17)^2 + (x-21)^2)/60) + 0.15sin(1.7y + 0.8x)
         for y in 1:33, x in 1:41]
    b = [exp(-((y-17)^2 + (x-21)^2)/60) + 0.03cos(0.4y - 0.9x)
         for y in 1:33, x in 1:41]
    # Missing samples, different direction textures, borders and subsampling
    # exercise the native preprocessing path without a raw reader or optimizer.
    f[15, 17] = NaN
    b[19, 23] = Inf
    original_f, original_b = copy(f), copy(b)
    for flatten in ("none", "plane", "rows", "plane+rows"), stride in (1, 2), radius in (0, 1, 2)
        cfg = PatternConfig(flatten=flatten, stride=stride, smooth_radius_px=radius,
            no_plot=true)
        img = fusion_image(f, b)
        fused = GaussianFit2D._fused_roi_data(img, cfg)
        swapped = GaussianFit2D._fused_roi_data(fusion_image(b, f), cfg)
        @test fused == swapped
        @test all(isfinite, fused[3])
        @test isequal(f, original_f) && isequal(b, original_b)

        # Alternate-channel fitting already specifies an unsmoothed mean;
        # the two entry points must agree on the same physical samples.
        other = GaussianFit2D._channel_roi_data(img, cfg, fused[4], "Z")
        @test fused == other

        pf = preprocess_channel(img, get_channel(img, "Z"; direction="fwd"), cfg)
        pb = preprocess_channel(img, get_channel(img, "Z"; direction="bwd"), cfg)
        # The historical ROI and noise estimator are unchanged by this fix.
        _, _, old_mask = GaussianFit2D.molecule_roi_mask_fused(img, cfg, (pf[5] + pb[5])/2)
        count(old_mask) == 0 && (old_mask = trues(size(old_mask)))
        @test fused[4] == old_mask
        @test fused[8] == max(pf[7], pb[7], GaussianFit2D.EPS)
    end

    @testset "Identical views retain the unsmoothed signal" begin
        a = [exp(-((y-17)^2 + (x-21)^2)/60) + 0.1*(-1.0)^(x+y)
             for y in 1:33, x in 1:41]
        cfg = PatternConfig(flatten="none", stride=1, smooth_radius_px=1, no_plot=true)
        img = fusion_image(a, copy(a))
        fused = GaussianFit2D._fused_roi_data(img, cfg)
        @test fused[3] ≈ a .- quantile(a[fused[4]], 0.05)
        smooth = preprocess_channel(img, img.channels[1], cfg)[5]
        center(v) = v .- mean(v)
        @test maximum(abs, center(fused[3]) - center(smooth)) > 0.05
        # Unit conversion happens independently before either view is fused.
        meters = GaussianFit2D._fused_roi_data(
            fusion_image(a .* 1e-9, copy(a); forward_unit="m"), cfg)
        @test meters[3] ≈ fused[3]
        @test meters[4] == fused[4]
    end
end

@testset "Forward/backward alignment on the molecule (fwd_bwd_align)" begin
    # A chain of three lobes; the backward view is the forward one moved 6 px along x.
    base(x0) = [exp(-((y-25)^2)/40) * sum(exp(-((x - x0 - 12k)^2)/20) for k in 0:2) for y in 1:49, x in 1:81]
    f = base(25.0); b = base(31.0)
    cfg0 = PatternConfig(flatten="none", stride=1, smooth_radius_px=0, no_plot=true)
    cfg1 = PatternConfig(flatten="none", stride=1, smooth_radius_px=0, no_plot=true, fwd_bwd_align="molecule_x_shift")
    xs = collect(range(0, 4.0, length=81))
    @test isapprox(GaussianFit2D._fwd_bwd_molecule_lag(f, b, (f .+ b) ./ 2, xs), 6.0; atol=0.3)
    @test GaussianFit2D._shift_rows(f, 6.0)[:, 40:60] ≈ b[:, 40:60]
    plain = GaussianFit2D._fused_roi_data(fusion_image(f, b), cfg0)[3]
    aligned = GaussianFit2D._fused_roi_data(fusion_image(f, b), cfg1)[3]
    # Aligned fusion keeps the lobes sharp: higher peak-to-valley contrast along the chain.
    contrast(z) = (r = z[25, :]; maximum(r[25:55]) - minimum(r[30:46]))
    @test contrast(aligned) > 1.5 * contrast(plain)
    # Swapping the directions gives the same fused image, and identical views are unchanged.
    @test GaussianFit2D._fused_roi_data(fusion_image(b, f), cfg1)[3] ≈ aligned
    @test GaussianFit2D._fused_roi_data(fusion_image(f, copy(f)), cfg1)[3] ≈ GaussianFit2D._fused_roi_data(fusion_image(f, copy(f)), cfg0)[3]
    @test_throws ErrorException GaussianFit2D._fused_roi_data(fusion_image(f, b),
        PatternConfig(flatten="none", stride=1, smooth_radius_px=0, no_plot=true, fwd_bwd_align="other"))
end
