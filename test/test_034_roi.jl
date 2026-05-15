#!/usr/bin/env julia
# Compare ROI detection between 003 and 034.
using GaussianFit2D, Printf, Statistics

const DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"

for fn in ["240817_003.sxm", "240817_034.sxm"]
    fp = joinpath(DIR, fn)
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
    
    img = GaussianFit2D.read_sxm(fp)
    xs, ys, raw, z, z_smooth, unit, noise = GaussianFit2D.preprocess_channel(img,
        GaussianFit2D.get_channel(img, "Z"; direction="fwd"), pcfg)
    
    # ROI detection
    mask, roi_info = GaussianFit2D.molecule_roi_mask(img, pcfg)
    
    signal = z_smooth .- minimum(z_smooth)
    threshold = pcfg.roi_threshold_fraction * maximum(signal)
    raw_mask = signal .> threshold
    
    println("\n=== $fn ===")
    @printf("  Image: %d × %d px\n", size(z_smooth)...)
    @printf("  Signal range: [%.4f, %.4f] nm\n", minimum(signal), maximum(signal))
    @printf("  roi_threshold_fraction=%.2f → threshold=%.4f nm\n", pcfg.roi_threshold_fraction, threshold)
    @printf("  Raw mask (before dilation): %d px (%.1f%%)\n", count(raw_mask), count(raw_mask)/length(raw_mask)*100)
    @printf("  Final mask (after dilation=%d px): %d px (%.1f%%)\n", pcfg.roi_dilate_px, count(mask), count(mask)/length(mask)*100)
    
    # Chain fit
    results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, GaussianFit2D.ChainSweepConfig(
        n_min=6, n_max=6, multistart=1, cv_method="gcv", intelligent_sweep=false,
        chain_circular_sigmas=true, chain_tilted_baseline=true,
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        sigma_parallel_min_nm=0.45/2.355, sigma_parallel_max_nm=1.20/2.355,
        sigma_perp_min_nm=0.45/2.355, sigma_perp_max_nm=1.20/2.355,
        support_threshold_fraction=0.25, support_noise_k=2.5, support_padding_nm=0.05,
        max_overlap=0.6, global_maxtime=10.0))
    
    fit_mask = ctx.mask
    @printf("  Chain fit mask: %d px (%.1f%% of ROI)\n", count(fit_mask), count(fit_mask)/count(mask)*100)
end
