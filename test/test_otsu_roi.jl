#!/usr/bin/env julia
# Test Otsu-adaptive ROI on 003, 034, 035.
using GaussianFit2D, Printf

const DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"
const SIGMA_MIN = 0.45 / 2.355
const SIGMA_MAX = 1.20 / 2.355

# Warmup
println("Warmup...")
fp0 = joinpath(DIR, "240817_003.sxm")
pcfg0 = GaussianFit2D.PatternConfig(filepath=fp0, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
img0 = GaussianFit2D.read_sxm(fp0)
GaussianFit2D.chain_gaussian_sweep(img0, pcfg0, GaussianFit2D.ChainSweepConfig(
    n_min=6, n_max=6, multistart=1, cv_method="gcv", intelligent_sweep=false,
    chain_circular_sigmas=true, chain_tilted_baseline=true,
    spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
    sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
    sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
    support_threshold_fraction=0.25, support_noise_k=2.5, support_padding_nm=0.05,
    max_overlap=0.6, global_maxtime=10.0))

for fn in ["240817_003.sxm", "240817_034.sxm", "240817_035.sxm"]
    fp = joinpath(DIR, fn)
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
    img = GaussianFit2D.read_sxm(fp)
    
    xs, ys, roi_mask = GaussianFit2D.molecule_roi_mask(img, pcfg)
    
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14, multistart=1, cv_method="gcv",
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        support_threshold_fraction=0.25, support_noise_k=2.5, support_padding_nm=0.05,
        max_overlap=0.6, global_maxtime=10.0,
        sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
        sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
        intelligent_sweep=false, chain_circular_sigmas=true, chain_tilted_baseline=true)
    
    results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    
    n_roi = count(roi_mask)
    n_fit = count(ctx.mask)
    n_img = prod(size(ctx.zimg))
    
    println("\n=== $fn ===")
    @printf("  ROI: %d px (%.1f%% of image), fit_mask: %d px\n", n_roi, n_roi/n_img*100, n_fit)
    
    println("  N    BIC          nll/pt     valid")
    for r in sort(results; by=r->r.n)
        r.success || continue
        @printf("  %-2d  %11.1f  %.5f    %s\n", r.n, r.bic, r.train_nll, r.valid)
    end
    @printf("  → Best: N=%d  BIC=%.1f\n", best.n, best.bic)
end
