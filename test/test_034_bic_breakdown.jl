#!/usr/bin/env julia
# BIC breakdown: n_eff, NLL, penalty for 003 vs 034.
using GaussianFit2D, Printf

const DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"
const FWHM_SIGMA = 2.355
const SIGMA_MIN = 0.45 / FWHM_SIGMA
const SIGMA_MAX = 1.20 / FWHM_SIGMA

for fn in ["240817_003.sxm", "240817_034.sxm"]
    fp = joinpath(DIR, fn)
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=6, n_max=6, multistart=1, cv_method="gcv",
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        support_threshold_fraction=0.25, support_noise_k=2.5, support_padding_nm=0.05,
        max_overlap=0.6, global_maxtime=10.0,
        sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
        sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
        intelligent_sweep=false, chain_circular_sigmas=true, chain_tilted_baseline=true)
    
    img = GaussianFit2D.read_sxm(fp)
    results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    
    # Get fit data
    xs, ys, zimg, mask = ctx.xs, ctx.ys, ctx.zimg, ctx.mask
    n_roi = count(mask)
    n_eff = max(10, n_roi ÷ 9)
    pcount = GaussianFit2D._chain_nparams(6, ccfg)
    
    @printf("\n=== %s ===\n", fn)
    @printf("  ROI pixels: %d  n_eff: %d\n", n_roi, n_eff)
    @printf("  noise: %.5f\n", ctx.noise)
    @printf("  BIC: %.1f = 2×NLL(%.1f) + %d×log(%d) = %.1f + %.1f\n",
        best.bic, best.train_nll * n_roi, pcount, n_eff,
        2 * best.train_nll * n_roi, pcount * log(n_eff))
    @printf("  train_nll/point: %.5f\n", best.train_nll)
    @printf("  RSS: %.4f  chi2_red: %.2f\n", best.rss, best.chi2_reduced)
end
