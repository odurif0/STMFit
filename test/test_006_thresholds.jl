#!/usr/bin/env julia
# Show noise-based support detection behavior for 006.
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_006.sxm"
const SIGMA_MIN = 0.45 / 2.355
const SIGMA_MAX = 1.20 / 2.355

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir="/tmp", no_plot=true)
img = GaussianFit2D.read_sxm(FILE)

for (label, noise_k, pad) in [("inspect-style", 2.5, 0.20), ("tighter-noise", 3.0, 0.05)]
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=6, n_max=6, multistart=1, cv_method="gcv",
        intelligent_sweep=false, chain_circular_sigmas=true, chain_tilted_baseline=true,
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        support_noise_k=noise_k, support_padding_nm=pad,
        sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
        sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
        max_overlap=0.6, global_maxtime=10.0)
    
    results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    sm = ctx.support_meta
    
    @printf("\n=== %s (noise_k=%.2f, pad=%.2f) ===\n", label, noise_k, pad)
    @printf("  threshold_noise=%.4f  effective=%.4f\n",
        sm.threshold_noise, sm.threshold)
    @printf("  baseline=%.4f  peak=%.4f  noise=%.4f\n", sm.baseline, sm.peak, sm.noise_sigma_profile)
    @printf("  Support: %.2f nm  padding: %.2f nm\n",
        get(sm, :final_support_length_nm, NaN), sm.padding_nm)
end
