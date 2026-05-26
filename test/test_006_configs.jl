#!/usr/bin/env julia
# Compare batch vs inspect configs on 006.
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_006.sxm"
const OUT = "/tmp/test_006_cfg"
const SIGMA_MIN = 0.45 / 2.355
const SIGMA_MAX = 1.20 / 2.355

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUT, no_plot=true)

# Warmup
println("Warmup...")
img = GaussianFit2D.read_sxm(FILE)
GaussianFit2D.chain_gaussian_sweep(img, pcfg, GaussianFit2D.ChainSweepConfig(
    n_min=6, n_max=6, multistart=1, cv_method="gcv", intelligent_sweep=false,
    chain_circular_sigmas=true, chain_tilted_baseline=true,
    spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
    sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
    sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
    support_noise_k=2.5, support_padding_nm=0.20,
    max_overlap=0.6, global_maxtime=10.0))

function test_config(label, img, pcfg; support_noise_k=2.5, support_padding_nm=0.20, chain_circular_sigmas=true)
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14, multistart=1, cv_method="gcv",
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        support_noise_k=support_noise_k,
        support_padding_nm=support_padding_nm,
        max_overlap=0.6, global_maxtime=10.0, global_maxiter=10000,
        sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
        sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
        intelligent_sweep=true, chain_circular_sigmas=chain_circular_sigmas,
        chain_tilted_baseline=true, fuse_z_bwd=true)
    
    results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    
    supp_len = hasproperty(ctx, :support_meta) ? get(ctx.support_meta, :final_support_length_nm, NaN) : NaN
    n_roi = count(ctx.mask)
    
    println("\n=== $label ===")
    @printf("  support_noise_k=%.2f  supp_pad=%.2f  circ=%s\n",
        support_noise_k, support_padding_nm, chain_circular_sigmas)
    @printf("  ROI: %d px, support: %.2f nm\n", n_roi, isnan(supp_len) ? ctx.axisctx.tmax - ctx.axisctx.tmin : supp_len)
    
    for r in sort(results; by=r->r.n)
        r.success || continue
        marker = r.n == best.n ? " <-- BEST" : ""
        @printf("  N=%-2d  BIC=%10.1f  nll/pt=%.5f  valid=%s%s\n",
            r.n, r.bic, r.train_nll, r.valid, marker)
    end
end

# inspect config
test_config("inspect-style (noise_k=2.5, pad=0.20)", img, pcfg;
    support_noise_k=2.5, support_padding_nm=0.20)

# batch config
test_config("legacy batch padding (noise_k=2.5, pad=0.05)", img, pcfg;
    support_noise_k=2.5, support_padding_nm=0.05)

# Isolate: only change padding
test_config("noise_k=3.0, pad=0.20", img, pcfg;
    support_noise_k=3.0, support_padding_nm=0.20)

# Isolate: only change noise multiplier
test_config("noise_k=2.0, pad=0.05", img, pcfg;
    support_noise_k=2.0, support_padding_nm=0.05)
