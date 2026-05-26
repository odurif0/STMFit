#!/usr/bin/env julia
# Show BIC for ALL N values on 034 (circular sweep, multistart=1)
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_034.sxm"
const OUT = "results/test_034_allN"
mkpath(OUT)

const FWHM_SIGMA = 2.355
const SIGMA_MIN = 0.45 / FWHM_SIGMA
const SIGMA_MAX = 1.20 / FWHM_SIGMA

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUT, no_plot=true)

ccfg = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14, multistart=1,
    spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
    support_noise_k=2.5, support_padding_nm=0.25,
    max_overlap=0.6, global_maxtime=10.0, global_maxiter=10000, cv_folds=5,
    cv_method="gcv",
    sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
    sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
    intelligent_sweep=false,  # force ALL N from 2 to 14
    chain_circular_sigmas=true, chain_tilted_baseline=true)

# Warmup
println("Warmup...")
img = GaussianFit2D.read_sxm(FILE)
ccfg_w = deepcopy(ccfg); ccfg_w.n_min=6; ccfg_w.n_max=6
GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_w)

println("\n=== Full sweep N=2..14 (circular) ===")
results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)

println("\n  N    BIC         GCV         valid  reason")
for r in sort(results; by=r->r.n)
    r.success || continue
    @printf("  %-2d  %10.1f  %11.6f  %-5s  %s\n", r.n, r.bic, r.gcv, r.valid, r.reason)
end
@printf("\n  Best by BIC: N=%d (BIC=%.1f)\n", best.n, best.bic)

# Now elliptical
ccfg_ell = deepcopy(ccfg)
ccfg_ell.chain_circular_sigmas = false
println("\n=== Full sweep N=2..14 (elliptical) ===")
results_e, best_e, ctx_e = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_ell)

println("\n  N    BIC         GCV         valid  reason")
for r in sort(results_e; by=r->r.n)
    r.success || continue
    @printf("  %-2d  %10.1f  %11.6f  %-5s  %s\n", r.n, r.bic, r.gcv, r.valid, r.reason)
end
@printf("\n  Best by BIC: N=%d (BIC=%.1f)\n", best_e.n, best_e.bic)
