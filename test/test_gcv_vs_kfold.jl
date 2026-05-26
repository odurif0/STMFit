#!/usr/bin/env julia
# Compare GCV vs kfold: scores per N + timing.
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_003.sxm"
const OUT = "results/gcv_vs_kfold"
mkpath(OUT)

const FWHM_SIGMA = 2.355
const SIGMA_MIN = 0.45 / FWHM_SIGMA
const SIGMA_MAX = 1.20 / FWHM_SIGMA

function make_ccfg(; cv_method="gcv")
    GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        support_noise_k=2.5, support_padding_nm=0.25,
        max_overlap=0.6, global_maxtime=10.0, global_maxiter=10000, cv_folds=5,
        cv_method=cv_method,
        sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
        sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
        intelligent_sweep=true, chain_circular_sigmas=true, chain_tilted_baseline=true)
end

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUT, no_plot=true)

# ── Warmup ──
println("=== WARMUP ===")
@time begin
    img = GaussianFit2D.read_sxm(FILE)
    GaussianFit2D.chain_gaussian_sweep(img, pcfg, make_ccfg(cv_method="gcv"))
end

# ── GCV sweep ──
println("\n=== GCV sweep ===")
t_gcv = @elapsed begin
    results_gcv, best_gcv, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, make_ccfg(cv_method="gcv"))
end
@printf("  Time: %.2f s  (best N=%d)\n", t_gcv, best_gcv.n)

# ── kfold sweep ──
println("\n=== kfold sweep ===")
t_kf = @elapsed begin
    results_kf, best_kf, _ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, make_ccfg(cv_method="kfold"))
end
@printf("  Time: %.2f s  (best N=%d)\n", t_kf, best_kf.n)

# ── Comparison ──
@printf("\n=== Speedup: %.1fx  (saved %.0f s) ===\n", t_kf / max(t_gcv, 0.01), t_kf - t_gcv)
println("\nPer-N comparison:")
println("  N   BIC        GCV_score   kfold_nll   GCV_sel   kfold_sel")
gcv_best_n = argmin([r.success ? r.gcv : Inf for r in results_gcv])
kf_best_n  = argmin([r.success ? r.cv_nll_mean : Inf for r in results_kf])
bic_best_n = argmin([r.success ? r.bic : Inf for r in results_gcv])

for (i, (rg, rk)) in enumerate(zip(results_gcv, results_kf))
    rg.success || continue
    n = rg.n
    gcv_sel = rg.gcv == minimum(r.gcv for r in results_gcv if r.success) ? " <--" : ""
    kf_sel  = rk.success && rk.cv_nll_mean == minimum(r.cv_nll_mean for r in results_kf if r.success) ? " <--" : ""
    @printf("  %-2d  %9.1f  %11.4f  %11.4f  %s  %s\n",
        n, rg.bic, rg.gcv, rk.success ? rk.cv_nll_mean : NaN, gcv_sel, kf_sel)
end

# Also check: does the N=0 result exist?
has_n0_gcv = any(r -> r.success && r.n == 0, results_gcv)
has_n0_kf  = any(r -> r.success && r.n == 0, results_kf)
@printf("\n  N=0 in results: GCV=%s, kfold=%s\n", has_n0_gcv, has_n0_kf)
