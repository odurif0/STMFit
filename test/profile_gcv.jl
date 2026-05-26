#!/usr/bin/env julia
# Compare GCV vs kfold CV timing for a single 2D sweep.
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_003.sxm"
const OUT = "results/profile_gcv"
mkpath(OUT)

const FWHM_SIGMA = 2.355
const SIGMA_MIN = 0.45 / FWHM_SIGMA
const SIGMA_MAX = 1.20 / FWHM_SIGMA

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUT, no_plot=true)

function make_ccfg(; cv_method="gcv", circular=true)
    GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
        spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
        support_noise_k=2.5, support_padding_nm=0.25,
        max_overlap=0.6, global_maxtime=10.0, global_maxiter=10000, cv_folds=5,
        cv_method=cv_method,
        sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
        sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
        intelligent_sweep=true, chain_circular_sigmas=circular, chain_tilted_baseline=true)
end

println("=== WARMUP ===")
t_warm = @elapsed begin
    img = GaussianFit2D.read_sxm(FILE)
    results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, make_ccfg(cv_method="gcv"))
end
@printf("  Warmup (GCV): %.1f s  N=%d\n", t_warm, best.n)

println("\n=== TIMED: GCV ===")
t_gcv = @elapsed begin
    results_gcv, best_gcv, ctx_gcv = GaussianFit2D.chain_gaussian_sweep(img, pcfg, make_ccfg(cv_method="gcv"))
end
@printf("  GCV sweep: %.2f s  N=%d  (%d results)\n", t_gcv, best_gcv.n, length(results_gcv))

println("\n=== TIMED: kfold (5-fold) ===")
t_kfold = @elapsed begin
    results_kf, best_kf, ctx_kf = GaussianFit2D.chain_gaussian_sweep(img, pcfg, make_ccfg(cv_method="kfold"))
end
@printf("  kfold sweep: %.2f s  N=%d  (%d results)\n", t_kfold, best_kf.n, length(results_kf))

@printf("\n  Speedup: %.1fx\n", t_kfold / max(t_gcv, 0.01))
@printf("  Time saved per file: %.1f s\n", t_kfold - t_gcv)

# Compare CV scores for each N
println("\n=== CV score comparison ===")
println("  N    GCV_nll    kfold_nll   delta")
for r_g in results_gcv
    r_g.success || continue
    n = r_g.n
    r_k = findfirst(r -> r.success && r.n == n, results_kf)
    if r_k !== nothing
        @printf("  %-3d  %9.2f  %9.2f  %+.2f\n", n, r_g.cv_nll_mean, results_kf[r_k].cv_nll_mean,
            results_kf[r_k].cv_nll_mean - r_g.cv_nll_mean)
    end
end
