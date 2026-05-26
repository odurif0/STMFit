#!/usr/bin/env julia --project=.
using GaussianFit2D, Printf

fp = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"
img = GaussianFit2D.read_sxm(fp)
pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1)

# GCV sweep
ccfg_gcv = GaussianFit2D.ChainSweepConfig(
    n_min=2, n_max=14, spacing_min_nm=0.35, spacing_max_nm=0.75,
    fit_width_nm=0.15, support_noise_k=2.5,
    support_padding_nm=0.25, max_overlap=0.6,
    global_maxtime=10.0, global_maxiter=10000,
    sigma_parallel_min_nm=0.191, sigma_parallel_max_nm=0.509,
    sigma_perp_min_nm=0.10, sigma_perp_max_nm=0.55,
    intelligent_sweep=true, fuse_z_bwd=true,
    chain_tilted_baseline=true, cv_method="gcv",
    selection_criterion="gcv",
)

# k-fold sweep
ccfg_kfold = deepcopy(ccfg_gcv)
ccfg_kfold.cv_method = "kfold"
ccfg_kfold.cv_folds = 3

t_gcv = @elapsed results_gcv, best_gcv, ctx_gcv = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_gcv)
t_kfold = @elapsed results_kfold, best_kfold, ctx_kfold = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_kfold)

println("GCV:    N=$(best_gcv.n)  time=$(round(t_gcv,digits=1))s")
println("k-fold: N=$(best_kfold.n)  time=$(round(t_kfold,digits=1))s")
println("Speedup: $(round(t_kfold/t_gcv,digits=1))×")
println()

println("GCV values:")
for r in filter(r->r.success, results_gcv)
    println("  N=$(r.n): GCV=$(round(r.gcv,digits=6))  BIC=$(round(r.bic,digits=1))  RSS=$(round(r.rss,digits=6))")
end
