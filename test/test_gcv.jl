#!/usr/bin/env julia --project=.
using GaussianFit2D, STMMolecularFit, Printf

fp = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"
img = GaussianFit2D.read_sxm(fp)
pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1)

ccfg = GaussianFit2D.ChainSweepConfig(
    n_min=2, n_max=14, spacing_min_nm=0.35, spacing_max_nm=0.75,
    fit_width_nm=0.15, support_threshold_fraction=0.25, support_noise_k=2.5,
    support_padding_nm=0.05, max_overlap=0.6,
    global_maxtime=10.0, global_maxiter=10000,
    sigma_parallel_min_nm=0.191, sigma_parallel_max_nm=0.509,
    sigma_perp_min_nm=0.10, sigma_perp_max_nm=0.55,
    intelligent_sweep=true, fuse_z_bwd=true,
    chain_tilted_baseline=true, cv_method="gcv",
    selection_criterion="gcv",  # ← test GCV selection
)

println("=== circular sweep with GCV selection ===")
results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
println()

# Print comparison: GCV vs BIC vs AICc vs CV
println(rpad("N",5), rpad("GCV",12), rpad("BIC",10), rpad("AICc",10), rpad("CV",10),
      rpad("RSS",10), rpad("valid",8))
println(repeat("-", 65))
for r in filter(r->r.success, results)
    println(rpad(string(r.n),5),
            rpad(@sprintf("%.6f", r.gcv),12),
            rpad(@sprintf("%.1f", r.bic),10),
            rpad(@sprintf("%.1f", r.aicc),10),
            rpad(@sprintf("%.4f", r.cv_nll_mean),10),
            rpad(@sprintf("%.6f", r.rss),10),
            rpad(r.valid ? "✓" : "✗",8))
end

println("\nGCV best N = $(best.n), valid=$(best.valid)")
println("GCV = $(round(best.gcv, digits=6))")
println("BIC = $(round(best.bic, digits=1)), AICc = $(round(best.aicc, digits=1))")

# Also compute what BIC selection would give
bic_best = argmin(r->r.bic, filter(r->r.success, results))
println("\nBIC best N = $(bic_best.n)")
