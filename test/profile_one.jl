#!/usr/bin/env julia
# Profile a single file to identify bottlenecks.
using STMMolecularFit, GaussianFit2D, GaussianFit1D
using Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_003.sxm"
const OUT = "results/profile_timing"
mkpath(OUT)

const FWHM_SIGMA = 2.355
const SIGMA_MIN = 0.45 / FWHM_SIGMA
const SIGMA_MAX = 1.20 / FWHM_SIGMA

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUT, no_plot=true)

ccfg_circ = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
    spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
    support_threshold_fraction=0.25, support_noise_k=2.5, support_padding_nm=0.05,
    max_overlap=0.6, global_maxtime=10.0, global_maxiter=10000, cv_folds=5,
    sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
    sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
    intelligent_sweep=true, chain_circular_sigmas=true, chain_tilted_baseline=true)

println("=== WARMUP (compilation) ===")
t_warm = @elapsed begin
    img = GaussianFit2D.read_sxm(FILE)
end
@printf("  read_sxm: %.2f s\n", t_warm)

t_warm = @elapsed begin
    GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
end
@printf("  2D sweep (warmup): %.2f s\n", t_warm)

println("\n=== TIMED RUNS ===")

# 1D
t_1d = @elapsed begin
    scfg = STMMolecularFit.SlideConfig(width_nm=0.30, support_threshold_fraction=0.20,
        support_noise_k=2.5, support_padding_nm=0.20, output_dir=OUT, no_plot=true)
    fcfg = STMMolecularFit.FitSlideConfig(min_spacing=0.35, max_spacing=0.75, max_overlap=0.6, output_dir=OUT)
    img1d = STMMolecularFit.read_sxm(FILE)
    slide = STMMolecularFit.extract_slide(img1d, scfg)
    fit_1d = STMMolecularFit.fit_slide(slide, fcfg)
end
@printf("  1D total:         %.2f s\n", t_1d)

# 2D circular
t_circ = @elapsed begin
    results_circ, best_circ, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
end
@printf("  2D circ sweep:    %.2f s  (%d N-values)\n", t_circ, length(results_circ))

# 2D elliptical (full NLopt sweep)
ccfg_ell = deepcopy(ccfg_circ)
ccfg_ell.chain_circular_sigmas = false
t_ell = @elapsed begin
    results_ell, best_ell, ctx_ell = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_ell)
end
@printf("  2D ell sweep:     %.2f s  (%d N-values)\n", t_ell, length(results_ell))

# circ→ell refinement (LsqFit only, no NLopt)
t_ref = @elapsed begin
    # reuse _refine_circ_to_ell logic inline
    ccfg_refine = deepcopy(ccfg_ell)
    ccfg_refine.skip_global = true
    ccfg_refine.max_iter = 50
    ccfg_refine.multistart = 1
end
@printf("  circ→ell refine:  %.2f s\n", t_ref)

@printf("\n  TOTAL per file:   %.2f s  (warmup excluded)\n", t_1d + t_circ + t_ell + t_ref)
@printf("  If NLopt=10s×4N×2sweeps = ~80s just in NLopt\n")
