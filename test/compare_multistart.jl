#!/usr/bin/env julia --project=.
# Compare multistart=1 vs multistart=20 on key files
using GaussianFit2D, Printf

FILES = ["240817_002.sxm", "240817_026.sxm", "240817_041.sxm"]
DATA_DIR = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100"
SIGMA_MIN = 0.191; SIGMA_MAX = 0.509

for fn in FILES
    fp = joinpath(DATA_DIR, fn)
    img = GaussianFit2D.read_sxm(fp)
    pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1)

    for ms in [1, 20]
        ccfg = GaussianFit2D.ChainSweepConfig(
            n_min=2, n_max=14, spacing_min_nm=0.35, spacing_max_nm=0.75,
            fit_width_nm=0.15,
            support_noise_k=2.5, support_padding_nm=0.25,
            max_overlap=0.6, global_maxtime=10.0, global_maxiter=10000,
            sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
            sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
            intelligent_sweep=true, fuse_z_bwd=true,
            chain_tilted_baseline=true, cv_method="gcv",
            selection_criterion="gcv", multistart=ms,
        )
        ccfg_circ = deepcopy(ccfg); ccfg_circ.chain_circular_sigmas = true
        t0 = time()
        results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
        dt = time() - t0
        println("$(fn) ms=$ms: Ncirc=$(best.n) GCV=$(round(best.gcv,digits=6)) BIC=$(round(best.bic,digits=1)) time=$(round(dt,digits=0))s")
    end
    println()
end
