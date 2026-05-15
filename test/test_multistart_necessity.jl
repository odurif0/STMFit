#!/usr/bin/env julia
# Check if multistart>1 finds a better minimum than start=1.
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_003.sxm"
const OUT = "results/multistart_test"
mkpath(OUT)

const FWHM_SIGMA = 2.355
const SIGMA_MIN = 0.45 / FWHM_SIGMA
const SIGMA_MAX = 1.20 / FWHM_SIGMA

pcfg = GaussianFit2D.PatternConfig(filepath=FILE, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUT, no_plot=true)

# Warmup
println("Warmup...")
img = GaussianFit2D.read_sxm(FILE)
ccfg_w = GaussianFit2D.ChainSweepConfig(n_min=6, n_max=6, multistart=1, cv_method="gcv",
    global_maxtime=10.0, global_maxiter=10000, intelligent_sweep=false,
    chain_circular_sigmas=true, chain_tilted_baseline=true)
GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_w)

println("\n=== N=6, circular, comparing multistart values ===")
for ms in [1, 3, 5, 20]
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=6, n_max=6, multistart=ms, cv_method="gcv",
        global_maxtime=10.0, global_maxiter=10000, intelligent_sweep=false,
        chain_circular_sigmas=true, chain_tilted_baseline=true)
    t = @elapsed begin
        results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    end
    @printf("  ms=%2d  BIC=%9.1f  RSS=%.6f  time=%.1fs\n", ms, best.bic, best.rss, t)
end
