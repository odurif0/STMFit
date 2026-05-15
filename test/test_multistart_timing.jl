#!/usr/bin/env julia
# Time a single N fit with varying multistart.
using GaussianFit2D, Printf

const FILE = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_003.sxm"
const OUT = "results/multistart_timing"
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

for ms in [1, 5, 10, 20]
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=6, n_max=6, multistart=ms, cv_method="gcv",
        global_maxtime=10.0, global_maxiter=10000, intelligent_sweep=false,
        chain_circular_sigmas=true, chain_tilted_baseline=true)
    t = @elapsed begin
        results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
    end
    @printf("  multistart=%2d  time=%.1f s  per_start=%.1f s\n", ms, t, t/ms)
end
