#!/usr/bin/env julia --project=.
using GaussianFit2D, Printf

fp = "/home/durif/Rebecca/data/data/20240817_LHe_Cu100/240817_002.sxm"
img = GaussianFit2D.read_sxm(fp)
pcfg = GaussianFit2D.PatternConfig(filepath=fp, channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1)

ccfg = GaussianFit2D.ChainSweepConfig(
    n_min=2, n_max=14, spacing_min_nm=0.35, spacing_max_nm=0.75,
    fit_width_nm=0.15, support_noise_k=2.5,
    support_padding_nm=0.25, max_overlap=0.6,
    global_maxtime=10.0, global_maxiter=10000,
    sigma_parallel_min_nm=0.191, sigma_parallel_max_nm=0.509,
    sigma_perp_min_nm=0.10, sigma_perp_max_nm=0.55,
    intelligent_sweep=true, fuse_z_bwd=true,
    chain_tilted_baseline=true, cv_method="gcv",
)

results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)

r_c = results_circ[1]  # N=6
n = r_c.n
println("Circ params length: $(length(r_c.params)), n=$(n)")
println("Param layout: b0,bx,by,amps($n),t0,deltas($(n-1)),us($n),sigs($n)")

# Debug expansion
xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
ac_full = ctx_circ.axisctx_full
xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ac_full, ccfg)

n_prefix = 3; split_idx = n_prefix + n + 1 + (n - 1) + n
println("split_idx = $split_idx (expected for n=$n with tilted: 3 + $n + 1 + $(n-1) + $n = $(3+n+1+(n-1)+n))")
p_init = vcat(r_c.params[1:split_idx],
              r_c.params[(split_idx+1):end],
              r_c.params[(split_idx+1):end])
println("p_init length: $(length(p_init)) (expected: $(3 + n + 1 + (n-1) + n + n + n))")

ccfg_refine = deepcopy(ccfg)
ccfg_refine.skip_global = true
ccfg_refine.max_iter = 50
ccfg_refine.multistart = 1

try
    r_ref = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise,
        n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
    println("Refinement: success=$(r_ref.success), RSS=$(round(r_ref.rss,digits=6))")
catch e
    println("Refinement FAILED: $e")
end
