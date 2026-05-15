#!/usr/bin/env julia --project=.
using GaussianFit2D, Printf

function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    refined = GaussianFit2D.ChainModelResult[]
    isempty(results_circ) && return refined
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    ac_full = ctx_circ.axisctx_full
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ac_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    ccfg_refine = deepcopy(ccfg_ell)
    ccfg_refine.skip_global = true
    ccfg_refine.max_iter = 50
    ccfg_refine.multistart = 1
    for r_c in results_circ
        r_c.success || continue; n = r_c.n
        try
            n_prefix = 3; split_idx = n_prefix + n + 1 + (n - 1) + n
            p_init = vcat(r_c.params[1:split_idx],
                          r_c.params[(split_idx+1):end],
                          r_c.params[(split_idx+1):end])
            r_ref = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise,
                n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r_ref.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r_ref.params, n,
                    ac_fit, ccfg_refine; amp_min=r_ref.amp_min, amp_range=r_ref.amp_range)
                GaussianFit2D._finalize_chain_result!(r_ref, zfit, pred, noise,
                    n, n_eff, z, xs, ys, zimg, xfit, yfit, ac_fit, ccfg_refine)
                push!(refined, r_ref)
            end
        catch
        end
    end
    return refined
end

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
    selection_criterion="gcv",
)

results_circ, best_circ_raw, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx_circ)

println("N   GCV(circ)   GCV(ell)    RSS(circ)   RSS(ell)")
println(repeat("-", 65))
for n in sort(unique([r.n for r in results_circ if r.success]))
    rc = findfirst(r -> r.n == n && r.success, results_circ)
    re = findfirst(r -> r.n == n && r.success, results_ell)
    rc = rc === nothing ? nothing : results_circ[rc]
    re = re === nothing ? nothing : results_ell[re]
    println(rpad(string(n),4),
            rc !== nothing ? rpad(@sprintf("%.6f", rc.gcv), 12) : rpad("---",12),
            re !== nothing ? rpad(@sprintf("%.6f", re.gcv), 12) : rpad("---",12),
            rc !== nothing ? rpad(@sprintf("%.6f", rc.rss), 12) : rpad("---",12),
            re !== nothing ? rpad(@sprintf("%.6f", re.rss), 12) : rpad("---",12))
end
