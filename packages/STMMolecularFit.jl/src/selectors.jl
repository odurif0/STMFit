# ──────────────────────────────────────────────────────────────────────────────
# selectors.jl — Label-free count-selection helpers for STMFit chain models.
# Included by STMMolecularFit.jl; used by test/batch_full.jl.
#
# The robust-AICc guard rescoring an exhaustive elliptical candidate set, and
# the dispatcher that reports the primary count for the configured policy.
# All functions are module-internal (_ prefix).
# ──────────────────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ══════════════════════════════════════════════════════════════════════════════

function _load_chain_data(img, pcfg, ccfg)
    has_bwd = any(c -> lowercase(c.name) == lowercase(pcfg.roi_channel) && lowercase(c.direction) == "bwd", img.channels)
    return ccfg.fuse_z_bwd && has_bwd ? GaussianFit2D._fused_roi_data(img, pcfg) : GaussianFit2D._robust_roi_data(img, pcfg)
end

function _robust_rescore(r, xfit, yfit, zfit, noise, n_eff, axisctx, ccfg, nu::Real)
    r.success || return nothing
    try
        pred = GaussianFit2D._chain_model_values(
            xfit, yfit, r.params, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)
        resid = zfit .- pred
        total_nll = GaussianFit2D._student_nll(resid, noise, Float64(nu))
        pcount = GaussianFit2D._chain_nparams(r.n, ccfg)
        n_eff_safe = max(n_eff, pcount + 2)
        robust_aicc = 2 * total_nll + 2 * pcount +
            (2 * pcount * (pcount + 1)) / max(n_eff_safe - pcount - 1, 1)
        return (robust_aicc=robust_aicc, pcount=pcount)
    catch
        return nothing
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Robust AICc guard
# ══════════════════════════════════════════════════════════════════════════════

function _integrated_robust_aicc_n(img, pcfg, ccfg_ell; nu=8.0)
    # Robust guard uses an auxiliary exhaustive elliptical candidate set, matching
    # the validated robust-rescore audit. It is label-free, but intentionally
    # separate from the fast circ→ell effective selector because robust rescoring
    # the effective circ/refined set was empirically less stable.
    ccfg_guard = deepcopy(ccfg_ell)
    ccfg_guard.chain_circular_sigmas = false
    ccfg_guard.intelligent_sweep = false

    results_guard, _, _ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_guard)

    xs, ys, zimg, mask, x, y, z, noise = _load_chain_data(img, pcfg, ccfg_ell)
    axisctx_full = GaussianFit2D._weighted_roi_axis(x, y, z)

    xfit_ell, yfit_ell, zfit_ell, axisctx_ell, _, _ = GaussianFit2D._chain_fit_data(x, y, z, axisctx_full, ccfg_guard)
    n_eff_ell = max(10, length(zfit_ell) ÷ 9)

    best_n = 0
    best_score = Inf
    best_source = "NA"
    for r in results_guard
        r.success && r.valid || continue
        resc = _robust_rescore(r, xfit_ell, yfit_ell, zfit_ell, noise, n_eff_ell, axisctx_ell, ccfg_guard, nu)
        resc === nothing && continue
        if resc.robust_aicc < best_score
            best_n = r.n
            best_score = resc.robust_aicc
            best_source = "ell_robust_aicc"
        end
    end
    return best_n == 0 ? (nothing, "NA", NaN) : (best_n, best_source, best_score)
end

# ══════════════════════════════════════════════════════════════════════════════
# Primary selection dispatcher
# ══════════════════════════════════════════════════════════════════════════════

function _select_primary(n_eff::Int, eff_source::AbstractString, refined, policy::AbstractString)
    if policy == "gcv_with_robust_aicc_guard" && refined.robust_n != "NA"
        # Flag the selection as guard-driven for both down-only and
        # up-when-ambiguous moves (n_refined != n_eff).
        return refined.n_refined, policy, refined.n_refined != n_eff ? "robust_aicc_guard" : eff_source
    elseif policy == "adaptive_support_rescue" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    end
    return n_eff, "gcv", eff_source
end
