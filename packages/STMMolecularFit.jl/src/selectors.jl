# ──────────────────────────────────────────────────────────────────────────────
# selectors.jl — Experimental model-selection policies for STMFit chain models.
# Included by STMMolecularFit.jl.
#
# All functions are module-internal (_ prefix).  They depend on GaussianFit2D
# (already `using`'d in the parent module) and on each other.
# ──────────────────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════

const SUPPORT_MARG_PAD_GRID = [0.00, 0.10, 0.20, 0.25, 0.35, 0.50]
const SUPPORT_MARG_REGRET_MARGIN = 0.02
const STABILITY_COMPETITIVE_TOL = 0.01
const LOBE_EVIDENCE_SEP_SIGMA = 2.0
const LOBE_EVIDENCE_VALLEY_SNR = 3.0
const LOBE_EVIDENCE_VALLEY_FRAC = 0.2
const LAPLACE_FD_REL_STEP = 1e-4
const LAPLACE_FD_MIN_STEP = 1e-6
const LAPLACE_EIGEN_CAP = 1e12

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
# Spatial blocked cross-validation
# ══════════════════════════════════════════════════════════════════════════════

function _blocked_cv_score(xs, ys, zimg, x, y, z, noise, n::Int, axisctx,
                           ccfg::GaussianFit2D.ChainSweepConfig, folds::Int;
                           warm_start=nothing)
    npoints = length(x)
    npoints > folds || return (Inf, 0.0, 0)

    t = (x .- axisctx.origin[1]) .* axisctx.axis[1] .+ (y .- axisctx.origin[2]) .* axisctx.axis[2]
    sort_idx = sortperm(t)

    block_ends = Int[0]
    base = npoints ÷ folds
    rem = npoints % folds
    for k in 1:folds
        push!(block_ends, block_ends[end] + base + (k <= rem ? 1 : 0))
    end

    ccfg_cv = deepcopy(ccfg)
    ccfg_cv.skip_global = true
    ccfg_cv.max_iter = min(max(50, ccfg.max_iter ÷ 2), 150)
    ccfg_cv.multistart = 1

    scores = Float64[]
    all_idx = collect(1:npoints)
    for k in 1:folds
        val_range = (block_ends[k] + 1):block_ends[k + 1]
        isempty(val_range) && continue
        val_idx = sort_idx[val_range]
        train_idx = setdiff(all_idx, val_idx)
        length(train_idx) > 10 && length(val_idx) > 5 || continue

        r = GaussianFit2D._fit_chain_n(xs, ys, zimg,
            x[train_idx], y[train_idx], z[train_idx], noise,
            n, axisctx, ccfg_cv; starts=1, warm_start=warm_start)
        r.success || continue

        pred = GaussianFit2D._chain_model_values(
            x[val_idx], y[val_idx], r.params, n, axisctx, ccfg_cv;
            amp_min=r.amp_min, amp_range=r.amp_range)
        nll = GaussianFit2D._student_nll(z[val_idx] .- pred, noise, ccfg_cv.student_nu)
        push!(scores, nll / length(val_idx))
    end

    length(scores) >= 2 || return (Inf, 0.0, length(scores))
    return (mean(scores), std(scores) / sqrt(length(scores)), length(scores))
end

function _spatial_blocked_cv_selection(img, pcfg, ccfg_ell, ctx_circ, results_ell; folds=5)
    xs, ys, zimg, _, x, y, z, noise = _load_chain_data(img, pcfg, ccfg_ell)
    xfit, yfit, zfit, axisctx, _, _ = GaussianFit2D._chain_fit_data(
        x, y, z, ctx_circ.axisctx_full, ccfg_ell)

    by_n = Dict{Int,Any}()
    for r in results_ell
        r.success && r.valid || continue
        if !haskey(by_n, r.n) || r.gcv < by_n[r.n].gcv
            by_n[r.n] = r
        end
    end

    best_n = 0
    best_score = Inf
    best_ok = 0
    best_se = 0.0
    for n in sort(collect(keys(by_n)))
        r = by_n[n]
        score, se, ok = _blocked_cv_score(xs, ys, zimg, xfit, yfit, zfit, noise,
            n, axisctx, ccfg_ell, folds; warm_start=r.params)
        ok >= 2 && isfinite(score) || continue
        if score < best_score
            best_n = n
            best_score = score
            best_ok = ok
            best_se = se
        end
    end

    return best_n == 0 ? (nothing, NaN, 0, NaN) : (best_n, best_score, best_ok, best_se)
end

# ══════════════════════════════════════════════════════════════════════════════
# Support-marginalized GCV
# ══════════════════════════════════════════════════════════════════════════════

function _result_by_n_prefer_ell(results_ell, results_circ)
    # Track provenance explicitly (whether the current entry for each N came
    # from the elliptical set) instead of inferring it via value-equality `in`,
    # which silently misfires when a circ result and an ell result are
    # field-equal (e.g. two degenerate fits with gcv == Inf): the circ entry
    # would then look "in results_ell" and wrongly survive the ell preference.
    by_n = Dict{Int,Any}()
    from_ell = Dict{Int,Bool}()
    for r in results_circ
        r.success && r.valid || continue
        if !haskey(by_n, r.n) || r.gcv < by_n[r.n].gcv
            by_n[r.n] = r
            from_ell[r.n] = false
        end
    end
    for r in results_ell
        r.success && r.valid || continue
        if !haskey(by_n, r.n) || r.gcv < by_n[r.n].gcv || !get(from_ell, r.n, false)
            by_n[r.n] = r
            from_ell[r.n] = true
        end
    end
    return by_n
end

function _support_rescore_gcv(x_roi, y_roi, z_roi, noise, axisctx_full, ccfg_p, r)
    xfit, yfit, zfit, axisctx_p, _, _ = GaussianFit2D._chain_fit_data(x_roi, y_roi, z_roi, axisctx_full, ccfg_p)
    nd = length(zfit)
    pcount = GaussianFit2D._chain_nparams(r.n, ccfg_p)
    nd > pcount || return Inf
    GaussianFit2D._chain_can_fit_support(r.n, axisctx_p, ccfg_p) || return Inf
    pred = GaussianFit2D._chain_model_values(xfit, yfit, r.params, r.n, axisctx_p, ccfg_p;
        amp_min=r.amp_min, amp_range=r.amp_range)
    rss = sum(abs2, zfit .- pred)
    return nd / (nd - pcount)^2 * rss
end

function _support_marginalized_gcv_selection(img, pcfg, ccfg_ell, results_ell, results_circ;
                                             pad_grid=SUPPORT_MARG_PAD_GRID)
    _, _, _, _, x_roi, y_roi, z_roi, noise = _load_chain_data(img, pcfg, ccfg_ell)
    axisctx_full = GaussianFit2D._weighted_roi_axis(x_roi, y_roi, z_roi)
    by_n = _result_by_n_prefer_ell(results_ell, results_circ)
    isempty(by_n) && return (nothing, NaN, 0, 0, NaN, Dict{Int,Float64}(), Dict{Int,Float64}())

    regrets = Dict{Int,Vector{Float64}}()
    supports_ok = 0
    for pad in pad_grid
        ccfg_p = deepcopy(ccfg_ell)
        ccfg_p.support_padding_nm = Float64(pad)
        scores = Tuple{Int,Float64}[]
        for n in sort(collect(keys(by_n)))
            gcv = try
                _support_rescore_gcv(x_roi, y_roi, z_roi, noise, axisctx_full, ccfg_p, by_n[n])
            catch
                Inf
            end
            isfinite(gcv) && push!(scores, (n, gcv))
        end
        isempty(scores) && continue
        best = minimum(last, scores)
        isfinite(best) || continue
        for (n, s) in scores
            push!(get!(regrets, n, Float64[]), (s - best) / max(abs(best), eps(Float64)))
        end
        supports_ok += 1
    end

    candidates = [n for (n, rs) in regrets if length(rs) >= 3]
    isempty(candidates) && return (nothing, NaN, supports_ok, 0, NaN, Dict{Int,Float64}(), Dict{Int,Float64}())
    med = Dict(n => median(regrets[n]) for n in candidates)
    q75 = Dict(n => quantile(regrets[n], 0.75) for n in candidates)
    ordered = sort(candidates; by=n -> (med[n], q75[n], n))
    best_n = first(ordered)
    runner_delta = length(ordered) >= 2 ? med[ordered[2]] - med[best_n] : NaN
    return (best_n, med[best_n], supports_ok, length(candidates), runner_delta, med, q75)
end

# ══════════════════════════════════════════════════════════════════════════════
# Stability selection
# ══════════════════════════════════════════════════════════════════════════════

function _stability_selection(img, pcfg, ccfg_ell, results_ell, results_circ;
                              pad_grid=SUPPORT_MARG_PAD_GRID,
                              competitive_tol=STABILITY_COMPETITIVE_TOL)
    _, _, _, _, x_roi, y_roi, z_roi, noise = _load_chain_data(img, pcfg, ccfg_ell)
    axisctx_full = GaussianFit2D._weighted_roi_axis(x_roi, y_roi, z_roi)
    by_n = _result_by_n_prefer_ell(results_ell, results_circ)
    isempty(by_n) && return nothing

    competitive = Dict{Int,Int}()
    feasible = Dict{Int,Int}()
    supports_ok = 0
    for pad in pad_grid
        ccfg_p = deepcopy(ccfg_ell)
        ccfg_p.support_padding_nm = Float64(pad)
        scores = Tuple{Int,Float64}[]
        for n in sort(collect(keys(by_n)))
            gcv = try
                _support_rescore_gcv(x_roi, y_roi, z_roi, noise, axisctx_full, ccfg_p, by_n[n])
            catch
                Inf
            end
            isfinite(gcv) && push!(scores, (n, gcv))
        end
        isempty(scores) && continue
        best_gcv = minimum(last, scores)
        isfinite(best_gcv) || continue
        supports_ok += 1
        for (n, gcv) in scores
            feasible[n] = get(feasible, n, 0) + 1
            if gcv <= best_gcv * (1.0 + competitive_tol)
                competitive[n] = get(competitive, n, 0) + 1
            end
        end
    end

    supports_ok >= 3 || return nothing
    candidates = [n for n in keys(feasible) if feasible[n] >= 3]
    isempty(candidates) && return nothing
    ordered = sort(candidates; by=n -> (-get(competitive, n, 0), -get(feasible, n, 0), n))
    best_n = first(ordered)
    win_pct = 100 * get(competitive, best_n, 0) / supports_ok
    return (best_n, win_pct, supports_ok, get(competitive, best_n, 0), get(feasible, best_n, 0))
end

# ══════════════════════════════════════════════════════════════════════════════
# Slope-heuristic MDL
# ══════════════════════════════════════════════════════════════════════════════

function _theil_sen_slope(xs::Vector{Float64}, ys::Vector{Float64})
    slopes = Float64[]
    for i in eachindex(xs), j in (i + 1):length(xs)
        dx = xs[j] - xs[i]
        abs(dx) <= eps(Float64) && continue
        push!(slopes, (ys[j] - ys[i]) / dx)
    end
    isempty(slopes) && return NaN
    return median(slopes)
end

function _slope_heuristic_mdl_selection(results_ell, ccfg_ell, results_circ, ccfg_circ, n_raw::Int)
    n_eff = max(10, n_raw ÷ 9)
    candidates = NamedTuple[]
    for r in results_ell
        r.success && r.valid && isfinite(r.rss) && r.rss > 0 || continue
        d = GaussianFit2D._chain_nparams(r.n, ccfg_ell)
        contrast = n_eff * log(max(r.rss, eps(Float64)) / max(n_raw, 1))
        push!(candidates, (n=r.n, source="ell", d=Float64(d), contrast=contrast, rss=r.rss, result=r))
    end
    for r in results_circ
        r.success && r.valid && isfinite(r.rss) && r.rss > 0 || continue
        d = GaussianFit2D._chain_nparams(r.n, ccfg_circ)
        contrast = n_eff * log(max(r.rss, eps(Float64)) / max(n_raw, 1))
        push!(candidates, (n=r.n, source="circ", d=Float64(d), contrast=contrast, rss=r.rss, result=r))
    end
    length(candidates) >= 4 || return (nothing, NaN, 0, NaN)

    dmed = median([c.d for c in candidates])
    high = [c for c in candidates if c.d >= dmed]
    length(high) >= 4 || return (nothing, NaN, length(candidates), NaN)
    slope = _theil_sen_slope([c.d for c in high], [c.contrast for c in high])
    isfinite(slope) || return (nothing, NaN, length(candidates), NaN)
    alpha = max(0.0, -slope)
    alpha > 0 || return (nothing, NaN, length(candidates), alpha)

    scored = [(c=c, score=c.contrast + 2.0 * alpha * c.d) for c in candidates]
    sort!(scored; by=x -> (x.score, x.c.d, x.c.rss, x.c.n))
    best = first(scored)
    return (best.c.n, best.score, length(candidates), alpha)
end

# ══════════════════════════════════════════════════════════════════════════════
# Laplace evidence (Fisher information approximation)
# ══════════════════════════════════════════════════════════════════════════════

function _chain_param_bounds_like_fit(n::Int, ccfg)
    np = GaussianFit2D._chain_nparams(n, ccfg)
    lower = fill(-10.0, np)
    upper = fill(10.0, np)
    lower[1] = -5.0; upper[1] = 5.0
    j = 2
    if ccfg.chain_tilted_baseline
        lower[j] = -1.0; upper[j] = 1.0; j += 1
        lower[j] = -1.0; upper[j] = 1.0; j += 1
    end
    for _ in 1:n
        lower[j] = -5.0; upper[j] = 5.0; j += 1
    end
    n_spacing = GaussianFit2D._chain_spacing_param_count(n, ccfg)
    lower[j] = -4.0; upper[j] = 4.0; j += 1
    for _ in 1:(n_spacing - 1)
        lower[j] = -5.0; upper[j] = 5.0; j += 1
    end
    for _ in 1:n
        lower[j] = -3.0; upper[j] = 3.0; j += 1
    end
    n_sigma_types = GaussianFit2D._chain_sigma_param_count(n, ccfg)
    if ccfg.chain_circular_sigmas
        for _ in 1:n_sigma_types
            lower[j] = -5.0; upper[j] = 5.0; j += 1
        end
    else
        for _ in 1:(2n_sigma_types)
            lower[j] = -5.0; upper[j] = 5.0; j += 1
        end
    end
    return lower, upper
end

function _laplace_fd_jacobian(r, xfit, yfit, axisctx, ccfg, lower, upper)
    p = collect(Float64, r.params)
    base = GaussianFit2D._chain_model_values(
        xfit, yfit, p, r.n, axisctx, ccfg;
        amp_min=r.amp_min, amp_range=r.amp_range)
    nd = length(base)
    d = length(p)
    Jq = Matrix{Float64}(undef, nd, d)
    prior_scale = max.((upper .- lower) ./ 2.0, LAPLACE_FD_MIN_STEP)
    for j in 1:d
        qj = p[j] / prior_scale[j]
        hq = max(LAPLACE_FD_REL_STEP * max(abs(qj), 1.0), LAPLACE_FD_MIN_STEP / prior_scale[j])
        dp = hq * prior_scale[j]
        lo = lower[j] + 1e-9
        hi = upper[j] - 1e-9
        p_minus_j = clamp(p[j] - dp, lo, hi)
        p_plus_j = clamp(p[j] + dp, lo, hi)
        if p_plus_j == p_minus_j
            Jq[:, j] .= 0.0
            continue
        end
        p_minus = copy(p); p_minus[j] = p_minus_j
        p_plus = copy(p); p_plus[j] = p_plus_j
        pred_minus = GaussianFit2D._chain_model_values(
            xfit, yfit, p_minus, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)
        pred_plus = GaussianFit2D._chain_model_values(
            xfit, yfit, p_plus, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)
        dq = (p_plus_j - p_minus_j) / prior_scale[j]
        Jq[:, j] .= (pred_plus .- pred_minus) ./ dq
    end
    return base, Jq
end

function _laplace_evidence_score(r, xfit, yfit, zfit, noise, axisctx, ccfg)
    r.success && r.valid || return nothing
    lower, upper = _chain_param_bounds_like_fit(r.n, ccfg)
    length(r.params) == length(lower) || return nothing
    pred, Jq = _laplace_fd_jacobian(r, xfit, yfit, axisctx, ccfg, lower, upper)
    resid = zfit .- pred
    nd = length(zfit)
    d = length(r.params)
    nd > d || return nothing
    n_eff = max(10, nd ÷ 9)
    eff_scale = n_eff / max(nd, 1)
    nu = ccfg.student_nu
    fit = 2.0 * eff_scale * GaussianFit2D._student_nll(resid, noise, nu)

    u = resid ./ max(noise, eps(Float64))
    weights = sqrt.(eff_scale .* (nu + 1.0) ./ (nu .+ u.^2)) ./ max(noise, eps(Float64))
    Jw = Jq .* reshape(weights, :, 1)
    lambdas = min.(svdvals(Jw).^2, LAPLACE_EIGEN_CAP)
    occam = sum(log1p, lambdas)
    d_eff = sum(lambdas ./ (1.0 .+ lambdas))
    sloppy_penalty = max(0.0, d - d_eff) * log(max(n_eff, 2))
    score = fit + occam + sloppy_penalty
    return (score=score, fit=fit, occam=occam, d_eff=d_eff,
            sloppy=sloppy_penalty, rss=sum(abs2, resid), d=d)
end

function _laplace_evidence_selection(img, pcfg, ccfg_ell, ccfg_circ, results_ell, results_circ)
    _, _, _, _, x_roi, y_roi, z_roi, noise = _load_chain_data(img, pcfg, ccfg_ell)
    axisctx_full = GaussianFit2D._weighted_roi_axis(x_roi, y_roi, z_roi)
    xfit_ell, yfit_ell, zfit_ell, axisctx_ell, _, _ = GaussianFit2D._chain_fit_data(
        x_roi, y_roi, z_roi, axisctx_full, ccfg_ell)
    xfit_circ, yfit_circ, zfit_circ, axisctx_circ, _, _ = GaussianFit2D._chain_fit_data(
        x_roi, y_roi, z_roi, axisctx_full, ccfg_circ)

    candidates = NamedTuple[]
    for r in results_ell
        resc = try
            _laplace_evidence_score(r, xfit_ell, yfit_ell, zfit_ell, noise, axisctx_ell, ccfg_ell)
        catch err
            nothing
        end
        resc === nothing || !isfinite(resc.score) || push!(candidates, merge((n=r.n, source="ell"), resc))
    end
    for r in results_circ
        resc = try
            _laplace_evidence_score(r, xfit_circ, yfit_circ, zfit_circ, noise, axisctx_circ, ccfg_circ)
        catch err
            nothing
        end
        resc === nothing || !isfinite(resc.score) || push!(candidates, merge((n=r.n, source="circ"), resc))
    end
    isempty(candidates) && return (nothing, NaN, 0, "NA", NaN, NaN, NaN)
    ordered = sort(candidates; by=c -> (c.score, c.d, c.rss, c.n))
    best = first(ordered)
    return (best.n, best.score, length(candidates), best.source, best.d_eff, best.occam, best.sloppy)
end

# ══════════════════════════════════════════════════════════════════════════════
# Resolved-lobe evidence (local geometrical resolvability)
# ══════════════════════════════════════════════════════════════════════════════

function _resolved_lobe_count(r, axisctx, ccfg, noise;
                              sep_thresh=LOBE_EVIDENCE_SEP_SIGMA,
                              valley_snr_thresh=LOBE_EVIDENCE_VALLEY_SNR,
                              valley_frac_thresh=LOBE_EVIDENCE_VALLEY_FRAC)
    r === nothing && return (0, 0)
    n = r.n
    n <= 1 && return (n, 0)

    _b, feats, _ts, _us, spars, _sperps = GaussianFit2D._decode_chain(
        r.params, n, axisctx, ccfg;
        amp_min=r.amp_min, amp_range=r.amp_range)

    parent = collect(1:n)
    function _find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    function _union!(a, b)
        ra, rb = _find(a), _find(b)
        ra != rb && (parent[rb] = ra)
    end

    n_unresolved = 0
    tiny = eps(Float64)
    for idx in 1:(n - 1)
        f_i = feats[idx]
        f_j = feats[idx + 1]
        d_nm = sqrt((f_i.x_nm - f_j.x_nm)^2 + (f_i.y_nm - f_j.y_nm)^2)
        sigma_pair = sqrt((spars[idx]^2 + spars[idx + 1]^2) / 2.0)
        sep_sigma = d_nm / max(sigma_pair, tiny)

        xs_line = collect(range(f_i.x_nm, f_j.x_nm, length=64))
        ys_line = collect(range(f_i.y_nm, f_j.y_nm, length=64))
        vals = GaussianFit2D._chain_model_values(
            xs_line, ys_line, r.params, n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)
        peak_i = GaussianFit2D._chain_model_values(
            [f_i.x_nm], [f_i.y_nm], r.params, n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)[1]
        peak_j = GaussianFit2D._chain_model_values(
            [f_j.x_nm], [f_j.y_nm], r.params, n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)[1]

        valley_min = minimum(vals[2:end-1])
        min_peak = min(peak_i, peak_j)
        valley_depth = min_peak - valley_min
        valley_snr = valley_depth / max(noise, tiny)
        valley_frac = valley_depth / max(abs(min_peak), tiny)

        if sep_sigma < sep_thresh &&
           (valley_snr < valley_snr_thresh || valley_frac < valley_frac_thresh)
            _union!(idx, idx + 1)
            n_unresolved += 1
        end
    end

    roots = Set{Int}()
    for k in 1:n
        push!(roots, _find(k))
    end
    return (length(roots), n_unresolved)
end

function _local_lobe_evidence_selection(n_eff::Int, results_ell, results_circ, axisctx, ccfg, noise)
    by_n = _result_by_n_prefer_ell(results_ell, results_circ)
    isempty(by_n) && return nothing

    eff_metrics = nothing
    for n_test in n_eff:-1:1
        r = get(by_n, n_test, nothing)
        r === nothing && continue
        n_resolved, n_unresolved = _resolved_lobe_count(r, axisctx, ccfg, noise)
        n_test == n_eff && (eff_metrics = (n_resolved, n_unresolved))
        # Down-only guard: accept the highest candidate with at most one
        # unresolved adjacent pair. This can reject redundant lobes, but it
        # cannot provide evidence for increasing N beyond the GCV choice.
        if n_resolved >= n_test - 1
            return (n_test, n_resolved, n_unresolved, true)
        end
    end
    # If no lower candidate passes the local resolvability rule, keep GCV.
    # The local diagnostic is then inconclusive rather than a selection failure.
    if eff_metrics !== nothing
        n_resolved, n_unresolved = eff_metrics
        return (n_eff, n_resolved, n_unresolved, false)
    end
    return nothing
end

# ══════════════════════════════════════════════════════════════════════════════
# Fwd/Bwd direction-consensus selector
# ══════════════════════════════════════════════════════════════════════════════

function _linear_recal_score(z_scan::Vector{Float64}, model::Vector{Float64})
    n = length(z_scan)
    n < 4 && return NaN, NaN
    A = hcat(model, ones(n))
    try
        coeff = A \ z_scan
        pred = coeff[1] .* model .+ coeff[2]
        rss = sum(abs2, z_scan .- pred)
        return coeff[1], rss
    catch
        return NaN, NaN
    end
end

_bilinear_val(xs, ys, zmat, x0, y0) = begin
    ix = clamp(searchsortedlast(xs, x0), 1, length(xs) - 1)
    iy = clamp(searchsortedlast(ys, y0), 1, length(ys) - 1)
    x1, x2 = xs[ix], xs[ix + 1]
    y1, y2 = ys[iy], ys[iy + 1]
    dx = (x0 - x1) / max(x2 - x1, eps(Float64))
    dy = (y0 - y1) / max(y2 - y1, eps(Float64))
    return (1 - dx) * (1 - dy) * zmat[iy, ix] +
           dx * (1 - dy) * zmat[iy, ix + 1] +
           (1 - dx) * dy * zmat[iy + 1, ix] +
           dx * dy * zmat[iy + 1, ix + 1]
end

function _fwd_bwd_consensus_selection(img, pcfg, ccfg_ell, results_ell, results_circ)
    # Use the fused ROI pixel coordinates (same for all channels)
    _, _, _, _, x_roi, y_roi, z_roi, noise = _load_chain_data(img, pcfg, ccfg_ell)
    axisctx_full = GaussianFit2D._weighted_roi_axis(x_roi, y_roi, z_roi)
    xfit, yfit, zfit, axisctx, _, _ = GaussianFit2D._chain_fit_data(
        x_roi, y_roi, z_roi, axisctx_full, ccfg_ell)

    # Preprocess fwd and bwd channels (same stride → same grid as fused)
    ch_fwd = GaussianFit2D.get_channel(img, "Z"; direction="fwd")
    xs_g, ys_g, _, z_fwd_g, _, _, _ = GaussianFit2D.preprocess_channel(img, ch_fwd, pcfg)
    ch_bwd = GaussianFit2D.get_channel(img, "Z"; direction="bwd")
    _, _, _, z_bwd_g, _, _, _ = GaussianFit2D.preprocess_channel(img, ch_bwd, pcfg)

    # Bilinear lookup: get fwd/bwd values at the same (xfit, yfit) positions
    z_fwd = [_bilinear_val(xs_g, ys_g, z_fwd_g, xfit[i], yfit[i]) for i in eachindex(xfit)]
    z_bwd = [_bilinear_val(xs_g, ys_g, z_bwd_g, xfit[i], yfit[i]) for i in eachindex(xfit)]

    # Get model candidates
    by_n = _result_by_n_prefer_ell(results_ell, results_circ)
    isempty(by_n) && return nothing

    scores = Tuple{Int,Float64}[]
    for n in sort(collect(keys(by_n)))
        r = by_n[n]
        # Model prediction at the same fit pixels
        pred = GaussianFit2D._chain_model_values(
            xfit, yfit, r.params, r.n, axisctx, ccfg_ell;
            amp_min=r.amp_min, amp_range=r.amp_range)
        good = isfinite.(z_fwd) .& isfinite.(z_bwd) .& isfinite.(pred)
        zf = z_fwd[good]; zb = z_bwd[good]; mp = pred[good]
        length(mp) < 20 && continue

        # Linear recalibration per scan
        a_fwd, rss_fwd = _linear_recal_score(zf, mp)
        a_bwd, rss_bwd = _linear_recal_score(zb, mp)
        (!isfinite(a_fwd) || !isfinite(a_bwd)) && continue

        nd_joint = 2 * length(mp)
        p_chain = GaussianFit2D._chain_nparams(r.n, ccfg_ell)
        p_eff = p_chain + 4
        nd_joint > p_eff || continue
        joint_gcv = nd_joint / (nd_joint - p_eff)^2 * (rss_fwd + rss_bwd)
        isfinite(joint_gcv) && push!(scores, (n, joint_gcv))
    end

    isempty(scores) && return nothing
    ordered = sort(scores; by=last)
    best = first(ordered)
    return (best[1], best[2], length(scores))
end

# ══════════════════════════════════════════════════════════════════════════════
# Primary selection dispatcher
# ══════════════════════════════════════════════════════════════════════════════

function _select_primary(n_eff::Int, eff_source::AbstractString, refined, policy::AbstractString)
    if policy == "gcv_with_robust_aicc_guard" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.n_refined < n_eff ? "robust_aicc_guard" : eff_source
    elseif policy == "spatial_blocked_cv" && refined.robust_n != "NA"
        return refined.n_refined, policy, "spatial_blocked_cv"
    elseif policy in ("support_marginalized_gcv", "support_marginalized_gcv_guard") && refined.robust_n != "NA"
        return refined.n_refined, policy, "support_marginalized_gcv"
    elseif policy == "slope_heuristic_mdl" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    elseif policy == "stability_selection" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    elseif policy == "local_lobe_evidence" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    elseif policy in ("laplace_evidence", "laplace_evidence_guard") && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    elseif policy == "fwd_bwd_consensus" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    elseif policy == "adaptive_support_rescue" && refined.robust_n != "NA"
        return refined.n_refined, policy, refined.source
    end
    return n_eff, "gcv", eff_source
end
