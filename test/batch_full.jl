#!/usr/bin/env julia
# Full 1D + 2D batch pipeline: fits, best-plots, and enriched summary.
# For each file: 1D slide fit + 2D elliptical chain + 2D circular chain + 6-panel plot.
# Standalone: discovers .sxm files directly if no triage TSV is present.
# Usage:
#   julia --project=. test/scripts/batch_full.jl [N_files]
#   julia --project=. test/scripts/batch_full.jl [N_files] --chunk i/n

using STMMolecularFit, GaussianFit2D, GaussianFit1D
using DelimitedFiles, Plots, Printf, Statistics

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
list_sxm_files(dir) = sort([f for f in readdir(dir) if endswith(lowercase(f), ".sxm")])
const TSV = "results/batch_triage_20240817_relaxed.tsv"
const OUTDIR = "results/best_plots"
const EXCLUDE = Set(["240817_001.sxm", "240817_008.sxm", "240817_009.sxm",
                      "240817_010.sxm", "240817_011.sxm", "240817_012.sxm",
                      "240817_013.sxm", "240817_014.sxm", "240817_016.sxm",
                      "240817_020.sxm", "240817_022.sxm", "240817_023.sxm",
                      "240817_025.sxm", "240817_028.sxm"])

function _parse_cli(args)
    n_files = 48
    chunk_idx = 1
    chunk_total = 1
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--chunk"
            i < length(args) || error("--chunk requires i/n, e.g. --chunk 2/4")
            chunk_arg = args[i + 1]
            i += 2
        elseif startswith(arg, "--chunk=")
            chunk_arg = split(arg, "=", limit=2)[2]
            i += 1
        elseif startswith(arg, "--")
            error("Unknown option: $arg")
        else
            n_files = parse(Int, arg)
            i += 1
            continue
        end

        parts = split(chunk_arg, "/")
        length(parts) == 2 || error("Invalid --chunk '$chunk_arg'; expected i/n, e.g. 2/4")
        chunk_idx = parse(Int, parts[1])
        chunk_total = parse(Int, parts[2])
    end
    chunk_total >= 1 || error("chunk total must be >= 1")
    1 <= chunk_idx <= chunk_total || error("chunk index must satisfy 1 <= i <= n")
    n_files >= 0 || error("N_files must be >= 0")
    return n_files, chunk_idx, chunk_total
end

const N_FILES, CHUNK_IDX, CHUNK_TOTAL = _parse_cli(ARGS)
mkpath(OUTDIR)

const FWHM_SIGMA = 2.355
const FWHM_MIN_1D_NM = 0.45
const FWHM_MAX_1D_NM = 1.20
const SIGMA_MIN_HARMONIZED_NM = FWHM_MIN_1D_NM / FWHM_SIGMA
const SIGMA_MAX_HARMONIZED_NM = FWHM_MAX_1D_NM / FWHM_SIGMA
const COLORMAP_RESID = cgrad([:blue, :lightgray, :red])
const SUMMARY_HEADER = [
    "filepath", "status", "classification",
    "N_ell", "N_circ", "N_1D", "N_eff", "eff_source",
    "dN_circ_ell", "dN_1D_ell", "dN_1D_circ",
    "BIC_ell", "BIC_circ", "sBIC_1D",
    "N_raw_ell", "BIC_raw_ell", "raw_valid_ell", "N_raw_circ", "BIC_raw_circ", "raw_valid_circ",
    "valid_ell", "valid_circ", "reason_ell", "reason_circ",
    "chi2_ell", "chi2_circ", "chi2_1D",
    "support_1D_nm", "support_2D_ell_nm", "support_2D_circ_nm", "support_mismatch_ell", "support_mismatch_circ",
    "Nset_ell_10", "Nset_circ_10", "Nset_1D_10", "common_N_10",
    "Nset_ell_hybrid", "Nset_circ_hybrid", "Nset_1D_hybrid", "common_N_hybrid",
    "threshold_ell_hybrid", "threshold_circ_hybrid", "threshold_1D_hybrid",
    "kappa_ell", "kappa_circ", "kappa_1D",
    "best_plot", "file_dir"
]

_bilinear_interp(xs, ys, z, x0, y0) = begin
    ix = clamp(searchsortedlast(xs, x0), 1, length(xs)-1)
    iy = clamp(searchsortedlast(ys, y0), 1, length(ys)-1)
    x1, x2 = xs[ix], xs[ix+1]; y1, y2 = ys[iy], ys[iy+1]
    tx = (x0-x1)/max(x2-x1, 1e-12); ty = (y0-y1)/max(y2-y1, 1e-12)
    return (1-tx)*(1-ty)*z[iy,ix] + tx*(1-ty)*z[iy,ix+1] + (1-tx)*ty*z[iy+1,ix] + tx*ty*z[iy+1,ix+1]
end

function _ellipse!(p, x0, y0, a, b, angle; color=:cyan, alpha=0.3, label="")
    θ = range(0, 2π, length=72)
    cosθ, sinθ = cos.(θ), sin.(θ)
    ca, sa = cos(angle), sin(angle)
    xe = x0 .+ a .* cosθ .* ca .- b .* sinθ .* sa
    ye = y0 .+ a .* cosθ .* sa .+ b .* sinθ .* ca
    plot!(p, xe, ye; color=color, alpha=alpha, label=label, linewidth=1.5)
end

_setstr(v) = isempty(v) ? "{}" : "{" * join(v, ",") * "}"

function _nset(results, best_score::Float64, threshold::Float64, scorefun)
    ns = Int[]
    for r in results
        getproperty(r, :success) || continue
        getproperty(r, :valid) || continue
        s = scorefun(r)
        isfinite(s) || continue
        s - best_score <= threshold + 1e-9 && push!(ns, getproperty(r, :n))
    end
    return sort(unique(ns))
end

function _best_valid_or_best(results, best_raw)
    valid = [r for r in results if getproperty(r, :success) && getproperty(r, :valid) && isfinite(getproperty(r, :bic))]
    isempty(valid) && return best_raw
    return sort(valid; by=r -> r.bic)[1]
end

function _select_effective_best(results_ell, results_circ; criterion="gcv")
    """Select best N using effective_criterion(N) = min(ell_criterion(N), circ_criterion(N)).
    Circular model is nested within elliptical → circ metric is a legitimate lower bound.
    Default criterion is GCV (analytical, no refit needed). Fallback to BIC for tiebreaking."""
    scorefun = criterion == "gcv" ? (r -> r.gcv) : (r -> r.bic)
    valid_ell = [r for r in results_ell if r.success && r.valid && isfinite(scorefun(r))]
    valid_circ = [r for r in results_circ if r.success && r.valid && isfinite(scorefun(r))]
    if isempty(valid_ell) && isempty(valid_circ)
        return nothing, 0
    end
    ell_by_n = Dict(r.n => r for r in valid_ell)
    circ_by_n = Dict(r.n => r for r in valid_circ)
    all_ns = sort(unique(vcat(collect(keys(ell_by_n)), collect(keys(circ_by_n)))))
    
    # Pass 1: find best N by effective criterion
    best_n = 0
    best_score = Inf
    for n in all_ns
        s_ell = haskey(ell_by_n, n) ? scorefun(ell_by_n[n]) : Inf
        s_circ = haskey(circ_by_n, n) ? scorefun(circ_by_n[n]) : Inf
        eff = min(s_ell, s_circ)
        if eff < best_score
            best_score = eff
            best_n = n
        end
    end
    final_n = best_n
    
    # CV tiebreaker: override if CV strongly disagrees
    best_n_cv = 0
    best_cv = Inf
    for n in all_ns
        for r in [get(ell_by_n, n, nothing), get(circ_by_n, n, nothing)]
            r === nothing && continue
            if isfinite(r.cv_nll_mean) && r.cv_nll_mean < best_cv
                best_cv = r.cv_nll_mean
                best_n_cv = n
            end
        end
    end
    
    if best_n_cv > 0 && best_n_cv != best_n && isfinite(best_cv)
        cv_at_best = Inf
        for r in [get(ell_by_n, best_n, nothing), get(circ_by_n, best_n, nothing)]
            r === nothing && continue
            if isfinite(r.cv_nll_mean); cv_at_best = min(cv_at_best, r.cv_nll_mean); end
        end
        cv_ratio = isfinite(cv_at_best) && best_cv > 0 ? cv_at_best / best_cv : 1.0
        if cv_ratio > 2.0 && best_n_cv < best_n
            final_n = best_n_cv
        end
    end
    
    r_ell = get(ell_by_n, final_n, nothing)
    r_circ = get(circ_by_n, final_n, nothing)
    if r_ell !== nothing && r_circ !== nothing
        best_result = r_ell.bic <= r_circ.bic + 10 ? r_ell : r_circ
    elseif r_ell !== nothing
        best_result = r_ell
    elseif r_circ !== nothing
        best_result = r_circ
    else
        return nothing, 0
    end
    return best_result, final_n
end

function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    """For each N fitted by the circular sweep, run a circ→ell LsqFit refinement.
    Returns a vector of elliptical ChainModelResult (one per N).
    NLopt is intentionally skipped — it always diverges from the isotropic start."""
    refined = GaussianFit2D.ChainModelResult[]
    isempty(results_circ) && return refined
    
    # Get fit data once
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    ac_full = ctx_circ.axisctx_full
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ac_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    
    # Config for LsqFit-only refinement
    ccfg_refine = deepcopy(ccfg_ell)
    ccfg_refine.skip_global = true
    ccfg_refine.max_iter = 50
    ccfg_refine.multistart = 1
    
    for r_c in results_circ
        r_c.success || continue
        n = r_c.n
        try
            # Expand circ params to elliptical format
            n_prefix = 3  # b0, bx, by
            split_idx = n_prefix + n + 1 + (n - 1) + n
            p_init = vcat(r_c.params[1:split_idx],
                          r_c.params[(split_idx+1):end],
                          r_c.params[(split_idx+1):end])
            # Run LsqFit-only from circ warm-start
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

function _nset_1d(results, best_score::Float64, threshold::Float64)
    ns = Int[]
    for r in results
        r.success || continue
        s = r.student_bic
        isfinite(s) || continue
        s - best_score <= threshold + 1e-9 && push!(ns, r.n_peaks)
    end
    return sort(unique(ns))
end

_intersect3(a, b, c) = sort(collect(intersect(Set(a), Set(b), Set(c))))

function _support_length(ctx)
    hasproperty(ctx, :support_meta) || return ctx.axisctx.tmax - ctx.axisctx.tmin
    return get(ctx.support_meta, :final_support_length_nm, ctx.axisctx.tmax - ctx.axisctx.tmin)
end

function _support_mismatch(l1d::Real, l2d::Real)
    return isfinite(l1d) && isfinite(l2d) && abs(l2d) > 1e-12 ? abs(l1d - l2d) / abs(l2d) : NaN
end

function _classification(best_ell, best_circ, best1d, common10, commonhybrid)
    !best_ell.valid && return "PROBLEMATIC"
    !best_circ.valid && return "PROBLEMATIC"
    ns = [best_ell.n, best_circ.n, best1d.n_peaks]
    length(unique(ns)) == 1 && return "ROBUST"
    !isempty(common10) && return "TOLERANT_10"
    !isempty(commonhybrid) && return "TOLERANT_HYBRID"
    maximum(ns) - minimum(ns) <= 1 && return "AMBIGU_MINOR"
    return "PROBLEMATIC"
end

function _write_scores(path::String, results, score_label::String, scorefun, selected)
    open(path, "w") do io
        println(io, join(["N", score_label, "delta_raw", "delta_valid", "valid", "raw_best", "selected", "reason", "chi2", "cv_nll_mean", "cv_nll_std"], '\t'))
        vals = [scorefun(r) for r in results if getproperty(r, :success) && isfinite(scorefun(r))]
        best_raw = isempty(vals) ? Inf : minimum(vals)
        valid_vals = [scorefun(r) for r in results if getproperty(r, :success) && getproperty(r, :valid) && isfinite(scorefun(r))]
        best_valid = isempty(valid_vals) ? best_raw : minimum(valid_vals)
        for r in sort(results; by=r -> getproperty(r, :n))
            getproperty(r, :success) || continue
            s = scorefun(r)
            isfinite(s) || continue
            println(io, join([r.n, s, s - best_raw, s - best_valid, getproperty(r, :valid), abs(s - best_raw) <= 1e-9, r === selected, getproperty(r, :reason), getproperty(r, :chi2_reduced), getproperty(r, :cv_nll_mean), getproperty(r, :cv_nll_std)], '\t'))
        end
    end
end

function _write_scores_1d(path::String, results)
    open(path, "w") do io
        println(io, join(["N", "sBIC", "delta", "competitive", "chi2", "kappa"], '\t'))
        vals = [r.student_bic for r in results if r.success]
        best = isempty(vals) ? Inf : minimum(vals)
        for r in sort(results; by=r -> r.n_peaks)
            r.success || continue
            println(io, join([r.n_peaks, r.student_bic, r.student_bic - best, r.competitive, r.chi2_red, r.kappa_max_adj], '\t'))
        end
    end
end

function make_best_plot(best_ell, ctx_ell, ccfg_ell, best_circ, ctx_circ, ccfg_circ, best_1d, x_1d, y_1d, cfg_1d, slide_mode, outpath)
    """6-panel (2×3) comparison: ell 2D | circ 2D | 1D  |  ell residuals | circ residuals | 1D residuals."""
    n_ell, n_circ, n1d = best_ell.n, best_circ.n, best_1d.n_peaks

    # ── Common helper: add overlays to a 2D panel ──
    function _add_2d_overlays!(p, best, ctx, ccfg, n)
        ax_ctx = ctx.axisctx
        ox, oy = ax_ctx.origin; ax, ay = ax_ctx.axis
        t_all = (ctx.x .- ox) .* ax .+ (ctx.y .- oy) .* ay
        # ROI contour (white)
        contour!(p, ctx.xs, ctx.ys, Float64.(ctx.mask); levels=[0.5], color=:white, linewidth=1.2)
        # Axis line (yellow)
        plot!(p, [ox+minimum(t_all)*ax, ox+maximum(t_all)*ax],
              [oy+minimum(t_all)*ay, oy+maximum(t_all)*ay];
              color=:yellow, linewidth=1.5, label="")
        # Ridge path (green): search ±0.35 nm perpendicular, then smooth
        perp = [-ay, ax]
        t_ridge = range(minimum(t_all), maximum(t_all), length=120)
        hw_ridge = 0.35
        rx, ry = Float64[], Float64[]
        for t_val in t_ridge
            best_v, bx, by = -Inf, 0.0, 0.0
            for u_r in range(-hw_ridge, hw_ridge, step=0.02)
                x0 = ox + t_val*ax + u_r*perp[1]
                y0 = oy + t_val*ay + u_r*perp[2]
                if minimum(ctx.xs) <= x0 <= maximum(ctx.xs) && minimum(ctx.ys) <= y0 <= maximum(ctx.ys)
                    v = _bilinear_interp(ctx.xs, ctx.ys, ctx.zimg, x0, y0)
                    if v > best_v; best_v = v; bx = x0; by = y0; end
                end
            end
            push!(rx, bx); push!(ry, by)
        end
        plot!(p, rx, ry; color=:white, linewidth=2, alpha=0.9, label="")
        # FWHM ellipses/circles (cyan)
        if n > 0 && best.success
            _, feats, _, _, _, _ = GaussianFit2D._decode_chain(best.params, n, ax_ctx, ccfg;
                amp_min=best.amp_min, amp_range=best.amp_range)
            axis_angle = atan(ax, ay)
            xs_c = Float64[]; ys_c = Float64[]
            for f in feats
                a_ell = f.sigma_x_nm * FWHM_SIGMA / 2
                b_ell = f.sigma_y_nm * FWHM_SIGMA / 2
                _ellipse!(p, f.x_nm, f.y_nm, a_ell, b_ell, axis_angle; color=:cyan, alpha=0.4, label="")
                push!(xs_c, f.x_nm); push!(ys_c, f.y_nm)
            end
            scatter!(p, xs_c, ys_c; marker=:cross, markersize=10, color=:red, linewidth=2.5, label="")
        end
    end

    # ── Shared ROI bounds (use ctx_ell as reference) ──
    xs_ell, ys_ell, zimg_ell, mask_ell = ctx_ell.xs, ctx_ell.ys, ctx_ell.zimg, ctx_ell.mask
    roi_rows = [iy for iy in eachindex(ys_ell) if any(mask_ell[iy, :])]
    roi_cols = [ix for ix in eachindex(xs_ell) if any(mask_ell[:, ix])]
    rx = isempty(roi_cols) ? (minimum(xs_ell)-0.3, maximum(xs_ell)+0.3) : (xs_ell[minimum(roi_cols)]-0.3, xs_ell[maximum(roi_cols)]+0.3)
    ry = isempty(roi_rows) ? (minimum(ys_ell)-0.3, maximum(ys_ell)+0.3) : (ys_ell[minimum(roi_rows)]-0.3, ys_ell[maximum(roi_rows)]+0.3)
    xmin, xmax = rx[1], rx[2]; ymin, ymax = ry[1], ry[2]
    if xmin >= xmax; xmin = minimum(xs_ell); xmax = maximum(xs_ell); end
    if ymin >= ymax; ymin = minimum(ys_ell); ymax = maximum(ys_ell); end

    z_clims = (quantile(vec(zimg_ell), 0.10), quantile(vec(zimg_ell), 0.995))

    # ═══════════════════════════════════════════════════════
    # Row 1 (top): heatmaps + fit
    # ═══════════════════════════════════════════════════════

    # ── Panel (1,1): 2D elliptical heatmap ──
    p_ell = heatmap(xs_ell, ys_ell, zimg_ell; aspect_ratio=:equal, colormap=:thermal, clims=z_clims,
                    title="2D ell N=$n_ell", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_ell, (xmin, xmax)); ylims!(p_ell, (ymin, ymax))
    _add_2d_overlays!(p_ell, best_ell, ctx_ell, ccfg_ell, n_ell)

    # ── Panel (1,2): 2D circular heatmap ──
    xs_circ, ys_circ, zimg_circ, mask_circ = ctx_circ.xs, ctx_circ.ys, ctx_circ.zimg, ctx_circ.mask
    p_circ = heatmap(xs_circ, ys_circ, zimg_circ; aspect_ratio=:equal, colormap=:thermal, clims=z_clims,
                     title="2D circ N=$n_circ", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_circ, (xmin, xmax)); ylims!(p_circ, (ymin, ymax))
    _add_2d_overlays!(p_circ, best_circ, ctx_circ, ccfg_circ, n_circ)

    # ── Panel (1,3): 1D data + fit + components ──
    y1d_pred = GaussianFit1D.predict_fit(x_1d, best_1d, cfg_1d)
    y1d_resid = y_1d .- y1d_pred

    p_1d = plot(x_1d, y_1d; color=:gray, alpha=0.7, label="data", linewidth=1)
    plot!(p_1d, x_1d, y1d_pred; color=:red, linewidth=2, label="fit N=$n1d")

    centers = GaussianFit1D._params_to_centers(best_1d.popt, n1d)
    comp_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :lime, :teal, :gold]
    asymmetric = cfg_1d.asymmetric_edges && n1d >= 2
    y0 = best_1d.popt[1]
    for (i, c) in enumerate(centers)
        idx = i - 1
        A = GaussianFit1D._get_amplitude(best_1d.popt, idx)
        σ_in = GaussianFit1D._get_sigma(best_1d.popt, idx)
        if asymmetric && (idx == 0 || idx == n1d - 1)
            σ_out = idx == 0 ? best_1d.popt[end-1] : best_1d.popt[end]
            z = x_1d .- c
            s = idx == 0 ? (z .< 0) .* σ_out .+ (z .>= 0) .* σ_in :
                           (z .< 0) .* σ_in .+ (z .>= 0) .* σ_out
            y_comp = y0 .+ A .* exp.(-0.5 .* (z ./ s).^2)
        else
            y_comp = y0 .+ A .* exp.(-0.5 .* ((x_1d .- c) ./ max(σ_in, 1e-9)).^2)
        end
        col = comp_colors[mod1(i, length(comp_colors))]
        plot!(p_1d, x_1d, y_comp; color=col, alpha=0.55, linestyle=:dash, linewidth=1.5, label="")
    end
    xlabel!(p_1d, "position (nm)"); ylabel!(p_1d, "intensity")
    title!(p_1d, "1D N=$n1d ΔsBIC=0 (sBIC=$(round(best_1d.student_bic, digits=0)))")

    # ═══════════════════════════════════════════════════════
    # Row 2 (bottom): residuals
    # ═══════════════════════════════════════════════════════

    # ── Panel (2,1): 2D elliptical residuals ──
    pred_img_ell = zeros(size(zimg_ell))
    ax_ell = ctx_ell.axisctx
    for iy in eachindex(ys_ell), ix in eachindex(xs_ell)
        pred_img_ell[iy, ix] = GaussianFit2D._chain_model_values([xs_ell[ix]], [ys_ell[iy]], best_ell.params, n_ell, ax_ell, ccfg_ell; amp_min=best_ell.amp_min, amp_range=best_ell.amp_range)[1]
    end
    noise_ell = ctx_ell.noise
    resid_ell = (zimg_ell .- pred_img_ell) .* Float64.(mask_ell) ./ max(noise_ell, 1e-12)
    p_res_ell = heatmap(xs_ell, ys_ell, resid_ell; aspect_ratio=:equal,
                        colormap=COLORMAP_RESID, clims=(-3, 3),
                        title="2D ell residuals", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_res_ell, (xmin, xmax)); ylims!(p_res_ell, (ymin, ymax))

    # ── Panel (2,2): 2D circular residuals ──
    pred_img_circ = zeros(size(zimg_circ))
    ax_circ = ctx_circ.axisctx
    for iy in eachindex(ys_circ), ix in eachindex(xs_circ)
        pred_img_circ[iy, ix] = GaussianFit2D._chain_model_values([xs_circ[ix]], [ys_circ[iy]], best_circ.params, n_circ, ax_circ, ccfg_circ; amp_min=best_circ.amp_min, amp_range=best_circ.amp_range)[1]
    end
    noise_circ = ctx_circ.noise
    resid_circ = (zimg_circ .- pred_img_circ) .* Float64.(mask_circ) ./ max(noise_circ, 1e-12)
    p_res_circ = heatmap(xs_circ, ys_circ, resid_circ; aspect_ratio=:equal,
                         colormap=COLORMAP_RESID, clims=(-3, 3),
                         title="2D circ residuals", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_res_circ, (xmin, xmax)); ylims!(p_res_circ, (ymin, ymax))

    # ── Panel (2,3): 1D residuals ──
    p_res_1d = plot(x_1d, y1d_resid; color=:red, label="", linewidth=1)
    hline!(p_res_1d, [0]; color=:gray, linestyle=:dash, label="")
    xlabel!(p_res_1d, "position (nm)"); ylabel!(p_res_1d, "residual")
    title!(p_res_1d, "1D residuals  σ=$(round(std(y1d_resid), digits=5))")

    # ── Global title ──
    title_str = "elliptical β=$(round(best_ell.bic, digits=0)) vs circular β=$(round(best_circ.bic, digits=0)) vs 1D β=$(round(best_1d.student_bic, digits=0)) | slide: $slide_mode"

    l = @layout grid(2, 3, heights=[0.65, 0.35])
    fig = plot(p_ell, p_circ, p_1d, p_res_ell, p_res_circ, p_res_1d;
               layout=l, size=(2400, 800),
               plot_title=title_str, plot_titlefontsize=8,
               left_margin=0Plots.mm, right_margin=0Plots.mm,
               top_margin=1Plots.mm, bottom_margin=0Plots.mm,
               margin=0.5Plots.mm)
    savefig(fig, outpath)
end

# ── Load candidate files: prefer triage TSV when present, otherwise standalone discovery ──
cands = Tuple{String,Int,Float64}[]
if isfile(TSV)
    data = readdlm(TSV)
    for i in 2:size(data, 1)
        data[i, 5] isa Bool && data[i, 5] || continue
        push!(cands, (string(data[i, 1]), Int(data[i, 3]), Float64(data[i, 4])))
    end
    sort!(cands, by = x -> x[3])
else
    @warn "Triage TSV not found: $TSV; discovering SXM files directly from $DATA_DIR"
    cands = [(fn, 0, Inf) for fn in list_sxm_files(DATA_DIR)]
end
# Exclude files that are not chitosan chains
cands = [(fn, n, bic) for (fn, n, bic) in cands if !(fn in EXCLUDE)]
to_process_base = cands[1:min(N_FILES, length(cands))]
to_process_all = CHUNK_TOTAL == 1 ? to_process_base : [f for (i, f) in enumerate(to_process_base) if mod1(i, CHUNK_TOTAL) == CHUNK_IDX]
if CHUNK_TOTAL > 1
    @printf("Chunk %d/%d: %d of %d selected files\n", CHUNK_IDX, CHUNK_TOTAL, length(to_process_all), length(to_process_base))
end

# ── Skip already-processed files (plot exists + enriched summary has ok row) ──
already_done = Set{String}()
summary_name = CHUNK_TOTAL == 1 ? "summary_overlap060_hard.tsv" : @sprintf("summary_overlap060_hard_chunk%02dof%02d.tsv", CHUNK_IDX, CHUNK_TOTAL)
summary_file = joinpath(OUTDIR, summary_name)
if isfile(summary_file)
    for (i, line) in enumerate(eachline(summary_file))
        i == 1 && continue
        isempty(strip(line)) && continue
        parts = split(line, '\t')
        length(parts) >= 2 || continue
        fn = parts[1]
        status = parts[2]
        png = joinpath(OUTDIR, replace(fn, r"\.sxm$"i => "_best.png"))
        if status == "ok" && isfile(png)
            push!(already_done, fn)
        end
    end
end
to_process = [f for f in to_process_all if !(f[1] in already_done)]
if length(already_done) > 0
    println("  $(length(already_done)) files already done, $(length(to_process)) remaining")
end

# ── 2D config (relaxed) ──
pcfg = GaussianFit2D.PatternConfig(filepath="", channel="Z", direction="fwd",
    stride=1, flatten="plane+rows", smooth_radius_px=1, output_dir=OUTDIR, no_plot=false)
ccfg = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
    spacing_min_nm=0.35, spacing_max_nm=0.75, fit_width_nm=0.15,
    support_threshold_fraction=0.25, support_noise_k=2.5, support_padding_nm=0.20,
    max_overlap=0.6,
    global_maxtime=10.0, global_maxiter=10000, cv_folds=5,
    multistart=1,      # bootstrap only (fast; default is 1)
    sigma_parallel_min_nm=SIGMA_MIN_HARMONIZED_NM,
    sigma_parallel_max_nm=SIGMA_MAX_HARMONIZED_NM,
    sigma_perp_min_nm=SIGMA_MIN_HARMONIZED_NM,
    sigma_perp_max_nm=SIGMA_MAX_HARMONIZED_NM,
    intelligent_sweep=true, fuse_z_bwd=true,
    chain_tilted_baseline=true)
# Circular 2D config (same settings, circular sigmas)
ccfg_circ = deepcopy(ccfg)
ccfg_circ.chain_circular_sigmas = true

# ── 1D config ──
scfg = STMMolecularFit.SlideConfig(width_nm=0.70, support_threshold_fraction=0.20,
    support_noise_k=2.5, support_padding_nm=0.20, output_dir=OUTDIR, no_plot=true)
fcfg = STMMolecularFit.FitSlideConfig(min_spacing=0.35, max_spacing=0.75, max_overlap=0.6, output_dir=OUTDIR)

# ── Process ──
# Write header once if summary file is new
if !isfile(summary_file)
    open(summary_file, "w") do io
        println(io, join(SUMMARY_HEADER, '\t'))
    end
end

println("Generating best plots in $OUTDIR/ ...")
ntot = length(to_process)
t0 = time()
summary_lock = ReentrantLock()
n_parallel = min(4, Threads.nthreads(), ntot)
println("Processing $ntot files ($n_parallel parallel workers)...")
flush(stdout)

processed = Threads.Atomic{Int}(0)

Threads.@threads for idx in 1:ntot
    local fn, n2d_tsv, bic2d_tsv = to_process[idx]
    local fp = joinpath(DATA_DIR, fn)
    local file_base = splitext(fn)[1]
    local file_dir = joinpath(OUTDIR, file_base)
    mkpath(file_dir)
    local outpath = joinpath(OUTDIR, replace(fn, r"\.sxm$"i => "_best.png"))
    
    Threads.atomic_add!(processed, 1)
    local done = processed[]
    local pct = done / ntot * 100
    local elapsed = time() - t0
    local eta = done > 0 ? elapsed / done * (ntot - done) : 0
    local eta_min = eta > 0 ? div(round(Int, eta), 60) : 0
    local eta_sec = eta > 0 ? round(Int, eta) % 60 : 0
    @printf("[%2d/%2d %3.0f%%  ETA %d:%02d] %-24s\n", done, ntot, pct, eta_min, eta_sec, fn)
    flush(stdout)
    try
        # 1D fit
        img = STMMolecularFit.read_sxm(fp)
        scfg_file = deepcopy(scfg); scfg_file.output_dir = file_dir
        fcfg_file = deepcopy(fcfg); fcfg_file.output_dir = file_dir
        slide = STMMolecularFit.extract_slide(img, scfg_file)
        if isempty(slide.x) || isempty(slide.y)
            @warn("$fn: empty slide profile (no molecule found?), skipping")
            lock(summary_lock) do
                open(summary_file, "a") do io
                    println(io, join([fn, "no_molecule", "0", "", "", "", "", "", "", "", ""], '\t'))
                end
            end
        else
        fit_1d = STMMolecularFit.fit_slide(slide, fcfg_file)
        best1d = GaussianFit1D.best_result(fit_1d.fit_run)
        x_1d, y_1d = fit_1d.fit_run.x, fit_1d.fit_run.y
        cfg_1d = fit_1d.fit_run.cfg

        # 2D circular sweep (primary: always converges)
        img2d = GaussianFit2D.read_sxm(fp)
        pcfg_file = deepcopy(pcfg); pcfg_file.filepath = fp; pcfg_file.output_dir = file_dir
        results_circ, best_circ_raw, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img2d, pcfg_file, ccfg_circ)
        best_circ_sweep = _best_valid_or_best(results_circ, best_circ_raw)

        # circ→ell LsqFit refinement at each N (replaces NLopt elliptical sweep)
        # NLopt is intentionally excluded — it always diverges in 33D sigma space.
        results_ell = _refine_circ_to_ell(results_circ, img2d, pcfg_file, ccfg, ctx_circ)
        ctx_ell = ctx_circ  # circ→ell refinement shares the circular context (same data/axis)
        best_ell_raw = isempty(results_ell) ? best_circ_raw : results_ell[1]
        best_ell_sweep = isempty(results_ell) ? best_circ_sweep : _best_valid_or_best(results_ell, best_ell_raw)

        # Effective best: use min(ell_GCV, circ_GCV) per N for model selection (default GCV)
        best_eff, best_n_eff = _select_effective_best(results_ell, results_circ; criterion="gcv")
        if best_eff === nothing
            best_eff = best_ell_sweep; best_n_eff = best_ell_sweep.n
        end

        # 6-panel combined plot
        make_best_plot(best_ell_sweep, ctx_circ, ccfg, best_circ_sweep, ctx_circ, ccfg_circ,
                       best1d, x_1d, y_1d, cfg_1d, string(scfg.slide_mode), outpath)

        dn = best_eff.n - best1d.n_peaks
        println("Neff=$(best_n_eff) Nell=$(best_ell_sweep.n) Ncirc=$(best_circ_sweep.n) N1D=$(best1d.n_peaks) Δ1D-eff=$(best1d.n_peaks - best_n_eff) ✓")

        _write_scores(joinpath(file_dir, "ell_scores.tsv"), results_ell, "BIC", r -> r.bic, best_ell_sweep)
        _write_scores(joinpath(file_dir, "circ_scores.tsv"), results_circ, "BIC", r -> r.bic, best_circ_sweep)
        _write_scores_1d(joinpath(file_dir, "fit_1d_scores.tsv"), fit_1d.fit_run.all_results)

        thr_ell = max(10.0, 0.01 * abs(best_ell_sweep.bic))
        thr_circ = max(10.0, 0.01 * abs(best_circ_sweep.bic))
        thr_1d = max(10.0, 0.01 * abs(best1d.student_bic))
        nset_ell_10 = _nset(results_ell, best_ell_sweep.bic, 10.0, r -> r.bic)
        nset_circ_10 = _nset(results_circ, best_circ_sweep.bic, 10.0, r -> r.bic)
        nset_1d_10 = _nset_1d(fit_1d.fit_run.all_results, best1d.student_bic, 10.0)
        common10 = _intersect3(nset_ell_10, nset_circ_10, nset_1d_10)
        nset_ell_h = _nset(results_ell, best_ell_sweep.bic, thr_ell, r -> r.bic)
        nset_circ_h = _nset(results_circ, best_circ_sweep.bic, thr_circ, r -> r.bic)
        nset_1d_h = _nset_1d(fit_1d.fit_run.all_results, best1d.student_bic, thr_1d)
        commonh = _intersect3(nset_ell_h, nset_circ_h, nset_1d_h)
        classif = _classification(best_eff, best_circ_sweep, best1d, common10, commonh)
        support_1d = slide.support_length_nm
        support_ell = _support_length(ctx_ell)
        support_circ = _support_length(ctx_circ)
        mismatch_ell = _support_mismatch(support_1d, support_ell)
        mismatch_circ = _support_mismatch(support_1d, support_circ)

        # Write summary (thread-safe)
        eff_source = "refined"
        lock(summary_lock) do
            open(summary_file, "a") do io
            println(io, join([fn, "ok", classif,
                              best_ell_sweep.n, best_circ_sweep.n, best1d.n_peaks,
                              best_n_eff, eff_source,
                              best_circ_sweep.n - best_ell_sweep.n, best1d.n_peaks - best_ell_sweep.n, best1d.n_peaks - best_circ_sweep.n,
                              round(best_ell_sweep.bic, digits=3), round(best_circ_sweep.bic, digits=3), round(best1d.student_bic, digits=3),
                              best_ell_raw.n, round(best_ell_raw.bic, digits=3), best_ell_raw.valid,
                              best_circ_raw.n, round(best_circ_raw.bic, digits=3), best_circ_raw.valid,
                              best_ell_sweep.valid, best_circ_sweep.valid, best_ell_sweep.reason, best_circ_sweep.reason,
                              best_ell_sweep.chi2_reduced, best_circ_sweep.chi2_reduced, best1d.chi2_red,
                              support_1d, support_ell, support_circ, mismatch_ell, mismatch_circ,
                              _setstr(nset_ell_10), _setstr(nset_circ_10), _setstr(nset_1d_10), _setstr(common10),
                              _setstr(nset_ell_h), _setstr(nset_circ_h), _setstr(nset_1d_h), _setstr(commonh),
                              thr_ell, thr_circ, thr_1d,
                              best_ell_sweep.kappa_max_adj, best_circ_sweep.kappa_max_adj, best1d.kappa_max_adj,
                              outpath, file_dir], '\t'))
        end
        end  # lock
        end  # else
    catch e
        msg = sprint(showerror, e)
        msg_short = length(msg) <= 80 ? msg : msg[1:min(end, 80)] * "..."
        println("FAILED: $msg_short")
        lock(summary_lock) do
            open(summary_file, "a") do io
            row = fill("ERR", length(SUMMARY_HEADER))
            row[1] = fn
            row[2] = "error"
            row[3] = "PROBLEMATIC"
            println(io, join(row, '\t'))
        end
        end  # lock
    end
end

# ── Count matches ──
n_done = length(already_done) + length(to_process)
n_success = n_done  # rough estimate; for precise count, re-read summary
println("\nDone! Processed $(length(to_process)) new files in $OUTDIR/")
println("Summary (appended): $summary_file")
