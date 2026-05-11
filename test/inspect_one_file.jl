#!/usr/bin/env julia
# Deep single-file comparison: 1D slide vs 2D elliptical vs 2D circular chain fits.
# Usage: julia --project=. test/scripts/inspect_one_file.jl <filepath.sxm> [output_dir]

using STMMolecularFit, GaussianFit2D, GaussianFit1D
using Plots, Printf, Statistics

const FWHM_SIGMA = 2.355
const FWHM_MIN_1D_NM = 0.45
const FWHM_MAX_1D_NM = 1.20
const SIGMA_MIN = FWHM_MIN_1D_NM / FWHM_SIGMA
const SIGMA_MAX = FWHM_MAX_1D_NM / FWHM_SIGMA
const SPACING_MIN = 0.35
const SPACING_MAX = 0.75
const OVERLAP = 0.60
const COLORMAP_RESID = cgrad([:blue, :lightgray, :red])

function _ellipse!(p, x0, y0, a, b, angle; color=:cyan, alpha=0.3, label="")
    θ = range(0, 2π, length=72)
    cosθ, sinθ = cos.(θ), sin.(θ)
    ca, sa = cos(angle), sin(angle)
    xe = x0 .+ a .* cosθ .* ca .- b .* sinθ .* sa
    ye = y0 .+ a .* cosθ .* sa .+ b .* sinθ .* ca
    plot!(p, xe, ye; color=color, alpha=alpha, label=label, linewidth=1.5)
end

function main()
    length(ARGS) >= 1 || error("Usage: julia --project=. test/scripts/inspect_one_file.jl <filepath.sxm> [output_dir]")
    filepath = ARGS[1]
    output_dir = length(ARGS) >= 2 ? ARGS[2] : "results/circ_vs_ell"
    mkpath(output_dir)

    println("━━━ Loading image: $filepath ━━━")

    # ── Shared 2D PatternConfig (from centralize_best_plots.jl) ──
    pcfg = GaussianFit2D.PatternConfig(
        filepath=filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1,
        output_dir=output_dir, no_plot=false)

    # ── 2D config base (widths harmonized with 1D FWHM bounds 0.45..1.20 nm) ──
    function make_ccfg(; circular=false)
        return GaussianFit2D.ChainSweepConfig(
            n_min=2, n_max=14,
            spacing_min_nm=SPACING_MIN, spacing_max_nm=SPACING_MAX,
            max_overlap=OVERLAP, fit_width_nm=0.15,
            support_threshold_fraction=0.20, support_noise_k=2.5, support_padding_nm=0.20,
            global_maxtime=10.0, global_maxiter=10000, cv_folds=3,
            sigma_parallel_min_nm=SIGMA_MIN, sigma_parallel_max_nm=SIGMA_MAX,
            sigma_perp_min_nm=SIGMA_MIN, sigma_perp_max_nm=SIGMA_MAX,
            intelligent_sweep=true, fuse_z_bwd=true,
            chain_circular_sigmas=circular)
    end

    # ── 1D config (from centralize_best_plots.jl) ──
    scfg = STMMolecularFit.SlideConfig(
        width_nm=0.30, support_threshold_fraction=0.20,
        support_noise_k=2.5, support_padding_nm=0.20,
        output_dir=output_dir, no_plot=true)
    fcfg = STMMolecularFit.FitSlideConfig(min_spacing=0.35, max_spacing=0.75, max_overlap=0.6, output_dir=output_dir)

    # ═══════════════════════════════════════════════════════════
    # 1. 2D elliptical sweep (chain_circular_sigmas=false)
    # ═══════════════════════════════════════════════════════════
    println("\n━━━ 2D elliptical sweep (circular=false) ━━━")
    ccfg_ell = make_ccfg(circular=false)
    img_ell = GaussianFit2D.read_sxm(filepath)
    results_ell, best_ell, ctx_ell = GaussianFit2D.chain_gaussian_sweep(img_ell, pcfg, ccfg_ell)
    println("  → N_ell=$(best_ell.n)  BIC_ell=$(round(best_ell.bic, digits=1))  spar=$(round(best_ell.sigma_parallel_nm, digits=3))  sperp=$(round(best_ell.sigma_perp_nm, digits=3))")

    # ═══════════════════════════════════════════════════════════
    # 2. 2D circular sweep (chain_circular_sigmas=true)
    # ═══════════════════════════════════════════════════════════
    println("\n━━━ 2D circular sweep (circular=true) ━━━")
    ccfg_circ = make_ccfg(circular=true)
    img_circ = GaussianFit2D.read_sxm(filepath)
    results_circ, best_circ, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img_circ, pcfg, ccfg_circ)
    println("  → N_circ=$(best_circ.n)  BIC_circ=$(round(best_circ.bic, digits=1))  sigma=$(round(best_circ.sigma_parallel_nm, digits=3))")

    # ═══════════════════════════════════════════════════════════
    # 3. 1D slide extraction + fit
    # ═══════════════════════════════════════════════════════════
    println("\n━━━ 1D slide fit ━━━")
    img_1d = STMMolecularFit.read_sxm(filepath)
    slide = STMMolecularFit.extract_slide(img_1d, scfg)
    fit_1d = STMMolecularFit.fit_slide(slide, fcfg)
    best1d = GaussianFit1D.best_result(fit_1d.fit_run)
    x_1d, y_1d = fit_1d.fit_run.x, fit_1d.fit_run.y
    cfg_1d = fit_1d.fit_run.cfg
    n1d = best1d.n_peaks
    println("  → N_1d=$(n1d)  sBIC_1d=$(round(best1d.student_bic, digits=1))")

    # ═══════════════════════════════════════════════════════════
    # 4. Build 3×2 comparison plot
    # ═══════════════════════════════════════════════════════════

    # helper: compute ROI bounds
    function roi_bounds(xs, ys, mask)
        roi_rows = [iy for iy in eachindex(ys) if any(mask[iy, :])]
        roi_cols = [ix for ix in eachindex(xs) if any(mask[:, ix])]
        rxmin = isempty(roi_cols) ? minimum(xs) : xs[minimum(roi_cols)] - 0.5
        rxmax = isempty(roi_cols) ? maximum(xs) : xs[maximum(roi_cols)] + 0.5
        rymin = isempty(roi_rows) ? minimum(ys) : ys[minimum(roi_rows)] - 0.5
        rymax = isempty(roi_rows) ? maximum(ys) : ys[maximum(roi_rows)] + 0.5
        return rxmin, rxmax, rymin, rymax
    end

    # helper: build prediction image and decode features
    function build_prediction(xs, ys, best, ctx, ccfg)
        n = best.n
        pred_img = zeros(size(ctx.zimg))
        for iy in eachindex(ys), ix in eachindex(xs)
            pred_img[iy, ix] = best.n == 0 ? best.params[1] :
                GaussianFit2D._chain_model_values([xs[ix]], [ys[iy]], best.params, n, ctx.axisctx, ccfg; amp_min=best.amp_min, amp_range=best.amp_range)[1]
        end
        feats = nothing
        if n > 0 && best.success
            _, feats, _, _, _, _ = GaussianFit2D._decode_chain(best.params, n, ctx.axisctx, ccfg;
                amp_min=best.amp_min, amp_range=best.amp_range)
        end
        return pred_img, feats
    end

    # helper: draw 2D heatmap + mask + axis + ellipses
    function make_2d_heatmap(xs, ys, zimg, mask, best, ctx, ccfg, axctx; title_str="", colorbar=false)
        n = best.n
        pred_img, feats = build_prediction(xs, ys, best, ctx, ccfg)
        rxmin, rxmax, rymin, rymax = roi_bounds(xs, ys, mask)

        z_clims = (quantile(vec(zimg), 0.10), quantile(vec(zimg), 0.995))
        p = heatmap(xs, ys, zimg; aspect_ratio=:equal, title=title_str,
                    xlabel="x (nm)", ylabel="y (nm)", colorbar=colorbar,
                    colormap=:thermal, clims=z_clims)
        xlims!(p, rxmin, rxmax); ylims!(p, rymin, rymax)
        contour!(p, xs, ys, Float64.(mask); levels=[0.5], color=:white, linewidth=1.5)
        ox, oy = axctx.origin; ax, ay = axctx.axis
        t_all = (ctx.x .- ox) .* ax .+ (ctx.y .- oy) .* ay
        plot!(p, [ox + minimum(t_all)*ax, ox + maximum(t_all)*ax],
              [oy + minimum(t_all)*ay, oy + maximum(t_all)*ay];
              color=:yellow, linewidth=2, label="")

        if feats !== nothing
            axis_angle = atan(ax, ay)
            for f in feats
                if ccfg.chain_circular_sigmas
                    # Circular sigma: same a and b
                    a_ell = f.sigma_x_nm * FWHM_SIGMA / 2
                    b_ell = a_ell
                else
                    a_ell = f.sigma_x_nm * FWHM_SIGMA / 2
                    b_ell = f.sigma_y_nm * FWHM_SIGMA / 2
                end
                _ellipse!(p, f.x_nm, f.y_nm, a_ell, b_ell, axis_angle; color=:cyan, alpha=0.5, label="")
            end
        end

        return p, pred_img
    end

    # helper: 2D residuals heatmap
    function make_2d_residuals(xs, ys, zimg, pred_img, mask, noise, axctx; title_str="", colorbar=false)
        resid_img = (zimg .- pred_img) .* Float64.(mask) ./ max(noise, 1e-12)
        rxmin, rxmax, rymin, rymax = roi_bounds(xs, ys, mask)
        p = heatmap(xs, ys, resid_img; aspect_ratio=:equal, title=title_str,
                    xlabel="x (nm)", ylabel="y (nm)", colorbar=colorbar,
                    colormap=COLORMAP_RESID, clims=(-3, 3))
        xlims!(p, rxmin, rxmax); ylims!(p, rymin, rymax)
        return p
    end

    # ── Compute 1D t-shift relative to elliptical as reference ──
    n_ell = best_ell.n
    n_circ = best_circ.n

    t_shift_ell = 0.0
    if n_ell > 0 && best_ell.success
        _, _, ts, _, _, _ = GaussianFit2D._decode_chain(best_ell.params, n_ell, ctx_ell.axisctx, ccfg_ell;
            amp_min=best_ell.amp_min, amp_range=best_ell.amp_range)
        t_shift_ell = ts[1]
    end

    t_shift_circ = 0.0
    if n_circ > 0 && best_circ.success
        _, _, ts, _, _, _ = GaussianFit2D._decode_chain(best_circ.params, n_circ, ctx_circ.axisctx, ccfg_circ;
            amp_min=best_circ.amp_min, amp_range=best_circ.amp_range)
        t_shift_circ = ts[1]
    end

    # Use elliptical t_shift as reference for 1D alignment
    t_shift = t_shift_ell
    x_1d_t = x_1d .+ t_shift

    # ── 1D fit prediction ──
    y1d_pred = GaussianFit1D.predict_fit(x_1d, best1d, cfg_1d)
    y1d_resid = y_1d .- y1d_pred

    # ── Top-left: elliptical 2D heatmap + ellipses ──
    p_tl, pred_ell = make_2d_heatmap(
        ctx_ell.xs, ctx_ell.ys, ctx_ell.zimg, ctx_ell.mask,
        best_ell, ctx_ell, ccfg_ell, ctx_ell.axisctx;
        title_str=@sprintf("2D elliptical N=%d", n_ell))

    # ── Top-center: circular 2D heatmap + circular ellipses ──
    p_tc, pred_circ = make_2d_heatmap(
        ctx_circ.xs, ctx_circ.ys, ctx_circ.zimg, ctx_circ.mask,
        best_circ, ctx_circ, ccfg_circ, ctx_circ.axisctx;
        title_str=@sprintf("2D circular N=%d", n_circ))

    # ── Top-right: 1D fit with multi-Gaussian components ──
    p_tr = plot(x_1d_t, y_1d; color=:gray, alpha=0.7, label="1D data", linewidth=1)
    plot!(p_tr, x_1d_t, y1d_pred; color=:red, label="1D fit N=$n1d", linewidth=2)
    centers = GaussianFit1D._params_to_centers(best1d.popt, n1d)
    comp_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :lime, :teal, :gold]
    asymmetric = cfg_1d.asymmetric_edges && n1d >= 2
    y0 = best1d.popt[1]
    for (i, c) in enumerate(centers)
        idx = i - 1
        A = GaussianFit1D._get_amplitude(best1d.popt, idx)
        σ_in = GaussianFit1D._get_sigma(best1d.popt, idx)
        if asymmetric && (idx == 0 || idx == n1d - 1)
            σ_out = idx == 0 ? best1d.popt[end-1] : best1d.popt[end]
            z = x_1d .- c
            s = idx == 0 ? (z .< 0) .* σ_out .+ (z .>= 0) .* σ_in :
                           (z .< 0) .* σ_in .+ (z .>= 0) .* σ_out
            y_comp = y0 .+ A .* exp.(-0.5 .* (z ./ s).^2)
        else
            y_comp = y0 .+ A .* exp.(-0.5 .* ((x_1d .- c) ./ max(σ_in, 1e-9)).^2)
        end
        col = comp_colors[mod1(i, length(comp_colors))]
        plot!(p_tr, x_1d_t, y_comp; color=col, alpha=0.35, linestyle=:dash, linewidth=1, label="")
    end
    xlabel!("position (nm)"); ylabel!(p_tr, "intensity")
    title!(p_tr, @sprintf("1D fit N=%d  sBIC=%.0f", n1d, best1d.student_bic))

    # ── Bottom-left: elliptical residuals ──
    p_bl = make_2d_residuals(
        ctx_ell.xs, ctx_ell.ys, ctx_ell.zimg, pred_ell, ctx_ell.mask, ctx_ell.noise, ctx_ell.axisctx;
        title_str=@sprintf("Elliptical residuals / noise (σ_∥=%.3f σ_⊥=%.3f)", best_ell.sigma_parallel_nm, best_ell.sigma_perp_nm))

    # ── Bottom-center: circular residuals ──
    p_bc = make_2d_residuals(
        ctx_circ.xs, ctx_circ.ys, ctx_circ.zimg, pred_circ, ctx_circ.mask, ctx_circ.noise, ctx_circ.axisctx;
        title_str=@sprintf("Circular residuals / noise (σ=%.3f)", best_circ.sigma_parallel_nm))

    # ── Bottom-right: 1D residuals ──
    p_br = plot(x_1d_t, y1d_resid; color=:red, label="1D resid", linewidth=1)
    hline!(p_br, [0]; color=:gray, linestyle=:dash, label="")
    xlabel!("position (nm)"); ylabel!(p_br, "residual")
    title!(p_br, @sprintf("1D residuals  σ=%.5f", std(y1d_resid)))

    # ── Global title ──
    title_str = @sprintf("N_ell=%d (BIC=%.0f)  vs  N_circ=%d (BIC=%.0f)  vs  N_1d=%d (sBIC=%.0f)",
                         n_ell, best_ell.bic, n_circ, best_circ.bic, n1d, best1d.student_bic)

    # ── Assemble 3×2 layout ──
    fig = plot(p_tl, p_tc, p_tr, p_bl, p_bc, p_br;
               layout=(2, 3), size=(2400, 1400),
               plot_title=title_str, plot_titlefontsize=12)

    outpath = joinpath(output_dir, "compare_circular_elliptical.png")
    savefig(fig, outpath)
    println("\n━━━ Saved: $outpath ━━━")

    # ── Print summary ──
    println()
    println("="^70)
    println("SUMMARY")
    println("="^70)
    println("Elliptical:  N=$n_ell  BIC=$(round(best_ell.bic, digits=1))  spar=$(round(best_ell.sigma_parallel_nm, digits=3))  sperp=$(round(best_ell.sigma_perp_nm, digits=3))")
    println("Circular:    N=$n_circ  BIC=$(round(best_circ.bic, digits=1))  sigma=$(round(best_circ.sigma_parallel_nm, digits=3))")
    println("1D slide:    N=$n1d  sBIC=$(round(best1d.student_bic, digits=1))")
    dbic = best_circ.bic - best_ell.bic
    println("ΔBIC(circ - ell) = $(round(dbic, digits=1))")
    if dbic < -10
        println("→ Circular model is strongly favoured (ΔBIC < -10)")
    elseif dbic < -2
        println("→ Circular model is moderately favoured (ΔBIC < -2)")
    elseif dbic > 10
        println("→ Elliptical model is strongly favoured (ΔBIC > 10)")
    elseif dbic > 2
        println("→ Elliptical model is moderately favoured (ΔBIC > 2)")
    else
        println("→ Models are comparable (-2 ≤ ΔBIC ≤ 2)")
    end
    println()
    println("Output: $outpath")
end

main()
