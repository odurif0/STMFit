"""
Visualization module for multi-Gaussian fitting results.

Uses Plots.jl with GR backend for publication-quality 3-panel figures.
"""

using Plots
gr()  # Fast, non-interactive GR backend

export plot_results

function _plot_single_model(x, y, result::FitResult, all_results, cfg::FitConfig, suffix)
    """Three-panel figure: fit, residuals, BIC comparison."""
    n_peaks = result.n_peaks
    popt = result.popt
    x_unit = cfg.x_unit
    asymmetric = cfg.asymmetric_edges && n_peaks >= 2

    n_colors = max(n_peaks, 1)
    peak_colors = Plots.palette(:tab10, n_colors)

    xx = range(minimum(x), maximum(x), length=cfg.fine_grid_points)
    y_total = predict_fit(xx, result, cfg)

    y0 = popt[1]
    centers = _params_to_centers(popt, n_peaks)

    # Panel 1: Data + Fit
    p1 = plot(x, y; seriestype=:scatter, label="Data",
              markercolor=:gray, markersize=3, markerstrokewidth=0,
              alpha=0.6, legend=:topright, legendfontsize=7)
    plot!(p1, xx, y_total; linewidth=2, color=:black, label="Total fit")
    hline!(p1, [y0]; linestyle=:dot, linewidth=1, color=:gray, alpha=0.5,
           label="Baseline y0=$(round(y0, digits=3))")

    for i in 0:(n_peaks - 1)
        A = _get_amplitude(popt, i; use_log=cfg.use_log_amplitude)
        sigma = _get_sigma(popt, i)
        mu = centers[i+1]
        fwhm = FWHM_TO_SIGMA * sigma

        if asymmetric && (i == 0 || i == n_peaks - 1)
            sigma_outer = i == 0 ? popt[end-1] : popt[end]
            z = xx .- mu
            s = i == 0 ? (z .< 0) .* sigma_outer .+ (z .>= 0) .* sigma :
                         (z .< 0) .* sigma .+ (z .>= 0) .* sigma_outer
            g = y0 .+ A .* exp.(-0.5 .* (z ./ s).^2)
        else
            g = y0 .+ A .* exp.(-0.5 .* ((xx .- mu) ./ sigma).^2)
        end

        color = peak_colors[i+1]
        plot!(p1, xx, g; linestyle=:dash, linewidth=1.2, color=color,
              label="Peak $(i+1): $(round(mu, digits=2)) $x_unit, FWHM=$(round(fwhm, digits=3))")
        plot!(p1, xx, g; fillrange=y0, fillalpha=0.12, color=color, label="")
    end

    xlabel!(p1, "Distance ($x_unit)")
    ylabel!(p1, "Intensity")

    # Panel 2: Residuals
    residuals = y - predict_fit(x, result, cfg)
    p2 = plot(x, residuals; seriestype=:scatter, label="Residuals",
              markercolor=:gray, markersize=3, markerstrokewidth=0, alpha=0.6)
    hline!(p2, [0.0]; linewidth=0.5, color=:black, label="")
    xlabel!(p2, "Distance ($x_unit)")
    ylabel!(p2, "Residuals")

    # Panel 3: BIC vs n_peaks
    n_list = [r.n_peaks for r in all_results]
    bic_label = cfg.use_student_bic ? "sBIC" : "BIC"
    bic_list = cfg.use_student_bic ? [r.student_bic for r in all_results] : [r.bic for r in all_results]
    best_idx = argmin(bic_list)
    best_bic_val = bic_list[best_idx]
    bic_threshold = cfg.bic_competition_threshold

    # Color-code
    bar_colors = []
    for (i, bic_val) in enumerate(bic_list)
        delta = bic_val - best_bic_val
        if i == best_idx
            push!(bar_colors, :red)
        elseif delta <= bic_threshold
            push!(bar_colors, :gold)
        else
            push!(bar_colors, :lightgray)
        end
    end

    p3 = bar(n_list, bic_list; color=bar_colors, label=bic_label,
             legend=:topleft, legendfontsize=8, bar_width=0.7)
    this_bic = cfg.use_student_bic ? result.student_bic : result.bic
    scatter!(p3, [n_list[best_idx]], [bic_list[best_idx]];
             markershape=:star, markersize=12, color=:red,
             label="Best: $(n_list[best_idx]) peaks")
    scatter!(p3, [n_peaks], [this_bic];
             markershape=:diamond, markersize=8, color=:orange,
             label="This: $n_peaks peaks")

    # Threshold line
    hline!(p3, [best_bic_val + bic_threshold];
           linestyle=:dash, linewidth=1, color=:orange,
           label="Δ$(bic_label) ≤ $(round(bic_threshold, digits=0))")

    xlabel!(p3, "Number of Gaussians")
    ylabel!(p3, bic_label)

    # Combine into a single figure
    bic_val = cfg.use_student_bic ? result.student_bic : result.bic
    title_text = "Multi-Gaussian Fit ($n_peaks peaks, $bic_label = $(round(bic_val, digits=1)), " *
                 "R² = $(round(result.r_squared, digits=4)))"

    fig = plot(p1, p2, p3;
               layout=(3, 1), size=(cfg.fig_width * 100, cfg.fig_height * 100),
               plot_title=title_text, plot_titlefontsize=12)

    plot_file = output_path(cfg, suffix)
    savefig(fig, plot_file)
    return plot_file
end

function plot_results(x, y, best_r, all_results, cfg)
    """Generate plots for all models."""
    plot_files = String[]
    overall_best = cfg.use_student_bic ? argmin(r -> r.student_bic, all_results) : argmin(r -> r.bic, all_results)

    for result in all_results
        n = result.n_peaks
        tag = result === overall_best ? "_best" : ""
        suffix = "_n$(n)$(tag).png"
        path = _plot_single_model(x, y, result, all_results, cfg, suffix)
        result.plot_file = path
        push!(plot_files, path)
        label = result === overall_best ? " (BEST)" : ""
        println("  Plot n=$n saved to: $path$label")
    end

    println("\n  $(length(plot_files)) plots generated.")
    return plot_files
end
