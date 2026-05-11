#!/usr/bin/env julia
# Compare 1D slide extraction modes on a single file:
#   :mean (strip), :ridge (crest max), :ridge_mean (crest + disc avg)
# Output: 4-panel plot (2D map + 3 slide profiles) + summary.
# Usage: julia --project=. scripts/compare_slide_modes.jl <filepath> [output_dir]

using STMMolecularFit, GaussianFit1D, GaussianFit2D
using Plots, Printf, Statistics

const FP = length(ARGS) >= 1 ? ARGS[1] : error("Usage: compare_slide_modes.jl <filepath> [output_dir]")
const OUTDIR = length(ARGS) >= 2 ? ARGS[2] : "results/slide_modes/$(replace(basename(FP), r"\.sxm$"i=>""))"
mkpath(OUTDIR)

# Helper
_my_bilinear(xs, ys, z, x0, y0) = begin
    ix = clamp(searchsortedlast(xs, x0), 1, length(xs)-1)
    iy = clamp(searchsortedlast(ys, y0), 1, length(ys)-1)
    x1, x2 = xs[ix], xs[ix+1]; y1, y2 = ys[iy], ys[iy+1]
    tx = (x0-x1)/max(x2-x1, 1e-12); ty = (y0-y1)/max(y2-y1, 1e-12)
    (1-tx)*(1-ty)*z[iy,ix] + tx*(1-ty)*z[iy,ix+1] + (1-tx)*ty*z[iy+1,ix] + tx*ty*z[iy+1,ix+1]
end

const MODES = [(:mean,       0.30, "mean (±0.15 nm)"),
               (:ridge,      0.70, "ridge  (±0.35 nm, max)"),
               (:ridge_mean, 0.70, "ridge mean (±0.35 nm, avg)")]
const COLORS = Dict(:mean => :steelblue, :ridge => :darkgreen, :ridge_mean => :darkorange)

img = STMMolecularFit.read_sxm(FP)
img2 = GaussianFit2D.read_sxm(FP)
results_1d = []

for (mode, width, label) in MODES
    println("Extracting slide: $(mode) width=$(width)...")
    scfg = STMMolecularFit.SlideConfig(width_nm=width, slide_mode=mode,
        support_threshold_fraction=0.20, support_noise_k=2.5, support_padding_nm=0.20,
        output_dir=OUTDIR, no_plot=true)
    slide = STMMolecularFit.extract_slide(img, scfg)
    fcfg = STMMolecularFit.FitSlideConfig(min_spacing=0.35, max_spacing=0.75, output_dir=OUTDIR)
    fit = STMMolecularFit.fit_slide(slide, fcfg)
    best = GaussianFit1D.best_result(fit.fit_run)
    push!(results_1d, (mode=mode, label=label, slide=slide, fit=fit, best=best))
    @printf("  → N=%d  ΔsBIC=0 (raw=%.0f)  support=%.2f nm\n", best.n_peaks, best.student_bic, slide.support_length_nm)
end

# ── 2D heatmap panel (top-left) ──
pcfg2d = GaussianFit2D.PatternConfig(filepath=FP, channel="Z", direction="fwd", stride=2,
    flatten="plane+rows", smooth_radius_px=2, output_dir=OUTDIR, no_plot=false)
xs2d, ys2d, zimg2d, mask2d, _, _, _, _ = GaussianFit2D._fused_roi_data(img2, pcfg2d)
z_clims = (Statistics.quantile(vec(zimg2d), 0.10), Statistics.quantile(vec(zimg2d), 0.995))
p_2d = heatmap(xs2d, ys2d, zimg2d; aspect_ratio=:equal, colormap=:thermal, clims=z_clims,
               title="2D image", xlabel="x (nm)", ylabel="y (nm)", colorbar=false, legend=:topright, legendfontsize=6)
contour!(p_2d, xs2d, ys2d, Float64.(mask2d); levels=[0.5], color=:white, linewidth=2, label="2D ROI mask")
# Slide mask: contour line + scatter (both for visibility)
xs_s, ys_s, slide_msk = STMMolecularFit.molecule_roi_mask(img, STMMolecularFit.PreprocessConfig(roi_threshold_fraction=0.35, stride=2))
sm = Float64.(slide_msk)
plot!(p_2d, xs_s, ys_s, sm; seriestype=:contour, levels=[0.5],
      color=:orangered, linewidth=1.5, label="slide mask")
bx2, by2 = Float64[], Float64[]
for iy in 2:size(slide_msk,1)-1, ix in 2:size(slide_msk,2)-1
    if slide_msk[iy,ix] && (!slide_msk[iy-1,ix] || !slide_msk[iy+1,ix] || !slide_msk[iy,ix-1] || !slide_msk[iy,ix+1])
        push!(bx2, xs_s[ix]); push!(by2, ys_s[iy])
    end
end
scatter!(p_2d, bx2, by2; markersize=2, color=:orangered, alpha=0.6,
         markerstrokewidth=0, label="")

# Draw straight axis and ridge path (use full line bounds, not cropped slide support)
ox, oy = results_1d[1].slide.origin
ax, ay = results_1d[1].slide.axis
perp = results_1d[1].slide.perp
tmin_img, tmax_img = STMMolecularFit._line_bounds(xs2d, ys2d, (ox,oy), (ax,ay))
# Straight axis (yellow)
t_axis = range(tmin_img, tmax_img, length=200)
plot!(p_2d, [ox + t*ax for t in t_axis], [oy + t*ay for t in t_axis];
      color=:yellow, linewidth=2, label="straight axis")

# Ridge path (green, crest tracing with smoothing)
halfw_ridge = 0.35
rx, ry = Float64[], Float64[]
for t in range(tmin_img, tmax_img, length=200)
    bv, bx, by = -Inf, ox + t*ax, oy + t*ay
    for u in range(-halfw_ridge, halfw_ridge, step=0.02)
        x0 = ox + t*ax + u*perp[1]; y0 = oy + t*ay + u*perp[2]
        minimum(xs2d) <= x0 <= maximum(xs2d) && minimum(ys2d) <= y0 <= maximum(ys2d) || continue
        v = _my_bilinear(xs2d, ys2d, zimg2d, x0, y0)
        if v > bv; bv = v; bx = x0; by = y0; end
    end
    push!(rx, bx); push!(ry, by)
end
# Smooth ridge to prevent zigzag
rad = max(2, length(rx) ÷ 25)
srx = [mean(rx[max(1,i-rad):min(end,i+rad)]) for i in 1:length(rx)]
sry = [mean(ry[max(1,i-rad):min(end,i+rad)]) for i in 1:length(ry)]
plot!(p_2d, srx, sry; color=:lime, linewidth=2.5, label="ridge path")

# ── 1D profile panels (3 rows) ──
plots_1d = []
all_ymax = max([maximum(r.fit.fit_run.y) for r in results_1d]...)
all_ymin = min([minimum(r.fit.fit_run.y) for r in results_1d]...)
yl_common = (all_ymin * 0.9, all_ymax * 1.15)

for r in results_1d
    x1d, y1d = r.fit.fit_run.x, r.fit.fit_run.y
    y_pred = GaussianFit1D.predict_fit(x1d, r.best, r.fit.fit_run.cfg)
    color = COLORS[r.mode]
    p = plot(x1d, y1d; color=:gray, alpha=0.6, label="data", linewidth=1.2)
    plot!(p, x1d, y_pred; color=color, linewidth=2.5,
          label="fit N=$(r.best.n_peaks) ΔsBIC=0 (raw=$(round(r.best.student_bic,digits=0)))")
    n_comp = min(r.best.n_peaks, 12)
    centers = GaussianFit1D._params_to_centers(r.best.popt, r.best.n_peaks)
    cfg1d = r.fit.fit_run.cfg
    asymmetric = cfg1d.asymmetric_edges && n_comp >= 2
    y0 = r.best.popt[1]
    for j in 1:n_comp
        idx = j - 1
        A = GaussianFit1D._get_amplitude(r.best.popt, idx)
        σ_in = GaussianFit1D._get_sigma(r.best.popt, idx)
        if asymmetric && (idx == 0 || idx == n_comp - 1)
            σ_out = idx == 0 ? r.best.popt[end-1] : r.best.popt[end]
            z = x1d .- centers[j]
            s = idx == 0 ? (z .< 0) .* σ_out .+ (z .>= 0) .* σ_in :
                           (z .< 0) .* σ_in .+ (z .>= 0) .* σ_out
            y_comp = y0 .+ A .* exp.(-0.5 .* (z ./ s).^2)
        else
            y_comp = y0 .+ A .* exp.(-0.5 .* ((x1d .- centers[j]) ./ max(σ_in, 1e-9)).^2)
        end
        plot!(p, x1d, y_comp; color=color, alpha=0.2, linestyle=:dash, linewidth=0.8, label="")
    end
    ylims!(p, yl_common...)
    title!(p, "$(r.label) — support=$(round(r.slide.support_length_nm,digits=2)) nm")
    xlabel!("position (nm)"); ylabel!(p, "intensity")
    push!(plots_1d, p)
end

# ── Assemble 4-panel layout ──
# ── Layout: 2D left (1 col × full height), 1D right (3 stacked rows) ──
right_panel = plot(plots_1d...; layout=(3,1))
l = @layout [a{0.55w} b{0.45w}]
fig = plot(p_2d, right_panel; layout=l, size=(1700, 950),
           plot_title="$(basename(FP)) — slide extraction modes")
savefig(fig, joinpath(OUTDIR, "slide_modes_comparison.png"))

# ── Summary ──
println("\n" * "="^50)
println("Summary — $(basename(FP))")
println("="^50)
for r in results_1d
    @printf("  %-6s  N=%d  ΔsBIC=0 (raw=%.0f)  support=%.2f nm  noise=%.5f\n",
        string(r.mode), r.best.n_peaks, r.best.student_bic, r.slide.support_length_nm, r.slide.noise_1d)
end
println("\nPlot: $(joinpath(OUTDIR, "slide_modes_comparison.png"))")
