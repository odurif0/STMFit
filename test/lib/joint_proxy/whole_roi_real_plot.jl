#!/usr/bin/env julia

using Printf
using Statistics
using Plots

const _REAL_RESIDUAL_CMAP = cgrad([:blue, :lightgray, :red])

function _real_masked(values, mask)
    result = fill(NaN, size(values))
    result[mask] .= values[mask]
    return result
end

function _real_heat(xs, ys, values, title; color=:thermal, clims=nothing)
    kwargs = (; aspect_ratio=:equal, color, title, xlabel="x (nm)", ylabel="y (nm)",
        colorbar=false)
    return clims === nothing ? heatmap(xs, ys, values; kwargs...) :
        heatmap(xs, ys, values; kwargs..., clims)
end

function _real_overlay!(panel, centers, theta)
    scatter!(panel, first.(centers), last.(centers); marker=:cross, markersize=7,
        markerstrokewidth=2, color=:cyan, label="")
    origin = (sum(first, centers) / length(centers),
              sum(last, centers) / length(centers))
    axis = (cos(theta), sin(theta))
    positions = [((x-origin[1])*axis[1] + (y-origin[2])*axis[2]) for (x,y) in centers]
    lo, hi = extrema(positions)
    plot!(panel, origin[1] .+ [lo,hi].*axis[1], origin[2] .+ [lo,hi].*axis[2];
        color=:yellow, linewidth=1.5, label="")
    return panel
end

function plot_frozen_real_diagnostic(filename, case, set, bundle, outpath;
                                     half_nm::Real, step_nm::Real)
    reg = bundle.nominal.registration
    score = bundle.nominal.scoring.best
    frozen = frozen_scorer_input(reg; theta_base=case.theta_base)
    mask = frozen.mask
    common_image = Matrix(reg.common_image)
    registered_common = reg.fit.background .+
        reg.fit.coefficients[1] .* case.backbone .+
        reg.fit.coefficients[2] .* common_image
    common_component = reg.fit.coefficients[2] .* common_image
    common_residual = case.observed .- registered_common
    contrast = Main._assemble_contrast_frozen(set, score.sequence, case.centers,
        case.xs, case.ys; half_nm=Float64(half_nm), step_nm=Float64(step_nm),
        theta_total=frozen.theta_total, direction=frozen.direction,
        phase=frozen.phase, mirror=frozen.mirror,
        shift_t=frozen.shift_t_nm, shift_u=frozen.shift_u_nm)
    contrast = Main._frozen_blur(contrast, frozen.blur_sigma_nm,
        case.xs[2] - case.xs[1])
    contrast_component = score.contrast_coef .* contrast
    final_model = score.background .+ score.backbone_coef .* case.backbone .+
        score.common_coef .* common_image .+ contrast_component
    final_residual = case.observed .- final_model
    improvement = common_residual.^2 .- final_residual.^2

    finite_observed = case.observed[mask]
    data_clims = quantile(finite_observed, [0.01, 0.99])
    residual_limit = max(quantile(abs.(vcat(common_residual[mask], final_residual[mask])), 0.99), eps())
    improvement_limit = max(quantile(abs.(improvement[mask]), 0.99), eps())
    panels = [
        _real_heat(case.xs, case.ys, _real_masked(case.observed, mask), "Observed";
            clims=(data_clims[1], data_clims[2])),
        _real_heat(case.xs, case.ys, _real_masked(case.backbone, mask), "Backbone"),
        _real_heat(case.xs, case.ys, _real_masked(common_component, mask), "Common contribution"; color=:viridis),
        _real_heat(case.xs, case.ys, _real_masked(registered_common, mask), "Registered common fit";
            clims=(data_clims[1], data_clims[2])),
        _real_heat(case.xs, case.ys, _real_masked(common_residual, mask), "Common residual";
            color=_REAL_RESIDUAL_CMAP, clims=(-residual_limit, residual_limit)),
        _real_heat(case.xs, case.ys, _real_masked(contrast_component, mask), "Scaled contrast"; color=:viridis),
        _real_heat(case.xs, case.ys, _real_masked(final_residual, mask), "Final residual";
            color=_REAL_RESIDUAL_CMAP, clims=(-residual_limit, residual_limit)),
        _real_heat(case.xs, case.ys, _real_masked(improvement, mask), "Per-pixel contrast improvement";
            color=_REAL_RESIDUAL_CMAP, clims=(-improvement_limit, improvement_limit)),
    ]
    foreach(panel -> _real_overlay!(panel, case.centers, case.theta_base), panels[1:4])
    sequence = bundle.final_abstain ? "?" : join(score.sequence)
    title = @sprintf("%s — fused fit-mask — seq=%s state=%d/%d/%d gain=%.5g abstain=%s",
        filename, sequence, reg.direction, reg.phase, reg.mirror,
        bundle.nominal.scoring.incremental_contrast_gain, bundle.final_abstain)
    figure = plot(panels...; layout=(2,4), size=(2200,1000), dpi=170,
        plot_title=title, plot_titlefontsize=9, margin=2Plots.mm)
    savefig(figure, outpath)
    return outpath
end
