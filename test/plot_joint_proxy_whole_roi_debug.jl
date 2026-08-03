#!/usr/bin/env julia

# Diagnostic-only plots for real whole-ROI DFT-mold transfer. This script does
# not alter fitting, N selection, calibration, or production unit assignment.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using Printf
using Statistics
using Plots
using STMSXMIO
using GaussianFit2D

include(joinpath(@__DIR__, "diagnose_joint_proxy_whole_roi_batch.jl"))

const DEBUG_FILES = ["240817_007.sxm", "240817_050.sxm"]
const DEFAULT_SUMMARY = joinpath(@__DIR__, "..", "results",
    "best_plots_240817_primary_rerun", "summary_overlap060_hard.tsv")
const DEFAULT_MOLDS = "/tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv"
const RESIDUAL_COLORMAP = cgrad([:blue, :lightgray, :red])

function _parse_debug_args(args)
    outdir = "/tmp/opencode/whole_roi_debug_plots"
    summary = DEFAULT_SUMMARY
    molds = DEFAULT_MOLDS
    files = String[]
    i = 1
    while i <= length(args)
        if args[i] == "--out-dir"
            outdir = args[i + 1]; i += 2
        elseif args[i] == "--selection-summary"
            summary = args[i + 1]; i += 2
        elseif args[i] == "--molds"
            molds = args[i + 1]; i += 2
        else
            push!(files, String(args[i])); i += 1
        end
    end
    return (; outdir=abspath(outdir), summary=abspath(summary), molds=abspath(molds),
        files=isempty(files) ? DEBUG_FILES : files)
end

function _fit_details(y, backbone, mold)
    predictors = [backbone, mold]
    candidates = NamedTuple[]
    for active in (Int[], [1], [2], [1, 2])
        design = ones(Float64, length(y), length(active) + 1)
        for (column, predictor) in enumerate(active)
            design[:, column + 1] .= predictors[predictor]
        end
        coefficients = design \ y
        any(coefficients[2:end] .< -sqrt(eps(Float64))) && continue
        coefficients[2:end] .= max.(coefficients[2:end], 0.0)
        full = zeros(2)
        for (column, predictor) in enumerate(active)
            full[predictor] = coefficients[column + 1]
        end
        residual = y .- design * coefficients
        push!(candidates, (; sse=sum(abs2, residual), background=coefficients[1],
            backbone_coefficient=full[1], mold_coefficient=full[2]))
    end
    return candidates[argmin(getproperty.(candidates, :sse))]
end

function _null_details(y, backbone)
    zero_mold = zeros(length(y))
    fit = _fit_details(y, backbone, zero_mold)
    return (; sse=fit.sse, background=fit.background,
        backbone_coefficient=fit.backbone_coefficient)
end

function _search_variant(entry, observed, backbone, mask, centers, xs, ys, theta;
                         half_nm=0.80, step_nm=0.08)
    good = mask .& isfinite.(observed) .& isfinite.(backbone)
    count(good) >= 20 || error("diagnostic mask has fewer than 20 finite pixels")
    y, bv = observed[good], backbone[good]
    null_fit = _null_details(y, bv)
    best = nothing
    sequence_sse = Dict{String,Float64}()
    n = length(centers)
    for direction in (0, 1), phase in (0, 1), mirror in (0, 1)
        images = precompute_lobe_images(entry, centers, n, xs, ys, theta,
            half_nm, step_nm, direction, phase, mirror)
        for bits in 0:(2^n - 1)
            mold = zeros(Float64, size(observed))
            for lobe in 1:n
                mold .+= ((bits >> (lobe - 1)) & 1 == 0) ? images[lobe][1] : images[lobe][2]
            end
            fit = _fit_details(y, bv, mold[good])
            sequence = join([Int((bits >> (lobe - 1)) & 1) for lobe in 1:n], "")
            sequence_sse[sequence] = min(get(sequence_sse, sequence, Inf), fit.sse)
            if best === nothing || fit.sse < best.fit.sse
                best = (; sequence, direction, phase, mirror, fit, mold)
            end
        end
    end
    centered = best.mold .- mean(best.mold[good])
    centered_fit = _fit_details(y, bv, centered[good])
    ranked = sort(collect(sequence_sse); by=last)
    sequence_gap_pct = length(ranked) > 1 ?
        100 * (last(ranked[2]) - last(ranked[1])) / null_fit.sse : Inf
    return (; good, null_fit, best, centered_fit,
        sse_reduction_pct=100 * (1 - best.fit.sse / null_fit.sse), sequence_gap_pct)
end

function _masked(values, mask)
    result = fill(NaN, size(values))
    result[mask] .= values[mask]
    return result
end

function _overlay!(panel, centers, axisctx)
    scatter!(panel, first.(centers), last.(centers); marker=:cross, markersize=7,
        markerstrokewidth=2, color=:cyan, label="")
    ts = [((x - axisctx.origin[1]) * axisctx.axis[1] +
           (y - axisctx.origin[2]) * axisctx.axis[2]) for (x, y) in centers]
    lo, hi = extrema(ts)
    xline = axisctx.origin[1] .+ [lo, hi] .* axisctx.axis[1]
    yline = axisctx.origin[2] .+ [lo, hi] .* axisctx.axis[2]
    plot!(panel, xline, yline; color=:yellow, linewidth=1.5, label="")
    return panel
end

function _heat(xs, ys, values, title; color=:thermal, clims=nothing)
    kwargs = (; aspect_ratio=:equal, color, title, xlabel="x (nm)", ylabel="y (nm)",
        colorbar=false)
    return clims === nothing ? heatmap(xs, ys, values; kwargs...) :
        heatmap(xs, ys, values; kwargs..., clims)
end

function _plot_variant(filename, variant_name, observed, backbone, result,
                       centers, xs, ys, axisctx, outpath)
    good, best, null_fit = result.good, result.best, result.null_fit
    null_model = null_fit.background .+ null_fit.backbone_coefficient .* backbone
    mold_component = best.fit.mold_coefficient .* best.mold
    combined = best.fit.background .+ best.fit.backbone_coefficient .* backbone .+ mold_component
    null_residual = observed .- null_model
    combined_residual = observed .- combined
    improvement = null_residual .^ 2 .- combined_residual .^ 2
    finite_observed = observed[good]
    data_clims = quantile(finite_observed, [0.01, 0.99])
    residual_limit = max(quantile(abs.(vcat(null_residual[good], combined_residual[good])), 0.99), eps())
    improvement_limit = max(quantile(abs.(improvement[good]), 0.99), eps())

    top = [
        _heat(xs, ys, _masked(observed, good), "Observed"; clims=(data_clims[1], data_clims[2])),
        _heat(xs, ys, _masked(null_model, good), "Backbone-only fit"; clims=(data_clims[1], data_clims[2])),
        _heat(xs, ys, _masked(mold_component, good), "Scaled mold contribution"; color=:viridis),
        _heat(xs, ys, _masked(combined, good), "Combined fit"; clims=(data_clims[1], data_clims[2])),
    ]
    foreach(panel -> _overlay!(panel, centers, axisctx), top)
    bottom = [
        _heat(xs, ys, _masked(null_residual, good), "Null residual";
            color=RESIDUAL_COLORMAP, clims=(-residual_limit, residual_limit)),
        _heat(xs, ys, _masked(combined_residual, good), "Combined residual";
            color=RESIDUAL_COLORMAP, clims=(-residual_limit, residual_limit)),
        _heat(xs, ys, _masked(best.mold, good), "Raw assembled mold"; color=:viridis),
        _heat(xs, ys, _masked(improvement, good), "Per-pixel SSE improvement";
            color=RESIDUAL_COLORMAP, clims=(-improvement_limit, improvement_limit)),
    ]
    title = @sprintf("%s — %s — seq=%s state=%d/%d/%d SSE↓=%.3f%% mold=%.5g",
        filename, variant_name, best.sequence, best.direction, best.phase, best.mirror,
        result.sse_reduction_pct, best.fit.mold_coefficient)
    figure = plot(top..., bottom...; layout=(2, 4), size=(2200, 1000), dpi=170,
        plot_title=title, plot_titlefontsize=9, margin=2Plots.mm)
    savefig(figure, outpath)
    return outpath
end

function _prepare_file(filename, selected_n, data_dir, entry, outdir)
    path = joinpath(data_dir, filename)
    isfile(path) || error("SXM file not found: $path")
    image = STMSXMIO.read_sxm(path)
    pattern = Inference._default_pattern_config(path, outdir)
    chain = Inference._default_chain_config()
    results, _, context = redirect_stdout(devnull) do
        GaussianFit2D.chain_gaussian_sweep(image, pattern, chain)
    end
    views = Inference._build_views(image, pattern,
        Main.JointProxyCandidateViews.build_view_data)
    report = Main.JointProxyCandidateViews.extract_candidate_views(
        results, context, chain, views; patch_half_nm=0.32, patch_step_nm=0.08)
    candidate_index = findfirst(item -> item.n == selected_n, report.candidates)
    candidate_index === nothing && error("selected N=$selected_n missing for $filename")
    candidate = report.candidates[candidate_index]
    best_index = findfirst(r -> r.n == selected_n && r.success && r.valid, results)
    best_index === nothing && error("valid selected fit missing for $filename")
    best_result = results[best_index]
    fwd = only(v for v in views if v.label == "fwd" && v.present)
    fwd.xs == context.xs && fwd.ys == context.ys || error("forward/fused grids differ")
    centers = [(lobe.x_nm, lobe.y_nm) for lobe in candidate.lobes]
    margin = 0.90
    ix = findall(x -> minimum(first, centers)-margin <= x <= maximum(first, centers)+margin, context.xs)
    iy = findall(y -> minimum(last, centers)-margin <= y <= maximum(last, centers)+margin, context.ys)
    xs, ys = context.xs[ix], context.ys[iy]
    forward = fwd.z[iy, ix]
    fused = context.zimg[iy, ix]
    model_values = GaussianFit2D._chain_model_values(
        vec(repeat(reshape(xs, 1, :), length(ys), 1)),
        vec(repeat(reshape(ys, :, 1), 1, length(xs))),
        best_result.params, selected_n, context.axisctx, chain;
        amp_min=best_result.amp_min, amp_range=best_result.amp_range)
    backbone = reshape(model_values, length(ys), length(xs))
    finite_mask = isfinite.(forward) .& isfinite.(fused) .& isfinite.(backbone)
    xgrid = repeat(reshape(xs, 1, :), length(ys), 1)
    ygrid = repeat(reshape(ys, :, 1), 1, length(xs))
    t, u = GaussianFit2D._chain_coordinates(xgrid, ygrid, context.axisctx)
    fit_mask = finite_mask .& (abs.(u) .<= context.fit_width_nm) .&
        (t .>= context.axisctx.tmin) .& (t .<= context.axisctx.tmax)
    theta = atan(context.axisctx.axis[2], context.axisctx.axis[1])
    variants = [("forward_full", forward, finite_mask),
                ("forward_fitmask", forward, fit_mask),
                ("fused_full", fused, finite_mask),
                ("fused_fitmask", fused, fit_mask)]
    rows = NamedTuple[]
    stem = splitext(filename)[1]
    for (name, observed, mask) in variants
        result = _search_variant(entry, observed, backbone, mask, centers, xs, ys, theta)
        png = joinpath(outdir, "$(stem)_$(name).png")
        _plot_variant(filename, name, observed, backbone, result, centers, xs, ys,
            context.axisctx, png)
        centered_delta = result.centered_fit.sse - result.best.fit.sse
        push!(rows, (; file=filename, variant=name, n=selected_n, pixels=count(result.good),
            reduction=result.sse_reduction_pct, sequence=result.best.sequence,
            sequence_gap=result.sequence_gap_pct,
            direction=result.best.direction, phase=result.best.phase, mirror=result.best.mirror,
            backbone_coeff=result.best.fit.backbone_coefficient,
            mold_coeff=result.best.fit.mold_coefficient,
            mold_mean=mean(result.best.mold[result.good]),
            mold_std=std(result.best.mold[result.good]), centered_sse_delta=centered_delta,
            png))
        println("  ", filename, " ", name, ": SSE↓=",
            @sprintf("%.3f%%", result.sse_reduction_pct), " seq=", result.best.sequence,
            " mold=", @sprintf("%.5g", result.best.fit.mold_coefficient))
    end
    return rows
end

function main(args)
    opts = _parse_debug_args(args)
    mkpath(opts.outdir)
    selected = read_selected_counts(opts.summary)
    entry = load_converged_entry(opts.molds)
    data_dir = something(get(ENV, "STMFIT_DATA_DIR", nothing), ".")
    rows = NamedTuple[]
    for file in opts.files
        haskey(selected, file) || error("file missing from selection summary: $file")
        append!(rows, _prepare_file(file, selected[file], data_dir, entry, opts.outdir))
    end
    table = joinpath(opts.outdir, "comparison.tsv")
    open(table, "w") do io
        println(io, "file\tvariant\tn\tpixels\tsse_reduction_pct\tsequence_gap_pct\tbest_sequence\tdirection\tphase\tmirror\tbackbone_coeff\tmold_coeff\tmold_mean\tmold_std\tcentered_sse_delta\tpng")
        for r in rows
            println(io, join((r.file, r.variant, r.n, r.pixels,
                @sprintf("%.8f", r.reduction), @sprintf("%.8f", r.sequence_gap),
                r.sequence, r.direction, r.phase, r.mirror,
                @sprintf("%.8g", r.backbone_coeff), @sprintf("%.8g", r.mold_coeff),
                @sprintf("%.8g", r.mold_mean), @sprintf("%.8g", r.mold_std),
                @sprintf("%.8g", r.centered_sse_delta), r.png), '\t'))
        end
    end
    println("wrote ", table)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
