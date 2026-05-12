# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

function _write_parameters(result::PatternFitResult, cfg::PatternConfig)
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "molecular_features.tsv")
    open(out, "w") do io
        println(io, "type\tfeature\tamplitude\tx_nm\ty_nm\tsigma_x_nm\tsigma_y_nm\tscore")
        for (i, f) in enumerate(result.raw_features)
            println(io, "detected\t$(i)\t$(f.amplitude)\t$(f.x_nm)\t$(f.y_nm)\t$(f.sigma_x_nm)\t$(f.sigma_y_nm)\t$(f.score)")
        end
        for (i, f) in enumerate(result.features)
            println(io, "fitted\t$(i)\t$(f.amplitude)\t$(f.x_nm)\t$(f.y_nm)\t$(f.sigma_x_nm)\t$(f.sigma_y_nm)\t$(f.score)")
        end
    end
    result.parameter_file = out
    estimate_out = joinpath(cfg.output_dir, "chain_estimate.tsv")
    open(estimate_out, "w") do io
        println(io, "roi_length_nm\trepeat_spacing_nm\tend_extension_nm\testimated_repeats\testimated_range_min\testimated_range_max\tvisible_chain_features")
        println(io, "$(result.roi_length_nm)\t$(cfg.repeat_spacing_nm)\t$(cfg.end_extension_nm)\t$(result.estimated_repeats)\t$(result.estimated_repeat_range[1])\t$(result.estimated_repeat_range[2])\t$(isempty(result.chains) ? 0 : length(result.chains[1].features))")
    end
    chain_out = joinpath(cfg.output_dir, "molecular_chains.tsv")
    open(chain_out, "w") do io
        println(io, "chain_id\tindex_in_chain\tx_nm\ty_nm\tscore\tmean_spacing_nm\tspacing_cv\tmax_turn_angle_deg\tchain_score")
        for chain in result.chains
            for (i, f) in enumerate(chain.features)
                println(io, "$(chain.id)\t$(i)\t$(f.x_nm)\t$(f.y_nm)\t$(f.score)\t$(chain.mean_spacing_nm)\t$(chain.spacing_cv)\t$(chain.max_turn_angle_deg)\t$(chain.score)")
            end
        end
    end
    return out
end

_safe_name(s::String) = replace(lowercase(strip(s)), r"[^a-z0-9]+" => "_")

function _write_all_images_summary(rows, cfg::PatternConfig)
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "all_images_summary.tsv")
    open(out, "w") do io
        println(io, "channel\tdirection\tsuccess\tn_detected\tn_fitted\tr_squared\trss\toutput_dir")
        for row in rows
            println(io, "$(row.channel)\t$(row.direction)\t$(row.success)\t$(row.n_detected)\t$(row.n_fitted)\t$(row.r_squared)\t$(row.rss)\t$(row.output_dir)")
        end
    end
    return out
end

function _plot_result(img::SXMImage, cfg::PatternConfig, result::PatternFitResult)
    ch = get_channel(img, cfg.channel; direction=cfg.direction)
    xs, ys, raw, z, z_smooth, unit, _noise = preprocess_channel(img, ch, cfg)
    xflat = repeat(xs, inner=length(ys))
    yflat = repeat(ys, outer=length(xs))
    pred = reshape(_model_values(xflat, yflat, result.params_unconstrained, img, cfg), length(ys), length(xs))

    p1 = heatmap(xs, ys, raw; aspect_ratio=:equal, title="Raw $(ch.name) $(ch.direction) [$unit]", xlabel="x (nm)", ylabel="y (nm)")
    evidence = result.evidence_map
    if cfg.fusion && !isempty(evidence)
        p2 = heatmap(xs, ys, evidence; aspect_ratio=:equal, title="Fused multi-view evidence + chains", xlabel="x (nm)", ylabel="y (nm)")
    else
        p2 = heatmap(xs, ys, z; aspect_ratio=:equal, title="Preprocessed + chains", xlabel="x (nm)", ylabel="y (nm)")
    end
    for f in result.raw_features
        scatter!(p2, [f.x_nm], [f.y_nm]; marker=:circle, markersize=3, label=false, color=:gray)
    end
    chain_colors = [:red, :cyan, :yellow, :lime, :magenta, :orange, :white, :blue]
    for (ci, chain) in enumerate(result.chains)
        color = chain_colors[mod1(ci, length(chain_colors))]
        xs_chain = [f.x_nm for f in chain.features]
        ys_chain = [f.y_nm for f in chain.features]
        plot!(p2, xs_chain, ys_chain; linewidth=2.5, color=color, label=false)
        scatter!(p2, xs_chain, ys_chain; marker=:cross, markersize=10, color=color, linewidth=2.5, label=false)
    end
    p3 = heatmap(xs, ys, pred; aspect_ratio=:equal, title="Chain-constrained Gaussian refinement ($(length(result.features)) fitted)", xlabel="x (nm)", ylabel="y (nm)")
    for f in result.features
        scatter!(p3, [f.x_nm], [f.y_nm]; marker=:cross, markersize=5, label=false, color=:red)
    end
    p4 = heatmap(xs, ys, z .- pred; aspect_ratio=:equal, title="Residuals", xlabel="x (nm)", ylabel="y (nm)")
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1200, 950))
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "molecular_pattern_overview.png")
    savefig(fig, out)
    result.plot_file = out
    return out
end

function _plot_all_images(img::SXMImage, cfg::PatternConfig)
    mkpath(cfg.output_dir)
    panels = Any[]
    for ch in img.channels
        local_cfg = deepcopy(cfg)
        local_cfg.channel = ch.name
        local_cfg.direction = ch.direction
        xs, ys, raw, z, z_smooth, unit, noise = preprocess_channel(img, ch, local_cfg)
        candidates = detect_blobs(z_smooth, xs, ys, local_cfg, noise)
        p_raw = heatmap(xs, ys, raw; aspect_ratio=:equal,
                        title="Raw $(ch.name) $(ch.direction) [$unit]",
                        xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
        p_pre = heatmap(xs, ys, z; aspect_ratio=:equal,
                        title="Preprocessed + candidates ($(length(candidates)))",
                        xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
        for f in candidates
            scatter!(p_pre, [f.x_nm], [f.y_nm]; marker=:circle, markersize=3,
                     label=false, color=:white)
        end
        push!(panels, p_raw, p_pre)
    end
    fig = plot(panels...; layout=(length(img.channels), 2), size=(1100, 280length(img.channels)))
    out = joinpath(cfg.output_dir, "all_images_overview.png")
    savefig(fig, out)
    return out
end

function _run_all_images(img::SXMImage, cfg::PatternConfig)
    rows = NamedTuple[]
    for ch in img.channels
        local_cfg = deepcopy(cfg)
        local_cfg.channel = ch.name
        local_cfg.direction = ch.direction
        local_cfg.output_dir = joinpath(cfg.output_dir, "$(_safe_name(ch.name))_$(ch.direction)")
        result = fit_molecular_pattern(img, local_cfg)
        _write_parameters(result, local_cfg)
        !cfg.no_plot && _plot_result(img, local_cfg, result)
        push!(rows, (channel=ch.name, direction=ch.direction, success=result.success,
                     n_detected=length(result.raw_features), n_fitted=length(result.features), r_squared=result.r_squared,
                     rss=result.rss, output_dir=local_cfg.output_dir))
        println("$(ch.name) $(ch.direction): detected=$(length(result.raw_features)) fitted=$(length(result.features)) R²=$(@sprintf("%.5f", result.r_squared)) dir=$(local_cfg.output_dir)")
    end
    summary = _write_all_images_summary(rows, cfg)
    overview = cfg.no_plot ? nothing : _plot_all_images(img, cfg)
    println("all_images_summary: $summary")
    overview !== nothing && println("all_images_overview: $overview")
end

function _parse_sweep(spec::String)
    parts = split(spec, ':')
    length(parts) == 3 || error("threshold sweep must be start:step:stop, e.g. 2.0:0.25:4.0")
    a, step, b = parse.(Float64, parts)
    step > 0 || error("threshold sweep step must be > 0")
    return collect(a:step:b)
end

function _run_threshold_sweep(img::SXMImage, cfg::PatternConfig, spec::String)
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "threshold_sweep.tsv")
    open(out, "w") do io
        println(io, "threshold_sigma\tn_candidates\tn_chains\tn_chain_features\tbest_chain_length\tbest_chain_score\troi_length_nm\testimated_repeats\testimated_range")
        for thr in _parse_sweep(spec)
            local_cfg = deepcopy(cfg)
            local_cfg.threshold_sigma = thr
            _xs, _ys, _evidence, candidates, chains, accepted, _roi, _axis, roi_len, est_n, est_range = detect_molecular_chains(img, local_cfg)
            best_len = isempty(chains) ? 0 : maximum(length(c.features) for c in chains)
            best_score = isempty(chains) ? -Inf : maximum(c.score for c in chains)
            println(io, "$(thr)\t$(length(candidates))\t$(length(chains))\t$(length(accepted))\t$(best_len)\t$(best_score)\t$(roi_len)\t$(est_n)\t$(est_range[1])-$(est_range[2])")
            println("threshold=$(thr): candidates=$(length(candidates)) chains=$(length(chains)) best_length=$(best_len) estimated_repeats=$(est_n)")
        end
    end
    println("threshold_sweep: $out")
    return out
end

function _write_chain_sweep(results::Vector{ChainModelResult}, best::ChainModelResult, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig)
    mkpath(cfg.output_dir)
    summary = joinpath(cfg.output_dir, "chain_model_selection.tsv")
    open(summary, "w") do io
        println(io, "N\tn_params\tsuccess\tvalid\ttrain_nll\tcv_nll_mean\tcv_nll_std\tbic\taicc\tbic_per_pixel\taicc_per_pixel\trss\tchi2_reduced\tmad\tresidual_peak_snr\tmean_spacing_nm\tspacing_cv\tmax_lateral_nm\tsigma_parallel_nm\tsigma_perp_nm\toverlap\tkappa_max_adj\tendpoint_overrun_nm\tbound_like\treason\twarnings\tselected")
        for r in results
            npx = max(1, length(ctx.z))
            println(io, "$(r.n)\t$(_chain_nparams(r.n, ccfg))\t$(r.success)\t$(r.valid)\t$(r.train_nll)\t$(r.cv_nll_mean)\t$(r.cv_nll_std)\t$(r.bic)\t$(r.aicc)\t$(r.bic/npx)\t$(r.aicc/npx)\t$(r.rss)\t$(r.chi2_reduced)\t$(r.mad)\t$(r.residual_peak_snr)\t$(r.mean_spacing_nm)\t$(r.spacing_cv)\t$(r.max_lateral_nm)\t$(r.sigma_parallel_nm)\t$(r.sigma_perp_nm)\t$(r.overlap)\t$(r.kappa_max_adj)\t$(r.endpoint_overrun_nm)\t$(r.bound_like)\t$(r.reason)\t$(r.valid ? "" : r.reason)\t$(r === best)")
        end
    end
    params = joinpath(cfg.output_dir, "chain_selected_lobes.tsv")
    open(params, "w") do io
        println(io, "lobe\tamplitude\tx_nm\ty_nm\tt_nm\tu_nm\tsigma_parallel_nm\tsigma_perp_nm\tspacing_prev_nm")
        if best.n > 0 && best.success
            _b, feats, ts, us, spars, sperps = _decode_chain(best.params, best.n, ctx.axisctx, ccfg;
                                                              amp_min=best.amp_min, amp_range=best.amp_range)
            for (i, f) in enumerate(feats)
                spacing_prev = i == 1 ? NaN : ts[i] - ts[i-1]
                println(io, "$(i)\t$(f.amplitude)\t$(f.x_nm)\t$(f.y_nm)\t$(ts[i])\t$(us[i])\t$(spars[i])\t$(sperps[i])\t$(spacing_prev)")
            end
        end
    end
    axisfile = joinpath(cfg.output_dir, "chain_axis.tsv")
    open(axisfile, "w") do io
        println(io, "origin_x_nm\torigin_y_nm\taxis_x\taxis_y\tperp_x\tperp_y\ttmin_nm\ttmax_nm")
        println(io, "$(ctx.axisctx.origin[1])\t$(ctx.axisctx.origin[2])\t$(ctx.axisctx.axis[1])\t$(ctx.axisctx.axis[2])\t$(ctx.axisctx.perp[1])\t$(ctx.axisctx.perp[2])\t$(ctx.axisctx.tmin)\t$(ctx.axisctx.tmax)")
    end
    supportfile = joinpath(cfg.output_dir, "chain_support_metadata.tsv")
    open(supportfile, "w") do io
        println(io, "key\tvalue")
        println(io, "pipeline\tauto ROI -> weighted PCA axis -> tube axial profile -> robust support -> adaptive n_max -> deterministic global+local fit -> BIC selection")
        println(io, "n_min_config\t$(ccfg.n_min)")
        println(io, "n_max_config\t$(ccfg.n_max)")
        println(io, "spacing_min_nm\t$(ccfg.spacing_min_nm)")
        println(io, "spacing_min_effective_nm\t$(_effective_spacing_min_nm(ccfg))")
        println(io, "spacing_min_effective_source\tmax(spacing_min_nm, sqrt(-2log(max_overlap))*sigma_max)")
        println(io, "spacing_max_nm\t$(ccfg.spacing_max_nm)")
        println(io, "max_overlap\t$(ccfg.max_overlap)")
        println(io, "kappa_max\t$(ccfg.kappa_max)")
        println(io, "sigma_parallel_min_nm\t$(ccfg.sigma_parallel_min_nm)")
        println(io, "sigma_parallel_max_nm\t$(ccfg.sigma_parallel_max_nm)")
        println(io, "sigma_perp_min_nm\t$(ccfg.sigma_perp_min_nm)")
        println(io, "sigma_perp_max_nm\t$(ccfg.sigma_perp_max_nm)")
        for k in keys(ctx.support_meta)
            println(io, "$(k)\t$(getfield(ctx.support_meta, k))")
        end
    end
    profile = joinpath(cfg.output_dir, "chain_axis_profile.tsv")
    pred = best.success ? (best.n == 0 ? fill(best.params[1], length(ctx.z)) :
        _chain_model_values(ctx.x, ctx.y, best.params, best.n, ctx.axisctx, ccfg;
                            amp_min=best.amp_min, amp_range=best.amp_range)) : fill(NaN, length(ctx.z))
    t = (ctx.x .- ctx.axisctx.origin[1]) .* ctx.axisctx.axis[1] .+ (ctx.y .- ctx.axisctx.origin[2]) .* ctx.axisctx.axis[2]
    nb = max(20, min(200, Int(ceil((maximum(t) - minimum(t)) / max(median(diff(ctx.xs)), EPS)))))
    sdata = zeros(nb); smodel = zeros(nb); sresid = zeros(nb); counts = zeros(Int, nb)
    for i in eachindex(t)
        b = clamp(Int(floor((t[i] - minimum(t)) / max(maximum(t) - minimum(t), EPS) * (nb - 1))) + 1, 1, nb)
        sdata[b] += ctx.z[i]; smodel[b] += pred[i]; sresid[b] += ctx.z[i] - pred[i]; counts[b] += 1
    end
    open(profile, "w") do io
        println(io, "t_nm\tdata_mean\tmodel_mean\tresidual_mean\tcount")
        for b in 1:nb
            tb = minimum(t) + (b - 0.5) / nb * (maximum(t) - minimum(t))
            c = max(counts[b], 1)
            println(io, "$(tb)\t$(sdata[b]/c)\t$(smodel[b]/c)\t$(sresid[b]/c)\t$(counts[b])")
        end
    end
    return summary, params, axisfile, profile, supportfile
end

function _plot_chain_sweep(results::Vector{ChainModelResult}, best::ChainModelResult, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig)
    xs, ys, zimg, mask = ctx.xs, ctx.ys, ctx.zimg, ctx.mask
    pred_img = zeros(size(zimg))
    if best.success
        for iy in eachindex(ys), ix in eachindex(xs)
            pred_img[iy, ix] = best.n == 0 ? best.params[1] :
                _chain_model_values([xs[ix]], [ys[iy]], best.params, best.n, ctx.axisctx, ccfg;
                                    amp_min=best.amp_min, amp_range=best.amp_range)[1]
        end
    end
    p1 = heatmap(xs, ys, zimg; aspect_ratio=:equal, title="Chain ROI data", xlabel="x (nm)", ylabel="y (nm)")
    contour!(p1, xs, ys, Float64.(mask); levels=[0.5], color=:white, linewidth=2)
    ox, oy = ctx.axisctx.origin; ax, ay = ctx.axisctx.axis
    plot!(p1, [ox + ctx.axisctx.tmin*ax, ox + ctx.axisctx.tmax*ax], [oy + ctx.axisctx.tmin*ay, oy + ctx.axisctx.tmax*ay]; color=:yellow, linewidth=2, label="axis")
    if best.n > 0 && best.success
        _b, feats, _ts, _us, _spars, _sperps = _decode_chain(best.params, best.n, ctx.axisctx, ccfg;
                                                              amp_min=best.amp_min, amp_range=best.amp_range)
        scatter!(p1, [f.x_nm for f in feats], [f.y_nm for f in feats]; color=:red, marker=:cross, markersize=10, linewidth=2.5, label="lobes")
    end
    p2 = heatmap(xs, ys, pred_img; aspect_ratio=:equal, title="Ordered chain model N=$(best.n)", xlabel="x (nm)", ylabel="y (nm)")
    p3 = heatmap(xs, ys, (zimg .- pred_img) .* mask ./ max(ctx.noise, EPS); aspect_ratio=:equal, title="Residuals / noise", xlabel="x (nm)", ylabel="y (nm)")
    Ns = [r.n for r in results]
    p4 = plot(Ns, [r.cv_nll_mean for r in results]; marker=:circle, label="CV NLL", xlabel="N", ylabel="score", title="Chain model selection")
    plot!(p4, Ns, [r.bic / max(1, length(ctx.z)) for r in results]; marker=:diamond, label="BIC / n")
    plot!(p4, Ns, [r.residual_peak_snr for r in results]; marker=:utriangle, label="resid SNR")
    vline!(p4, [best.n]; label="selected", color=:red)
    fig = plot(p1, p2, p3, p4; layout=(2,2), size=(1200, 950))
    out = joinpath(cfg.output_dir, "chain_model_selection.png")
    savefig(fig, out)
    return out
end

function _plot_chain_model_grid(r::ChainModelResult, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig, outdir::String)
    r.success || return nothing
    xs, ys, zimg, mask = ctx.xs, ctx.ys, ctx.zimg, ctx.mask
    pred_img = zeros(size(zimg))
    for iy in eachindex(ys), ix in eachindex(xs)
        pred_img[iy, ix] = r.n == 0 ? r.params[1] :
            _chain_model_values([xs[ix]], [ys[iy]], r.params, r.n, ctx.axisctx, ccfg;
                                amp_min=r.amp_min, amp_range=r.amp_range)[1]
    end
    resid = (zimg .- pred_img) .* mask ./ max(ctx.noise, EPS)
    p1 = heatmap(xs, ys, zimg; aspect_ratio=:equal, title="Data ROI", xlabel="x (nm)", ylabel="y (nm)")
    contour!(p1, xs, ys, Float64.(mask); levels=[0.5], color=:white, linewidth=2)
    ox, oy = ctx.axisctx.origin; ax, ay = ctx.axisctx.axis
    plot!(p1, [ox + ctx.axisctx.tmin*ax, ox + ctx.axisctx.tmax*ax], [oy + ctx.axisctx.tmin*ay, oy + ctx.axisctx.tmax*ay]; color=:yellow, linewidth=2, label=false)
    if r.n > 0
        _b, feats, _ts, _us, _spars, _sperps = _decode_chain(r.params, r.n, ctx.axisctx, ccfg;
                                                              amp_min=r.amp_min, amp_range=r.amp_range)
        scatter!(p1, [f.x_nm for f in feats], [f.y_nm for f in feats]; color=:red, marker=:cross, markersize=10, linewidth=2.5, label=false)
    end
    p2 = heatmap(xs, ys, pred_img; aspect_ratio=:equal, title="Chain model N=$(r.n)", xlabel="x (nm)", ylabel="y (nm)")
    p3 = heatmap(xs, ys, resid; aspect_ratio=:equal, title="Residual / noise", xlabel="x (nm)", ylabel="y (nm)")
    p4 = plot(title="scores", legend=false, framestyle=:box)
    annotate!(p4, 0.05, 0.82, text("valid=$(r.valid)", :left, 10))
    annotate!(p4, 0.05, 0.68, text(@sprintf("CV %.4g ± %.3g", r.cv_nll_mean, r.cv_nll_std), :left, 10))
    annotate!(p4, 0.05, 0.54, text(@sprintf("BIC %.4g", r.bic), :left, 10))
    annotate!(p4, 0.05, 0.40, text(@sprintf("resid SNR %.3g", r.residual_peak_snr), :left, 10))
    annotate!(p4, 0.05, 0.26, text(@sprintf("spacing %.3g CV %.3g", r.mean_spacing_nm, r.spacing_cv), :left, 10))
    annotate!(p4, 0.05, 0.12, text(@sprintf("σ∥ %.3g σ⊥ %.3g", r.sigma_parallel_nm, r.sigma_perp_nm), :left, 10))
    xlims!(p4, 0, 1); ylims!(p4, 0, 1)
    fig = plot(p1, p2, p3, p4; layout=(2,2), size=(1200, 950))
    out = joinpath(outdir, @sprintf("chain_model_N%02d.png", r.n))
    savefig(fig, out)
    return out
end

function _plot_chain_models_by_n(results::Vector{ChainModelResult}, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig)
    outdir = joinpath(cfg.output_dir, "models_by_N")
    mkpath(outdir)
    paths = String[]
    for r in results
        path = _plot_chain_model_grid(r, ctx, cfg, ccfg, outdir)
        path !== nothing && push!(paths, path)
    end
    return paths
end

function _run_chain_sweep(img::SXMImage, cfg::PatternConfig, ccfg::ChainSweepConfig)
    results, best, ctx = chain_gaussian_sweep(img, cfg, ccfg)
    summary, params, axisfile, profile, supportfile = _write_chain_sweep(results, best, ctx, cfg, ccfg)
    plot_path = cfg.no_plot ? nothing : _plot_chain_sweep(results, best, ctx, cfg, ccfg)
    model_plots = cfg.no_plot ? String[] : _plot_chain_models_by_n(results, ctx, cfg, ccfg)
    println("chain_selected_N: $(best.n)")
    println("chain_selected_valid: $(best.valid)")
    println("chain_selected_cv_nll: $(best.cv_nll_mean) ± $(best.cv_nll_std)")
    println("chain_selected_bic: $(best.bic)")
    println("chain_residual_peak_snr: $(best.residual_peak_snr)")
    println("chain_mean_spacing_nm: $(best.mean_spacing_nm)")
    println("chain_summary: $summary")
    println("chain_lobes: $params")
    println("chain_axis: $axisfile")
    println("chain_axis_profile: $profile")
    println("chain_support_metadata: $supportfile")
    plot_path !== nothing && println("chain_plot: $plot_path")
    !isempty(model_plots) && println("chain_models_by_N_dir: $(joinpath(cfg.output_dir, "models_by_N"))")
end

