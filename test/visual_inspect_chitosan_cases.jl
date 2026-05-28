#!/usr/bin/env julia

using GaussianFit2D, Printf, TOML, Plots

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = ["240817_017.sxm", "240817_019.sxm", "240817_043.sxm", "240817_058.sxm"]

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    outdir = "results/chitosan_visual_inspection"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            config_file = args[i + 1]; i += 2
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--files"
            files = split(args[i + 1], ","); i += 2
        elseif startswith(arg, "--files=")
            files = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--outdir"
            outdir = args[i + 1]; i += 2
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]; i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, outdir
end

function _configs(model, preproc, output_dir)
    pcfg = GaussianFit2D.PatternConfig(filepath="", channel="Z", direction="fwd",
        stride=get(preproc, "stride", 1),
        flatten=get(preproc, "flatten", "plane+rows"),
        smooth_radius_px=get(preproc, "smooth_radius_px", 1),
        output_dir=output_dir, no_plot=true)
    ccfg = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
        spacing_min_nm=model["spacing_min_nm"], spacing_max_nm=model["spacing_max_nm"],
        fit_width_nm=model["fit_width_nm"],
        support_noise_k=model["support_noise_k"],
        support_padding_nm=model["support_padding_nm"],
        support_min_length_nm=get(model, "support_min_length_nm", 1.0),
        support_baseline_quantile=get(model, "support_baseline_quantile", 0.10),
        max_overlap=model["max_overlap"],
        global_maxtime=model["global_maxtime"], global_maxiter=model["global_maxiter"],
        max_iter=get(model, "max_iter", 300), multistart=get(model, "multistart", 1),
        cv_folds=get(model, "cv_folds", 5), cv_method=get(model, "cv_method", "gcv"),
        selection_criterion=get(model, "selection_criterion", "gcv"),
        sigma_parallel_min_nm=model["sigma_parallel_min_nm"],
        sigma_parallel_max_nm=model["sigma_parallel_max_nm"],
        sigma_perp_min_nm=model["sigma_parallel_min_nm"],
        sigma_perp_max_nm=model["sigma_parallel_max_nm"],
        kappa_max=get(model, "kappa_max", 10.0),
        kappa_weight=get(model, "kappa_weight", 1.0),
        min_amplitude_fraction=get(model, "min_amplitude_fraction", 0.3),
        shared_sigma_types=get(model, "shared_sigma_types", 0),
        chain_spacing_model=get(model, "chain_spacing_model", "free"),
        chain_tilted_baseline=get(model, "chain_tilted_baseline", true),
        intelligent_sweep=false, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg); ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

_score(r, criterion) = lowercase(String(criterion)) == "gcv" ? r.gcv : lowercase(String(criterion)) == "aicc" ? r.aicc : lowercase(String(criterion)) == "cv" ? r.cv_nll_mean : r.bic
_valid(r, criterion) = r.success && r.valid && isfinite(_score(r, criterion))

function _best_by_n(results, criterion)
    by_n = Dict{Int,Any}()
    for r in results
        _valid(r, criterion) || continue
        if !haskey(by_n, r.n) || _score(r, criterion) < _score(by_n[r.n], criterion)
            by_n[r.n] = r
        end
    end
    return by_n
end

function _effective_best(by_ell, by_circ, criterion)
    best_n, best_r, best_source, best_s = 0, nothing, "NA", Inf
    for n in sort(unique(vcat(collect(keys(by_ell)), collect(keys(by_circ)))))
        r_ell = get(by_ell, n, nothing); r_circ = get(by_circ, n, nothing)
        s_ell = r_ell === nothing ? Inf : _score(r_ell, criterion)
        s_circ = r_circ === nothing ? Inf : _score(r_circ, criterion)
        if min(s_ell, s_circ) < best_s
            best_n, best_s = n, min(s_ell, s_circ)
            best_r, best_source = s_ell <= s_circ ? (r_ell, "ell") : (r_circ, "circ")
        end
    end
    return best_n, best_r, best_source
end

function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    refined = GaussianFit2D.ChainModelResult[]
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ctx_circ.axisctx_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    ccfg_refine = deepcopy(ccfg_ell); ccfg_refine.skip_global = true; ccfg_refine.max_iter = 50; ccfg_refine.multistart = 1
    for r_c in results_circ
        r_c.success || continue
        n = r_c.n
        n_prefix = 1 + (ccfg_refine.chain_tilted_baseline ? 2 : 0)
        split_idx = n_prefix + n + GaussianFit2D._chain_spacing_param_count(n, ccfg_refine) + n
        p_init = vcat(r_c.params[1:split_idx], r_c.params[(split_idx+1):end], r_c.params[(split_idx+1):end])
        try
            r = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r.params, n, ac_fit, ccfg_refine; amp_min=r.amp_min, amp_range=r.amp_range)
                GaussianFit2D._finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, ac_fit, ccfg_refine)
                push!(refined, r)
            end
        catch err
            @warn "elliptical refinement failed" n exception=(err, catch_backtrace())
        end
    end
    return refined
end

function _prediction_image(ctx, r, ccfg)
    pred = fill(NaN, size(ctx.zimg))
    for iy in eachindex(ctx.ys), ix in eachindex(ctx.xs)
        ctx.mask[iy, ix] || continue
        pred[iy, ix] = GaussianFit2D._chain_model_values([ctx.xs[ix]], [ctx.ys[iy]], r.params, r.n, ctx.axisctx, ccfg; amp_min=r.amp_min, amp_range=r.amp_range)[1]
    end
    return pred
end

function _overlay_features!(p, ctx, r, ccfg; color=:cyan)
    _b, feats, _ts, _us, _spars, _sperps = GaussianFit2D._decode_chain(r.params, r.n, ctx.axisctx, ccfg; amp_min=r.amp_min, amp_range=r.amp_range)
    scatter!(p, [f.x_nm for f in feats], [f.y_nm for f in feats]; color=color, markersize=4, label="")
    for (i, f) in enumerate(feats)
        annotate!(p, f.x_nm, f.y_nm, text(string(i), 8, color))
    end
end

function _panel(ctx, r, ccfg, title)
    pred = _prediction_image(ctx, r, ccfg)
    resid = (ctx.zimg .- pred) ./ max(ctx.noise, eps(Float64))
    p1 = heatmap(ctx.xs, ctx.ys, ctx.zimg; aspect_ratio=:equal, colorbar=false, title=title * " data")
    _overlay_features!(p1, ctx, r, ccfg)
    p2 = heatmap(ctx.xs, ctx.ys, pred; aspect_ratio=:equal, colorbar=false, title="model")
    _overlay_features!(p2, ctx, r, ccfg)
    p3 = heatmap(ctx.xs, ctx.ys, resid; aspect_ratio=:equal, colorbar=false, clims=(-3, 3), color=:balance, title="residual/noise")
    _overlay_features!(p3, ctx, r, ccfg; color=:black)
    return p1, p2, p3
end

function main()
    config_file, files, outdir = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    model, preproc = cfg["model"], cfg["preprocessing"]
    criterion = get(model, "selection_criterion", "gcv")
    mkpath(outdir)
    summary_tsv = joinpath(outdir, "visual_inspection_summary.tsv")
    features_tsv = joinpath(outdir, "visual_inspection_features.tsv")
    open(summary_tsv, "w") do io
        println(io, join(["file", "model", "N", "source", "gcv", "rss", "chi2", "residual_peak_snr", "mean_spacing_nm", "spacing_cv", "overlap", "kappa", "support_nm", "plot"], '\t'))
    end
    open(features_tsv, "w") do io
        println(io, join(["file", "model", "N", "source", "lobe", "t_nm", "u_nm", "amplitude", "amp_rel", "sigma_parallel_nm", "sigma_perp_nm", "spacing_prev_nm"], '\t'))
    end
    for fn in files
        @printf("Visual inspection %s\n", fn); flush(stdout)
        pcfg, ccfg, ccfg_circ = _configs(model, preproc, joinpath(outdir, splitext(fn)[1]))
        pcfg.filepath = joinpath(DATA_DIR, fn)
        img = GaussianFit2D.read_sxm(pcfg.filepath)
        results_circ, _, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
        results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx)
        by_ell, by_circ = _best_by_n(results_ell, criterion), _best_by_n(results_circ, criterion)
        selected_n, selected_r, selected_source = _effective_best(by_ell, by_circ, criterion)
        r6 = get(by_ell, 6, get(by_circ, 6, nothing))
        r6 === nothing && (@warn "No N=6 result" file=fn; continue)
        ccfg_sel = selected_source == "circ" ? ccfg_circ : ccfg
        p = plot(_panel(ctx, r6, ccfg, @sprintf("%s N=6 GCV=%.3g", fn, r6.gcv))...,
                 _panel(ctx, selected_r, ccfg_sel, @sprintf("selected N=%d %s GCV=%.3g", selected_n, selected_source, selected_r.gcv))...;
                 layout=(2, 3), size=(1500, 900))
        out = joinpath(outdir, splitext(fn)[1] * "_N6_vs_selected.png")
        savefig(p, out)
        support_nm = ctx.axisctx.tmax - ctx.axisctx.tmin
        for (label, r, source, cfg_use) in [("N6", r6, "ell", ccfg), ("selected", selected_r, selected_source, ccfg_sel)]
            open(summary_tsv, "a") do io
                println(io, join([fn, label, r.n, source, r.gcv, r.rss, r.chi2_reduced,
                                  r.residual_peak_snr, r.mean_spacing_nm, r.spacing_cv,
                                  r.overlap, r.kappa_max_adj, support_nm, out], '\t'))
            end
            _b, feats, ts, us, spars, sperps = GaussianFit2D._decode_chain(r.params, r.n, ctx.axisctx, cfg_use; amp_min=r.amp_min, amp_range=r.amp_range)
            amax = maximum([f.amplitude for f in feats])
            open(features_tsv, "a") do io
                for i in 1:length(feats)
                    println(io, join([fn, label, r.n, source, i, ts[i], us[i], feats[i].amplitude,
                                      feats[i].amplitude / max(amax, eps(Float64)), spars[i], sperps[i],
                                      i == 1 ? "" : ts[i] - ts[i-1]], '\t'))
                end
            end
        end
        println("  wrote $out")
    end
    println("Wrote $summary_tsv")
    println("Wrote $features_tsv")
end

main()
