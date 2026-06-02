#!/usr/bin/env julia

using GaussianFit2D, Printf, TOML, Statistics, LinearAlgebra, Plots

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = [
    "240817_017.sxm", "240817_019.sxm", "240817_043.sxm", "240817_058.sxm",
    "240817_018.sxm", "240817_021.sxm", "240817_060.sxm",
    "240817_029.sxm", "240817_030.sxm", "240817_031.sxm", "240817_032.sxm",
    "240817_034.sxm", "240817_035.sxm", "240817_037.sxm", "240817_038.sxm", "240817_051.sxm",
]

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    outdir = "results/stm_artifact_audit"
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
        intelligent_sweep=true, fuse_z_bwd=true)
    return pcfg, ccfg
end

_mad(v) = 1.4826 * median(abs.(v .- median(v)))

function _ncc_shift(a::AbstractMatrix, b::AbstractMatrix, dy::Int, dx::Int)
    ny, nx = size(a)
    ya = max(1, 1 + dy):min(ny, ny + dy)
    yb = max(1, 1 - dy):min(ny, ny - dy)
    xa = max(1, 1 + dx):min(nx, nx + dx)
    xb = max(1, 1 - dx):min(nx, nx - dx)
    (length(ya) < 8 || length(xa) < 8) && return NaN
    av = vec(@view a[ya, xa]); bv = vec(@view b[yb, xb])
    am = av .- mean(av); bm = bv .- mean(bv)
    den = norm(am) * norm(bm)
    den <= eps(Float64) && return NaN
    return dot(am, bm) / den
end

function _best_shift_corr(a, b, px_nm, py_nm; min_shift_nm=0.45, max_shift_nm=1.80)
    stride = max(1, Int(ceil(maximum(size(a)) / 96)))
    if stride > 1
        a = a[1:stride:end, 1:stride:end]
        b = b[1:stride:end, 1:stride:end]
        px_nm *= stride
        py_nm *= stride
    end
    max_dx = max(1, Int(round(max_shift_nm / px_nm)))
    max_dy = max(1, Int(round(max_shift_nm / py_nm)))
    min_r2 = min_shift_nm^2
    best = (corr=-Inf, dx=0, dy=0, dist_nm=NaN)
    for dy in -max_dy:max_dy, dx in -max_dx:max_dx
        dx == 0 && dy == 0 && continue
        dist2 = (dx * px_nm)^2 + (dy * py_nm)^2
        dist2 >= min_r2 || continue
        c = _ncc_shift(a, b, dy, dx)
        if isfinite(c) && c > best.corr
            best = (corr=c, dx=dx, dy=dy, dist_nm=sqrt(dist2))
        end
    end
    return best
end

function _model_image(xs, ys, best, axisctx, ccfg)
    xmat = repeat(reshape(xs, 1, :), length(ys), 1)
    ymat = repeat(reshape(ys, :, 1), 1, length(xs))
    pred = GaussianFit2D._chain_model_values(vec(xmat), vec(ymat), best.params, best.n, axisctx, ccfg;
        amp_min=best.amp_min, amp_range=best.amp_range)
    return reshape(pred, length(ys), length(xs))
end

function _plot_diagnostics(path, z, model, resid, title)
    p1 = heatmap(z, aspect_ratio=:equal, colorbar=false, title="image")
    p2 = heatmap(model, aspect_ratio=:equal, colorbar=false, title="model")
    p3 = heatmap(resid, aspect_ratio=:equal, colorbar=false, title="residual")
    savefig(plot(p1, p2, p3, layout=(1,3), size=(1200, 360), plot_title=title), path)
end

function _fmt(x)
    x isa AbstractString && return x
    x isa Integer && return string(x)
    x isa Bool && return string(x)
    x === nothing && return ""
    try
        isfinite(x) ? @sprintf("%.8g", x) : ""
    catch
        string(x)
    end
end

function main()
    config_file, files, outdir = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    mkpath(outdir)
    out_tsv = joinpath(outdir, "artifact_metrics.tsv")
    header = ["file", "N", "score_gcv", "fwd_bwd_corr", "fwd_bwd_nrmse", "double_tip_corr", "double_tip_dx_nm", "double_tip_dy_nm", "double_tip_dist_nm",
        "ghost_corr", "ghost_dx_nm", "ghost_dy_nm", "ghost_dist_nm", "line_discontinuity", "stripe_fft", "residual_peak_snr", "overlap", "kappa"]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end
    model = cfg["model"]; preproc = cfg["preprocessing"]
    for fn in files
        @printf("Artifact audit %s\n", fn)
        file_out = joinpath(outdir, splitext(fn)[1]); mkpath(file_out)
        pcfg, ccfg = _configs(model, preproc, file_out)
        pcfg.filepath = joinpath(DATA_DIR, fn)
        img = GaussianFit2D.read_sxm(pcfg.filepath)
        art = GaussianFit2D.compute_image_artifact_diagnostics(img, pcfg)
        results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
        xs, ys, zimg = ctx.xs, ctx.ys, ctx.zimg
        px_nm = median(diff(xs)); py_nm = median(diff(ys))
        z0 = zimg .- median(vec(zimg))
        dt = _best_shift_corr(z0, z0, px_nm, py_nm)
        line_score = art.line_discontinuity
        stripe_score = art.stripe_periodicity
        if best !== nothing && best.success
            model_img = _model_image(xs, ys, best, ctx.axisctx, ccfg)
            resid = zimg .- model_img
            ghost = _best_shift_corr(resid .- median(vec(resid)), model_img .- median(vec(model_img)), px_nm, py_nm)
            _plot_diagnostics(joinpath(file_out, "artifact_diagnostics.png"), zimg, model_img, resid,
                @sprintf("%s N=%d", fn, best.n))
            row = Any[fn, best.n, best.gcv, art.fwd_bwd_corr, art.fwd_bwd_nrmse, dt.corr, dt.dx * px_nm, dt.dy * py_nm, dt.dist_nm,
                ghost.corr, ghost.dx * px_nm, ghost.dy * py_nm, ghost.dist_nm, line_score, stripe_score,
                best.residual_peak_snr, best.overlap, best.kappa_max_adj]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(row), '\t'))
            end
        else
            row = Any[fn, "", "", art.fwd_bwd_corr, art.fwd_bwd_nrmse, dt.corr, dt.dx * px_nm, dt.dy * py_nm, dt.dist_nm,
                "", "", "", "", line_score, stripe_score, "", "", ""]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(row), '\t'))
            end
        end
    end
    println("Wrote $out_tsv")
end

main()
