#!/usr/bin/env julia
#
# Experimental: Joint fwd/bwd refit audit for chain model selection.
#
# For each file and each N:
#   1. Warm-start from the fused-fit result parameters.
#   2. Refine chain geometry and linear recalibration jointly via LsqFit,
#      minimizing (z_fwd - a_fwd·model - b_fwd)² + (z_bwd - a_bwd·model - b_bwd)².
#   3. If refit fails, fall back to per-scan OLS rescoring (no geometry change).
#   4. Compute joint metrics and compare source vs. joint selection.
#
# Usage:
#   julia --project=. test/audit_fwd_bwd_joint_refit.jl \
#     --config config/chitosan.toml \
#     --files 240817_019.sxm \
#     --out results/fwd_bwd_joint_refit_smoke/joint_refit_metrics.tsv
#
# Notes:
#   - Does NOT modify any configuration file or core selection logic.
#   - The joint GCV is nd/(nd - p_eff)² · (rss_fwd + rss_bwd)
#     with p_eff = _chain_nparams(N) + 4, nd = 2 · npixels_in_mask.
#   - `selected_by_source` = source GCV pick; `selected_by_joint_refit` = refit joint GCV pick.

using GaussianFit2D, Printf, TOML, Statistics, LinearAlgebra

# ---------------------------------------------------------------------------
# Constants & CLI
# ---------------------------------------------------------------------------

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = ["240817_017.sxm", "240817_019.sxm"]

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    out_tsv = "results/fwd_bwd_joint_refit_audit/joint_refit_metrics.tsv"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i < length(args) || error("--config requires a file path")
            config_file = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--files"
            i < length(args) || error("--files requires comma-separated names")
            files = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--files=")
            files = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--out"
            i < length(args) || error("--out requires a TSV path")
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, out_tsv
end

# ---------------------------------------------------------------------------
# Config builders (mirrors audit scripts, intelligent_sweep=false)
# ---------------------------------------------------------------------------

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
        intelligent_sweep=false,
        fuse_z_bwd=true)
    return pcfg, ccfg
end

# ---------------------------------------------------------------------------
# Chain parameter bounds (mirrors _fit_chain_n in core.jl)
# ---------------------------------------------------------------------------

function _chain_param_bounds(n::Int, ccfg::GaussianFit2D.ChainSweepConfig)
    np = GaussianFit2D._chain_nparams(n, ccfg)
    lower = fill(-10.0, np); upper = fill(10.0, np)
    lower[1] = -5.0; upper[1] = 5.0   # b0
    j = 2
    if ccfg.chain_tilted_baseline
        lower[j] = -1.0; upper[j] = 1.0; j += 1  # bx
        lower[j] = -1.0; upper[j] = 1.0; j += 1  # by
    end
    for _ in 1:n; lower[j] = -5.0; upper[j] = 5.0; j += 1 end   # amps
    n_spacing = GaussianFit2D._chain_spacing_param_count(n, ccfg)
    lower[j] = -4.0; upper[j] = 4.0; j += 1                     # t0
    for _ in 1:(n_spacing - 1); lower[j] = -5.0; upper[j] = 5.0; j += 1 end  # deltas
    for _ in 1:n; lower[j] = -3.0; upper[j] = 3.0; j += 1 end   # us (laterals)
    n_sigma_types = GaussianFit2D._chain_sigma_param_count(n, ccfg)
    if ccfg.chain_circular_sigmas
        for _ in 1:n_sigma_types; lower[j] = -5.0; upper[j] = 5.0; j += 1 end
    else
        for _ in 1:(2n_sigma_types); lower[j] = -5.0; upper[j] = 5.0; j += 1 end
    end
    return lower, upper
end

# ---------------------------------------------------------------------------
# Per-scan OLS recalibration (fallback / rescore)
# ---------------------------------------------------------------------------

function _linear_recalibrate(z_scan::Vector{Float64}, model::Vector{Float64})
    n = length(z_scan)
    n < 4 && return NaN, NaN, NaN, Float64[]
    A = hcat(model, ones(n))
    try
        coeff = A \ z_scan
        a, b = coeff[1], coeff[2]
        pred = a .* model .+ b
        rss = sum(abs2, z_scan .- pred)
        return a, b, rss, pred
    catch
        return NaN, NaN, NaN, Float64[]
    end
end

# ---------------------------------------------------------------------------
# Joint refit via LsqFit (uses subsampled data for speed)
# ---------------------------------------------------------------------------

const MAX_REFIT_PIXELS = 600  # subsample to this many pixels for Jacobian

function _subsample_indices(n::Int, max_n::Int)
    n <= max_n && return 1:n
    step = n ÷ max_n
    return 1:step:n
end

function _joint_refit(xdata, zdata, p0, n_val, axisctx, ccfg, amp_min, amp_range,
                      npix, lower, upper)
    # Build closure for the joint model.
    # xdata is 2×(2*npix): first npix cols = fwd coords, last npix = bwd coords.
    function model_f(xdata, p)
        chain_p = @view p[1:end-4]
        a_fwd = p[end-3]; b_fwd = p[end-2]
        a_bwd = p[end-1]; b_bwd = p[end]
        x = @view xdata[1, :]
        y = @view xdata[2, :]
        # Compute chain model once for all points (fwd and bwd share coords)
        chain_val = GaussianFit2D._chain_model_values(x, y, chain_p, n_val,
            axisctx, ccfg; amp_min=amp_min, amp_range=amp_range)
        result = similar(chain_val)
        @. result[1:npix] = a_fwd * chain_val[1:npix] + b_fwd
        @. result[npix+1:end] = a_bwd * chain_val[npix+1:end] + b_bwd
        return result
    end

    try
        fit = GaussianFit2D.LsqFit.curve_fit(model_f, xdata, zdata, p0;
                                              lower=lower, upper=upper,
                                              maxIter=80, autodiff=:finite)
        p_final = fit.param
        pred = model_f(xdata, p_final)
        rss = sum(abs2, zdata .- pred)
        return p_final, rss, pred, true
    catch
        return p0, NaN, Float64[], false
    end
end

# ---------------------------------------------------------------------------
# Metrics computation (shared between refit and fallback)
# ---------------------------------------------------------------------------

function _compute_joint_metrics(z_fwd, z_bwd, pred_fwd, pred_bwd, n, ccfg)
    npix = length(z_fwd)
    nd_joint = 2 * npix

    rss_fwd = sum(abs2, z_fwd .- pred_fwd)
    rss_bwd = sum(abs2, z_bwd .- pred_bwd)

    resid_fwd = z_fwd .- pred_fwd
    resid_bwd = z_bwd .- pred_bwd

    # Residual correlation
    rf = resid_fwd .- mean(resid_fwd)
    rb = resid_bwd .- mean(resid_bwd)
    denom = norm(rf) * norm(rb)
    resid_corr = denom > eps(Float64) ? dot(rf, rb) / denom : NaN

    # Joint GCV
    p_chain = GaussianFit2D._chain_nparams(n, ccfg)
    p_eff = p_chain + 4
    joint_gcv = nd_joint > p_eff ?
        nd_joint / (nd_joint - p_eff)^2 * (rss_fwd + rss_bwd) : Inf

    # Joint NRMSE
    var_total = var(z_fwd) + var(z_bwd)
    joint_nrmse = (var_total > eps(Float64) && isfinite(rss_fwd + rss_bwd)) ?
        sqrt((rss_fwd + rss_bwd) / (2 * npix)) / sqrt(var_total) : NaN

    return rss_fwd, rss_bwd, joint_gcv, joint_nrmse, resid_corr
end

# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

function _fmt(x)
    x === nothing && return ""
    x isa AbstractString && return x
    x isa Bool && return string(x)
    x isa Integer && return string(x)
    x isa Real && return isfinite(x) ? @sprintf("%.8g", x) : ""
    return string(x)
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    config_file, files, out_tsv = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    model = cfg["model"]
    preproc = cfg["preprocessing"]
    mkpath(dirname(out_tsv))

    header = ["file", "N", "source_gcv", "source_valid",
              "refit_success",
              "joint_refit_gcv", "joint_refit_nrmse",
              "joint_refit_rss_fwd", "joint_refit_rss_bwd",
              "a_fwd", "b_fwd", "a_bwd", "b_bwd",
              "residual_corr_fwd_bwd",
              "selected_by_source", "selected_by_joint_refit"]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end

    for fn in files
        @printf("Joint fwd/bwd refit audit: %s\n", fn)
        file_out = joinpath(dirname(out_tsv), splitext(fn)[1])
        mkpath(file_out)
        pcfg, ccfg = _configs(model, preproc, file_out)
        pcfg.filepath = joinpath(DATA_DIR, fn)

        # ── Read SXM ──
        img = GaussianFit2D.read_sxm(pcfg.filepath)

        # ── Fused sweep (intelligent_sweep=false → exhaustive) ──
        results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
        xs, ys = ctx.xs, ctx.ys
        mask = ctx.mask
        axisctx = ctx.axisctx

        # Validate mask
        if mask === nothing || isempty(mask)
            mask = trues(length(ys), length(xs))
        end
        if size(mask) != (length(ys), length(xs))
            @warn("$fn: mask size mismatch; fallback full grid")
            mask = trues(length(ys), length(xs))
        end
        mask_inds = findall(mask)
        npix = length(mask_inds)

        # ── Separate fwd / bwd preprocessing ──
        ch_fwd = GaussianFit2D.get_channel(img, "Z"; direction="fwd")
        _xsf, _ysf, _raw_fwd, z_fwd, _zs_fwd, _ufwd, _nfwd =
            GaussianFit2D.preprocess_channel(img, ch_fwd, pcfg)
        ch_bwd = GaussianFit2D.get_channel(img, "Z"; direction="bwd")
        _xsb, _ysb, _raw_bwd, z_bwd, _zs_bwd, _ubwd, _nbwd =
            GaussianFit2D.preprocess_channel(img, ch_bwd, pcfg)

        # Extract masked vectors: build coordinate grids then index by mask
        x_grid = repeat(reshape(xs, 1, :), length(ys), 1)
        y_grid = repeat(reshape(ys, :, 1), 1, length(xs))
        z_fwd_masked = vec(z_fwd[mask_inds])
        z_bwd_masked = vec(z_bwd[mask_inds])
        x_masked = vec(x_grid[mask_inds])
        y_masked = vec(y_grid[mask_inds])

        # Remove non-finite points from all arrays
        good = isfinite.(z_fwd_masked) .& isfinite.(z_bwd_masked)
        z_fwd_masked = z_fwd_masked[good]
        z_bwd_masked = z_bwd_masked[good]
        x_masked = x_masked[good]
        y_masked = y_masked[good]
        npix = length(x_masked)
        npix >= 20 || (@printf("  %s: too few valid mask pixels (%d), skip\n", fn, npix); continue)

        # ── Per-N loop ──
        rows = []
        joint_refit_scores = Dict{Int,Float64}()

        source_selected_n = (best !== nothing && best.success) ? best.n : 0

        # Precompute chain param bounds per N (cache)
        chain_bounds_cache = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

        for r in results
            (r.success && isfinite(r.gcv)) || continue
            n = r.n
            haskey(chain_bounds_cache, n) || (chain_bounds_cache[n] = _chain_param_bounds(n, ccfg))
            lb_chain, ub_chain = chain_bounds_cache[n]

            # ── Compute model on full grid (once, for both refit and fallback) ──
            xmat = repeat(reshape(xs, 1, :), length(ys), 1)
            ymat = repeat(reshape(ys, :, 1), 1, length(xs))
            model_full = GaussianFit2D._chain_model_values(
                vec(xmat), vec(ymat), r.params, n, axisctx, ccfg;
                amp_min=r.amp_min, amp_range=r.amp_range)
            model_img = reshape(model_full, length(ys), length(xs))
            model_masked_all = vec(model_img[mask_inds])[good]

            # Subsample for refit (speed), metrics on all data
            refit_idx = _subsample_indices(npix, MAX_REFIT_PIXELS)

            x_sub = x_masked[refit_idx]
            y_sub = y_masked[refit_idx]
            z_fwd_sub = z_fwd_masked[refit_idx]
            z_bwd_sub = z_bwd_masked[refit_idx]
            nsub = length(refit_idx)

            xdata_sub = hcat(vcat(x_sub', y_sub'), vcat(x_sub', y_sub'))
            zdata_sub = vcat(z_fwd_sub, z_bwd_sub)

            # Initialise joint params
            p0 = vcat(r.params, [1.0, 0.0, 1.0, 0.0])  # a_fwd,b_fwd, a_bwd,b_bwd
            lb = vcat(lb_chain, [-5.0, -2.0, -5.0, -2.0])
            ub = vcat(ub_chain, [ 5.0,  2.0,  5.0,  2.0])

            # ── Attempt joint refit (subsampled) ──
            p_refit, rss_joint_sub, pred_joint_sub, _ =
                _joint_refit(xdata_sub, zdata_sub, p0, n, axisctx, ccfg,
                             r.amp_min, r.amp_range, nsub, lb, ub)

            refit_success = isfinite(rss_joint_sub)
            rss_fwd_r = NaN; rss_bwd_r = NaN
            joint_gcv = Inf; joint_nrmse = NaN; resid_corr = NaN
            a_fwd_r = NaN; b_fwd_r = NaN; a_bwd_r = NaN; b_bwd_r = NaN

            if refit_success
                # Compute metrics on ALL data using refit params
                p_refit_chain = @view p_refit[1:end-4]
                a_fwd_r = p_refit[end-3]; b_fwd_r = p_refit[end-2]
                a_bwd_r = p_refit[end-1]; b_bwd_r = p_refit[end]

                # Full-resolution model with refit params
                model_refit_full = GaussianFit2D._chain_model_values(
                    vec(xmat), vec(ymat), p_refit_chain, n, axisctx, ccfg;
                    amp_min=r.amp_min, amp_range=r.amp_range)
                model_refit_img = reshape(model_refit_full, length(ys), length(xs))
                model_refit_masked = vec(model_refit_img[mask_inds])[good]

                pred_fwd = a_fwd_r .* model_refit_masked .+ b_fwd_r
                pred_bwd = a_bwd_r .* model_refit_masked .+ b_bwd_r

                rss_fwd_r, rss_bwd_r, joint_gcv, joint_nrmse, resid_corr =
                    _compute_joint_metrics(z_fwd_masked, z_bwd_masked,
                                           pred_fwd, pred_bwd, n, ccfg)
            else
                # ── Fallback: OLS rescore (no geometry change) ──
                a_fwd_r, b_fwd_r, rss_fwd_r, pred_fwd = _linear_recalibrate(z_fwd_masked, model_masked_all)
                a_bwd_r, b_bwd_r, rss_bwd_r, pred_bwd = _linear_recalibrate(z_bwd_masked, model_masked_all)

                if isfinite(a_fwd_r) && isfinite(a_bwd_r)
                    rss_fwd_r, rss_bwd_r, joint_gcv, joint_nrmse, resid_corr =
                        _compute_joint_metrics(z_fwd_masked, z_bwd_masked,
                                               pred_fwd, pred_bwd, n, ccfg)
                end
            end

            joint_refit_scores[n] = joint_gcv
            selected_by_source = (best !== nothing && best.success && best.n == n) ? 1 : 0

            push!(rows, (n=n, source_gcv=r.gcv, source_valid=r.valid,
                         refit_success=refit_success,
                         joint_gcv=joint_gcv, joint_nrmse=joint_nrmse,
                         rss_fwd=rss_fwd_r, rss_bwd=rss_bwd_r,
                         a_fwd=a_fwd_r, b_fwd=b_fwd_r,
                         a_bwd=a_bwd_r, b_bwd=b_bwd_r,
                         resid_corr=resid_corr,
                         selected_by_source=selected_by_source))
        end

        # ── Joint refit selection ──
        best_joint_n = -1
        if !isempty(joint_refit_scores)
            best_joint_n = argmin(joint_refit_scores)
        end

        # ── Write rows ──
        for row in rows
            sel_joint = (row.n == best_joint_n) ? 1 : 0
            out_row = Any[fn, row.n,
                          row.source_gcv, row.source_valid ? 1 : 0,
                          row.refit_success ? 1 : 0,
                          row.joint_gcv, row.joint_nrmse,
                          row.rss_fwd, row.rss_bwd,
                          row.a_fwd, row.b_fwd, row.a_bwd, row.b_bwd,
                          row.resid_corr,
                          row.selected_by_source, sel_joint]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(out_row), '\t'))
            end
        end

        if isempty(rows)
            open(out_tsv, "a") do io
                println(io, "$fn\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
            end
        end
    end

    println("Wrote $out_tsv")
end

main()
