#!/usr/bin/env julia
#
# Experimental: Joint fwd/bwd consistency audit for chain model selection.
#
# For each file:
#   1. Run fused 2D chain sweep (intelligent_sweep=false) to get per-N results.
#   2. For each successful N, compute the model image on the full grid.
#   3. Preprocess fwd and bwd scans separately.
#   4. Within the fused ROI mask, fit a linear recalibration
#      z_scan ≈ a * model + b  (by OLS, per scan).
#   5. Compute joint metrics: joint GCV, joint NRMSE, per-scan RSS,
#      recalibration coefficients, and residual correlation.
#   6. Mark the source-selected N and the joint-selected N.
#
# Usage:
#   julia --project=. test/audit_fwd_bwd_joint.jl \
#     --config config/chitosan.toml \
#     --files 240817_017.sxm,240817_019.sxm \
#     --out results/fwd_bwd_joint_audit/joint_metrics.tsv
#
# Notes:
#   - Does NOT modify any configuration file or core selection logic.
#   - The joint GCV is: nd / (nd - pcount_eff)^2 * (rss_fwd + rss_bwd)
#     where pcount_eff = _chain_nparams(N) + 4, nd = 2 * npixels_in_mask.
#   - `selected_by_source` = 1 for the N that the fused sweep picks (lowest gcv).
#   - `selected_by_joint`  = 1 for the N with the lowest joint_gcv.
#   - Metrics use the fused ROI mask region (documented below).

using GaussianFit2D, Printf, TOML, Statistics, LinearAlgebra

# ---------------------------------------------------------------------------
# Constants & CLI
# ---------------------------------------------------------------------------

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = ["240817_017.sxm", "240817_019.sxm"]

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    out_tsv = "results/fwd_bwd_joint_audit/joint_metrics.tsv"
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
# Config builders (mirrors audit_chitosan_cases.jl but with intelligent_sweep=false)
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
        intelligent_sweep=false,   # exhaustive linear sweep → per-N results
        fuse_z_bwd=true)
    return pcfg, ccfg
end

# ---------------------------------------------------------------------------
# Linear recalibration:  z_scan ≈ a * model + b  (OLS)
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

    header = ["file", "N", "source_score_gcv", "source_valid",
              "joint_gcv", "joint_nrmse", "rss_fwd", "rss_bwd",
              "fit_a_fwd", "fit_b_fwd", "fit_a_bwd", "fit_b_bwd",
              "residual_corr_fwd_bwd", "selected_by_source", "selected_by_joint"]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end

    for fn in files
        @printf("Joint fwd/bwd audit: %s\n", fn)
        file_out = joinpath(dirname(out_tsv), splitext(fn)[1])
        mkpath(file_out)
        pcfg, ccfg = _configs(model, preproc, file_out)
        pcfg.filepath = joinpath(DATA_DIR, fn)

        # ── Read SXM ──
        img = GaussianFit2D.read_sxm(pcfg.filepath)

        # ── Fused sweep (intelligent_sweep=false → linear exhaustive sweep) ──
        results, best, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)
        xs, ys = ctx.xs, ctx.ys
        mask = ctx.mask
        axisctx = ctx.axisctx

        # Validate mask
        if mask === nothing || isempty(mask)
            mask = trues(length(ys), length(xs))
        end
        if size(mask) != (length(ys), length(xs))
            @warn("$fn: mask size $(size(mask)) != grid ($(length(ys)), $(length(xs))); fallback full grid")
            mask = trues(length(ys), length(xs))
        end
        mask_inds = findall(mask)

        # ── Separate fwd / bwd preprocessing ──
        ch_fwd = GaussianFit2D.get_channel(img, "Z"; direction="fwd")
        _xsf, _ysf, _raw_fwd, z_fwd, _zs_fwd, _ufwd, _nfwd =
            GaussianFit2D.preprocess_channel(img, ch_fwd, pcfg)
        ch_bwd = GaussianFit2D.get_channel(img, "Z"; direction="bwd")
        _xsb, _ysb, _raw_bwd, z_bwd, _zs_bwd, _ubwd, _nbwd =
            GaussianFit2D.preprocess_channel(img, ch_bwd, pcfg)

        # ── Per-N evaluation ──
        rows = []
        joint_scores = Dict{Int,Float64}()  # n → joint_gcv

        source_selected_n = (best !== nothing && best.success) ? best.n : 0

        for r in results
            (r.success && isfinite(r.gcv)) || continue

            # Compute model on full grid
            xmat = repeat(reshape(xs, 1, :), length(ys), 1)
            ymat = repeat(reshape(ys, :, 1), 1, length(xs))
            model_full = GaussianFit2D._chain_model_values(
                vec(xmat), vec(ymat), r.params, r.n, axisctx, ccfg;
                amp_min=r.amp_min, amp_range=r.amp_range)
            model_img = reshape(model_full, length(ys), length(xs))

            # Extract masked pixels
            z_fwd_masked = vec(z_fwd[mask_inds])
            z_bwd_masked = vec(z_bwd[mask_inds])
            model_masked = vec(model_img[mask_inds])

            # Remove non-finite points
            good = isfinite.(z_fwd_masked) .& isfinite.(z_bwd_masked) .& isfinite.(model_masked)
            z_fwd_masked = z_fwd_masked[good]
            z_bwd_masked = z_bwd_masked[good]
            model_masked = model_masked[good]

            length(model_masked) >= 20 || continue

            # ── Linear recalibration per scan ──
            a_fwd, b_fwd, rss_fwd, pred_fwd = _linear_recalibrate(z_fwd_masked, model_masked)
            a_bwd, b_bwd, rss_bwd, pred_bwd = _linear_recalibrate(z_bwd_masked, model_masked)
            (isfinite(a_fwd) && isfinite(a_bwd)) || continue

            # ── Residuals ──
            resid_fwd = z_fwd_masked .- pred_fwd
            resid_bwd = z_bwd_masked .- pred_bwd

            # ── Residual correlation (Pearson) ──
            rf = resid_fwd .- mean(resid_fwd)
            rb = resid_bwd .- mean(resid_bwd)
            denom = norm(rf) * norm(rb)
            resid_corr = denom > eps(Float64) ? dot(rf, rb) / denom : NaN

            # ── Joint GCV ──
            #   nd   = 2 * N_pixels (fwd + bwd jointly)
            #   p_eff = n_chain_params + 4  (a_fwd, b_fwd, a_bwd, b_bwd)
            nd_joint = 2 * length(model_masked)
            p_chain = GaussianFit2D._chain_nparams(r.n, ccfg)
            p_eff = p_chain + 4
            joint_gcv = nd_joint > p_eff ?
                nd_joint / (nd_joint - p_eff)^2 * (rss_fwd + rss_bwd) : Inf

            # ── Joint NRMSE (normalized by total variance) ──
            var_total = var(z_fwd_masked) + var(z_bwd_masked)
            joint_nrmse = (var_total > eps(Float64) && isfinite(rss_fwd + rss_bwd)) ?
                sqrt((rss_fwd + rss_bwd) / (2 * length(model_masked))) / sqrt(var_total) : NaN

            joint_scores[r.n] = joint_gcv

            selected_by_source = (best !== nothing && best.success && best.n == r.n) ? 1 : 0

            push!(rows, (n=r.n, gcv=r.gcv, valid=r.valid,
                         joint_gcv=joint_gcv, joint_nrmse=joint_nrmse,
                         rss_fwd=rss_fwd, rss_bwd=rss_bwd,
                         a_fwd=a_fwd, b_fwd=b_fwd, a_bwd=a_bwd, b_bwd=b_bwd,
                         resid_corr=resid_corr,
                         selected_by_source=selected_by_source))
        end

        # ── Joint selection (lowest joint_gcv) ──
        best_joint_n = -1
        if !isempty(joint_scores)
            best_joint_n = argmin(joint_scores)
        end

        # ── Write rows ──
        for row in rows
            sel_joint = (row.n == best_joint_n) ? 1 : 0
            out_row = Any[fn, row.n, row.gcv, row.valid ? 1 : 0,
                          row.joint_gcv, row.joint_nrmse,
                          row.rss_fwd, row.rss_bwd,
                          row.a_fwd, row.b_fwd, row.a_bwd, row.b_bwd,
                          row.resid_corr, row.selected_by_source, sel_joint]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(out_row), '\t'))
            end
        end

        if isempty(rows)
            # Placeholder row if nothing succeeded
            open(out_tsv, "a") do io
                println(io, "$fn\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
            end
        end
    end

    println("Wrote $out_tsv")
end

main()
