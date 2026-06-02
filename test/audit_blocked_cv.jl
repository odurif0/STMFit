#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# audit_blocked_cv.jl — External spatial blocked cross-validation of 2D
#                        Gaussian chain fits for chitosan.
#
# This script evaluates chain model complexity (N lobes) by spatially blocked
# cross-validation along the molecular axis.  The validation folds are formed
# from CONTIGUOUS blocks of the axial coordinate t, NOT interleaved cyclic
# splits.  This is a more realistic generalisation test: the model must predict
# entire unseen spatial segments.
#
# IMPORTANT: External audit only — this is NOT a fit-time selector.  Results
# should be compared to the GCV/BIC-selected N for diagnostic purposes.
#
# Usage:
#   julia --project=. test/audit_blocked_cv.jl \
#       --config config/chitosan.toml \
#       --files 240817_017.sxm,240817_019.sxm,240817_043.sxm,240817_058.sxm \
#       --out results/blocked_cv_audit/blocked_cv.tsv \
#       --folds 5 \
#       [--n-min 4] [--n-max 10] \
#       [--cv-global]
#       [--set key=value]
# ──────────────────────────────────────────────────────────────────────────────

using STMMolecularFit, GaussianFit2D
using Printf, TOML, LinearAlgebra, Statistics

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = [
    "240817_017.sxm", "240817_019.sxm", "240817_043.sxm", "240817_058.sxm",
    "240817_002.sxm", "240817_003.sxm", "240817_018.sxm", "240817_039.sxm",
    "240817_060.sxm",
]
const OUTDIR = "results/blocked_cv_audit"
const EPS = 1e-12

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    out_tsv = joinpath(OUTDIR, "blocked_cv.tsv")
    folds = 5
    n_min = nothing         # use config defaults
    n_max = nothing
    overrides = Pair{String,Any}[]
    cv_global = false
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
        elseif arg == "--folds"
            i < length(args) || error("--folds requires an integer")
            folds = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--folds=")
            folds = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--n-min"
            i < length(args) || error("--n-min requires an integer")
            n_min = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--n-min=")
            n_min = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--n-max"
            i < length(args) || error("--n-max requires an integer")
            n_max = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--n-max=")
            n_max = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--set"
            i < length(args) || error("--set requires key=value")
            push!(overrides, _parse_override(args[i+1])); i += 2
        elseif startswith(arg, "--set=")
            push!(overrides, _parse_override(split(arg, "=", limit=2)[2])); i += 1
        elseif arg == "--cv-global"
            cv_global = true; i += 1
        elseif arg in ("-h", "--help")
            _print_help()
            exit(0)
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, out_tsv, folds, n_min, n_max, overrides, cv_global
end

function _parse_override(spec::AbstractString)
    parts = split(spec, "=", limit=2)
    length(parts) == 2 || error("Invalid --set '$spec'; expected key=value")
    key, raw = parts
    value = if lowercase(raw) in ("true", "false")
        lowercase(raw) == "true"
    else
        try
            occursin(r"[\.eE]", raw) ? parse(Float64, raw) : parse(Int, raw)
        catch
            raw
        end
    end
    return key => value
end

function _print_help()
    println("""
    audit_blocked_cv.jl — Spatial blocked cross-validation for chitosan fits

    Options:
      --config PATH     Config TOML [config/chitosan.toml]
      --files F1,F2,... Comma-separated .sxm files [default targets + controls]
      --out PATH        Output TSV [$(OUTDIR)/blocked_cv.tsv]
      --folds INT       Number of contiguous blocks [5]
      --n-min INT       Minimum N to evaluate [from config]
      --n-max INT       Maximum N to evaluate [from config]
      --cv-global       Use global+local optimization inside each CV fold
                        (slower; default is local-only for smoke audits)
      --set KEY=VAL     Override model config
      -h, --help        This message
    """)
end

# ══════════════════════════════════════════════════════════════════════════════
# Config construction (same pattern as audit_chitosan_cases.jl)
# ══════════════════════════════════════════════════════════════════════════════

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
        max_iter=get(model, "max_iter", 300),
        multistart=get(model, "multistart", 1),
        cv_folds=get(model, "cv_folds", 5),
        cv_method=get(model, "cv_method", "gcv"),
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
        chain_circular_sigmas=false,  # elliptical (full 2D) sigmas
        intelligent_sweep=false,      # exhaustive linear sweep
        fuse_z_bwd=true)
    return pcfg, ccfg
end

# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

function _score(r, criterion::AbstractString)
    c = lowercase(String(criterion))
    c == "gcv" && return r.gcv
    c == "aicc" && return r.aicc
    c == "cv" && return r.cv_nll_mean
    return r.bic
end

_valid(r, criterion) = r !== nothing && r.success && r.valid && isfinite(_score(r, criterion))

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

function _fmt(x)
    x === nothing && return ""
    x isa AbstractString && return x
    x isa Bool && return string(x)
    x isa Integer && return string(x)
    x isa Real && return isfinite(x) ? @sprintf("%.8g", x) : string(x)
    return string(x)
end

# ══════════════════════════════════════════════════════════════════════════════
# Contiguous-block cross-validation along axial coordinate t
# ══════════════════════════════════════════════════════════════════════════════

"""
    _blocked_cv(xs, ys, zimg, x, y, z, noise, n, axisctx, ccfg, folds)

Perform K-fold block CV along the axial coordinate t.

1. Compute t for each point in (x, y).
2. Sort by t.
3. Split sorted indices into K contiguous blocks of equal size (±1).
4. For each fold k:
   - Train on all blocks except k.
   - Validation on block k.
   - Fit chain model on training data (`_fit_chain_n` with reduced settings).
   - Predict on held-out block, compute mean Student-t NLL per pixel.
5. Return (mean_cv_nll, stderr_cv_nll, n_folds_ok).

Returns (Inf, 0.0, 0) if fewer than 2 folds succeed.
"""
function _blocked_cv(xs, ys, zimg, x, y, z, noise, n::Int, axisctx, ccfg::GaussianFit2D.ChainSweepConfig, folds::Int; cv_global::Bool=false)
    npoints = length(x)
    # Compute axial coordinate t for each point
    t = (x .- axisctx.origin[1]) .* axisctx.axis[1] .+ (y .- axisctx.origin[2]) .* axisctx.axis[2]
    sort_idx = sortperm(t)

    # Contiguous block boundaries
    # Each block has either n÷folds or n÷folds+1 points
    block_ends = Int[0]
    base = npoints ÷ folds
    rem = npoints % folds
    for k in 1:folds
        block_size = base + (k <= rem ? 1 : 0)
        push!(block_ends, block_ends[end] + block_size)
    end
    # block_ends[k-1]+1 : block_ends[k] = indices of block k in sorted order

    # CV config: reduced iterations; global optimization is optional because it
    # is much slower but less dependent on the full-data optimum.
    ccfg_cv = deepcopy(ccfg)
    ccfg_cv.skip_global = !cv_global
    ccfg_cv.max_iter = max(50, ccfg.max_iter ÷ 2)
    ccfg_cv.multistart = 1

    scores = Float64[]
    for k in 1:folds
        # Validation = block k (contiguous in sorted t order)
        val_range = (block_ends[k]+1):block_ends[k+1]
        val_indices = sort_idx[val_range]
        # Train = everything else
        train_indices = setdiff(1:npoints, val_indices)

        length(train_indices) > 10 && length(val_indices) > 5 || continue

        xtrain = x[train_indices]; ytrain = y[train_indices]; ztrain = z[train_indices]
        xval   = x[val_indices];   yval   = y[val_indices];   zval   = z[val_indices]

        # Fit on training data
        r = GaussianFit2D._fit_chain_n(xs, ys, zimg, xtrain, ytrain, ztrain,
                                        noise, n, axisctx, ccfg_cv; starts=1)
        r.success || continue

        # Predict on held-out
        pred = GaussianFit2D._chain_model_values(xval, yval, r.params, n, axisctx, ccfg_cv;
                                                  amp_min=r.amp_min, amp_range=r.amp_range)
        nll = GaussianFit2D._student_nll(zval .- pred, noise, ccfg_cv.student_nu)
        push!(scores, nll / length(xval))
    end

    if length(scores) < 2
        return Inf, 0.0, length(scores)
    end
    mu = mean(scores)
    se = std(scores) / sqrt(length(scores))
    return mu, se, length(scores)
end

# ══════════════════════════════════════════════════════════════════════════════
# Data loading helper (matches chain_gaussian_sweep path)
# ══════════════════════════════════════════════════════════════════════════════

function _load_data(img, cfg, ccfg)
    """Load image data following chain_gaussian_sweep's decision logic."""
    has_bwd = any(c -> lowercase(c.name) == lowercase(cfg.roi_channel) && lowercase(c.direction) == "bwd", img.channels)
    use_fusion = ccfg.fuse_z_bwd && has_bwd
    if use_fusion
        return GaussianFit2D._fused_roi_data(img, cfg)
    else
        return GaussianFit2D._robust_roi_data(img, cfg)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    config_file, files, out_tsv, folds, n_min_opt, n_max_opt, overrides, cv_global = _parse_cli(ARGS)

    cfg = TOML.parsefile(config_file)
    model = cfg["model"]
    preproc = cfg["preprocessing"]
    for (key, value) in overrides
        model[key] = value
    end
    criterion = get(model, "selection_criterion", "gcv")

    @printf("Blocked-CV audit\n")
    @printf("  config:    %s\n", config_file)
    @printf("  folds:     %d (contiguous axial blocks)\n", folds)
    @printf("  files:     %s\n", join(files, ", "))
    @printf("  criterion: %s\n", criterion)
    @printf("  cv_global: %s\n", cv_global)
    @printf("  output:    %s\n", out_tsv)

    mkpath(dirname(out_tsv))

    # ── Write header ──
    header = [
        "file", "N", "is_selected",
        "blocked_cv_nll_mean", "blocked_cv_nll_se", "folds_ok",
        "success", "reason", "gcv", "bic", "aicc", "rss", "valid",
        "baseline_selected_N",
    ]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end

    for fn in files
        @printf("\n━━━ %s ━━━\n", fn)
        file_out = joinpath(dirname(out_tsv), splitext(fn)[1])
        mkpath(file_out)
        pcfg, ccfg = _configs(model, preproc, file_out)
        pcfg.filepath = joinpath(DATA_DIR, fn)

        # ── Full-data fitting (baseline selected_N) ──
        @printf("  Baseline full sweep...\n")
        img = GaussianFit2D.read_sxm(pcfg.filepath)

        # Get data once, reuse for both baseline sweep and CV
        xs, ys, zimg, mask, x, y, z, noise = _load_data(img, pcfg, ccfg)
        axisctx_full = GaussianFit2D._weighted_roi_axis(x, y, z)
        xfit, yfit, zfit, axisctx, fit_keep, support_meta = GaussianFit2D._chain_fit_data(x, y, z, axisctx_full, ccfg)
        axis_length = axisctx.tmax - axisctx.tmin
        spacing_min_eff = GaussianFit2D._effective_spacing_min_nm(ccfg)

        # Baseline: run exhaustive sweep to get selected N
        # Reuse data via override
        baseline_results, baseline_best, _ = GaussianFit2D.chain_gaussian_sweep(
            img, pcfg, ccfg;
            override_data=(xs, ys, zimg, mask, x, y, z, noise),
            override_axisctx=axisctx_full)
        baseline_selected_N = baseline_best.n

        @printf("  Baseline: N=%d (by %s)  support=%.2f nm  noise=%.4g\n",
                baseline_selected_N, criterion, axis_length, noise)

        # ── Feasible N range ──
        n_max_data = max(1, Int(floor(axis_length / max(spacing_min_eff, EPS))) + 1)
        n_min_data = max(2, Int(floor(axis_length / max(ccfg.spacing_max_nm, EPS))))
        n_min_eff = max(ccfg.n_min, n_min_data)
        n_max_eff = min(ccfg.n_max, n_max_data)
        if n_min_opt !== nothing
            n_min_eff = max(n_min_eff, n_min_opt)
        end
        if n_max_opt !== nothing
            n_max_eff = min(n_max_eff, n_max_opt)
        end
        n_min_eff = min(n_min_eff, n_max_eff)

        @printf("  N sweep: %d..%d\n", n_min_eff, n_max_eff)

        # ── Blocked CV for each N ──
        cv_results = Vector{NamedTuple}()
        for n in n_min_eff:n_max_eff
            @printf("  N=%d blocked CV...", n)
            cv_mean, cv_se, folds_ok = _blocked_cv(
                xs, ys, zimg, xfit, yfit, zfit, noise, n, axisctx, ccfg, folds;
                cv_global=cv_global)
            @printf("  mean=%.6g  se=%.6g  folds=%d\n", cv_mean, cv_se, folds_ok)

            # Also try to get the full-fit result for model diagnostics
            # Look up from baseline results by N
            bres = filter(r -> r.n == n, baseline_results)
            full_r = isempty(bres) ? nothing : bres[1]
            full_success = full_r !== nothing && full_r.success
            full_reason = full_r !== nothing ? full_r.reason : ""
            full_gcv = full_r !== nothing ? full_r.gcv : Inf
            full_bic = full_r !== nothing ? full_r.bic : Inf
            full_aicc = full_r !== nothing ? full_r.aicc : Inf
            full_rss = full_r !== nothing ? full_r.rss : Inf
            full_valid = full_r !== nothing && full_r.valid

            push!(cv_results, (n=n, cv_mean=cv_mean, cv_se=cv_se, folds_ok=folds_ok,
                               success=full_success, reason=full_reason,
                               gcv=full_gcv, bic=full_bic, aicc=full_aicc,
                               rss=full_rss, valid=full_valid))
        end

        # ── Select by minimum blocked CV ──
        ok = filter(r -> isfinite(r.cv_mean) && r.folds_ok >= 2, cv_results)
        selected_n = if isempty(ok)
            @printf("  WARNING: no N with valid CV; defaulting to baseline\n")
            baseline_selected_N
        else
            best = argmin(r -> r.cv_mean, ok)
            best.n
        end
        @printf("  CV-selected: N=%d (baseline: N=%d)\n", selected_n, baseline_selected_N)

        # ── Write TSV ──
        for r in cv_results
            is_sel = r.n == selected_n
            row = Any[
                fn, r.n, is_sel,
                r.cv_mean, r.cv_se, r.folds_ok,
                r.success, r.reason, r.gcv, r.bic, r.aicc, r.rss, r.valid,
                baseline_selected_N,
            ]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(row), '\t'))
            end
        end
    end

    @printf("\nWrote %s\n", out_tsv)
end

main()
