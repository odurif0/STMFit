#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# audit_resolved_lobes.jl — Physical resolvability of adjacent Gaussian lobes
#
# For a fitted chain model with N Gaussian lobes, compute how many adjacent
# lobes are physically resolved. If two adjacent lobes are too close relative
# to their longitudinal widths AND the model has no significant valley between
# their maxima relative to noise, treat them as unresolved and merge them for
# N_resolved.
#
# This is a generic audit: no expected N, no target-specific thresholds, no
# N=6 prior.  It may map any N to any lower N.
#
# Usage:
#   julia test/audit_resolved_lobes.jl \
#       --config config/chitosan.toml \
#       --files 240817_019.sxm,240817_058.sxm \
#       --out results/resolved_lobes_audit/resolved_lobes.tsv \
#       [--sep-threshold 2.0] \
#       [--valley-snr-threshold 3.0] \
#       [--valley-frac-threshold 0.2]
# ──────────────────────────────────────────────────────────────────────────────

using STMMolecularFit, GaussianFit2D
using Printf, TOML, LinearAlgebra

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = ["240817_019.sxm", "240817_058.sxm"]
const OUTDIR = "results/resolved_lobes_audit"
const EPS = 1e-12

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    out_tsv = joinpath(OUTDIR, "resolved_lobes.tsv")
    sep_thresh = 2.0
    valley_snr_thresh = 3.0
    valley_frac_thresh = 0.2
    overrides = Pair{String,Any}[]
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
        elseif arg == "--sep-threshold"
            i < length(args) || error("--sep-threshold requires a float")
            sep_thresh = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--sep-threshold=")
            sep_thresh = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--valley-snr-threshold"
            i < length(args) || error("--valley-snr-threshold requires a float")
            valley_snr_thresh = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--valley-snr-threshold=")
            valley_snr_thresh = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--valley-frac-threshold"
            i < length(args) || error("--valley-frac-threshold requires a float")
            valley_frac_thresh = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--valley-frac-threshold=")
            valley_frac_thresh = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--set"
            i < length(args) || error("--set requires key=value")
            push!(overrides, _parse_override(args[i+1])); i += 2
        elseif startswith(arg, "--set=")
            push!(overrides, _parse_override(split(arg, "=", limit=2)[2])); i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, out_tsv, sep_thresh, valley_snr_thresh, valley_frac_thresh, overrides
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
        intelligent_sweep=false, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg)
    ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

# ══════════════════════════════════════════════════════════════════════════════
# Helpers from audit_chitosan_cases.jl
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

function _effective_best(by_ell, by_circ, criterion)
    best_n, best_score, best_source = 0, Inf, "NA"
    for n in sort(unique(vcat(collect(keys(by_ell)), collect(keys(by_circ)))))
        s_ell = haskey(by_ell, n) ? _score(by_ell[n], criterion) : Inf
        s_circ = haskey(by_circ, n) ? _score(by_circ[n], criterion) : Inf
        source = s_ell <= s_circ ? "ell" : "circ"
        score = min(s_ell, s_circ)
        if score < best_score
            best_n, best_score, best_source = n, score, source
        end
    end
    return best_n, best_source, best_score
end

function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    refined = GaussianFit2D.ChainModelResult[]
    isempty(results_circ) && return refined
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    ac_full = ctx_circ.axisctx_full
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ac_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    ccfg_refine = deepcopy(ccfg_ell)
    ccfg_refine.skip_global = true
    ccfg_refine.max_iter = 50
    ccfg_refine.multistart = 1
    for r_c in results_circ
        r_c.success || continue
        n = r_c.n
        try
            n_prefix = 1 + (ccfg_refine.chain_tilted_baseline ? 2 : 0)
            split_idx = n_prefix + n + GaussianFit2D._chain_spacing_param_count(n, ccfg_refine) + n
            p_init = vcat(r_c.params[1:split_idx], r_c.params[(split_idx+1):end], r_c.params[(split_idx+1):end])
            r_ref = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise,
                n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r_ref.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r_ref.params, n,
                    ac_fit, ccfg_refine; amp_min=r_ref.amp_min, amp_range=r_ref.amp_range)
                GaussianFit2D._finalize_chain_result!(r_ref, zfit, pred, noise,
                    n, n_eff, z, xs, ys, zimg, xfit, yfit, ac_fit, ccfg_refine)
                push!(refined, r_ref)
            end
        catch err
            @warn "elliptical refinement failed" n exception=(err, catch_backtrace())
        end
    end
    return refined
end

function _support_length(ctx)
    hasproperty(ctx, :support_meta) || return ctx.axisctx.tmax - ctx.axisctx.tmin
    return get(ctx.support_meta, :final_support_length_nm, ctx.axisctx.tmax - ctx.axisctx.tmin)
end

# ══════════════════════════════════════════════════════════════════════════════
# Pair-based resolvability metrics
# ══════════════════════════════════════════════════════════════════════════════

"""
    _pair_metrics(r, axisctx, ccfg, noise, sep_thresh, valley_snr_thresh, valley_frac_thresh)

Compute per-adjacent-pair metrics for a successful chain model result.
Returns a vector of named tuples, one per adjacent pair (n-1 pairs for N lobes).
Also returns N_resolved after merging unresolved adjacent pairs.
"""
function _pair_metrics(r, axisctx, ccfg, noise, sep_thresh, valley_snr_thresh, valley_frac_thresh)
    # Default empty return
    empty_pair = (d_nm=NaN, sigma_pair_nm=NaN, sep_sigma=NaN,
                  valley_depth=NaN, valley_snr=NaN, valley_frac=NaN,
                  amp_ratio=NaN, unresolved=false, i=0, j=0)
    if r === nothing || !r.success || r.n <= 1
        return [empty_pair], r.n, "", Float64[], Float64[], Float64[]
    end

    # Decode chain
    _b, feats, ts, _us, spars, _sperps = GaussianFit2D._decode_chain(r.params, r.n, axisctx, ccfg;
        amp_min=r.amp_min, amp_range=r.amp_range)

    # Sort by axial coordinate (ts) — already sorted from _decode_chain
    n = r.n
    pairs = Vector{NamedTuple}()
    unresolved_edges = Tuple{Int,Int}[]  # (i,j) pairs that are unresolved
    min_sep_sigma_arr = Float64[]
    min_valley_snr_arr = Float64[]
    min_valley_frac_arr = Float64[]

    for idx in 1:(n-1)
        i, j = idx, idx+1
        f_i = feats[i]
        f_j = feats[j]

        # Distance between centers
        d_nm = sqrt((f_i.x_nm - f_j.x_nm)^2 + (f_i.y_nm - f_j.y_nm)^2)

        # Mean sigma parallel (longitudinal width)
        sp_i = spars[i]
        sp_j = spars[j]
        sigma_pair_nm = sqrt((sp_i^2 + sp_j^2) / 2.0)
        sep_sigma = d_nm / max(sigma_pair_nm, EPS)

        # Amplitude ratio
        amp_i = f_i.amplitude
        amp_j = f_j.amplitude
        amp_ratio = min(amp_i, amp_j) / max(max(amp_i, amp_j), EPS)

        # --- Valley computation: evaluate model along line between centers ---
        # Sample n_samples points along the line
        n_samples = 64
        xs_line = range(f_i.x_nm, f_j.x_nm, length=n_samples)
        ys_line = range(f_i.y_nm, f_j.y_nm, length=n_samples)

        # Full model values at each sample point
        vals = GaussianFit2D._chain_model_values(
            collect(xs_line), collect(ys_line),
            r.params, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)

        # Peak heights at lobe centers
        peak_i = GaussianFit2D._chain_model_values(
            [f_i.x_nm], [f_i.y_nm],
            r.params, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)[1]
        peak_j = GaussianFit2D._chain_model_values(
            [f_j.x_nm], [f_j.y_nm],
            r.params, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)[1]

        # Valley is the minimum between the two peaks (excluding the endpoints)
        valley_min = minimum(vals[2:end-1])
        min_peak = min(peak_i, peak_j)
        valley_depth = min_peak - valley_min

        valley_snr = valley_depth / max(noise, EPS)
        valley_frac = valley_depth / max(min_peak, EPS)

        # Determine unresolved status
        unresolved = (sep_sigma < sep_thresh) && (valley_snr < valley_snr_thresh || valley_frac < valley_frac_thresh)

        push!(pairs, (d_nm=d_nm, sigma_pair_nm=sigma_pair_nm, sep_sigma=sep_sigma,
                      valley_depth=valley_depth, valley_snr=valley_snr,
                      valley_frac=valley_frac, amp_ratio=amp_ratio,
                      unresolved=unresolved, i=i, j=j))
        if unresolved
            push!(unresolved_edges, (i, j))
            push!(min_sep_sigma_arr, sep_sigma)
            push!(min_valley_snr_arr, valley_snr)
            push!(min_valley_frac_arr, valley_frac)
        end
    end

    # Compute N_resolved via connected components (union-find)
    parent = collect(1:n)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    function union!(a, b)
        ra, rb = find(a), find(b)
        ra != rb && (parent[rb] = ra)
    end
    for (a, b) in unresolved_edges
        union!(a, b)
    end
    roots = Set{Int}()
    for k in 1:n
        push!(roots, find(k))
    end
    n_resolved = length(roots)

    return pairs, n_resolved, unresolved_edges, min_sep_sigma_arr, min_valley_snr_arr, min_valley_frac_arr
end

"""
    _format_pair_details(pairs)

Format pair details as a compact string for TSV output.
Format: "i-j:d/sep/vsnr/vfrac/ampratio/U" for each pair,
  where U=1 if unresolved, 0 otherwise. Separated by ";".
"""
function _format_pair_details(pairs)
    parts = String[]
    for p in pairs
        unresolved_flag = p.unresolved ? 1 : 0
        push!(parts, @sprintf("%d-%d:%.4g/%.4g/%.4g/%.4g/%.4g/%d",
            p.i, p.j, p.d_nm, p.sep_sigma, p.valley_snr, p.valley_frac,
            p.amp_ratio, unresolved_flag))
    end
    return join(parts, ";")
end

function _format_unresolved_pairs(unresolved_edges)
    return join(["$a-$b" for (a, b) in unresolved_edges], ";")
end

# ══════════════════════════════════════════════════════════════════════════════
# Row formatting helpers
# ══════════════════════════════════════════════════════════════════════════════

function _fmt(x)
    x === nothing && return ""
    x isa AbstractString && return x
    x isa Bool && return string(x)
    x isa Integer && return string(x)
    x isa Real && return isfinite(x) ? @sprintf("%.8g", x) : string(x)
    return string(x)
end

function _row_metrics(r, axisctx, ccfg)
    r === nothing && return fill("", 22)
    return Any[r.success, r.valid, r.reason, r.gcv, r.cv_nll_mean, r.bic, r.aicc,
        r.rss, r.chi2_reduced, r.residual_peak_snr, r.overlap, r.kappa_max_adj,
        r.endpoint_overrun_nm, r.mean_spacing_nm, r.spacing_cv, r.sigma_parallel_nm,
        r.sigma_perp_nm]
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    config_file, files, out_tsv, sep_thresh, valley_snr_thresh, valley_frac_thresh, overrides = _parse_cli(ARGS)

    cfg = TOML.parsefile(config_file)
    model = cfg["model"]
    preproc = cfg["preprocessing"]
    for (key, value) in overrides
        model[key] = value
    end
    criterion = get(model, "selection_criterion", "gcv")

    @printf("Resolved-lobe audit\n")
    @printf("  sep_threshold       = %.4f\n", sep_thresh)
    @printf("  valley_snr_threshold = %.4f\n", valley_snr_thresh)
    @printf("  valley_frac_threshold = %.4f\n", valley_frac_thresh)
    @printf("  selection_criterion  = %s\n", criterion)
    @printf("  output               = %s\n", out_tsv)

    mkpath(dirname(out_tsv))

    # ── Write header ──
    header = [
        # File-level identification
        "file", "selected_N", "N_resolved", "source", "criterion", "score", "noise",
        "support_length_nm",
        # Per-pair summary stats
        "min_sep_sigma", "min_valley_snr", "min_valley_frac",
        # Unresolved pair info
        "unresolved_pairs", "n_unresolved_pairs",
        # Pair details compact
        "pair_details",
        # Per-N details
        "N", "is_selected",
        # Model result diagnostics
        "success", "valid", "reason", "gcv", "cv_nll_mean", "bic", "aicc",
        "rss", "chi2_reduced", "residual_peak_snr", "overlap", "kappa_max_adj",
        "endpoint_overrun_nm", "mean_spacing_nm", "spacing_cv",
        "sigma_parallel_nm", "sigma_perp_nm",
    ]

    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end

    for fn in files
        @printf("\n━━━ Auditing %s (resolved lobes) ━━━\n", fn)
        file_out = joinpath(dirname(out_tsv), splitext(fn)[1])
        mkpath(file_out)
        pcfg, ccfg, ccfg_circ = _configs(model, preproc, file_out)
        pcfg.filepath = joinpath(DATA_DIR, fn)

        img = GaussianFit2D.read_sxm(pcfg.filepath)
        results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
        results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx_circ)

        by_ell = _best_by_n(results_ell, criterion)
        by_circ = _best_by_n(results_circ, criterion)
        selected_n, eff_source, eff_score = _effective_best(by_ell, by_circ, criterion)
        noise = ctx_circ.noise
        support_len = _support_length(ctx_circ)

        ns = sort(unique(vcat(collect(keys(by_ell)), collect(keys(by_circ)))))
        @printf("  Selected: N=%d source=%s score=%.6g noise=%.6g support=%.2f\n",
                selected_n, eff_source, eff_score, noise, support_len)

        for n in ns
            r_ell = get(by_ell, n, nothing)
            r_circ = get(by_circ, n, nothing)
            s_ell = r_ell === nothing ? Inf : _score(r_ell, criterion)
            s_circ = r_circ === nothing ? Inf : _score(r_circ, criterion)
            source = s_ell <= s_circ ? "ell" : "circ"
            r = source == "ell" ? r_ell : r_circ
            score = min(s_ell, s_circ)
            is_selected = (n == selected_n)

            # Pick the correct axisctx and ccfg for this result
            result_axisctx = ctx_circ.axisctx
            result_ccfg = source == "ell" ? ccfg : ccfg_circ

            # Compute resolvability metrics
            pairs, n_resolved, unresolved_edges,
                min_sep_arr, min_vsnr_arr, min_vfrac_arr = _pair_metrics(
                    r, result_axisctx, result_ccfg, noise,
                    sep_thresh, valley_snr_thresh, valley_frac_thresh)

            min_sep = isempty(min_sep_arr) ? NaN : minimum(min_sep_arr)
            min_vsnr = isempty(min_vsnr_arr) ? NaN : minimum(min_vsnr_arr)
            min_vfrac = isempty(min_vfrac_arr) ? NaN : minimum(min_vfrac_arr)
            unresolved_pairs_str = _format_unresolved_pairs(unresolved_edges)
            n_unresolved = length(unresolved_edges)
            pair_details_str = _format_pair_details(pairs)

            row = Any[
                fn, selected_n, n_resolved, source, criterion, score, noise,
                support_len,
                min_sep, min_vsnr, min_vfrac,
                unresolved_pairs_str, n_unresolved,
                pair_details_str,
                n, is_selected,
                _row_metrics(r, result_axisctx, result_ccfg)...,
            ]

            open(out_tsv, "a") do io
                println(io, join(_fmt.(row), '\t'))
            end

            if is_selected
                @printf("  N=%d -> N_resolved=%d  unresolved: %s\n",
                        n, n_resolved, isempty(unresolved_pairs_str) ? "none" : unresolved_pairs_str)
            end
        end
    end

    @printf("\nWrote %s\n", out_tsv)
end

main()
