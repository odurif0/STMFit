#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# audit_robust_rescore.jl — External robust-likelihood rescoring of 2D Gaussian
#                            chain fits for chitosan.
#
# This script fits chain models (via one of two pipelines), then rescores each
# successful candidate with Student-t negative log-likelihood over a grid of
# nu (degrees of freedom) values.  For each nu it computes robust IC penalties
# (AICc, BIC) to test whether robust scoring improves N selection without an
# N=6 prior.
#
# Pipelines:
#   ell        — exhaustive elliptical 2D sweep (chain_circular_sigmas=false)
#   effective  — circular sweep + elliptical refinement (matches batch_full
#                candidate space via audit_chitosan_cases.jl patterns)
#
# IMPORTANT: External audit only — this is NOT a fit-time selector.  The fits
# themselves are standard RSS-optimised models; only the scoring changes.
#
# Usage:
#   julia --project=. test/audit_robust_rescore.jl \
#       --config config/chitosan.toml \
#       --pipeline effective \
#       --files 240817_017.sxm,240817_026.sxm \
#       --out results/robust_rescore_audit/robust_rescore.tsv \
#       --nu-grid 1.5,2,3,4,8,1000000 \
#       --criterion robust_aicc
# ──────────────────────────────────────────────────────────────────────────────

using STMMolecularFit, GaussianFit2D
using Printf, TOML, LinearAlgebra, Statistics

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = [
    "240817_017.sxm", "240817_019.sxm", "240817_043.sxm", "240817_058.sxm",
    "240817_002.sxm", "240817_003.sxm", "240817_018.sxm", "240817_039.sxm",
    "240817_060.sxm",
]
const OUTDIR = "results/robust_rescore_audit"
const DEFAULT_NU_GRID = [1.5, 2.0, 3.0, 4.0, 8.0, 1e6]
const EPS = 1e-12

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    out_tsv = joinpath(OUTDIR, "robust_rescore.tsv")
    nu_grid = copy(DEFAULT_NU_GRID)
    criterion = "robust_aicc"
    pipeline = "ell"
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
        elseif arg == "--nu-grid"
            i < length(args) || error("--nu-grid requires comma-separated floats")
            nu_grid = _parse_csv_floats(args[i+1]); i += 2
        elseif startswith(arg, "--nu-grid=")
            nu_grid = _parse_csv_floats(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--criterion"
            i < length(args) || error("--criterion requires a name")
            criterion = args[i+1]; i += 2
        elseif startswith(arg, "--criterion=")
            criterion = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--pipeline"
            i < length(args) || error("--pipeline requires ell or effective")
            pipeline = lowercase(args[i+1]); i += 2
        elseif startswith(arg, "--pipeline=")
            pipeline = lowercase(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--set"
            i < length(args) || error("--set requires key=value")
            push!(overrides, _parse_override(args[i+1])); i += 2
        elseif startswith(arg, "--set=")
            push!(overrides, _parse_override(split(arg, "=", limit=2)[2])); i += 1
        elseif arg in ("-h", "--help")
            _print_help()
            exit(0)
        else
            error("Unknown option: $arg")
        end
    end
    pipeline in ("ell", "effective") || error("--pipeline must be 'ell' or 'effective'")
    return config_file, files, out_tsv, nu_grid, criterion, pipeline, overrides
end

function _parse_csv_floats(s::AbstractString)
    parts = split(s, ',')
    vals = Float64[]
    for p in parts
        v = tryparse(Float64, strip(p))
        v === nothing && error("Invalid float in grid: '$p'")
        push!(vals, v)
    end
    return vals
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
    audit_robust_rescore.jl — Robust-likelihood rescoring for chitosan fits

    Options:
      --config PATH       Config TOML [config/chitosan.toml]
      --files F1,F2,...   Comma-separated .sxm files [default targets + controls]
      --out PATH          Output TSV [$(OUTDIR)/robust_rescore.tsv]
      --nu-grid V,...     Comma-separated nu values [$(join(string.(DEFAULT_NU_GRID), ","))]
      --criterion NAME    Selection: robust_nll | robust_aicc | robust_bic [robust_aicc]
      --pipeline NAME     Candidate pipeline: ell | effective [ell]
      --set KEY=VAL       Override model config
      -h, --help          This message
    """)
end

# ══════════════════════════════════════════════════════════════════════════════
# Config construction — returns (pcfg, ccfg_ell, ccfg_circ)
# Matches audit_chitosan_cases.jl pattern.
# ══════════════════════════════════════════════════════════════════════════════

function _configs(model, preproc, output_dir)
    pcfg = GaussianFit2D.PatternConfig(filepath="", channel="Z", direction="fwd",
        stride=get(preproc, "stride", 1),
        flatten=get(preproc, "flatten", "plane+rows"),
        smooth_radius_px=get(preproc, "smooth_radius_px", 1),
        output_dir=output_dir, no_plot=true)
    ccfg_ell = GaussianFit2D.ChainSweepConfig(n_min=2, n_max=14,
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
        chain_circular_sigmas=false,
        intelligent_sweep=false,
        fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg_ell)
    ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg_ell, ccfg_circ
end

# ══════════════════════════════════════════════════════════════════════════════
# Helpers ported from audit_chitosan_cases.jl
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
    """Refine circular-fit results to elliptical. Returns (refined_results, xfit_ell, yfit_ell, zfit_ell, axisctx_ell)."""
    refined = GaussianFit2D.ChainModelResult[]
    isempty(results_circ) && return refined, Float64[], Float64[], Float64[], nothing
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
    return refined, xfit, yfit, zfit, ac_fit
end

# ══════════════════════════════════════════════════════════════════════════════
# Data loading helper
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
# Rescoring function
# ══════════════════════════════════════════════════════════════════════════════

"""
    _rescore(r, xfit, yfit, zfit, noise, n_eff, axisctx, ccfg, nu)

Compute robust NLL/AICc/BIC for a fitted result at a given nu.
Returns NamedTuple or nothing on failure.
Formulas:
  total_nll = _student_nll(resid, noise, nu)           # sum, not mean
  robust_nll_mean = total_nll / length(zfit)
  pcount = _chain_nparams(n, ccfg)
  robust_aicc = 2*total_nll + 2*pcount + 2*pcount*(pcount+1) / max(n_eff - pcount - 1, 1)
  robust_bic  = 2*total_nll + pcount * log(n_eff)
"""
function _rescore(r, xfit, yfit, zfit, noise, n_eff, axisctx, ccfg, nu)
    r.success || return nothing
    try
        pred = GaussianFit2D._chain_model_values(
            xfit, yfit, r.params, r.n, axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)
        resid = zfit .- pred
        total_nll = GaussianFit2D._student_nll(resid, noise, nu)
        robust_nll_mean = total_nll / length(zfit)
        pcount = GaussianFit2D._chain_nparams(r.n, ccfg)
        n_eff_safe = max(n_eff, pcount + 2)
        robust_aicc = 2*total_nll + 2*pcount +
            (2*pcount*(pcount+1)) / max(n_eff_safe - pcount - 1, 1)
        robust_bic = 2*total_nll + pcount * log(n_eff_safe)
        return (robust_nll_mean=robust_nll_mean,
                robust_aicc=robust_aicc,
                robust_bic=robust_bic,
                pcount=pcount)
    catch
        return nothing
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Formatting
# ══════════════════════════════════════════════════════════════════════════════

function _fmt(x)
    x === nothing && return ""
    x isa AbstractString && return x
    x isa Bool && return string(x)
    x isa Integer && return string(x)
    x isa Real && return isfinite(x) ? @sprintf("%.8g", x) : string(x)
    return string(x)
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    args = _parse_cli(ARGS)
    config_file, files, out_tsv, nu_grid, criterion, pipeline, overrides = args

    cfg = TOML.parsefile(config_file)
    model = cfg["model"]
    preproc = cfg["preprocessing"]
    for (key, value) in overrides
        model[key] = value
    end

    @printf("Robust-rescore audit\n")
    @printf("  config:    %s\n", config_file)
    @printf("  files:     %s\n", join(files, ", "))
    @printf("  pipeline:  %s\n", pipeline)
    @printf("  nu-grid:   %s\n", join([_fmt(v) for v in nu_grid], ", "))
    @printf("  criterion: %s\n", criterion)
    @printf("  output:    %s\n", out_tsv)
    @printf("  NOTE: Rescoring fitted models only — no fits are modified.\n")

    criterion in ("robust_nll", "robust_aicc", "robust_bic") ||
        error("Unknown criterion '$criterion'; use robust_nll, robust_aicc, or robust_bic")

    mkpath(dirname(out_tsv))

    # ── Write header ──
    header = [
        "file", "nu", "criterion", "source", "N", "is_selected",
        "robust_nll_mean", "robust_aicc", "robust_bic",
        "pcount", "n_eff",
        "success", "valid", "reason",
        "gcv", "bic", "aicc", "rss",
        "baseline_selected_N",
    ]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end

    for fn in files
        @printf("\n━━━ %s ━━━\n", fn)
        file_out = joinpath(dirname(out_tsv), splitext(fn)[1])
        mkpath(file_out)
        pcfg, ccfg_ell, ccfg_circ = _configs(model, preproc, file_out)
        pcfg.filepath = joinpath(DATA_DIR, fn)

        # ── Load image and shared data ──
        @printf("  Loading...\n")
        img = GaussianFit2D.read_sxm(pcfg.filepath)
        xs, ys, zimg, mask, x, y, z, noise = _load_data(img, pcfg, ccfg_ell)
        axisctx_full = GaussianFit2D._weighted_roi_axis(x, y, z)

        # Precompute elliptical fit data (needed for both pipelines)
        xfit_ell, yfit_ell, zfit_ell, axisctx_ell, _, _ = GaussianFit2D._chain_fit_data(x, y, z, axisctx_full, ccfg_ell)
        n_eff_ell = max(10, length(zfit_ell) ÷ 9)

        # ── Pipeline-specific candidate construction ──
        if pipeline == "ell"
            # ---- Exhaustive elliptical sweep ----
            @printf("  Elliptical sweep...\n")
            results_all, baseline_best, ctx = GaussianFit2D.chain_gaussian_sweep(
                img, pcfg, ccfg_ell;
                override_data=(xs, ys, zimg, mask, x, y, z, noise),
                override_axisctx=axisctx_full)
            baseline_selected_N = baseline_best.n

            # Build single-source candidate list
            candidates = Vector{NamedTuple}()
            for r in results_all
                r.success || continue
                push!(candidates, (n=r.n, source="ell",
                    result=r, ccfg=ccfg_ell, axisctx=axisctx_ell,
                    xfit=xfit_ell, yfit=yfit_ell, zfit=zfit_ell, n_eff=n_eff_ell))
            end

        else  # pipeline == "effective"
            # ---- Circular sweep + elliptical refinement ----
            @printf("  Circular sweep...\n")
            results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(
                img, pcfg, ccfg_circ;
                override_data=(xs, ys, zimg, mask, x, y, z, noise),
                override_axisctx=axisctx_full)

            # Also need circular fit data for rescoring circ candidates
            xfit_circ, yfit_circ, zfit_circ, axisctx_circ, _, _ = GaussianFit2D._chain_fit_data(x, y, z, axisctx_full, ccfg_circ)
            n_eff_circ = max(10, length(zfit_circ) ÷ 9)

            @printf("  Refining circ→ell (%d circ results)...\n", count(r -> r.success, results_circ))
            results_ell_ref, xfit_ref, yfit_ref, zfit_ref, axisctx_ref = _refine_circ_to_ell(
                results_circ, img, pcfg, ccfg_ell, ctx_circ)

            # Baseline: effective best (circ or refined ell by config criterion)
            sel_criterion = get(model, "selection_criterion", "gcv")
            by_circ = _best_by_n(results_circ, sel_criterion)
            by_ell  = _best_by_n(results_ell_ref, sel_criterion)
            baseline_selected_N, eff_source, eff_score = _effective_best(by_ell, by_circ, sel_criterion)

            # Build candidate list from both sources
            candidates = Vector{NamedTuple}()
            for r in results_circ
                r.success || continue
                push!(candidates, (n=r.n, source="circ",
                    result=r, ccfg=ccfg_circ, axisctx=axisctx_circ,
                    xfit=xfit_circ, yfit=yfit_circ, zfit=zfit_circ, n_eff=n_eff_circ))
            end
            for r in results_ell_ref
                r.success || continue
                push!(candidates, (n=r.n, source="ell_refined",
                    result=r, ccfg=ccfg_ell, axisctx=axisctx_ref,
                    xfit=xfit_ref, yfit=yfit_ref, zfit=zfit_ref, n_eff=n_eff_ell))
            end
        end

        @printf("  Baseline selected N=%d (by %s)\n", baseline_selected_N,
                get(model, "selection_criterion", "gcv"))
        @printf("  %d candidates to rescore\n", length(candidates))

        # ── Rescore all candidates across all nu ──
        cv_results = Vector{NamedTuple}()
        for cand in candidates
            r = cand.result
            for nu in nu_grid
                resc = _rescore(r, cand.xfit, cand.yfit, cand.zfit, noise,
                                cand.n_eff, cand.axisctx, cand.ccfg, nu)
                if resc === nothing
                    push!(cv_results, (nu=nu, n=cand.n, source=cand.source,
                        robust_nll_mean=NaN, robust_aicc=NaN, robust_bic=NaN,
                        pcount=0, n_eff=Int(cand.n_eff),
                        success=false, valid=false, reason="rescore_failed",
                        gcv=NaN, bic=NaN, aicc=NaN, rss=NaN))
                else
                    push!(cv_results, (nu=nu, n=cand.n, source=cand.source,
                        robust_nll_mean=resc.robust_nll_mean,
                        robust_aicc=resc.robust_aicc,
                        robust_bic=resc.robust_bic,
                        pcount=resc.pcount, n_eff=Int(cand.n_eff),
                        success=r.success, valid=r.valid, reason=r.reason,
                        gcv=r.gcv, bic=r.bic, aicc=r.aicc, rss=r.rss))
                end
            end
        end

        # ── Select best candidate per nu (across all sources) ──
        nu_keys = sort(unique([r.nu for r in cv_results]))
        selected_by_nu = Dict{Float64,Tuple{Int,String}}()
        for nu in nu_keys
            subset = filter(r -> r.nu == nu && isfinite(r.robust_nll_mean), cv_results)
            isempty(subset) && continue
            if criterion == "robust_nll"
                best = argmin(r -> r.robust_nll_mean, subset)
            elseif criterion == "robust_aicc"
                best = argmin(r -> r.robust_aicc, subset)
            else
                best = argmin(r -> r.robust_bic, subset)
            end
            selected_by_nu[nu] = (best.n, best.source)
        end

        # ── Write TSV rows ──
        for r in cv_results
            sel_info = get(selected_by_nu, r.nu, (-1, ""))
            is_sel = (r.n == sel_info[1] && r.source == sel_info[2])
            row = Any[
                fn, r.nu, criterion, r.source, r.n, is_sel,
                r.robust_nll_mean, r.robust_aicc, r.robust_bic,
                r.pcount, r.n_eff,
                r.success, r.valid, r.reason,
                r.gcv, r.bic, r.aicc, r.rss,
                baseline_selected_N,
            ]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(row), '\t'))
            end
        end

        # Print summary
        @printf("  Selection by %s:\n", criterion)
        for nu in nu_keys
            sel_n, sel_src = get(selected_by_nu, nu, (-1, ""))
            @printf("    nu=%-10.4g → N=%d source=%s\n", nu, sel_n, sel_src)
        end
    end

    @printf("\nWrote %s\n", out_tsv)
end

main()
