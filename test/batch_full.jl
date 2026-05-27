#!/usr/bin/env julia
# Full 1D + 2D batch pipeline: fits, best-plots, and enriched summary.
# For each file: 1D slide fit + 2D elliptical chain + 2D circular chain + 6-panel plot.
# Configuration loaded from config/chitosan.toml (or --config path/to/calibration.toml).
# Usage:
#   julia --project=. test/batch_full.jl [N_files]
#   julia --project=. test/batch_full.jl [N_files] --config config/my_system.toml
#   julia --project=. test/batch_full.jl [N_files] --chunk i/n

using STMMolecularFit, GaussianFit2D, GaussianFit1D
using DelimitedFiles, Plots, Printf, Statistics, TOML

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
list_sxm_files(dir) = sort([f for f in readdir(dir) if endswith(lowercase(f), ".sxm")])
const TSV = "results/batch_triage_20240817_relaxed.tsv"
const OUTDIR = "results/best_plots"
const EXCLUDE = Set(["240817_001.sxm", "240817_008.sxm", "240817_009.sxm",
                      "240817_010.sxm", "240817_011.sxm", "240817_012.sxm",
                      "240817_013.sxm", "240817_014.sxm", "240817_016.sxm",
                      "240817_020.sxm", "240817_022.sxm", "240817_023.sxm",
                      "240817_025.sxm", "240817_028.sxm", "240817_056.sxm",
                      "240817_057.sxm", "240817_063.sxm",
                      "240817_015.sxm", "240817_027.sxm"])

function _parse_cli(args)
    n_files = 48
    chunk_idx = 1
    chunk_total = 1
    config_file = "config/chitosan.toml"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--chunk"
            i < length(args) || error("--chunk requires i/n, e.g. --chunk 2/4")
            chunk_arg = args[i + 1]
            i += 2
        elseif startswith(arg, "--chunk=")
            chunk_arg = split(arg, "=", limit=2)[2]
            i += 1
        elseif arg == "--config"
            i < length(args) || error("--config requires a file path")
            config_file = args[i + 1]
            i += 2
            continue
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]
            i += 1
            continue
        elseif startswith(arg, "--")
            error("Unknown option: $arg")
        else
            n_files = parse(Int, arg)
            i += 1
            continue
        end

        parts = split(chunk_arg, "/")
        length(parts) == 2 || error("Invalid --chunk '$chunk_arg'; expected i/n, e.g. 2/4")
        chunk_idx = parse(Int, parts[1])
        chunk_total = parse(Int, parts[2])
    end
    chunk_total >= 1 || error("chunk total must be >= 1")
    1 <= chunk_idx <= chunk_total || error("chunk index must satisfy 1 <= i <= n")
    n_files >= 0 || error("N_files must be >= 0")
    return n_files, chunk_idx, chunk_total, config_file
end

const N_FILES, CHUNK_IDX, CHUNK_TOTAL, CONFIG_FILE = _parse_cli(ARGS)
mkpath(OUTDIR)

const FWHM_SIGMA = 2.355
const COLORMAP_RESID = cgrad([:blue, :lightgray, :red])
const SUMMARY_HEADER = [
    "filepath", "status", "classification",
    "N_ell", "N_circ", "N_1D", "N_eff", "eff_source",
    "dN_circ_ell", "dN_1D_ell", "dN_1D_circ",
    "BIC_ell", "BIC_circ", "sBIC_1D",
    "N_raw_ell", "BIC_raw_ell", "raw_valid_ell", "N_raw_circ", "BIC_raw_circ", "raw_valid_circ",
    "valid_ell", "valid_circ", "reason_ell", "reason_circ",
    "chi2_ell", "chi2_circ", "chi2_1D",
    "support_1D_nm", "support_2D_ell_nm", "support_2D_circ_nm", "support_mismatch_ell", "support_mismatch_circ",
    "Nset_ell_10", "Nset_circ_10", "Nset_1D_10", "common_N_10",
    "Nset_ell_hybrid", "Nset_circ_hybrid", "Nset_1D_hybrid", "common_N_hybrid",
    "threshold_ell_hybrid", "threshold_circ_hybrid", "threshold_1D_hybrid",
    "ambiguous_ell", "runnerup_N_ell", "delta_GCV_ell", "delta_GCV_rel_ell",
    "ambiguous_eff", "runnerup_N_eff", "delta_GCV_eff", "delta_GCV_rel_eff",
    "kappa_ell", "kappa_circ", "kappa_1D",
    "best_plot", "file_dir"
]

const GCV_AMBIGUITY_REL_THRESHOLD = 0.05

_bilinear_interp(xs, ys, z, x0, y0) = begin
    ix = clamp(searchsortedlast(xs, x0), 1, length(xs)-1)
    iy = clamp(searchsortedlast(ys, y0), 1, length(ys)-1)
    x1, x2 = xs[ix], xs[ix+1]; y1, y2 = ys[iy], ys[iy+1]
    tx = (x0-x1)/max(x2-x1, 1e-12); ty = (y0-y1)/max(y2-y1, 1e-12)
    return (1-tx)*(1-ty)*z[iy,ix] + tx*(1-ty)*z[iy,ix+1] + (1-tx)*ty*z[iy+1,ix] + tx*ty*z[iy+1,ix+1]
end

function _ellipse!(p, x0, y0, a, b, angle; color=:cyan, alpha=0.3, label="")
    θ = range(0, 2π, length=72)
    cosθ, sinθ = cos.(θ), sin.(θ)
    ca, sa = cos(angle), sin(angle)
    xe = x0 .+ a .* cosθ .* ca .- b .* sinθ .* sa
    ye = y0 .+ a .* cosθ .* sa .+ b .* sinθ .* ca
    plot!(p, xe, ye; color=color, alpha=alpha, label=label, linewidth=1.5)
end

_setstr(v) = isempty(v) ? "{}" : "{" * join(v, ",") * "}"

function _write_summary_row(io, row)
    length(row) == length(SUMMARY_HEADER) || error("summary row has $(length(row)) columns; expected $(length(SUMMARY_HEADER))")
    println(io, join(row, '\t'))
end

function _error_summary_row(fn::AbstractString)
    row = fill("ERR", length(SUMMARY_HEADER))
    row[1] = fn
    row[2] = "error"
    row[3] = "PROBLEMATIC"
    return row
end

function _ok_summary_row(fn, classif, best_ell_sweep, best_circ_sweep, best1d,
                         best_n_eff, eff_source, best_ell_raw, best_circ_raw,
                         qc, ambiguity, outpath, file_dir)
    n1d = best1d === nothing ? "NA" : best1d.n_peaks
    sBIC_1d = best1d === nothing ? "NA" : round(best1d.student_bic, digits=3)
    chi2_1d = best1d === nothing ? "NA" : best1d.chi2_red
    kappa_1d = best1d === nothing ? "NA" : best1d.kappa_max_adj
    dN_1d_ell = best1d === nothing ? "NA" : best1d.n_peaks - best_ell_sweep.n
    dN_1d_circ = best1d === nothing ? "NA" : best1d.n_peaks - best_circ_sweep.n

    return [fn, "ok", classif,
            best_ell_sweep.n, best_circ_sweep.n, n1d,
            best_n_eff, eff_source,
            best_circ_sweep.n - best_ell_sweep.n,
            dN_1d_ell,
            dN_1d_circ,
            round(best_ell_sweep.bic, digits=3), round(best_circ_sweep.bic, digits=3), sBIC_1d,
            best_ell_raw.n, round(best_ell_raw.bic, digits=3), best_ell_raw.valid,
            best_circ_raw.n, round(best_circ_raw.bic, digits=3), best_circ_raw.valid,
            best_ell_sweep.valid, best_circ_sweep.valid, best_ell_sweep.reason, best_circ_sweep.reason,
            best_ell_sweep.chi2_reduced, best_circ_sweep.chi2_reduced, chi2_1d,
            qc.support_1d, qc.support_ell, qc.support_circ, qc.mismatch_ell, qc.mismatch_circ,
            _setstr(qc.nset_ell_10), _setstr(qc.nset_circ_10), _setstr(qc.nset_1d_10), _setstr(qc.common10),
            _setstr(qc.nset_ell_h), _setstr(qc.nset_circ_h), _setstr(qc.nset_1d_h), _setstr(qc.commonh),
            qc.thr_ell, qc.thr_circ, qc.thr_1d,
            ambiguity.amb_ell, ambiguity.runner_ell, ambiguity.dgcv_ell, ambiguity.dgcv_rel_ell,
            ambiguity.amb_eff, ambiguity.runner_eff, ambiguity.dgcv_eff, ambiguity.dgcv_rel_eff,
            best_ell_sweep.kappa_max_adj, best_circ_sweep.kappa_max_adj, kappa_1d,
            outpath, file_dir]
end

function _nset(results, best_score::Float64, threshold::Float64, scorefun)
    ns = Int[]
    for r in results
        getproperty(r, :success) || continue
        getproperty(r, :valid) || continue
        s = scorefun(r)
        isfinite(s) || continue
        s - best_score <= threshold + 1e-9 && push!(ns, getproperty(r, :n))
    end
    return sort(unique(ns))
end

function _scorefun(criterion::AbstractString)
    c = lowercase(String(criterion))
    c == "gcv" && return r -> getproperty(r, :gcv)
    c == "aicc" && return r -> getproperty(r, :aicc)
    c == "cv" && return r -> getproperty(r, :cv_nll_mean)
    return r -> getproperty(r, :bic)
end

_valid_scored(results, scorefun) = [r for r in results if getproperty(r, :success) && getproperty(r, :valid) && isfinite(scorefun(r))]

function _best_by_n(results, scorefun)
    by_n = Dict{Int,Any}()
    for r in _valid_scored(results, scorefun)
        n = getproperty(r, :n)
        if !haskey(by_n, n) || scorefun(r) < scorefun(by_n[n])
            by_n[n] = r
        end
    end
    return by_n
end

function _best_valid_or_best(results, best_raw; criterion="gcv")
    scorefun = _scorefun(criterion)
    valid = _valid_scored(results, scorefun)
    isempty(valid) && return best_raw
    return sort(valid; by=scorefun)[1]
end

function _ambiguity_stats(results, selected_n::Int; criterion="gcv", rel_threshold=GCV_AMBIGUITY_REL_THRESHOLD)
    lowercase(String(criterion)) == "gcv" || return (false, "NA", NaN, NaN)
    by_n = _best_by_n(results, r -> r.gcv)
    selected = get(by_n, selected_n, nothing)
    selected === nothing && return (false, "NA", NaN, NaN)
    others = sort([r for (n, r) in by_n if n != selected_n]; by=r -> r.gcv)
    isempty(others) && return (false, "NA", NaN, NaN)
    runner = first(others)
    delta = runner.gcv - selected.gcv
    rel = delta / max(abs(selected.gcv), eps(Float64))
    return (rel <= rel_threshold, runner.n, delta, rel)
end

function _effective_score_by_n(results_ell, results_circ, scorefun)
    by_n = Dict{Int,Float64}()
    for r in vcat(_valid_scored(results_ell, scorefun), _valid_scored(results_circ, scorefun))
        by_n[r.n] = min(get(by_n, r.n, Inf), scorefun(r))
    end
    return by_n
end

function _effective_ambiguity_stats(results_ell, results_circ, selected_n::Int; rel_threshold=GCV_AMBIGUITY_REL_THRESHOLD)
    by_n = _effective_score_by_n(results_ell, results_circ, r -> r.gcv)
    selected_score = get(by_n, selected_n, Inf)
    isfinite(selected_score) || return (false, "NA", NaN, NaN)
    others = sort([(n, s) for (n, s) in by_n if n != selected_n && isfinite(s)]; by=x -> x[2])
    isempty(others) && return (false, "NA", NaN, NaN)
    runner_n, runner_score = first(others)
    delta = runner_score - selected_score
    rel = delta / max(abs(selected_score), eps(Float64))
    return (rel <= rel_threshold, runner_n, delta, rel)
end

function _select_effective_best(results_ell, results_circ; criterion="gcv")
    """Select best N using effective_criterion(N) = min(ell_criterion(N), circ_criterion(N)).
    Circular model is nested within elliptical → circ metric is a legitimate lower bound.
    Default criterion is GCV (analytical, no refit needed)."""
    scorefun = _scorefun(criterion)
    valid_ell = _valid_scored(results_ell, scorefun)
    valid_circ = _valid_scored(results_circ, scorefun)
    if isempty(valid_ell) && isempty(valid_circ)
        return nothing, 0, "NA"
    end
    ell_by_n = Dict(r.n => r for r in valid_ell)
    circ_by_n = Dict(r.n => r for r in valid_circ)
    all_ns = sort(unique(vcat(collect(keys(ell_by_n)), collect(keys(circ_by_n)))))

    # Find best N by the configured effective criterion only:
    # effective_score(N) = min(score_ell(N), score_circ(N)).
    best_n = 0
    best_score = Inf
    best_source = "NA"
    for n in all_ns
        s_ell = haskey(ell_by_n, n) ? scorefun(ell_by_n[n]) : Inf
        s_circ = haskey(circ_by_n, n) ? scorefun(circ_by_n[n]) : Inf
        eff = min(s_ell, s_circ)
        if eff < best_score
            best_score = eff
            best_n = n
            best_source = s_ell <= s_circ ? "ell" : "circ"
        end
    end

    r_ell = get(ell_by_n, best_n, nothing)
    r_circ = get(circ_by_n, best_n, nothing)
    if best_source == "ell" && r_ell !== nothing
        best_result = r_ell
    elseif best_source == "circ" && r_circ !== nothing
        best_result = r_circ
    elseif r_ell !== nothing
        best_result = r_ell; best_source = "ell"
    elseif r_circ !== nothing
        best_result = r_circ; best_source = "circ"
    else
        return nothing, 0, "NA"
    end
    return best_result, best_n, best_source
end

function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    """For each N fitted by the circular sweep, run a circ→ell LsqFit refinement.
    Returns a vector of elliptical ChainModelResult (one per N).
    NLopt is intentionally skipped — it always diverges from the isotropic start."""
    refined = GaussianFit2D.ChainModelResult[]
    isempty(results_circ) && return refined
    
    # Get fit data once
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    ac_full = ctx_circ.axisctx_full
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ac_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    
    # Config for LsqFit-only refinement
    ccfg_refine = deepcopy(ccfg_ell)
    ccfg_refine.skip_global = true
    ccfg_refine.max_iter = 50
    ccfg_refine.multistart = 1
    
    for r_c in results_circ
        r_c.success || continue
        n = r_c.n
        try
            # Expand circular sigma params to elliptical format.  Works for both
            # per-lobe sigmas and shared K-type sigmas because the sigma block is
            # always the final block in the parameter vector.
            n_prefix = 1 + (ccfg_refine.chain_tilted_baseline ? 2 : 0)
            split_idx = n_prefix + n + GaussianFit2D._chain_spacing_param_count(n, ccfg_refine) + n
            p_init = vcat(r_c.params[1:split_idx],
                          r_c.params[(split_idx+1):end],
                          r_c.params[(split_idx+1):end])
            # Run LsqFit-only from circ warm-start
            r_ref = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise,
                n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r_ref.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r_ref.params, n,
                    ac_fit, ccfg_refine; amp_min=r_ref.amp_min, amp_range=r_ref.amp_range)
                GaussianFit2D._finalize_chain_result!(r_ref, zfit, pred, noise,
                    n, n_eff, z, xs, ys, zimg, xfit, yfit, ac_fit, ccfg_refine)
                push!(refined, r_ref)
            end
        catch
        end
    end
    return refined
end

function _nset_1d(results, best_score::Float64, threshold::Float64)
    ns = Int[]
    for r in results
        r.success || continue
        s = r.student_bic
        isfinite(s) || continue
        s - best_score <= threshold + 1e-9 && push!(ns, r.n_peaks)
    end
    return sort(unique(ns))
end

_intersect3(a, b, c) = sort(collect(intersect(Set(a), Set(b), Set(c))))

function _support_length(ctx)
    hasproperty(ctx, :support_meta) || return ctx.axisctx.tmax - ctx.axisctx.tmin
    return get(ctx.support_meta, :final_support_length_nm, ctx.axisctx.tmax - ctx.axisctx.tmin)
end

function _support_mismatch(l1d::Real, l2d::Real)
    return isfinite(l1d) && isfinite(l2d) && abs(l2d) > 1e-12 ? abs(l1d - l2d) / abs(l2d) : NaN
end

function _classification(best_ell, best_circ, best1d, common10, commonhybrid)
    !best_ell.valid && return "PROBLEMATIC"
    !best_circ.valid && return "PROBLEMATIC"
    best1d === nothing && return best_ell.n == best_circ.n ? "ROBUST_2D_ONLY" : "PROBLEMATIC_1D_FAIL"
    ns = [best_ell.n, best_circ.n, best1d.n_peaks]
    length(unique(ns)) == 1 && return "ROBUST"
    !isempty(common10) && return "TOLERANT_10"
    !isempty(commonhybrid) && return "TOLERANT_HYBRID"
    maximum(ns) - minimum(ns) <= 1 && return "AMBIGU_MINOR"
    return "PROBLEMATIC"
end

function _write_scores(path::String, results, score_label::String, scorefun, selected)
    open(path, "w") do io
        println(io, join(["N", score_label, "delta_raw", "delta_valid", "BIC", "GCV", "AICc", "valid", "raw_best", "selected", "reason", "chi2", "cv_nll_mean", "cv_nll_std"], '\t'))
        vals = [scorefun(r) for r in results if getproperty(r, :success) && isfinite(scorefun(r))]
        best_raw = isempty(vals) ? Inf : minimum(vals)
        valid_vals = [scorefun(r) for r in results if getproperty(r, :success) && getproperty(r, :valid) && isfinite(scorefun(r))]
        best_valid = isempty(valid_vals) ? best_raw : minimum(valid_vals)
        for r in sort(results; by=r -> getproperty(r, :n))
            getproperty(r, :success) || continue
            s = scorefun(r)
            isfinite(s) || continue
            println(io, join([r.n, s, s - best_raw, s - best_valid,
                              getproperty(r, :bic), getproperty(r, :gcv), getproperty(r, :aicc),
                              getproperty(r, :valid), abs(s - best_raw) <= 1e-9, r === selected,
                              getproperty(r, :reason), getproperty(r, :chi2_reduced),
                              getproperty(r, :cv_nll_mean), getproperty(r, :cv_nll_std)], '\t'))
        end
    end
end

function _write_scores_1d(path::String, results)
    open(path, "w") do io
        println(io, join(["N", "sBIC", "delta", "competitive", "chi2", "kappa"], '\t'))
        vals = [r.student_bic for r in results if r.success]
        best = isempty(vals) ? Inf : minimum(vals)
        for r in sort(results; by=r -> r.n_peaks)
            r.success || continue
            println(io, join([r.n_peaks, r.student_bic, r.student_bic - best, r.competitive, r.chi2_red, r.kappa_max_adj], '\t'))
        end
    end
end

function _quality_warnings(n_eff::Integer, n_1d, support_1d::Real, support_2d::Real)
    warnings = String[]
    if n_1d === nothing
        push!(warnings, "1D fit failed; 2D result only")
        mismatch = _support_mismatch(support_1d, support_2d)
        if isfinite(mismatch) && mismatch > 1.0
            push!(warnings, @sprintf("support mismatch: 1D %.1fnm vs 2D %.1fnm", support_1d, support_2d))
        end
        return warnings
    end
    n_delta = n_1d - n_eff
    mismatch = _support_mismatch(support_1d, support_2d)
    if isfinite(mismatch) && mismatch > 1.0
        push!(warnings, @sprintf("support mismatch: 1D %.1fnm vs 2D %.1fnm", support_1d, support_2d))
    end
    if n_delta >= 4
        push!(warnings, "1D overcounts by +$n_delta peaks")
    elseif n_delta <= -4
        push!(warnings, "1D undercounts by $(abs(n_delta)) peaks")
    end
    return warnings
end

function _selection_diagnostics(best_ell_sweep, best_circ_sweep, best1d, best_eff,
                                results_ell, results_circ, fit_1d, slide, ctx_ell, ctx_circ)
    thr_ell = max(10.0, 0.01 * abs(best_ell_sweep.bic))
    thr_circ = max(10.0, 0.01 * abs(best_circ_sweep.bic))
    thr_1d = best1d === nothing ? NaN : max(10.0, 0.01 * abs(best1d.student_bic))

    nset_ell_10 = _nset(results_ell, best_ell_sweep.bic, 10.0, r -> r.bic)
    nset_circ_10 = _nset(results_circ, best_circ_sweep.bic, 10.0, r -> r.bic)
    nset_1d_10 = best1d === nothing ? Int[] : _nset_1d(fit_1d.fit_run.all_results, best1d.student_bic, 10.0)
    common10 = _intersect3(nset_ell_10, nset_circ_10, nset_1d_10)

    nset_ell_h = _nset(results_ell, best_ell_sweep.bic, thr_ell, r -> r.bic)
    nset_circ_h = _nset(results_circ, best_circ_sweep.bic, thr_circ, r -> r.bic)
    nset_1d_h = best1d === nothing ? Int[] : _nset_1d(fit_1d.fit_run.all_results, best1d.student_bic, thr_1d)
    commonh = _intersect3(nset_ell_h, nset_circ_h, nset_1d_h)

    support_1d = slide.support_length_nm
    support_ell = _support_length(ctx_ell)
    support_circ = _support_length(ctx_circ)

    return (;
        thr_ell, thr_circ, thr_1d,
        nset_ell_10, nset_circ_10, nset_1d_10, common10,
        nset_ell_h, nset_circ_h, nset_1d_h, commonh,
        support_1d, support_ell, support_circ,
        mismatch_ell=_support_mismatch(support_1d, support_ell),
        mismatch_circ=_support_mismatch(support_1d, support_circ),
        classif=_classification(best_eff, best_circ_sweep, best1d, common10, commonh),
    )
end

function make_best_plot(best_ell, ctx_ell, ccfg_ell, best_circ, ctx_circ, ccfg_circ, best_1d, x_1d, y_1d, cfg_1d, slide_mode, outpath; warnings=String[])
    """6-panel (2×3) comparison: ell 2D | circ 2D | 1D  |  ell residuals | circ residuals | 1D residuals."""
    n_ell, n_circ = best_ell.n, best_circ.n
    n1d = best_1d === nothing ? 0 : best_1d.n_peaks

    function _add_roi_contour!(p, ctx; color=:white, linewidth=1.6)
        contour!(p, ctx.xs, ctx.ys, Float64.(ctx.mask);
                 levels=[0.5], color=color, linewidth=linewidth,
                 colorbar=false, label="")
    end

    # ── Common helper: add overlays to a 2D panel ──
    function _add_2d_overlays!(p, best, ctx, ccfg, n)
        ax_ctx = ctx.axisctx
        ox, oy = ax_ctx.origin; ax, ay = ax_ctx.axis
        t_all = (ctx.x .- ox) .* ax .+ (ctx.y .- oy) .* ay
        # ROI contour (white)
        _add_roi_contour!(p, ctx; color=:white, linewidth=1.8)
        # Axis line (yellow)
        plot!(p, [ox+minimum(t_all)*ax, ox+maximum(t_all)*ax],
              [oy+minimum(t_all)*ay, oy+maximum(t_all)*ay];
              color=:yellow, linewidth=1.5, label="")
        # Ridge path (green): search ±0.35 nm perpendicular, then smooth
        perp = [-ay, ax]
        t_ridge = range(minimum(t_all), maximum(t_all), length=120)
        hw_ridge = 0.35
        rx, ry = Float64[], Float64[]
        for t_val in t_ridge
            best_v, bx, by = -Inf, 0.0, 0.0
            for u_r in range(-hw_ridge, hw_ridge, step=0.02)
                x0 = ox + t_val*ax + u_r*perp[1]
                y0 = oy + t_val*ay + u_r*perp[2]
                if minimum(ctx.xs) <= x0 <= maximum(ctx.xs) && minimum(ctx.ys) <= y0 <= maximum(ctx.ys)
                    v = _bilinear_interp(ctx.xs, ctx.ys, ctx.zimg, x0, y0)
                    if v > best_v; best_v = v; bx = x0; by = y0; end
                end
            end
            push!(rx, bx); push!(ry, by)
        end
        plot!(p, rx, ry; color=:white, linewidth=2, alpha=0.9, label="")
        # FWHM ellipses/circles (cyan)
        if n > 0 && best.success
            _, feats, _, _, _, _ = GaussianFit2D._decode_chain(best.params, n, ax_ctx, ccfg;
                amp_min=best.amp_min, amp_range=best.amp_range)
            axis_angle = atan(ax, ay)
            xs_c = Float64[]; ys_c = Float64[]
            for f in feats
                a_ell = f.sigma_x_nm * FWHM_SIGMA / 2
                b_ell = f.sigma_y_nm * FWHM_SIGMA / 2
                _ellipse!(p, f.x_nm, f.y_nm, a_ell, b_ell, axis_angle; color=:cyan, alpha=0.4, label="")
                push!(xs_c, f.x_nm); push!(ys_c, f.y_nm)
            end
            scatter!(p, xs_c, ys_c; marker=:cross, markersize=10, color=:red, linewidth=2.5, label="")
        end
    end

    # ── Shared ROI bounds (use ctx_ell as reference) ──
    xs_ell, ys_ell, zimg_ell, mask_ell = ctx_ell.xs, ctx_ell.ys, ctx_ell.zimg, ctx_ell.mask
    roi_rows = [iy for iy in eachindex(ys_ell) if any(mask_ell[iy, :])]
    roi_cols = [ix for ix in eachindex(xs_ell) if any(mask_ell[:, ix])]
    rx = isempty(roi_cols) ? (minimum(xs_ell)-0.3, maximum(xs_ell)+0.3) : (xs_ell[minimum(roi_cols)]-0.3, xs_ell[maximum(roi_cols)]+0.3)
    ry = isempty(roi_rows) ? (minimum(ys_ell)-0.3, maximum(ys_ell)+0.3) : (ys_ell[minimum(roi_rows)]-0.3, ys_ell[maximum(roi_rows)]+0.3)
    xmin, xmax = rx[1], rx[2]; ymin, ymax = ry[1], ry[2]
    if xmin >= xmax; xmin = minimum(xs_ell); xmax = maximum(xs_ell); end
    if ymin >= ymax; ymin = minimum(ys_ell); ymax = maximum(ys_ell); end

    z_clims = (quantile(vec(zimg_ell), 0.10), quantile(vec(zimg_ell), 0.995))

    # ═══════════════════════════════════════════════════════
    # Row 1 (top): heatmaps + fit
    # ═══════════════════════════════════════════════════════

    # ── Panel (1,1): 2D elliptical heatmap ──
    p_ell = heatmap(xs_ell, ys_ell, zimg_ell; aspect_ratio=:equal, colormap=:thermal, clims=z_clims,
                    title="2D ell N=$n_ell", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_ell, (xmin, xmax)); ylims!(p_ell, (ymin, ymax))
    _add_2d_overlays!(p_ell, best_ell, ctx_ell, ccfg_ell, n_ell)

    # ── Panel (1,2): 2D circular heatmap ──
    xs_circ, ys_circ, zimg_circ, mask_circ = ctx_circ.xs, ctx_circ.ys, ctx_circ.zimg, ctx_circ.mask
    p_circ = heatmap(xs_circ, ys_circ, zimg_circ; aspect_ratio=:equal, colormap=:thermal, clims=z_clims,
                     title="2D circ N=$n_circ", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_circ, (xmin, xmax)); ylims!(p_circ, (ymin, ymax))
    _add_2d_overlays!(p_circ, best_circ, ctx_circ, ccfg_circ, n_circ)

    # ── Panel (1,3): 1D data + fit + components ──
    p_1d = plot(x_1d, y_1d; color=:gray, alpha=0.7, label="data", linewidth=1)
    y1d_resid = y_1d .- mean(y_1d)
    if best_1d !== nothing
        y1d_pred = GaussianFit1D.predict_fit(x_1d, best_1d, cfg_1d)
        y1d_resid = y_1d .- y1d_pred
        plot!(p_1d, x_1d, y1d_pred; color=:red, linewidth=2, label="fit N=$n1d")

        centers = GaussianFit1D._params_to_centers(best_1d.popt, n1d)
        comp_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :lime, :teal, :gold]
        asymmetric = cfg_1d.asymmetric_edges && n1d >= 2
        y0 = best_1d.popt[1]
        for (i, c) in enumerate(centers)
            idx = i - 1
            A = GaussianFit1D._get_amplitude(best_1d.popt, idx)
            σ_in = GaussianFit1D._get_sigma(best_1d.popt, idx)
            if asymmetric && (idx == 0 || idx == n1d - 1)
                σ_out = idx == 0 ? best_1d.popt[end-1] : best_1d.popt[end]
                z = x_1d .- c
                s = idx == 0 ? (z .< 0) .* σ_out .+ (z .>= 0) .* σ_in :
                               (z .< 0) .* σ_in .+ (z .>= 0) .* σ_out
                y_comp = y0 .+ A .* exp.(-0.5 .* (z ./ s).^2)
            else
                y_comp = y0 .+ A .* exp.(-0.5 .* ((x_1d .- c) ./ max(σ_in, 1e-9)).^2)
            end
            col = comp_colors[mod1(i, length(comp_colors))]
            plot!(p_1d, x_1d, y_comp; color=col, alpha=0.55, linestyle=:dash, linewidth=1.5, label="")
        end
        title!(p_1d, "1D N=$n1d ΔsBIC=0 (sBIC=$(round(best_1d.student_bic, digits=0)))")
    else
        title!(p_1d, "1D fit failed")
    end
    xlabel!(p_1d, "position (nm)"); ylabel!(p_1d, "intensity")

    # ═══════════════════════════════════════════════════════
    # Row 2 (bottom): residuals
    # ═══════════════════════════════════════════════════════

    # ── Panel (2,1): 2D elliptical residuals ──
    pred_img_ell = zeros(size(zimg_ell))
    ax_ell = ctx_ell.axisctx
    for iy in eachindex(ys_ell), ix in eachindex(xs_ell)
        pred_img_ell[iy, ix] = GaussianFit2D._chain_model_values([xs_ell[ix]], [ys_ell[iy]], best_ell.params, n_ell, ax_ell, ccfg_ell; amp_min=best_ell.amp_min, amp_range=best_ell.amp_range)[1]
    end
    noise_ell = ctx_ell.noise
    resid_ell = (zimg_ell .- pred_img_ell) .* Float64.(mask_ell) ./ max(noise_ell, 1e-12)
    p_res_ell = heatmap(xs_ell, ys_ell, resid_ell; aspect_ratio=:equal,
                        colormap=COLORMAP_RESID, clims=(-3, 3),
                        title="2D ell residuals", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_res_ell, (xmin, xmax)); ylims!(p_res_ell, (ymin, ymax))
    _add_roi_contour!(p_res_ell, ctx_ell; color=:black, linewidth=1.4)

    # ── Panel (2,2): 2D circular residuals ──
    pred_img_circ = zeros(size(zimg_circ))
    ax_circ = ctx_circ.axisctx
    for iy in eachindex(ys_circ), ix in eachindex(xs_circ)
        pred_img_circ[iy, ix] = GaussianFit2D._chain_model_values([xs_circ[ix]], [ys_circ[iy]], best_circ.params, n_circ, ax_circ, ccfg_circ; amp_min=best_circ.amp_min, amp_range=best_circ.amp_range)[1]
    end
    noise_circ = ctx_circ.noise
    resid_circ = (zimg_circ .- pred_img_circ) .* Float64.(mask_circ) ./ max(noise_circ, 1e-12)
    p_res_circ = heatmap(xs_circ, ys_circ, resid_circ; aspect_ratio=:equal,
                         colormap=COLORMAP_RESID, clims=(-3, 3),
                         title="2D circ residuals", xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
    xlims!(p_res_circ, (xmin, xmax)); ylims!(p_res_circ, (ymin, ymax))
    _add_roi_contour!(p_res_circ, ctx_circ; color=:black, linewidth=1.4)

    # ── Panel (2,3): 1D residuals ──
    p_res_1d = plot(x_1d, y1d_resid; color=:red, label="", linewidth=1)
    hline!(p_res_1d, [0]; color=:gray, linestyle=:dash, label="")
    xlabel!(p_res_1d, "position (nm)"); ylabel!(p_res_1d, "residual")
    title!(p_res_1d, "1D residuals  σ=$(round(std(y1d_resid), digits=5))")

    # ── Global title ──
    bic_1d_str = best_1d === nothing ? "NA" : string(round(best_1d.student_bic, digits=0))
    title_str = "elliptical β=$(round(best_ell.bic, digits=0)) vs circular β=$(round(best_circ.bic, digits=0)) vs 1D β=$bic_1d_str | slide: $slide_mode"
    if !isempty(warnings)
        title_str *= " | ⚠ " * join(warnings, "; ")
    end

    l = @layout grid(2, 3, heights=[0.65, 0.35])
    fig = plot(p_ell, p_circ, p_1d, p_res_ell, p_res_circ, p_res_1d;
               layout=l, size=(2400, 800),
               plot_title=title_str, plot_titlefontsize=8,
               left_margin=0Plots.mm, right_margin=0Plots.mm,
               top_margin=1Plots.mm, bottom_margin=0Plots.mm,
               margin=0.5Plots.mm)
    savefig(fig, outpath)
end

# ── Load candidate files: prefer triage TSV when present, otherwise standalone discovery ──
cands = Tuple{String,Int,Float64}[]
if isfile(TSV)
    data = readdlm(TSV)
    for i in 2:size(data, 1)
        data[i, 5] isa Bool && data[i, 5] || continue
        push!(cands, (string(data[i, 1]), Int(data[i, 3]), Float64(data[i, 4])))
    end
    sort!(cands, by = x -> x[3])
else
    @warn "Triage TSV not found: $TSV; discovering SXM files directly from $DATA_DIR"
    cands = [(fn, 0, Inf) for fn in list_sxm_files(DATA_DIR)]
end
# Exclude files that are not chitosan chains
cands = [(fn, n, bic) for (fn, n, bic) in cands if !(fn in EXCLUDE)]
to_process_base = cands[1:min(N_FILES, length(cands))]
to_process_all = CHUNK_TOTAL == 1 ? to_process_base : [f for (i, f) in enumerate(to_process_base) if mod1(i, CHUNK_TOTAL) == CHUNK_IDX]
if CHUNK_TOTAL > 1
    @printf("Chunk %d/%d: %d of %d selected files\n", CHUNK_IDX, CHUNK_TOTAL, length(to_process_all), length(to_process_base))
end

# ── Skip already-processed files (plot exists + enriched summary has ok row) ──
already_done = Set{String}()
summary_name = CHUNK_TOTAL == 1 ? "summary_overlap060_hard.tsv" : @sprintf("summary_overlap060_hard_chunk%02dof%02d.tsv", CHUNK_IDX, CHUNK_TOTAL)
summary_file = joinpath(OUTDIR, summary_name)
if isfile(summary_file)
    for (i, line) in enumerate(eachline(summary_file))
        i == 1 && continue
        isempty(strip(line)) && continue
        parts = split(line, '\t')
        length(parts) >= 2 || continue
        fn = parts[1]
        status = parts[2]
        png = joinpath(OUTDIR, replace(fn, r"\.sxm$"i => "_best.png"))
        if status == "ok" && isfile(png)
            push!(already_done, fn)
        end
    end
end
to_process = [f for f in to_process_all if !(f[1] in already_done)]
if length(already_done) > 0
    println("  $(length(already_done)) files already done, $(length(to_process)) remaining")
end

# ── Load calibration ──
@info "Loading calibration from $CONFIG_FILE"
cfg_toml = TOML.parsefile(CONFIG_FILE)
model = cfg_toml["model"]
preproc = cfg_toml["preprocessing"]

const SIGMA_MIN_HARMONIZED_NM = model["sigma_parallel_min_nm"]
const SIGMA_MAX_HARMONIZED_NM = model["sigma_parallel_max_nm"]

# ── 2D config (from calibration) ──
pcfg = GaussianFit2D.PatternConfig(filepath="", channel="Z", direction="fwd",
    stride=get(preproc, "stride", 1),
    flatten=get(preproc, "flatten", "plane+rows"),
    smooth_radius_px=get(preproc, "smooth_radius_px", 1),
    output_dir=OUTDIR, no_plot=true)
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
    sigma_parallel_min_nm=SIGMA_MIN_HARMONIZED_NM,
    sigma_parallel_max_nm=SIGMA_MAX_HARMONIZED_NM,
    sigma_perp_min_nm=model["sigma_parallel_min_nm"],
    sigma_perp_max_nm=model["sigma_parallel_max_nm"],
    kappa_max=get(model, "kappa_max", 10.0),
    kappa_weight=get(model, "kappa_weight", 1.0),
    min_amplitude_fraction=get(model, "min_amplitude_fraction", 0.3),
    shared_sigma_types=get(model, "shared_sigma_types", 0),
    chain_spacing_model=get(model, "chain_spacing_model", "free"),
    chain_tilted_baseline=get(model, "chain_tilted_baseline", true),
    intelligent_sweep=true, fuse_z_bwd=true)
# Circular 2D config (same settings, circular sigmas)
ccfg_circ = deepcopy(ccfg)
ccfg_circ.chain_circular_sigmas = true

# ── 1D config (from calibration) ──
scfg = STMMolecularFit.SlideConfig(
    width_nm=get(model, "fit_width_nm", 0.70),
    support_noise_k=model["support_noise_k"],
    support_padding_nm=model["support_padding_nm"],
    output_dir=OUTDIR, no_plot=true)
fcfg = STMMolecularFit.FitSlideConfig(
    min_spacing=model["spacing_min_nm"], max_spacing=model["spacing_max_nm"],
    max_overlap=model["max_overlap"], output_dir=OUTDIR, no_plot=true)

# ── Process ──
# Write header once if summary file is new
if !isfile(summary_file)
    open(summary_file, "w") do io
        println(io, join(SUMMARY_HEADER, '\t'))
    end
end

println("Generating best plots in $OUTDIR/ ...")
ntot = length(to_process)
t0 = time()
summary_lock = ReentrantLock()
plot_lock = ReentrantLock()
n_parallel = min(4, Threads.nthreads(), ntot)
println("Processing $ntot files ($n_parallel parallel workers)...")
flush(stdout)

processed = Threads.Atomic{Int}(0)

Threads.@threads for idx in 1:ntot
    local fn, n2d_tsv, bic2d_tsv = to_process[idx]
    local fp = joinpath(DATA_DIR, fn)
    local file_base = splitext(fn)[1]
    local file_dir = joinpath(OUTDIR, file_base)
    mkpath(file_dir)
    local outpath = joinpath(OUTDIR, replace(fn, r"\.sxm$"i => "_best.png"))
    
    Threads.atomic_add!(processed, 1)
    local done = processed[]
    local pct = done / ntot * 100
    local elapsed = time() - t0
    local eta = done > 0 ? elapsed / done * (ntot - done) : 0
    local eta_min = eta > 0 ? div(round(Int, eta), 60) : 0
    local eta_sec = eta > 0 ? round(Int, eta) % 60 : 0
    @printf("[%2d/%2d %3.0f%%  ETA %d:%02d] %-24s\n", done, ntot, pct, eta_min, eta_sec, fn)
    flush(stdout)
    try
        # 1D fit
        img = STMMolecularFit.read_sxm(fp)
        scfg_file = deepcopy(scfg); scfg_file.output_dir = file_dir
        fcfg_file = deepcopy(fcfg); fcfg_file.output_dir = file_dir
        slide = STMMolecularFit.extract_slide(img, scfg_file)
        if isempty(slide.x) || isempty(slide.y)
            @warn("$fn: empty slide profile (no molecule found?), skipping")
            lock(summary_lock) do
                open(summary_file, "a") do io
                    println(io, join([fn, "no_molecule", "0", "", "", "", "", "", "", "", ""], '\t'))
                end
            end
        else
        fit_1d = nothing
        best1d = nothing
        x_1d, y_1d = slide.x, slide.y
        cfg_1d = fcfg_file
        try
            fit_1d = STMMolecularFit.fit_slide(slide, fcfg_file)
            best1d = GaussianFit1D.best_result(fit_1d.fit_run)
            x_1d, y_1d = fit_1d.fit_run.x, fit_1d.fit_run.y
            cfg_1d = fit_1d.fit_run.cfg
        catch e1d
            @warn "$fn: 1D fit failed; continuing with 2D-only result" reason=sprint(showerror, e1d)
        end

        # 2D circular sweep (primary: always converges)
        img2d = GaussianFit2D.read_sxm(fp)
        pcfg_file = deepcopy(pcfg); pcfg_file.filepath = fp; pcfg_file.output_dir = file_dir
        results_circ, best_circ_raw, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img2d, pcfg_file, ccfg_circ)
        criterion = get(model, "selection_criterion", "gcv")
        scorefun = _scorefun(criterion)
        score_label = uppercase(String(criterion))
        best_circ_sweep = _best_valid_or_best(results_circ, best_circ_raw; criterion=criterion)

        # circ→ell LsqFit refinement at each N (replaces NLopt elliptical sweep)
        # NLopt is intentionally excluded — it always diverges in 33D sigma space.
        results_ell = _refine_circ_to_ell(results_circ, img2d, pcfg_file, ccfg, ctx_circ)
        ctx_ell = ctx_circ  # circ→ell refinement shares the circular context (same data/axis)
        best_ell_raw = isempty(results_ell) ? best_circ_raw : results_ell[1]
        best_ell_sweep = isempty(results_ell) ? best_circ_sweep : _best_valid_or_best(results_ell, best_ell_raw; criterion=criterion)

        # ── Selection: choose best valid models by configured criterion (default GCV) ──
        # Effective best uses min(ell_score, circ_score) per N.
        best_eff, best_n_eff, eff_source = _select_effective_best(results_ell, results_circ; criterion=criterion)
        if best_eff === nothing
            best_eff = best_ell_sweep; best_n_eff = best_ell_sweep.n; eff_source = "ell"
        end

        # ── QC diagnostics only: ambiguity/support/tolerant N sets do not alter selection ──
        amb_ell, runner_ell, dgcv_ell, dgcv_rel_ell = _ambiguity_stats(results_ell, best_ell_sweep.n; criterion=criterion)
        amb_eff, runner_eff, dgcv_eff, dgcv_rel_eff = _effective_ambiguity_stats(results_ell, results_circ, best_n_eff)
        ambiguity = (; amb_ell, runner_ell, dgcv_ell, dgcv_rel_ell,
                       amb_eff, runner_eff, dgcv_eff, dgcv_rel_eff)

        n1d_val = best1d === nothing ? "NA" : best1d.n_peaks
        d1d_eff_val = best1d === nothing ? "NA" : best1d.n_peaks - best_n_eff
        println("Neff=$(best_n_eff) Nell=$(best_ell_sweep.n) Ncirc=$(best_circ_sweep.n) N1D=$(n1d_val) Δ1D-eff=$(d1d_eff_val) ✓")

        _write_scores(joinpath(file_dir, "ell_scores.tsv"), results_ell, score_label, scorefun, best_ell_sweep)
        _write_scores(joinpath(file_dir, "circ_scores.tsv"), results_circ, score_label, scorefun, best_circ_sweep)
        if fit_1d !== nothing
            _write_scores_1d(joinpath(file_dir, "fit_1d_scores.tsv"), fit_1d.fit_run.all_results)
        end

        qc = _selection_diagnostics(best_ell_sweep, best_circ_sweep, best1d, best_eff,
                                    results_ell, results_circ, fit_1d, slide, ctx_ell, ctx_circ)
        classif = qc.classif
        plot_warnings = _quality_warnings(best_n_eff, best1d === nothing ? nothing : best1d.n_peaks, qc.support_1d, qc.support_ell)
        amb_ell && push!(plot_warnings, @sprintf("ambiguous ell GCV: selected N=%d; second best N=%s (ΔGCV=%.1f%%)", best_ell_sweep.n, string(runner_ell), 100 * dgcv_rel_ell))
        amb_eff && push!(plot_warnings, @sprintf("ambiguous eff GCV: selected N=%d; second best N=%s (ΔGCV=%.1f%%)", best_n_eff, string(runner_eff), 100 * dgcv_rel_eff))

        # 6-panel combined plot. GR/Plots is not thread-safe, so serialize savefig.
        lock(plot_lock) do
            make_best_plot(best_ell_sweep, ctx_circ, ccfg, best_circ_sweep, ctx_circ, ccfg_circ,
                           best1d, x_1d, y_1d, cfg_1d, string(scfg.slide_mode), outpath;
                           warnings=plot_warnings)
        end

        # Write summary (thread-safe)
        row = _ok_summary_row(fn, classif, best_ell_sweep, best_circ_sweep, best1d,
                              best_n_eff, eff_source, best_ell_raw, best_circ_raw,
                              qc, ambiguity, outpath, file_dir)
        lock(summary_lock) do
            open(summary_file, "a") do io
                _write_summary_row(io, row)
        end
        end  # lock
        end  # else
    catch e
        msg = sprint(showerror, e)
        msg_short = length(msg) <= 80 ? msg : msg[1:min(end, 80)] * "..."
        println("FAILED: $msg_short")
        lock(summary_lock) do
            open(summary_file, "a") do io
                _write_summary_row(io, _error_summary_row(fn))
        end
        end  # lock
    end
end

# ── Count matches ──
n_done = length(already_done) + length(to_process)
n_success = n_done  # rough estimate; for precise count, re-read summary
println("\nDone! Processed $(length(to_process)) new files in $OUTDIR/")
println("Summary (appended): $summary_file")
