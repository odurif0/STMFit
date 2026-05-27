#!/usr/bin/env julia

using STMMolecularFit, GaussianFit2D
using Printf, TOML

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const DEFAULT_FILES = ["240817_017.sxm", "240817_019.sxm", "240817_043.sxm", "240817_058.sxm"]
const OUTDIR = "results/chitosan_case_audit"

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = copy(DEFAULT_FILES)
    out_tsv = joinpath(OUTDIR, "audit_exhaustive.tsv")
    overrides = Pair{String,Any}[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i < length(args) || error("--config requires a file path")
            config_file = args[i + 1]; i += 2
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--files"
            i < length(args) || error("--files requires comma-separated names")
            files = split(args[i + 1], ","); i += 2
        elseif startswith(arg, "--files=")
            files = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--out"
            i < length(args) || error("--out requires a TSV path")
            out_tsv = args[i + 1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--set"
            i < length(args) || error("--set requires key=value")
            push!(overrides, _parse_override(args[i + 1])); i += 2
        elseif startswith(arg, "--set=")
            push!(overrides, _parse_override(split(arg, "=", limit=2)[2])); i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, out_tsv, overrides
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

function _configs(model, preproc, output_dir; exhaustive=true)
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
        intelligent_sweep=!exhaustive, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg)
    ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

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

function _fmt(x)
    x === nothing && return ""
    x isa AbstractString && return x
    x isa Bool && return string(x)
    x isa Integer && return string(x)
    x isa Real && return isfinite(x) ? @sprintf("%.8g", x) : string(x)
    return string(x)
end

function _shape_metrics(r, axisctx, ccfg)
    if r === nothing || !r.success || isempty(r.params) || r.n <= 0
        return fill("", 5)
    end
    _b, feats, ts, _us, _spars, _sperps = GaussianFit2D._decode_chain(r.params, r.n, axisctx, ccfg;
        amp_min=r.amp_min, amp_range=r.amp_range)
    ds = diff(ts)
    amps = [f.amplitude for f in feats]
    return Any[isempty(ds) ? NaN : minimum(ds), isempty(ds) ? NaN : maximum(ds),
        isempty(amps) ? NaN : minimum(amps), isempty(amps) ? NaN : maximum(amps),
        isempty(amps) ? NaN : minimum(amps) / max(maximum(amps), eps(Float64))]
end

function _row_metrics(r, axisctx, ccfg)
    r === nothing && return fill("", 22)
    return Any[r.success, r.valid, r.reason, r.gcv, r.cv_nll_mean, r.bic, r.aicc,
        r.rss, r.chi2_reduced, r.residual_peak_snr, r.overlap, r.kappa_max_adj,
        r.endpoint_overrun_nm, r.mean_spacing_nm, r.spacing_cv, r.sigma_parallel_nm,
        r.sigma_perp_nm, _shape_metrics(r, axisctx, ccfg)...]
end

function main()
    config_file, files, out_tsv, overrides = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    model = cfg["model"]
    preproc = cfg["preprocessing"]
    for (key, value) in overrides
        model[key] = value
    end
    criterion = get(model, "selection_criterion", "gcv")
    mkpath(dirname(out_tsv))
    header = ["file", "criterion", "support_2d_nm", "selected_N_eff", "eff_source", "eff_score", "N",
              "eff_score_N", "winner_N_source",
              "ell_success", "ell_valid", "ell_reason", "ell_gcv", "ell_student_gcv", "ell_bic", "ell_aicc", "ell_rss", "ell_chi2", "ell_resid_snr", "ell_overlap", "ell_kappa", "ell_endpoint_overrun_nm", "ell_mean_spacing_nm", "ell_spacing_cv", "ell_sigma_parallel_nm", "ell_sigma_perp_nm", "ell_min_spacing_nm", "ell_max_spacing_nm", "ell_min_amp", "ell_max_amp", "ell_min_amp_rel",
              "circ_success", "circ_valid", "circ_reason", "circ_gcv", "circ_student_gcv", "circ_bic", "circ_aicc", "circ_rss", "circ_chi2", "circ_resid_snr", "circ_overlap", "circ_kappa", "circ_endpoint_overrun_nm", "circ_mean_spacing_nm", "circ_spacing_cv", "circ_sigma_parallel_nm", "circ_sigma_perp_nm", "circ_min_spacing_nm", "circ_max_spacing_nm", "circ_min_amp", "circ_max_amp", "circ_min_amp_rel"]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end
    for fn in files
        @printf("Auditing %s (exhaustive N sweep)...\n", fn)
        file_out = joinpath(dirname(out_tsv), splitext(fn)[1])
        mkpath(file_out)
        pcfg, ccfg, ccfg_circ = _configs(model, preproc, file_out; exhaustive=true)
        pcfg.filepath = joinpath(DATA_DIR, fn)
        img = GaussianFit2D.read_sxm(pcfg.filepath)
        results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
        results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx_circ)
        by_ell = _best_by_n(results_ell, criterion)
        by_circ = _best_by_n(results_circ, criterion)
        selected_n, eff_source, eff_score = _effective_best(by_ell, by_circ, criterion)
        ns = sort(unique(vcat(collect(keys(by_ell)), collect(keys(by_circ)))))
        for n in ns
            r_ell = get(by_ell, n, nothing)
            r_circ = get(by_circ, n, nothing)
            s_ell = r_ell === nothing ? Inf : _score(r_ell, criterion)
            s_circ = r_circ === nothing ? Inf : _score(r_circ, criterion)
            source = s_ell <= s_circ ? "ell" : "circ"
            row = Any[fn, criterion, _support_length(ctx_circ), selected_n, eff_source, eff_score,
                      n, min(s_ell, s_circ), source,
                      _row_metrics(r_ell, ctx_circ.axisctx, ccfg)...,
                      _row_metrics(r_circ, ctx_circ.axisctx, ccfg_circ)...]
            open(out_tsv, "a") do io
                println(io, join(_fmt.(row), '\t'))
            end
        end
    end
    println("Wrote $out_tsv")
end

main()
