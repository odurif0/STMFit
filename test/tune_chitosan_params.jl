#!/usr/bin/env julia

using STMMolecularFit, GaussianFit2D
using Printf, TOML, DelimitedFiles

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const OUTDIR = "results/param_tuning"

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    # Calibration/tuning defaults intentionally exclude suspected artefact or
    # poor-quality files: 029,030,031,032,034,035,037,038,051.
    files = ["240817_003.sxm", "240817_017.sxm", "240817_019.sxm", "240817_043.sxm",
             "240817_058.sxm", "240817_060.sxm", "240817_061.sxm"]
    out_tsv = joinpath(OUTDIR, "tuning_summary.tsv")
    variants_filter = String[]
    global_maxtime = NaN
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
        elseif arg == "--variants"
            i < length(args) || error("--variants requires comma-separated variant names")
            variants_filter = split(args[i + 1], ","); i += 2
        elseif startswith(arg, "--variants=")
            variants_filter = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--global-maxtime"
            i < length(args) || error("--global-maxtime requires seconds")
            global_maxtime = parse(Float64, args[i + 1]); i += 2
        elseif startswith(arg, "--global-maxtime=")
            global_maxtime = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, out_tsv, variants_filter, global_maxtime
end

function _variant_specs()
    return [
        ("baseline", Dict{String,Any}()),
        ("support_k_2p0", Dict("support_noise_k" => 2.0)),
        ("support_k_3p0", Dict("support_noise_k" => 3.0)),
        ("support_k_3p5", Dict("support_noise_k" => 3.5)),
        ("fit_width_0p10", Dict("fit_width_nm" => 0.10)),
        ("fit_width_0p11", Dict("fit_width_nm" => 0.11)),
        ("fit_width_0p12", Dict("fit_width_nm" => 0.12)),
        ("fit_width_0p13", Dict("fit_width_nm" => 0.13)),
        ("fit_width_0p14", Dict("fit_width_nm" => 0.14)),
        ("fit_width_0p16", Dict("fit_width_nm" => 0.16)),
        ("fit_width_0p18", Dict("fit_width_nm" => 0.18)),
        ("fit_width_0p20", Dict("fit_width_nm" => 0.20)),
        ("fit_width_0p25", Dict("fit_width_nm" => 0.25)),
        ("sigma_min_0p15", Dict("sigma_parallel_min_nm" => 0.15)),
        ("sigma_min_0p22", Dict("sigma_parallel_min_nm" => 0.22)),
        ("sigma_max_0p45", Dict("sigma_parallel_max_nm" => 0.45)),
        ("sigma_max_0p55", Dict("sigma_parallel_max_nm" => 0.55)),
        ("sigma_0p17_0p45", Dict("sigma_parallel_min_nm" => 0.17, "sigma_parallel_max_nm" => 0.45)),
        ("sigma_0p22_0p45", Dict("sigma_parallel_min_nm" => 0.22, "sigma_parallel_max_nm" => 0.45)),
        ("sigma_0p17_0p55", Dict("sigma_parallel_min_nm" => 0.17, "sigma_parallel_max_nm" => 0.55)),
        ("spacing_0p40_0p75", Dict("spacing_min_nm" => 0.40, "spacing_max_nm" => 0.75)),
        ("spacing_0p45_0p75", Dict("spacing_min_nm" => 0.45, "spacing_max_nm" => 0.75)),
        ("spacing_0p35_0p70", Dict("spacing_min_nm" => 0.35, "spacing_max_nm" => 0.70)),
        ("spacing_0p40_0p70", Dict("spacing_min_nm" => 0.40, "spacing_max_nm" => 0.70)),
        ("spacing_0p45_0p70", Dict("spacing_min_nm" => 0.45, "spacing_max_nm" => 0.70)),
        ("spacing_0p30_0p80", Dict("spacing_min_nm" => 0.30, "spacing_max_nm" => 0.80)),
        ("amp_min_0p20", Dict("min_amplitude_fraction" => 0.20)),
        ("amp_min_0p35", Dict("min_amplitude_fraction" => 0.35)),
        ("amp_min_0p40", Dict("min_amplitude_fraction" => 0.40)),
        ("overlap_0p50", Dict("max_overlap" => 0.50)),
        ("overlap_0p55", Dict("max_overlap" => 0.55)),
        ("overlap_0p70", Dict("max_overlap" => 0.70)),
        ("amp035_width012", Dict("min_amplitude_fraction" => 0.35, "fit_width_nm" => 0.12)),
        ("amp035_width013", Dict("min_amplitude_fraction" => 0.35, "fit_width_nm" => 0.13)),
        ("amp035_overlap055", Dict("min_amplitude_fraction" => 0.35, "max_overlap" => 0.55)),
        ("amp035_sigma022", Dict("min_amplitude_fraction" => 0.35, "sigma_parallel_min_nm" => 0.22)),
        ("amp035_k3p0", Dict("min_amplitude_fraction" => 0.35, "support_noise_k" => 3.0)),
        ("width012_overlap055", Dict("fit_width_nm" => 0.12, "max_overlap" => 0.55)),
        ("amp035_width012_overlap055", Dict("min_amplitude_fraction" => 0.35, "fit_width_nm" => 0.12, "max_overlap" => 0.55)),
        ("amp040_sigma022", Dict("min_amplitude_fraction" => 0.40, "sigma_parallel_min_nm" => 0.22)),
        ("amp040_overlap050", Dict("min_amplitude_fraction" => 0.40, "max_overlap" => 0.50)),
        ("amp040_overlap055", Dict("min_amplitude_fraction" => 0.40, "max_overlap" => 0.55)),
        ("amp040_width010", Dict("min_amplitude_fraction" => 0.40, "fit_width_nm" => 0.10)),
        ("amp035_sigma022_overlap055", Dict("min_amplitude_fraction" => 0.35, "sigma_parallel_min_nm" => 0.22, "max_overlap" => 0.55)),
        ("sigma022_overlap050", Dict("sigma_parallel_min_nm" => 0.22, "max_overlap" => 0.50)),
        ("sigma022_overlap055", Dict("sigma_parallel_min_nm" => 0.22, "max_overlap" => 0.55)),
        ("kappa_6", Dict("kappa_max" => 6.0)),
        ("kappa_4", Dict("kappa_max" => 4.0)),
        ("kappa_9", Dict("kappa_max" => 9.0)),
        ("kappa_9p25", Dict("kappa_max" => 9.25)),
        ("kappa_9p5", Dict("kappa_max" => 9.5)),
        ("kappa_9p75", Dict("kappa_max" => 9.75)),
        ("kappa_9p85", Dict("kappa_max" => 9.85)),
        ("kappa_9p90", Dict("kappa_max" => 9.90)),
        ("kappa_9p95", Dict("kappa_max" => 9.95)),
        ("kappa_10", Dict("kappa_max" => 10.0)),
        ("kappa_11", Dict("kappa_max" => 11.0)),
        ("kappa_12", Dict("kappa_max" => 12.0)),
        ("kappa_off", Dict("kappa_max" => 0.0)),
        ("amp_min_0p31", Dict("min_amplitude_fraction" => 0.31)),
        ("amp_min_0p32", Dict("min_amplitude_fraction" => 0.32)),
        ("amp_min_0p33", Dict("min_amplitude_fraction" => 0.33)),
        ("spacing_0p37_0p75", Dict("spacing_min_nm" => 0.37, "spacing_max_nm" => 0.75)),
        ("spacing_0p38_0p75", Dict("spacing_min_nm" => 0.38, "spacing_max_nm" => 0.75)),
        ("sigma_min_0p20", Dict("sigma_parallel_min_nm" => 0.20)),
        ("sigma_min_0p205", Dict("sigma_parallel_min_nm" => 0.205)),
        ("pad_0p10", Dict("support_padding_nm" => 0.10)),
        ("pad_0p15", Dict("support_padding_nm" => 0.15)),
        ("pad_0p18", Dict("support_padding_nm" => 0.18)),
        ("pad_0p20", Dict("support_padding_nm" => 0.20)),
        ("pad_0p22", Dict("support_padding_nm" => 0.22)),
        ("pad_0p25", Dict("support_padding_nm" => 0.25)),
        ("pad_0p30", Dict("support_padding_nm" => 0.30)),
        ("pad_0p40", Dict("support_padding_nm" => 0.40)),
        ("support_k_2p55", Dict("support_noise_k" => 2.55)),
        ("support_k_2p60", Dict("support_noise_k" => 2.60)),
        ("support_k_2p65", Dict("support_noise_k" => 2.65)),
        ("support_k_2p70", Dict("support_noise_k" => 2.70)),
        ("support_k_2p75", Dict("support_noise_k" => 2.75)),
        ("support_k_3p0", Dict("support_noise_k" => 3.0)),
        ("support_k_3p25", Dict("support_noise_k" => 3.25)),
        ("baseline_q_0p05", Dict("support_baseline_quantile" => 0.05)),
        ("baseline_q_0p08", Dict("support_baseline_quantile" => 0.08)),
        ("baseline_q_0p12", Dict("support_baseline_quantile" => 0.12)),
        ("baseline_q_0p15", Dict("support_baseline_quantile" => 0.15)),
        # Physically reasonable, mildly constrained candidates for chitosan/Cu(100)
        ("phys_sigma_0p17_0p50", Dict("sigma_parallel_min_nm" => 0.17, "sigma_parallel_max_nm" => 0.50)),
        ("phys_sigma_0p18_0p48", Dict("sigma_parallel_min_nm" => 0.18, "sigma_parallel_max_nm" => 0.48)),
        ("phys_sigma_0p20_0p48", Dict("sigma_parallel_min_nm" => 0.20, "sigma_parallel_max_nm" => 0.48)),
        ("phys_spacing_0p40_0p75", Dict("spacing_min_nm" => 0.40, "spacing_max_nm" => 0.75)),
        ("phys_spacing_0p40_0p72", Dict("spacing_min_nm" => 0.40, "spacing_max_nm" => 0.72)),
        ("phys_spacing_0p38_0p72", Dict("spacing_min_nm" => 0.38, "spacing_max_nm" => 0.72)),
        ("phys_pad025_sigma018048", Dict("support_padding_nm" => 0.25, "sigma_parallel_min_nm" => 0.18, "sigma_parallel_max_nm" => 0.48)),
        ("phys_space040075_sigma018048", Dict("spacing_min_nm" => 0.40, "spacing_max_nm" => 0.75, "sigma_parallel_min_nm" => 0.18, "sigma_parallel_max_nm" => 0.48)),
        ("phys_space040072_sigma018048", Dict("spacing_min_nm" => 0.40, "spacing_max_nm" => 0.72, "sigma_parallel_min_nm" => 0.18, "sigma_parallel_max_nm" => 0.48)),
        ("phys_space038072_sigma018048", Dict("spacing_min_nm" => 0.38, "spacing_max_nm" => 0.72, "sigma_parallel_min_nm" => 0.18, "sigma_parallel_max_nm" => 0.48)),
        ("phys_pad025_space040075", Dict("support_padding_nm" => 0.25, "spacing_min_nm" => 0.40, "spacing_max_nm" => 0.75)),
    ]
end

function _merge_model(base_model, overrides)
    model = deepcopy(base_model)
    for (k, v) in overrides
        model[k] = v
    end
    return model
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
        intelligent_sweep=true, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg)
    ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

function _score(r, criterion::String)
    criterion == "gcv" && return r.gcv
    criterion == "aicc" && return r.aicc
    criterion == "cv" && return r.cv_nll_mean
    return r.bic
end

function _best_valid(results, criterion::String)
    valid = [r for r in results if r.success && r.valid && isfinite(_score(r, criterion))]
    isempty(valid) && return nothing
    return sort(valid; by=r -> _score(r, criterion))[1]
end

function _select_effective(results_ell, results_circ, criterion::String)
    valid_ell = [r for r in results_ell if r.success && r.valid && isfinite(_score(r, criterion))]
    valid_circ = [r for r in results_circ if r.success && r.valid && isfinite(_score(r, criterion))]
    isempty(valid_ell) && isempty(valid_circ) && return nothing
    by_ell = Dict(r.n => r for r in valid_ell)
    by_circ = Dict(r.n => r for r in valid_circ)
    best_n, best_s = 0, Inf
    for n in sort(unique(vcat(collect(keys(by_ell)), collect(keys(by_circ)))))
        s = min(haskey(by_ell, n) ? _score(by_ell[n], criterion) : Inf,
                haskey(by_circ, n) ? _score(by_circ[n], criterion) : Inf)
        if s < best_s
            best_n, best_s = n, s
        end
    end
    r_ell = get(by_ell, best_n, nothing)
    r_circ = get(by_circ, best_n, nothing)
    r_ell === nothing && return r_circ
    r_circ === nothing && return r_ell
    return _score(r_ell, criterion) <= _score(r_circ, criterion) ? r_ell : r_circ
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

function main()
    config_file, files, out_tsv, variants_filter, global_maxtime = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    base_model = cfg["model"]
    preproc = cfg["preprocessing"]
    mkpath(dirname(out_tsv))
    header = ["variant", "file", "status", "criterion", "N_ell", "N_circ", "N_eff",
              "score_ell", "score_circ", "score_eff", "bic_ell", "bic_circ",
              "gcv_ell", "gcv_circ", "chi2_ell", "chi2_circ", "support_2d_nm",
              "support_noise_k", "fit_width_nm", "sigma_min", "sigma_max",
              "min_amp_frac", "max_overlap", "kappa_max", "support_padding_nm"]
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
    end
    variants = _variant_specs()
    if !isempty(variants_filter)
        wanted = Set(variants_filter)
        variants = [(name, spec) for (name, spec) in variants if name in wanted]
        isempty(variants) && error("No variants matched --variants=$(join(variants_filter, ','))")
    end
    total = length(variants) * length(files)
    done = 0
    for (vname, overrides) in variants
        model = _merge_model(base_model, overrides)
        if isfinite(global_maxtime)
            model["global_maxtime"] = global_maxtime
        end
        criterion = get(model, "selection_criterion", "gcv")
        for fn in files
            done += 1
            @printf("[%3d/%3d] %-16s %s\n", done, total, vname, fn)
            flush(stdout)
            fp = joinpath(DATA_DIR, fn)
            file_out = joinpath(OUTDIR, vname, splitext(fn)[1])
            mkpath(file_out)
            pcfg, ccfg, ccfg_circ = _configs(model, preproc, file_out)
            pcfg.filepath = fp
            try
                img = GaussianFit2D.read_sxm(fp)
                results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
                results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx_circ)
                best_ell = _best_valid(results_ell, criterion)
                best_circ = _best_valid(results_circ, criterion)
                best_eff = _select_effective(results_ell, results_circ, criterion)
                row = Any[vname, fn, "ok", criterion,
                    best_ell === nothing ? "" : best_ell.n,
                    best_circ === nothing ? "" : best_circ.n,
                    best_eff === nothing ? "" : best_eff.n,
                    best_ell === nothing ? "" : _score(best_ell, criterion),
                    best_circ === nothing ? "" : _score(best_circ, criterion),
                    best_eff === nothing ? "" : _score(best_eff, criterion),
                    best_ell === nothing ? "" : best_ell.bic,
                    best_circ === nothing ? "" : best_circ.bic,
                    best_ell === nothing ? "" : best_ell.gcv,
                    best_circ === nothing ? "" : best_circ.gcv,
                    best_ell === nothing ? "" : best_ell.chi2_reduced,
                    best_circ === nothing ? "" : best_circ.chi2_reduced,
                    _support_length(ctx_circ),
                    model["support_noise_k"], model["fit_width_nm"], model["sigma_parallel_min_nm"],
                    model["sigma_parallel_max_nm"], get(model, "min_amplitude_fraction", 0.3),
                    model["max_overlap"], get(model, "kappa_max", 10.0), model["support_padding_nm"]]
                open(out_tsv, "a") do io
                    println(io, join(row, '\t'))
                end
            catch err
                @warn "variant/file failed" variant=vname file=fn exception=(err, catch_backtrace())
                row = Any[vname, fn, "error", criterion, fill("", length(header) - 4)...]
                open(out_tsv, "a") do io
                    println(io, join(row, '\t'))
                end
            end
        end
    end
    println("\nWrote $out_tsv")
end

main()
