#!/usr/bin/env julia

using GaussianFit2D, TOML, Printf, Statistics

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")
const OUTDIR = "results/param_tuning"

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = ["240817_017.sxm", "240817_019.sxm", "240817_049.sxm", "240817_058.sxm"]
    paddings = [0.10, 0.15, 0.25]
    out_tsv = joinpath(OUTDIR, "support_candidate_audit.tsv")
    global_maxtime = NaN
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            config_file = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--files"
            files = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--files=")
            files = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--paddings"
            paddings = parse.(Float64, split(args[i+1], ",")); i += 2
        elseif startswith(arg, "--paddings=")
            paddings = parse.(Float64, split(split(arg, "=", limit=2)[2], ",")); i += 1
        elseif arg == "--out"
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--global-maxtime"
            global_maxtime = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--global-maxtime=")
            global_maxtime = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, paddings, out_tsv, global_maxtime
end

function _configs(model, preproc, filepath, outdir)
    pcfg = GaussianFit2D.PatternConfig(filepath=filepath, channel="Z", direction="fwd",
        stride=get(preproc, "stride", 1),
        flatten=get(preproc, "flatten", "plane+rows"),
        smooth_radius_px=get(preproc, "smooth_radius_px", 1),
        output_dir=outdir, no_plot=true)
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
        chain_spacing_model=get(model, "chain_spacing_model", "free"),
        chain_tilted_baseline=get(model, "chain_tilted_baseline", true),
        intelligent_sweep=true, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg)
    ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

function _score(r, criterion)
    criterion == "gcv" && return r.gcv
    criterion == "aicc" && return r.aicc
    criterion == "cv" && return r.cv_nll_mean
    return r.bic
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
            n_prefix = ccfg_refine.chain_tilted_baseline ? 3 : 1
            split_idx = n_prefix + n + GaussianFit2D._chain_spacing_param_count(n, ccfg_refine) + n
            p_init = vcat(r_c.params[1:split_idx], r_c.params[(split_idx+1):end], r_c.params[(split_idx+1):end])
            r_ref = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise,
                n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r_ref.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r_ref.params, n, ac_fit, ccfg_refine;
                    amp_min=r_ref.amp_min, amp_range=r_ref.amp_range)
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

function _audit_gcv(r, audit_ctx, fit_axisctx, ccfg)
    pred = GaussianFit2D._chain_model_values(audit_ctx.x, audit_ctx.y, r.params, r.n, fit_axisctx, ccfg;
        amp_min=r.amp_min, amp_range=r.amp_range)
    n_params = GaussianFit2D._chain_nparams(r.n, ccfg)
    gcv, _ = GaussianFit2D._chain_gcv_score(audit_ctx.z, pred, audit_ctx.noise, n_params, ccfg.student_nu)
    return gcv
end

function main()
    config_file, files, paddings, out_tsv, global_maxtime = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    base_model = cfg["model"]
    preproc = cfg["preprocessing"]
    criterion = get(base_model, "selection_criterion", "gcv")
    mkpath(dirname(out_tsv))
    open(out_tsv, "w") do io
        println(io, join(["file", "status", "selected_N", "selected_padding", "selected_audit_gcv",
                          "single_best_N", "single_best_padding", "single_best_score",
                          "n_candidates", "support_audit_nm"], '\t'))
    end
    for fn in files
        @printf("Audit support candidates: %s\n", fn); flush(stdout)
        fp = joinpath(DATA_DIR, fn)
        try
            img = GaussianFit2D.read_sxm(fp)
            # Audit on the largest support candidate, so smaller supports must predict held-out edges.
            audit_pad = maximum(paddings)
            audit_model = deepcopy(base_model); audit_model["support_padding_nm"] = audit_pad
            if isfinite(global_maxtime); audit_model["global_maxtime"] = global_maxtime; end
            audit_pcfg, _audit_ccfg, audit_circ = _configs(audit_model, preproc, fp, joinpath(OUTDIR, "support_audit", "audit", splitext(fn)[1]))
            _, _, audit_ctx = GaussianFit2D.chain_gaussian_sweep(img, audit_pcfg, audit_circ)

            candidates = []
            single_best = nothing
            single_best_pad = NaN
            for pad in paddings
                model = deepcopy(base_model); model["support_padding_nm"] = pad
                if isfinite(global_maxtime); model["global_maxtime"] = global_maxtime; end
                pcfg, ccfg, ccfg_circ = _configs(model, preproc, fp, joinpath(OUTDIR, "support_audit", @sprintf("pad_%0.2f", pad), splitext(fn)[1]))
                results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
                results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx_circ)
                valid = [r for r in results_ell if r.success && r.valid && isfinite(_score(r, criterion))]
                if !isempty(valid)
                    local_best = sort(valid; by=r -> _score(r, criterion))[1]
                    if single_best === nothing || _score(local_best, criterion) < _score(single_best, criterion)
                        single_best = local_best
                        single_best_pad = pad
                    end
                    for r in valid
                        agcv = _audit_gcv(r, audit_ctx, ctx_circ.axisctx, ccfg)
                        push!(candidates, (result=r, padding=pad, audit_gcv=agcv, score=_score(r, criterion), support=ctx_circ.axisctx.tmax-ctx_circ.axisctx.tmin))
                    end
                end
            end
            if isempty(candidates)
                open(out_tsv, "a") do io
                    println(io, join([fn, "no_candidates", "", "", "", "", "", "", 0, audit_ctx.axisctx.tmax-audit_ctx.axisctx.tmin], '\t'))
                end
                continue
            end
            best = sort(candidates; by=x -> x.audit_gcv)[1]
            open(out_tsv, "a") do io
                println(io, join([fn, "ok", best.result.n, best.padding, best.audit_gcv,
                                  single_best === nothing ? "" : single_best.n, single_best_pad,
                                  single_best === nothing ? "" : _score(single_best, criterion),
                                  length(candidates), audit_ctx.axisctx.tmax-audit_ctx.axisctx.tmin], '\t'))
            end
        catch err
            @warn "support audit failed" file=fn exception=(err, catch_backtrace())
            open(out_tsv, "a") do io
                println(io, join([fn, "error", "", "", "", "", "", "", 0, ""], '\t'))
            end
        end
    end
    println("Wrote $out_tsv")
end

main()
