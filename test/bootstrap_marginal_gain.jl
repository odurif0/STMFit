#!/usr/bin/env julia

using GaussianFit2D, Random, Printf, TOML, Statistics

const DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "/home/durif/Rebecca/data/data/20240817_LHe_Cu100")

function _parse_cli(args)
    config_file = "config/chitosan.toml"
    files = ["240817_017.sxm", "240817_019.sxm", "240817_058.sxm"]
    out_tsv = "results/chitosan_case_audit/bootstrap_marginal_gain.tsv"
    n0 = 6
    n1 = 7
    reps = 30
    block_nm = 0.45
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
        elseif arg == "--out"
            out_tsv = args[i + 1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--n0"
            n0 = parse(Int, args[i + 1]); i += 2
        elseif startswith(arg, "--n0=")
            n0 = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--n1"
            n1 = parse(Int, args[i + 1]); i += 2
        elseif startswith(arg, "--n1=")
            n1 = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--reps"
            reps = parse(Int, args[i + 1]); i += 2
        elseif startswith(arg, "--reps=")
            reps = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--block-nm"
            block_nm = parse(Float64, args[i + 1]); i += 2
        elseif startswith(arg, "--block-nm=")
            block_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        else
            error("Unknown option: $arg")
        end
    end
    return config_file, files, out_tsv, n0, n1, reps, block_nm
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
        chain_tilted_baseline=get(model, "chain_tilted_baseline", true),
        intelligent_sweep=false, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg); ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    refined = Dict{Int,GaussianFit2D.ChainModelResult}()
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ctx_circ.axisctx_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    ccfg_refine = deepcopy(ccfg_ell); ccfg_refine.skip_global = true; ccfg_refine.max_iter = 50; ccfg_refine.multistart = 1
    for r_c in results_circ
        r_c.success || continue
        n = r_c.n
        n_prefix = 1 + (ccfg_refine.chain_tilted_baseline ? 2 : 0)
        split_idx = n_prefix + n + 1 + (n - 1) + n
        p_init = vcat(r_c.params[1:split_idx], r_c.params[(split_idx+1):end], r_c.params[(split_idx+1):end])
        try
            r = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r.params, n, ac_fit, ccfg_refine; amp_min=r.amp_min, amp_range=r.amp_range)
                GaussianFit2D._finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, ac_fit, ccfg_refine)
                refined[n] = r
            end
        catch err
            @warn "refinement failed" n exception=(err, catch_backtrace())
        end
    end
    return refined, (xs=xs, ys=ys, zimg=zimg, xfit=xfit, yfit=yfit, zfit=zfit, noise=noise, axisctx=ac_fit)
end

function _block_wild_residual(resid, t, block_nm, rng)
    out = similar(resid)
    order = sortperm(t)
    t_sorted = t[order]
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && t_sorted[j + 1] - t_sorted[i] <= block_nm
            j += 1
        end
        sign = rand(rng, Bool) ? 1.0 : -1.0
        out[order[i:j]] .= sign .* resid[order[i:j]]
        i = j + 1
    end
    return out
end

function main()
    config_file, files, out_tsv, n0, n1, reps, block_nm = _parse_cli(ARGS)
    cfg = TOML.parsefile(config_file)
    model, preproc = cfg["model"], cfg["preprocessing"]
    mkpath(dirname(out_tsv))
    open(out_tsv, "w") do io
        println(io, join(["file", "n0", "n1", "reps", "block_nm", "rss_n0", "rss_n1", "T_obs", "p_boot", "T_boot_median", "T_boot_q95"], '\t'))
    end
    rng = MersenneTwister(4321)
    for fn in files
        @printf("Bootstrap marginal gain %s N=%d→%d B=%d\n", fn, n0, n1, reps); flush(stdout)
        pcfg, ccfg, ccfg_circ = _configs(model, preproc, joinpath("results/chitosan_case_audit/bootstrap", splitext(fn)[1]))
        pcfg.filepath = joinpath(DATA_DIR, fn)
        img = GaussianFit2D.read_sxm(pcfg.filepath)
        results_circ, _, ctx_circ = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg_circ)
        refined, dat = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx_circ)
        haskey(refined, n0) && haskey(refined, n1) || error("Missing fitted N=$n0 or N=$n1 for $fn")
        r0, r1 = refined[n0], refined[n1]
        pred0 = GaussianFit2D._chain_model_values(dat.xfit, dat.yfit, r0.params, n0, dat.axisctx, ccfg; amp_min=r0.amp_min, amp_range=r0.amp_range)
        rss0 = sum(abs2, dat.zfit .- pred0)
        pred1 = GaussianFit2D._chain_model_values(dat.xfit, dat.yfit, r1.params, n1, dat.axisctx, ccfg; amp_min=r1.amp_min, amp_range=r1.amp_range)
        rss1 = sum(abs2, dat.zfit .- pred1)
        Tobs = max(rss0 - rss1, 0.0)
        t, _u = GaussianFit2D._chain_coordinates(dat.xfit, dat.yfit, dat.axisctx)
        resid0 = dat.zfit .- pred0
        ccfg_boot = deepcopy(ccfg); ccfg_boot.skip_global = true; ccfg_boot.max_iter = 40; ccfg_boot.multistart = 1
        Tboot = Float64[]
        for _ in 1:reps
            zstar = pred0 .+ _block_wild_residual(resid0, t, block_nm, rng)
            b0 = GaussianFit2D._fit_chain_n(dat.xs, dat.ys, dat.zimg, dat.xfit, dat.yfit, zstar, dat.noise, n0, dat.axisctx, ccfg_boot; starts=1, warm_start=r0.params)
            b1 = GaussianFit2D._fit_chain_n(dat.xs, dat.ys, dat.zimg, dat.xfit, dat.yfit, zstar, dat.noise, n1, dat.axisctx, ccfg_boot; starts=1, warm_start=r1.params)
            b0.success && b1.success || continue
            p0 = GaussianFit2D._chain_model_values(dat.xfit, dat.yfit, b0.params, n0, dat.axisctx, ccfg_boot; amp_min=b0.amp_min, amp_range=b0.amp_range)
            p1 = GaussianFit2D._chain_model_values(dat.xfit, dat.yfit, b1.params, n1, dat.axisctx, ccfg_boot; amp_min=b1.amp_min, amp_range=b1.amp_range)
            push!(Tboot, max(sum(abs2, zstar .- p0) - sum(abs2, zstar .- p1), 0.0))
        end
        pboot = isempty(Tboot) ? NaN : count(>=(Tobs), Tboot) / length(Tboot)
        q95 = isempty(Tboot) ? NaN : quantile(Tboot, 0.95)
        med = isempty(Tboot) ? NaN : median(Tboot)
        open(out_tsv, "a") do io
            println(io, join([fn, n0, n1, length(Tboot), block_nm, rss0, rss1, Tobs, pboot, med, q95], '\t'))
        end
    end
    println("Wrote $out_tsv")
end

main()
