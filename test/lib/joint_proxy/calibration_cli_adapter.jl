_ensure_real_adapter_dependencies_loaded() = nothing

function _real_pattern_config(case::SyntheticCase, opt::CliOptions)
    return GaussianFit2D.PatternConfig(filepath=case.case_id, channel="Z", direction="fwd",
        roi_channel="Z", stride=1, flatten="plane+rows", smooth_radius_px=1,
        output_dir="/tmp/opencode/joint_proxy_calibration", no_plot=true,
        fusion=true, fusion_channels="Z")
end

function _real_chain_config(opt::CliOptions)
    return GaussianFit2D.ChainSweepConfig(n_min=opt.n_min, n_max=opt.n_max, multistart=1,
        intelligent_sweep=false, skip_global=false, fuse_z_bwd=true,
        global_maxtime=1.0, global_maxiter=1500, max_iter=90,
        spacing_min_nm=0.35, spacing_max_nm=0.80, fit_width_nm=0.35,
        support_noise_k=2.5, support_padding_nm=0.40,
        sigma_parallel_min_nm=0.10, sigma_parallel_max_nm=0.35,
        sigma_perp_min_nm=0.10, sigma_perp_max_nm=0.35,
        max_overlap=0.75, kappa_max=10.0, residual_peak_snr_threshold=Inf,
        selection_criterion="gcv")
end

function _real_candidate_report(case::SyntheticCase, opt::CliOptions,
                                bundle::CalibrationBundle;
                                fit_runner=nothing, view_builder=nothing,
                                extractor=nothing, load_deps::Bool=true,
                                pcfg=nothing, ccfg=nothing)
    if load_deps
        _ensure_real_adapter_dependencies_loaded()
        fit_runner === nothing && (fit_runner = GaussianFit2D.chain_gaussian_sweep)
        view_builder === nothing && (view_builder = JointProxyCandidateViews.build_view_data)
        extractor === nothing && (extractor = JointProxyCandidateViews.extract_candidate_views)
        pcfg === nothing && (pcfg = _real_pattern_config(case, opt))
        ccfg === nothing && (ccfg = _real_chain_config(opt))
    end
    fit_runner === nothing && error("fit_runner is required when load_deps=false")
    view_builder === nothing && error("view_builder is required when load_deps=false")
    extractor === nothing && error("extractor is required when load_deps=false")
    results, _best, ctx = redirect_stdout(devnull) do
        fit_runner(case.img, pcfg, ccfg)
    end
    views = JointProxyCandidateViews.ViewData[]
    push!(views, view_builder(case.img, "Z", pcfg; direction="fwd"))
    if any(ch -> lowercase(ch.direction) == "bwd", case.img.channels)
        push!(views, view_builder(case.img, "Z", pcfg; direction="bwd"))
    end
    return extractor(results, ctx, ccfg, views;
        patch_half_nm=bundle.registry.grid_half_nm,
        patch_step_nm=bundle.registry.grid_step_nm)
end

function _real_count_observations(cases::Vector{SyntheticCase}, opt::CliOptions,
                                  bundle::CalibrationBundle; fit_runner=nothing,
                                  view_builder=nothing, extractor=nothing,
                                  load_deps::Bool=true)
    reports = [_real_candidate_report(case, opt, bundle; fit_runner=fit_runner,
        view_builder=view_builder, extractor=extractor, load_deps=load_deps) for case in cases]
    return (count_cases=_count_cases_from_reports(cases, reports, opt, bundle, opt.cases ÷ 2),
        observation_mode="fitted_candidate_views", candidate_reports=reports)
end
