_ensure_real_adapter_dependencies_loaded() = nothing

function _default_pattern_config(file::String, outdir::String)
    return GaussianFit2D.PatternConfig(filepath=file, channel="Z", direction="fwd", output_dir=outdir,
        stride=1, flatten="plane+rows", smooth_radius_px=1, fusion=true, fusion_channels="Z", no_plot=true)
end

function _default_chain_config()
    return GaussianFit2D.ChainSweepConfig(n_min=2, n_max=12, multistart=1, cv_folds=5, cv_method="gcv",
        intelligent_sweep=false, skip_global=false, fuse_z_bwd=true, global_maxtime=1.0, global_maxiter=2000,
        max_iter=90, spacing_min_nm=0.35, spacing_max_nm=0.80, fit_width_nm=0.35, support_noise_k=2.5,
        support_padding_nm=0.40, sigma_parallel_min_nm=0.10, sigma_parallel_max_nm=0.35,
        sigma_perp_min_nm=0.10, sigma_perp_max_nm=0.35, max_overlap=0.75, kappa_max=10.0,
        selection_criterion="gcv")
end

function _calibrate_probability(p::Real, temperature::Real)
    t = Float64(temperature)
    isfinite(t) && t > 0 || return Float64(p)
    q = clamp(Float64(p), eps(Float64), 1.0 - eps(Float64))
    return 1.0 / (1.0 + exp(-log(q / (1.0 - q)) / t))
end

function _row_has_bwd(img)
    any(ch -> lowercase(ch.direction) == "bwd", img.channels)
end

function _build_views(img, pcfg, view_builder)
    views = [view_builder(img, "Z", pcfg; direction="fwd")]
    _row_has_bwd(img) && push!(views, view_builder(img, "Z", pcfg; direction="bwd"))
    return views
end

function _default_count_result(report, calibration::InferenceCalibration)
    try
        return score_count_report(report, calibration.count_temperature; threshold=calibration.count_threshold)
    catch err
        err isa CountCalibrationError || rethrow()
        return CountCasePosterior("", -1, CountPosteriorRow[], missing, 0.0, true, 0, NaN)
    end
end

function _default_type_result(lobes, ensemble, candidate, calibration::Union{Nothing,InferenceCalibration})
    raw = try
        infer_type_posterior(lobes, ensemble; rho=candidate.residual_corr,
            effective_factor=candidate.effective_view_factor)
    catch err
        err isa ArgumentError || rethrow()
        TypePosteriorResult(fill(0.5, length(lobes), 2), fill(0, length(lobes)),
            0.0, TypePosteriorGlobalState[], Pair{String,Float64}[])
    end
    temp = calibration === nothing ? 1.0 : calibration.type_temperature
    rows = NamedTuple[]
    for (i, lobe) in enumerate(candidate.lobes)
        p1 = _calibrate_probability(raw.lobe_marginals[i, 2], temp)
        p0 = 1.0 - p1
        push!(rows, (; p0, p1, raw_p0=raw.lobe_marginals[i, 1], raw_p1=raw.lobe_marginals[i, 2],
            predicted=(calibration !== nothing && max(p0, p1) >= calibration.type_threshold) ? (p1 >= p0 ? 1 : 0) : missing,
            confidence=max(p0, p1)))
    end
    return raw, rows
end

function _candidate_rows(file::String, report, count_result)
    rows = NamedTuple[]
    best_n = count_result.predicted_n === missing ? nothing : count_result.predicted_n
    cand_by_n = Dict(c.n => c for c in report.candidates)
    for r in count_result.rows
        cand = cand_by_n[r.n]
        push!(rows, (; file=file, n=r.n, joint_gcv=r.joint_gcv, delta_rel=r.delta_rel,
            probability=r.probability, rank=r.rank, is_map=best_n === nothing ? false : r.n == best_n,
            view_count=cand.view_count, bwd_missing=cand.bwd_missing, residual_corr=cand.residual_corr,
            effective_view_factor=cand.effective_view_factor, effective_n_views=cand.effective_n_views,
            p_eff=cand.p_eff, nd_joint=cand.nd_joint, source_gcv=cand.source_gcv,
            joint_nrmse=cand.joint_nrmse))
    end
    return rows
end

function _prediction_rows(file::String, candidate, count_result, type_rows, opts::CliOptions, calibration)
    n_prediction = opts.uncalibrated || count_result.abstained ? "?" : count_result.predicted_n
    n_probability = count_result.confidence
    n_confidence = count_result.confidence
    rows = NamedTuple[]
    for (i, lobe) in enumerate(candidate.lobes)
        type_row = type_rows[i]
        pred = opts.uncalibrated || count_result.abstained ? "?" : (ismissing(type_row.predicted) ? "?" : type_row.predicted)
        conf = type_row.confidence
        push!(rows, (; file=file, lobe=lobe.index, N_prediction=n_prediction, N_probability=n_probability,
            N_confidence=n_confidence, predicted=pred, confidence=conf,
            type_probability_0=type_row.p0, type_probability_1=type_row.p1))
    end
    return rows
end

function _summary_row(file::String, report, count_result)
    return (; file=file, N_prediction=count_result.abstained ? "?" : count_result.predicted_n,
        N_probability=count_result.confidence, N_confidence=count_result.confidence,
        N_abstained=count_result.abstained, candidate_count=length(report.candidates),
        lobe_count=sum(length(c.lobes) for c in report.candidates), view_count=report.view_count,
        bwd_missing=report.bwd_missing)
end

function _abstained_summary_row(file::String, report)
    return (; file=file, N_prediction="?", N_probability=0.0, N_confidence=0.0,
        N_abstained=true, candidate_count=0, lobe_count=0,
        view_count=report.view_count, bwd_missing=report.bwd_missing)
end

function _manifest(bundle::InferenceBundle, opts::CliOptions, files::Vector{String}, art_rows)
    return (; provenance=(; config_sha256=bundle.config_hash, source_sha256=bundle.source_hash,
        payload_sha256=bundle.payload_hash, uncalibrated=opts.uncalibrated,
        calibration_path=opts.uncalibrated ? "none" : something(opts.calibration, "none")),
        inputs=(; config=opts.config, data_dir=opts.data_dir, files=files,
            chunk=opts.chunk === nothing ? "none" : "$(opts.chunk[1])/$(opts.chunk[2])"),
        outputs=(; candidate_n_rows=length(art_rows[1]), candidate_lobes_rows=length(art_rows[2]),
            predictions_rows=length(art_rows[3]), chain_summary_rows=length(art_rows[4])))
end

function _real_inference_artifacts(opts::CliOptions, bundle::InferenceBundle, files::Vector{String};
    read_sxm_fn=nothing, fit_runner=nothing, view_builder=nothing, extractor=nothing,
    count_posterior_fn=nothing, type_posterior_fn=nothing, load_deps::Bool=true,
    pcfg=nothing, ccfg=nothing)

    load_deps && _ensure_real_adapter_dependencies_loaded()
    read_sxm_fn === nothing && (read_sxm_fn = STMSXMIO.read_sxm)
    fit_runner === nothing && (fit_runner = GaussianFit2D.chain_gaussian_sweep)
    view_builder === nothing && (view_builder = build_view_data)
    extractor === nothing && (extractor = extract_candidate_views)
    count_posterior_fn === nothing && (count_posterior_fn = _default_count_result)
    type_posterior_fn === nothing && (type_posterior_fn = _default_type_result)

    cn_rows = NamedTuple[]
    cl_rows = NamedTuple[]
    pr_rows = NamedTuple[]
    cs_rows = NamedTuple[]
    for file in files
        path = joinpath(opts.data_dir, file)
        img = read_sxm_fn(path)
        local_pcfg = pcfg === nothing ? _default_pattern_config(path, opts.outdir) : pcfg
        local_ccfg = ccfg === nothing ? _default_chain_config() : ccfg
        results, _best, ctx = fit_runner(img, local_pcfg, local_ccfg)
        views = _build_views(img, local_pcfg, view_builder)
        report = extractor(results, ctx, local_ccfg, views; patch_half_nm=bundle.registry.grid_half_nm,
            patch_step_nm=bundle.registry.grid_step_nm)
        count_result = count_posterior_fn(report, bundle.calibration === nothing ?
            InferenceCalibration(1.0, 1.1, 1.0, 1.1, bundle.config_hash, bundle.source_hash, bundle.payload_hash, "none") : bundle.calibration)
        if !hasproperty(count_result, :rows)
            error("count_posterior_fn must return a result with rows")
        end
        append!(cn_rows, _candidate_rows(file, report, count_result))
        if isempty(count_result.rows)
            push!(cs_rows, _abstained_summary_row(file, report))
            continue
        end
        selected_n = count_result.abstained ? nothing : count_result.predicted_n
        selected_idx = selected_n === nothing ? nothing : findfirst(c -> c.n == selected_n, report.candidates)
        selected_candidate = selected_idx === nothing ? first(report.candidates) : report.candidates[selected_idx]
        _, type_rows = type_posterior_fn([TypePosteriorLobeEvidence(l.residual_patches) for l in selected_candidate.lobes],
            bundle.registry, selected_candidate, bundle.calibration)
        for cand in report.candidates
            _, cand_rows = type_posterior_fn([TypePosteriorLobeEvidence(l.residual_patches) for l in cand.lobes], bundle.registry,
                cand, bundle.calibration)
            for (lobe, row) in zip(cand.lobes, cand_rows)
                pred = row.predicted === missing ? "?" : row.predicted
                push!(cl_rows, (; file=file, n=cand.n, lobe=lobe.index, x_nm=lobe.x_nm, y_nm=lobe.y_nm,
                    t_nm=lobe.t_nm, u_nm=lobe.u_nm, amplitude=lobe.amplitude, sigma_parallel_nm=lobe.sigma_parallel_nm,
                    sigma_perp_nm=lobe.sigma_perp_nm, skew_ratio=lobe.skew_ratio, type_probability_0=row.p0,
                    type_probability_1=row.p1, predicted_type=pred, confidence=row.confidence))
            end
        end
        append!(pr_rows, _prediction_rows(file, selected_candidate, count_result, type_rows, opts, bundle.calibration))
        push!(cs_rows, _summary_row(file, report, count_result))
    end
    artifacts = InferenceArtifacts(vcat(cn_rows...), vcat(cl_rows...), vcat(pr_rows...), cs_rows,
        _manifest(bundle, opts, files, (cn_rows, cl_rows, pr_rows, cs_rows)))
    return artifacts
end
