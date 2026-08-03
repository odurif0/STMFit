function _control_type(control::ControlType)
    control in (CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS) && return :null
    control in (CONTROL_SWAPPED_TYPES, CONTROL_CORRUPTED_VIEW) && return :family_conflict
    return :ok
end

_case_bias(seed::UInt64, n::Int) = 1e-6 * Float64(mod(seed + UInt64(13 * n), UInt64(97)))

_effective_view_factor(rho) = isnan(rho) ? 1.0 : 1.0 / (1.0 + clamp(Float64(rho), 0.0, 0.95))

function _make_views(case::SyntheticCase)
    view_count = case.control == CONTROL_MISSING_BWD ? 1 : 2
    residual_corr = case.control == CONTROL_IDENTICAL_MOLDS ? 0.995 : case.control == CONTROL_CORRUPTED_VIEW ? 0.10 : case.control == CONTROL_MISSING_BWD ? NaN : 0.35
    eff = _effective_view_factor(residual_corr)
    return view_count, residual_corr, eff
end

function _candidate_report(case::SyntheticCase, opt::CliOptions)
    view_count, rho, eff = _make_views(case)
    bwd_missing = view_count == 1
    cands = SyntheticCountCandidate[]
    for n in opt.n_min:opt.n_max
        base = abs(n - case.truth.N)
        ctrl = case.control == CONTROL_NO_MOLECULE ? 0.20 : case.control == CONTROL_CORRUPTED_VIEW ? 0.25 : case.control == CONTROL_IDENTICAL_MOLDS ? 0.15 : 0.0
        score = (base + 1)^2 + ctrl + _case_bias(case.case_seed, n)
        push!(cands, SyntheticCountCandidate(n, true, score))
    end
    push!(cands, SyntheticCountCandidate(opt.n_min - 1, false, Inf))
    return SyntheticCountReport(cands, view_count, bwd_missing, rho)
end

function _type_probs(true_types::Vector{Int}, control::Symbol, family::String)
    n = length(true_types)
    probs = zeros(Float64, n, 2)
    nullish = control == :null
    conflict = control == :family_conflict
    flip = conflict && family == "stm"
    for (i, t) in enumerate(true_types)
        if nullish
            p1 = 0.5
        else
            p1 = family == "geom" ? 0.90 : 0.84
            p1 -= 0.01 * mod(i, 3)
            p1 = clamp(p1, 0.55, 0.95)
        end
        pred = flip ? 1 - t : t
        probs[i, pred + 1] = p1
        probs[i, 2 - pred] = 1.0 - p1
    end
    return probs
end

function _type_result(case::SyntheticCase, control::Symbol, family::String)
    probs = _type_probs(case.truth.sequence, control, family)
    mapseq = [row[2] >= row[1] ? 1 : 0 for row in eachrow(probs)]
    sens = [family => mean(probs[:, 2])]
    evidence = sum(log(max(probs[i, case.truth.sequence[i] + 1], eps(Float64))) for i in eachindex(case.truth.sequence))
    return TypePosteriorResult(probs, mapseq, evidence, TypePosteriorGlobalState[], sens)
end

function _count_cases_from_reports(cases::Vector{SyntheticCase}, reports, opt::CliOptions,
                                   bundle::CalibrationBundle, n_cal::Int)
    out = CountCalibrationCase[]
    for (i, case) in enumerate(cases)
        candidates = hasproperty(reports[i], :candidates) ? getproperty(reports[i], :candidates) : Any[]
        has_true_candidate = any(candidates) do candidate
            hasproperty(candidate, :n) && hasproperty(candidate, :valid) &&
                hasproperty(candidate, :joint_gcv) && candidate.n == case.truth.N &&
                candidate.valid && isfinite(candidate.joint_gcv)
        end
        has_true_candidate || continue
        push!(out, CountCalibrationCase(case.case_id, i <= n_cal ? :calibration : :heldout,
            case.truth.N, reports[i], UInt64(opt.seed), bundle.config_hash))
    end
    return out
end

function _fast_count_observations(cases::Vector{SyntheticCase}, opt::CliOptions,
                                  bundle::CalibrationBundle)
    reports = [_candidate_report(case, opt) for case in cases]
    return (count_cases=_count_cases_from_reports(cases, reports, opt, bundle, opt.cases ÷ 2),
        observation_mode="lightweight_fast", candidate_reports=reports)
end

function _type_cases(cases::Vector{SyntheticCase}, opt::CliOptions, n_cal::Int, bundle::CalibrationBundle)
    out = TypeCalibrationCase[]
    for (i, case) in enumerate(cases)
        split = i <= n_cal ? :calibration : :heldout
        ctl = i <= n_cal ? :ok : _control_type(case.control)
        push!(out, TypeCalibrationCase(case_id=case.case_id, split=split, true_types=copy(case.truth.sequence),
            family_reports=["geom" => _type_result(case, case.control == CONTROL_SWAPPED_TYPES ? :family_conflict : ctl, "geom"),
                            "stm" => _type_result(case, ctl, "stm")],
            seed=UInt64(opt.seed), config_hash=bundle.config_hash, manifest_hash=bundle.payload_hash,
            source_hash=bundle.source_hash, control=ctl))
    end
    return out
end

function _report_total_candidates(reports)
    total = 0
    for r in reports
        hasproperty(r, :candidates) || continue
        total += length(getproperty(r, :candidates))
    end
    return total
end

function _write_section(io::IO, name::String, pairs::Vector{Pair{String,Any}})
    println(io, "[$name]")
    for (k, v) in pairs
        println(io, "$k = $(_fmt(v))")
    end
    println(io)
end

function _write_output(path::String, opt::CliOptions, bundle::CalibrationBundle,
                       cases::Vector{SyntheticCase}, count_model, type_model,
                       count_cases::Vector{CountCalibrationCase}, type_cases::Vector{TypeCalibrationCase},
                       reports, observation_mode::String)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# synthetic-only joint proxy calibration")
        println(io)
        _write_section(io, "provenance", Pair{String,Any}[
            "seed" => opt.seed, "cases" => opt.cases, "fast" => opt.fast,
            "config_sha256" => bundle.config_hash, "source_sha256" => bundle.source_hash,
            "payload_sha256" => bundle.payload_hash, "observation_mode" => observation_mode])
        bundle.code_revision === nothing || _write_section(io, "provenance.extra", Pair{String,Any}["code_revision" => bundle.code_revision])
        _write_section(io, "simulator", Pair{String,Any}[
            "n_min" => opt.n_min, "n_max" => opt.n_max,
            "width" => bundle.sim_cfg.width, "height" => bundle.sim_cfg.height,
            "range_nm" => [bundle.sim_cfg.range_nm[1], bundle.sim_cfg.range_nm[2]],
            "spacing_nm" => bundle.sim_cfg.spacing_nm, "sigma_par_nm" => bundle.sim_cfg.sigma_par_nm,
            "sigma_perp_nm" => bundle.sim_cfg.sigma_perp_nm, "amplitude_nm" => bundle.sim_cfg.amplitude_nm,
            "noise_sigma" => bundle.sim_cfg.noise_sigma, "correlated_noise_frac" => bundle.sim_cfg.correlated_noise_frac,
            "blur_sigma_px" => bundle.sim_cfg.blur_sigma_px, "affine_scale_jitter" => bundle.sim_cfg.affine_scale_jitter,
            "drift_strength" => bundle.sim_cfg.drift_strength, "row_offset_sigma" => bundle.sim_cfg.row_offset_sigma])
        count_cal = count(c -> c.split == :calibration, count_cases)
        type_cal = count(c -> c.split == :calibration, type_cases)
        count_hold = count(c -> c.split == :heldout, count_cases)
        _write_section(io, "splits", Pair{String,Any}[
            "cases" => opt.cases, "count_calibration" => count_cal, "count_heldout" => count_hold,
            "type_calibration" => type_cal, "type_heldout" => opt.cases - type_cal])
        _write_section(io, "count", Pair{String,Any}[
            "temperature" => count_model.temperature, "confidence_threshold" => count_model.confidence_threshold,
            "calibration_nll" => count_model.calibration_nll, "heldout_nll" => count_model.heldout_nll,
            "baseline_temp1_nll" => count_model.baseline_temp1_nll])
        _write_section(io, "count.diagnostics", Pair{String,Any}[
            "temperature_grid" => count_model.diagnostics.temperature_grid,
            "calibration_nll_grid" => count_model.diagnostics.calibration_nll_grid,
            "selected_grid_index" => count_model.diagnostics.selected_grid_index,
            "threshold_source" => count_model.diagnostics.threshold_source,
            "heldout_cases" => count_model.diagnostics.heldout_cases,
            "calibration_cases" => count_model.diagnostics.calibration_cases])
        _write_section(io, "type", Pair{String,Any}[
            "temperature" => type_model.temperature, "confidence_threshold" => type_model.confidence_threshold,
            "calibration_nll" => type_model.calibration_nll, "heldout_nll" => type_model.heldout_nll,
            "baseline_temp1_nll" => type_model.baseline_temp1_nll])
        _write_section(io, "type.diagnostics", Pair{String,Any}[
            "temperature_grid" => type_model.diagnostics.temperature_grid,
            "calibration_nll_grid" => type_model.diagnostics.calibration_nll_grid,
            "selected_grid_index" => type_model.diagnostics.selected_grid_index,
            "threshold_source" => type_model.diagnostics.threshold_source,
            "heldout_cases" => type_model.diagnostics.heldout_cases,
            "calibration_cases" => type_model.diagnostics.calibration_cases])
        null_cases = count(c -> c.control == :null, type_cases)
        conflict_cases = count(c -> c.control == :family_conflict, type_cases)
        _write_section(io, "diagnostics", Pair{String,Any}[
            "candidate_recovered" => count(r -> !r.bwd_missing, reports),
            "candidate_missing_bwd" => count(r -> r.bwd_missing, reports),
            "candidate_report_count" => length(reports),
            "count_candidate_recovered" => length(count_cases),
            "count_candidate_missing_true" => length(reports) - length(count_cases),
            "candidate_observation_total" => _report_total_candidates(reports),
            "null_cases" => null_cases, "family_conflict_cases" => conflict_cases,
            "count_threshold_source" => count_model.diagnostics.threshold_source,
            "type_threshold_source" => type_model.diagnostics.threshold_source])
    end
end

function run_cli(args::AbstractVector{<:AbstractString}=String[]; count_adapter=nothing)
    opt = parse_cli(args)
    bundle = _load_bundle(opt)
    n_cal = opt.cases ÷ 2
    controls = [i <= n_cal ? CONTROL_NORMAL : CONTROL_CYCLE[mod1(i - n_cal, length(CONTROL_CYCLE))] for i in 1:opt.cases]
    cases = generate_batch(MersenneTwister(opt.seed), bundle.sim_cfg, default_proxy_ensemble(); n_cases=opt.cases, controls=controls)
    adapter = count_adapter === nothing ? (opt.fast ? _fast_count_observations : _real_count_observations) : count_adapter
    count_observations = adapter(cases, opt, bundle)
    reports = count_observations.candidate_reports
    count_cases = count_observations.count_cases
    type_cases = _type_cases(cases, opt, n_cal, bundle)
    count_model = calibrate_count_temperature(count_cases; seed=opt.seed, config_hash=bundle.config_hash)
    type_model = calibrate_type_temperature(type_cases; seed=opt.seed, config_hash=bundle.config_hash,
        manifest_hash=bundle.payload_hash, source_hash=bundle.source_hash)
    _write_output(opt.out, opt, bundle, cases, count_model, type_model, count_cases, type_cases, reports,
        count_observations.observation_mode)
    return opt.out
end
