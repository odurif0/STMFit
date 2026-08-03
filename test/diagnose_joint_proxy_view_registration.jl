#!/usr/bin/env julia

include(joinpath(@__DIR__, "diagnose_joint_proxy_view_row_drift.jl"))

const REGISTRATION_FIELDS = [:file, :stage, :ablation, :estimated_lag_px,
    :estimated_lag_nm, :n, :lobe, :n_pairs, :corr, :slope, :affine_nrmse,
    :fwd_p1, :bwd_p1, :fwd_argmax, :bwd_argmax, :argmax_concordant,
    :fwd_confidence, :bwd_confidence]

parse_registration_options(args) = parse_row_drift_options(args)

function shift_view_x(view::Main.JointProxyCandidateViews.ViewData, lag_px::Int)
    shifted = fill(NaN, size(view.z))
    for ix in axes(view.z, 2)
        source = ix + lag_px
        firstindex(view.z, 2) <= source <= lastindex(view.z, 2) || continue
        shifted[:, ix] .= view.z[:, source]
    end
    return Main.JointProxyCandidateViews.ViewData(view.label, copy(view.xs), copy(view.ys),
        shifted, view.noise, view.present)
end

function robust_global_lag(rows)
    lags = [Int(row.best_lag_px) for row in rows if row.status == "ok"]
    isempty(lags) && error("no finite row lag is available for registration")
    return round(Int, median(lags))
end

_argmax_type(p0::Real, p1::Real) = p1 > p0 ? 1 : p0 > p1 ? 0 : -1

function _registered_views(views, lag_px::Int)
    return [view.label == "bwd" ? shift_view_x(view, lag_px) : view for view in views]
end

function _registration_stage_rows(file::String, stage::String, lag_px::Int, lag_nm::Float64,
    candidate, evidence, registry::ProxyEnsemble)
    ensembles = Dict(family_ensembles(registry))
    rows = NamedTuple[]
    for ablation in ("stm_dft_v1", "combined")
        ensemble = ensembles[ablation]
        fwd_patches = [Dict("fwd" => lobe.view_patches["fwd"]) for lobe in evidence]
        bwd_patches = [Dict("bwd" => lobe.view_patches["bwd"]) for lobe in evidence]
        fwd_result = generative_type_posterior(candidate, fwd_patches, ensemble;
            effective_factor=1.0)
        bwd_result = generative_type_posterior(candidate, bwd_patches, ensemble;
            effective_factor=1.0)
        fwd_result === nothing && error("non-finite fwd generative evidence after registration")
        bwd_result === nothing && error("non-finite bwd generative evidence after registration")
        for (i, (lobe, lobe_evidence)) in enumerate(zip(candidate.lobes, evidence))
            transfer = _affine_patch_metrics(lobe_evidence.view_patches["fwd"],
                lobe_evidence.view_patches["bwd"])
            fp0, fp1 = fwd_result.lobe_marginals[i, 1], fwd_result.lobe_marginals[i, 2]
            bp0, bp1 = bwd_result.lobe_marginals[i, 1], bwd_result.lobe_marginals[i, 2]
            farg, barg = _argmax_type(fp0, fp1), _argmax_type(bp0, bp1)
            push!(rows, (; file, stage, ablation, estimated_lag_px=lag_px,
                estimated_lag_nm=lag_nm, n=candidate.n, lobe=lobe.index,
                transfer.n_pairs, transfer.corr, transfer.slope, transfer.affine_nrmse,
                fwd_p1=fp1, bwd_p1=bp1, fwd_argmax=farg, bwd_argmax=barg,
                argmax_concordant=farg == barg,
                fwd_confidence=max(fp0, fp1), bwd_confidence=max(bp0, bp1)))
        end
    end
    return rows
end

function run_registration_diagnostic(options::RowDriftOptions)
    case = fixed_geometry_case(options.config, options.data_dir, options.file,
        dirname(options.out), options.seed)
    lag_rows = compute_row_drift_rows(case, options.file, options.flatten_mode,
        options.max_lag_px, options.min_pairs)
    lag_px = robust_global_lag(lag_rows)

    pattern = pattern_with_flatten(case.base_pattern, options.flatten_mode)
    views = Inference._build_views(case.image, pattern,
        Main.JointProxyCandidateViews.build_view_data)
    has_bwd = any(view -> view.present && view.label == "bwd", views)
    has_bwd || error("registration requires a bwd view")
    x_step_nm = length(first(views).xs) > 1 ? median(diff(first(views).xs)) : NaN
    lag_nm = lag_px * x_step_nm

    before_evidence = image_lobe_evidence(case.candidate, views)
    registered = _registered_views(views, lag_px)
    after_evidence = image_lobe_evidence(case.candidate, registered)
    rows = vcat(
        _registration_stage_rows(options.file, "before", lag_px, lag_nm,
            case.candidate, before_evidence, case.registry),
        _registration_stage_rows(options.file, "registered", lag_px, lag_nm,
            case.candidate, after_evidence, case.registry),
    )
    write_tsv(options.out, rows; fields=REGISTRATION_FIELDS)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_registration_diagnostic(parse_registration_options(ARGS))
end
