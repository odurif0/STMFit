# OLS / joint scoring helpers for Todo 2.

"""
    extract_candidate_views(results, ctx, ccfg, views; kwargs...) -> CandidateViewReport

Retain every `success && valid && isfinite(gcv)` candidate from `results`,
decode each with its own N, and build paired-view joint metrics + aligned
residual patches.
"""
function extract_candidate_views(
        results::AbstractVector{ChainModelResult},
        ctx::NamedTuple,
        ccfg::ChainSweepConfig,
        views::AbstractVector{ViewData};
        patch_half_nm::Float64=0.32,
        patch_step_nm::Float64=0.08,
        min_pixels::Int=20)::CandidateViewReport

    present_views = [v for v in views if v.present]
    isempty(present_views) && throw(ArgumentError(
        "no present views: at least one ViewData with present=true is required"))

    view_count = length(present_views)
    bwd_missing = !any(v -> lowercase(v.label) == "bwd", present_views)

    xs = ctx.xs::Vector{Float64}
    ys = ctx.ys::Vector{Float64}
    nx, ny = length(xs), length(ys)

    # Validate each present view shares the sweep grid (mask alignment contract).
    for v in present_views
        if length(v.xs) != nx || length(v.ys) != ny
            throw(ArgumentError(
                "view '$(v.label)' grid ($(length(v.xs))×$(length(v.ys))) " *
                "does not match sweep grid ($nx×$ny); cannot align mask"))
        end
    end

    mask = _resolve_mask(ctx, nx, ny, min_pixels)
    mask_inds = findall(mask)

    patch_tu = _patch_offset_grid(patch_half_nm, patch_step_nm)
    ax, ay = ctx.axisctx.axis
    px, py = ctx.axisctx.perp

    xmat = repeat(reshape(xs, 1, :), ny, 1)
    ymat = repeat(reshape(ys, :, 1), 1, nx)
    vec_x = vec(xmat)
    vec_y = vec(ymat)

    candidates = CandidateNView[]
    skipped = SkippedCandidate[]
    audit = String[]
    push!(audit, "view_count=$view_count bwd_missing=$bwd_missing " *
                 "mask_pixels=$(length(mask_inds)) patch_grid=$(length(patch_tu))")

    for r in results
        if !r.success
            push!(skipped, SkippedCandidate(r.n, "not_success: $(r.reason)"))
            continue
        end
        if !r.valid
            push!(skipped, SkippedCandidate(r.n, "not_valid: $(r.reason)"))
            continue
        end
        if !isfinite(r.gcv)
            push!(skipped, SkippedCandidate(r.n, "non_finite_gcv"))
            continue
        end

        model_vec = GaussianFit2D._chain_model_values(
            vec_x, vec_y, r.params, r.n, ctx.axisctx, ccfg;
            amp_min=r.amp_min, amp_range=r.amp_range)
        model_img = reshape(model_vec, ny, nx)

        good = trues(ny, nx)
        @inbounds for idx in mask_inds
            m = model_img[idx]
            ok = isfinite(m)
            for v in present_views
                ok = ok && isfinite(v.z[idx])
            end
            good[idx] = ok
        end
        n_good = count(good)
        if n_good < min_pixels
            push!(skipped, SkippedCandidate(r.n,
                "insufficient_finite_pixels: $n_good < $min_pixels"))
            continue
        end
        good_inds = findall(good)
        model_masked = [model_img[i] for i in good_inds]

        recalibrations = ViewRecalibration[]
        resid_imgs = Dict{String,Matrix{Float64}}()
        total_rss = 0.0
        var_total = 0.0
        for v in present_views
            z_masked = [v.z[i] for i in good_inds]
            rec = _ols_recalibrate(v.label, z_masked, model_masked)
            push!(recalibrations, rec)
            if rec.valid
                total_rss += rec.rss
                var_total += _var_skipnan(z_masked)
                resid_imgs[v.label] = v.z .- (rec.a .* model_img .+ rec.b)
            else
                resid_imgs[v.label] = fill(NaN, ny, nx)
            end
        end

        p_chain = GaussianFit2D._chain_nparams(r.n, ccfg)
        p_eff = p_chain + 2 * view_count
        nd_joint = view_count * n_good
        joint_gcv = nd_joint > p_eff ? nd_joint / (nd_joint - p_eff)^2 * total_rss : Inf
        joint_nrmse = (var_total > eps(Float64) && isfinite(total_rss)) ?
            sqrt(total_rss / (view_count * n_good)) / sqrt(var_total) : NaN

        residual_corr = if view_count >= 2 &&
                          recalibrations[1].valid && recalibrations[2].valid
            _pearson(recalibrations[1].resid, recalibrations[2].resid)
        else
            NaN
        end

        eff_factor = effective_view_factor(residual_corr)
        eff_n_views = view_count * eff_factor

        lobes = _decode_and_patch(r, ccfg, ctx.axisctx, patch_tu,
                                  ax, ay, px, py, present_views, resid_imgs,
                                  model_img)

        push!(candidates, CandidateNView(
            r.n, r.gcv, r.success, r.valid, r.reason,
            view_count, bwd_missing, recalibrations,
            joint_gcv, joint_nrmse, residual_corr,
            eff_factor, eff_n_views, p_eff, nd_joint,
            lobes, patch_tu))
        push!(audit, "retained N=$(r.n) gcv=$(r.gcv) joint_gcv=$joint_gcv " *
                     "rho=$residual_corr eff_views=$eff_n_views lobes=$(length(lobes))")
    end

    return CandidateViewReport(candidates, skipped, view_count, bwd_missing, audit)
end

function _resolve_mask(ctx::NamedTuple, nx::Int, ny::Int, min_pixels::Int)
    hasproperty(ctx, :mask) || return trues(ny, nx)
    mask = ctx.mask
    if mask === nothing || isempty(mask)
        return trues(ny, nx)
    end
    if size(mask) != (ny, nx)
        throw(ArgumentError(
            "mask size $(size(mask)) != grid ($ny, $nx); refusing to emit " *
            "misaligned rows — pass a mask matching the sweep grid or nothing"))
    end
    if !(eltype(mask) <: Bool)
        throw(ArgumentError(
            "mask eltype $(eltype(mask)) is not Bool; non-finite/typed masks " *
            "are rejected explicitly"))
    end
    n_true = count(mask)
    if n_true < min_pixels
        throw(ArgumentError(
            "undersized mask: $n_true true pixels < $min_pixels minimum"))
    end
    return mask
end

function _ols_recalibrate(label::String, z_scan::Vector{Float64}, model::Vector{Float64})
    n = length(z_scan)
    n < 4 && return ViewRecalibration(label, NaN, NaN, NaN, Float64[], n, false,
                                      "insufficient_pixels: $n < 4")
    A = hcat(model, ones(n))
    try
        coeff = A \ z_scan
        a, b = coeff[1], coeff[2]
        if !isfinite(a) || !isfinite(b)
            return ViewRecalibration(label, NaN, NaN, NaN, Float64[], n, false,
                                     "non_finite_ols")
        end
        pred = a .* model .+ b
        resid = z_scan .- pred
        rss = sum(abs2, resid)
        return ViewRecalibration(label, a, b, rss, resid, n, true, "ok")
    catch
        return ViewRecalibration(label, NaN, NaN, NaN, Float64[], n, false,
                                 "ols_singular")
    end
end

function _pearson(a::Vector{Float64}, b::Vector{Float64})
    n = length(a)
    n < 2 && return NaN
    af = a .- mean(a)
    bf = b .- mean(b)
    denom = norm(af) * norm(bf)
    denom <= eps(Float64) && return NaN
    return dot(af, bf) / denom
end

_var_skipnan(v::Vector{Float64}) = (g = filter(isfinite, v); isempty(g) ? 0.0 : var(g))
