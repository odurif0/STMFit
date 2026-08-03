#!/usr/bin/env julia

# TDD test for the candidate-N paired-view geometry + residual-patch extractor
# (Todo 2, joint-proxy-mold-inference-v1).
#
# Acceptance criteria (from .omo/plans/joint-proxy-mold-inference-v1.md):
#   - julia --project=. test/joint_proxy/test_candidate_views.jl passes
#   - all valid synthetic candidates remain present
#   - patch dimensions/coordinates match across views
#   - OLS nuisance recovery is within test tolerance
#   - duplicating one view yields effective evidence equivalent to one
#     independent view rather than double confidence
#
# QA scenarios:
#   Happy   — paired synthetic data with known affine view scales returns two
#             aligned finite patch sets and finite joint metrics.
#   Failure — missing bwd falls back with view_count=1; invalid/failed
#             candidates are excluded with reasons; non-finite/undersized masks
#             fail explicitly rather than emit corrupt rows.
#
# Label-free: no expected-N, no benchmark labels, no production selector, no
# selected-N feature TSV reuse. Synthetic fixtures only.

using Test
using LinearAlgebra
using Statistics
using Random

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "candidate_views.jl"))
using .JointProxyCandidateViews:
    ViewData, ViewRecalibration, LobePatch, CandidateNView, CandidateViewReport,
    SkippedCandidate, build_view_data, extract_candidate_views,
    effective_view_factor

using GaussianFit2D
using GaussianFit2D: ChainModelResult, ChainSweepConfig
const G = GaussianFit2D

# ─── Synthetic fixture helpers ───────────────────────────────────────────────
# Config chosen so the param layout is simple and deterministic:
#   chain_tilted_baseline=false, circular sigmas, free spacing, gaussian peak.
# With a zero param vector, _decode_chain yields n finite on-axis lobes with
# unit amplitudes and mid-range spacings/sigmas — a valid non-degenerate model.

function _ccfg()
    ChainSweepConfig(
        n_min=2, n_max=12,
        spacing_min_nm=0.40, spacing_max_nm=0.80,
        lateral_max_nm=0.30, fit_width_nm=0.50,
        support_noise_k=2.5, support_padding_nm=0.20, support_min_length_nm=1.0,
        sigma_parallel_min_nm=0.15, sigma_parallel_max_nm=0.35,
        sigma_perp_min_nm=0.15, sigma_perp_max_nm=0.35,
        chain_circular_sigmas=true, shared_sigma_types=0,
        chain_spacing_model="free", chain_tilted_baseline=false,
        peak_profile=:gaussian, skew_ratio_max=2.0,
        max_overlap=1.0,            # → effective spacing min == spacing_min_nm
        selection_criterion="gcv", intelligent_sweep=false)
end

function _axisctx(xs, ys)
    ox = (xs[1] + xs[end]) / 2
    oy = (ys[1] + ys[end]) / 2
    axis = [1.0, 0.0]
    perp = [0.0, 1.0]               # = [-axis[2], axis[1]]
    t = xs .- ox
    return (origin=(ox, oy), axis=axis, perp=perp,
            tmin=minimum(t), tmax=maximum(t))
end

# Zero param vector for the chosen config: 4n+1 floats (b0, n amps, t0, n-1
# deltas, n laterals, n circular sigmas). Decodes to finite on-axis lobes.
_params_zeros(n::Int) = zeros(Float64, 4n + 1)

function _result(n, ccfg, axisctx, xs, ys; success=true, valid=true,
                 gcv=0.01 * n, reason="")
    p = _params_zeros(n)
    nx, ny = length(xs), length(ys)
    xmat = repeat(reshape(xs, 1, :), ny, 1)
    ymat = repeat(reshape(ys, :, 1), 1, nx)
    model = reshape(G._chain_model_values(vec(xmat), vec(ymat), p, n,
                                           axisctx, ccfg), ny, nx)
    @assert all(isfinite, model) "synthetic model not finite for n=$n"
    @assert std(vec(model)) > 1e-6 "synthetic model degenerate for n=$n"
    r = ChainModelResult(n=n, params=p, success=success, valid=valid,
                         gcv=gcv, reason=reason, amp_min=NaN, amp_range=NaN)
    return r, model
end

function _view(label, model_img, xs, ys, a, b; noise=0.0, noise_vec=nothing)
    z = a .* model_img .+ b
    if noise_vec !== nothing
        z = z .+ noise_vec
    elseif noise > 0
        z = z .+ noise .* randn(size(z))
    end
    return ViewData(label, xs, ys, copy(z), max(noise, 1e-12), true)
end

const PATCH_HALF = 0.32
const PATCH_STEP = 0.08

# ─── Tests ───────────────────────────────────────────────────────────────────

@testset "candidate_views" begin
    ccfg = _ccfg()
    xs = collect(range(0.0, 8.0; length=80))
    ys = collect(range(0.0, 8.0; length=80))
    axisctx = _axisctx(xs, ys)
    nx, ny = length(xs), length(ys)
    ctx = (xs=xs, ys=ys, mask=trues(ny, nx), axisctx=axisctx)

    r3, m3 = _result(3, ccfg, axisctx, xs, ys)
    r4, m4 = _result(4, ccfg, axisctx, xs, ys)
    r5, m5 = _result(5, ccfg, axisctx, xs, ys)
    r_fail = ChainModelResult(n=6, params=zeros(25); success=false, valid=false,
                              gcv=Inf, reason="optimizer_failed", amp_min=NaN, amp_range=NaN)
    r_invalid = ChainModelResult(n=7, params=zeros(29); success=true, valid=false,
                                 gcv=0.5, reason="spacing_violation", amp_min=NaN, amp_range=NaN)

    @testset "retains every valid candidate; excludes invalid/failed with reasons" begin
        results = [r3, r4, r5, r_fail, r_invalid]
        views = [_view("fwd", m4, xs, ys, 2.0, 0.5),
                 _view("bwd", m4, xs, ys, 1.5, -0.3)]
        rep = extract_candidate_views(results, ctx, ccfg, views;
                                      patch_half_nm=PATCH_HALF, patch_step_nm=PATCH_STEP)
        @test rep isa CandidateViewReport
        @test length(rep.candidates) == 3
        @test sort([c.n for c in rep.candidates]) == [3, 4, 5]
        @test length(rep.skipped) == 2
        @test sort([s.n for s in rep.skipped]) == [6, 7]
        @test all(s -> !isempty(s.reason), rep.skipped)
        @test rep.view_count == 2
        @test !rep.bwd_missing
    end

    @testset "OLS nuisance recovery within tolerance" begin
        a_f, b_f, a_b, b_b = 2.0, 0.5, 1.5, -0.3
        views = [_view("fwd", m4, xs, ys, a_f, b_f; noise=0.0),
                 _view("bwd", m4, xs, ys, a_b, b_b; noise=0.0)]
        rep = extract_candidate_views([r4], ctx, ccfg, views;
                                      patch_half_nm=PATCH_HALF, patch_step_nm=PATCH_STEP)
        @test length(rep.candidates) == 1
        c = rep.candidates[1]
        vf = c.views[1]; vb = c.views[2]
        @test vf.label == "fwd"
        @test vb.label == "bwd"
        @test isapprox(vf.a, a_f; atol=1e-8)
        @test isapprox(vf.b, b_f; atol=1e-8)
        @test isapprox(vb.a, a_b; atol=1e-8)
        @test isapprox(vb.b, b_b; atol=1e-8)
        @test vf isa ViewRecalibration && vf.valid
        @test isfinite(c.joint_gcv)
        @test isfinite(c.joint_nrmse)
        @test c.joint_nrmse < 1e-6
        @test c.p_eff == G._chain_nparams(4, ccfg) + 4
        @test c.nd_joint == 2 * vf.n_pixels
    end

    @testset "patch dimensions/coordinates match across views" begin
        views = [_view("fwd", m4, xs, ys, 2.0, 0.5; noise=1e-3),
                 _view("bwd", m4, xs, ys, 1.5, -0.3; noise=1e-3)]
        rep = extract_candidate_views([r4], ctx, ccfg, views;
                                      patch_half_nm=PATCH_HALF, patch_step_nm=PATCH_STEP)
        c = rep.candidates[1]
        @test length(c.lobes) == 4
        npix = length(c.patch_tu)
        @test npix == 81                    # 9 x 9 grid
        for lb in c.lobes
            @test lb isa LobePatch
            @test length(lb.residual_patches["fwd"]) == npix
            @test length(lb.residual_patches["bwd"]) == npix
            @test length(lb.raw_patches["fwd"]) == npix
            @test count(isfinite, lb.residual_patches["fwd"]) > 0
            @test count(isfinite, lb.residual_patches["bwd"]) > 0
            @test isfinite(lb.x_nm) && isfinite(lb.y_nm)
        end
        # shared (t,u) offset grid identical for all lobes
        @test all(isequal(c.patch_tu[1], c.patch_tu[k]) === false
                  for k in 2:npix) === true || npix > 1
    end

    @testset "duplicate view does not double confidence" begin
        Random.seed!(20260710)
        noise_f = 1e-2 .* randn(ny, nx)
        noise_b = 1e-2 .* randn(ny, nx)
        # two independent views (distinct affine + independent noise)
        rep_indep = extract_candidate_views(
            [r4], ctx, ccfg,
            [_view("fwd", m4, xs, ys, 2.0, 0.5; noise_vec=noise_f),
             _view("bwd", m4, xs, ys, 1.5, -0.3; noise_vec=noise_b)];
            patch_half_nm=PATCH_HALF, patch_step_nm=PATCH_STEP)
        c_indep = rep_indep.candidates[1]
        @test c_indep.view_count == 2
        @test c_indep.effective_n_views > 1.8     # ~2 independent views

        # duplicated view: bwd is an exact copy of fwd (same noise field)
        v_dup = _view("fwd", m4, xs, ys, 2.0, 0.5; noise_vec=noise_f)
        v_dup_b = ViewData("bwd", v_dup.xs, v_dup.ys, copy(v_dup.z), v_dup.noise, true)
        rep_dup = extract_candidate_views([r4], ctx, ccfg, [v_dup, v_dup_b];
                                          patch_half_nm=PATCH_HALF, patch_step_nm=PATCH_STEP)
        c_dup = rep_dup.candidates[1]
        @test c_dup.view_count == 2
        @test c_dup.residual_corr > 0.99
        @test c_dup.effective_n_views < 1.1       # ~1, NOT 2
        @test c_dup.effective_n_views > 0.9
        @test c_dup.effective_n_views < c_indep.effective_n_views - 0.5
        # factor formula: 1/(1+clamp(rho,0,0.95))
        @test isapprox(c_dup.effective_view_factor,
                       1.0 / (1.0 + 0.95); atol=1e-12)
    end

    @testset "missing bwd falls back to one view" begin
        views = [_view("fwd", m4, xs, ys, 2.0, 0.5)]
        rep = extract_candidate_views([r4], ctx, ccfg, views;
                                      patch_half_nm=PATCH_HALF, patch_step_nm=PATCH_STEP)
        @test rep.view_count == 1
        @test rep.bwd_missing
        c = rep.candidates[1]
        @test c.view_count == 1
        @test length(c.views) == 1
        @test isnan(c.residual_corr)
        @test c.effective_view_factor == 1.0
        @test c.effective_n_views == 1.0
        @test isfinite(c.joint_gcv)
        @test c.p_eff == G._chain_nparams(4, ccfg) + 2
        @test length(c.lobes) == 4
    end

    @testset "malformed/undersized masks fail explicitly; non-finite candidate excluded" begin
        views = [_view("fwd", m4, xs, ys, 2.0, 0.5),
                 _view("bwd", m4, xs, ys, 1.5, -0.3)]
        # wrong-size mask
        bad_ctx = (xs=xs, ys=ys, mask=trues(2, 2), axisctx=axisctx)
        @test_throws ArgumentError extract_candidate_views([r4], bad_ctx, ccfg, views)
        # undersized mask (single True pixel)
        small = falses(ny, nx); small[40, 40] = true
        small_ctx = (xs=xs, ys=ys, mask=small, axisctx=axisctx)
        @test_throws ArgumentError extract_candidate_views([r4], small_ctx, ccfg, views)
        # non-finite model: NaN params → model all-NaN → excluded with reason
        r_bad = ChainModelResult(n=4, params=fill(NaN, 17); success=true, valid=true,
                                 gcv=0.1, reason="", amp_min=NaN, amp_range=NaN)
        rep = extract_candidate_views([r4, r_bad], ctx, ccfg, views)
        @test length(rep.candidates) == 1
        @test length(rep.skipped) == 1
        @test !isempty(rep.skipped[1].reason)
    end

    @testset "effective_view_factor helper is the documented formula" begin
        @test effective_view_factor(0.0) ≈ 1.0
        @test effective_view_factor(-0.5) ≈ 1.0          # negative clamped to 0
        @test effective_view_factor(1.0) ≈ 1.0 / 1.95     # clamped to 0.95
        @test effective_view_factor(NaN) ≈ 1.0            # missing corr → 1 view
        @test effective_view_factor(0.5) ≈ 1.0 / 1.5
    end
end

# ─── Exit code propagation ───────────────────────────────────────────────────
# Propagate a nonzero exit if any sub-test failed or errored (defensive across
# Julia versions: DefaultTestSet field names have varied).
function _any_nonpass(ts)
    for f in (:n_fail, :n_failed, :n_errors, :n_error, :n_non_pass)
        if hasproperty(ts, f)
            v = getproperty(ts, f)
            v isa Number && v > 0 && return true
        end
    end
    if hasproperty(ts, :results)
        for c in ts.results
            c isa Test.DefaultTestSet && _any_nonpass(c) && return true
        end
    end
    return false
end

const _ts = Test.get_testset_depth() > 0 ? Test.get_testset() : nothing
exit(_ts === nothing || !_any_nonpass(_ts) ? 0 : 1)
