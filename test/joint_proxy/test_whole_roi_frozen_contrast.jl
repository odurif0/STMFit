#!/usr/bin/env julia

# Red-first focused tests for test/lib/joint_proxy/whole_roi_frozen_contrast.jl.
# Diagnostic-only: this file does not enter production fitting, selection,
# calibration, or unit assignment. No benchmark label, expected N, NKNNKN,
# class count, or composition is read anywhere below.
#
# Contract verified here (Todo 4 of common-template-frozen-contrast-v2):
#   - the scorer accepts ONLY the immutable frozen input and verifies its
#     common_image / mask hashes and grid identity before scoring;
#   - it enumerates ONLY binary sequences under the already chosen
#     direction/phase/mirror, with no path to registration/global-state helpers;
#   - for each sequence it fits background + backbone + frozen common + signed
#     contrast and ranks all 2^n unique sequences by SSE with a deterministic
#     tie-break;
#   - it reports best, runner-up, bitwise-complement SSE, incremental contrast
#     gain, best/runner and complement margins, and the contrast coefficient;
#   - it abstains with a NAMED reason for zero contrast gain, zero contrast
#     coefficient, nonfinite evidence, or best/runner/complement ties within the
#     configured scale-aware numerical tolerance — and never for any arbitrary
#     chemistry-confidence floor;
#   - common-only / identical-template data abstain; exact typed data recovers
#     the injected sequence; swapping type metadata inverts the recovered
#     sequence while the frozen common registration stays byte-identical;
#   - the scorer rejects a record whose frozen hashes or grid identity were
#     changed, and a record whose observation mask drifted from the frozen mask.

using Test
using LinearAlgebra
using Random
using Statistics

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "whole_roi_frozen_contrast.jl"))

const WR_REG = Main.JointProxyRegistry
const DC = Main.WholeRoiDiagnosticConfig

# --- synthetic entry builders (match the Todo 2/3 test convention) -----------

function build_entry(M0_pixels::Vector{Float64}, M1_pixels::Vector{Float64})
    length(M0_pixels) == length(M1_pixels) || throw(ArgumentError("M0/M1 length mismatch"))
    templates = WR_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        pix = typ == 0 ? copy(M0_pixels) : copy(M1_pixels)
        push!(templates, WR_REG.ProxyTemplate(typ, parity, mirror, pix))
    end
    src = WR_REG.ProxySource("synthetic", "frozen_contrast_test", "", "", 0.0, 1.0, true)
    return WR_REG.ProxyEntry(src, templates)
end

function build_swapped_entry(M0_pixels::Vector{Float64}, M1_pixels::Vector{Float64})
    # type 0 carries M1 pixels and vice-versa: a pure metadata swap.
    return build_entry(M1_pixels, M0_pixels)
end

const HALF, STEP, GRID_N = 0.2, 0.2, 3

function baseline_M0()
    p = zeros(Float64, 9); p[5] = 1.0; return p
end
function baseline_M1()
    p = zeros(Float64, 9); p[8] = 1.0; p[5] = 0.2; return p
end

# --- frozen-input builder (identity transform, no registration call) ---------

"""Build a FrozenScorerInput at the identity transform (zero shift, zero
rotation, zero blur) for the given direction/phase/mirror, computing the frozen
common image and the chemistry-blind common-only baseline SSE directly. This
does NOT call register_common: it constructs the contract the scorer requires,
so the test stays isolated from the registration lane."""
function identity_frozen_input(set::CommonContrastSet, centers, xs, ys,
                               observed, backbone; direction::Int=0,
                               phase::Int=0, mirror::Int=0)
    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP, theta=0.0,
        direction=direction, phase=phase, mirror=mirror)
    common_only = active_set_nnls(observed, [backbone, common_img])
    return FrozenScorerInput(direction, phase, mirror, 0.0, 0.0, 0.0, 0.0,
        common_img, Base.hash(common_img),
        isfinite.(observed) .& isfinite.(common_img) .& isfinite.(backbone),
        Base.hash(isfinite.(observed) .& isfinite.(common_img) .& isfinite.(backbone)),
        set.grid_n, common_only.sse)
end

"""Config matching config/joint_proxy_whole_roi_diagnostic.toml tie tolerances."""
function tie_config()
    return DC.DiagnosticConfig(
        DC.DiagnosticParamBounds(DC.DiagnosticRange(0.0, 0.0, 0.01),
                                 DC.DiagnosticRange(0.0, 0.0, 0.01)),
        DC.DiagnosticParamBounds(DC.DiagnosticRange(0.0, 0.0, 0.01),
                                 DC.DiagnosticRange(0.0, 0.0, 0.01)),
        DC.DiagnosticParamBounds(DC.DiagnosticRange(0.0, 0.0, 0.01),
                                 DC.DiagnosticRange(0.0, 0.0, 0.01)),
        DC.DiagnosticParamBounds(DC.DiagnosticRange(0.0, 0.0, 0.01),
                                 DC.DiagnosticRange(0.0, 0.0, 0.01)),
        1.0e-6, 1.0e-9)
end

function make_backbone(xs::Vector{Float64}, ys::Vector{Float64}, sigma::Float64)
    return [exp(-(x^2 + y^2) / (2.0 * sigma^2)) for y in ys, x in xs]
end

# Common geometry for most tests.
const XS = collect(-0.6:0.1:0.6)
const YS = collect(-0.4:0.1:0.4)
const CENTERS = [(-0.3, 0.0), (-0.1, 0.0), (0.1, 0.0), (0.3, 0.0)]
const SEQ_INJECTED = Int[1, 0, 1, 0]

# --- integrity / tamper tests ------------------------------------------------

@testset "tampered common_image hash is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    # Corrupt the stored hash but keep the image: scorer must detect mismatch.
    bad = FrozenScorerInput(input.direction, input.phase, input.mirror,
        input.theta_total, input.shift_t_nm, input.shift_u_nm, input.blur_sigma_nm,
        input.common_image, 0xdeadbeef % UInt64, input.mask, input.mask_hash,
        input.grid_n, input.common_only_sse)
    @test_throws ArgumentError score_frozen_contrast(
        bad, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

@testset "tampered common_image pixels are detected via hash" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    tampered_image = copy(input.common_image)
    tampered_image[1, 1] += 1.0e-3
    bad = FrozenScorerInput(input.direction, input.phase, input.mirror,
        input.theta_total, input.shift_t_nm, input.shift_u_nm, input.blur_sigma_nm,
        tampered_image, input.common_image_hash, input.mask, input.mask_hash,
        input.grid_n, input.common_only_sse)
    @test_throws ArgumentError score_frozen_contrast(
        bad, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

@testset "tampered mask hash is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    bad = FrozenScorerInput(input.direction, input.phase, input.mirror,
        input.theta_total, input.shift_t_nm, input.shift_u_nm, input.blur_sigma_nm,
        input.common_image, input.common_image_hash, input.mask, 0x0badbad % UInt64,
        input.grid_n, input.common_only_sse)
    @test_throws ArgumentError score_frozen_contrast(
        bad, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

@testset "grid identity mismatch is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    bad = FrozenScorerInput(input.direction, input.phase, input.mirror,
        input.theta_total, input.shift_t_nm, input.shift_u_nm, input.blur_sigma_nm,
        input.common_image, input.common_image_hash, input.mask, input.mask_hash,
        99, input.common_only_sse)   # wrong grid_n
    @test_throws ArgumentError score_frozen_contrast(
        bad, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

@testset "ROI shape mismatch is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    wrong_observed = zeros(Float64, 5, 5)
    wrong_backbone = ones(Float64, 5, 5)
    @test_throws ArgumentError score_frozen_contrast(
        input, set, wrong_observed, wrong_backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

@testset "frozen state must be 0/1" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    for bad_direction in (2, -1)
        bad = FrozenScorerInput(bad_direction, input.phase, input.mirror,
            input.theta_total, input.shift_t_nm, input.shift_u_nm, input.blur_sigma_nm,
            input.common_image, input.common_image_hash, input.mask, input.mask_hash,
            input.grid_n, input.common_only_sse)
        @test_throws ArgumentError score_frozen_contrast(
            bad, set, observed, backbone, CENTERS, XS, YS;
            half_nm=HALF, step_nm=STEP, config=tie_config())
    end
end

@testset "stale grid (half/step implying wrong grid_n) is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    # half=0.3 step=0.1 -> grid_n = 7 != GRID_N. The scorer must reject this
    # before enumerating, because the frozen mold grid is stale for these params.
    @test_throws ArgumentError score_frozen_contrast(
        input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=0.3, step_nm=0.1, config=tie_config())
end

@testset "max_n exceeded is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    @test_throws ArgumentError score_frozen_contrast(
        input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config(), max_n=2)
end

@testset "n must be >= 1" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    @test_throws ArgumentError score_frozen_contrast(
        input, set, observed, backbone, Float64[], XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

# --- exact recovery ----------------------------------------------------------

@testset "exact typed data recovers the injected sequence (identity transform)" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, direction=0, phase=0, mirror=0)
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, direction=0, phase=0, mirror=0)
    # Independent physical synthesis: bg + backbone + common + signed contrast.
    bg, bb_c, common_c, contrast_c = 0.37, 1.4, 0.9, 0.55
    observed = bg .+ bb_c .* backbone .+ common_c .* common_img .+ contrast_c .* contrast_img

    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone;
        direction=0, phase=0, mirror=0)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())

    @test result.best.sequence == SEQ_INJECTED
    @test result.best.contrast_coef ≈ contrast_c atol = 1e-8
    @test result.best.common_coef ≈ common_c atol = 1e-8
    @test result.best.backbone_coef ≈ bb_c atol = 1e-8
    @test result.best.background ≈ bg atol = 1e-8
    @test result.best.sse ≤ 1e-16
    @test result.incremental_contrast_gain > 0
    @test !result.abstain
    @test isempty(result.abstention_reasons)
    # Full ranking must contain every unique binary sequence once.
    @test length(result.ranking) == 2^length(CENTERS)
    @test issorted([(s.sse, s.bits) for s in result.ranking]; by=identity)
    @test first(result.ranking).sequence == SEQ_INJECTED
    # Runner-up must be a different sequence and strictly worse within tol.
    @test result.runner_up !== nothing
    @test result.runner_up.sequence != result.best.sequence
    @test result.best_runner_margin > 0
end

@testset "complement of best is reported and its margin is positive (non-abstaining)" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.37 .+ 1.4 .* backbone .+ 0.9 .* common_img .+ 0.55 .* contrast_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())

    expected_complement = Int[1 - b for b in SEQ_INJECTED]
    @test result.complement.sequence == expected_complement
    @test result.complement_margin > 0
    # Complement SSE must be (approximately) the common-only baseline: under the
    # NNLS sign constraint the complement's contrast coefficient is clamped to 0.
    @test result.complement.contrast_coef == 0.0
    @test result.complement.sse ≈ input.common_only_sse rtol = 1e-6
end

# --- null / abstention -------------------------------------------------------

@testset "common-only data abstains with zero contrast gain" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    # Pure common-only observation: no contrast signal of any kind.
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())

    @test result.abstain
    @test "zero_contrast_gain" in result.abstention_reasons
    @test result.incremental_contrast_gain ≥ -1e-12
    # No sequence should be reported as a decisive chemical assignment here.
    @test result.best.contrast_coef == 0.0
end

@testset "identical type templates -> exactly zero contrast, full tie, abstention" begin
    # M0 == M1: contrast is exactly zero for every (parity, mirror). Every
    # sequence produces the same contrast image (zero) and the same SSE, so the
    # best/runner and complement margins are exactly zero and multiple named
    # reasons fire. This is the canonical exact numerical tie.
    entry = build_entry(baseline_M0(), baseline_M0())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())

    @test result.abstain
    @test "zero_contrast_gain" in result.abstention_reasons
    @test "complement_tie" in result.abstention_reasons
    @test "best_runner_tie" in result.abstention_reasons
    @test result.incremental_contrast_gain == 0.0
    @test result.complement_margin == 0.0
    @test result.best_runner_margin == 0.0
    # Every sequence genuinely ties at the common-only SSE.
    @test all(s -> isapprox(s.sse, result.best.sse; atol=1e-12), result.ranking)
end

@testset "weak contrast does NOT abstain (no arbitrary confidence floor)" begin
    # A small but clearly nonzero contrast signal (well above the numerical tie
    # tolerance) with a unique winner must NOT abstain: the only abstention
    # triggers are structural (zero/tie/nonfinite), not an arbitrary floor.
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    weak_c = 0.05   # small but >> numerical_tie (rtol=1e-6, atol=1e-9)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img .+ weak_c .* contrast_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())

    @test !result.abstain
    @test result.best.sequence == SEQ_INJECTED
    @test result.best.contrast_coef ≈ weak_c atol = 1e-8
    @test result.incremental_contrast_gain > 0
end

# --- swap invariance / inversion ---------------------------------------------

@testset "swapping type metadata inverts the recovered sequence" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    swapped_entry = build_swapped_entry(baseline_M0(), baseline_M1())
    set_swapped = derive_common_contrast(swapped_entry)

    backbone = make_backbone(XS, YS, 0.5)
    # Frozen common image: identical for both sets because common is the exact
    # average and is invariant under a type-0/type-1 metadata swap.
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    common_img_swapped = assemble_common_image(set_swapped, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    @test common_img == common_img_swapped

    # Synthesize observed using the ORIGINAL set and the injected sequence.
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img .+ 0.55 .* contrast_img

    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    # The frozen contract is identical for both scorings (common is swap-invariant
    # and common_only_sse is the same); only the contrast templates differ.
    result_orig = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    result_swap = score_frozen_contrast(input, set_swapped, observed, backbone,
        CENTERS, XS, YS; half_nm=HALF, step_nm=STEP, config=tie_config())

    @test result_orig.best.sequence == SEQ_INJECTED
    @test result_swap.best.sequence == Int[1 - b for b in SEQ_INJECTED]
    @test !result_orig.abstain && !result_swap.abstain
    # The recovered contrast coefficient flips sign in magnitude parity: the
    # swapped set's contrast is -1x the original, so the coefficient matches.
    @test result_swap.best.contrast_coef ≈ result_orig.best.contrast_coef atol = 1e-8
end

@testset "frozen common registration is byte-identical under swap (sentinel)" begin
    # Chemistry-blind registration output is invariant under a metadata swap by
    # construction; the scorer must rely on that property and not re-register.
    # Here we assert the property the scorer relies on: assemble_common_image is
    # swap-invariant for every (direction, phase, mirror).
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    swapped_set = derive_common_contrast(build_swapped_entry(baseline_M0(), baseline_M1()))
    for direction in (0, 1), phase in (0, 1), mirror in (0, 1)
        a = assemble_common_image(set, CENTERS, XS, YS; half_nm=HALF, step_nm=STEP,
            direction=direction, phase=phase, mirror=mirror)
        b = assemble_common_image(swapped_set, CENTERS, XS, YS; half_nm=HALF, step_nm=STEP,
            direction=direction, phase=phase, mirror=mirror)
        @test a == b
    end
end

# --- determinism / immutability ----------------------------------------------

@testset "scoring is deterministic across repeated calls" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img .+ 0.55 .* contrast_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    r1 = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    r2 = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    @test r1.best.sequence == r2.best.sequence
    @test r1.best.sse == r2.best.sse
    @test (r1.best.background, r1.best.backbone_coef, r1.best.common_coef,
           r1.best.contrast_coef) ==
          (r2.best.background, r2.best.backbone_coef, r2.best.common_coef,
           r2.best.contrast_coef)
    # Structs have no automatic value-==; compare the full ranking by a tuple
    # signature (vectors/tuples of plain numbers compare by value).
    sig(rk) = [(s.sequence, s.bits, s.sse, s.background, s.backbone_coef,
                s.common_coef, s.contrast_coef, s.n_pixels) for s in rk]
    @test sig(r1.ranking) == sig(r2.ranking)
    @test r1.abstention_reasons == r2.abstention_reasons
    @test r1.complement_margin == r2.complement_margin
end

@testset "frozen input is immutable and scorer does not mutate it" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    @test isimmutable(input)
    snapshot_common = copy(input.common_image)
    snapshot_mask = copy(input.mask)
    score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    @test input.common_image == snapshot_common
    @test input.mask == snapshot_mask
end

# --- nonfinite / dirty -------------------------------------------------------

@testset "NaN in observed is masked consistently with frozen mask" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img .+ 0.55 .* contrast_img
    observed[1, 1] = NaN
    observed[3, 5] = NaN
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    @test isfinite(result.best.sse)
    @test result.best.sequence == SEQ_INJECTED
    @test !result.abstain
end

@testset "mask drift (new NaN after freezing) is rejected" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    # Introduce a NEW nonfinite pixel AFTER freezing. The frozen mask no longer
    # matches the observed finite-pixel set, so the scorer must reject it rather
    # than silently scoring against a drifted baseline.
    drifted = copy(observed)
    drifted[2, 2] = NaN
    @test_throws ArgumentError score_frozen_contrast(
        input, set, drifted, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
end

@testset "nonfinite common_only_sse propagates as abstention (not a crash)" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    input = identity_frozen_input(set, CENTERS, XS, YS, observed, backbone)
    bad = FrozenScorerInput(input.direction, input.phase, input.mirror,
        input.theta_total, input.shift_t_nm, input.shift_u_nm, input.blur_sigma_nm,
        input.common_image, input.common_image_hash, input.mask, input.mask_hash,
        input.grid_n, NaN)
    result = score_frozen_contrast(bad, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    @test result.abstain
    @test "nonfinite_evidence" in result.abstention_reasons
end

# --- single-lobe edge case ---------------------------------------------------

@testset "single lobe: two complementary sequences only" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    centers = [(0.0, 0.0)]
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, centers, XS, YS;
        half_nm=HALF, step_nm=STEP)
    contrast_img = assemble_contrast_image(set, Int[1], centers, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img .+ 0.5 .* contrast_img
    common_only = active_set_nnls(observed, [backbone, common_img])
    input = FrozenScorerInput(0, 0, 0, 0.0, 0.0, 0.0, 0.0, common_img,
        Base.hash(common_img), isfinite.(observed) .& isfinite.(common_img) .& isfinite.(backbone),
        Base.hash(isfinite.(observed) .& isfinite.(common_img) .& isfinite.(backbone)),
        set.grid_n, common_only.sse)
    result = score_frozen_contrast(input, set, observed, backbone, centers, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    @test length(result.ranking) == 2
    @test result.best.sequence == [1]
    @test result.complement.sequence == [0]
    @test !result.abstain
end

# --- adapter compatibility with Todo 3 FrozenRegistration --------------------

@testset "adapter consumes a Todo 3 FrozenRegistration (integration smoke)" begin
    # Build a real FrozenRegistration via register_common (identity transform)
    # and confirm the adapter + scorer round-trip recovers the injected sequence.
    if !isdefined(Main, :FrozenRegistration)
        include(joinpath(dirname(@__DIR__), "lib", "joint_proxy",
                         "whole_roi_common_registration.jl"))
    end
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    backbone = make_backbone(XS, YS, 0.5)
    common_img = assemble_common_image(set, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    contrast_img = assemble_contrast_image(set, SEQ_INJECTED, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img .+ 0.55 .* contrast_img

    zero_range = DC.DiagnosticRange(0.0, 0.0, 0.01)
    reg_config = DC.DiagnosticConfig(
        DC.DiagnosticParamBounds(zero_range, zero_range),
        DC.DiagnosticParamBounds(zero_range, zero_range),
        DC.DiagnosticParamBounds(zero_range, zero_range),
        DC.DiagnosticParamBounds(zero_range, zero_range),
        1e-6, 1e-9)
    reg = register_common(set, observed, backbone, CENTERS, XS, YS;
        config=reg_config, half_nm=HALF, step_nm=STEP, theta_base=0.0)
    input = frozen_scorer_input(reg; theta_base=0.0)
    result = score_frozen_contrast(input, set, observed, backbone, CENTERS, XS, YS;
        half_nm=HALF, step_nm=STEP, config=tie_config())
    @test result.best.sequence == SEQ_INJECTED
    @test !result.abstain
end

# --- source guards -----------------------------------------------------------

@testset "API signature: no sequence/registration parameter" begin
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_frozen_contrast.jl"))
    src = read(helper, String)
    m = match(r"function score_frozen_contrast\(([^)]+)\)", src)
    @test m !== nothing
    sig = lowercase(m.captures[1])
    # The scorer ENUMERATES sequences internally; it must not receive one, and it
    # must not take any registration/contrast-template handle.
    forbidden = ["sequence::", "seq::", "register", "expected_n", "truth"]
    found = filter(p -> occursin(p, sig), forbidden)
    @test isempty(found)
end

@testset "no registration/global-state helper is called" begin
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_frozen_contrast.jl"))
    src = read(helper, String)
    no_comments = replace(src, r"#.*" => "")
    no_strings = replace(no_comments, r"\"[^\"]*\"" => "\"\"")
    lower = lowercase(no_strings)
    # The scorer must not reach registration or re-derive global state. Pure
    # coordinate/template/parity/blur helpers (_sample_template, _whole_roi_parity,
    # _gaussian_blur) ARE permitted: they are not registration.
    forbidden_calls = ["register_common", "_evaluate_point", "_assemble_common_shifted",
                       "register_global", "search_whole_roi_sequences"]
    leaked = filter(t -> occursin(t, lower), forbidden_calls)
    @test isempty(leaked)
end

@testset "runtime source guard: no executable benchmark-truth/composition read" begin
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_frozen_contrast.jl"))
    src = read(helper, String)
    no_comments = replace(src, r"#.*" => "")
    no_strings = replace(no_comments, r"\"[^\"]*\"" => "\"\"")
    lower = lowercase(no_strings)
    forbidden = ["nknnkn", "010010", "101101", "ground_truth", "ground truth",
                 "benchmark", "composition", "topk", "top_k", "expected_n",
                 "expectedn", "confidence_floor"]
    leaked = filter(tok -> occursin(tok, lower), forbidden)
    @test isempty(leaked)
    # Positive control: the strip is non-vacuous. The header comment documents
    # the prohibition using these tokens, so the RAW source contains them.
    raw_lower = lowercase(src)
    @test any(occursin(tok, raw_lower) for tok in ("nknnkn", "benchmark", "composition"))
    @test !any(occursin(tok, lower) for tok in ("nknnkn", "benchmark", "composition"))
end

@testset "no arbitrary confidence floor in abstention logic" begin
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_frozen_contrast.jl"))
    src = read(helper, String)
    no_comments = replace(src, r"#.*" => "")
    no_strings = replace(no_comments, r"\"[^\"]*\"" => "\"\"")
    # Abstention must be gated ONLY by the named structural checks and the
    # configured numerical tie tolerance, never by a hardcoded margin threshold.
    @test !occursin(r"gain\s*<\s*0\.[0-9]+", no_strings)
    @test !occursin(r"margin\s*<\s*0\.[0-9]+", no_strings)
end
