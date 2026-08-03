#!/usr/bin/env julia

# Red-first focused tests for test/lib/joint_proxy/whole_roi_common_contrast.jl.
# Diagnostic-only: this file does not enter production fitting, selection,
# calibration, or unit assignment. No benchmark label, expected N, NKNNKN,
# class count, or composition is read anywhere below.
#
# Contract verified here (Todo 2 of common-template-frozen-contrast-v2):
#   - exact M0 = common - contrast, M1 = common + contrast per (parity,mirror);
#   - chemistry-blind common assembly and signed typed-contrast assembly that
#     preserve the existing pixel ordering (u outer, t inner) and parity formula;
#   - generalized exhaustive active-set NNLS with a free intercept that never
#     returns a negative predictor slope;
#   - loud errors on incomplete (type,parity,mirror) state, stale/malformed
#     entries, and mismatched grid parameters.

using Test
using LinearAlgebra
using Random
using Statistics

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "whole_roi_common_contrast.jl"))

const WR_REG = Main.JointProxyRegistry

# --- synthetic entry builders ------------------------------------------------

"""Build an entry with identical M0/M1 pixel vectors for every (parity,mirror)."""
function build_entry(M0_pixels::Vector{Float64}, M1_pixels::Vector{Float64})
    length(M0_pixels) == length(M1_pixels) || throw(ArgumentError("M0/M1 length mismatch"))
    templates = WR_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        pix = typ == 0 ? copy(M0_pixels) : copy(M1_pixels)
        push!(templates, WR_REG.ProxyTemplate(typ, parity, mirror, pix))
    end
    src = WR_REG.ProxySource("synthetic", "whole_roi_common_contrast_test",
        "", "", 0.0, 1.0, true)
    return WR_REG.ProxyEntry(src, templates)
end

"""Build an entry from a per-(type,parity,mirror) pixel function (allows gaps)."""
function build_entry(pix_fn)
    templates = WR_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        pix = pix_fn(typ, parity, mirror)
        pix === nothing && continue
        push!(templates, WR_REG.ProxyTemplate(typ, parity, mirror, copy(pix)))
    end
    src = WR_REG.ProxySource("synthetic", "whole_roi_common_contrast_test",
        "", "", 0.0, 1.0, true)
    return WR_REG.ProxyEntry(src, templates)
end

# pixel vector convention matches _template_matrix: index = t + (u-1)*grid_n,
# coords = collect(-half_nm:step_nm:half_nm). grid_n = 3 -> coords [-0.2,0,0.2].
const HALF, STEP, GRID_N = 0.2, 0.2, 3

# M0: single center peak (t=0,u=0) -> pixel[5]. M1: shifted peak (t=0,u=0.2) ->
# pixel[8] plus a smaller center bump. M0 != M1 -> nonzero contrast.
function baseline_M0()
    p = zeros(Float64, 9); p[5] = 1.0; return p
end
function baseline_M1()
    p = zeros(Float64, 9); p[8] = 1.0; p[5] = 0.2; return p
end

# --- tests -------------------------------------------------------------------

@testset "common/contrast reconstruction (all 8 states)" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    @test set.grid_n == GRID_N
    @test length(set.molds) == 4  # one per (parity, mirror)
    for parity in (0, 1), mirror in (0, 1)
        mold = common_mold(set, parity, mirror)
        @test mold.parity == parity && mold.mirror == mirror
        @test mold.grid_n == GRID_N
        # Independently rebuild M0/M1 matrices in the same [iu,it] convention.
        M0 = Main._template_matrix(baseline_M0(), GRID_N)
        M1 = Main._template_matrix(baseline_M1(), GRID_N)
        @test maximum(abs.(mold.common .- mold.contrast .- M0)) ≤ 1e-12
        @test maximum(abs.(mold.common .+ mold.contrast .- M1)) ≤ 1e-12
        # common is exactly the average; contrast the signed half-difference.
        @test mold.common ≈ (M0 .+ M1) ./ 2 atol = 1e-12
        @test mold.contrast ≈ (M1 .- M0) ./ 2 atol = 1e-12
        # Never all-type-0: common must differ from M0 wherever M0 != M1.
        @test maximum(abs.(mold.common .- M0)) > 0
    end
end

@testset "identical type templates produce exactly zero contrast" begin
    entry = build_entry(baseline_M0(), baseline_M0())  # M0 == M1
    set = derive_common_contrast(entry)
    for parity in (0, 1), mirror in (0, 1)
        mold = common_mold(set, parity, mirror)
        @test all(iszero, mold.contrast)
        @test mold.common == Main._template_matrix(baseline_M0(), GRID_N)
    end
end

@testset "asymmetric peak ordering preserved (t vs u)" begin
    # M0 peak at (t=0.2, u=0.0) -> pixel[6]; M1 peak at (t=0.0, u=0.2) -> pixel[8].
    M0p = zeros(Float64, 9); M0p[6] = 1.0   # t-asymmetric
    M1p = zeros(Float64, 9); M1p[8] = 1.0   # u-asymmetric
    entry = build_entry(M0p, M1p)
    set = derive_common_contrast(entry)

    xs = collect(0.0:0.1:1.0)
    ys = collect(0.0:0.1:1.0)
    cx, cy = 0.5, 0.5

    # theta = 0: t maps to x, u maps to y. M0 peak (t=0.2) -> x=cx+0.2;
    # M1 peak (u=0.2) -> y=cy+0.2.
    common_img = assemble_common_image(set, [(cx, cy)], xs, ys;
        half_nm=HALF, step_nm=STEP, theta=0.0, direction=0, phase=0, mirror=0)
    # Two equal 0.5 peaks: at (x=cx+0.2, y=cy) and (x=cx, y=cy+0.2).
    @test common_img[findfirst(==(cy), ys), findfirst(==(cx + 0.2), xs)] ≈ 0.5 atol = 1e-12
    @test common_img[findfirst(==(cy + 0.2), ys), findfirst(==(cx), xs)] ≈ 0.5 atol = 1e-12

    # Signed contrast for sequence [1]: +0.5 at M1 site, -0.5 at M0 site.
    contrast_img = assemble_contrast_image(set, [1], [(cx, cy)], xs, ys;
        half_nm=HALF, step_nm=STEP, theta=0.0, direction=0, phase=0, mirror=0)
    (maxval, imax) = findmax(contrast_img)
    (minval, imin) = findmin(contrast_img)
    iy_max, ix_max = Tuple(imax)
    iy_min, ix_min = Tuple(imin)
    @test xs[ix_max] == cx && ys[iy_max] == cy + 0.2   # M1 -> +u direction
    @test xs[ix_min] == cx + 0.2 && ys[iy_min] == cy   # M0 -> +t direction
    @test maxval ≈ 0.5 atol = 1e-12
    @test minval ≈ -0.5 atol = 1e-12

    # theta = pi/2: rotation R sends (t,u) -> (x,y) with x = cos*t - sin*u...
    # Here dx = cos*t - sin*u? No: assemble uses t = dx*cos+dy*sin,
    # u = -dx*sin+dy*cos. For theta=pi/2, cos=0,sin=1 -> t=dy, u=-dx.
    # So feature at template (t=0.2,u=0) [M0] appears at dy=0.2,-dx=0 -> dx=0,dy=0.2.
    # Feature at template (t=0,u=0.2) [M1] appears at dy=0, -dx=0.2 -> dx=-0.2,dy=0.
    common_rot = assemble_common_image(set, [(cx, cy)], xs, ys;
        half_nm=HALF, step_nm=STEP, theta=pi / 2, direction=0, phase=0, mirror=0)
    @test common_rot[findfirst(==(cy + 0.2), ys), findfirst(==(cx), xs)] ≈ 0.5 atol = 1e-12
    @test common_rot[findfirst(==(cy), ys), findfirst(==(cx - 0.2), xs)] ≈ 0.5 atol = 1e-12
end

@testset "assembly equivalence: common + signed contrast == reference" begin
    # The chemistry-blind common plus the signed typed contrast must reproduce
    # the existing assemble_whole_roi reference for every sequence, parity
    # assignment, and mirror. This guards parity-formula and pixel-ordering
    # divergence without duplicating the reference.
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    xs = collect(-0.6:0.1:0.8)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.3, 0.05), (0.0, -0.05), (0.3, 0.05), (0.6, -0.05)]
    seqs = [[0, 0, 0, 0], [1, 1, 1, 1], [0, 1, 0, 1], [1, 0, 1, 0]]
    for direction in (0, 1), phase in (0, 1), mirror in (0, 1), seq in seqs
        ref = Main.assemble_whole_roi(entry, seq, centers, xs, ys;
            half_nm=HALF, step_nm=STEP, theta=0.0,
            direction=direction, phase=phase, mirror=mirror)
        common_img = assemble_common_image(set, centers, xs, ys;
            half_nm=HALF, step_nm=STEP, theta=0.0,
            direction=direction, phase=phase, mirror=mirror)
        contrast_img = assemble_contrast_image(set, seq, centers, xs, ys;
            half_nm=HALF, step_nm=STEP, theta=0.0,
            direction=direction, phase=phase, mirror=mirror)
        @test maximum(abs.((common_img .+ contrast_img) .- ref)) ≤ 1e-12
        # Common component must be sequence-independent.
        @test common_img == assemble_common_image(set, centers, xs, ys;
            half_nm=HALF, step_nm=STEP, theta=0.0,
            direction=direction, phase=phase, mirror=mirror)
    end
end

@testset "analytic active-set NNLS: exact positive recovery" begin
    # y = 3 + 2*p1 + 0.5*p2 exactly, both coefficients positive -> full active set.
    p1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    p2 = [2.0, 1.0, 0.5, 0.25, 0.1]
    y = 3.0 .+ 2.0 .* p1 .+ 0.5 .* p2
    fit = active_set_nnls(y, [p1, p2])
    @test fit.background ≈ 3.0 atol = 1e-10
    @test fit.coefficients ≈ [2.0, 0.5] atol = 1e-10
    @test fit.sse ≤ 1e-20
    @test sort(fit.active) == [1, 2]

    # Single-predictor analytic case.
    y1 = 1.5 .+ 0.7 .* p1
    fit1 = active_set_nnls(y1, [p1])
    @test fit1.background ≈ 1.5 atol = 1e-10
    @test fit1.coefficients ≈ [0.7] atol = 1e-10
    @test fit1.sse ≤ 1e-20
end

@testset "analytic active-set NNLS: negative coefficient forced to zero" begin
    # True coefficient of p1 is negative; the active set containing p1 is
    # rejected, so p1 must be dropped (coefficient exactly zero). The intercept
    # absorbs the mean of y.
    p1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = 3.0 .- 1.0 .* p1            # exact negative slope on p1
    fit = active_set_nnls(y, [p1])
    @test fit.coefficients[1] == 0.0
    @test fit.active == Int[]       # only the free intercept remains
    @test fit.background ≈ mean(y) atol = 1e-10
    # With a genuinely useful positive predictor that is orthogonal to p1 and to
    # the constant, dropping p1 leaves the positive slope unbiased.
    q1 = [-2.0, -1.0, 0.0, 1.0, 2.0]            # mean 0
    q2 = [1.0, -2.0, 0.0, 2.0, -1.0]            # mean 0, orthogonal to q1
    @test dot(q1, q2) ≈ 0.0 atol = 1e-12
    y2 = 1.0 .+ 0.8 .* q2 .- 1.0 .* q1
    fit2 = active_set_nnls(y2, [q1, q2])
    @test fit2.coefficients[1] == 0.0          # q1 dropped
    @test fit2.coefficients[2] ≈ 0.8 atol = 1e-10
    @test fit2.background ≈ 1.0 atol = 1e-10
end

@testset "active-set never returns a negative slope (randomized)" begin
    rng = MersenneTwister(20240725)
    for _ in 1:200
        n = 4
        A = randn(rng, n, 3)
        p1, p2, p3 = A[:, 1], A[:, 2], A[:, 3]
        true_coef = randn(rng, 3)
        y = 0.5 .+ A * true_coef
        fit = active_set_nnls(y, [p1, p2, p3])
        @test all(>=(0.0), fit.coefficients)
        @test fit.sse ≥ 0.0
        @test isfinite(fit.background)
    end
end

@testset "active-set is deterministic (flaky probe)" begin
    rng = MersenneTwister(7)
    A = randn(rng, 8, 3)
    y = randn(rng, 8)
    predictors = [A[:, 1], A[:, 2], A[:, 3]]
    a = active_set_nnls(y, predictors)
    b = active_set_nnls(y, predictors)
    @test a.sse == b.sse
    @test a.coefficients == b.coefficients
    @test a.active == b.active
    @test a.background == b.background
end

@testset "matrix-form active-set matches fit_whole_roi_nuisance" begin
    # The generalized solver with 2 predictors must reproduce the existing
    # 2-predictor exhaustive solver exactly on identical finite inputs.
    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    entry = build_entry(baseline_M0(), baseline_M1())
    backbone = [exp(-((x + 0.1)^2 + y^2) / 0.35) for y in ys, x in xs]
    mold = Main.assemble_whole_roi(entry, [0, 1], [(-0.2, 0.0), (0.2, 0.0)], xs, ys;
        half_nm=HALF, step_nm=STEP, direction=0, phase=0, mirror=0)
    observed = 0.7 .+ 1.4 .* backbone .+ 0.55 .* mold

    ref = Main.fit_whole_roi_nuisance(observed, backbone, mold)
    gen = active_set_nnls(observed, [backbone, mold])
    @test gen.sse ≈ ref.sse atol = 1e-9
    @test gen.background ≈ ref.background atol = 1e-9
    @test gen.coefficients[1] ≈ ref.backbone_coefficient atol = 1e-9
    @test gen.coefficients[2] ≈ ref.mold_coefficient atol = 1e-9
    @test gen.sse ≤ 1e-18
end

@testset "incomplete (type,parity,mirror) state -> loud error" begin
    # Missing type 1 for (parity=1, mirror=0): a naive all-type-0 fallback must
    # NOT happen; derive must throw an ArgumentError naming the missing state.
    full_M0, full_M1 = baseline_M0(), baseline_M1()
    incomplete = build_entry() do typ, parity, mirror
        (typ == 1 && parity == 1 && mirror == 0) ? nothing :
        typ == 0 ? full_M0 : full_M1
    end
    @test_throws ArgumentError derive_common_contrast(incomplete)

    # Missing an entire type for one (parity,mirror) is caught via the same path.
    one_type_only = build_entry() do typ, parity, mirror
        typ == 0 ? full_M0 : nothing
    end
    @test_throws ArgumentError derive_common_contrast(one_type_only)

    # Asking for a (parity,mirror) pair that was never built must also error.
    set = derive_common_contrast(build_entry(full_M0, full_M1))
    @test_throws ArgumentError common_mold(set, 2, 0)
    @test_throws KeyError common_mold(set, 0, 3)
end

@testset "malformed entries -> loud error" begin
    # Duplicate (type,parity,mirror) states: _whole_roi_template throws because
    # exactly one match is expected.
    M0, M1 = baseline_M0(), baseline_M1()
    dup_templates = WR_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        push!(dup_templates, WR_REG.ProxyTemplate(typ, parity, mirror,
            typ == 0 ? copy(M0) : copy(M1)))
    end
    push!(dup_templates, WR_REG.ProxyTemplate(0, 0, 0, copy(M0)))  # duplicate
    src = WR_REG.ProxySource("synthetic", "dup", "", "", 0.0, 1.0, true)
    @test_throws ArgumentError derive_common_contrast(WR_REG.ProxyEntry(src, dup_templates))

    # Wrong pixel length (not a perfect square): _template_matrix throws.
    bad_M0 = zeros(Float64, 8)
    bad_M1 = zeros(Float64, 8)
    @test_throws ArgumentError derive_common_contrast(build_entry(bad_M0, bad_M1))

    # Inconsistent pixel length across states.
    mixed = build_entry() do typ, parity, mirror
        (typ == 0 && parity == 0 && mirror == 0) ? zeros(Float64, 16) :
        typ == 0 ? copy(M0) : copy(M1)
    end
    @test_throws ArgumentError derive_common_contrast(mixed)
end

@testset "stale: mismatched grid parameters at assembly -> error" begin
    set = derive_common_contrast(build_entry(baseline_M0(), baseline_M1()))
    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    # half/step implying a different grid_n than the templates (3) must error.
    @test_throws ArgumentError assemble_common_image(set, [(0.0, 0.0)], xs, ys;
        half_nm=0.3, step_nm=0.1)  # grid_n = length(-0.3:0.1:0.3) = 7
    @test_throws ArgumentError assemble_contrast_image(set, [0], [(0.0, 0.0)], xs, ys;
        half_nm=0.3, step_nm=0.1)
    # Correct grid_n does not error.
    img = assemble_common_image(set, [(0.0, 0.0)], xs, ys;
        half_nm=HALF, step_nm=STEP)
    @test size(img) == (length(ys), length(xs))
end

@testset "dirty: NaN in templates propagates to assembly (documented)" begin
    # NaN pixels are not silently masked during pure template arithmetic; the
    # assembled image carries NaN where the NaN template contributes. This is
    # the documented behavior for the diagnostic representation itself (the
    # finite-pixel mask is a fit-time concern, handled by the matrix-form NNLS).
    nan_M0 = zeros(Float64, 9); nan_M0[5] = NaN
    M1 = baseline_M1()
    set = derive_common_contrast(build_entry(nan_M0, M1))
    xs = collect(-0.2:0.2:0.2)
    ys = collect(-0.2:0.2:0.2)
    img = assemble_common_image(set, [(0.0, 0.0)], xs, ys;
        half_nm=HALF, step_nm=STEP, direction=0, phase=0, mirror=0)
    @test any(isnan, img)

    # The matrix-form NNLS must tolerate NaN/Inf via the finite-pixel good mask
    # (matching fit_whole_roi_nuisance) provided enough finite pixels remain.
    good_y = ones(Float64, 4, 4)
    good_y[1, 1] = NaN
    p = reshape(1.0:16.0, 4, 4)
    fit = active_set_nnls(good_y, [p])
    @test isfinite(fit.sse)
    @test all(>=(0.0), fit.coefficients)
end

@testset "real registry smoke: derive works on the DFT entry" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    reg_path = joinpath(root, "config", "joint_proxy_molds.toml")
    isfile(reg_path) || return  # skip when the untracked config is absent
    registry = WR_REG.load_registry(reg_path)
    entries = [e for e in registry.entries if e.source.family == "stm_dft_v1"]
    isempty(entries) && return
    entry = first(entries)
    set = derive_common_contrast(entry)
    @test set.grid_n == registry.grid_n
    @test length(set.molds) == 4
    for mold in values(set.molds)
        @test all(isfinite, mold.common)
        @test all(isfinite, mold.contrast)
        # common and contrast reconstruct both typed templates.
        M0 = Main._template_matrix(
            Main._whole_roi_template(entry, 0, mold.parity, mold.mirror), set.grid_n)
        M1 = Main._template_matrix(
            Main._whole_roi_template(entry, 1, mold.parity, mold.mirror), set.grid_n)
        @test maximum(abs.(mold.common .- mold.contrast .- M0)) ≤ 1e-9
        @test maximum(abs.(mold.common .+ mold.contrast .- M1)) ≤ 1e-9
    end
end

# --- AdversarialVerify follow-up regressions ----------------------------------

@testset "rank-deficient candidate never throws (collinear with intercept)" begin
    # A constant predictor is collinear with the free intercept -> the [1]
    # design is rank-deficient. The solver must NOT throw SingularException; it
    # must skip the unidentifiable candidate and return a deterministic feasible
    # parsimonious solution (intercept-only here).
    fit = active_set_nnls([2.0, 2.0], [[1.0, 1.0]])
    @test fit.active == Int[]
    @test fit.coefficients == [0.0]
    @test fit.background ≈ 2.0 atol = 1e-9
    @test isfinite(fit.sse)

    # A useful predictor alongside a collinear one: only the useful predictor's
    # active sets are admissible; the collinear predictor is dropped to zero.
    p_const = [1.0, 1.0, 1.0, 1.0]            # collinear with intercept
    p_useful = [0.0, 1.0, 2.0, 3.0]           # full-rank with intercept
    y = 1.0 .+ 0.5 .* p_useful
    fit2 = active_set_nnls(y, [p_const, p_useful])
    @test fit2.coefficients[1] == 0.0          # collinear predictor dropped
    @test fit2.coefficients[2] ≈ 0.5 atol = 1e-9
    @test fit2.background ≈ 1.0 atol = 1e-9
    @test 2 in fit2.active && !(1 in fit2.active)

    # Two mutually-collinear constant predictors, constant y -> intercept-only.
    fit3 = active_set_nnls([5.0, 5.0], [[1.0, 1.0], [2.0, 2.0]])
    @test fit3.active == Int[]
    @test fit3.coefficients == [0.0, 0.0]
    @test fit3.background ≈ 5.0 atol = 1e-9   # mean([5,5]) == 5.0
end

@testset "deterministic tie-break across distinct equal-SSE active sets" begin
    # [] and [1] both fit y perfectly; the documented parsimony rule must pick
    # the smallest active set, deterministically, under BOTH default and tight
    # tolerances. (Previously the default tolerance let floating-point noise
    # select [1] while atol=1e-12 selected [].)
    for tol_kw in (() , (atol=1e-12,))
        fit = active_set_nnls([2.0, 2.0], [[-1.0, 1.0]]; tol_kw...)
        @test fit.active == Int[]
        @test fit.coefficients == [0.0]
    end

    # Two distinct predictors each fit y exactly; [1] and [2] tie and both beat
    # []. The tie-break must select [1] (popcount-equal, lexicographically first).
    y = [1.0, 3.0]
    p1 = [0.0, 1.0]                            # bg=1, c=2 -> y exact
    p2 = [0.0, 2.0]                            # bg=1, c=1 -> y exact
    fit = active_set_nnls(y, [p1, p2])
    @test fit.active == [1]
    @test fit.coefficients ≈ [2.0, 0.0] atol = 1e-9
    @test fit.background ≈ 1.0 atol = 1e-9
    # Same input must produce byte-identical results across repeated calls.
    again = active_set_nnls(y, [p1, p2])
    @test again.active == fit.active && again.coefficients == fit.coefficients &&
          again.sse == fit.sse && again.background == fit.background
end

@testset "runtime source guard: no executable benchmark-truth/composition read" begin
    # A runtime (not just grep-by-reviewer) guard: read the helper source,
    # strip Julia comments and string literals, and assert that no forbidden
    # token survives in the executable region. AGENTS.md forbids reading
    # expected N, NKNNKN, class counts, benchmark labels, or composition here.
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_common_contrast.jl"))
    src = read(helper, String)
    no_comments = replace(src, r"#.*" => "")                       # strip line comments
    no_strings = replace(no_comments, r"\"[^\"]*\"" => "\"\"")      # strip string literals
    lower = lowercase(no_strings)
    forbidden = ["nknnkn", "010010", "101101", "ground_truth", "ground truth",
                 "benchmark", "composition", "topk", "top_k", "expected_n",
                 "expectedn"]
    leaked = filter(tok -> occursin(tok, lower), forbidden)
    @test isempty(leaked)
    # Positive control: the guard is not vacuous. The helper's header comment
    # documents the prohibition using these tokens, so the RAW source contains
    # them; the strip must remove them. This proves the mechanism catches a
    # deliberate injection in executable code.
    raw_lower = lowercase(src)
    @test any(occursin(tok, raw_lower) for tok in ("nknnkn", "benchmark", "composition"))
    @test !any(occursin(tok, lower) for tok in ("nknnkn", "benchmark", "composition"))
end
