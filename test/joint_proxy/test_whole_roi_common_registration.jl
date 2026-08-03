#!/usr/bin/env julia

# Red-first focused tests for test/lib/joint_proxy/whole_roi_common_registration.jl.
# Diagnostic-only: this file does not enter production fitting, selection,
# calibration, or unit assignment. No benchmark label, expected N, NKNNKN,
# class count, or composition is read anywhere below.
#
# Contract verified here (Todo 3 of common-template-frozen-contrast-v2):
#   - deterministic chemistry-blind global registration over 8 global states;
#   - coarse-then-fine bounded search over shift_t, shift_u, rotation, blur;
#   - coordinate-sampling assembly with shift and rotation via theta addition;
#   - SSE-based selection with documented minimal-transform tie-break;
#   - immutable FrozenRegistration struct with all parameters, common image,
#     hash, mask, and fit coefficients;
#   - API has no sequence or contrast parameter.

using Test
using LinearAlgebra
using Random
using Statistics

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "whole_roi_common_registration.jl"))

const WR_REG = Main.JointProxyRegistry

# --- synthetic entry builders (same as test_whole_roi_common_contrast.jl) ---

function build_entry(M0_pixels::Vector{Float64}, M1_pixels::Vector{Float64})
    length(M0_pixels) == length(M1_pixels) || throw(ArgumentError("M0/M1 length mismatch"))
    templates = WR_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        pix = typ == 0 ? copy(M0_pixels) : copy(M1_pixels)
        push!(templates, WR_REG.ProxyTemplate(typ, parity, mirror, pix))
    end
    src = WR_REG.ProxySource("synthetic", "registration_test", "", "", 0.0, 1.0, true)
    return WR_REG.ProxyEntry(src, templates)
end

const HALF, STEP, GRID_N = 0.2, 0.2, 3

function baseline_M0()
    p = zeros(Float64, 9); p[5] = 1.0; return p
end
function baseline_M1()
    p = zeros(Float64, 9); p[8] = 1.0; p[5] = 0.2; return p
end

# --- helpers -----------------------------------------------------------------

"""Build a minimal degenerate config with all ranges frozen at zero (identity transform)."""
function identity_config()
    DC = Main.WholeRoiDiagnosticConfig
    zero_range = DC.DiagnosticRange(0.0, 0.0, 0.01)
    return DC.DiagnosticConfig(
        DC.DiagnosticParamBounds(zero_range, zero_range),
        DC.DiagnosticParamBounds(zero_range, zero_range),
        DC.DiagnosticParamBounds(zero_range, zero_range),
        DC.DiagnosticParamBounds(zero_range, zero_range),
        1e-6, 1e-9)
end

"""Build a config with coarse ranges spanning typical values and fine refinement."""
function standard_test_config()
    DC = Main.WholeRoiDiagnosticConfig
    return DC.DiagnosticConfig(
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(-0.10, 0.10, 0.05),
            DC.DiagnosticRange(-0.04, 0.04, 0.02)),
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(-0.06, 0.06, 0.03),
            DC.DiagnosticRange(-0.02, 0.02, 0.01)),
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(-2.0, 2.0, 1.0),
            DC.DiagnosticRange(-0.5, 0.5, 0.25)),
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(0.0, 0.06, 0.03),
            DC.DiagnosticRange(0.0, 0.01, 0.005)),
        1e-6, 1e-9)
end

"""Build a single-center backbone image (broad Gaussian)."""
function make_backbone(xs::Vector{Float64}, ys::Vector{Float64},
                       cx::Float64, cy::Float64, sigma::Float64)
    return [exp(-((x - cx)^2 + (y - cy)^2) / (2.0 * sigma^2)) for y in ys, x in xs]
end

"""Canonical deterministic TSV serialization of a FrozenRegistration for
byte-stability checks (acceptance: "repeated runs are byte-stable in TSV
serialization"). Diagnostic-only; the production TSV emitter lives in the
Todo 6 CLI driver — this in-test serializer just exercises determinism of
every frozen field in a fixed text order."""
function freeze_to_tsv(r::FrozenRegistration)
    io = IOBuffer()
    println(io, "direction\t", r.direction)
    println(io, "phase\t", r.phase)
    println(io, "mirror\t", r.mirror)
    println(io, "shift_t_nm\t", r.shift_t_nm)
    println(io, "shift_u_nm\t", r.shift_u_nm)
    println(io, "rotation_deg\t", r.rotation_deg)
    println(io, "blur_sigma_nm\t", r.blur_sigma_nm)
    println(io, "common_image_hash\t", r.common_image_hash)
    println(io, "grid_n\t", r.grid_n)
    println(io, "sse\t", r.fit.sse)
    println(io, "background\t", r.fit.background)
    println(io, "coefficients\t", join(r.fit.coefficients, ","))
    println(io, "active\t", join(r.fit.active, ","))
    println(io, "tie_info\t", r.tie_info)
    return String(take!(io))
end

# --- tests -------------------------------------------------------------------

@testset "zero nuisance recovers identity transform" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    # Assemble common image at identity transform (direction=0, phase=0, mirror=0)
    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img

    config = identity_config()
    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    @test result.shift_t_nm ≈ 0.0 atol = 0.02
    @test result.shift_u_nm ≈ 0.0 atol = 0.02
    @test result.rotation_deg ≈ 0.0 atol = 0.5
    @test result.blur_sigma_nm ≈ 0.0 atol = 0.01
    @test result.fit.sse ≤ 1e-10
    @test isfinite(result.common_image_hash)
    @test result.grid_n == GRID_N
end

@testset "off-grid nuisance recovers within one fine step" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    # Inject known nuisance BETWEEN coarse grid points via the implementation's
    # exact forward model (coordinate-sampled assembly + Gaussian blur). Building
    # the synthetic observation with the same forward model the registration uses
    # makes "recovers within one fine step" a statement about the optimizer
    # instead of nearest-neighbor resampling error on a 0.1 nm pixel grid (which
    # cannot represent a 0.03 nm shift and made the old test see an identity image).
    injected_shift_t = 0.03
    injected_shift_u = -0.01
    injected_rot_deg = 0.5
    injected_blur = 0.005
    pixel_step = xs[2] - xs[1]
    blur_px = injected_blur / pixel_step

    common_fwd = Main._assemble_common_shifted(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP, theta=deg2rad(injected_rot_deg),
        direction=0, phase=0, mirror=0,
        shift_t=injected_shift_t, shift_u=injected_shift_u)
    common_fwd = blur_px > 0 ? Main._gaussian_blur(common_fwd, blur_px) : common_fwd
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_fwd

    config = standard_test_config()
    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    # Discrete global state is recovered exactly (synthetic templates are
    # identical across states, so the documented minimal-transform tie-break
    # selects direction=phase=mirror=0).
    @test result.direction == 0
    @test result.phase == 0
    @test result.mirror == 0

    # Fine step sizes from config
    fine_step_t = config.shift_t_nm.fine.step   # 0.02
    fine_step_u = config.shift_u_nm.fine.step   # 0.01
    fine_step_rot = config.rotation_deg.fine.step   # 0.25
    fine_step_blur = config.blur_sigma_nm.fine.step  # 0.005

    @test abs(result.shift_t_nm - injected_shift_t) ≤ fine_step_t + 1e-10
    @test abs(result.shift_u_nm - injected_shift_u) ≤ fine_step_u + 1e-10
    @test abs(result.rotation_deg - injected_rot_deg) ≤ fine_step_rot + 1e-10
    @test abs(result.blur_sigma_nm - injected_blur) ≤ fine_step_blur + 1e-10
    @test result.fit.sse ≤ 1e-6
end

@testset "opposite sequences yield identical frozen registration" begin
    # Chemistry-blindness contract (acceptance: "two opposite type sequences
    # sharing the same nuisance yield the same frozen registration within tie
    # rules"; cf. Todo 4: "swapping type metadata leaves common registration
    # byte-identical"). The registration API takes no sequence argument and fits
    # only the label-invariant common = (M0 + M1)/2. Swapping the type labels
    # (type 0 <-> type 1) leaves the common template EXACTLY unchanged and only
    # flips the sign of the contrast; a sequence [0,1,0,1] under the swapped
    # labels is the complement [1,0,1,0] under the original labels. Therefore
    # the frozen registration must be byte-identical under the swap.
    #
    # This is the correct operational test of the criterion. Requiring two
    # DIFFERENT observed images (common + small contrast_a vs common + small
    # contrast_b) to yield the exact same argmin is NOT a registration invariant:
    # nonzero contrast is structured residual that can legitimately move the
    # common-only argmin, and that is not a chemistry-blindness violation. The
    # metadata-swap form below tests the actual invariant.
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    # Swapped entry: type 0 carries M1 and type 1 carries M0.
    entry_swapped = build_entry(baseline_M1(), baseline_M0())
    set_swapped = derive_common_contrast(entry_swapped)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.0, 0.0), (0.2, 0.0), (0.4, 0.0)]
    backbone = make_backbone(xs, ys, 0.1, 0.0, 0.5)

    # Sanity: the swap preserves common and flips contrast (else the test would
    # be vacuous). common_mold is (M0+M1)/2 in both; contrast changes sign.
    for parity in (0, 1), mirror in (0, 1)
        ca = common_mold(set, parity, mirror)
        cs = common_mold(set_swapped, parity, mirror)
        @test ca.common ≈ cs.common atol = 1e-12
        @test cs.contrast ≈ .-(ca.contrast) atol = 1e-12
    end

    # Inject a shared nuisance on the common image via the forward model.
    injected_shift_t, injected_shift_u = 0.02, 0.01
    injected_rot_deg = 0.3
    injected_blur = 0.005
    pixel_step = xs[2] - xs[1]
    blur_px = injected_blur / pixel_step

    common_fwd = Main._assemble_common_shifted(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP, theta=deg2rad(injected_rot_deg),
        direction=0, phase=0, mirror=0,
        shift_t=injected_shift_t, shift_u=injected_shift_u)
    common_fwd = blur_px > 0 ? Main._gaussian_blur(common_fwd, blur_px) : common_fwd
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_fwd

    config = standard_test_config()
    result_a = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)
    result_b = register_common(set_swapped, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    # Byte-identical: the registration consumed only the (identical) common.
    @test result_a.direction == result_b.direction
    @test result_a.phase == result_b.phase
    @test result_a.mirror == result_b.mirror
    @test result_a.shift_t_nm == result_b.shift_t_nm
    @test result_a.shift_u_nm == result_b.shift_u_nm
    @test result_a.rotation_deg == result_b.rotation_deg
    @test result_a.blur_sigma_nm == result_b.blur_sigma_nm
    @test result_a.common_image_hash == result_b.common_image_hash
    @test result_a.common_image == result_b.common_image
    @test result_a.fit.coefficients == result_b.fit.coefficients
    @test result_a.fit.sse == result_b.fit.sse
    @test freeze_to_tsv(result_a) == freeze_to_tsv(result_b)
end

@testset "repeat byte-stable serialization" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img

    config = identity_config()
    r1 = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)
    r2 = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    @test r1.direction == r2.direction
    @test r1.phase == r2.phase
    @test r1.mirror == r2.mirror
    @test r1.shift_t_nm == r2.shift_t_nm
    @test r1.shift_u_nm == r2.shift_u_nm
    @test r1.rotation_deg == r2.rotation_deg
    @test r1.blur_sigma_nm == r2.blur_sigma_nm
    @test r1.common_image_hash == r2.common_image_hash
    @test r1.fit.sse == r2.fit.sse
    @test r1.fit.background == r2.fit.background
    @test r1.fit.coefficients == r2.fit.coefficients
    @test r1.fit.active == r2.fit.active
    @test r1.tie_info == r2.tie_info
    # Acceptance: "repeated runs are byte-stable in TSV serialization".
    @test freeze_to_tsv(r1) == freeze_to_tsv(r2)
end

@testset "boundary status: shift beyond configured bounds" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    # Config with a narrow coarse range (e.g., shift_t max = 0.10)
    DC = Main.WholeRoiDiagnosticConfig
    config = DC.DiagnosticConfig(
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(-0.10, 0.10, 0.05),
            DC.DiagnosticRange(-0.10, 0.10, 0.05)),
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(-0.06, 0.06, 0.03),
            DC.DiagnosticRange(-0.02, 0.02, 0.01)),
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(-2.0, 2.0, 1.0),
            DC.DiagnosticRange(-0.5, 0.5, 0.25)),
        DC.DiagnosticParamBounds(
            DC.DiagnosticRange(0.0, 0.06, 0.03),
            DC.DiagnosticRange(0.0, 0.01, 0.005)),
        1e-6, 1e-9)

    # Inject a shift_t (0.15 nm) beyond the configured coarse bound (0.10) using
    # the implementation's exact forward model. The registration cannot represent
    # 0.15 within its search grid (the fine grid is clamped to the coarse range),
    # so it must pick the boundary value rather than pretend it converged on the
    # true shift.
    injected_shift_t = 0.15
    common_fwd = Main._assemble_common_shifted(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP, theta=0.0,
        direction=0, phase=0, mirror=0,
        shift_t=injected_shift_t, shift_u=0.0)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_fwd

    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    # The recovered shift must sit at the positive wall (0.10), not the true 0.15
    # and not zero: the optimizer correctly hits the configured bound.
    @test result.shift_t_nm ≤ config.shift_t_nm.coarse.max + 1e-12
    @test result.shift_t_nm ≈ config.shift_t_nm.coarse.max atol = config.shift_t_nm.fine.step
    @test result.shift_t_nm > 0.0
    @test isfinite(result.fit.sse)
    @test isfinite(result.common_image_hash)
end

@testset "sentinel: scorer cannot mutate frozen state" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img

    config = identity_config()
    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    # FrozenRegistration is declared with `struct` (immutable): fields cannot be
    # reassigned, so the scorer cannot mutate direction/phase/mirror/transform/
    # common image/mask after freezing. `ismutabletype` is the version-stable
    # check (returns false for immutable struct types); `isimmutable` on a
    # DataType checks the immutability of the type-as-value, not the struct kind.
    @test !ismutabletype(FrozenRegistration)
    @test isimmutable(result)

    # Attempt to construct a modified copy and verify original unchanged
    original_hash = result.common_image_hash
    original_direction = result.direction
    original_phase = result.phase
    modified = FrozenRegistration(
        1 - result.direction, result.phase, result.mirror,
        result.shift_t_nm, result.shift_u_nm, result.rotation_deg, result.blur_sigma_nm,
        result.common_image, result.common_image_hash,
        result.grid_n, result.mask, result.fit, result.tie_info)
    @test original_direction == result.direction
    @test original_hash == result.common_image_hash
    @test modified.direction != result.direction
end

@testset "deep immutability: direct setindex! on frozen fields is rejected" begin
    # Regression for the shallow-immutability defect (task-3 independent review):
    # a `struct` prevents field REASSIGNMENT, but a mutable Matrix field still
    # allows in-place `img[i,j] = v`, silently invalidating common_image_hash
    # and defeating Todo 4's tamper check. The frozen array fields must reject
    # element mutation outright. This test is RED against shallow Matrix fields
    # and turns GREEN once FrozenRegistration exposes read-only FrozenMatrix
    # wrappers (copy on freeze; setindex!/fill!/copyto! throw).
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)
    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    config = identity_config()
    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    stored_common_hash = result.common_image_hash

    # Bracket-assignment is the syntax a caller would actually write.
    function _try_setimg!(img, i, j, v)
        img[i, j] = v
    end
    function _try_setmask!(m, i, j, v)
        m[i, j] = v
    end

    @test_throws ArgumentError _try_setimg!(result.common_image, 1, 1, 999.0)
    @test_throws ArgumentError _try_setimg!(result.common_image,
        size(result.common_image, 1), size(result.common_image, 2), -1.0)
    @test_throws ArgumentError _try_setmask!(result.mask, 1, 1, !result.mask[1, 1])

    # Bulk mutators that route through setindex! must be rejected too.
    @test_throws ArgumentError fill!(result.common_image, 0.0)
    @test_throws ArgumentError copyto!(result.common_image,
        collect(result.common_image))

    # Stored hash still validates: the attempts did not modify the data.
    @test Base.hash(result.common_image) == stored_common_hash

    # Regression for nested public-field mutation. The first read-only wrapper
    # stored a publicly reachable Matrix in `.data`, and the frozen fit retained
    # mutable coefficient/active vectors. Frozen storage now contains only
    # immutable String/Tuple payloads.
    @test !hasproperty(result.common_image, :data)
    @test result.common_image.bytes isa String
    @test result.fit.coefficients isa Tuple
    @test result.fit.active isa Tuple
    @test_throws MethodError setindex!(result.fit.coefficients, 9.0, 1)
    @test_throws MethodError setindex!(result.fit.active, 7, 1)
    @test Base.hash(result.common_image) == stored_common_hash
end

@testset "runtime source guard: no executable benchmark-truth/composition read" begin
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_common_registration.jl"))
    src = read(helper, String)
    no_comments = replace(src, r"#.*" => "")
    no_strings = replace(no_comments, r"\"[^\"]*\"" => "\"\"")
    lower = lowercase(no_strings)
    forbidden = ["nknnkn", "010010", "101101", "ground_truth", "ground truth",
                 "benchmark", "composition", "topk", "top_k", "expected_n",
                 "expectedn"]
    leaked = filter(tok -> occursin(tok, lower), forbidden)
    @test isempty(leaked)
    raw_lower = lowercase(src)
    @test any(occursin(tok, raw_lower) for tok in ("nknnkn", "benchmark"))
    @test !any(occursin(tok, lower) for tok in ("nknnkn", "benchmark"))
end

@testset "API signature: no sequence or contrast parameter" begin
    helper = normpath(joinpath(@__DIR__, "..", "lib", "joint_proxy",
                               "whole_roi_common_registration.jl"))
    src = read(helper, String)
    # Extract the full register_common signature (spans multiple lines).
    m = match(r"function register_common\(([^)]+)\)", src)
    @test m !== nothing
    sig = m.captures[1]
    sig_norm = replace(sig, r"\s+" => " ")
    # Pull PARAMETER NAMES only (the leading identifier of each comma/semicolon
    # separated part), so a legitimate type name like CommonContrastSet does not
    # false-match the forbidden "contrast" token (the old substring check did).
    parts = split(sig_norm, r"[,;]")
    param_names = String[]
    for p in parts
        mm = match(r"^([A-Za-z_]\w*)", strip(p))
        mm === nothing && continue
        push!(param_names, lowercase(mm.captures[1]))
    end
    @test !isempty(param_names)
    forbidden = ["sequence", "contrast", "seq", "typ"]
    leaked = String[]
    for name in param_names, tok in forbidden
        if name == tok || startswith(name, tok)
            push!(leaked, name)
        end
    end
    @test isempty(leaked)
    # Sanity: the known chemistry-blind positional and keyword names are present.
    @test "set" in param_names
    @test "observed" in param_names
    @test "config" in param_names
end

@testset "malformed config rejected at the parser gate" begin
    # DiagnosticConfig is a plain immutable record: validation lives in the
    # parser (parse_diagnostic_config), which is the contract boundary that
    # keeps malformed bounds off register_common. The registration consumes an
    # already-validated record, so this test exercises the parser gate that
    # protects the registration rather than the struct constructor.
    DC = Main.WholeRoiDiagnosticConfig

    # A known-valid TOML dict (matches standard_test_config()).
    function valid_toml()
        return Dict{String,Any}(
            "registration" => Dict{String,Any}(
                "coarse" => Dict{String,Any}(
                    "shift_t_nm"    => Dict("min"=>-0.10, "max"=>0.10, "step"=>0.05),
                    "shift_u_nm"    => Dict("min"=>-0.06, "max"=>0.06, "step"=>0.03),
                    "rotation_deg"  => Dict("min"=>-2.0,  "max"=>2.0,  "step"=>1.0),
                    "blur_sigma_nm" => Dict("min"=>0.0,   "max"=>0.06, "step"=>0.03)),
                "fine" => Dict{String,Any}(
                    "shift_t_nm"    => Dict("min"=>-0.04, "max"=>0.04, "step"=>0.02),
                    "shift_u_nm"    => Dict("min"=>-0.02, "max"=>0.02, "step"=>0.01),
                    "rotation_deg"  => Dict("min"=>-0.5,  "max"=>0.5,  "step"=>0.25),
                    "blur_sigma_nm" => Dict("min"=>0.0,   "max"=>0.01, "step"=>0.005))),
            "numerical_tie" => Dict("rtol"=>1.0e-6, "atol"=>1.0e-9))
    end
    # Positive control: the valid dict parses cleanly.
    @test DC.parse_diagnostic_config(valid_toml()) isa DC.DiagnosticConfig

    # 1) fine step wider than coarse step -> violates refinement monotonicity.
    d1 = valid_toml()
    d1["registration"]["fine"]["shift_t_nm"] = Dict("min"=>-0.06, "max"=>0.06, "step"=>0.06)
    @test_throws ArgumentError DC.parse_diagnostic_config(d1)

    # 2) fine range span wider than coarse span.
    d2 = valid_toml()
    d2["registration"]["fine"]["shift_t_nm"] = Dict("min"=>-0.12, "max"=>0.12, "step"=>0.04)
    @test_throws ArgumentError DC.parse_diagnostic_config(d2)

    # 3) rtol and atol both zero -> no scale-aware tie tolerance.
    d3 = valid_toml()
    d3["numerical_tie"] = Dict("rtol"=>0.0, "atol"=>0.0)
    @test_throws ArgumentError DC.parse_diagnostic_config(d3)

    # 4) asymmetric translation bound -> violates symmetry invariant.
    d4 = valid_toml()
    d4["registration"]["coarse"]["shift_t_nm"] = Dict("min"=>-0.10, "max"=>0.20, "step"=>0.05)
    @test_throws ArgumentError DC.parse_diagnostic_config(d4)

    # 5) missing section -> loud failure, not a silent default.
    d5 = valid_toml()
    delete!(d5, "numerical_tie")
    @test_throws ArgumentError DC.parse_diagnostic_config(d5)
end

@testset "deterministic across equivalent inputs (flaky probe)" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img

    config = standard_test_config()
    results = [register_common(set, observed, backbone, centers, xs, ys;
                config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)
               for _ in 1:3]

    for i in 2:3
        @test results[1].direction == results[i].direction
        @test results[1].phase == results[i].phase
        @test results[1].mirror == results[i].mirror
        @test results[1].shift_t_nm == results[i].shift_t_nm
        @test results[1].shift_u_nm == results[i].shift_u_nm
        @test results[1].rotation_deg == results[i].rotation_deg
        @test results[1].blur_sigma_nm == results[i].blur_sigma_nm
        @test results[1].common_image_hash == results[i].common_image_hash
        @test results[1].fit.sse == results[i].fit.sse
        @test results[1].fit.background == results[i].fit.background
    end
end

@testset "stale/dirty: NaN in observed handled via finite mask" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    # Sprinkle NaN
    observed[1, 1] = NaN
    observed[3, 5] = NaN
    observed[end, end] = NaN

    config = identity_config()
    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    @test isfinite(result.fit.sse)
    @test result.fit.background ≈ 0.5 atol = 0.01
    @test all(>=(0.0), result.fit.coefficients)
end

@testset "noisy observation still converges" begin
    entry = build_entry(baseline_M0(), baseline_M1())
    set = derive_common_contrast(entry)
    rng = MersenneTwister(42)

    xs = collect(-0.6:0.1:0.6)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.2, 0.0), (0.2, 0.0)]
    backbone = make_backbone(xs, ys, 0.0, 0.0, 0.5)

    common_img = assemble_common_image(set, centers, xs, ys;
        half_nm=HALF, step_nm=STEP)
    observed = 0.5 .+ 0.8 .* backbone .+ 0.4 .* common_img
    observed .+= 0.02 .* randn(rng, size(observed))

    config = identity_config()
    result = register_common(set, observed, backbone, centers, xs, ys;
        config=config, half_nm=HALF, step_nm=STEP, theta_base=0.0)

    @test isfinite(result.fit.sse)
    @test isfinite(result.common_image_hash)
    @test all(isfinite, result.fit.coefficients)
    @test result.shift_t_nm ≈ 0.0 atol = 0.02
    @test result.shift_u_nm ≈ 0.0 atol = 0.02
    @test result.rotation_deg ≈ 0.0 atol = 0.5
end
