#!/usr/bin/env julia

# Synthetic end-to-end gate for common-only registration followed by frozen
# type-contrast scoring. Diagnostic-only; no benchmark labels or expected
# sequence enter registration/scoring.

using Test
using Random

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "whole_roi_end_to_end.jl"))

const E2E_REG = Main.JointProxyRegistry
const E2E_DC = Main.WholeRoiDiagnosticConfig
const E2E_HALF = 0.4
const E2E_STEP = 0.2
const E2E_XS = collect(-1.0:0.1:1.0)
const E2E_YS = collect(-0.6:0.1:0.6)
const E2E_CENTERS = [(-0.45, 0.0), (-0.15, 0.0), (0.15, 0.0), (0.45, 0.0)]
const E2E_SEQUENCE = Int[1, 0, 1, 0]

function e2e_pixels(matrix::Matrix{Float64})
    return vec(permutedims(matrix))
end

function e2e_entry(M0::Matrix{Float64}, M1::Matrix{Float64})
    templates = E2E_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        pixels = typ == 0 ? e2e_pixels(M0) : e2e_pixels(M1)
        push!(templates, E2E_REG.ProxyTemplate(typ, parity, mirror, copy(pixels)))
    end
    source = E2E_REG.ProxySource(
        "synthetic", "whole_roi_end_to_end", "", "", 0.0, 1.0, true)
    return E2E_REG.ProxyEntry(source, templates)
end

function e2e_matrices()
    M0 = zeros(5, 5)
    M1 = zeros(5, 5)
    M0[3, 2] = 1.0
    M0[2, 3] = 0.35
    M1[3, 4] = 1.2
    M1[4, 3] = 0.20
    return M0, M1
end

function e2e_backbone()
    return [sum(exp(-((x - cx)^2 / 0.22^2 + (y - cy)^2 / 0.15^2) / 2)
                for (cx, cy) in E2E_CENTERS)
            for y in E2E_YS, x in E2E_XS]
end

function e2e_config(; shift=0.0, rotation=0.0, blur=0.0)
    range(limit, step) = E2E_DC.DiagnosticRange(-limit, limit, step)
    positive_range(limit, step) = E2E_DC.DiagnosticRange(0.0, limit, step)
    shift_coarse = range(shift, shift == 0 ? 0.01 : shift)
    shift_fine = range(shift == 0 ? 0.0 : shift / 2, shift == 0 ? 0.01 : shift / 2)
    rot_coarse = range(rotation, rotation == 0 ? 0.1 : rotation)
    rot_fine = range(rotation == 0 ? 0.0 : rotation / 2,
                     rotation == 0 ? 0.1 : rotation / 2)
    blur_coarse = positive_range(blur, blur == 0 ? 0.01 : blur)
    blur_fine = positive_range(blur == 0 ? 0.0 : blur / 2,
                               blur == 0 ? 0.01 : blur / 2)
    return E2E_DC.DiagnosticConfig(
        E2E_DC.DiagnosticParamBounds(shift_coarse, shift_fine),
        E2E_DC.DiagnosticParamBounds(shift_coarse, shift_fine),
        E2E_DC.DiagnosticParamBounds(rot_coarse, rot_fine),
        E2E_DC.DiagnosticParamBounds(blur_coarse, blur_fine),
        1e-6, 1e-9)
end

function e2e_case(set; sequence=E2E_SEQUENCE, shift_t=0.0, shift_u=0.0,
                  rotation_deg=0.0, blur_sigma=0.0, common_amp=1.0,
                  contrast_amp=0.3, noise_sigma=0.0, seed=20260725)
    return synthesize_whole_roi_observation(set, sequence, e2e_backbone(),
        E2E_CENTERS, E2E_XS, E2E_YS;
        half_nm=E2E_HALF, step_nm=E2E_STEP,
        shift_t_nm=shift_t, shift_u_nm=shift_u,
        rotation_deg=rotation_deg, blur_sigma_nm=blur_sigma,
        background=0.25, backbone_amp=0.7, common_amp=common_amp,
        contrast_amp=contrast_amp, noise_sigma=noise_sigma, seed=seed)
end

function e2e_run(case, set, config)
    return run_whole_roi_common_contrast(case, set;
        config=config, half_nm=E2E_HALF, step_nm=E2E_STEP)
end

@testset "truth-withheld exact and common-only controls" begin
    M0, M1 = e2e_matrices()
    set = derive_common_contrast(e2e_entry(M0, M1))
    exact_case = e2e_case(set)
    @test !any(name -> name in (:sequence, :truth, :direction, :phase, :mirror),
               propertynames(exact_case))
    exact = e2e_run(exact_case, set, e2e_config())
    @test !exact.scoring.abstain
    @test exact.scoring.best.sequence == E2E_SEQUENCE
    @test exact.registration.direction == 0
    @test exact.registration.phase == 0
    @test exact.registration.mirror == 0

    common_only = e2e_run(e2e_case(set; contrast_amp=0.0), set, e2e_config())
    @test common_only.scoring.abstain
    @test "best_runner_tie" in common_only.scoring.abstention_reasons
    @test "complement_tie" in common_only.scoring.abstention_reasons
end

@testset "metadata swap and exact complement tie" begin
    M0, M1 = e2e_matrices()
    set = derive_common_contrast(e2e_entry(M0, M1))
    swapped = derive_common_contrast(e2e_entry(M1, M0))
    case = e2e_case(set)
    nominal = e2e_run(case, set, e2e_config())
    inverted = e2e_run(case, swapped, e2e_config())
    @test freeze_registration_tsv(nominal.registration) ==
          freeze_registration_tsv(inverted.registration)
    @test inverted.scoring.best.sequence == 1 .- nominal.scoring.best.sequence

    identical = derive_common_contrast(e2e_entry(M0, M0))
    tied = e2e_run(e2e_case(identical), identical, e2e_config())
    @test tied.scoring.abstain
    @test "zero_contrast_gain" in tied.scoring.abstention_reasons
    @test "complement_tie" in tied.scoring.abstention_reasons
end

@testset "frozen nuisance, noise, boundary, and transpose controls" begin
    M0, M1 = e2e_matrices()
    set = derive_common_contrast(e2e_entry(M0, M1))
    config = e2e_config(shift=0.04, rotation=0.5, blur=0.02)
    nuisance = e2e_case(set; shift_t=0.02, shift_u=-0.02,
        rotation_deg=0.25, blur_sigma=0.01, contrast_amp=0.15)
    recovered = e2e_run(nuisance, set, config)
    @test abs(recovered.registration.shift_t_nm - 0.02) <= 0.02 + 1e-12
    @test abs(recovered.registration.shift_u_nm + 0.02) <= 0.02 + 1e-12
    @test abs(recovered.registration.rotation_deg - 0.25) <= 0.25 + 1e-12
    @test abs(recovered.registration.blur_sigma_nm - 0.01) <= 0.01 + 1e-12
    @test recovered.scoring.best.sequence == E2E_SEQUENCE

    clean = e2e_run(e2e_case(set), set, e2e_config())
    noisy_runs = [e2e_run(e2e_case(set; noise_sigma=0.02, seed=seed),
                          set, e2e_config()) for seed in 1:4]
    @test all(run -> isfinite(run.scoring.best.sse), noisy_runs)
    @test all(run -> all(score -> isfinite(score.sse), run.scoring.ranking),
              noisy_runs)
    @test all(run -> run.scoring.best_runner_margin <=
                     clean.scoring.best_runner_margin, noisy_runs)

    outside = e2e_run(e2e_case(set; shift_t=0.08, contrast_amp=0.1), set, config)
    @test :shift_t_nm in outside.boundary_parameters

    transposed = derive_common_contrast(e2e_entry(permutedims(M0), permutedims(M1)))
    orientation_case = e2e_case(set; contrast_amp=0.0)
    bad = e2e_run(orientation_case, transposed, e2e_config())
    good = e2e_run(orientation_case, set, e2e_config())
    require_exact_common(result) = result.registration.fit.sse <= 1e-8 ||
        throw(AssertionError("transposed asymmetric template failed exact common fit"))
    @test require_exact_common(good)
    @test_throws AssertionError require_exact_common(bad)
end
