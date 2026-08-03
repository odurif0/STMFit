#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# test/joint_proxy/test_simulator.jl
# TDD tests for the physics-proxy synthetic chain simulator (Todo 3 of
# joint-proxy-mold-inference-v1).
#
# Run: julia --project=. test/joint_proxy/test_simulator.jl
# ──────────────────────────────────────────────────────────────────────────────

using Test
using Random
using STMSXMIO
using SHA

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "simulator.jl"))
include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "proxy_registry.jl"))
include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "type_posterior.jl"))
using .JointProxySimulator

const REG = JointProxyRegistry
const TP = JointProxyTypePosterior

# ── Shared fixtures ──────────────────────────────────────────────────────────

const DEFAULT_CFG = SimulatorConfig(seed=20260710)
const ENSEMBLE = default_proxy_ensemble()
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _exact_posterior_lobes(case::SyntheticCase)
    templates = REG.render_geometric_templates(joinpath(ROOT, "templates", "chitosan_geometric_sites.tsv");
                                               half_nm=0.32, step_nm=0.08,
                                               parity_flip="t", mirror_flip="u",
                                               normalize="zscore")
    by_key = Dict((t.type, t.parity, t.mirror) => t for t in templates)
    lobes = TP.TypePosteriorLobeEvidence[]
    for i in 1:case.truth.N
        parity = JointProxySimulator._parity_for_lobe(i, case.truth.N,
                                                      case.truth.direction, case.truth.phase)
        templ = by_key[(case.truth.sequence[i], parity, case.truth.mirror)]
        push!(lobes, TP.TypePosteriorLobeEvidence(Dict("fwd" => copy(templ.pixels),
                                                      "bwd" => copy(templ.pixels))))
    end
    reg = REG.ProxyEnsemble([
        REG.ProxyEntry(REG.ProxySource("geometric", "geometric",
                                       joinpath(ROOT, "templates", "chitosan_geometric_sites.tsv"),
                                       "sha_geometric", 0.0, 1.0, true), templates)
    ], 0.32, 0.08, 9, 81, "payload", String[])
    return lobes, reg
end

# ── Tests ────────────────────────────────────────────────────────────────────

@testset "JointProxySimulator — Todo 3" begin

    # ════════════════════════════════════════════════════════════════════════
    # 1. Config validation — invalid ranges throw BEFORE allocation
    # ════════════════════════════════════════════════════════════════════════
    @testset "Config validation throws before allocation" begin
        # Invalid N range
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(n_min=0, n_max=5), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(n_min=5, n_max=3), ENSEMBLE)
        # Impossible image geometry
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(width=0), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(height=0), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(range_nm=(0.0, 4.0)), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(range_nm=(8.0, -1.0)), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(spacing_nm=0.0), ENSEMBLE)
        # Negative noise / blur
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(noise_sigma=-0.1), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(blur_sigma_px=-1.0), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(row_offset_sigma=-0.01), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(affine_scale_jitter=-0.1), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(drift_strength=-0.01), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(contrast_strength=-0.01), ENSEMBLE)
        # correlated_noise_frac out of [0,1]
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(correlated_noise_frac=1.5), ENSEMBLE)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), SimulatorConfig(correlated_noise_frac=-0.1), ENSEMBLE)
        # Empty proxy ensemble
        empty_ens = ProxyEnsemble(ProxySite[], ProxySite[])
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), DEFAULT_CFG, empty_ens)
        # Only type 0
        half_ens = ProxyEnsemble(ENSEMBLE.sites0, ProxySite[])
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), DEFAULT_CFG, half_ens)
        # Only type 1
        half_ens2 = ProxyEnsemble(ProxySite[], ENSEMBLE.sites1)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), DEFAULT_CFG, half_ens2)
        # Unknown control
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), DEFAULT_CFG, ENSEMBLE; control=:unknown_control)
    end

    # ════════════════════════════════════════════════════════════════════════
    # 2. Determinism — same seed = byte-identical, different seeds = different
    # ════════════════════════════════════════════════════════════════════════
    @testset "Determinism" begin
        case1 = generate_case(MersenneTwister(42), DEFAULT_CFG, ENSEMBLE; case_id="det")
        case2 = generate_case(MersenneTwister(42), DEFAULT_CFG, ENSEMBLE; case_id="det")
        @test case_checksum(case1) == case_checksum(case2)
        @test case1.case_seed == case2.case_seed

        case3 = generate_case(MersenneTwister(99), DEFAULT_CFG, ENSEMBLE; case_id="det")
        @test case_checksum(case1) != case_checksum(case3)
    end

    @testset "Replay from recorded seed" begin
        case = generate_case(MersenneTwister(42), DEFAULT_CFG, ENSEMBLE; case_id="replay")
        replayed = replay_case(case.case_seed, DEFAULT_CFG, ENSEMBLE;
                               control=case.control, case_id=case.case_id)
        @test case == replayed
        @test case_checksum(case) == case_checksum(replayed)
        @test image_data_bytes(case) == image_data_bytes(replayed)
    end

    @testset "Recorded view nuisance truth" begin
        case = generate_case(MersenneTwister(123), DEFAULT_CFG, ENSEMBLE; case_id="nuisance")
        for nuisance in (case.truth.fwd_nuisance, case.truth.bwd_nuisance)
            @test all(isfinite, (nuisance.scale_x, nuisance.scale_y, nuisance.shear, nuisance.blur_sigma_px))
            @test (1.0 - DEFAULT_CFG.affine_scale_jitter) <= nuisance.scale_x <= (1.0 + DEFAULT_CFG.affine_scale_jitter)
            @test (1.0 - DEFAULT_CFG.affine_scale_jitter) <= nuisance.scale_y <= (1.0 + DEFAULT_CFG.affine_scale_jitter)
            @test -DEFAULT_CFG.drift_strength <= nuisance.shear <= DEFAULT_CFG.drift_strength
            @test (0.5 * DEFAULT_CFG.blur_sigma_px) <= nuisance.blur_sigma_px <= (1.5 * DEFAULT_CFG.blur_sigma_px)
            xs, ys = JointProxySimulator._inverse_affine_coords(nuisance, 1.25, 2.5)
            @test xs ≈ 1.25 / nuisance.scale_x
            @test ys ≈ (2.5 - nuisance.shear * 1.25) / nuisance.scale_y
        end
        @test case.truth.fwd_nuisance != case.truth.bwd_nuisance
        @test case_checksum(case) == case_checksum(replay_case(case.case_seed, DEFAULT_CFG, ENSEMBLE;
                                                                 control=case.control, case_id=case.case_id))
    end

    @testset "Direction / phase / mirror transforms" begin
        case = nothing
        for seed in 1:200
            candidate = generate_case(MersenneTwister(seed), DEFAULT_CFG, ENSEMBLE;
                                      control=CONTROL_NORMAL, case_id="transform")
            if candidate.truth.N >= 4
                case = candidate
                break
            end
        end
        @test case !== nothing
        @test case.truth.direction in (0, 1)
        @test case.truth.phase in (0, 1)
        @test case.truth.mirror in (0, 1)

        parity_seq = [JointProxySimulator._parity_for_lobe(i, case.truth.N,
                                                           case.truth.direction,
                                                           case.truth.phase)
                      for i in 1:case.truth.N]
        @test all(p in (0, 1) for p in parity_seq)

        for i in 1:case.truth.N
            sites = JointProxySimulator.lobe_proxy_sites(case, ENSEMBLE, i)
            expected = JointProxySimulator._transformed_sites(
                case.truth.sequence[i] == 0 ? ENSEMBLE.sites0 : ENSEMBLE.sites1,
                parity_seq[i], case.truth.mirror)
            @test sites == expected
        end
    end

    # ════════════════════════════════════════════════════════════════════════
    # 3. Sequence length equals sampled N
    # ════════════════════════════════════════════════════════════════════════
    @testset "Sequence length equals N" begin
        for seed in 1:20
            case = generate_case(MersenneTwister(seed), DEFAULT_CFG, ENSEMBLE)
            @test length(case.truth.sequence) == case.truth.N
            @test all(s in (0, 1) for s in case.truth.sequence)
        end
    end

    # ════════════════════════════════════════════════════════════════════════
    # 4. Both identities occur across a fixed batch
    # ════════════════════════════════════════════════════════════════════════
    @testset "Both identities in batch" begin
        batch = generate_batch(MersenneTwister(20260710), DEFAULT_CFG, ENSEMBLE; n_cases=32)
        all_seq = Int[]
        for c in batch
            append!(all_seq, c.truth.sequence)
        end
        @test 0 in all_seq
        @test 1 in all_seq
    end

    # ════════════════════════════════════════════════════════════════════════
    # 5. N is in configured range
    # ════════════════════════════════════════════════════════════════════════
    @testset "N in range" begin
        cfg = SimulatorConfig(n_min=3, n_max=8, seed=1)
        rng = MersenneTwister(1)
        for _ in 1:50
            case = generate_case(rng, cfg, ENSEMBLE)
            @test cfg.n_min <= case.truth.N <= cfg.n_max
        end
    end

    # ════════════════════════════════════════════════════════════════════════
    # 6. fwd/bwd share geometry but differ in nuisance
    # ════════════════════════════════════════════════════════════════════════
    @testset "Shared geometry, different nuisance" begin
        # With zero nuisance, fwd == bwd (shared geometry, no per-view difference)
        zero_nuisance = SimulatorConfig(
            noise_sigma=0.0, blur_sigma_px=0.0,
            affine_scale_jitter=0.0, drift_strength=0.0,
            row_offset_sigma=0.0, seed=42)
        case = generate_case(MersenneTwister(42), zero_nuisance, ENSEMBLE)
        fwd = get_channel(case.img, "Z"; direction="fwd").data
        bwd = get_channel(case.img, "Z"; direction="bwd").data
        @test fwd == bwd

        # With non-zero nuisance, fwd != bwd
        case2 = generate_case(MersenneTwister(42), DEFAULT_CFG, ENSEMBLE)
        fwd2 = get_channel(case2.img, "Z"; direction="fwd").data
        bwd2 = get_channel(case2.img, "Z"; direction="bwd").data
        @test fwd2 != bwd2
    end

    # ════════════════════════════════════════════════════════════════════════
    # 7. Named controls produce their declared condition
    # ════════════════════════════════════════════════════════════════════════
    @testset "Named controls" begin
        @testset "no_molecule" begin
            zero_cfg = SimulatorConfig(noise_sigma=0.0, blur_sigma_px=0.0,
                                       affine_scale_jitter=0.0, drift_strength=0.0,
                                       row_offset_sigma=0.0, seed=1)
            case = generate_case(MersenneTwister(1), zero_cfg, ENSEMBLE;
                                  control=CONTROL_NO_MOLECULE)
            fwd = get_channel(case.img, "Z"; direction="fwd").data
            @test all(iszero, fwd)
            @test case.control == CONTROL_NO_MOLECULE
            @test case.truth.control == CONTROL_NO_MOLECULE

            bg_case = generate_case(MersenneTwister(2), DEFAULT_CFG, ENSEMBLE;
                                    control=CONTROL_NO_MOLECULE)
            bg_fwd = get_channel(bg_case.img, "Z"; direction="fwd").data
            @test !all(iszero, bg_fwd)
            @test all(isfinite, bg_fwd)
        end

        @testset "identical_molds" begin
            # With identical molds and zero nuisance, a case with all type-0
            # sequence and a case with all type-1 sequence should produce
            # the same image (no type-discriminating contrast).
            cfg = SimulatorConfig(noise_sigma=0.0, blur_sigma_px=0.0,
                                   affine_scale_jitter=0.0, drift_strength=0.0,
                                   row_offset_sigma=0.0, seed=7)
            case = generate_case(MersenneTwister(7), cfg, ENSEMBLE;
                                  control=CONTROL_IDENTICAL_MOLDS)
            @test case.control == CONTROL_IDENTICAL_MOLDS
            @test case.truth.control == CONTROL_IDENTICAL_MOLDS
            # Verify the image is non-trivial (backbone present)
            fwd = get_channel(case.img, "Z"; direction="fwd").data
            @test maximum(abs, fwd) > 0
        end

        @testset "swapped_types" begin
            case = generate_case(MersenneTwister(3), DEFAULT_CFG, ENSEMBLE;
                                  control=CONTROL_SWAPPED_TYPES)
            @test case.control == CONTROL_SWAPPED_TYPES
            @test case.truth.control == CONTROL_SWAPPED_TYPES
            # Image should still be finite and non-trivial
            fwd = get_channel(case.img, "Z"; direction="fwd").data
            @test all(isfinite, fwd)
            @test maximum(abs, fwd) > 0
        end

        @testset "missing_bwd" begin
            case = generate_case(MersenneTwister(5), DEFAULT_CFG, ENSEMBLE;
                                  control=CONTROL_MISSING_BWD)
            @test case.control == CONTROL_MISSING_BWD
            @test case.truth.control == CONTROL_MISSING_BWD
            dirs = [c.direction for c in case.img.channels]
            @test "fwd" in dirs
            @test !("bwd" in dirs)
            @test length(case.img.channels) == 1
        end

        @testset "corrupted_view" begin
            case = generate_case(MersenneTwister(9), DEFAULT_CFG, ENSEMBLE;
                                  control=CONTROL_CORRUPTED_VIEW)
            @test case.control == CONTROL_CORRUPTED_VIEW
            @test case.truth.control == CONTROL_CORRUPTED_VIEW
            bwd = get_channel(case.img, "Z"; direction="bwd").data
            @test any(isnan, bwd)
            # fwd should still be finite
            fwd = get_channel(case.img, "Z"; direction="fwd").data
            @test all(isfinite, fwd)
        end
    end

    # ════════════════════════════════════════════════════════════════════════
    # 8. Finite dimensions and ranges
    # ════════════════════════════════════════════════════════════════════════
    @testset "Finite dimensions and ranges" begin
        batch = generate_batch(MersenneTwister(123), DEFAULT_CFG, ENSEMBLE; n_cases=16)
        for c in batch
            @test c.img.width == DEFAULT_CFG.width
            @test c.img.height == DEFAULT_CFG.height
            @test c.img.range_nm == DEFAULT_CFG.range_nm
            # fwd is always present and finite
            fwd = get_channel(c.img, "Z"; direction="fwd").data
            @test size(fwd) == (DEFAULT_CFG.height, DEFAULT_CFG.width)
            if c.control != CONTROL_CORRUPTED_VIEW
                @test all(isfinite, fwd)
            end
        end
    end

    # ════════════════════════════════════════════════════════════════════════
    # 9. 32-case checksum fixture — reproducible
    # ════════════════════════════════════════════════════════════════════════
    @testset "32-case checksum fixture" begin
        batch1 = generate_batch(MersenneTwister(20260710), DEFAULT_CFG, ENSEMBLE; n_cases=32)
        batch2 = generate_batch(MersenneTwister(20260710), DEFAULT_CFG, ENSEMBLE; n_cases=32)
        cs1 = [case_checksum(c) for c in batch1]
        cs2 = [case_checksum(c) for c in batch2]
        @test cs1 == cs2
        # All 32 checksums are unique (different cases produce different images)
        @test length(unique(cs1)) == 32
        # All checksums are 64-char hex strings
        @test all(length(cs) == 64 for cs in cs1)
        @test all(all(c in "0123456789abcdef" for c in cs) for cs in cs1)
    end

    # ════════════════════════════════════════════════════════════════════════
    # 10. Truth isolation — SXMImage has no truth fields
    # ════════════════════════════════════════════════════════════════════════
    @testset "Truth isolation" begin
        case = generate_case(MersenneTwister(1), DEFAULT_CFG, ENSEMBLE)
        # SXMImage fields: filepath, header, width, height, range_nm, offset_nm, channels
        @test :N ∉ propertynames(case.img)
        @test :sequence ∉ propertynames(case.img)
        @test :truth ∉ propertynames(case.img)
        # Channels are plain SXMChannel (name, unit, direction, data)
        for ch in case.img.channels
            @test :N ∉ propertynames(ch)
            @test :sequence ∉ propertynames(ch)
        end
        # Truth lives only in SyntheticCase.truth
        @test case.truth.N > 0
        @test length(case.truth.sequence) == case.truth.N
    end

    # ════════════════════════════════════════════════════════════════════════
    # 11. load_proxy_ensemble from tracked TSV
    # ════════════════════════════════════════════════════════════════════════
    @testset "load_proxy_ensemble" begin
        tsv = joinpath(dirname(dirname(@__DIR__)), "templates",
                        "chitosan_geometric_sites.tsv")
        if isfile(tsv)
            ens = load_proxy_ensemble(tsv)
            @test !isempty(ens.sites0)
            @test !isempty(ens.sites1)
            @test ens == ENSEMBLE
        end
    end

    @testset "Exact-backbone / true-center posterior compatibility" begin
        zero_cfg = SimulatorConfig(noise_sigma=0.0, blur_sigma_px=0.0,
                                   affine_scale_jitter=0.0, drift_strength=0.0,
                                   row_offset_sigma=0.0, seed=11)
        case = generate_case(MersenneTwister(11), zero_cfg, ENSEMBLE;
                             control=CONTROL_NORMAL, case_id="exact_backbone")
        lobes, reg = _exact_posterior_lobes(case)
        res = TP.infer_type_posterior(lobes, reg; rho=0.0)
        @test all(isapprox.(sum(res.lobe_marginals, dims=2), ones(case.truth.N, 1); atol=1e-12))
        @test res.map_sequence == case.truth.sequence
        @test isapprox(sum(last.(res.proxy_family_sensitivity)), 1.0; atol=1e-12)
    end

    # ════════════════════════════════════════════════════════════════════════
    # 12. Adversarial: stale proxy state (ensemble mutated after creation)
    # ════════════════════════════════════════════════════════════════════════
    @testset "Adversarial: stale proxy state" begin
        ens = default_proxy_ensemble()
        # Generate a valid case first
        case = generate_case(MersenneTwister(1), DEFAULT_CFG, ens)
        @test case.truth.N > 0
        # Now mutate: clear sites1 — should throw on next generate
        empty!(ens.sites1)
        @test_throws ArgumentError generate_case(
            MersenneTwister(1), DEFAULT_CFG, ens)
    end

    # ════════════════════════════════════════════════════════════════════════
    # 13. All controls cycle through a batch
    # ════════════════════════════════════════════════════════════════════════
    @testset "All controls in batch" begin
        controls = collect(ALL_CONTROLS)
        batch = generate_batch(MersenneTwister(55), DEFAULT_CFG, ENSEMBLE;
                                n_cases=length(controls), controls=controls)
        @test length(batch) == length(ALL_CONTROLS)
        for (i, expected) in enumerate(ALL_CONTROLS)
            @test batch[i].control == expected
        end
    end
end

println("\n=== JointProxySimulator tests complete ===")
