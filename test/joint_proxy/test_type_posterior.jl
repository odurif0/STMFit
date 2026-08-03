#!/usr/bin/env julia

using Test
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "proxy_registry.jl"))
include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "type_posterior.jl"))

using .JointProxyRegistry: ProxySource, ProxyTemplate, ProxyEntry, ProxyEnsemble,
                           render_geometric_templates
using .JointProxyTypePosterior: TypePosteriorPriors, TypePosteriorResult,
                                TypePosteriorLobeEvidence, effective_view_factor,
                                score_ncc_cost, binary_chain_forward_backward,
                                binary_chain_viterbi, infer_type_posterior

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "type_posterior_test_support.jl"))


@testset "JointProxyTypePosterior" begin
    @testset "NCC scorer / Viterbi baseline" begin
        a = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        @test score_ncc_cost(a, a) ≈ -1.0
        @test score_ncc_cost(a, b) ≈ 0.0
        unary = [0.0 1.0; 2.0 -1.0; -0.5 -0.25]
        seq, val = binary_chain_viterbi(unary)
        @test seq == [1, 0, 1]
        @test isfinite(val)
    end

    priors = TypePosteriorPriors(type=[1.0, 1.0], direction=[1.0, 1.0], phase=[1.0, 1.0],
                                 mirror=[1.0, 1.0], source=nothing)
    ens = _ensemble()

    @testset "exact brute force N=1..6" begin
        for n in 1:6
            lobes, truth, rho = _evidence(n)
            got = infer_type_posterior(lobes, ens; priors=priors, rho=rho)
            logz, probs, seqs, _ = _bruteforce(lobes, ens, priors; rho=rho)
            @test got isa TypePosteriorResult
            @test isapprox(got.log_evidence, logz; atol=1e-10, rtol=0)
            @test isapprox(sum(got.lobe_marginals, dims=2), ones(n, 1); atol=1e-10, rtol=0)
            @test all(0.0 .<= got.lobe_marginals .<= 1.0)
            map_idx = argmax(probs)
            @test got.map_sequence == seqs[map_idx]
        end
    end

    @testset "actual geometric order and shuffled vector remain identical" begin
        lobes, rho = _geometric_evidence(4)
        ordered = _actual_ensemble()
        shuffled = _actual_ensemble(shuffle=true)
        got_ordered = infer_type_posterior(lobes, ordered; priors=priors, rho=rho)
        got_shuffled = infer_type_posterior(lobes, shuffled; priors=priors, rho=rho)
        @test got_ordered.lobe_marginals ≈ got_shuffled.lobe_marginals atol=1e-10 rtol=0
        @test got_ordered.map_sequence == got_shuffled.map_sequence
        @test got_ordered.log_evidence ≈ got_shuffled.log_evidence atol=1e-10 rtol=0
    end

    @testset "metadata swap inverts the physical MAP" begin
        lobes, rho = _geometric_evidence(4)
        ordered = _actual_ensemble()
        swapped = _actual_ensemble(swap_metadata=true)
        got_ordered = infer_type_posterior(lobes, ordered; priors=priors, rho=rho)
        got_swapped = infer_type_posterior(lobes, swapped; priors=priors, rho=rho)
        @test got_swapped.map_sequence == (1 .- got_ordered.map_sequence)
        @test got_ordered.lobe_marginals[:, 1] ≈ got_swapped.lobe_marginals[:, 2] atol=1e-10 rtol=0
    end

    @testset "reverse symmetry" begin
        lobes, rho = _symmetric_evidence(4)
        got = infer_type_posterior(lobes, ens; priors=priors, rho=rho)
        revlobes, _ = _symmetric_evidence(4; reversed=true)
        got_rev = infer_type_posterior(revlobes, ens; priors=priors, rho=rho)
        @test got.lobe_marginals ≈ reverse(got_rev.lobe_marginals; dims=1) atol=1e-10 rtol=0
    end

    @testset "identical molds => 0.5" begin
        idens = ProxyEnsemble([
            ProxyEntry(e.source, [ProxyTemplate(t.type, t.parity, t.mirror, fill(1.0, length(t.pixels)))
                                  for t in e.templates])
            for e in ens.entries
        ], ens.grid_half_nm, ens.grid_step_nm, ens.grid_n, ens.npix, ens.payload_sha256, String[])
        lobes, _, rho = _evidence(3)
        got = infer_type_posterior(lobes, idens; priors=priors, rho=rho)
        @test all(isapprox.(got.lobe_marginals[:, 1], 0.5; atol=1e-10, rtol=0))
        @test all(isapprox.(got.lobe_marginals[:, 2], 0.5; atol=1e-10, rtol=0))
    end

    @testset "invalid priors / NaN / duplicate-view correction" begin
        lobes, _, rho = _evidence(4; duplicate=true)
        @test_throws ArgumentError infer_type_posterior(lobes, ens;
            priors=TypePosteriorPriors(type=[NaN, 1.0], direction=[1.0, 1.0], phase=[1.0, 1.0], mirror=[1.0, 1.0]))
        bad = ProxyEnsemble([
            i == 1 ? ProxyEntry(e.source, [j == 1 ? ProxyTemplate(t.type, t.parity, t.mirror, fill(NaN, length(t.pixels))) : t for (j, t) in enumerate(e.templates)]) : e
            for (i, e) in enumerate(ens.entries)
        ], ens.grid_half_nm, ens.grid_step_nm, ens.grid_n, ens.npix, ens.payload_sha256, String[])
        @test_throws ArgumentError infer_type_posterior(lobes, bad; priors=priors, rho=rho)
        got = infer_type_posterior(lobes, ens; priors=priors, rho=rho)
        logz, probs, _, _ = _bruteforce(lobes, ens, priors; rho=rho)
        @test isapprox(got.log_evidence, logz; atol=1e-10, rtol=0)
        @test got.map_sequence == _bruteforce(lobes, ens, priors; rho=rho)[3][argmax(probs)]
        @test_throws ArgumentError infer_type_posterior(lobes, ens; priors=priors, rho=rho, pair_bond_evidence=[1.0])
    end

    @testset "missing and duplicate template states fail" begin
        lobes, _, rho = _evidence(2)
        missing = ProxyEnsemble([
            ProxyEntry(ens.entries[1].source, ens.entries[1].templates[2:end])
        ], ens.grid_half_nm, ens.grid_step_nm, ens.grid_n, ens.npix, ens.payload_sha256, String[])
        dup = ProxyEnsemble([
            ProxyEntry(ens.entries[1].source, vcat(ens.entries[1].templates, ens.entries[1].templates[1:1]))
        ], ens.grid_half_nm, ens.grid_step_nm, ens.grid_n, ens.npix, ens.payload_sha256, String[])
        @test_throws ArgumentError infer_type_posterior(lobes, missing; priors=priors, rho=rho)
        @test_throws ArgumentError infer_type_posterior(lobes, dup; priors=priors, rho=rho)
    end

    @testset "proxy-family sensitivity normalizes" begin
        lobes, _, rho = _evidence(2)
        got = infer_type_posterior(lobes, ens; priors=priors, rho=rho)
        @test isapprox(sum(last.(got.proxy_family_sensitivity)), 1.0; atol=1e-10, rtol=0)
        @test length(got.global_state_posterior) == 16
    end
end
