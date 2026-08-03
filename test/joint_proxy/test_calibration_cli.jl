#!/usr/bin/env julia

# TDD tests for the deterministic synthetic-only joint-proxy calibration CLI.
#
# Scope guardrails:
# - synthetic registry/simulator/calibration path only
# - no real-data, manifest, benchmark, truth or grading imports
# - only the Todo 7 CLI surface is exercised

using Test
using TOML

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "calibration_cli.jl"))
using .JointProxyCalibrationCLI

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

function _read_toml(path::String)
    return TOML.parsefile(path)
end

@testset "calibrate_joint_proxy_molds CLI" begin
    @testset "red-first: helper module is callable via include" begin
        @test isdefined(JointProxyCalibrationCLI, :parse_cli)
        @test isdefined(JointProxyCalibrationCLI, :run_cli)
        @test isdefined(JointProxyCalibrationCLI, :GaussianFit2D)
    end

    @testset "happy path emits deterministic TOML without forbidden fields" begin
        mktempdir() do dir
            out = joinpath(dir, "joint_proxy_calibration.toml")
            args = ["--config", CONFIG, "--seed", "20260710", "--cases", "12",
                    "--n-min", "3", "--n-max", "8", "--fast", "--out", out]
            JointProxyCalibrationCLI.run_cli(args)
            bytes_a = read(out)
            parsed = _read_toml(out)
            @test haskey(parsed, "provenance")
            @test haskey(parsed, "simulator")
            @test haskey(parsed, "splits")
            @test haskey(parsed, "count")
            @test haskey(parsed, "type")
            @test parsed["provenance"]["observation_mode"] == "lightweight_fast"
            @test parsed["diagnostics"]["candidate_report_count"] == 12
            @test !haskey(parsed, "truth")
            @test !haskey(parsed, "benchmark")
            @test !haskey(parsed, "expected_N")
            @test !haskey(parsed, "target_N")

            JointProxyCalibrationCLI.run_cli(args)
            bytes_b = read(out)
            @test bytes_a == bytes_b

            @test isfinite(parsed["count"]["temperature"])
            @test isfinite(parsed["type"]["temperature"])
            @test parsed["splits"]["cases"] == 12
        end
    end

    @testset "non-fast invokes adapter and can be unit-tested with injected fitter/view builders" begin
        mktempdir() do dir
            out = joinpath(dir, "joint_proxy_calibration.toml")
            called = Ref(false)
            fake = function (cases, opt, bundle)
                called[] = true
                return JointProxyCalibrationCLI._fast_count_observations(cases, opt, bundle)
            end
            args = ["--config", CONFIG, "--seed", "20260710", "--cases", "12",
                    "--n-min", "2", "--n-max", "4", "--out", out]
            JointProxyCalibrationCLI.run_cli(args; count_adapter=fake)
            @test called[]
            parsed = _read_toml(out)
            @test parsed["provenance"]["observation_mode"] == "lightweight_fast"

            fit_calls = Ref(0)
            view_calls = Ref(0)
            extract_calls = Ref(0)
            fake_fit = function (img, pcfg, ccfg)
                fit_calls[] += 1
                return (Any[], nothing, (xs=[0.0, 1.0], ys=[0.0, 1.0], mask=trues(2, 2), axisctx=(axis=[1.0, 0.0], perp=[0.0, 1.0])))
            end
            fake_view = function (img, channel, pcfg; direction="fwd")
                view_calls[] += 1
                return JointProxyCalibrationCLI.JointProxyCandidateViews.ViewData(
                    direction, [0.0, 1.0], [0.0, 1.0], [1.0 2.0; 3.0 4.0], 0.1, true)
            end
            fake_extract = function (results, ctx, ccfg,
                    views::AbstractVector{JointProxyCalibrationCLI.JointProxyCandidateViews.ViewData};
                    patch_half_nm, patch_step_nm)
                extract_calls[] += 1
                return (; candidates=[(; n=3, valid=true, joint_gcv=1.0)], bwd_missing=false, audit=["fake"])
            end
            opt = JointProxyCalibrationCLI.parse_cli(["--config", CONFIG, "--seed", "20260710",
                "--cases", "12", "--n-min", "2", "--n-max", "4", "--out", out])
            @test isinf(JointProxyCalibrationCLI._real_chain_config(opt).residual_peak_snr_threshold)
            bundle = JointProxyCalibrationCLI._load_bundle(opt)
            case = JointProxyCalibrationCLI.generate_case(JointProxyCalibrationCLI.MersenneTwister(1),
                JointProxyCalibrationCLI.SimulatorConfig(n_min=2, n_max=3, seed=1),
                JointProxyCalibrationCLI.default_proxy_ensemble(); case_id="smoke")
            obs = JointProxyCalibrationCLI._real_count_observations([case], opt, bundle;
                fit_runner=fake_fit, view_builder=fake_view, extractor=fake_extract,
                load_deps=false)
            @test fit_calls[] == 1
            @test view_calls[] >= 1
            @test extract_calls[] == 1
            @test length(obs.count_cases) == 1
            @test obs.count_cases[1].report.candidates[1].joint_gcv == 1.0

            missing_true = (; candidates=[(; n=case.truth.N + 1, valid=true, joint_gcv=1.0)],
                bwd_missing=false, audit=["missing true candidate"])
            filtered = JointProxyCalibrationCLI._count_cases_from_reports(
                [case], [missing_true], opt, bundle, 1)
            @test isempty(filtered)
        end
    end

    @testset "missing true candidates remain explicit diagnostics" begin
        mktempdir() do dir
            out = joinpath(dir, "joint_proxy_calibration.toml")
            args = ["--config", CONFIG, "--seed", "20260710", "--cases", "12",
                    "--n-min", "2", "--n-max", "4", "--out", out]
            filtered_adapter = function (cases, opt, bundle)
                observations = JointProxyCalibrationCLI._fast_count_observations(cases, opt, bundle)
                reports = copy(observations.candidate_reports)
                report = reports[1]
                candidates = filter(candidate -> candidate.n != cases[1].truth.N, report.candidates)
                reports[1] = JointProxyCalibrationCLI.SyntheticCountReport(
                    candidates, report.view_count, report.bwd_missing, report.residual_corr)
                return (
                    count_cases=JointProxyCalibrationCLI._count_cases_from_reports(
                        cases, reports, opt, bundle, opt.cases ÷ 2),
                    observation_mode="filtered_test_fixture",
                    candidate_reports=reports,
                )
            end

            JointProxyCalibrationCLI.run_cli(args; count_adapter=filtered_adapter)
            parsed = _read_toml(out)
            @test parsed["splits"]["count_calibration"] == 5
            @test parsed["splits"]["count_heldout"] == 6
            @test parsed["diagnostics"]["count_candidate_recovered"] == 11
            @test parsed["diagnostics"]["count_candidate_missing_true"] == 1
        end
    end

    @testset "failure modes are explicit" begin
        mktempdir() do dir
            out = joinpath(dir, "joint_proxy_calibration.toml")
            @test_throws Exception JointProxyCalibrationCLI.run_cli([
                "--config", CONFIG, "--seed", "20260710", "--cases", "3",
                "--n-min", "3", "--n-max", "8", "--fast", "--out", out])

            mkpath(joinpath(dir, "locked"))
            chmod(joinpath(dir, "locked"), 0o500)
            @test_throws Exception JointProxyCalibrationCLI.run_cli([
                "--config", CONFIG, "--seed", "20260710", "--cases", "12",
                "--n-min", "3", "--n-max", "8", "--fast",
                "--out", joinpath(dir, "locked", "out.toml")])

            @test_throws Exception JointProxyCalibrationCLI.run_cli([
                "--config", CONFIG, "--seed", "20260710", "--cases", "12",
                "--n-min", "9", "--n-max", "8", "--fast", "--out", out])

            @test_throws Exception JointProxyCalibrationCLI.run_cli([
                "--config", CONFIG, "--seed", "20260710", "--cases", "12",
                "--n-min", "3", "--n-max", "8", "--fast", "--out", out,
                "--benchmark"])
        end
    end

    @testset "hash validation rejects altered provenance" begin
        mktempdir() do dir
            cfg = joinpath(dir, "cfg.toml")
            cp(CONFIG, cfg; force=true)
            txt = read(cfg, String) * "\n[provenance]\nexpected_proxy_payload_sha256 = \"deadbeef\"\n"
            open(cfg, "w") do io
                write(io, txt)
            end
            @test_throws Exception JointProxyCalibrationCLI.run_cli([
                "--config", cfg, "--seed", "20260710", "--cases", "12",
                "--n-min", "3", "--n-max", "8", "--fast", "--out", joinpath(dir, "out.toml")])
        end
    end

    @testset "config parsing rejects forbidden truth/benchmark fields" begin
        mktempdir() do dir
            cfg = joinpath(dir, "bad.toml")
            write(cfg, "[model]\nfoo = 1\n[benchmark]\ntruth = 1\n")
            @test_throws Exception JointProxyCalibrationCLI.run_cli([
                "--config", cfg, "--seed", "20260710", "--cases", "12",
                "--n-min", "3", "--n-max", "8", "--fast", "--out", joinpath(dir, "out.toml")])
        end
    end
end
