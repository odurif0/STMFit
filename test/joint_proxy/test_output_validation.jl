#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "validate_joint_proxy_predictions.jl"))
using .JointProxyOutputValidation

include(joinpath(@__DIR__, "..", "merge_joint_proxy_shards.jl"))
using .JointProxyShardMerge

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "output_validation_fixtures.jl"))
using .JointProxyOutputValidationFixtures

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

_err(f) = try f(); nothing catch e; e end

@testset "joint proxy output validation" begin
    @testset "red-first: validation and merge entrypoints load" begin
        @test isdefined(JointProxyOutputValidation, :run_cli)
        @test isdefined(JointProxyShardMerge, :run_cli)
    end

    @testset "merge sorts two-digit candidate and lobe indices numerically" begin
        candidate_rows = Any[Dict("file" => "a.sxm", "n" => "10"), Dict("file" => "a.sxm", "n" => "2")]
        JointProxyShardMerge._stable_sort!(candidate_rows, ["file", "n"])
        @test [row["n"] for row in candidate_rows] == ["2", "10"]

        lobe_rows = Any[Dict("file" => "a.sxm", "n" => "2", "lobe" => "10"),
            Dict("file" => "a.sxm", "n" => "2", "lobe" => "2")]
        JointProxyShardMerge._stable_sort!(lobe_rows, ["file", "n", "lobe"])
        @test [row["lobe"] for row in lobe_rows] == ["2", "10"]
    end

    @testset "summary-only abstention is valid when no count candidate survives" begin
        summary = Dict(
            "file" => "bad.sxm", "N_prediction" => "?", "N_probability" => "0.0",
            "N_confidence" => "0.0", "N_abstained" => "true", "candidate_count" => "0",
            "lobe_count" => "0", "view_count" => "2", "bwd_missing" => "false",
        )
        @test JointProxyOutputValidation._validate_summaries(
            [summary], "chain_summary.tsv", Dict(), Dict(), Dict()) === nothing
        invalid = copy(summary); invalid["N_abstained"] = "false"
        @test_throws JointProxyOutputValidation.ValidationError JointProxyOutputValidation._validate_summaries(
            [invalid], "chain_summary.tsv", Dict(), Dict(), Dict())
    end

    @testset "happy path validates and merge is byte-identical to unchunked fixture" begin
        mktempdir() do dir
            cfg_hash = config_hash(CONFIG)
            expected = write_fixture(joinpath(dir, "expected"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, chunk="none", data_dir=joinpath(dir, "data"))
            shard1 = write_fixture(joinpath(dir, "shard1"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, files=["a.sxm"], chunk="1/3", data_dir=joinpath(dir, "data"))
            shard2 = write_fixture(joinpath(dir, "shard2"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, files=["b.sxm"], chunk="2/3", data_dir=joinpath(dir, "data"))
            shard3 = write_fixture(joinpath(dir, "shard3"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, files=["c.sxm"], chunk="3/3", data_dir=joinpath(dir, "data"))
            @test JointProxyOutputValidation.run_cli(["--artifacts", expected, "--calibration", joinpath(expected, "joint_proxy_calibration.toml")]) isa ValidationReport
            out = joinpath(dir, "merged")
            @test JointProxyShardMerge.run_cli(["--shard", shard1, "--shard", shard2, "--shard", shard3, "--outdir", out]) == out
            for name in ("candidate_n.tsv", "candidate_lobes.tsv", "predictions.tsv", "chain_summary.tsv", "run_manifest.toml")
                @test read_bytes(joinpath(out, name)) == read_bytes(joinpath(expected, name))
            end
        end
    end

    @testset "validator rejects corruptions with named diagnostics" begin
        mktempdir() do dir
            cfg_hash = config_hash(CONFIG)
            base = write_fixture(joinpath(dir, "base"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, data_dir=joinpath(dir, "data"))

            dup = copy_fixture(base, joinpath(dir, "dup"))
            write_tsv(joinpath(dup, "predictions.tsv"), JointProxyOutputValidation.PREDICTION_FIELDS, vcat(fixture_rows()[3], fixture_rows()[3][1]))
            @test occursin("duplicate key", sprint(showerror, _err(() -> JointProxyOutputValidation.run_cli(["--artifacts", dup, "--calibration", joinpath(dup, "joint_proxy_calibration.toml")]))))

            drift = copy_fixture(base, joinpath(dir, "drift"))
            rows = collect(fixture_rows()[3]); rows[1] = (; rows[1]..., N_probability=0.2); rows[2] = (; rows[2]..., N_probability=0.2)
            write_tsv(joinpath(drift, "predictions.tsv"), JointProxyOutputValidation.PREDICTION_FIELDS, rows)
            @test occursin("probability drift", sprint(showerror, _err(() -> JointProxyOutputValidation.run_cli(["--artifacts", drift, "--calibration", joinpath(drift, "joint_proxy_calibration.toml")]))))

            missing_lobe = copy_fixture(base, joinpath(dir, "missing_lobe"))
            cand_lobes = collect(fixture_rows()[2]); filter!(r -> !(r.file == "b.sxm" && r.n == 4 && r.lobe == 4), cand_lobes)
            write_tsv(joinpath(missing_lobe, "candidate_lobes.tsv"), JointProxyOutputValidation.CANDIDATE_LOBE_FIELDS, cand_lobes)
            @test occursin("missing candidate lobe", sprint(showerror, _err(() -> JointProxyOutputValidation.run_cli(["--artifacts", missing_lobe, "--calibration", joinpath(missing_lobe, "joint_proxy_calibration.toml")]))))

            forced = copy_fixture(base, joinpath(dir, "forced"))
            preds = collect(fixture_rows()[3]); preds[4] = (; preds[4]..., predicted=0)
            write_tsv(joinpath(forced, "predictions.tsv"), JointProxyOutputValidation.PREDICTION_FIELDS, preds)
            @test occursin("forced final label under uncertain N", sprint(showerror, _err(() -> JointProxyOutputValidation.run_cli(["--artifacts", forced, "--calibration", joinpath(forced, "joint_proxy_calibration.toml")]))))

            forbidden = write_fixture(joinpath(dir, "forbidden"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, forbidden_column="truth")
            @test occursin("forbidden expected_N/sequence/truth column", sprint(showerror, _err(() -> JointProxyOutputValidation.run_cli(["--artifacts", forbidden, "--calibration", joinpath(forbidden, "joint_proxy_calibration.toml")]))))

            bad_src = write_fixture(joinpath(dir, "bad_src"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, bad_source_hash=true)
            @test occursin("source hash mismatch", sprint(showerror, _err(() -> JointProxyOutputValidation.run_cli(["--artifacts", bad_src, "--calibration", joinpath(bad_src, "joint_proxy_calibration.toml")]))))
        end
    end

    @testset "merge rejects overlap by default" begin
        mktempdir() do dir
            cfg_hash = config_hash(CONFIG)
            shard1 = write_fixture(joinpath(dir, "shard1"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, files=["a.sxm"], chunk="1/2", data_dir=joinpath(dir, "data"))
            shard2 = write_fixture(joinpath(dir, "shard2"); config_hash=cfg_hash, source_hash=SOURCE_HASH,
                payload_hash=PAYLOAD_HASH, config_path=CONFIG, files=["a.sxm"], chunk="2/2", data_dir=joinpath(dir, "data"))
            @test occursin("overlapping shard", sprint(showerror, _err(() -> JointProxyShardMerge.run_cli(["--shard", shard1, "--shard", shard2, "--outdir", joinpath(dir, "out")]))))
        end
    end
end
