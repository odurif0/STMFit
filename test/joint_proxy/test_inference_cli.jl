#!/usr/bin/env julia

using Test
using SHA
using TOML

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "inference_cli.jl"))
using .JointProxyInferenceCLI
using STMSXMIO

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG = joinpath(ROOT, "config", "joint_proxy_molds.toml")

function _write(path::AbstractString, text::AbstractString)
    open(path, "w") do io
        write(io, text)
    end
end

function _make_calibration(path::AbstractString; config_hash::AbstractString, source_hash::AbstractString,
                           payload_hash::AbstractString, count_temp::Real=1.0,
                           count_threshold::Real=0.6, type_temp::Real=1.0,
                           type_threshold::Real=0.5)
    open(path, "w") do io
        print(io, "[provenance]\n")
        print(io, "config_sha256 = \"$config_hash\"\n")
        print(io, "source_sha256 = \"$source_hash\"\n")
        print(io, "payload_sha256 = \"$payload_hash\"\n\n")
        print(io, "[count]\n")
        print(io, "temperature = $count_temp\n")
        print(io, "confidence_threshold = $count_threshold\n\n")
        print(io, "[type]\n")
        print(io, "temperature = $type_temp\n")
        print(io, "confidence_threshold = $type_threshold\n")
    end
end

function _fake_artifacts(opts, bundle, files)
    sorted = sort(files)
    candidate_n_rows = NamedTuple[]
    candidate_lobes_rows = NamedTuple[]
    predictions_rows = NamedTuple[]
    chain_summary_rows = NamedTuple[]
    for file in sorted
        idx = 1 + mod(sum(Int(c) for c in codeunits(file)), 3)
        hi_prob = idx == 1 ? 0.65 : idx == 2 ? 0.75 : 0.85
        lo_prob = 1.0 - hi_prob
        push!(candidate_n_rows, (; file=file, n=3, joint_gcv=2.0 + idx, delta_rel=0.40, probability=lo_prob, rank=2,
            is_map=false, view_count=2, bwd_missing=false, residual_corr=0.2, effective_view_factor=0.83,
            effective_n_views=1.66, p_eff=8, nd_joint=16, source_gcv=2.8 + idx, joint_nrmse=0.20))
        push!(candidate_n_rows, (; file=file, n=2, joint_gcv=1.0 + idx, delta_rel=0.0, probability=hi_prob, rank=1,
            is_map=true, view_count=2, bwd_missing=false, residual_corr=0.1, effective_view_factor=0.91,
            effective_n_views=1.82, p_eff=6, nd_joint=12, source_gcv=1.5 + idx, joint_nrmse=0.10))
        for n in (2, 3)
            for lobe in 1:2
                p1 = n == 2 ? (lobe == 1 ? 0.20 : 0.70) : (lobe == 1 ? 0.35 : 0.75)
                pred = p1 >= 0.5 ? 1 : 0
                push!(candidate_lobes_rows, (; file=file, n=n, lobe=lobe, x_nm=0.5 * lobe, y_nm=1.0 + lobe,
                    t_nm=0.1 * lobe, u_nm=0.2 * lobe, amplitude=0.4 + 0.1 * lobe, sigma_parallel_nm=0.25,
                    sigma_perp_nm=0.15, skew_ratio=1.0 + 0.1 * (n - 2), type_probability_0=1.0 - p1,
                    type_probability_1=p1, predicted_type=opts.uncalibrated ? "?" : pred, confidence=max(1.0 - p1, p1)))
            end
        end
        for lobe in 1:2
            p1 = lobe == 1 ? 0.20 : 0.70
            pred = opts.uncalibrated ? "?" : (p1 >= 0.5 ? 1 : 0)
            push!(predictions_rows, (; file=file, lobe=lobe, n_prediction=opts.uncalibrated ? "?" : 2,
                n_probability=hi_prob, n_confidence=hi_prob, type_probability_0=1.0 - p1, type_probability_1=p1,
                predicted=pred, confidence=max(1.0 - p1, p1)))
        end
        push!(chain_summary_rows, (; file=file, candidate_count=2, lobe_count=4, n_prediction=opts.uncalibrated ? "?" : 2,
            n_probability=hi_prob, n_confidence=hi_prob, n_abstained=opts.uncalibrated, view_count=2, bwd_missing=false))
    end
    manifest = (; provenance=(; config_sha256=bundle.config_hash, source_sha256=bundle.source_hash,
                              payload_sha256=bundle.payload_hash, uncalibrated=opts.uncalibrated),
                 inputs=(; config=opts.config, data_dir=opts.data_dir, files=sorted,
                    chunk=opts.chunk === nothing ? "none" : "$(opts.chunk[1])/$(opts.chunk[2])"),
                 outputs=(; rows=length(predictions_rows)))
    return JointProxyInferenceCLI.InferenceArtifacts(candidate_n_rows, candidate_lobes_rows,
        predictions_rows, chain_summary_rows, manifest)
end

@testset "infer_joint_proxy_molds CLI" begin
    @test isdefined(JointProxyInferenceCLI, :parse_cli)
    @test isdefined(JointProxyInferenceCLI, :run_cli)
    @test isdefined(JointProxyInferenceCLI, :STMSXMIO)
    @test isdefined(JointProxyInferenceCLI, :GaussianFit2D)
    calibration_sources = JointProxyInferenceCLI._calibration_source_files()
    @test any(endswith(path, "calibration_cli_adapter.jl") for path in calibration_sources)
    @test all(isfile, calibration_sources)
    @test length(unique(calibration_sources)) == length(calibration_sources)

    @testset "parses valid file selection and chunking" begin
        opt = JointProxyInferenceCLI.parse_cli(["--config", CONFIG, "--data-dir", "data",
            "--files", "b.sxm,a.sxm", "--calibration", "cal.toml", "--outdir", "out",
            "--chunk", "2/3"])
        @test opt.files == ["b.sxm", "a.sxm"]
        @test opt.chunk == (2, 3)
        @test !opt.uncalibrated
    end

    @testset "missing option values raise typed CLI errors" begin
        for flag in ("--config", "--data-dir", "--files", "--files-from", "--calibration", "--outdir", "--chunk")
            @test_throws JointProxyInferenceCLI.InferenceCliError JointProxyInferenceCLI.parse_cli([flag])
        end
    end

    @testset "rejects malicious or duplicate file lists before loading data" begin
        mktempdir() do dir
            list = joinpath(dir, "files.txt")
            _write(list, "# comment\nfoo.sxm\tbar\n")
            saw_load = Ref(false)
            @test_throws Exception JointProxyInferenceCLI.run_cli([
                "--config", CONFIG, "--data-dir", dir, "--files-from", list,
                "--uncalibrated", "--outdir", joinpath(dir, "out")],
                adapter=(opts, bundle, files) -> (saw_load[] = true; _fake_artifacts(opts, bundle, files)))
            @test !saw_load[]
        end
    end

    @testset "writes deterministic artifacts and stable ordering" begin
        mktempdir() do dir
            data_dir = joinpath(dir, "data")
            mkpath(data_dir)
            for f in ("a.sxm", "b.sxm")
                _write(joinpath(data_dir, f), "stub")
            end
            cfg_hash = bytes2hex(sha256(read(CONFIG)))
            payload_hash = JointProxyRegistry.load_registry(CONFIG).payload_sha256
            source_hash = JointProxyInferenceCLI._calibration_source_sha256()
            cal = joinpath(dir, "calibration.toml")
            _make_calibration(cal; config_hash=cfg_hash, source_hash=source_hash, payload_hash=payload_hash)
            out = joinpath(dir, "out")
            JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "b.sxm,a.sxm", "--calibration", cal, "--outdir", out,
                "--chunk", "1/1"], adapter=_fake_artifacts)

            @test isfile(joinpath(out, "candidate_n.tsv"))
            @test isfile(joinpath(out, "candidate_lobes.tsv"))
            @test isfile(joinpath(out, "predictions.tsv"))
            @test isfile(joinpath(out, "chain_summary.tsv"))
            @test isfile(joinpath(out, "run_manifest.toml"))

            cand = JointProxyInferenceCLI._read_tsv(joinpath(out, "candidate_n.tsv"))
            @test [r["file"] for r in cand] == ["a.sxm", "a.sxm", "b.sxm", "b.sxm"]
            @test isapprox(sum(parse.(Float64, [r["probability"] for r in cand if r["file"] == "a.sxm"])), 1.0; atol=1e-12)
            lobes = JointProxyInferenceCLI._read_tsv(joinpath(out, "candidate_lobes.tsv"))
            @test length(unique((r["file"], r["n"], r["lobe"]) for r in lobes)) == length(lobes)
            pred = JointProxyInferenceCLI._read_tsv(joinpath(out, "predictions.tsv"))
            @test all(r["N_prediction"] == "2" for r in pred)
            mani = TOML.parsefile(joinpath(out, "run_manifest.toml"))
            @test mani["provenance"]["config_sha256"] == cfg_hash
            @test mani["provenance"]["payload_sha256"] == payload_hash
        end
    end

    @testset "uncalibrated mode keeps probabilities but forces hard labels to question marks" begin
        mktempdir() do dir
            data_dir = joinpath(dir, "data")
            mkpath(data_dir)
            _write(joinpath(data_dir, "a.sxm"), "stub")
            out = joinpath(dir, "uncalibrated")
            JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm", "--uncalibrated", "--outdir", out], adapter=_fake_artifacts)
            pred = JointProxyInferenceCLI._read_tsv(joinpath(out, "predictions.tsv"))
            @test all(r["predicted"] == "?" for r in pred)
            @test all(parse(Float64, r["N_probability"]) > 0 for r in pred)
        end
    end

    @testset "chunked output unions match unchunked after sort" begin
        mktempdir() do dir
            data_dir = joinpath(dir, "data")
            mkpath(data_dir)
            for f in ("a.sxm", "b.sxm", "c.sxm")
                _write(joinpath(data_dir, f), "stub")
            end
            cfg_hash = bytes2hex(sha256(read(CONFIG)))
            payload_hash = JointProxyRegistry.load_registry(CONFIG).payload_sha256
            source_hash = JointProxyInferenceCLI._calibration_source_sha256()
            cal = joinpath(dir, "calibration.toml")
            _make_calibration(cal; config_hash=cfg_hash, source_hash=source_hash, payload_hash=payload_hash)
            out_all = joinpath(dir, "all")
            out_1 = joinpath(dir, "c1")
            out_2 = joinpath(dir, "c2")
            JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm,b.sxm,c.sxm", "--calibration", cal, "--outdir", out_all], adapter=_fake_artifacts)
            JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm,b.sxm,c.sxm", "--calibration", cal, "--outdir", out_1,
                "--chunk", "1/2"], adapter=_fake_artifacts)
            JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm,b.sxm,c.sxm", "--calibration", cal, "--outdir", out_2,
                "--chunk", "2/2"], adapter=_fake_artifacts)
            all_rows = JointProxyInferenceCLI._read_tsv(joinpath(out_all, "candidate_n.tsv"))
            chunk_rows = vcat(JointProxyInferenceCLI._read_tsv(joinpath(out_1, "candidate_n.tsv")),
                              JointProxyInferenceCLI._read_tsv(joinpath(out_2, "candidate_n.tsv")))
            sort_key(r) = (r["file"], parse(Int, r["n"]))
            @test sort(all_rows; by=sort_key) == sort(chunk_rows; by=sort_key)
        end
    end

    @testset "real adapter dispatches fit, view, extract, count and type hooks" begin
        mktempdir() do dir
            img = Main.STMSXMIO.SXMImage("toy.sxm", Dict{String,String}(), 2, 2, (1.0, 1.0), (0.0, 0.0), [
                Main.STMSXMIO.SXMChannel("Z", "nm", "fwd", [1.0 2.0; 3.0 4.0]),
                Main.STMSXMIO.SXMChannel("Z", "nm", "bwd", [4.0 3.0; 2.0 1.0]),
            ])
            fit_calls = Ref(0)
            view_calls = Ref(0)
            extract_calls = Ref(0)
            count_calls = Ref(0)
            type_calls = Ref(0)
            fake_fit = function (image, pcfg, ccfg)
                fit_calls[] += 1
                return (Any[1], nothing, (xs=[0.0, 1.0], ys=[0.0, 1.0], mask=trues(2, 2), axisctx=(axis=[1.0, 0.0], perp=[0.0, 1.0])))
            end
            fake_view = function (image, channel, pcfg; direction="fwd")
                view_calls[] += 1
                return Main.JointProxyCandidateViews.ViewData(direction, [0.0, 1.0], [0.0, 1.0], [1.0 2.0; 3.0 4.0], 0.1, true)
            end
            fake_extract = function (results, ctx, ccfg, views; patch_half_nm, patch_step_nm, min_pixels=20)
                extract_calls[] += 1
                lobe = Main.JointProxyCandidateViews.LobePatch(1, 0.5, 0.5, 0.1, 0.2, 1.0, 0.3, 0.2, 1.0,
                    Dict("fwd" => [0.1, 0.2], "bwd" => [0.2, 0.1]), Dict("fwd" => [0.3, 0.4], "bwd" => [0.4, 0.3]))
                cand = Main.JointProxyCandidateViews.CandidateNView(2, 1.0, true, true, "ok", 2, false,
                    Any[], 1.0, 0.1, 0.2, 0.8, 1.6, 6, 12, [lobe], [(t=0.0, u=0.0)])
                return Main.JointProxyCandidateViews.CandidateViewReport([cand], Main.JointProxyCandidateViews.SkippedCandidate[], 2, false, ["ok"])
            end
            fake_count = function (report, calibration)
                count_calls[] += 1
                return (; rows=[(; n=2, joint_gcv=1.0, delta_rel=0.0, probability=0.9, rank=1)], predicted_n=2,
                    probability=0.9, abstained=false, confidence=0.9)
            end
            fake_type = function (lobes, ensemble, candidate, calibration)
                type_calls[] += 1
                return (Main.JointProxyTypePosterior.TypePosteriorResult([0.75 0.25], [0], 0.0,
                    Main.JointProxyTypePosterior.TypePosteriorGlobalState[], ["geom" => 1.0]),
                    [(; p0=0.75, p1=0.25, raw_p0=0.75, raw_p1=0.25, predicted=0, confidence=0.75)])
            end
            bundle = JointProxyInferenceCLI._load_bundle_for_test(CONFIG)
            opts = JointProxyInferenceCLI.CliOptions(config=CONFIG, data_dir=dir, files=["toy.sxm"], files_from=nothing,
                calibration=nothing, outdir=dir, chunk=nothing, uncalibrated=true)
            artifacts = JointProxyInferenceCLI._real_inference_artifacts(opts, bundle, ["toy.sxm"],
                read_sxm_fn=_ -> img, fit_runner=fake_fit, view_builder=fake_view, extractor=fake_extract,
                count_posterior_fn=fake_count, type_posterior_fn=fake_type, load_deps=false,
                pcfg=(; label="pcfg"), ccfg=(; label="ccfg"))
            @test fit_calls[] == 1
            @test view_calls[] == 2
            @test extract_calls[] == 1
            @test count_calls[] == 1
            @test type_calls[] == 2
            @test length(artifacts.candidate_n_rows) == 1
            @test artifacts.candidate_n_rows[1].probability == 0.9
        end
    end

    @testset "real adapter converts unusable count and type evidence to abstention" begin
        invalid_candidate = (; valid=false, joint_gcv=Inf, n=2)
        report = (; candidates=[invalid_candidate])
        calibration = JointProxyInferenceCLI.InferenceCalibration(1.0, 0.6, 1.0, 0.6,
            "config", "source", "payload", "calibration")
        count_result = JointProxyInferenceCLI._default_count_result(report, calibration)
        @test count_result.abstained
        @test ismissing(count_result.predicted_n)
        @test isempty(count_result.rows)

        ensemble = JointProxyRegistry.load_registry(CONFIG)
        lobes = [Main.JointProxyTypePosterior.TypePosteriorLobeEvidence(
            Dict("fwd" => fill(NaN, ensemble.npix)))]
        candidate = (; residual_corr=0.0, effective_view_factor=1.0, lobes=[nothing])
        _, type_rows = JointProxyInferenceCLI._default_type_result(lobes, ensemble, candidate, calibration)
        @test ismissing(type_rows[1].predicted)
        @test type_rows[1].p0 == 0.5
        @test type_rows[1].p1 == 0.5
    end

    @testset "stale calibration and chunk errors fail" begin
        mktempdir() do dir
            data_dir = joinpath(dir, "data")
            mkpath(data_dir)
            _write(joinpath(data_dir, "a.sxm"), "stub")
            bad_cal = joinpath(dir, "bad.toml")
            _make_calibration(bad_cal; config_hash="deadbeef", source_hash="deadbeef", payload_hash="deadbeef")
            @test_throws Exception JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm", "--calibration", bad_cal, "--outdir", joinpath(dir, "out")],
                adapter=_fake_artifacts)
            @test_throws Exception JointProxyInferenceCLI.parse_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm", "--calibration", bad_cal, "--outdir", joinpath(dir, "out"), "--chunk", "0/2"])
            @test_throws Exception JointProxyInferenceCLI.parse_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "a.sxm,a.sxm", "--uncalibrated", "--outdir", joinpath(dir, "dup")])
            @test_throws Exception JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", joinpath(dir, "missing"),
                "--files", "a.sxm", "--uncalibrated", "--outdir", joinpath(dir, "out")], adapter=_fake_artifacts)
            @test_throws Exception JointProxyInferenceCLI.run_cli(["--config", CONFIG, "--data-dir", data_dir,
                "--files", "missing.sxm", "--uncalibrated", "--outdir", joinpath(dir, "out")], adapter=_fake_artifacts)
        end
    end

    @testset "uncalibrated mode rejects calibration requests" begin
        @test_throws Exception JointProxyInferenceCLI.parse_cli(["--config", CONFIG, "--data-dir", "data",
            "--files", "a.sxm", "--uncalibrated", "--calibration", "x.toml", "--outdir", "out"])
    end
end
