#!/usr/bin/env julia

# T6 integration tests for the unknown-production hierarchical profile,
# executable leakage firewall, and hierarchical QC abstention reasons.

using Test

const ROOT = dirname(@__DIR__)
const RUNNER = joinpath(@__DIR__, "run_unknown_unit_assignment.jl")
const VALIDATOR = joinpath(@__DIR__, "validate_unit_predictions.jl")
const QC = joinpath(@__DIR__, "summarize_unknown_unit_qc.jl")

function _run(cmd::Cmd)
    out = IOBuffer()
    err = IOBuffer()
    process = run(pipeline(cmd; stdout=out, stderr=err), wait=false)
    wait(process)
    return process.exitcode, String(take!(out)), String(take!(err))
end

_julia(script::AbstractString) = `$(Base.julia_cmd()) --project=$ROOT $script`

function _write_features(path::AbstractString)
    open(path, "w") do io
        println(io, "file\tN\tlobe\tamplitude\tamp_prominence\tamp_neighbor_ratio\tintegrated_prominence\tamp_rel")
        for file in ("scan_A.sxm", "scan_B.sxm"), lobe in 1:6
            high = lobe > 3
            amplitude = high ? 0.10 + 0.001lobe : 0.05 + 0.001lobe
            feature = high ? 1.0 + 0.02lobe : -1.0 + 0.02lobe
            println(io, join((file, 6, lobe, amplitude, feature, feature / 2,
                              feature * 0.8, feature * 0.6), '\t'))
        end
    end
    return path
end

function _rows(path::AbstractString)
    lines = readlines(path)
    header = split(lines[1], '\t'; keepempty=true)
    return [Dict(String(header[i]) => String(value)
                 for (i, value) in enumerate(split(line, '\t'; keepempty=true)))
            for line in lines[2:end]]
end

@testset "T6 unknown hierarchical integration" begin
    mktempdir() do tmp
        features = _write_features(joinpath(tmp, "features.tsv"))

        @testset "hierarchical happy path through production wrapper" begin
            outdir = joinpath(tmp, "hierarchical")
            cmd = `$(_julia(RUNNER)) --features $features --outdir $outdir --profile hierarchical_equalprior --seeds 2`
            code, out, err = _run(cmd)
            @test code == 0
            @test isempty(err)
            @test occursin("profiles: hierarchical_equalprior", out)
            @test sort(readdir(outdir)) == [
                "manifest.tsv",
                "predictions_hierarchical_equalprior.tsv",
                "summary.tsv",
                "validation_hierarchical_equalprior.log",
            ]
            predictions = joinpath(outdir, "predictions_hierarchical_equalprior.tsv")
            rows = _rows(predictions)
            @test !isempty(rows)
            @test all(row["model"] == "hierarchical_equalprior" for row in rows)
            @test all(length(row["provenance_sha256"]) == 64 for row in rows)
            @test occursin("hierarchical_equalprior", read(joinpath(outdir, "summary.tsv"), String))
            manifest = read(joinpath(outdir, "manifest.tsv"), String)
            @test occursin("hierarchical_equalprior", manifest)
            @test occursin("base_local", manifest)
            @test !occursin("base_bwd_consensus", manifest)

            repeat_outdir = joinpath(tmp, "hierarchical_repeat")
            repeat_cmd = `$(_julia(RUNNER)) --features $features --outdir $repeat_outdir --profile hierarchical_equalprior --seeds 2`
            repeat_code, _, repeat_err = _run(repeat_cmd)
            @test repeat_code == 0
            @test isempty(repeat_err)
            @test read(predictions) == read(joinpath(repeat_outdir, "predictions_hierarchical_equalprior.tsv"))
        end

        @testset "invalid profile is rejected" begin
            outdir = joinpath(tmp, "invalid_profile")
            code, _, err = _run(`$(_julia(RUNNER)) --features $features --outdir $outdir --profile not_a_profile`)
            @test code != 0
            @test occursin("profile", lowercase(err))
            @test !isdir(outdir)
        end

        @testset "validator rejects forbidden columns" begin
            predictions = joinpath(tmp, "forbidden_column.tsv")
            write(predictions, "file\tlobe\tpredicted\texpected_N\nscan_A.sxm\t1\t0\t6\n")
            code, _, err = _run(`$(_julia(VALIDATOR)) --predictions $predictions`)
            @test code != 0
            @test occursin("forbidden", lowercase(err))
        end

        @testset "validator rejects forbidden path-valued cells but permits SXM basenames" begin
            bad = joinpath(tmp, "forbidden_path.tsv")
            write(bad, "file\tlobe\tpredicted\tmodel\nscan_A.sxm\t1\t0\tresults/benchmark_grades/output.tsv\n")
            code_bad, _, err_bad = _run(`$(_julia(VALIDATOR)) --predictions $bad`)
            @test code_bad != 0
            @test occursin("path", lowercase(err_bad))

            good = joinpath(tmp, "ordinary_file.tsv")
            write(good, "file\tlobe\tpredicted\tmodel\nscan_grade_report.sxm\t1\t0\thierarchical_equalprior\n")
            code_good, out_good, err_good = _run(`$(_julia(VALIDATOR)) --predictions $good`)
            @test code_good == 0
            @test isempty(err_good)
            @test occursin("status:      ok", out_good)

            for (index, value) in enumerate(("benchmark:truth", "BeNcHmArK:TrUtH"))
                uri = joinpath(tmp, "forbidden_uri_$index.tsv")
                write(uri, "file\tlobe\tpredicted\tmodel\nscan_A.sxm\t1\t0\t$value\n")
                code_uri, _, err_uri = _run(`$(_julia(VALIDATOR)) --predictions $uri`)
                @test code_uri != 0
                @test occursin("path", lowercase(err_uri))
            end

            for (index, value) in enumerate(("file:benchmark_truth.sxm", "FiLe:BeNcHmArK_TrUtH.sxm"))
                uri_file = joinpath(tmp, "forbidden_uri_file_$index.tsv")
                write(uri_file, "file\tlobe\tpredicted\n$value\t1\t0\n")
                code_uri_file, _, err_uri_file = _run(`$(_julia(VALIDATOR)) --predictions $uri_file`)
                @test code_uri_file != 0
                @test occursin("path", lowercase(err_uri_file))
            end
        end

        @testset "QC surfaces hierarchical abstention reasons" begin
            predictions = joinpath(tmp, "hierarchical_abstentions.tsv")
            open(predictions, "w") do io
                println(io, "file\tlobe\tpredicted\tconfidence\tprobability_1\tviews_used\tinvalid_reason\tmodel\tmodel_version\tprovenance_sha256")
                sha = repeat("a", 64)
                for (file, reason) in (("missing.sxm", "missing_view"),
                                       ("unstable.sxm", "degenerate_view"),
                                       ("one.sxm", "one_component_evidence")), lobe in 1:2
                    println(io, join((file, lobe, "?", "0.00000000", "NA", 0,
                                      reason, "hierarchical_equalprior", 1, sha), '\t'))
                end
                for lobe in 1:2
                    println(io, join(("partial.sxm", lobe, "1", "0.75000000",
                                      "0.75000000", 1, "ok_partial_views",
                                      "hierarchical_equalprior", 1, sha), '\t'))
                end
            end
            queue = joinpath(tmp, "review_queue.tsv")
            code, _, err = _run(`$(_julia(QC)) --predictions $predictions --out $queue --profile hierarchical_equalprior`)
            @test code == 0
            @test isempty(err)
            rows = Dict(row["file"] => row for row in _rows(queue))
            @test occursin("hierarchical_missing_view_abstention", rows["missing.sxm"]["review_reasons"])
            @test occursin("hierarchical_unstable_or_degenerate_model_abstention", rows["unstable.sxm"]["review_reasons"])
            @test occursin("hierarchical_one_component_abstention", rows["one.sxm"]["review_reasons"])
            @test occursin("hierarchical_partial_view_prediction", rows["partial.sxm"]["review_reasons"])
            @test !occursin("abstention", rows["partial.sxm"]["review_reasons"])
            @test all(row["review_status"] == "review" for row in values(rows))
        end
    end
end

exit(0)
