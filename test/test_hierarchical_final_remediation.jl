#!/usr/bin/env julia

using Test

const ROOT = dirname(@__DIR__)
const CLI = joinpath(@__DIR__, "build_hierarchical_unit_predictions.jl")

include(joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl"))
using .HierarchicalUnitAssignment
include(joinpath(@__DIR__, "evaluate_hierarchical_unit_assignment.jl"))
const Evaluation = HierarchicalUnitAssignmentEvaluation

function _run(cmd::Cmd)
    stdout = IOBuffer()
    stderr = IOBuffer()
    process = run(pipeline(cmd; stdout=stdout, stderr=stderr), wait=false)
    wait(process)
    return process.exitcode, String(take!(stdout)), String(take!(stderr))
end

function _prediction_rows(path::AbstractString)
    lines = readlines(path)
    header = split(first(lines), '\t'; keepempty=true)
    return [Dict(String(name) => String(value)
                 for (name, value) in zip(header, split(line, '\t'; keepempty=true)))
            for line in lines[2:end]]
end

function _prediction_provenance(path::AbstractString)
    values = unique(row["provenance_sha256"] for row in _prediction_rows(path))
    @test length(values) == 1
    return only(values)
end

function _status_fixture(; converged::Bool, monotone::Bool)
    fit2 = TwoComponentFit(reshape([-1.0, 1.0], 2, 1), ones(2, 1),
                           copy(EQUAL_PRIOR_WEIGHTS), -1.0, [-1.0],
                           converged, monotone)
    fit1 = OneComponentFit([0.0], [1.0], -2.0)
    return ViewFit("fixture", ["amplitude"], fit2, fit1, [1, 2], 2,
                   false, false, [0.1, 0.9], [-2.0, 2.0])
end

@testset "hierarchical final remediation" begin
    @testset "focused hierarchical implementation files stay below 250 pure LOC" begin
        implementation_dir = joinpath(@__DIR__, "lib", "hierarchical")
        implementation_files = sort(filter(path -> endswith(path, ".jl"),
                                           readdir(implementation_dir; join=true)))
        for path in implementation_files
            pure_loc = count(readlines(path)) do line
                stripped = strip(line)
                !isempty(stripped) && !startswith(stripped, '#')
            end
            @test pure_loc < 250
        end
    end

    @testset "public EM parameters reject invalid values before fitting" begin
        function assert_parameter_error(call, parameter::AbstractString)
            err = try
                call()
                nothing
            catch caught
                caught
            end
            @test err isa ArgumentError
            if err isa ArgumentError
                @test occursin(parameter, sprint(showerror, err))
            end
        end

        # One finite row is deliberately insufficient for fitting. Parameter
        # errors must win over data checks, initialization, and start indexing.
        insufficient_X = [0.0]
        for value in (0, -1, 1.5, true)
            assert_parameter_error("n_starts") do
                fit_em_two_component(insufficient_X, 0; n_starts=value)
            end
        end
        for value in (-1, 1.5, true)
            assert_parameter_error("first_seed") do
                fit_em_two_component(insufficient_X, value)
            end
        end
        for value in (0.0, -1.0, Inf, -Inf, NaN, true)
            assert_parameter_error("tol") do
                fit_em_two_component(insufficient_X, 0; tol=value)
            end
        end
        for value in (0.0, -1.0, Inf, -Inf, NaN, true)
            assert_parameter_error("cov_floor") do
                fit_em_two_component(insufficient_X, 0; cov_floor=value)
            end
        end

        # The start offset is an index into the declared deterministic start
        # grid, not a request to silently regenerate or alias another start.
        # Both integer parameters are validated before this relational check,
        # and the relational check precedes all data conversion and fitting.
        poisoned_X = Any["must not be converted"]
        for value in (2, 3, typemax(Int))
            assert_parameter_error("first_seed") do
                fit_em_two_component(poisoned_X, value; n_starts=2,
                                     tol=0.0, cov_floor=0.0)
            end
        end
        assert_parameter_error("n_starts") do
            fit_em_two_component(poisoned_X, 2; n_starts=0,
                                 tol=0.0, cov_floor=0.0)
        end
        assert_parameter_error("first_seed") do
            fit_em_two_component(poisoned_X, -1; n_starts=2,
                                 tol=0.0, cov_floor=0.0)
        end
    end

    @testset "fewer than two occupied hard clusters are one-component evidence" begin
        Z = reshape([-2.0, -1.0, 1.0, 2.0], :, 1)
        records = [LobeRecord("fixture.sxm", i, amplitude,
                              Dict("amplitude" => amplitude))
                   for (i, amplitude) in enumerate((-2.0, -1.0, 1.0, 2.0))]
        fit2 = TwoComponentFit(reshape([-3.0, 3.0], 2, 1), ones(2, 1),
                               copy(EQUAL_PRIOR_WEIGHTS), 100.0, [100.0],
                               true, true)
        fit1 = OneComponentFit([0.0], [1.0], -100.0)

        # Separation, likelihood improvement, and raw-amplitude spread all
        # favor two components; occupancy must nevertheless dominate.
        @test HierarchicalUnitAssignment._one_component_evidence(
            fit2, fit1, Z, records, fill(1, length(records)))
        @test HierarchicalUnitAssignment._one_component_evidence(
            fit2, fit1, Z, records, fill(2, length(records)))
    end

    @testset "unstable and nonmonotone views abstain explicitly" begin
        unstable = combine_views([_status_fixture(converged=false, monotone=true)], 2)
        @test all(isnan, unstable.probability_1)
        @test unstable.views_used == [0, 0]
        @test unstable.invalid_reason == ["unstable_model", "unstable_model"]

        nonmonotone = combine_views([_status_fixture(converged=true, monotone=false)], 2)
        @test all(isnan, nonmonotone.probability_1)
        @test nonmonotone.views_used == [0, 0]
        @test nonmonotone.invalid_reason == ["nonmonotone_model", "nonmonotone_model"]

        mktempdir() do tmp
            records = [LobeRecord("fixture.sxm", i, 0.1i,
                                  Dict("amplitude" => 0.1i)) for i in 1:2]
            for (name, combined) in (("unstable", unstable),
                                     ("nonmonotone", nonmonotone))
                path = joinpath(tmp, "$name.tsv")
                write_predictions(path, records, combined, repeat("a", 64))
                rows = _prediction_rows(path)
                @test all(row["predicted"] == "?" for row in rows)
                @test all(row["invalid_reason"] == "$(name)_model" for row in rows)
            end
        end
    end

    @testset "held-out evaluator rejects unstable fits before scoring" begin
        unstable = Evaluation.HierarchicalUnitAssignment.TwoComponentFit(
            reshape([-1.0, 1.0], 2, 1), ones(2, 1), [0.5, 0.5],
            -1.0, [-1.0], false, true)
        nonmonotone = Evaluation.HierarchicalUnitAssignment.TwoComponentFit(
            reshape([-1.0, 1.0], 2, 1), ones(2, 1), [0.5, 0.5],
            -1.0, [-1.0], true, false)
        @test_throws Evaluation.EvaluationError Evaluation._require_stable_fit(unstable, "fixture")
        @test_throws Evaluation.EvaluationError Evaluation._require_stable_fit(nonmonotone, "fixture")
    end

    @testset "prediction output replaces symlinks and commits atomically" begin
        mktempdir() do tmp
            records = [LobeRecord("fixture.sxm", i, 0.1i,
                                  Dict("amplitude" => 0.1i)) for i in 1:2]
            combined = combine_views([_status_fixture(converged=true, monotone=true)], 2)
            target = joinpath(tmp, "target.tsv")
            destination = joinpath(tmp, "predictions.tsv")
            write(target, "target sentinel\n")
            symlink(target, destination)

            write_predictions(destination, records, combined, repeat("a", 64))

            @test !islink(destination)
            @test startswith(read(destination, String),
                             join(HierarchicalUnitAssignment.PREDICTION_COLUMNS, '\t'))
            @test read(target, String) == "target sentinel\n"
        end

        mktempdir() do tmp
            destination = joinpath(tmp, "stable.tsv")
            write(destination, "committed bytes\n")

            @test_throws ErrorException HierarchicalUnitAssignment._atomic_write(destination) do io
                write(io, "partial replacement\n")
                error("injected writer failure")
            end

            @test read(destination, String) == "committed bytes\n"
            @test readdir(tmp) == ["stable.tsv"]
        end
    end

    @testset "CLI consumes merged split and backward descriptors" begin
        mktempdir() do tmp
            features = joinpath(tmp, "features.tsv")
            split_features = joinpath(tmp, "split.tsv")
            patches = joinpath(tmp, "patches.tsv")
            predictions = joinpath(tmp, "predictions.tsv")
            repeated_predictions = joinpath(tmp, "predictions-repeat.tsv")
            split_changed_predictions = joinpath(tmp, "predictions-split-changed.tsv")
            patch_changed_predictions = joinpath(tmp, "predictions-patch-changed.tsv")

            open(features, "w") do io
                println(io, "file\tlobe\tamplitude")
                for file in ("scan_a.sxm", "scan_b.sxm"), lobe in 1:6
                    println(io, join((file, lobe, lobe <= 3 ? 0.05 : 0.10), '\t'))
                end
            end
            open(split_features, "w") do io
                println(io, "file\tlobe\tskew_ratio")
                for file in ("scan_a.sxm", "scan_b.sxm"), lobe in 1:6
                    println(io, join((file, lobe, lobe <= 3 ? 0.5 : 2.0), '\t'))
                end
            end
            patch_columns = ["bwd_res_p" * lpad(i, 3, '0') for i in 1:9]
            open(patches, "w") do io
                println(io, join(["file", "lobe", patch_columns...], '\t'))
                for file in ("scan_a.sxm", "scan_b.sxm"), lobe in 1:6
                    values = zeros(9)
                    values[lobe <= 3 ? 1 : 3] = -1.0
                    println(io, join((file, lobe, values...), '\t'))
                end
            end

            original_split = read(split_features, String)
            original_patches = read(patches, String)

            cmd = `$(Base.julia_cmd()) --project=$ROOT $CLI --features $features --split-features $split_features --patches $patches --out $predictions --view split=split_log_skew --view backward=bwd_neg_com_t --em-starts 2`
            code, stdout, stderr = _run(cmd)
            @test code == 0
            @test isempty(stderr)
            @test occursin("view split", stdout)
            @test occursin("view backward", stdout)
            rows = _prediction_rows(predictions)
            @test !isempty(rows)
            @test all(row["views_used"] == "2" for row in rows)
            @test all(row["predicted"] in ("0", "1") for row in rows)

            baseline_provenance = _prediction_provenance(predictions)

            repeat_cmd = `$(Base.julia_cmd()) --project=$ROOT $CLI --features $features --split-features $split_features --patches $patches --out $repeated_predictions --view split=split_log_skew --view backward=bwd_neg_com_t --em-starts 2`
            repeat_code, _, repeat_stderr = _run(repeat_cmd)
            @test repeat_code == 0
            @test isempty(repeat_stderr)
            @test read(predictions) == read(repeated_predictions)

            write(split_features, replace(original_split, "0.5" => "0.6"; count=1))
            split_cmd = `$(Base.julia_cmd()) --project=$ROOT $CLI --features $features --split-features $split_features --patches $patches --out $split_changed_predictions --view split=split_log_skew --view backward=bwd_neg_com_t --em-starts 2`
            split_code, _, split_stderr = _run(split_cmd)
            @test split_code == 0
            @test isempty(split_stderr)
            @test _prediction_provenance(split_changed_predictions) != baseline_provenance

            write(split_features, original_split)
            write(patches, replace(original_patches, "-1.0" => "-0.9"; count=1))
            patch_cmd = `$(Base.julia_cmd()) --project=$ROOT $CLI --features $features --split-features $split_features --patches $patches --out $patch_changed_predictions --view split=split_log_skew --view backward=bwd_neg_com_t --em-starts 2`
            patch_code, _, patch_stderr = _run(patch_cmd)
            @test patch_code == 0
            @test isempty(patch_stderr)
            @test _prediction_provenance(patch_changed_predictions) != baseline_provenance
        end
    end
end
