#!/usr/bin/env julia

using SHA
using Test

const ROOT = dirname(@__DIR__)
const EVALUATOR = joinpath(@__DIR__, "evaluate_hierarchical_unit_assignment.jl")
const MERGER = joinpath(@__DIR__, "merge_hierarchical_unit_assignment_shards.jl")
const ARRAY = joinpath(ROOT, "hpc", "hierarchical_unit_assignment_array.sbatch")
const LAUNCHER = joinpath(ROOT, "hpc", "launch_hierarchical_unit_assignment.sh")

@test isfile(EVALUATOR)
@test isfile(MERGER)
@test isfile(ARRAY)
@test isfile(LAUNCHER)

if isfile(EVALUATOR)
    include(EVALUATOR)
end
if isfile(MERGER)
    include(MERGER)
end

function _write_shard(path, rows; header=HierarchicalUnitAssignmentMerge.SHARD_COLUMNS)
    open(path, "w") do io
        println(io, join(header, '\t'))
        for row in rows
            println(io, join((string(row[c]) for c in header), '\t'))
        end
    end
end

function _row(; fold="20240101", kind="unbootstrapped", seed=-1,
              block="unbootstrapped", scan="20240101/a.sxm", improvement=0.5,
              occupancy=0.5, config_hash=repeat("a", 64),
              input_hash=repeat("b", 64), source_hash=repeat("c", 64),
              expected_rows=1)
    return Dict(
        "schema_version" => "1", "fold_date" => fold,
        "bootstrap_kind" => kind, "bootstrap_seed" => seed,
        "bootstrap_block" => block, "scan_key" => scan,
        "scan_date" => fold, "training_scan_count" => 2,
        "heldout_scan_count" => 1, "lobe_count" => 6,
        "ll_one" => -20.0, "ll_two" => -17.0,
        "improvement" => 3.0, "improvement_per_lobe" => improvement,
        "forward_backward_agreement" => 1.0,
        "perturbation_min_agreement" => 1.0,
        "occupancy_1" => occupancy, "config_sha256" => config_hash,
        "input_sha256" => input_hash, "source_sha256" => source_hash,
        "resample_unit" => "whole_scan_within_training_dates",
        "row_complete" => "true", "shard_complete" => "true",
        "expected_rows" => expected_rows, "expected_seed_count" => 500,
    )
end

function _complete_rows(; improvement=0.5, occupancy=0.5)
    rows = [_row(improvement=improvement, occupancy=occupancy)]
    for seed in 0:499
        push!(rows, _row(kind="bootstrap", seed=seed,
                         block="$(50 * (seed ÷ 50)):$(50 * (seed ÷ 50) + 49)",
                         improvement=improvement, occupancy=occupancy))
    end
    for row in rows
        row["expected_rows"] = length(rows)
    end
    return rows
end

@testset "leading date and folds" begin
    if isdefined(Main, :HierarchicalUnitAssignmentEvaluation)
        E = HierarchicalUnitAssignmentEvaluation
        @test E.parse_leading_date("20240102/scan.sxm") == "20240102"
        @test E.parse_leading_date("folder/20240102_scan.sxm") == "20240102"
        @test_throws E.EvaluationError E.parse_leading_date("folder/scan.sxm")
        @test_throws E.EvaluationError E.parse_leading_date("20240102/20240103_scan.sxm")
        @test_throws E.EvaluationError E.parse_leading_date("20241340_scan.sxm")

        paths = ["20240103/c.sxm", "20240101/a.sxm", "20240102/b.sxm"]
        folds1 = E.build_date_folds(paths)
        folds2 = E.build_date_folds(reverse(paths))
        @test folds1 == folds2
        @test [f.fold_date for f in folds1] == ["20240101", "20240102", "20240103"]
        @test folds1[1].heldout_scans == ["20240101/a.sxm"]
        @test folds1[1].training_scans == ["20240102/b.sxm", "20240103/c.sxm"]
    end
end

@testset "whole-scan bootstrap" begin
    if isdefined(Main, :HierarchicalUnitAssignmentEvaluation)
        E = HierarchicalUnitAssignmentEvaluation
        scans = ["20240101/a.sxm", "20240101/b.sxm",
                 "20240102/c.sxm", "20240102/d.sxm"]
        a = E.whole_scan_bootstrap(scans, 17)
        b = E.whole_scan_bootstrap(reverse(scans), 17)
        @test a == b
        @test length(a) == length(scans)
        @test all(s in scans for s in a)
        @test [E.parse_leading_date(s) for s in a] ==
              sort([E.parse_leading_date(s) for s in scans])
        @test E.BOOTSTRAP_SEEDS == collect(0:499)
        @test_throws E.EvaluationError E.bootstrap_training_units(scans, 0; resample_unit="lobe")
    end
end

@testset "scientific gate is strict and occupancy-free" begin
    if isdefined(Main, :HierarchicalUnitAssignmentMerge)
        M = HierarchicalUnitAssignmentMerge
        pass_rows = _complete_rows(improvement=0.25, occupancy=0.01)
        g1 = M.evaluate_gate(pass_rows, ["20240101"])
        @test g1.passed
        pass_rows2 = _complete_rows(improvement=0.25, occupancy=0.99)
        g2 = M.evaluate_gate(pass_rows2, ["20240101"])
        @test g2 == g1

        zero_fold = _complete_rows(improvement=0.25)
        zero_fold[1]["improvement_per_lobe"] = 0.0
        @test !M.evaluate_gate(zero_fold, ["20240101"]).passed

        zero_lcb = _complete_rows(improvement=0.0)
        zero_lcb[1]["improvement_per_lobe"] = 0.1
        @test M.evaluate_gate(zero_lcb, ["20240101"]).lower_bound == 0.0
        @test !M.evaluate_gate(zero_lcb, ["20240101"]).passed
    end
end

@testset "strict shard merge corruptions" begin
    if isdefined(Main, :HierarchicalUnitAssignmentMerge)
        M = HierarchicalUnitAssignmentMerge
        mktempdir() do tmp
            expected = M.MergeContract(
                ["20240101"], repeat("a", 64), repeat("b", 64), repeat("c", 64))
            good = joinpath(tmp, "good.tsv")
            _write_shard(good, _complete_rows())
            result = M.merge_shards([good], joinpath(tmp, "merged.tsv"), expected)
            @test result.gate.passed
            @test length(result.rows) == 501

            cases = Dict{String,Vector{Dict{String,Any}}}()
            rows = _complete_rows(); deleteat!(rows, 2); cases["missing_seed"] = rows
            rows = _complete_rows(); push!(rows, copy(rows[2])); cases["duplicate"] = rows
            rows = _complete_rows(); rows[2]["bootstrap_seed"] = 500; cases["wrong_seed"] = rows
            rows = _complete_rows(); rows[1]["fold_date"] = "20240202"; cases["wrong_fold"] = rows
            rows = _complete_rows(); rows[1]["config_sha256"] = repeat("d", 64); cases["wrong_config_hash"] = rows
            rows = _complete_rows(); rows[1]["input_sha256"] = repeat("d", 64); cases["wrong_input_hash"] = rows
            rows = _complete_rows(); rows[1]["source_sha256"] = repeat("d", 64); cases["wrong_source_hash"] = rows
            rows = _complete_rows(); rows[1]["shard_complete"] = "false"; cases["partial"] = rows
            for (name, badrows) in cases
                path = joinpath(tmp, "$name.tsv")
                _write_shard(path, badrows)
                @test_throws M.MergeError M.merge_shards([path], joinpath(tmp, "$name-out.tsv"), expected)
            end
            malformed = joinpath(tmp, "malformed.tsv")
            _write_shard(malformed, _complete_rows(); header=M.SHARD_COLUMNS[1:end-1])
            @test_throws M.MergeError M.merge_shards([malformed], joinpath(tmp, "malformed-out.tsv"), expected)
        end
    end
end

@testset "merged TSV and gate publish as a recoverable set" begin
    if isdefined(Main, :HierarchicalUnitAssignmentMerge)
        M = HierarchicalUnitAssignmentMerge
        mktempdir() do tmp
            expected = M.MergeContract(
                ["20240101"], repeat("a", 64), repeat("b", 64), repeat("c", 64))
            shard = joinpath(tmp, "complete.tsv")
            merged = joinpath(tmp, "merged.tsv")
            gate = string(merged, ".gate.tsv")
            _write_shard(shard, _complete_rows())
            write(merged, "old merged\n")
            write(gate, "old gate\n")

            @test_throws ErrorException M.merge_shards(
                [shard], merged, expected; failpoint=:before_gate)
            @test read(merged, String) == "old merged\n"
            @test read(gate, String) == "old gate\n"

            @test_throws M.MergeTransactionInterruption M.merge_shards(
                [shard], merged, expected; interruptpoint=:after_install_1)
            @test startswith(read(merged, String), join(M.SHARD_COLUMNS, '\t'))
            @test !ispath(gate)

            result = M.merge_shards([shard], merged, expected)
            @test result.gate.passed
            @test startswith(read(merged, String), join(M.SHARD_COLUMNS, '\t'))
            @test startswith(read(gate, String), "terminal_status\t")
            @test isempty(filter(name -> occursin("stmfit-txn", name), readdir(tmp)))
        end
    end
end

@testset "merge transaction recovery is deterministic" begin
    if isdefined(Main, :HierarchicalUnitAssignmentMerge)
        M = HierarchicalUnitAssignmentMerge
        for (interruptpoint, visible_generation) in (
                (:after_marker, :old),
                (:after_backup_1, :old),
                (:after_backup_2, :old),
                (:after_install_1, :old),
                (:after_gate, :old),
                (:after_commit, :new))
            mktempdir() do tmp
                merged = joinpath(tmp, "merged.tsv")
                gate = string(merged, ".gate.tsv")
                destinations = [merged, gate]
                write(merged, "old merged\n")
                write(gate, "old gate\n")

                @test_throws M.MergeTransactionInterruption M._with_merge_transaction(
                        destinations; gate, interruptpoint) do staged
                    write(staged[abspath(merged)], "new merged\n")
                    write(staged[abspath(gate)], "new gate\n")
                end
                if interruptpoint == :after_install_1
                    @test read(merged, String) == "new merged\n"
                    @test !ispath(gate)
                end

                observed = Ref{Symbol}(:unset)
                M._with_merge_transaction(destinations; gate) do staged
                    observed[] = startswith(read(merged, String), "old") ? :old : :new
                    write(staged[abspath(merged)], "rerun merged\n")
                    write(staged[abspath(gate)], "rerun gate\n")
                end
                @test observed[] == visible_generation
                @test read(merged, String) == "rerun merged\n"
                @test read(gate, String) == "rerun gate\n"
                @test isempty(filter(name -> occursin("stmfit-txn", name), readdir(tmp)))
            end
        end
    end
end

@testset "merge transaction refuses unsafe markers" begin
    if isdefined(Main, :HierarchicalUnitAssignmentMerge)
        M = HierarchicalUnitAssignmentMerge
        for marker_kind in (:malformed, :symlink)
            mktempdir() do tmp
                merged = joinpath(tmp, "merged.tsv")
                gate = string(merged, ".gate.tsv")
                marker = M._txn_marker(gate)
                write(merged, "old merged\n")
                write(gate, "old gate\n")
                if marker_kind == :malformed
                    write(marker, "not a transaction marker\n")
                else
                    victim = joinpath(tmp, "marker-victim")
                    write(victim, "protected marker target\n")
                    symlink(victim, marker)
                end

                @test_throws M.MergeError M._with_merge_transaction([merged, gate]; gate) do staged
                    write(staged[abspath(merged)], "new merged\n")
                    write(staged[abspath(gate)], "new gate\n")
                end
                @test read(merged, String) == "old merged\n"
                @test read(gate, String) == "old gate\n"
                @test marker_kind != :symlink || read(readlink(marker), String) == "protected marker target\n"
            end
        end
    end
end

@testset "merge transaction replaces destination symlinks without following them" begin
    if isdefined(Main, :HierarchicalUnitAssignmentMerge)
        M = HierarchicalUnitAssignmentMerge
        for failpoint in (:before_gate, nothing)
            mktempdir() do tmp
                merged = joinpath(tmp, "merged.tsv")
                gate = string(merged, ".gate.tsv")
                merged_target = joinpath(tmp, "merged-target.tsv")
                gate_target = joinpath(tmp, "gate-target.tsv")
                write(merged_target, "protected merged target\n")
                write(gate_target, "protected gate target\n")
                symlink(merged_target, merged)
                symlink(gate_target, gate)

                publish = () -> M._with_merge_transaction([merged, gate]; gate, failpoint) do staged
                    write(staged[abspath(merged)], "new merged\n")
                    write(staged[abspath(gate)], "new gate\n")
                end
                if failpoint === nothing
                    publish()
                    @test !islink(merged)
                    @test !islink(gate)
                    @test read(merged, String) == "new merged\n"
                    @test read(gate, String) == "new gate\n"
                else
                    @test_throws ErrorException publish()
                    @test islink(merged)
                    @test islink(gate)
                    @test readlink(merged) == merged_target
                    @test readlink(gate) == gate_target
                end
                @test read(merged_target, String) == "protected merged target\n"
                @test read(gate_target, String) == "protected gate target\n"
                @test isempty(filter(name -> occursin("stmfit-txn", name), readdir(tmp)))
            end
        end
    end
end

@testset "hierarchical shard and gate outputs replace symlinks atomically" begin
    if isdefined(Main, :HierarchicalUnitAssignmentEvaluation) &&
       isdefined(Main, :HierarchicalUnitAssignmentMerge)
        E = HierarchicalUnitAssignmentEvaluation
        M = HierarchicalUnitAssignmentMerge
        writers = [
            ("evaluator", (path, rows) -> E._write_rows(path, rows)),
            ("merged", (path, rows) -> M._write_rows(path, rows)),
        ]
        for (name, writer) in writers
            mktempdir() do tmp
                target = joinpath(tmp, "$name-target.tsv")
                destination = joinpath(tmp, "$name.tsv")
                write(target, "target sentinel\n")
                symlink(target, destination)

                writer(destination, [_row()])

                @test !islink(destination)
                @test startswith(read(destination, String), join(M.SHARD_COLUMNS, '\t'))
                @test read(target, String) == "target sentinel\n"
            end
        end

        mktempdir() do tmp
            target = joinpath(tmp, "gate-target.tsv")
            destination = joinpath(tmp, "gate.tsv")
            write(target, "target sentinel\n")
            symlink(target, destination)
            gate = M.GateResult(true, true, 0.25, ["20240101" => 0.5], [0.25])

            M._write_gate(destination, gate)

            @test !islink(destination)
            @test startswith(read(destination, String), "terminal_status\t")
            @test read(target, String) == "target sentinel\n"
        end

        for module_under_test in (E, M)
            mktempdir() do tmp
                destination = joinpath(tmp, "stable.tsv")
                write(destination, "committed bytes\n")

                @test_throws ErrorException module_under_test._atomic_write(destination) do io
                    write(io, "partial replacement\n")
                    error("injected writer failure")
                end

                @test read(destination, String) == "committed bytes\n"
                @test readdir(tmp) == ["stable.tsv"]
            end
        end
    end
end

@testset "evaluator to complete 500-seed merge round trip" begin
    if isdefined(Main, :HierarchicalUnitAssignmentEvaluation) &&
       isdefined(Main, :HierarchicalUnitAssignmentMerge)
        E = HierarchicalUnitAssignmentEvaluation
        M = HierarchicalUnitAssignmentMerge
        mktempdir() do tmp
            features = joinpath(tmp, "label_free_features.tsv")
            open(features, "w") do io
                println(io, "file\tlobe\tamplitude\tamp_prominence\tamp_neighbor_ratio\tintegrated_prominence\tamp_rel")
                values = [-2.0, -1.8, -1.6, -1.4, 1.4, 1.6, 1.8, 2.0]
                scans = ["20240101/a.sxm", "20240101/b.sxm",
                         "20240102/c.sxm", "20240102/d.sxm"]
                for (scan_index, scan) in enumerate(scans), (lobe, value) in enumerate(values)
                    amplitude = value < 0 ? 0.05 + scan_index / 1000 : 0.10 + scan_index / 1000
                    println(io, join((scan, lobe, amplitude, value,
                                      value * 0.5, value * 0.8, value * 0.6), '\t'))
                end
            end
            config = joinpath(ROOT, "config", "unit_assignment_candidate.toml")
            shards = String[]
            blocks = ["unbootstrapped"; ["$start:$(start + 49)" for start in 0:50:450]]
            for fold in ("20240101", "20240102"), block in blocks
                rows = E.evaluate_fold(features, config, fold, block;
                                       view_specs=["base" => ["amp_prominence"]])
                path = joinpath(tmp, "fold_$(fold)_$(replace(block, ':' => '-')).tsv")
                E._write_rows(path, rows)
                if fold == "20240101" && block == "0:49"
                    repeat_rows = E.evaluate_fold(features, config, fold, block;
                                                  view_specs=["base" => ["amp_prominence"]])
                    repeat_path = joinpath(tmp, "deterministic-repeat.tsv")
                    E._write_rows(repeat_path, repeat_rows)
                    @test read(path) == read(repeat_path)
                end
                push!(shards, path)
            end
            contract = M.MergeContract(
                ["20240101", "20240102"], bytes2hex(sha256(read(config))),
                bytes2hex(sha256(read(features))), E._source_sha256())
            result = M.merge_shards(shards, joinpath(tmp, "merged.tsv"), contract)
            result_repeat = M.merge_shards(shards, joinpath(tmp, "merged-repeat.tsv"), contract)
            @test length(result.rows) == 2 * 2 * 501
            @test Set(row["bootstrap_seed"] for row in result.rows
                      if row["bootstrap_kind"] == "bootstrap") == Set(0:499)
            @test isfile(result.gate_out)
            @test read(result.out) == read(result_repeat.out)
            @test read(result.gate_out) == read(result_repeat.gate_out)
        end
    end
end

@testset "HPC launcher and array contract" begin
    if all(isfile, (ARRAY, LAUNCHER))
        array_text = read(ARRAY, String)
        launcher_text = read(LAUNCHER, String)
        @test occursin("--cpus-per-task=1", array_text)
        @test occursin("JULIA_NUM_THREADS=1", array_text)
        @test occursin("bootstrap_block", lowercase(array_text))
        @test occursin("date_fold", lowercase(array_text))
        @test occursin("--array", launcher_text)
        @test occursin("%8", launcher_text)
        @test occursin("--config-sha256", launcher_text)
        @test occursin("--input-sha256", launcher_text)
        @test occursin("--source-sha256", launcher_text)
        @test occursin("--dry-run", launcher_text)
        @test occursin("--submit", launcher_text)
        @test occursin("--watch", launcher_text)
        @test occursin("instantiate", lowercase(launcher_text))
        @test occursin("merge_hierarchical_unit_assignment_shards.jl", launcher_text)
        @test occursin("rsync", launcher_text)
        @test occursin("JobIDRaw,State,ExitCode", launcher_text)
        @test !occursin("ArrayTaskID", launcher_text)
        @test !occursin("rm -rf", launcher_text)
        @test success(`bash -n $LAUNCHER`)
        @test success(`bash -n $ARRAY`)
        @test !success(pipeline(`$LAUNCHER`; stdout=devnull, stderr=devnull))
        @test !success(pipeline(`$LAUNCHER --invalid`; stdout=devnull, stderr=devnull))
    end
end
