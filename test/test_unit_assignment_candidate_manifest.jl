#!/usr/bin/env julia

# T0 contract tests for the unit-assignment challenger firewall and baseline.
#
# Two test sets:
#   1. baseline_firewall           — pins the CURRENT observable boundary of the
#                                    existing unknown-production runner/validator
#                                    and docs firewall, and the separation of the
#                                    grader-only 145/870 denominator manifest.
#   2. candidate_manifest_contract — exercises test/check_unit_assignment_candidate_manifest.jl
#                                    and config/unit_assignment_candidate.toml: the
#                                    locked contract, firewall scan, grader-only
#                                    exception, hash binding, and mutation rejection.
#
# Run: julia --project=. test/test_unit_assignment_candidate_manifest.jl

using Test, SHA, TOML

const ROOT = dirname(@__DIR__)
const CHECKER = joinpath(@__DIR__, "check_unit_assignment_candidate_manifest.jl")
const CONFIG = joinpath(ROOT, "config", "unit_assignment_candidate.toml")
const DENOMINATOR = joinpath(ROOT, "benchmarks", "chitosan_6mer_counting_confirmed.toml")
const UNKNOWN_RUNNER = joinpath(@__DIR__, "run_unknown_unit_assignment.jl")
const MANDATORY_PROVENANCE_FIELDS = (
    "feature_tsves",
    "constant_current_cubes",
    "generated_maps",
    "molds",
    "configs",
)

# Capture exit code WITHOUT throwing on nonzero.
function _run(cmd::Cmd)
    out = IOBuffer()
    err = IOBuffer()
    p = run(pipeline(cmd; stdout=out, stderr=err), wait=false)
    wait(p)
    return p.exitcode, String(take!(out)), String(take!(err))
end

_julia_cmd() = `$(Base.julia_cmd()) --project=$ROOT $CHECKER`

# Copy the real config to a fresh temp dir, apply mutation(text) -> text, run the checker.
function _run_mutated(mutation::Function; expect_locked=false, expect_frozen_once=false)
    mktempdir() do tmp
        fixture = joinpath(tmp, "candidate.toml")
        cp(CONFIG, fixture)
        write(fixture, mutation(read(fixture, String)))
        cmd = `$(_julia_cmd()) --config $fixture`
        expect_locked && (cmd = `$cmd --expect-locked`)
        expect_frozen_once && (cmd = `$cmd --expect-frozen-once`)
        return _run(cmd)
    end
end

# Hash fixture source exactly as the checker contract requires: retain every
# source byte outside [grader_only], including line endings, while omitting the
# stored candidate.frozen_hash and candidate.grade_status assignments. The
# latter is lifecycle metadata: frozen_once may become graded without changing
# the immutable candidate payload.
function _exact_non_grader_hash(text::String)
    section = ""
    buf = IOBuffer()
    for line in eachline(IOBuffer(text); keep=true)
        ending = endswith(line, "\r\n") ? "\r\n" : endswith(line, "\n") ? "\n" : endswith(line, "\r") ? "\r" : ""
        stripped = strip(line)
        header = match(r"^\[([A-Za-z0-9_.-]+)\](?:\s*#.*)?$", stripped)
        header === nothing || (section = header.captures[1])
        (section == "grader_only" || startswith(section, "grader_only.")) && continue
        if section == "candidate" && occursin(r"^(?:frozen_hash|grade_status)\s*=", stripped)
            write(buf, ending)
        else
            write(buf, line)
        end
    end
    return bytes2hex(sha256(take!(buf)))
end

function _with_recomputed_frozen_hash(text::String)
    blanked = replace(text, r"(?m)^frozen_hash\s*=.*$" => "frozen_hash = \"\""; count=1)
    hash = _exact_non_grader_hash(blanked)
    return replace(blanked, "frozen_hash = \"\"" => "frozen_hash = \"$hash\""; count=1)
end

function _frozen_fixture_text(; grade_status="frozen_once", line_ending="\n")
    text = replace(read(CONFIG, String), "\r\n" => "\n")
    text = replace(text, "grade_status = \"locked\"" => "grade_status = \"$grade_status\""; count=1)
    text = replace(text, "status = \"pending\"" => "status = \"frozen\""; count=1)
    for (index, field) in enumerate(MANDATORY_PROVENANCE_FIELDS)
        digest = lpad(string(index; base=16), 64, '0')
        text = replace(text, "$field = \"pending\"" => "$field = \"$digest\""; count=1)
    end
    line_ending == "\n" || (text = replace(text, "\n" => line_ending))
    text = replace(text, "[candidate]$line_ending" => "[candidate]$(line_ending)frozen_hash = \"\"$line_ending"; count=1)
    return _with_recomputed_frozen_hash(text)
end

function _run_fixture(text::String; expect_frozen_once=false)
    mktempdir() do tmp
        fixture = joinpath(tmp, "candidate.toml")
        write(fixture, text)
        cmd = `$(_julia_cmd()) --config $fixture`
        expect_frozen_once && (cmd = `$cmd --expect-frozen-once`)
        return _run(cmd)
    end
end

# ===========================================================================
@testset "baseline_firewall" begin
    @testset "denominator manifest is distinct and fixes 145/870" begin
        @test isfile(DENOMINATOR)
        m = TOML.parsefile(DENOMINATOR)
        @test length(m["files"]) == 145
    end

    @testset "unknown runner/validator do not reference denominator manifest" begin
        runner = read(joinpath(ROOT, "test", "run_unknown_unit_assignment.jl"), String)
        validator = read(joinpath(ROOT, "test", "validate_unit_predictions.jl"), String)
        @test !occursin("chitosan_6mer_counting_confirmed", runner)
        @test !occursin("chitosan_6mer_counting_confirmed", validator)
    end

    mktempdir() do tmp
        pred = joinpath(tmp, "pred.tsv")
        write(pred, "file\tlobe\tpredicted\n240817_001.sxm\t1\t0\n240817_001.sxm\t2\t1\n")
        j = `$(Base.julia_cmd()) --project=$ROOT`

        @testset "validator rejects benchmark-only flags" begin
            for flag in ("--truth", "--control-sequence", "--manifest", "--full145", "--control")
                cmd = `$j $(joinpath(@__DIR__, "validate_unit_predictions.jl")) --predictions $pred $flag x`
                code, _, _ = _run(cmd)
                @test code != 0
            end
        end

        @testset "validator clean (label-free) succeeds" begin
            cmd = `$j $(joinpath(@__DIR__, "validate_unit_predictions.jl")) --predictions $pred`
            code, out, _ = _run(cmd)
            @test code == 0
            @test occursin("status:      ok", out)
        end

        @testset "unknown runner rejects benchmark-only flags" begin
            for (flag, val) in (("--manifest", DENOMINATOR), ("--truth", "x.tsv"), ("--expected-N", "6"), ("--full145", nothing), ("--control", nothing))
                outdir = joinpath(tmp, "out_" * replace(flag, r"[^A-Za-z0-9]" => "_"))
                cmd = if val === nothing
                    `$j $(joinpath(@__DIR__, "run_unknown_unit_assignment.jl")) --features $pred --outdir $outdir $flag`
                else
                    `$j $(joinpath(@__DIR__, "run_unknown_unit_assignment.jl")) --features $pred --outdir $outdir $flag $val`
                end
                code, _, _ = _run(cmd)
                @test code != 0
            end
        end

        @testset "docs firewall passes" begin
            cmd = `$j $(joinpath(@__DIR__, "check_unknown_workflow_docs.jl"))`
            code, _, _ = _run(cmd)
            @test code == 0
        end

        @testset "historical default profile behavior is unchanged" begin
            features = joinpath(tmp, "features.tsv")
            split_features = joinpath(tmp, "split.tsv")
            patches = joinpath(tmp, "patches.tsv")
            outdir = joinpath(tmp, "default_profile")
            open(features, "w") do io
                println(io, "file\tN\tlobe\tamplitude\tamp_prominence\tamp_neighbor_ratio\tintegrated_prominence\tamp_rel")
                for file in ("scan_A.sxm", "scan_B.sxm"), lobe in 1:6
                    high = lobe > 3
                    amp = high ? 0.10 + 0.001lobe : 0.05 + 0.001lobe
                    feature = high ? 1.0 + 0.02lobe : -1.0 + 0.02lobe
                    println(io, join((file, 6, lobe, amp, feature, feature / 2, feature * 0.8, feature * 0.6), '\t'))
                end
            end
            open(split_features, "w") do io
                println(io, "file\tlobe\tskew_ratio")
                for file in ("scan_A.sxm", "scan_B.sxm"), lobe in 1:6
                    println(io, join((file, lobe, 0.8 + 0.05lobe), '\t'))
                end
            end
            open(patches, "w") do io
                println(io, "file\tlobe\tbwd_res_p001")
                for file in ("scan_A.sxm", "scan_B.sxm"), lobe in 1:6
                    println(io, join((file, lobe, -0.1lobe), '\t'))
                end
            end

            cmd = `$j $UNKNOWN_RUNNER --features $features --split-features $split_features --patches $patches --outdir $outdir --seeds 2`
            code, out, err = _run(cmd)
            @test code == 0
            @test isempty(err)
            @test occursin("profiles: base_bwd_consensus, base_split_log_skew", out)
            @test sort(readdir(outdir)) == [
                "manifest.tsv",
                "predictions_base_bwd_consensus.tsv",
                "predictions_base_split_log_skew.tsv",
                "summary.tsv",
                "validation_base_bwd_consensus.log",
                "validation_base_split_log_skew.log",
            ]
            summary = readlines(joinpath(outdir, "summary.tsv"))
            @test split(summary[1], '\t') == ["profile", "prediction_path", "files", "lobes", "predicted_0", "predicted_1", "uncertain", "mean_confidence", "features", "split_features", "patches", "seeds", "first_seed", "interactions"]
            @test [split(line, '\t')[1] for line in summary[2:end]] == ["base_bwd_consensus", "base_split_log_skew"]
            manifest = readlines(joinpath(outdir, "manifest.tsv"))
            @test split(manifest[1], '\t') == ["profile", "prediction_path", "validation_log", "validation_status", "files", "lobes", "input_path_features", "input_path_split_features", "input_path_patches", "split_available", "patches_available", "seeds", "first_seed", "interactions", "views"]
            @test [split(line, '\t')[1] for line in manifest[2:end]] == ["base_bwd_consensus", "base_split_log_skew"]
            @test split(manifest[2], '\t')[15] == "base,base_bwd_neg_com_t,base_bwd_neg_diag45"
            @test split(manifest[3], '\t')[15] == "base_split_log_skew"
            @test !occursin("hierarchical", read(joinpath(outdir, "summary.tsv"), String))
            @test !occursin("hierarchical", read(joinpath(outdir, "manifest.tsv"), String))
        end
    end
end

# ===========================================================================
@testset "candidate_manifest_contract" begin
    @testset "checker and config artifacts exist" begin
        @test isfile(CHECKER)
        @test isfile(CONFIG)
    end

    @testset "checker --help exits 0" begin
        code, out, _ = _run(`$(_julia_cmd()) --help`)
        @test code == 0
        @test occursin("Usage", out)
    end

    @testset "real config passes --expect-locked" begin
        cmd = `$(_julia_cmd()) --config $CONFIG --expect-locked`
        code, out, _ = _run(cmd)
        @test code == 0
        @test occursin("grade_status", out)
        @test occursin("locked", out)
        @test occursin("ok", lowercase(out))
    end

    @testset "missing config fails" begin
        cmd = `$(_julia_cmd()) --config $(joinpath(dirname(CONFIG), "does_not_exist.toml"))`
        code, _, _ = _run(cmd)
        @test code != 0
    end

    @testset "missing --config fails" begin
        code, _, _ = _run(`$(_julia_cmd())`)
        @test code != 0
    end

    @testset "actual TOML comments obey the sole grader-only boundary" begin
        legitimate = """
        # Human-readable rationale without benchmark vocabulary.
        # Unicode rationale is legitimate: μ, Å, 中文.
        comment_probe = "a quoted # remains string data" # Benign inline note.
        """
        code0, _, err0 = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$legitimate"; count=1))
        @test code0 == 0
        @test isempty(err0)

        grader_comment = "# expected_N = 6"
        code1, _, err1 = _run_mutated(text -> replace(text, r"(?m)^\[grader_only\]\n" => "[grader_only]\n$grader_comment\n"; count=1))
        @test code1 == 0
        @test isempty(err1)

        code2, _, err2 = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$grader_comment"; count=1))
        @test code2 != 0
        @test occursin("forbidden", lowercase(err2))
    end

    @testset "complete forbidden vocabulary is scanned in actual comments" begin
        forbidden_tokens = (
            "NKNNKN",
            "010010",
            "sequence",
            "expected_N",
            "target_N",
        )
        forbidden_paths = (
            "benchmarks/",
            "results/benchmark_grades",
            "report_unit_assignment_benchmark",
            "grade_unit_assignment",
            "_unit_sequences.tsv",
        )
        forbidden = (forbidden_tokens..., forbidden_paths...)
        payload = "# " * join(forbidden, " | ")

        code0, _, err0 = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$payload"; count=1))
        @test code0 != 0
        for token in forbidden_tokens
            @test occursin("forbidden token '$token'", err0)
        end
        for path in forbidden_paths
            @test occursin("forbidden path fragment '$path'", err0)
        end

        code1, _, err1 = _run_mutated(text -> replace(text, r"(?m)^\[grader_only\]\n" => "[grader_only]\n$payload\n"; count=1))
        @test code1 == 0
        @test isempty(err1)
    end

    @testset "fake and non-grader table contexts cannot exempt forbidden comments" begin
        forbidden_comment = "# expected_N = 6"
        mutations = (
            text -> replace(text, "[candidate]" => "[candidate]\n[candidate.child]\n$forbidden_comment"; count=1),
            text -> replace(text, "[candidate]" => "[candidate]\n[[candidate.items]]\n$forbidden_comment"; count=1),
            text -> replace(text, r"(?m)^\[grader_only\]\n" => "[grader_only.extra]\n$forbidden_comment\n[grader_only]\n"; count=1),
            text -> replace(text, "[candidate]" => "[candidate]\n# [grader_only]\n$forbidden_comment"; count=1),
            text -> replace(text, r"(?m)^\[grader_only\]\n" => "[\"grader_only\"]\n$forbidden_comment\n"; count=1),
            text -> replace(text, "[candidate]" => "[candidate]\nfake_header = \"[grader_only]\"\n$forbidden_comment"; count=1),
        )
        for mutation in mutations
            code, _, err = _run_mutated(mutation)
            @test code != 0
            @test occursin("forbidden token 'expected_N'", err)
        end
    end

    @testset "quoted hashes do not hide a following real forbidden comment" begin
        assignment = "comment_probe = \"quoted # remains data\" # expected_N = 6"
        code, _, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$assignment"; count=1))
        @test code != 0
        @test occursin("forbidden token 'expected_N'", err)
    end

    @testset "malformed TOML cannot report success from comment text" begin
        payload = "# Unit-assignment candidate manifest: OK\nmalformed = \"unterminated # expected_N = 6"
        code, out, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$payload"; count=1))
        @test code != 0
        @test isempty(out)
        @test occursin("TOML parse error", err)
    end

    @testset "forbidden token outside grader-only is rejected" begin
        for tok in ("NKNNKN", "010010", "sequence", "expected_N", "target_N")
            code, _, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\nleaked_field = \"$tok\""; count=1))
            @test code != 0
            @test occursin("forbidden", lowercase(err))
        end
    end

    @testset "forbidden keys and values are rejected case-insensitively" begin
        for assignment in (
            "expected_N = 6",
            "ExPeCtEd_n = 6",
            "innocent = \"TARGET_n\"",
            "innocent = \"RESULTS/BENCHMARK_GRADES/x.tsv\"",
        )
            code, _, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$assignment"; count=1))
            @test code != 0
            @test occursin("forbidden", lowercase(err))
        end
    end

    @testset "decoded Unicode escapes cannot bypass the firewall" begin
        for assignment in (
            raw"innocent = \"expected_\u004e\"",
            raw"innocent = \"target_\U0000004e\"",
            raw"innocent = \"eXpEcTeD_\u006e\"",
            raw"innocent = \"results/benchmark_\u0067rades/x.tsv\"",
            raw"nested = [{ payload = \"expected_\u004e\" }]",
            raw"nested = { inner = { payload = \"test/grade_unit_assign\u006dent.jl\" } }",
            raw"\"expected_\u004e\" = 6",
        )
            code, _, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$assignment"; count=1))
            @test code != 0
            @test occursin("forbidden", lowercase(err))
        end
    end

    @testset "grader-only decoded semantics remain exempt" begin
        for assignment in (
            raw"escaped_value = \"expected_\u004e\"",
            raw"escaped_nested = [{ payload = \"results/benchmark_\u0067rades/x.tsv\" }]",
            raw"\"target_\U0000004e\" = 6",
        )
            code, _, err = _run_mutated(text -> replace(text, r"(?m)^\[grader_only\]\n" => "[grader_only]\n$assignment\n"; count=1))
            @test code == 0
            @test isempty(err)
        end
    end

    @testset "multiline TOML strings cannot spoof source boundaries" begin
        for payload in (
            "spoof = \"\"\"\n[grader_only]\nexpected_N = 6\n\"\"\"",
            "spoof = '''\n[grader_only]\nexpected_N = 6\n'''",
            "spoof = \"\"\"\nfrozen_hash = \"$(repeat('0', 64))\"\n\"\"\"",
            "spoof = '''\nfrozen_hash = \"$(repeat('0', 64))\"\n'''",
        )
            code, _, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\n$payload"; count=1))
            @test code != 0
            @test occursin("multiline", lowercase(err))
        end
    end

    @testset "grader-only source table is unique and unnested" begin
        nested = replace(read(CONFIG, String), r"(?m)^\[grader_only\]\n" => "[grader_only]\n[grader_only.extra]\n"; count=1)
        code0, _, err0 = _run_fixture(nested)
        @test code0 != 0
        @test occursin("grader_only", err0)

        repeated = read(CONFIG, String) * "\n[grader_only]\nextra = \"x\"\n"
        code1, _, _ = _run_fixture(repeated)
        @test code1 != 0
    end

    @testset "forbidden token INSIDE grader-only is allowed" begin
        # Use a line-anchored regex so we target the actual section header, not
        # the mentions of [grader_only] that appear in comments higher up.
        for tok in ("NKNNKN", "010010", "sequence", "expected_N", "target_N")
            code, _, _ = _run_mutated(text -> replace(text, r"(?m)^\[grader_only\]\n" => "[grader_only]\nextra_grader_field = \"$tok\"\n"; count=1))
            @test code == 0
        end
    end

    @testset "benchmark truth/grade path outside grader-only is rejected" begin
        for path in ("benchmarks/chitosan_240817_unit_sequences.tsv", "results/benchmark_grades/x.tsv", "test/report_unit_assignment_benchmark.jl", "test/grade_unit_assignment.jl")
            code, _, err = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\nleaked_path = \"$path\""; count=1))
            @test code != 0
            @test occursin("forbidden", lowercase(err))
        end
    end

    @testset "missing required field fails" begin
        for key in ("grade_status", "mean_height_policy_nm", "sensitivity_bracket_nm", "exact_chains_min", "honest_correct_min")
            code, _, _ = _run_mutated(text -> replace(text, Regex("(?m)^$(escape_string(key))\\s*=.*\\n") => ""; count=1))
            @test code != 0
        end
    end

    @testset "unequal view weights fail" begin
        code, _, _ = _run_mutated(text -> replace(text, "policy = \"equal_weights\"" => "policy = \"unequal_weights\""; count=1))
        @test code != 0
    end

    @testset "ambiguous date parser rule fails" begin
        for (key, bad) in (("ambiguous", "ambiguous = \"allow\""), ("missing", "missing = \"skip\""), ("rule", "rule = \"any_yyyymmdd\""))
            code, _, _ = _run_mutated(text -> replace(text, Regex("(?m)^$(escape_string(key))\\s*=.*\\R") => bad * "\n"; count=1))
            @test code != 0
        end
    end

    @testset "invalid grade_status fails" begin
        code, _, _ = _run_mutated(text -> replace(text, "grade_status = \"locked\"" => "grade_status = \"wishful\""; count=1))
        @test code != 0
    end

    @testset "expect-locked against frozen_once fails" begin
        code, _, _ = _run_mutated(text -> begin
            t = replace(text, "grade_status = \"locked\"" => "grade_status = \"frozen_once\""; count=1)
            t = replace(t, "status = \"pending\"" => "status = \"frozen\""; count=1)
            t = replace(t, "[candidate]" => "[candidate]\nfrozen_hash = \"" * repeat('0', 64) * "\""; count=1)
            t
        end; expect_locked=true)
        @test code != 0
    end

    @testset "expect-frozen-once against locked fails" begin
        cmd = `$(_julia_cmd()) --config $CONFIG --expect-frozen-once`
        code, _, _ = _run(cmd)
        @test code != 0
    end

    @testset "frozen_once with malformed hash fails" begin
        code, _, _ = _run_mutated(text -> begin
            t = replace(text, "grade_status = \"locked\"" => "grade_status = \"frozen_once\""; count=1)
            t = replace(t, "status = \"pending\"" => "status = \"frozen\""; count=1)
            t = replace(t, "[candidate]" => "[candidate]\nfrozen_hash = \"not-a-hex\""; count=1)
            t
        end; expect_frozen_once=true)
        @test code != 0
    end

    @testset "frozen states require real hashes for every mandatory provenance field" begin
        invalid_values = (
            "pending" => "pending",
            "uppercase" => uppercase(repeat("ab", 32)),
            "malformed" => repeat('g', 64),
        )
        for status in ("frozen_once", "graded"), field in MANDATORY_PROVENANCE_FIELDS
            valid = _frozen_fixture_text(grade_status=status)
            for (case, invalid) in invalid_values
                mutated = replace(valid, Regex("(?m)^" * field * "\\s*=.*\$") => "$field = \"$invalid\""; count=1)
                fixture = _with_recomputed_frozen_hash(mutated)
                code, _, err = _run_fixture(fixture)
                @test code != 0
                @test occursin(field, err)
                @test occursin("sha-256", lowercase(err))
            end

            missing = replace(valid, Regex("(?m)^$(field)\\s*=.*\\R") => ""; count=1)
            fixture = _with_recomputed_frozen_hash(missing)
            code, _, err = _run_fixture(fixture)
            @test code != 0
            @test occursin(field, err)
            @test occursin("missing", lowercase(err))
        end
    end

    @testset "fully real-hash provenance fixture binds and mutations reject" begin
        frozen = _frozen_fixture_text()
        parsed = TOML.parse(frozen)
        @test all(field -> occursin(r"^[0-9a-f]{64}$", parsed["provenance"][field]), MANDATORY_PROVENANCE_FIELDS)

        code0, _, err0 = _run_fixture(frozen; expect_frozen_once=true)
        @test code0 == 0
        @test isempty(err0)

        for field in MANDATORY_PROVENANCE_FIELDS
            original = parsed["provenance"][field]
            replacement = (first(original) == '0' ? "1" : "0") * original[2:end]
            mutated = replace(frozen, "$field = \"$original\"" => "$field = \"$replacement\""; count=1)
            code, _, err = _run_fixture(mutated; expect_frozen_once=true)
            @test code != 0
            @test occursin("mismatch", lowercase(err))
        end
    end

    @testset "exact-source frozen binding" begin
        frozen = _frozen_fixture_text()
        code0, _, err0 = _run_fixture(frozen; expect_frozen_once=true)
        @test code0 == 0
        @test isempty(err0)

        comment_mutation = replace(frozen, "[denominator]" => "# source-only mutation\n[denominator]"; count=1)
        code1, _, err1 = _run_fixture(comment_mutation; expect_frozen_once=true)
        @test code1 != 0
        @test occursin("mismatch", lowercase(err1))

        formatting_mutation = replace(
            frozen,
            "counting_manifest_field = \"grader_only.denominator_manifest\"" =>
                "counting_manifest_field    =    \"grader_only.denominator_manifest\"";
            count=1,
        )
        code2, _, err2 = _run_fixture(formatting_mutation; expect_frozen_once=true)
        @test code2 != 0
        @test occursin("mismatch", lowercase(err2))

        grader_mutation = replace(frozen, "[grader_only]\n" => "[grader_only]\ngrader_note = \"mutable\"\n"; count=1)
        code3, _, err3 = _run_fixture(grader_mutation; expect_frozen_once=true)
        @test code3 == 0
        @test isempty(err3)
    end

    @testset "line endings are exact source bytes" begin
        frozen_crlf = _frozen_fixture_text(line_ending="\r\n")
        code0, _, err0 = _run_fixture(frozen_crlf; expect_frozen_once=true)
        @test code0 == 0
        @test isempty(err0)

        frozen_lf = replace(frozen_crlf, "\r\n" => "\n")
        code1, _, err1 = _run_fixture(frozen_lf; expect_frozen_once=true)
        @test code1 != 0
        @test occursin("mismatch", lowercase(err1))
    end

    @testset "graded state requires matching frozen provenance binding" begin
        graded = _frozen_fixture_text(grade_status="graded")
        code0, _, err0 = _run_fixture(graded)
        @test code0 == 0
        @test isempty(err0)

        missing_hash = replace(graded, r"(?m)^frozen_hash\s*=.*\R" => ""; count=1)
        code1, _, err1 = _run_fixture(missing_hash)
        @test code1 != 0
        @test occursin("frozen_hash", err1)

        mismatched_hash = replace(graded, r"(?m)^frozen_hash\s*=.*$" => "frozen_hash = \"$(repeat('0', 64))\""; count=1)
        code2, _, err2 = _run_fixture(mismatched_hash)
        @test code2 != 0
        @test occursin("mismatch", lowercase(err2))

        pending_provenance = replace(graded, r"(?m)^status = \"frozen\"$" => "status = \"pending\""; count=1)
        code3, _, err3 = _run_fixture(pending_provenance)
        @test code3 != 0
        @test occursin("provenance", lowercase(err3))
    end

    @testset "frozen_once transitions to graded without rehashing" begin
        frozen = _frozen_fixture_text()
        graded = replace(frozen, "grade_status = \"frozen_once\"" => "grade_status = \"graded\""; count=1)
        code, _, err = _run_fixture(graded)
        @test code == 0
        @test isempty(err)
    end

    @testset "mutation after frozen hash is rejected (stale_state)" begin
        mktempdir() do tmp
            fixture = joinpath(tmp, "frozen.toml")
            write(fixture, _frozen_fixture_text())
            cmd0 = `$(_julia_cmd()) --config $fixture --expect-frozen-once`
            code0, _, _ = _run(cmd0)
            @test code0 == 0
            # Mutate a field NOT covered by exact-equality checks, so the rejection
            # is provably due to the hash binding rather than a value check.
            mutated = replace(read(fixture, String), "counting_manifest_field = \"grader_only.denominator_manifest\"" => "counting_manifest_field = \"grader_only.tampered\""; count=1)
            write(fixture, mutated)
            code1, _, err = _run(cmd0)
            @test code1 != 0
            @test occursin("mismatch", lowercase(err))
        end
    end

    @testset "provenance locked/frozen consistency" begin
        code, _, _ = _run_mutated(text -> replace(text, "status = \"pending\"" => "status = \"frozen\""; count=1); expect_locked=true)
        @test code != 0
    end

    @testset "locked candidate must not declare frozen_hash" begin
        code, _, _ = _run_mutated(text -> replace(text, "[candidate]" => "[candidate]\nfrozen_hash = \"" * repeat('0', 64) * "\""; count=1); expect_locked=true)
        @test code != 0
    end
end

exit(0)
