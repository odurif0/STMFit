#!/usr/bin/env julia

using SHA
using TOML
using Test

const ROOT = dirname(@__DIR__)
const EVIDENCE_SOURCE = joinpath(@__DIR__, "lib", "structured_assignment", "evidence.jl")

include(EVIDENCE_SOURCE)
using .StructuredEvidence

file_sha256(path::String) = bytes2hex(sha256(read(path)))

function fixture_artifacts(root::String, name::String)
    directory = joinpath(root, name)
    mkdir(directory)
    declarations = [
        ("config", "model.toml", "[model]\nname = \"fixture\"\n"),
        ("source", "source-b.jl", "source-b-v1\n"),
        ("source", "source-a.jl", "source-a-v1\n"),
        ("input", "input-b.tsv", "key\tvalue\nb\t2\n"),
        ("input", "input-a.tsv", "key\tvalue\na\t1\n"),
        ("output", "output-b.tsv", "file\tvalue\nb\t20\n"),
        ("output", "output-a.tsv", "file\tvalue\na\t10\n"),
        ("expected_rows", "expected.tsv", "file\nscan-a\nscan-b\n"),
        ("log", "command.log", "completed without a terminal claim\n"),
    ]
    specs = ArtifactSpec[]
    for (role, filename, contents) in declarations
        path = joinpath(directory, filename)
        write(path, contents)
        push!(specs, ArtifactSpec(role, relpath(path, root)))
    end
    return specs
end

function artifact_path(root::String, spec::ArtifactSpec)
    return joinpath(root, split(spec.path, '/')...)
end

function publish_fixture(
    root::String,
    evidence_directory::String,
    specs::Vector{ArtifactSpec};
    campaign_id::String="campaign-fixture",
    model_id::String="structured-v2",
    terminal_status::String="PASS",
    reason_codes::Vector{String}=["ok"],
    reasons::Vector{String}=["all declared evidence is complete"],
    generated_at::String="2026-08-08T12:00:00Z",
    exact_command::String="julia --project=$root test/fake_gate.jl --input fixture",
    failpoint::Union{Nothing,Symbol}=nothing,
)
    return publish_evidence_bundle(
        root,
        evidence_directory;
        campaign_id,
        model_id,
        artifacts=specs,
        exact_command,
        generated_at,
        terminal_status,
        reason_codes,
        reasons,
        failpoint,
    )
end

function assert_blocked(result, code::String)
    @test !result.valid
    @test !result.published
    @test result.terminal_status == "BLOCKED"
    @test code in result.reason_codes
end

function independent_command_identity_sha256(root::String, exact_command::String)
    normalized = replace(exact_command, realpath(root) => "\$ROOT")
    material = "stmfit-structured-command-identity-v1\ncommand=$normalized\n"
    return bytes2hex(sha256(codeunits(material)))
end

function concrete_validation(
    root::String,
    receipt_path::String;
    expected_campaign_id::String,
    expected_model_id::String,
    expected_receipt_sha256::String,
    expected_command_identity_sha256::String,
)
    try
        return validate_evidence_bundle(
            root,
            receipt_path;
            expected_campaign_id,
            expected_model_id,
            expected_receipt_sha256,
            expected_command_identity_sha256,
        )
    catch error
        if error isa MethodError || error isa UndefKeywordError
            return validate_evidence_bundle(root, receipt_path)
        end
        rethrow()
    end
end

function validate_published(
    root::String,
    result::EvidenceResult;
    receipt_path::String=result.receipt_path,
    expected_receipt_sha256::String=result.receipt_sha256,
)
    return concrete_validation(
        root,
        receipt_path;
        expected_campaign_id=result.campaign_id,
        expected_model_id=result.model_id,
        expected_receipt_sha256,
        expected_command_identity_sha256=result.command_identity_sha256,
    )
end

function validate_malformed_receipt(root::String, receipt_path::String)
    return concrete_validation(
        root,
        receipt_path;
        expected_campaign_id="campaign-fixture",
        expected_model_id="structured-v2",
        expected_receipt_sha256=repeat('0', 64),
        expected_command_identity_sha256=repeat('0', 64),
    )
end

function rebuilt_index(
    base::ArtifactIndex;
    campaign_id::String=base.campaign_id,
    model_id::String=base.model_id,
    metadata::Dict{String,String}=copy(base.metadata),
)
    return ArtifactIndex(
        base.schema,
        base.schema_version,
        campaign_id,
        model_id,
        base.config_sha256,
        base.source_sha256,
        base.input_sha256,
        base.output_sha256,
        base.expected_rows_sha256,
        base.terminal_status,
        base.selection_eligible,
        base.diagnostic_only,
        copy(base.reason_codes),
        copy(base.reasons),
        copy(base.artifacts),
        metadata,
    )
end

function receipt_for_index(root::String, index_path::String, index::ArtifactIndex)
    index_bytes = artifact_index_bytes(index)
    index_sha256 = bytes2hex(sha256(index_bytes))
    relative_index = replace(relpath(index_path, root), '\\' => '/')
    provisional = GateReceipt(
        GATE_RECEIPT_SCHEMA,
        1,
        index.campaign_id,
        index.model_id,
        index.config_sha256,
        index.source_sha256,
        index.input_sha256,
        index.output_sha256,
        index.expected_rows_sha256,
        index.terminal_status,
        index.selection_eligible,
        index.diagnostic_only,
        copy(index.reason_codes),
        copy(index.reasons),
        relative_index,
        index_sha256,
        artifact_index_canonical_sha256(index),
        repeat('0', 64),
        copy(index.metadata),
    )
    return GateReceipt(
        provisional.schema,
        provisional.schema_version,
        provisional.campaign_id,
        provisional.model_id,
        provisional.config_sha256,
        provisional.source_sha256,
        provisional.input_sha256,
        provisional.output_sha256,
        provisional.expected_rows_sha256,
        provisional.terminal_status,
        provisional.selection_eligible,
        provisional.diagnostic_only,
        provisional.reason_codes,
        provisional.reasons,
        provisional.artifact_index_path,
        provisional.artifact_index_sha256,
        provisional.artifact_index_canonical_sha256,
        bytes2hex(sha256(canonical_gate_receipt_bytes(provisional))),
        provisional.metadata,
    )
end

function install_variant(root::String, evidence::String, index::ArtifactIndex)
    mkdir(evidence)
    index_bytes = artifact_index_bytes(index)
    index_sha256 = bytes2hex(sha256(index_bytes))
    index_path = joinpath(evidence, "artifact-index-$index_sha256.tsv")
    write(index_path, index_bytes)
    receipt_path = joinpath(
        evidence,
        "gate-receipt-$(index.campaign_id)-$(index.model_id).toml",
    )
    receipt = receipt_for_index(root, index_path, index)
    write(receipt_path, gate_receipt_bytes(receipt))
    return index_path, receipt_path, receipt
end

@testset "structured evidence contracts" begin
    @testset "deterministic artifact index and canonical metadata exclusion" begin
        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            kwargs = (
                campaign_id="campaign-fixture",
                model_id="structured-v2",
                exact_command="julia --project=$root test/fake_gate.jl --input fixture",
                generated_at="2026-08-08T12:00:00Z",
                terminal_status="PASS",
                reason_codes=["ok"],
                reasons=["complete fixture"],
            )
            index_a = build_artifact_index(root, specs; kwargs...)
            index_b = build_artifact_index(root, reverse(specs); kwargs...)
            bytes_a = artifact_index_bytes(index_a)
            bytes_b = artifact_index_bytes(index_b)

            @test bytes_a == bytes_b
            @test endswith(String(copy(bytes_a)), "\n")
            @test artifact_index_bytes(parse_artifact_index(bytes_a)) == bytes_a
            @test canonical_artifact_index_bytes(parse_artifact_index(bytes_a)) ==
                  canonical_artifact_index_bytes(index_a)
            @test artifact_index_canonical_sha256(index_a) ==
                  bytes2hex(sha256(canonical_artifact_index_bytes(index_a)))
            @test index_a.config_sha256 ==
                  file_sha256(artifact_path(root, only(filter(s -> s.role == "config", specs))))
            @test index_a.expected_rows_sha256 ==
                  file_sha256(artifact_path(root, only(filter(s -> s.role == "expected_rows", specs))))
            @test index_a.source_sha256 == role_bundle_sha256(index_a, "source")
            @test index_a.input_sha256 == role_bundle_sha256(index_a, "input")
            @test index_a.output_sha256 == role_bundle_sha256(index_a, "output")
            @test [record.path for record in index_a.artifacts] ==
                  sort([spec.path for spec in specs]; by=path ->
                      (only(filter(s -> s.path == path, specs)).role, path))

            changed = with_noncanonical_metadata(
                index_a;
                repository_root="/different/machine/STMFit",
                generated_at="2030-01-02T03:04:05Z",
                exact_command="julia --project=/different/machine/STMFit test/fake_gate.jl --input fixture",
            )
            @test artifact_index_bytes(changed) != bytes_a
            @test canonical_artifact_index_bytes(changed) ==
                  canonical_artifact_index_bytes(index_a)
            @test artifact_index_canonical_sha256(changed) ==
                  artifact_index_canonical_sha256(index_a)
            @test artifact_index_bytes(parse_artifact_index(artifact_index_bytes(changed))) ==
                  artifact_index_bytes(changed)
            @test changed.metadata["repository_root"] == "/different/machine/STMFit"
            @test changed.metadata["generated_at"] == "2030-01-02T03:04:05Z"
            @test changed.metadata["exact_command"] ==
                  "julia --project=/different/machine/STMFit test/fake_gate.jl --input fixture"
        end
    end

    @testset "complete publication round-trips and binds every field" begin
        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            result = publish_fixture(root, evidence, specs)

            @test result.valid
            @test result.published
            @test result.terminal_status == "PASS"
            @test result.reason_codes == ["ok"]
            @test is_v2_selection_eligible(result)
            @test isfile(result.artifact_index_path)
            @test isfile(result.receipt_path)
            @test !islink(result.artifact_index_path)
            @test !islink(result.receipt_path)
            @test basename(result.artifact_index_path) ==
                  "artifact-index-$(result.artifact_index_sha256).tsv"
            @test file_sha256(result.artifact_index_path) == result.artifact_index_sha256
            @test file_sha256(result.receipt_path) == result.receipt_sha256

            index_bytes = read(result.artifact_index_path)
            receipt_bytes = read(result.receipt_path)
            index = parse_artifact_index(index_bytes)
            receipt = parse_gate_receipt(receipt_bytes)
            @test artifact_index_bytes(index) == index_bytes
            @test gate_receipt_bytes(receipt) == receipt_bytes
            @test receipt.campaign_id == index.campaign_id == "campaign-fixture"
            @test receipt.model_id == index.model_id == "structured-v2"
            @test receipt.config_sha256 == index.config_sha256
            @test receipt.source_sha256 == index.source_sha256
            @test receipt.input_sha256 == index.input_sha256
            @test receipt.output_sha256 == index.output_sha256
            @test receipt.expected_rows_sha256 == index.expected_rows_sha256
            @test receipt.terminal_status == index.terminal_status == "PASS"
            @test receipt.reason_codes == index.reason_codes == ["ok"]
            @test receipt.reasons == index.reasons == ["all declared evidence is complete"]
            @test receipt.metadata == index.metadata
            @test receipt.metadata["exact_command"] ==
                  "julia --project=$root test/fake_gate.jl --input fixture"
            @test receipt.artifact_index_sha256 == file_sha256(result.artifact_index_path)
            @test receipt.artifact_index_canonical_sha256 ==
                  artifact_index_canonical_sha256(index)
            @test receipt.canonical_sha256 == gate_receipt_canonical_sha256(receipt)
            @test gate_receipt_canonical_sha256(receipt) ==
                  bytes2hex(sha256(canonical_gate_receipt_bytes(receipt)))

            validation = validate_published(root, result)
            @test validation.valid
            @test !validation.published
            @test validation.terminal_status == "PASS"
            @test validation.artifact_index_path == result.artifact_index_path
            @test validation.artifact_index_sha256 == result.artifact_index_sha256
            @test validation.receipt_sha256 == result.receipt_sha256
            @test is_v2_selection_eligible(validation)

            repeat_result = publish_fixture(root, evidence, reverse(specs))
            @test repeat_result.valid
            @test repeat_result.published
            @test repeat_result.artifact_index_path == result.artifact_index_path
            @test repeat_result.artifact_index_sha256 == result.artifact_index_sha256
            @test read(repeat_result.receipt_path) == receipt_bytes
            @test length(filter(name -> startswith(name, "artifact-index-"), readdir(evidence))) == 1

            changed_receipt = with_noncanonical_metadata(
                receipt;
                repository_root="/other/host/STMFit",
                generated_at="2031-02-03T04:05:06Z",
                exact_command="julia --project=/other/host/STMFit test/fake_gate.jl --input fixture",
            )
            @test gate_receipt_bytes(changed_receipt) != receipt_bytes
            @test canonical_gate_receipt_bytes(changed_receipt) ==
                  canonical_gate_receipt_bytes(receipt)
            @test gate_receipt_canonical_sha256(changed_receipt) ==
                  gate_receipt_canonical_sha256(receipt)
            @test gate_receipt_bytes(parse_gate_receipt(gate_receipt_bytes(changed_receipt))) ==
                  gate_receipt_bytes(changed_receipt)
        end
    end

    @testset "terminal statuses are exact and follow-up is never selectable" begin
        @test Set(TERMINAL_STATUSES) ==
              Set(["PASS", "FAIL", "BLOCKED", "SKIPPED", "FOLLOWUP_REQUIRED"])
        for terminal in TERMINAL_STATUSES
            mktempdir() do root
                specs = fixture_artifacts(root, "artifacts")
                evidence = joinpath(root, "evidence")
                mkdir(evidence)
                result = publish_fixture(
                    root,
                    evidence,
                    specs;
                    terminal_status=terminal,
                    reason_codes=[terminal == "PASS" ? "ok" : "fixture_terminal"],
                    reasons=["terminal fixture for $terminal"],
                )
                @test result.valid
                @test result.published
                @test result.terminal_status == terminal
                receipt = parse_gate_receipt(read(result.receipt_path))
                index = parse_artifact_index(read(result.artifact_index_path))
                @test receipt.diagnostic_only == (terminal == "FOLLOWUP_REQUIRED")
                @test index.diagnostic_only == (terminal == "FOLLOWUP_REQUIRED")
                @test receipt.selection_eligible == (terminal == "PASS")
                @test index.selection_eligible == (terminal == "PASS")
                @test is_v2_selection_eligible(result) == (terminal == "PASS")
                @test !is_v2_selection_eligible(receipt)
                @test !is_v2_selection_eligible(index)
            end
        end

        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            assert_blocked(
                publish_fixture(root, evidence, specs; terminal_status="pass"),
                "invalid_terminal_status",
            )
            @test isempty(readdir(evidence))
            assert_blocked(
                publish_fixture(root, evidence, specs; reason_codes=String[]),
                "invalid_reasons",
            )
            assert_blocked(
                publish_fixture(root, evidence, specs; reasons=["one", "two"]),
                "invalid_reasons",
            )
            assert_blocked(
                publish_fixture(root, evidence, specs; exact_command="julia\nfalse"),
                "invalid_text",
            )
        end
    end

    @testset "current bytes are authoritative for every artifact" begin
        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            published = publish_fixture(root, evidence, specs)
            @test published.valid
            receipt_before = read(published.receipt_path)
            index_before = read(published.artifact_index_path)
            originals = Dict(spec.path => read(artifact_path(root, spec)) for spec in specs)

            for spec in specs
                path = artifact_path(root, spec)
                rm(path)
                assert_blocked(validate_published(root, published),
                    "missing_artifact")
                write(path, originals[spec.path])
                @test validate_published(root, published).valid

                write(path, vcat(originals[spec.path], UInt8('X')))
                assert_blocked(validate_published(root, published),
                    "artifact_size_mismatch")
                write(path, originals[spec.path])

                replacement = copy(originals[spec.path])
                replacement[1] = xor(replacement[1], 0x01)
                write(path, replacement)
                assert_blocked(validate_published(root, published),
                    "artifact_hash_mismatch")
                write(path, originals[spec.path])
                @test validate_published(root, published).valid
            end

            regular = artifact_path(root, first(specs))
            rm(regular)
            mkdir(regular)
            assert_blocked(validate_published(root, published),
                "artifact_not_regular")
            rm(regular)
            write(regular, originals[first(specs).path])

            index_mutation = copy(index_before)
            index_mutation[end - 1] = xor(index_mutation[end - 1], 0x01)
            write(published.artifact_index_path, index_mutation)
            assert_blocked(validate_published(root, published),
                "artifact_index_hash_mismatch")
            write(published.artifact_index_path, index_before)
            @test validate_published(root, published).valid
            @test read(published.receipt_path) == receipt_before
        end
    end

    @testset "path, symlink, duplicate, and schema evasions block publication" begin
        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)

            assert_blocked(
                publish_fixture(root, evidence, vcat(specs, [first(specs)])),
                "duplicate_artifact",
            )
            assert_blocked(
                publish_fixture(root, evidence,
                    vcat(specs[2:end], [ArtifactSpec("config", "artifacts/./model.toml")])),
                "noncanonical_path",
            )
            outside = joinpath(dirname(root), "$(basename(root))-outside.tsv")
            write(outside, "outside\n")
            try
                assert_blocked(
                    publish_fixture(root, evidence,
                        vcat(specs[2:end], [ArtifactSpec("config", "../$(basename(outside))")])),
                    "path_escape",
                )
            finally
                rm(outside; force=true)
            end

            target = artifact_path(root, first(specs))
            link = joinpath(root, "artifact-link.toml")
            symlink(target, link)
            assert_blocked(
                publish_fixture(root, evidence,
                    vcat(specs[2:end], [ArtifactSpec("config", relpath(link, root))])),
                "symlink_rejected",
            )
            rm(link)

            linked_parent = joinpath(root, "linked-parent")
            symlink(dirname(target), linked_parent)
            assert_blocked(
                publish_fixture(root, evidence,
                    vcat(specs[2:end], [ArtifactSpec("config", "linked-parent/model.toml")])),
                "symlink_rejected",
            )
            rm(linked_parent)
            @test isempty(readdir(evidence))
        end

        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            published = publish_fixture(root, evidence, specs)
            @test published.valid
            valid_index = read(published.artifact_index_path, String)
            lines = split(chomp(valid_index), '\n')

            duplicate_path = joinpath(root, "duplicate-index.tsv")
            write(duplicate_path, join(vcat(lines, [lines[2]]), "\n") * "\n")
            assert_blocked(validate_artifact_index(root, duplicate_path), "duplicate_artifact")

            malformed_path = joinpath(root, "malformed-index.tsv")
            malformed = replace(valid_index, "schema\tschema_version" =>
                "wrong_schema\tschema_version"; count=1)
            write(malformed_path, malformed)
            assert_blocked(validate_artifact_index(root, malformed_path), "schema_invalid")

            fake_log = joinpath(root, "fabricated-pass.log")
            write(fake_log, "terminal_status = \"PASS\"\nall checks passed\n")
            assert_blocked(validate_malformed_receipt(root, fake_log), "schema_invalid")

            receipt_text = read(published.receipt_path, String)
            schema_invalid = replace(receipt_text,
                "schema_version = 1" => "schema_version = 2"; count=1)
            write(published.receipt_path, schema_invalid)
            assert_blocked(validate_published(root, published),
                "schema_invalid")
            write(published.receipt_path, receipt_text)

            canonical_invalid = replace(receipt_text,
                r"canonical_sha256 = \"[0-9a-f]{64}\"" =>
                    "canonical_sha256 = \"$(repeat('0', 64))\"";
                count=1)
            write(published.receipt_path, canonical_invalid)
            assert_blocked(validate_published(root, published),
                "canonical_hash_mismatch")
            write(published.receipt_path, receipt_text)
            @test validate_published(root, published).valid
        end
    end

    @testset "publication is content-addressed, atomic, and interruption-safe" begin
        mktempdir() do root
            old_specs = fixture_artifacts(root, "old-artifacts")
            new_specs = fixture_artifacts(root, "new-artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            old = publish_fixture(root, evidence, old_specs)
            @test old.valid
            old_receipt = read(old.receipt_path)
            old_index = read(old.artifact_index_path)
            old_names = Set(readdir(evidence))

            for failpoint in (:before_index_install, :after_index_install,
                              :before_receipt_install)
                interrupted = publish_fixture(
                    root,
                    evidence,
                    new_specs;
                    generated_at="2026-08-08T12:00:01Z",
                    failpoint,
                )
                assert_blocked(interrupted, "interrupted_publication")
                @test read(old.receipt_path) == old_receipt
                @test read(old.artifact_index_path) == old_index
                @test Set(readdir(evidence)) == old_names
                @test validate_published(root, old).valid
                @test !ispath(publication_lock_path(old.receipt_path))
            end

            lock_path = publication_lock_path(old.receipt_path)
            mkdir(lock_path)
            concurrent = publish_fixture(root, evidence, new_specs)
            assert_blocked(concurrent, "concurrent_publication")
            @test read(old.receipt_path) == old_receipt
            @test read(old.artifact_index_path) == old_index
            rm(lock_path)

            outside_target = joinpath(root, "outside-receipt-target.toml")
            write(outside_target, "unchanged\n")
            rm(old.receipt_path)
            symlink(outside_target, old.receipt_path)
            linked = publish_fixture(root, evidence, new_specs)
            assert_blocked(linked, "symlink_rejected")
            @test read(outside_target, String) == "unchanged\n"
            @test islink(old.receipt_path)
        end

        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            interrupted = publish_fixture(
                root,
                evidence,
                specs;
                failpoint=:before_receipt_install,
            )
            assert_blocked(interrupted, "interrupted_publication")
            @test isempty(readdir(evidence))
        end
    end

    @testset "receipt substitution and stale state never inherit PASS text" begin
        mktempdir() do root
            specs_a = fixture_artifacts(root, "artifacts-a")
            specs_b = fixture_artifacts(root, "artifacts-b")
            evidence_a = joinpath(root, "evidence-a")
            evidence_b = joinpath(root, "evidence-b")
            mkdir(evidence_a)
            mkdir(evidence_b)
            result_a = publish_fixture(root, evidence_a, specs_a; campaign_id="campaign-a")
            result_b = publish_fixture(root, evidence_b, specs_b; campaign_id="campaign-b")
            @test result_a.valid
            @test result_b.valid

            receipt_a = read(result_a.receipt_path)
            write(result_a.receipt_path, read(result_b.receipt_path))
            assert_blocked(validate_published(root, result_a),
                "receipt_path_mismatch")
            write(result_a.receipt_path, receipt_a)

            output = artifact_path(root, first(filter(spec -> spec.role == "output", specs_a)))
            output_bytes = read(output)
            write(output, vcat(output_bytes, UInt8('!')))
            success_log = joinpath(evidence_a, "success.log")
            write(success_log, "PASS: downstream publication succeeded\n")
            assert_blocked(validate_published(root, result_a),
                "artifact_size_mismatch")
            assert_blocked(validate_malformed_receipt(root, success_log), "schema_invalid")
            write(output, output_bytes)
            @test validate_published(root, result_a).valid
        end
    end

    @testset "correction regressions bind concrete identity and reconcile commit" begin
        mktempdir() do root
            specs = fixture_artifacts(root, "artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            initial = publish_fixture(root, evidence, specs)
            @test initial.valid
            index = parse_artifact_index(read(initial.artifact_index_path))
            receipt = parse_gate_receipt(read(initial.receipt_path))
            initial_receipt_sha256 = file_sha256(initial.receipt_path)
            expected_command_identity = independent_command_identity_sha256(
                root,
                receipt.metadata["exact_command"],
            )

            exact = concrete_validation(
                root,
                initial.receipt_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=initial_receipt_sha256,
                expected_command_identity_sha256=expected_command_identity,
            )
            @test exact.valid
            @test hasproperty(exact, :identity_verified) && exact.identity_verified
            @test hasproperty(exact, :campaign_id) && exact.campaign_id == "campaign-fixture"
            @test hasproperty(exact, :model_id) && exact.model_id == "structured-v2"
            @test hasproperty(exact, :command_identity_sha256) &&
                  exact.command_identity_sha256 == expected_command_identity
            @test is_v2_selection_eligible(exact)

            missing = validate_evidence_bundle(root, initial.receipt_path)
            @test !missing.valid
            @test "missing_expected_identity" in missing.reason_codes
            @test !is_v2_selection_eligible(missing)
            @test !is_v2_selection_eligible(index)
            @test !is_v2_selection_eligible(receipt)

            changed_command = with_noncanonical_metadata(
                index;
                repository_root=root,
                generated_at=index.metadata["generated_at"],
                exact_command="julia --project=$root attacker-command.jl",
            )
            _, command_path, command_receipt = install_variant(
                root,
                joinpath(root, "command-evidence"),
                changed_command,
            )
            @test command_receipt.canonical_sha256 == receipt.canonical_sha256
            command_full_mismatch = concrete_validation(
                root,
                command_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=initial_receipt_sha256,
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !command_full_mismatch.valid
            @test "receipt_sha256_mismatch" in command_full_mismatch.reason_codes
            command_identity_mismatch = concrete_validation(
                root,
                command_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=file_sha256(command_path),
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !command_identity_mismatch.valid
            @test "command_identity_mismatch" in command_identity_mismatch.reason_codes
            @test !is_v2_selection_eligible(command_identity_mismatch)

            changed_root = with_noncanonical_metadata(
                index;
                repository_root="/attacker/other-root",
                generated_at=index.metadata["generated_at"],
                exact_command=index.metadata["exact_command"],
            )
            _, root_path, root_receipt = install_variant(
                root,
                joinpath(root, "root-evidence"),
                changed_root,
            )
            @test root_receipt.canonical_sha256 == receipt.canonical_sha256
            root_mismatch = concrete_validation(
                root,
                root_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=file_sha256(root_path),
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !root_mismatch.valid
            @test "repository_root_mismatch" in root_mismatch.reason_codes

            changed_time = with_noncanonical_metadata(
                index;
                repository_root=root,
                generated_at="2035-01-02T03:04:05Z",
                exact_command=index.metadata["exact_command"],
            )
            _, time_path, time_receipt = install_variant(
                root,
                joinpath(root, "time-evidence"),
                changed_time,
            )
            @test time_receipt.canonical_sha256 == receipt.canonical_sha256
            time_mismatch = concrete_validation(
                root,
                time_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=initial_receipt_sha256,
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !time_mismatch.valid
            @test "receipt_sha256_mismatch" in time_mismatch.reason_codes

            changed_campaign = rebuilt_index(index; campaign_id="campaign-other")
            _, campaign_path, _ = install_variant(
                root,
                joinpath(root, "campaign-evidence"),
                changed_campaign,
            )
            campaign_mismatch = concrete_validation(
                root,
                campaign_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=file_sha256(campaign_path),
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !campaign_mismatch.valid
            @test "campaign_identity_mismatch" in campaign_mismatch.reason_codes
            @test hasproperty(campaign_mismatch, :campaign_id) &&
                  campaign_mismatch.campaign_id == "campaign-other"

            changed_model = rebuilt_index(index; model_id="structured-v2-other")
            _, model_path, _ = install_variant(
                root,
                joinpath(root, "model-evidence"),
                changed_model,
            )
            model_mismatch = concrete_validation(
                root,
                model_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=file_sha256(model_path),
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !model_mismatch.valid
            @test "model_identity_mismatch" in model_mismatch.reason_codes
            @test hasproperty(model_mismatch, :model_id) &&
                  model_mismatch.model_id == "structured-v2-other"

            wrong_full_hash = concrete_validation(
                root,
                initial.receipt_path;
                expected_campaign_id="campaign-fixture",
                expected_model_id="structured-v2",
                expected_receipt_sha256=repeat('0', 64),
                expected_command_identity_sha256=expected_command_identity,
            )
            @test !wrong_full_hash.valid
            @test "receipt_sha256_mismatch" in wrong_full_hash.reason_codes
        end

        mktempdir() do root_a
            specs_a = fixture_artifacts(root_a, "artifacts")
            evidence_a = joinpath(root_a, "evidence")
            mkdir(evidence_a)
            result_a = publish_fixture(root_a, evidence_a, specs_a)
            receipt_a = parse_gate_receipt(read(result_a.receipt_path))
            mktempdir() do root_b
                specs_b = fixture_artifacts(root_b, "artifacts")
                evidence_b = joinpath(root_b, "evidence")
                mkdir(evidence_b)
                result_b = publish_fixture(root_b, evidence_b, specs_b)
                receipt_b = parse_gate_receipt(read(result_b.receipt_path))
                @test receipt_a.canonical_sha256 == receipt_b.canonical_sha256
                command_a = independent_command_identity_sha256(
                    root_a,
                    receipt_a.metadata["exact_command"],
                )
                command_b = independent_command_identity_sha256(
                    root_b,
                    receipt_b.metadata["exact_command"],
                )
                @test command_a == command_b
                validation_a = concrete_validation(
                    root_a,
                    result_a.receipt_path;
                    expected_campaign_id="campaign-fixture",
                    expected_model_id="structured-v2",
                    expected_receipt_sha256=file_sha256(result_a.receipt_path),
                    expected_command_identity_sha256=command_a,
                )
                validation_b = concrete_validation(
                    root_b,
                    result_b.receipt_path;
                    expected_campaign_id="campaign-fixture",
                    expected_model_id="structured-v2",
                    expected_receipt_sha256=file_sha256(result_b.receipt_path),
                    expected_command_identity_sha256=command_b,
                )
                @test hasproperty(validation_a, :identity_verified) && validation_a.identity_verified
                @test hasproperty(validation_b, :identity_verified) && validation_b.identity_verified
            end
        end

        mktempdir() do root
            old_specs = fixture_artifacts(root, "old-artifacts")
            new_specs = fixture_artifacts(root, "new-artifacts")
            evidence = joinpath(root, "evidence")
            mkdir(evidence)
            old = publish_fixture(root, evidence, old_specs)
            old_receipt = read(old.receipt_path)
            old_index = read(old.artifact_index_path)

            precommit = publish_fixture(
                root,
                evidence,
                new_specs;
                exact_command="julia --project=$root replacement.jl",
                generated_at="2026-08-08T12:00:01Z",
                failpoint=:index_directory_fsync,
            )
            @test !precommit.valid
            @test "fsync_failed" in precommit.reason_codes
            @test read(old.receipt_path) == old_receipt
            @test read(old.artifact_index_path) == old_index

            postcommit = publish_fixture(
                root,
                evidence,
                new_specs;
                exact_command="julia --project=$root replacement.jl",
                generated_at="2026-08-08T12:00:01Z",
                failpoint=:receipt_directory_fsync,
            )
            @test postcommit.valid && postcommit.published &&
                  postcommit.terminal_status == "PASS"
            @test hasproperty(postcommit, :identity_verified) && postcommit.identity_verified
            @test hasproperty(postcommit, :durability_status) &&
                  postcommit.durability_status == "WARNING"
            @test hasproperty(postcommit, :warnings) && !isempty(postcommit.warnings)
            @test read(old.receipt_path) != old_receipt
            @test read(old.artifact_index_path) == old_index
            if postcommit.valid && hasproperty(postcommit, :command_identity_sha256)
                current = concrete_validation(
                    root,
                    postcommit.receipt_path;
                    expected_campaign_id=postcommit.campaign_id,
                    expected_model_id=postcommit.model_id,
                    expected_receipt_sha256=postcommit.receipt_sha256,
                    expected_command_identity_sha256=postcommit.command_identity_sha256,
                )
                @test current.valid && current.identity_verified &&
                      is_v2_selection_eligible(current)
            else
                @test false
            end
            @test !ispath(publication_lock_path(old.receipt_path))
        end
    end
end
