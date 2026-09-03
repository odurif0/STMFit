#!/usr/bin/env julia

using Printf
using SHA
using Test

const ROOT = dirname(@__DIR__)
const ADAPTER_SOURCE = joinpath(
    @__DIR__, "lib", "structured_assignment", "champion_adapter.jl")
const FROZEN_SOURCE = joinpath(
    ROOT, "results", "unit_assignment", "best_labelfree_cc_soft_20260802.tsv")
const FROZEN_SHA256 =
    "86c7b0985413cf4ceb61ee903ed9bd27e375b5cc2ed314fa1b15b9648dd4ecf5"
const FROZEN_BYTES = 27_451
const FROZEN_ROWS = 892
const FROZEN_KEY_ORDER_SHA256 =
    "ca21d9004b3ba95711d143359999b223fcbf6e69601454f5f9d1ccc04df6e1ba"
const PROTECTED_HASHES = Dict(
    "docs/src/journal.md" =>
        "8ec2bd10703ccad097ddb378d2b705248fdfa1a42ca2e2a8889e1472dc27ab29",
    "hpc/v24_pair_sweep.sbatch" =>
        "d2d91f2218439b54c25b141a9ffdf42b27e7ffd5424dad50f51f9bcad15a62c5",
    "test/detect_missed_lobes.py" =>
        "b6a46945ffb58c81e47ff9ca72238e2a5022aba71dd545fd6938507d7a8ffc4b",
    "test/refit_missed_lobes.py" =>
        "e37ba297e57753b748a9818b48ce65c03d63ea615a1c1f14a01006f9c8cce45e",
)

sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))
file_sha256(path::AbstractString) = sha256_hex(read(path))

if !isfile(ADAPTER_SOURCE)
    @testset "structured champion adapter RED: implementation is required" begin
        @test isfile(ADAPTER_SOURCE)
    end
    exit(1)
end

module ChampionAdapterUnderTest
include(Main.ADAPTER_SOURCE)
end

const CA = ChampionAdapterUnderTest

function rejection_code(f::Function)
    try
        f()
    catch err
        @test err isa CA.ChampionAdapterError
        return err isa CA.ChampionAdapterError ? err.code : :wrong_exception
    end
    @test false
    return :not_rejected
end

function frozen_lines()
    text = String(read(FROZEN_SOURCE))
    @test isvalid(text)
    @test endswith(text, "\n")
    return String.(split(text[1:prevind(text, lastindex(text))], '\n'; keepempty=true))
end

function lines_bytes(lines::AbstractVector{<:AbstractString})
    return collect(codeunits(join(lines, "\n") * "\n"))
end

function replace_field(
    lines::AbstractVector{<:AbstractString},
    row::Int,
    column::Int,
    value::String,
)
    changed = String.(lines)
    fields = split(changed[row], '\t'; keepempty=true)
    fields[column] = value
    changed[row] = join(fields, '\t')
    return changed
end

function stream_digest(stream)
    io = IOBuffer()
    for record in stream
        write(io, record.key.file, '\t', string(record.key.lobe), '\t',
              string(record.unary.predicted), '\t', record.unary.confidence_text,
              '\t', bitstring(record.unary.probability_1), '\n')
    end
    return sha256_hex(take!(io))
end

function child_process_snapshot()
    path = "/proc/$(getpid())/task/$(getpid())/children"
    return isfile(path) ? strip(read(path, String)) : "unavailable"
end

@testset "immutable frozen-champion unary adapter" begin
    @testset "typed real-source stream and exact round trip" begin
        @test file_sha256(FROZEN_SOURCE) == FROZEN_SHA256
        children_before = child_process_snapshot()
        stream = CA.load_frozen_champion_unaries()
        children_after = child_process_snapshot()

        @test stream isa CA.BaselineUnaryStream
        @test stream isa AbstractVector{CA.BaselineUnaryRecord}
        @test !ismutabletype(typeof(stream))
        @test !ismutabletype(CA.FrozenChampionSource)
        @test !ismutabletype(CA.BaselineUnaryKey)
        @test !ismutabletype(CA.BaselineUnaryEvidence)
        @test !ismutabletype(CA.BaselineUnaryRecord)
        @test length(stream) == FROZEN_ROWS
        @test size(stream) == (FROZEN_ROWS,)
        @test axes(stream) == (Base.OneTo(FROZEN_ROWS),)
        @test stream.source.path == realpath(FROZEN_SOURCE)
        @test stream.source.sha256 == FROZEN_SHA256
        @test stream.source.byte_count == FROZEN_BYTES
        @test stream.source.row_count == FROZEN_ROWS
        @test stream.source.key_order_sha256 == FROZEN_KEY_ORDER_SHA256
        @test children_after == children_before

        lines = frozen_lines()
        @test length(lines) == FROZEN_ROWS + 1
        @test lines[1] == "file\tlobe\tpredicted\tconfidence"
        seen = Set{Tuple{String,Int}}()
        zero_confidence_predictions = Int[]
        for (index, source_line) in enumerate(lines[2:end])
            fields = split(source_line, '\t'; keepempty=true)
            @test length(fields) == 4
            record = stream[index]
            source_key = (fields[1], parse(Int, fields[2]))
            @test (record.key.file, record.key.lobe) == source_key
            @test !(source_key in seen)
            push!(seen, source_key)
            @test string(record.unary.predicted) == fields[3]
            @test record.unary.confidence_text == fields[4]
            @test @sprintf("%.8f", record.unary.confidence) == fields[4]
            @test record.unary.confidence == parse(Float64, fields[4])
            probability = record.unary.predicted == 1 ?
                0.5 + record.unary.confidence / 2 :
                0.5 - record.unary.confidence / 2
            @test isequal(record.unary.probability_1, probability)
            @test 0.0 <= record.unary.probability_1 <= 1.0
            if iszero(record.unary.confidence)
                push!(zero_confidence_predictions, record.unary.predicted)
                @test record.unary.probability_1 == 0.5
            end
        end
        @test length(seen) == FROZEN_ROWS
        @test !isempty(zero_confidence_predictions)
        @test all(predicted -> predicted in (0, 1), zero_confidence_predictions)
        @test file_sha256(FROZEN_SOURCE) == FROZEN_SHA256
    end

    @testset "join metadata remains separate from unary evidence" begin
        @test fieldnames(CA.BaselineUnaryKey) == (:file, :lobe)
        @test fieldnames(CA.BaselineUnaryEvidence) ==
              (:predicted, :confidence, :confidence_text, :probability_1)
        @test fieldnames(CA.BaselineUnaryRecord) == (:key, :unary)
        for confidence in (0.0, 0.25, 0.5, 1.0)
            text = @sprintf("%.8f", confidence)
            for predicted in (0, 1)
                first_unary = CA._unary_evidence(predicted, confidence, text)
                second_unary = CA._unary_evidence(predicted, confidence, text)
                @test first_unary == second_unary
                @test first_unary.probability_1 ==
                      (predicted == 1 ? 0.5 + confidence / 2 :
                                        0.5 - confidence / 2)
            end
        end
        @test !applicable(CA._unary_evidence, "scan.sxm", 1, 0, 1.0, "1.00000000")
    end

    @testset "deterministic replay" begin
        first_stream = CA.load_frozen_champion_unaries()
        second_stream = CA.load_frozen_champion_unaries(FROZEN_SOURCE)
        @test length(first_stream) == length(second_stream) == FROZEN_ROWS
        @test stream_digest(first_stream) == stream_digest(second_stream)
        @test first_stream.source == second_stream.source
        @test collect(first_stream) == collect(second_stream)
    end

    @testset "strict row and key integrity" begin
        lines = frozen_lines()
        @test length(CA._parse_frozen_rows(lines_bytes(lines))) == FROZEN_ROWS

        duplicate = copy(lines)
        duplicate[3] = duplicate[2]
        @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(duplicate))) ==
              :duplicate_key

        missing = copy(lines)
        deleteat!(missing, 3)
        @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(missing))) ==
              :row_count_mismatch

        extra = copy(lines)
        push!(extra, "fixture_scan.sxm\t1\t0\t1.00000000")
        @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(extra))) ==
              :row_count_mismatch

        reordered = copy(lines)
        reordered[2], reordered[3] = reordered[3], reordered[2]
        @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(reordered))) ==
              :key_order_mismatch

        for header in (
            "lobe\tfile\tpredicted\tconfidence",
            "file\tlobe\tpredicted",
            "file\tlobe\tpredicted\tconfidence\textra",
            "file\tlobe\tpredicted\tpredicted",
        )
            malformed = copy(lines)
            malformed[1] = header
            @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(malformed))) ==
                  :invalid_columns
        end

        missing_cell = copy(lines)
        missing_cell[2] = join(split(missing_cell[2], '\t')[1:3], '\t')
        @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(missing_cell))) ==
              :invalid_row_columns
        extra_cell = copy(lines)
        extra_cell[2] *= "\textra"
        @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(extra_cell))) ==
              :invalid_row_columns
    end

    @testset "strict prediction and confidence semantics" begin
        lines = frozen_lines()
        for value in ("2", "?", "-1", "", "01", " 1")
            fixture = replace_field(lines, 2, 3, value)
            @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(fixture))) ==
                  :invalid_prediction
        end
        for value in ("NaN", "Inf", "+Inf", "-Inf")
            fixture = replace_field(lines, 2, 4, value)
            @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(fixture))) ==
                  :nonfinite_confidence
        end
        for value in ("1.00000001", "-0.00000000", "-0.10000000", "2.00000000")
            fixture = replace_field(lines, 2, 4, value)
            @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(fixture))) ==
                  :confidence_out_of_range
        end
        for value in ("", "abc", "0.5", ".50000000", " 0.50000000", "0.50000000 ")
            fixture = replace_field(lines, 2, 4, value)
            @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(fixture))) ==
                  :malformed_confidence
        end

        invalid_utf8 = lines_bytes(lines)
        invalid_utf8[findfirst(==(UInt8('2')), invalid_utf8)] = 0xff
        @test rejection_code(() -> CA._parse_frozen_rows(invalid_utf8)) ==
              :invalid_encoding
    end

    @testset "forbidden declarations fail before row materialization" begin
        lines = frozen_lines()
        for name in ("benchmark", "truth", "grade", "grader")
            fixture = copy(lines)
            fixture[1] = "file\tlobe\tpredicted\tconfidence\t$name"
            @test rejection_code(() -> CA._parse_frozen_rows(lines_bytes(fixture))) ==
                  :forbidden_declaration
        end
    end

    @testset "path identity, symlinks, and source mutation" begin
        mktempdir() do temporary
            copy_path = joinpath(temporary, "copy.tsv")
            write(copy_path, read(FROZEN_SOURCE))
            @test rejection_code(() -> CA.load_frozen_champion_unaries(copy_path)) ==
                  :source_path_substitution

            link_path = joinpath(temporary, "link.tsv")
            symlink(FROZEN_SOURCE, link_path)
            @test rejection_code(() -> CA.load_frozen_champion_unaries(link_path)) ==
                  :symlink_rejected

            real_parent = joinpath(temporary, "real-parent")
            mkdir(real_parent)
            parent_copy = joinpath(real_parent, "source.tsv")
            write(parent_copy, read(FROZEN_SOURCE))
            linked_parent = joinpath(temporary, "linked-parent")
            symlink(real_parent, linked_parent)
            @test rejection_code(() -> CA.load_frozen_champion_unaries(
                joinpath(linked_parent, "source.tsv"))) == :symlink_rejected

            dotdot_path = joinpath(ROOT, "results", "unit_assignment", "..",
                                   "unit_assignment", basename(FROZEN_SOURCE))
            @test rejection_code(() -> CA.load_frozen_champion_unaries(dotdot_path)) ==
                  :source_path_substitution

            snapshot = CA._read_source_snapshot(copy_path, FROZEN_SHA256)
            changed = read(copy_path)
            changed[end - 1] = xor(changed[end - 1], 0x01)
            write(copy_path, changed)
            before_verify_paths = sort(readdir(temporary))
            @test rejection_code(() -> CA._verify_source_snapshot(snapshot)) ==
                  :source_changed
            @test sort(readdir(temporary)) == before_verify_paths
        end
        @test file_sha256(FROZEN_SOURCE) == FROZEN_SHA256
    end

    @testset "hash gate precedes parsing and return gate rechecks bytes" begin
        mktempdir() do temporary
            path = joinpath(temporary, "malformed.tsv")
            lines = replace_field(frozen_lines(), 2, 3, "2")
            bytes = lines_bytes(lines)
            write(path, bytes)
            @test rejection_code(() -> CA._load_source_unaries(
                path, path, FROZEN_SHA256)) == :source_hash_mismatch
            @test rejection_code(() -> CA._load_source_unaries(
                path, path, sha256_hex(bytes))) == :invalid_prediction
            @test readdir(temporary) == ["malformed.tsv"]
        end
    end

    @testset "interruption cannot leave partial output or adapter state" begin
        mktempdir() do temporary
            paths_before = readdir(temporary)
            interrupted = try
                CA.load_frozen_champion_unaries()
                throw(InterruptException())
            catch err
                err
            end
            @test interrupted isa InterruptException
            @test readdir(temporary) == paths_before
            @test file_sha256(FROZEN_SOURCE) == FROZEN_SHA256
        end
    end

    @testset "stdlib-only source and no child-process surface" begin
        source_text = read(ADAPTER_SOURCE, String)
        lowered = lowercase(source_text)
        for forbidden in (
            "benchmark", "truth", "grade", "grader", "labels", "label_path",
            "known_label", "python",
            "expected_n", "target_n", "class_count", "composition",
            "top_k", "nknnkn", "010010", "101101", "build_cc_soft_champion",
        )
            @test !occursin(forbidden, lowered)
        end
        for process_api in ("run(", "cmd(", "pipeline(", "base.julia_cmd", "`julia")
            @test !occursin(process_api, lowered)
        end
        for write_api in ("write(", "mkpath(", "mkdir(", "rm(", "mv(", "cp(")
            @test !occursin(write_api, lowered)
        end
        imports = [strip(match.captures[2]) for match in eachmatch(
            r"(?m)^\s*(using|import)\s+([^\n]+)$", source_text)]
        @test Set(imports) == Set(["Printf", "SHA"])
    end

    @testset "protected baseline and frozen source remain exact" begin
        for (relative_path, digest) in PROTECTED_HASHES
            @test file_sha256(joinpath(ROOT, relative_path)) == digest
        end
        @test file_sha256(FROZEN_SOURCE) == FROZEN_SHA256
        @test filesize(FROZEN_SOURCE) == FROZEN_BYTES
    end
end
