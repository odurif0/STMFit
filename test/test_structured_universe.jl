#!/usr/bin/env julia

using SHA
using Statistics
using Test
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const UNIVERSE_SOURCE = joinpath(@__DIR__, "lib", "structured_assignment", "universe.jl")
const UNIVERSE_CLI = joinpath(@__DIR__, "build_structured_scan_universe.jl")
const FIREWALL_SOURCE = joinpath(@__DIR__, "lib", "structured_assignment", "firewall.jl")

if !isfile(UNIVERSE_SOURCE)
    @testset "structured universe contract source" begin
        @test isfile(UNIVERSE_SOURCE)
    end
    exit(1)
end

include(UNIVERSE_SOURCE)
using .StructuredUniverse

const REPLAY_SHA256 = Ref("")

_sha(bytes) = bytes2hex(sha256(bytes))
_sha_file(path) = _sha(read(path))

function _write(path::String, text::String)
    mkpath(dirname(path))
    write(path, text)
    return path
end

function _feature_text(; rows=nothing)
    body = rows === nothing ? [
        "20240229_b.sxm\t2\t2\t2.2",
        "20240101_b.sxm\t1\t1\t1.1",
        "20240101_a.sxm\t3\t3\t1.3",
        "20240229_b.sxm\t2\t1\t2.1",
        "20240101_a.sxm\t3\t1\t1.1",
        "20240101_a.sxm\t3\t2\t1.2",
    ] : rows
    return "file\tN\tlobe\tamplitude\n" * join(body, "\n") * "\n"
end

function _fixture()
    root = mktempdir()
    for relative in (
        "config/unit_assignment_structured_candidate.toml",
        "config/unit_assignment_structured_model.toml",
        "config/chitosan.toml",
        "test/extract_lobe_patches.jl",
        "test/extract_lobe_patches_bwd.jl",
    )
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(joinpath(REPOSITORY_ROOT, relative), destination)
    end
    features = _write(joinpath(root, "inputs", "features.tsv"), _feature_text())
    return (
        root=root,
        features=features,
        feature_sha256=_sha_file(features),
        candidate="config/unit_assignment_structured_candidate.toml",
        model="config/unit_assignment_structured_model.toml",
    )
end

function _freeze(fixture; kwargs...)
    return freeze_universe(
        fixture.root;
        features=relpath(fixture.features, fixture.root),
        feature_sha256=fixture.feature_sha256,
        candidate_config=fixture.candidate,
        model_config=fixture.model,
        source_paths=[UNIVERSE_SOURCE, UNIVERSE_CLI, FIREWALL_SOURCE],
        kwargs...,
    )
end

function _code(action)
    try
        action()
    catch error
        error isa UniverseError || rethrow()
        return error.code
    end
    return :accepted
end

function _table(path::String)
    lines = split(chomp(read(path, String)), '\n')
    header = split(first(lines), '\t'; keepempty=true)
    rows = [Dict(header[index] => value for (index, value) in enumerate(
                split(line, '\t'; keepempty=true))) for line in lines[2:end]]
    return header, rows
end

function _bundle_digest(path::String)
    io = IOBuffer()
    for name in sort(readdir(path))
        bytes = read(joinpath(path, name))
        println(io, name, '\t', _sha(bytes), '\t', length(bytes))
    end
    return _sha(take!(io))
end

function _patch_text(keys; value="0.25")
    rows = ["$(key.file)\t$(key.lobe)\t$value" for key in keys]
    return "file\tlobe\tres_p001\n" * join(rows, "\n") * "\n"
end

function _run_cli(arguments, stdout_path::String, stderr_path::String)
    command = `$(Base.julia_cmd()) --project=$(REPOSITORY_ROOT) $(UNIVERSE_CLI) $arguments`
    process = open(stdout_path, "w") do stdout_io
        open(stderr_path, "w") do stderr_io
            run(pipeline(command; stdout=stdout_io, stderr=stderr_io); wait=false)
        end
    end
    wait(process)
    return process.exitcode
end

@testset "structured immutable universe" begin
    @testset "dates, rows, and immutable keys" begin
        fixture = _fixture()
        frozen = _freeze(fixture)

        @test parse_scan_date("20240229_scan.sxm") == "20240229"
        @test _code(() -> parse_scan_date("scan.sxm")) == :invalid_scan_name
        @test _code(() -> parse_scan_date("20230229_scan.sxm")) == :invalid_scan_date
        @test _code(() -> parse_scan_date("20240101_20240202_scan.sxm")) == :ambiguous_scan_date
        @test _code(() -> parse_scan_date("folder/20240101_scan.sxm")) == :noncanonical_scan_name
        @test _code(() -> parse_scan_date("20240101_SCAN.SXM")) == :noncanonical_scan_name

        @test Tuple(row.file for row in frozen.scans) ==
              ("20240101_a.sxm", "20240101_b.sxm", "20240229_b.sxm")
        @test Tuple(row.date for row in frozen.scans) ==
              ("20240101", "20240101", "20240229")
        @test Tuple(row.lobe_count for row in frozen.scans) == (3, 1, 2)
        @test frozen.keys == (
            LobeKey("20240101_a.sxm", 1),
            LobeKey("20240101_a.sxm", 2),
            LobeKey("20240101_a.sxm", 3),
            LobeKey("20240101_b.sxm", 1),
            LobeKey("20240229_b.sxm", 1),
            LobeKey("20240229_b.sxm", 2),
        )
        @test frozen.feature_sha256 == fixture.feature_sha256
        @test frozen.candidate_config_sha256 == CANDIDATE_CONFIG_SHA256
        @test frozen.model_config_sha256 == MODEL_CONFIG_SHA256
        @test frozen.forward_producer_sha256 == FORWARD_PRODUCER_SHA256
        @test frozen.backward_producer_sha256 == BACKWARD_PRODUCER_SHA256
        @test frozen.grid_sha256 == GRID_SHA256
        @test validate_binding(frozen) === nothing
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "strict feature table failures" begin
        cases = Dict{Symbol,String}(
            :duplicate_row => _feature_text(rows=[
                "20240101_a.sxm\t1\t1\t1.0",
                "20240101_a.sxm\t1\t1\t1.0",
                "20240202_b.sxm\t1\t1\t2.0",
            ]),
            :duplicate_key => _feature_text(rows=[
                "20240101_a.sxm\t1\t1\t1.0",
                "20240101_a.sxm\t1\t1\t1.1",
                "20240202_b.sxm\t1\t1\t2.0",
            ]),
            :lobe_key_mismatch => _feature_text(rows=[
                "20240101_a.sxm\t2\t1\t1.0",
                "20240202_b.sxm\t1\t1\t2.0",
            ]),
            :lobe_count_mismatch => _feature_text(rows=[
                "20240101_a.sxm\t2\t1\t1.0",
                "20240101_a.sxm\t3\t2\t1.1",
                "20240202_b.sxm\t1\t1\t2.0",
            ]),
            :invalid_integer => replace(_feature_text(), "\t3\t3\t" => "\t3.0\t3\t"),
            :noncanonical_integer => replace(_feature_text(), "\t3\t3\t" => "\t03\t3\t"),
            :row_width_mismatch => replace(_feature_text(), "20240101_b.sxm\t1\t1\t1.1" => "20240101_b.sxm\t1\t1"),
            :duplicate_column => replace(_feature_text(), "file\tN\tlobe\tamplitude" => "file\tN\tlobe\tlobe"),
            :missing_column => replace(_feature_text(), "file\tN\tlobe\tamplitude" => "file\tN\tunit\tamplitude"),
            :missing_final_lf => chop(_feature_text()),
            :carriage_return => replace(_feature_text(), "\n" => "\r\n"),
        )
        for (wanted, text) in cases
            fixture = _fixture()
            write(fixture.features, text)
            current = merge(fixture, (feature_sha256=_sha_file(fixture.features),))
            @test _code(() -> _freeze(current)) == wanted
            rm(fixture.root; recursive=true, force=true)
        end

        fixture = _fixture()
        @test _code(() -> _freeze(fixture; feature_sha256=repeat("0", 64))) ==
              :feature_hash_mismatch
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "explicit negative text fixtures" begin
        negative_values = [
            "benchmark", "grader", "truth", "expected_count", "target_count",
            "composition", "top_k", "class_count", "label",
        ]
        for value in negative_values
            fixture = _fixture()
            text = replace(_feature_text(), "amplitude" => value)
            write(fixture.features, text)
            current = merge(fixture, (feature_sha256=_sha_file(fixture.features),))
            @test _code(() -> _freeze(current)) == :forbidden_input
            rm(fixture.root; recursive=true, force=true)
        end
    end

    @testset "mutation and declared hash binding" begin
        for mutation in (
            text -> replace(text, "20240101_a.sxm" => "20240101_c.sxm"),
            text -> replace(text, "20240101_b.sxm\t1\t1\t1.1\n" => ""),
            text -> text * "20240303_c.sxm\t1\t1\t3.0\n",
            text -> replace(text, "20240101_a.sxm" => "20240101_20240202_a.sxm"),
            text -> replace(text, "1.3" => "1.4"),
        )
            fixture = _fixture()
            frozen = _freeze(fixture)
            write(fixture.features, mutation(read(fixture.features, String)))
            @test _code(() -> validate_binding(frozen)) == :dependency_changed
            out = joinpath(fixture.root, "out")
            @test _code(() -> publish_bundle(frozen, out)) == :dependency_changed
            @test !ispath(out)
            rm(fixture.root; recursive=true, force=true)
        end

        fixture = _fixture()
        @test _code(() -> _freeze(fixture; model_config_sha256=repeat("0", 64))) ==
              :config_hash_mismatch
        rm(fixture.root; recursive=true, force=true)

        fixture = _fixture()
        producer = joinpath(fixture.root, "test", "extract_lobe_patches.jl")
        write(producer, read(producer, String) * "\n")
        @test _code(() -> _freeze(fixture)) == :source_hash_mismatch
        rm(fixture.root; recursive=true, force=true)

        fixture = _fixture()
        frozen = _freeze(fixture)
        model = joinpath(fixture.root, fixture.model)
        write(model, read(model, String) * "\n")
        @test _code(() -> validate_binding(frozen)) == :dependency_changed
        rm(fixture.root; recursive=true, force=true)

        fixture = _fixture()
        candidate = joinpath(fixture.root, fixture.candidate)
        changed = replace(read(candidate, String), FORWARD_PRODUCER_SHA256 => repeat("0", 64))
        write(candidate, changed)
        @test _code(() -> _freeze(fixture;
            candidate_config_sha256=_sha_file(candidate))) == :declaration_mismatch
        rm(fixture.root; recursive=true, force=true)

        fixture = _fixture()
        model = joinpath(fixture.root, fixture.model)
        changed_model = replace(read(model, String), GRID_SHA256 => repeat("0", 64))
        write(model, changed_model)
        candidate = joinpath(fixture.root, fixture.candidate)
        changed_candidate = replace(
            read(candidate, String), MODEL_CONFIG_SHA256 => _sha_file(model))
        write(candidate, changed_candidate)
        @test _code(() -> _freeze(fixture;
            candidate_config_sha256=_sha_file(candidate),
            model_config_sha256=_sha_file(model))) == :declaration_mismatch
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "canonical bundle, work rows, and replay" begin
        fixture = _fixture()
        frozen = _freeze(fixture)
        first_out = joinpath(fixture.root, "bundle-one")
        second_out = joinpath(fixture.root, "bundle-two")
        first_hashes = publish_bundle(frozen, first_out)
        second_hashes = publish_bundle(frozen, second_out)
        names = [
            "folds.tsv", "patch_keys.tsv", "perturbations.tsv", "receipt.toml",
            "scan_universe.tsv", "seeds.tsv", "shards.tsv",
        ]
        @test sort(readdir(first_out)) == names
        @test sort(readdir(second_out)) == names
        @test first_hashes == second_hashes
        @test all(name -> read(joinpath(first_out, name)) == read(joinpath(second_out, name)), names)
        @test validate_bundle(frozen, first_out) == first_hashes
        @test publish_bundle(frozen, first_out) == first_hashes

        _, scan_rows = _table(joinpath(first_out, "scan_universe.tsv"))
        _, key_rows = _table(joinpath(first_out, "patch_keys.tsv"))
        _, fold_rows = _table(joinpath(first_out, "folds.tsv"))
        _, seed_rows = _table(joinpath(first_out, "seeds.tsv"))
        _, perturbation_rows = _table(joinpath(first_out, "perturbations.tsv"))
        _, shard_rows = _table(joinpath(first_out, "shards.tsv"))
        dates = unique(row.date for row in frozen.scans)
        @test length(scan_rows) == length(frozen.scans)
        @test length(key_rows) == length(frozen.keys)
        @test length(fold_rows) == length(dates) * length(frozen.scans)
        seed_count_per_fold = 2length(frozen.bootstrap_seeds) +
                              length(frozen.shuffle_seeds) + length(frozen.restart_seeds)
        @test length(seed_rows) == length(dates) * seed_count_per_fold
        perturbation_count_per_scan = 1 + length(frozen.views) + 4 + 2 +
                                      length(frozen.bootstrap_seeds) +
                                      length(frozen.restart_seeds)
        @test length(perturbation_rows) == length(frozen.scans) * perturbation_count_per_scan
        @test length(shard_rows) == length(fold_rows) + length(seed_rows) +
                                    length(perturbation_rows)
        @test length(unique(row["shard_id"] for row in shard_rows)) == length(shard_rows)
        @test all(row -> row["feature_sha256"] == frozen.feature_sha256, scan_rows)
        @test all(row -> row["grid_sha256"] == GRID_SHA256, scan_rows)
        @test count(row -> row["perturbation"] == "noise", perturbation_rows) ==
              length(frozen.scans) * length(frozen.bootstrap_seeds)
        @test count(row -> row["perturbation"] == "seed_restart", perturbation_rows) ==
              length(frozen.scans) * length(frozen.restart_seeds)
        @test any(row -> row["perturbation"] == "shift_t_plus" &&
                         row["delta_t_nm"] == "0.02", perturbation_rows)
        @test any(row -> row["perturbation"] == "contrast_low" &&
                         row["factor"] == "0.95", perturbation_rows)

        receipt = TOML.parsefile(joinpath(first_out, "receipt.toml"))
        @test receipt["feature_sha256"] == frozen.feature_sha256
        @test receipt["candidate_config_sha256"] == CANDIDATE_CONFIG_SHA256
        @test receipt["model_config_sha256"] == MODEL_CONFIG_SHA256
        @test receipt["source_sha256"] == frozen.source_sha256
        @test receipt["forward_producer_sha256"] == FORWARD_PRODUCER_SHA256
        @test receipt["backward_producer_sha256"] == BACKWARD_PRODUCER_SHA256
        @test receipt["grid_sha256"] == GRID_SHA256
        @test Set(keys(receipt["artifacts"])) == Set(filter(!=("receipt.toml"), names))

        negative_output_values = [
            "benchmark", "grader", "truth", "expected_count", "target_count",
            "composition", "top_k", "class_count", "label",
        ]
        for name in names
            text = lowercase(read(joinpath(first_out, name), String))
            @test all(value -> !occursin(value, text), negative_output_values)
        end

        REPLAY_SHA256[] = _bundle_digest(first_out)
        @test REPLAY_SHA256[] == _bundle_digest(second_out)

        interrupted = joinpath(fixture.root, "interrupted")
        @test _code(() -> publish_bundle(frozen, interrupted; failpoint=:after_first_stage)) ==
              :interrupted_publication
        @test !ispath(interrupted)
        @test isempty(filter(name -> startswith(name, ".structured-universe-stage-"),
                             readdir(fixture.root)))
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "patch-key independence and all-or-nothing receipt" begin
        fixture = _fixture()
        frozen = _freeze(fixture)
        forward = _write(joinpath(fixture.root, "patches", "forward.tsv"),
                         _patch_text(frozen.keys))
        backward = _write(joinpath(fixture.root, "patches", "backward.tsv"),
                          _patch_text(frozen.keys; value="-0.5"))
        @test validate_patch_keys(frozen, forward).row_count == length(frozen.keys)
        @test validate_patch_keys(frozen, backward).keys_sha256 == frozen.keys_sha256

        variants = [
            _patch_text(frozen.keys[1:end-1]),
            _patch_text(vcat(collect(frozen.keys), [LobeKey("20240303_x.sxm", 1)])),
        ]
        for text in variants
            path = _write(joinpath(fixture.root, "patches", "variant-$(_sha(codeunits(text))).tsv"), text)
            @test _code(() -> validate_patch_keys(frozen, path)) == :patch_key_mismatch
        end

        duplicate = collect(frozen.keys)
        push!(duplicate, first(frozen.keys))
        duplicate_path = _write(joinpath(fixture.root, "patches", "duplicate.tsv"),
                                _patch_text(duplicate))
        @test _code(() -> validate_patch_keys(frozen, duplicate_path)) == :duplicate_patch_key

        reordered_path = _write(joinpath(fixture.root, "patches", "reordered.tsv"),
                                _patch_text(reverse(collect(frozen.keys))))
        @test _code(() -> validate_patch_keys(frozen, reordered_path)) == :patch_key_order

        substituted = collect(frozen.keys)
        substituted[1] = LobeKey("20240303_x.sxm", 1)
        substituted_path = _write(joinpath(fixture.root, "patches", "substituted.tsv"),
                                  _patch_text(substituted))
        @test _code(() -> validate_patch_keys(frozen, substituted_path)) == :patch_key_mismatch

        malformed = _write(joinpath(fixture.root, "patches", "malformed.tsv"),
                           replace(_patch_text(frozen.keys), "res_p001" => "lobe"))
        @test _code(() -> validate_patch_keys(frozen, malformed)) == :duplicate_column
        nonfinite = _write(joinpath(fixture.root, "patches", "nonfinite.tsv"),
                           _patch_text(frozen.keys; value="Inf"))
        @test _code(() -> validate_patch_keys(frozen, nonfinite)) == :invalid_patch_value

        receipt = joinpath(fixture.root, "patches", "pair-receipt.toml")
        pair_hash = publish_patch_pair_receipt(frozen, forward, backward, receipt)
        prior = read(receipt)
        @test pair_hash == _sha(prior)
        @test publish_patch_pair_receipt(frozen, forward, backward, receipt) == pair_hash
        write(backward, _patch_text(frozen.keys[1:end-1]))
        @test _code(() -> publish_patch_pair_receipt(frozen, forward, backward, receipt)) ==
              :patch_key_mismatch
        @test read(receipt) == prior
        write(backward, _patch_text(frozen.keys; value="-0.5"))

        interrupted_receipt = joinpath(fixture.root, "patches", "interrupted.toml")
        @test _code(() -> publish_patch_pair_receipt(
            frozen, forward, backward, interrupted_receipt; failpoint=:before_install)) ==
              :interrupted_publication
        @test !ispath(interrupted_receipt)
        @test isempty(filter(name -> occursin(".stage-", name),
                             readdir(dirname(interrupted_receipt))))

        for (name, path, changed_text) in (
            ("forward", forward, _patch_text(frozen.keys; value="0.75")),
            ("backward", backward, _patch_text(frozen.keys; value="-0.75")),
        )
            original = read(path)
            changed_receipt = joinpath(fixture.root, "patches", "$name-changed.toml")
            code = _code(() -> publish_patch_pair_receipt(
                frozen,
                forward,
                backward,
                changed_receipt;
                _before_patch_reconciliation=() -> write(path, changed_text),
            ))
            @test code == :patch_content_changed
            @test !ispath(changed_receipt)
            @test isempty(filter(item -> occursin(".stage-", item),
                                 readdir(dirname(changed_receipt))))
            write(path, original)
            @test validate_patch_keys(frozen, path).sha256 == _sha(original)
        end

        authoritative = read(receipt)
        backward_original = read(backward)
        existing_code = _code(() -> publish_patch_pair_receipt(
            frozen,
            forward,
            backward,
            receipt;
            _before_patch_reconciliation=() -> write(
                backward,
                _patch_text(frozen.keys; value="-0.875"),
            ),
        ))
        @test existing_code == :patch_content_changed
        @test read(receipt) == authoritative
        @test isempty(filter(item -> occursin(".stage-", item),
                             readdir(dirname(receipt))))
        write(backward, backward_original)
        @test publish_patch_pair_receipt(frozen, forward, backward, receipt) == pair_hash

        symlink_path = joinpath(fixture.root, "patches", "linked.tsv")
        symlink(forward, symlink_path)
        @test _code(() -> validate_patch_keys(frozen, symlink_path)) == :symlink_rejected
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "fixed patch transforms" begin
        coordinates = collect(GRID_COORDINATES)
        patch = Matrix{Float64}(undef, 17, 17)
        for (u_index, u) in enumerate(coordinates), (t_index, t) in enumerate(coordinates)
            patch[u_index, t_index] = t + 10u
        end

        shifted_t = shift_patch(patch, 0.02, 0.0)
        shifted_u = shift_patch(patch, 0.0, -0.02)
        @test shifted_t[9, 9] ≈ -0.02 atol=1e-15
        @test shifted_t[9, 1] ≈ -0.32 atol=1e-15
        @test shifted_u[9, 9] ≈ 0.2 atol=1e-15
        @test shifted_u[17, 9] ≈ 3.2 atol=1e-15
        @test shift_patch(patch, -0.02, 0.0)[9, 17] ≈ 0.32 atol=1e-15
        @test shift_patch(patch, 0.0, 0.02)[1, 9] ≈ -3.2 atol=1e-15
        @test _code(() -> shift_patch(patch, 0.04, 0.0)) == :invalid_shift

        marked = Matrix{Union{Missing,Float64}}(undef, 2, 3)
        marked[1, :] = Union{Missing,Float64}[2.0, missing, Inf]
        marked[2, :] = Union{Missing,Float64}[-4.0, NaN, -Inf]
        contrasted = contrast_patch(marked, 0.95)
        @test contrasted[1, 1] === 1.9
        @test ismissing(contrasted[1, 2])
        @test contrasted[1, 3] == Inf
        @test contrasted[2, 1] === -3.8
        @test isnan(contrasted[2, 2])
        @test contrasted[2, 3] == -Inf
        @test _code(() -> contrast_patch(marked, 1.0)) == :invalid_contrast

        dropped = drop_channel(marked)
        @test size(dropped) == size(marked)
        @test all(ismissing, dropped)
        @test !any(value -> value === 0.0, dropped)
    end

    @testset "training-only frozen noise scale" begin
        fixture = _fixture()
        frozen = _freeze(fixture)
        base = [Float64(row + 3column) for row in 1:17, column in 1:17]
        patches = Dict{LobeKey,Matrix{Float64}}()
        for key in frozen.keys
            multiplier = startswith(key.file, "20240101") ? 1.0 : 20.0
            patches[key] = multiplier .* base .+ key.lobe
        end
        scale = freeze_noise_scale(frozen, "20240229", "forward", patches)
        changed_heldout = deepcopy(patches)
        for key in frozen.keys
            startswith(key.file, "20240229") || continue
            changed_heldout[key] .= 1.0e12 .* base
        end
        replayed = freeze_noise_scale(frozen, "20240229", "forward", changed_heldout)
        @test scale.target_mad == replayed.target_mad
        @test scale.binding_sha256 == replayed.binding_sha256
        @test scale.training_files == ("20240101_a.sxm", "20240101_b.sxm")

        changed_training = deepcopy(patches)
        changed_training[LobeKey("20240101_a.sxm", 1)] .*= 3.0
        altered = freeze_noise_scale(frozen, "20240229", "forward", changed_training)
        @test altered.target_mad != scale.target_mad
        @test altered.binding_sha256 != scale.binding_sha256

        unit_noise = correlated_noise(77; target_mad=1.0)
        @test size(unit_noise) == (17, 17)
        @test median(vec(unit_noise)) == 0.0
        @test median_absolute_deviation(vec(unit_noise)) ≈ 1.0 atol=2e-15
        @test correlated_noise(77; target_mad=1.0) == unit_noise
        @test correlated_noise(78; target_mad=1.0) != unit_noise

        zero_patch = zeros(17, 17)
        noisy = apply_noise_patch(zero_patch, scale, 19)
        @test median_absolute_deviation(vec(noisy)) ≈ scale.target_mad rtol=1e-14
        marked = Matrix{Union{Missing,Float64}}(zero_patch)
        marked[1, 1] = missing
        marked[1, 2] = NaN
        marked[1, 3] = Inf
        noisy_marked = apply_noise_patch(marked, scale, 19)
        @test ismissing(noisy_marked[1, 1])
        @test isnan(noisy_marked[1, 2])
        @test noisy_marked[1, 3] == Inf
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "path confinement and command surface" begin
        fixture = _fixture()
        @test _code(() -> freeze_universe(
            fixture.root;
            features="../features.tsv",
            feature_sha256=fixture.feature_sha256,
            candidate_config=fixture.candidate,
            model_config=fixture.model,
            source_paths=[UNIVERSE_SOURCE, UNIVERSE_CLI, FIREWALL_SOURCE],
        )) == :path_escape

        link = joinpath(fixture.root, "inputs", "features-link.tsv")
        symlink(fixture.features, link)
        linked = merge(fixture, (features=link, feature_sha256=_sha_file(fixture.features)))
        @test _code(() -> _freeze(linked)) == :symlink_rejected

        frozen = _freeze(fixture)
        real_parent = joinpath(fixture.root, "real-parent")
        mkpath(real_parent)
        linked_parent = joinpath(fixture.root, "linked-parent")
        symlink(real_parent, linked_parent)
        @test _code(() -> publish_bundle(frozen, joinpath(linked_parent, "out"))) ==
              :symlink_rejected

        stdout_one = joinpath(fixture.root, "cli-one.stdout")
        stderr_one = joinpath(fixture.root, "cli-one.stderr")
        arguments_one = [
            "--root", fixture.root,
            "--features", relpath(fixture.features, fixture.root),
            "--feature-sha256", fixture.feature_sha256,
            "--candidate-config", fixture.candidate,
            "--model-config", fixture.model,
            "--out-dir", "cli-one",
        ]
        @test _run_cli(arguments_one, stdout_one, stderr_one) == 0
        @test isdir(joinpath(fixture.root, "cli-one"))

        stdout_two = joinpath(fixture.root, "cli-two.stdout")
        stderr_two = joinpath(fixture.root, "cli-two.stderr")
        arguments_two = copy(arguments_one)
        arguments_two[end] = "cli-two"
        @test _run_cli(arguments_two, stdout_two, stderr_two) == 0
        @test _bundle_digest(joinpath(fixture.root, "cli-one")) ==
              _bundle_digest(joinpath(fixture.root, "cli-two"))

        bad_stdout = joinpath(fixture.root, "cli-bad.stdout")
        bad_stderr = joinpath(fixture.root, "cli-bad.stderr")
        bad_arguments = copy(arguments_one)
        bad_arguments[6] = repeat("0", 64)
        bad_arguments[end] = "cli-bad"
        @test _run_cli(bad_arguments, bad_stdout, bad_stderr) != 0
        @test !ispath(joinpath(fixture.root, "cli-bad"))
        rm(fixture.root; recursive=true, force=true)
    end

    @testset "new product text boundary" begin
        explicit_negative_terms = [
            "benchmark", "grader", "truth", "expected_count", "target_count",
            "composition", "top_k", "class_count", "label",
        ]
        for path in (UNIVERSE_SOURCE, UNIVERSE_CLI)
            source = lowercase(read(path, String))
            @test all(term -> !occursin(term, source), explicit_negative_terms)
            @test !occursin(r"\b145\b|\b870\b", source)
        end
    end
end

println("byte_replay_sha256=", REPLAY_SHA256[])
