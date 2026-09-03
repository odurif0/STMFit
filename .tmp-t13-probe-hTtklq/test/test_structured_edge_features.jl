#!/usr/bin/env julia

using Printf
using SHA
using Test
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const EDGE_SOURCE = joinpath(@__DIR__, "lib", "structured_assignment", "edge_features.jl")
const EDGE_CLI = joinpath(@__DIR__, "build_label_free_edge_features.jl")

if !isfile(EDGE_SOURCE) || !isfile(EDGE_CLI)
    @testset "structured edge features RED: frozen products are required" begin
        @test isfile(EDGE_SOURCE)
        @test isfile(EDGE_CLI)
    end
    exit(1)
end

include(EDGE_SOURCE)
using .StructuredEdgeFeatures
const SEF = StructuredEdgeFeatures

module EdgeCLIContract
include(Main.EDGE_CLI)
end
const EDGE_CLI_CONTRACT = EdgeCLIContract

const REPLAY_SHA256 = Ref("")

_sha(bytes) = bytes2hex(sha256(bytes))
_sha_file(path::String) = _sha(read(path))
_quote(value::String) = repr(value)

function _write(path::String, text::String)
    mkpath(dirname(path))
    write(path, text)
    return path
end

function _copy_dependency(root::String, relative::String)
    destination = joinpath(root, split(relative, '/')...)
    mkpath(dirname(destination))
    cp(joinpath(REPOSITORY_ROOT, split(relative, '/')...), destination; force=true)
    return destination
end

function _default_specs()
    return [
        (file="20260101_alpha.sxm", n=4, lobe=1, t=0.0, u=0.0, amplitude=1.0),
        (file="20260101_alpha.sxm", n=4, lobe=2, t=0.5, u=0.0, amplitude=1.1),
        (file="20260101_alpha.sxm", n=4, lobe=3, t=1.0, u=0.0, amplitude=1.2),
        (file="20260101_alpha.sxm", n=4, lobe=4, t=1.5, u=0.0, amplitude=1.3),
        (file="20260102_single.sxm", n=1, lobe=1, t=0.2, u=0.0, amplitude=0.8),
    ]
end

function _feature_header()
    return [
        "file", "N", "lobe", "amplitude", "x_nm", "y_nm", "t_nm", "u_nm",
        "sigma_parallel_nm", "sigma_perp_nm", "spacing_prev_nm", "amp_rel",
        "skew_ratio", "axis_x", "axis_y", "origin_x_nm", "origin_y_nm",
        "baseline", "tilt_x", "tilt_y", "gcv", "source",
    ]
end

function _feature_rows(specs)
    rows = Vector{Vector{String}}()
    by_file = Dict{String,Vector{Any}}()
    for spec in specs
        push!(get!(by_file, spec.file, Any[]), spec)
    end
    for file in sort!(collect(keys(by_file)))
        prior_t = nothing
        for spec in sort!(by_file[file]; by=item -> item.lobe)
            spacing = prior_t === nothing ? "NA" : @sprintf("%.6f", spec.t - prior_t)
            push!(rows, [
                spec.file,
                string(spec.n),
                string(spec.lobe),
                @sprintf("%.8e", spec.amplitude),
                @sprintf("%.6f", spec.t),
                @sprintf("%.6f", spec.u),
                @sprintf("%.6f", spec.t),
                @sprintf("%.6f", spec.u),
                "0.250000",
                "0.180000",
                spacing,
                @sprintf("%.6f", spec.amplitude / 1.3),
                "1.000000",
                "1.00000000",
                "0.00000000",
                "0.000000",
                "0.000000",
                "0.00000000e+00",
                "0.00000000e+00",
                "0.00000000e+00",
                "1.00000000e+00",
                "ell",
            ])
            prior_t = spec.t
        end
    end
    return rows
end

function _feature_text(specs)
    rows = _feature_rows(specs)
    return join(_feature_header(), '\t') * "\n" *
           join((join(row, '\t') for row in rows), "\n") * "\n"
end

function _sparse_vector(entries::Pair{Int,Float64}...)
    values = zeros(Float64, 289)
    for (index, value) in entries
        values[index] = value
    end
    return values
end

function _default_vectors(specs)
    forward_basis = [
        _sparse_vector(1 => 1.0, 2 => -1.0),
        7.0 .+ 3.0 .* _sparse_vector(1 => 1.0, 2 => -1.0),
        _sparse_vector(1 => 1.0, 2 => 1.0, 3 => -1.0, 4 => -1.0),
        _sparse_vector(1 => -1.0, 2 => -1.0, 3 => 1.0, 4 => 1.0),
    ]
    backward_basis = [
        _sparse_vector(7 => 1.0, 8 => -1.0),
        _sparse_vector(7 => -1.0, 8 => 1.0),
        _sparse_vector(7 => 1.0, 8 => -1.0),
        _sparse_vector(7 => 1.0, 8 => -1.0),
    ]
    forward = Dict{Tuple{String,Int},Vector{Float64}}()
    backward = Dict{Tuple{String,Int},Vector{Float64}}()
    for spec in specs
        slot = mod1(spec.lobe, length(forward_basis))
        forward[(spec.file, spec.lobe)] = copy(forward_basis[slot])
        backward[(spec.file, spec.lobe)] = copy(backward_basis[slot])
    end
    return forward, backward
end

_pixel_text(value::Float64) = isfinite(value) ? @sprintf("%.7g", value) : "NA"

function _metadata_by_key(specs)
    result = Dict{Tuple{String,Int},NTuple{3,String}}()
    for spec in specs
        result[(spec.file, spec.lobe)] = (
            @sprintf("%.6f", spec.t),
            @sprintf("%.6f", spec.u),
            @sprintf("%.8e", spec.amplitude),
        )
    end
    return result
end

function _patch_text(channel::Symbol, specs, vectors)
    header = channel == :forward ? SEF.FORWARD_PATCH_HEADER : SEF.BACKWARD_PATCH_HEADER
    metadata = _metadata_by_key(specs)
    rows = String[]
    for key in sort!(collect(keys(metadata)))
        t_text, u_text, amplitude_text = metadata[key]
        residual = _pixel_text.(vectors[key])
        payload = if channel == :forward
            vcat(fill("0", 289), residual)
        else
            vcat(fill("0", 289), residual, fill("0", 289), fill("0", 289))
        end
        push!(rows, join(vcat([key[1], string(key[2]), t_text, u_text, amplitude_text], payload), '\t'))
    end
    return join(header, '\t') * "\n" * join(rows, "\n") * "\n"
end

function _producer_receipt_text(
    fixture,
    channel::Symbol,
    frozen,
    command::String,
)
    output = channel == :forward ? fixture.forward : fixture.backward
    source = channel == :forward ?
        joinpath(fixture.root, "test", "extract_lobe_patches.jl") :
        joinpath(fixture.root, "test", "extract_lobe_patches_bwd.jl")
    header = channel == :forward ? SEF.FORWARD_PATCH_HEADER : SEF.BACKWARD_PATCH_HEADER
    fields = [
        "schema = \"stmfit-structured-patch-producer-receipt-v1\"",
        "schema_version = 1",
        "status = \"PASS\"",
        "channel = $(_quote(String(channel)))",
        "execution_site = \"Viper_Slurm_compute_node\"",
        "julia_threads = 1",
        "root = $(_quote(realpath(fixture.root)))",
        "cwd = $(_quote(realpath(fixture.root)))",
        "features_path = $(_quote(realpath(fixture.features)))",
        "data_path = $(_quote(realpath(fixture.data)))",
        "output_path = $(_quote(realpath(output)))",
        "source_path = $(_quote(realpath(source)))",
        "source_sha256 = $(_quote(_sha_file(source)))",
        "feature_sha256 = $(_quote(fixture.feature_sha256))",
        "output_sha256 = $(_quote(_sha_file(output)))",
        "candidate_config_sha256 = $(_quote(SEF.CANDIDATE_CONFIG_SHA256))",
        "model_config_sha256 = $(_quote(SEF.MODEL_CONFIG_SHA256))",
        "grid_sha256 = $(_quote(SEF.GRID_SHA256))",
        "keys_sha256 = $(_quote(frozen.keys_sha256))",
        "table_header_sha256 = $(_quote(_sha(codeunits(join(header, '\t') * "\n"))))",
        "row_count = $(length(frozen.keys))",
        "command = $(_quote(command))",
    ]
    return join(fields, "\n") * "\n"
end

function _make_fixture(; specs=_default_specs(), vector_transform=nothing)
    root = mktempdir(; prefix="stmfit-edge-features-")
    for relative in (
        "config/unit_assignment_structured_candidate.toml",
        "config/unit_assignment_structured_model.toml",
        "config/chitosan.toml",
        "test/extract_lobe_patches.jl",
        "test/extract_lobe_patches_bwd.jl",
    )
        _copy_dependency(root, relative)
    end
    data = joinpath(root, "synthetic-data")
    mkpath(data)
    features = _write(joinpath(root, "inputs", "features.tsv"), _feature_text(specs))
    feature_sha256 = _sha_file(features)
    forward_vectors, backward_vectors = _default_vectors(specs)
    if vector_transform !== nothing
        forward_vectors, backward_vectors = vector_transform(forward_vectors, backward_vectors)
    end
    forward = _write(joinpath(root, "patches", "forward.tsv"),
                     _patch_text(:forward, specs, forward_vectors))
    backward = _write(joinpath(root, "patches", "backward.tsv"),
                      _patch_text(:backward, specs, backward_vectors))
    fixture = (
        root=root,
        data=data,
        features=features,
        feature_sha256=feature_sha256,
        forward=forward,
        backward=backward,
        forward_receipt=joinpath(root, "patches", "forward-receipt.toml"),
        backward_receipt=joinpath(root, "patches", "backward-receipt.toml"),
        universe=joinpath(root, "universe"),
        specs=specs,
        forward_vectors=forward_vectors,
        backward_vectors=backward_vectors,
    )
    frozen = SEF.StructuredUniverse.freeze_universe(
        root;
        features=relpath(features, root),
        feature_sha256=feature_sha256,
        candidate_config="config/unit_assignment_structured_candidate.toml",
        model_config="config/unit_assignment_structured_model.toml",
        source_paths=[
            SEF.UNIVERSE_SOURCE_PATH,
            SEF.UNIVERSE_CLI_PATH,
            SEF.FIREWALL_SOURCE_PATH,
        ],
    )
    SEF.StructuredUniverse.publish_bundle(frozen, fixture.universe)
    commands = SEF.expected_producer_commands(
        root;
        features=features,
        data_dir=data,
        forward_patches=forward,
        backward_patches=backward,
    )
    _write(fixture.forward_receipt,
           _producer_receipt_text(fixture, :forward, frozen, commands.forward))
    _write(fixture.backward_receipt,
           _producer_receipt_text(fixture, :backward, frozen, commands.backward))
    return merge(fixture, (frozen=frozen, commands=commands))
end

function _input_keywords(fixture)
    return (
        features=relpath(fixture.features, fixture.root),
        feature_sha256=fixture.feature_sha256,
        candidate_config="config/unit_assignment_structured_candidate.toml",
        model_config="config/unit_assignment_structured_model.toml",
        universe_dir=relpath(fixture.universe, fixture.root),
        forward_patches=relpath(fixture.forward, fixture.root),
        backward_patches=relpath(fixture.backward, fixture.root),
        forward_receipt=relpath(fixture.forward_receipt, fixture.root),
        backward_receipt=relpath(fixture.backward_receipt, fixture.root),
    )
end

function _build(fixture)
    return SEF.build_edge_bundle(fixture.root; _input_keywords(fixture)...)
end

function _publish(fixture, output::String)
    return SEF.publish_edge_bundle(fixture.root, output; _input_keywords(fixture)...)
end

function _code(action::Function)
    try
        action()
    catch error
        if error isa SEF.EdgeFeatureError || error isa SEF.StructuredUniverse.UniverseError
            return error.code
        end
        rethrow()
    end
    return :accepted
end

function _race_publish(
    install_competitor::Function,
    fixture,
    output::String,
)
    hook_calls = Ref(0)
    previous_hook = SEF._PUBLICATION_INSTALL_HOOK[]
    SEF._PUBLICATION_INSTALL_HOOK[] = function(stage, destination)
        @test startswith(basename(stage), ".structured-edge-stage-")
        @test destination == output
        hook_calls[] += 1
        install_competitor(destination)
        return nothing
    end
    code = try
        _code(() -> _publish(fixture, output))
    finally
        SEF._PUBLICATION_INSTALL_HOOK[] = previous_hook
    end
    return code, hook_calls[]
end

function _parse_tsv_bytes(bytes::Vector{UInt8})
    lines = split(chomp(String(copy(bytes))), '\n')
    header = split(first(lines), '\t'; keepempty=true)
    rows = [Dict(header[index] => value for (index, value) in enumerate(
                 split(line, '\t'; keepempty=true))) for line in lines[2:end]]
    return header, rows
end

function _bundle_digest(bundle)
    io = IOBuffer()
    for name in sort!(collect(keys(bundle.files)))
        println(io, name, '\t', _sha(bundle.files[name]), '\t', length(bundle.files[name]))
    end
    return _sha(take!(io))
end

function _table_replay_digest(bundle)
    material = "edge-observations=" *
               _sha(bundle.files["edge_observations.tsv"]) * "\n" *
               "node-segments=" *
               _sha(bundle.files["node_segments.tsv"]) * "\n"
    return _sha(codeunits(material))
end

function _centered_cosine(left::Vector{Float64}, right::Vector{Float64})
    left_centered = left .- sum(left) / 289
    right_centered = right .- sum(right) / 289
    return sum(left_centered .* right_centered) /
           sqrt(sum(abs2, left_centered) * sum(abs2, right_centered))
end

function _refresh_receipt!(fixture, channel::Symbol)
    command = channel == :forward ? fixture.commands.forward : fixture.commands.backward
    path = channel == :forward ? fixture.forward_receipt : fixture.backward_receipt
    write(path, _producer_receipt_text(fixture, channel, fixture.frozen, command))
    return nothing
end

function _mutate_header!(path::String, transform::Function)
    lines = split(chomp(read(path, String)), '\n')
    header = split(first(lines), '\t'; keepempty=true)
    lines[1] = join(transform(header), '\t')
    write(path, join(lines, "\n") * "\n")
    return nothing
end

function _mutate_rows!(path::String, transform::Function)
    lines = split(chomp(read(path, String)), '\n')
    header = first(lines)
    rows = transform(copy(lines[2:end]))
    write(path, header * "\n" * join(rows, "\n") * "\n")
    return nothing
end

function _with_patch_mutation(
    fixture,
    channel::Symbol,
    wanted::Symbol,
    mutation::Function,
)
    path = channel == :forward ? fixture.forward : fixture.backward
    receipt = channel == :forward ? fixture.forward_receipt : fixture.backward_receipt
    original_path = read(path)
    original_receipt = read(receipt)
    try
        mutation(path)
        _refresh_receipt!(fixture, channel)
        @test _code(() -> _build(fixture)) == wanted
    finally
        write(path, original_path)
        write(receipt, original_receipt)
    end
    return nothing
end

function _run_cli(arguments)
    stdout = IOBuffer()
    stderr = IOBuffer()
    command = `$(Base.julia_cmd()) --project=$(REPOSITORY_ROOT) $(EDGE_CLI) $arguments`
    process = run(pipeline(command; stdout=stdout, stderr=stderr); wait=false)
    wait(process)
    return process.exitcode, String(take!(stdout)), String(take!(stderr))
end

function _cli_parse_code(arguments::Vector{String})
    try
        EDGE_CLI_CONTRACT.parse_arguments(arguments)
    catch error
        return Symbol(nameof(typeof(error)))
    end
    return :accepted
end

@testset "provenance-bound observed-edge features" begin
    @testset "frozen grid bytes, command capture, and public surface" begin
        expected_grid = "version=1\n" *
                        "half_t_nm=0.32\n" *
                        "half_u_nm=0.32\n" *
                        "step_nm=0.04\n" *
                        "order=u_outer,t_inner\n" *
                        "coordinates_nm=-0.32,-0.28,-0.24,-0.20,-0.16,-0.12,-0.08,-0.04,0.00,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32\n" *
                        "forward_prefix=res_p\n" *
                        "backward_prefix=bwd_res_p\n" *
                        "pixels=289\n"
        @test String(SEF.grid_contract_bytes()) == expected_grid
        @test _sha(SEF.grid_contract_bytes()) == SEF.GRID_SHA256
        @test SEF.GRID_SHA256 == "d281d836d4bc5a1657762a46bc2ee0ff51ff9b3a4aff6068dee55632797d950e"
        @test SEF.CANDIDATE_CONFIG_SHA256 == "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
        @test SEF.MODEL_CONFIG_SHA256 == "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
        @test _sha_file(joinpath(REPOSITORY_ROOT, "config", "unit_assignment_structured_candidate.toml")) ==
              SEF.CANDIDATE_CONFIG_SHA256
        @test _sha_file(joinpath(REPOSITORY_ROOT, "config", "unit_assignment_structured_model.toml")) ==
              SEF.MODEL_CONFIG_SHA256
        @test _sha_file(SEF.T2_REVIEW_PATH) == SEF.T2_REVIEW_SHA256
        @test _sha_file(SEF.T3_REVIEW_PATH) == SEF.T3_REVIEW_SHA256

        fixture = _make_fixture()
        try
            root = realpath(fixture.root)
            @test fixture.commands.forward ==
                  "julia --project=$root $root/test/extract_lobe_patches.jl --features $(realpath(fixture.features)) --data-dir $(realpath(fixture.data)) --half-nm 0.32 --half-u-nm 0.32 --step-nm 0.04 --out $(realpath(fixture.forward))"
            @test fixture.commands.backward ==
                  "julia --project=$root $root/test/extract_lobe_patches_bwd.jl --features $(realpath(fixture.features)) --data-dir $(realpath(fixture.data)) --half-nm 0.32 --step-nm 0.04 --out $(realpath(fixture.backward))"
            @test !occursin("--step-nm 0.08", fixture.commands.forward)
            @test !occursin("--step-nm 0.08", fixture.commands.backward)
        finally
            rm(fixture.root; recursive=true, force=true)
        end

        exported = Set(names(SEF; all=false, imported=false))
        forbidden_public = Set([
            :build_from_matrix, :edge_matrix, :alternate_observation,
            :pair_prior, :transition, :switch_penalty, :run_length,
        ])
        @test isempty(intersect(exported, forbidden_public))
        @test !occursin(r"\brun\s*\(", read(EDGE_SOURCE, String))
        @test !occursin(r"\brun\s*\(", read(EDGE_CLI, String))
    end

    @testset "centered cosine, positive affine invariance, and exact schemas" begin
        fixture = _make_fixture()
        try
            bundle = _build(fixture)
            repeat_bundle = _build(fixture)
            @test bundle.files == repeat_bundle.files
            @test _bundle_digest(bundle) == _bundle_digest(repeat_bundle)
            REPLAY_SHA256[] = _table_replay_digest(bundle)

            edge_header, edge_rows = _parse_tsv_bytes(bundle.files["edge_observations.tsv"])
            node_header, node_rows = _parse_tsv_bytes(bundle.files["node_segments.tsv"])
            @test edge_header == SEF.EDGE_HEADER
            @test node_header == SEF.NODE_HEADER
            @test length(edge_rows) == 3
            @test length(node_rows) == 5
            @test all(row -> row["edge_status"] == "eligible", edge_rows)
            @test all(row -> row["split_reason"] == "none", edge_rows)
            @test all(row -> Set([row["corr_fwd"], row["corr_bwd"]]) != Set(["NA"]), edge_rows)
            @test Set(edge_header) == Set([
                "file", "left_lobe", "right_lobe", "left_t_nm", "right_t_nm",
                "gap_nm", "left_segment_id", "right_segment_id", "edge_status",
                "split_reason", "corr_fwd", "corr_bwd", "grid_hash",
                "forward_patch_hash", "backward_patch_hash", "feature_hash",
                "config_hash", "source_hash",
            ])
            @test isempty(filter(name -> occursin("raw", name) || occursin("diff", name), edge_header))

            first_row = first(edge_rows)
            key_left = (first_row["file"], parse(Int, first_row["left_lobe"]))
            key_right = (first_row["file"], parse(Int, first_row["right_lobe"]))
            expected_forward = _centered_cosine(
                fixture.forward_vectors[key_left], fixture.forward_vectors[key_right])
            expected_backward = _centered_cosine(
                fixture.backward_vectors[key_left], fixture.backward_vectors[key_right])
            @test parse(Float64, first_row["corr_fwd"]) ≈ expected_forward atol=1e-15
            @test parse(Float64, first_row["corr_bwd"]) ≈ expected_backward atol=1e-15
            @test first_row["corr_fwd"] ==
                  @sprintf("%.17g", parse(Float64, first_row["corr_fwd"]))
            @test first_row["corr_bwd"] ==
                  @sprintf("%.17g", parse(Float64, first_row["corr_bwd"]))
            @test expected_forward ≈ 1.0 atol=1e-15

            receipt = TOML.parse(String(bundle.files["receipt.toml"]))
            @test receipt["schema"] == "stmfit-structured-edge-feature-receipt-v1"
            @test receipt["edge_observations_sha256"] == _sha(bundle.files["edge_observations.tsv"])
            @test receipt["node_segments_sha256"] == _sha(bundle.files["node_segments.tsv"])
            @test receipt["forward_patch_sha256"] == _sha_file(fixture.forward)
            @test receipt["backward_patch_sha256"] == _sha_file(fixture.backward)
            @test receipt["feature_sha256"] == fixture.feature_sha256
            @test receipt["candidate_config_sha256"] == SEF.CANDIDATE_CONFIG_SHA256
            @test receipt["model_config_sha256"] == SEF.MODEL_CONFIG_SHA256
            @test receipt["grid_sha256"] == SEF.GRID_SHA256
            @test receipt["forward_command"] == fixture.commands.forward
            @test receipt["backward_command"] == fixture.commands.backward

            alpha_segments = unique(row["segment_id"] for row in node_rows
                                    if row["file"] == "20260101_alpha.sxm")
            @test length(alpha_segments) == 1
            expected_segment = _sha(codeunits(
                "segment-v1\n" *
                "file=20260101_alpha.sxm\n" *
                "first_lobe=1\n" *
                "last_lobe=4\n" *
                "grid_hash=$(SEF.GRID_SHA256)\n",
            ))[1:16]
            @test only(alpha_segments) == expected_segment
            singleton = only(filter(row -> row["file"] == "20260102_single.sxm", node_rows))
            @test singleton["node_status"] == "isolated"
            @test singleton["left_boundary_reason"] == "start"
            @test singleton["right_boundary_reason"] == "end"
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end

    @testset "exact one-based reversal and topology equivariance" begin
        reversal_specs = [
            (file="20260101_alpha.sxm", n=4, lobe=1, t=0.0, u=0.0, amplitude=1.0),
            (file="20260101_alpha.sxm", n=4, lobe=2, t=0.5, u=0.0, amplitude=1.1),
            (file="20260101_alpha.sxm", n=4, lobe=3, t=1.4, u=0.0, amplitude=1.2),
            (file="20260101_alpha.sxm", n=4, lobe=4, t=1.9, u=0.0, amplitude=1.3),
            (file="20260102_single.sxm", n=1, lobe=1, t=0.2, u=0.0, amplitude=0.8),
        ]
        original = _make_fixture(specs=reversal_specs)
        reverse_specs = [merge(spec, (t=-spec.t,)) for spec in original.specs]
        reversed_fixture = _make_fixture(
            specs=reverse_specs,
            vector_transform=(forward, backward) -> (
                Dict(key => SEF._reverse_patch_values(value) for (key, value) in forward),
                Dict(key => SEF._reverse_patch_values(value) for (key, value) in backward),
            ),
        )
        try
            marker = collect(1.0:289.0)
            transformed = SEF._reverse_patch_values(marker)
            for u_index in 1:17, t_index in 1:17
                k = 17 * (u_index - 1) + t_index
                reversed_index = 17 * (u_index - 1) + (18 - t_index)
                @test transformed[k] == marker[reversed_index]
            end
            @test transformed != reverse(marker)
            @test transformed != marker

            original_bundle = _build(original)
            reversed_bundle = _build(reversed_fixture)
            _, original_edges = _parse_tsv_bytes(original_bundle.files["edge_observations.tsv"])
            _, reversed_edges = _parse_tsv_bytes(reversed_bundle.files["edge_observations.tsv"])
            canonical(rows) = Dict(
                (row["file"], min(parse(Int, row["left_lobe"]), parse(Int, row["right_lobe"])),
                 max(parse(Int, row["left_lobe"]), parse(Int, row["right_lobe"]))) => (
                    corr_fwd=row["corr_fwd"], corr_bwd=row["corr_bwd"],
                    status=row["edge_status"], reason=row["split_reason"],
                    segments=Set([row["left_segment_id"], row["right_segment_id"]]),
                 ) for row in rows)
            @test canonical(original_edges) == canonical(reversed_edges)

            _, original_nodes = _parse_tsv_bytes(original_bundle.files["node_segments.tsv"])
            _, reversed_nodes = _parse_tsv_bytes(reversed_bundle.files["node_segments.tsv"])
            original_segments = Dict((row["file"], row["lobe"]) => row["segment_id"]
                                     for row in original_nodes)
            reversed_segments = Dict((row["file"], row["lobe"]) => row["segment_id"]
                                     for row in reversed_nodes)
            @test original_segments == reversed_segments

            left = original.forward_vectors[("20260101_alpha.sxm", 1)]
            right = original.forward_vectors[("20260101_alpha.sxm", 2)]
            @test (right, left) != (SEF._reverse_patch_values(right),
                                    SEF._reverse_patch_values(left))
        finally
            rm(original.root; recursive=true, force=true)
            rm(reversed_fixture.root; recursive=true, force=true)
        end
    end

    @testset "gap splits, per-channel failures, singleton preservation, and precedence" begin
        specs = [
            (file="20260101_split.sxm", n=4, lobe=1, t=0.0, u=0.0, amplitude=1.0),
            (file="20260101_split.sxm", n=4, lobe=2, t=0.5, u=0.0, amplitude=1.1),
            (file="20260101_split.sxm", n=4, lobe=3, t=1.4, u=0.0, amplitude=1.2),
            (file="20260101_split.sxm", n=4, lobe=4, t=1.9, u=0.0, amplitude=1.3),
            (file="20260102_single.sxm", n=1, lobe=1, t=0.0, u=0.0, amplitude=0.8),
        ]
        fixture = _make_fixture(specs=specs,
            vector_transform=(forward, backward) -> begin
                backward[("20260101_split.sxm", 4)] = zeros(289)
                return forward, backward
            end)
        try
            bundle = _build(fixture)
            _, edges = _parse_tsv_bytes(bundle.files["edge_observations.tsv"])
            @test [(row["edge_status"], row["split_reason"]) for row in edges] == [
                ("eligible", "none"),
                ("split", "gap_out_of_range"),
                ("split", "zero_variance_bwd"),
            ]
            @test edges[2]["corr_fwd"] == "NA"
            @test edges[2]["corr_bwd"] == "NA"
            @test edges[2]["left_segment_id"] != edges[2]["right_segment_id"]
            _, nodes = _parse_tsv_bytes(bundle.files["node_segments.tsv"])
            split_nodes = filter(row -> row["file"] == "20260101_split.sxm", nodes)
            @test [row["node_status"] for row in split_nodes] ==
                  ["connected", "connected", "isolated", "isolated"]
            @test split_nodes[2]["right_boundary_reason"] == "gap_out_of_range"
            @test split_nodes[3]["left_boundary_reason"] == "gap_out_of_range"
            @test split_nodes[3]["right_boundary_reason"] == "zero_variance_bwd"
            @test split_nodes[4]["left_boundary_reason"] == "zero_variance_bwd"

            @test SEF._select_split_reason(false, :missing, :missing, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :gap_out_of_range
            @test SEF._select_split_reason(true, :missing, :missing, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :missing_forward_key
            @test SEF._select_split_reason(true, :valid, :missing, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :missing_backward_key
            @test SEF._select_split_reason(true, :nonfinite, :valid, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :nonfinite_fwd
            @test SEF._select_split_reason(true, :valid, :nonfinite, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :nonfinite_bwd
            @test SEF._select_split_reason(true, :zero_variance, :valid, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :zero_variance_fwd
            @test SEF._select_split_reason(true, :valid, :zero_variance, :nonfinite,
                                          :nonfinite, :out_of_range, :out_of_range) ==
                  :zero_variance_bwd
            @test SEF._select_split_reason(true, :valid, :valid, :out_of_range,
                                          :valid, :out_of_range, :out_of_range) ==
                  :corr_out_of_range_fwd
            @test SEF._select_split_reason(true, :valid, :valid, :valid,
                                          :out_of_range, :out_of_range, :out_of_range) ==
                  :corr_out_of_range_bwd
            @test SEF._guard_dot(1.0 + 5.0e-13) == (status=:valid, value=1.0)
            @test SEF._guard_dot(-1.0 - 5.0e-13) == (status=:valid, value=-1.0)
            @test SEF._guard_dot(1.0 + 2.0e-12).status == :out_of_range
            @test SEF._guard_dot(NaN).status == :out_of_range
        finally
            rm(fixture.root; recursive=true, force=true)
        end

        inserted_specs = [
            (file="20260101_inserted.sxm", n=3, lobe=1, t=0.0, u=0.0, amplitude=1.0),
            (file="20260101_inserted.sxm", n=3, lobe=2, t=1.0, u=0.0, amplitude=1.1),
            (file="20260101_inserted.sxm", n=3, lobe=3, t=0.5, u=0.0, amplitude=1.2),
            (file="20260102_single.sxm", n=1, lobe=1, t=0.0, u=0.0, amplitude=0.8),
        ]
        inserted = _make_fixture(specs=inserted_specs)
        try
            bundle = _build(inserted)
            _, edges = _parse_tsv_bytes(bundle.files["edge_observations.tsv"])
            @test [(row["left_lobe"], row["right_lobe"], row["gap_nm"]) for row in edges] ==
                  [("1", "3", "0.5"), ("3", "2", "0.5")]
            @test all(row -> row["edge_status"] == "eligible", edges)
            _, nodes = _parse_tsv_bytes(bundle.files["node_segments.tsv"])
            inserted_segment = only(unique(
                row["segment_id"] for row in nodes
                if row["file"] == "20260101_inserted.sxm"
            ))
            expected_inserted_segment = _sha(codeunits(
                "segment-v1\n" *
                "file=20260101_inserted.sxm\n" *
                "first_lobe=1\n" *
                "last_lobe=2\n" *
                "grid_hash=$(SEF.GRID_SHA256)\n",
            ))[1:16]
            @test inserted_segment == expected_inserted_segment
        finally
            rm(inserted.root; recursive=true, force=true)
        end
    end

    @testset "strict schema, keys, producer receipts, and nonfinite split paths" begin
        fixture = _make_fixture()
        try
            _with_patch_mutation(fixture, :forward, :missing_residual_column,
                path -> _mutate_header!(path, header -> filter(!=("res_p010"), header)))
            _with_patch_mutation(fixture, :forward, :extra_residual_column,
                path -> _mutate_header!(path, header -> vcat(header, ["res_p290"])))
            _with_patch_mutation(fixture, :forward, :duplicate_column,
                path -> _mutate_header!(path, header -> vcat(header[1:20], [header[20]], header[21:end])))
            _with_patch_mutation(fixture, :forward, :reordered_residual_columns,
                path -> _mutate_header!(path, header -> begin
                    changed = copy(header)
                    first_index = findfirst(==("res_p010"), changed)
                    second_index = findfirst(==("res_p011"), changed)
                    changed[first_index], changed[second_index] = changed[second_index], changed[first_index]
                    changed
                end))
            _with_patch_mutation(fixture, :backward, :patch_schema_mismatch,
                path -> _mutate_header!(path, header -> vcat(header, ["alternate_pair_observation"])))

            _with_patch_mutation(fixture, :forward, :patch_key_mismatch,
                path -> _mutate_rows!(path, rows -> rows[1:end-1]))
            _with_patch_mutation(fixture, :forward, :patch_key_mismatch,
                path -> _mutate_rows!(path, rows -> vcat(rows, [replace(first(rows),
                    "20260101_alpha.sxm" => "20260103_extra.sxm")])))
            _with_patch_mutation(fixture, :forward, :duplicate_patch_key,
                path -> _mutate_rows!(path, rows -> vcat(rows, [first(rows)])))
            _with_patch_mutation(fixture, :forward, :patch_key_order,
                path -> _mutate_rows!(path, reverse))

            original_forward = read(fixture.forward)
            original_receipt = read(fixture.forward_receipt)
            try
                lines = split(chomp(read(fixture.forward, String)), '\n')
                header = split(lines[1], '\t'; keepempty=true)
                residual_index = findfirst(==("res_p001"), header)
                fields = split(lines[2], '\t'; keepempty=true)
                fields[residual_index] = "NA"
                lines[2] = join(fields, '\t')
                write(fixture.forward, join(lines, "\n") * "\n")
                _refresh_receipt!(fixture, :forward)
                bundle = _build(fixture)
                _, edges = _parse_tsv_bytes(bundle.files["edge_observations.tsv"])
                @test first(edges)["split_reason"] == "nonfinite_fwd"
            finally
                write(fixture.forward, original_forward)
                write(fixture.forward_receipt, original_receipt)
            end

            original_forward = read(fixture.forward)
            original_receipt = read(fixture.forward_receipt)
            try
                lines = split(chomp(read(fixture.forward, String)), '\n')
                header = split(lines[1], '\t'; keepempty=true)
                residual_indices = [findfirst(==(@sprintf("res_p%03d", index)), header)
                                    for index in 1:289]
                fields = split(lines[2], '\t'; keepempty=true)
                fields[residual_indices] .= "2.5"
                lines[2] = join(fields, '\t')
                write(fixture.forward, join(lines, "\n") * "\n")
                _refresh_receipt!(fixture, :forward)
                bundle = _build(fixture)
                _, edges = _parse_tsv_bytes(bundle.files["edge_observations.tsv"])
                @test first(edges)["split_reason"] == "zero_variance_fwd"
            finally
                write(fixture.forward, original_forward)
                write(fixture.forward_receipt, original_receipt)
            end

            forward_receipt = read(fixture.forward_receipt, String)
            cases = [
                (:producer_command_mismatch,
                 replace(forward_receipt, "--step-nm 0.04" => "--step-nm 0.08")),
                (:producer_source_hash_mismatch,
                 replace(forward_receipt, SEF.FORWARD_PRODUCER_SHA256 => repeat("0", 64))),
                (:grid_hash_mismatch,
                 replace(forward_receipt, SEF.GRID_SHA256 => repeat("1", 64))),
                (:patch_hash_mismatch,
                 replace(forward_receipt, _sha_file(fixture.forward) => repeat("2", 64))),
                (:patch_receipt_schema,
                 forward_receipt * "unexpected = \"field\"\n"),
            ]
            for (wanted, text) in cases
                write(fixture.forward_receipt, text)
                @test _code(() -> _build(fixture)) == wanted
            end
            write(fixture.forward_receipt, forward_receipt)

            alternate_data = joinpath(fixture.root, "alternate-synthetic-data")
            mkpath(alternate_data)
            backward_receipt = read(fixture.backward_receipt, String)
            write(
                fixture.backward_receipt,
                replace(
                    backward_receipt,
                    realpath(fixture.data) => realpath(alternate_data),
                ),
            )
            @test _code(() -> _build(fixture)) == :producer_path_mismatch
            write(fixture.backward_receipt, backward_receipt)

            candidate = joinpath(
                fixture.root,
                "config",
                "unit_assignment_structured_candidate.toml",
            )
            candidate_bytes = read(candidate)
            write(candidate, vcat(candidate_bytes, UInt8('\n')))
            @test _code(() -> _build(fixture)) == :config_hash_mismatch
            write(candidate, candidate_bytes)

            model = joinpath(
                fixture.root,
                "config",
                "unit_assignment_structured_model.toml",
            )
            model_bytes = read(model)
            write(model, replace(String(copy(model_bytes)), "step_nm = 0.04" => "step_nm = 0.08"))
            @test _code(() -> _build(fixture)) == :config_hash_mismatch
            write(model, model_bytes)

            producer = joinpath(fixture.root, "test", "extract_lobe_patches.jl")
            producer_bytes = read(producer)
            write(producer, vcat(producer_bytes, UInt8('\n')))
            @test _code(() -> _build(fixture)) == :source_hash_mismatch
            write(producer, producer_bytes)
            @test _build(fixture) isa SEF.EdgeBundle
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end

    @testset "atomic publication, replay validation, and CLI firewall" begin
        fixture = _make_fixture()
        try
            output_one = joinpath(fixture.root, "edge-one")
            output_two = joinpath(fixture.root, "edge-two")
            first_hashes = _publish(fixture, output_one)
            second_hashes = _publish(fixture, output_two)
            @test first_hashes == second_hashes
            @test sort(readdir(output_one)) ==
                  ["edge_observations.tsv", "node_segments.tsv", "receipt.toml"]
            @test all(name -> read(joinpath(output_one, name)) == read(joinpath(output_two, name)),
                      readdir(output_one))
            @test SEF.validate_edge_bundle(fixture.root, output_one;
                _input_keywords(fixture)...) == first_hashes
            @test _publish(fixture, output_one) == first_hashes

            for attempt in 1:3
                interrupted = joinpath(fixture.root, "interrupted-$(attempt)")
                @test _code(() -> SEF._publish_edge_bundle(
                    fixture.root, interrupted; _input_keywords(fixture)...,
                    failpoint=:before_install)) == :interrupted_publication
                @test !ispath(interrupted)
                @test isempty(filter(name -> startswith(name, ".structured-edge-stage-"),
                                     readdir(fixture.root)))
            end

            raced_empty = joinpath(fixture.root, "raced-empty")
            competitor_inode = Ref{UInt64}(0)
            race_code, hook_calls = _race_publish(fixture, raced_empty) do destination
                mkdir(destination)
                competitor_inode[] = lstat(destination).inode
            end
            @test race_code == :publication_collision
            @test hook_calls == 1
            @test isdir(raced_empty)
            @test !islink(raced_empty)
            @test lstat(raced_empty).inode == competitor_inode[]
            @test isempty(readdir(raced_empty))
            @test isempty(filter(name -> startswith(name, ".structured-edge-stage-"),
                                 readdir(fixture.root)))

            raced_file = joinpath(fixture.root, "raced-file")
            file_bytes = collect(codeunits("competitor-file\n"))
            file_inode = Ref{UInt64}(0)
            file_code, file_hook_calls = _race_publish(fixture, raced_file) do destination
                write(destination, file_bytes)
                file_inode[] = lstat(destination).inode
            end
            @test file_code == :publication_collision
            @test file_hook_calls == 1
            @test isfile(raced_file)
            @test !islink(raced_file)
            @test lstat(raced_file).inode == file_inode[]
            @test read(raced_file) == file_bytes

            symlink_target = joinpath(fixture.root, "raced-symlink-target")
            symlink_target_bytes = collect(codeunits("competitor-target\n"))
            write(symlink_target, symlink_target_bytes)
            raced_symlink = joinpath(fixture.root, "raced-symlink")
            symlink_inode = Ref{UInt64}(0)
            symlink_code, symlink_hook_calls = _race_publish(
                fixture,
                raced_symlink,
            ) do destination
                symlink(symlink_target, destination)
                symlink_inode[] = lstat(destination).inode
            end
            @test symlink_code == :publication_collision
            @test symlink_hook_calls == 1
            @test islink(raced_symlink)
            @test lstat(raced_symlink).inode == symlink_inode[]
            @test readlink(raced_symlink) == symlink_target
            @test read(symlink_target) == symlink_target_bytes

            raced_nonempty = joinpath(fixture.root, "raced-nonempty")
            marker_name = "owner.marker"
            marker_bytes = collect(codeunits("competitor-directory\n"))
            directory_inode = Ref{UInt64}(0)
            directory_code, directory_hook_calls = _race_publish(
                fixture,
                raced_nonempty,
            ) do destination
                mkdir(destination)
                write(joinpath(destination, marker_name), marker_bytes)
                directory_inode[] = lstat(destination).inode
            end
            @test directory_code == :publication_collision
            @test directory_hook_calls == 1
            @test isdir(raced_nonempty)
            @test !islink(raced_nonempty)
            @test lstat(raced_nonempty).inode == directory_inode[]
            @test readdir(raced_nonempty) == [marker_name]
            @test read(joinpath(raced_nonempty, marker_name)) == marker_bytes
            @test isempty(filter(name -> startswith(name, ".structured-edge-stage-"),
                                 readdir(fixture.root)))

            edge_path = joinpath(output_one, "edge_observations.tsv")
            edge_bytes = read(edge_path)
            write(edge_path, vcat(edge_bytes, UInt8('\n')))
            @test _code(() -> SEF.validate_edge_bundle(
                fixture.root, output_one; _input_keywords(fixture)...)) == :bundle_mismatch
            write(edge_path, edge_bytes)

            cli_output = joinpath(fixture.root, "cli-output")
            arguments = [
                "--root", fixture.root,
                "--features", relpath(fixture.features, fixture.root),
                "--feature-sha256", fixture.feature_sha256,
                "--candidate-config", "config/unit_assignment_structured_candidate.toml",
                "--model-config", "config/unit_assignment_structured_model.toml",
                "--universe-dir", relpath(fixture.universe, fixture.root),
                "--forward-patches", relpath(fixture.forward, fixture.root),
                "--backward-patches", relpath(fixture.backward, fixture.root),
                "--forward-receipt", relpath(fixture.forward_receipt, fixture.root),
                "--backward-receipt", relpath(fixture.backward_receipt, fixture.root),
                "--out-dir", relpath(cli_output, fixture.root),
            ]
            code, stdout, stderr = _run_cli(arguments)
            @test code == 0
            @test isdir(cli_output)
            @test occursin("status=PASS", stdout)
            @test isempty(stderr)

            forbidden_controls = (
                ["--matrix", "matrix.tsv"],
                ["--alternate-edge-observation", "x"],
                ["--pair-prior", "0.25"],
                ["--transition", "0.1"],
                ["--switch-penalty", "1"],
                ["--run-length", "2"],
                ["--composition", "0.5"],
                ["--N", "6"],
            )
            for extra in forbidden_controls
                @test _cli_parse_code(extra) != :accepted
            end

            rejected_output = joinpath(fixture.root, "rejected")
            rejected_arguments = vcat(
                arguments[1:end-2],
                ["--out-dir", rejected_output],
                first(forbidden_controls),
            )
            rejected_code, _, rejected_stderr = _run_cli(rejected_arguments)
            @test rejected_code != 0
            @test !ispath(rejected_output)
            @test !isempty(rejected_stderr)

            help_code, help_stdout, help_stderr = _run_cli(["--help"])
            @test help_code == 0
            @test isempty(help_stderr)
            for token in ("matrix", "alternate", "pair-prior", "transition", "switch",
                          "run-length", "composition", "benchmark", "grader")
                @test !occursin(token, lowercase(help_stdout))
            end
        finally
            rm(fixture.root; recursive=true, force=true)
        end
    end
end

println("edge_tables_sha256=", REPLAY_SHA256[])
