#!/usr/bin/env julia

using Test

const ROOT = dirname(@__DIR__)
const FIREWALL_SOURCE = joinpath(@__DIR__, "lib", "structured_assignment", "firewall.jl")

if !isfile(FIREWALL_SOURCE)
    @testset "structured firewall RED: source is required" begin
        @test isfile(FIREWALL_SOURCE)
    end
    exit(1)
end

module StructuredFirewallUnderTest
include(Main.FIREWALL_SOURCE)
end

const SAF = StructuredFirewallUnderTest
const PHYSICAL = [
    "amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated",
    "amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
    "bwd_neg_com_t", "bwd_neg_diag45", "split_log_skew", "corr_fwd", "corr_bwd",
]

function nested_percent(text::String, additional_levels::Int)
    result = text
    for _ in 1:additional_levels
        result = replace(result, "%" => "%25")
    end
    return result
end

percent_encode(text::String) = join(
    "%" * uppercase(string(byte; base=16, pad=2)) for byte in codeunits(text))

function rejection_code(f::Function)
    try
        f()
    catch err
        @test err isa SAF.StructuredFirewallError
        return err isa SAF.StructuredFirewallError ? err.code : :wrong_exception
    end
    @test false
    return :not_rejected
end

function run_captured(cmd::Cmd)
    stdout = IOBuffer()
    stderr = IOBuffer()
    process = run(pipeline(cmd; stdout=stdout, stderr=stderr), wait=false)
    wait(process)
    return process.exitcode, String(take!(stdout)), String(take!(stderr))
end

function option_values(arguments::Vector{String})
    values = Dict{String,Vector{String}}()
    for argument in arguments
        startswith(argument, "--") || continue
        pair = split(argument, '='; limit=2)
        length(pair) == 2 || continue
        push!(get!(values, pair[1], String[]), pair[2])
    end
    return values
end

function probe_main(arguments::Vector{String})
    allowed = Set(["--out", "--feature", "--column", "--value", "--features",
                   "--interrupt-before-output"])
    valued = setdiff(allowed, Set(["--interrupt-before-output"]))
    try
        SAF.validate_cli_arguments(arguments;
            allowed_flags=allowed,
            value_flags=valued,
            repeatable_flags=Set(["--feature"]))
        options = option_values(arguments)
        out = only(get(options, "--out", String[]))
        selected = get(options, "--feature", ["amp_prominence"])
        column = only(get(options, "--column", ["amp_prominence"]))
        value = only(get(options, "--value", ["1.25"]))
        row = Dict{String,Any}("file" => "clean_scan.sxm", "lobe" => 1,
                               column => value)
        SAF.build_feature_batch([row], selected, ["amp_prominence"];
            required_metadata=["file", "lobe"])
        "--interrupt-before-output" in arguments && return 75
        open(out, "w") do io
            write(io, "validated\n")
        end
        return 0
    catch err
        showerror(stderr, err)
        println(stderr)
        return 2
    end
end

if !isempty(ARGS) && first(ARGS) == "--probe"
    exit(probe_main(ARGS[2:end]))
end

@testset "structured-assignment firewall" begin
    @testset "separate exact allowlists" begin
        @test SAF.JOIN_TOPOLOGY_METADATA == Set(["file", "lobe", "left_lobe", "right_lobe"])
        @test SAF.NUMERICAL_MODEL_FEATURES == Set(PHYSICAL)
        @test isempty(intersect(SAF.JOIN_TOPOLOGY_METADATA,
                                SAF.NUMERICAL_MODEL_FEATURES))
        @test SAF.validate_declared_features(PHYSICAL) == PHYSICAL
        for feature in PHYSICAL
            @test !SAF.is_forbidden_feature(feature)
        end
    end

    @testset "metadata never becomes numerical evidence" begin
        for metadata in sort!(collect(SAF.JOIN_TOPOLOGY_METADATA))
            @test SAF.is_join_metadata(metadata)
            @test SAF.is_forbidden_feature(metadata)
            @test rejection_code(() -> SAF.validate_feature_selection(
                [metadata], PHYSICAL)) == :metadata_as_feature
            @test rejection_code(() -> SAF.validate_declared_features(
                [metadata])) == :metadata_as_feature
        end
    end

    @testset "forbidden concepts survive spelling evasions" begin
        forbidden = [
            "benchmark_score", "truth_label", "grade_path", "grader_output",
            "expected_N", "target-n", "class_count", "composition_fraction",
            "top_k", "normalized_position", "ordinal_position", "terminal_pos",
            "centered_pos", "lobe_order", "lobe_index", "run_length",
            "state_frequency", "transition_prob", "switch_penalty", "sequence",
            "control_sequence", "nknnkn", "010010", "101101",
            "TrUtH_Label", "tr%75th_label", "tr%2575th_label",
            "tr\\u0075th_label", "tr&#x75;th_label", "ｔｒｕｔｈ_label",
            "trυth_label",
        ]
        for feature in forbidden
            @test SAF.is_forbidden_feature(feature)
            @test rejection_code(() -> SAF.validate_feature_selection(
                [feature], PHYSICAL)) == :forbidden_feature
        end
        @test rejection_code(() -> SAF.validate_feature_selection(
            ["amp_prominence", "amp_prominence"], PHYSICAL)) == :duplicate_feature
        @test rejection_code(() -> SAF.validate_feature_selection(
            ["harmless_but_undeclared"], PHYSICAL)) == :unknown_feature
    end

    @testset "CLI rejects forbidden and unknown controls" begin
        allowed = Set(["--config", "--features", "--out", "--view"])
        value_flags = copy(allowed)
        repeatable = Set(["--view"])
        @test isnothing(SAF.validate_cli_arguments(
            ["--config", "config/model.toml", "--features=inputs/features.tsv",
             "--view=base=amp_prominence,amp_rel", "--out", "out.tsv"];
            allowed_flags=allowed, value_flags=value_flags,
            repeatable_flags=repeatable))
        for arguments in (
            ["--truth=labels.tsv"],
            ["--TrUtH=labels.tsv"],
            ["--tr%75th=labels.tsv"],
            ["--ｔｒｕｔｈ=labels.tsv"],
            ["--features=/tmp/BENCHMARK.tsv"],
            ["--features", "/tmp/gr%61de.tsv"],
            ["--features=/tmp/trυth.tsv"],
            ["--view=base=state_frequency"],
        )
            @test rejection_code(() -> SAF.validate_cli_arguments(arguments;
                allowed_flags=allowed, value_flags=value_flags,
                repeatable_flags=repeatable)) == :forbidden_cli
        end
        @test rejection_code(() -> SAF.validate_cli_arguments(["--mystery=1"];
            allowed_flags=allowed, value_flags=value_flags,
            repeatable_flags=repeatable)) == :unknown_flag
        @test rejection_code(() -> SAF.validate_cli_arguments(
            ["--out=a.tsv", "--out=b.tsv"];
            allowed_flags=allowed, value_flags=value_flags,
            repeatable_flags=repeatable)) == :duplicate_flag
    end

    @testset "review correction regressions" begin
        @test !applicable(SAF.validate_feature_selection, ["amplitude"])

        @test rejection_code(() -> SAF.validate_cli_arguments(["--control"];
            allowed_flags=Set(["--control"]), value_flags=Set{String}())) ==
              :forbidden_cli
        @test rejection_code(() -> SAF.validate_cli_arguments(
            ["--features=/tmp/control.tsv"];
            allowed_flags=Set(["--features"]),
            value_flags=Set(["--features"]))) == :forbidden_cli

        for name in ("control", "expected_unit_count", "target_unit_count")
            @test SAF.is_forbidden_feature(name)
            @test rejection_code(() -> SAF.validate_table_header(
                ["file", "lobe", name], ["amplitude"])) == :forbidden_column

            metadata_row = Dict{String,Any}(
                "file" => "$(name)_scan.sxm", "lobe" => 1, "amplitude" => 0.1,
            )
            @test rejection_code(() -> SAF.build_feature_batch([metadata_row],
                ["amplitude"], ["amplitude"];
                required_metadata=["file", "lobe"])) == :forbidden_value

            provenance_row = Dict{String,Any}(
                "file" => "scan.sxm", "lobe" => 1, "amplitude" => 0.1,
                "source_path" => "/tmp/$name.tsv",
            )
            @test rejection_code(() -> SAF.build_feature_batch([provenance_row],
                ["amplitude"], ["amplitude"];
                required_metadata=["file", "lobe"],
                provenance_fields=["source_path"])) == :forbidden_value

            numerical_row = Dict{String,Any}(
                "file" => "scan.sxm", "lobe" => 1, "amplitude" => name,
            )
            @test rejection_code(() -> SAF.build_feature_batch([numerical_row],
                ["amplitude"], ["amplitude"];
                required_metadata=["file", "lobe"])) == :forbidden_value
        end

        five_level_truth = "tr" * nested_percent("%75", 4) * "th"
        six_level_control = nested_percent("%63ontrol", 5)
        mixed_backslash = "tr" * nested_percent("%5cu0075", 5) * "th"
        mixed_html = "tr" * nested_percent("%26%23117%3B", 5) * "th"
        for encoded in (five_level_truth, six_level_control,
                        mixed_backslash, mixed_html)
            @test SAF.is_forbidden_feature(encoded)
            @test rejection_code(() -> SAF.validate_cli_arguments(
                ["--features=/tmp/$encoded.tsv"];
                allowed_flags=Set(["--features"]),
                value_flags=Set(["--features"]))) == :forbidden_cli
        end

        @test isdefined(SAF, :_decode_fixed_point)
        if isdefined(SAF, :_decode_fixed_point)
            deeply_encoded = "tr" * nested_percent("%75", 6) * "th"
            @test rejection_code(() -> SAF._decode_fixed_point(
                deeply_encoded, 1)) == :decode_not_converged
        end

        clean_provenance = Dict{String,Any}(
            "file" => "scan.sxm", "lobe" => 1, "amplitude" => 0.1,
            "source_path" => "/tmp/features.tsv",
        )
        batch = SAF.build_feature_batch([clean_provenance],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"],
            provenance_fields=["source_path"])
        @test batch.matrix == reshape([0.1], 1, 1)
    end

    @testset "invalid percent UTF-8 fails closed" begin
        concept_payloads = [
            percent_encode("truth") * "%FF",
            percent_encode("control") * "%FF",
            percent_encode("expected_unit_count") * "%FF",
            percent_encode("target_unit_count") * "%FF",
        ]
        malformed_payloads = [
            "%80",
            "%C2",
            "%E2%82",
            "%C0%AF",
            "%ED%A0%80",
            "%F4%90%80%80",
        ]
        for payload in vcat(concept_payloads, malformed_payloads)
            @test rejection_code(() -> SAF._canonical_text(payload)) ==
                  :invalid_encoding
            @test rejection_code(() -> SAF.is_forbidden_feature(payload)) ==
                  :invalid_encoding
            @test rejection_code(() -> SAF.validate_cli_arguments(
                ["--features=/tmp/$payload.tsv"];
                allowed_flags=Set(["--features"]),
                value_flags=Set(["--features"]))) == :invalid_encoding
        end

        @test SAF._percent_decode_once("features%zz.tsv") == "features%zz.tsv"
        @test SAF._canonical_text("features%zz.tsv") == "features%zz.tsv"
        @test isnothing(SAF.validate_cli_arguments(
            ["--features=/tmp/features%zz.tsv"];
            allowed_flags=Set(["--features"]),
            value_flags=Set(["--features"])))
        @test SAF._canonical_text("features%2Etsv") == "features.tsv"

        valid_nested = "tr" * nested_percent("%75", 8) * "th"
        valid_mixed = "tr" * nested_percent("%26%23117%3B", 8) * "th"
        @test SAF._canonical_text(valid_nested) == "truth"
        @test SAF._canonical_text(valid_mixed) == "truth"
    end

    @testset "clean node and edge fixtures keep metadata out of matrices" begin
        nodes = [
            Dict{String,Any}("file" => "scan_A.sxm", "lobe" => 1,
                             "amplitude" => "0.10", "source_path" => "/tmp/features.tsv"),
            Dict{String,Any}("file" => "scan_A.sxm", "lobe" => 2,
                             "amplitude" => 0.20, "source_path" => "/tmp/features.tsv"),
        ]
        node_batch = SAF.build_feature_batch(nodes, ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"], provenance_fields=["source_path"])
        @test node_batch.feature_names == ["amplitude"]
        @test node_batch.matrix == reshape([0.10, 0.20], 2, 1)
        @test all(Set(keys(row)) == Set(["file", "lobe"])
                  for row in node_batch.metadata)

        edges = [Dict{String,Any}(
            "file" => "scan_A.sxm", "lobe" => 99,
            "left_lobe" => 1, "right_lobe" => 2,
            "corr_fwd" => 0.25, "corr_bwd" => -0.50,
        )]
        edge_batch = SAF.build_feature_batch(edges, ["corr_fwd", "corr_bwd"],
            ["corr_fwd", "corr_bwd"];
            required_metadata=["file", "lobe", "left_lobe", "right_lobe"])
        @test edge_batch.feature_names == ["corr_fwd", "corr_bwd"]
        @test edge_batch.matrix == reshape([0.25, -0.50], 1, 2)
        node_by_key = Dict((row["file"], parse(Int, row["lobe"])) =>
                           node_batch.matrix[index, 1]
                           for (index, row) in enumerate(node_batch.metadata))
        edge_meta = only(edge_batch.metadata)
        joined = (node_by_key[(edge_meta["file"], parse(Int, edge_meta["left_lobe"]))],
                  node_by_key[(edge_meta["file"], parse(Int, edge_meta["right_lobe"]))])
        @test joined == (0.10, 0.20)
        @test size(edge_batch.matrix, 2) == 2
    end

    @testset "columns, cells, and paths are fail closed" begin
        @test rejection_code(() -> SAF.validate_table_header(
            ["file", "lobe", "amplitude", "amplitude"], ["amplitude"])) ==
              :duplicate_column
        @test rejection_code(() -> SAF.validate_table_header(
            ["file", "lobe", "truth_label"], ["amplitude"])) == :forbidden_column
        @test rejection_code(() -> SAF.validate_table_header(
            ["file", "lobe", "undeclared"], ["amplitude"])) == :unknown_column

        path_row = Dict{String,Any}("file" => "scan.sxm", "lobe" => 1,
                                    "amp_prominence" => "../secret/features.tsv")
        @test rejection_code(() -> SAF.build_feature_batch([path_row],
            ["amp_prominence"], ["amp_prominence"];
            required_metadata=["file", "lobe"])) == :path_value
        basename_path = Dict{String,Any}("file" => "scan.sxm", "lobe" => 1,
                                         "amp_prominence" => "secret.tsv")
        @test rejection_code(() -> SAF.build_feature_batch([basename_path],
            ["amp_prominence"], ["amp_prominence"];
            required_metadata=["file", "lobe"])) == :path_value
        escaped_file = Dict{String,Any}("file" => "../scan.sxm", "lobe" => 1,
                                        "amplitude" => 0.1)
        @test rejection_code(() -> SAF.build_feature_batch([escaped_file],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"])) == :path_value
        encoded_file = Dict{String,Any}("file" => "%2e%2e%2fscan.sxm", "lobe" => 1,
                                        "amplitude" => 0.1)
        @test rejection_code(() -> SAF.build_feature_batch([encoded_file],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"])) == :path_value
        truth_path = Dict{String,Any}("file" => "scan.sxm", "lobe" => 1,
                                      "amplitude" => 0.1,
                                      "source_path" => "/tmp/truth.tsv")
        @test rejection_code(() -> SAF.build_feature_batch([truth_path],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"],
            provenance_fields=["source_path"])) == :forbidden_value
        unicode_truth_path = Dict{String,Any}("file" => "scan.sxm", "lobe" => 1,
                                              "amplitude" => 0.1,
                                              "source_path" => "/tmp/trυth.tsv")
        @test rejection_code(() -> SAF.build_feature_batch([unicode_truth_path],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"],
            provenance_fields=["source_path"])) == :forbidden_value
        @test rejection_code(() -> SAF.build_feature_batch([truth_path],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"],
            provenance_fields=["truth_path"])) == :forbidden_provenance_field
        nonfinite = Dict{String,Any}("file" => "scan.sxm", "lobe" => 1,
                                     "amplitude" => "Inf")
        @test rejection_code(() -> SAF.build_feature_batch([nonfinite],
            ["amplitude"], ["amplitude"];
            required_metadata=["file", "lobe"])) == :nonfinite_feature
    end

    @testset "bad CLI/data exits before output creation" begin
        mktempdir() do temporary
            script = @__FILE__
            julia = `$(Base.julia_cmd()) --project=$ROOT $script --probe`
            cases = [
                ["--feature=lobe"],
                ["--tr%75th=x.tsv"],
                ["--ｔｒｕｔｈ=x.tsv"],
                ["--features=/tmp/GrAdE.tsv"],
                ["--features=/tmp/control.tsv"],
                ["--features=/tmp/tr%2525252575th.tsv"],
                ["--features=/tmp/expected_unit_count.tsv"],
                ["--features=/tmp/target_unit_count.tsv"],
                ["--features=/tmp/$(nested_percent("%63ontrol", 5)).tsv"],
                ["--features=/tmp/$("tr" * nested_percent("%5cu0075", 5) * "th").tsv"],
                ["--features=/tmp/$("tr" * nested_percent("%26%23117%3B", 5) * "th").tsv"],
                ["--features=/tmp/$(percent_encode("truth"))%FF.tsv"],
                ["--features=/tmp/$(percent_encode("control"))%FF.tsv"],
                ["--features=/tmp/$(percent_encode("expected_unit_count"))%FF.tsv"],
                ["--features=/tmp/$(percent_encode("target_unit_count"))%FF.tsv"],
                ["--features=/tmp/%80.tsv"],
                ["--features=/tmp/%C2.tsv"],
                ["--features=/tmp/%E2%82.tsv"],
                ["--features=/tmp/%C0%AF.tsv"],
                ["--features=/tmp/%ED%A0%80.tsv"],
                ["--features=/tmp/%F4%90%80%80.tsv"],
                ["--column=state_frequency", "--value=0.5"],
                ["--column=amp_prominence", "--value=../escape.tsv"],
            ]
            for (index, extra) in enumerate(cases)
                output = joinpath(temporary, "forbidden-$index.tsv")
                command = `$julia --out=$output $extra`
                code, _, error_text = run_captured(command)
                @test code != 0
                @test !ispath(output)
                @test !isempty(error_text)
            end
            output = joinpath(temporary, "clean.tsv")
            code, _, error_text = run_captured(
                `$julia --out=$output --feature=amp_prominence`)
            @test code == 0
            @test isfile(output)
            @test read(output, String) == "validated\n"
            @test isempty(error_text)

            for replay in 1:2
                interrupted = joinpath(temporary, "interrupted-$replay.tsv")
                code, _, _ = run_captured(
                    `$julia --out=$interrupted --interrupt-before-output`)
                @test code == 75
                @test !ispath(interrupted)
            end
        end
    end
end
