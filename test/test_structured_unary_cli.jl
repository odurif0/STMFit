#!/usr/bin/env julia

using Printf
using Random
using SHA
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
const LIB_SOURCE = joinpath(@__DIR__, "lib", "structured_unit_assignment.jl")
const CLI_SOURCE = joinpath(@__DIR__, "build_structured_unit_predictions.jl")
const MODEL_CONFIG = joinpath(ROOT, "config", "unit_assignment_structured_model.toml")
const CHAMPION_SOURCE = joinpath(
    ROOT, "results", "unit_assignment", "best_labelfree_cc_soft_20260802.tsv")

const OUTPUT_COLUMNS = [
    "file", "lobe", "predicted", "probability_1", "confidence", "model",
    "status", "invalid_reason", "views_used", "segment_id", "edge_status",
    "edge_admission_hash", "trigger_hash", "source_hash", "config_hash",
    "input_hash",
]
const VIEW_FEATURES = [
    "amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel",
    "bwd_neg_com_t", "bwd_neg_diag45", "split_log_skew",
]
const INPUT_COLUMNS = vcat(["file", "lobe", "amplitude"], VIEW_FEATURES)
const MODES = ["gaussian_one", "gaussian_two", "student_t_two"]
const EXPECTED_IDENTITIES = Dict(
    "config" => "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad",
    "candidate" => "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd",
    "champion_review" => "6ff61299a9a9145df0ed21aab14f4d57ed3ff60f762158dc30c4a87fa667dd7f",
    "robust_review" => "8b0fb2c1a324a345de19d4e718278a52a6d568c0e474939cfb6bfe643e061635",
    "robust_source" => "0ca4a3863a4b9ed13aee5cd3a0229df4f2764a8524b764eeb0fd27238f14f241",
    "champion_source" => "86c7b0985413cf4ceb61ee903ed9bd27e375b5cc2ed314fa1b15b9648dd4ecf5",
)
const UNCHANGED_SOURCES = Dict(
    "test/run_unknown_unit_assignment.jl" =>
        "62637365bb0483117c04183343cc9958d99ee3feb46444f3f8ede3aea5738f82",
    "test/validate_unit_predictions.jl" =>
        "13a4fca86db4420e659fda708f91e4808030abf3819b035e94a5c36db2667b48",
    "test/summarize_unknown_unit_qc.jl" =>
        "cf063140cabd0f5c34032b01b3308cff06c2278802e0d80ac186ff57b534ddfe",
    "test/lib/hierarchical_unit_assignment.jl" =>
        "4024d47b82077436a1c22a9511930effc5bca110f4c9eba5daef487b923513a8",
    "test/lib/hierarchical/loading.jl" =>
        "41171af123d489f8e3d4e20ba0ad2e55641a3cbb5768b3debc711c3f8262422a",
    "test/lib/hierarchical/nuisance.jl" =>
        "f256eb82966350544b55a23413417aa6eacda7cf6669c049f8dc90feb4fd2203",
    "test/lib/hierarchical/views.jl" =>
        "40380fb586a7866e718225eac3fd610e4eba9a0e0c3bddb924f108655e5fdfcf",
    "test/build_hierarchical_unit_predictions.jl" =>
        "5d810dc3fcd70d0363ba31961604062d0f20f101abdfd7f258dd853bfedcb451",
)

sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))
file_sha256(path::AbstractString) = sha256_hex(read(path))

missing_products = filter(path -> !isfile(path), [LIB_SOURCE, CLI_SOURCE])
if !isempty(missing_products)
    @testset "structured unary CLI RED: product sources are required" begin
        @test isfile(LIB_SOURCE)
        @test isfile(CLI_SOURCE)
    end
    exit(1)
end

include(LIB_SOURCE)
using .StructuredUnitAssignment
const SUA = StructuredUnitAssignment

include(CLI_SOURCE)
const SUCLI = StructuredUnaryPredictionCLI

function _row(file::String, lobe::Int, amplitude::Float64, latent::Float64;
              missing_base::Bool=false, missing_extras::Bool=false)
    jitter = 0.025 * sin(0.73 * lobe) + 0.018 * cos(0.41 * lobe)
    values = Dict{String,String}(
        "file" => file,
        "lobe" => string(lobe),
        "amplitude" => @sprintf("%.17g", amplitude),
        "amp_prominence" => @sprintf("%.17g", latent + jitter),
        "amp_neighbor_ratio" => @sprintf("%.17g", 0.61 * latent - 0.4 * jitter),
        "integrated_prominence" => @sprintf("%.17g", 1.17 * latent + 0.7 * jitter),
        "amp_rel" => @sprintf("%.17g", 0.83 * latent - 0.2 * jitter),
        "bwd_neg_com_t" => @sprintf("%.17g", 0.47 * latent + 0.9 * jitter),
        "bwd_neg_diag45" => @sprintf("%.17g", -0.39 * latent + 0.8 * jitter),
        "split_log_skew" => @sprintf("%.17g", 0.29 * latent - 1.1 * jitter),
    )
    if missing_base
        for name in VIEW_FEATURES[1:4]
            values[name] = "NA"
        end
    end
    if missing_extras
        for name in VIEW_FEATURES[5:7]
            values[name] = "NA"
        end
    end
    return values
end

function separated_rows(; amplitude_mode::Symbol=:normal,
                        missing_base::Bool=false,
                        missing_extras::Bool=false)
    rows = Dict{String,String}[]
    for (scan_index, file) in enumerate(("scan_alpha.sxm", "scan_beta.sxm"))
        for lobe in 1:24
            high = lobe > 12
            latent = (high ? 2.4 : -2.2) + 0.07 * sin(0.31 * lobe + scan_index)
            low_amplitude = 1.0 + 0.025 * sin(lobe + scan_index)
            high_amplitude = 3.0 + 0.025 * cos(lobe + scan_index)
            amplitude = if amplitude_mode == :normal
                high ? high_amplitude : low_amplitude
            elseif amplitude_mode == :reversed
                high ? low_amplitude : high_amplitude
            elseif amplitude_mode == :tied
                2.0
            else
                error("unsupported fixture amplitude mode")
            end
            push!(rows, _row(file, lobe, amplitude, latent;
                             missing_base=missing_base,
                             missing_extras=missing_extras))
        end
    end
    return rows
end

function rows_with_constant_required_view()
    rows = separated_rows()
    for row in rows
        row["bwd_neg_com_t"] = "0"
    end
    return rows
end

function rows_with_missing_required_view()
    rows = separated_rows()
    for row in rows
        row["bwd_neg_com_t"] = "NA"
    end
    return rows
end

function one_population_rows()
    rows = Dict{String,String}[]
    rng = MersenneTwister(808)
    for (scan_index, file) in enumerate(("scan_gamma.sxm", "scan_delta.sxm"))
        for lobe in 1:90
            row = Dict{String,String}(
                "file" => file,
                "lobe" => string(lobe),
                "amplitude" => @sprintf("%.17g", 2.0 + 0.15 * randn(rng)),
            )
            for feature in VIEW_FEATURES
                row[feature] = @sprintf("%.17g", randn(rng))
            end
            push!(rows, row)
        end
    end
    return rows
end

function write_features(path::AbstractString, rows;
                        columns::Vector{String}=copy(INPUT_COLUMNS),
                        extra_cell::Bool=false,
                        missing_cell::Bool=false)
    open(path, "w") do io
        println(io, join(columns, '\t'))
        for (row_index, row) in enumerate(rows)
            values = [get(row, column, "") for column in columns]
            row_index == 1 && extra_cell && push!(values, "unexpected")
            row_index == 1 && missing_cell && pop!(values)
            println(io, join(values, '\t'))
        end
    end
    return path
end

function read_table(path::AbstractString)
    lines = readlines(path)
    @test !isempty(lines)
    header = String.(split(first(lines), '\t'; keepempty=true))
    rows = Dict{String,String}[]
    for line in lines[2:end]
        values = String.(split(line, '\t'; keepempty=true))
        @test length(values) == length(header)
        push!(rows, Dict(header[index] => values[index] for index in eachindex(header)))
    end
    return header, rows
end

function run_captured(command::Cmd)
    stdout = IOBuffer()
    stderr = IOBuffer()
    process = run(pipeline(command; stdout=stdout, stderr=stderr), wait=false)
    wait(process)
    return process.exitcode, String(take!(stdout)), String(take!(stderr))
end

function cli_command(arguments::Vector{String})
    base = `$(Base.julia_cmd()) --project=$ROOT $CLI_SOURCE`
    return Cmd(vcat(base.exec, arguments))
end

function run_entrypoint(arguments::Vector{String})
    code = redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            SUCLI.entrypoint(arguments)
        end
    end
    return code, "", ""
end

function run_mode(features::String, output::String, mode::String)
    return run_entrypoint([
        "--config", MODEL_CONFIG,
        "--features", features,
        "--out", output,
        "--model", mode,
    ])
end

function captured_error(action::Function)
    try
        action()
        return nothing
    catch error
        return error
    end
end

function owner_snapshot(path::String)
    if islink(path)
        return (:symlink, lstat(path).inode, readlink(path))
    elseif isdir(path)
        return (:directory, stat(path).inode, sort(readdir(path)))
    elseif isfile(path)
        return (:file, stat(path).inode, read(path))
    end
    return (:absent, UInt(0), nothing)
end

function create_owner(path::String, kind::Symbol, token::String)
    if kind == :file
        write(path, "OWNER-$token\n")
    elseif kind == :symlink
        target = joinpath(dirname(path), "target-$token.txt")
        write(target, "TARGET-$token\n")
        symlink(target, path)
    elseif kind == :directory
        mkdir(path)
        write(joinpath(path, "marker.txt"), "DIRECTORY-$token\n")
    else
        error("unsupported owner kind")
    end
    return owner_snapshot(path)
end

function rows_by_key(rows)
    return Dict((row["file"], parse(Int, row["lobe"])) => row for row in rows)
end

function parsed_probability(row)
    return row["probability_1"] == "NA" ? nothing : parse(Float64, row["probability_1"])
end

function assert_output_contract(path::String, mode::String, input_hash::String)
    header, rows = read_table(path)
    @test header == OUTPUT_COLUMNS
    @test length(rows) == 48
    @test [(row["file"], parse(Int, row["lobe"])) for row in rows] ==
          sort([(row["file"], parse(Int, row["lobe"])) for row in rows])
    @test length(Set((row["file"], row["lobe"]) for row in rows)) == length(rows)
    for row in rows
        @test !isempty(row["file"])
        @test parse(Int, row["lobe"]) > 0
        @test row["predicted"] in ("0", "1", "?")
        probability = parsed_probability(row)
        @test probability === nothing || (probability isa Float64 && 0.0 <= probability <= 1.0)
        confidence = parse(Float64, row["confidence"])
        @test isfinite(confidence)
        @test 0.0 <= confidence <= 1.0
        @test row["model"] == mode
        @test row["status"] in ("ok", "abstained", "failed")
        @test !isempty(row["invalid_reason"])
        @test 0 <= parse(Int, row["views_used"]) <= 4
        @test row["segment_id"] == "pending"
        @test row["edge_status"] == "pending"
        @test row["edge_admission_hash"] == "pending"
        @test occursin(r"^[0-9a-f]{64}$", row["trigger_hash"])
        @test occursin(r"^[0-9a-f]{64}$", row["source_hash"])
        @test row["config_hash"] == EXPECTED_IDENTITIES["config"]
        @test row["input_hash"] == input_hash
        if row["status"] == "ok"
            @test row["predicted"] in ("0", "1")
            @test probability isa Float64
            @test confidence == abs(2.0 * probability - 1.0)
        elseif row["status"] == "abstained"
            @test row["predicted"] == "?"
            @test probability == 0.5
            @test iszero(confidence)
        else
            @test row["predicted"] == "?"
            @test probability === nothing
            @test iszero(confidence)
        end
    end
    @test length(unique(row["trigger_hash"] for row in rows)) == 1
    @test length(unique(row["source_hash"] for row in rows)) == 1
    return rows
end

function champion_stream_hash(stream)
    io = IOBuffer()
    for record in stream
        print(io, record.key.file, '\t', record.key.lobe, '\t',
              record.unary.predicted, '\t', record.unary.confidence_text, '\t',
              bitstring(record.unary.probability_1), '\n')
    end
    return sha256_hex(take!(io))
end

const REPLAY_SHA256 = Ref("")

@testset "structured unary facade and C0/C1/C2 CLI" begin
    @testset "immutable dependency identities and exact reused APIs" begin
        @test file_sha256(MODEL_CONFIG) == EXPECTED_IDENTITIES["config"]
        @test file_sha256(joinpath(ROOT, "config", "unit_assignment_structured_candidate.toml")) ==
              EXPECTED_IDENTITIES["candidate"]
        @test file_sha256(joinpath(
            ROOT, ".omo", "evidence", "structured-label-free-unit-assignment", "t6",
            "review", "AdversarialVerify.json")) == EXPECTED_IDENTITIES["champion_review"]
        @test file_sha256(joinpath(
            ROOT, ".omo", "evidence", "structured-label-free-unit-assignment", "t7",
            "correction", "review", "AdversarialVerify.json")) ==
              EXPECTED_IDENTITIES["robust_review"]
        @test file_sha256(joinpath(
            @__DIR__, "lib", "structured_assignment", "robust_emissions.jl")) ==
              EXPECTED_IDENTITIES["robust_source"]
        @test Tuple(SUA.SUPPORTED_MODELS) == Tuple(MODES)
        @test SUA.OUTPUT_COLUMNS == OUTPUT_COLUMNS

        adapter_source = read(joinpath(
            @__DIR__, "lib", "structured_assignment", "champion_adapter.jl"), String)
        robust_source = read(joinpath(
            @__DIR__, "lib", "structured_assignment", "robust_emissions.jl"), String)
        facade_source = read(LIB_SOURCE, String)
        @test occursin("champion_adapter.jl", facade_source)
        @test occursin("robust_emissions.jl", facade_source)
        @test !occursin("function load_frozen_champion_unaries", facade_source)
        @test !occursin("function log_student_t_diag", facade_source)
        @test occursin("function load_frozen_champion_unaries", adapter_source)
        @test occursin("function log_student_t_diag", robust_source)
    end

    @testset "Todo 6 stream round-trips through the facade" begin
        @test file_sha256(CHAMPION_SOURCE) == EXPECTED_IDENTITIES["champion_source"]
        stream = SUA.load_champion_unaries()
        @test length(stream) == 892
        @test stream.source.sha256 == EXPECTED_IDENTITIES["champion_source"]
        @test champion_stream_hash(stream) ==
              "5f8aca367cf5fd71e014dda1ad4b1324a5d0307d24dbfcfb78f0217fac938b4c"
        @test SUA.load_champion_unaries() == stream
    end

    @testset "all three modes, exact schema, deterministic replay, and equal fusion" begin
        mktempdir() do temporary
            features = write_features(joinpath(temporary, "features.tsv"), separated_rows())
            input_hash = file_sha256(features)
            bytes_by_mode = Dict{String,Vector{UInt8}}()
            parsed_by_mode = Dict{String,Vector{Dict{String,String}}}()
            for mode in MODES
                first_output = joinpath(temporary, "$(mode)-first.tsv")
                second_output = joinpath(temporary, "$(mode)-second.tsv")
                first_code, _, first_stderr = run_mode(features, first_output, mode)
                second_code, _, second_stderr = run_mode(features, second_output, mode)
                @test first_code == 0
                @test second_code == 0
                @test isempty(first_stderr)
                @test isempty(second_stderr)
                @test read(first_output) == read(second_output)
                bytes_by_mode[mode] = read(first_output)
                parsed_by_mode[mode] = assert_output_contract(first_output, mode, input_hash)
            end

            c0 = parsed_by_mode["gaussian_one"]
            @test all(row["status"] == "abstained" for row in c0)
            @test all(row["predicted"] == "?" for row in c0)
            @test all(row["probability_1"] == "0.5" for row in c0)
            @test all(row["invalid_reason"] == "one_component_evidence" for row in c0)
            @test all(row["views_used"] == "0" for row in c0)

            for mode in ("gaussian_two", "student_t_two")
                rows = parsed_by_mode[mode]
                @test all(row["status"] == "ok" for row in rows)
                @test all(row["views_used"] == "4" for row in rows)
                @test Set(row["predicted"] for row in rows) == Set(["0", "1"])
            end

            contract = SUA.load_contract(MODEL_CONFIG)
            records = SUA.HierarchicalUnitAssignment.load_records(features)
            direct_fits = [SUA.HierarchicalUnitAssignment.fit_view(
                records, spec.second, spec.first;
                first_seed=0,
                n_starts=contract.n_starts,
                cov_floor=contract.covariance_floor,
                max_iter=contract.max_iter,
                tol=contract.tolerance,
            ) for spec in contract.views]
            direct = SUA.HierarchicalUnitAssignment.combine_views(direct_fits, length(records))
            c1 = parsed_by_mode["gaussian_two"]
            for (index, row) in enumerate(c1)
                @test row["file"] == records[index].file
                @test parse(Int, row["lobe"]) == records[index].lobe
                @test parse(Float64, row["probability_1"]) == direct.probability_1[index]
                logits = [fit.log_odds_for_label1[index] for fit in direct_fits]
                @test all(isfinite, logits)
                expected = 1.0 / (1.0 + exp(-sum(logits) / 4))
                @test parse(Float64, row["probability_1"]) == expected
            end

            replay_material = IOBuffer()
            for mode in MODES
                println(replay_material, mode)
                write(replay_material, bytes_by_mode[mode])
            end
            REPLAY_SHA256[] = sha256_hex(take!(replay_material))
            @test occursin(r"^[0-9a-f]{64}$", REPLAY_SHA256[])
        end
    end

    @testset "amplitude-only orientation and corrected exact-tie abstention" begin
        mktempdir() do temporary
            ordinary_path = write_features(
                joinpath(temporary, "ordinary.tsv"), separated_rows())
            reversed_path = write_features(
                joinpath(temporary, "reversed.tsv"),
                separated_rows(amplitude_mode=:reversed))
            tied_path = write_features(
                joinpath(temporary, "tied.tsv"), separated_rows(amplitude_mode=:tied))
            for mode in ("gaussian_two", "student_t_two")
                ordinary = SUA.build_predictions(MODEL_CONFIG, ordinary_path, mode)
                reversed = SUA.build_predictions(MODEL_CONFIG, reversed_path, mode)
                @test [(row.file, row.lobe) for row in ordinary.rows] ==
                      [(row.file, row.lobe) for row in reversed.rows]
                for (left, right) in zip(ordinary.rows, reversed.rows)
                    @test left.status == right.status == "ok"
                    @test left.predicted != right.predicted
                    @test left.probability_1 + right.probability_1 ≈ 1.0 atol=2e-14 rtol=0
                end

                tied = SUA.build_predictions(MODEL_CONFIG, tied_path, mode)
                @test all(row.status == "abstained" for row in tied.rows)
                @test all(row.predicted == "?" for row in tied.rows)
                @test all(row.probability_1 == 0.5 for row in tied.rows)
                @test all(row.confidence == 0.0 for row in tied.rows)
                @test all(row.invalid_reason == "amplitude_orientation_undefined"
                          for row in tied.rows)
            end
        end
    end

    @testset "partial, absent, and one-component view evidence fails closed" begin
        mktempdir() do temporary
            partial_path = write_features(
                joinpath(temporary, "partial.tsv"),
                separated_rows(missing_extras=true))
            absent_path = write_features(
                joinpath(temporary, "absent.tsv"),
                separated_rows(missing_base=true, missing_extras=true))
            one_path = write_features(
                joinpath(temporary, "one.tsv"), one_population_rows())
            for mode in ("gaussian_two", "student_t_two")
                partial = SUA.build_predictions(MODEL_CONFIG, partial_path, mode)
                @test all(row.status == "ok" for row in partial.rows)
                @test all(row.views_used == 1 for row in partial.rows)
                @test all(row.invalid_reason == "ok_partial_views" for row in partial.rows)

                absent = SUA.build_predictions(MODEL_CONFIG, absent_path, mode)
                @test all(row.status == "abstained" for row in absent.rows)
                @test all(row.predicted == "?" for row in absent.rows)
                @test all(row.probability_1 == 0.5 for row in absent.rows)
                @test all(row.views_used == 0 for row in absent.rows)
                @test all(row.invalid_reason in ("degenerate_view", "missing_view")
                          for row in absent.rows)

                one = SUA.build_predictions(MODEL_CONFIG, one_path, mode)
                @test all(row.status == "abstained" for row in one.rows)
                @test all(row.predicted == "?" for row in one.rows)
                @test all(row.probability_1 == 0.5 for row in one.rows)
                @test all(occursin("one_component", row.invalid_reason) ||
                          row.invalid_reason == "amplitude_orientation_undefined"
                          for row in one.rows)
            end
        end
    end

    @testset "missing views marginalize but invalid available required views abstain" begin
        mktempdir() do temporary
            missing_path = write_features(
                joinpath(temporary, "missing-required.tsv"),
                rows_with_missing_required_view())
            invalid_path = write_features(
                joinpath(temporary, "invalid-required.tsv"),
                rows_with_constant_required_view())
            for mode in ("gaussian_two", "student_t_two")
                missing = SUA.build_predictions(MODEL_CONFIG, missing_path, mode)
                @test all(row.status == "ok" for row in missing.rows)
                @test all(row.views_used == 3 for row in missing.rows)
                @test all(row.invalid_reason == "ok_partial_views" for row in missing.rows)

                invalid = SUA.build_predictions(MODEL_CONFIG, invalid_path, mode)
                @test all(row.status == "abstained" for row in invalid.rows)
                @test all(row.predicted == "?" for row in invalid.rows)
                @test all(row.probability_1 == 0.5 for row in invalid.rows)
                @test all(row.confidence == 0.0 for row in invalid.rows)
                @test all(row.views_used == 0 for row in invalid.rows)
                @test all(row.invalid_reason == "invalid_required_view_normalization"
                          for row in invalid.rows)

                output = joinpath(temporary, "invalid-required-$mode.tsv")
                code, _, stderr = run_mode(invalid_path, output, mode)
                @test code == 0
                @test isempty(stderr)
                _, cli_rows = read_table(output)
                @test length(cli_rows) == length(invalid.rows)
                @test all(row["status"] == "abstained" for row in cli_rows)
                @test all(row["invalid_reason"] ==
                          "invalid_required_view_normalization" for row in cli_rows)
                @test all(row["views_used"] == "0" for row in cli_rows)
            end
        end
    end

    @testset "strict feature input failures happen before output" begin
        mktempdir() do temporary
            rows = separated_rows()
            function rejected_fixture(writer::Function, name::String)
                features = joinpath(temporary, "$name.tsv")
                writer(features)
                output = joinpath(temporary, "$name-output.tsv")
                sentinel = collect(codeunits("preexisting-output\n"))
                write(output, sentinel)
                code, _, _ = run_mode(features, output, "gaussian_two")
                @test code != 0
                @test read(output) == sentinel
            end

            rejected_fixture("missing-column") do path
                write_features(path, rows; columns=INPUT_COLUMNS[1:end-1])
            end
            rejected_fixture("duplicate-column") do path
                write_features(path, rows; columns=vcat(INPUT_COLUMNS, [last(INPUT_COLUMNS)]))
            end
            rejected_fixture("extra-column") do path
                changed = deepcopy(rows)
                for row in changed
                    row["sigma_parallel_nm"] = "0.4"
                end
                write_features(path, changed; columns=vcat(INPUT_COLUMNS, ["sigma_parallel_nm"]))
            end
            rejected_fixture("forbidden-column") do path
                forbidden = join(("class", "count"), "_")
                changed = deepcopy(rows)
                for row in changed
                    row[forbidden] = "1"
                end
                write_features(path, changed; columns=vcat(INPUT_COLUMNS, [forbidden]))
            end
            rejected_fixture("duplicate-key") do path
                changed = deepcopy(rows)
                changed[2]["file"] = changed[1]["file"]
                changed[2]["lobe"] = changed[1]["lobe"]
                write_features(path, changed)
            end
            rejected_fixture("nonfinite") do path
                changed = deepcopy(rows)
                changed[1]["amp_prominence"] = "Inf"
                write_features(path, changed)
            end
            rejected_fixture("extra-cell") do path
                write_features(path, rows; extra_cell=true)
            end
            rejected_fixture("missing-cell") do path
                write_features(path, rows; missing_cell=true)
            end

            missing_features = joinpath(temporary, "missing.tsv")
            absent_output = joinpath(temporary, "absent-output.tsv")
            code, _, _ = run_mode(missing_features, absent_output, "gaussian_two")
            @test code != 0
            @test !ispath(absent_output)

            good_features = write_features(joinpath(temporary, "good.tsv"), rows)
            missing_config = joinpath(temporary, "missing.toml")
            config_output = joinpath(temporary, "config-output.tsv")
            code, _, _ = run_entrypoint([
                "--config", missing_config,
                "--features", good_features,
                "--out", config_output,
                "--model", "gaussian_two",
            ])
            @test code != 0
            @test !ispath(config_output)
        end
    end

    @testset "prediction publication is atomic no-replace for every destination owner" begin
        mktempdir() do temporary
            features = write_features(joinpath(temporary, "features.tsv"), separated_rows())
            bundle = SUA.build_predictions(MODEL_CONFIG, features, "gaussian_two")

            for kind in (:file, :symlink, :directory)
                destination = joinpath(temporary, "existing-$kind.tsv")
                before = create_owner(destination, kind, "existing-$kind")
                directory_before = Set(readdir(temporary))
                error = captured_error() do
                    SUA.write_prediction_tsv(destination, bundle)
                end
                @test error isa SUA.StructuredUnaryError
                if error isa SUA.StructuredUnaryError
                    @test error.code == :output_collision
                end
                @test owner_snapshot(destination) == before
                @test Set(readdir(temporary)) == directory_before
            end

            cli_destination = joinpath(temporary, "existing-cli.tsv")
            cli_before = create_owner(cli_destination, :file, "existing-cli")
            cli_entries = Set(readdir(temporary))
            code, _, _ = run_mode(features, cli_destination, "gaussian_two")
            @test code == 2
            @test owner_snapshot(cli_destination) == cli_before
            @test Set(readdir(temporary)) == cli_entries

            @test isdefined(SUA, :_before_noreplace_hook)
            if isdefined(SUA, :_before_noreplace_hook)
                for kind in (:file, :symlink, :directory)
                    destination = joinpath(temporary, "racing-$kind.tsv")
                    before_entries = Set(readdir(temporary))
                    competitor = Ref{Any}(nothing)
                    SUA._before_noreplace_hook[] = (_stage, path) -> begin
                        competitor[] = create_owner(path, kind, "racing-$kind")
                    end
                    error = try
                        captured_error() do
                            SUA.write_prediction_tsv(destination, bundle)
                        end
                    finally
                        SUA._before_noreplace_hook[] = nothing
                    end
                    @test error isa SUA.StructuredUnaryError
                    if error isa SUA.StructuredUnaryError
                        @test error.code == :output_collision
                    end
                    @test competitor[] !== nothing
                    if competitor[] !== nothing
                        @test owner_snapshot(destination) == competitor[]
                    end
                    additions = kind == :symlink ?
                        Set([basename(destination), "target-racing-symlink.txt"]) :
                        Set([basename(destination)])
                    expected_entries = union(before_entries, additions)
                    @test Set(readdir(temporary)) == expected_entries
                end
            end

            @test isdefined(SUCLI.SUA, :_before_noreplace_hook)
            if isdefined(SUCLI.SUA, :_before_noreplace_hook)
                race_destination = joinpath(temporary, "racing-cli.tsv")
                race_before = Ref{Any}(nothing)
                before_entries = Set(readdir(temporary))
                SUCLI.SUA._before_noreplace_hook[] = (_stage, path) -> begin
                    race_before[] = create_owner(path, :file, "racing-cli")
                end
                code = try
                    first(run_mode(features, race_destination, "student_t_two"))
                finally
                    SUCLI.SUA._before_noreplace_hook[] = nothing
                end
                @test code == 2
                @test race_before[] !== nothing
                if race_before[] !== nothing
                    @test owner_snapshot(race_destination) == race_before[]
                end
                expected_entries = union(before_entries, Set([basename(race_destination)]))
                @test Set(readdir(temporary)) == expected_entries
            end
        end
    end

    @testset "CLI has only the frozen controls and rejects forbidden controls" begin
        @test_throws Exception SUCLI.parse_cli(String[])
        @test_throws Exception SUCLI.parse_cli([
            "--config", MODEL_CONFIG, "--features", "features.tsv", "--out", "out.tsv",
        ])
        @test_throws Exception SUCLI.parse_cli([
            "--config", MODEL_CONFIG, "--features", "features.tsv", "--out", "out.tsv",
            "--model", "alternate",
        ])

        mktempdir() do temporary
            features = write_features(joinpath(temporary, "features.tsv"), separated_rows())
            controls = [
                "--" * "gra" * "der",
                "--" * "bench" * "mark",
                "--" * "class" * "-count",
                "--" * "transition",
                "--" * "run" * "-length",
                "--" * "composition",
                "--" * "pair" * "-prior",
                "--" * "observation",
                "--view",
            ]
            for (index, flag) in enumerate(controls)
                output = joinpath(temporary, "rejected-$index.tsv")
                sentinel = collect(codeunits("unchanged\n"))
                write(output, sentinel)
                arguments = [
                    "--config", MODEL_CONFIG,
                    "--features", features,
                    "--out", output,
                    "--model", "gaussian_two",
                    flag, "value",
                ]
                code, _, _ = run_entrypoint(arguments)
                @test code != 0
                @test read(output) == sentinel
            end
        end
    end

    @testset "help is explicit and label-free" begin
        code, stdout, stderr = run_captured(cli_command(["--help"]))
        @test code == 0
        @test isempty(stderr)
        help = lowercase(stdout)
        @test occursin("label-free", help)
        for flag in ("--config", "--features", "--out", "--model")
            @test occursin(flag, help)
        end
        for mode in MODES
            @test occursin(mode, help)
        end
        for fragment in (
            "gra" * "der", "bench" * "mark", "class count", "transition",
            "run-length", "composition", "pair-prior", "observation", "--view",
        )
            @test !occursin(fragment, help)
        end
    end

    @testset "historical default and v1 implementation sources remain exact" begin
        for (relative, expected) in UNCHANGED_SOURCES
            @test file_sha256(joinpath(ROOT, relative)) == expected
        end
    end
end

println("T8_REPLAY_SHA256=", REPLAY_SHA256[])
println("T8_CONFIG_SHA256=", file_sha256(MODEL_CONFIG))
println("T8_SOURCE_SHA256=", SUA.source_hash())
