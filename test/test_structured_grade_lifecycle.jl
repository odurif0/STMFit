#!/usr/bin/env julia

using Printf
using SHA
using TOML
using Test

const ROOT = realpath(dirname(@__DIR__))
const FREEZE_WRAPPER = joinpath(@__DIR__, "freeze_structured_unit_assignment_candidate.jl")
const GRADE_WRAPPER = joinpath(@__DIR__, "grade_frozen_structured_candidate.jl")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_structured_candidate.toml")
const MODEL = joinpath(ROOT, "config", "unit_assignment_structured_model.toml")
const T5_EVIDENCE = joinpath(
    ROOT, ".omo", "evidence", "structured-label-free-unit-assignment", "t5")

@testset "structured freeze and one-shot grade contract" begin
    @test isfile(FREEZE_WRAPPER)
    @test isfile(GRADE_WRAPPER)

    if isfile(FREEZE_WRAPPER) && isfile(GRADE_WRAPPER)
        include(FREEZE_WRAPPER)
        include(GRADE_WRAPPER)

        @test isdefined(Main, :StructuredFreezeLifecycle)
        @test isdefined(Main, :StructuredGradeLifecycle)
        @test isdefined(Main.StructuredFreezeLifecycle, :freeze_candidate)
        @test isdefined(Main.StructuredFreezeLifecycle, :derive_lifecycle)
        @test isdefined(Main.StructuredGradeLifecycle, :grade_frozen_candidate)
        @test isdefined(Main.StructuredGradeLifecycle, :recover_reserved_grade)
        @test isdefined(Main.StructuredGradeLifecycle, :parse_grade_summary)
    end
end

const FreezeLifecycle = Main.StructuredFreezeLifecycle
const GradeLifecycle = Main.StructuredGradeLifecycle
const StructuredEvidence = Main.StructuredFreezeLifecycle.StructuredEvidence
const FIXTURE_COUNTER = Ref(0)

file_sha256(path::String) = bytes2hex(sha256(read(path)))
relative_path(path::String) = replace(relpath(path, ROOT), '\\' => '/')

function write_fixture_file(path::String, contents)
    mkpath(dirname(path))
    write(path, contents)
    return path
end

function fake_grader_source()
    return raw"""
    #!/usr/bin/env julia
    using TOML

    function parse_options(arguments)
        values = Dict{String,String}()
        index = 1
        while index <= length(arguments)
            flag = arguments[index]
            index < length(arguments) || error("missing value for $flag")
            values[flag] = arguments[index + 1]
            index += 2
        end
        return values
    end

    options = parse_options(copy(ARGS))
    prediction = options["--prediction"]
    destination = options["--destination"]
    reservation = options["--reservation"]
    counter = options["--counter"]
    mode = options["--mode"]
    sleep_ms = parse(Int, options["--sleep-ms"])

    isfile(reservation) || exit(91)
    receipt = TOML.parsefile(reservation)
    get(receipt, "state", "") == "grade_reserved" || exit(92)
    get(receipt, "report_destination", "") == destination || exit(93)
    get(receipt, "prediction_path", "") == prediction || exit(94)

    open(counter, "a") do io
        println(io, "reservation_preceded_side_effect=true")
        flush(io)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 || exit(95)
    end
    env_keys = sort!(collect(keys(ENV)))
    write(joinpath(destination, "fake-environment.keys"), join(env_keys, "\n") * "\n")
    any(key -> startswith(key, "STMFIT_"), env_keys) && exit(96)
    sleep_ms > 0 && sleep(sleep_ms / 1000)
    mode == "error" && exit(7)

    grade_tsv = joinpath(destination, "grades", "structured_v2.tsv")
    write(grade_tsv, "file\tN_truth\tN_pred\nfixture.sxm\t6\t6\n")
    write(joinpath(destination, "report.md"), "# Generated fake report\n")
    write(joinpath(destination, "lobe_position_errors.tsv"),
        "profile\tlobe\tpossible\tclassified\tcorrect\twrong\tuncertain\tclassified_accuracy\thonest_correct_frac\n")

    summary = joinpath(destination, "summary.tsv")
    if mode == "malformed"
        write(summary, "PASS: fake grader says success\n")
        exit(0)
    end
    header = [
        "profile", "files", "lobes_possible", "lobes_classified", "coverage",
        "prediction_lobes", "missing_truth_lobes", "extra_prediction_lobes",
        "files_short_N", "files_extra_N", "physical_correct", "physical_errors_emitted",
        "physical_accuracy_classified", "honest_correct", "honest_uncertain",
        "honest_correct_frac", "sequence_exact", "oracle_correct",
        "oracle_accuracy_classified", "prediction_path", "grade_tsv", "note",
    ]
    row = [
        "structured_v2", "145", "870", "854", "98.2%", "854", "16", "0",
        "16", "0", "677", "177", "79.3%", "677", "193", "77.8%",
        "36", "700", "82.0%", prediction, grade_tsv, "generated fake grader",
    ]
    write(summary, join(header, '\t') * "\n" * join(row, '\t') * "\n")
    println("generated fake grade")
    exit(0)
    """
end

struct Fixture
    campaign::String
    campaign_id::String
    attempt_id::String
    prediction::String
    universe::String
    source::String
    freeze_manifest::String
    eligibility_receipt::String
    eligibility_result
    gates::Vector{Any}
    fake_grader::String
    invocation_counter::String
end

function artifact_specs(directory::String, stem::String)
    config = write_fixture_file(joinpath(directory, "$stem-config.toml"),
        "[model]\nname = \"$stem\"\n")
    source = write_fixture_file(joinpath(directory, "$stem-source.jl"), "$stem-source\n")
    input = write_fixture_file(joinpath(directory, "$stem-input.tsv"), "key\tvalue\na\t1\n")
    output = write_fixture_file(joinpath(directory, "$stem-output.tsv"), "key\tvalue\na\t2\n")
    expected = write_fixture_file(joinpath(directory, "$stem-expected.tsv"), "key\na\n")
    return StructuredEvidence.ArtifactSpec[
        StructuredEvidence.ArtifactSpec("config", relative_path(config)),
        StructuredEvidence.ArtifactSpec("source", relative_path(source)),
        StructuredEvidence.ArtifactSpec("input", relative_path(input)),
        StructuredEvidence.ArtifactSpec("output", relative_path(output)),
        StructuredEvidence.ArtifactSpec("expected_rows", relative_path(expected)),
    ]
end

function publish_gate(campaign::String, campaign_id::String, name::String, terminal::String)
    model_id = "gate-$name"
    artifacts = artifact_specs(joinpath(campaign, "gate-artifacts", name), name)
    directory = joinpath(campaign, "evidence", "gates", name)
    mkpath(directory)
    command = "julia --project=$ROOT test/fake-$name-gate.jl --campaign $campaign_id"
    result = StructuredEvidence.publish_evidence_bundle(
        ROOT,
        directory;
        campaign_id,
        model_id,
        artifacts,
        exact_command=command,
        generated_at="2026-08-08T00:00:00Z",
        terminal_status=terminal,
        reason_codes=[terminal == "PASS" ? "ok" : "gate_terminal"],
        reasons=["generated mandatory gate $name terminal $terminal"],
    )
    @test result.valid
    @test result.published
    return (name=name, campaign_id=campaign_id, model_id=model_id, result=result,
            command_identity=StructuredEvidence.command_identity_sha256(ROOT, command))
end

function freeze_manifest_text(
    campaign_id::String,
    prediction::String,
    universe::String,
    eligibility_command_identity::String,
    gates,
)
    io = IOBuffer()
    println(io, "schema = \"structured_freeze_manifest_v1\"")
    println(io, "schema_version = 1")
    println(io, "campaign_id = \"$campaign_id\"")
    println(io, "model_id = \"structured-v2\"")
    println(io, "prediction_path = \"$(relative_path(prediction))\"")
    println(io, "universe_path = \"$(relative_path(universe))\"")
    println(io, "eligibility_command_identity_sha256 = \"$eligibility_command_identity\"")
    for gate in sort(gates; by=item -> item.name)
        println(io, "\n[[mandatory_gates]]")
        println(io, "name = \"$(gate.name)\"")
        println(io, "campaign_id = \"$(gate.campaign_id)\"")
        println(io, "model_id = \"$(gate.model_id)\"")
        println(io, "receipt_path = \"$(relative_path(gate.result.receipt_path))\"")
        println(io, "receipt_sha256 = \"$(gate.result.receipt_sha256)\"")
        println(io, "command_identity_sha256 = \"$(gate.command_identity)\"")
    end
    return String(take!(io))
end

function build_fixture(;
    gate_statuses::Vector{String}=["PASS", "PASS"],
    primary_status::String="PASS",
)
    FIXTURE_COUNTER[] += 1
    campaign_id = "campaign-t5-$(FIXTURE_COUNTER[])"
    attempt_id = "attempt-1"
    campaign = mktempdir(T5_EVIDENCE; prefix="lifecycle-fixture-")
    prediction = write_fixture_file(joinpath(campaign, "predictions.tsv"),
        "file\tlobe\tpredicted\nfixture.sxm\t1\t0\n")
    universe = write_fixture_file(joinpath(campaign, "evidence", "universe.tsv"),
        "file\nfixture.sxm\n")
    source = write_fixture_file(joinpath(campaign, "evidence", "structured-source.jl"),
        "generated structured source\n")
    config = write_fixture_file(joinpath(campaign, "evidence", "model.toml"),
        "[model]\nname = \"structured-v2\"\n")
    expected = write_fixture_file(joinpath(campaign, "evidence", "expected.tsv"),
        "file\nfixture.sxm\n")
    gates = Any[
        publish_gate(campaign, campaign_id, "gate-$(index)", terminal)
        for (index, terminal) in enumerate(gate_statuses)
    ]
    eligibility_command =
        "julia --project=$ROOT test/fake-eligibility.jl --campaign $campaign_id"
    eligibility_command_identity =
        StructuredEvidence.command_identity_sha256(ROOT, eligibility_command)
    manifest = write_fixture_file(
        joinpath(campaign, "evidence", "freeze-manifest.toml"),
        freeze_manifest_text(
            campaign_id,
            prediction,
            universe,
            eligibility_command_identity,
            gates,
        ),
    )
    specs = StructuredEvidence.ArtifactSpec[
        StructuredEvidence.ArtifactSpec("config", relative_path(config)),
        StructuredEvidence.ArtifactSpec("source", relative_path(source)),
        StructuredEvidence.ArtifactSpec("input", relative_path(universe)),
        StructuredEvidence.ArtifactSpec("output", relative_path(prediction)),
        StructuredEvidence.ArtifactSpec("expected_rows", relative_path(expected)),
        StructuredEvidence.ArtifactSpec("freeze_manifest", relative_path(manifest)),
    ]
    append!(specs, [
        StructuredEvidence.ArtifactSpec("gate_receipt", relative_path(gate.result.receipt_path))
        for gate in gates
    ])
    eligibility_directory = joinpath(campaign, "evidence", "eligibility")
    mkpath(eligibility_directory)
    eligibility = StructuredEvidence.publish_evidence_bundle(
        ROOT,
        eligibility_directory;
        campaign_id,
        model_id="structured-v2",
        artifacts=specs,
        exact_command=eligibility_command,
        generated_at="2026-08-08T00:00:00Z",
        terminal_status=primary_status,
        reason_codes=[primary_status == "PASS" ? "ok" : "eligibility_terminal"],
        reasons=["generated eligibility terminal $primary_status"],
    )
    @test eligibility.valid
    @test eligibility.published
    fake_grader = write_fixture_file(joinpath(campaign, "fake_grader.jl"), fake_grader_source())
    counter = joinpath(campaign, "fake-invocations.log")
    return Fixture(
        campaign,
        campaign_id,
        attempt_id,
        prediction,
        universe,
        source,
        manifest,
        eligibility.receipt_path,
        eligibility,
        gates,
        fake_grader,
        counter,
    )
end

function freeze_request(fixture::Fixture)
    return FreezeLifecycle.FreezeRequest(
        fixture.campaign,
        fixture.campaign_id,
        fixture.attempt_id,
        "structured-v2",
        fixture.eligibility_receipt,
        fixture.eligibility_result.receipt_sha256,
        fixture.eligibility_result.command_identity_sha256,
    )
end

freeze_fixture(fixture::Fixture) = FreezeLifecycle.freeze_candidate(freeze_request(fixture))

function fake_launcher(fixture::Fixture; mode::String="success", sleep_ms::Int=0)
    return function(arguments, environment, cwd, stdout_path, stderr_path)
        destination = joinpath(fixture.campaign, "grade")
        expected = GradeLifecycle.command_arguments(
            FreezeLifecycle.production_context(), fixture.prediction, destination)
        arguments == expected ||
            throw(FreezeLifecycle.LifecycleError("fake_command_mismatch", "reserved command changed"))
        cwd == ROOT ||
            throw(FreezeLifecycle.LifecycleError("fake_cwd_mismatch", "launch cwd changed"))
        Set(keys(environment)) <= Set(FreezeLifecycle.ALLOWED_ENVIRONMENT_KEYS) ||
            throw(FreezeLifecycle.LifecycleError("fake_environment_mismatch", "forbidden environment key"))
        any(key -> startswith(key, "STMFIT_"), keys(environment)) &&
            throw(FreezeLifecycle.LifecycleError("fake_environment_mismatch", "STMFIT key leaked"))
        command = Cmd(Cmd(vcat(
            collect(Base.julia_cmd()),
            [
                "--project=$ROOT",
                fixture.fake_grader,
                "--prediction", fixture.prediction,
                "--destination", destination,
                "--reservation", joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"),
                "--counter", fixture.invocation_counter,
                "--mode", mode,
                "--sleep-ms", string(sleep_ms),
            ],
        )); dir=cwd)
        command = setenv(command, environment)
        open(stdout_path, "a") do stdout_io
            open(stderr_path, "a") do stderr_io
                process = run(pipeline(ignorestatus(command); stdout=stdout_io, stderr=stderr_io); wait=false)
                wait(process)
                return process.exitcode
            end
        end
    end
end

function invocation_count(fixture::Fixture)
    !isfile(fixture.invocation_counter) && return 0
    return length(split(chomp(read(fixture.invocation_counter, String)), '\n'; keepempty=false))
end

function run_capture(command::Cmd)
    stdout = IOBuffer()
    stderr = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout, stderr); wait=false)
    wait(process)
    return process.exitcode, String(take!(stdout)), String(take!(stderr))
end

function lifecycle_checker(fixture::Fixture, expectation::String)
    checker = joinpath(ROOT, "test", "check_unit_assignment_candidate_manifest.jl")
    return run_capture(`$(Base.julia_cmd()) --project=$ROOT $checker --config $CANDIDATE --model $MODEL --lifecycle-dir $(joinpath(fixture.campaign, "lifecycle")) $expectation`)
end

function assert_no_lifecycle_debris(fixture::Fixture)
    lifecycle = joinpath(fixture.campaign, "lifecycle")
    !isdir(lifecycle) && return
    debris = filter(readdir(lifecycle)) do name
        occursin("lock", name) || occursin("stage", name) || occursin("partial", name) ||
            occursin("tmp", name) || startswith(name, "jl_")
    end
    @test isempty(debris)
end

function cleanup_fixture(fixture::Fixture)
    ispath(fixture.campaign) && rm(fixture.campaign; recursive=true, force=true)
    @test !ispath(fixture.campaign)
end

@testset "generated evidence freezes exactly one immutable candidate" begin
    candidate_before = read(CANDIDATE)
    model_before = read(MODEL)
    fixture = build_fixture()
    try
        frozen = freeze_fixture(fixture)
        @test frozen.document["state"] == "frozen_once"
        @test frozen.document["campaign_id"] == fixture.campaign_id
        @test frozen.document["prediction_path"] == fixture.prediction
        @test frozen.document["prediction_sha256"] == file_sha256(fixture.prediction)
        @test frozen.document["source_bundle_sha256"] ==
              StructuredEvidence.parse_artifact_index(
                  read(fixture.eligibility_result.artifact_index_path)).source_sha256
        @test frozen.document["universe_sha256"] == file_sha256(fixture.universe)
        @test frozen.document["evidence_index_sha256"] ==
              file_sha256(fixture.eligibility_result.artifact_index_path)
        @test frozen.document["eligibility_receipt_sha256"] ==
              fixture.eligibility_result.receipt_sha256
        @test frozen.document["grader_source_sha256"] ==
              file_sha256(joinpath(ROOT, "test", "report_unit_assignment_benchmark.jl"))
        @test stat(frozen.path).mode & 0o222 == 0
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
        checker_exit, checker_out, checker_err = lifecycle_checker(fixture, "--expect-frozen-once")
        @test checker_exit == 0
        @test isempty(checker_err)
        @test occursin("lifecycle_state: frozen_once", checker_out)
        @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
        @test read(CANDIDATE) == candidate_before
        @test read(MODEL) == model_before
        @test FreezeLifecycle.serialize_freeze(frozen.document) ==
              FreezeLifecycle.serialize_freeze(frozen.document)
        assert_no_lifecycle_debris(fixture)
    finally
        cleanup_fixture(fixture)
    end
end

@testset "successful fake grade is reserved before one side effect" begin
    candidate_before = read(CANDIDATE)
    model_before = read(MODEL)
    fixture = build_fixture()
    try
        freeze_fixture(fixture)
        evidence_index_before = read(fixture.eligibility_result.artifact_index_path)
        eligibility_receipt_before = read(fixture.eligibility_receipt)
        previous = get(ENV, "STMFIT_FAKE_SECRET", nothing)
        ENV["STMFIT_FAKE_SECRET"] = "must-not-leak"
        outcome = try
            GradeLifecycle.grade_frozen_candidate(
                fixture.campaign;
                launcher=fake_launcher(fixture),
            )
        finally
            if previous === nothing
                delete!(ENV, "STMFIT_FAKE_SECRET")
            else
                ENV["STMFIT_FAKE_SECRET"] = previous
            end
        end
        @test outcome.outcome == "success"
        @test outcome.process_exit_code == 0
        @test outcome.invocation_count == 1
        @test invocation_count(fixture) == 1
        @test read(fixture.invocation_counter, String) ==
              "reservation_preceded_side_effect=true\n"
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
        reservation_path = joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml")
        reservation = TOML.parsefile(reservation_path)
        expected_command =
            "julia --project=$ROOT $(joinpath(ROOT, "test", "report_unit_assignment_benchmark.jl")) " *
            "--full145-own-n --profile structured_v2=$(fixture.prediction) " *
            "--outdir $(joinpath(fixture.campaign, "grade"))"
        @test reservation["repository_root"] == ROOT
        @test reservation["project_path"] == joinpath(ROOT, "Project.toml")
        @test reservation["project_sha256"] == file_sha256(joinpath(ROOT, "Project.toml"))
        @test reservation["grade_wrapper_path"] == GRADE_WRAPPER
        @test reservation["grade_wrapper_sha256"] == file_sha256(GRADE_WRAPPER)
        @test reservation["report_script_path"] ==
              joinpath(ROOT, "test", "report_unit_assignment_benchmark.jl")
        @test reservation["report_script_sha256"] ==
              file_sha256(joinpath(ROOT, "test", "report_unit_assignment_benchmark.jl"))
        @test reservation["prediction_path"] == fixture.prediction
        @test reservation["report_destination"] == joinpath(fixture.campaign, "grade")
        @test reservation["exact_command"] == expected_command
        @test reservation["working_directory"] == ROOT
        @test reservation["environment_keys"] ==
              [key for key in FreezeLifecycle.ALLOWED_ENVIRONMENT_KEYS if haskey(ENV, key)]
        @test stat(reservation_path).mode & 0o222 == 0
        checker_exit, checker_out, checker_err = lifecycle_checker(fixture, "--expect-graded")
        @test checker_exit == 0
        @test isempty(checker_err)
        @test occursin("lifecycle_state: graded", checker_out)
        final = TOML.parsefile(outcome.final_path)
        @test final["state"] == "graded"
        @test final["outcome"] == "success"
        @test final["summary_tsv_sha256"] ==
              file_sha256(joinpath(fixture.campaign, "grade", "summary.tsv"))
        @test final["report_md_sha256"] ==
              file_sha256(joinpath(fixture.campaign, "grade", "report.md"))
        @test final["lobe_position_errors_tsv_sha256"] ==
              file_sha256(joinpath(fixture.campaign, "grade", "lobe_position_errors.tsv"))
        @test final["grade_tsv_sha256"] ==
              file_sha256(joinpath(fixture.campaign, "grade", "grades", "structured_v2.tsv"))
        @test final["stdout_sha256"] ==
              file_sha256(joinpath(fixture.campaign, "grade", "stdout.log"))
        @test final["stderr_sha256"] ==
              file_sha256(joinpath(fixture.campaign, "grade", "stderr.log"))
        @test stat(outcome.final_path).mode & 0o222 == 0
        metrics = GradeLifecycle.parse_grade_summary(
            joinpath(fixture.campaign, "grade", "summary.tsv");
            frozen_prediction=fixture.prediction,
            reserved_grade_directory=joinpath(fixture.campaign, "grade"),
        )
        @test length(metrics.integers) == 14
        @test metrics.integers["physical_correct"] == 677
        @test metrics.fractions["coverage"] == 0.982
        @test metrics.fractions["physical_accuracy_classified"] ≈ 0.793 atol=1.0e-15
        @test metrics.fractions["honest_correct_frac"] ≈ 0.778 atol=1.0e-15
        @test metrics.fractions["oracle_accuracy_classified"] == 0.82
        environment_keys = split(chomp(read(
            joinpath(fixture.campaign, "grade", "fake-environment.keys"), String)), '\n')
        @test !any(key -> startswith(key, "STMFIT_"), environment_keys)
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        @test invocation_count(fixture) == 1
        @test read(CANDIDATE) == candidate_before
        @test read(MODEL) == model_before
        @test read(fixture.eligibility_result.artifact_index_path) == evidence_index_before
        @test read(fixture.eligibility_receipt) == eligibility_receipt_before
        assert_no_lifecycle_debris(fixture)
    finally
        cleanup_fixture(fixture)
    end
end

function valid_summary_row(prediction::String, grade_tsv::String)
    return Dict(
        "profile" => "structured_v2",
        "files" => "145",
        "lobes_possible" => "870",
        "lobes_classified" => "854",
        "coverage" => "98.2%",
        "prediction_lobes" => "854",
        "missing_truth_lobes" => "16",
        "extra_prediction_lobes" => "0",
        "files_short_N" => "16",
        "files_extra_N" => "0",
        "physical_correct" => "677",
        "physical_errors_emitted" => "177",
        "physical_accuracy_classified" => "79.3%",
        "honest_correct" => "677",
        "honest_uncertain" => "193",
        "honest_correct_frac" => "77.8%",
        "sequence_exact" => "36",
        "oracle_correct" => "700",
        "oracle_accuracy_classified" => "82.0%",
        "prediction_path" => prediction,
        "grade_tsv" => grade_tsv,
        "note" => "generated parser fixture",
    )
end

function write_summary_fixture(
    path::String,
    row::Dict{String,String};
    header::Vector{String}=copy(GradeLifecycle.GRADE_HEADER),
    rows::Int=1,
    trailing_lf::Bool=true,
)
    io = IOBuffer()
    println(io, join(header, '\t'))
    for _ in 1:rows
        println(io, join([get(row, field, "") for field in header], '\t'))
    end
    bytes = take!(io)
    !trailing_lf && !isempty(bytes) && pop!(bytes)
    write(path, bytes)
    return path
end

function with_parser_fixture(function_body::Function)
    campaign = mktempdir(T5_EVIDENCE; prefix="parser-fixture-")
    try
        prediction = write_fixture_file(joinpath(campaign, "predictions.tsv"),
            "file\tlobe\tpredicted\nfixture.sxm\t1\t0\n")
        destination = joinpath(campaign, "grade")
        mkpath(joinpath(destination, "grades"))
        grade_tsv = write_fixture_file(joinpath(destination, "grades", "structured_v2.tsv"),
            "file\tN_truth\tN_pred\nfixture.sxm\t6\t6\n")
        summary = joinpath(destination, "summary.tsv")
        row = valid_summary_row(prediction, grade_tsv)
        write_summary_fixture(summary, row)
        return function_body(campaign, prediction, destination, grade_tsv, summary, row)
    finally
        rm(campaign; recursive=true, force=true)
        @test !ispath(campaign)
    end
end

function parser_call(summary::String, prediction::String, destination::String)
    return GradeLifecycle.parse_grade_summary(
        summary;
        frozen_prediction=prediction,
        reserved_grade_directory=destination,
    )
end

@testset "strict exact 22-column grade parser adversarial matrix" begin
    with_parser_fixture() do campaign, prediction, destination, grade_tsv, summary, row
        valid_bytes = read(summary)
        metrics = parser_call(summary, prediction, destination)
        @test length(metrics.integers) == 14
        @test length(metrics.fractions) == 4

        malformed_headers = [
            vcat(GradeLifecycle.GRADE_HEADER, ["extra"]),
            GradeLifecycle.GRADE_HEADER[1:end-1],
            vcat(["profile", "profile"], GradeLifecycle.GRADE_HEADER[3:end]),
            vcat(["PROFILE"], GradeLifecycle.GRADE_HEADER[2:end]),
        ]
        for header in malformed_headers
            write_summary_fixture(summary, row; header)
            @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        end

        write_summary_fixture(summary, row; rows=2)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        write(summary, join(GradeLifecycle.GRADE_HEADER, '\t') * "\n")
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        write_summary_fixture(summary, row; trailing_lf=false)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)

        for field in GradeLifecycle.INTEGER_FIELDS
            changed = copy(row)
            changed[field] = "-1"
            write_summary_fixture(summary, changed)
            @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
            changed[field] = "01"
            write_summary_fixture(summary, changed)
            @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        end

        for field in GradeLifecycle.PERCENT_FIELDS
            for invalid in ("NA", "79%", "079.3%", "79.30%", "-1.0%", "101.0%")
                changed = copy(row)
                changed[field] = invalid
                write_summary_fixture(summary, changed)
                @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
            end
        end

        for (field, wrong) in (
            "coverage" => "98.1%",
            "physical_accuracy_classified" => "79.2%",
            "honest_correct_frac" => "77.9%",
            "oracle_accuracy_classified" => "81.9%",
        )
            changed = copy(row)
            changed[field] = wrong
            write_summary_fixture(summary, changed)
            @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        end

        for (field, wrong) in (
            "files" => "144",
            "lobes_possible" => "869",
            "physical_errors_emitted" => "176",
            "honest_correct" => "676",
            "honest_uncertain" => "192",
        )
            changed = copy(row)
            changed[field] = wrong
            write_summary_fixture(summary, changed)
            @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        end

        changed = copy(row)
        changed["profile"] = "PASS"
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)

        other_prediction = write_fixture_file(joinpath(campaign, "other-predictions.tsv"),
            "file\tlobe\tpredicted\nfixture.sxm\t1\t1\n")
        changed = copy(row)
        changed["prediction_path"] = other_prediction
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)

        changed = copy(row)
        changed["prediction_path"] = joinpath(campaign, "missing-predictions.tsv")
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)

        changed = copy(row)
        changed["grade_tsv"] = prediction
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)

        changed = copy(row)
        changed["grade_tsv"] = string(joinpath(destination, "grades"), "/../grades/structured_v2.tsv")
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)

        link = joinpath(destination, "grades", "linked.tsv")
        symlink(grade_tsv, link)
        changed = copy(row)
        changed["grade_tsv"] = link
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        rm(link)

        fifo = joinpath(destination, "grades", "special.tsv")
        run(`mkfifo $fifo`)
        changed["grade_tsv"] = fifo
        write_summary_fixture(summary, changed)
        @test_throws FreezeLifecycle.LifecycleError parser_call(summary, prediction, destination)
        rm(fifo)

        write(summary, valid_bytes)
        @test parser_call(summary, prediction, destination).integers["files"] == 145
    end
end

function with_fixture(function_body::Function; kwargs...)
    fixture = build_fixture(; kwargs...)
    try
        return function_body(fixture)
    finally
        cleanup_fixture(fixture)
    end
end

@testset "freeze rejects absent, ineligible, stale, and misleading receipts" begin
    for primary_status in ("FAIL", "BLOCKED", "SKIPPED", "FOLLOWUP_REQUIRED")
        with_fixture(; primary_status) do fixture
            @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "freeze_receipt.toml"))
        end
    end

    for gate_status in ("FAIL", "BLOCKED", "SKIPPED", "FOLLOWUP_REQUIRED")
        with_fixture(; gate_statuses=[gate_status]) do fixture
            @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "freeze_receipt.toml"))
        end
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        rm(fixture.eligibility_receipt)
        @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        gate_path = fixture.gates[1].result.receipt_path
        rm(gate_path)
        @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        chmod(fixture.eligibility_receipt, 0o600)
        write(fixture.eligibility_receipt, "PASS: all mandatory gates passed\n")
        @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
    end

    with_fixture(; gate_statuses=["PASS", "PASS"]) do fixture
        first_gate = fixture.gates[1].result.receipt_path
        chmod(first_gate, 0o600)
        write(first_gate, read(fixture.gates[2].result.receipt_path))
        @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        index_path = fixture.eligibility_result.artifact_index_path
        bytes = read(index_path)
        bytes[end - 1] = xor(bytes[end - 1], 0x01)
        chmod(index_path, 0o600)
        write(index_path, bytes)
        @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        request = freeze_request(fixture)
        wrong_receipt = FreezeLifecycle.FreezeRequest(
            request.campaign_dir, request.campaign_id, request.attempt_id, request.model_id,
            request.eligibility_receipt, repeat('0', 64),
            request.expected_command_identity_sha256,
        )
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.freeze_candidate(wrong_receipt)
        wrong_command = FreezeLifecycle.FreezeRequest(
            request.campaign_dir, request.campaign_id, request.attempt_id, request.model_id,
            request.eligibility_receipt, request.expected_receipt_sha256, repeat('0', 64),
        )
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.freeze_candidate(wrong_command)
        wrong_model = FreezeLifecycle.FreezeRequest(
            request.campaign_dir, request.campaign_id, request.attempt_id, "structured-v2-other",
            request.eligibility_receipt, request.expected_receipt_sha256,
            request.expected_command_identity_sha256,
        )
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.freeze_candidate(wrong_model)
        @test !ispath(joinpath(fixture.campaign, "lifecycle", "freeze_receipt.toml"))
    end
end

function replaced_toml_field(text::String, key::String, replacement::String)
    pattern = Regex("(?m)^" * key * " = \"[^\"]*\"\$")
    length(findall(pattern, text)) == 1 || error("expected one $key field")
    return replace(text, pattern => "$key = \"$replacement\""; count=1)
end

function assert_no_transient_debris(fixture::Fixture)
    debris = String[]
    for (directory, directories, files) in walkdir(fixture.campaign)
        for name in vcat(directories, files)
            lowered = lowercase(name)
            if occursin(".lock", lowered) || occursin(".stage", lowered) ||
               occursin(".partial", lowered) || endswith(lowered, ".tmp") ||
               startswith(name, "jl_")
                push!(debris, joinpath(directory, name))
            end
        end
    end
    @test isempty(debris)
end

@testset "root, project, prediction, source, grader, and path identities fail pre-reservation" begin
    candidate_before = read(CANDIDATE)
    model_before = read(MODEL)
    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        context = FreezeLifecycle.production_context()
        fake_project = write_fixture_file(joinpath(fixture.campaign, "OtherProject.toml"),
            "name = \"other\"\n")
        fake_source = write_fixture_file(joinpath(fixture.campaign, "other_report.jl"),
            "println(\"not grader\")\n")
        fake_wrapper = write_fixture_file(joinpath(fixture.campaign, "other_wrapper.jl"),
            "println(\"not wrapper\")\n")
        changed_contexts = [
            FreezeLifecycle.LifecycleContext(
                dirname(ROOT), context.candidate_path, context.model_path, context.project_path,
                context.grade_wrapper_path, context.report_script_path),
            FreezeLifecycle.LifecycleContext(
                ROOT, context.candidate_path, context.model_path, fake_project,
                context.grade_wrapper_path, context.report_script_path),
            FreezeLifecycle.LifecycleContext(
                ROOT, context.candidate_path, context.model_path, context.project_path,
                fake_wrapper, context.report_script_path),
            FreezeLifecycle.LifecycleContext(
                ROOT, context.candidate_path, context.model_path, context.project_path,
                context.grade_wrapper_path, fake_source),
        ]
        for changed_context in changed_contexts
            @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                fixture.campaign;
                context=changed_context,
                launcher=fake_launcher(fixture),
            )
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))
        end

        prediction_bytes = read(fixture.prediction)
        write(fixture.prediction, vcat(prediction_bytes, UInt8('X')))
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))
        write(fixture.prediction, prediction_bytes)
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"

        source_bytes = read(fixture.source)
        write(fixture.source, vcat(source_bytes, UInt8('X')))
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))
        write(fixture.source, source_bytes)
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"

        grade_path = joinpath(fixture.campaign, "grade")
        target = joinpath(fixture.campaign, "grade-target")
        mkdir(target)
        symlink(target, grade_path)
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        rm(grade_path)
        rm(target)

        write(grade_path, "regular destination\n")
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        rm(grade_path)

        mkdir(grade_path)
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        rm(grade_path)

        run(`mkfifo $grade_path`)
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        rm(grade_path)
        @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))

        lock_path = joinpath(
            fixture.campaign, "lifecycle", ".grade-reservation.atomic-create.lock")
        write(lock_path, "foreign lock\n")
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        @test isfile(lock_path)
        rm(lock_path)

        linked_campaign = joinpath(T5_EVIDENCE, "linked-campaign-$(fixture.campaign_id)")
        symlink(fixture.campaign, linked_campaign)
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.derive_lifecycle(linked_campaign)
        rm(linked_campaign)

        outside = mktempdir(; prefix="outside-campaign-")
        try
            @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.derive_lifecycle(outside)
        finally
            rm(outside; recursive=true, force=true)
        end
        @test read(CANDIDATE) == candidate_before
        @test read(MODEL) == model_before
        assert_no_transient_debris(fixture)
    end
end

@testset "receipt substitution and current artifact staleness invalidate lifecycle" begin
    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        freeze_path = joinpath(fixture.campaign, "lifecycle", "freeze_receipt.toml")
        freeze_bytes = read(freeze_path)
        freeze_text = String(copy(freeze_bytes))
        mutated = replaced_toml_field(freeze_text, "prediction_sha256", repeat('0', 64))
        chmod(freeze_path, 0o600)
        write(freeze_path, mutated)
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.derive_lifecycle(fixture.campaign)
        write(freeze_path, freeze_bytes)
        chmod(freeze_path, 0o444)
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"

        misleading = joinpath(fixture.campaign, "lifecycle", "misleading.toml")
        write(misleading, "state = \"PASS\"\n")
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.derive_lifecycle(fixture.campaign)
        rm(misleading)

        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign;
            launcher=fake_launcher(fixture),
            failpoint=:after_reservation,
            leave_reserved=true,
        )
        reservation_path = joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml")
        reservation_bytes = read(reservation_path)
        reservation_text = String(copy(reservation_bytes))
        context = FreezeLifecycle.production_context()
        other_prediction = write_fixture_file(joinpath(fixture.campaign, "other.tsv"), "other\n")
        other_project = write_fixture_file(joinpath(fixture.campaign, "OtherProject.toml"), "other\n")
        other_source = write_fixture_file(joinpath(fixture.campaign, "other-source.jl"), "other\n")
        substitutions = [
            "repository_root" => fixture.campaign,
            "project_path" => other_project,
            "project_sha256" => file_sha256(other_project),
            "prediction_path" => other_prediction,
            "prediction_sha256" => file_sha256(other_prediction),
            "grade_wrapper_path" => other_source,
            "grade_wrapper_sha256" => file_sha256(other_source),
            "report_script_path" => other_source,
            "report_script_sha256" => file_sha256(other_source),
            "exact_command" => "true",
        ]
        for (field, value) in substitutions
            chmod(reservation_path, 0o600)
            write(reservation_path, replaced_toml_field(reservation_text, field, value))
            @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.derive_lifecycle(fixture.campaign)
            write(reservation_path, reservation_bytes)
            chmod(reservation_path, 0o444)
        end
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "grade_reserved"
        recovered = GradeLifecycle.recover_reserved_grade(fixture.campaign)
        @test recovered.outcome == "error"
        @test recovered.invocation_count == 0
        checker_exit, _, checker_err = lifecycle_checker(fixture, "--expect-graded")
        @test checker_exit == 0
        @test isempty(checker_err)
        assert_no_transient_debris(fixture)
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        outcome = GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture))
        @test outcome.outcome == "success"
        grade_tsv = joinpath(fixture.campaign, "grade", "grades", "structured_v2.tsv")
        original = read(grade_tsv)
        write(grade_tsv, vcat(original, UInt8('X')))
        @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.derive_lifecycle(fixture.campaign)
        write(grade_tsv, original)
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
        assert_no_transient_debris(fixture)
    end
end

@testset "pre-reservation crashes leave no attempt and permit one retry" begin
    for phase in (:before_lock, :after_lock, :before_reservation)
        @testset "$phase" begin
            with_fixture(; gate_statuses=["PASS"]) do fixture
                freeze_fixture(fixture)
                @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign;
                    launcher=fake_launcher(fixture),
                    failpoint=phase,
                )
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
                @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))
                @test !ispath(joinpath(fixture.campaign, "grade"))
                @test invocation_count(fixture) == 0
                assert_no_transient_debris(fixture)

                retry = GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
                @test retry.outcome == "success"
                @test retry.invocation_count == 1
                @test invocation_count(fixture) == 1
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
                assert_no_transient_debris(fixture)
            end
        end
    end
end

@testset "every post-reservation crash consumes attempt and recovers only ERROR" begin
    phases = [
        :after_reservation,
        :after_grade_directory,
        :before_launch_intent,
        :after_launch_intent,
        :after_process,
        :after_artifact_validation,
        :before_final,
    ]
    for phase in phases
        @testset "$phase" begin
            with_fixture(; gate_statuses=["PASS"]) do fixture
                freeze_fixture(fixture)
                launcher = fake_launcher(fixture)
                @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign;
                    launcher,
                    failpoint=phase,
                    leave_reserved=true,
                )
                lifecycle = FreezeLifecycle.derive_lifecycle(fixture.campaign)
                @test lifecycle.state == "grade_reserved"
                count_before_retry = invocation_count(fixture)
                @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher)
                @test invocation_count(fixture) == count_before_retry
                recovered = GradeLifecycle.recover_reserved_grade(fixture.campaign)
                @test recovered.outcome == "error"
                @test recovered.error_reason != ""
                expected_invocations = phase in (
                    :after_reservation, :after_grade_directory, :before_launch_intent) ? 0 : 1
                @test recovered.invocation_count == expected_invocations
                @test recovered.process_exit_code == (expected_invocations == 0 ? -1 :
                    (phase in (:after_process, :after_artifact_validation, :before_final) ? 0 : 255))
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
                @test TOML.parsefile(recovered.final_path)["outcome"] == "error"
                @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher)
                @test invocation_count(fixture) == count_before_retry
                checker_exit, _, checker_err = lifecycle_checker(fixture, "--expect-graded")
                @test checker_exit == 0
                @test isempty(checker_err)
                assert_no_transient_debris(fixture)
            end
        end
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        committed = GradeLifecycle.grade_frozen_candidate(
            fixture.campaign;
            launcher=fake_launcher(fixture),
            failpoint=:after_final,
            leave_reserved=true,
        )
        @test committed.outcome == "success"
        @test committed.invocation_count == 1
        @test invocation_count(fixture) == 1
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.recover_reserved_grade(
            fixture.campaign)
        assert_no_transient_debris(fixture)
    end
end

@testset "process and artifact errors finalize once without relaunch" begin
    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        outcome = GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture; mode="error"))
        @test outcome.outcome == "error"
        @test outcome.error_reason == "grader_process_failed"
        @test outcome.process_exit_code == 7
        @test outcome.invocation_count == 1
        @test invocation_count(fixture) == 1
        @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture; mode="success"))
        @test invocation_count(fixture) == 1
        checker_exit, _, checker_err = lifecycle_checker(fixture, "--expect-graded")
        @test checker_exit == 0
        @test isempty(checker_err)
        assert_no_transient_debris(fixture)
    end

    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        outcome = GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture; mode="malformed"))
        @test outcome.outcome == "error"
        @test startswith(outcome.error_reason, "invalid_grade_artifacts_")
        @test outcome.process_exit_code == 0
        @test outcome.invocation_count == 1
        @test invocation_count(fixture) == 1
        @test read(joinpath(fixture.campaign, "grade", "summary.tsv"), String) ==
              "PASS: fake grader says success\n"
        @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
            fixture.campaign; launcher=fake_launcher(fixture; mode="success"))
        @test invocation_count(fixture) == 1
        assert_no_transient_debris(fixture)
    end
end

function race_worker_source()
    return raw"""
    #!/usr/bin/env julia
    grade_wrapper, campaign, fake_grader, prediction, counter = ARGS
    include(grade_wrapper)
    using .StructuredGradeLifecycle

    function launcher(arguments, environment, cwd, stdout_path, stderr_path)
        destination = joinpath(campaign, "grade")
        command = Cmd(Cmd(vcat(
            collect(Base.julia_cmd()),
            [
                "--project=$cwd",
                fake_grader,
                "--prediction", prediction,
                "--destination", destination,
                "--reservation", joinpath(campaign, "lifecycle", "grade_reservation.toml"),
                "--counter", counter,
                "--mode", "success",
                "--sleep-ms", "1000",
            ],
        )); dir=cwd)
        command = setenv(command, environment)
        open(stdout_path, "a") do stdout_io
            open(stderr_path, "a") do stderr_io
                process = run(pipeline(ignorestatus(command); stdout=stdout_io, stderr=stderr_io); wait=false)
                wait(process)
                return process.exitcode
            end
        end
    end

    try
        outcome = StructuredGradeLifecycle.grade_frozen_candidate(campaign; launcher)
        println("outcome=$(outcome.outcome)")
        exit(outcome.outcome == "success" ? 0 : 1)
    catch error
        println(stderr, "blocked: ", sprint(showerror, error))
        exit(2)
    end
    """
end

@testset "real process race permits exactly one durable reservation and invocation" begin
    candidate_before = read(CANDIDATE)
    model_before = read(MODEL)
    with_fixture(; gate_statuses=["PASS"]) do fixture
        freeze_fixture(fixture)
        worker = write_fixture_file(joinpath(fixture.campaign, "race-worker.jl"),
            race_worker_source())
        stdout_paths = [joinpath(fixture.campaign, "worker-$index.stdout") for index in 1:2]
        stderr_paths = [joinpath(fixture.campaign, "worker-$index.stderr") for index in 1:2]
        stdout_ios = [open(path, "w") for path in stdout_paths]
        stderr_ios = [open(path, "w") for path in stderr_paths]
        processes = Any[]
        try
            for index in 1:2
                command = `$(Base.julia_cmd()) --project=$ROOT $worker $GRADE_WRAPPER $(fixture.campaign) $(fixture.fake_grader) $(fixture.prediction) $(fixture.invocation_counter)`
                push!(processes, run(pipeline(ignorestatus(command);
                    stdout=stdout_ios[index], stderr=stderr_ios[index]); wait=false))
            end
            foreach(wait, processes)
        finally
            foreach(io -> isopen(io) && close(io), stdout_ios)
            foreach(io -> isopen(io) && close(io), stderr_ios)
        end
        @test sort([process.exitcode for process in processes]) == [0, 2]
        @test invocation_count(fixture) == 1
        @test read(fixture.invocation_counter, String) ==
              "reservation_preceded_side_effect=true\n"
        lifecycle = FreezeLifecycle.derive_lifecycle(fixture.campaign)
        @test lifecycle.state == "graded"
        @test lifecycle.final.document["outcome"] == "success"
        @test lifecycle.final.document["invocation_count"] == 1
        @test count(path -> occursin("outcome=success", read(path, String)), stdout_paths) == 1
        @test count(path -> occursin("blocked:", read(path, String)), stderr_paths) == 1
        checker_exit, checker_out, checker_err = lifecycle_checker(fixture, "--expect-graded")
        @test checker_exit == 0
        @test isempty(checker_err)
        @test occursin("lifecycle_state: graded", checker_out)
        @test read(CANDIDATE) == candidate_before
        @test read(MODEL) == model_before
        assert_no_transient_debris(fixture)
    end
end

function capture_durability_fault(
    function_body::Function,
    site::Symbol;
    occurrence::Int=1,
    before_throw::Function=(_path::String) -> nothing,
)
    hits = Ref(0)
    hook = function(actual_site::Symbol, path::String)
        actual_site == site || return nothing
        hits[] += 1
        hits[] == occurrence || return nothing
        before_throw(path)
        throw(FreezeLifecycle.LifecycleError(
            "injected_durability_fault",
            "injected durability failure at $site for $path",
        ))
    end
    value = nothing
    captured = nothing
    try
        value = FreezeLifecycle.with_durability_fault_hook(hook) do
            function_body()
        end
    catch error
        captured = error
    end
    return (value=value, error=captured, hits=hits[])
end

function assert_durability_warning(
    result,
    expected_path::String;
    code::String="publication_durability_uncertain",
)
    @test result.durability == "uncertain"
    @test length(result.durability_warnings) == 1
    warning = only(result.durability_warnings)
    @test warning.code == code
    @test warning.path == expected_path
    warning_io = IOBuffer()
    FreezeLifecycle.write_durability_warnings(warning_io, result.durability_warnings)
    warning_text = String(take!(warning_io))
    @test occursin("WARNING [$code]", warning_text)
    @test occursin(expected_path, warning_text)
    @test !occursin("BLOCKED", warning_text)
end

@testset "primitive durability faults reconcile exact authority and retryability" begin
    required_api = (
        :DurabilityWarning,
        :PublicationOutcome,
        :with_durability_fault_hook,
        :write_durability_warnings,
    )
    durability_api_available = all(name -> isdefined(FreezeLifecycle, name), required_api)
    @test durability_api_available

    if durability_api_available
        lock_names = (
            ".freeze.atomic-create.lock",
            ".grade-reservation.atomic-create.lock",
            ".grade-final.atomic-create.lock",
        )
        for site in (:lock_file_fsync, :lock_parent_fsync), lock_name in lock_names
            parent = mktempdir(T5_EVIDENCE; prefix="lock-fault-")
            try
                lock_path = joinpath(parent, lock_name)
                probe = capture_durability_fault(site) do
                    FreezeLifecycle.acquire_create_lock(parent, lock_name)
                end
                @test probe.error isa FreezeLifecycle.LifecycleError
                if probe.error isa FreezeLifecycle.LifecycleError
                    @test probe.error.code == "injected_durability_fault"
                end
                @test probe.hits == 1
                @test !ispath(lock_path)
                retry_lock = FreezeLifecycle.acquire_create_lock(parent, lock_name)
                @test isfile(lock_path)
                FreezeLifecycle.release_create_lock(retry_lock)
                @test !ispath(lock_path)
                @test isempty(readdir(parent))
            finally
                rm(parent; recursive=true, force=true)
            end
        end

        parent = mktempdir(T5_EVIDENCE; prefix="lock-cleanup-fault-")
        try
            lock_path = joinpath(parent, ".freeze.atomic-create.lock")
            hook = function(site::Symbol, _path::String)
                site in (:lock_file_fsync, :lock_cleanup_parent_fsync) || return nothing
                throw(FreezeLifecycle.LifecycleError(
                    "injected_durability_fault", "injected lock cleanup uncertainty"))
            end
            captured = try
                FreezeLifecycle.with_durability_fault_hook(hook) do
                    FreezeLifecycle.acquire_create_lock(parent, basename(lock_path))
                end
                nothing
            catch error
                error
            end
            @test captured isa FreezeLifecycle.LifecycleError
            if captured isa FreezeLifecycle.LifecycleError
                @test captured.code == "stale_lock_cleanup_unresolved"
            end
            @test !ispath(lock_path)
        finally
            rm(parent; recursive=true, force=true)
        end

        for lock_name in lock_names
            parent = mktempdir(T5_EVIDENCE; prefix="lock-release-fault-")
            try
                lock_path = joinpath(parent, lock_name)
                lock = FreezeLifecycle.acquire_create_lock(parent, lock_name)
                release = capture_durability_fault(:lock_release_parent_fsync) do
                    FreezeLifecycle.release_create_lock(lock)
                end
                @test release.error === nothing
                @test release.hits == 1
                @test !ispath(lock_path)
                @test length(release.value) == 1
                warning = only(release.value)
                @test warning.code == "lock_removal_durability_uncertain"
                @test warning.path == lock_path
                retry_lock = FreezeLifecycle.acquire_create_lock(parent, lock_name)
                FreezeLifecycle.release_create_lock(retry_lock)
                @test isempty(readdir(parent))
            finally
                rm(parent; recursive=true, force=true)
            end
        end

        parent = mktempdir(T5_EVIDENCE; prefix="lock-ownership-fault-")
        try
            lock_path = joinpath(parent, ".freeze.atomic-create.lock")
            replace_owned_lock = function(created_path::String)
                rm(created_path)
                write(created_path, "competing creator\n")
            end
            replaced = capture_durability_fault(
                :lock_file_fsync; before_throw=replace_owned_lock) do
                FreezeLifecycle.acquire_create_lock(parent, basename(lock_path))
            end
            @test replaced.error isa FreezeLifecycle.LifecycleError
            if replaced.error isa FreezeLifecycle.LifecycleError
                @test replaced.error.code == "stale_lock_cleanup_unresolved"
            end
            @test read(lock_path, String) == "competing creator\n"
            rm(lock_path)
        finally
            rm(parent; recursive=true, force=true)
        end

        receipt_names = ("freeze_receipt.toml", "grade_reservation.toml", "grade_final.toml")
        for receipt_name in receipt_names
            parent = mktempdir(T5_EVIDENCE; prefix="receipt-fault-")
            try
                bytes = collect(codeunits("kind = \"$receipt_name\"\n"))
                preinstall_path = joinpath(parent, "pre-$receipt_name")
                preinstall = capture_durability_fault(:receipt_preinstall_fsync) do
                    FreezeLifecycle.atomic_create_file(preinstall_path, bytes)
                end
                @test preinstall.error isa FreezeLifecycle.LifecycleError
                if preinstall.error isa FreezeLifecycle.LifecycleError
                    @test preinstall.error.code == "injected_durability_fault"
                end
                @test !ispath(preinstall_path)
                retry = FreezeLifecycle.atomic_create_file(preinstall_path, bytes)
                @test retry.durability == "confirmed"
                @test isempty(retry.durability_warnings)
                @test read(preinstall_path) == bytes

                postinstall_path = joinpath(parent, "post-$receipt_name")
                postinstall = capture_durability_fault(:receipt_postinstall_parent_fsync) do
                    FreezeLifecycle.atomic_create_file(postinstall_path, bytes)
                end
                @test postinstall.error === nothing
                @test postinstall.hits == 1
                @test read(postinstall_path) == bytes
                @test stat(postinstall_path).mode & 0o222 == 0
                assert_durability_warning(postinstall.value, postinstall_path)
                @test_throws FreezeLifecycle.LifecycleError FreezeLifecycle.atomic_create_file(
                    postinstall_path, bytes)
                @test read(postinstall_path) == bytes
                @test isempty(filter(name -> startswith(name, "jl_"), readdir(parent)))
            finally
                rm(parent; recursive=true, force=true)
            end
        end

        parent = mktempdir(T5_EVIDENCE; prefix="receipt-authority-fault-")
        try
            path = joinpath(parent, "freeze_receipt.toml")
            intended = collect(codeunits("state = \"frozen_once\"\n"))
            tamper = function(installed_path::String)
                chmod(installed_path, 0o600)
                write(installed_path, "state = \"tampered\"\n")
                chmod(installed_path, 0o444)
            end
            unresolved = capture_durability_fault(
                :receipt_postinstall_parent_fsync; before_throw=tamper) do
                FreezeLifecycle.atomic_create_file(path, intended)
            end
            @test unresolved.error isa FreezeLifecycle.LifecycleError
            if unresolved.error isa FreezeLifecycle.LifecycleError
                @test unresolved.error.code == "publication_authority_unresolved"
            end
            @test isfile(path)
            @test read(path) != intended
        finally
            rm(parent; recursive=true, force=true)
        end
    end
end

@testset "freeze and grade durability faults preserve one authoritative lifecycle" begin
    durability_api_available = isdefined(FreezeLifecycle, :with_durability_fault_hook)
    @test durability_api_available

    if durability_api_available
        for site in (:lock_file_fsync, :lock_parent_fsync)
            with_fixture(; gate_statuses=["PASS"]) do fixture
                failed_lock = capture_durability_fault(site) do
                    freeze_fixture(fixture)
                end
                @test failed_lock.error isa FreezeLifecycle.LifecycleError
                if failed_lock.error isa FreezeLifecycle.LifecycleError
                    @test failed_lock.error.code == "injected_durability_fault"
                end
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "unfrozen"
                assert_no_transient_debris(fixture)
                retry = freeze_fixture(fixture)
                @test retry.durability == "confirmed"
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
            end

            with_fixture(; gate_statuses=["PASS"]) do fixture
                freeze_fixture(fixture)
                failed_lock = capture_durability_fault(site) do
                    GradeLifecycle.grade_frozen_candidate(
                        fixture.campaign; launcher=fake_launcher(fixture))
                end
                @test failed_lock.error isa FreezeLifecycle.LifecycleError
                if failed_lock.error isa FreezeLifecycle.LifecycleError
                    @test failed_lock.error.code == "injected_durability_fault"
                end
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
                @test invocation_count(fixture) == 0
                assert_no_transient_debris(fixture)
                retry = GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
                @test retry.outcome == "success"
                @test invocation_count(fixture) == 1
            end

            with_fixture(; gate_statuses=["PASS"]) do fixture
                freeze_fixture(fixture)
                failed_final_lock = capture_durability_fault(site; occurrence=2) do
                    GradeLifecycle.grade_frozen_candidate(
                        fixture.campaign; launcher=fake_launcher(fixture))
                end
                @test failed_final_lock.error === nothing
                @test failed_final_lock.value.outcome == "error"
                @test failed_final_lock.value.invocation_count == 1
                @test invocation_count(fixture) == 1
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
                @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
                @test invocation_count(fixture) == 1
                assert_no_transient_debris(fixture)
            end
        end

        with_fixture(; gate_statuses=["PASS"]) do fixture
            released = capture_durability_fault(:lock_release_parent_fsync) do
                freeze_fixture(fixture)
            end
            @test released.error === nothing
            lock_path = joinpath(fixture.campaign, "lifecycle", ".freeze.atomic-create.lock")
            assert_durability_warning(
                released.value,
                lock_path;
                code="lock_removal_durability_uncertain",
            )
            @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
            assert_no_transient_debris(fixture)
        end

        for occurrence in 1:2
            with_fixture(; gate_statuses=["PASS"]) do fixture
                freeze_fixture(fixture)
                released = capture_durability_fault(
                    :lock_release_parent_fsync; occurrence) do
                    GradeLifecycle.grade_frozen_candidate(
                        fixture.campaign; launcher=fake_launcher(fixture))
                end
                @test released.error === nothing
                @test released.value.outcome == "success"
                @test invocation_count(fixture) == 1
                lock_name = occurrence == 1 ?
                    ".grade-reservation.atomic-create.lock" : ".grade-final.atomic-create.lock"
                assert_durability_warning(
                    released.value,
                    joinpath(fixture.campaign, "lifecycle", lock_name);
                    code="lock_removal_durability_uncertain",
                )
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
                assert_no_transient_debris(fixture)
            end
        end

        with_fixture(; gate_statuses=["PASS"]) do fixture
            preinstall = capture_durability_fault(:receipt_preinstall_fsync) do
                freeze_fixture(fixture)
            end
            @test preinstall.error isa FreezeLifecycle.LifecycleError
            @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "unfrozen"
            assert_no_transient_debris(fixture)
            retry = freeze_fixture(fixture)
            @test retry.durability == "confirmed"
            @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
        end

        with_fixture(; gate_statuses=["PASS"]) do fixture
            postinstall = capture_durability_fault(:receipt_postinstall_parent_fsync) do
                freeze_fixture(fixture)
            end
            @test postinstall.error === nothing
            freeze_path = joinpath(fixture.campaign, "lifecycle", "freeze_receipt.toml")
            assert_durability_warning(postinstall.value, freeze_path)
            @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
            @test_throws FreezeLifecycle.LifecycleError freeze_fixture(fixture)
            assert_no_transient_debris(fixture)
        end

        with_fixture(; gate_statuses=["PASS"]) do fixture
            freeze_fixture(fixture)
            preinstall = capture_durability_fault(:receipt_preinstall_fsync) do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            @test preinstall.error isa FreezeLifecycle.LifecycleError
            @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "frozen_once"
            @test invocation_count(fixture) == 0
            assert_no_transient_debris(fixture)
            retry = GradeLifecycle.grade_frozen_candidate(
                fixture.campaign; launcher=fake_launcher(fixture))
            @test retry.outcome == "success"
            @test invocation_count(fixture) == 1
        end

        for occurrence in 1:2
            with_fixture(; gate_statuses=["PASS"]) do fixture
                freeze_fixture(fixture)
                postinstall = capture_durability_fault(
                    :receipt_postinstall_parent_fsync; occurrence) do
                    GradeLifecycle.grade_frozen_candidate(
                        fixture.campaign; launcher=fake_launcher(fixture))
                end
                @test postinstall.error === nothing
                outcome = postinstall.value
                @test outcome.outcome == "success"
                @test outcome.invocation_count == 1
                @test invocation_count(fixture) == 1
                expected_path = joinpath(
                    fixture.campaign,
                    "lifecycle",
                    occurrence == 1 ? "grade_reservation.toml" : "grade_final.toml",
                )
                assert_durability_warning(outcome, expected_path)
                @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
                @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
                @test invocation_count(fixture) == 1
                assert_no_transient_debris(fixture)
            end
        end

        with_fixture(; gate_statuses=["PASS"]) do fixture
            freeze_fixture(fixture)
            final_preinstall = capture_durability_fault(
                :receipt_preinstall_fsync; occurrence=2) do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            @test final_preinstall.error === nothing
            @test final_preinstall.value.outcome == "error"
            @test final_preinstall.value.invocation_count == 1
            @test invocation_count(fixture) == 1
            @test FreezeLifecycle.derive_lifecycle(fixture.campaign).state == "graded"
            @test_throws FreezeLifecycle.LifecycleError GradeLifecycle.grade_frozen_candidate(
                fixture.campaign; launcher=fake_launcher(fixture))
            @test invocation_count(fixture) == 1
            assert_no_transient_debris(fixture)
        end
    end
end

function capture_lifecycle_call(function_body::Function)
    value = nothing
    captured = nothing
    try
        value = function_body()
    catch error
        captured = error
    end
    return (value=value, error=captured)
end

function assert_publication_authority_unresolved(result)
    @test result.error isa FreezeLifecycle.LifecycleError
    if result.error isa FreezeLifecycle.LifecycleError
        @test result.error.code == "publication_authority_unresolved"
    end
    @test result.value === nothing
    if result.value !== nothing && hasproperty(result.value, :durability)
        @test result.value.durability != "confirmed"
    end
end

make_receipt_writable(path::String) = chmod(path, 0o600)

@testset "writable exact lifecycle receipts never establish authority" begin
    @testset "mode-0600 freeze fails faulting and later grade invocations" begin
        with_fixture(; gate_statuses=["PASS"]) do fixture
            faulting = capture_durability_fault(
                :receipt_postinstall_parent_fsync;
                before_throw=make_receipt_writable,
            ) do
                freeze_fixture(fixture)
            end
            assert_publication_authority_unresolved(faulting)
            freeze_path = joinpath(fixture.campaign, "lifecycle", "freeze_receipt.toml")
            @test isfile(freeze_path)
            @test stat(freeze_path).mode & 0o222 != 0
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))
            @test invocation_count(fixture) == 0

            later_derive = capture_lifecycle_call() do
                FreezeLifecycle.derive_lifecycle(fixture.campaign)
            end
            assert_publication_authority_unresolved(later_derive)
            later_grade = capture_lifecycle_call() do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            assert_publication_authority_unresolved(later_grade)
            @test invocation_count(fixture) == 0
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_reservation.toml"))
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_final.toml"))
            assert_no_transient_debris(fixture)
        end
    end

    @testset "mode-0600 reservation never consumes or finalizes the attempt" begin
        with_fixture(; gate_statuses=["PASS"]) do fixture
            freeze_fixture(fixture)
            faulting = capture_durability_fault(
                :receipt_postinstall_parent_fsync;
                before_throw=make_receipt_writable,
            ) do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            assert_publication_authority_unresolved(faulting)
            reservation_path = joinpath(
                fixture.campaign, "lifecycle", "grade_reservation.toml")
            @test isfile(reservation_path)
            @test stat(reservation_path).mode & 0o222 != 0
            @test invocation_count(fixture) == 0
            @test !ispath(joinpath(fixture.campaign, "grade"))
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_final.toml"))

            later_derive = capture_lifecycle_call() do
                FreezeLifecycle.derive_lifecycle(fixture.campaign)
            end
            assert_publication_authority_unresolved(later_derive)
            later_grade = capture_lifecycle_call() do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            assert_publication_authority_unresolved(later_grade)
            @test invocation_count(fixture) == 0
            @test !ispath(joinpath(fixture.campaign, "lifecycle", "grade_final.toml"))
            assert_no_transient_debris(fixture)
        end
    end

    @testset "mode-0600 SUCCESS final never returns confirmed authority" begin
        with_fixture(; gate_statuses=["PASS"]) do fixture
            freeze_fixture(fixture)
            faulting = capture_durability_fault(
                :receipt_postinstall_parent_fsync;
                occurrence=2,
                before_throw=make_receipt_writable,
            ) do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            assert_publication_authority_unresolved(faulting)
            final_path = joinpath(fixture.campaign, "lifecycle", "grade_final.toml")
            @test isfile(final_path)
            @test stat(final_path).mode & 0o222 != 0
            @test TOML.parsefile(final_path)["outcome"] == "success"
            @test invocation_count(fixture) == 1

            later_derive = capture_lifecycle_call() do
                FreezeLifecycle.derive_lifecycle(fixture.campaign)
            end
            assert_publication_authority_unresolved(later_derive)
            later_grade = capture_lifecycle_call() do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            assert_publication_authority_unresolved(later_grade)
            @test invocation_count(fixture) == 1
            @test stat(final_path).mode & 0o222 != 0
            assert_no_transient_debris(fixture)
        end
    end

    @testset "mode-0600 ERROR final never returns confirmed authority" begin
        with_fixture(; gate_statuses=["PASS"]) do fixture
            freeze_fixture(fixture)
            faulting = capture_durability_fault(
                :receipt_postinstall_parent_fsync;
                occurrence=2,
                before_throw=make_receipt_writable,
            ) do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture; mode="error"))
            end
            assert_publication_authority_unresolved(faulting)
            final_path = joinpath(fixture.campaign, "lifecycle", "grade_final.toml")
            @test isfile(final_path)
            @test stat(final_path).mode & 0o222 != 0
            @test TOML.parsefile(final_path)["outcome"] == "error"
            @test invocation_count(fixture) == 1

            later_derive = capture_lifecycle_call() do
                FreezeLifecycle.derive_lifecycle(fixture.campaign)
            end
            assert_publication_authority_unresolved(later_derive)
            later_grade = capture_lifecycle_call() do
                GradeLifecycle.grade_frozen_candidate(
                    fixture.campaign; launcher=fake_launcher(fixture))
            end
            assert_publication_authority_unresolved(later_grade)
            @test invocation_count(fixture) == 1
            @test stat(final_path).mode & 0o222 != 0
            assert_no_transient_debris(fixture)
        end
    end
end
