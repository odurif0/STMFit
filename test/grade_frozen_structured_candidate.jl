#!/usr/bin/env julia

if !isdefined(Main, :StructuredFreezeLifecycle)
    include(joinpath(@__DIR__, "freeze_structured_unit_assignment_candidate.jl"))
end

module StructuredGradeLifecycle

using Printf
using SHA
using TOML

using Main.StructuredFreezeLifecycle
const Freeze = Main.StructuredFreezeLifecycle

const GRADE_HEADER = [
    "profile", "files", "lobes_possible", "lobes_classified", "coverage",
    "prediction_lobes", "missing_truth_lobes", "extra_prediction_lobes",
    "files_short_N", "files_extra_N", "physical_correct", "physical_errors_emitted",
    "physical_accuracy_classified", "honest_correct", "honest_uncertain",
    "honest_correct_frac", "sequence_exact", "oracle_correct",
    "oracle_accuracy_classified", "prediction_path", "grade_tsv", "note",
]
const INTEGER_FIELDS = [
    "files", "lobes_possible", "lobes_classified", "prediction_lobes",
    "missing_truth_lobes", "extra_prediction_lobes", "files_short_N", "files_extra_N",
    "physical_correct", "physical_errors_emitted", "honest_correct", "honest_uncertain",
    "sequence_exact", "oracle_correct",
]
const PERCENT_FIELDS = [
    "coverage", "physical_accuracy_classified", "honest_correct_frac",
    "oracle_accuracy_classified",
]
const FAILPOINTS = Set([
    :before_lock,
    :after_lock,
    :before_reservation,
    :after_reservation,
    :after_grade_directory,
    :before_launch_intent,
    :after_launch_intent,
    :after_process,
    :after_artifact_validation,
    :before_final,
    :after_final,
])

export GRADE_HEADER,
       GradeMetrics,
       GradeOutcome,
       command_arguments,
       grade_frozen_candidate,
       parse_grade_summary,
       production_launcher,
       recover_reserved_grade

struct GradeMetrics
    integers::Dict{String,Int}
    fractions::Dict{String,Float64}
    prediction_path::String
    grade_tsv_path::String
end

struct GradeOutcome
    state::String
    outcome::String
    process_exit_code::Int
    invocation_count::Int
    error_reason::String
    final_path::String
    final_sha256::String
    durability::String
    durability_warnings::Vector{DurabilityWarning}
end

durability_state(warnings::Vector{DurabilityWarning}) =
    isempty(warnings) ? "confirmed" : "uncertain"

function grade_outcome_from_final(final, warnings::Vector{DurabilityWarning})::GradeOutcome
    return GradeOutcome(
        "graded",
        String(final.document["outcome"]),
        Int(final.document["process_exit_code"]),
        Int(final.document["invocation_count"]),
        String(final.document["error_reason"]),
        final.path,
        final.sha256,
        durability_state(warnings),
        copy(warnings),
    )
end

function inject_failpoint(actual::Union{Nothing,Symbol}, expected::Symbol)
    actual == expected && throw(LifecycleError(
        "simulated_crash_$expected",
        "simulated crash at lifecycle phase $expected",
    ))
    return nothing
end

function validate_failpoint(failpoint::Union{Nothing,Symbol})
    failpoint === nothing || failpoint in FAILPOINTS ||
        throw(LifecycleError("invalid_failpoint", "unsupported lifecycle failpoint: $failpoint"))
end

function command_arguments(context::LifecycleContext, prediction::String, destination::String)
    exact_grade_command(context, prediction, destination)
    return [
        "julia",
        "--project=$(context.root)",
        context.report_script_path,
        "--full145-own-n",
        "--profile",
        "structured_v2=$prediction",
        "--outdir",
        destination,
    ]
end

function destination_absent(context::LifecycleContext, campaign::String)::String
    destination = joinpath(campaign, "grade")
    Freeze.reject_symlink_ancestors(context.root, destination, "grade destination")
    islink(destination) &&
        throw(LifecycleError("symlink_rejected", "grade destination is a symlink"))
    if ispath(destination)
        if isdir(destination)
            throw(LifecycleError("grade_destination_exists", "grade destination already exists"))
        elseif isfile(destination)
            throw(LifecycleError("grade_destination_exists", "grade destination is a regular file"))
        end
        throw(LifecycleError("special_file_rejected", "grade destination is a special file"))
    end
    stat(campaign).device == stat(dirname(destination)).device ||
        throw(LifecycleError("cross_filesystem_publication", "grade destination parent changed filesystem"))
    return destination
end

function reservation_document(context::LifecycleContext, freeze, environment::Dict{String,String})
    destination = joinpath(freeze.campaign, "grade")
    prediction = freeze.document["prediction_path"]
    return Dict{String,Any}(
        "schema" => "structured_grade_reservation_v1",
        "schema_version" => 1,
        "state" => "grade_reserved",
        "campaign_id" => freeze.document["campaign_id"],
        "attempt_id" => freeze.document["attempt_id"],
        "candidate_sha256" => freeze.document["candidate_sha256"],
        "model_sha256" => freeze.document["model_sha256"],
        "prediction_path" => prediction,
        "prediction_sha256" => freeze.document["prediction_sha256"],
        "repository_root" => context.root,
        "project_path" => context.project_path,
        "project_sha256" => sha256_file(context.project_path),
        "grade_wrapper_path" => context.grade_wrapper_path,
        "grade_wrapper_sha256" => sha256_file(context.grade_wrapper_path),
        "report_script_path" => context.report_script_path,
        "report_script_sha256" => sha256_file(context.report_script_path),
        "exact_command" => exact_grade_command(context, prediction, destination),
        "working_directory" => context.root,
        "report_destination" => destination,
        "environment_keys" => [key for key in ALLOWED_ENVIRONMENT_KEYS if haskey(environment, key)],
        "previous_receipt" => "freeze_receipt.toml",
        "previous_receipt_sha256" => freeze.sha256,
    )
end

function create_empty_file(path::String)
    parent = dirname(path)
    isdir(parent) && !islink(parent) ||
        throw(LifecycleError("grade_artifact_setup_failed", "artifact parent is not a directory"))
    islink(path) && throw(LifecycleError("symlink_rejected", "grade artifact is a symlink: $path"))
    if ispath(path)
        isfile(path) ||
            throw(LifecycleError("special_file_rejected", "grade artifact is a special file: $path"))
        return false
    end
    flags = Cint(0x0001 | 0x0040 | 0x0080 | 0x80000)
    descriptor = ccall(:open, Cint, (Cstring, Cint, Cint), path, flags, Cint(0o600))
    descriptor >= 0 || throw(LifecycleError(
        ispath(path) ? "grade_artifact_collision" : "grade_artifact_setup_failed",
        "exclusive grade artifact creation failed: $path",
    ))
    try
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            throw(LifecycleError("fsync_failed", "empty grade artifact fsync failed"))
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
    isfile(path) && !islink(path) ||
        throw(LifecycleError("grade_artifact_setup_failed", "grade artifact is not regular"))
    Freeze.fsync_directory(parent)
    return true
end

function ensure_grade_artifacts(reservation)
    destination = reservation.destination
    islink(destination) &&
        throw(LifecycleError("symlink_rejected", "grade destination is a symlink"))
    if !ispath(destination)
        mkdir(destination; mode=0o700)
        Freeze.fsync_directory(dirname(destination))
    end
    isdir(destination) && !islink(destination) ||
        throw(LifecycleError("special_file_rejected", "grade destination is not a real directory"))
    stat(destination).device == stat(reservation.campaign).device ||
        throw(LifecycleError("cross_filesystem_publication", "grade destination changed filesystem"))
    grades = joinpath(destination, "grades")
    islink(grades) && throw(LifecycleError("symlink_rejected", "grades directory is a symlink"))
    if !ispath(grades)
        mkdir(grades; mode=0o700)
        Freeze.fsync_directory(destination)
    end
    isdir(grades) && !islink(grades) ||
        throw(LifecycleError("special_file_rejected", "grades path is not a real directory"))
    paths = grade_artifact_paths(destination)
    for path in values(paths)
        create_empty_file(path)
    end
    return paths
end

function append_durable(path::String, text::String)
    Freeze.validate_text(text, "process evidence line")
    open(path, "a") do io
        write(io, text, '\n')
        flush(io)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 ||
            throw(LifecycleError("fsync_failed", "process evidence fsync failed"))
    end
    return nothing
end

function launch_intent_line(reservation)::String
    command_digest = bytes2hex(sha256(codeunits(reservation.document["exact_command"])))
    return "stmfit-structured-launch-intent-v1 reservation_sha256=$(reservation.sha256) command_sha256=$command_digest"
end

process_exit_line(code::Int) = "stmfit-structured-process-exit-v1 code=$code"

function process_evidence(reservation, stderr_path::String)
    bytes = snapshot_bytes(stderr_path, "grade stderr evidence")
    text = try
        String(copy(bytes))
    catch error
        throw(LifecycleError("process_evidence_invalid", sprint(showerror, error)))
    end
    lines = split(text, '\n'; keepempty=false)
    intent_lines = filter(line -> startswith(line, "stmfit-structured-launch-intent-v1 "), lines)
    exit_lines = filter(line -> startswith(line, "stmfit-structured-process-exit-v1 "), lines)
    length(intent_lines) <= 1 ||
        throw(LifecycleError("process_evidence_invalid", "duplicate launch intent markers"))
    length(exit_lines) <= 1 ||
        throw(LifecycleError("process_evidence_invalid", "duplicate process exit markers"))
    if isempty(intent_lines)
        isempty(exit_lines) ||
            throw(LifecycleError("process_evidence_invalid", "process exit exists without launch intent"))
        return (invocation_count=0, process_exit_code=-1)
    end
    only(intent_lines) == launch_intent_line(reservation) ||
        throw(LifecycleError("process_evidence_invalid", "launch intent identity mismatch"))
    if isempty(exit_lines)
        return (invocation_count=1, process_exit_code=255)
    end
    matched = match(r"^stmfit-structured-process-exit-v1 code=(-?[0-9]+)$", only(exit_lines))
    matched === nothing &&
        throw(LifecycleError("process_evidence_invalid", "process exit marker is malformed"))
    code = tryparse(Int, matched.captures[1])
    code === nothing &&
        throw(LifecycleError("process_evidence_invalid", "process exit code is outside Int range"))
    return (invocation_count=1, process_exit_code=code)
end

function production_launcher(
    arguments::Vector{String},
    environment::Dict{String,String},
    working_directory::String,
    stdout_path::String,
    stderr_path::String,
)::Int
    Set(keys(environment)) <= Set(ALLOWED_ENVIRONMENT_KEYS) ||
        throw(LifecycleError("environment_firewall_failed", "launch environment contains forbidden keys"))
    any(key -> startswith(key, "STMFIT_"), keys(environment)) &&
        throw(LifecycleError("environment_firewall_failed", "launch environment contains STMFIT variables"))
    command = setenv(Cmd(Cmd(arguments); dir=working_directory), environment)
    open(stdout_path, "a") do stdout_io
        open(stderr_path, "a") do stderr_io
            process = run(pipeline(ignorestatus(command); stdout=stdout_io, stderr=stderr_io); wait=false)
            wait(process)
            flush(stdout_io)
            flush(stderr_io)
            ccall(:fsync, Cint, (Cint,), fd(stdout_io)) == 0 ||
                throw(LifecycleError("fsync_failed", "stdout fsync failed"))
            ccall(:fsync, Cint, (Cint,), fd(stderr_io)) == 0 ||
                throw(LifecycleError("fsync_failed", "stderr fsync failed"))
            return process.exitcode
        end
    end
end

function parse_nonnegative_integer(value::String, field::String)::Int
    occursin(r"^(0|[1-9][0-9]*)$", value) ||
        throw(LifecycleError("grade_integer_invalid", "$field must be a canonical nonnegative integer"))
    parsed = tryparse(Int, value)
    parsed === nothing &&
        throw(LifecycleError("grade_integer_invalid", "$field is outside Int range"))
    return parsed
end

function parse_percentage(value::String, field::String)::Float64
    occursin(r"^(0|[1-9][0-9]?|100)\.[0-9]%$", value) ||
        throw(LifecycleError("grade_percentage_invalid", "$field has invalid percentage bytes"))
    parsed = tryparse(Float64, chop(value; tail=1))
    parsed === nothing &&
        throw(LifecycleError("grade_percentage_invalid", "$field is not numeric"))
    fraction = parsed / 100.0
    0.0 <= fraction <= 1.0 ||
        throw(LifecycleError("grade_percentage_invalid", "$field is outside [0,1]"))
    return fraction
end

function exact_percentage(numerator::Int, denominator::Int, field::String)::String
    denominator > 0 ||
        throw(LifecycleError("grade_invariant_failed", "$field denominator must be positive"))
    return @sprintf("%.1f%%", 100 * numerator / denominator)
end

function parse_grade_summary(
    summary_path::String;
    frozen_prediction::String,
    reserved_grade_directory::String,
)::GradeMetrics
    summary_path == joinpath(reserved_grade_directory, "summary.tsv") ||
        throw(LifecycleError("grade_summary_path_mismatch", "summary path differs from reservation"))
    summary = Freeze.existing_file(ROOT, summary_path, "grade summary")
    prediction = Freeze.existing_file(ROOT, frozen_prediction, "frozen prediction")
    destination = Freeze.existing_directory(ROOT, reserved_grade_directory, "reserved grade directory")
    bytes = snapshot_bytes(summary, "grade summary")
    text = try
        String(copy(bytes))
    catch error
        throw(LifecycleError("grade_summary_invalid", sprint(showerror, error)))
    end
    startswith(text, "\ufeff") &&
        throw(LifecycleError("grade_summary_invalid", "summary contains a BOM"))
    occursin('\r', text) &&
        throw(LifecycleError("grade_summary_invalid", "summary must use LF line endings"))
    endswith(text, "\n") ||
        throw(LifecycleError("grade_summary_invalid", "summary must be LF terminated"))
    lines = split(text, '\n'; keepempty=true)
    length(lines) == 3 && lines[3] == "" ||
        throw(LifecycleError("grade_row_count_invalid", "summary must contain exactly one row"))
    header = split(lines[1], '\t'; keepempty=true)
    header == GRADE_HEADER ||
        throw(LifecycleError("grade_header_invalid", "summary header differs from exact 22-column contract"))
    fields = split(lines[2], '\t'; keepempty=true)
    length(fields) == length(GRADE_HEADER) ||
        throw(LifecycleError("grade_row_invalid", "summary row has the wrong field count"))
    row = Dict(GRADE_HEADER[index] => String(fields[index]) for index in eachindex(GRADE_HEADER))
    row["profile"] == "structured_v2" ||
        throw(LifecycleError("grade_identity_mismatch", "summary profile is not structured_v2"))
    row_prediction = Freeze.existing_file(ROOT, row["prediction_path"], "summary prediction path")
    row_prediction == prediction ||
        throw(LifecycleError("grade_identity_mismatch", "summary prediction path differs from frozen path"))
    expected_grade = joinpath(destination, "grades", "structured_v2.tsv")
    row_grade = Freeze.existing_file(ROOT, row["grade_tsv"], "summary grade TSV")
    row_grade == expected_grade ||
        throw(LifecycleError("grade_path_escape", "summary grade_tsv differs from reserved grade path"))
    Freeze.path_is_within(row_grade, destination) ||
        throw(LifecycleError("grade_path_escape", "summary grade_tsv escapes reserved directory"))

    integers = Dict(field => parse_nonnegative_integer(row[field], field) for field in INTEGER_FIELDS)
    fractions = Dict(field => parse_percentage(row[field], field) for field in PERCENT_FIELDS)
    integers["files"] == 145 ||
        throw(LifecycleError("grade_invariant_failed", "files must equal 145"))
    integers["lobes_possible"] == 870 ||
        throw(LifecycleError("grade_invariant_failed", "lobes_possible must equal 870"))
    integers["physical_correct"] + integers["physical_errors_emitted"] ==
        integers["lobes_classified"] ||
        throw(LifecycleError("grade_invariant_failed", "physical totals do not equal classified lobes"))
    integers["honest_correct"] == integers["physical_correct"] ||
        throw(LifecycleError("grade_invariant_failed", "honest_correct differs from physical_correct"))
    integers["honest_correct"] + integers["honest_uncertain"] ==
        integers["lobes_possible"] ||
        throw(LifecycleError("grade_invariant_failed", "honest totals do not equal possible lobes"))
    expected_percentages = Dict(
        "coverage" => exact_percentage(integers["lobes_classified"], integers["lobes_possible"], "coverage"),
        "physical_accuracy_classified" => exact_percentage(
            integers["physical_correct"], integers["lobes_classified"],
            "physical_accuracy_classified",
        ),
        "honest_correct_frac" => exact_percentage(
            integers["honest_correct"], integers["lobes_possible"], "honest_correct_frac",
        ),
        "oracle_accuracy_classified" => exact_percentage(
            integers["oracle_correct"], integers["lobes_classified"],
            "oracle_accuracy_classified",
        ),
    )
    for (field, expected) in expected_percentages
        row[field] == expected ||
            throw(LifecycleError("grade_rounding_mismatch", "$field expected $expected, got $(row[field])"))
    end
    return GradeMetrics(integers, fractions, prediction, row_grade)
end

function artifact_hashes(paths::Dict{String,String})::Dict{String,String}
    hashes = Dict{String,String}()
    for (path_key, path) in paths
        Freeze.existing_file(ROOT, path, "grade artifact")
        hashes[replace(path_key, "_path" => "_sha256")] = sha256_file(path)
    end
    return hashes
end

function final_document(
    reservation,
    paths::Dict{String,String},
    hashes::Dict{String,String};
    outcome::String,
    error_reason::String,
    process_exit_code::Int,
    invocation_count::Int,
)
    document = Dict{String,Any}(
        "schema" => "structured_grade_final_v1",
        "schema_version" => 1,
        "state" => "graded",
        "outcome" => outcome,
        "error_reason" => error_reason,
        "campaign_id" => reservation.document["campaign_id"],
        "attempt_id" => reservation.document["attempt_id"],
        "candidate_sha256" => reservation.document["candidate_sha256"],
        "model_sha256" => reservation.document["model_sha256"],
        "profile" => "structured_v2",
        "prediction_path" => reservation.document["prediction_path"],
        "prediction_sha256" => reservation.document["prediction_sha256"],
        "process_exit_code" => process_exit_code,
        "invocation_count" => invocation_count,
        "previous_receipt" => "grade_reservation.toml",
        "previous_receipt_sha256" => reservation.sha256,
    )
    for (key, value) in paths
        document[key] = value
    end
    for (key, value) in hashes
        document[key] = value
    end
    return document
end

function finalize_grade(
    campaign::String;
    context::LifecycleContext,
    outcome::String,
    error_reason::String,
    process_exit_code::Int,
    invocation_count::Int,
    durability_warnings::Vector{DurabilityWarning}=DurabilityWarning[],
)
    outcome in ("success", "error") ||
        throw(LifecycleError("internal_error", "invalid final outcome"))
    reservation = validate_reservation_chain(context, campaign)
    lock = acquire_create_lock(reservation.lifecycle, ".grade-final.atomic-create.lock")
    warnings = copy(durability_warnings)
    final = nothing
    try
        final_path = joinpath(reservation.lifecycle, "grade_final.toml")
        if ispath(final_path) || islink(final_path)
            final = Freeze.validate_final_chain(context, campaign)
        else
            reservation = validate_reservation_chain(context, campaign)
            paths = ensure_grade_artifacts(reservation)
            if outcome == "success"
                parse_grade_summary(
                    paths["summary_tsv_path"];
                    frozen_prediction=reservation.document["prediction_path"],
                    reserved_grade_directory=reservation.destination,
                )
                isempty(error_reason) ||
                    throw(LifecycleError("internal_error", "success cannot record an error reason"))
                process_exit_code == 0 && invocation_count == 1 ||
                    throw(LifecycleError("internal_error", "success requires exit 0 and one invocation"))
            elseif invocation_count == 0
                error_reason == "crash_after_reservation_before_launch" ||
                    throw(LifecycleError("internal_error", "prelaunch error reason is not canonical"))
                process_exit_code == -1 ||
                    throw(LifecycleError("internal_error", "prelaunch error exit is not canonical"))
            else
                invocation_count == 1 ||
                    throw(LifecycleError("internal_error", "launched error requires one invocation"))
                process_exit_code != -1 ||
                    throw(LifecycleError("internal_error", "launched error uses prelaunch sentinel"))
                isempty(error_reason) &&
                    throw(LifecycleError("internal_error", "launched error requires a reason"))
            end
            hashes = artifact_hashes(paths)
            document = final_document(
                reservation,
                paths,
                hashes;
                outcome,
                error_reason,
                process_exit_code,
                invocation_count,
            )
            bytes = serialize_final(document)
            TOML.parse(String(copy(bytes)))
            publication = atomic_create_file(final_path, bytes)
            append!(warnings, publication.durability_warnings)
            final = Freeze.validate_final_chain(context, campaign)
            final.sha256 == publication.sha256 || throw(LifecycleError(
                "publication_failed", "final receipt hash changed after publication"))
            outcome == "success" && parse_grade_summary(
                paths["summary_tsv_path"];
                frozen_prediction=reservation.document["prediction_path"],
                reserved_grade_directory=reservation.destination,
            )
        end
    finally
        append!(warnings, release_create_lock(lock))
    end
    final === nothing && throw(LifecycleError(
        "publication_authority_unresolved", "final receipt authority was not established"))
    return grade_outcome_from_final(final, warnings)
end

function outcome_from_existing_final(
    context::LifecycleContext,
    campaign::String;
    durability_warnings::Vector{DurabilityWarning}=DurabilityWarning[],
    cause=nothing,
)::GradeOutcome
    final = Freeze.validate_final_chain(context, campaign)
    warnings = copy(durability_warnings)
    if cause isa LifecycleError && cause.code == "publication_authority_unresolved"
        push!(warnings, Freeze.publication_warning(final.path, cause))
    end
    return grade_outcome_from_final(final, warnings)
end

function recover_reserved_grade(
    campaign_dir::String;
    context::LifecycleContext=production_context(),
)::GradeOutcome
    lifecycle = derive_lifecycle(campaign_dir; context)
    lifecycle.state == "grade_reserved" ||
        throw(LifecycleError("recovery_not_permitted", "recovery requires grade_reserved state"))
    reservation = lifecycle.reservation
    paths = ensure_grade_artifacts(reservation)
    evidence = process_evidence(reservation, paths["stderr_path"])
    reason = evidence.invocation_count == 0 ?
        "crash_after_reservation_before_launch" : "recovery_after_reserved_launch"
    return finalize_grade(
        lifecycle.campaign;
        context,
        outcome="error",
        error_reason=reason,
        process_exit_code=evidence.process_exit_code,
        invocation_count=evidence.invocation_count,
    )
end

function grade_frozen_candidate(
    campaign_dir::String;
    context::LifecycleContext=production_context(),
    launcher::Function=production_launcher,
    failpoint::Union{Nothing,Symbol}=nothing,
    leave_reserved::Bool=false,
)::GradeOutcome
    validate_failpoint(failpoint)
    reservation_created = false
    reservation_commit_started = false
    reservation_expected_bytes = UInt8[]
    campaign = ""
    lifecycle_directory = ""
    durability_warnings = DurabilityWarning[]
    try
        Freeze.validate_context(context)
        lifecycle = derive_lifecycle(campaign_dir; context)
        lifecycle.state == "frozen_once" ||
            throw(LifecycleError("grade_not_permitted", "grading requires frozen_once state"))
        campaign = lifecycle.campaign
        freeze = lifecycle.freeze
        destination_absent(context, campaign)
        environment = grade_environment()
        document = reservation_document(context, freeze, environment)
        bytes = serialize_reservation(document)
        reservation_expected_bytes = bytes
        exact_grade_command(context, document["prediction_path"], document["report_destination"]) ==
            document["exact_command"] ||
            throw(LifecycleError("command_identity_mismatch", "reservation command is not exact"))
        lifecycle_directory = lifecycle.lifecycle

        inject_failpoint(failpoint, :before_lock)
        lock = acquire_create_lock(lifecycle_directory, ".grade-reservation.atomic-create.lock")
        try
            inject_failpoint(failpoint, :after_lock)
            locked = derive_lifecycle(campaign; context)
            locked.state == "frozen_once" ||
                throw(LifecycleError("grade_not_permitted", "another attempt already reserved grading"))
            destination_absent(context, campaign)
            Freeze.validate_freeze_chain(context, campaign)
            inject_failpoint(failpoint, :before_reservation)
            reservation_commit_started = true
            publication = atomic_create_file(
                joinpath(lifecycle_directory, "grade_reservation.toml"), bytes)
            append!(durability_warnings, publication.durability_warnings)
            reservation_created = true
        finally
            append!(durability_warnings, release_create_lock(lock))
        end

        reservation = validate_reservation_chain(context, campaign)
        inject_failpoint(failpoint, :after_reservation)
        paths = ensure_grade_artifacts(reservation)
        inject_failpoint(failpoint, :after_grade_directory)
        inject_failpoint(failpoint, :before_launch_intent)
        append_durable(paths["stderr_path"], launch_intent_line(reservation))
        inject_failpoint(failpoint, :after_launch_intent)

        arguments = command_arguments(
            context,
            reservation.document["prediction_path"],
            reservation.destination,
        )
        join(arguments, ' ') == reservation.document["exact_command"] ||
            throw(LifecycleError("command_identity_mismatch", "launch arguments differ from reservation"))
        process_exit_code = launcher(
            arguments,
            environment,
            context.root,
            paths["stdout_path"],
            paths["stderr_path"],
        )
        append_durable(paths["stderr_path"], process_exit_line(process_exit_code))
        inject_failpoint(failpoint, :after_process)
        if process_exit_code != 0
            return finalize_grade(
                campaign;
                context,
                outcome="error",
                error_reason="grader_process_failed",
                process_exit_code,
                invocation_count=1,
                durability_warnings,
            )
        end
        try
            parse_grade_summary(
                paths["summary_tsv_path"];
                frozen_prediction=reservation.document["prediction_path"],
                reserved_grade_directory=reservation.destination,
            )
        catch error
            code = error isa LifecycleError ? error.code : "internal_error"
            return finalize_grade(
                campaign;
                context,
                outcome="error",
                error_reason="invalid_grade_artifacts_$code",
                process_exit_code=0,
                invocation_count=1,
                durability_warnings,
            )
        end
        inject_failpoint(failpoint, :after_artifact_validation)
        inject_failpoint(failpoint, :before_final)
        outcome = finalize_grade(
            campaign;
            context,
            outcome="success",
            error_reason="",
            process_exit_code=0,
            invocation_count=1,
            durability_warnings,
        )
        inject_failpoint(failpoint, :after_final)
        return outcome
    catch error
        final_path = isempty(campaign) ? "" :
            joinpath(campaign, "lifecycle", "grade_final.toml")
        if !isempty(final_path) && (ispath(final_path) || islink(final_path))
            return outcome_from_existing_final(
                context,
                campaign;
                durability_warnings,
                cause=error,
            )
        end
        reservation_path = isempty(lifecycle_directory) ? "" :
            joinpath(lifecycle_directory, "grade_reservation.toml")
        if !reservation_created && reservation_commit_started &&
           !isempty(reservation_path) && (ispath(reservation_path) || islink(reservation_path))
            exact, reason = Freeze.reconcile_exact_publication(
                reservation_path, reservation_expected_bytes)
            if !exact
                error isa LifecycleError && error.code == "publication_authority_unresolved" &&
                    throw(error)
                throw(LifecycleError(
                    "publication_authority_unresolved",
                    "reservation authority is unresolved: $reason; " *
                    "original_error=$(Freeze.one_line_error(error))",
                ))
            end
            reservation = validate_reservation_chain(context, campaign)
            reservation.sha256 == bytes2hex(sha256(reservation_expected_bytes)) ||
                throw(LifecycleError(
                    "publication_authority_unresolved",
                    "reservation authority differs from the intended receipt",
                ))
            reservation_created = true
            push!(
                durability_warnings,
                Freeze.publication_warning(reservation_path, error),
            )
        end
        if reservation_created
            leave_reserved && rethrow(error)
            reservation = validate_reservation_chain(context, campaign)
            paths = ensure_grade_artifacts(reservation)
            failure_code = error isa LifecycleError ? error.code : "internal_error"
            failure_message = replace(sprint(showerror, error), r"[\t\r\n]+" => " ")
            append_durable(
                paths["stderr_path"],
                "stmfit-structured-wrapper-error-v1 code=$failure_code message=$(failure_message)",
            )
            evidence = process_evidence(reservation, paths["stderr_path"])
            reason = evidence.invocation_count == 0 ?
                "crash_after_reservation_before_launch" :
                "post_reservation_failure_$failure_code"
            return finalize_grade(
                campaign;
                context,
                outcome="error",
                error_reason=reason,
                process_exit_code=evidence.process_exit_code,
                invocation_count=evidence.invocation_count,
                durability_warnings,
            )
        end
        rethrow(error)
    end
end

function parse_cli(arguments::Vector{String})
    campaign = nothing
    recover = false
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--campaign-dir"
            campaign === nothing ||
                throw(LifecycleError("cli_error", "--campaign-dir may be supplied only once"))
            index < length(arguments) ||
                throw(LifecycleError("cli_error", "--campaign-dir requires a path"))
            campaign = arguments[index + 1]
            index += 2
        elseif startswith(argument, "--campaign-dir=")
            campaign === nothing ||
                throw(LifecycleError("cli_error", "--campaign-dir may be supplied only once"))
            campaign = split(argument, '='; limit=2)[2]
            index += 1
        elseif argument == "--recover-error"
            recover = true
            index += 1
        elseif argument in ("-h", "--help")
            println("""
            Usage: julia --project=. test/grade_frozen_structured_candidate.jl \\
              --campaign-dir PATH [--recover-error]

            Without --recover-error, atomically reserves and launches the one exact
            structured_v2 benchmark command. Recovery never relaunches; it converts an
            orphaned durable reservation into an immutable ERROR final receipt.
            """)
            return nothing
        else
            throw(LifecycleError("cli_error", "unknown argument: $argument"))
        end
    end
    campaign === nothing &&
        throw(LifecycleError("cli_error", "--campaign-dir is required"))
    return (campaign=String(campaign), recover=recover)
end

function main(arguments::Vector{String}=copy(ARGS))::Int
    try
        options = parse_cli(arguments)
        options === nothing && return 0
        outcome = options.recover ? recover_reserved_grade(options.campaign) :
            grade_frozen_candidate(options.campaign)
        write_durability_warnings(stderr, outcome.durability_warnings)
        println("GRADED state=$(outcome.state) outcome=$(outcome.outcome) receipt=$(outcome.final_path)")
        return outcome.outcome == "success" ? 0 : 1
    catch error
        failure = error isa LifecycleError ? error :
            LifecycleError("internal_error", sprint(showerror, error, catch_backtrace()))
        println(stderr, "BLOCKED [$(failure.code)]: $(failure.message)")
        return 2
    end
end

end # module StructuredGradeLifecycle

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StructuredGradeLifecycle.main())
end
