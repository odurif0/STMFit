#!/usr/bin/env julia

# Focused, non-authoritative suite for the inactive T13 activation candidate.
#
# Driver mode builds disposable bound roots and spawns focused worker
# subprocesses that load the copied evaluator from the bound root with
# `--project=<bound-root>` and `--startup-file=no --history-file=no`.  This is
# the only way to exercise the restored `realpath(@__FILE__)` root boundary.
# All records are synthetic, disposable, non-authoritative fixtures; the live
# repository Boulder, config, Project/Manifest and sealed evidence are read-only
# copy sources and are never mutated.

module StructuredEvaluatorActivationTests

using Test
using TOML
using SHA
using LinearAlgebra

include(joinpath(@__DIR__, "evaluate_structured_unit_assignment.jl"))
const E = StructuredUnitAssignmentEvaluator

const RUNTIME_V4_EVALUATOR_CONFIG =
    "config/unit_assignment_structured_evaluator_runtime_v4.toml"
const COVERAGE_MAP = joinpath(@__DIR__, "structured_evaluator_gate5_coverage.toml")
const ACTIVATION_ENVIRONMENT = E._T13_SYNTHETIC_ENVIRONMENT
const IDENTITY_REASONS = (:runtime_version_mismatch, :runtime_executable_mismatch,
                          :runtime_sysimage_mismatch)
const BOULDER_RELATIVE = ".omo/boulder.json"
const JULIA_BINARY = joinpath(Sys.BINDIR, "julia")
const WORKER_LOG_COUNTER = Ref(0)
const CONTROL_MODE = 0o444
const SOURCE_REVIEW_RELATIVE = "records/source-review.json"
const SOURCE_REPLACEMENT_BACKUP = "records/retained-original-evaluator-source.jl"
const DIRECTORY_REPLACEMENT_RELATIVE = "universe"
const DIRECTORY_REPLACEMENT_BACKUP = "records/retained-original-universe"

ENV[ACTIVATION_ENVIRONMENT] = "1"

sha_file(path::AbstractString) = bytes2hex(sha256(read(path)))
mode_text(path::AbstractString) =
    lpad(string(UInt(stat(path).mode) & UInt(0o777); base=8), 4, '0')
repo_root() = realpath(joinpath(@__DIR__, ".."))

# Parent-reserved neutral canonical scratch for CLI-through-main fixtures.
# The stage TMPDIR pathname itself can contain a forbidden firewall concept, so
# the full CLI route is exercised from this fixed scratch instead.  It is used
# only for disposable mock test roots, never an operational ledger.
neutral_cli_scratch() =
    joinpath(repo_root(), ".omo", "run-continuation", "t13-mc-v1")

function with_neutral_cli_fixture(f::Function)
    scratch = neutral_cli_scratch()
    isdir(scratch) || error("neutral CLI scratch is absent: $scratch")
    canonical = realpath(scratch)
    previous = get(ENV, "TMPDIR", nothing)
    try
        ENV["TMPDIR"] = canonical
        return f(finalize_activation_fixture(activation_fixture()))
    finally
        previous === nothing ? delete!(ENV, "TMPDIR") : (ENV["TMPDIR"] = previous)
    end
end

function write_toml(path::AbstractString, document)
    open(path, "w") do io
        TOML.print(io, document)
    end
    return path
end

function write_control_toml(path::AbstractString, document)
    isfile(path) && chmod(path, 0o644)
    write_toml(path, document)
    chmod(path, CONTROL_MODE)
    return path
end

function run_worker(root::String, op::String; foreign::Union{Nothing,String}=nothing,
                    extra::Union{Nothing,String}=nothing,
                    extra_env::Union{Nothing,Dict{String,String}}=nothing)
    arguments = String[@__FILE__, "--worker", op, root]
    foreign === nothing || push!(arguments, foreign)
    extra === nothing || push!(arguments, extra)
    base = Cmd(vcat(JULIA_BINARY, "--startup-file=no", "--history-file=no",
                    "--project=$root", arguments))
    environment = Dict{String,String}(String(key) => String(value)
                                      for (key, value) in ENV)
    if extra_env !== nothing
        for (key, value) in extra_env
            environment[String(key)] = String(value)
        end
    end
    cmd = Cmd(base; dir=root, env=environment)
    log_root = joinpath(get(ENV, "TMPDIR", root), "worker-logs")
    mkpath(log_root)
    # Diagnostic retention: worker logs stay under the bounded TMPDIR for
    # post-run inspection instead of an atexit cleanup.
    log_dir = mktempdir(log_root; prefix="$(op)-", cleanup=false)
    write(joinpath(log_dir, "argv.txt"),
          join(vcat(JULIA_BINARY, "--startup-file=no", "--history-file=no",
                    "--project=$root", arguments), "\n") * "\n")
    write(joinpath(log_dir, "cwd.txt"), root * "\n")
    allowlist = ("PATH", "HOME", "LANG", "LC_ALL", "TMPDIR", "TMP", "TEMP",
                 "JULIA_DEPOT_PATH", "JULIA_LOAD_PATH", "JULIA_PROJECT",
                 "JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                 "OMP_NUM_THREADS", "JULIA_PKG_OFFLINE",
                 "STMFIT_T13_SYNTHETIC_EXECUTION_ONLY", "STMFIT_HARNESS_PROBE",
                 "STMFIT_DATA_DIR", "STMFIT_T13_PUBLIC_REVALIDATION",
                 "JULIA_DEBUG")
    write(joinpath(log_dir, "env.txt"),
          join(["$(key)=$(get(environment, key, ""))" for key in allowlist],
               "\n") * "\n")
    stdout_path = joinpath(log_dir, "stdout.txt")
    stderr_path = joinpath(log_dir, "stderr.txt")
    process = run(pipeline(ignorestatus(cmd); stdout=stdout_path,
                           stderr=stderr_path); wait=true)
    exitcode = process.exitcode
    write(joinpath(log_dir, "exit-code.txt"), string(exitcode) * "\n")
    return (success=exitcode == 0, exitcode=exitcode,
            stdout=read(stdout_path, String), stderr=read(stderr_path, String),
            log_dir=log_dir)
end

function _load_fixture_module(root::String)
    sandbox = Module(gensym(:FixtureEvaluatorHost))
    Base.include(sandbox, joinpath(root, "test", "evaluate_structured_unit_assignment.jl"))
    return getfield(sandbox, :StructuredUnitAssignmentEvaluator)
end

function _worker_capture(root::String)
    document = Dict{String,Any}(
        "active_project" => something(Base.active_project(), ""),
        "cwd" => pwd(),
        "startup_file" => Base.JLOptions().startupfile == 2,
        "history_file" => Base.JLOptions().historyfile == 0,
        "threads" => Base.Threads.nthreads(),
        "blas_threads" => LinearAlgebra.BLAS.get_num_threads(),
        "depot_path" => first(DEPOT_PATH),
        "load_path" => String.(LOAD_PATH),
        "offline" => get(ENV, "JULIA_PKG_OFFLINE", "") in ("true", "1"),
        "forbidden_environment" => collect(E._T13_FORBIDDEN_ENVIRONMENT),
    )
    write_toml(joinpath(root, "captured-launch.toml"), document)
    return nothing
end

function _worker_options(root::String)
    options = Dict{String,String}()
    for (key, value) in TOML.parsefile(joinpath(root, "options.toml"))
        options[String(key)] = String(value)
    end
    return options
end

function _worker_authorize(root::String, foreign::Union{Nothing,String})
    fixture = _load_fixture_module(foreign === nothing ? root : foreign)
    return Base.invokelatest(_worker_authorize_impl, fixture, root)
end

function _worker_authorize_impl(fixture::Module, root::String)
    options = _worker_options(root)
    try
        fixture._authorize_t13_activation(options)
        println("WORKER_STATUS=ok")
    catch error
        error isa fixture.EvaluatorError || rethrow()
        println("WORKER_STATUS=blocked")
        println("WORKER_REASON=", error.reason)
        println("WORKER_MESSAGE=", repr(sprint(showerror, error)))
    end
    return nothing
end

function _worker_toctou(root::String, relative::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_toctou_impl, fixture, root, relative)
end

function _worker_toctou_impl(fixture::Module, root::String, relative::String)
    options = _worker_options(root)
    try
        authorization = fixture._authorize_t13_activation(options)
        write(joinpath(root, relative), "post-auth dependency mutation\n")
        context = fixture._ProductionContext(
            root, authorization.snapshots, authorization.inventories, "",
            authorization.bindings)
        try
            fixture._verify_context(context)
            println("WORKER_STATUS=ok")
        catch error
            error isa fixture.EvaluatorError || rethrow()
            println("WORKER_STATUS=blocked")
            println("WORKER_REASON=", error.reason)
        end
    catch error
        error isa fixture.EvaluatorError || rethrow()
        println("WORKER_STATUS=blocked")
        println("WORKER_REASON=", error.reason)
    end
    return nothing
end

function _worker_expect_blocked(fixture::Module, label::String, thunk::Function)
    try
        thunk()
        println("WORKER_", label, "_REASON=ok")
        println("WORKER_", label, "_MESSAGE=")
        return "ok"
    catch error
        error isa fixture.EvaluatorError || rethrow()
        println("WORKER_", label, "_REASON=", error.reason)
        println("WORKER_", label, "_MESSAGE=", repr(sprint(showerror, error)))
        return String(error.reason)
    end
end

function _worker_expect_publication_rejection(fixture::Module, label::String,
                                              thunk::Function)
    try
        thunk()
        println("WORKER_", label, "_STATE=committed")
        println("WORKER_", label, "_REASON=none")
        println("WORKER_", label, "_MESSAGE=")
    catch error
        error isa fixture._PublicationError || rethrow()
        println("WORKER_", label, "_STATE=", error.state)
        println("WORKER_", label, "_REASON=", error.reason)
        println("WORKER_", label, "_MESSAGE=", repr(sprint(showerror, error)))
    end
    return nothing
end

function _worker_ordinary_authority(fixture::Module, root::String)
    # Explicit non-authoritative ordinary authority fixture: only the approved
    # runtime-v4 config snapshot is collected, through the real `_snapshot`
    # helper.  This is enough for `_assemble_production_context`'s merge and
    # `_verify_context` revalidation and for the activated blocker serializer's
    # config digest, but it does NOT exercise the full ordinary authority
    # collector, the ordinary input inventory collector, or the heavy positive
    # T8/T11/T12 loaders.
    seen = Dict{Tuple{UInt64,UInt64},String}()
    config = fixture._snapshot(root, RUNTIME_V4_EVALUATOR_CONFIG,
                              fixture._RUNTIME_V4_CONFIG_SHA256,
                              "ordinary authority fixture config", seen)
    return (snapshots=[config],
            inventories=Tuple{String,Vector{String}}[])
end

function _worker_context_boundary(root::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_context_boundary_impl, fixture, root)
end

function _worker_context_boundary_impl(fixture::Module, root::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    # Unit-level boundary: the real receiving-options validator is the
    # production gate that runs before any context assembly.
    fixture._validate_activation_receiving_options(authorization, root, options)
    println("WORKER_RECEIVING_VALID=ok")
    _worker_expect_blocked(fixture, "RECEIVING_ROOT", () ->
        fixture._validate_activation_receiving_options(
            authorization, root * "-unbound-root", options))
    mutated_options = copy(options)
    mutated_options["--out-dir"] = "other-output"
    _worker_expect_blocked(fixture, "RECEIVING_COMMAND", () ->
        fixture._validate_activation_receiving_options(
            authorization, root, mutated_options))
    # Assemble through the actual shared production-context builder.  The
    # returned context must still carry the original activation bindings and
    # every original protected path, mode, identity and snapshot.
    ordinary = _worker_ordinary_authority(fixture, root)
    context = fixture._assemble_production_context(root, options, ordinary;
                                                   activation=authorization)
    retained = context.root == root &&
               context.activation === authorization.bindings &&
               context.activation.destination ==
                   authorization.bindings.destination
    modes_complete = keys(context.activation.modes) ==
                     keys(authorization.bindings.modes) &&
                     all(context.activation.modes[path] ==
                         authorization.bindings.modes[path]
                         for path in keys(authorization.bindings.modes))
    identities_complete = context.activation.directory_identities ==
                          authorization.bindings.directory_identities
    context_by_path = Dict(snapshot.path => snapshot
                           for snapshot in context.snapshots)
    snapshots_complete = length(context.snapshots) ==
                         length(authorization.snapshots) + 1 &&
                         all(begin
                                 match = get(context_by_path, snapshot.path, nothing)
                                 match !== nothing &&
                                     match.sha256 == snapshot.sha256 &&
                                     (match.device, match.inode, match.nlink) ==
                                     (snapshot.device, snapshot.inode,
                                      snapshot.nlink) &&
                                     match.bytes == snapshot.bytes
                             end for snapshot in authorization.snapshots)
    inventories_complete = context.inventories == authorization.inventories
    println("WORKER_CONTEXT_ACTIVATION=", retained ? "retained" : "lost")
    println("WORKER_CONTEXT_MODES=", modes_complete ? "complete" : "mismatch")
    println("WORKER_CONTEXT_IDENTITIES=", identities_complete ? "complete" : "mismatch")
    println("WORKER_CONTEXT_SNAPSHOTS=", snapshots_complete ? "complete" : "mismatch")
    println("WORKER_CONTEXT_INVENTORIES=", inventories_complete ? "complete" : "mismatch")
    # Activation-only record drift after context construction.  The byte change
    # requires restoring the captured read-only mode so the identity/mode
    # guards are not the first boundary that trips.
    record = joinpath(root, SOURCE_REVIEW_RELATIVE)
    record_mode = UInt(stat(record).mode) & UInt(0o777)
    chmod(record, 0o644)
    write(record, "mutated activation-only record after context construction\n")
    chmod(record, record_mode)
    println("WORKER_RECORD_MODE_RESTORED=",
            (UInt(stat(record).mode) & UInt(0o777)) == record_mode ? "true" : "false")
    _worker_expect_blocked(fixture, "VERIFY", () ->
        fixture._verify_context(context))
    # Simulated downstream loader error (no heavy loader is invoked): the real
    # activated blocker serializer and publisher call shape must consume the
    # assembled context and reject the drifted authorization.
    loader_error = fixture.EvaluatorError(
        :BLOCKED, :simulated_downstream_loader_failure,
        "simulated downstream loader failure before publication")
    blocker_files = fixture._blocker_files(loader_error; context=context)
    blocker = TOML.parse(String(blocker_files["receipt.toml"]))
    println("WORKER_BLOCKER_SCHEMA=", blocker["schema"])
    println("WORKER_BLOCKER_SCOPE=", blocker["activation_scope"])
    println("WORKER_BLOCKER_RECEIPT_SHA=", blocker["activation_receipt_sha256"])
    println("WORKER_BLOCKER_REASON=", blocker["reason"])
    _worker_expect_publication_rejection(fixture, "PUBLISH", () ->
        fixture._publish_atomic(root, options["--out-dir"], blocker_files;
                                context=context))
    destination = joinpath(root, options["--out-dir"])
    println("WORKER_DESTINATION=", ispath(destination) ? "present" : "absent")
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_replace_source(root::String)
    fixture = _load_fixture_module(root)
    relative = first(fixture._T13_SOURCE_CLOSURE)
    return Base.invokelatest(_worker_replace_file_impl, fixture, root, relative,
                             SOURCE_REPLACEMENT_BACKUP)
end

function _worker_main_root(root::String)
    path = joinpath(root, "main-root.txt")
    isfile(path) || error("main-root.txt is absent")
    # `strip(read(path, String))` returns a SubString; the production helpers
    # accept `String` only, so convert once here.
    return String(strip(read(path, String)))
end

function _worker_control_report(root::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_control_report_impl, fixture, root)
end

function _worker_control_report_impl(fixture::Module, root::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    control = authorization.bindings.control
    main = _worker_main_root(root)
    decoy = joinpath(root, BOULDER_RELATIVE)
    println("WORKER_CONTROL_ROOT=", control.root)
    println("WORKER_CONTROL_FILE=", control.file.path)
    println("WORKER_CONTROL_SHA=", control.file.sha256)
    println("WORKER_DECOY_SHA=", sha_file(decoy))
    println("WORKER_CONTROL_ROOT_IS_MAIN=",
            control.root == main ? "true" : "false")
    println("WORKER_CONTROL_FILE_FIXED=",
            control.file.path == joinpath(main, E._T13_CONTROL_FILENAME) ?
                "true" : "false")
    println("WORKER_CONTROL_SNAPSHOT_ORDINARY=",
            any(snapshot.path == control.file.path
                for snapshot in authorization.snapshots) ? "true" : "false")
    ordinary = _worker_ordinary_authority(fixture, root)
    context = fixture._assemble_production_context(root, options, ordinary;
                                                   activation=authorization)
    records = fixture._manifest_records(
        root, context.snapshots, context.inventories,
        [fixture._T13_ACTIVATION_AUTHORITY_SCHEMA];
        control=context.activation.control)
    tagged = [record for record in records
              if startswith(record, "main_control\t")]
    file_records = [record for record in records if startswith(record, "file\t")]
    println("WORKER_MANIFEST_TAGGED=", length(tagged))
    println("WORKER_MANIFEST_CONTROL_RECORD=", isempty(tagged) ? "" : tagged[1])
    println("WORKER_MANIFEST_ESCAPES=",
            any(record -> occursin("../", record), records) ? "true" : "false")
    println("WORKER_MANIFEST_FILES_UNDER_ROOT=",
            all(record -> !startswith(split(record, '\t')[2], ".."),
                file_records) ? "true" : "false")
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_control_preauthorize(root::String, mutation::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_control_preauthorize_impl, fixture, root,
                             mutation)
end

function _worker_control_preauthorize_impl(fixture::Module, root::String,
                                           mutation::String)
    options = _worker_options(root)
    main = _worker_main_root(root)
    parent = joinpath(main, E._T13_CONTROL_PARENT)
    control_file = joinpath(main, E._T13_CONTROL_FILENAME)
    decoy = joinpath(root, BOULDER_RELATIVE)
    original_mode = UInt(stat(control_file).mode) & UInt(0o777)
    original_parent_mode = UInt(stat(parent).mode) & UInt(0o777)
    original_root_mode = UInt(stat(main).mode) & UInt(0o777)
    original_file_identity =
        (UInt64(stat(control_file).device), UInt64(stat(control_file).inode))
    original_parent_identity =
        (UInt64(stat(parent).device), UInt64(stat(parent).inode))
    original_root_identity =
        (UInt64(stat(main).device), UInt64(stat(main).inode))
    backup_root = joinpath(root, "records", "retained-main-control")
    mkpath(backup_root)
    if mutation in ("missing-file", "wrong-hash", "unreadable")
        # A valid R-local decoy with the reviewed bytes: a fallback to R would
        # find "right" bytes and must still be rejected.
        chmod(decoy, 0o644)
        write(decoy, read(control_file))
        chmod(decoy, original_mode)
    end
    if mutation == "missing-file"
        rm(control_file; force=true)
    elseif mutation == "missing-root"
        rm(main; recursive=true, force=true)
    elseif mutation == "wrong-hash"
        chmod(control_file, 0o644)
        write(control_file, "mutated Main control bytes\n")
        chmod(control_file, original_mode)
    elseif mutation == "unreadable"
        chmod(control_file, 0o000)
    elseif mutation == "mode-file"
        chmod(control_file, 0o644)
    elseif mutation == "mode-parent"
        chmod(parent, original_parent_mode == 0o755 ? 0o700 : 0o755)
    elseif mutation == "mode-root"
        chmod(main, original_root_mode == 0o700 ? 0o755 : 0o700)
    elseif mutation == "parent-replace-initial"
        backup = joinpath(backup_root, "omo")
        mv(parent, backup)
        mkdir(parent)
        chmod(parent, original_parent_mode)
        mv(joinpath(backup, "boulder.json"), control_file)
        new_parent_identity =
            (UInt64(stat(parent).device), UInt64(stat(parent).inode))
        new_file_identity =
            (UInt64(stat(control_file).device), UInt64(stat(control_file).inode))
        println("WORKER_PARENT_INODE_CHANGED=",
                new_parent_identity != original_parent_identity ? "true" : "false")
        println("WORKER_FILE_INODE_PRESERVED=",
                new_file_identity == original_file_identity ? "true" : "false")
        println("WORKER_BACKUP_PARENT_INODE=", UInt64(stat(backup).inode))
        println("WORKER_ORIGINAL_PARENT_INODE=", original_parent_identity[2])
    elseif mutation == "root-replace-initial"
        backup = joinpath(backup_root, "main")
        mv(main, backup)
        mkdir(main)
        chmod(main, original_root_mode)
        mv(joinpath(backup, ".omo"), parent)
        new_root_identity =
            (UInt64(stat(main).device), UInt64(stat(main).inode))
        new_parent_identity =
            (UInt64(stat(parent).device), UInt64(stat(parent).inode))
        println("WORKER_ROOT_INODE_CHANGED=",
                new_root_identity != original_root_identity ? "true" : "false")
        println("WORKER_PARENT_INODE_PRESERVED=",
                new_parent_identity == original_parent_identity ? "true" : "false")
        println("WORKER_BACKUP_ROOT_INODE=", UInt64(stat(backup).inode))
        println("WORKER_ORIGINAL_ROOT_INODE=", original_root_identity[2])
    end
    println("WORKER_MUTATION=", mutation)
    if mutation == "r-decoy"
        chmod(decoy, 0o644)
        write(decoy, "mutated R-local decoy before authorization\n")
        chmod(decoy, 0o644)
        authorization = fixture._authorize_t13_activation(options)
        println("WORKER_CONTROL_SHA=", authorization.bindings.control.file.sha256)
        println("WORKER_DECOY_SHA=", sha_file(decoy))
        println("WORKER_AUTHORIZE_REASON=ok")
    elseif mutation in ("parent-replace-initial", "root-replace-initial")
        # Initial acquisition: the replacement happened before the
        # authorization capture, so the new directory identity must become the
        # captured one (and no revalidation may reject it).
        authorization = fixture._authorize_t13_activation(options)
        println("WORKER_AUTHORIZE_REASON=ok")
        println("WORKER_CAPTURED_ROOT_INODE=",
                authorization.bindings.control.root_identity[2])
        println("WORKER_CAPTURED_PARENT_INODE=",
                authorization.bindings.control.parent_identity[2])
        println("WORKER_CURRENT_ROOT_INODE=", UInt64(stat(main).inode))
        println("WORKER_CURRENT_PARENT_INODE=", UInt64(stat(parent).inode))
        println("WORKER_CAPTURED_FILE_INODE=",
                authorization.bindings.control.file.inode)
        println("WORKER_ORIGINAL_FILE_INODE=", original_file_identity[2])
    else
        _worker_expect_blocked(fixture, "AUTHORIZE", () ->
            fixture._authorize_t13_activation(options))
    end
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_control_replacement(root::String, mutation::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_control_replacement_impl, fixture, root,
                             mutation)
end

function _worker_control_replacement_impl(fixture::Module, root::String,
                                          mutation::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    ordinary = _worker_ordinary_authority(fixture, root)
    context = fixture._assemble_production_context(root, options, ordinary;
                                                   activation=authorization)
    println("WORKER_CONTEXT=assembled")
    control = authorization.bindings.control
    main = _worker_main_root(root)
    parent = joinpath(main, E._T13_CONTROL_PARENT)
    control_file = joinpath(main, E._T13_CONTROL_FILENAME)
    backup_root = joinpath(root, "records", "retained-main-control")
    mkpath(backup_root)
    original_sha = sha_file(control_file)
    original_mode = UInt(stat(control_file).mode) & UInt(0o777)
    original_parent_mode = UInt(stat(parent).mode) & UInt(0o777)
    original_root_mode = UInt(stat(main).mode) & UInt(0o777)
    original_file_identity =
        (UInt64(stat(control_file).device), UInt64(stat(control_file).inode))
    original_parent_identity =
        (UInt64(stat(parent).device), UInt64(stat(parent).inode))
    original_root_identity =
        (UInt64(stat(main).device), UInt64(stat(main).inode))
    if mutation == "r-decoy-mutate"
        decoy = joinpath(root, BOULDER_RELATIVE)
        chmod(decoy, 0o644)
        write(decoy, "mutated R-local decoy after context assembly\n")
        chmod(decoy, 0o644)
        println("WORKER_CONTROL_SHA=", control.file.sha256)
        println("WORKER_DECOY_SHA=", sha_file(decoy))
        fixture._verify_context(context)
        println("WORKER_VERIFY=ok")
        println("WORKER_STATUS=ok")
        return nothing
    elseif mutation == "file-inode"
        backup = joinpath(backup_root, "boulder.json")
        mv(control_file, backup)
        cp(backup, control_file; force=true)
        chmod(control_file, original_mode)
        replacement_identity =
            (UInt64(stat(control_file).device), UInt64(stat(control_file).inode))
        println("WORKER_FILE_INODE_CHANGED=",
                replacement_identity != original_file_identity ? "true" : "false")
        println("WORKER_BYTES_EQUAL=",
                sha_file(control_file) == original_sha ? "true" : "false")
        println("WORKER_MODE_EQUAL=",
                (UInt(stat(control_file).mode) & UInt(0o777)) == original_mode ?
                    "true" : "false")
        println("WORKER_BACKUP_INODE=", UInt64(stat(backup).inode))
        println("WORKER_ORIGINAL_INODE=", original_file_identity[2])
        println("WORKER_REPLACEMENT_INODE=", replacement_identity[2])
        println("WORKER_AUTH_INODE=", control.file.inode)
        # The assembled context must retain the dedicated Main-control object
        # (with its original captured inode), not an ordinary snapshot merge.
        println("WORKER_CONTEXT_RETAINS_CONTROL=",
                context.activation.control === control ? "true" : "false")
        # A fresh observation through the real helper is logged but never
        # adopted or merged into the retained authorization.
        fresh = fixture._snapshot(main, E._T13_CONTROL_FILENAME, original_sha,
                                  "fresh control snapshot",
                                  Dict{Tuple{UInt64,UInt64},String}())
        println("WORKER_FRESH_INODE=", fresh.inode)
        println("WORKER_FRESH_ADOPTED=false")
        # Reassembly with the original authorization must reject: the retained
        # control object still holds the original captured inode.
        _worker_expect_blocked(fixture, "ASSEMBLER", () ->
            fixture._assemble_production_context(root, options, ordinary;
                                                 activation=authorization))
    elseif mutation == "file-bytes"
        chmod(control_file, 0o644)
        write(control_file, "mutated Main control bytes\n")
        chmod(control_file, original_mode)
        println("WORKER_BYTES_EQUAL=",
                sha_file(control_file) == original_sha ? "true" : "false")
    elseif mutation == "file-mode"
        chmod(control_file, 0o644)
        println("WORKER_MODE_EQUAL=",
                (UInt(stat(control_file).mode) & UInt(0o777)) == original_mode ?
                    "true" : "false")
    elseif mutation == "parent-mode"
        chmod(parent, original_parent_mode == 0o755 ? 0o700 : 0o755)
    elseif mutation == "root-mode"
        chmod(main, original_root_mode == 0o700 ? 0o755 : 0o700)
    elseif mutation == "parent-inode"
        backup = joinpath(backup_root, "omo")
        mv(parent, backup)
        mkdir(parent)
        chmod(parent, original_parent_mode)
        mv(joinpath(backup, "boulder.json"), control_file)
        new_parent_identity =
            (UInt64(stat(parent).device), UInt64(stat(parent).inode))
        new_file_identity =
            (UInt64(stat(control_file).device), UInt64(stat(control_file).inode))
        println("WORKER_PARENT_INODE_CHANGED=",
                new_parent_identity != original_parent_identity ? "true" : "false")
        println("WORKER_FILE_INODE_PRESERVED=",
                new_file_identity == original_file_identity ? "true" : "false")
        println("WORKER_BACKUP_PARENT_INODE=", UInt64(stat(backup).inode))
        println("WORKER_ORIGINAL_PARENT_INODE=", original_parent_identity[2])
    elseif mutation == "root-inode"
        backup = joinpath(backup_root, "main")
        mv(main, backup)
        mkdir(main)
        chmod(main, original_root_mode)
        mv(joinpath(backup, ".omo"), parent)
        new_root_identity =
            (UInt64(stat(main).device), UInt64(stat(main).inode))
        new_parent_identity =
            (UInt64(stat(parent).device), UInt64(stat(parent).inode))
        println("WORKER_ROOT_INODE_CHANGED=",
                new_root_identity != original_root_identity ? "true" : "false")
        println("WORKER_PARENT_INODE_PRESERVED=",
                new_parent_identity == original_parent_identity ? "true" : "false")
        println("WORKER_BACKUP_ROOT_INODE=", UInt64(stat(backup).inode))
        println("WORKER_ORIGINAL_ROOT_INODE=", original_root_identity[2])
    else
        error("unknown control replacement mutation: $mutation")
    end
    _worker_expect_blocked(fixture, "VERIFY", () ->
        fixture._verify_context(context))
    loader_error = fixture.EvaluatorError(
        :BLOCKED, :simulated_downstream_loader_failure,
        "simulated downstream loader failure before publication")
    blocker_files = fixture._blocker_files(loader_error; context=context)
    _worker_expect_publication_rejection(fixture, "PUBLISH", () ->
        fixture._publish_atomic(root, options["--out-dir"], blocker_files;
                                context=context))
    destination = joinpath(root, options["--out-dir"])
    println("WORKER_DESTINATION=", ispath(destination) ? "present" : "absent")
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_publication_hooks(root::String, mutation::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_publication_hooks_impl, fixture, root,
                             mutation)
end

function _worker_publication_hooks_impl(fixture::Module, root::String,
                                        mutation::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    ordinary = _worker_ordinary_authority(fixture, root)
    context = fixture._assemble_production_context(root, options, ordinary;
                                                   activation=authorization)
    main = _worker_main_root(root)
    control_file = joinpath(main, E._T13_CONTROL_FILENAME)
    control_mode = UInt(stat(control_file).mode) & UInt(0o777)
    mutate_control = function ()
        chmod(control_file, 0o644)
        write(control_file, "publication-hook Main control mutation\n")
        chmod(control_file, control_mode)
        return nothing
    end
    files = Dict{String,Vector{UInt8}}(name => Vector{UInt8}(name * "\n")
                                       for name in E._FILES)
    destination = joinpath(root, options["--out-dir"])
    if mutation == "existing-fast-path"
        mkpath(destination)
        chmod(destination, 0o700)
        for (name, bytes) in files
            write(joinpath(destination, name), bytes)
            chmod(joinpath(destination, name), 0o644)
        end
        mutate_control()
        _worker_expect_publication_rejection(fixture, "PUBLISH", () ->
            fixture._publish_atomic(root, options["--out-dir"], files;
                                    context=context))
        println("WORKER_DESTINATION=",
                ispath(destination) ? "present" : "absent")
        println("WORKER_DESTINATION_FILES=",
                isdir(destination) ? length(readdir(destination)) : -1)
    elseif mutation == "preferred-rename-destination-exists"
        hook = function (event, args...)
            event === :rename_noreplace || return nothing
            dest = String(args[2])
            mkpath(dest)
            chmod(dest, 0o700)
            for (name, bytes) in files
                write(joinpath(dest, name), bytes)
                chmod(joinpath(dest, name), 0o644)
            end
            mutate_control()
            return 17
        end
        fixture._PUBLICATION_HOOK[] = hook
        try
            _worker_expect_publication_rejection(fixture, "PUBLISH", () ->
                fixture._publish_atomic(root, options["--out-dir"], files;
                                        context=context))
        finally
            fixture._PUBLICATION_HOOK[] = nothing
        end
        println("WORKER_DESTINATION=",
                ispath(destination) ? "present" : "absent")
        println("WORKER_DESTINATION_FILES=",
                isdir(destination) ? length(readdir(destination)) : -1)
    elseif mutation == "postcommit"
        # The actual commit event differs by filesystem: the preferred
        # rename path emits `after_preferred_rename`, while the unsupported-
        # rename fallback emits `after_receipt_link`.  Subscribe to both and
        # mutate exactly once at whichever event really fires.
        mutation_count = Ref(0)
        observed_events = String[]
        control_sha_before = sha_file(control_file)
        hook = function (event, args...)
            event in (:after_preferred_rename, :after_receipt_link) || return nothing
            mutation_count[] == 0 || return nothing
            mutation_count[] += 1
            push!(observed_events, String(event))
            mutate_control()
            return nothing
        end
        fixture._PUBLICATION_HOOK[] = hook
        try
            _worker_expect_publication_rejection(fixture, "PUBLISH", () ->
                fixture._publish_atomic(root, options["--out-dir"], files;
                                        context=context))
        finally
            fixture._PUBLICATION_HOOK[] = nothing
        end
        println("WORKER_POSTCOMMIT_EVENTS=", join(observed_events, ","))
        println("WORKER_POSTCOMMIT_MUTATIONS=", mutation_count[])
        println("WORKER_CONTROL_SHA_BEFORE=", control_sha_before)
        println("WORKER_CONTROL_SHA_AFTER=", sha_file(control_file))
        println("WORKER_CONTROL_MODE_RESTORED=",
                (UInt(stat(control_file).mode) & UInt(0o777)) == control_mode ?
                    "true" : "false")
        println("WORKER_DESTINATION=",
                ispath(destination) ? "present" : "absent")
        println("WORKER_DESTINATION_FILES=",
                isdir(destination) ? length(readdir(destination)) : -1)
    else
        error("unknown publication hook mutation: $mutation")
    end
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_replace_file_impl(fixture::Module, root::String, relative::String,
                                   backup_relative::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    original = joinpath(root, relative)
    backup = joinpath(root, backup_relative)
    mkpath(dirname(backup))
    original_sha = sha_file(original)
    original_mode = UInt(stat(original).mode) & UInt(0o777)
    original_identity = (UInt64(stat(original).device), UInt64(stat(original).inode))
    # Actual filesystem replacement: rename the original pathname to a retained
    # backup (which keeps the original inode allocated) and create a new file
    # with identical bytes and mode at the original pathname.
    mv(original, backup)
    cp(backup, original; force=true)
    chmod(original, original_mode)
    replacement_identity = (UInt64(stat(original).device), UInt64(stat(original).inode))
    replacement_mode = UInt(stat(original).mode) & UInt(0o777)
    bound = [snapshot for snapshot in authorization.snapshots
             if snapshot.path == original]
    length(bound) == 1 ||
        error("authorization does not bind exactly one snapshot for $relative")
    println("WORKER_RELATIVE=", relative)
    println("WORKER_INODE_CHANGED=",
            replacement_identity != original_identity ? "true" : "false")
    println("WORKER_BYTES_EQUAL=", sha_file(original) == original_sha ? "true" : "false")
    println("WORKER_MODE_EQUAL=", replacement_mode == original_mode ? "true" : "false")
    println("WORKER_ORIGINAL_INODE=", original_identity[2])
    println("WORKER_REPLACEMENT_INODE=", replacement_identity[2])
    println("WORKER_AUTH_INODE=", only(bound).inode)
    # The retained backup still holds the original inode, so the replacement
    # inode cannot be a reuse of the original allocation.
    backup_identity = (UInt64(stat(backup).device), UInt64(stat(backup).inode))
    println("WORKER_BACKUP_INODE=", backup_identity[2])
    println("WORKER_BACKUP_BYTES_EQUAL=",
            sha_file(backup) == original_sha ? "true" : "false")
    println("WORKER_BACKUP_MODE_EQUAL=",
            (UInt(stat(backup).mode) & UInt(0o777)) == original_mode ? "true" : "false")
    # Fresh ordinary snapshot through the real helper, then the actual merge and
    # assembler: the original authorization identity must not be rebased.
    seen = Dict{Tuple{UInt64,UInt64},String}()
    fresh = fixture._snapshot(root, relative, original_sha,
                              "fresh ordinary snapshot for $relative", seen)
    println("WORKER_FRESH_INODE=", fresh.inode)
    _worker_expect_blocked(fixture, "MERGE", () ->
        fixture._merge_snapshots(authorization.snapshots,
                                 [fresh], "activation binding"))
    _worker_expect_blocked(fixture, "ASSEMBLER", () ->
        fixture._assemble_production_context(
            root, options,
            (snapshots=[fresh],
             inventories=Tuple{String,Vector{String}}[]);
            activation=authorization))
    destination = joinpath(root, options["--out-dir"])
    println("WORKER_DESTINATION=", ispath(destination) ? "present" : "absent")
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_replace_directory(root::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_replace_directory_impl, fixture, root)
end

function _worker_replace_directory_impl(fixture::Module, root::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    ordinary = _worker_ordinary_authority(fixture, root)
    context = fixture._assemble_production_context(root, options, ordinary;
                                                   activation=authorization)
    println("WORKER_CONTEXT=assembled")
    original = joinpath(root, DIRECTORY_REPLACEMENT_RELATIVE)
    backup = joinpath(root, DIRECTORY_REPLACEMENT_BACKUP)
    mkpath(dirname(backup))
    original_mode = UInt(stat(original).mode) & UInt(0o777)
    original_directory_identity =
        (UInt64(stat(original).device), UInt64(stat(original).inode))
    names = sort(readdir(original))
    member_state = Dict{String,NamedTuple}()
    for name in names
        path = joinpath(original, name)
        member_state[name] = (
            identity=(UInt64(stat(path).device), UInt64(stat(path).inode)),
            sha256=sha_file(path),
            mode=UInt(stat(path).mode) & UInt(0o777))
    end
    # Actual directory replacement with preserved member files: plain renames
    # keep every member dev/inode/bytes/name/mode, while the directory at the
    # original pathname gets a new inode.  No symlink or hardlink alias is used.
    mv(original, backup)
    mkdir(original)
    chmod(original, original_mode)
    for name in names
        mv(joinpath(backup, name), joinpath(original, name))
    end
    replacement_directory_identity =
        (UInt64(stat(original).device), UInt64(stat(original).inode))
    names_after = sort(readdir(original))
    member_identities_equal = true
    member_bytes_equal = true
    member_modes_equal = true
    for name in names_after
        path = joinpath(original, name)
        state = member_state[name]
        member_identities_equal &=
            (UInt64(stat(path).device), UInt64(stat(path).inode)) == state.identity
        member_bytes_equal &= sha_file(path) == state.sha256
        member_modes_equal &= (UInt(stat(path).mode) & UInt(0o777)) == state.mode
    end
    println("WORKER_DIRECTORY=", DIRECTORY_REPLACEMENT_RELATIVE)
    println("WORKER_ORIGINAL_DIRECTORY_INODE=", original_directory_identity[2])
    println("WORKER_DIRECTORY_INODE_CHANGED=",
            replacement_directory_identity != original_directory_identity ?
                "true" : "false")
    println("WORKER_MEMBER_NAMES_EQUAL=",
            names_after == names ? "true" : "false")
    println("WORKER_MEMBER_IDENTITIES_EQUAL=",
            member_identities_equal ? "true" : "false")
    println("WORKER_MEMBER_BYTES_EQUAL=", member_bytes_equal ? "true" : "false")
    println("WORKER_MEMBER_MODES_EQUAL=", member_modes_equal ? "true" : "false")
    # The retained sibling backup still holds the original directory inode, so
    # the replacement directory inode cannot be a reuse of the original.
    backup_directory_identity =
        (UInt64(stat(backup).device), UInt64(stat(backup).inode))
    println("WORKER_BACKUP_DIRECTORY_INODE=", backup_directory_identity[2])
    _worker_expect_blocked(fixture, "VERIFY", () ->
        fixture._verify_context(context))
    loader_error = fixture.EvaluatorError(
        :BLOCKED, :simulated_downstream_loader_failure,
        "simulated downstream loader failure before publication")
    blocker_files = fixture._blocker_files(loader_error; context=context)
    _worker_expect_publication_rejection(fixture, "PUBLISH", () ->
        fixture._publish_atomic(root, options["--out-dir"], blocker_files;
                                context=context))
    destination = joinpath(root, options["--out-dir"])
    println("WORKER_DESTINATION=", ispath(destination) ? "present" : "absent")
    println("WORKER_STATUS=ok")
    return nothing
end

function _worker_bypass(root::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_bypass_impl, fixture, root)
end

function _worker_bypass_impl(fixture::Module, root::String)
    options = _worker_options(root)
    authorization = fixture._authorize_t13_activation(options)
    plain_options = Dict{String,String}(
        "--root" => root, "--evaluator-config" => RUNTIME_V4_EVALUATOR_CONFIG)
    plain = try
        fixture._preflight_config_context(plain_options)
        :ok
    catch error
        error isa fixture.EvaluatorError ? error.reason : :unexpected
    end
    config_path = joinpath(root, RUNTIME_V4_EVALUATOR_CONFIG)
    info = stat(config_path)
    config_snapshot = fixture._Snapshot(
        config_path, Vector{UInt8}(read(config_path)), fixture._RUNTIME_V4_CONFIG_SHA256,
        UInt64(info.device), UInt64(info.inode), UInt64(info.nlink))
    activated = try
        fixture._collect_authority_snapshots(
            root, RUNTIME_V4_EVALUATOR_CONFIG, :runtime_v4, config_snapshot,
            Dict{Tuple{UInt64,UInt64},String}(); activation=authorization)
        :ok
    catch error
        error isa fixture.EvaluatorError ? error.reason : :unexpected
    end
    mutated = copy(options)
    mutated["--out-dir"] = "other-output"
    receiving = try
        fixture._preflight_config_context(mutated; activation=authorization)
        :ok
    catch error
        error isa fixture.EvaluatorError ? error.reason : :unexpected
    end
    println("WORKER_STATUS=ok")
    println("WORKER_PLAIN=", plain)
    println("WORKER_ACTIVATED=", activated)
    println("WORKER_RECEIVING=", receiving)
    return nothing
end

function _worker_main(root::String)
    fixture = _load_fixture_module(root)
    return Base.invokelatest(_worker_main_impl, fixture, root)
end

function _worker_main_impl(fixture::Module, root::String)
    document = TOML.parsefile(joinpath(root, "main-args.toml"))
    arguments = String[String(entry) for entry in document["arguments"]]
    status = try
        fixture.main(arguments)
    catch error
        println("WORKER_STATUS=unexpected")
        println("WORKER_MESSAGE=", sprint(showerror, error))
        exit(2)
    end
    println("WORKER_STATUS=ok")
    println("WORKER_MAIN_EXIT=", status)
    return nothing
end

function _worker_main_cli(arguments::Vector{String})
    op = arguments[2]
    root = arguments[3]
    if op == "capture"
        _worker_capture(root)
        println("WORKER_STATUS=ok")
    elseif op == "authorize"
        _worker_authorize(root, nothing)
    elseif op == "authorize-foreign"
        _worker_authorize(root, arguments[4])
    elseif op == "bypass"
        _worker_bypass(root)
    elseif op == "toctou"
        _worker_toctou(root, arguments[4])
    elseif op == "context-boundary"
        _worker_context_boundary(root)
    elseif op == "replace-source"
        _worker_replace_source(root)
    elseif op == "replace-directory"
        _worker_replace_directory(root)
    elseif op == "control-report"
        _worker_control_report(root)
    elseif op == "control-preauthorize"
        _worker_control_preauthorize(root, arguments[4])
    elseif op == "control-replacement"
        _worker_control_replacement(root, arguments[4])
    elseif op == "publication-hooks"
        _worker_publication_hooks(root, arguments[4])
    elseif op == "env-probe"
        println("WORKER_ENV=", get(ENV, "STMFIT_HARNESS_PROBE", "unset"))
    elseif op == "main"
        _worker_main(root)
    else
        println("WORKER_STATUS=unknown-op")
        exit(2)
    end
    return nothing
end

if !isempty(ARGS) && first(ARGS) == "--worker"
    _worker_main_cli(copy(ARGS))
    exit(0)
end

# ---------------------------------------------------------------------------
# Driver helpers
# ---------------------------------------------------------------------------

function activation_fixture()
    repo = repo_root()
    scratch = get(ENV, "TMPDIR", "")
    isempty(scratch) && error("TMPDIR must point at the fixture scratch")
    scratch = realpath(scratch)
    root = realpath(mktempdir(scratch; prefix="t13-activation-root-",
                              cleanup=false))
    for relative in E._T13_SOURCE_CLOSURE
        source = joinpath(repo, relative)
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(source, destination; force=true)
        chmod(destination, UInt(stat(source).mode) & UInt(0o777))
    end
    for relative in ("Project.toml", "Manifest.toml", RUNTIME_V4_EVALUATOR_CONFIG)
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(joinpath(repo, relative), destination; force=true)
    end
    for relative in vcat(String.(collect(E._T13_FIXED_PROVENANCE_FILES)),
                         [script.path for script in E._T13_PRODUCER_SCRIPTS])
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(joinpath(repo, relative), destination; force=true)
    end
    for stem in ("gate4_expected_boulder_postimage", "gate4_transition_receipt",
                 "gate4_post_review", "gate4_oracle_final_review")
        relative = E._t13_authority_binding(stem).path
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(joinpath(repo, relative), destination; force=true)
        mode = E._runtime_v4_expected_mode(destination)
        mode === nothing || chmod(destination, mode)
    end
    for (relative, mode) in ((E._T13_GATE5_PUBLICATION_PATH, 0o444),
                             (E._T13_GATE2_ORACLE_REVIEW_PATH, 0o444))
        destination = joinpath(root, relative)
        mkpath(dirname(destination))
        cp(joinpath(repo, relative), destination; force=true)
        chmod(destination, mode)
    end
    # Two genuine disposable roots: the execution root R (above) and the fake
    # Main-control root M (below).  Neither is the real checkout and neither is
    # an operational ledger.  The single live control file lives only under M;
    # R carries a non-consumed decoy at the same relative name to prove the
    # activation branch never falls back to R for control bytes.
    main = realpath(mktempdir(scratch; prefix="t13-main-control-root-",
                              cleanup=false))
    main_parent = joinpath(main, ".omo")
    mkpath(main_parent)
    boulder = joinpath(main_parent, "boulder.json")
    write(boulder, "reviewed Main-control postimage (synthetic fixture)\n")
    chmod(boulder, CONTROL_MODE)
    decoy = joinpath(root, BOULDER_RELATIVE)
    mkpath(dirname(decoy))
    write(decoy, "R-local non-consumed Boulder decoy (synthetic fixture)\n")
    chmod(decoy, 0o644)
    write(joinpath(root, "main-root.txt"), main * "\n")
    records = joinpath(root, "records")
    mkpath(records)
    source_review = joinpath(records, "source-review.json")
    transition_review = joinpath(records, "transition-review.json")
    publication = joinpath(records, "transition-publication.json")
    write(source_review, "ACT-2 source review record (synthetic, non-authoritative)\n")
    write(transition_review, "ACT-3 transition review record (synthetic, non-authoritative)\n")
    write(publication, "transition publication record (synthetic, non-authoritative)\n")
    inputs_dir = joinpath(root, "inputs")
    mkpath(inputs_dir)
    input_paths = Dict{String,String}(
        "--features" => joinpath(inputs_dir, "features.tsv"),
        "--candidate-config" => joinpath(inputs_dir, "candidate.toml"),
        "--model-config" => joinpath(inputs_dir, "model.toml"),
        "--forward-receipt" => joinpath(inputs_dir, "forward.toml"),
        "--backward-receipt" => joinpath(inputs_dir, "backward.toml"),
    )
    for (flag, path) in input_paths
        write(path, flag * " synthetic fixture\n")
    end
    directory_paths = Dict{String,String}(
        "--universe-dir" => joinpath(root, "universe"),
        "--edge-dir" => joinpath(root, "edge"),
        "--admission-dir" => joinpath(root, "admission"),
    )
    for (flag, path) in directory_paths
        mkpath(path)
        write(joinpath(path, "member-one.tsv"), flag * " one\n")
        write(joinpath(path, "member-two.tsv"), flag * " two\n")
    end
    data_dir = joinpath(root, "data")
    mkpath(data_dir)
    feature_table = input_paths["--features"]
    write(feature_table, "feature synthetic fixture\n")
    forward_patch = joinpath(inputs_dir, "forward_patch.tsv")
    backward_patch = joinpath(inputs_dir, "backward_patch.tsv")
    write(forward_patch, "forward patch synthetic fixture\n")
    write(backward_patch, "backward patch synthetic fixture\n")
    features_hash = sha_file(feature_table)
    forward_hash = sha_file(forward_patch)
    backward_hash = sha_file(backward_patch)
    for (orientation, output, script, script_hash) in (
            ("forward", forward_patch, E._T13_PRODUCER_SCRIPTS[1].path,
             E._T13_PRODUCER_SCRIPTS[1].sha256),
            ("backward", backward_patch, E._T13_PRODUCER_SCRIPTS[2].path,
             E._T13_PRODUCER_SCRIPTS[2].sha256))
        write_toml(input_paths["--" * orientation * "-receipt"],
                   Dict{String,Any}("root" => root, "cwd" => root,
                                    "features_path" => feature_table,
                                    "data_path" => data_dir,
                                    "output_path" => output,
                                    "source_path" => script,
                                    "source_sha256" => script_hash,
                                    "feature_sha256" => features_hash))
    end
    forward_receipt_hash = sha_file(input_paths["--forward-receipt"])
    backward_receipt_hash = sha_file(input_paths["--backward-receipt"])
    write_toml(joinpath(directory_paths["--universe-dir"], "receipt.toml"),
               Dict{String,Any}("feature_sha256" => features_hash))
    write_toml(joinpath(directory_paths["--edge-dir"], "receipt.toml"),
               Dict{String,Any}("feature_sha256" => features_hash,
                                "forward_patch_sha256" => forward_hash,
                                "backward_patch_sha256" => backward_hash,
                                "forward_patch_receipt_sha256" => forward_receipt_hash,
                                "backward_patch_receipt_sha256" => backward_receipt_hash))
    return (repo=repo, root=root, main=main, main_parent=main_parent,
            input_paths=input_paths,
            directory_paths=directory_paths, records=records,
            source_review=source_review, transition_review=transition_review,
            publication=publication, boulder=boulder, decoy=decoy,
            data_dir=data_dir)
end

function finalize_activation_fixture(fixture)
    root = fixture.root
    capture = run_worker(root, "capture")
    capture.success || error("worker capture failed: $(capture.stderr)")
    launch = TOML.parsefile(joinpath(root, "captured-launch.toml"))
    destination = joinpath(root, "output")
    receipt_path = joinpath(fixture.records, "receipt.toml")
    spec_path = joinpath(fixture.records, "execution-spec.toml")
    options = Dict{String,String}(
        "--root" => root,
        "--evaluator-config" => RUNTIME_V4_EVALUATOR_CONFIG,
        "--features" => relpath(fixture.input_paths["--features"], root),
        "--candidate-config" => relpath(fixture.input_paths["--candidate-config"], root),
        "--model-config" => relpath(fixture.input_paths["--model-config"], root),
        "--universe-dir" => relpath(fixture.directory_paths["--universe-dir"], root),
        "--edge-dir" => relpath(fixture.directory_paths["--edge-dir"], root),
        "--forward-receipt" => relpath(fixture.input_paths["--forward-receipt"], root),
        "--backward-receipt" => relpath(fixture.input_paths["--backward-receipt"], root),
        "--admission-dir" => relpath(fixture.directory_paths["--admission-dir"], root),
        "--out-dir" => relpath(destination, root),
        "--activation-receipt" => relpath(receipt_path, root),
        "--activation-sha256" => "0"^64,
    )
    source_files = [Dict{String,Any}(
                        "path" => relative,
                        "sha256" => sha_file(joinpath(root, relative)),
                        "mode" => mode_text(joinpath(root, relative)),
                    ) for relative in E._T13_SOURCE_CLOSURE]
    input_entries = [Dict{String,Any}(
                         "path" => relpath(fixture.input_paths[flag], root),
                         "sha256" => sha_file(fixture.input_paths[flag]),
                         "mode" => mode_text(fixture.input_paths[flag]),
                     ) for flag in ("--features", "--candidate-config",
                                    "--model-config", "--forward-receipt",
                                    "--backward-receipt")]
    directory_entries = [Dict{String,Any}(
                             "path" => relpath(fixture.directory_paths[flag], root),
                             "mode" => mode_text(fixture.directory_paths[flag]),
                             "members" => [Dict{String,Any}(
                                               "name" => name,
                                               "sha256" => sha_file(
                                                   joinpath(fixture.directory_paths[flag], name)),
                                               "mode" => mode_text(
                                                   joinpath(fixture.directory_paths[flag], name)),
                                           ) for name in sort(readdir(fixture.directory_paths[flag]))],
                         ) for flag in ("--universe-dir", "--edge-dir",
                                        "--admission-dir")]
    forward_receipt_document = TOML.parsefile(fixture.input_paths["--forward-receipt"])
    backward_receipt_document = TOML.parsefile(fixture.input_paths["--backward-receipt"])
    features_rel = relpath(String(forward_receipt_document["features_path"]), root)
    forward_output_rel = relpath(String(forward_receipt_document["output_path"]), root)
    backward_output_rel = relpath(String(backward_receipt_document["output_path"]), root)
    dependency_paths = sort(unique(vcat(
        String.(collect(E._T13_FIXED_PROVENANCE_FILES)),
        [features_rel, forward_output_rel, backward_output_rel])))
    dependency_entries = [Dict{String,Any}(
                              "path" => relative,
                              "sha256" => sha_file(joinpath(root, relative)),
                              "mode" => mode_text(joinpath(root, relative)),
                          ) for relative in dependency_paths]
    hierarchical = joinpath(root, "test/lib/hierarchical")
    dependency_directories = [Dict{String,Any}(
        "path" => "test/lib/hierarchical",
        "mode" => mode_text(hierarchical),
        "members" => [Dict{String,Any}(
                          "name" => name,
                          "sha256" => sha_file(joinpath(hierarchical, name)),
                          "mode" => mode_text(joinpath(hierarchical, name)),
                      ) for name in sort(readdir(hierarchical))])]
    identity_entries = [Dict{String,Any}("path" => script.path, "kind" => "file")
                        for script in E._T13_PRODUCER_SCRIPTS]
    push!(identity_entries,
          Dict{String,Any}(
              "path" => relpath(String(forward_receipt_document["data_path"]), root),
              "kind" => "directory"))
    dependencies = Dict{String,Any}("files" => dependency_entries,
                                    "directories" => dependency_directories,
                                    "identity_paths" => identity_entries)
    spec = Dict{String,Any}(
        "schema" => E._T13_EXECUTION_SPEC_SCHEMA,
        "scope" => E._T13_SYNTHETIC_SCOPE,
        "root" => Dict{String,Any}("path" => root),
        "control" => Dict{String,Any}(
            "role" => E._T13_CONTROL_ROLE,
            "root" => fixture.main,
            "root_mode" => mode_text(fixture.main),
            "parent_mode" => mode_text(fixture.main_parent),
            "file_mode" => mode_text(fixture.boulder)),
        "entrypoint" => Dict{String,Any}(
            "path" => "test/evaluate_structured_unit_assignment.jl",
            "sha256" => sha_file(joinpath(root,
                                          "test/evaluate_structured_unit_assignment.jl"))),
        "config" => Dict{String,Any}(
            "path" => RUNTIME_V4_EVALUATOR_CONFIG,
            "sha256" => E._RUNTIME_V4_CONFIG_SHA256),
        "project" => Dict{String,Any}(
            "path" => "Project.toml",
            "sha256" => sha_file(joinpath(root, "Project.toml"))),
        "manifest" => Dict{String,Any}(
            "path" => "Manifest.toml",
            "sha256" => sha_file(joinpath(root, "Manifest.toml"))),
        "source" => Dict{String,Any}("files" => source_files),
        "command" => Dict{String,Any}(
            "body" => E._t13_canonical_command(options),
            "environment" => [Dict{String,Any}("name" => ACTIVATION_ENVIRONMENT,
                                               "value" => "1")]),
        "launch" => launch,
        "destination" => Dict{String,Any}(
            "path" => relpath(destination, root),
            "reports" => collect(E._FILES)),
        "inputs" => Dict{String,Any}("files" => input_entries,
                                     "directories" => directory_entries),
        "dependencies" => dependencies,
        "predecessors" => predecessor_table(),
        "permissions" => Dict{String,Any}(
            "real_campaign" => false, "public_producer" => false, "t14" => false,
            "downstream" => false, "grading" => false,
            "runtime_equivalence" => false),
    )
    write_control_toml(spec_path, spec)
    receipt = Dict{String,Any}(
        "schema" => E._T13_ACTIVATION_SCHEMA,
        "scope" => E._T13_SYNTHETIC_SCOPE,
        "spec" => Dict{String,Any}("path" => relpath(spec_path, root),
                                   "sha256" => sha_file(spec_path)),
        "source_review" => Dict{String,Any}(
            "path" => relpath(fixture.source_review, root),
            "sha256" => sha_file(fixture.source_review)),
        "transition_review" => Dict{String,Any}(
            "path" => relpath(fixture.transition_review, root),
            "sha256" => sha_file(fixture.transition_review)),
        "boulder" => Dict{String,Any}(
            "control_root" => fixture.main,
            "preimage_path" => BOULDER_RELATIVE,
            "preimage_sha256" => E._RUNTIME_V4_BOULDER_SHA256,
            "postimage_path" => BOULDER_RELATIVE,
            "postimage_sha256" => sha_file(fixture.boulder)),
        "publication" => Dict{String,Any}(
            "path" => relpath(fixture.publication, root),
            "sha256" => sha_file(fixture.publication)),
    )
    write_control_toml(receipt_path, receipt)
    for path in (fixture.source_review, fixture.transition_review, fixture.publication)
        chmod(path, CONTROL_MODE)
    end
    options["--activation-sha256"] = sha_file(receipt_path)
    write_toml(joinpath(root, "options.toml"), options)
    return merge(fixture, (options=options, receipt_path=receipt_path,
                           spec_path=spec_path, destination=destination,
                           launch=launch))
end

function predecessor_table()
    table = Dict{String,Any}()
    for stem in ("gate4_expected_boulder_postimage", "gate4_transition_receipt",
                 "gate4_post_review", "gate4_oracle_final_review")
        binding = E._t13_authority_binding(stem)
        table[stem * "_path"] = binding.path
        table[stem * "_sha256"] = binding.sha256
    end
    table["gate5_publication_receipt_path"] = E._T13_GATE5_PUBLICATION_PATH
    table["gate5_publication_receipt_sha256"] = E._T13_GATE5_PUBLICATION_SHA256
    table["gate2_oracle_review_path"] = E._T13_GATE2_ORACLE_REVIEW_PATH
    table["gate2_oracle_review_sha256"] = E._T13_GATE2_ORACLE_REVIEW_SHA256
    return table
end

function with_fixture(f::Function; finalize::Bool=true)
    fixture = activation_fixture()
    # Diagnostic retention: fixture roots stay under the bounded TMPDIR for
    # post-run inspection; the marker records the deliberate non-deletion.
    write(joinpath(fixture.root, "T13_RETAINED_DIAGNOSTIC_FIXTURE"),
          "retained diagnostic fixture root; safe to remove with TMPDIR\n")
    return f(finalize ? finalize_activation_fixture(fixture) : fixture)
end

function rewrite_receipt!(mutator::Function, fixture)
    document = TOML.parsefile(fixture.receipt_path)
    mutator(document)
    write_control_toml(fixture.receipt_path, document)
    fixture.options["--activation-sha256"] = sha_file(fixture.receipt_path)
    write_toml(joinpath(fixture.root, "options.toml"), fixture.options)
    return document
end

function rewrite_spec!(mutator::Function, fixture)
    document = TOML.parsefile(fixture.spec_path)
    mutator(document)
    write_control_toml(fixture.spec_path, document)
    receipt = TOML.parsefile(fixture.receipt_path)
    receipt["spec"]["sha256"] = sha_file(fixture.spec_path)
    write_control_toml(fixture.receipt_path, receipt)
    fixture.options["--activation-sha256"] = sha_file(fixture.receipt_path)
    write_toml(joinpath(fixture.root, "options.toml"), fixture.options)
    return document
end

# Malicious M == R control-as-input fixture: the R-local decoy becomes the
# single Main control file and one approved input file (the candidate config)
# is the same pathname, so the dedicated control file would also be captured as
# an ordinary snapshot.  Authorization must reject the overlap.
function malicious_control_as_input!(fixture)
    chmod(fixture.decoy, CONTROL_MODE)
    decoy_sha = sha_file(fixture.decoy)
    decoy_mode = mode_text(fixture.decoy)
    original_candidate_relative =
        relpath(fixture.input_paths["--candidate-config"], fixture.root)
    fixture.options["--candidate-config"] = BOULDER_RELATIVE
    rewrite_spec!(fixture) do document
        document["control"]["root"] = fixture.root
        document["control"]["root_mode"] = mode_text(fixture.root)
        document["control"]["parent_mode"] =
            mode_text(joinpath(fixture.root, ".omo"))
        document["control"]["file_mode"] = decoy_mode
        for entry in document["inputs"]["files"]
            if entry["path"] == original_candidate_relative
                entry["path"] = BOULDER_RELATIVE
                entry["sha256"] = decoy_sha
                entry["mode"] = decoy_mode
            end
        end
        document["command"]["body"] = E._t13_canonical_command(fixture.options)
    end
    rewrite_receipt!(fixture) do document
        document["boulder"]["control_root"] = fixture.root
        document["boulder"]["postimage_sha256"] = decoy_sha
    end
    return fixture
end

function worker_reason(result)
    for line in split(result.stdout, '\n')
        startswith(line, "WORKER_REASON=") && return String(line[15:end])
    end
    return ""
end

function worker_value(result, prefix::String)
    for line in split(result.stdout, '\n')
        startswith(line, prefix) && return String(line[length(prefix) + 1:end])
    end
    return ""
end

function main_arguments(options::Dict{String,String})
    arguments = String[]
    for flag in E._T13_ACTIVATION_COMMAND_FLAGS
        push!(arguments, flag)
        push!(arguments, options[flag])
    end
    push!(arguments, "--activation-sha256", options["--activation-sha256"])
    return arguments
end

function dummy_report()
    return E.EvaluatorReport(
        :PASS, :ok, E.NamedTuple[], E.NamedTuple[], E.NamedTuple[],
        E.NamedTuple[], E.NamedTuple[], E.NamedTuple[], NamedTuple(),
        "a"^64, "b"^64, "c"^64, Float64[], "d"^64, "e"^64)
end

function dummy_control(root::String)
    main = joinpath(root, "dummy-main-control")
    parent = joinpath(main, ".omo")
    mkpath(parent)
    path = joinpath(parent, "boulder.json")
    isfile(path) || write(path, "dummy Main-control fixture\n")
    chmod(path, CONTROL_MODE)
    bytes = Vector{UInt8}(read(path))
    info = stat(path)
    snapshot = E._Snapshot(path, bytes, bytes2hex(sha256(bytes)),
                           UInt64(info.device), UInt64(info.inode),
                           UInt64(info.nlink))
    return E._MainControl(
        main, snapshot,
        UInt(stat(main).mode) & UInt(0o777),
        UInt(stat(parent).mode) & UInt(0o777),
        UInt(stat(path).mode) & UInt(0o777),
        (UInt64(stat(main).device), UInt64(stat(main).inode)),
        (UInt64(stat(parent).device), UInt64(stat(parent).inode)),
        E._RUNTIME_V4_BOULDER_SHA256, snapshot.sha256)
end

function dummy_bindings(root::String)
    modes = Dict{String,UInt}()
    identities = Dict{String,Tuple{UInt64,UInt64}}()
    return E._ActivationBindings(
        E._T13_SYNTHETIC_SCOPE,
        "receipt.toml", "1"^64, "spec.toml", "2"^64,
        "source-review.json", "3"^64, "transition-review.json", "4"^64,
        "publication.json", "5"^64,
        BOULDER_RELATIVE, E._RUNTIME_V4_BOULDER_SHA256,
        joinpath(root, "dummy-main-control", ".omo", "boulder.json"), "6"^64,
        "test/evaluate_structured_unit_assignment.jl", "7"^64,
        "output", modes, root, String[], identities,
        dummy_control(root))
end

# ---------------------------------------------------------------------------
# Driver tests
# ---------------------------------------------------------------------------

function _harness_fail(label::String, result)
    println("HARNESS_CHECK=failed")
    println("HARNESS_LABEL=", label)
    println("HARNESS_EXIT=", result.exitcode)
    println("HARNESS_STDOUT=", result.stdout)
    println("HARNESS_STDERR=", result.stderr)
    exit(1)
end

function _harness_check()
    scratch = get(ENV, "TMPDIR", "")
    isempty(scratch) && error("TMPDIR is required for the harness check")
    scratch = realpath(scratch)
    root = realpath(mktempdir(scratch; prefix="t13-harness-root-",
                              cleanup=false))
    cp(joinpath(repo_root(), "Project.toml"), joinpath(root, "Project.toml");
       force=true)
    cp(joinpath(repo_root(), "Manifest.toml"), joinpath(root, "Manifest.toml");
       force=true)
    capture = run_worker(root, "capture")
    (capture.success && occursin("WORKER_STATUS=ok", capture.stdout)) ||
        _harness_fail("capture probe", capture)
    unknown = run_worker(root, "unknown-op")
    (!unknown.success && unknown.exitcode != 0 &&
     occursin("WORKER_STATUS=unknown-op", unknown.stdout)) ||
        _harness_fail("nonzero probe", unknown)
    probe = run_worker(root, "env-probe";
                       extra_env=Dict("STMFIT_HARNESS_PROBE" => "override"))
    (probe.success && occursin("WORKER_ENV=override", probe.stdout)) ||
        _harness_fail("environment override probe", probe)
    println("HARNESS_CHECK=ok")
    println("HARNESS_CAPTURE_LOG=", capture.log_dir)
    println("HARNESS_NONZERO_LOG=", unknown.log_dir)
    println("HARNESS_ENV_LOG=", probe.log_dir)
    return nothing
end


function _fixture_tree_members(root::String)
    found = String[]
    for (directory, _, files) in walkdir(root)
        for name in files
            push!(found, relpath(joinpath(directory, name), root))
        end
    end
    return sort(found)
end

function _copy_fixture_tree(source::String, destination::String)
    isdir(source) || error("fixture source is absent: $source")
    isdir(destination) && error("fixture destination exists: $destination")
    mkpath(destination)
    for name in sort(readdir(source))
        src = joinpath(source, name)
        dst = joinpath(destination, name)
        if islink(src)
            error("fixture copy refuses symlink: $name")
        elseif isdir(src)
            _copy_fixture_tree(src, dst)
        elseif isfile(src)
            stat(src).nlink == 1 || error("fixture copy refuses hardlink: $name")
            cp(src, dst; force=true)
            chmod(dst, UInt(stat(src).mode) & UInt(0o777))
        else
            error("fixture copy refuses special entry: $name")
        end
    end
    return destination
end

function _verify_fixture_copy(source::String, destination::String)
    source_members = _fixture_tree_members(source)
    destination_members = _fixture_tree_members(destination)
    source_members == destination_members ||
        error("fixture copy membership differs")
    for relative in source_members
        sha_file(joinpath(source, relative)) ==
            sha_file(joinpath(destination, relative)) ||
            error("fixture copy bytes differ: $relative")
        (UInt(stat(joinpath(source, relative)).mode) & UInt(0o777)) ==
            (UInt(stat(joinpath(destination, relative)).mode) & UInt(0o777)) ||
            error("fixture copy mode differs: $relative")
    end
    return destination_members
end

function _smoke_fail(label::String, result)
    println("DISPATCH_SMOKE=failed")
    println("SMOKE_CASE=", label)
    println("SMOKE_EXIT=", result.exitcode)
    println("SMOKE_STDOUT=", result.stdout)
    println("SMOKE_STDERR=", result.stderr)
    println("SMOKE_LOG=", result.log_dir)
    exit(1)
end

function _dispatch_smoke()
    scratch = get(ENV, "TMPDIR", "")
    isempty(scratch) && error("TMPDIR is required for the dispatch smoke")
    scratch = realpath(scratch)
    fixture = finalize_activation_fixture(activation_fixture())
    authorize = run_worker(fixture.root, "authorize")
    (authorize.success && occursin("WORKER_STATUS=ok", authorize.stdout) &&
     !ispath(fixture.destination)) || _smoke_fail("authorize", authorize)
    foreign = realpath(mktempdir(scratch; prefix="t13-smoke-foreign-",
                                 cleanup=false))
    rm(foreign)
    _copy_fixture_tree(fixture.root, foreign)
    foreign_members = _verify_fixture_copy(fixture.root, foreign)
    foreign_result = run_worker(fixture.root, "authorize-foreign";
                                foreign=foreign)
    (foreign_result.success &&
     occursin("WORKER_STATUS=blocked", foreign_result.stdout) &&
     worker_reason(foreign_result) == "activation_source_mismatch" &&
     occursin("entrypoint", foreign_result.stdout) &&
     !ispath(fixture.destination)) ||
        _smoke_fail("authorize-foreign", foreign_result)
    println("DISPATCH_SMOKE=ok")
    println("SMOKE_FOREIGN_MEMBERS=", length(foreign_members))
    println("SMOKE_FIXTURE_ROOT=", fixture.root)
    println("SMOKE_FOREIGN_ROOT=", foreign)
    println("SMOKE_AUTHORIZE_LOG=", authorize.log_dir)
    println("SMOKE_FOREIGN_LOG=", foreign_result.log_dir)
    return nothing
end

# Focused assembled-context/source/directory context testsets.  The exact
# bodies are factored into this helper so the default full suite and the
# explicit test-only `--context-smoke` branch execute the same code once.
# Main-control replacement/mutation coverage lives in
# `_mc_remediation_smoke_tests` below.
function _context_smoke_tests()
    # Unit-level shared-production-context boundary.  The worker obtains a
    # genuine fixture authorization and calls the actual
    # `_assemble_production_context` with an explicit non-authoritative
    # ordinary authority fixture (only the approved runtime-v4 config snapshot
    # collected through the real `_snapshot` helper).  This does not exercise
    # the full ordinary authority collector or the heavy positive T8/T11/T12
    # loaders, and it is not a main/end-to-end integration claim.
    first_result = @testset "assembled production context retains activation and rejects record drift" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "context-boundary")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_RECEIVING_VALID=") == "ok"
            @test worker_value(result, "WORKER_RECEIVING_ROOT_REASON=") ==
                  "activation_binding_mismatch"
            @test worker_value(result, "WORKER_RECEIVING_COMMAND_REASON=") ==
                  "activation_command_mismatch"
            @test worker_value(result, "WORKER_CONTEXT_ACTIVATION=") == "retained"
            @test worker_value(result, "WORKER_CONTEXT_MODES=") == "complete"
            @test worker_value(result, "WORKER_CONTEXT_IDENTITIES=") == "complete"
            @test worker_value(result, "WORKER_CONTEXT_SNAPSHOTS=") == "complete"
            @test worker_value(result, "WORKER_CONTEXT_INVENTORIES=") == "complete"
            @test worker_value(result, "WORKER_RECORD_MODE_RESTORED=") == "true"
            @test worker_value(result, "WORKER_VERIFY_REASON=") ==
                  "authority_snapshot_changed"
            @test worker_value(result, "WORKER_BLOCKER_SCHEMA=") ==
                  E._T13_ACTIVATED_BLOCKER_SCHEMA
            @test worker_value(result, "WORKER_BLOCKER_SCOPE=") ==
                  E._T13_SYNTHETIC_SCOPE
            @test worker_value(result, "WORKER_BLOCKER_RECEIPT_SHA=") ==
                  fixture.options["--activation-sha256"]
            @test worker_value(result, "WORKER_BLOCKER_REASON=") ==
                  "simulated_downstream_loader_failure"
            @test worker_value(result, "WORKER_PUBLISH_STATE=") == "not_committed"
            @test worker_value(result, "WORKER_PUBLISH_REASON=") ==
                  "authority_snapshot_changed"
            @test worker_value(result, "WORKER_DESTINATION=") == "absent"
            @test !ispath(fixture.destination)
        end
    end

    # Actual filesystem replacement of an activation source file: same bytes
    # and mode, new inode, original inode retained at a backup path.  The
    # original authorization identity must not be rebased by a fresh ordinary
    # snapshot.  The Main-control replacement cases live below.
    second_result = @testset "actual source pathname replacement is not rebased" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "replace-source")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            relative = worker_value(result, "WORKER_RELATIVE=")
            @test relative in E._T13_SOURCE_CLOSURE
            @test worker_value(result, "WORKER_INODE_CHANGED=") == "true"
            @test worker_value(result, "WORKER_BYTES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_MODE_EQUAL=") == "true"
            @test worker_value(result, "WORKER_AUTH_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_INODE=")
            @test worker_value(result, "WORKER_BACKUP_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_INODE=")
            @test worker_value(result, "WORKER_BACKUP_BYTES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_BACKUP_MODE_EQUAL=") == "true"
            @test worker_value(result, "WORKER_REPLACEMENT_INODE=") !=
                  worker_value(result, "WORKER_ORIGINAL_INODE=")
            @test worker_value(result, "WORKER_FRESH_INODE=") ==
                  worker_value(result, "WORKER_REPLACEMENT_INODE=")
            @test worker_value(result, "WORKER_MERGE_REASON=") ==
                  "activation_binding_mismatch"
            @test worker_value(result, "WORKER_ASSEMBLER_REASON=") ==
                  "activation_binding_mismatch"
            @test worker_value(result, "WORKER_DESTINATION=") == "absent"
            @test !ispath(fixture.destination)
        end
    end

    # Actual directory replacement: an approved nonempty input directory is
    # renamed to a retained sibling and recreated at the original pathname with
    # the original member files moved in, so member dev/inode/bytes/names/modes
    # stay identical while the directory inode changes.  The captured directory
    # identity is the boundary that must reject verification and publication.
    third_result = @testset "actual approved input directory replacement is rejected" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "replace-directory")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_CONTEXT=") == "assembled"
            @test worker_value(result, "WORKER_DIRECTORY=") ==
                  DIRECTORY_REPLACEMENT_RELATIVE
            @test worker_value(result, "WORKER_DIRECTORY_INODE_CHANGED=") == "true"
            @test worker_value(result, "WORKER_BACKUP_DIRECTORY_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_DIRECTORY_INODE=")
            @test worker_value(result, "WORKER_MEMBER_NAMES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_MEMBER_IDENTITIES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_MEMBER_BYTES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_MEMBER_MODES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_VERIFY_REASON=") ==
                  "activation_snapshot_changed"
            @test worker_value(result, "WORKER_PUBLISH_STATE=") == "not_committed"
            @test worker_value(result, "WORKER_PUBLISH_REASON=") ==
                  "activation_snapshot_changed"
            @test worker_value(result, "WORKER_DESTINATION=") == "absent"
            @test !ispath(fixture.destination)
        end
    end
    return (first_result, second_result, third_result)
end

function _activation_test_counts(result)
    counts = Test.get_test_counts(result)
    return (passes=counts.passes + counts.cumulative_passes,
            fails=counts.fails + counts.cumulative_fails,
            errors=counts.errors + counts.cumulative_errors,
            broken=counts.broken + counts.cumulative_broken)
end

function _print_context_smoke_summary(smoke_results)
    total_passes = 0
    total_fails = 0
    total_errors = 0
    total_broken = 0
    for smoke_result in smoke_results
        counts = _activation_test_counts(smoke_result)
        total_passes += counts.passes
        total_fails += counts.fails
        total_errors += counts.errors
        total_broken += counts.broken
        println("CONTEXT_SMOKE_TESTSET=", smoke_result.description,
                " passes=", counts.passes,
                " fails=", counts.fails,
                " errors=", counts.errors,
                " broken=", counts.broken)
    end
    println("CONTEXT_SMOKE_TOTAL=passes=", total_passes,
            " fails=", total_fails, " errors=", total_errors,
            " broken=", total_broken)
    return nothing
end

# ---------------------------------------------------------------------------
# MC-2 remediation smoke helper
# ---------------------------------------------------------------------------
#
# The exact repaired test bodies for the MC-2 failure groups and the new
# production-gap checks are factored here so the parent can run the seven
# targeted groups through `--mc-remediation-smoke` and the default full suite
# executes the same code once.  Group 7's source checks assert the actual
# post-read/stage helpers and the precise wiring/order in the frozen evaluator
# source; they do not invent assertion counts.

function _evaluator_source_text()
    return read(joinpath(repo_root(), "test",
                         "evaluate_structured_unit_assignment.jl"), String)
end

function _source_function_body(text::String, signature::String)
    start = findfirst(signature, text)
    start === nothing && error("source function is absent: $signature")
    stop = findnext("\nend\n", text, last(start))
    stop === nothing && error("source function is unterminated: $signature")
    return text[first(start):last(stop)]
end

function _source_ordered(body::String, markers::Vector{String})
    positions = Int[]
    cursor = firstindex(body)
    for marker in markers
        index = findnext(marker, body, cursor)
        index === nothing && return nothing
        push!(positions, first(index))
        cursor = nextind(body, last(index))
    end
    return issorted(positions) ? positions : nothing
end

function _mc_remediation_smoke_tests()
    # Group 1: neutral parser-only pairing with a complete literal argv, the
    # firewall's continued rejection of a value containing a forbidden concept,
    # and the full CLI-through-main route from the neutral canonical scratch
    # (the stage TMPDIR pathname itself is not firewall-clean).
    cli_result = @testset "MC remediation: neutral CLI pairing, firewall value, and CLI-through-main" begin
        neutral = String[
            "--root", "/neutral/root",
            "--evaluator-config", "config/unit_assignment_structured_evaluator_runtime_v4.toml",
            "--features", "inputs/features.tsv",
            "--candidate-config", "inputs/candidate.toml",
            "--model-config", "inputs/model.toml",
            "--universe-dir", "universe",
            "--edge-dir", "edge",
            "--forward-receipt", "inputs/forward.toml",
            "--backward-receipt", "inputs/backward.toml",
            "--admission-dir", "admission",
            "--out-dir", "output",
            "--activation-receipt", "records/receipt.toml",
            "--activation-sha256", "0"^64,
        ]
        @test E.parse_cli(neutral)["--activation-receipt"] == "records/receipt.toml"
        digest_index = findlast(==("--activation-sha256"), neutral)
        partial = copy(neutral)
        deleteat!(partial, digest_index + 1)
        deleteat!(partial, digest_index)
        error = try E.parse_cli(partial); nothing catch caught; caught end
        @test error isa E.EvaluatorError && error.reason == :cli_error
        @test occursin("must be supplied together", error.message)
        receipt_index = findlast(==("--activation-receipt"), neutral)
        digest_only = copy(neutral)
        deleteat!(digest_only, receipt_index + 1)
        deleteat!(digest_only, receipt_index)
        error = try E.parse_cli(digest_only); nothing catch caught; caught end
        @test error isa E.EvaluatorError && error.reason == :cli_error
        @test occursin("must be supplied together", error.message)
        duplicate = vcat(neutral, ["--activation-sha256", "0"^64])
        error = try E.parse_cli(duplicate); nothing catch caught; caught end
        @test error isa E.EvaluatorError && error.reason == :cli_error
        unknown = vcat(partial, ["--activation-sha", "0"^64])
        error = try E.parse_cli(unknown); nothing catch caught; caught end
        @test error isa E.EvaluatorError && error.reason == :cli_error
        # Separate firewall rejection: the parser values above are neutral, but
        # a value containing a forbidden concept is still rejected.  The root
        # is not exempted from the firewall.
        forbidden = copy(neutral)
        forbidden[findfirst(==("--features"), forbidden) + 1] = "control-value.tsv"
        error = try E.parse_cli(forbidden); nothing catch caught; caught end
        @test error isa E.EvaluatorError && error.reason == :cli_error
        @test occursin("forbidden", error.message)
        with_neutral_cli_fixture() do fixture
            @test startswith(fixture.root, realpath(neutral_cli_scratch()))
            arguments = main_arguments(fixture.options)
            digest_index = findlast(==("--activation-sha256"), arguments)
            arguments[digest_index + 1] = "0"^64
            write_toml(joinpath(fixture.root, "main-args.toml"),
                       Dict{String,Any}("arguments" => arguments))
            result = run_worker(fixture.root, "main")
            @test result.success
            @test occursin("WORKER_MAIN_EXIT=2", result.stdout)
            # The intended downstream reason, not merely a firewall rejection:
            # the receipt digest does not match the reviewed receipt.
            @test occursin("authority_hash_mismatch", result.stderr)
            @test occursin("activation receipt", result.stderr)
            @test !occursin("forbidden_cli", result.stderr)
            @test !ispath(fixture.destination)
        end
    end

    # Group 2: the real two-root M/R positive, the tagged manifest, and the
    # missing-M/wrong-bytes/decoy non-consumption rejections.
    two_root_result = @testset "MC remediation: real M/R positive and missing-M decoy" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "control-report")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_CONTROL_ROOT=") == fixture.main
            @test worker_value(result, "WORKER_CONTROL_ROOT_IS_MAIN=") == "true"
            @test worker_value(result, "WORKER_CONTROL_FILE_FIXED=") == "true"
            @test worker_value(result, "WORKER_CONTROL_FILE=") ==
                  joinpath(fixture.main, BOULDER_RELATIVE)
            @test worker_value(result, "WORKER_CONTROL_SHA=") ==
                  sha_file(fixture.boulder)
            @test worker_value(result, "WORKER_CONTROL_SHA=") !=
                  worker_value(result, "WORKER_DECOY_SHA=")
            @test worker_value(result, "WORKER_CONTROL_SNAPSHOT_ORDINARY=") ==
                  "false"
            @test worker_value(result, "WORKER_MANIFEST_TAGGED=") == "1"
            @test worker_value(result, "WORKER_MANIFEST_ESCAPES=") == "false"
            @test worker_value(result, "WORKER_MANIFEST_FILES_UNDER_ROOT=") ==
                  "true"
            record = worker_value(result, "WORKER_MANIFEST_CONTROL_RECORD=")
            @test occursin("role=main_boulder", record)
            @test occursin("root=$(fixture.main)", record)
            @test occursin("file=.omo/boulder.json", record)
            @test occursin("sha256=$(sha_file(fixture.boulder))", record)
            @test occursin("mode=0444", record)
            @test !occursin("../", record)
            @test !ispath(fixture.destination)
        end
        # M == R is allowed only when that root is the reviewed Main.  The
        # control record must still be unique and the control path must not
        # appear as an ordinary snapshot alias.
        with_fixture() do fixture
            chmod(fixture.decoy, CONTROL_MODE)
            rewrite_spec!(fixture) do document
                document["control"]["root"] = fixture.root
                document["control"]["root_mode"] = mode_text(fixture.root)
                document["control"]["parent_mode"] =
                    mode_text(joinpath(fixture.root, ".omo"))
                document["control"]["file_mode"] = mode_text(fixture.decoy)
            end
            rewrite_receipt!(fixture) do document
                document["boulder"]["control_root"] = fixture.root
                document["boulder"]["postimage_sha256"] = sha_file(fixture.decoy)
            end
            result = run_worker(fixture.root, "control-report")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_CONTROL_ROOT=") == fixture.root
            @test worker_value(result, "WORKER_CONTROL_FILE=") ==
                  joinpath(fixture.root, BOULDER_RELATIVE)
            @test worker_value(result, "WORKER_CONTROL_SHA=") ==
                  sha_file(fixture.decoy)
            @test worker_value(result, "WORKER_CONTROL_SNAPSHOT_ORDINARY=") ==
                  "false"
            @test worker_value(result, "WORKER_MANIFEST_TAGGED=") == "1"
            @test worker_value(result, "WORKER_MANIFEST_ESCAPES=") == "false"
            @test worker_value(result, "WORKER_MANIFEST_FILES_UNDER_ROOT=") ==
                  "true"
            @test !ispath(fixture.destination)
        end
        for (mutation, expected) in (
                ("missing-file", "activation_control_mismatch"),
                ("missing-root", "activation_control_mismatch"),
                ("wrong-hash", "authority_hash_mismatch"),
                ("unreadable", "authority_path_mismatch"),
                ("mode-file", "authority_path_mismatch"),
                ("mode-parent", "authority_path_mismatch"),
                ("mode-root", "authority_path_mismatch"))
            with_fixture() do fixture
                result = run_worker(fixture.root, "control-preauthorize";
                                    extra=mutation)
                @test result.success
                @test occursin("WORKER_STATUS=ok", result.stdout)
                @test worker_value(result, "WORKER_MUTATION=") == mutation
                @test worker_value(result, "WORKER_AUTHORIZE_REASON=") == expected
                @test !ispath(fixture.destination)
            end
        end
        with_fixture() do fixture
            result = run_worker(fixture.root, "control-preauthorize";
                                extra="r-decoy")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_AUTHORIZE_REASON=") == "ok"
            @test worker_value(result, "WORKER_CONTROL_SHA=") ==
                  sha_file(fixture.boulder)
            @test worker_value(result, "WORKER_DECOY_SHA=") !=
                  worker_value(result, "WORKER_CONTROL_SHA=")
        end
        with_fixture() do fixture
            result = run_worker(fixture.root, "control-replacement";
                                extra="r-decoy-mutate")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_VERIFY=") == "ok"
            @test worker_value(result, "WORKER_CONTROL_SHA=") ==
                  sha_file(fixture.boulder)
            @test worker_value(result, "WORKER_DECOY_SHA=") !=
                  worker_value(result, "WORKER_CONTROL_SHA=")
        end
    end

    # Group 3: replacements performed before the authorization capture are
    # valid initial acquisitions with the new captured identities; the same
    # replacements after the capture stay negative cases.
    initial_result = @testset "MC remediation: initial acquisition vs post-capture identity replacement" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "control-preauthorize";
                                extra="parent-replace-initial")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_MUTATION=") ==
                  "parent-replace-initial"
            @test worker_value(result, "WORKER_AUTHORIZE_REASON=") == "ok"
            @test worker_value(result, "WORKER_PARENT_INODE_CHANGED=") == "true"
            @test worker_value(result, "WORKER_FILE_INODE_PRESERVED=") == "true"
            @test worker_value(result, "WORKER_BACKUP_PARENT_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_PARENT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_PARENT_INODE=") ==
                  worker_value(result, "WORKER_CURRENT_PARENT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_PARENT_INODE=") !=
                  worker_value(result, "WORKER_ORIGINAL_PARENT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_ROOT_INODE=") ==
                  worker_value(result, "WORKER_CURRENT_ROOT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_FILE_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_FILE_INODE=")
            @test !ispath(fixture.destination)
        end
        with_fixture() do fixture
            result = run_worker(fixture.root, "control-preauthorize";
                                extra="root-replace-initial")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_MUTATION=") ==
                  "root-replace-initial"
            @test worker_value(result, "WORKER_AUTHORIZE_REASON=") == "ok"
            @test worker_value(result, "WORKER_ROOT_INODE_CHANGED=") == "true"
            @test worker_value(result, "WORKER_PARENT_INODE_PRESERVED=") == "true"
            @test worker_value(result, "WORKER_BACKUP_ROOT_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_ROOT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_ROOT_INODE=") ==
                  worker_value(result, "WORKER_CURRENT_ROOT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_ROOT_INODE=") !=
                  worker_value(result, "WORKER_ORIGINAL_ROOT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_PARENT_INODE=") ==
                  worker_value(result, "WORKER_CURRENT_PARENT_INODE=")
            @test worker_value(result, "WORKER_CAPTURED_FILE_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_FILE_INODE=")
            @test !ispath(fixture.destination)
        end
        for mutation in ("parent-inode", "root-inode")
            with_fixture() do fixture
                result = run_worker(fixture.root, "control-replacement";
                                    extra=mutation)
                @test result.success
                @test occursin("WORKER_STATUS=ok", result.stdout)
                @test worker_value(result, "WORKER_CONTEXT=") == "assembled"
                @test worker_value(result, "WORKER_VERIFY_REASON=") ==
                      "activation_snapshot_changed"
                @test worker_value(result, "WORKER_PUBLISH_STATE=") ==
                      "not_committed"
                @test worker_value(result, "WORKER_PUBLISH_REASON=") ==
                      "activation_snapshot_changed"
                @test worker_value(result, "WORKER_DESTINATION=") == "absent"
                @test !ispath(fixture.destination)
                if mutation == "parent-inode"
                    @test worker_value(result, "WORKER_PARENT_INODE_CHANGED=") ==
                          "true"
                    @test worker_value(result, "WORKER_FILE_INODE_PRESERVED=") ==
                          "true"
                else
                    @test worker_value(result, "WORKER_ROOT_INODE_CHANGED=") ==
                          "true"
                    @test worker_value(result, "WORKER_PARENT_INODE_PRESERVED=") ==
                          "true"
                end
            end
        end
    end

    # Group 4: the retained dedicated Main-control object, the post-read
    # file-mode/byte rechecks, and the rejection of a same-bytes new-inode
    # replacement (fresh observation is logged but never adopted).
    inode_result = @testset "MC remediation: Main file inode keeps the dedicated control object" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "control-replacement";
                                extra="file-inode")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_CONTEXT=") == "assembled"
            @test worker_value(result, "WORKER_FILE_INODE_CHANGED=") == "true"
            @test worker_value(result, "WORKER_BYTES_EQUAL=") == "true"
            @test worker_value(result, "WORKER_MODE_EQUAL=") == "true"
            @test worker_value(result, "WORKER_AUTH_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_INODE=")
            @test worker_value(result, "WORKER_BACKUP_INODE=") ==
                  worker_value(result, "WORKER_ORIGINAL_INODE=")
            @test worker_value(result, "WORKER_REPLACEMENT_INODE=") !=
                  worker_value(result, "WORKER_ORIGINAL_INODE=")
            @test worker_value(result, "WORKER_FRESH_INODE=") ==
                  worker_value(result, "WORKER_REPLACEMENT_INODE=")
            @test worker_value(result, "WORKER_FRESH_ADOPTED=") == "false"
            @test worker_value(result, "WORKER_CONTEXT_RETAINS_CONTROL=") == "true"
            @test worker_value(result, "WORKER_ASSEMBLER_REASON=") ==
                  "authority_snapshot_changed"
            @test worker_value(result, "WORKER_VERIFY_REASON=") ==
                  "authority_snapshot_changed"
            @test worker_value(result, "WORKER_PUBLISH_STATE=") == "not_committed"
            @test worker_value(result, "WORKER_PUBLISH_REASON=") ==
                  "authority_snapshot_changed"
            @test worker_value(result, "WORKER_DESTINATION=") == "absent"
            @test !ispath(fixture.destination)
        end
        for (mutation, expected) in (
                ("file-bytes", "authority_snapshot_changed"),
                ("file-mode", "authority_path_mismatch"),
                ("parent-mode", "authority_path_mismatch"),
                ("root-mode", "authority_path_mismatch"))
            with_fixture() do fixture
                result = run_worker(fixture.root, "control-replacement";
                                    extra=mutation)
                @test result.success
                @test occursin("WORKER_STATUS=ok", result.stdout)
                @test worker_value(result, "WORKER_VERIFY_REASON=") == expected
                @test worker_value(result, "WORKER_PUBLISH_STATE=") ==
                      "not_committed"
                @test worker_value(result, "WORKER_PUBLISH_REASON=") == expected
                @test worker_value(result, "WORKER_DESTINATION=") == "absent"
                @test !ispath(fixture.destination)
            end
        end
    end

    # Group 5: M == R with the R-local control file also declared as an
    # approved input must be rejected as an overlap (never silently dropped),
    # at authorization, at context revalidation and at manifest construction.
    alias_result = @testset "MC remediation: M==R control-as-input overlap is rejected" begin
        with_fixture() do fixture
            @test fixture.main != fixture.root
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            malicious_control_as_input!(fixture)
            @test E._t13_relative_from_option(
                      fixture.root,
                      fixture.options["--candidate-config"]) == BOULDER_RELATIVE
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "activation_binding_mismatch"
            @test !ispath(fixture.destination)
        end
        sandbox = mktempdir()
        try
            bindings = dummy_bindings(sandbox)
            control = bindings.control
            snapshots = [control.file]
            context = E._ProductionContext(
                sandbox, snapshots, Tuple{String,Vector{String}}[], "", bindings)
            error = try E._verify_context(context); nothing catch caught; caught end
            @test error isa E.EvaluatorError
            @test error.reason == :activation_binding_mismatch
            error = try
                E._manifest_records(sandbox, snapshots,
                                    Tuple{String,Vector{String}}[], String[];
                                    control=control)
                nothing
            catch caught
                caught
            end
            @test error isa E.EvaluatorError
            @test error.reason == :activation_binding_mismatch
        finally
            rm(sandbox; recursive=true, force=true)
        end
    end

    # Group 6: the retained existing-destination fast path and forced-EEXIST
    # cases, plus a real postcommit hook that fires exactly once on whichever
    # commit event the filesystem actually takes.
    hooks_result = @testset "MC remediation: existing/EEXIST and real postcommit hooks" begin
        for mutation in ("existing-fast-path",
                         "preferred-rename-destination-exists")
            with_fixture() do fixture
                result = run_worker(fixture.root, "publication-hooks";
                                    extra=mutation)
                @test result.success
                @test occursin("WORKER_STATUS=ok", result.stdout)
                @test worker_value(result, "WORKER_PUBLISH_STATE=") ==
                      "not_committed"
                @test worker_value(result, "WORKER_PUBLISH_REASON=") ==
                      "authority_snapshot_changed"
                @test worker_value(result, "WORKER_DESTINATION=") == "present"
                @test worker_value(result, "WORKER_DESTINATION_FILES=") ==
                      string(length(E._FILES))
                @test ispath(fixture.destination)
            end
        end
        with_fixture() do fixture
            result = run_worker(fixture.root, "publication-hooks";
                                extra="postcommit")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test worker_value(result, "WORKER_POSTCOMMIT_MUTATIONS=") == "1"
            @test worker_value(result, "WORKER_POSTCOMMIT_EVENTS=") in
                  ("after_preferred_rename", "after_receipt_link")
            @test worker_value(result, "WORKER_CONTROL_SHA_BEFORE=") !=
                  worker_value(result, "WORKER_CONTROL_SHA_AFTER=")
            @test worker_value(result, "WORKER_CONTROL_MODE_RESTORED=") == "true"
            @test worker_value(result, "WORKER_PUBLISH_STATE=") == "ambiguous"
            @test worker_value(result, "WORKER_PUBLISH_REASON=") ==
                  "authority_snapshot_changed"
            @test worker_value(result, "WORKER_DESTINATION=") == "present"
            @test worker_value(result, "WORKER_DESTINATION_FILES=") ==
                  string(length(E._FILES))
            @test ispath(fixture.destination)
        end
    end

    # Group 7: actual post-read helper behavior plus precise source wiring and
    # stage order.  The static checks use the real frozen source text; no
    # deterministic mid-read race is claimed here.
    static_result = @testset "MC remediation: post-read recheck and loader stage boundaries" begin
        sandbox = mktempdir()
        try
            control = dummy_control(sandbox)
            @test E._verify_main_control(control) === nothing
            file = control.file.path
            original_mode = UInt(stat(file).mode) & UInt(0o777)
            chmod(file, original_mode == 0o444 ? 0o644 : 0o444)
            error = try E._verify_main_control(control); nothing catch caught; caught end
            @test error isa E.EvaluatorError
            @test error.reason == :authority_path_mismatch
            chmod(file, original_mode)
            @test E._verify_main_control(control) === nothing
            backup = joinpath(sandbox, "retained-boulder.json")
            mv(file, backup)
            cp(backup, file; force=true)
            chmod(file, original_mode)
            @test (UInt64(stat(file).device), UInt64(stat(file).inode)) !=
                  (control.file.device, control.file.inode)
            error = try E._verify_main_control(control); nothing catch caught; caught end
            @test error isa E.EvaluatorError
            @test error.reason == :authority_snapshot_changed
        finally
            rm(sandbox; recursive=true, force=true)
        end
        source = _evaluator_source_text()
        verify_body = _source_function_body(source, "function _verify_main_control(")
        @test _source_ordered(verify_body, [
            "_main_control_directory_recheck(control.root",
            "_main_control_directory_recheck(parent",
            "_verify_snapshot(control.file",
            "_main_control_file_recheck(control.file.path",
            "_main_control_directory_recheck(control.root",
            "_main_control_directory_recheck(parent",
        ]) !== nothing
        read_body = _source_function_body(source, "function _read_main_control(")
        @test _source_ordered(read_body, [
            "_snapshot(root, _T13_CONTROL_FILENAME",
            "_main_control_file_recheck(file.path",
            "_main_control_directory_recheck(root",
            "_main_control_directory_recheck(parent",
        ]) !== nothing
        loader_body = _source_function_body(source, "function _load_production_input(")
        @test _source_ordered(loader_body, [
            "T8.load_contract",
            "_verify_context(context)",
            "T11.load_admission_data",
            "_verify_context(context)",
            "T11.evaluate_admission",
            "_verify_context(context)",
            "T11.report_files",
            "_verify_context(context)",
            "T11.normalize_admission_data",
            "_verify_context(context)",
            "T12.infer_structured_chains",
            "_verify_context(context)",
            "return input",
        ]) !== nothing
    end

    return (cli_result, two_root_result, initial_result, inode_result,
            alias_result, hooks_result, static_result)
end

function _print_mc_remediation_smoke_summary(smoke_results)
    total_passes = 0
    total_fails = 0
    total_errors = 0
    total_broken = 0
    for smoke_result in smoke_results
        counts = _activation_test_counts(smoke_result)
        total_passes += counts.passes
        total_fails += counts.fails
        total_errors += counts.errors
        total_broken += counts.broken
        println("MC_REMEDIATION_SMOKE_TESTSET=", smoke_result.description,
                " passes=", counts.passes,
                " fails=", counts.fails,
                " errors=", counts.errors,
                " broken=", counts.broken)
    end
    println("MC_REMEDIATION_SMOKE_TOTAL=passes=", total_passes,
            " fails=", total_fails, " errors=", total_errors,
            " broken=", total_broken)
    return nothing
end

if "--mc-remediation-smoke" in ARGS
    _print_mc_remediation_smoke_summary(_mc_remediation_smoke_tests())
    exit(0)
end

if "--context-smoke" in ARGS
    _print_context_smoke_summary(_context_smoke_tests())
    exit(0)
end

if "--dispatch-smoke" in ARGS
    _dispatch_smoke()
    exit(0)
end

_harness_check()

if "--harness-check" in ARGS
    exit(0)
end

const RESULTS = @testset "T13 activation inactive candidate" verbose=true begin

    @testset "predecessor lookup is stem-to-path exact" begin
        @test E._t13_authority_binding("gate4_transition_receipt").path ==
              E._RUNTIME_V4_AUTHORITY_BINDINGS[findfirst(
                  entry -> entry[1] == "gate4_transition_receipt_path",
                  E._RUNTIME_V4_AUTHORITY_BINDINGS)][2]
        @test E._t13_authority_binding("gate4_transition_receipt").sha256 ==
              "87818b63b132d8856dec22cb389cc936283bb0827cfa4884cbfd6e8c169de092"
        unknown = try
            E._t13_authority_binding("gate4_transition")
            nothing
        catch error
            error
        end
        @test unknown isa E.EvaluatorError
        @test unknown.reason == :activation_predecessor_mismatch
        @test E._T13_GATE2_ORACLE_REVIEW_MODE == 0o444
        @test E._hash_bytes(read(joinpath(repo_root(), E._T13_GATE2_ORACLE_REVIEW_PATH))) ==
              E._T13_GATE2_ORACLE_REVIEW_SHA256
    end

    @testset "consumed include closure matches source independently" begin
        found = Set{String}()
        pending = String["test/evaluate_structured_unit_assignment.jl"]
        while !isempty(pending)
            relative = pop!(pending)
            relative in found && continue
            push!(found, relative)
            text = read(joinpath(repo_root(), relative), String)
            include_root = dirname(relative)
            dir_match = match(r"const _DIR = joinpath\(@__DIR__, \"([^\"]+)\"\)",
                              text)
            if dir_match !== nothing
                include_root = joinpath(include_root, dir_match.captures[1])
            end
            for match in eachmatch(r"include\(joinpath\(([^)]*)\)\)", text)
                literals = [entry.captures[1] for entry in
                            eachmatch(r"\"([^\"]+)\"", match.captures[1])]
                isempty(literals) && continue
                candidate = normpath(joinpath(include_root, joinpath(literals...)))
                isfile(joinpath(repo_root(), candidate)) || continue
                push!(pending, candidate)
            end
        end
        @test found == Set(String.(collect(E._T13_SOURCE_CLOSURE)))
    end

    @testset "positive real-context authorization in the bound root" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            @test !ispath(fixture.destination)
            @test E._runtime_v3_mode(
                fixture.receipt_path, UInt(CONTROL_MODE), "receipt") === nothing
        end
    end

    @testset "foreign-root byte-identical copy is rejected" begin
        with_fixture() do fixture
            foreign_scratch = realpath(get(ENV, "TMPDIR", fixture.root))
            foreign = realpath(mktempdir(foreign_scratch;
                                         prefix="t13-activation-foreign-",
                                         cleanup=false))
            rm(foreign)
            # Diagnostic retention: keep the foreign copy under the bounded
            # TMPDIR instead of deleting it in a finally block.
            _copy_fixture_tree(fixture.root, foreign)
            _verify_fixture_copy(fixture.root, foreign)
            result = run_worker(fixture.root, "authorize-foreign";
                                foreign=foreign)
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "activation_source_mismatch"
        end
    end

    @testset "guarded activation route bypasses only the default gate" begin
        with_fixture() do fixture
            result = run_worker(fixture.root, "bypass")
            @test result.success
            @test occursin("WORKER_STATUS=ok", result.stdout)
            plain = ""
            activated = ""
            receiving = ""
            for line in split(result.stdout, '\n')
                startswith(line, "WORKER_PLAIN=") && (plain = String(line[14:end]))
                startswith(line, "WORKER_ACTIVATED=") &&
                    (activated = String(line[18:end]))
                startswith(line, "WORKER_RECEIVING=") &&
                    (receiving = String(line[18:end]))
            end
            @test plain == "todo_state_mismatch" || plain in String.(IDENTITY_REASONS)
            @test activated != "todo_state_mismatch"
            # The receiving-options validator compares the full canonical
            # command before destination scope, so an --out-dir mismatch is a
            # command mismatch; destination checks stay intact.
            @test receiving == "activation_command_mismatch"
            @test !ispath(fixture.destination)
            @test !ispath(joinpath(fixture.root, "other-output"))
        end
    end

    # MC-2 remediation groups 1-7 (neutral CLI pairing and CLI-through-main,
    # real M/R and missing-M decoy, initial vs post-capture replacement, Main
    # inode, M==R control-as-input, existing/EEXIST and real postcommit hook,
    # and the post-read/stage source checks).  The same helper is exposed as
    # `--mc-remediation-smoke` for targeted runs before the full suite.
    _mc_remediation_smoke_tests()

    @testset "receipt and native Main-control schema boundaries" begin
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["boulder"]["postimage_path"] = "records/detached-copy.json"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_state_mismatch"
        end
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["boulder"]["preimage_path"] = "records/other.json"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_state_mismatch"
        end
        with_fixture() do fixture
            detached = joinpath(fixture.root, "records", "detached-copy.json")
            write(detached, "detached correctly hashed postimage copy\n")
            rewrite_receipt!(fixture) do document
                document["boulder"]["postimage_sha256"] = sha_file(detached)
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "authority_hash_mismatch"
        end
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["scope"] = "T13_REAL_CAMPAIGN"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_scope_mismatch"
        end
        # The receipt control_root is a canonical absolute path, but it must
        # agree with the spec's Main-control root: no arbitrary external
        # control directory and no execution-root substitution is accepted.
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["boulder"]["control_root"] = fixture.records
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_binding_mismatch"
        end
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["boulder"]["control_root"] = fixture.root
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_binding_mismatch"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["control"]["root"] = fixture.root
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_binding_mismatch"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["control"]["role"] = "secondary_boulder"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_binding_mismatch"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["control"]["extra_mode"] = "0444"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_schema_mismatch"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                delete!(document["control"], "parent_mode")
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_schema_mismatch"
        end
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["boulder"]["extra_control"] = "x"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_schema_mismatch"
        end
        # Arbitrary external review/data paths are rejected: every referenced
        # activation record and input stays repository-relative.
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["source_review"]["path"] = "/tmp/external-review.json"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_path_invalid"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["inputs"]["files"][1]["path"] = "/tmp/external-input.tsv"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_path_invalid"
        end
    end

    @testset "v1 activation schemas cannot enable the successor route" begin
        with_fixture() do fixture
            rewrite_receipt!(fixture) do document
                document["schema"] = "stmfit_t13_activation_receipt_v1"
            end
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "activation_schema_mismatch"
            @test !ispath(fixture.destination)
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["schema"] = "stmfit_t13_execution_spec_v1"
            end
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "activation_schema_mismatch"
            @test !ispath(fixture.destination)
        end
        @test E._T13_ACTIVATION_SCHEMA == "stmfit_t13_activation_receipt_v2"
        @test E._T13_EXECUTION_SPEC_SCHEMA == "stmfit_t13_execution_spec_v2"
        @test E._T13_ACTIVATION_AUTHORITY_SCHEMA ==
              "schema=structured-evaluator-authority-v8-t13-main-control"
        @test E._T13_CONTROL_ROLE == "main_boulder"
        @test E._T13_CONTROL_FILENAME == ".omo/boulder.json"
        @test E._T13_CONTROL_KEYS ==
              ("role", "root", "root_mode", "parent_mode", "file_mode")
    end

    @testset "launch environment bindings are enforced" begin
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["launch"]["cwd"] = "/tmp"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_launch_mismatch"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["launch"]["threads"] = document["launch"]["threads"] + 1
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_launch_mismatch"
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                document["launch"]["active_project"] = "/tmp/other/Project.toml"
            end
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "activation_launch_mismatch"
        end
        with_fixture() do fixture
            result = run_worker(fixture.root, "authorize";
                                extra_env=Dict("STMFIT_DATA_DIR" => "/tmp/forbidden"))
            @test worker_reason(result) == "activation_launch_mismatch"
        end
    end

    @testset "writable control record is rejected" begin
        with_fixture() do fixture
            chmod(fixture.receipt_path, 0o644)
            result = run_worker(fixture.root, "authorize")
            @test worker_reason(result) == "authority_path_mismatch"
        end
    end

    # Lower-level fabricated identity mutations: the snapshot is constructed
    # with `inode + 1` and the directory is recreated with `mkdir`.  These are
    # not actual filesystem replacements; the actual-replacement cases live in
    # the source, Main-control and approved-directory testsets.
    @testset "lower-level fabricated identity mutations are rejected (not filesystem replacement)" begin
        sandbox = mktempdir()
        try
            path = joinpath(sandbox, "snapshot.bin")
            write(path, "snapshot\n")
            info = stat(path)
            snapshot = E._Snapshot(path, Vector{UInt8}(read(path)), sha_file(path),
                                   UInt64(info.device), UInt64(info.inode),
                                   UInt64(info.nlink))
            duplicate = E._Snapshot(path, snapshot.bytes, snapshot.sha256,
                                    snapshot.device, snapshot.inode + 1,
                                    snapshot.nlink)
            error = try E._merge_snapshots([snapshot], [duplicate], "test"); nothing
                    catch caught; caught end
            @test error isa E.EvaluatorError
            @test error.reason == :activation_binding_mismatch

            directory = joinpath(sandbox, "watched")
            mkdir(directory)
            info = stat(directory)
            identities = Dict{String,Tuple{UInt64,UInt64}}(
                directory => (UInt64(info.device), UInt64(info.inode)))
            base = dummy_bindings(sandbox)
            bindings = E._ActivationBindings(
                base.scope, base.receipt_path, base.receipt_sha256,
                base.spec_path, base.spec_sha256,
                base.source_review_path, base.source_review_sha256,
                base.transition_review_path, base.transition_review_sha256,
                base.publication_path, base.publication_sha256,
                base.boulder_preimage_path, base.boulder_preimage_sha256,
                base.boulder_postimage_path, base.boulder_postimage_sha256,
                base.entrypoint_path, base.entrypoint_sha256,
                base.destination, base.modes, base.root, base.command_body,
                identities, base.control)
            rm(directory; recursive=true, force=true)
            mkdir(directory)
            error = try E._verify_activation_bindings(bindings); nothing
                    catch caught; caught end
            @test error isa E.EvaluatorError
            @test error.reason == :activation_snapshot_changed
        finally
            rm(sandbox; recursive=true, force=true)
        end
    end

    @testset "activated receipt schema and legacy equivalence" begin
        report = dummy_report()
        tsv = Dict{String,Vector{UInt8}}("bootstrap.tsv" => Vector{UInt8}("x\n"))
        legacy = String(E._receipt_bytes(report, tsv))
        legacy_direct = String(E._legacy_receipt_bytes(report, tsv))
        @test legacy == legacy_direct
        @test occursin("gateclosure_sha256", legacy)
        @test occursin("closure_v1_claim_sha256", legacy)
        @test !occursin("activation_schema", legacy)
        @test !occursin("activation_control", legacy)
        @test !occursin("main_control", legacy)
        sandbox = mktempdir()
        try
            bindings = dummy_bindings(sandbox)
            activated = String(E._receipt_bytes(report, tsv; activation=bindings))
            @test occursin(E._T13_ACTIVATED_REPORT_SCHEMA, activated)
            @test occursin("schema_version = 2", activated)
            blocker = String(E._activated_blocker_files(
                E.EvaluatorError(:BLOCKED, :activation_test_reason, "blocked"),
                "e"^64, bindings)["receipt.toml"])
            @test occursin(E._T13_ACTIVATED_BLOCKER_SCHEMA, blocker)
            @test occursin("schema_version = 2", blocker)
            @test occursin("activation_control_sha256 = " *
                           repr(bindings.control.file.sha256), blocker)
            @test occursin("activation_schema = \"$(E._T13_ACTIVATION_SCHEMA)\"",
                           activated)
            @test occursin("activation_authority_schema = " *
                           repr(E._T13_ACTIVATION_AUTHORITY_SCHEMA), activated)
            @test occursin("activation_boulder_postimage_sha256", activated)
            @test occursin("activation_control_role = " *
                           repr(E._T13_CONTROL_ROLE), activated)
            @test occursin("activation_control_root = " *
                           repr(bindings.control.root), activated)
            @test occursin("activation_control_file = " *
                           repr(E._T13_CONTROL_FILENAME), activated)
            @test occursin("activation_control_sha256 = " *
                           repr(bindings.control.file.sha256), activated)
            @test occursin("activation_control_file_mode = \"0444\"", activated)
            @test !occursin("gateclosure_sha256", activated)
            @test !occursin("closure_v1_claim_sha256", activated)
            @test !occursin("live_plan_boulder_runtime_authority", activated)
        finally
            rm(sandbox; recursive=true, force=true)
        end
    end

    @testset "lower-level publisher nine-file contract (not activation integration)" begin
        files = Dict{String,Vector{UInt8}}(name => Vector{UInt8}(name * "\n")
                                           for name in E._FILES)
        root = mktempdir()
        try
            @test E._publish_atomic(root, "out", files) == :committed_verified
            @test E._publish_atomic(root, "out", files) == :committed_verified
            @test sort(readdir(joinpath(root, "out"))) == sort(collect(E._FILES))
        finally
            rm(root; recursive=true, force=true)
        end
    end

    @testset "default gate and public producer stay blocked" begin
        document = TOML.parse(String(read(joinpath(repo_root(),
                                                   RUNTIME_V4_EVALUATOR_CONFIG))))
        denied = try E._validate_runtime_v4_execution_gate(document); nothing
                 catch caught; caught end
        @test denied isa E.EvaluatorError && denied.reason == :todo_state_mismatch
        error = try
            E.produce_unary_selection_evidence(repo_root();
                evaluator_config=RUNTIME_V4_EVALUATOR_CONFIG,
                features=".tmp-t13-activation-missing/features.tsv",
                candidate_config="config/unit_assignment_structured_candidate.toml",
                model_config="config/unit_assignment_structured_model.toml",
                universe_dir=".tmp-t13-activation-missing/universe",
                edge_dir=".tmp-t13-activation-missing/edges",
                forward_receipt=".tmp-t13-activation-missing/forward.toml",
                backward_receipt=".tmp-t13-activation-missing/backward.toml",
                admission_dir=".tmp-t13-activation-missing/admission")
            nothing
        catch caught
            caught
        end
        @test error isa E.EvaluatorError && error.status == :BLOCKED
        @test error.reason == :todo_state_mismatch || error.reason in IDENTITY_REASONS
    end

    @testset "finite F3 dependencies are bound and revalidated" begin
        with_fixture() do fixture
            target = joinpath(fixture.root, first(E._T13_FIXED_PROVENANCE_FILES))
            write(target, "mutated transitive dependency\n")
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "authority_hash_mismatch"
        end
        with_fixture() do fixture
            # Deliberate duplicate: the replacement path is already declared,
            # so spec validation rejects the duplicated dependency path before
            # the derived-set comparison can raise activation_dependency_mismatch.
            rewrite_spec!(fixture) do document
                entries = document["dependencies"]["files"]
                entry = first(candidate for candidate in entries
                              if candidate["path"] ==
                                 "test/build_structured_unit_predictions.jl")
                entry["path"] = "test/build_label_free_edge_features.jl"
            end
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "activation_binding_mismatch"
        end
        with_fixture() do fixture
            # Real unique replacement: an unlisted synthetic provenance file
            # preserves cardinality and path uniqueness, so the derived
            # closed-set comparison is the first check that fails.
            unapproved_relative = "test/activation-unapproved-provenance.txt"
            unapproved = joinpath(fixture.root, unapproved_relative)
            write(unapproved,
                  "# Synthetic unapproved provenance; never included or executed.\n")
            rewrite_spec!(fixture) do document
                entries = document["dependencies"]["files"]
                @test all(candidate -> candidate["path"] != unapproved_relative,
                          entries)
                entry = first(candidate for candidate in entries
                              if candidate["path"] ==
                                 "test/build_structured_unit_predictions.jl")
                entry["path"] = unapproved_relative
                entry["sha256"] = sha_file(unapproved)
                entry["mode"] = mode_text(unapproved)
            end
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "activation_dependency_mismatch"
            @test !ispath(fixture.destination)
        end
        with_fixture() do fixture
            rewrite_spec!(fixture) do document
                files = document["dependencies"]["files"]
                index = findfirst(candidate -> candidate["path"] ==
                                  "test/build_structured_unit_predictions.jl", files)
                deleteat!(files, index)
            end
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test worker_reason(result) == "activation_dependency_mismatch"
        end
        with_fixture() do fixture
            extra_path = joinpath(fixture.root, "test",
                                  "test_structured_evaluator.jl")
            isfile(extra_path) || write(
                extra_path, "# Synthetic extra dependency; never included/executed.\n")
            rewrite_spec!(fixture) do document
                push!(document["dependencies"]["files"],
                      Dict{String,Any}(
                          "path" => "test/test_structured_evaluator.jl",
                          "sha256" => sha_file(joinpath(fixture.root,
                                                        "test/test_structured_evaluator.jl")),
                          "mode" => mode_text(joinpath(fixture.root,
                                                       "test/test_structured_evaluator.jl"))))
            end
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test worker_reason(result) == "activation_dependency_mismatch"
        end
        with_fixture() do fixture
            write(joinpath(fixture.root, "test/lib/hierarchical", "intruder.jl"),
                  "unapproved hierarchical member\n")
            result = run_worker(fixture.root, "authorize")
            @test result.success
            @test worker_reason(result) == "activation_binding_mismatch"
        end
        with_fixture() do fixture
            result = run_worker(fixture.root, "toctou";
                                extra=first(E._T13_FIXED_PROVENANCE_FILES))
            @test result.success
            @test occursin("WORKER_STATUS=blocked", result.stdout)
            @test worker_reason(result) == "authority_snapshot_changed"
        end
    end

    _context_smoke_tests()

    @testset "activation coverage map successor entries" begin
        coverage = TOML.parsefile(COVERAGE_MAP)
        @test coverage["schema"] == "structured_evaluator_gate5_coverage_v1"
        ids = Set(String(record["id"]) for record in coverage["category"])
        for id in ("activation_default_denial", "activation_cli_pairing",
                   "activation_receipt_schema_and_digest_anchor",
                   "activation_scope_and_permissions",
                   "activation_boulder_preimage_postimage",
                   "activation_source_closure_and_entrypoint",
                   "activation_command_environment_destination",
                   "activation_input_inventory_binding",
                   "activation_toctou_mode_identity_membership",
                   "activation_receipt_provenance_fields",
                   "activation_public_producer_denial",
                   "activation_publication_equivalence",
                   "activation_assembled_context_boundary",
                   "activation_actual_source_boulder_replacement",
                   "activation_actual_directory_replacement",
                   "activation_main_control_two_root_positive",
                   "activation_main_control_schema_boundaries",
                   "activation_main_control_v1_rejection",
                   "activation_main_control_preauthorization",
                   "activation_main_control_replacement_mutation",
                   "activation_main_control_publication_hooks",
                   "activation_main_control_manifest_tagging",
                   "activation_main_control_r_decoy",
                   "activation_cli_through_main_neutral",
                   "activation_main_control_initial_acquisition",
                   "activation_main_control_postread_recheck",
                   "activation_main_control_input_overlap",
                   "activation_loader_stage_context_revalidation")
            @test id in ids
        end
        retired = [record for record in coverage["category"]
                   if record["status"] == "intentionally_retired"]
        @test length(retired) == 1
        @test retired[1]["id"] == "integrated_historical_v0_v3_positive"
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    counts = _activation_test_counts(RESULTS)
    println("assertions=", counts.passes,
            " fails=", counts.fails,
            " errors=", counts.errors,
            " broken=", counts.broken)
end

end
