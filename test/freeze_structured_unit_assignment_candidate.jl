#!/usr/bin/env julia

module StructuredFreezeLifecycle

using SHA
using TOML

const ROOT = realpath(dirname(@__DIR__))
const CANDIDATE_SHA256 = "09bf73577bdfbcc2fd6a2643c1f80c872bd14c21da38a2b68c3c056c8b7f69fd"
const MODEL_SHA256 = "b3bac29d7dbecb0a9a46ec4b81a283c6b6cd4dda586c639b29d8ea105ecbd5ad"
const EVIDENCE_MODEL_ID = "structured-v2"
const ALLOWED_ENVIRONMENT_KEYS = [
    "PATH", "HOME", "USER", "JULIA_DEPOT_PATH", "JULIA_LOAD_PATH",
    "LANG", "LC_ALL", "GKSwstype",
]
const FREEZE_KEYS = [
    "schema", "schema_version", "state", "campaign_id", "attempt_id",
    "candidate_path", "candidate_sha256", "model_path", "model_sha256",
    "prediction_path", "prediction_sha256", "source_bundle_sha256",
    "universe_sha256", "evidence_index_sha256", "eligibility_receipt_sha256",
    "grader_source_sha256", "previous_receipt", "previous_receipt_sha256",
]
const RESERVATION_KEYS = [
    "schema", "schema_version", "state", "campaign_id", "attempt_id",
    "candidate_sha256", "model_sha256", "prediction_path", "prediction_sha256",
    "repository_root", "project_path", "project_sha256", "grade_wrapper_path",
    "grade_wrapper_sha256", "report_script_path", "report_script_sha256",
    "exact_command", "working_directory", "report_destination", "environment_keys",
    "previous_receipt", "previous_receipt_sha256",
]
const FINAL_KEYS = [
    "schema", "schema_version", "state", "outcome", "error_reason",
    "campaign_id", "attempt_id", "candidate_sha256", "model_sha256", "profile",
    "prediction_path", "prediction_sha256", "process_exit_code", "invocation_count",
    "summary_tsv_path", "summary_tsv_sha256", "report_md_path", "report_md_sha256",
    "lobe_position_errors_tsv_path", "lobe_position_errors_tsv_sha256",
    "grade_tsv_path", "grade_tsv_sha256", "stdout_path", "stdout_sha256",
    "stderr_path", "stderr_sha256", "previous_receipt", "previous_receipt_sha256",
]

include(joinpath(@__DIR__, "lib", "structured_assignment", "evidence.jl"))
using .StructuredEvidence

export ALLOWED_ENVIRONMENT_KEYS,
       CANDIDATE_SHA256,
       EVIDENCE_MODEL_ID,
       FINAL_KEYS,
       FREEZE_KEYS,
       DurabilityWarning,
       LifecycleContext,
       LifecycleError,
       MODEL_SHA256,
       PublicationOutcome,
       RESERVATION_KEYS,
       ROOT,
       FreezeRequest,
       acquire_create_lock,
       atomic_create_file,
       derive_lifecycle,
       exact_grade_command,
       freeze_candidate,
       grade_artifact_paths,
       grade_environment,
       production_context,
       read_exact_final,
       read_exact_freeze,
       read_exact_reservation,
       release_create_lock,
       serialize_final,
       serialize_freeze,
       serialize_reservation,
       sha256_file,
       snapshot_bytes,
       validate_campaign_directory,
       validate_freeze_chain,
       validate_reservation_chain,
       with_durability_fault_hook,
       write_durability_warnings

struct LifecycleError <: Exception
    code::String
    message::String
end

Base.showerror(io::IO, error::LifecycleError) =
    print(io, error.code, ": ", error.message)

struct LifecycleContext
    root::String
    candidate_path::String
    model_path::String
    project_path::String
    grade_wrapper_path::String
    report_script_path::String
end

struct FreezeRequest
    campaign_dir::String
    campaign_id::String
    attempt_id::String
    model_id::String
    eligibility_receipt::String
    expected_receipt_sha256::String
    expected_command_identity_sha256::String
end

struct GateExpectation
    name::String
    campaign_id::String
    model_id::String
    receipt_path::String
    receipt_sha256::String
    command_identity_sha256::String
end

struct FreezeManifest
    campaign_id::String
    model_id::String
    prediction_path::String
    universe_path::String
    eligibility_command_identity_sha256::String
    mandatory_gates::Vector{GateExpectation}
end

struct CreateLock
    path::String
    device::UInt64
    inode::UInt64
    acquired::Bool
end

struct DurabilityWarning
    code::String
    operation::String
    path::String
    message::String
end

struct PublicationOutcome
    path::String
    sha256::String
    durability::String
    durability_warnings::Vector{DurabilityWarning}
end

const DURABILITY_FAULT_HOOK = Ref{Union{Nothing,Function}}(nothing)

function with_durability_fault_hook(function_body::Function, hook::Function)
    previous = DURABILITY_FAULT_HOOK[]
    DURABILITY_FAULT_HOOK[] = hook
    try
        return function_body()
    finally
        DURABILITY_FAULT_HOOK[] = previous
    end
end

function inject_durability_fault(site::Symbol, path::String)
    hook = DURABILITY_FAULT_HOOK[]
    hook === nothing || hook(site, path)
    return nothing
end

one_line_error(error)::String = replace(sprint(showerror, error), r"[\t\r\n]+" => " ")

function write_durability_warnings(io::IO, warnings::Vector{DurabilityWarning})
    for warning in warnings
        println(io,
            "WARNING [$(warning.code)] operation=$(warning.operation) " *
            "path=$(toml_quote(warning.path)) message=$(toml_quote(warning.message))")
    end
    return nothing
end

fail(code::String, message::String) = throw(LifecycleError(code, message))
sha256_hex(bytes) = bytes2hex(sha256(bytes))

function validate_sha256(value, context::String)::String
    value isa String && occursin(r"^[0-9a-f]{64}$", value) ||
        fail("schema_invalid", "$context must be a lowercase SHA-256")
    return value
end

function validate_identifier(value, context::String)::String
    value isa String && occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", value) ||
        fail("invalid_identifier", "$context must be a nonempty portable identifier")
    return value
end

function validate_text(value, context::String; allow_empty::Bool=false)::String
    value isa String || fail("schema_invalid", "$context must be a string")
    !allow_empty && isempty(value) && fail("invalid_text", "$context must not be empty")
    occursin('\0', value) && fail("invalid_text", "$context contains NUL")
    any(character -> occursin(character, value), ('\t', '\n', '\r')) &&
        fail("invalid_text", "$context must be one line")
    return value
end

function toml_quote(value::String)::String
    io = IOBuffer()
    write(io, UInt8('"'))
    for character in value
        codepoint = Int(character)
        if character == '"'
            write(io, "\\\"")
        elseif character == '\\'
            write(io, "\\\\")
        elseif character == '\b'
            write(io, "\\b")
        elseif character == '\t'
            write(io, "\\t")
        elseif character == '\n'
            write(io, "\\n")
        elseif character == '\f'
            write(io, "\\f")
        elseif character == '\r'
            write(io, "\\r")
        elseif codepoint < 0x20 || codepoint == 0x7f
            write(io, "\\u", uppercase(string(codepoint; base=16, pad=4)))
        else
            write(io, string(character))
        end
    end
    write(io, UInt8('"'))
    return String(take!(io))
end

toml_string_array(values::Vector{String}) =
    "[" * join((toml_quote(value) for value in values), ", ") * "]"

function snapshot_bytes(
    path::String,
    context::String="file";
    require_immutable::Bool=false,
)::Vector{UInt8}
    islink(path) && fail("symlink_rejected", "$context is a symlink: $path")
    ispath(path) || fail("missing_file", "$context is missing: $path")
    isfile(path) || fail("special_file_rejected", "$context is not a regular file: $path")
    before = stat(path)
    require_immutable && before.mode & 0o222 != 0 && fail(
        "publication_authority_unresolved", "$context is writable: $path")
    bytes = read(path)
    islink(path) && fail("symlink_rejected", "$context became a symlink: $path")
    isfile(path) || fail("special_file_rejected", "$context stopped being regular: $path")
    after = stat(path)
    stable = before.device == after.device && before.inode == after.inode &&
              before.size == after.size && before.mtime == after.mtime &&
              before.ctime == after.ctime && before.mode == after.mode
    stable || fail("file_changed_during_read", "$context changed while read: $path")
    require_immutable && after.mode & 0o222 != 0 && fail(
        "publication_authority_unresolved", "$context is writable: $path")
    length(bytes) == after.size ||
        fail("file_changed_during_read", "$context size changed while read: $path")
    return bytes
end

sha256_file(path::String) = sha256_hex(snapshot_bytes(path, "hash source"))

function path_is_within(path::String, directory::String)::Bool
    relative = relpath(path, directory)
    parts = splitpath(relative)
    return relative != ".." && (isempty(parts) || first(parts) != "..")
end

function lexical_absolute(root::String, supplied::String, context::String)::String
    validate_text(supplied, context)
    occursin('\\', supplied) &&
        fail("noncanonical_path", "$context must use native '/' separators")
    components = split(supplied, '/'; keepempty=true)
    any(component -> component == "..", components) &&
        fail("path_escape", "$context contains parent traversal")
    any(component -> component == ".", components) &&
        fail("noncanonical_path", "$context contains a dot component")
    absolute = isabspath(supplied) ? supplied : joinpath(root, supplied)
    absolute = abspath(absolute)
    normpath(absolute) == absolute ||
        fail("noncanonical_path", "$context is not lexically canonical")
    path_is_within(absolute, root) ||
        fail("path_escape", "$context escapes repository root")
    return absolute
end

function reject_symlink_ancestors(root::String, absolute::String, context::String)
    relative = relpath(absolute, root)
    cursor = root
    for component in splitpath(relative)
        component == "." && continue
        cursor = joinpath(cursor, component)
        islink(cursor) &&
            fail("symlink_rejected", "$context contains symlink component: $cursor")
    end
    return nothing
end

function existing_file(root::String, supplied::String, context::String)::String
    absolute = lexical_absolute(root, supplied, context)
    reject_symlink_ancestors(root, absolute, context)
    ispath(absolute) || fail("missing_file", "$context is missing: $absolute")
    isfile(absolute) ||
        fail("special_file_rejected", "$context is not a regular file: $absolute")
    realpath(absolute) == absolute ||
        fail("noncanonical_path", "$context is not its canonical real path")
    return absolute
end

function existing_directory(root::String, supplied::String, context::String)::String
    absolute = lexical_absolute(root, supplied, context)
    reject_symlink_ancestors(root, absolute, context)
    ispath(absolute) || fail("missing_directory", "$context is missing: $absolute")
    isdir(absolute) ||
        fail("special_file_rejected", "$context is not a directory: $absolute")
    realpath(absolute) == absolute ||
        fail("noncanonical_path", "$context is not its canonical real path")
    return absolute
end

function repository_relative(root::String, absolute::String)::String
    path_is_within(absolute, root) || fail("path_escape", "path escapes repository root")
    return join(splitpath(relpath(absolute, root)), "/")
end

function production_context()::LifecycleContext
    root = ROOT
    return LifecycleContext(
        root,
        joinpath(root, "config", "unit_assignment_structured_candidate.toml"),
        joinpath(root, "config", "unit_assignment_structured_model.toml"),
        joinpath(root, "Project.toml"),
        joinpath(root, "test", "grade_frozen_structured_candidate.jl"),
        joinpath(root, "test", "report_unit_assignment_benchmark.jl"),
    )
end

function validate_context(context::LifecycleContext; require_grade_wrapper::Bool=true)
    context.root == ROOT || fail("root_identity_mismatch", "ROOT must equal $ROOT")
    realpath(context.root) == ROOT ||
        fail("root_identity_mismatch", "repository root changed")
    expected = production_context()
    for field in (:candidate_path, :model_path, :project_path, :report_script_path)
        getfield(context, field) == getfield(expected, field) ||
            fail("context_identity_mismatch", "$field differs from the fixed repository path")
    end
    require_grade_wrapper && context.grade_wrapper_path != expected.grade_wrapper_path &&
        fail("context_identity_mismatch", "grade_wrapper_path differs from the fixed path")
    candidate = existing_file(ROOT, context.candidate_path, "v2 candidate TOML")
    model = existing_file(ROOT, context.model_path, "v2 model TOML")
    project = existing_file(ROOT, context.project_path, "root Project.toml")
    report = existing_file(ROOT, context.report_script_path, "benchmark report source")
    sha256_file(candidate) == CANDIDATE_SHA256 ||
        fail("candidate_hash_mismatch", "v2 candidate TOML bytes changed")
    sha256_file(model) == MODEL_SHA256 ||
        fail("model_hash_mismatch", "v2 model TOML bytes changed")
    parsed = try
        TOML.parse(String(snapshot_bytes(candidate, "v2 candidate TOML")))
    catch error
        fail("candidate_schema_invalid", sprint(showerror, error))
    end
    grader = get(parsed, "grader_only", nothing)
    grader isa AbstractDict || fail("candidate_schema_invalid", "missing [grader_only]")
    declared = get(grader, "grade_script_sha256", nothing)
    validate_sha256(declared, "grader_only.grade_script_sha256")
    sha256_file(report) == declared ||
        fail("grader_source_hash_mismatch", "benchmark report source differs from candidate binding")
    if require_grade_wrapper
        existing_file(ROOT, context.grade_wrapper_path, "grade wrapper source")
    end
    return (candidate=candidate, model=model, project=project, report=report)
end

function validate_campaign_directory(context::LifecycleContext, supplied::String)::String
    campaign = existing_directory(context.root, supplied, "campaign directory")
    campaign == context.root && fail("path_escape", "campaign directory may not be repository root")
    return campaign
end

function ensure_lifecycle_directory(campaign::String)::Tuple{String,Bool}
    lifecycle = joinpath(campaign, "lifecycle")
    islink(lifecycle) && fail("symlink_rejected", "lifecycle directory is a symlink")
    if ispath(lifecycle)
        isdir(lifecycle) ||
            fail("special_file_rejected", "lifecycle path is not a directory")
        return lifecycle, false
    end
    mkdir(lifecycle; mode=0o700)
    isdir(lifecycle) && !islink(lifecycle) ||
        fail("publication_failed", "could not create lifecycle directory")
    stat(lifecycle).device == stat(campaign).device ||
        fail("cross_filesystem_publication", "lifecycle directory is on another filesystem")
    fsync_directory(campaign)
    return lifecycle, true
end

function fsync_directory(
    path::String;
    fault_site::Union{Nothing,Symbol}=nothing,
    fault_path::String=path,
)
    descriptor = ccall(:open, Cint, (Cstring, Cint), path, 0)
    descriptor >= 0 || fail("fsync_failed", "could not open directory for fsync: $path")
    try
        fault_site === nothing || inject_durability_fault(fault_site, fault_path)
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            fail("fsync_failed", "directory fsync failed: $path")
    finally
        ccall(:close, Cint, (Cint,), descriptor)
    end
    return nothing
end

function lock_matches(path::String, device::UInt64, inode::UInt64)::Bool
    (ispath(path) || islink(path)) || return false
    islink(path) && return false
    isfile(path) || return false
    current = stat(path)
    return current.device == device && current.inode == inode
end

function cleanup_owned_create_lock(path::String, device::UInt64, inode::UInt64)
    if ispath(path) || islink(path)
        lock_matches(path, device, inode) || fail(
            "stale_lock_cleanup_unresolved",
            "exclusive lock ownership changed before cleanup: $path",
        )
        rm(path; force=true)
    end
    (ispath(path) || islink(path)) && fail(
        "stale_lock_cleanup_unresolved",
        "owned exclusive lock remains after cleanup: $path",
    )
    try
        fsync_directory(
            dirname(path);
            fault_site=:lock_cleanup_parent_fsync,
            fault_path=path,
        )
    catch error
        fail(
            "stale_lock_cleanup_unresolved",
            "owned lock is absent but removal durability is unresolved: " * one_line_error(error),
        )
    end
    return nothing
end

function acquire_create_lock(parent::String, name::String)::CreateLock
    validate_text(name, "lock name")
    occursin('/', name) && fail("invalid_lock", "lock name must be a leaf")
    lock_path = joinpath(parent, name)
    islink(lock_path) && fail("symlink_rejected", "lock path is a symlink")
    ispath(lock_path) && fail("concurrent_lifecycle_operation", "lock already exists: $lock_path")
    flags = Cint(0x0001 | 0x0040 | 0x0080 | 0x80000) # O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC
    descriptor = ccall(:open, Cint, (Cstring, Cint, Cint), lock_path, flags, Cint(0o600))
    if descriptor < 0
        ispath(lock_path) &&
            fail("concurrent_lifecycle_operation", "lock was created concurrently")
        fail("lock_creation_failed", "exclusive lock creation failed")
    end
    identity = nothing
    descriptor_open = true
    try
        identity = stat(RawFD(descriptor))
        inject_durability_fault(:lock_file_fsync, lock_path)
        ccall(:fsync, Cint, (Cint,), descriptor) == 0 ||
            fail("fsync_failed", "lock fsync failed")
        ccall(:close, Cint, (Cint,), descriptor)
        descriptor_open = false
        lock_matches(lock_path, identity.device, identity.inode) ||
            fail("lock_creation_failed", "exclusive lock identity changed")
        stat(lock_path).device == stat(parent).device ||
            fail("cross_filesystem_publication", "lock is on another filesystem")
        fsync_directory(
            parent;
            fault_site=:lock_parent_fsync,
            fault_path=lock_path,
        )
        return CreateLock(lock_path, identity.device, identity.inode, true)
    catch error
        descriptor_open && ccall(:close, Cint, (Cint,), descriptor)
        identity === nothing && fail(
            "stale_lock_cleanup_unresolved",
            "exclusive lock was created but its ownership identity is unavailable: $lock_path",
        )
        cleanup_owned_create_lock(lock_path, identity.device, identity.inode)
        throw(error)
    end
end

function release_create_lock(lock::CreateLock)
    warnings = DurabilityWarning[]
    if lock.acquired && (ispath(lock.path) || islink(lock.path))
        lock_matches(lock.path, lock.device, lock.inode) ||
            fail("lock_ownership_lost", "refusing to remove a lock not owned by this invocation")
        rm(lock.path; force=true)
        try
            fsync_directory(
                dirname(lock.path);
                fault_site=:lock_release_parent_fsync,
                fault_path=lock.path,
            )
        catch error
            (ispath(lock.path) || islink(lock.path)) && fail(
                "lock_release_authority_unresolved",
                "lock removal failed after release durability error: $(one_line_error(error))",
            )
            push!(warnings, DurabilityWarning(
                "lock_removal_durability_uncertain",
                "release_create_lock",
                lock.path,
                "owned lock is absent, but parent-directory durability was not confirmed after: " *
                one_line_error(error),
            ))
        end
    end
    return warnings
end

function reconcile_exact_publication(path::String, bytes::Vector{UInt8})::Tuple{Bool,String}
    try
        snapshot_bytes(path, "publication reconciliation"; require_immutable=true) == bytes ||
            return false, "destination bytes differ from the intended immutable receipt"
        return true, "exact intended immutable bytes are installed"
    catch error
        return false, "destination could not be validated: $(one_line_error(error))"
    end
end

function publication_warning(path::String, error)::DurabilityWarning
    return DurabilityWarning(
        "publication_durability_uncertain",
        "atomic_create_file",
        path,
        "exact intended immutable bytes are authoritative, but parent-directory " *
        "durability was not confirmed after: $(one_line_error(error))",
    )
end

function atomic_create_file(path::String, bytes::Vector{UInt8})::PublicationOutcome
    parent = dirname(path)
    isdir(parent) && !islink(parent) ||
        fail("publication_failed", "publication parent is not a real directory")
    islink(path) && fail("symlink_rejected", "publication destination is a symlink")
    ispath(path) && fail("immutable_receipt_exists", "refusing to replace existing receipt: $path")
    stage, io = mktemp(parent; cleanup=false)
    installed = false
    pending_error = nothing
    try
        stat(stage).device == stat(parent).device ||
            fail("cross_filesystem_publication", "stage is on another filesystem")
        write(io, bytes)
        flush(io)
        inject_durability_fault(:receipt_preinstall_fsync, path)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 ||
            fail("fsync_failed", "receipt stage fsync failed")
        close(io)
        chmod(stage, 0o444)
        open(stage, "r") do verify
            ccall(:fsync, Cint, (Cint,), fd(verify)) == 0 ||
                fail("fsync_failed", "receipt metadata fsync failed")
        end
        snapshot_bytes(stage, "receipt stage") == bytes ||
            fail("publication_failed", "receipt stage bytes changed")
        ispath(path) && fail("immutable_receipt_exists", "receipt appeared concurrently")
        linked = ccall(:link, Cint, (Cstring, Cstring), stage, path)
        if linked != 0
            ispath(path) && fail("immutable_receipt_exists", "receipt appeared concurrently")
            fail("publication_failed", "atomic no-replace receipt link failed")
        end
        installed = true
        rm(stage; force=true)
        fsync_directory(
            parent;
            fault_site=:receipt_postinstall_parent_fsync,
            fault_path=path,
        )
        snapshot_bytes(path, "published receipt") == bytes ||
            fail("publication_failed", "published receipt bytes changed")
        stat(path).mode & 0o222 == 0 ||
            fail("publication_failed", "published receipt is not immutable")
    catch error
        pending_error = error
    finally
        isopen(io) && close(io)
        if ispath(stage)
            try
                rm(stage; force=true)
            catch cleanup_error
                pending_error === nothing && (pending_error = cleanup_error)
            end
        end
    end
    digest = sha256_hex(bytes)
    if pending_error !== nothing
        if installed
            exact, reason = reconcile_exact_publication(path, bytes)
            if exact
                warning = publication_warning(path, pending_error)
                return PublicationOutcome(path, digest, "uncertain", [warning])
            end
            fail(
                "publication_authority_unresolved",
                "receipt install reached its commit point but authority is unresolved: $reason; " *
                "original_error=$(one_line_error(pending_error))",
            )
        end
        throw(pending_error)
    end
    return PublicationOutcome(path, digest, "confirmed", DurabilityWarning[])
end

function require_keys(document::AbstractDict, expected::Vector{String}, context::String)
    actual = Set(String(key) for key in keys(document))
    wanted = Set(expected)
    actual == wanted ||
        fail("schema_invalid", "$context keys differ; expected $(sort(expected)), got $(sort!(collect(actual)))")
    return nothing
end

expect_string(document, key::String, context::String; allow_empty::Bool=false) =
    validate_text(get(document, key, nothing), "$context.$key"; allow_empty)

function expect_integer(document, key::String, context::String)::Int
    value = get(document, key, nothing)
    value isa Integer && !(value isa Bool) ||
        fail("schema_invalid", "$context.$key must be an integer")
    return try
        Int(value)
    catch
        fail("schema_invalid", "$context.$key is outside Int range")
    end
end

function expect_string_array(document, key::String, context::String)::Vector{String}
    value = get(document, key, nothing)
    value isa Vector && all(item -> item isa String, value) ||
        fail("schema_invalid", "$context.$key must be an array of strings")
    return String[String(item) for item in value]
end

function serialize_freeze(document::AbstractDict)::Vector{UInt8}
    require_keys(document, FREEZE_KEYS, "freeze receipt")
    io = IOBuffer()
    for key in FREEZE_KEYS
        value = get(document, key, nothing)
        if key == "schema_version"
            value isa Integer && !(value isa Bool) || fail("schema_invalid", "$key must be integer")
            write(io, key, " = ", string(value), '\n')
        else
            value isa String || fail("schema_invalid", "$key must be string")
            write(io, key, " = ", toml_quote(value), '\n')
        end
    end
    return take!(io)
end

function serialize_reservation(document::AbstractDict)::Vector{UInt8}
    require_keys(document, RESERVATION_KEYS, "grade reservation")
    io = IOBuffer()
    for key in RESERVATION_KEYS
        value = get(document, key, nothing)
        if key == "schema_version"
            value isa Integer && !(value isa Bool) || fail("schema_invalid", "$key must be integer")
            write(io, key, " = ", string(value), '\n')
        elseif key == "environment_keys"
            value isa Vector && all(item -> item isa String, value) ||
                fail("schema_invalid", "environment_keys must be strings")
            write(io, key, " = ", toml_string_array(String[String(item) for item in value]), '\n')
        else
            value isa String || fail("schema_invalid", "$key must be string")
            write(io, key, " = ", toml_quote(value), '\n')
        end
    end
    return take!(io)
end

function serialize_final(document::AbstractDict)::Vector{UInt8}
    require_keys(document, FINAL_KEYS, "grade final")
    io = IOBuffer()
    integer_keys = Set(["schema_version", "process_exit_code", "invocation_count"])
    for key in FINAL_KEYS
        value = get(document, key, nothing)
        if key in integer_keys
            value isa Integer && !(value isa Bool) || fail("schema_invalid", "$key must be integer")
            write(io, key, " = ", string(value), '\n')
        else
            value isa String || fail("schema_invalid", "$key must be string")
            write(io, key, " = ", toml_quote(value), '\n')
        end
    end
    return take!(io)
end

function parse_exact_receipt(path::String, keys::Vector{String}, serializer::Function, context::String)
    bytes = snapshot_bytes(path, context; require_immutable=true)
    text = try
        String(copy(bytes))
    catch error
        fail("schema_invalid", "$context is not UTF-8: $(sprint(showerror, error))")
    end
    startswith(text, "\ufeff") && fail("schema_invalid", "$context contains a BOM")
    occursin('\r', text) && fail("schema_invalid", "$context must use LF line endings")
    endswith(text, "\n") || fail("schema_invalid", "$context must be LF terminated")
    document = try
        TOML.parse(text)
    catch error
        fail("schema_invalid", "$context TOML parse failed: $(sprint(showerror, error))")
    end
    require_keys(document, keys, context)
    serializer(document) == bytes ||
        fail("schema_invalid", "$context is not in deterministic serialization")
    return document, bytes
end

read_exact_freeze(path::String) =
    parse_exact_receipt(path, FREEZE_KEYS, serialize_freeze, "freeze receipt")
read_exact_reservation(path::String) =
    parse_exact_receipt(path, RESERVATION_KEYS, serialize_reservation, "grade reservation")
read_exact_final(path::String) =
    parse_exact_receipt(path, FINAL_KEYS, serialize_final, "grade final")

function parse_freeze_manifest(path::String, root::String)::FreezeManifest
    bytes = snapshot_bytes(path, "freeze manifest")
    document = try
        TOML.parse(String(copy(bytes)))
    catch error
        fail("freeze_manifest_invalid", sprint(showerror, error))
    end
    require_keys(document, [
        "schema", "schema_version", "campaign_id", "model_id", "prediction_path",
        "universe_path", "eligibility_command_identity_sha256", "mandatory_gates",
    ], "freeze manifest")
    get(document, "schema", nothing) == "structured_freeze_manifest_v1" ||
        fail("freeze_manifest_invalid", "freeze manifest schema is invalid")
    expect_integer(document, "schema_version", "freeze manifest") == 1 ||
        fail("freeze_manifest_invalid", "freeze manifest schema_version must equal 1")
    campaign_id = validate_identifier(get(document, "campaign_id", nothing), "freeze manifest campaign")
    model_id = validate_identifier(get(document, "model_id", nothing), "freeze manifest model")
    prediction = existing_file(root, expect_string(document, "prediction_path", "freeze manifest"),
        "frozen prediction")
    universe = existing_file(root, expect_string(document, "universe_path", "freeze manifest"),
        "frozen universe")
    command_identity = validate_sha256(
        get(document, "eligibility_command_identity_sha256", nothing),
        "freeze manifest eligibility command identity",
    )
    gates_value = get(document, "mandatory_gates", nothing)
    gates_value isa Vector ||
        fail("freeze_manifest_invalid", "mandatory_gates must be an array of tables")
    isempty(gates_value) &&
        fail("freeze_manifest_invalid", "at least one mandatory gate receipt is required")
    gates = GateExpectation[]
    for (index, value) in enumerate(gates_value)
        value isa AbstractDict ||
            fail("freeze_manifest_invalid", "mandatory_gates[$index] must be a table")
        require_keys(value, [
            "name", "campaign_id", "model_id", "receipt_path", "receipt_sha256",
            "command_identity_sha256",
        ], "mandatory_gates[$index]")
        push!(gates, GateExpectation(
            validate_identifier(get(value, "name", nothing), "mandatory gate name"),
            validate_identifier(get(value, "campaign_id", nothing), "mandatory gate campaign"),
            validate_identifier(get(value, "model_id", nothing), "mandatory gate model"),
            existing_file(root, expect_string(value, "receipt_path", "mandatory gate"),
                "mandatory gate receipt"),
            validate_sha256(get(value, "receipt_sha256", nothing), "mandatory gate receipt hash"),
            validate_sha256(get(value, "command_identity_sha256", nothing),
                "mandatory gate command identity"),
        ))
    end
    names = [gate.name for gate in gates]
    length(unique(names)) == length(names) ||
        fail("freeze_manifest_invalid", "mandatory gate names must be unique")
    names == sort(names) ||
        fail("freeze_manifest_invalid", "mandatory gates must be sorted by name")
    return FreezeManifest(campaign_id, model_id, prediction, universe, command_identity, gates)
end

function conventional_eligibility_receipt(campaign::String, campaign_id::String)::String
    return joinpath(
        campaign,
        "evidence",
        "eligibility",
        "gate-receipt-$campaign_id-$EVIDENCE_MODEL_ID.toml",
    )
end

function artifact_absolute(root::String, record)::String
    return joinpath(root, split(record.path, '/')...)
end

function exactly_one_record(index, role::String, context::String)
    records = filter(record -> record.role == role, index.artifacts)
    length(records) == 1 ||
        fail("evidence_contract_mismatch", "$context requires exactly one role=$role artifact")
    return only(records)
end

function validate_gate_receipts(context::LifecycleContext, campaign::String,
                                manifest::FreezeManifest, primary_index)
    registered = Dict(record.path => record for record in primary_index.artifacts
                      if record.role == "gate_receipt")
    length(registered) == length(manifest.mandatory_gates) ||
        fail("mandatory_gate_mismatch", "primary evidence must bind every and only mandatory gate receipt")
    for gate in manifest.mandatory_gates
        path_is_within(gate.receipt_path, campaign) ||
            fail("path_escape", "mandatory gate receipt is outside campaign")
        relative = repository_relative(context.root, gate.receipt_path)
        haskey(registered, relative) ||
            fail("mandatory_gate_mismatch", "mandatory gate $(gate.name) is absent from primary index")
        registered[relative].sha256 == gate.receipt_sha256 ||
            fail("mandatory_gate_mismatch", "mandatory gate $(gate.name) hash differs from primary index")
        result = validate_evidence_bundle(
            context.root,
            gate.receipt_path;
            expected_campaign_id=gate.campaign_id,
            expected_model_id=gate.model_id,
            expected_receipt_sha256=gate.receipt_sha256,
            expected_command_identity_sha256=gate.command_identity_sha256,
        )
        result.valid || fail("mandatory_gate_invalid",
            "mandatory gate $(gate.name) failed validation: $(join(result.reason_codes, ','))")
        result.terminal_status == "PASS" && is_v2_selection_eligible(result) ||
            fail("mandatory_gate_ineligible",
                "mandatory gate $(gate.name) terminal=$(result.terminal_status) is not eligible")
    end
    return nothing
end

function validate_primary_bundle(
    context::LifecycleContext,
    campaign::String,
    campaign_id::String,
    model_id::String,
    receipt_path::String,
    expected_receipt_sha256::String,
    expected_command_identity_sha256::String,
)
    model_id == EVIDENCE_MODEL_ID ||
        fail("model_identity_mismatch", "freeze expects model $EVIDENCE_MODEL_ID")
    expected_path = conventional_eligibility_receipt(campaign, campaign_id)
    receipt_path == expected_path ||
        fail("receipt_path_mismatch", "eligibility receipt must be $expected_path")
    path_is_within(receipt_path, campaign) ||
        fail("path_escape", "eligibility receipt is outside campaign")
    result = validate_evidence_bundle(
        context.root,
        receipt_path;
        expected_campaign_id=campaign_id,
        expected_model_id=model_id,
        expected_receipt_sha256,
        expected_command_identity_sha256,
    )
    result.valid || fail("eligibility_receipt_invalid",
        "eligibility receipt validation failed: $(join(result.reason_codes, ','))")
    result.terminal_status == "PASS" && is_v2_selection_eligible(result) ||
        fail("eligibility_receipt_ineligible",
            "eligibility receipt terminal=$(result.terminal_status) is not selection eligible")
    index_bytes = snapshot_bytes(result.artifact_index_path, "eligibility artifact index")
    index = StructuredEvidence.parse_artifact_index(index_bytes)
    freeze_record = exactly_one_record(index, "freeze_manifest", "eligibility evidence")
    freeze_path = artifact_absolute(context.root, freeze_record)
    manifest = parse_freeze_manifest(freeze_path, context.root)
    manifest.campaign_id == campaign_id ||
        fail("campaign_identity_mismatch", "freeze manifest campaign differs")
    manifest.model_id == model_id ||
        fail("model_identity_mismatch", "freeze manifest model differs")
    manifest.eligibility_command_identity_sha256 == expected_command_identity_sha256 ||
        fail("command_identity_mismatch", "freeze manifest eligibility command identity differs")
    dirname(manifest.prediction_path) == campaign ||
        fail("prediction_path_mismatch", "frozen prediction must be a direct campaign artifact")
    path_is_within(manifest.universe_path, campaign) ||
        fail("path_escape", "frozen universe is outside campaign")
    prediction_records = filter(
        record -> record.role == "output" &&
                  artifact_absolute(context.root, record) == manifest.prediction_path,
        index.artifacts,
    )
    length(prediction_records) == 1 ||
        fail("prediction_path_mismatch", "primary evidence must bind the frozen prediction as one output")
    prediction_record = only(prediction_records)
    prediction_record.sha256 == sha256_file(manifest.prediction_path) ||
        fail("prediction_hash_mismatch", "frozen prediction current bytes differ")
    universe_records = filter(
        record -> record.role == "input" &&
                  artifact_absolute(context.root, record) == manifest.universe_path,
        index.artifacts,
    )
    length(universe_records) == 1 ||
        fail("universe_path_mismatch", "primary evidence must bind the frozen universe as one input")
    universe_record = only(universe_records)
    universe_record.sha256 == sha256_file(manifest.universe_path) ||
        fail("universe_hash_mismatch", "frozen universe current bytes differ")
    validate_gate_receipts(context, campaign, manifest, index)
    return (
        result=result,
        index=index,
        index_sha256=sha256_hex(index_bytes),
        manifest=manifest,
        prediction_sha256=prediction_record.sha256,
        universe_sha256=universe_record.sha256,
    )
end

function observed_receipt_command_identity(context::LifecycleContext, receipt_path::String)::String
    receipt = StructuredEvidence.parse_gate_receipt(snapshot_bytes(receipt_path, "eligibility receipt"))
    return command_identity_sha256(context.root, receipt.metadata["exact_command"])
end

function freeze_document(context::LifecycleContext, request::FreezeRequest, bundle)::Dict{String,Any}
    return Dict{String,Any}(
        "schema" => "structured_freeze_receipt_v1",
        "schema_version" => 1,
        "state" => "frozen_once",
        "campaign_id" => request.campaign_id,
        "attempt_id" => request.attempt_id,
        "candidate_path" => context.candidate_path,
        "candidate_sha256" => CANDIDATE_SHA256,
        "model_path" => context.model_path,
        "model_sha256" => MODEL_SHA256,
        "prediction_path" => bundle.manifest.prediction_path,
        "prediction_sha256" => bundle.prediction_sha256,
        "source_bundle_sha256" => bundle.index.source_sha256,
        "universe_sha256" => bundle.universe_sha256,
        "evidence_index_sha256" => bundle.index_sha256,
        "eligibility_receipt_sha256" => request.expected_receipt_sha256,
        "grader_source_sha256" => sha256_file(context.report_script_path),
        "previous_receipt" => "absence",
        "previous_receipt_sha256" => "absence",
    )
end

function validate_freeze_document(document, context::LifecycleContext, campaign::String, bytes::Vector{UInt8})
    get(document, "schema", nothing) == "structured_freeze_receipt_v1" ||
        fail("schema_invalid", "freeze receipt schema is invalid")
    expect_integer(document, "schema_version", "freeze receipt") == 1 ||
        fail("schema_invalid", "freeze receipt schema_version must equal 1")
    get(document, "state", nothing) == "frozen_once" ||
        fail("schema_invalid", "freeze receipt state must be frozen_once")
    campaign_id = validate_identifier(get(document, "campaign_id", nothing), "freeze campaign")
    validate_identifier(get(document, "attempt_id", nothing), "freeze attempt")
    get(document, "candidate_path", nothing) == context.candidate_path ||
        fail("candidate_path_mismatch", "freeze receipt candidate path changed")
    get(document, "model_path", nothing) == context.model_path ||
        fail("model_path_mismatch", "freeze receipt model path changed")
    get(document, "candidate_sha256", nothing) == CANDIDATE_SHA256 ||
        fail("candidate_hash_mismatch", "freeze receipt candidate hash changed")
    get(document, "model_sha256", nothing) == MODEL_SHA256 ||
        fail("model_hash_mismatch", "freeze receipt model hash changed")
    get(document, "previous_receipt", nothing) == "absence" &&
        get(document, "previous_receipt_sha256", nothing) == "absence" ||
        fail("receipt_chain_mismatch", "freeze receipt must follow absence")
    for key in (
        "prediction_sha256", "source_bundle_sha256", "universe_sha256",
        "evidence_index_sha256", "eligibility_receipt_sha256", "grader_source_sha256",
    )
        validate_sha256(get(document, key, nothing), "freeze receipt $key")
    end
    prediction = existing_file(context.root,
        expect_string(document, "prediction_path", "freeze receipt"), "frozen prediction")
    dirname(prediction) == campaign ||
        fail("prediction_path_mismatch", "freeze receipt prediction is not a direct campaign artifact")
    sha256_file(prediction) == document["prediction_sha256"] ||
        fail("prediction_hash_mismatch", "frozen prediction bytes changed")
    sha256_file(context.report_script_path) == document["grader_source_sha256"] ||
        fail("grader_source_hash_mismatch", "grader/report source bytes changed")
    eligibility_path = conventional_eligibility_receipt(campaign, campaign_id)
    eligibility_path = existing_file(context.root, eligibility_path, "eligibility receipt")
    sha256_file(eligibility_path) == document["eligibility_receipt_sha256"] ||
        fail("eligibility_receipt_hash_mismatch", "eligibility receipt bytes changed")
    observed_command = observed_receipt_command_identity(context, eligibility_path)
    bundle = validate_primary_bundle(
        context,
        campaign,
        campaign_id,
        EVIDENCE_MODEL_ID,
        eligibility_path,
        document["eligibility_receipt_sha256"],
        observed_command,
    )
    bundle.manifest.eligibility_command_identity_sha256 == observed_command ||
        fail("command_identity_mismatch", "eligibility command identity changed")
    expected = Dict(
        "prediction_path" => bundle.manifest.prediction_path,
        "prediction_sha256" => bundle.prediction_sha256,
        "source_bundle_sha256" => bundle.index.source_sha256,
        "universe_sha256" => bundle.universe_sha256,
        "evidence_index_sha256" => bundle.index_sha256,
        "grader_source_sha256" => sha256_file(context.report_script_path),
    )
    for (key, value) in expected
        get(document, key, nothing) == value ||
            fail("freeze_binding_mismatch", "freeze receipt $key differs from current validated chain")
    end
    serialize_freeze(document) == bytes ||
        fail("schema_invalid", "freeze receipt serialization changed")
    return (document=document, bundle=bundle, sha256=sha256_hex(bytes))
end

function validate_freeze_chain(context::LifecycleContext, campaign::String)
    validate_context(context)
    campaign = validate_campaign_directory(context, campaign)
    lifecycle = joinpath(campaign, "lifecycle")
    isdir(lifecycle) && !islink(lifecycle) ||
        fail("missing_freeze_receipt", "lifecycle directory is absent")
    path = existing_file(context.root, joinpath(lifecycle, "freeze_receipt.toml"), "freeze receipt")
    document, bytes = read_exact_freeze(path)
    return merge(validate_freeze_document(document, context, campaign, bytes),
        (campaign=campaign, lifecycle=lifecycle, path=path))
end

function exact_grade_command(context::LifecycleContext, prediction::String, destination::String)::String
    for (value, label) in ((context.root, "ROOT"), (context.report_script_path, "report script"),
                           (prediction, "prediction"), (destination, "grade destination"))
        occursin(r"\s", value) && fail("unsupported_path", "$label contains whitespace")
    end
    return "julia --project=$(context.root) $(context.report_script_path) --full145-own-n " *
           "--profile structured_v2=$prediction --outdir $destination"
end

grade_environment() = Dict(key => ENV[key] for key in ALLOWED_ENVIRONMENT_KEYS if haskey(ENV, key))

function grade_artifact_paths(destination::String)::Dict{String,String}
    return Dict(
        "summary_tsv_path" => joinpath(destination, "summary.tsv"),
        "report_md_path" => joinpath(destination, "report.md"),
        "lobe_position_errors_tsv_path" => joinpath(destination, "lobe_position_errors.tsv"),
        "grade_tsv_path" => joinpath(destination, "grades", "structured_v2.tsv"),
        "stdout_path" => joinpath(destination, "stdout.log"),
        "stderr_path" => joinpath(destination, "stderr.log"),
    )
end

function validate_reservation_document(document, context::LifecycleContext, freeze, bytes::Vector{UInt8})
    get(document, "schema", nothing) == "structured_grade_reservation_v1" ||
        fail("schema_invalid", "grade reservation schema is invalid")
    expect_integer(document, "schema_version", "grade reservation") == 1 ||
        fail("schema_invalid", "grade reservation schema_version must equal 1")
    get(document, "state", nothing) == "grade_reserved" ||
        fail("schema_invalid", "grade reservation state must be grade_reserved")
    for key in ("campaign_id", "attempt_id", "prediction_path", "prediction_sha256",
                "candidate_sha256", "model_sha256")
        get(document, key, nothing) == get(freeze.document, key, nothing) ||
            fail("receipt_chain_mismatch", "grade reservation $key differs from freeze receipt")
    end
    get(document, "repository_root", nothing) == context.root ||
        fail("root_identity_mismatch", "grade reservation ROOT changed")
    get(document, "project_path", nothing) == context.project_path ||
        fail("project_identity_mismatch", "grade reservation project path changed")
    get(document, "grade_wrapper_path", nothing) == context.grade_wrapper_path ||
        fail("grader_source_mismatch", "grade wrapper path changed")
    get(document, "report_script_path", nothing) == context.report_script_path ||
        fail("grader_source_mismatch", "report source path changed")
    get(document, "working_directory", nothing) == context.root ||
        fail("root_identity_mismatch", "grade working directory changed")
    destination = joinpath(freeze.campaign, "grade")
    get(document, "report_destination", nothing) == destination ||
        fail("grade_destination_mismatch", "reserved destination changed")
    command = exact_grade_command(context, freeze.document["prediction_path"], destination)
    get(document, "exact_command", nothing) == command ||
        fail("command_identity_mismatch", "reserved exact command changed")
    get(document, "project_sha256", nothing) == sha256_file(context.project_path) ||
        fail("project_identity_mismatch", "root Project.toml bytes changed")
    get(document, "grade_wrapper_sha256", nothing) == sha256_file(context.grade_wrapper_path) ||
        fail("grader_source_hash_mismatch", "grade wrapper bytes changed")
    get(document, "report_script_sha256", nothing) == sha256_file(context.report_script_path) ||
        fail("grader_source_hash_mismatch", "report source bytes changed")
    get(document, "previous_receipt", nothing) == "freeze_receipt.toml" ||
        fail("receipt_chain_mismatch", "reservation previous receipt name changed")
    get(document, "previous_receipt_sha256", nothing) == freeze.sha256 ||
        fail("receipt_chain_mismatch", "reservation freeze hash changed")
    environment_keys = expect_string_array(document, "environment_keys", "grade reservation")
    length(unique(environment_keys)) == length(environment_keys) ||
        fail("schema_invalid", "environment_keys contains duplicates")
    all(key -> key in ALLOWED_ENVIRONMENT_KEYS, environment_keys) ||
        fail("schema_invalid", "environment_keys contains a forbidden key")
    serialize_reservation(document) == bytes ||
        fail("schema_invalid", "grade reservation serialization changed")
    return (document=document, sha256=sha256_hex(bytes), destination=destination)
end

function validate_reservation_chain(context::LifecycleContext, campaign::String)
    freeze = validate_freeze_chain(context, campaign)
    path = existing_file(context.root,
        joinpath(freeze.lifecycle, "grade_reservation.toml"), "grade reservation")
    document, bytes = read_exact_reservation(path)
    reservation = validate_reservation_document(document, context, freeze, bytes)
    return merge(reservation, (freeze=freeze, path=path, lifecycle=freeze.lifecycle,
                               campaign=freeze.campaign))
end

function validate_final_document(document, context::LifecycleContext, reservation, bytes::Vector{UInt8})
    get(document, "schema", nothing) == "structured_grade_final_v1" ||
        fail("schema_invalid", "grade final schema is invalid")
    expect_integer(document, "schema_version", "grade final") == 1 ||
        fail("schema_invalid", "grade final schema_version must equal 1")
    get(document, "state", nothing) == "graded" ||
        fail("schema_invalid", "grade final state must be graded")
    get(document, "profile", nothing) == "structured_v2" ||
        fail("schema_invalid", "grade final profile must be structured_v2")
    for key in ("campaign_id", "attempt_id", "candidate_sha256", "model_sha256",
                "prediction_path", "prediction_sha256")
        get(document, key, nothing) == get(reservation.document, key, nothing) ||
            fail("receipt_chain_mismatch", "grade final $key differs from reservation")
    end
    get(document, "previous_receipt", nothing) == "grade_reservation.toml" ||
        fail("receipt_chain_mismatch", "final previous receipt name changed")
    get(document, "previous_receipt_sha256", nothing) == reservation.sha256 ||
        fail("receipt_chain_mismatch", "final reservation hash changed")
    outcome = get(document, "outcome", nothing)
    outcome in ("success", "error") ||
        fail("schema_invalid", "grade final outcome must be success or error")
    reason = expect_string(document, "error_reason", "grade final"; allow_empty=true)
    exit_code = expect_integer(document, "process_exit_code", "grade final")
    invocation_count = expect_integer(document, "invocation_count", "grade final")
    if outcome == "success"
        isempty(reason) || fail("schema_invalid", "successful final must have empty error_reason")
        exit_code == 0 || fail("schema_invalid", "successful final must have exit code zero")
        invocation_count == 1 || fail("schema_invalid", "successful final must have one invocation")
    else
        isempty(reason) && fail("schema_invalid", "error final requires an error_reason")
        if invocation_count == 0
            reason == "crash_after_reservation_before_launch" ||
                fail("schema_invalid", "zero-invocation error must use the prelaunch reason")
            exit_code == -1 ||
                fail("schema_invalid", "zero-invocation error must use exit sentinel -1")
        elseif invocation_count == 1
            reason != "crash_after_reservation_before_launch" ||
                fail("schema_invalid", "launched error uses the prelaunch-only reason")
            exit_code != -1 ||
                fail("schema_invalid", "launched error uses the prelaunch-only exit sentinel")
        else
            fail("schema_invalid", "error final invocation_count must be zero or one")
        end
    end
    expected_paths = grade_artifact_paths(reservation.destination)
    for (path_key, expected_path) in expected_paths
        get(document, path_key, nothing) == expected_path ||
            fail("grade_artifact_path_mismatch", "grade final $path_key changed")
        actual_path = existing_file(context.root, expected_path, "grade final artifact")
        hash_key = replace(path_key, "_path" => "_sha256")
        validate_sha256(get(document, hash_key, nothing), "grade final $hash_key")
        sha256_file(actual_path) == document[hash_key] ||
            fail("grade_artifact_hash_mismatch", "grade final $hash_key is stale")
    end
    serialize_final(document) == bytes ||
        fail("schema_invalid", "grade final serialization changed")
    return (document=document, sha256=sha256_hex(bytes), outcome=outcome)
end

function validate_final_chain(context::LifecycleContext, campaign::String)
    reservation = validate_reservation_chain(context, campaign)
    path = existing_file(context.root,
        joinpath(reservation.lifecycle, "grade_final.toml"), "grade final")
    document, bytes = read_exact_final(path)
    final = validate_final_document(document, context, reservation, bytes)
    return merge(final, (reservation=reservation, path=path, lifecycle=reservation.lifecycle,
                         campaign=reservation.campaign))
end

function lifecycle_toml_names(lifecycle::String)::Vector{String}
    return sort!([name for name in readdir(lifecycle) if endswith(name, ".toml")])
end

function derive_lifecycle(campaign_dir::String; context::LifecycleContext=production_context())
    validate_context(context)
    campaign = validate_campaign_directory(context, campaign_dir)
    lifecycle = joinpath(campaign, "lifecycle")
    islink(lifecycle) && fail("symlink_rejected", "lifecycle directory is a symlink")
    !ispath(lifecycle) && return (state="unfrozen", campaign=campaign)
    isdir(lifecycle) || fail("special_file_rejected", "lifecycle path is not a directory")
    allowed = Set(["freeze_receipt.toml", "grade_reservation.toml", "grade_final.toml"])
    names = lifecycle_toml_names(lifecycle)
    unknown = setdiff(Set(names), allowed)
    isempty(unknown) ||
        fail("unknown_lifecycle_receipt", "unregistered lifecycle TOML: $(join(sort!(collect(unknown)), ','))")
    freeze_path = joinpath(lifecycle, "freeze_receipt.toml")
    reservation_path = joinpath(lifecycle, "grade_reservation.toml")
    final_path = joinpath(lifecycle, "grade_final.toml")
    has_freeze = ispath(freeze_path) || islink(freeze_path)
    has_reservation = ispath(reservation_path) || islink(reservation_path)
    has_final = ispath(final_path) || islink(final_path)
    !has_freeze && (has_reservation || has_final) &&
        fail("receipt_chain_mismatch", "reservation/final exists without freeze")
    has_final && !has_reservation &&
        fail("receipt_chain_mismatch", "final exists without reservation")
    !has_freeze && return (state="unfrozen", campaign=campaign, lifecycle=lifecycle)
    freeze = validate_freeze_chain(context, campaign)
    !has_reservation && return (state="frozen_once", campaign=campaign,
                                lifecycle=lifecycle, freeze=freeze)
    reservation = validate_reservation_chain(context, campaign)
    !has_final && return (state="grade_reserved", campaign=campaign,
                          lifecycle=lifecycle, freeze=freeze, reservation=reservation)
    final = validate_final_chain(context, campaign)
    return (state="graded", campaign=campaign, lifecycle=lifecycle, freeze=freeze,
            reservation=reservation, final=final)
end

function assert_no_lifecycle_receipts(lifecycle::String)
    for name in ("freeze_receipt.toml", "grade_reservation.toml", "grade_final.toml")
        path = joinpath(lifecycle, name)
        islink(path) && fail("symlink_rejected", "$name is a symlink")
        ispath(path) && fail("already_frozen", "lifecycle already contains $name")
    end
    isempty(lifecycle_toml_names(lifecycle)) ||
        fail("unknown_lifecycle_receipt", "lifecycle contains unregistered TOML")
    return nothing
end

function freeze_candidate(request::FreezeRequest; context::LifecycleContext=production_context())
    validate_context(context)
    campaign = validate_campaign_directory(context, request.campaign_dir)
    validate_identifier(request.campaign_id, "campaign_id")
    validate_identifier(request.attempt_id, "attempt_id")
    validate_identifier(request.model_id, "model_id")
    validate_sha256(request.expected_receipt_sha256, "expected eligibility receipt hash")
    validate_sha256(request.expected_command_identity_sha256,
        "expected eligibility command identity")
    receipt = existing_file(context.root, request.eligibility_receipt, "eligibility receipt")
    receipt == conventional_eligibility_receipt(campaign, request.campaign_id) ||
        fail("receipt_path_mismatch", "eligibility receipt path is not conventional")

    bundle = validate_primary_bundle(
        context,
        campaign,
        request.campaign_id,
        request.model_id,
        receipt,
        request.expected_receipt_sha256,
        request.expected_command_identity_sha256,
    )
    lifecycle, created_lifecycle = ensure_lifecycle_directory(campaign)
    lock = nothing
    committed = false
    publication = nothing
    frozen = nothing
    lock_warnings = DurabilityWarning[]
    try
        lock = acquire_create_lock(lifecycle, ".freeze.atomic-create.lock")
        assert_no_lifecycle_receipts(lifecycle)
        bundle = validate_primary_bundle(
            context,
            campaign,
            request.campaign_id,
            request.model_id,
            receipt,
            request.expected_receipt_sha256,
            request.expected_command_identity_sha256,
        )
        document = freeze_document(context, request, bundle)
        bytes = serialize_freeze(document)
        TOML.parse(String(copy(bytes)))
        path = joinpath(lifecycle, "freeze_receipt.toml")
        publication = atomic_create_file(path, bytes)
        committed = true
        frozen = validate_freeze_chain(context, campaign)
        frozen.sha256 == publication.sha256 ||
            fail("publication_failed", "freeze receipt did not validate after publication")
    finally
        lock !== nothing && append!(lock_warnings, release_create_lock(lock))
        if created_lifecycle && !committed && isdir(lifecycle) && isempty(readdir(lifecycle))
            rm(lifecycle)
            fsync_directory(campaign)
        end
    end
    warnings = vcat(publication.durability_warnings, lock_warnings)
    return merge(frozen, (
        durability=isempty(warnings) ? "confirmed" : "uncertain",
        durability_warnings=warnings,
    ))
end

function option_value(arguments::Vector{String}, index::Int, flag::String)
    index < length(arguments) || fail("cli_error", "$flag requires a value")
    value = arguments[index + 1]
    startswith(value, "--") && fail("cli_error", "$flag requires a value")
    return value, index + 2
end

function parse_cli(arguments::Vector{String})
    values = Dict{String,String}()
    flags = Dict(
        "--campaign-dir" => "campaign_dir",
        "--campaign-id" => "campaign_id",
        "--attempt-id" => "attempt_id",
        "--model-id" => "model_id",
        "--eligibility-receipt" => "eligibility_receipt",
        "--eligibility-receipt-sha256" => "expected_receipt_sha256",
        "--eligibility-command-identity-sha256" => "expected_command_identity_sha256",
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument in ("-h", "--help")
            println("""
            Usage: julia --project=. test/freeze_structured_unit_assignment_candidate.jl \\
              --campaign-dir PATH --campaign-id ID --attempt-id ID --model-id structured-v2 \\
              --eligibility-receipt PATH --eligibility-receipt-sha256 HEX \\
              --eligibility-command-identity-sha256 HEX

            Validates the concrete Todo 4 eligibility bundle, its hash-bound freeze
            manifest, and every mandatory gate receipt before atomically creating the
            immutable v2 freeze receipt. It never reads benchmark labels or grades.
            """)
            return nothing
        elseif haskey(flags, argument)
            key = flags[argument]
            haskey(values, key) && fail("cli_error", "$argument may be supplied only once")
            value, index = option_value(arguments, index, argument)
            values[key] = value
        else
            fail("cli_error", "unknown argument: $argument")
        end
    end
    missing = sort!([key for key in values(flags) if !haskey(values, key)])
    isempty(missing) || fail("cli_error", "missing required options: $(join(missing, ','))")
    return FreezeRequest(
        values["campaign_dir"],
        values["campaign_id"],
        values["attempt_id"],
        values["model_id"],
        values["eligibility_receipt"],
        values["expected_receipt_sha256"],
        values["expected_command_identity_sha256"],
    )
end

function main(arguments::Vector{String}=copy(ARGS))::Int
    try
        request = parse_cli(arguments)
        request === nothing && return 0
        frozen = freeze_candidate(request)
        write_durability_warnings(stderr, frozen.durability_warnings)
        println("FROZEN state=frozen_once receipt=$(frozen.path) sha256=$(frozen.sha256)")
        return 0
    catch error
        failure = error isa LifecycleError ? error :
            LifecycleError("internal_error", sprint(showerror, error, catch_backtrace()))
        println(stderr, "BLOCKED [$(failure.code)]: $(failure.message)")
        return 1
    end
end

end # module StructuredFreezeLifecycle

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StructuredFreezeLifecycle.main())
end
