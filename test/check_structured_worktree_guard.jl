#!/usr/bin/env julia

using SHA
using TOML
using Test

const BASELINE_SCHEMA = "stmfit-structured-worktree-baseline-v1"
const RECEIPT_SCHEMA = "stmfit-structured-worktree-receipt-v1"
const PROTECTED_PATHS = [
    "docs/src/journal.md",
    "hpc/v24_pair_sweep.sbatch",
    "test/detect_missed_lobes.py",
    "test/refit_missed_lobes.py",
]
const JOURNAL_PATH = first(PROTECTED_PATHS)

struct GuardFailure <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, err::GuardFailure) = print(io, err.code, ": ", err.message)

struct RepositoryIdentity
    root::String
    git_toplevel::String
    git_dir::String
    sha256::String
end

struct ProtectedSnapshot
    path::String
    exists::Bool
    tracking::String
    status_entry::String
    byte_count::Int
    sha256::String
end

struct WorktreeBaseline
    repository::RepositoryIdentity
    status_entries::Vector{String}
    status_sha256::String
    journal_prefix_byte_count::Int
    journal_prefix_sha256::String
    protected::Vector{ProtectedSnapshot}
end

struct CliOptions
    mode::Symbol
    target::Union{Nothing,String}
    baseline::Union{Nothing,String}
    receipt::Union{Nothing,String}
end

sha256_hex(bytes::AbstractVector{UInt8}) = bytes2hex(SHA.sha256(bytes))

function string_bytes(value::String)::Vector{UInt8}
    return collect(codeunits(value))
end

function toml_quote(value::String)::String
    io = IOBuffer()
    write(io, UInt8('"'))
    for char in value
        codepoint = Int(char)
        if char == '"'
            write(io, "\\\"")
        elseif char == '\\'
            write(io, "\\\\")
        elseif char == '\b'
            write(io, "\\b")
        elseif char == '\t'
            write(io, "\\t")
        elseif char == '\n'
            write(io, "\\n")
        elseif char == '\f'
            write(io, "\\f")
        elseif char == '\r'
            write(io, "\\r")
        elseif codepoint < 0x20 || codepoint == 0x7f
            if codepoint <= 0xffff
                write(io, "\\u", uppercase(string(codepoint; base=16, pad=4)))
            else
                write(io, "\\U", uppercase(string(codepoint; base=16, pad=8)))
            end
        else
            write(io, string(char))
        end
    end
    write(io, UInt8('"'))
    return String(take!(io))
end

function write_string_array(io::IO, name::String, values::Vector{String})
    if isempty(values)
        write(io, name, " = []\n")
        return
    end
    write(io, name, " = [\n")
    for value in values
        write(io, "  ", toml_quote(value), ",\n")
    end
    write(io, "]\n")
end

function canonical_status_bytes(entries::Vector{String})::Vector{UInt8}
    isempty(entries) && return UInt8[]
    return string_bytes(join(entries, "\n") * "\n")
end

function git_read(root::String, arguments::Vector{String})::String
    command = Cmd(vcat(["git", "-C", root], arguments))
    try
        return read(command, String)
    catch err
        throw(GuardFailure(:git_command_failed,
            "git $(join(arguments, ' ')) failed: $(sprint(showerror, err))"))
    end
end

function git_success(root::String, arguments::Vector{String})::Bool
    command = Cmd(vcat(["git", "-C", root], arguments))
    return success(pipeline(command; stdout=devnull, stderr=devnull))
end

function status_lines(text::String)::Vector{String}
    isempty(text) && return String[]
    lines = split(chomp(text), '\n'; keepempty=false)
    result = sort!(String[String(line) for line in lines])
    length(result) == length(unique(result)) ||
        throw(GuardFailure(:git_status_ambiguous, "git status contains duplicate entries"))
    return result
end

function git_status_entries(root::String)::Vector{String}
    return status_lines(git_read(root, ["status", "--short"]))
end

function git_path_status(root::String, relative_path::String)::String
    entries = status_lines(git_read(root, ["status", "--short", "--", relative_path]))
    return join(entries, "\n")
end

function git_tracking(root::String, relative_path::String)::String
    tracked = git_success(root, ["ls-files", "--error-unmatch", "--", relative_path])
    return tracked ? "tracked" : "untracked"
end

function repository_identity(root::String)::RepositoryIdentity
    canonical_root = realpath(root)
    git_toplevel = realpath(chomp(git_read(canonical_root, ["rev-parse", "--show-toplevel"])))
    git_toplevel == canonical_root || throw(GuardFailure(:repository_root_mismatch,
        "script root $canonical_root is not git top-level $git_toplevel"))
    git_dir = realpath(chomp(git_read(canonical_root, ["rev-parse", "--absolute-git-dir"])))
    material = "stmfit-structured-worktree-root-v1\n" *
               "root=$(canonical_root)\n" *
               "git_toplevel=$(git_toplevel)\n" *
               "git_dir=$(git_dir)\n"
    return RepositoryIdentity(canonical_root, git_toplevel, git_dir,
        sha256_hex(string_bytes(material)))
end

function has_parent_component(path::String)::Bool
    normalized = replace(path, '\\' => '/')
    return any(component -> component == "..", split(normalized, '/'; keepempty=true))
end

function resolve_repository_path(
    root::String,
    supplied::String;
    must_exist::Bool=false,
    output::Bool=false,
)::Tuple{String,String}
    isempty(supplied) && throw(GuardFailure(:invalid_path, "path must not be empty"))
    occursin('\0', supplied) && throw(GuardFailure(:invalid_path, "path contains NUL"))
    has_parent_component(supplied) &&
        throw(GuardFailure(:path_escape, "parent traversal is forbidden: $supplied"))

    canonical_root = realpath(root)
    candidate = normpath(isabspath(supplied) ? supplied : joinpath(canonical_root, supplied))
    relative_native = relpath(candidate, canonical_root)
    relative_parts = splitpath(relative_native)
    (relative_native == "." || (!isempty(relative_parts) && first(relative_parts) == "..")) &&
        throw(GuardFailure(:path_escape, "path is outside repository root: $supplied"))
    relative_path = join(relative_parts, "/")

    cursor = canonical_root
    for (index, component) in enumerate(relative_parts)
        cursor = joinpath(cursor, component)
        islink(cursor) && throw(GuardFailure(:symlink_rejected,
            "symlink component is forbidden: $relative_path"))
        if index < length(relative_parts)
            isdir(cursor) || throw(GuardFailure(:invalid_path,
                "parent directory does not exist: $(dirname(relative_path))"))
        end
    end

    if must_exist
        isfile(candidate) || throw(GuardFailure(:missing_path,
            "required regular file does not exist: $relative_path"))
    elseif output
        isdir(dirname(candidate)) || throw(GuardFailure(:invalid_path,
            "output parent directory does not exist: $(dirname(relative_path))"))
        if ispath(candidate) && !isfile(candidate)
            throw(GuardFailure(:invalid_path, "output is not a regular file: $relative_path"))
        end
    end

    return candidate, relative_path
end

function protected_absolute_path(root::String, relative_path::String)::String
    relative_path in PROTECTED_PATHS ||
        throw(GuardFailure(:invalid_protected_path, "unexpected protected path: $relative_path"))
    absolute_path = joinpath(root, split(relative_path, '/')...)
    cursor = root
    for component in split(relative_path, '/')
        cursor = joinpath(cursor, component)
        islink(cursor) && throw(GuardFailure(:protected_symlink,
            "protected path contains a symlink: $relative_path"))
    end
    return absolute_path
end

function snapshot_protected(
    root::String,
    relative_path::String;
    require_exists::Bool,
)::ProtectedSnapshot
    absolute_path = protected_absolute_path(root, relative_path)
    exists = isfile(absolute_path)
    if ispath(absolute_path) && !exists
        throw(GuardFailure(:protected_not_regular,
            "protected path is not a regular file: $relative_path"))
    end
    require_exists && !exists && throw(GuardFailure(:protected_missing,
        "protected path does not exist: $relative_path"))

    tracking = git_tracking(root, relative_path)
    status_entry = git_path_status(root, relative_path)
    if !exists
        return ProtectedSnapshot(relative_path, false, tracking, status_entry, 0, "")
    end
    bytes = read(absolute_path)
    return ProtectedSnapshot(relative_path, true, tracking, status_entry,
        length(bytes), sha256_hex(bytes))
end

function all_protected_snapshots(root::String; require_exists::Bool)::Vector{ProtectedSnapshot}
    return [snapshot_protected(root, path; require_exists=require_exists)
            for path in PROTECTED_PATHS]
end

function snapshot_signature(snapshot::ProtectedSnapshot)
    return (snapshot.path, snapshot.exists, snapshot.tracking, snapshot.status_entry,
        snapshot.byte_count, snapshot.sha256)
end

function snapshots_equal(left::Vector{ProtectedSnapshot}, right::Vector{ProtectedSnapshot})::Bool
    length(left) == length(right) || return false
    return all(snapshot_signature(a) == snapshot_signature(b) for (a, b) in zip(left, right))
end

function atomic_write_bytes(
    path::String,
    bytes::Vector{UInt8};
    before_rename::Function=() -> nothing,
)
    parent = dirname(path)
    isdir(parent) || throw(GuardFailure(:invalid_path, "atomic-write parent is absent: $parent"))
    islink(path) && throw(GuardFailure(:symlink_rejected, "refusing to replace symlink: $path"))
    temporary_path, io = mktemp(parent; cleanup=false)
    try
        write(io, bytes)
        flush(io)
        result = ccall(:fsync, Cint, (Cint,), fd(io))
        result == 0 || throw(GuardFailure(:fsync_failed,
            "fsync failed for temporary publication file"))
        close(io)
        before_rename()
        islink(path) && throw(GuardFailure(:symlink_rejected,
            "destination became a symlink before publication: $path"))
        Base.Filesystem.rename(temporary_path, path)
    finally
        isopen(io) && close(io)
        ispath(temporary_path) && rm(temporary_path; force=true)
    end
    return nothing
end

function serialize_baseline(baseline::WorktreeBaseline)::Vector{UInt8}
    io = IOBuffer()
    write(io, "kind = \"protected_worktree_baseline\"\n")
    write(io, "schema = ", toml_quote(BASELINE_SCHEMA), "\n")
    write(io, "schema_version = 1\n\n")

    write(io, "[repository]\n")
    write(io, "git_dir = ", toml_quote(baseline.repository.git_dir), "\n")
    write(io, "git_toplevel = ", toml_quote(baseline.repository.git_toplevel), "\n")
    write(io, "identity_sha256 = ", toml_quote(baseline.repository.sha256), "\n")
    write(io, "root = ", toml_quote(baseline.repository.root), "\n\n")

    write(io, "[git_status]\n")
    write(io, "command = \"git status --short\"\n")
    write_string_array(io, "entries", baseline.status_entries)
    write(io, "sha256 = ", toml_quote(baseline.status_sha256), "\n\n")

    write(io, "[journal_prefix]\n")
    write(io, "byte_count = $(baseline.journal_prefix_byte_count)\n")
    write(io, "path = ", toml_quote(JOURNAL_PATH), "\n")
    write(io, "sha256 = ", toml_quote(baseline.journal_prefix_sha256), "\n")

    for snapshot in sort(baseline.protected; by=item -> item.path)
        write(io, "\n[[protected]]\n")
        write(io, "byte_count = $(snapshot.byte_count)\n")
        write(io, "exists = ", snapshot.exists ? "true\n" : "false\n")
        write(io, "path = ", toml_quote(snapshot.path), "\n")
        write(io, "sha256 = ", toml_quote(snapshot.sha256), "\n")
        write(io, "status_entry = ", toml_quote(snapshot.status_entry), "\n")
        write(io, "tracking = ", toml_quote(snapshot.tracking), "\n")
    end
    return take!(io)
end

function capture_baseline(root::String, destination::String)::String
    destination_path, relative_destination = resolve_repository_path(root, destination; output=true)
    relative_destination in PROTECTED_PATHS && throw(GuardFailure(:protected_destination,
        "baseline destination overlaps a protected path: $relative_destination"))

    identity = repository_identity(root)
    status_before = git_status_entries(root)
    protected_before = all_protected_snapshots(root; require_exists=true)
    status_after = git_status_entries(root)
    protected_after = all_protected_snapshots(root; require_exists=true)
    status_before == status_after || throw(GuardFailure(:capture_race,
        "git status changed while the baseline was captured"))
    snapshots_equal(protected_before, protected_after) || throw(GuardFailure(:capture_race,
        "protected bytes or statuses changed while the baseline was captured"))

    journal = only(filter(item -> item.path == JOURNAL_PATH, protected_before))
    baseline = WorktreeBaseline(
        identity,
        status_before,
        sha256_hex(canonical_status_bytes(status_before)),
        journal.byte_count,
        journal.sha256,
        protected_before,
    )
    atomic_write_bytes(destination_path, serialize_baseline(baseline))
    return relative_destination
end

function expect_keys(table::AbstractDict, expected::Vector{String}, context::String)
    actual = Set(String(key) for key in keys(table))
    wanted = Set(expected)
    actual == wanted || throw(GuardFailure(:malformed_baseline,
        "$context keys differ; expected $(sort(expected)), got $(sort!(collect(actual)))"))
end

function expect_table(table::AbstractDict, key::String, context::String)::AbstractDict
    value = get(table, key, nothing)
    value isa AbstractDict || throw(GuardFailure(:malformed_baseline,
        "$context.$key must be a TOML table"))
    return value
end

function expect_string(table::AbstractDict, key::String, context::String)::String
    value = get(table, key, nothing)
    value isa String || throw(GuardFailure(:malformed_baseline,
        "$context.$key must be a string"))
    return value
end

function expect_bool(table::AbstractDict, key::String, context::String)::Bool
    value = get(table, key, nothing)
    value isa Bool || throw(GuardFailure(:malformed_baseline,
        "$context.$key must be a boolean"))
    return value
end

function expect_integer(table::AbstractDict, key::String, context::String)::Int
    value = get(table, key, nothing)
    (value isa Integer && !(value isa Bool)) || throw(GuardFailure(:malformed_baseline,
        "$context.$key must be an integer"))
    try
        return Int(value)
    catch
        throw(GuardFailure(:malformed_baseline, "$context.$key is outside Int range"))
    end
end

function expect_sha256(value::String, context::String; allow_empty::Bool=false)::String
    allow_empty && isempty(value) && return value
    occursin(r"^[0-9a-f]{64}$", value) || throw(GuardFailure(:malformed_baseline,
        "$context must be a lowercase SHA-256"))
    return value
end

function parse_baseline(bytes::Vector{UInt8})::WorktreeBaseline
    document = try
        TOML.parse(String(copy(bytes)))
    catch err
        throw(GuardFailure(:malformed_baseline, "TOML parse failed: $(sprint(showerror, err))"))
    end
    expect_keys(document,
        ["git_status", "journal_prefix", "kind", "protected", "repository", "schema", "schema_version"],
        "baseline")
    expect_string(document, "kind", "baseline") == "protected_worktree_baseline" ||
        throw(GuardFailure(:malformed_baseline, "baseline.kind is invalid"))
    expect_string(document, "schema", "baseline") == BASELINE_SCHEMA ||
        throw(GuardFailure(:malformed_baseline, "baseline.schema is invalid"))
    expect_integer(document, "schema_version", "baseline") == 1 ||
        throw(GuardFailure(:malformed_baseline, "baseline.schema_version is invalid"))

    repository = expect_table(document, "repository", "baseline")
    expect_keys(repository, ["git_dir", "git_toplevel", "identity_sha256", "root"], "repository")
    identity = RepositoryIdentity(
        expect_string(repository, "root", "repository"),
        expect_string(repository, "git_toplevel", "repository"),
        expect_string(repository, "git_dir", "repository"),
        expect_sha256(expect_string(repository, "identity_sha256", "repository"),
            "repository.identity_sha256"),
    )
    identity_material = "stmfit-structured-worktree-root-v1\n" *
                        "root=$(identity.root)\n" *
                        "git_toplevel=$(identity.git_toplevel)\n" *
                        "git_dir=$(identity.git_dir)\n"
    sha256_hex(string_bytes(identity_material)) == identity.sha256 ||
        throw(GuardFailure(:malformed_baseline, "repository identity digest is inconsistent"))

    status = expect_table(document, "git_status", "baseline")
    expect_keys(status, ["command", "entries", "sha256"], "git_status")
    expect_string(status, "command", "git_status") == "git status --short" ||
        throw(GuardFailure(:malformed_baseline, "git_status.command is invalid"))
    status_value = get(status, "entries", nothing)
    status_value isa Vector || throw(GuardFailure(:malformed_baseline,
        "git_status.entries must be an array"))
    all(item -> item isa String, status_value) || throw(GuardFailure(:malformed_baseline,
        "git_status.entries must contain only strings"))
    status_entries = String[String(item) for item in status_value]
    status_entries == sort(unique(status_entries)) || throw(GuardFailure(:malformed_baseline,
        "git_status.entries must be sorted and unique"))
    status_sha256 = expect_sha256(expect_string(status, "sha256", "git_status"),
        "git_status.sha256")
    sha256_hex(canonical_status_bytes(status_entries)) == status_sha256 ||
        throw(GuardFailure(:malformed_baseline, "git status digest is inconsistent"))

    prefix = expect_table(document, "journal_prefix", "baseline")
    expect_keys(prefix, ["byte_count", "path", "sha256"], "journal_prefix")
    prefix_path = expect_string(prefix, "path", "journal_prefix")
    prefix_path == JOURNAL_PATH || throw(GuardFailure(:malformed_baseline,
        "journal_prefix.path must be $JOURNAL_PATH"))
    prefix_byte_count = expect_integer(prefix, "byte_count", "journal_prefix")
    prefix_byte_count >= 0 || throw(GuardFailure(:malformed_baseline,
        "journal_prefix.byte_count must be nonnegative"))
    prefix_sha256 = expect_sha256(expect_string(prefix, "sha256", "journal_prefix"),
        "journal_prefix.sha256")

    protected_value = get(document, "protected", nothing)
    protected_value isa Vector || throw(GuardFailure(:malformed_baseline,
        "protected must be an array of tables"))
    length(protected_value) == length(PROTECTED_PATHS) ||
        throw(GuardFailure(:malformed_baseline,
            "baseline must name exactly $(length(PROTECTED_PATHS)) protected paths"))
    snapshots = ProtectedSnapshot[]
    for (index, item) in enumerate(protected_value)
        item isa AbstractDict || throw(GuardFailure(:malformed_baseline,
            "protected[$index] must be a table"))
        expect_keys(item,
            ["byte_count", "exists", "path", "sha256", "status_entry", "tracking"],
            "protected[$index]")
        path = expect_string(item, "path", "protected[$index]")
        exists = expect_bool(item, "exists", "protected[$index]")
        exists || throw(GuardFailure(:malformed_baseline,
            "captured protected path must exist: $path"))
        tracking = expect_string(item, "tracking", "protected[$index]")
        tracking in ("tracked", "untracked") || throw(GuardFailure(:malformed_baseline,
            "protected[$index].tracking is invalid"))
        byte_count = expect_integer(item, "byte_count", "protected[$index]")
        byte_count >= 0 || throw(GuardFailure(:malformed_baseline,
            "protected[$index].byte_count must be nonnegative"))
        digest = expect_sha256(expect_string(item, "sha256", "protected[$index]"),
            "protected[$index].sha256")
        push!(snapshots, ProtectedSnapshot(path, true, tracking,
            expect_string(item, "status_entry", "protected[$index]"), byte_count, digest))
    end
    sort!(snapshots; by=item -> item.path)
    [item.path for item in snapshots] == PROTECTED_PATHS ||
        throw(GuardFailure(:malformed_baseline,
            "baseline must name exactly the frozen protected paths without duplicates"))
    journal = only(filter(item -> item.path == JOURNAL_PATH, snapshots))
    (journal.byte_count == prefix_byte_count && journal.sha256 == prefix_sha256) ||
        throw(GuardFailure(:malformed_baseline,
            "journal prefix must equal the complete captured journal"))

    return WorktreeBaseline(identity, status_entries, status_sha256,
        prefix_byte_count, prefix_sha256, snapshots)
end

function assert_repository_identity(captured::RepositoryIdentity, current::RepositoryIdentity)
    captured.root == current.root || throw(GuardFailure(:repository_identity_changed,
        "repository root differs from the captured baseline"))
    captured.git_toplevel == current.git_toplevel ||
        throw(GuardFailure(:repository_identity_changed, "git top-level differs from baseline"))
    captured.git_dir == current.git_dir ||
        throw(GuardFailure(:repository_identity_changed, "git directory differs from baseline"))
    captured.sha256 == current.sha256 ||
        throw(GuardFailure(:repository_identity_changed, "repository identity digest differs"))
end

function assert_status_and_tracking(captured::ProtectedSnapshot, current::ProtectedSnapshot)
    captured.tracking == current.tracking || throw(GuardFailure(:dirty_status_changed,
        "tracking state changed for $(captured.path): $(captured.tracking) -> $(current.tracking)"))
    captured.status_entry == current.status_entry || throw(GuardFailure(:dirty_status_changed,
        "git status entry changed for $(captured.path): $(repr(captured.status_entry)) -> $(repr(current.status_entry))"))
end

function assert_exact_bytes(captured::ProtectedSnapshot, current::ProtectedSnapshot)
    current.exists || throw(GuardFailure(:protected_missing,
        "protected path was deleted: $(captured.path)"))
    captured.byte_count == current.byte_count || throw(GuardFailure(:protected_changed,
        "byte count changed for $(captured.path): $(captured.byte_count) -> $(current.byte_count)"))
    captured.sha256 == current.sha256 || throw(GuardFailure(:protected_changed,
        "SHA-256 changed for $(captured.path)"))
end

function verify_baseline(
    root::String,
    baseline::WorktreeBaseline,
    mode::Symbol,
)::Vector{ProtectedSnapshot}
    mode in (:pre_write, :allow_journal_append) ||
        throw(GuardFailure(:invalid_mode, "unsupported verification mode: $mode"))
    assert_repository_identity(baseline.repository, repository_identity(root))
    current = all_protected_snapshots(root; require_exists=false)
    for (captured, observed) in zip(baseline.protected, current)
        captured.path == observed.path || throw(GuardFailure(:internal_error,
            "protected path ordering changed"))
        assert_status_and_tracking(captured, observed)
        if mode == :pre_write || captured.path != JOURNAL_PATH
            assert_exact_bytes(captured, observed)
        else
            observed.exists || throw(GuardFailure(:protected_missing,
                "journal was deleted"))
            journal_bytes = read(protected_absolute_path(root, JOURNAL_PATH))
            length(journal_bytes) >= baseline.journal_prefix_byte_count ||
                throw(GuardFailure(:journal_truncated,
                    "journal is shorter than its captured prefix"))
            prefix = baseline.journal_prefix_byte_count == 0 ? UInt8[] :
                journal_bytes[1:baseline.journal_prefix_byte_count]
            sha256_hex(prefix) == baseline.journal_prefix_sha256 ||
                throw(GuardFailure(:journal_prefix_changed,
                    "journal captured prefix bytes changed"))
        end
    end
    return current
end

function mode_name(mode::Symbol)::String
    mode == :pre_write && return "pre-write"
    mode == :allow_journal_append && return "allow-journal-append"
    return String(mode)
end

function serialize_receipt(
    identity::RepositoryIdentity,
    mode::Symbol,
    baseline_path::String,
    baseline_sha256::String,
    decision::String,
    reason_code::Symbol,
    reason::String,
    status_entries::Vector{String},
)::Vector{UInt8}
    io = IOBuffer()
    write(io, "baseline_path = ", toml_quote(baseline_path), "\n")
    write(io, "baseline_sha256 = ", toml_quote(baseline_sha256), "\n")
    write(io, "decision = ", toml_quote(decision), "\n")
    write(io, "kind = \"protected_worktree_guard_receipt\"\n")
    write(io, "mode = ", toml_quote(mode_name(mode)), "\n")
    write(io, "reason = ", toml_quote(reason), "\n")
    write(io, "reason_code = ", toml_quote(String(reason_code)), "\n")
    write(io, "repository_identity_sha256 = ", toml_quote(identity.sha256), "\n")
    write(io, "repository_root = ", toml_quote(identity.root), "\n")
    write(io, "schema = ", toml_quote(RECEIPT_SCHEMA), "\n")
    write(io, "schema_version = 1\n")
    write_string_array(io, "status_entries", status_entries)
    write(io, "status_sha256 = ",
        toml_quote(sha256_hex(canonical_status_bytes(status_entries))), "\n")
    write_string_array(io, "verified_paths", copy(PROTECTED_PATHS))
    return take!(io)
end

function terminal_receipt_decision(path::String)::Union{Nothing,String}
    (isfile(path) && !islink(path)) || return nothing
    document = try
        TOML.parsefile(path)
    catch
        return nothing
    end
    get(document, "kind", nothing) == "protected_worktree_guard_receipt" || return nothing
    get(document, "schema", nothing) == RECEIPT_SCHEMA || return nothing
    get(document, "schema_version", nothing) == 1 || return nothing
    decision = get(document, "decision", nothing)
    decision in ("PASS", "BLOCKED") || return nothing
    mode = get(document, "mode", nothing)
    mode in ("pre-write", "allow-journal-append") || return nothing
    baseline_sha = get(document, "baseline_sha256", nothing)
    baseline_sha isa String || return nothing
    decision == "PASS" && !occursin(r"^[0-9a-f]{64}$", baseline_sha) && return nothing
    return String(decision)
end

function guard_failure(err)::GuardFailure
    err isa GuardFailure && return err
    return GuardFailure(:internal_error, sprint(showerror, err))
end

function safe_status_entries(root::String)::Vector{String}
    try
        return git_status_entries(root)
    catch
        return String[]
    end
end

function verify_with_receipt(
    root::String,
    mode::Symbol,
    baseline_argument::String,
    receipt_argument::String,
)
    identity = repository_identity(root)
    receipt_path, receipt_relative = resolve_repository_path(root, receipt_argument; output=true)
    receipt_relative in PROTECTED_PATHS && throw(GuardFailure(:protected_destination,
        "receipt destination overlaps a protected path: $receipt_relative"))

    baseline_relative = baseline_argument
    baseline_digest = ""
    baseline_path = ""
    try
        baseline_path, baseline_relative = resolve_repository_path(root, baseline_argument;
            must_exist=true)
        baseline_path == receipt_path && throw(GuardFailure(:path_collision,
            "baseline and receipt paths must differ"))
        baseline_bytes = read(baseline_path)
        baseline_digest = sha256_hex(baseline_bytes)
        baseline = parse_baseline(baseline_bytes)
        verify_baseline(root, baseline, mode)

        reread_baseline = read(baseline_path)
        reread_baseline == baseline_bytes || throw(GuardFailure(:baseline_changed_during_check,
            "baseline bytes changed during verification"))
        verify_baseline(root, baseline, mode)
        status_entries = safe_status_entries(root)
        receipt_bytes = serialize_receipt(identity, mode, baseline_relative, baseline_digest,
            "PASS", :ok, "all protected invariants satisfied", status_entries)
        atomic_write_bytes(receipt_path, receipt_bytes)
        terminal_receipt_decision(receipt_path) == "PASS" ||
            throw(GuardFailure(:receipt_publication_failed,
                "published PASS receipt is not terminal and parseable"))
        return (ok=true, code=:ok, message="all protected invariants satisfied",
            receipt=receipt_relative)
    catch err
        failure = guard_failure(err)
        if !isempty(baseline_path) && baseline_path == receipt_path
            throw(failure)
        end
        status_entries = safe_status_entries(root)
        receipt_bytes = serialize_receipt(identity, mode, baseline_relative, baseline_digest,
            "BLOCKED", failure.code, failure.message, status_entries)
        atomic_write_bytes(receipt_path, receipt_bytes)
        terminal_receipt_decision(receipt_path) == "BLOCKED" ||
            throw(GuardFailure(:receipt_publication_failed,
                "published BLOCKED receipt is not terminal and parseable"))
        return (ok=false, code=failure.code, message=failure.message,
            receipt=receipt_relative)
    end
end

function assign_once(current::Union{Nothing,String}, value::String, option::String)::String
    current === nothing || throw(GuardFailure(:cli_error, "$option may be supplied only once"))
    return value
end

function option_value(arguments::Vector{String}, index::Int, option::String)::Tuple{String,Int}
    index < length(arguments) || throw(GuardFailure(:cli_error, "$option requires a path"))
    value = arguments[index + 1]
    startswith(value, "--") && throw(GuardFailure(:cli_error, "$option requires a path"))
    isempty(value) && throw(GuardFailure(:cli_error, "$option requires a nonempty path"))
    return value, index + 2
end

function parse_cli(arguments::Vector{String})::CliOptions
    mode::Union{Nothing,Symbol} = nothing
    target::Union{Nothing,String} = nothing
    baseline::Union{Nothing,String} = nothing
    receipt::Union{Nothing,String} = nothing
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--capture"
            mode === nothing || throw(GuardFailure(:cli_error, "exactly one mode is required"))
            mode = :capture
            value, index = option_value(arguments, index, argument)
            target = assign_once(target, value, argument)
        elseif argument == "--pre-write"
            mode === nothing || throw(GuardFailure(:cli_error, "exactly one mode is required"))
            mode = :pre_write
            index += 1
        elseif argument == "--allow-journal-append"
            mode === nothing || throw(GuardFailure(:cli_error, "exactly one mode is required"))
            mode = :allow_journal_append
            index += 1
        elseif argument == "--self-test"
            mode === nothing || throw(GuardFailure(:cli_error, "exactly one mode is required"))
            mode = :self_test
            index += 1
        elseif argument == "--baseline"
            value, index = option_value(arguments, index, argument)
            baseline = assign_once(baseline, value, argument)
        elseif argument == "--receipt"
            value, index = option_value(arguments, index, argument)
            receipt = assign_once(receipt, value, argument)
        else
            throw(GuardFailure(:cli_error, "unknown argument: $argument"))
        end
    end

    mode === nothing && throw(GuardFailure(:cli_error, "one mode is required"))
    if mode == :capture
        target === nothing && throw(GuardFailure(:cli_error, "--capture requires a path"))
        (baseline === nothing && receipt === nothing) || throw(GuardFailure(:cli_error,
            "--capture does not accept --baseline or --receipt"))
    elseif mode in (:pre_write, :allow_journal_append)
        baseline === nothing && throw(GuardFailure(:cli_error, "--baseline is required"))
        receipt === nothing && throw(GuardFailure(:cli_error, "--receipt is required"))
        target === nothing || throw(GuardFailure(:cli_error, "verification mode does not accept --capture"))
    else
        (target === nothing && baseline === nothing && receipt === nothing) ||
            throw(GuardFailure(:cli_error, "--self-test accepts no other options"))
    end
    return CliOptions(mode, target, baseline, receipt)
end

function mutate_byte(bytes::Vector{UInt8}, index::Int)::Vector{UInt8}
    isempty(bytes) && error("fixture unexpectedly empty")
    mutated = copy(bytes)
    mutated[index] = xor(mutated[index], 0x01)
    return mutated
end

function with_fixture(function_body::Function, source_root::String, cleanup_roots::Vector{String})
    fixture_root = mktempdir(; prefix="stmfit-guard-selftest-")
    push!(cleanup_roots, fixture_root)
    try
        for relative_path in PROTECTED_PATHS
            source = protected_absolute_path(source_root, relative_path)
            destination = joinpath(fixture_root, split(relative_path, '/')...)
            mkpath(dirname(destination))
            write(destination, read(source))
        end
        write(joinpath(fixture_root, ".gitignore"), ".omo/\n")
        run(Cmd(["git", "-C", fixture_root, "init", "-q"]))
        mkpath(joinpath(fixture_root, ".omo", "evidence"))
        return function_body(fixture_root)
    finally
        rm(fixture_root; recursive=true, force=true)
    end
end

function baseline_and_receipt_paths(root::String)
    evidence = joinpath(root, ".omo", "evidence")
    return joinpath(evidence, "baseline.toml"), joinpath(evidence, "receipt.toml")
end

function receipt_reason_code(path::String)::String
    document = TOML.parsefile(path)
    return String(document["reason_code"])
end

function run_self_tests(source_root::String)
    source_hashes_before = Dict(path => sha256_hex(read(protected_absolute_path(source_root, path)))
                                for path in PROTECTED_PATHS)
    cleanup_roots = String[]

    @testset "structured protected-worktree guard" begin
        @testset "deterministic capture and exact pre-write" begin
            with_fixture(source_root, cleanup_roots) do root
                baseline, receipt = baseline_and_receipt_paths(root)
                repeat_baseline = joinpath(dirname(baseline), "baseline-repeat.toml")
                capture_baseline(root, baseline)
                capture_baseline(root, repeat_baseline)
                @test read(baseline) == read(repeat_baseline)
                parsed = parse_baseline(read(baseline))
                @test [item.path for item in parsed.protected] == PROTECTED_PATHS
                @test parsed.journal_prefix_byte_count ==
                      length(read(protected_absolute_path(root, JOURNAL_PATH)))
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test outcome.ok
                @test terminal_receipt_decision(receipt) == "PASS"
                first_receipt = read(receipt)
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test outcome.ok
                @test read(receipt) == first_receipt
            end
        end

        @testset "replacement deletion and byte mutation matrix" begin
            with_fixture(source_root, cleanup_roots) do root
                baseline, receipt = baseline_and_receipt_paths(root)
                capture_baseline(root, baseline)
                originals = Dict(path => read(protected_absolute_path(root, path))
                                 for path in PROTECTED_PATHS)

                for relative_path in PROTECTED_PATHS
                    destination = protected_absolute_path(root, relative_path)
                    for position in (1, length(originals[relative_path]))
                        write(destination, mutate_byte(originals[relative_path], position))
                        outcome = verify_with_receipt(root, :allow_journal_append,
                            baseline, receipt)
                        @test !outcome.ok
                        @test terminal_receipt_decision(receipt) == "BLOCKED"
                        write(destination, originals[relative_path])
                    end
                end

                for relative_path in PROTECTED_PATHS
                    destination = protected_absolute_path(root, relative_path)
                    replacement = fill(UInt8('X'), length(originals[relative_path]))
                    replacement == originals[relative_path] && (replacement[1] = UInt8('Y'))
                    temporary = joinpath(dirname(destination), ".replacement")
                    write(temporary, replacement)
                    Base.Filesystem.rename(temporary, destination)
                    outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                    @test !outcome.ok
                    @test receipt_reason_code(receipt) == "protected_changed"
                    write(destination, originals[relative_path])
                end

                for relative_path in PROTECTED_PATHS
                    destination = protected_absolute_path(root, relative_path)
                    rm(destination)
                    outcome = verify_with_receipt(root, :allow_journal_append,
                        baseline, receipt)
                    @test !outcome.ok
                    @test terminal_receipt_decision(receipt) == "BLOCKED"
                    write(destination, originals[relative_path])
                end

                journal = protected_absolute_path(root, JOURNAL_PATH)
                write(journal, originals[JOURNAL_PATH][1:end-1])
                outcome = verify_with_receipt(root, :allow_journal_append, baseline, receipt)
                @test !outcome.ok
                @test receipt_reason_code(receipt) == "journal_truncated"
                write(journal, originals[JOURNAL_PATH])

                prefix_edit = copy(originals[JOURNAL_PATH])
                middle = max(1, length(prefix_edit) ÷ 2)
                prefix_edit[middle] = xor(prefix_edit[middle], 0x01)
                write(journal, prefix_edit)
                outcome = verify_with_receipt(root, :allow_journal_append, baseline, receipt)
                @test !outcome.ok
                @test receipt_reason_code(receipt) == "journal_prefix_changed"
                write(journal, originals[JOURNAL_PATH])

                open(journal, "a") do io
                    write(io, "\nself-test journal suffix\n")
                end
                strict = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test !strict.ok
                append_only = verify_with_receipt(root, :allow_journal_append, baseline, receipt)
                @test append_only.ok
                @test terminal_receipt_decision(receipt) == "PASS"
                write(journal, originals[JOURNAL_PATH])

                final_check = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test final_check.ok
            end
        end

        @testset "CLI and malformed baseline rejection" begin
            for arguments in (
                String[],
                ["--unknown"],
                ["--capture"],
                ["--pre-write"],
                ["--pre-write", "--baseline", "x"],
                ["--self-test", "--receipt", "x"],
                ["--self-test", "--pre-write"],
            )
                @test_throws GuardFailure parse_cli(arguments)
            end
            parsed = parse_cli(["--pre-write", "--baseline", "a", "--receipt", "b"])
            @test parsed.mode == :pre_write

            with_fixture(source_root, cleanup_roots) do root
                baseline, receipt = baseline_and_receipt_paths(root)
                capture_baseline(root, baseline)
                valid_bytes = read(baseline)

                write(baseline, "broken = [")
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test !outcome.ok
                @test receipt_reason_code(receipt) == "malformed_baseline"

                valid_text = String(copy(valid_bytes))
                last_header = findlast("[[protected]]", valid_text)
                @test last_header !== nothing
                missing_text = valid_text[1:first(last_header)-1]
                write(baseline, missing_text)
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test !outcome.ok
                @test receipt_reason_code(receipt) == "malformed_baseline"

                duplicate_block = valid_text[first(last_header):end]
                write(baseline, valid_text * "\n" * duplicate_block)
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test !outcome.ok
                @test receipt_reason_code(receipt) == "malformed_baseline"

                write(baseline, valid_bytes)
                stale_path = protected_absolute_path(root, "test/detect_missed_lobes.py")
                stale_bytes = read(stale_path)
                write(stale_path, mutate_byte(stale_bytes, 1))
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test !outcome.ok
                @test terminal_receipt_decision(receipt) == "BLOCKED"
            end
        end

        @testset "path confinement and symlink rejection" begin
            with_fixture(source_root, cleanup_roots) do root
                baseline, receipt = baseline_and_receipt_paths(root)
                @test_throws GuardFailure capture_baseline(root, "../escaped.toml")
                capture_baseline(root, baseline)

                baseline_link = joinpath(dirname(baseline), "baseline-link.toml")
                symlink(baseline, baseline_link)
                outcome = verify_with_receipt(root, :pre_write, baseline_link, receipt)
                @test !outcome.ok
                @test receipt_reason_code(receipt) == "symlink_rejected"

                outside = mktempdir(; prefix="stmfit-guard-outside-")
                try
                    symlink(outside, joinpath(root, "outside-link"))
                    @test_throws GuardFailure capture_baseline(root,
                        "outside-link/escaped-baseline.toml")

                    target = joinpath(outside, "receipt-target.toml")
                    write(target, "unchanged")
                    receipt_link = joinpath(dirname(receipt), "receipt-link.toml")
                    symlink(target, receipt_link)
                    @test_throws GuardFailure verify_with_receipt(root, :pre_write,
                        baseline, receipt_link)
                    @test read(target, String) == "unchanged"
                finally
                    rm(outside; recursive=true, force=true)
                end
            end
        end

        @testset "interrupted and misleading receipt state" begin
            with_fixture(source_root, cleanup_roots) do root
                baseline, receipt = baseline_and_receipt_paths(root)
                capture_baseline(root, baseline)
                outcome = verify_with_receipt(root, :pre_write, baseline, receipt)
                @test outcome.ok
                valid_receipt = read(receipt)

                rm(receipt)
                partial = receipt * ".partial"
                write(partial, valid_receipt)
                @test terminal_receipt_decision(receipt) === nothing
                rm(partial)
                @test !ispath(partial)

                fake_log = joinpath(dirname(receipt), "success.log")
                write(fake_log, "PASS: everything succeeded\n")
                @test terminal_receipt_decision(fake_log) === nothing

                old_bytes = string_bytes("old terminal bytes\n")
                write(receipt, old_bytes)
                directory_before = Set(readdir(dirname(receipt)))
                @test_throws ErrorException atomic_write_bytes(receipt,
                    string_bytes("new terminal bytes\n");
                    before_rename=() -> error("simulated interruption"))
                @test read(receipt) == old_bytes
                @test Set(readdir(dirname(receipt))) == directory_before
            end
        end
    end

    source_hashes_after = Dict(path => sha256_hex(read(protected_absolute_path(source_root, path)))
                               for path in PROTECTED_PATHS)
    @test source_hashes_after == source_hashes_before
    @test all(path -> !ispath(path), cleanup_roots)
    return nothing
end

function main(arguments::Vector{String})::Int
    root = realpath(joinpath(@__DIR__, ".."))
    options = try
        parse_cli(arguments)
    catch err
        failure = guard_failure(err)
        println(stderr, "BLOCKED [$(failure.code)]: $(failure.message)")
        return 2
    end

    try
        if options.mode == :self_test
            run_self_tests(root)
            println("SELF-TEST PASS")
            return 0
        elseif options.mode == :capture
            relative = capture_baseline(root, something(options.target))
            println("CAPTURED $relative")
            return 0
        else
            outcome = verify_with_receipt(root, options.mode,
                something(options.baseline), something(options.receipt))
            if outcome.ok
                println("PASS receipt=$(outcome.receipt)")
                return 0
            end
            println(stderr, "BLOCKED [$(outcome.code)]: $(outcome.message)")
            return 1
        end
    catch err
        failure = guard_failure(err)
        println(stderr, "BLOCKED [$(failure.code)]: $(failure.message)")
        return 1
    end
end

exit(main(copy(ARGS)))
