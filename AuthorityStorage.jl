module AuthorityStorage

"""
Storage-v2's small, deliberately boring filesystem boundary.

The module does not discover a checkout, infer a project, or read the process
environment.  Callers provide every authority binding explicitly.  The
builders use this module for staging, manifests, publication, and the
execution-image move/restore contract.
"""

using SHA
using TOML

export ManifestEntry, Manifest, ExecutionImage, AmbiguousPublicationError,
       COMMON_CATEGORIES, IMAGE_LABELS, receipt_path, sha256_file,
       tree_digest, snapshot, write_manifest!, read_manifest, verify_manifest,
       verify_source_provenance, publish_staged!, write_external_receipt!,
       execution_image, relocate_no_replace, verify_relocated,
       restore_relocated!, verify_image_set, build_common_authority,
       build_execution_images, sanitize_environment, canonical_command,
       materialize_common_payload!

const COMMON_CATEGORIES = ("config", "package", "test", "audit", "harness",
                           "source", "closure")
const IMAGE_LABELS = ("t10", "t11", "t12", "reverse")
const FORMAT = "stmfit-storage-v2"
const SELF_SHA = "self-reference"

struct ManifestEntry
    path::String
    kind::Symbol
    sha256::Union{Nothing,String}
    bytes::Int64
    mode::UInt64
    nlink::UInt64
    device::UInt64
    inode::UInt64
    children::Vector{String}
end

struct Manifest
    entries::Vector{ManifestEntry}
    digest::String
end

mutable struct ExecutionImage
    root::String
    original_root::String
    label::String
    common_receipt_sha256::String
    prestate_slot::String
    fixed_lexical_path::String
    manifest::Manifest
end

struct AmbiguousPublicationError <: Exception
    root::String
    receipt::String
    cause
end

Base.showerror(io::IO, e::AmbiguousPublicationError) = print(
    io, "ambiguous publication (root was renamed but receipt was not " *
    "published): ", e.root, "; receipt: ", e.receipt, "; cause: ", e.cause)

struct StorageError <: Exception
    message::String
end
Base.showerror(io::IO, e::StorageError) = print(io, e.message)

_err(msg) = throw(StorageError(String(msg)))

"Return the external receipt name.  It is never inside the authority root."
receipt_path(root::AbstractString) = String(root) * ".receipt.toml"

_normroot(p::AbstractString) = normpath(abspath(String(p)))
_rel(p) = replace(String(p), '\\' => '/')

function _permission_mode(mode::Integer)
    # Keep file type bits: they are part of the identity we verify, while the
    # permission bits remain directly visible in the resulting manifest.
    UInt64(mode)
end

function _lstat(path::AbstractString)
    try
        Base.Filesystem.lstat(String(path))
    catch e
        _err("cannot lstat $(path): $(sprint(showerror, e))")
    end
end

function _kind(st)
    if Base.Filesystem.isfile(st)
        return :file
    elseif Base.Filesystem.isdir(st)
        return :dir
    end
    _err("unsupported filesystem object (symlink, device, fifo, or socket): " *
         "object mode $(st.mode)")
end

function _sha_bytes(bytes::Vector{UInt8})
    bytes2hex(SHA.sha256(bytes))
end

function sha256_file(path::AbstractString)
    st = _lstat(path)
    _kind(st) === :file || _err("expected a regular file: $(path)")
    open(String(path), "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _digest_tree(path::AbstractString)
    st = _lstat(path)
    kind = _kind(st)
    kind === :file && return sha256_file(path)
    rows = String[]
    for name in sort(readdir(String(path)))
        child = joinpath(String(path), name)
        cst = _lstat(child)
        ck = _kind(cst)
        if ck === :file
            push!(rows, "file\t$(_rel(name))\t$(sha256_file(child))\t$(cst.size)")
        else
            push!(rows, "dir\t$(_rel(name))\t$(_digest_tree(child))")
        end
    end
    _sha_bytes(Vector{UInt8}(codeunits(join(rows, "\n"))))
end

tree_digest(path::AbstractString) = _digest_tree(path)

function _walk(path::AbstractString, rel::String=".")
    st = _lstat(path)
    kind = _kind(st)
    out = Tuple{String,String,Any}[(rel, String(kind), st)]
    if kind === :dir
        for name in sort(readdir(String(path)))
            childrel = rel == "." ? name : joinpath(rel, name)
            append!(out, _walk(joinpath(String(path), name), _rel(childrel)))
        end
    end
    out
end

function _direct_children(path::AbstractString)
    sort!([_rel(name) for name in readdir(String(path))])
end

function snapshot(root::AbstractString; manifest_name::AbstractString="manifest.toml")
    root = _normroot(root)
    ispath(root) || _err("root does not exist: $(root)")
    rows = _walk(root)
    entries = ManifestEntry[]
    for (rel, kind_string, st) in rows
        kind = Symbol(kind_string)
        if kind === :file
            # A manifest cannot contain its own cryptographic hash.  It does
            # contain all of its other metadata and a typed self-reference;
            # the external receipt authenticates the complete manifest bytes.
            digest = rel == manifest_name ? SELF_SHA : sha256_file(joinpath(root, rel))
            push!(entries, ManifestEntry(rel, :file, digest, Int64(st.size),
                                         _permission_mode(st.mode), UInt64(st.nlink),
                                         UInt64(st.device), UInt64(st.inode), String[]))
        else
            push!(entries, ManifestEntry(rel, :dir, nothing, Int64(0),
                                         _permission_mode(st.mode), UInt64(st.nlink),
                                         UInt64(st.device), UInt64(st.inode),
                                         _direct_children(joinpath(root, rel == "." ? "" : rel))))
        end
    end
    sort!(entries, by=e -> e.path)
    Manifest(entries, "")
end

function _entry_dict(e::ManifestEntry)
    d = Dict{String,Any}(
        "path" => e.path,
        "kind" => String(e.kind),
        "bytes" => e.bytes,
        "mode" => Int(e.mode),
        "nlink" => Int(e.nlink),
        "device" => Int(e.device),
        "inode" => Int(e.inode),
        "children" => e.children,
    )
    d["sha256"] = e.sha256 === nothing ? "" : e.sha256
    d
end

function _manifest_data(m::Manifest)
    Dict{String,Any}(
        "format" => FORMAT,
        "entry_count" => length(m.entries),
        "entries" => [_entry_dict(e) for e in m.entries],
    )
end

function _manifest_bytes(m::Manifest)
    io = IOBuffer()
    TOML.print(io, _manifest_data(m); sorted=true)
    Vector{UInt8}(take!(io))
end

function _with_digest(m::Manifest, bytes::Vector{UInt8})
    Manifest(m.entries, _sha_bytes(bytes))
end

function _write_bytes(path::AbstractString, bytes::Vector{UInt8}; mode::Union{Nothing,Integer}=nothing)
    open(String(path), "w") do io
        write(io, bytes)
    end
    mode === nothing || chmod(String(path), mode)
    String(path)
end

function _manifest_entry_with_actual_self(m::Manifest, root::String)
    # The manifest's inode is allocated by the first write.  Its size is
    # solved to a fixed point because the size itself is recorded.
    entries = copy(m.entries)
    idx = findfirst(e -> e.path == "manifest.toml", entries)
    idx === nothing && _err("snapshot did not include manifest.toml")
    provisional = Manifest(entries, "")
    _write_bytes(joinpath(root, "manifest.toml"), _manifest_bytes(provisional); mode=0o444)
    for _ in 1:8
        st = _lstat(joinpath(root, "manifest.toml"))
        old = entries[idx]
        entries[idx] = ManifestEntry(old.path, old.kind, SELF_SHA, Int64(st.size),
                                     _permission_mode(st.mode), UInt64(st.nlink),
                                     UInt64(st.device), UInt64(st.inode), old.children)
        current = Manifest(entries, "")
        bytes = _manifest_bytes(current)
        _write_bytes(joinpath(root, "manifest.toml"), bytes; mode=0o444)
        nst = _lstat(joinpath(root, "manifest.toml"))
        if Int64(nst.size) == Int64(length(bytes))
            # Capture all metadata once more; rewriting keeps inode and link
            # count stable, but this also protects against unusual filesystems.
            entries[idx] = ManifestEntry(old.path, old.kind, SELF_SHA,
                                         Int64(nst.size), _permission_mode(nst.mode),
                                         UInt64(nst.nlink), UInt64(nst.device),
                                         UInt64(nst.inode), old.children)
            final = Manifest(entries, "")
            final_bytes = _manifest_bytes(final)
            if length(final_bytes) == Int(nst.size)
                _write_bytes(joinpath(root, "manifest.toml"), final_bytes; mode=0o444)
                return _with_digest(final, final_bytes)
            end
        end
    end
    _err("manifest self-metadata did not reach a fixed point")
end

function write_manifest!(root::AbstractString)
    root = _normroot(root)
    m = snapshot(root)
    _manifest_entry_with_actual_self(m, root)
end

function read_manifest(path::AbstractString)
    data = TOML.parsefile(String(path))
    get(data, "format", nothing) == FORMAT || _err("wrong manifest format: $(path)")
    raw = get(data, "entries", nothing)
    raw isa Vector || _err("manifest entries are not an array: $(path)")
    entries = ManifestEntry[]
    for x in raw
        kind = Symbol(String(x["kind"]))
        sha = String(get(x, "sha256", ""))
        push!(entries, ManifestEntry(String(x["path"]), kind,
                                     isempty(sha) ? nothing : sha,
                                     Int64(x["bytes"]), UInt64(x["mode"]),
                                     UInt64(x["nlink"]), UInt64(x["device"]),
                                     UInt64(x["inode"]),
                                     String.(get(x, "children", String[]))))
    end
    bytes = read(String(path))
    Manifest(entries, _sha_bytes(bytes))
end

function _same_entry(expected::ManifestEntry, actual::ManifestEntry; check_sha=true)
    expected.path == actual.path || return false
    expected.kind == actual.kind || return false
    expected.bytes == actual.bytes || return false
    expected.mode == actual.mode || return false
    expected.nlink == actual.nlink || return false
    expected.device == actual.device || return false
    expected.inode == actual.inode || return false
    expected.children == actual.children || return false
    (!check_sha || expected.sha256 == actual.sha256) || return false
    true
end

function verify_manifest(root::AbstractString, expected::Union{Nothing,Manifest}=nothing)
    root = _normroot(root)
    path = joinpath(root, "manifest.toml")
    ispath(path) || _err("manifest missing: $(path)")
    recorded = read_manifest(path)
    current = snapshot(root)
    # The current snapshot has the self-reference marker too.  The manifest
    # bytes are authenticated by the receipt or by `recorded.digest` supplied
    # by a caller.
    length(recorded.entries) == length(current.entries) || _err("manifest entry count mismatch")
    bypath = Dict(e.path => e for e in current.entries)
    for e in recorded.entries
        haskey(bypath, e.path) || _err("manifest missing path $(e.path)")
        _same_entry(e, bypath[e.path]; check_sha=e.path != "manifest.toml") ||
            _err("manifest metadata mismatch at $(e.path)")
    end
    expected === nothing || begin
        recorded.digest == expected.digest || _err("manifest digest mismatch")
        length(expected.entries) == length(recorded.entries) || _err("expected manifest entry count mismatch")
        for e in expected.entries
            found = findfirst(x -> x.path == e.path, recorded.entries)
            found === nothing && _err("expected manifest path missing $(e.path)")
            _same_entry(e, recorded.entries[found]; check_sha=true) ||
                _err("expected manifest changed at $(e.path)")
        end
    end
    recorded
end

function _provenance_rows(root::String)
    source = joinpath(root, "source")
    isdir(source) || _err("source category is missing")
    rows = String[]
    for (rel, kind_string, st) in _walk(source, ".")
        kind_string == "file" || continue
        rel == "provenance.toml" && continue
        p = joinpath(source, rel == "." ? "" : rel)
        push!(rows, "path=$(repr(_rel(rel))) sha256=$(repr(sha256_file(p))) bytes=$(st.size)")
    end
    sort!(rows)
end

function _write_source_provenance!(root::String)
    rows = _provenance_rows(root)
    isempty(rows) && _err("source category has no regular files")
    text = "format = \"stmfit-source-provenance-v2\"\n" * join(rows, "\n") * "\n"
    _write_bytes(joinpath(root, "source-provenance.toml"), Vector{UInt8}(codeunits(text)); mode=0o444)
end

function verify_source_provenance(root::AbstractString)
    root = _normroot(root)
    p = joinpath(root, "source-provenance.toml")
    isfile(p) || _err("source provenance is missing")
    expected = Set{String}()
    for line in split(read(p, String), '\n')
        startswith(line, "path=") || continue
        m = match(r"^path=\"(.*)\" sha256=\"([0-9a-f]+)\" bytes=([0-9]+)$", line)
        m === nothing && _err("malformed source provenance row")
        rel, digest, bytes = m.captures
        full = joinpath(root, "source", rel)
        isfile(full) || _err("source provenance file missing: $(rel)")
        sha256_file(full) == digest || _err("source provenance hash mismatch: $(rel)")
        stat(full).size == parse(Int, bytes) || _err("source provenance size mismatch: $(rel)")
        push!(expected, rel)
    end
    actual = Set{String}()
    for (rel, kind, _) in _walk(joinpath(root, "source"), ".")
        kind == "file" && rel != "provenance.toml" && push!(actual, rel)
    end
    expected == actual || _err("source provenance inventory mismatch")
    true
end

function _fsync_path(path::AbstractString; directory=false)
    flags = directory ? Int(Base.Filesystem.JL_O_RDONLY | Base.Filesystem.JL_O_DIRECTORY) :
                        Int(Base.Filesystem.JL_O_RDONLY)
    fd = ccall(:open, Cint, (Cstring, Cint), String(path), flags)
    fd < 0 && _err("open for fsync failed: $(path)")
    try
        rc = ccall(:fsync, Cint, (Cint,), fd)
        rc == 0 || _err("fsync failed: $(path)")
    finally
        ccall(:close, Cint, (Cint,), fd)
    end
    nothing
end

function _fsync_tree(root::String)
    for (rel, kind, _) in reverse(_walk(root))
        p = rel == "." ? root : joinpath(root, rel)
        _fsync_path(p; directory=kind == "dir")
    end
end

function _rename_noreplace(src::String, dst::String)
    # Linux's renameat2 is the only operation here that provides the required
    # no-replace guarantee for directories.  There is intentionally no
    # check-then-rename fallback.
    parent = dirname(dst)
    isdir(parent) || _err("destination parent is absent: $(parent)")
    flags = 1 # RENAME_NOREPLACE
    rc = ccall(:renameat2, Cint,
               (Cint, Cstring, Cint, Cstring, Cuint), -100, src, -100, dst, flags)
    rc == 0 || begin
        errno = Base.Libc.errno()
        _err("no-replace rename failed ($(errno)): $(src) -> $(dst)")
    end
    dst
end

function _mkdir_stage(parent::String, base::String)
    for _ in 1:100
        candidate = joinpath(parent, "." * base * ".stage-" * string(rand(UInt128), base=16))
        try
            mkdir(candidate, mode=0o700)
            return candidate
        catch e
            ispath(candidate) || rethrow(e)
        end
    end
    _err("could not allocate a hidden sibling stage")
end

function _remove_stage(path::String)
    ispath(path) && rm(path; recursive=true, force=true)
    nothing
end

function _o_excl_write(path::String, bytes::Vector{UInt8}; mode::Integer=0o444)
    flags = Int(Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT |
                Base.Filesystem.JL_O_EXCL | Base.Filesystem.JL_O_CLOEXEC)
    fd = ccall(:open, Cint, (Cstring, Cint, Cuint), path, flags, UInt32(mode))
    fd < 0 && _err("exclusive receipt create failed: $(path) (errno $(Base.Libc.errno()))")
    io = Base.fdio(fd, true)
    try
        write(io, bytes)
        flush(io)
    finally
        close(io)
    end
    path
end

function write_external_receipt!(root::AbstractString, data::AbstractDict)
    root = _normroot(root)
    rp = receipt_path(root)
    ispath(rp) && _err("receipt already exists: $(rp)")
    io = IOBuffer()
    TOML.print(io, Dict{String,Any}(String(k) => v for (k, v) in pairs(data)); sorted=true)
    _o_excl_write(rp, Vector{UInt8}(take!(io)))
end

function publish_staged!(stage::AbstractString, final_root::AbstractString,
                         receipt::AbstractDict; fault=nothing)
    stage = _normroot(stage)
    final_root = _normroot(final_root)
    isdir(stage) || _err("stage is not a directory")
    parent = dirname(final_root)
    isdir(parent) || _err("fixed lexical parent is absent: $(parent)")
    ispath(final_root) && _err("final root already exists")
    rp = receipt_path(final_root)
    ispath(rp) && _err("final receipt already exists")
    _fsync_tree(stage)
    fault === nothing || fault(:before_rename, stage, final_root)
    fault === nothing || fault(:during_rename, stage, final_root)
    renamed = false
    try
        _rename_noreplace(stage, final_root)
        renamed = true
        _fsync_path(parent; directory=true)
        fault === nothing || fault(:after_rename_before_receipt, final_root, rp)
        write_external_receipt!(final_root, receipt)
    catch e
        if renamed
            throw(AmbiguousPublicationError(final_root, rp, e))
        end
        rethrow()
    end
    final_root
end

function sanitize_environment(values)
    pairs_ = values isa AbstractDict ? collect(pairs(values)) : values
    out = Dict{String,String}()
    for item in pairs_
        k, v = item isa Pair ? (String(first(item)), String(last(item))) : begin
            s = String(item)
            occursin('=', s) || _err("environment binding must be KEY=VALUE")
            split(s, "="; limit=2)
        end
        occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", k) || _err("invalid environment key: $(k)")
        occursin(r"(?i)(token|secret|password|credential|cookie|private[_-]?key|auth)", k) &&
            _err("secret-bearing environment key is not allowed: $(k)")
        haskey(out, k) && _err("duplicate environment key: $(k)")
        out[k] = v
    end
    isempty(out) && _err("sanitized environment must be nonempty")
    Dict{String,String}(k => out[k] for k in sort(collect(keys(out))))
end

function canonical_command(command)
    if command isa AbstractString
        isempty(strip(command)) && _err("bootstrap command is empty")
        return String(command)
    end
    vals = String.(collect(command))
    isempty(vals) && _err("bootstrap command is empty")
    join(vals, " ")
end

function _require_binding(name::String, path)
    path === nothing && _err("missing authority binding: $(name)")
    s = _normroot(String(path))
    ispath(s) || _err("authority binding $(name) is absent: $(s)")
    st = _lstat(s)
    _kind(st)
    s
end

function _binding_record(name::String, path::String)
    st = _lstat(path)
    Dict{String,Any}("name" => name, "path" => path,
                     "kind" => String(_kind(st)),
                     "sha256" => tree_digest(path),
                     "bytes" => Int(st.size), "mode" => Int(st.mode),
                     "nlink" => Int(st.nlink), "device" => Int(st.device),
                     "inode" => Int(st.inode))
end

function _bindings_file!(stage::String, b::Dict{String,String}, bootstrap,
                         bootstrap_version, bootstrap_sha, env, command)
    data = Dict{String,Any}(
        "format" => "stmfit-authority-bindings-v2",
        "bindings" => [_binding_record(k, b[k]) for k in sort(collect(keys(b)))],
        "bootstrap" => Dict{String,Any}(
            "executable" => bootstrap,
            "version" => bootstrap_version,
            "sha256" => bootstrap_sha,
        ),
        "environment" => env,
        "command" => command,
    )
    io = IOBuffer()
    TOML.print(io, data; sorted=true)
    _write_bytes(joinpath(stage, "authority-bindings.toml"), Vector{UInt8}(take!(io)); mode=0o444)
end

function _ensure_categories(stage::String)
    for category in COMMON_CATEGORIES
        p = joinpath(stage, category)
        isdir(p) || _err("missing authority category: $(category)")
        isempty(readdir(p)) && _err("authority category is empty: $(category)")
    end
end

function _reject_runs(stage::String)
    ispath(joinpath(stage, "runs")) && _err("common authority must not contain runs")
end

function _verify_tree_objects(root::String)
    for (_, _, st) in _walk(root)
        _kind(st) # deliberately reject links and special files
    end
end

function _chmod_for_authority!(stage::String)
    chmod(stage, 0o555)
    _fsync_tree(stage)
end

"The production Lanes3 seam.  It is intentionally not a public CLI option."
function materialize_common_payload!(; kwargs...)
    lanes = if isdefined(Main, :Lanes3)
        getfield(Main, :Lanes3)
    else
        try
            Base.require(Main, :Lanes3)
            getfield(Main, :Lanes3)
        catch
            _err("Lanes3 is required by the production common builder")
        end
    end
    isdefined(lanes, :materialize_common_payload!) ||
        _err("Lanes3.materialize_common_payload! is not exported by Lanes3")
    getfield(lanes, :materialize_common_payload!)(; kwargs...)
end

function _run_materializer(materializer, stage::String, bindings; production::Bool)
    kwargs = (; root=stage,
              transitive=bindings["transitive"],
              project=bindings["project"],
              manifest=bindings["manifest"],
              fixture=bindings["fixture"],
              source=bindings["source"],
              config=bindings["config"],
              package=bindings["package"],
              test=bindings["test"],
              audit=bindings["audit"],
              harness=bindings["harness"],
              closure=bindings["closure"],
              materialize_synthetic_inputs=true)
    if production
        # Keep this call visibly and mechanically separate: production cannot
        # choose an injected function through a CLI flag.
        materialize_common_payload!(; kwargs...)
    else
        materializer === nothing && _err("test materializer was not supplied")
        materializer(; kwargs...)
    end
end

function _validate_bootstrap(executable, version, digest)
    exe = _require_binding("bootstrap_executable", executable)
    _kind(_lstat(exe)) === :file || _err("bootstrap executable is not a file")
    isempty(strip(String(version))) && _err("bootstrap version is empty")
    actual = sha256_file(exe)
    lowercase(String(digest)) == actual ||
        _err("bootstrap executable SHA-256 mismatch: expected $(digest), got $(actual)")
    exe, String(version), actual
end

function _build_common(; output_root, transitive, project, manifest, fixture, source,
                        config, package, test, audit, harness, closure,
                        bootstrap_executable, bootstrap_version,
                        bootstrap_sha256, env, command, materializer=nothing,
                        production=false, dry_run=false, fault=nothing)
    output_root = _normroot(output_root)
    parent = dirname(output_root)
    isdir(parent) || _err("fixed lexical parent is absent: $(parent)")
    ispath(output_root) && _err("final common root already exists: $(output_root)")
    ispath(receipt_path(output_root)) && _err("final common receipt already exists")
    env = sanitize_environment(env)
    command = canonical_command(command)
    exe, version, sha = _validate_bootstrap(bootstrap_executable, bootstrap_version, bootstrap_sha256)
    names = Dict{String,Any}("transitive" => transitive, "project" => project,
                             "manifest" => manifest, "fixture" => fixture,
                             "source" => source, "config" => config,
                             "package" => package, "test" => test,
                             "audit" => audit, "harness" => harness,
                             "closure" => closure)
    bindings = Dict{String,String}()
    for (name, value) in names
        # `nothing` is a validated default only for optional categories.  The
        # production CLI supplies all eleven names explicitly.
        value === nothing && (name in ("harness", "closure") ? continue :
                              _err("missing authority binding: $(name)"))
        bindings[name] = _require_binding(name, value)
    end
    dry_run && return (dry_run=true, success=false, root=output_root,
                       receipt=receipt_path(output_root), schema="A1a-storage-v2")

    stage = _mkdir_stage(parent, basename(output_root))
    renamed = false
    try
        _run_materializer(materializer, stage, bindings; production=production)
        _verify_tree_objects(stage)
        _reject_runs(stage)
        _ensure_categories(stage)
        _write_bindings_file!(stage, bindings, exe, version, sha, env, command)
        _write_source_provenance!(stage)
        _verify_tree_objects(stage)
        manifest_ = write_manifest!(stage)
        verify_source_provenance(stage)
        verify_manifest(stage, manifest_)
        _chmod_for_authority!(stage)
        # Recheck after chmod and fsync: the stage identity and every recorded
        # field must be the one that is about to cross the rename boundary.
        verify_manifest(stage, manifest_)
        fault === nothing || fault(:before_receipt, stage, output_root)
        receipt = Dict{String,Any}(
            "format" => "stmfit-common-receipt-v2",
            "root" => output_root,
            "manifest_sha256" => manifest_.digest,
            "entry_count" => length(manifest_.entries),
            "labels" => collect(IMAGE_LABELS),
            "authority_sha256" => sha,
            "command" => command,
            "environment" => env,
        )
        result = publish_staged!(stage, output_root, receipt; fault=fault)
        renamed = true
        (root=result, receipt=receipt_path(result), manifest=manifest_, receipt_data=receipt)
    catch e
        # After a successful rename the root is intentionally preserved.  A
        # missing/partial receipt is an ambiguous state, never a repair case.
        if !renamed && ispath(stage)
            _remove_stage(stage)
        end
        rethrow()
    end
end

function build_common_authority(; output_root, transitive, project, manifest, fixture,
                                source, config, package, test, audit,
                                harness=nothing, closure=nothing,
                                bootstrap_executable, bootstrap_version,
                                bootstrap_sha256, env, command,
                                materializer=nothing, dry_run=false, fault=nothing)
    _build_common(; output_root, transitive, project, manifest, fixture, source,
                  config, package, test, audit, harness, closure,
                  bootstrap_executable, bootstrap_version, bootstrap_sha256,
                  env, command, materializer, production=false, dry_run, fault)
end

function _copy_tree(src::String, dst::String; root_mode=nothing)
    st = _lstat(src)
    kind = _kind(st)
    if kind === :file
        mkpath(dirname(dst))
        _write_bytes(dst, read(src); mode=Int(st.mode & 0o7777))
        return
    end
    mkdir(dst; mode=root_mode === nothing ? Int(st.mode & 0o7777) : root_mode)
    for name in sort(readdir(src))
        _copy_tree(joinpath(src, name), joinpath(dst, name))
    end
    chmod(dst, Int(st.mode & 0o7777))
end

function _copy_authority(common::String, image::String)
    for name in sort(readdir(common))
        _copy_tree(joinpath(common, name), joinpath(image, name))
    end
end

function _inode_set(root::String)
    Set{Tuple{UInt64,UInt64}}((UInt64(st.device), UInt64(st.inode)) for (_, _, st) in _walk(root))
end

function _assert_disjoint(roots::Vector{String})
    seen = Set{Tuple{UInt64,UInt64}}()
    for root in roots
        current = _inode_set(root)
        !isempty(intersect(seen, current)) && _err("shared inode between authority copies")
        union!(seen, current)
    end
    true
end

function _image_receipt!(image::String, label::String, common_receipt_sha::String,
                         prestate_slot::String, fixed_path::String, common_manifest_sha::String)
    data = Dict{String,Any}(
        "format" => "stmfit-image-receipt-v2",
        "label" => label,
        "common_receipt_sha256" => common_receipt_sha,
        "common_manifest_sha256" => common_manifest_sha,
        "prestate_slot" => prestate_slot,
        "fixed_lexical_path" => fixed_path,
    )
    io = IOBuffer()
    TOML.print(io, data; sorted=true)
    _write_bytes(joinpath(image, "image-receipt.toml"), Vector{UInt8}(take!(io)); mode=0o444)
end

function _read_receipt(path::String)
    isfile(path) || _err("receipt is missing: $(path)")
    TOML.parsefile(path)
end

function _image_manifest!(image::String)
    m = write_manifest!(image)
    verify_manifest(image, m)
    m
end

function execution_image(root::AbstractString; label::AbstractString="",
                         common_receipt_sha256::AbstractString="",
                         prestate_slot::AbstractString="",
                         fixed_lexical_path::AbstractString=String(root),
                         manifest::Union{Nothing,Manifest}=nothing)
    root = _normroot(root)
    m = manifest === nothing ? verify_manifest(root) : manifest
    ExecutionImage(root, root, String(label), String(common_receipt_sha256),
                   String(prestate_slot), _rel(fixed_lexical_path), m)
end

function relocate_no_replace(image::ExecutionImage, destination::AbstractString)
    src = _normroot(image.root)
    dst = _normroot(destination)
    isdir(src) || _err("execution image source is absent")
    ispath(dst) && _err("fixed-path collision: $(dst)")
    isdir(dirname(dst)) || _err("relocation parent is absent")
    _lstat(src).device == _lstat(dirname(dst)).device || _err("cross-device relocation")
    _rename_noreplace(src, dst)
    image.root = dst
    dst
end

function verify_relocated(image::ExecutionImage)
    isdir(image.root) || _err("relocated image is absent")
    ispath(image.original_root) && image.original_root != image.root &&
        _err("relocation residue remains at original path")
    verify_manifest(image.root, image.manifest)
    true
end

function restore_relocated!(image::ExecutionImage)
    image.root == image.original_root && _err("image is not relocated")
    ispath(image.original_root) && _err("restore destination collision")
    isdir(dirname(image.original_root)) || _err("restore parent is absent")
    _lstat(image.root).device == _lstat(dirname(image.original_root)).device ||
        _err("cross-device restore")
    _rename_noreplace(image.root, image.original_root)
    image.root = image.original_root
    verify_manifest(image.root, image.manifest)
    image.root
end

function _verify_image_root(image::String, label::String, common_sha::String,
                            prestate::String, fixed_path::String, common_manifest_sha::String)
    isdir(image) || _err("missing image root: $(label)")
    st = _lstat(image)
    (st.mode & 0o7777) == 0o555 || _err("wrong image root mode: $(label)")
    r = _read_receipt(joinpath(image, "image-receipt.toml"))
    r["label"] == label || _err("wrong image label receipt")
    r["common_receipt_sha256"] == common_sha || _err("wrong common receipt binding")
    r["common_manifest_sha256"] == common_manifest_sha || _err("wrong common manifest binding")
    r["prestate_slot"] == prestate || _err("wrong prestate slot")
    r["fixed_lexical_path"] == _rel(fixed_path) || _err("wrong fixed lexical path")
    runs = joinpath(image, "runs")
    isdir(runs) || _err("image runs directory is missing")
    ( _lstat(runs).mode & 0o7777) == 0o700 || _err("wrong runs directory mode")
    isempty(readdir(runs)) || _err("image runs directory is not empty")
    verify_manifest(image)
    true
end

function verify_image_set(set_root::AbstractString; common_receipt_sha256::AbstractString,
                          prestate_slot::AbstractString, future_run_root::AbstractString,
                          common_manifest_sha256::AbstractString)
    set_root = _normroot(set_root)
    isdir(set_root) || _err("image set is absent")
    labels = sort(readdir(set_root))
    labels == sort(collect(IMAGE_LABELS)) || _err("extra or missing image-set member")
    roots = String[]
    for label in IMAGE_LABELS
        image = joinpath(set_root, label)
        fixed = joinpath(String(future_run_root), label)
        _verify_image_root(image, label, String(common_receipt_sha256),
                           String(prestate_slot), fixed, String(common_manifest_sha256))
        push!(roots, image)
    end
    _assert_disjoint(roots)
    true
end

function _build_images(; common_root, output_parent, set_name, prestate_slot,
                        future_run_root, dry_run=false, fault=nothing)
    common_root = _normroot(common_root)
    output_parent = _normroot(output_parent)
    set_root = joinpath(output_parent, String(set_name))
    isdir(output_parent) || _err("image-set parent is absent")
    ( _lstat(output_parent).mode & 0o7777) == 0o700 ||
        _err("image-set parent must remain mode 0700")
    ispath(set_root) && _err("final image set already exists")
    ispath(receipt_path(set_root)) && _err("final image-set receipt already exists")
    isdir(common_root) || _err("common authority root is absent")
    common_receipt = receipt_path(common_root)
    common_data = _read_receipt(common_receipt)
    common_manifest = verify_manifest(common_root)
    verify_source_provenance(common_root)
    String(common_data["manifest_sha256"]) == common_manifest.digest ||
        _err("common receipt does not bind the common manifest")
    future = _rel(abspath(String(future_run_root)))
    ispath(String(future_run_root)) && _err("future run root must be absent")
    isempty(strip(String(prestate_slot))) && _err("prestate slot is empty")
    dry_run && return (dry_run=true, success=false, root=set_root,
                       receipt=receipt_path(set_root), schema="A1a-storage-v2-images")
    # The lexical parent is fixed and is never created or chmod'd here.
    _lstat(output_parent).device == _lstat(dirname(set_root)).device ||
        _err("image stage would cross device boundary")
    stage = _mkdir_stage(output_parent, String(set_name))
    renamed = false
    try
        image_roots = String[]
        for label in IMAGE_LABELS
            image = joinpath(stage, label)
            mkdir(image; mode=0o700)
            _copy_authority(common_root, image)
            mkdir(joinpath(image, "runs"); mode=0o700)
            _write_bytes(joinpath(image, "common-receipt.toml"), read(common_receipt); mode=0o444)
            _image_receipt!(image, label, sha256_file(common_receipt), String(prestate_slot),
                            joinpath(String(future_run_root), label), common_manifest.digest)
            im = _image_manifest!(image)
            chmod(image, 0o555)
            _fsync_tree(image)
            verify_manifest(image, im)
            push!(image_roots, image)
        end
        _assert_disjoint(image_roots)
        for label in IMAGE_LABELS
            image = joinpath(stage, label)
            _lstat(image).device == _lstat(output_parent).device || _err("image device mismatch")
        end
        chmod(stage, 0o555)
        _fsync_tree(stage)
        verify_image_set(stage; common_receipt_sha256=sha256_file(common_receipt),
                         prestate_slot=String(prestate_slot), future_run_root=future,
                         common_manifest_sha256=common_manifest.digest)
        fault === nothing || fault(:before_receipt, stage, set_root)
        members = [Dict{String,Any}("label" => label,
                                    "manifest_sha256" => read_manifest(joinpath(stage, label, "manifest.toml")).digest,
                                    "path" => joinpath(String(future_run_root), label))
                   for label in IMAGE_LABELS]
        receipt = Dict{String,Any}(
            "format" => "stmfit-image-set-receipt-v2",
            "root" => set_root,
            "common_receipt_sha256" => sha256_file(common_receipt),
            "common_manifest_sha256" => common_manifest.digest,
            "prestate_slot" => String(prestate_slot),
            "future_run_root" => future,
            "labels" => collect(IMAGE_LABELS),
            "members" => members,
        )
        result = publish_staged!(stage, set_root, receipt; fault=fault)
        renamed = true
        (root=result, receipt=receipt_path(result), receipt_data=receipt)
    catch e
        !renamed && ispath(stage) && _remove_stage(stage)
        rethrow()
    end
end

function build_execution_images(; common_root, output_parent, set_name,
                                prestate_slot, future_run_root,
                                dry_run=false, fault=nothing)
    _build_images(; common_root, output_parent, set_name, prestate_slot,
                  future_run_root, dry_run, fault)
end

end # module AuthorityStorage
