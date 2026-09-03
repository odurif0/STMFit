module AuthorityIdentity

"""
The small, deliberately boring, runtime identity boundary used by the
integration harness.  This file has no knowledge of the comparator or of any
scientific result.  In particular, a successful identity capture is not a
Gate-2 or a PASS claim.
"""

using Dates
using Libdl
using LinearAlgebra
using Printf
using SHA
using TOML
using UUIDs

export IdentityError, DSORecord, DSOContract,
       parse_dso_contract, write_dso_contract, sha256_file,
       capture_identity, validate_identity, validate_runtime_identity,
       publish_identity, identity_summary

const IDENTITY_SCHEMA = "stmfit.authority-identity-v2"
const DSO_SCHEMA = "runtime-dso-contract-v1"
const MANIFEST_SCHEMA = "stmfit.authority-identity-manifest-v2"
const RECEIPT_SCHEMA = "stmfit.authority-identity-receipt-v2"
const _BAD_WORDS = Set(["", "UNSET", "UNAVAILABLE", "UNKNOWN", "NONE"])

struct IdentityError <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, e::IdentityError) = print(io, "[", e.code, "] ", e.message)

_fail(code::Symbol, message) = throw(IdentityError(code, String(message)))

function _text(x, field::AbstractString)
    x === nothing && _fail(:unavailable, "$field is unavailable")
    s = strip(String(x))
    upper = uppercase(String(s))
    (upper in _BAD_WORDS || startswith(upper, "UNSET") || startswith(upper, "UNAVAILABLE")) &&
        _fail(:unavailable, "$field is unset or unavailable")
    return String(s)
end

function _hex(x, field::AbstractString)
    s = lowercase(_text(x, field))
    (length(s) == 64 && all(c -> c in "0123456789abcdef", s)) ||
        _fail(:invalid_hash, "$field is not a lowercase SHA-256 digest")
    return s
end

_sha256_bytes(x) = bytes2hex(SHA.sha256(x))
sha256_file(path::AbstractString) = _sha256_bytes(read(path))

function _canonical(path::AbstractString)
    p = isabspath(path) ? String(path) : abspath(String(path))
    ispath(p) || _fail(:missing_path, "path does not exist: $p")
    return realpath(p)
end

function _default_file_info(path::AbstractString)
    p = _canonical(path)
    isfile(p) || _fail(:not_regular_file, "expected a regular file: $p")
    st = stat(p)
    return Dict{String,Any}(
        "path" => p,
        "size" => Int64(st.size),
        "sha256" => sha256_file(p),
    )
end

function _file_mode(path::AbstractString)
    return Int(stat(path).mode & 0o7777)
end

function _fsync_fd(fd)
    # Julia 1.10 exposes an Int64 here; newer Julia exposes a four-byte
    # RawFD.  Normalize both without relying on an internal type name.
    number = fd isa Integer ? Int32(fd) : reinterpret(Int32, fd)
    rc = ccall(:fsync, Cint, (Cint,), number)
    rc == 0 || _fail(:fsync_failed, "fsync failed")
    return nothing
end

function _fsync_file(path::AbstractString)
    io = open(path, "r")
    try
        flush(io)
        _fsync_fd(Base.fd(io))
    finally
        close(io)
    end
    return nothing
end

function _fsync_directory(path::AbstractString)
    flags = Base.Filesystem.JL_O_RDONLY | Base.Filesystem.JL_O_DIRECTORY
    f = Base.Filesystem.open(path, flags)
    try
        _fsync_fd(Base.fd(f))
    finally
        close(f)
    end
    return nothing
end

function _write_toml(path::AbstractString, value::Dict{String,Any}; mode::Integer=0o644)
    open(path, "w") do io
        TOML.print(io, value)
        flush(io)
        _fsync_fd(Base.fd(io))
    end
    chmod(path, mode)
    _fsync_file(path)
    return path
end

function _mode_string(mode::Integer)
    return @sprintf("%04o", mode & 0o7777)
end

# --------------------------------------------------------------------------
# Gate-1 DSO contract

struct DSORecord
    name::String
    path::String
    size::Int64
    sha256::String
    runtime_label::Union{Nothing,String}
    runtime_version::Union{Nothing,String}
end

function DSORecord(; name, path, size, sha256, runtime_label=nothing,
                   runtime_version=nothing)
    n = _text(name, "DSO name")
    p = _canonical(String(path))
    z = try Int64(size) catch; _fail(:dso_schema, "DSO size is not an integer") end
    z >= 0 || _fail(:dso_schema, "DSO size is negative")
    h = _hex(sha256, "DSO sha256")
    rl = runtime_label === nothing ? nothing : _text(runtime_label, "DSO runtime label")
    rv = runtime_version === nothing ? nothing : _text(runtime_version, "DSO runtime version")
    basename(p) == n || _fail(:dso_schema, "DSO name does not match its basename: $n")
    return DSORecord(n, p, z, h, rl, rv)
end

struct DSOContract
    schema::String
    runtime_label::String
    runtime_version::String
    common::Vector{DSORecord}
    versioned::Vector{DSORecord}
end

function DSOContract(; schema=DSO_SCHEMA, runtime_label, runtime_version,
                     common=DSORecord[], versioned=DSORecord[])
    s = _text(schema, "DSO contract schema")
    s == DSO_SCHEMA || _fail(:dso_schema, "unsupported DSO contract schema: $s")
    rl = _text(runtime_label, "DSO contract runtime label")
    rv = _text(runtime_version, "DSO contract runtime version")
    c = collect(common)
    v = collect(versioned)
    all(r -> r isa DSORecord, c) || _fail(:dso_schema, "common providers are not DSO records")
    all(r -> r isa DSORecord, v) || _fail(:dso_schema, "versioned components are not DSO records")
    length(unique(getfield.(c, :path))) == length(c) ||
        _fail(:dso_schema, "duplicate common DSO path")
    length(unique(getfield.(v, :path))) == length(v) ||
        _fail(:dso_schema, "duplicate versioned DSO path")
    isempty(intersect(Set(getfield.(c, :path)), Set(getfield.(v, :path)))) ||
        _fail(:dso_schema, "a DSO occurs in both common and versioned sets")
    for r in c
        r.runtime_label === nothing && r.runtime_version === nothing ||
            _fail(:dso_schema, "common DSO has runtime-specific fields: $(r.name)")
    end
    for r in v
        r.runtime_label == rl && r.runtime_version == rv ||
            _fail(:dso_runtime, "versioned DSO is for a different runtime: $(r.name)")
    end
    return DSOContract(s, rl, rv, c, v)
end

function _exact_keys(d::AbstractDict, expected, where::AbstractString)
    got = sort!(String.(collect(keys(d))))
    want = sort!(String.(collect(expected)))
    got == want || _fail(:dso_schema, "$where has keys $(got), expected $(want)")
    return nothing
end

function _dict(x, where)
    x isa AbstractDict || _fail(:dso_schema, "$where must be an object/table")
    return x
end

function _records(x, where::AbstractString, versioned::Bool, rl::String, rv::String)
    x isa AbstractVector || _fail(:dso_schema, "$where must be an array")
    records = DSORecord[]
    for (i, item) in enumerate(x)
        d = _dict(item, "$where[$i]")
        keys = versioned ? ["name", "path", "size", "sha256", "runtime_label", "runtime_version"] :
                          ["name", "path", "size", "sha256"]
        _exact_keys(d, keys, "$where[$i]")
        for key in ("name", "path", "sha256")
            d[key] isa AbstractString || _fail(:dso_schema, "$where[$i].$key must be a string")
        end
        d["size"] isa Integer || _fail(:dso_schema, "$where[$i].size must be an integer")
        d["size"] >= 0 || _fail(:dso_schema, "$where[$i].size is negative")
        if versioned
            d["runtime_label"] isa AbstractString || _fail(:dso_schema, "$where[$i].runtime_label must be a string")
            d["runtime_version"] isa AbstractString || _fail(:dso_schema, "$where[$i].runtime_version must be a string")
        end
        r = versioned ? DSORecord(name=d["name"], path=d["path"], size=d["size"],
                                  sha256=d["sha256"], runtime_label=d["runtime_label"],
                                  runtime_version=d["runtime_version"]) :
                        DSORecord(name=d["name"], path=d["path"], size=d["size"],
                                  sha256=d["sha256"])
        versioned && (r.runtime_label == rl && r.runtime_version == rv ||
                      _fail(:dso_runtime, "$where[$i] targets another runtime"))
        push!(records, r)
    end
    return records
end

function _contract_from_dict(root::AbstractDict)
    _exact_keys(root, ["schema", "runtime_label", "runtime_version", "common", "versioned"],
                "DSO contract")
    root["schema"] isa AbstractString || _fail(:dso_schema, "DSO contract schema must be a string")
    root["runtime_label"] isa AbstractString || _fail(:dso_schema, "DSO contract runtime label must be a string")
    root["runtime_version"] isa AbstractString || _fail(:dso_schema, "DSO contract runtime version must be a string")
    schema = _text(root["schema"], "DSO contract schema")
    rl = _text(root["runtime_label"], "DSO contract runtime label")
    rv = _text(root["runtime_version"], "DSO contract runtime version")
    common = _dict(root["common"], "common")
    versioned = _dict(root["versioned"], "versioned")
    _exact_keys(common, ["count", "providers"], "common")
    _exact_keys(versioned, ["count", "components"], "versioned")
    common["count"] isa Integer || _fail(:dso_schema, "common.count must be an integer")
    versioned["count"] isa Integer || _fail(:dso_schema, "versioned.count must be an integer")
    common["count"] >= 0 || _fail(:dso_schema, "common.count is negative")
    versioned["count"] >= 0 || _fail(:dso_schema, "versioned.count is negative")
    cp = _records(common["providers"], "common.providers", false, rl, rv)
    vp = _records(versioned["components"], "versioned.components", true, rl, rv)
    Int(common["count"]) == length(cp) || _fail(:dso_count, "common.count does not match providers")
    Int(versioned["count"]) == length(vp) || _fail(:dso_count, "versioned.count does not match components")
    return DSOContract(schema=schema, runtime_label=rl, runtime_version=rv,
                       common=cp, versioned=vp)
end

# A minimal JSON reader is kept here so the Gate-1 contract does not acquire a
# package dependency merely to read a six-key, content-addressed contract.
mutable struct _JSONParser
    text::String
    index::Int
end

function _jskip!(p::_JSONParser)
    while p.index <= lastindex(p.text) && p.text[p.index] in (' ', '\n', '\r', '\t')
        p.index = nextind(p.text, p.index)
    end
end

function _jchar!(p::_JSONParser)
    p.index > lastindex(p.text) && _fail(:dso_schema, "truncated JSON string")
    c = p.text[p.index]
    p.index = nextind(p.text, p.index)
    return c
end

function _jstring!(p::_JSONParser)
    _jchar!(p) == '"' || _fail(:dso_schema, "expected JSON string")
    io = IOBuffer()
    while true
        p.index > lastindex(p.text) && _fail(:dso_schema, "unterminated JSON string")
        c = _jchar!(p)
        c == '"' && return String(take!(io))
        c == '\\' || (write(io, c); continue)
        e = _jchar!(p)
        if e == '"' || e == '\\' || e == '/'
            write(io, e)
        elseif e == 'b'
            write(io, '\b')
        elseif e == 'f'
            write(io, '\f')
        elseif e == 'n'
            write(io, '\n')
        elseif e == 'r'
            write(io, '\r')
        elseif e == 't'
            write(io, '\t')
        elseif e == 'u'
            p.index + 3 <= lastindex(p.text) || _fail(:dso_schema, "short JSON unicode escape")
            h = p.text[p.index:nextind(p.text, nextind(p.text, nextind(p.text, p.index)))]
            all(c -> c in "0123456789abcdefABCDEF", h) || _fail(:dso_schema, "bad JSON unicode escape")
            write(io, Char(parse(UInt32, h; base=16)))
            for _ in 1:4
                p.index = nextind(p.text, p.index)
            end
        else
            _fail(:dso_schema, "bad JSON escape")
        end
    end
end

function _jvalue!(p::_JSONParser)
    _jskip!(p)
    p.index > lastindex(p.text) && _fail(:dso_schema, "truncated JSON")
    c = p.text[p.index]
    c == '{' && return _jobject!(p)
    c == '[' && return _jarray!(p)
    c == '"' && return _jstring!(p)
    rest = p.text[p.index:end]
    startswith(rest, "true") && (p.index += 4; return true)
    startswith(rest, "false") && (p.index += 5; return false)
    startswith(rest, "null") && (p.index += 4; return nothing)
    m = match(r"^-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?", rest)
    m === nothing && _fail(:dso_schema, "invalid JSON value")
    token = m.match
    p.index += ncodeunits(token)
    return occursin(r"[.eE]", token) ? parse(Float64, token) : parse(Int, token)
end

function _jobject!(p::_JSONParser)
    _jchar!(p) == '{' || _fail(:dso_schema, "expected JSON object")
    d = Dict{String,Any}()
    _jskip!(p)
    if p.index <= lastindex(p.text) && p.text[p.index] == '}'
        p.index = nextind(p.text, p.index)
        return d
    end
    while true
        _jskip!(p)
        k = _jstring!(p)
        haskey(d, k) && _fail(:dso_schema, "duplicate JSON key: $k")
        _jskip!(p)
        _jchar!(p) == ':' || _fail(:dso_schema, "expected ':' after JSON key")
        d[k] = _jvalue!(p)
        _jskip!(p)
        c = _jchar!(p)
        c == '}' && return d
        c == ',' || _fail(:dso_schema, "expected ',' in JSON object")
    end
end

function _jarray!(p::_JSONParser)
    _jchar!(p) == '[' || _fail(:dso_schema, "expected JSON array")
    a = Any[]
    _jskip!(p)
    if p.index <= lastindex(p.text) && p.text[p.index] == ']'
        p.index = nextind(p.text, p.index)
        return a
    end
    while true
        push!(a, _jvalue!(p))
        _jskip!(p)
        c = _jchar!(p)
        c == ']' && return a
        c == ',' || _fail(:dso_schema, "expected ',' in JSON array")
    end
end

function _parse_json(text::String)
    p = _JSONParser(text, firstindex(text))
    value = _jvalue!(p)
    _jskip!(p)
    p.index > lastindex(p.text) || _fail(:dso_schema, "trailing JSON content")
    return _dict(value, "DSO contract")
end

function parse_dso_contract(path::AbstractString)
    p = _canonical(path)
    root = endswith(lowercase(p), ".json") ? _parse_json(String(read(p))) : TOML.parsefile(p)
    return _contract_from_dict(root)
end

function _dso_dict(r::DSORecord; versioned=(r.runtime_label !== nothing))
    d = Dict{String,Any}("name" => r.name, "path" => r.path,
                         "size" => r.size, "sha256" => r.sha256)
    versioned && begin
        d["runtime_label"] = r.runtime_label
        d["runtime_version"] = r.runtime_version
    end
    return d
end

function _contract_dict(c::DSOContract)
    return Dict{String,Any}(
        "schema" => c.schema,
        "runtime_label" => c.runtime_label,
        "runtime_version" => c.runtime_version,
        "common" => Dict{String,Any}("count" => length(c.common),
                                     "providers" => [_dso_dict(r) for r in c.common]),
        "versioned" => Dict{String,Any}("count" => length(c.versioned),
                                         "components" => [_dso_dict(r; versioned=true) for r in c.versioned]),
    )
end

function _json_quote(s::AbstractString)
    io = IOBuffer(); write(io, '"')
    for c in s
        if c == '"'; write(io, "\\\"")
        elseif c == '\\'; write(io, "\\\\")
        elseif c == '\n'; write(io, "\\n")
        elseif c == '\r'; write(io, "\\r")
        elseif c == '\t'; write(io, "\\t")
        else; write(io, c)
        end
    end
    write(io, '"'); return String(take!(io))
end

function _json_write(io::IO, x)
    if x isa AbstractDict
        print(io, '{'); first = true
        for k in sort!(String.(collect(keys(x))))
            first || print(io, ','); first = false
            print(io, _json_quote(k), ':'); _json_write(io, x[k])
        end
        print(io, '}')
    elseif x isa AbstractVector
        print(io, '[')
        for (i, v) in enumerate(x); i > 1 && print(io, ','); _json_write(io, v); end
        print(io, ']')
    elseif x isa AbstractString
        print(io, _json_quote(x))
    elseif x === nothing
        print(io, "null")
    elseif x isa Bool || x isa Integer || x isa AbstractFloat
        print(io, x)
    else
        _fail(:dso_schema, "cannot encode JSON value of type $(typeof(x))")
    end
end

function write_dso_contract(path::AbstractString, contract::DSOContract)
    p = String(path)
    data = _contract_dict(contract)
    if endswith(lowercase(p), ".json")
        open(p, "w") do io
            _json_write(io, data); write(io, '\n'); flush(io); _fsync_fd(Base.fd(io))
        end
    else
        _write_toml(p, data)
    end
    return p
end

function _validate_dso!(contract::DSOContract, actual_paths, actual_info,
                        runtime_label::String, runtime_version::String)
    contract.runtime_label == runtime_label || _fail(:dso_runtime, "DSO contract runtime label differs")
    contract.runtime_version == runtime_version || _fail(:dso_runtime, "DSO contract runtime version differs")
    expected = vcat(contract.common, contract.versioned)
    ep = Set(r.path for r in expected)
    ap = Set(actual_paths)
    isempty(setdiff(ep, ap)) || _fail(:dso_missing, "expected loaded DSO is missing")
    isempty(setdiff(ap, ep)) || _fail(:dso_extra, "an extra relevant loaded DSO is present")
    for r in expected
        a = actual_info[r.path]
        a["path"] == r.path || _fail(:dso_path, "loaded DSO canonical path differs: $(r.name)")
        Int64(a["size"]) == r.size || _fail(:dso_hash, "loaded DSO size differs: $(r.name)")
        lowercase(String(a["sha256"])) == r.sha256 || _fail(:dso_hash, "loaded DSO hash differs: $(r.name)")
    end
    return Dict{String,Any}(
        "schema" => contract.schema,
        "runtime_label" => runtime_label,
        "runtime_version" => runtime_version,
        "common_count" => length(contract.common),
        "versioned_count" => length(contract.versioned),
        "actual_count" => length(actual_paths),
        "common_providers" => [_dso_dict(r) for r in contract.common],
        "versioned_julia_components" => [_dso_dict(r; versioned=true) for r in contract.versioned],
    )
end

# --------------------------------------------------------------------------
# Private probe seam.  It is intentionally not exported and the CLI never
# accepts one.  Tests use it to replace OS observations without replacing the
# production capture path.

_as_probe_function(x) = x isa Function ? x : (() -> x)
_default_command() = Base.julia_cmd().exec[1]
_default_running_executable() = Sys.isunix() && ispath("/proc/self/exe") ? "/proc/self/exe" : _default_command()
_default_sysimage() = String(Base.JLOptions().image_file)
_default_project() = Base.active_project()
_default_manifest() = joinpath(dirname(String(Base.active_project())), "Manifest.toml")
_default_environment() = Dict{String,String}(String(k) => String(v) for (k, v) in ENV)
_default_kernel() = readchomp(`uname -srvm`)
_default_cpu() = read("/proc/cpuinfo", String)

function _default_microcode()
    for p in ("/sys/devices/system/cpu/microcode/version", "/sys/devices/system/cpu/cpu0/microcode/version")
        isfile(p) && return readchomp(p)
    end
    isfile("/proc/cpuinfo") || return nothing
    for line in split(read("/proc/cpuinfo", String), '\n')
        startswith(lowercase(strip(line)), "microcode") && return strip(split(line, ':'; limit=2)[end])
    end
    return nothing
end

function _default_lscpu()
    path = Sys.which("lscpu")
    path === nothing && return nothing
    return read(`$path --json`)
end

function _default_slurm()
    e = _default_environment()
    return Dict{String,Any}(
        "job_id" => get(e, "SLURM_JOB_ID", nothing),
        "node" => get(e, "SLURMD_NODENAME", nothing),
        "task_id" => get(e, "SLURM_PROCID", nothing),
        "array_job_id" => get(e, "SLURM_ARRAY_JOB_ID", nothing),
        "array_task_id" => get(e, "SLURM_ARRAY_TASK_ID", nothing),
    )
end

struct _SystemProbe
    version::Function
    julia_command::Function
    sys_bindir::Function
    running_executable::Function
    canonical_path::Function
    file_info::Function
    force_numerics::Function
    dllist::Function
    julia_threads::Function
    blas_threads::Function
    sysimage::Function
    project::Function
    manifest::Function
    environment::Function
    hostname::Function
    kernel::Function
    cpu::Function
    microcode::Function
    lscpu::Function
    slurm::Function
end

function _SystemProbe(; version=()->VERSION, julia_command=_default_command,
    sys_bindir=()->Sys.BINDIR, running_executable=_default_running_executable,
    canonical_path=_canonical, file_info=_default_file_info,
    force_numerics=()->nothing, dllist=Libdl.dllist,
    julia_threads=()->Threads.nthreads(), blas_threads=()->BLAS.get_num_threads(),
    sysimage=_default_sysimage, project=_default_project, manifest=_default_manifest,
    environment=_default_environment, hostname=()->gethostname(), kernel=_default_kernel,
    cpu=_default_cpu, microcode=_default_microcode, lscpu=_default_lscpu,
    slurm=_default_slurm)
    return _SystemProbe(map(_as_probe_function,
        (version, julia_command, sys_bindir, running_executable, canonical_path,
         file_info, force_numerics, dllist, julia_threads, blas_threads, sysimage,
         project, manifest, environment, hostname, kernel, cpu, microcode, lscpu, slurm))...)
end

function _call(f::Function, label::String)
    try
        return f()
    catch e
        e isa IdentityError && rethrow()
        _fail(:unavailable, "$label is unavailable: $(sprint(showerror, e))")
    end
end

function _payload_text(x, label)
    x === nothing && _fail(:unavailable, "$label is unavailable")
    x isa AbstractVector{UInt8} && return String(x)
    return _text(x, label)
end

function _probe_info(probe::_SystemProbe, path::AbstractString)
    path = String(path)
    x = try probe.file_info(path) catch e
        _fail(:unavailable, "file information is unavailable for $path: $(sprint(showerror, e))")
    end
    d = x isa AbstractDict ? Dict{String,Any}(String(k) => v for (k, v) in x) :
        (x isa NamedTuple ? Dict{String,Any}(String(k) => getfield(x, k) for k in keys(x)) :
         _fail(:probe, "file_info must return a dictionary or named tuple"))
    _exact_keys(d, ["path", "size", "sha256"], "file information")
    d["path"] isa AbstractString || _fail(:probe, "file_info path must be a string")
    d["size"] isa Integer || _fail(:probe, "file_info size must be an integer")
    d["sha256"] isa AbstractString || _fail(:probe, "file_info sha256 must be a string")
    d["path"] = String(d["path"])
    d["size"] = Int64(d["size"])
    d["sha256"] = _hex(d["sha256"], "file sha256")
    return d
end

function _probe_path(probe::_SystemProbe, path, label)
    p = _text(path, label)
    q = try probe.canonical_path(p) catch e
        _fail(:missing_path, "$label cannot be resolved: $(sprint(showerror, e))")
    end
    return _text(q, "$label canonical path")
end

function _executable_identity(probe::_SystemProbe, supplied)
    cmd = _call(probe.julia_command, "Base.julia_cmd().exec[1]")
    cmd isa Cmd && (cmd = cmd.exec[1])
    cmd isa AbstractVector && (cmd = cmd[1])
    bindir = _call(probe.sys_bindir, "Sys.BINDIR")
    cmdtext = _text(cmd, "Base.julia_cmd().exec[1]")
    b = _probe_path(probe, bindir, "Sys.BINDIR")
    cmdpath = isabspath(cmdtext) ? cmdtext : joinpath(b, cmdtext)
    runpath = _probe_path(probe, _call(probe.running_executable, "/proc/self/exe"), "running executable")
    candidates = unique([_probe_path(probe, cmdpath, "Base.julia_cmd executable"),
                         _probe_path(probe, joinpath(b, basename(cmdpath)), "Sys.BINDIR executable"),
                         runpath])
    length(candidates) == 1 || _fail(:runtime_executable, "runtime executable paths disagree: $candidates")
    path = only(candidates)
    binding = _probe_info(probe, path)
    String(binding["path"]) == path || _fail(:runtime_executable, "executable file_info path differs")
    if supplied !== nothing
        sp = _probe_path(probe, supplied, "supplied executable")
        sp == path || _fail(:runtime_executable, "supplied executable path differs")
        sb = _probe_info(probe, sp)
        sb["size"] == binding["size"] && sb["sha256"] == binding["sha256"] ||
            _fail(:runtime_executable, "supplied executable bytes differ")
    end
    return binding
end

function validate_runtime_identity(; current_version=VERSION, runtime_label,
                                   runtime_version, executable=nothing,
                                   detected_executable=nothing)
    actual_version = string(VERSION)
    string(current_version) == actual_version ||
        _fail(:runtime_version, "current VERSION is $actual_version, not $current_version")
    string(runtime_version) == actual_version ||
        _fail(:runtime_version, "supplied runtime version differs from VERSION")
    expected_label = "julia-" * actual_version
    String(runtime_label) == expected_label ||
        _fail(:runtime_label, "supplied runtime label must be $expected_label")
    if executable !== nothing && detected_executable !== nothing
        a = _canonical(String(executable)); b = _canonical(String(detected_executable))
        a == b || _fail(:runtime_executable, "supplied executable path differs")
        sha256_file(a) == sha256_file(b) || _fail(:runtime_executable, "supplied executable hash differs")
        stat(a).size == stat(b).size || _fail(:runtime_executable, "supplied executable size differs")
    else
        # Unlike a caller-supplied version string, the executable is always
        # checked against the live command/bindir/proc identity when this
        # public validator is used directly.
        _executable_identity(_SystemProbe(), executable)
    end
    return true
end

function _runtime_identity(probe::_SystemProbe, supplied_executable, supplied_label, supplied_version)
    v = string(_call(probe.version, "current VERSION"))
    v == string(VERSION) || _fail(:runtime_version, "probe VERSION differs from current VERSION")
    rv = _text(supplied_version, "supplied runtime version")
    rv == v || _fail(:runtime_version, "supplied runtime version differs from current VERSION")
    rl = _text(supplied_label, "supplied runtime label")
    rl == "julia-" * v || _fail(:runtime_label, "supplied runtime label differs from current VERSION")
    exe = _executable_identity(probe, supplied_executable)
    return Dict{String,Any}(
        "current_version" => v,
        "runtime_label" => rl,
        "runtime_version" => rv,
        "executable" => exe,
    )
end

function _actual_dsos(probe::_SystemProbe)
    raw = _call(probe.dllist, "Libdl.dllist()")
    raw isa AbstractVector || _fail(:dso_probe, "Libdl.dllist() did not return a vector")
    paths = String[]
    for item in raw
        item === nothing && continue
        text = String(item)
        isfile(text) || continue # linux-vdso and similar loader handles are not DSOs
        p = _probe_path(probe, text, "loaded DSO")
        p in paths || push!(paths, p)
    end
    sort!(paths)
    info = Dict{String,Any}(p => _probe_info(probe, p) for p in paths)
    return paths, info
end

function _environment_contract(probe::_SystemProbe, project::String, nthreads::Int)
    raw = _call(probe.environment, "environment")
    raw isa AbstractDict || _fail(:environment, "environment probe did not return a dictionary")
    env = Dict{String,String}(String(k) => String(v) for (k, v) in raw)
    depot = get(env, "JULIA_DEPOT_PATH", join(String.(Base.DEPOT_PATH), Sys.iswindows() ? ';' : ':'))
    load = get(env, "JULIA_LOAD_PATH", join(String.(LOAD_PATH), Sys.iswindows() ? ';' : ':'))
    project_env = get(env, "JULIA_PROJECT", project)
    history = get(env, "JULIA_HISTORY", "default")
    startup = get(env, "JULIA_STARTUP_FILE", "default")
    thread_env = get(env, "JULIA_NUM_THREADS", string(nthreads))
    temp = get(env, "TMPDIR", tempdir())
    compiled = get(env, "JULIA_COMPILED_MODULES", string(Base.JLOptions().use_compiled_modules))
    fields = Dict{String,Any}(
        "depot_path" => _text(depot, "JULIA_DEPOT_PATH"),
        "load_path" => _text(load, "JULIA_LOAD_PATH"),
        "project" => _text(project_env, "JULIA_PROJECT"),
        "history" => _text(history, "JULIA_HISTORY"),
        "startup" => _text(startup, "JULIA_STARTUP_FILE"),
        "threads" => _text(thread_env, "JULIA_NUM_THREADS"),
        "temp" => _text(temp, "temporary directory"),
        "compiled_modules" => _text(compiled, "JULIA_COMPILED_MODULES"),
    )
    fields["threads"] == "1" || _fail(:threads, "JULIA_NUM_THREADS is not exactly one")
    fields["project"] == project || _fail(:environment, "JULIA_PROJECT differs from active project")
    fields["compiled_modules"] in ("0", "1") || _fail(:environment, "invalid compiled-modules setting")
    return fields
end

function _host_contract(probe::_SystemProbe)
    hostname = _text(_call(probe.hostname, "hostname"), "hostname")
    kernel = _payload_text(_call(probe.kernel, "kernel"), "kernel")
    cpu = _payload_text(_call(probe.cpu, "CPU information"), "CPU information")
    microcode = _payload_text(_call(probe.microcode, "microcode"), "microcode")
    lscpu = _call(probe.lscpu, "lscpu")
    lscpu === nothing && _fail(:unavailable, "lscpu is unavailable")
    lbytes = lscpu isa AbstractVector{UInt8} ? Vector{UInt8}(lscpu) : Vector{UInt8}(codeunits(String(lscpu)))
    isempty(lbytes) && _fail(:unavailable, "lscpu output is empty")
    return Dict{String,Any}(
        "hostname" => hostname, "hostname_sha256" => _sha256_bytes(codeunits(hostname)),
        "kernel" => kernel, "kernel_sha256" => _sha256_bytes(codeunits(kernel)),
        "cpu" => cpu, "cpu_sha256" => _sha256_bytes(codeunits(cpu)),
        "microcode" => microcode, "microcode_sha256" => _sha256_bytes(codeunits(microcode)),
        "lscpu_sha256" => _sha256_bytes(lbytes), "lscpu_size" => length(lbytes),
    )
end

function _slurm_contract(probe::_SystemProbe, hostname::String)
    x = _call(probe.slurm, "Slurm state")
    x isa AbstractDict || _fail(:slurm, "Slurm probe did not return a dictionary")
    fields = Dict{String,Any}()
    for key in ("job_id", "node", "task_id", "array_job_id", "array_task_id")
        haskey(x, key) || _fail(:slurm, "Slurm field is missing: $key")
        fields[key] = _text(x[key], "Slurm $key")
    end
    fields["node"] == hostname || _fail(:slurm_node, "Slurm node does not match hostname")
    return fields
end

function _binding_from_path(probe::_SystemProbe, supplied, actual::String, label)
    p = _probe_path(probe, supplied, label)
    p == actual || _fail(:path, "$label differs from the active path")
    return p
end

function _capture_with_probe(probe::_SystemProbe; dso_contract=nothing,
    dso_contract_path=nothing, runtime_label=nothing, runtime_version=nothing,
    executable=nothing, project_path=nothing, manifest_path=nothing,
    sysimage_path=nothing)
    contract = dso_contract === nothing ? parse_dso_contract(dso_contract_path) : dso_contract
    contract isa DSOContract || _fail(:dso_schema, "dso_contract is not a DSOContract")
    rl = runtime_label === nothing ? contract.runtime_label : runtime_label
    rv = runtime_version === nothing ? contract.runtime_version : runtime_version
    runtime = _runtime_identity(probe, executable, rl, rv)
    # Numerical operations are intentionally before dllist(), not merely an
    # import-time observation.
    _call(probe.force_numerics, "forced numerical operations")
    paths, info = _actual_dsos(probe)
    dso = _validate_dso!(contract, paths, info, runtime["runtime_label"], runtime["runtime_version"])

    nt = Int(_call(probe.julia_threads, "Julia thread count"))
    nt == 1 || _fail(:threads, "Julia thread count is $nt, expected exactly one")
    bt = Int(_call(probe.blas_threads, "BLAS thread count"))
    bt == 1 || _fail(:blas_threads, "BLAS thread count is $bt, expected exactly one")

    actual_project = _probe_path(probe, _call(probe.project, "active project"), "active project")
    actual_manifest = _probe_path(probe, _call(probe.manifest, "Manifest.toml"), "Manifest.toml")
    project_path !== nothing && _binding_from_path(probe, project_path, actual_project, "supplied project")
    manifest_path !== nothing && _binding_from_path(probe, manifest_path, actual_manifest, "supplied manifest")
    sysimage = _probe_path(probe, _call(probe.sysimage, "sysimage"), "sysimage")
    sysimage_path !== nothing && _binding_from_path(probe, sysimage_path, sysimage, "supplied sysimage")

    project_binding = _probe_info(probe, actual_project)
    manifest_binding = _probe_info(probe, actual_manifest)
    sysimage_binding = _probe_info(probe, sysimage)
    host = _host_contract(probe)
    environment = _environment_contract(probe, actual_project, nt)
    slurm = _slurm_contract(probe, host["hostname"])
    report = Dict{String,Any}(
        "schema" => IDENTITY_SCHEMA,
        "runtime" => runtime,
        "dso" => dso,
        "threads" => Dict{String,Any}("julia" => nt, "blas" => bt),
        "sysimage" => sysimage_binding,
        "project" => project_binding,
        "manifest" => manifest_binding,
        "environment" => environment,
        "host" => host,
        "slurm" => slurm,
    )
    validate_identity(report)
    return report
end

"""Capture using only the current process and operating-system observations."""
function capture_identity(; dso_contract=nothing, dso_contract_path=nothing,
    contract_path=nothing, runtime_label=nothing, runtime_version=nothing,
    executable=nothing, project_path=nothing, manifest_path=nothing, sysimage_path=nothing)
    path = dso_contract_path === nothing ? contract_path : dso_contract_path
    path === nothing && dso_contract === nothing && _fail(:arguments, "a DSO contract path is required")
    runtime_label === nothing && _fail(:arguments, "runtime_label is required")
    runtime_version === nothing && _fail(:arguments, "runtime_version is required")
    executable === nothing && _fail(:arguments, "executable is required")
    # The public production entrypoint deliberately has no probe keyword.
    return _capture_with_probe(_SystemProbe(); dso_contract=dso_contract,
        dso_contract_path=path, runtime_label=runtime_label,
        runtime_version=runtime_version, executable=executable,
        project_path=project_path, manifest_path=manifest_path,
        sysimage_path=sysimage_path)
end

# --------------------------------------------------------------------------
# Report validation and publication

const _REPORT_KEYS = ["schema", "runtime", "dso", "threads", "sysimage", "project",
                      "manifest", "environment", "host", "slurm"]

function _validate_binding_dict(d, label)
    d isa AbstractDict || _fail(:identity_schema, "$label is not a table")
    _exact_keys(d, ["path", "size", "sha256"], label)
    d["path"] isa AbstractString || _fail(:identity_schema, "$label path is not a string")
    d["size"] isa Integer || _fail(:identity_schema, "$label size is not an integer")
    d["sha256"] isa AbstractString || _fail(:identity_schema, "$label sha256 is not a string")
    p = _canonical(String(d["path"]))
    p == String(d["path"]) || _fail(:identity_path, "$label is not canonical")
    st = stat(p)
    Int64(d["size"]) == st.size || _fail(:identity_hash, "$label size changed")
    lowercase(String(d["sha256"])) == sha256_file(p) || _fail(:identity_hash, "$label hash changed")
    _hex(d["sha256"], "$label sha256")
    return true
end

function validate_identity(report::AbstractDict; runtime_label=nothing,
    runtime_version=nothing, executable=nothing, project_path=nothing,
    manifest_path=nothing, sysimage_path=nothing)
    _exact_keys(report, _REPORT_KEYS, "identity report")
    report["schema"] == IDENTITY_SCHEMA || _fail(:identity_schema, "unsupported identity report schema")
    rt = report["runtime"]
    rt isa AbstractDict || _fail(:identity_schema, "runtime is not a table")
    _exact_keys(rt, ["current_version", "runtime_label", "runtime_version", "executable"], "runtime")
    v = _text(rt["current_version"], "current VERSION")
    _text(rt["runtime_version"], "runtime version") == v || _fail(:runtime_version, "runtime report version mismatch")
    _text(rt["runtime_label"], "runtime label") == "julia-" * v || _fail(:runtime_label, "runtime report label mismatch")
    runtime_label !== nothing && String(runtime_label) == rt["runtime_label"] || runtime_label === nothing ||
        _fail(:runtime_label, "supplied runtime label differs")
    runtime_version !== nothing && String(runtime_version) == rt["runtime_version"] || runtime_version === nothing ||
        _fail(:runtime_version, "supplied runtime version differs")
    _validate_binding_dict(rt["executable"], "runtime executable")
    executable !== nothing && _canonical(String(executable)) == rt["executable"]["path"] || executable === nothing ||
        _fail(:runtime_executable, "supplied executable differs")
    for key in ("sysimage", "project", "manifest")
        _validate_binding_dict(report[key], key)
    end
    project_path !== nothing && _canonical(String(project_path)) == report["project"]["path"] || project_path === nothing ||
        _fail(:path, "supplied project differs")
    manifest_path !== nothing && _canonical(String(manifest_path)) == report["manifest"]["path"] || manifest_path === nothing ||
        _fail(:path, "supplied manifest differs")
    sysimage_path !== nothing && _canonical(String(sysimage_path)) == report["sysimage"]["path"] || sysimage_path === nothing ||
        _fail(:path, "supplied sysimage differs")
    th = report["threads"]
    th isa AbstractDict || _fail(:identity_schema, "threads is not a table")
    _exact_keys(th, ["julia", "blas"], "threads")
    th["julia"] == 1 || _fail(:threads, "identity does not contain exactly one Julia thread")
    th["blas"] == 1 || _fail(:blas_threads, "identity does not contain exactly one BLAS thread")
    env = report["environment"]
    env isa AbstractDict || _fail(:identity_schema, "environment is not a table")
    _exact_keys(env, ["depot_path", "load_path", "project", "history", "startup", "threads", "temp", "compiled_modules"], "environment")
    for (k, value) in env
        _text(value, "environment.$k")
    end
    env["threads"] == "1" || _fail(:threads, "environment thread contract is not one")
    env["project"] == report["project"]["path"] || _fail(:environment, "environment project differs")
    env["compiled_modules"] in ("0", "1") || _fail(:environment, "invalid compiled module contract")
    host = report["host"]
    host isa AbstractDict || _fail(:identity_schema, "host is not a table")
    _exact_keys(host, ["hostname", "hostname_sha256", "kernel", "kernel_sha256", "cpu", "cpu_sha256", "microcode", "microcode_sha256", "lscpu_sha256", "lscpu_size"], "host")
    for key in ("hostname", "kernel", "cpu", "microcode")
        _text(host[key], "host.$key")
        _hex(host[key * "_sha256"], "host.$key sha256")
    end
    _hex(host["lscpu_sha256"], "host lscpu sha256")
    Int(host["lscpu_size"]) > 0 || _fail(:unavailable, "lscpu output is empty")
    sl = report["slurm"]
    sl isa AbstractDict || _fail(:identity_schema, "slurm is not a table")
    _exact_keys(sl, ["job_id", "node", "task_id", "array_job_id", "array_task_id"], "slurm")
    for key in keys(sl); _text(sl[key], "slurm.$key"); end
    sl["node"] == host["hostname"] || _fail(:slurm_node, "identity Slurm node differs from hostname")
    dso = report["dso"]
    dso isa AbstractDict || _fail(:identity_schema, "dso is not a table")
    _exact_keys(dso, ["schema", "runtime_label", "runtime_version", "common_count",
                      "versioned_count", "actual_count", "common_providers",
                      "versioned_julia_components"], "dso")
    for key in ("runtime_label", "runtime_version", "schema")
        _text(dso[key], "dso.$key")
    end
    dso["runtime_label"] == rt["runtime_label"] && dso["runtime_version"] == rt["runtime_version"] ||
        _fail(:dso_runtime, "identity DSO runtime differs")
    Int(dso["actual_count"]) == Int(dso["common_count"]) + Int(dso["versioned_count"]) ||
        _fail(:dso_count, "identity DSO counts do not add up")
    dso["common_providers"] isa AbstractVector || _fail(:identity_schema, "common providers are not an array")
    dso["versioned_julia_components"] isa AbstractVector || _fail(:identity_schema, "versioned components are not an array")
    length(dso["common_providers"]) == Int(dso["common_count"]) || _fail(:dso_count, "common provider count differs")
    length(dso["versioned_julia_components"]) == Int(dso["versioned_count"]) || _fail(:dso_count, "versioned component count differs")
    paths = String[]
    for (i, item) in enumerate(dso["common_providers"])
        item isa AbstractDict || _fail(:identity_schema, "common provider $i is not a table")
        _exact_keys(item, ["name", "path", "size", "sha256"], "dso.common_providers[$i]")
        item["name"] == basename(String(item["path"])) || _fail(:dso_schema, "common provider basename differs")
        _validate_binding_dict(Dict{String,Any}("path"=>item["path"], "size"=>item["size"], "sha256"=>item["sha256"]), "common provider $i")
        push!(paths, String(item["path"]))
    end
    for (i, item) in enumerate(dso["versioned_julia_components"])
        item isa AbstractDict || _fail(:identity_schema, "versioned component $i is not a table")
        _exact_keys(item, ["name", "path", "size", "sha256", "runtime_label", "runtime_version"], "dso.versioned_julia_components[$i]")
        item["runtime_label"] == rt["runtime_label"] && item["runtime_version"] == rt["runtime_version"] ||
            _fail(:dso_runtime, "versioned component runtime differs")
        item["name"] == basename(String(item["path"])) || _fail(:dso_schema, "versioned component basename differs")
        _validate_binding_dict(Dict{String,Any}("path"=>item["path"], "size"=>item["size"], "sha256"=>item["sha256"]), "versioned component $i")
        push!(paths, String(item["path"]))
    end
    length(unique(paths)) == length(paths) || _fail(:dso_schema, "identity DSO paths are duplicated")
    return true
end

function identity_summary(report::AbstractDict)
    validate_identity(report)
    return Dict{String,Any}(
        "runtime_version" => report["runtime"]["runtime_version"],
        "runtime_executable_sha256" => report["runtime"]["executable"]["sha256"],
        "sysimage_sha256" => report["sysimage"]["sha256"],
        "project_sha256" => report["project"]["sha256"],
        "manifest_sha256" => report["manifest"]["sha256"],
        "dso_common_count" => report["dso"]["common_count"],
        "dso_versioned_count" => report["dso"]["versioned_count"],
        "dso_actual_count" => report["dso"]["actual_count"],
    )
end

function _manifest(identity_path::String; root_mode=_mode_string(_file_mode(dirname(identity_path))))
    st = stat(identity_path)
    return Dict{String,Any}(
        "schema" => MANIFEST_SCHEMA,
        "files" => [Dict{String,Any}(
            "path" => basename(identity_path), "size" => Int64(st.size),
            "sha256" => sha256_file(identity_path), "mode" => _mode_string(_file_mode(identity_path)),
        )],
        "root_mode" => root_mode,
    )
end

function _rename_noreplace(src::String, dst::String)
    ispath(dst) && _fail(:root_collision, "publication root already exists: $dst")
    if isdefined(Base.Filesystem, :_mv_noreplace)
        try
            Base.Filesystem._mv_noreplace(src, dst)
        catch e
            ispath(dst) && _fail(:root_collision, "publication root appeared during publication: $dst")
            _fail(:root_publish, "atomic no-replace root publish failed: $(sprint(showerror, e))")
        end
    else
        # This branch is only for old/non-Unix Julia builds.  The pre-check is
        # retained as a defensive fallback; supported Linux uses _mv_noreplace.
        ispath(dst) && _fail(:root_collision, "publication root already exists: $dst")
        Base.Filesystem.rename(src, dst)
    end
    return nothing
end

function _write_exclusive(path::String, value::Dict{String,Any})
    io = IOBuffer(); TOML.print(io, value); write(io, '\n')
    payload = take!(io)
    flags = Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL
    f = try
        Base.Filesystem.open(path, flags, 0o644)
    catch e
        _fail(:receipt_collision, "external receipt was not created exclusively: $(sprint(showerror, e))")
    end
    try
        write(f, payload)
        _fsync_fd(Base.fd(f))
    finally
        close(f)
    end
    _fsync_directory(dirname(path))
    return path
end

function _fault_point(fault, point::Symbol)
    fault === nothing && return nothing
    result = fault(point)
    result === nothing && return nothing
    result isa Exception && throw(result)
    _fail(:injected_fault, "injected publication fault at $point")
end

"""
Publish an already validated report using the receipt-last protocol.  The
root is assembled and fsynced under a sibling staging name, moved with
no-replace semantics, and only then is the external O_EXCL receipt created.
There is intentionally no operation on the root after receipt creation.
"""
function publish_identity(report::AbstractDict, root::AbstractString;
                          receipt=nothing, receipt_path=nothing, _fault=nothing)
    validate_identity(report)
    target = abspath(String(root))
    parent = dirname(target)
    isdir(parent) || _fail(:publication_parent, "publication parent does not exist: $parent")
    ispath(target) && _fail(:root_collision, "publication root already exists: $target")
    rpath = receipt === nothing ? receipt_path : receipt
    rpath = rpath === nothing ? joinpath(parent, basename(target) * ".receipt.toml") : abspath(String(rpath))
    startswith(rpath, target * "/") && _fail(:receipt_path, "receipt must be external to publication root")
    stage = joinpath(parent, "." * basename(target) * ".stage-" * string(uuid4()))
    mkdir(stage, mode=0o755)
    root_published = false
    try
        identity_path = joinpath(stage, "identity.toml")
        _fault_point(_fault, :before_identity)
        _write_toml(identity_path, Dict{String,Any}(String(k) => v for (k, v) in report))
        chmod(identity_path, 0o644)
        manifest_path = joinpath(stage, "manifest.toml")
        manifest = _manifest(identity_path; root_mode=_mode_string(_file_mode(stage)))
        _fault_point(_fault, :before_manifest)
        _write_toml(manifest_path, manifest)
        chmod(manifest_path, 0o644)
        _fsync_file(identity_path); _fsync_file(manifest_path); _fsync_directory(stage)
        _fault_point(_fault, :before_root_rename)
        _rename_noreplace(stage, target)
        root_published = true
        _fsync_directory(parent)
        receipt_data = Dict{String,Any}(
            "schema" => RECEIPT_SCHEMA,
            "root" => target,
            "identity_sha256" => sha256_file(joinpath(target, "identity.toml")),
            "manifest_sha256" => sha256_file(joinpath(target, "manifest.toml")),
            "file_count" => 2,
            "root_mode" => _mode_string(_file_mode(target)),
            "published_utc" => string(now(UTC)),
        )
        # This is the terminal filesystem action.  No root file is touched
        # after this call.
        _fault_point(_fault, :before_receipt)
        _write_exclusive(rpath, receipt_data)
        return Dict{String,Any}("root" => target, "receipt" => rpath,
                                "identity_sha256" => receipt_data["identity_sha256"],
                                "manifest_sha256" => receipt_data["manifest_sha256"],
                                "file_count" => 2)
    catch e
        !root_published && ispath(stage) && rm(stage; recursive=true, force=true)
        e isa IdentityError && rethrow()
        _fail(:publication, sprint(showerror, e))
    end
end

end # module AuthorityIdentity
