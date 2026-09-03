#!/usr/bin/env julia

# This is intentionally an integration entrypoint, not a package wrapper.
# AuthorityIdentity is included directly and the private test probe seam is
# not reachable through these arguments.
include(joinpath(@__DIR__, "AuthorityIdentity.jl"))
using .AuthorityIdentity

const _OPTIONS = Set([
    "--dso-contract", "--contract", "--output", "--receipt",
    "--runtime-label", "--runtime-version", "--executable",
    "--project", "--manifest", "--sysimage", "--dry-run",
])

function _usage(io::IO=stdout)
    println(io, "Usage: julia integration/identity_probe.jl --dso-contract PATH --output ROOT")
    println(io, "       --runtime-label julia-VERSION --runtime-version VERSION --executable PATH")
    println(io, "       [--project PATH --manifest PATH --sysimage PATH --receipt PATH]")
    println(io, "       --dry-run [schema|full]")
    println(io, "\nThe normal path captures only the current process/OS identity and publishes")
    println(io, "a receipt-last identity root. It never accepts an injected system probe.")
end

function _parse(args)
    out = Dict{String,String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        arg == "--help" && return Dict("help" => "true")
        startswith(arg, "--") || error("unexpected positional argument: $arg")
        key, value = if occursin('=', arg)
            split(arg, "="; limit=2)
        else
            arg, nothing
        end
        key in _OPTIONS || error("unknown option: $key")
        if key == "--dry-run" && value === nothing
            if i < length(args) && !startswith(args[i + 1], "--")
                i += 1; value = args[i]
            else
                value = "schema"
            end
        elseif value === nothing
            i == length(args) && error("missing value for $key")
            i += 1
            value = args[i]
            startswith(value, "--") && error("missing value for $key")
        end
        key == "--contract" && (key = "--dso-contract")
        haskey(out, key) && error("duplicate option: $key")
        out[key] = value
        i += 1
    end
    return out
end

function _required(opts, key)
    haskey(opts, key) || error("missing required option: $key")
    isempty(strip(opts[key])) && error("empty required option: $key")
    return opts[key]
end

function _schema_dry_run()
    println(stdout, "schema=", AuthorityIdentity.IDENTITY_SCHEMA)
    println(stdout, "dso_schema=", AuthorityIdentity.DSO_SCHEMA)
    println(stdout, "manifest_schema=", AuthorityIdentity.MANIFEST_SCHEMA)
    println(stdout, "receipt_schema=", AuthorityIdentity.RECEIPT_SCHEMA)
    println(stdout, "api=capture_identity,validate_identity,publish_identity,parse_dso_contract")
    return 0
end

function _full_dry_run(opts)
    contract_path = _required(opts, "--dso-contract")
    label = _required(opts, "--runtime-label")
    version = _required(opts, "--runtime-version")
    executable = _required(opts, "--executable")
    contract = parse_dso_contract(contract_path)
    contract.runtime_label == label || error("runtime label does not match the DSO contract")
    contract.runtime_version == version || error("runtime version does not match the DSO contract")
    println(stdout, "dry_run=full")
    println(stdout, "runtime_label=", label)
    println(stdout, "runtime_version=", version)
    println(stdout, "executable=", executable)
    println(stdout, "dso_common_count=", length(contract.common))
    println(stdout, "dso_versioned_count=", length(contract.versioned))
    println(stdout, "writes=none")
    return 0
end

function main(args=ARGS)
    opts = try
        _parse(args)
    catch e
        println(stderr, "BLOCKED arguments: ", e)
        return 2
    end
    haskey(opts, "help") && (_usage(); return 0)
    if haskey(opts, "--dry-run")
        mode = lowercase(opts["--dry-run"])
        mode == "schema" && return _schema_dry_run()
        mode == "full" && return try _full_dry_run(opts) catch e
            println(stderr, "BLOCKED dry-run: ", e); 2
        end
        println(stderr, "BLOCKED dry-run mode must be schema or full")
        return 2
    end

    try
        contract = _required(opts, "--dso-contract")
        output = _required(opts, "--output")
        label = _required(opts, "--runtime-label")
        version = _required(opts, "--runtime-version")
        executable = _required(opts, "--executable")
        report = capture_identity(
            dso_contract_path=contract,
            runtime_label=label,
            runtime_version=version,
            executable=executable,
            project_path=get(opts, "--project", nothing),
            manifest_path=get(opts, "--manifest", nothing),
            sysimage_path=get(opts, "--sysimage", nothing),
        )
        publication = publish_identity(report, output; receipt=get(opts, "--receipt", nothing))
        summary = identity_summary(report)
        println(stdout, "root=", publication["root"])
        println(stdout, "receipt=", publication["receipt"])
        println(stdout, "runtime_executable_sha256=", summary["runtime_executable_sha256"])
        println(stdout, "sysimage_sha256=", summary["sysimage_sha256"])
        println(stdout, "project_sha256=", summary["project_sha256"])
        println(stdout, "manifest_sha256=", summary["manifest_sha256"])
        println(stdout, "dso_counts=", summary["dso_common_count"], "/",
                summary["dso_versioned_count"], "/", summary["dso_actual_count"])
        return 0
    catch e
        println(stderr, "BLOCKED identity: ", e)
        return 2
    end
end

exit(main())
