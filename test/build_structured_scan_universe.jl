#!/usr/bin/env julia

using SHA

include(joinpath(@__DIR__, "lib", "structured_assignment", "universe.jl"))
using .StructuredUniverse

const VALUE_FLAGS = Set([
    "--root",
    "--features",
    "--feature-sha256",
    "--candidate-config",
    "--model-config",
    "--out-dir",
])
const ALLOWED_FLAGS = union(VALUE_FLAGS, Set(["--help"]))

function _usage(io::IO=stdout)
    println(io, "Usage: julia --project=. test/build_structured_scan_universe.jl \\")
    println(io, "  --features PATH --feature-sha256 HEX \\")
    println(io, "  --candidate-config PATH --model-config PATH --out-dir PATH [--root PATH]")
    return nothing
end

function _arguments(arguments)
    StructuredUniverse.InputBoundary.validate_cli_arguments(
        arguments;
        allowed_flags=ALLOWED_FLAGS,
        value_flags=VALUE_FLAGS,
    )
    values = Dict{String,String}()
    index = 1
    while index <= length(arguments)
        argument = String(arguments[index])
        pair = split(argument, '='; limit=2)
        flag = pair[1]
        if flag == "--help"
            values[flag] = "true"
            index += 1
            continue
        end
        if length(pair) == 2
            values[flag] = pair[2]
            index += 1
        else
            values[flag] = String(arguments[index + 1])
            index += 2
        end
    end
    return values
end

function main(arguments=ARGS)::Int
    values = _arguments(arguments)
    if haskey(values, "--help")
        _usage()
        return 0
    end
    for flag in ("--features", "--feature-sha256", "--candidate-config", "--model-config", "--out-dir")
        haskey(values, flag) || throw(UniverseError(:missing_argument, "required option is absent: $flag"))
    end
    root = get(values, "--root", normpath(joinpath(@__DIR__, "..")))
    source = joinpath(@__DIR__, "lib", "structured_assignment", "universe.jl")
    boundary = joinpath(@__DIR__, "lib", "structured_assignment", "firewall.jl")
    frozen = freeze_universe(
        root;
        features=values["--features"],
        feature_sha256=values["--feature-sha256"],
        candidate_config=values["--candidate-config"],
        model_config=values["--model-config"],
        source_paths=[source, abspath(@__FILE__), boundary],
    )
    hashes = publish_bundle(frozen, values["--out-dir"])
    println("status=PASS")
    println("files=", length(hashes))
    println("bundle_sha256=", bytes2hex(SHA.sha256(codeunits(join(
        ["$name=$(hashes[name])" for name in sort!(collect(keys(hashes)))], "\n") * "\n"))))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(main())
    catch error
        if error isa UniverseError ||
           error isa StructuredUniverse.InputBoundary.StructuredFirewallError
            showerror(stderr, error)
            println(stderr)
            exit(2)
        end
        rethrow()
    end
end
