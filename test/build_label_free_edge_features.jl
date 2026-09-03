#!/usr/bin/env julia

using SHA

include(joinpath(@__DIR__, "lib", "structured_assignment", "edge_features.jl"))
using .StructuredEdgeFeatures

const VALUE_FLAGS = Set([
    "--root",
    "--features",
    "--feature-sha256",
    "--candidate-config",
    "--model-config",
    "--universe-dir",
    "--forward-patches",
    "--backward-patches",
    "--forward-receipt",
    "--backward-receipt",
    "--out-dir",
])
const ALLOWED_FLAGS = union(VALUE_FLAGS, Set(["--help"]))
const REQUIRED_FLAGS = sort!(collect(VALUE_FLAGS))

function usage(io::IO=stdout)
    println(io, "Usage: julia --project=. test/build_label_free_edge_features.jl \\")
    println(io, "  --root PATH --features PATH --feature-sha256 HEX \\")
    println(io, "  --candidate-config PATH --model-config PATH --universe-dir PATH \\")
    println(io, "  --forward-patches PATH --backward-patches PATH \\")
    println(io, "  --forward-receipt PATH --backward-receipt PATH --out-dir PATH")
    return nothing
end

function parse_arguments(arguments::Vector{String})::Dict{String,String}
    StructuredEdgeFeatures.StructuredUniverse.InputBoundary.validate_cli_arguments(
        arguments;
        allowed_flags=ALLOWED_FLAGS,
        value_flags=VALUE_FLAGS,
    )
    values = Dict{String,String}()
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        pair = split(argument, '='; limit=2)
        flag = pair[1]
        if flag == "--help"
            values[flag] = "true"
            index += 1
        elseif length(pair) == 2
            values[flag] = pair[2]
            index += 1
        else
            values[flag] = arguments[index + 1]
            index += 2
        end
    end
    if haskey(values, "--help")
        length(values) == 1 || throw(EdgeFeatureError(
            :cli_error,
            "--help accepts no construction inputs",
        ))
        return values
    end
    for flag in REQUIRED_FLAGS
        haskey(values, flag) || throw(EdgeFeatureError(
            :missing_argument,
            "required option is absent: $flag",
        ))
    end
    return values
end

function main(arguments::Vector{String}=copy(ARGS))::Int
    values = parse_arguments(arguments)
    if haskey(values, "--help")
        usage()
        return 0
    end
    hashes = publish_edge_bundle(
        values["--root"],
        values["--out-dir"];
        features=values["--features"],
        feature_sha256=values["--feature-sha256"],
        candidate_config=values["--candidate-config"],
        model_config=values["--model-config"],
        universe_dir=values["--universe-dir"],
        forward_patches=values["--forward-patches"],
        backward_patches=values["--backward-patches"],
        forward_receipt=values["--forward-receipt"],
        backward_receipt=values["--backward-receipt"],
    )
    canonical = join(["$name=$(hashes[name])" for name in sort!(collect(keys(hashes)))], "\n") * "\n"
    println("status=PASS")
    println("files=", length(hashes))
    println("bundle_sha256=", bytes2hex(sha256(codeunits(canonical))))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(main())
    catch error
        if error isa EdgeFeatureError ||
           error isa StructuredEdgeFeatures.StructuredUniverse.UniverseError ||
           error isa StructuredEdgeFeatures.StructuredUniverse.InputBoundary.StructuredFirewallError
            showerror(stderr, error)
            println(stderr)
            exit(2)
        end
        rethrow()
    end
end
