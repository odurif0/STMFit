#!/usr/bin/env julia

module StructuredUnaryPredictionCLI

include(joinpath(@__DIR__, "lib", "structured_unit_assignment.jl"))
using .StructuredUnitAssignment
const SUA = StructuredUnitAssignment

export Options, parse_cli, show_help, main, entrypoint

struct Options
    config::String
    features::String
    out::String
    model::String
end

const ALLOWED_FLAGS = Set(["--config", "--features", "--out", "--model", "--help"])
const VALUE_FLAGS = Set(["--config", "--features", "--out", "--model"])

function _option_values(arguments::Vector{String})
    values = Dict{String,String}()
    index = 1
    while index <= length(arguments)
        pair = split(arguments[index], '='; limit=2)
        flag = pair[1]
        if flag == "--help"
            index += 1
        elseif length(pair) == 2
            values[flag] = pair[2]
            index += 1
        else
            values[flag] = arguments[index + 1]
            index += 2
        end
    end
    return values
end

function parse_cli(arguments::Vector{String})
    SUA.StructuredFirewall.validate_cli_arguments(
        arguments;
        allowed_flags=ALLOWED_FLAGS,
        value_flags=VALUE_FLAGS,
    )
    "--help" in arguments &&
        throw(ArgumentError("--help must be used by itself"))
    values = _option_values(arguments)
    for flag in ("--config", "--features", "--out", "--model")
        haskey(values, flag) || throw(ArgumentError("$flag is required"))
    end
    model = values["--model"]
    model in SUA.SUPPORTED_MODELS ||
        throw(ArgumentError("--model must name a configured unary mode"))
    return Options(values["--config"], values["--features"], values["--out"], model)
end

function show_help(io::IO=stdout)
    println(io, """
    Usage: julia --project=. test/build_structured_unit_predictions.jl \\
      --config PATH --features PATH --out PATH --model MODE

    Label-free structured unary predictor with fixed equal view fusion.

    Required:
      --config PATH    Immutable structured model TOML
      --features PATH  Per-lobe physical feature TSV
      --out PATH       Exact-schema prediction TSV
      --model MODE     gaussian_one | gaussian_two | student_t_two

    Optional:
      --help           Show this text
    """)
    return nothing
end

function main(arguments::Vector{String}=copy(ARGS))
    if arguments == ["--help"]
        show_help()
        return 0
    end
    options = parse_cli(arguments)
    bundle = SUA.build_predictions(options.config, options.features, options.model)
    destination = SUA.write_prediction_tsv(options.out, bundle)
    println("structured unary status=ok")
    println("model=", options.model)
    println("rows=", length(bundle.rows))
    println("out=", destination)
    return 0
end

function entrypoint(arguments::Vector{String}=copy(ARGS))
    try
        return main(arguments)
    catch error
        print(stderr, "BLOCKED: ")
        showerror(stderr, error)
        println(stderr)
        return 2
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StructuredUnaryPredictionCLI.entrypoint(copy(ARGS)))
end
