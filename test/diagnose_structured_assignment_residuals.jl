#!/usr/bin/env julia

module StructuredDiagnosticsCLI

include(joinpath(@__DIR__, "lib", "structured_assignment", "diagnostics.jl"))
module Firewall
include(joinpath(@__DIR__, "lib", "structured_assignment", "firewall.jl"))
end

const SAD = StructuredAssignmentDiagnostics
const ALLOWED_FLAGS = Set([
    "--config", "--features", "--folds", "--receipt", "--out", "--help",
])
const VALUE_FLAGS = Set(["--config", "--features", "--folds", "--receipt", "--out"])

export Options, parse_cli, show_help, main, entrypoint

struct Options
    config::String
    features::String
    folds::String
    receipt::String
    out::String
end

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
    Firewall.validate_cli_arguments(
        arguments; allowed_flags=ALLOWED_FLAGS, value_flags=VALUE_FLAGS)
    "--help" in arguments && throw(ArgumentError("--help must be used by itself"))
    values = _option_values(arguments)
    for flag in ("--config", "--features", "--folds", "--receipt", "--out")
        haskey(values, flag) || throw(ArgumentError("$flag is required"))
    end
    return Options(values["--config"], values["--features"], values["--folds"],
                   values["--receipt"], values["--out"])
end

function show_help(io::IO=stdout)
    println(io, """
    Usage: julia --project=. test/diagnose_structured_assignment_residuals.jl \\
      --config PATH --features PATH --folds PATH --receipt PATH --out PATH

    Label-free, follow-up-only structured residual diagnostics.

    Required:
      --config PATH   Immutable structured model TOML
      --features PATH Exact configured physical-feature TSV
      --folds PATH    Complete nested date-fold TSV
      --receipt PATH  Hash-bound edge-feasibility input receipt
      --out PATH      Deterministic diagnostic evidence TSV

    Optional:
      --help          Show this text
    """)
    return nothing
end

function main(options::Options)
    run = SAD.run_diagnostics(
        options.config, options.features, options.folds, options.receipt)
    bytes = SAD.report_tsv_bytes(run)
    destination = SAD.write_report(
        options.out, bytes; snapshots=run.snapshots,
        expected_source_hash=run.source_hash)
    println("structured diagnostic status=", run.status)
    println("reason=", run.reason)
    println("null_campaigns=", run.null_campaigns.campaigns)
    println("null_false_triggers=", run.null_campaigns.failures)
    println("null_clopper_pearson_upper=", run.null_campaigns.clopper_pearson_upper)
    println("run_sha256=", run.run_sha256)
    println("output_sha256=", SAD.sha256_hex(bytes))
    println("out=", destination)
    return run.status == "BLOCKED" ? 2 : 0
end

function main(arguments::Vector{String}=copy(ARGS))
    if arguments == ["--help"]
        show_help()
        return 0
    end
    return main(parse_cli(arguments))
end

function entrypoint(arguments::Vector{String}=copy(ARGS))
    if arguments == ["--help"]
        show_help()
        return 0
    end
    options = try
        parse_cli(arguments)
    catch error
        print(stderr, "BLOCKED: ")
        showerror(stderr, error)
        println(stderr)
        return 2
    end
    try
        return main(options)
    catch error
        digest_or_na(path::String) = try
            isfile(path) ? SAD.sha256_hex(read(path)) : "NA"
        catch
            "NA"
        end
        source_or_na = try
            SAD.model_source_hash()
        catch
            "NA"
        end
        blocked_bytes = SAD.blocked_report_bytes(
            error;
            config_hash=digest_or_na(options.config),
            feature_hash=digest_or_na(options.features),
            fold_hash=digest_or_na(options.folds),
            input_receipt_hash=digest_or_na(options.receipt),
            unary_review_hash=digest_or_na(SAD.UNARY_REVIEW_PATH),
            source_hash=source_or_na,
        )
        publication_error = try
            expected_source = occursin(r"^[0-9a-f]{64}$", source_or_na) ?
                              source_or_na : nothing
            SAD.write_report(
                options.out, blocked_bytes; expected_source_hash=expected_source)
            nothing
        catch caught
            caught
        end
        print(stderr, "BLOCKED: ")
        showerror(stderr, error)
        if publication_error !== nothing
            print(stderr, "; publication=")
            showerror(stderr, publication_error)
        end
        println(stderr)
        return 2
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StructuredDiagnosticsCLI.entrypoint(copy(ARGS)))
end
