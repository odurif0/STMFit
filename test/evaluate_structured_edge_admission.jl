#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "structured_assignment", "edge_admission.jl"))
using .StructuredEdgeAdmission

const SEA = StructuredEdgeAdmission
const VALUE_FLAGS = Set([
    "--root",
    "--features",
    "--candidate-config",
    "--model-config",
    "--universe-dir",
    "--edge-dir",
    "--forward-receipt",
    "--backward-receipt",
    "--out-dir",
])
const ALLOWED_FLAGS = union(VALUE_FLAGS, Set(["--help"]))
const REQUIRED_FLAGS = sort!(collect(VALUE_FLAGS))

function usage(io::IO=stdout)
    println(io, "Usage: julia --project=. test/evaluate_structured_edge_admission.jl \\")
    println(io, "  --root PATH --features PATH --candidate-config PATH --model-config PATH \\")
    println(io, "  --universe-dir PATH --edge-dir PATH \\")
    println(io, "  --forward-receipt PATH --backward-receipt PATH --out-dir PATH")
    return nothing
end

function parse_arguments(arguments::Vector{String})::Dict{String,String}
    SEA.StructuredEdgeFeatures.StructuredUniverse.InputBoundary.validate_cli_arguments(
        arguments;
        allowed_flags=ALLOWED_FLAGS,
        value_flags=VALUE_FLAGS,
    )
    values = Dict{String,String}()
    index = 1
    while index <= length(arguments)
        pair = split(arguments[index], '='; limit=2)
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
        length(values) == 1 ||
            throw(EdgeAdmissionError(:cli_error,
                                     "--help accepts no evaluation inputs"))
        return values
    end
    for flag in REQUIRED_FLAGS
        haskey(values, flag) ||
            throw(EdgeAdmissionError(:missing_argument,
                                     "required option is absent: $flag"))
    end
    return values
end

function _as_admission_error(error)
    error isa EdgeAdmissionError && return error
    if hasproperty(error, :code)
        return EdgeAdmissionError(Symbol(getproperty(error, :code)),
                                  sprint(showerror, error))
    end
    return EdgeAdmissionError(:evaluation_error, sprint(showerror, error))
end

function _safe_hashes(values::Dict{String,String})
    hashes = Dict{String,String}()
    haskey(values, "--root") || return hashes
    root = values["--root"]
    isdir(root) || return hashes
    for (name, flag) in (
        "input_sha256" => "--features",
        "candidate_config_sha256" => "--candidate-config",
        "model_config_sha256" => "--model-config",
    )
        haskey(values, flag) || continue
        supplied = values[flag]
        path = isabspath(supplied) ? normpath(supplied) : normpath(joinpath(root, supplied))
        isfile(path) && !islink(path) || continue
        hashes[name] = SEA._sha256_file(path)
    end
    for (name, flag) in (
        "universe_receipt_sha256" => "--universe-dir",
        "edge_receipt_sha256" => "--edge-dir",
    )
        haskey(values, flag) || continue
        supplied = values[flag]
        directory = isabspath(supplied) ? normpath(supplied) : normpath(joinpath(root, supplied))
        path = joinpath(directory, "receipt.toml")
        isfile(path) && !islink(path) || continue
        hashes[name] = SEA._sha256_file(path)
    end
    return hashes
end

function run_evaluation(values::Dict{String,String})
    data = load_admission_data(
        values["--root"];
        features=values["--features"],
        candidate_config=values["--candidate-config"],
        model_config=values["--model-config"],
        universe_dir=values["--universe-dir"],
        edge_dir=values["--edge-dir"],
        forward_receipt=values["--forward-receipt"],
        backward_receipt=values["--backward-receipt"],
    )
    report = evaluate_admission(data)
    hashes = publish_report(values["--root"], values["--out-dir"], report)
    return report, hashes
end

function main(arguments::Vector{String}=copy(ARGS))::Int
    values = parse_arguments(arguments)
    if haskey(values, "--help")
        usage()
        return 0
    end
    try
        report, hashes = run_evaluation(values)
        println("status=", report.status)
        println("result_sha256=", report.result_sha256)
        println("files=", length(hashes))
        println("receipt_sha256=", hashes["receipt.toml"])
        return 0
    catch error
        admission_error = _as_admission_error(error)
        try
            publish_blocker_receipt(
                values["--root"],
                values["--out-dir"],
                admission_error;
                hashes=_safe_hashes(values),
            )
        catch publication_error
            showerror(stderr, _as_admission_error(publication_error))
            println(stderr)
        end
        throw(admission_error)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(main())
    catch error
        showerror(stderr, _as_admission_error(error))
        println(stderr)
        exit(2)
    end
end
