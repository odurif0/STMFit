#!/usr/bin/env julia
# Saved prediction ablation only. No fitting, benchmark membership or truth input.
module SavedAssignmentComponents

using SHA, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: lobe_table, require_same_keys, write_table

const FIELDS = ["file", "lobe", "predicted", "confidence", "probability_1", "invalid_reason", "model"]
const REQUIRED = ["predicted", "confidence", "probability_1", "invalid_reason"]
const ENDPOINTS = ("kmeans_endpoint_common_validity", "gmm_endpoint_common_validity")
sha(path) = open(io -> bytes2hex(sha256(io)), path, "r")
number(value) = something(tryparse(Float64, value), NaN)

function available(row)
    row["predicted"] in ("0", "1", "?") || error("Invalid predicted label")
    p = number(row["probability_1"])
    confidence = number(row["confidence"])
    isfinite(confidence) && 0 <= confidence <= 1 || error("Invalid confidence field")
    isfinite(p) && !(0 <= p <= 1) && error("Vote frequency outside [0,1]")
    return row["predicted"] != "?" && isfinite(p)
end

function export_components(control_path, kmeans_path, gmm_path, destination)
    v"1.13" <= VERSION < v"1.14" || error("Julia 1.13 is required")
    paths = abspath.((control_path, kmeans_path, gmm_path))
    out = abspath(destination)
    (ispath(out) || islink(out)) && error("Output already exists: $out")
    isdir(dirname(out)) || error("Output parent must exist")
    hashes = sha.(paths)
    control = last(lobe_table(paths[1]; required=REQUIRED))
    components = [last(lobe_table(p; required=REQUIRED)) for p in paths[2:3]]
    for (name, rows) in zip(ENDPOINTS, components)
        require_same_keys(control, rows, name)
    end
    keys_sorted = sort(collect(keys(control)))
    eligible = Dict{Tuple{String,Int},Bool}()
    for key in keys_sorted
        c = control[key]
        flags = [available(rows[key]) for rows in components]
        joint = all(flags)
        available(c) == joint || error("Control/component availability mismatch: $key")
        (c["predicted"] != "?") == joint || error("Control must preserve unavailable keys: $key")
        if !joint
            !isfinite(number(c["probability_1"])) || error("Unavailable control has a finite vote")
        end
        eligible[key] = joint
    end
    # No membership filter. Every scientific key, including extra-N rows, survives.
    mkdir(out)
    cp(paths[1], joinpath(out, "control.tsv"); force=false)
    counts = Dict{String,Any}()
    for (name, rows) in zip(ENDPOINTS, components)
        exported = Dict{String,String}[]
        for key in keys_sorted
            chosen = eligible[key] ? rows[key] : control[key]
            r = Dict(field => chosen[field] for field in REQUIRED)
            r["file"] = key[1]; r["lobe"] = string(key[2]); r["model"] = name
            # Copy the original head decision, even if printed p rounds to 0.5.
            # Neither a new threshold nor a new class mapping is introduced.
            push!(exported, r)
        end
        write_table(joinpath(out, name * ".tsv"), FIELDS, exported)
        counts[name] = Dict("rows" => length(exported), "available" => count(values(eligible)),
                            "unavailable" => count(!, values(eligible)))
    end
    sha.(paths) == hashes || error("An input changed during export")
    sha(joinpath(out, "control.tsv")) == hashes[1] || error("Control copy changed")
    metadata = Dict{String,Any}(
        "schema" => "saved_assignment_component_endpoints_v1", "julia_version" => string(VERSION),
        "inputs" => Dict(name => Dict("path" => p, "sha256" => h) for (name,p,h) in zip(("control","kmeans","gmm"),paths,hashes)),
        "sources" => Dict("exporter" => sha(@__FILE__), "table_helper" => sha(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))),
        "outputs" => Dict(file => sha(joinpath(out,file)) for file in ("control.tsv", (name * ".tsv" for name in ENDPOINTS)...)),
        "counts" => counts,
        "policy" => "Keep every key; use exact saved head decisions on the common component-validity mask; otherwise retain control abstention and reason.",
        "limits" => "No refit, labels, benchmark filter, calibration, threshold change, confidence-as-posterior claim, or production promotion. Counts are the unchanged cached own-N lineage.")
    open(joinpath(out, "metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return metadata
end

function main(args=ARGS)
    if args == ["--help"]
        println("Usage: julia --startup-file=no --threads=1 --project=. test/export_saved_assignment_components.jl CONTROL KMEANS GMM NEW_OUT")
        println("Only saved label-free predictions. Preserves common availability and all keys; never loads benchmark membership/truth or fits a model.")
        return
    end
    length(args) == 4 || error("Expected CONTROL KMEANS GMM NEW_OUT; use --help")
    metadata = export_components(args...)
    println("Frozen saved-component endpoints exported; no recognition gain is claimed before external grading.")
    println(metadata["counts"])
end

end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && SavedAssignmentComponents.main()
