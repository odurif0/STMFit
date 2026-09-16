#!/usr/bin/env julia
# External comparison only. Never called by the production runner.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using Printf, TOML

function compare_predictions(predictions, reference, outdir; benchmark_manifest="")
    pheader, pred = lobe_table(predictions; required=["predicted", "confidence"])
    _, ref = lobe_table(reference; required=["predicted", "confidence"])
    for table in (pred, ref), row in values(table)
        row["predicted"] in ("0", "1", "?") || error("Invalid assignment in comparison input")
    end
    (ispath(outdir) || islink(outdir)) && error("Comparison output already exists")
    mkpath(outdir)
    overlap = intersect(Set(keys(pred)), Set(keys(ref)))
    missing = setdiff(Set(keys(ref)), Set(keys(pred)))
    extra = setdiff(Set(keys(pred)), Set(keys(ref)))
    rows = Dict{String,String}[]
    agrees = 0
    deltas = Float64[]
    for key in sort(collect(union(Set(keys(pred)), Set(keys(ref)))))
        p, r = get(pred, key, nothing), get(ref, key, nothing)
        status = p === nothing ? "missing_prediction" : r === nothing ? "extra_prediction" :
                 p["predicted"] == r["predicted"] ? "same_assignment" : "different_assignment"
        delta = NaN
        if p !== nothing && r !== nothing
            p["predicted"] == r["predicted"] && (agrees += 1)
            pc = something(tryparse(Float64, p["confidence"]), NaN)
            rc = something(tryparse(Float64, r["confidence"]), NaN)
            delta = abs(pc - rc)
            isfinite(delta) && push!(deltas, delta)
        end
        push!(rows, Dict("file" => key[1], "lobe" => string(key[2]), "status" => status,
            "prediction" => p === nothing ? "NA" : p["predicted"],
            "reference_prediction" => r === nothing ? "NA" : r["predicted"],
            "confidence_abs_difference" => isfinite(delta) ? @sprintf("%.10g", delta) : "NA"))
    end
    write_table(joinpath(outdir, "lobe_comparison.tsv"),
        ["file", "lobe", "status", "prediction", "reference_prediction", "confidence_abs_difference"], rows)
    same_keys = isempty(missing) && isempty(extra)
    same_assignments = same_keys && agrees == length(ref)
    same_confidence = same_assignments && length(deltas) == length(ref) && all(iszero, deltas)
    summary = Dict("prediction_lobes" => string(length(pred)), "reference_lobes" => string(length(ref)),
        "overlap_lobes" => string(length(overlap)), "same_assignment" => string(agrees),
        "missing_lobes" => string(length(missing)), "extra_lobes" => string(length(extra)),
        "identical_key_set" => string(same_keys), "identical_assignments" => string(same_assignments),
        "identical_reported_confidences" => string(same_confidence),
        "max_confidence_abs_difference" => isempty(deltas) ? "NA" : @sprintf("%.10g", maximum(deltas)))
    write_table(joinpath(outdir, "summary.tsv"),
        ["prediction_lobes", "reference_lobes", "overlap_lobes", "same_assignment", "missing_lobes", "extra_lobes",
         "identical_key_set", "identical_assignments", "identical_reported_confidences", "max_confidence_abs_difference"], [summary])
    if !isempty(benchmark_manifest)
        # Membership filtering belongs here, outside all prediction code.
        manifest = TOML.parsefile(benchmark_manifest)
        files = Set(String.(keys(manifest["files"])))
        candidate_files = Set(first.(collect(keys(pred))))
        isempty(setdiff(files, candidate_files)) || error("Predictions miss benchmark files; no partial grade emitted")
        graded_rows = [pred[k] for k in sort(collect(keys(pred))) if first(k) in files]
        write_table(joinpath(outdir, "predictions_for_external_grade.tsv"), pheader, graded_rows)
    end
    println("Reference agreement: $agrees/$(length(overlap)) overlapping lobes; $(length(missing)) missing, $(length(extra)) extra")
    println("Identical key set/assignments/reported confidences: $same_keys / $same_assignments / $same_confidence")
    println("This is output agreement on these inputs, not proof of descriptor identity or unknown-chain accuracy.")
    return summary
end

function main(args=ARGS)
    if "--help" in args || "-h" in args
        println("Usage: julia --project=. test/compare_reconstructed_champion.jl --predictions PATH --reference PATH --outdir NEW_DIR [--benchmark-manifest PATH]")
        println("External report only. Optional manifest selects rows for a separate existing grader; no fitting or parameter changes.")
        return
    end
    iseven(length(args)) || error("Expected --option VALUE pairs")
    opts = Dict{String,String}()
    for i in 1:2:length(args)
        args[i] in ("--predictions", "--reference", "--outdir", "--benchmark-manifest") || error("Unknown option: $(args[i])")
        haskey(opts, args[i]) && error("Repeated option: $(args[i])")
        opts[args[i]] = args[i+1]
    end
    all(k -> haskey(opts, k), ("--predictions", "--reference", "--outdir")) || error("Missing required options")
    compare_predictions(opts["--predictions"], opts["--reference"], opts["--outdir"];
                        benchmark_manifest=get(opts, "--benchmark-manifest", ""))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
