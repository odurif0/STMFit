#!/usr/bin/env julia
# Saved-table arithmetic only. Requires Julia 1.13, one thread, explicit inputs.
module ReconstructedRepresentationCLI
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
using .ReconstructedRepresentationDiagnostics

const USAGE = """
Usage: julia --startup-file=no --threads=1 --project=. test/diagnose_reconstructed_representation.jl \\
  --patches-bwd9 PATH --patches-bwd17 PATH --patches-fwd17 PATH \\
  --descriptor PATH --predictor PATH --fisher PATH \\
  --production-config PATH --identity-atol FLOAT --descriptor-atol FLOAT --outdir NEW_DIR

All options are required. Inputs are the six complete saved native tables, not
predictions, labels, raw scans or reference data. All keys and invalid rows remain.
The output directory must not exist (including empty directories and symlinks).
Grid and zero-L1 settings come from the unchanged production config. Arithmetic
comparison tolerances are explicit, not fitted and not classification thresholds.
No clustering, Fisher refit, new real score, replacement feature or grading runs.
"""

function main(args=ARGS)
    args in (["--help"], ["-h"]) && (print(USAGE); return 0)
    required = ["patches-bwd9", "patches-bwd17", "patches-fwd17", "descriptor", "predictor", "fisher",
                "production-config", "identity-atol", "descriptor-atol", "outdir"]
    iseven(length(args)) || error("Expected option/value pairs; use --help")
    options = Dict{String,String}()
    for i in 1:2:length(args)
        startswith(args[i], "--") || error("Expected named option: $(args[i])")
        name, value = args[i][3:end], args[i+1]
        name in required || error("Unknown option: $(args[i])")
        haskey(options, name) && error("Duplicate option: $(args[i])")
        isempty(value) && error("Empty option: $(args[i])")
        options[name] = value
    end
    all(haskey(options, name) for name in required) || error("Missing required option; use --help")
    paths = Dict(name => options[name] for name in required[1:6])
    run_audit(paths; production_config=options["production-config"], outdir=options["outdir"],
              identity_atol=parse(Float64, options["identity-atol"]),
              descriptor_atol=parse(Float64, options["descriptor-atol"]))
    return 0
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(ReconstructedRepresentationCLI.main())
    catch err
        err isa InterruptException && rethrow()
        print(stderr, "Saved representation diagnostic: ")
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end
