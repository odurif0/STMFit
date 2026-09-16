#!/usr/bin/env julia
# Explicit label-free descriptor reconstruction; no benchmark inputs.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment

function main(args=ARGS)
    if "--help" in args || "-h" in args
        println("Usage: julia --project=. test/build_reconstructed_descriptor.jl --features PATH --patches PATH --config PATH --out PATH")
        println("Adds patch_u_asym_reconstructed = sum(sign(u)*p)/sum(abs(p)) on normalized backward residual 9x9 patches.")
        println("Pixel order: u outer, t inner. Invalid/zero-signal patches remain NA with descriptor_reason. No benchmark inputs.")
        return
    end
    opts = Dict{String,String}()
    allowed = Set(["--features", "--patches", "--config", "--out"])
    iseven(length(args)) || error("Expected --option VALUE pairs")
    for i in 1:2:length(args)
        args[i] in allowed || error("Unknown or forbidden option: $(args[i])")
        haskey(opts, args[i]) && error("Repeated option: $(args[i])")
        opts[args[i]] = args[i+1]
    end
    all(k -> haskey(opts, k), allowed) || error("--features, --patches, --config and --out are required")
    augment_descriptor(opts["--features"], opts["--patches"], opts["--out"], opts["--config"])
end

abspath(PROGRAM_FILE) == @__FILE__ && main()
