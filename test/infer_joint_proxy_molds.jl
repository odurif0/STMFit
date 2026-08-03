#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "joint_proxy", "inference_cli.jl"))
using .JointProxyInferenceCLI

if abspath(PROGRAM_FILE) == @__FILE__
    run_cli(ARGS)
end
