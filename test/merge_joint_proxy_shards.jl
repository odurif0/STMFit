#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "joint_proxy", "output_merge.jl"))
using .JointProxyShardMerge

if abspath(PROGRAM_FILE) == @__FILE__
    run_cli(ARGS)
end
