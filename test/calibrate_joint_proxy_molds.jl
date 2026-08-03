#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "joint_proxy", "calibration_cli.jl"))
using .JointProxyCalibrationCLI

function _main()
    try
        run_cli(ARGS)
    catch err
        println(stderr, err)
        exit(1)
    end
end

abspath(PROGRAM_FILE) == abspath(joinpath(@__DIR__, "calibrate_joint_proxy_molds.jl")) && _main()
