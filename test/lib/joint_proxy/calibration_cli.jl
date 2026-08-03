include("proxy_registry.jl")
include("simulator.jl")
include("count_calibration.jl")
include("type_posterior.jl")
include("type_calibration.jl")

module JointProxyCalibrationCLI

using GaussianFit2D

include("calibration_source_files.jl")
include("candidate_views.jl")
include("calibration_cli_types.jl")
include("calibration_cli_core.jl")
include("calibration_cli_adapter.jl")

export CalibrationCliError, CliOptions, parse_cli, run_cli

end # module
