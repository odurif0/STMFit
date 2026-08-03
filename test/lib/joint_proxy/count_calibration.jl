module JointProxyCountCalibration

using Statistics

include("count_calibration_types.jl")
include("count_calibration_scoring.jl")
include("count_calibration_fitting.jl")

export CountCalibrationError, CountCalibrationCase, CountPosteriorRow,
       CountCasePosterior, CountReliabilityBin, CountCalibrationDiagnostics,
       CountCalibrationModel, score_count_report, predict_count,
       calibrate_count_temperature

end # module
