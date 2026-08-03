module JointProxyTypeCalibration

using Statistics

using Main.JointProxyRegistry: ProxySource, ProxyTemplate, ProxyEntry, ProxyEnsemble
using Main.JointProxyTypePosterior: TypePosteriorPriors, TypePosteriorResult,
    TypePosteriorLobeEvidence, infer_type_posterior

include("type_calibration_types.jl")
include("type_calibration_scoring.jl")
include("type_calibration_fitting.jl")

export TypeCalibrationError, TypeCalibrationCase, TypePosteriorRow, TypeCasePosterior,
       TypeReliabilityBin, TypeCalibrationDiagnostics, TypeCalibrationModel,
       score_type_case, predict_type, calibrate_type_temperature

end # module
