include("proxy_registry.jl")
include("candidate_views.jl")
include("count_calibration.jl")
include("type_posterior.jl")

module JointProxyInferenceCLI

include(joinpath(@__DIR__, "..", "script_utils.jl"))
using .ScriptUtils: _ensure_parent

using SHA
using Statistics
using Printf
using TOML
using STMSXMIO
using GaussianFit2D

using Main.JointProxyRegistry: load_registry, sha256_file
using Main.JointProxyCandidateViews: ViewData, ViewRecalibration, LobePatch,
    CandidateNView, CandidateViewReport, SkippedCandidate,
    build_view_data, extract_candidate_views, effective_view_factor
using Main.JointProxyCountCalibration: CountPosteriorRow, CountCasePosterior,
    CountCalibrationError, score_count_report
using Main.JointProxyTypePosterior: TypePosteriorPriors, TypePosteriorGlobalState,
    TypePosteriorLobeEvidence, TypePosteriorResult, infer_type_posterior

include("calibration_source_files.jl")
include("inference_cli_types.jl")
include("inference_cli_core.jl")
include("inference_cli_real_adapter.jl")

export InferenceCliError, CliOptions, InferenceArtifacts, parse_cli, run_cli

end # module
