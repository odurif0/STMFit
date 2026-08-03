module JointProxySyntheticValidation

using Statistics
using Printf
using SHA
using STMSXMIO

isdefined(Main, :JointProxySimulator) || Base.include(Main, joinpath(@__DIR__, "simulator.jl"))
isdefined(Main, :JointProxyRegistry) || Base.include(Main, joinpath(@__DIR__, "proxy_registry.jl"))
isdefined(Main, :JointProxyCandidateViews) || Base.include(Main, joinpath(@__DIR__, "candidate_views.jl"))
isdefined(Main, :JointProxyCountCalibration) || Base.include(Main, joinpath(@__DIR__, "count_calibration.jl"))
isdefined(Main, :JointProxyTypePosterior) || Base.include(Main, joinpath(@__DIR__, "type_posterior.jl"))
isdefined(Main, :JointProxyTypeCalibration) || Base.include(Main, joinpath(@__DIR__, "type_calibration.jl"))
isdefined(Main, :JointProxyOutputValidation) || Base.include(Main, joinpath(@__DIR__, "output_validation.jl"))

using Main.JointProxySimulator: SyntheticCase
using Main.JointProxyRegistry: ProxyEnsemble, ProxyEntry, ProxyTemplate
using Main.JointProxyCandidateViews: ViewData, ViewRecalibration, LobePatch, CandidateNView,
    CandidateViewReport, SkippedCandidate, effective_view_factor
using Main.JointProxyTypePosterior: TypePosteriorLobeEvidence, infer_type_posterior

export swap_registry, oracle_geometry_adapter_name, oracle_geometry_type_result,
    physical_candidate_report, selected_candidate, type_result_for_candidate,
    selected_type_result, calibration_observation_mode, candidate_report_summary, write_calibration_toml, write_metrics_tsv

include(joinpath(@__DIR__, "synthetic_validation_oracle.jl"))
include(joinpath(@__DIR__, "synthetic_validation_count.jl"))
include(joinpath(@__DIR__, "synthetic_validation_report.jl"))

end # module
