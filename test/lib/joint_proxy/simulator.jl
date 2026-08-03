module JointProxySimulator

using Random
using SHA
using Printf
using STMSXMIO: SXMImage, SXMChannel

include("simulator_types.jl")
include("simulator_proxy.jl")
include("simulator_nuisance.jl")
include("simulator_generation.jl")

export ProxySite, ProxyEnsemble, SimulatorConfig, ViewNuisanceTruth, SyntheticTruth,
       SyntheticCase, ControlType,
       CONTROL_NORMAL, CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS,
       CONTROL_SWAPPED_TYPES, CONTROL_MISSING_BWD, CONTROL_CORRUPTED_VIEW,
       ALL_CONTROLS,
       default_proxy_ensemble, load_proxy_ensemble,
       validate_config, replay_case, generate_case, generate_batch,
       case_checksum, image_data_bytes

end # module JointProxySimulator
