module JointProxyTypePosterior

using LinearAlgebra
using Statistics

using Main.JointProxyRegistry: ProxySource, ProxyTemplate, ProxyEntry, ProxyEnsemble

include("type_posterior_types.jl")
include("type_posterior_scoring.jl")
include("type_posterior_inference.jl")

export TypePosteriorPriors, TypePosteriorGlobalState, TypePosteriorLobeEvidence, TypePosteriorResult,
       score_ncc_cost, effective_view_factor, binary_chain_forward_backward, binary_chain_viterbi,
       infer_type_posterior

end # module
