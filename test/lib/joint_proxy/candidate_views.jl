module JointProxyCandidateViews

using LinearAlgebra
using Statistics

using GaussianFit2D
using GaussianFit2D: ChainModelResult, ChainSweepConfig

export ViewData, ViewRecalibration, LobePatch, CandidateNView,
       CandidateViewReport, SkippedCandidate,
       build_view_data, extract_candidate_views, effective_view_factor

include("candidate_views_types.jl")
include("candidate_views_patches.jl")
include("candidate_views_scoring.jl")

end # module
