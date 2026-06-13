module GaussianFit2D

using ArgParse
using LinearAlgebra
using Random
using LsqFit
using Optimization
using OptimizationNLopt
using Plots
using Printf
using Statistics
using STMFitCore: effective_spacing_min, endpoint_overrun, overlap_condition_number, kappa_penalty, adjacent_kappa_max, ResidualDiagnostics, compute_residual_diagnostics

export SXMImage, SXMChannel, PatternConfig, PatternFitResult, MolecularFeature,
       MolecularChain, ImageArtifactDiagnostics
export read_sxm, channel_names, get_channel, preprocess_channel, detect_blobs,
       fit_molecular_pattern, detect_molecular_chains, main_cli
export chain_gaussian_sweep, chain_direct_fit, fit_chain_consensus, fit_chain_batch
export compute_image_artifact_diagnostics

const EPS = 1e-12

include("types.jl")
include("core.jl")
include("plot.jl")
include("cli.jl")

end # module
