"""
GaussianFit1D.jl — Multi-Gaussian fitting of STM/AFM line profiles.

Fits a sum of Gaussians (with physical constraints) to a 1D intensity profile,
sweeps multiple values of n_peaks, and selects the best model via BIC.

Uses the most performant Julia packages available:
- Optimization.jl + NLopt for global optimization
- LsqFit.jl for Levenberg-Marquardt (local refinement + covariance)
- Plots.jl for publication-quality visualization
"""
module GaussianFit1D

using DelimitedFiles
using Statistics
using Printf
using LinearAlgebra
using Optimization
using OptimizationNLopt
using LsqFit
using JLD2
using Base.Threads  # for parallel sweep
using STMFitCore: FWHM_TO_SIGMA, effective_spacing_min, sigma_from_fwhm, adjacent_kappa_max, kappa_penalty, ResidualDiagnostics, compute_residual_diagnostics

# Types and configuration (must come first)
include("types.jl")

# Core computation
include("core.jl")

# Plotting and CLI
include("plot.jl")
include("cli.jl")

# GUI separated into its own repository: STMMolecularFitGUI.jl

# ===========================================================================
# Exports
# ===========================================================================

# Types
export FitConfig, FitResult, FitMetrics, FitRunResult, DEFAULT_CONFIG,
       FWHM_TO_SIGMA, MGF_VERSION

# Config
export build_config, output_path, user_data_dir

# Data & Model
export load_data, multi_gaussian,
       _params_to_centers, _get_amplitude, _get_sigma

# Fitting
export build_bounds, fit_model, run_model_comparison, compute_metrics, predict_fit

# API
export run_fit, select_model_result, best_result, result_for_n,
       update_model_rankings

# Serialization & Export
export save_results, load_results, export_results

# CLI
export parse_args, config_from_args, main_cli

# Plotting
export plot_results

end # module
