# Composition-free hierarchical equal-prior unit-assignment model (T4).
#
# This is the entry-point module. The implementation is split across focused
# files under `test/lib/hierarchical/`, each below 250 pure LOC:
#
#   firewall.jl    - feature-name and CLI-flag firewall
#   io_helpers.jl  - self-contained TSV/path helpers
#   loading.jl     - LobeRecord, deterministic (file,lobe)-keyed loading
#   nuisance.jl    - per-scan median+MAD profiling (no class assignments)
#   emission_math.jl - diagonal Gaussian and log-likelihood helpers
#   emissions.jl     - shared 2-component EM (fixed 0.5/0.5 priors, covariance
#                      floors, monotone checks, deterministic multi-start) plus
#                      one-component companion
#   views.jl       - per-view fit, physical higher-amplitude mapping,
#                    in-sample one-component identifiability check,
#                    equal-weight view log-likelihood combination
#   pipeline.jl    - prediction writer (existing schema + model/version/
#                    provenance) and the top-level pipeline used by the CLI
#
# SCIENTIFIC CONSTRAINTS. Labels, expected N/NKNNKN, lobe index/order,
# transitions, run lengths, per-chain class counts, top-k rules, occupancy
# priors, and benchmark/truth columns may not enter fitting, features,
# selection, confidence, abstention, or calibration. Equal priors are
# immutable. The feature-name firewall and CLI flag firewall enforce this at
# the parse boundary.
#
# This library is consumed by test/build_hierarchical_unit_predictions.jl.

module HierarchicalUnitAssignment

using Printf
using Random
using SHA
using Statistics

# Immutable model identity.
const MODEL_NAME = "hierarchical_equalprior"
const MODEL_VERSION = 1
const EQUAL_PRIOR_WEIGHTS = [0.5, 0.5]

# Include focused sub-files (each < 250 pure LOC). All code lands in this
# module's scope, so cross-references need no qualification.
const _DIR = joinpath(@__DIR__, "hierarchical")
include(joinpath(_DIR, "io_helpers.jl"))
include(joinpath(_DIR, "firewall.jl"))
include(joinpath(_DIR, "loading.jl"))
include(joinpath(_DIR, "nuisance.jl"))
include(joinpath(_DIR, "emission_math.jl"))
include(joinpath(_DIR, "emissions.jl"))
include(joinpath(_DIR, "views.jl"))
include(joinpath(_DIR, "pipeline.jl"))

export MODEL_NAME, MODEL_VERSION, EQUAL_PRIOR_WEIGHTS,
       LobeRecord, ScanProfile, TwoComponentFit, OneComponentFit, ViewFit,
       is_forbidden_feature, is_forbidden_flag,
       load_records, feature_matrix, profile_scan, normalize_per_scan,
       fit_em_two_component, fit_one_component, responsibilities,
       log_likelihood_under,
       physical_high_amplitude_cluster,
       fit_view, combine_views,
       write_predictions, run_pipeline,
       cov_floor_default, default_n_starts

end # module HierarchicalUnitAssignment
