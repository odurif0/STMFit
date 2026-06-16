# Configuration Reference

## User-facing TOML calibration

Batch runs are configured from TOML files. The default calibration is
`config/chitosan.toml` and `test/batch_full.jl` accepts an override:

```bash
julia --project=. test/batch_full.jl --config config/my_system.toml
```

The chitosan calibration currently uses a noise-only support rule:

```toml
[model]
fit_width_nm = 0.16
support_noise_k = 2.5
support_padding_nm = 0.25
kappa_max = 10.0
selection_criterion = "gcv"
cv_method = "gcv"
selection_policy = "gcv_with_robust_aicc_guard"
```

The chitosan default batch selection policy is the integrated robust overfit
guard, configured in the TOML calibration.  It writes `N_selected` as the
guarded primary result while preserving `N_eff` for comparison.  The raw
GCV/effective baseline remains available as an explicit command-line override:

```bash
julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy gcv
```

For manually labelled folders, pass the benchmark manifest to suppress best-plot
generation for non-chitosan/excluded files without using labels in model
selection:

```bash
julia --project=. test/batch_full.jl 28 \
  --data-dir /home/durif/Rebecca/data/data/20240818_LHe_Cu100 \
  --outdir results/best_plots_20240818 \
  --config config/chitosan.toml \
  --plot-manifest benchmarks/chitosan_manual_20240814_20240818.toml
```

By default this skips plots for `quality = "excluded"`; override with
`--skip-plot-quality excluded,ambiguous` if ambiguous files should be hidden too.

## Generic adaptive-support workflow

The benchmark-validated generic workflow is `adaptive_support_rescue`: standard
support first, objective support-rescue only if the support appears truncated,
then the same robust-AICc down-only guard on the active support.  Benchmark
labels are used only for external grading, never inside fitting or selection.

Short-chain benchmark-style pass:

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 39 \
  --data-dir /home/durif/Rebecca/data/data/20240817_LHe_Cu100 \
  --outdir results/best_plots_240817_adaptive_support_rescue \
  --tsv results/best_plots_240817_adaptive_support_rescue/primary_files.tsv \
  --config config/chitosan_adaptive_support_rescue.toml
```

For curated long-chain 10–20mer analyses, use the same workflow with only the
allowed N range extended to `n_max = 24`.

10–20mer adaptive pass:

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_adaptive_support_rescue \
  --tsv results/10_20mer_analysis_adaptive_support_rescue/triage_unused.tsv \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --skip-1d
```

The older standard/rescue/aggressive passes remain useful for comparison and
audit, but should not be treated as ground truth when they disagree with the
generic adaptive workflow.

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_rescue \
  --tsv results/10_20mer_analysis_rescue/triage_unused.tsv \
  --config config/chitosan_10_20mer_rescue.toml \
  --skip-1d

JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_rescue_aggressive \
  --tsv results/10_20mer_analysis_rescue_aggressive/triage_unused.tsv \
  --config config/chitosan_10_20mer_rescue_aggressive.toml \
  --skip-1d
```

Optional legacy guard-audit pass for comparison:

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_guard_audit \
  --tsv results/10_20mer_analysis_guard_audit/triage_unused.tsv \
  --config config/chitosan_10_20mer.toml \
  --selection-policy gcv_with_robust_aicc_guard \
  --skip-1d
```

For legacy comparisons, build the consolidated table and annotated plots:

```bash
python3 test/finalize_10_20mer_results.py \
  --output-dir results/10_20mer_analysis_final
```

Outputs:

- `results/10_20mer_analysis_final/final_results.tsv`
- `results/10_20mer_analysis_final/final_results.md`
- `results/10_20mer_analysis_final/plots/*.png`

Each legacy final plot keeps the original fit panels intact and adds a footer
showing `N final`, selected pass, confidence, standard/rescue/aggressive GCV
results, and whether the robust guard would change the final result.  `review`
is a QC confidence label for support sensitivity or diagnostic disagreement; it
is not an exclusion flag and does not change `N_final`.

A more diagnostic spatial blocked-CV selector is also available:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy spatial_blocked_cv \
  --cv-folds 3
```

It is more directly objectivable as a predictive-risk estimate, but current
smoke tests show it is not stable enough for default use; see
[Model Selection](selection.md#experimental-spatial-blocked-cv-selector).

A cheap support-sensitivity diagnostic can be enabled with:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy support_marginalized_gcv
```

or with a one-lobe capped overfit guard:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy support_marginalized_gcv_guard
```

It rescores fitted candidates across a fixed support-padding grid and selects by
median relative GCV regret.  It is useful for support ambiguity audits, but is
not recommended as the default selector; see
[Model Selection](selection.md#experimental-support-marginalized-gcv-selector).

A file-adaptive slope-heuristic MDL selector can also be enabled:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy slope_heuristic_mdl
```

It estimates the model-complexity penalty from the file's own
contrast–dimension curve.  This is statistically principled, but current smoke
tests make it diagnostic rather than default; see
[Model Selection](selection.md#experimental-slope-heuristic-mdl-selector).

A support-perturbation stability selector can be enabled with:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy stability_selection
```

It chooses the `N` that is most often within 1% of the best GCV across the fixed
support-padding grid.  This is useful for stability audits, but current smoke
tests make it diagnostic rather than default; see
[Model Selection](selection.md#experimental-stability-selection-selector).

A local lobe-resolvability guard can be enabled with:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy local_lobe_evidence
```

It checks whether adjacent fitted lobes are locally separated by a valley and is
strictly down-only.  Current smoke tests make it a separability diagnostic rather
than a recommended primary selector; see
[Model Selection](selection.md#experimental-local-lobe-evidence-guard).

An approximate Laplace-evidence selector and safer guard can be enabled with:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy laplace_evidence_guard
```

It scores fitted candidates with a local Gauss–Newton/Laplace evidence
approximation.  The direct selector is currently too parsimonious on some files;
the guard only permits one-lobe downshifts from `N_eff`.  See
[Model Selection](selection.md#experimental-laplace-evidence-selector).

A fwd/bwd direction-consensus selector can be enabled with:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy fwd_bwd_consensus
```

It exploits forward/backward scan replication by evaluating the fused model on
separate fwd and bwd channels.  Current smoke tests make it a replication
diagnostic rather than a recommended primary selector; see
[Model Selection](selection.md#experimental-fwdbwd-direction-consensus-selector).

The former contrast-fraction support threshold has been removed from the
program. Support is defined from the axial profile as
`baseline + support_noise_k * noise`, then expanded by `support_padding_nm`.
This avoids coupling support detection to the brightest lobe or to occasional
artefacts.

## ChainSweepConfig (GaussianFit2D)

Full configuration for the 2D chain model sweep.

```julia
GaussianFit2D.ChainSweepConfig(
    # ── Sweep range ──
    n_min              = 2,       # Minimum N (safety bound)
    n_max              = 14,      # Maximum N (safety bound)
    intelligent_sweep  = true,    # Adaptive range from support length
    early_stop_patience = 3,      # Consecutive BIC increases before stop
    early_stop_dbic    = 100.0,   # BIC increase threshold for early stop

    # ── Physical constraints ──
    spacing_min_nm     = 0.35,    # Minimum inter-lobe spacing (nm)
    spacing_max_nm     = 0.75,    # Maximum inter-lobe spacing (nm)
    max_overlap        = 0.60,    # Maximum lobe overlap fraction
    fit_width_nm       = 0.45,    # Struct default; chitosan.toml overrides to 0.16

    # ── Sigma bounds ──
    sigma_parallel_min_nm = 0.191, # Min axial sigma (FWHM 0.45 nm)
    sigma_parallel_max_nm = 0.509, # Max axial sigma (FWHM 1.20 nm)
    sigma_perp_min_nm   = 0.10,   # Min perpendicular sigma
    sigma_perp_max_nm   = 0.55,   # Max perpendicular sigma

    # ── Model variants ──
    chain_circular_sigmas = false, # σ∥=σ⟂ per lobe (simpler, more robust)
    chain_tilted_baseline = true,  # Linear tilt in 2D background

    # ── Optimization ──
    global_maxtime     = 10.0,    # NLopt timeout per N
    global_maxiter     = 10000,   # NLopt max iterations per N
    global_tol         = 1e-5,    # NLopt tolerance
    max_iter           = 300,     # LsqFit max iterations
    multistart         = 1,       # Number of random starts

    # ── Support detection ──
    support_noise_k    = 2.5,
    support_padding_nm = 0.20,    # Struct default; chitosan.toml uses 0.25
    support_min_length_nm = 1.0,
    support_baseline_quantile = 0.10,

    # ── Penalties ──
    kappa_max          = 10.0,    # Condition number penalty threshold
    kappa_weight       = 1.0,     # Condition number penalty strength
    peak_profile       = :gaussian, # :gaussian (2D: only :gaussian supported)
    min_amplitude_fraction = 0.3, # Min lobe amplitude (fraction of max data)

    # ── Cross-validation ──
    cv_folds           = 5,       # Number of CV folds (kfold only)
    cv_method          = "gcv",   # "gcv" (analytical, free) | "kfold" (refit per fold)
    student_nu         = 4.0,     # Student-t degrees of freedom
    residual_peak_snr_threshold = 3.5,

    # ── Selection ──
    selection_criterion = "gcv",  # "gcv" | "bic" | "aicc" | "cv"
)
```

`selection_criterion` controls the score used inside each model sweep.  The
batch-level `selection_policy` / `--selection-policy` is separate: it controls
whether the final reported primary result is the standard `N_eff` or a guarded
`N_selected` such as the chitosan default robust-AICc guard.

Experimental support rescue is available via
`config/chitosan_adaptive_support_rescue.toml` or
`--selection-policy adaptive_support_rescue`.  It runs the standard support
first and only tries a permissive support pass when the selected `N_eff` sits at
the objective support-feasibility ceiling.  Rescue acceptance is label-free and
requires a larger support, higher selected `N`, and circ/ell coherence.  The
robust-AICc guard is then applied down-only on the active support.
`adaptive_robust_guard_max_drop` is available only as a non-default diagnostic
cap on automatic robust-AICc downshifts; it is not used by the common
benchmark-aligned workflow.

## PatternConfig (GaussianFit2D)

Image preprocessing and blob detection configuration.

```julia
GaussianFit2D.PatternConfig(
    filepath   = "",          # SXM file path
    channel    = "Z",         # Channel name
    direction  = "fwd",       # Scan direction
    stride     = 2,           # Struct default; chitosan.toml uses 1
    flatten    = "plane+rows",# Background flattening
    smooth_radius_px = 1,     # Preprocessing smoothing
    threshold_sigma = 2.5,    # Blob detection threshold
    min_distance_px = 10,     # Minimum blob separation
    fusion     = true,        # Fuse Z fwd+bwd channels
    fuse_z_bwd = true,        # Fuse Z forward+backward scans
)
```

## FitSlideConfig (STMMolecularFit)

1D slide profile fitting configuration.

```julia
STMMolecularFit.FitSlideConfig(
    min_spacing    = 0.35,    # Minimum peak spacing (nm)
    max_spacing    = 0.75,    # Maximum peak spacing (nm)
    fwhm_min       = 0.45,    # Minimum FWHM (nm)
    fwhm_max       = 1.20,    # Maximum FWHM (nm)
    max_overlap    = 0.60,    # Maximum peak overlap
    kappa_max      = 10.0,    # Condition number threshold
    peak_profile   = :gaussian,  # :gaussian | :lorentzian | :pseudo_voigt
    amplitude_min_fraction = 0.3,
    global_maxtime = 8.0,     # NLopt timeout (s)
    global_maxiter = 5000,    # NLopt max iterations
)
```
