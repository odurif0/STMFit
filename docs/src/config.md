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
selection_policy = "support_midpoint_hybrid"
```

The chitosan default batch selection policy is now `support_midpoint_hybrid`,
configured in the TOML calibration. It first computes the integrated robust
overfit guard, then can move the guarded result by at most one lobe toward the
midpoint of the measured 2D support's physical feasible-N interval. The raw
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

### Sections

A calibration TOML has three sections; only `[model]` is mandatory (the others
have built-in defaults):

- **`[model]`** — physical calibration: lobe width (`sigma_parallel_*`), repeat
  spacing, ROI/support geometry, optimizer budget, overlap floor. These are the
  values to re-derive for a new molecule.
- **`[selection]`** — model-selection thresholds (label-free). Keep the defaults
  unless a sensitivity check (see `test/sensitivity_thresholds.jl`) shows your
  molecule's GCV curve needs a different ambiguity band:
  - `gcv_ambiguity_rel_threshold` (default `0.05`): relative GCV gap below which
    two N are considered indistinguishable. Feeds the up-when-ambiguous branch
    and the `ambiguous_eff` summary column. Overridable per-run with
    `--gcv-ambiguity-rel-threshold`.
  - `robust_guard_nu` (default `8.0`): Student-t degrees of freedom for the
    robust-AICc guard. Overridable per-run with `--robust-guard-nu`.
  - `support_midpoint_up_gcv_rel_threshold` (default `0.30`): relative effective
    GCV gap allowed for one-step support-midpoint upshifts under
    `support_midpoint_hybrid`.
- **`[preprocessing]`** — SXM channel name/direction, stride, flatten, smoothing.

### Structured Todo 9 diagnostic policy

The label-free structured diagnostics use an immutable policy block in
`config/unit_assignment_structured_model.toml`:

```toml
[selection.diagnostics.policy]
schema = "structured_diagnostics_policy_v1"
residual_model = "fixed_nu8_diagonal_student_t_two"
residual_formula = "x_minus_posterior_weighted_component_mean"
residual_fallback = "matched_one_component_student_t_mean"
residual_failure = "BLOCKED"
icc_estimator = "ICC_1_1_oneway_random_unbalanced"
icc_date_centering = "row_weighted_per_date_per_feature"
icc_feature_pooling = "equal_feature_summed_squares"
icc_unbalanced_group_size = "n0=(N-sum(n_s^2)/N)/(K-1)"
view_score = "mean_lobes(log_predictive_density/dimension)"
view_contrast = "mean(bwd_com,bwd_diag)-mean(base,split)"
channel_drop_features = ["bwd_neg_com_t", "bwd_neg_diag45"]
channel_drop_refit = "fresh_inner_training_base_local"
channel_drop_statistic = "mean(full_bwd_scores)-refit_base_score"
channel_drop_same_sign = "strict_nonzero_each_inner_date"
quantile_method = "Hyndman_Fan_type_7_explicit"
permutation_tail = "upper_inclusive_plus_one"
view_tail = "two_sided_zero_inclusive_plus_one"
threshold_equality = "SKIPPED"
holm_order = ["mfa_q1", "scan_effects", "view_asymmetry"]
holm_reject_equality = true
```

Missing or differing policy keys are `BLOCKED`. Diagnostics remain follow-up
only and are recomputed inside the outer/inner training partitions. Channel
dropout removes both backward descriptors, freshly refits `base_local` on the
same inner-training rows, and rescores the untouched held-out rows; it is not
the original backward-versus-non-backward view contrast.

### Calibrating a new molecule

Start from the annotated template:

```bash
cp config/template.toml config/my_molecule.toml
```

The template comments explain how to derive each value from a few representative
scans (FWHM -> sigma, observed pitch -> spacing, etc.). Only `[model]` values
must change for a molecule on the same STM. The `[selection]` defaults are the
current chitosan defaults; re-validate them before promoting the
support-midpoint hybrid to a new molecule.

To exclude non-target files (noise, test images, other molecules), pass an
exclusion list instead of relying on hard-coded defaults:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/my_molecule.toml \
  --exclude-from results/my_molecule_exclude.txt
```

The exclusion file is one `.sxm` filename per line (`#` comments allowed).

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
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml
```

> `--skip-1d` is now the default (the 1D fit is diagnostic-only and never
> affects `N_selected`). Add `--no-skip-1d` to compute the `N_1D` columns.

The older standard/rescue/aggressive passes remain useful for comparison and
audit, but should not be treated as ground truth when they disagree with the
generic adaptive workflow.

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_rescue \
  --tsv results/10_20mer_analysis_rescue/triage_unused.tsv \
  --config config/chitosan_10_20mer_rescue.toml \

JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_rescue_aggressive \
  --tsv results/10_20mer_analysis_rescue_aggressive/triage_unused.tsv \
  --config config/chitosan_10_20mer_rescue_aggressive.toml \
```

Optional legacy guard-audit pass for comparison:

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_guard_audit \
  --tsv results/10_20mer_analysis_guard_audit/triage_unused.tsv \
  --config config/chitosan_10_20mer.toml \
  --selection-policy gcv_with_robust_aicc_guard \
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
`N_selected` such as the chitosan default support-midpoint hybrid.

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

## Structured evaluator-v1 policy (correction3 pending review)

`config/unit_assignment_structured_evaluator.toml` is the correction3
policy/evidence-only candidate. It is authoritative only after parent
acceptance, a fresh independent Oracle PASS, and reviewer-owned `GateClosure`;
until then Todo 13 remains
blocked. The replaced predecessor config is preserved as an exact historical
preimage in correction2 evidence. The live config binds canonical
repository-relative paths, exact bytes, source-bundle members, Julia
1.12/1.12.6, the checked Todo 12/Todo 13 markers, and immutable T8, T11, and
T12 authorities. Symlinks, hardlinks, path escapes, stale hashes, duplicate or
substituted members, and runtime mismatches fail closed before scoring work.

The frozen implementation semantics are explicit: `j_mivc=log(0.5)+log
f_mvc(x_iv)`, `a_miv=logsumexp_c(j_mivc)`, `pi_mivc=exp(j-a)`, clipping only
state 1 to `[1e-12,1-1e-12]`, shared `A_i`, and
`eta_mic=(1/|A_i|)sum_v(a_miv+log(pi_clipped_mivc))`. Then
`U_mi=logsumexp_c(eta_mic)` and `q_mic=exp(eta-U)`. Structural absence is
omitted for both models; partial views use `1/|A_i|`; empty `A_i` is
`U=0,q=(0.5,0.5),?`, and an eligible incident edge is `BLOCKED`.

The exact status event/reason/consequence tables, edge-null lifecycle, graph
reference selection and per-scan `logZ` rule are in the config. Graph output
uses T12 marginals, never Viterbi or report-wide totals. Scoring fixes
`E_s`, `M_s=nodes_s+|E_s|`,
`L_C1=sum_i(U_C1_i)+sum_e(N_e)`, and
`L_meta=sum_i(U_selected_i)+sum_e(N_e)+I_graph*logZ_T12_s`; T12 is only a
relative graph normalizer. The gate uses one Mersenne Twister per seed 0–499,
paired whole-scan resampling, Type 7 interpolation
`0.525*x_(13)+0.475*x_(14)`, and exhaustive inclusive sign masks with no `+1`.

The validator's results are policy/static/synthetic evidence only. They are not
physical calibration, application processing, benchmark validation, or a
permission to use labels; any external labels or grading are post-gate only.
The correction2 validator snapshots authorities descriptor-relatively,
revalidates identities and bytes before return, and reports structured status,
reason, and truthful work counters. Exact selected-model/unary-fit/T11/T12
reference keys and cardinalities are required; metric populations are pooled
over their explicitly frozen node/pair scopes.

### Correction3 integration addendum

Correction3 integration is a worker/static/synthetic evidence phase only. It
does not create, import, execute, or authorize either Todo13 evaluator product.
The parent-reproduced policy and authority lanes are bound losslessly: 680
policy/semantic mutations plus 163 authority mutations produce 843 projected
rows, but the integrator does **not** independently reimplement those 843
mutations. The integration validator snapshots the live configuration and all
consumed lane/predecessor inputs descriptor-relatively, validates the
source-authored 202-key policy, binds 20 roles, 32 bundle members, 42
claim/review checks, 13 structural checks, and revalidates every snapshot after
semantic fixtures before its final Julia runtime check. Its evidence is not a
benchmark or 10–20mer application claim; labels remain external reporting only.

### Correction4 no-replace provenance successor

Correction3 is technically green but remains historically blocked because its
canonical paths were replaced while correcting the combined-row projection.
Correction4 is a provenance-only successor: it freshly regenerates the same six
canonical static/synthetic files twice from the immutable correction3 validator,
then publishes the run-1 bytes once with exclusive no-replace creation. It does
not change science, configuration, thresholds, T8/T11/T12, GCV, `n_eff`, labels,
benchmarks, application claims, or Todo behavior. The live evaluator config
remains a candidate; Todo13 products remain absent and blocked pending parent
acceptance, a fresh Oracle, and reviewer-owned GateClosure.

### Correction5 terminalization caveat

Correction4's six canonical `O_EXCL` bytes remain valid and are referenced
without republishing. Its `DoneClaim` is non-authoritative because it carries
the predecessor-hash terminalization defect and its finalizer failed before the
final Boulder closure. Correction5 preserves those failure artifacts and cleans
only the captured staging residue. No science, configuration, threshold, GCV,
`n_eff`, T8/T11/T12, label, benchmark, application, or Todo behavior changes;
Todo13 remains blocked pending parent acceptance, a fresh Oracle PASS, and
reviewer-owned GateClosure.

### Correction6 parent-owned atomic publication boundary

Correction5 cleanup remains valid, but Correction5 is non-authoritative because
its evidence closure omitted the cleanup receipt and its terminal replay omitted
all six required per-path canonical bindings. Correction6 references the
existing Correction4 canonical bytes without republishing them. Its complete
19-file payload is sealed in a hidden same-parent stage, and publication is
established only when the parent supplies and accepts the external checkpoint,
the descriptor-bound publisher performs exactly one
`renameat2(RENAME_NOREPLACE)`, and the parent records the external publication
receipt. The publisher reports both the transient basename-manifest namespace
and the final repository-relative manifest namespace; it writes neither
manifest nor receipt into the repository.

Correction6 never authorizes Todo13 itself. Todo13 remains blocked pending
parent acceptance of that receipt, a fresh independent Oracle PASS, and
reviewer-owned `GateClosure`. This is administrative provenance only: it
changes no policy, configuration, calibration, threshold, GCV, `n_eff`, T8,
T11, T12, label, benchmark, application, or Todo behavior.
