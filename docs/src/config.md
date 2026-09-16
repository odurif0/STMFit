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

> **Historical/current status note (2026-09-07).** The correction3--6 material
> below, including runtime-v3, records historical Julia 1.12.6-era authorities
> and is not the current runtime status. Historical evaluator configs v0, v1,
> v2, and v3 remain byte-identical records with their historical SHA-256 values.
> A live first-failure check of each historical config stops at stale T12 claim
> bytes; no old T12 bytes were reconstructed. The additive runtime-v4 Gate5
> candidate is documented separately at the end of this section.

## Historical structured evaluator-v1 policy (correction3)

`config/unit_assignment_structured_evaluator.toml` was the correction3
policy/evidence-only candidate. It was authoritative only after parent
acceptance, a fresh independent Oracle PASS, and reviewer-owned `GateClosure`;
until then Todo 13 remained blocked. The replaced predecessor config is
preserved as an exact historical preimage in correction2 evidence. That
historical config bound canonical repository-relative paths, exact bytes,
source-bundle members, Julia 1.12/1.12.6, the checked Todo 12/Todo 13
markers, and immutable T8, T11, and T12 authorities. Symlinks, hardlinks, path
escapes, stale hashes, duplicate or substituted members, and runtime mismatches
failed closed before scoring work.

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

### Gate5 runtime-v4 administrative candidate (2026-09-07)

`config/unit_assignment_structured_evaluator_runtime_v4.toml` is the additive
Julia 1.12.7 Gate5 implementation candidate. Its final SHA-256 is
`abe22ed5047c898f594067d4a54cbe8fcc99c15fef17b429fa366f64c255547b`.
The root `Project.toml` requires `julia = "=1.12.7"`. On Julia 1.12.7,
`Pkg.resolve()` made no change to `Manifest.toml`, and
`Pkg.Types.workspace_resolve_hash(ctx.env)` equals the recorded Manifest
`project_hash` `af56fc4d010e652f014ae3cec6d4e804d4722f8d`.

The scientific block from `[unary]` through EOF is byte-identical to certified
runtime-v3. No scientific, physical, model, calibration, threshold, GCV, or
`n_eff` parameter changed. Historical config paths and SHA-256 values remain
historical records; they are not rewritten to the runtime-v4 digest.

#### Administrative tables and fields

The fields before `[unary]` are administrative, not scientific model inputs:

| TOML table | Fields |
|---|---|
| `[evaluator]` | `schema`, `version`, `runtime_generation`, `allow_extra_keys`, `status_precedence`, `authority_state`, `authorizes_todo13` |
| `[prerequisite]` | `scope`, `julia`, `execution_used`, `real_data_used`, `labels_used`, `expected_N_used`, `NKNNKN_used`, `class_counts_used`, `composition_prior_used`, `equivalence_used`, `benchmark_or_grader_used` |
| `[history]` | `original_config_path`, `original_config_sha256`, `runtime_v1_config_path`, `runtime_v1_config_sha256`, `runtime_v2_config_path`, `runtime_v2_config_sha256`, `runtime_v3_config_path`, `runtime_v3_config_sha256`, `scope`, `runtime_equivalence` |
| `[runtime]` | `schema`, `version`, `julia`, `executable`, `executable_sha256`, `sysimage`, `sysimage_sha256`, `cross_runtime_compare`, `julia_1_12_6_equivalence`, `equivalence` |
| `[state]` | `gate4`, `gate5_effective`, `gate5_work_started`, `next`, `evaluator_execution_ready`, `evaluator_runtime_ready`, `runtime_v4_certified`, `todo13_authorized`, `t14_authorized`, `downstream_authorized` |
| `[outputs]` | `state`, `execution`, `expected_count`, `present_count`, `expected_names`, `absent_names`, `present_files`, `output_hashes_recorded` |
| `[authority]` | The proposal, plan, candidate, model, T8, fresh T11/T12 package, Gate4 postimage, transition receipt, post-review, and Oracle-review path/hash/mode bindings described below |

The shared `[authority]` fields are the `*_path`, `*_sha256`, and (where
sealed) `*_mode` triples for `proposal`, `plan`, `candidate`, `model`,
`t8_facade`, `t8_source_bundle`, `t8_claim`, and `t8_review`. The package
archive fields are exactly `t11_source_tar_path`,
`t11_source_tar_sha256`, `t11_source_tar_mode`, `t11_source_manifest_path`,
`t11_source_manifest_sha256`, `t11_source_manifest_mode`,
`t11_source_symlink_manifest_path`, `t11_source_symlink_manifest_sha256`,
`t11_source_symlink_manifest_mode`, and the corresponding `t12_` fields. The
package authority fields additionally bind the claim, review, publication
receipt, payload/evidence manifest, and T12 runtime-waiver path/hash/mode
triples. Gate4 uses `gate4_expected_boulder_postimage_*`,
`gate4_transition_receipt_*`, `gate4_post_review_*`, and
`gate4_oracle_final_review_*` path/hash/mode triples.

The v4 authority does **not** bind mutable tracked evaluator or test source
files. For each fresh T11 and T12 authority package it directly binds the
immutable `source.tar`, `source-manifest.tsv`, and
`source-symlink-manifest.tsv` members, including path, SHA-256, and mode. The
sealed T11 bindings are:

| Member | Path | SHA-256 | Mode |
|---|---|---|---|
| `source.tar` | `.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/v7-run/input/source.tar` | `aaf515357901c3276edf867c04ecfd0229b494f836ddd72ef3f0effed9f367cc` | `0444` |
| `source-manifest.tsv` | `.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/v7-run/input/source-manifest.tsv` | `2e7696824e208ba0c3ad960a08e452f81d79f892c776a1ebd0b28cbf1a196714` | `0444` |
| `source-symlink-manifest.tsv` | `.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t11-full-chain-julia-1.12.7-v1/v7-run/input/source-symlink-manifest.tsv` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `0444` |

The sealed T12 bindings are:

| Member | Path | SHA-256 | Mode |
|---|---|---|---|
| `source.tar` | `.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/v9-input/source.tar` | `e3538df88145753fc4eb65b1b5a0032e01391ef35f2ce54f16a1aede837b52b6` | `0400` |
| `source-manifest.tsv` | `.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/v9-input/source-manifest.tsv` | `ab9deb0c68f349a1908d9dba81a5cc5187d85332846ba9f0e9855b88bff9a3b2` | `0400` |
| `source-symlink-manifest.tsv` | `.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t12-full-chain-julia-1.12.7-v1/evidence/v9-input/source-symlink-manifest.tsv` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `0400` |

The remaining `[authority]` fields bind the T11/T12 claims, sidecars, reviews,
publication receipts, payload/evidence manifests, and the T12 runtime waiver;
the Gate4 fields bind the expected Boulder postimage, `TransitionReceipt`,
independent post-review, and Oracle final review. These are sealed provenance
records, not mutable source inputs.

#### Runtime guard and fail-closed state

The runtime guard checks Julia `1.12.7`, the exact runtime-v4 config bytes, the
binary and sysimage SHA-256 values, and canonical executable identity including
`/proc/self/exe`. The configured executable is
`/raven/u/system/soft/SLE_15/packages/x86_64/julia/1.12.7/bin/julia` with
SHA-256 `582718c20c563d824b88a08e2a9d0afcf5015b3053f425697cecc79ec4e9ea36`;
the configured sysimage is
`/raven/u/system/soft/SLE_15/packages/x86_64/julia/1.12.7/lib/julia/sys.so`
with SHA-256
`03041fa63edac21123a7d4ced06364243c5b9653f4ceb6ba618f69876eefa66f`.
The agents' `/u/...` namespace does not canonicalize to the configured
`/raven/u/...` path, so the production guard correctly blocks there. Positive
path/hash behavior is tested only through controlled seams; there is no
production bypass.

The state is fail-closed: `authorizes_todo13=false`,
`evaluator_execution_ready=false`, `evaluator_runtime_ready=false`,
`runtime_v4_certified=false`, `t14_authorized=false`, and
`downstream_authorized=false`. The `[outputs]` table records
`state = "ABSENT_AT_GATE5_CONFIG_PUBLICATION"`, nine expected evaluator
artifacts, `present_count = 0`, an empty `present_files` list, and
`output_hashes_recorded = false`. No evaluator output or product hash exists.
The candidate used no real data, labels, expected `N`, `NKNNKN`, class counts,
composition prior, benchmark/grader, or runtime-equivalence result.

The 2026-09-09 clean-checkout verification exposed six accidental `.omo` reads
in the earlier repository-local `82/82` synthetic run. After moving the 24
live-file assertions into the authority-backed suite, the isolated synthetic
contract passes `58/58` and the full suite passes `66,989/66,989` on two
consecutive runs, with identical stdout and no failures, errors or broken tests.
Source loading also passes. These test results do not establish hosted CI or
Gate2 approval. The machine-readable
`test/structured_evaluator_gate5_coverage.toml` records retained lower-level
coverage; only the impossible integrated-positive historical status is retired.
`Manifest.toml` is ignored by Git. The evidence archive must therefore include
it explicitly; a plain tracked-files export is insufficient for Pkg replay.
CI resolves its own dependencies, rejects changes to tracked `Project.toml`,
and hashes both environment files after bootstrap and after the synthetic test.
This detects test-time drift, not reproducibility against the local ignored
Manifest. The older closure observation retains its original input bindings;
it is not a review of the revised tests, coverage map or workflow.

Gate4 and the Gate1 transition are final, so Gate5 implementation may continue.
Gate2 Oracle integration review remains pending. Todo13 evaluator execution,
its outputs, T14, and downstream work are unauthorized and unstarted.

#### Inactive T13 activation successor contract v1 (2026-09-10) — historical, superseded by v2

**Historical record.** This v1 contract is superseded by the proposed
Main-control successor v2 below. It must not be read as the current candidate
or as an execution fallback; it is retained to document the old `ed1ac6b1`
bytes and the rejected execution-root-owned Boulder design.

The accepted T13 successor is an inactive candidate contract, not an active
grant. The entire
`config/unit_assignment_structured_evaluator_runtime_v4.toml` remains
unchanged, and the ordinary invocation, default evaluator CLI, and public
unary-selection producer remain fail-closed. The successor adds only a paired
optional CLI route: `--activation-receipt <repository-relative-path>` together
with `--activation-sha256 <independently approved expected digest>`. With both
absent, the baseline stays blocked; with either one alone, or with duplicate or
unknown flags, the invocation fails with `cli_error`.

Validation uses strict TOML records with the exact schemas
`stmfit_t13_execution_spec_v1` and `stmfit_t13_activation_receipt_v1`, plus the
synthetic scope string `T13_SYNTHETIC_EXECUTION_ONLY` and the environment
marker `STMFIT_T13_SYNTHETIC_EXECUTION_ONLY=1`. The activation receipt has
exactly the top-level keys `schema`, `scope`, `spec`, `source_review`,
`transition_review`, `boulder`, and `publication`; the execution spec has
exactly `schema`, `scope`, `root`, `entrypoint`, `config`, `project`,
`manifest`, `source`, `command`, `launch`, `destination`, `inputs`,
`dependencies`, `predecessors`, and `permissions`. Each path/sha256 pair is
repository-relative and normalized, and no referenced path may be the receipt
itself. These records bind the fixed 19-file consumed include closure
(`source.files`), the unchanged runtime-v4 config, `Project.toml`,
`Manifest.toml`, the canonical absolute root, the five file inputs, the three
input directories with complete member inventories, the exact canonical
command body and environment, the nine frozen report names, and the
predecessor and dependency records described below.

Boulder state is bound natively: both `boulder.preimage_path` and
`boulder.postimage_path` must be exactly `.omo/boulder.json` under the bound
root, the preimage digest is the fixed Gate4
`e509c919734c8110ab7ad4a7e2d475a552a0f2fca5f3533c9863683be1ac9133`, and the
postimage digest must differ. A detached or relocated copy is rejected even
when correctly hashed. Julia snapshots the live bound-root file and retains it
for revalidation; it does not parse or reinterpret the JSON state.

The execution spec binds the launch and process state it was reviewed under:
`launch.active_project` must be the bound-root `Project.toml`, `launch.cwd`
must equal the running cwd, `startup_file`/`history_file` must be false,
`threads`/`blas_threads` must equal the Julia and BLAS thread counts,
`depot_path`/`load_path` must equal the active depot and load path, `offline`
must equal the `JULIA_PKG_OFFLINE` state, and `forbidden_environment` must be
exactly `STMFIT_DATA_DIR`, `STMFIT_T13_PUBLIC_REVALIDATION`, and `JULIA_DEBUG`,
each absent or empty. These are observable-process-state integrity and scope
checks, not proof of process identity.

Predecessors are explicit: Gate4 expected Boulder postimage, Gate4 transition
receipt, Gate4 post-execution review, Gate4 Oracle final review, Gate5
publication receipt, and the completed Gate2 Oracle review
(`.omo/evidence/structured-label-free-unit-assignment/runtime_v4/t13-gate5-gate2-oracle-review-v1/Gate2OracleReview.json`,
SHA-256 `fe0b60720af63655f5415def50ef79a3096858840ed9d4f386e24569544e6fbe`).
The earlier Gate5 publication receipt alone is insufficient. The new control
records (activation receipt, execution spec, source review, transition review,
transition publication) require captured read-only mode `0444`; Julia compares
captured modes and never changes them or seals anything.

The finite additional read surface is bound as `dependencies`. `files` is the
closed set of provenance-only and receipt-derived content reads: the fixed
27-path provenance list (plan; T2 and T7 reviews; T7/T10 tests; T10 CLI; T8
build script; and 20 provenance-rebind claim, review, and canonical
source-bundle records, which include the T3 and T10 reviews), plus the producer
feature table and the forward and backward patch tables. `directories` is
`test/lib/hierarchical` with complete membership and per-member hashes and
modes. `identity_paths` covers the two
producer scripts (`kind = "file"`) and the producer `data_path` raw-data
directory (`kind = "directory"`); their bytes are not read. Cross-bindings are
checked against the already-bound universe and edge receipts (feature hash and
forward/backward patch hashes), and overlapping roles deduplicate only when the
normalized path, declared hash, mode, and captured identity agree. Missing,
extra, duplicate, substituted, or cross-bound-mismatched entries fail closed
before any publication-capable context exists. The T11 recomputation and T12
consumption are part of the consumed read surface but are not executed in this
inactive phase.

Activated output uses the successor report schema
`structured_label_free_unit_assignment_evaluator_v1_activation_receipt_v1`
(version 1) and the successor blocker schema
`structured_label_free_unit_assignment_evaluator_v1_activation_blocker_receipt_v1`;
the legacy v2 serializer remains byte-identical for non-activated runs. The
eight scientific TSV schemas, nine filenames, and publisher policy are
unchanged. The activated receipt carries the activation binding fields
(`activation_schema`, `activation_authority_schema`, `activation_scope`,
`activation_receipt_sha256`, `activation_spec_sha256`,
`activation_source_review_sha256`, `activation_transition_review_sha256`,
`activation_boulder_preimage_sha256`, `activation_boulder_postimage_sha256`,
`activation_publication_sha256`, `activation_entrypoint_sha256`,
`activation_destination`) and deliberately omits the legacy
`gateclosure_sha256`, `preclosure_*`, `review_sha256`, `closure_v1_*`, and
`live_plan_boulder_runtime_authority` fields, which do not certify an activated
run.

The trust boundary is explicit. The `--activation-sha256` digest is an
integrity pin; the provenance that makes it trustworthy is the external,
independently reviewed ACT-3 launch, not the digest itself. Julia compares the
receipt bytes to that supplied pin and checks scope, order, and bindings; it
performs no cryptographic signature or signer-identity verification, and no
caller-selected digest authorizes anything. The order is acyclic: successor
source/specification -> ACT-2 source review -> Boulder preimage, transition
patch, expected postimage -> ACT-3 transition review -> apply transition ->
earlier transition publication record -> final activation receipt last ->
independently pinned launch. The command body deliberately excludes
`--activation-sha256`, and no record contains its own or a future record's
hash.

Gate5 and Gate2 certify earlier bytes only; any successor requires new hashes
and new evidence. Boulder remains
`e509c919734c8110ab7ad4a7e2d475a552a0f2fca5f3533c9863683be1ac9133`, and no
production T13 run is authorized. The first future grant, if any, covers only
frozen synthetic T13 execution; real campaigns, the public T14-facing producer,
T14/downstream work, grading, and Julia 1.12.6 equivalence are excluded. The
successor reuses the existing nine-file publisher and the eight scientific TSV
schemas; its administrative provenance succeeds the historical closure-v1
observation without rewriting it. No scientific formula, parameter,
calibration, threshold, GCV, `n_eff`, or label-firewall behavior changes.

This section is inactive baseline metadata: it describes the frozen inactive
candidate, not a grant. ACT-2 remains BLOCKED, no activation receipt or live
Boulder transition exists, and the section authorizes no execution.

#### Proposed inactive Main-control successor v2 (2026-09-10)

The v2 candidate is an implementation candidate only: the code exists and its
remediated revision has passed fresh validation on 2026-09-11 (see "Tests and
limits" below), but MC-2 has not reviewed or accepted it. It is not deployed,
and the v1 section above is superseded historical design that must not be read
as the current contract or as an execution fallback.

**Two roots.** The execution/source/input root `R` owns the running code,
`Project.toml`/`Manifest.toml`, inputs, provenance, and outputs. The canonical
Main control root `M` owns the single live `.omo/boulder.json`. The v2 route
binds exactly that one file, rejects v1 spec/receipt schemas by equality with
no fallback, and never creates or consults an `R`-local `.omo/boulder.json` for
control bytes: an `R`-local file is a non-consumed decoy. No mirror, symlink,
hardlink, bind-mount, static-copy fallback, or automatic replication is
permitted, no operational `R` ledger is created, and a missing `M` fails closed
even when an `R` decoy holds the reviewed bytes. The unchanged runtime-v4
config, default denial, and public unary-selection producer denial still
govern.

**Schemas.** Execution spec `stmfit_t13_execution_spec_v2`; activation receipt
`stmfit_t13_activation_receipt_v2`. The spec top-level keys are exactly
`schema`, `scope`, `root`, `control`, `entrypoint`, `config`, `project`,
`manifest`, `source`, `command`, `launch`, `destination`, `inputs`,
`dependencies`, `predecessors`, and `permissions`. The new exact `[control]`
table has exactly `role` (`main_boulder`), `root` (canonical absolute `M`),
`root_mode`, `parent_mode`, and `file_mode`; the fixed control filename is
`.omo/boulder.json`, and there is no CLI control-root override, arbitrary
filename, multi-root whitelist, or fallback discovery. The receipt's `boulder`
table has exactly `control_root` (must agree with the spec `[control].root`),
`preimage_path`/`preimage_sha256` (fixed `.omo/boulder.json` and the Gate4
`e509c919…`), and `postimage_path`/`postimage_sha256` (same filename; reviewed
post digest that must differ from the preimage). The receipt carries the
reviewed post digest; the spec carries no future Boulder digest and the receipt
binds the spec hash, so the spec -> receipt -> externally reviewed launch pin
order stays acyclic. No receipt field contains its own hash.

**Dedicated control object.** A small `_MainControl` object holds the one fixed
`M/.omo/boulder.json` snapshot, the approved `root_mode`/`parent_mode`/
`file_mode`, the captured `M` and `M/.omo` device/inode identities, and the
preimage/postimage digests. It is deliberately separate from the `R`-relative
data snapshots. The required contract brackets the captured Main around the
read: complete pre-read checks, the read, then complete post-read revalidation
of the `M`/`M/.omo` directory modes and identities and of the file identity,
mode, and bytes. There is no `chmod` repair and no identity refresh.
`_verify_main_control` repeats those checks plus the byte check at aggregate
authorization return, at every context verification (preflight/assembly return,
before and after the expensive loader stages, before and after the manifest
rebuild, at loader return, and on loader/error paths), and at every publication
verification including the new preferred-rename destination-exists branch. The
same captured object is carried through authorization, authority collection,
context assembly, loading, blocker/error serialization, and publication; a
same-byte atomic replacement of `M` (new inode) invalidates, and no fresh
ordinary snapshot can rebase the captured control. The distinction is
initial capture versus change after capture: a replacement that happens before
the initial capture is a valid new baseline and the capture describes the
post-replacement state, while only a change after the capture is a revocation.
The initial MC-2 review found the post-read directory and file revalidation
incomplete. The bounded, non-scientific remediation completes it: the
post-read file identity/link/mode recheck and the captured `M`/`M/.omo`
directory identity and mode rechecks now run before `_read_main_control`
returns, `_verify_main_control` brackets the file-byte check with the captured
root/parent rechecks, and the original captured `_MainControl` object and
identities are retained with no rebasing, refresh, or `chmod` repair. Fresh
focused evidence exercises the initial-acquisition versus post-capture
replacement distinction and the retained inode; MC-2 attempt 2/3 review is
still pending. All other reads stay confined to `R` by the unchanged
path/snapshot confinement, and the running entrypoint, active project, cwd,
executable, and sysimage guards are unchanged.

**Authority manifest.** Activated contexts use
`schema=structured-evaluator-authority-v8-t13-main-control`. Ordinary records
stay strictly `R`-relative (an escaping path fails closed) and the manifest
appends exactly one fixed tagged `main_control` record carrying `role`, `root`,
`file`, the three modes, `sha256`, and `bytes`. No inode number and no `../`
path enters the deterministic content digest. The control file is exclusive
control state, not ordinary data: it may not also appear as an ordinary
`R`-relative snapshot or input, including when `M == R` is the actually
reviewed Main, in which case exactly one tagged control record is emitted and
the control path is excluded from ordinary records. The initial MC-2 review
found that exclusion not enforced; the bounded, non-scientific remediation now
rejects the overlap as `activation_binding_mismatch` before authorization is
usable, at every activation-context revalidation, and at tagged-manifest
construction, and it rejects the `M == R` control-as-input case before any
context exists. MC-2 attempt 2/3 review is still pending. The non-activated
runtime-v4 manifest keeps its historical v6 schema, and the legacy serializers
and the eight scientific TSV schemas are untouched.

**Activated output.** Activated report and blocker schemas are true successors,
`structured_label_free_unit_assignment_evaluator_v1_activation_receipt_v2` and
`structured_label_free_unit_assignment_evaluator_v1_activation_blocker_receipt_v2`,
with `schema_version = 2` and the Main-control provenance fields
`activation_control_role`, `activation_control_root`,
`activation_control_file`, `activation_control_sha256`,
`activation_control_root_mode`, `activation_control_parent_mode`, and
`activation_control_file_mode`. The legacy `gateclosure_sha256`,
`closure_v1_*`, plan, and preclosure Boulder fields remain absent from
activated receipts because they do not certify an activated run.

**Publisher.** One narrow addition is made to the otherwise unchanged
publisher: in `_publish_atomic`, the preferred-rename `destination_exists`
branch calls `_verify_context_for_publication(context, :not_committed)` before
returning `committed_verified` when the destination already holds byte-identical
files. A precommit Main-control invalidation therefore surfaces as
`_PublicationError(:not_committed, <existing reason>)`, never a false verified
publication. The rename/fallback algorithms, the nine filenames, receipt-last
ordering, and postcommit `:ambiguous` semantics are unchanged; postcommit
invalidations remain ambiguous/error with no automatic rollback. This is one
explicit publisher check addition, not a claim that the publisher is completely
unchanged. MC-2 confirmed this `destination_exists` addition correct, and no
further publisher algorithm change is planned.

**Scoped consultation and trust.** "Continuously consulted" means stage-boundary
consultation of the same live `M`: before and after the expensive T8/T11/T12
loader stages and the manifest rebuild, and at the authorization, error, and
publication boundaries. The remediation makes those interstage checks explicit
in `_load_production_input` at the T8/T11/T12 boundaries, including after
`T11.evaluate_admission`, after the admission file-set comparison, before T12
inference, and after `T12.infer_structured_chains`. It is not a continuous
instant poll, does not preempt inside an opaque numerical call, and does not
promise an atomic transaction between Boulder and outputs. A detected Main change prevents
verified publication; postcommit detection stays ambiguous/error; immediate
mid-call interruption still requires job cancellation. Any `M` byte or inode
change can invalidate, including unrelated administrative fields; no
semantic-revocation-only filter is claimed. No new error reason symbol is
introduced: Main-control failures reuse the existing vocabulary
(`activation_control_mismatch`, `activation_path_invalid`,
`activation_binding_mismatch`, `authority_path_mismatch`,
`activation_snapshot_changed`, `authority_hash_mismatch`,
`authority_hardlink`, and the existing publication classification). Trust stays
external: the reviewed spec and the actual parent-approved launch pin the real
Main path and filesystem view; a role label or caller-selected digest is not
authentication, and no signer, signature, reviewer identity, or full-process
attestation is claimed.

**Tests and limits.** The focused tests build two disposable roots under the
allocated `TMPDIR`: an execution root `R` with copied modules and a separately
allocated fake controller Main root `M`, with an `R`-local non-consumed Boulder
decoy. Neither root is the real checkout and neither is an operational ledger.
The frozen remediation revision (`cb5aacb8…` evaluator, `fb4faf31…` activation
tests, `939e17a3…` coverage map) passed static `Meta.parseall`/TOML parsing and
five fresh validation commands on 2026-09-11: targeted `227/227`, focused
`414/414`, synthetic `58/58`, existing full `67002/67002` with the synthetic
contract flag absent, and firewall `379/379`, all exit 0 with zero
fails/errors/broken and the twelve frozen source/control/docs hashes and
`git status` unchanged. The real postcommit hook fired once on
`after_receipt_link`, the Main control hash changed with the mode restored, the
outcome was ambiguous with `authority_snapshot_changed`, and all nine committed
destination files were preserved with no rollback; the CLI-through-Main
neutral-root case reached the actual activation receipt and failed with
`authority_hash_mismatch` rather than the earlier forbidden-CLI firewall. The
MC-2 attempt 1/3 failure (352 pass / 13 fail / 0 error / 0 broken, 365 total,
exit 1) is preserved as historical evidence for the pre-remediation bytes, no
failed case is relabeled as passing, and MC-2 attempt 2/3 review of the
remediated bytes remains pending with no MC-2 PASS claimed. The old
`ed1ac6b1` counts (164/58/66989/379) remain old evidence only. The known-heavy
positive T8/T11/T12 pipeline remains untested. Mock producer receipts with
literal `PASS` and Viper/Slurm fields are declared synthetic protocol
fixtures, not proof that the named extractors ran; only planning is authorized,
and no input generation, submission, ledger initialization, transition, or T13
execution is permitted. Historical T2/T7 files are not repaired; their
original bytes may later be copied into a separate compatibility tree under a
separately reviewed protocol.

This section documents an inactive remediated implementation candidate. It is
not an MC-2 acceptance, creates no ledger or grant, and authorizes no
execution.
