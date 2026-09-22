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

## Reconstructed Julia unit assignment

`config/unit_assignment_reconstructed.toml` is separate from the unchanged
counting/calibration config. It names `cc_soft_reconstructed_v1` and fixes the
new descriptor and numerical port before benchmark comparison:

- `[model] descriptor*`: backward-residual half-plane asymmetry on the specified
  9×9 grid; `descriptor_zero_l1` rejects zero/near-zero signal.
- `[model] split_skew_ratio_max`: bound for the split-width feature refit.
- `[model] mold_*`: template patch grid, target height, normal sampling interval
  (upper bound exclusive), and the original first-below-target isovalue scan.
- `[model] fisher_*`: PCA dimension, noise/GMM covariance regularization, native
  GMM iteration limit/tolerance and deterministic seed.
- `[selection]`: k-means/GMM seed counts, initial seed, GMM self-training,
  interaction features and soft-vote threshold. No composition prior.
- `[preprocessing]`: extractor normalization, `u_outer_t_inner` descriptor pixel
  order, the historical Fisher array layout, and explicit `patch_residual_filter`.
  `smooth_data_only` preserves `S(data) - model`; `smooth_residual` computes
  `S(data - model)` with the same native box smoother and count-config radius.
  Production assignment configs must state this field; existing custom configs
  need `smooth_data_only` to retain their previous behavior.

These are fixed settings, not knobs to search against the benchmark. Invalid
patches remain explicit missing inputs. Counting continues to use its own
existing `[model]`, `[selection]`, and `[preprocessing]` config unchanged.

The opt-in `config/unit_assignment_transverse_fisher.toml` differs only in the
method name (`cc_soft_transverse_fisher_v1`) and `[preprocessing] fisher_layout`.
`legacy_row_major` preserves the historical physical-t mirror;
`physical_u_outer_t_inner` reverses physical u in the extractor's u-outer/t-inner
serialization. Both keep the same 197 disk pixels and identical Fisher fitting.
The physical-u mode requires a reflection-closed disk (otherwise it errors,
rather than adding zero pixels). This is one scientific comparison, not a new
default, parameter sweep or change to the reconstructed half-plane descriptor.
The September 20 full145 comparison is negative (669 correct / 27 exact vs
671 / 28 for the symmetric-fusion control); keep the original layout for the
working candidate. The opt-in config records the unsuccessful experiment.

The separate opt-in `config/unit_assignment_matched_residual.toml` changes only
the method name (`cc_soft_matched_residual_v1`) and
`patch_residual_filter = "smooth_residual"` relative to the reconstructed config.
It keeps the legacy Fisher layout and every numerical setting. The pipeline
passes `--assignment-config` to both extractors and rejects all three patch-cache
options in matched mode, so legacy patches cannot silently be reused. Geometry,
split-width features and templates may still be reused. Standalone extractors
without `--assignment-config` retain legacy extraction for compatibility.
The completed fixed comparison gains 671/28 → **672/29** (correct positions /
exact chains), at unchanged coverage. This was the working candidate before the
descriptor comparison below, not a replacement for historical 677/36. No setting
is tuned after grading; the reconstructed config still selects the legacy mode.

Two further opt-in configs test only the descriptor on those same matched
residuals: `unit_assignment_transverse_moment.toml` and
`unit_assignment_affine_residual.toml`. Each differs from the matched-residual
control solely in `model.name` and `model.descriptor`; numerical settings,
column schema and all other preprocessing/selection settings are identical.

- `transverse_half_plane_asymmetry` remains the control: `sum(sign(u)*p)/sum(abs(p))`.
- `transverse_first_moment` uses `sum((u/max(abs(u)))*p)/sum(abs(p))`. It distinguishes
  transverse displacement within one half-plane and is dimensionless.
- `affine_residual_half_plane_asymmetry` first projects the complete normalized
  9×9 patch off `[1,t,u]` by least squares, then applies the original half-plane
  formula, including the projected patch's L1 denominator. A projected L1 mass
  at or below the unchanged `descriptor_zero_l1` becomes NA with
  `zero_affine_residual_mass`; no missing pixel is imputed and no row is dropped.

The storage column remains `patch_u_asym_reconstructed` for predictor-schema
compatibility. The explicit method name/config identifies its definition; neither
candidate is asserted to recover the lost historical producer. The plane
projection affects this one descriptor, not images, patches, CC/Fisher or
k-means inputs. These two predeclared comparisons are not promoted defaults.
Their completed full145 grades are **672/849, 28 exact** for the first moment
and **673/849, 33 exact** for affine residuals, versus control 672/849 and 29.
That comparison retained `unit_assignment_affine_residual.toml` as the opt-in working candidate;
the first moment is a recorded negative experiment. The affine candidate still
trails historical 677/854 and 36 exact. All sixteen final changes end at the
unchanged 0.5 vote tie, not at a newly calibrated confidence. No threshold or
other setting changes after this grade.

Two further independent candidates start from the 673/33 affine-descriptor
config. `unit_assignment_signed_mold.toml` changes only the method name and
`model.mold_margin_mode = "signed_cost_difference"`; the forward/backward
features become `cost_GlcN-cost_GlcNAc`, not their absolute value. Positive
means the GlcNAc physical template costs less under the existing alignment.
No decoded template label, benchmark label or expected composition enters it.
The explicit legacy mode, `absolute_cost_margin`, reads the saved `cost_margin`
without recomputing its rounding.

`unit_assignment_affine_fisher.toml` instead changes only the method name and
`preprocessing.fisher_patch_projection = "affine_disk"`. A fixed QR basis on
the configured Fisher disk removes `[1,t,u]` at training and held-out scoring,
including both unchanged mirror orientations. Original center-pixel amplitudes
still orient the learned groups. `model.fisher_projection_zero_l1 = 1e-12`
marks projected zero-mass rows NA with `zero_affine_patch_mass`; all keys remain.
Existing configs explicitly choose projection `none` (the floor is inactive).
These independent candidates are not combined, defaults or promoted methods.

The completed September 21 comparison is negative for both: signed CC gives
**668/849 (78.7%), 32 exact**, affine Fisher **672/849 (79.2%), 32 exact**, versus
the exactly replayed **673/849, 33 exact** control. All retain 849/870 coverage
and fixed counts. It kept `unit_assignment_affine_residual.toml` as the opt-in
working candidate at that stage; the two new configs document rejected experiments. No sign,
threshold, seed, combined arm or other parameter is retuned after this grade.


### Independent numerical candidates (2026-09-21)

Two opt-in configs copy the 673/33 affine-descriptor control and each change
only the method name and one explicit setting. All native assignment configs
now declare the historical settings, without changing their arithmetic:

| `[model]` field | Historical value | Independent experimental value |
|---|---|---|
| `fisher_score_center` | `"legacy_centered_mean"` | `"training_mean"` in `unit_assignment_centered_fisher.toml` |
| `gmm_final_covariance` | `"ridge"` | `"ledoit_wolf"` in `unit_assignment_shrunk_gmm.toml` |
| `gmm_covariance_ridge` | `1e-6` | unchanged |

Training-mean centering subtracts the opposite fold's mean patch at Fisher
scoring. Historical mode subtracts the near-zero mean of already centered
training patches. Neither changes the learned Fisher direction, amplitude
anchor, parity folds or mirror. It is not the equal-weight midpoint of clusters.

Ledoit-Wolf acts **only on the last hard self-training covariance**, after the
existing two reassignment iterations. EM, earlier hard iterations, means and
free mixture weights retain the historical ridge. `gmm_selftrain >= 1` is
required. With centered columns `x_i`, `S = XX'/n`, `T = tr(S)/p I`, use
`beta = max(0, (mean(norm(x_i)^4) - sum(abs2,S))/n)`,
`lambda = clamp(beta/sum(abs2,S-T),0,1)` and
`(1-lambda)S + lambda*T + gmm_covariance_ridge*I`.
If `S == T`, lambda is zero. No coefficient is fitted to benchmark results.
See [Ledoit and Wolf (2004)](https://doi.org/10.1016/S0047-259X(03)00096-4).
Dependent lobes and learned memberships preclude claiming iid optimality or
calibrated confidence. The GMM CLI accepts `--config`; omitting it explicitly
loads the historical reconstructed config. These candidates are not combined
or defaults. The completed comparison gives **675/849 (79.5%), 33/145 exact**
for training-mean centering versus control **673/849, 33 exact**. Retain
`unit_assignment_centered_fisher.toml` as the opt-in working candidate at that stage,
not a champion: historical 677/870 and 36/145 still lead. Final covariance
shrinkage gives **665/849 (78.3%), 34 exact** and is not retained as a working
replacement despite its better conditioning. All native arms keep 849/870
coverage, fixed counts and seven unavailable full-cohort keys. No combined arm
or parameter change follows grading.

### Independent patch-support and final-score candidates (2026-09-21)

Both configs copy the centered-Fisher candidate (675 correct / 33 exact), changing
only the method name and one explicit setting. They are not combined or defaults.

| Section / field | Existing value | Independent candidate |
|---|---|---|
| `[preprocessing] assignment_patch_support` | `"full_square"` | `"complete_disk_symmetric"` in `unit_assignment_patch_support.toml` |
| `[model] gmm_final_score` | `"mahalanobis"` | `"gaussian_density"` in `unit_assignment_gaussian_score.toml` |

All nine existing native assignment configs explicitly retain the old values.
Missing/unknown settings fail validation. The support mode leaves complete
patches unchanged. Fisher accepts a partial square only when its actual 197-pixel
scoring disk is complete; invalid values outside that disk are unused. Missing
schema columns still invalidate the row. The affine 9x9 backward descriptor
requires all 49 integer-grid disk pixels (`u_index^2+t_index^2 <= 4^2`). Outside
that disk it keeps only observed pixels with all their u/t reflection partners.
It fits `[1,u,t]` and computes the same half-plane/L1 statistic on this symmetric
observed support. No filling, coverage threshold or second normalization is
introduced. An incomplete central disk remains unavailable. This mode is defined
only for the affine-residual descriptor; existing CC pixel handling is untouched.

The final Gaussian score is `log(weight) - 0.5*(d^2 + logdet(Sigma) + p*log(2pi))`,
with the existing `1e-8*I` factorization guard included in Sigma. It acts only
after hard self-training and requires `gmm_selftrain >= 1`. EM, both hard updates,
ridge covariances, free weights and feature scaling are unchanged. Cluster naming
still uses mean raw amplitude of final assigned members, so that naming can
change when final assignments change. Votes remain hard seed votes, not calibrated
Gaussian probabilities. Neither new mode reads labels or imposes composition.

The completed independent comparison retains **`unit_assignment_patch_support.toml`**
at **676/870 correct, 676/852 (79.3%), 34/145 exact** as the latest opt-in working
candidate. Three missing predictions become usable, but two existing decisions
regress: net +1 correct and +1 exact chain against 675/33. Historical 677/36 and
854/870 coverage still lead. The Gaussian-score arm gives **666/849 (78.4%),
32 exact** and is rejected. No mode combination or post-grade tuning follows;
defaults, selected N and the separate counting/application claims are unchanged.

### Complete-patch training, admissible-patch prediction (2026-09-21)

`[selection] assignment_training_support` is required in every native assignment
config: `"all_admissible"` preserves the existing calculation; `"complete_patches"`
is the sole change (besides name) in `unit_assignment_complete_training.toml`,
copied from the current 676/34 patch-support candidate. This is one bounded
experiment, not a default or champion promotion.

Fisher fits each lobe-parity fold only on finite full forward 17x17 squares,
then scores every admissible opposite-fold disk. GMM fitting requires complete
forward 17x17, backward 17x17 and backward 9x9 patches, plus its usual finite-view
checks. Per-file feature mean/sample-standard-deviation moments use only finite
values from these complete-patch rows; the same frozen moments transform all
rows before unchanged pairwise products. A file without finite training values
has no usable transformed view; it does not borrow partial-row or global moments.
The existing zero/undefined-standard-deviation numerical guard remains one.

GMM EM, hard updates and free weights use training rows only. The high-raw-amplitude
group name also uses only assigned training members and is frozen before scoring
partial rows. Partial predictions can therefore neither move nor rename a group.
Insufficient training or inadmissible scoring retains keys with NA, without
imputation. K-means and its normalization are unchanged. There is no equal-class
prior, robust scaling, scan weighting, new threshold or post-grade tuning.

The runner writes `training_support.tsv` with observed-pixel counts, from the
three extracted patch tables. Fisher/GMM require `--training-support PATH` only
in the new mode and validate its exact key set and counts. Existing configs reject
that extra input. Older Fisher-attribution diagnostics reject this new mode
because they do not implement its training masks.

**Completed result: negative.** Complete training gives **674/870 correct,
674/850 (79.3%), 33/145 exact**, versus support control **676/852, 34 exact**.
The 893 complete training rows leave one file without moments: `240818_019`
has only partial patches, so two formerly available decisions become abstentions.
One other chain loses exactness. Keep `unit_assignment_patch_support.toml` at
676/34; the complete-training config remains an opt-in recorded experiment,
not the working replacement. No partial/global normalization fallback, different
eligibility rule or other post-grade change is introduced.

### Robust per-scan GMM normalization (2026-09-21)

Two explicit `[preprocessing]` fields now describe native GMM feature scaling:

- `gmm_feature_normalization = "mean_sample_std"` preserves all twelve earlier
  configs. The new `unit_assignment_robust_normalization.toml` copies the 676/34
  support candidate and changes only its name and this field to `"median_iqr"`.
- `gmm_scale_fallback = 1.0` exposes the existing degenerate-scale guard.
  Both fields are required; the fallback must be finite and strictly positive.

Each feature is centered/scaled separately within its file, using the same finite
values and training eligibility as before. The robust formula is
`(x - median(x)) / (Q75(x) - Q25(x))`, with Hyndman-Fan Type 7 quantiles
(`alpha=beta=1`). A zero or nonfinite IQR uses the declared fallback; a positive
IQR is not floored, even if small. Empty support remains NA. There is no clipping,
imputation, global fallback, Gaussian-consistency factor or quantile sweep.
The 28 pairwise products are formed after scaling the same eight descriptors.

This one candidate retains `all_admissible` training, including usable partial
patches. Fisher, k-means, pixel normalization, free GMM mixture weights, seeds,
hard updates, ridge, physical amplitude naming and voting are unchanged. The
completed experiment is **negative: 667/870 correct, 667/852 (78.3%), 10/145
exact**, versus 676/34 at the same coverage. No scale fallback occurs in any of
the 1,168 scan/feature pairs; the largest absolute expanded feature nevertheless
rises from 5.44 to 335.47. Reject this candidate and keep the support config at
676/34. No IQR floor, clipping, consistency factor or other post-grade adjustment
is added. The config records a negative experiment, not a new default.

### Equal-scan GMM training weights (2026-09-21)

Required `[selection] gmm_training_weighting` is `"equal_lobes"` in all thirteen
earlier native configs. The opt-in `unit_assignment_scan_weighting.toml` copies
the 676/34 support candidate, changing only its name and this field to
`"equal_scans"`. For each view, count `m_s` usable training rows per represented
scan after eligibility and finite-feature checks. With `n` rows and `S` scans,
the observation weight is `n/(S*m_s)`: mean one overall, equal total `n/S` per
scan. This is proportional to `1/m_s`; no expected count or class frequency is
used. An unavailable scan contributes nothing and retains its output keys.

Weights enter the two-component k-means++ draws, k-means centroids/objective,
initial means/covariances, EM objective/responsibilities, both hard moment updates
and raw-amplitude group naming. EM responsibilities and scoring are conditional
on a lobe, so weights multiply their training contribution, not their individual
classification score. Mixture weights are estimated freely from weighted mass;
the initial 1/2 values are not a composition constraint. The final ridge
covariance divides by weighted mass (maximum likelihood), not a degrees-of-
freedom correction. Weighted Ledoit-Wolf is not implemented and is rejected.

Per-file mean/sample-std scaling, finite per-feature normalization support,
partial-patch eligibility, interactions, Fisher, the separate k-means head,
ridge, iteration budgets, seed integers, Mahalanobis scoring and final soft vote
stay fixed. The weighted draws can select different initial centers with the
same RNG seeds. No scan resampling or robust scaling is part of this comparison.
The completed comparison gives a **tradeoff: 675/870 correct, 675/852 (79.2%),
36/145 exact**, versus the exactly replayed 676/34 support control at identical
coverage. All seven changed decisions land at the existing zero-margin vote tie.
Retain the `equal_lobes` support config as the primary working reference; this
opt-in config records a mixed result, not a promoted default. Historical 677/36
remains unexceeded. No weight formula, seed, threshold or other setting is
adjusted after the grade.

### Continuous GMM seed aggregation (2026-09-21)

Required `[selection] gmm_seed_aggregation` is `"hard_vote"` in all fourteen
earlier native configs. The opt-in `unit_assignment_continuous_vote.toml` copies
the 676/34 support candidate, changing only its name and this field to
`"mean_membership"`. Every seed fits exactly the same model and names its
high-amplitude component from the same hard assignments. Only its contribution
to the average changes: from `argmax(resp) == high_cluster` (0 or 1) to
`resp[high_cluster]`, where `resp` is the already-computed normalized exponential
of the two final scores.

The candidate keeps the existing `log(weight) - Mahalanobis_distance²/2`
scores, without adding a covariance-volume term, temperature or calibration.
These normalized memberships are not calibrated chemical probabilities.
EM, both hard updates, covariance ridge, free composition, raw-amplitude naming,
mean/std scaling, equal-lobe training weights, partial-patch support, seeds,
interactions, Fisher and the separate k-means head remain unchanged. GMM output
keeps the existing eight-decimal serialization; the final mean-of-two-heads vote
still uses `>=0.5`, with the same unavailable-input abstention rule.

This is one authorized comparison, not a change of default or a combined
variant. Its completed result is **negative: 671/870 correct, 671/852 (78.8%),
26/145 exact**, versus the exactly replayed 676/34 control at identical coverage.
All 24 final decision changes are 1→0 from old exact ties: GMM=1 and k-means=0
become GMM<1 and k-means=0. One exact chain is gained, nine lost. Retain the
`hard_vote` support config as the working reference; the new config records a
failed experiment, not a promoted default. No precision, threshold, temperature
or naming adjustment follows. Historical 677/36 remains unexceeded.

### Whole-scan GMM bootstrap (2026-09-21)

Required `[selection]` fields `gmm_resampling`, `gmm_bootstrap_replicates` and
`gmm_bootstrap_seed` explicitly preserve `"none"`, `1`, `0` in the fifteen
earlier native configs. The opt-in `unit_assignment_scan_bagging.toml` copies
the 676/34 support control, changing its name, mode to `"whole_scans"` and
replicate count to **20**. Seeds for resampling are 0–19; each replicate keeps
the existing ten GMM initialization seeds 0–9. There is no seed search.

For each view, draw S scans uniformly with replacement from the S scans with
usable rows. Copy all usable rows of each drawn scan with its multiplicity;
natural row counts are retained, not equalized. Per-scan scaling is computed
once before resampling; an unsampled scan uses only its own usual normalization
for prediction. EM, hard updates and mean-amplitude physical naming use the
duplicated training rows only. Predict all valid rows, including unsampled scans.
Average each replicate's valid binary seed votes, then average valid replicate
means equally. Unnamed seeds/replicates are omitted without retry; a row with
no valid result stays unavailable. No chemical class population is imposed.

This mode currently requires `equal_lobes`, `hard_vote`, `all_admissible`,
`mean_sample_std` and ridge covariance. The root pipeline writes
`gmm_scan_bootstrap.tsv` (view, replicate, draw seed, scan identity, usable rows,
multiplicity, training rows and valid/total seed counts). Direct GMM CLI use
requires `--bootstrap-audit NEW_PATH`; it rejects an existing audit and a path
equal to the prediction output. The historical representation diagnostic rejects
bootstrap rather than silently describing unresampled training.

Final voting, threshold, precision, unavailable-input handling, k-means, Fisher,
CC, geometry and N stay fixed. This is bagging, not out-of-bag validation or
probability calibration. The fixed comparison completes on Viper, job
**11925188**, **7m17s**, exit **0:0**, within its one-hour cap. All twenty bags
accept ten initialization seeds each. The result is **negative: 665/870 correct,
10/145 exact**, at unchanged **852/870 coverage**, versus the exactly replayed
676/34 control. It loses 24 exact chains and gains none; 57 of 62 final 1→0
changes leave old GMM=1 / k-means=0 ties. Retain the **676/34 support config**,
below historical 677/36. No replicate/seed search, fusion/threshold adjustment
or post-grade combination follows. See `results/scan_bagging_20260921/report.md`.

### Tied GMM covariance throughout learning (2026-09-22)

Required `[model] gmm_covariance_structure` is `"full"` in all sixteen earlier
native configs. The opt-in `unit_assignment_tied_covariance.toml` copies the
676/34 support control, changing only its name and this field to `"tied"`.
It is distinct from `gmm_final_covariance`, which remains `"ridge"`.

Tied mode shares one within-component covariance at k-means initialization,
every EM M-step and each of the two hard self-training updates. For n usable
training rows, covariance is the sum over rows and components of responsibility
times the outer product about that component's updated mean, divided by n,
plus the unchanged `1e-6` ridge. Initialization/hard updates use indicator
responsibilities. This pools scatter by observation mass, not equal component
weights and not the global covariance including between-component separation.
Initial mixture weights stay 1/2 as before; subsequent masses and both means
are learned freely. Two components do not impose a chemical composition.

This scoped mode requires `ridge`, `mahalanobis`, `equal_lobes`, `hard_vote`,
`none` resampling, `all_admissible` training and `mean_sample_std` scaling.
Final score arithmetic, naming by mean raw amplitude, seeds, interactions,
vote threshold, precision, missing-row policy, upstream signals and N are
unchanged. Historical representation diagnostics reject tied mode rather than
silently reporting the old model. One control and one candidate are authorized
on Viper, with four requested CPUs, 16 GB, a 30-minute limit and no requeue.
The completed result is **negative: 666/870 correct, 666/852 (78.2%), 6/145
exact**, versus the exactly replayed 676/34 support control at unchanged
coverage. It gains no exact chain and loses 28. All ten candidate seed fits
have identical component covariances and valid physical naming; masses remain
free. Reject this variant and retain **`unit_assignment_patch_support.toml`**,
still below historical 677/36. No other proposed lead, combination or post-grade
adjustment follows. See `results/tied_covariance_20260922/report.md`.

### Whole-scan Fisher and relative GMM naming (2026-09-22)

Three explicit `[selection]` fields preserve the seventeen earlier native
configs: `fisher_cv_scheme="lobe_parity"`, `fisher_scan_split_seed=0`, and
`gmm_cluster_naming="raw_amplitude"`. Two independent candidate configs copy
the support control; each changes its name and one mode only.

`unit_assignment_scan_fisher.toml` uses `fisher_cv_scheme="scan_hash_twofold"`.
Sort all unique scan basenames by SHA-256 of decimal seed + NUL + basename
(UTF-8), then basename for digest ties, and alternate ranks into groups 0/1.
Do this before invalid-row/training filtering. It balances numbers of scans,
not class proportions or lobe counts; seed zero is fixed without a search.
Every held-out scan uses only the other group's PCA/GMM/centering/naming/Fisher
fit. Training keys are sorted for deterministic order; missing or degenerate
training groups yield NA/reasons, never a self-trained fallback. Reordering
rows or reversing lobes does not change membership, but changing scan identities
or the cohort can. This does not make the downstream transductive classifier
or reused external benchmark independently validated.

`unit_assignment_relative_naming.toml` uses `gmm_cluster_naming="within_scan_z"`.
Only GMM group naming changes. On eligible, feature-valid training rows, raw
amplitude is centered and sample-std scaled within each scan, with the existing
`preprocessing.gmm_scale_fallback` for a constant or singleton scan. The global
group with larger mean standardized amplitude receives the existing class-1
name. Fitting, memberships, free masses, scoring, seeds and votes stay fixed;
there is no per-scan quota, forced second group or composition prior. The
scoped option rejects combinations outside the support-control policies.

Historical representation/attribution diagnostics reject unsupported grouping
or naming rather than silently describe the old model. The completed three-arm
comparison gives **668/870 correct, 34/145 exact** for whole-scan Fisher and
**676/870, 34/145** for relative naming, at unchanged **852/870 coverage**.
Both Fisher fits converge on disjoint 73-scan groups. Relative naming selects
the same high group for all ten seeds, with byte-identical fitted-parameter
hashes and GMM output; final scores/decisions stay unchanged (only the model
name differs). Neither config replaces **`unit_assignment_patch_support.toml`**
at 676/34; historical 677/36 remains the target. No combined arm, partition
search or post-grade tuning follows. See `results/scan_fisher_naming_20260922/report.md`.

### Fixed factor-analyzer and Student learning candidates (2026-09-22)

`unit_assignment_factor_analyzer.toml` and `unit_assignment_student_t.toml`
change only `[model] gmm_learning_family` (and the output model name) relative
to the support control. All native profiles now explicitly declare:

| Key in `[model]` | Fixed value / meaning |
|---|---|
| `gmm_learning_family` | `gaussian`, `factor_analyzer` or `student_t` |
| `gmm_factor_rank` | 4 latent factors per component; strictly below feature dimension |
| `gmm_student_df` | 5.0, shared fixed degrees of freedom, not fitted or searched |
| `gmm_learning_maxiter` | 200 alternative-family EM updates maximum |
| `gmm_learning_tolerance` | 1e-6 relative change of summed mixture log kernels |
| `gmm_learning_min_mass` | 1e-12 minimum component/precision mass; otherwise unavailable seed |
| `gmm_learning_cholesky_guard` | 1e-8, required to match the unchanged final-score guard |

These numerical settings govern the alternative families only; the Gaussian
route keeps its existing operations exactly. The family is retained through
the configured two hard updates. Both alternatives require separate component
matrices, ridge regularization, Mahalanobis final scores, parity Fisher,
raw-amplitude naming, equal-lobe weights, hard seed votes, no resampling,
all-admissible training and mean/sample-standard-deviation scaling. No
combination, rank selection, df estimation, new threshold or class-count prior
is implemented. The standalone builder requires `--selftrain >= 1`.

MFA stores `L*L' + Diagonal(psi)` with a separate positive noise diagonal for
each component, floored at the existing `gmm_covariance_ridge=1e-6`.
Student stores its **scale**, not its df-dependent covariance. Both use the
unchanged final `log(weight) - Mahalanobis²/2` vote; this is a learning-only
comparison, not Student posterior classification or chemical calibration.
Common density constants are omitted only inside alternative-family learning;
`log_kernel_sum` is not an across-family selection criterion. A finite fit
at the iteration cap is reported as nonconverged, not retried. Invalid seeds
are logged and excluded from the fixed-seed average; no usable seed means
`?`/`no_valid_view`, without dropping rows or falling back to Gaussian learning.
The historical representation diagnostic rejects both alternative families.

Completed comparison: **620/870 correct, 23/145 exact** for factors and
**631/870, 25/145** for Student, versus exactly replayed **676/34** support;
all retain **852/870 coverage**. Student converges for all ten seeds;
factors reach the predeclared 200-update cap for all ten. Retain support,
not either alternative. No setting changes after grading, including the cap;
the bounded factor result is not evidence about its fully converged optimum.
Source **e2205cb**, job **11935072**, report
`results/factor_student_mixtures_20260922/report.md`.

### All-update shrinkage and coherent Student decisions (2026-09-22)

Two isolated opt-in profiles extend the completed learning-only comparison:

| `[model]` key | Support / old profiles | `unit_assignment_em_shrinkage.toml` | `unit_assignment_student_density.toml` |
|---|---|---|---|
| `gmm_covariance_scope` | `final_only` | `all_updates` | `final_only` |
| `gmm_final_covariance` | unchanged | `ledoit_wolf` | `ridge` |
| `gmm_learning_family` | unchanged | `gaussian` | `student_t` |
| `gmm_hard_assignment` | `mahalanobis` | `mahalanobis` | `student_density` |
| `gmm_final_score` | unchanged | `mahalanobis` | `student_density` |

Both new keys are explicit in every native config. The legacy values preserve
the previous arithmetic. All-update shrinkage uses the existing spherical
Ledoit-Wolf estimator for the initial and hard groups, and its declared
fixed-responsibility plug-in extension at each EM M-step. Component covariances
remain separate. Its intensity is calculated from current observations, not
searched or chosen by recognition. It is numerical regularization, not an
iid-optimality or noise-calibration claim. See [the equations](calibration.md).

The Student profile requires matching hard/final `student_density` rules and
`student_t` learning; mismatched combinations are rejected. It uses the same
df-five scale model throughout, including component masses and matrix volumes.
The two hard updates, ten-seed hard vote and final k-means fusion stay in place;
this is not a continuous-vote or fusion-weight experiment.

The new Gaussian arm uses the explicit `gmm_learning_maxiter=200` and
`gmm_learning_tolerance=1e-6`, matching the previous Gaussian limits. These
settings now also reach the ordinary unresampled Gaussian builder; every
existing native config retains the identical values. Student retains its
existing limits, ridge and guard. Shrinkage is not a likelihood-maximizing
M-step: log-likelihood decreases and the largest decrease are reported, not
used to choose seeds. Finite capped fits are retained and identified as such.
No combination with tied covariance, alternative naming, weighting, resampling,
normalization or Fisher grouping is enabled. Selected N, features, thresholds,
DFT sources and unknown25 remain unchanged.

The completed comparison is negative: all-update shrinkage gives **637/870
correct, 24/145 exact**; coherent Student gives **627/870, 25/145**, both at
**852/870 coverage**, versus exactly replayed **676/34** support. All ten fits
per arm satisfy their stopping criterion. Keep the support profile; historical
677/36 remains unexceeded. No coefficient, df, threshold or fusion rule changes
after grading. Source **bb6943b**, job **11936229**; full evidence and losses:
`results/em_shrinkage_student_density_20260922/report.md`.

## Opt-in diagnostic exploration settings (2026-09-18)

`config/label_free_exploration.toml` is used only by standalone exploration tools;
it is not a production config. Pass the original molecule config separately.
Its empty `[model]`, `[selection]`, `[preprocessing]` tables make this separation
explicit and must not be interpreted as physical defaults.

- `[counting_variable_projection]`: `outer_maxeval`, `outer_maxtime_s`,
  `outer_xtol_rel`, `outer_ftol_rel` bound the profiled optimization;
  `linear_maxiter`, `linear_kkt_atol`, `linear_kkt_rtol`, `linear_bound_atol`,
  `linear_svd_rtol` control the bounded linear subproblem;
  `mapping_atol`, `mapping_rtol` check agreement with the native forward model;
  `native_elliptical_maxiter` records the native refinement budget explicitly.
  Full model parameter counts are retained in GCV.
- `[representation]`: `identity_atol` checks signed-mass/linear identities;
  `descriptor_atol` accounts for the saved descriptor's serialization precision.
  Neither changes descriptor values or production predictions.
- `[acquisition_noise]`: `channel`, `auxiliary_channel`; declared `lag_x_min_px`, `lag_x_max_px`,
  `lag_y_min_px`, `lag_y_max_px`, `lag_bound_basis`; registration sufficiency
  `min_registration_pixels`, `min_row_pixels`, `min_positive_correlation`,
  `ambiguity_correlation_gap`; molecular/background exclusion `footprint_sigma`,
  `background_guard_nm`, `min_background_pixels`, `min_background_fraction`;
  ACF reporting `acf_max_lag_px`, `min_acf_pairs`, `acf_threshold`; proposed
  block feasibility `block_length_multiplier`, `min_usable_blocks`,
  `block_min_occupancy`; local-view reporting `local_patch_radius_nm`,
  `min_patch_pixels`. An ACF that does not cross within its measured range is
  censored, not evidence of a precisely measured correlation length. Search
  boundaries and insufficient support are reported, not silently relaxed.
  A later resurgence in absolute ACF magnitude also leaves the length unresolved.
  `auxiliary_channel = "Current"` adds direct scale and cross-view coupling on
  the same Z support/lag, not a new likelihood, fit weight or nuisance projection.

The separate fixed Fisher replay uses the unchanged
`config/unit_assignment_reconstructed.toml` and the original forward17 `res`
patches. It adds no fitted settings and exports newly fitted fold coefficients,
not a recovered historical model.

These settings were declared without benchmark outcomes. No arbitrary `n_eff`
calibration, chemical prior, new production selection or automatic promotion is
introduced. See the dated research journal for the bounded scope and results.


## Masked robust background prototype (2026-09-18)

`config/masked_robust_preprocessing.toml` is an opt-in, **synthetic-first**
diagnostic config. Its empty `[model]`, `[selection]`, and `[preprocessing]`
tables do not supply a new production preprocessing policy. The standalone
prototype takes height values in nm, physical coordinate vectors, and an
explicit background mask. It never reads expected N, chemical labels or votes.

All controls in `[masked_robust_preprocessing]` must be supplied:

| Key | Initial value | Meaning |
|---|---:|---|
| `huber_delta` | 1.345 | Fixed dimensionless Huber cutoff; no benchmark tuning |
| `maxiter` | 50 | Maximum IRLS iterations; exhaustion remains nonconverged |
| `coefficient_rtol` | 1e-8 | Relative convergence tolerance in centered/scaled plane coordinates |
| `weight_atol` | 1e-8 | Absolute convergence tolerance for IRLS weights |
| `scale_floor_nm` | 1e-9 | Numerical zero-scale guard in nm, not a measured noise floor |
| `rank_rtol` | 1e-12 | Relative numerical rank tolerance for the plane design |
| `min_background_pixels` | 200 | Minimum finite supplied-background pixels for a plane |
| `min_background_fraction` | 0.05 | Minimum background count divided by all image pixels |
| `min_row_background_pixels` | 20 | Minimum observed background samples for a row offset |

The Huber scale is fixed to the normalized MAD of the initial background-only
OLS residuals, with the explicit numerical floor. It is not updated during
IRLS, not an independent noise estimate, and can be inflated by a contaminated
OLS initializer. A nonconverged or rank-deficient plane does not yield an
accepted corrected image. Final iteration diagnostics remain available.

The prototype fits a plane, then subtracts background-only row medians so that
usable background rows have zero median. This is sequential, not a simultaneous
identifiable fit of y tilt and arbitrary row offsets. Unsupported rows and raw
nonfinite pixels remain unavailable; no replacement values or smoothing are
introduced. The supplied exclusion is not relaxed to obtain a result. The
native reference retains its own global-median convention; comparisons must
report and account for that constant-level difference on the same observed
background support.

This config is not connected to a production driver. No N selection, chemical
threshold, GCV complexity, `n_eff`, physical bound or reference prediction changes.


Focused native checks (Julia 1.13):

```bash
julia --startup-file=no --threads=1 --project=. test/test_masked_robust_preprocessing.jl
julia --startup-file=no --threads=1 --project=. test/test_masked_preprocessing_signal.jl
# Optional fresh output directory: full synthetic arrays, metrics and used settings.
julia --startup-file=no --threads=1 --project=. test/test_masked_preprocessing_signal.jl \
  --outdir results/masked_preprocessing_example
```

The signal study reads this config, not a second hidden set of controls. The
optional output path must not exist, including as a dangling symlink. Saved
`initial_objective_nm2`, `objective_nm2` and `stationarity_inf_nm` describe the
plane stage on supplied finite background, before row-median subtraction; they
are not selection scores or noise estimates. Metrics retain explicit unavailable
anchor/support reasons and partial observed coverage. Fixed synthetic fixtures
are evaluation inputs, not learned physical calibration. See the journal for
verified results and adverse cases; there is no raw-data or HPC entrypoint yet.


### Four-scan masked-preprocessing comparison

`test/run_masked_preprocessing.jl` is a separate opt-in diagnostic driver. It
uses `--config` for the unchanged physical support settings,
`--acquisition-settings` for the frozen `config/label_free_exploration.toml`
registration/exclusion controls, and `--settings` for the unchanged
`config/masked_robust_preprocessing.toml`. There are no new numerical parameters
and no automatic choice of preprocessing method.

The four-file sample, original selected summary and base geometry are inputs,
not labels. They freeze the native support and the guarded background exclusion;
they do not independently validate that exclusion. The reused exclusion is the
original conservative square dilation by eight pixels, including more y margin
than the ±2-pixel registration window. Every method/direction remains in the
output, including unsupported rows and nonconvergence.

The comparison uses method-own, all-method common, and native-versus-each-method
common supports. Each support is valid over the entire unchanged 85-lag grid.
Physical level-aligned differences remove only one observable median per view
from the same shared off-footprint background pixels. These constants are fixed
across lags. The existing minimum background count/fraction controls determine
whether that level alignment is available; they do not calibrate uncertainty.
Background mean, MAD, standard deviation and changes versus native are post-fit
descriptions, not ground-truth errors or stationary noise estimates.

The per-case command `test/diagnose_masked_preprocessing_real.jl --help` documents
its array and table outputs. `arrays.jls` contains only built-in scalar, array
and dictionary data serialized by Julia 1.13. It retains observed raw nm arrays,
corrected images, masks and estimator diagnostics for saved-only checks. Load
only trusted files using Julia 1.13. No missing sample is turned into an observed
one. The existing native reader incidentally preprocesses auxiliary Current
arrays; the comparison discards them and computes no Current evidence.
