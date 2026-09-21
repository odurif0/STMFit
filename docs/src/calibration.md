# Calibration: deriving parameters objectively

The pipeline has ~25 calibration parameters. Most can be **measured** from a
single clean scan rather than hand-tuned, which makes the analysis generalizable
to a new molecule on the same STM. This page documents which parameters are
objective, which are principled choices, and which remain free.

## Auto-calibration

```bash
julia --project=. test/measure_calibration.jl path/to/clean_scan.sxm
```

This measures the objectivable quantities and emits a ready-to-use TOML.
Evaluate it externally on the benchmark *after* generating the TOML — do not
adjust measured parameters to recover benchmark labels:

```bash
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
    --config chitosan_auto.toml
```

## Parameter classification

### Measured from a single scan [objective]

| Parameter | Measurement method |
|---|---|
| `noise σ` | 1.4826·MAD of (raw − smoothed) high-frequency band, via the standard preprocessing pipeline |
| `pixel resolution` | `range_nm / width` (from the SXM header) |
| `FWHM range [lo, hi]` | Detect peaks in the chain-axis profile (weighted PCA → bright-pixel strip), fit half-max width per peak, take [25%, 95%] quantiles (25% excludes under-resolved outliers that would starve the fit) |
| `repeat spacing` | Median peak-to-peak distance along the chain axis |
| `spatial correlation range` | 2D isotropic autocorrelation on the full preprocessed image; first lag where ρ(h) drops to 1/e. Descriptive image correlation, not necessarily noise-only correlation. |

### Derived from a physical/numerical principle [principled, one fixed choice]

| Parameter | Derivation |
|---|---|
| `sigma_parallel_*` | `FWHM / 2.355` (Gaussian width relation) |
| `spacing_min/max` | `±30%` around the measured repeat spacing |
| `fit_width_nm` | `= 1.25 × σ_min` (tube half-width; the margin avoids lateral truncation of the narrowest lobe) |
| `support_min_length_nm` | `3 × spacing` (at least 3 repeats to call it a chain) |
| `n_max` | `longest_image_axis / spacing_min + 2` (generous cap; the chain may orient along either image axis) |
| `max_overlap` | 0.60 (Gaussian pair-overlap floor; sets the spacing lower bound) |
| `support_noise_k` | 2.5 (SNR threshold k·σ on the support envelope) |
| `support_padding_nm` | `= fit_width_nm` (pad by one tube half-width to avoid edge truncation) |
| `selection_criterion` | `gcv` (canonical practical score, with correlation/nonlinearity limits — see §Effective sample size) |
| `flatten` | `plane+rows` (STM scan-line + plane correction) |

### Free (not objectively measurable; left to default)

| Parameter | Why free |
|---|---|
| `global_maxtime`, `global_maxiter`, `max_iter` | Optimizer budget (numerical, not physical) |
| `chain_tilted_baseline`, `chain_circular_sigmas` | Model-form switches (domain choice) |
| `channel`, `direction` | Acquisition-dependent (Z topography by convention) |
| `selection_policy`, `gcv_ambiguity_rel_threshold`, `robust_guard_nu` | Selection-rule knobs (validated robust on [0.03, 0.06]) |

### Structured diagnostic policy

Todo 9's structured follow-up diagnostics are not calibrated from benchmark
labels. Their estimator semantics are frozen in the
`[selection.diagnostics.policy]` block of
`config/unit_assignment_structured_model.toml`; missing or differing keys
fail closed as `BLOCKED`. The policy fixes the fixed-ν=8 Student-t residual and
one-component fallback, date-centered unbalanced ICC(1,1), equal-feature
summed-squares pooling, the configured view contrast, explicit Hyndman–Fan
Type 7 quantiles, inclusive finite-sample tails, strict zero/equality handling,
and deterministic Holm ties.

Channel dropout is a separate follow-up statistic: it removes exactly
`bwd_neg_com_t` and `bwd_neg_diag45`, freshly refits `base_local` using only
the same inner-training partition, and rescored the untouched held-out rows.
It must never reuse or alias the original view contrast. The proposed
featurewise-standardized residual tests and hierarchical reliability model are
deferred to a separately preregistered v3 study and do not alter v2.

## Effective sample size — why GCV is the canonical criterion

The reported STM image/residual correlations are strong (ρ ≈ 0.9–0.95 at
lag 1; reported range 17–100 px). A narrow fit tube (~10 px across) may contain
too little independent spatial information to estimate a correlation model
reliably. Correlation measured over the molecular image also includes signal
structure: it is not automatically a background-noise covariance estimate.

The `n÷9` heuristic is a fixed placeholder, not a measured number of independent
observations. Replacing it with another ad hoc effective sample size can greatly
change absolute BIC/AICc values without improving the noise model. Their iid
absolute values must therefore not be interpreted as calibrated model evidence;
they remain secondary diagnostics/guards, not a replacement for canonical GCV.

**GCV** (`RSS·n/(n−p)²`) avoids inserting an arbitrary `n_eff` into the
per-candidate score and remains the canonical practical criterion. For a linear
smoother, generalized cross-validation approximates leave-one-out error by
replacing individual leverage corrections with their average; it is not the
exact leave-one-out identity in general. This code also uses a parameter-count
approximation in a nonlinear, constrained fit. Spatial correlation, active
bounds and model mismatch can therefore affect its predictive interpretation.
GCV does not by itself calibrate count uncertainty or establish chemical truth.
BIC/AICc remain secondary diagnostics/guards; the `n ÷ 9` placeholder is unchanged.
Numerically profiling amplitudes/background must not reduce the model parameter
count used in GCV.

## Calibrating a new molecule

1. Pick **one** clean, well-resolved scan of the new molecule.
2. Run `test/measure_calibration.jl <scan>` → produces `<scan>_calibration.toml`.
3. Inspect the measured values (especially FWHM and spacing — sanity-check
   against the visible structure).
4. Run the batch with the auto-calibrated TOML; spot-check N_selected on a few
   files visually.
5. If a parameter looks off (e.g. FWHM under-estimated on a noisy scan), measure
   on 2–3 scans and take the median.

Noise scale and correlation depend on acquisition, tip/feedback state and
preprocessing. They should not be assumed identical across sessions merely
because the instrument and flatten/smooth settings are unchanged. Estimates
from an image containing the molecule can also depend on its structure.
Recheck these quantities from suitable background or independent acquisition
evidence when transferring calibration, as well as molecule-specific FWHM,
spacing and `n_max`. This caution does not change the current production config.

## Structured evaluator-v1 is not physical calibration (correction3 pending review)

The evaluator config is a policy/static/synthetic correction3 prerequisite, not a
physical calibration and not an application or benchmark result. It becomes
authoritative only after parent acceptance, a fresh independent Oracle PASS, and
reviewer-owned `GateClosure`. Its
canonical path/hash/runtime/member checks must pass before any formula,
bootstrap, or graph work; Todo 13 remains blocked until that gate.

The exact graph reference is one partition/outer-score T12 block per held-out
date, with per-scan `logZ` summed only from that scan's selected blocks. T12
marginals, not Viterbi labels or report-wide sums, provide enabled probabilities;
disabled output uses selected unary `q`. Status rows distinguish `BLOCKED`,
`FAIL`, consumed `SKIPPED`, and `PASS`/`SKIPPED` final gates. Type 7 uses
`h=1+(n-1)p`, with 500-value lower quantiles exactly
`0.525*x_(13)+0.475*x_(14)`; sign masks use sorted dates, inclusive `>=`,
`count/2^K`, and no `+1`.

Evidence is limited to static policy checks and hand-computed/synthetic fixtures.
It creates no 10–20mer application claim and no benchmark claim. External labels
or grader outputs, if ever used, are permitted only after the independent policy
gate and outside this label-free prerequisite.

Correction2 distinguishes consumed diagnostic `SKIPPED` rows from terminal
statuses, resolves exactly one selected unary fit and matching T11/T12 partition
reference per held-out date, and computes entropy and view agreement over their
pooled frozen node/pair populations. Descriptor-based authority snapshots reject
path substitution, symlinks, hardlinks, identity collisions, and changed bytes
before any scientific work.

Correction3 integration binds worker-produced static and synthetic evidence; it
does not execute the absent Todo13 evaluator, use real data, use labels or
composition priors, or make a benchmark/10–20mer claim. The bound mutation
projection is 680 policy/semantic rows plus 163 authority rows (843 total),
not an integrator reimplementation. The validator uses the exact policy
bootstrap contract (500 Mersenne-Twister seeds, Type 7 lower quantile and 32
exhaustive sign masks), descriptor snapshots, final revalidation, and the
unchanged GCV/`n_eff` authority boundary. `n_eff` remains an authority-bound
policy statement rather than a newly calibrated physical quantity.

Correction4 is only publication-provenance remediation. It regenerates the
correction3 static/synthetic bytes twice and publishes them once with
descriptor-relative, exclusive no-replace creation. It introduces no physical
calibration, threshold, model, T8/T11/T12, GCV, `n_eff`, label, benchmark,
application, or Todo change; the live configuration remains a candidate and
Todo13 stays absent and blocked pending parent acceptance, a fresh Oracle, and
GateClosure.

Correction5 is a narrow administrative terminalization caveat. The six
correction4 canonical bytes remain valid by reference and are not republished;
the correction4 claim is non-authoritative because its predecessor hash and
finalizer/Boulder closure were defective. Correction5 records the failure and
cleans only the captured residue. It adds no physical calibration or science,
and changes no configuration, threshold, GCV, `n_eff`, T8/T11/T12, label,
benchmark, application, or Todo behavior. Todo13 remains blocked pending parent
acceptance, a fresh Oracle PASS, and reviewer-owned GateClosure.

Correction6 is an administrative publication boundary, not a calibration
change. Correction5 cleanup remains valid, but Correction5 is non-authoritative
because its evidence closure omitted the cleanup receipt and its terminal replay
omitted all six required per-path canonical bindings. Correction6 references the
existing Correction4 canonical bytes without republishing them. The parent owns
the external checkpoint and performs atomic no-replace directory publication;
publication is established only by the parent's external receipt after the
staged 19-file root and both manifest namespaces have been verified.

Correction6 never authorizes Todo13 itself. Todo13 remains blocked pending
parent acceptance of the external receipt, a fresh independent Oracle PASS, and
reviewer-owned `GateClosure`. No policy, configuration, calibration, threshold,
GCV, `n_eff`, T8, T11, T12, label, benchmark, application, or Todo behavior
changes.

## Reconstruction versus calibration (2026-09-16)

The user-approved Julia reconstruction has its own explicit numerical settings
in `config/unit_assignment_reconstructed.toml`; see the config and unit-assignment
references. The backward residual half-plane descriptor is a new definition,
not a recovered calibration or an inferred original formula. Its 9×9 grid and
normalization are fixed before comparison. The mold/Fisher settings port the
existing method; no benchmark labels, expected counts, or frozen predictions
are used to calibrate these quantities. The counting calibration is unchanged.

The September 20 opt-in `unit_assignment_transverse_fisher.toml` tests physical-u
instead of physical-t reflection when scoring the empirical Fisher mold. Its
`fisher_layout` is a coordinate convention, not a calibrated length, noise
estimate or fitted threshold. Patch sampling, normalization, disk support,
Fisher training, seeds and voting remain fixed. Neither reflection symmetry nor
a synthetic test alone establishes improved chemical recognition; the complete
external benchmark comparison is required. Original numerical settings are unchanged.
That fixed comparison regresses to 669/870 correct and 27/145 exact, versus
671/870 and 28/145 for the control at identical coverage; the transverse
variant is not retained. No recalibration or threshold adjustment follows.

The separate September 20 matched-residual comparison tests
`patch_residual_filter = "smooth_residual"`: `S(data - model)` rather than
`S(data) - model`. This makes an exact model image yield zero residual and
retains an injected residual as `S(signal)`, not as an unsmoothed signal.
It changes no smoothing radius, physical calibration, geometry, selected count,
normalization, classifier setting or missing-value treatment. It removes the
deterministic mismatch between smoothed data and unsmoothed model, not an
arbitrary background or transverse gradient. Synthetic consistency alone does
not establish better recognition; the full external comparison is required.
Both older assignment configs now state `smooth_data_only` explicitly, changing
their file hashes but not their effective settings.
The completed paired grade gives **672/849 (79.2%), 29 exact**, versus
671/849 and 28 for the control. This small development-benchmark gain supports
keeping the opt-in candidate, not a claim of calibrated chemistry or independent
generalization; historical 677/854 and 36 exact remains better. No threshold or
other calibration adjustment follows.

The subsequent two-descriptor comparison changes neither calibration nor
residual extraction. Its first moment uses the known transverse coordinate,
normalized by the existing patch half-width. Its affine-residual asymmetry uses
the orthogonal least-squares projection away from `[1,t,u]` on the full 9×9 grid.
This removes the affine component of molecular signal too: only the orthogonal
remainder is retained. It is not an estimated physical background, an acquisition
correction or a guarantee of preserving chemistry. Synthetic plane invariance
and agreement with an independent QR projection establish the algebra only.
All geometry, physical scales, missingness rules, classifier settings and voting
stay fixed; only the one descriptor supplied to GMM changes. Both variants are
declared before a single full-cohort comparison and external grade.
The completed comparison gives **672/849, 28 exact** for the first moment and
**673/849 (79.3%), 33 exact** for affine residuals, versus **672/849, 29 exact**
for the replayed control. Retain the affine config as a working candidate only;
history still has four more correct positions and three more exact chains.
The sixteen changed decisions all end at zero vote margin, and the benchmark
has been reused extensively. Neither this small gain nor the projection algebra
establishes a physical background correction or calibrated chemical confidence.

The independent signed-CC and affine-Fisher comparisons introduce no new
molecular calibration. Signed CC uses the existing physical template cost
ordering, not an experimental class or an oracle sign flip. Affine Fisher removes
the geometric `[1,t,u]` subspace from its 197-pixel disk at training and scoring;
the original unprojected center amplitude still names the learned groups.
Its explicit `fisher_projection_zero_l1 = 1e-12` numerical floor marks a vanished
projected signal unavailable, without missing-pixel imputation. This also removes
affine molecular signal. Templates, physical scales, mixture composition freedom
and count selection are unchanged. The completed September 21 comparison is
negative: signed CC gives **668/849 (78.7%), 32 exact** and affine Fisher
**672/849 (79.2%), 32 exact**, against the exactly replayed **673/849, 33 exact**
control, all at 849/870 coverage. Neither replaces the affine-descriptor working
candidate. Some incorrect new decisions have voting margin one; this is not
evidence of calibrated chemical confidence. No post-grade sign flip, threshold,
seed or hybrid is fitted to the benchmark.


The September 21 numerical comparison also introduces no molecular calibration.
`fisher_score_center` explicitly chooses the historical near-zero centered mean
or the original opposite-fold training mean; it does not impose equal cluster
populations. `gmm_final_covariance` chooses historical `ridge` or analytical
`ledoit_wolf` shrinkage on the final hard covariance only. The explicit
`gmm_covariance_ridge = 1e-6` is the unchanged numerical floor, not a measured
noise variance. No direction, sign, composition, threshold or seed is calibrated
from truth. Whole-scan withdrawal diagnoses sensitivity to training data, not
held-out recognition; Fisher features are not retrained within those halves.
Neither covariance conditioning nor stability calibrates chemical confidence.
The completed comparison improves **673 → 675 correct / 870**, with unchanged
**33/145 exact** and 849/870 coverage, using training-mean centering. This is the
latest opt-in working candidate, still below historical 677/36. Shrinkage gives
**665 correct / 34 exact** and is not retained as the replacement: its median
full-cohort covariance condition number falls from about 23,138 to 191, but
recognition loses eight correct positions. Whole-scan withdrawal changes 370/893
GMM decisions for ridge versus 350/893 for shrinkage: substantial sensitivity
remains despite modest improvement of that diagnostic. Better conditioning and
stability are not a substitute for the external recognition grade. No combined
variant, threshold, seed or parameter is tuned afterward.

### Patch support and covariance-volume score (2026-09-21)

`assignment_patch_support = "complete_disk_symmetric"` changes availability, not
the finite-pixel median/sample-standard-deviation patch normalization. The Fisher
disk must be completely observed. The affine backward descriptor additionally
uses reflection-closed observed support around a complete central disk; discarded
corner partners are not estimated or filled. Complete patches retain identical
arithmetic. Newly usable rows can change learned Fisher/GMM features throughout
the cohort, so this is a full-cohort experiment, not a local repair of only the
previously missing predictions. Missing disk pixels still require abstention.

`gmm_final_score = "gaussian_density"` adds the covariance-volume normalization
only at final assignment. It does not refit means, covariances or mixture weights,
nor change the two distance-based hard self-training iterations. The raw-amplitude
cluster naming convention and hard seed-vote aggregation remain, but their
resulting assignments may change. Dependent lobes, per-file standardization and
learned clusters mean that neither normalized density nor vote margin establishes
chemical confidence. No benchmark-derived sign, class ratio or threshold enters
either candidate. These are two independent, uncombined alternatives to 675/33;
no normalization or scan-weighting experiment is included.

## Opt-in label-free exploration (2026-09-18)

`config/label_free_exploration.toml` is diagnostic configuration, not a new
molecular calibration. Its empty `[model]`, `[selection]` and `[preprocessing]`
sections do not replace the explicit original molecule config. Native support,
amplitude/geometry bounds, and selection settings remain fixed. Unknown adaptive
inputs still require their original selected summary.

`[counting_variable_projection]` supplies bounded optimizer budgets, linear
solver/KKT/SVD tolerances, model-mapping tolerances and the explicit native
elliptical iteration budget. They are numerical controls, not fitted physical
quantities or externally graded hyperparameters. `[representation]` supplies
arithmetic/serialization comparison tolerances only. `[acquisition_noise]`
declares a small lag window, background exclusion, sample sufficiency checks and
correlation/block diagnostics. These are exploratory validity rules, not
calibrated uncertainty guarantees or replacement production thresholds.

An empirical supervised classifier score is not an information-theoretic upper
bound on the signal. Conversely, stable label-free clusters, lower residuals,
view agreement or reproducible votes do not prove chemical identity. Parameters
must not be selected by repeatedly reading external benchmark grades.


The implemented auxiliary `Current` diagnostic reports direct coupled channel
scales/correlations only; the optional nuisance projection was not implemented.
It does not calibrate chemistry, modify fit weights or supply an independent
likelihood. A fixed Fisher replay exports new fold weights under unchanged
reconstruction settings and checks saved scores; it is neither chemical
calibration nor recovery of unique historical coefficients. Raw-QC plates also
retain and flag nonfinite samples that native preprocessing imputes.


## Masked-background synthetic prototype (2026-09-18)

The user approved starting the first alternative after the completed label-free
exploration. The initial deliverable is an opt-in prototype and synthetic
signal-preservation evidence, not a new production calibration or a real-data
result. Controls are explicit in `config/masked_robust_preprocessing.toml`.

Native median imputation precedes plane/row flattening. Restoring missing masks
afterward cannot remove the imputation's influence on valid pixels. The prototype
instead estimates the background from finite supplied-background pixels only.
It does not remove a transverse component from each molecular patch by default.
A correct foreground exclusion can protect an injected localized signal, but
this conditional property does not prove that a real exclusion is correct.
Gaussian tails, unrecognized molecules or instrumental structure can remain in
the chosen background.

Sequential plane fitting followed by row medians also has a specific limitation:
unequal x support across rows can make row offsets bias the fitted x slope.
Separate y tilt and arbitrary row offsets are not uniquely identifiable.
Signal-preservation checks must therefore include absolute background error,
not only the difference between signal-present and signal-absent images.

The predeclared synthetic comparison separates four operations: native
imputation/reference flattening; finite-only OLS without foreground exclusion;
finite guarded OLS; and finite guarded Huber. No control is selected using
benchmark outcomes. Coverage loss and unsupported/nonconverged cases remain
visible alongside same-support errors. Huber's fixed initial-residual MAD and
its numerical floor are not calibrated STM noise, and improved agreement does
not establish chemical identity or correct counts.

A possible later real-data comparison would retain the original four diagnostic
cases and saved geometry/support. Such geometry already depends on the original
preprocessing; freezing it makes comparisons controlled, not independent. No
new real-data execution or cluster submission is part of this synthetic-first
step. Registration search boundaries, GCV and production defaults stay unchanged.


Synthetic verification is complete: 441 engine assertions, 1,116 signal-study
assertions and 12,641 independent saved-table checks pass. No Huber control was
retuned. The tests demonstrate conditional preservation with a correct supplied
mask, but also sequential x-slope bias, outlier-sensitive initial scale, signal
attenuation with a leaky exclusion and explicit loss of unsupported foreground
coverage. Huber is not uniformly better than guarded OLS. There is no new STM
calibration, real-data result, calibrated uncertainty or production default.


The authorized four-scan real-data pilot keeps the synthetic prototype's
controls unchanged. It measures background changes, observed coverage and
forward/backward agreement, not molecular accuracy. Comparisons use fixed common
observed supports over the original lag grid; dropping an unavailable method or
silently comparing different pixel populations is not allowed. A failed method
can leave the all-method comparison unavailable while native-versus-other-method
comparisons remain separately reportable.

The level convention is one median per view on the same observable guarded
background support, held fixed over all lags. It is not truth alignment or gain
calibration. Better post-fit background statistics or image agreement would not
show that molecular contrast was preserved on real images. The exclusion and
fit support still depend on the old native preprocessing and saved geometry.
No new count, chemical fit, noise covariance, independent-sample estimate,
benchmark grade or production-default change is part of this pilot.
