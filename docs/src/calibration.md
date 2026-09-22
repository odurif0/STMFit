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
opt-in working candidate at that stage, still below historical 677/36. Shrinkage gives
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

The completed comparison gives **676/870 correct and 34/145 exact** for support,
versus 675/33: three recovered correct predictions are offset by two regressions
among previously available rows. Coverage rises **849 → 852/870**, while
classified accuracy falls **79.5% → 79.3%**. The gain on the fixed denominator,
exact chains and coverage supports retaining it as the latest opt-in candidate,
not as a calibrated-confidence improvement or a new champion. Historical 677/36
still leads. The Gaussian-volume score gives **666/870 and 32/145**, with all
19 changed final decisions at zero vote margin, and is rejected. No threshold,
class count, per-file choice or parameter is adjusted after these results.

### Separate complete training from partial prediction (2026-09-21)

The sole new `[selection] assignment_training_support = "complete_patches"`
mode restricts Fisher/GMM learning by observed patch completeness, never by
predicted type or external correctness. Fisher requires its full forward square;
GMM requires all three patch families it consumes. The support candidate's
disk/symmetry rule still governs whether partial patches can be scored. This
does not change extraction, patch normalization, physical calibration, selected
N, DFT molds, vote thresholds or k-means.

The GMM's per-file feature-normalization moments, fits, weights and raw-amplitude
group naming are estimated from complete rows only, then frozen. Moments still
use finite values separately per feature; there is no robust estimator, global
fallback or scan reweighting. A file with no finite complete training value for
a feature cannot provide that view. Opposite-lobe-parity Fisher scoring remains;
this is not a new held-out-scan validation design. Complete training rows can
still be noisy and fewer rows can destabilize learning, so improvement is a
hypothesis, not a guarantee. The full 146-file output and all regressions must
be checked before the separate 145-file external grade is considered.

The completed test is **negative: 674/870 correct, 33/145 exact**, versus
676/34, with coverage falling **852 → 850/870**. All four fitted lobes of
`240818_019` have partial squares, so that file has no complete-row moments;
the prescribed abstention loses its two previously usable predictions. Elsewhere,
one gain is offset by one loss, including loss of exactness for `240818_020`.
This is evidence against this complete-training variant as a replacement, not
permission to add a benchmark-directed fallback. Keep the 676/34 support
candidate; normalization rules, counts and the frozen application stay unchanged.

### Robust GMM feature normalization (2026-09-21)

The newly authorized comparison changes only GMM per-file feature moments from
mean/sample-std to median/IQR, before pairwise products. This is not physical
pixel calibration and does not affect Fisher or k-means. The support candidate
at 676/34 is the control; complete-case training is not reused. All finite
admissible values remain eligible, including usable partial patches.

`[preprocessing] gmm_feature_normalization` selects the explicit formula, and
`gmm_scale_fallback = 1.0` records the degenerate-scale guard in every native
config. Q25/Q75 use Type 7 interpolation, with no normal-consistency rescaling.
An empty feature support stays unavailable; a zero/nonfinite IQR uses one.
Positive small IQRs are not floored. Median/IQR resists changes in extreme
observations when estimating moments, but does not bound their transformed
values; short scans or interactions can therefore become less stable, not more.
The completed comparison is negative: **667/870 correct, 667/852 (78.3%),
10/145 exact**, versus support control 676/34 at unchanged coverage. It loses
24 exact chains and gains none. No IQR fallback occurs in 1,168 scan/feature
pairs, but the expanded maximum grows from 5.44 to 335.47. This observed tail
amplification does not, by itself, establish the cause of the recognition loss.
Retain the mean/std support candidate; no clipping or alternative quantiles are
introduced after the grade.

All choices precede grading. No expected count, sequence, composition, benchmark
label, class balancing, scan weighting or benchmark-tuned fallback enters this
experiment. Counting, DFT and the frozen unknown-chain application are untouched.

### Equal-scan GMM weighting (2026-09-21)

The next authorized comparison leaves physical calibration and normalization
unchanged. It tests whether giving each scan equal total training influence
helps the 676/34 support candidate. `[selection] gmm_training_weighting` selects
`equal_lobes` (existing behavior) or `equal_scans` (one opt-in candidate).
Weights are proportional to inverse usable training-row count, normalized to
mean one to retain the objective's numerical scale. They affect GMM learning
and physical amplitude naming consistently, not the separate Fisher/k-means
heads or final vote. Within-scan finite feature moments stay unchanged.

This is not class balancing, an expected-N prior, or an effective-sample-size
estimate. Both chemical populations remain free. A scan with few usable lobes
can gain influence even if noisy; equal influence is a hypothesis, not a
guarantee of correct chemistry. The fixed comparison is one control
and one variant on Viper with a 30-minute ceiling. No count refit, bootstrap,
DFT change, unknown-chain rerun or post-grade adjustment is included.

The completed result is mixed: **675/870 correct, 675/852 (79.2%), 36/145 exact**,
versus the replayed support control **676/852, 34 exact**. It gains two exact
chains but loses one correct position; all seven changed decisions end at zero
vote margin. Retain the 676/34 support candidate as the primary reference, with
historical 677/36 still unexceeded. This result does not establish that short
or noisy scans caused the loss, nor justify a new physical calibration or
post-grade weight adjustment. The reused benchmark is development evidence.

### Continuous GMM seed aggregation (2026-09-21)

One authorized comparison retains the GMM's normalized per-seed scores instead
of replacing each by an argmax decision before averaging seeds. The explicit
`selection.gmm_seed_aggregation` modes are `hard_vote` (legacy) and
`mean_membership` (one opt-in candidate). Only aggregation changes, not the
learned groups, their amplitude-based physical names or the final threshold.
There is no temperature, label-fitted probability calibration or mixture-size
constraint. Scores may already be saturated, and smooth votes need not improve
chemical correctness. Their confidence remains an uncalibrated vote margin.

Physical calibration, selected N, missing-data handling, eight-decimal output
precision and the two-head vote stay fixed. Scope: support control plus one
candidate, one Viper allocation capped at 30 minutes; no scan weighting,
resampling, new naming rule, count refit, DFT change, unknown-chain rerun or
post-grade adjustment.

The completed result is **negative: 671/870 correct, 671/852 (78.8%), 26/145
exact**, versus the exactly replayed 676/34 support control at the same coverage.
All 24 changed decisions leave old exact ties (GMM=1, k-means=0) and become 0
as the continuous GMM score drops below 1. The GMM still emits class 1 for these
24 rows; the loss arises in the unchanged final mean-of-two-heads decision.
One exact chain is gained and nine lost. Reject this variant, retain the 676/34
support config and leave calibration, precision and thresholds unchanged.
Historical 677/36 remains the target. This is evidence about this specific
uncalibrated fusion, not a general rejection of continuous scores or independent
validation on unknown chains.

### Whole-scan GMM bootstrap (2026-09-21)

One fixed comparison resamples **whole scans**, not individual lobes or chemical
classes. Twenty replicates each draw S usable scans with replacement, using
seeds 0–19; all admissible rows of a drawn scan travel together. The ten existing
initialization seeds are reused within each replicate. Row multiplicity enters
learning and the unchanged raw-amplitude naming rule, not an imposed class prior.
Per-scan scaling and physical calibration remain unchanged. Unsampled scans
still receive predictions, so this is not held-out validation.

The explicit settings are `selection.gmm_resampling="whole_scans"`,
`gmm_bootstrap_replicates=20` and `gmm_bootstrap_seed=0` in the separately named
opt-in config. Earlier configs retain `none`, `1`, `0`. The candidate averages
hard seed votes within each replicate and then valid replicate means equally;
it does not select a best seed or tune confidence. Scores remain uncalibrated
frequencies. Final threshold, precision, N, unavailable-row handling, k-means,
Fisher, CC and all upstream observations remain fixed.

Scope is one support control and one candidate, one Viper job limited to one
hour, four requested CPUs and 16 GB. No retry, combined variant, new naming
rule, count refit, DFT change, unknown-chain rerun or post-grade adjustment.
Job **11925188** completes **0:0 in 7m17s**. All 200 seed fits are accepted;
bagging gives **665/870 correct, 665/852 (78.1%), 10/145 exact**, versus the
exactly replayed **676/852, 34 exact** support control. Coverage is unchanged;
24 exact chains are lost and none gained. Of 62 final 1→0 changes, 57 leave
old GMM=1 / k-means=0 ties under the unchanged fusion. Reject the variant;
retain the **676/34 support candidate**, still below historical 677/36. No
post-grade change to fusion, seeds, replicate count, naming, precision or
threshold follows. This negative result concerns this specified ensemble and
fusion, not all bagging methods; the reused grade is not independent validation.

### Tied GMM covariance throughout learning (2026-09-22)

The user authorizes one comparison of the 676/34 support control with a GMM
whose two components share covariance throughout learning. This is not the
previously rejected final-only Ledoit-Wolf estimate or final Gaussian-volume
score. The same features, per-scan scaling, interactions, ten seeds, two hard
updates, physical group-naming rule and final vote stay fixed. N, Fisher, CC,
support, acquisition preprocessing and unknown25 are not changed.

`model.gmm_covariance_structure="tied"` pools within-component scatter using
the learned responsibilities in every EM update and hard indicators during
initialization/self-training; divide by the number of usable rows and add the
unchanged ridge. Means and mixture masses stay free, with no class-count prior.
The covariance has 666 independent entries in the 36-dimensional view instead
of 1,332 across two separate matrices. This reduces model freedom, not a claim
of calibrated noise or chemical certainty. Better numerical conditioning alone
would not establish better recognition.

One Viper job is capped at 30 minutes, four requested CPUs and 16 GB, without
automatic retry. Exact control replay and all-key/N/upstream/vote checks precede
external grading. Labels, expected N, benchmark membership and historical
predictions are not learning or inference inputs. No naming, grouped-Fisher,
fusion, threshold, resampling or post-grade search accompanies this test.
The completed result is **negative: 666/870 correct, 666/852 (78.2%), 6/145
exact**, versus the exactly replayed 676/34 support control at the same
coverage. No exact chain is gained and 28 are lost. Reducing covariance freedom
does not improve recognition in this fixed representation/fusion. Learned
smaller-component masses are about 13.2–14.4%, not an imposed class proportion;
final class-1 emissions decrease from 223 to 156 across all 900 rows. These
post-hoc observations do not authorize a composition prior, new naming rule or
threshold adjustment. Reject this variant and retain **676/34**; historical
677/36 remains the target. This reused grade is not independent validation.

### Whole-scan Fisher and relative GMM naming (2026-09-22)

Two independent tests retain the 676/34 support control and all physical/count
calibration. The first changes only Fisher's cross-fitting unit: entire scans
are kept in one of two deterministic SHA-256-ranked groups, with an explicit
split seed zero. Every fitted transformation, center and naming operation
uses only the opposite group. Invalid rows do not choose the grouping, and
neither chemical labels nor a required number of lobes is used. This tests
sensitivity to shared scan structure; it is not proof that the parity feature
was chemically wrong or that scans are independent molecules.

The second changes only the global names of the already fitted GMM components.
`gmm_cluster_naming="within_scan_z"` centers/scales raw amplitudes within scans
using eligible, feature-valid training rows; the higher global group mean names
class 1. The existing explicit `gmm_scale_fallback` applies to constant or
singleton training scans. Positive scan gains/offsets then cancel for nonzero
scales, but no chemical accuracy or invariance under arbitrary tip changes is
claimed. No class quota, composition prior, per-scan forced group or DFT anchor
is introduced; the underlying higher-amplitude physical convention remains an
assumption. GMM means, covariances, memberships and free masses stay unchanged.

The two candidates are not combined. Labels belong only to the external grade
after fixed outputs; one Viper job is capped at 30 minutes with no retry.
The completed comparison gives **668/870 correct, 668/852 (78.4%), 34/145 exact**
for whole-scan Fisher: eight fewer correct positions than the replayed support
control, with four exact chains gained and four lost. Both disjoint 73-scan
fits converge. All 26 final flips enter existing zero-margin vote ties; no
threshold or split-seed adjustment follows. Relative naming gives **676/870,
676/852 (79.3%), 34/145**, with unchanged fitted parameters, component names,
scores and decisions across all ten seeds. These naming conventions agree on
this cohort, not necessarily under other acquisition conditions. Retain
**676/34**; historical 677/36 remains the target. This negative recognition
result does not establish that parity cross-fitting is independent validation,
and the reused benchmark remains development evidence.

### Fixed alternative mixture learning (2026-09-22)

Two independently tested profiles retain the support control's 8 features,
28 pairwise products, per-scan scaling, amplitude naming, seeds, selected N,
Fisher and final vote. Rank four (`gmm_factor_rank`) and Student df five
(`gmm_student_df`) are fixed before real inference, with no benchmark-driven
choice. The iteration, tolerance, minimum-mass and numerical-guard parameters
are explicit in the [configuration reference](config.md); Gaussian learning
is unchanged. Neither alternative estimates a chemical class proportion.

The factor-analyzer model initializes each component spectrally, using the
mean of the bottom `p-rank` covariance eigenvalues for isotropic residual
noise. EM jointly updates the component center and loadings using augmented
latent moments; separate diagonal noise is floored at the existing ridge.
This follows [Ghahramani and Hinton's MFA equations](https://www.cs.toronto.edu/~hinton/absps/tr-96-1.pdf),
including their component-specific-noise option. The two pure-Mahalanobis
hard steps each retain one such conditional parameter update.

For Student learning, latent precisions are `(df+p)/(df+distance²)`.
Centers use responsibility-times-precision weights, while each scale scatter
is divided by the **responsibility mass**, not the precision-weighted mass;
the existing ridge is then added. See
[Peel and McLachlan](https://people.smp.uq.edu.au/GeoffMcLachlan/pm_sc00.pdf).
The same update is retained during hard reassignment. The stored scale is
not multiplied by `df/(df-2)`. Final scoring deliberately remains the
control's log-mass minus half squared distance, not a Student density.
Thus neither scores nor seed votes are calibrated chemical probabilities.

No preprocessing, DFT calibration, count selection, abstention threshold or
unknown-chain claim changes. The repeatedly reused benchmark can compare
these frozen choices, but cannot provide independent validation.

The completed comparison is negative: factors **620/870 correct, 23/145
exact**, Student **631/870, 25/145**, versus exactly replayed **676/34** support,
at **852/870 coverage** throughout. All ten Student fits converge in 65–86
updates; all factor fits reach 200 updates without satisfying tolerance.
Their smallest noise variance is about 7.59e-5, above the 1e-6 floor, so no
floor saturation or collapsed seed is observed. The capped MFA result does
not characterize a fully converged optimum. False positives increase more
than false negatives decrease; almost every changed final decision enters a
zero-margin vote tie. Retain support without post-grade cap/rank/df,
threshold, naming or fusion adjustment. Source **e2205cb**, job **11935072**;
full losses: `results/factor_student_mixtures_20260922/report.md`.

### Covariance regularization during learning and Student density (2026-09-22)

`unit_assignment_em_shrinkage.toml` applies spherical covariance shrinkage at
initialization, each Gaussian EM M-step and both hard updates. For centered
vectors `y_i`, responsibility weights `a_i = r_i / sum(r)`, define
`S = sum(a_i y_i y_i')`, `T = tr(S)/p I`, and
`b = sum(a_i^2 ||y_i y_i' - S||_F^2)`. The fixed coefficient is
`alpha = clamp(b / ||S-T||_F^2, 0, 1)`, or zero if `S == T`. Use
`(1-alpha) S + alpha T + 1e-6 I`. Hard memberships reduce to the previously
tested equal-weight formula. Means and free mixture masses use unchanged
updates; neither group is forced to contain a specified number of lobes.

This is a **declared fixed-weight plug-in extension**, not a reproduced
implementation or optimality claim for
[Halbe, Bortman and Aladjem's regularized GMM](https://cris.bgu.ac.il/en/publications/regularized-mixture-density-estimation-with-an-analytical-setting/).
The latter motivates regularization during learning. The extension is checked
against direct weighted outer products and the hard-membership limit.
Learned responsibilities, dependent lobes and estimated centering invalidate
any automatic interpretation as independent measurement-error variance. No
`n_eff`, physical noise estimate or chemically calibrated confidence is inferred.
Shrinkage need not increase unpenalized likelihood monotonically; decreases
and bounded-fit convergence are reported without filtering seeds by grade.

`unit_assignment_student_density.toml` retains fixed df `nu=5` and the
[Student latent-precision updates](https://people.smp.uq.edu.au/GeoffMcLachlan/pm_sc00.pdf),
but uses the same log mass-density in the E-step, both hard assignments and
final scoring: `log(pi_c) - log(det(scale_c + guard I))/2 -
(nu+p)/2 * log1p(distance_c^2/nu)`. Only the common fixed-dimension/fixed-df
normalizing constant is omitted. Scale is not converted to covariance.
Hard assignment now maximizes this score, rather than minimizing distance
alone; the parameter update still uses responsibility-times-precision means
and responsibility-mass scale denominators. This is a classification-style
hard-update extension of the soft Student mixture, not pure soft EM throughout.

Both arms preserve raw-amplitude naming, ten hard seed votes, k-means fusion,
thresholds, support and selected N. Internal model memberships and final vote
fractions are not calibrated chemical probabilities. Benchmark labels are
used only in external grading.

Completed outcome: **637/870 correct, 24/145 exact** for all-update shrinkage,
**627/870, 25/145** for coherent Student, versus exactly replayed **676/34**
support, all at **852/870 coverage**. All ten fits per arm meet tolerance:
25 or 43 updates for shrinkage, 65–86 for Student. Shrinkage coefficients
range about 0.02993–0.06040; seeds 1 and 3 each show eleven ordinary-likelihood
decreases, with maximum drop 1.14377. None is excluded or retried. Student
latent precisions range about 0.003045–5.447092. Numerical convergence does
not establish good recognition or a global optimum. Reject these two fixed
versions, retain support and leave the historical 677/36 target open. No
post-grade tuning or unknown-chain claim follows. Source **bb6943b**, job
**11936229**, full report `results/em_shrinkage_student_density_20260922/report.md`.

### Fresh raw-GCV counts with unchanged calibration (2026-09-22)

Repeating the original direct-extractor GCV sweep after symmetric filtering
gives **661/870 correct, 32/145 exact**, versus fresh fixed-N control
**675/34** and saved support **676/34**. Exact counts fall **106 → 101/145**;
coverage falls **852 → 849/870**. No physical bound, GCV/guard, `n_eff`,
threshold, optimizer limit or DFT source changes. This is not the promoted
batch hybrid counting lineage or an exhaustive sweep.

The fresh control changes three decisions despite identical N/settings.
Same-N base/split centers vary by up to 0.00845/0.07637 nm between fresh arms;
their maximum relative GCV differences are 0.115/1.049%. The existing
time-limited global search is a plausible, unproven source of variability,
not evidence for a new calibrated noise model. Whole-cohort learning also
couples predictions across scans. Reject this version without retuning;
the result is end-to-end development evidence, not an isolated causal effect
of N or independent validation. See `results/gcv_reselection_20260922/report.md`.

## Opt-in label-free exploration (2026-09-18)

The September 22 paired-acquisition feasibility experiment preserves the native
Gaussian/split mean shape and all molecular bounds at saved N. It adds
equal/opposite view gain and background-plane terms, fixing mean gain to one
to remove the trivial gain/amplitude ambiguity. `paired_acquisition.toml` bounds
gains to `[0.5,1.5]`, the differential intercept to `[-5,5] nm` and differential
tilts to `[-1,1]`; these are exploratory numerical boxes, not a new physical
calibration. The matched fused and paired continuations share their start and
300-iteration budget. Native mean validity is never relaxed; paired residuals
must additionally pass the unchanged guard in both acquisitions. Stacked GCV
is a heuristic with four added coefficients, not an independence/noise claim.
The first real check is confined to the previously failed scan, without labels.

That check does not clear native validity: Gaussian elliptical mean maximum/
noise is **3.651 fused / 3.518 paired**, above 3.5, with capped nonconverged
continuations. Paired split passes mean validity but fails the forward-view
guard. This is not proof of a wrong N or inadequate physical widths; numerical
and model limitations remain unresolved. No calibration is changed and no
full-cohort recognition grade follows. Report:
`results/paired_acquisition_20260922/report.md`.

The newly authorized `registered_refit_original_support.toml` follow-up freezes
the original unregistered native ROI, axis/origin, tube and axial bounds for
both refit arms. It changes no physical constant and does not expand support
to accommodate saved N. The legacy preprocessing's imputation influence remains
in this original support; only actually observed samples enter either new
objective. The offset is estimated on the observed original ROI. Frame/cache
consistency is checked, with a saved-only replay against prior native geometry
before grading. No new N, noise estimate or independence claim. The three
previous span failures disappear, but `240817_006.sxm` fails the native Gaussian
residual guard after registration; 145/146 scans complete all fits. No partial
grade is produced, and no guard is relaxed. All original supports replay
exactly; retain support **676/34**, with recognition still inconclusive. See
`results/registered_original_support_20260922/report.md`.

The authorized September 22 registered-refit follow-up keeps physical bounds,
registration gates, selected N and assignment settings unchanged. It refits
both Gaussian and split profiles with native optimizers, selecting their
family by GCV at the saved N; it does not reselect N. The observed unsmoothed
mean supplies the fit, with native finite-only ROI statistics and no fabricated
fit pixels. Nearest-observed sampling is initialization-only. Masks cannot undo
earlier native flattening's imputation influence. The zero-shift arm follows
the identical rules and its fits are reused exactly when registration is zero.
ROIs/axes/support may differ for nonzero shifts; lower RSS across these arms is
not itself an improvement in a common objective. Native noise estimates and
the `n_eff` placeholder are unchanged; no independence or convergence claim is
made. Settings: `registered_refit.toml`. The result is **inconclusive**: two
zero-shift controls and one registered scan cannot fit their saved N within
the new observed support under the unchanged conservative parametrization.
143/146 scans complete all fits. This does not prove incorrect counts or worse
recognition; no candidate grade is available. No bound, support parameter or N
is adjusted afterward. Retain support **676/34**; full diagnostic:
`results/registered_refit_20260922/report.md`.

The September 22 acquisition-translation experiment changes neither molecular
calibration nor cached N/geometry. `acquisition_registration.toml` declares all
search and acceptance settings: up to 2 nm / one quarter width in x, four
contiguous bands, signed correlation ≥0.60, distant-peak gap ≥0.01 and 0.08 nm
peak-neighborhood/band-agreement distances. At least 256 common pixels / eight
rows per band and 16 pixels per row are required. These are experimental
identifiability gates, not a noise estimate or chemical validation. Forward is
the fixed reference; only backward observations move before subtracting the
saved Gaussian model. Native flattening stays unchanged; restored raw masks
cannot undo earlier imputation influence on that background. A separate
zero-shift mask control isolates this effect. No N, optimizer, GCV, `n_eff`,
class threshold, composition prior or missing-head fallback changes.
The completed comparison is **negative: 671/870 correct, 26/145 exact**, versus
both controls **676/34**, coverage **850 versus 852/870**. Stronger local
image agreement does not validate this assignment change. Keep saved support,
without post-grade calibration changes; report:
`results/acquisition_registration_20260922/report.md`.

The subsequent September 22 geometric-profile comparison keeps N, physical
calibration, preprocessing and assignment settings fixed. Joint and profiled
LN_BOBYQA arms share one fresh native initialization chosen by valid GCV at the
saved N. Both use the native κ-penalized RSS and native validity checks, unlike
the old raw-RSS-only prototype. Equal evaluation/time ceilings are not equal
computational cost and do not prove convergence; actual stops and costs are
saved. Split skew remains cached. No new noise covariance, effective sample
size, composition prior or class threshold is introduced. The completed external
comparison is **negative: 632/870 correct, 20/145 exact**, versus joint and saved
support **676/34**, at unchanged **852/870 coverage**. All 146 profiled pixel
objectives improve but recognition worsens; 103 searches are evaluation-capped.
Do not replace the saved support geometry or retune calibration from this result.
Settings remain in `geometry_profile_comparison.toml`; full report:
`results/geometry_profile_20260922/report.md`.

The authorized September 22 fixed-geometry amplitude comparison changes no
physical calibration. It solves only the Gaussian amplitudes and tilted plane
inside their native finite bounds, using the existing linear-profile solver.
The literal saved decimal geometry and N stay fixed; the split shape is reused.
GCV keeps every original model parameter and is diagnostic only. A lower RSS
or satisfied KKT condition does not establish better chemical recognition.
The completed comparison demonstrates that distinction: RSS decreases on all
146 scans (median **0.681%**), but recognition falls to **675/870 correct,
32/145 exact**, versus exactly replayed support **676/34**, at unchanged
**852/870 coverage**. N and every geometric field remain identical. Reject
this variant without retuning any bound, calibration, threshold or vote rule;
historical 677/36 remains unexceeded. See the dated [journal](journal.md) and
`results/frozen_amplitude_profile_20260922/report.md` for full diagnostics.

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
