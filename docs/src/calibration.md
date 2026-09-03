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
| `spatial correlation range` | 2D isotropic autocorrelation on the full preprocessed image; range = first lag where ρ(h) drops to 1/e |

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
| `selection_criterion` | `gcv` (valid under strong spatial correlation — see §Effective sample size) |
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

The STM residual field is strongly spatially correlated (ρ ≈ 0.9–0.95 at
lag 1; autocorrelation range 17–100 px). The fit window (~10 px) is **smaller**
than this correlation range, so the number of independent observations inside
the window is not meaningfully estimable: any `n_eff` (the `n÷9` heuristic, a
Durbin–Watson AR(1) estimate, or a variogram estimate) is an arbitrary choice
that changes the absolute scale of BIC/AICc by orders of magnitude.

BIC and AICc assume `n` independent observations. Because `n_eff` is undefined
here, **their absolute values are not interpretable** as model-selection scores;
they are retained only as shape diagnostics (how their *ranking* changes across
N, not their magnitudes).

**GCV** (`RSS·n/(n−p)²`) sidesteps the issue entirely: it is the analytical
leave-one-out cross-validation error of a linear smoother and does not require
choosing an `n_eff`. It is therefore the canonical practical criterion for
`N_selected`, and the selection by GCV is robust to the threshold and
reproducible across runs. BIC/AICc remain secondary diagnostics whose absolute
scale must not be trusted.

## Calibrating a new molecule

1. Pick **one** clean, well-resolved scan of the new molecule.
2. Run `test/measure_calibration.jl <scan>` → produces `<scan>_calibration.toml`.
3. Inspect the measured values (especially FWHM and spacing — sanity-check
   against the visible structure).
4. Run the batch with the auto-calibrated TOML; spot-check N_selected on a few
   files visually.
5. If a parameter looks off (e.g. FWHM under-estimated on a noisy scan), measure
   on 2–3 scans and take the median.

The correlation range and noise level are **instrument + preprocessing**
properties, not molecule properties: they are the same for any molecule on the
same STM with the same flatten/smooth settings. Only the molecule-specific
quantities (FWHM, spacing, n_max) need re-measurement.

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
