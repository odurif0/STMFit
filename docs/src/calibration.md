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
