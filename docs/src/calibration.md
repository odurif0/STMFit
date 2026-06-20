# Calibration: deriving parameters objectively

The pipeline has ~25 calibration parameters. Most can be **measured** from a
single clean scan rather than hand-tuned, which makes the analysis generalizable
to a new molecule on the same STM. This page documents which parameters are
objective, which are principled choices, and which remain free.

## Auto-calibration

```bash
julia --project=. test/measure_calibration.jl path/to/clean_scan.sxm
```

This measures the objectivable quantities and emits a ready-to-use TOML. Validate
it on the benchmark before adopting it:

```bash
julia --project=. test/batch_full.jl 48 --config chitosan_auto.toml ...
```

## Parameter classification

### Measured from a single scan [objective]

| Parameter | Measurement method |
|---|---|
| `noise σ` | 1.4826·MAD of (raw − smoothed) high-frequency band, via the standard preprocessing pipeline |
| `pixel resolution` | `range_nm / width` (from the SXM header) |
| `FWHM range [lo, hi]` | Detect peaks in the chain-axis profile (weighted PCA → bright-pixel strip), fit half-max width per peak, take [5%, 95%] quantiles |
| `repeat spacing` | Median peak-to-peak distance along the chain axis |
| `spatial correlation range` | 2D isotropic autocorrelation on the full preprocessed image; range = first lag where ρ(h) drops to 1/e |

### Derived from a physical/numerical principle [principled, one fixed choice]

| Parameter | Derivation |
|---|---|
| `sigma_parallel_*` | `FWHM / 2.355` (Gaussian width relation) |
| `spacing_min/max` | `±30%` around the measured repeat spacing |
| `fit_width_nm` | `= σ_min` (tube half-width ≈ narrowest lobe half-width) |
| `support_min_length_nm` | `3 × spacing` (at least 3 repeats to call it a chain) |
| `n_max` | `axis_length / spacing_min + 2` (generous cap from chain geometry) |
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

The STM produces residuals with strong spatial correlation (ρ ≈ 0.9–0.95 at
lag 1; autocorrelation range 17–100 px). The fit window (~10 px) is **smaller**
than the correlation range, so the number of statistically independent points
in the window is effectively zero:

> n_eff(window) = window_area / (π · range²) ≈ 100 / (π · 23²) ≈ 0.06

This means **BIC and AICc — which assume n independent observations — are not
well-defined** in this configuration. Their `n_eff` is an arbitrary constant
(the current `n÷9` heuristic or any DW/variogram estimate), and the absolute
values of BIC/AICc are therefore unreliable diagnostics.

**GCV** (`RSS·n/(n−p)²`) does not assume independence: it is the analytical
leave-one-out cross-validation error of a linear smoother, valid under spatial
correlation (smooth-spline theory). This is why the selection by GCV is robust
to the threshold and reproducible across runs. **GCV is the canonical selection
criterion**; BIC/AICc are retained only as secondary diagnostics whose absolute
scale is not interpretable.

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
