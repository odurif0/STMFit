# Pipeline Architecture

## Data Flow

```
SXM File (.sxm)
  │
  ├─→ STMMolecularFit: extract_slide()
  │     └─→ 1D profile along chain axis
  │           └─→ GaussianFit1D: fit_slide()
  │                 └─→ Independent QC/reference count
  │
  └─→ GaussianFit2D: chain_gaussian_sweep()
        │
        ├─→ Circular sweep (σ∥ = σ⟂)
        │     └─→ Deterministic 2D-only initialization + NLopt + LsqFit
        │
        ├─→ circ→ell LsqFit refinement (per N)
        │     └─→ Warm-start from circular, local only
        │
        └─→ Model selection
              └─→ configured criterion (default GCV), plus QC metrics
```

## Component Roles

### STMFitCore.jl
Shared mathematical utilities:
- `effective_spacing_min(spacing_min, spacing_max, sigma_max, max_overlap)`
- `kappa_penalty(κ; kappa_max, weight)` — condition number penalty
- `adjacent_kappa_max(deltas, sigmas)` — max adjacent condition number
- `endpoint_overrun(ts, tmin, tmax)` — support boundary check

### GaussianFit1D.jl
1D multi-Gaussian fitting on the axial slide profile:
- Sweeps N=2..max using NLopt + LsqFit
- Ghost peak filter (rejects models with ≥2 unconstrained edge peaks)
- `sBIC` (Student-t BIC) for model selection
- Produces an independent reference count (`N_1D`) and support length for QC.
- It is not used to initialize the standard 2D circular batch sweep.

### GaussianFit2D.jl
2D chain model with Gaussian lobes along a PCA-derived axis:
- `_weighted_roi_axis()` — intensity-weighted PCA via SVD
- `_active_t_support()` — adaptive support detection from axial profile
- `_chain_fit_data()` — extracts tube around axis, fits support bounds
- `_deterministic_chain_seed()` — autonomous 2D circular initialization from raw binned axial signal
- `_decode_chain()` — converts optimizer params → MolecularFeature list
- `_chain_model_values()` — evaluates 2D Gaussian model at grid points
- `_fit_chain_n()` — single-N optimizer (NLopt global + LsqFit local)
- `chain_gaussian_sweep()` — bidirectional N sweep with early stopping

### STMMolecularFit.jl
Orchestration and I/O:
- SXM file reading (Nanomics format)
- Slide profile extraction and arc-length correction
- Batch orchestration and 1D/2D QC comparison
- Plot generation and output file management

## Optimization Strategy

### Circular Model (anchor model)
```
deterministic raw-2D seed → NLopt (GN_DIRECT_L) → LsqFit (LM)
```
The circular model is initialized independently from the 1D fit. Candidate
centres are derived from the raw 2D axial profile (uniform, weighted quantile,
raw local maxima, and edge-aware seeds); the selected seed initializes the
single global/local optimization path for that N.

### circ→ell Refinement (elliptical model)
```
LsqFit (LM, local, 50 iter) — warm-started from circular solution
```
Finds the true elliptical minimum without global exploration.
This replaces the elliptical NLopt sweep entirely (see Research Journal §7, §11).

## Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `spacing_min_nm` | 0.35 | Minimum inter-lobe spacing |
| `spacing_max_nm` | 0.75 | Maximum inter-lobe spacing |
| `max_overlap` | 0.60 | Maximum lobe overlap fraction |
| `sigma_parallel_min_nm` | 0.191 | Minimum axial sigma (FWHM 0.45) |
| `sigma_parallel_max_nm` | 0.509 | Maximum axial sigma (FWHM 1.20) |
| `sigma_perp_min_nm` | 0.10 | Minimum perpendicular sigma |
| `sigma_perp_max_nm` | 0.55 | Maximum perpendicular sigma |
| `fit_width_nm` | 0.16 | Tube half-width around the molecular axis |
| `support_noise_k` | 2.5 | Threshold multiplier: baseline + k·noise |
| `support_padding_nm` | 0.25 | Chitosan calibrated support edge padding |
| `selection_criterion` | gcv | Primary criterion: gcv, bic, aicc, or cv |
| `cv_method` | gcv | Analytical GCV by default; kfold is slower |
| `cv_folds` | 5 | Cross-validation folds when `cv_method="kfold"` |
| `kappa_max` | 10.0 | Chitosan calibrated condition-number penalty threshold |
