# Pipeline Architecture

## Data Flow

```
SXM File (.sxm)
  │
  ├─→ STMMolecularFit: extract_slide()
  │     └─→ 1D profile along chain axis
  │           └─→ GaussianFit1D: fit_slide()
  │                 └─→ Peak centers & amplitudes (bootstrap for 2D)
  │
  └─→ GaussianFit2D: chain_gaussian_sweep()
        │
        ├─→ Circular sweep (σ∥ = σ⟂)
        │     └─→ Reliable convergence, isotropic gaussians
        │
        ├─→ Elliptical sweep (σ∥, σ⟂ independent)
        │     └─→ NLopt global + LsqFit local
        │
        ├─→ circ→ell LsqFit refinement (per N)
        │     └─→ Warm-start from circular, local only
        │
        └─→ Model selection
              └─→ min(BIC_circ, BIC_ell_refined) + CV tiebreaker
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
- Outputs peak centers and amplitudes as 2D bootstrap

### GaussianFit2D.jl
2D chain model with Gaussian lobes along a PCA-derived axis:
- `_weighted_roi_axis()` — intensity-weighted PCA via SVD
- `_active_t_support()` — adaptive support detection from axial profile
- `_chain_fit_data()` — extracts tube around axis, fits support bounds
- `_decode_chain()` — converts optimizer params → MolecularFeature list
- `_chain_model_values()` — evaluates 2D Gaussian model at grid points
- `_fit_chain_n()` — single-N optimizer (NLopt global + LsqFit local)
- `chain_gaussian_sweep()` — bidirectional N sweep with early stopping

### STMMolecularFit.jl
Orchestration and I/O:
- SXM file reading (Nanomics format)
- Slide profile extraction and arc-length correction
- 1D→2D bridge (bootstrap initialization)
- Plot generation and output file management

## Optimization Strategy

### Circular Model (N params, reliable)
```
NLopt (GN_DIRECT_L, global, 10s) → LsqFit (LM, local, 300 iter)
```
Converges reliably. Fewer parameters, better-conditioned landscape.

### Elliptical Model (N params + N extra sigmas, unstable)
```
NLopt (GN_DIRECT_L, global, 10s) → LsqFit (LM, local, 300 iter)
```
**Known issue**: NLopt diverges from isotropic solution in 33D space.
The global minimum exists but is unreachable by gradient-free optimization.

### circ→ell Refinement (current, optimal)
```
LsqFit (LM, local, 50 iter) — warm-started from circular solution
```
Finds the true elliptical minimum without global exploration.
Replaces the elliptical NLopt sweep entirely (see Research Journal §7, §11).

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
| `cv_folds` | 5 | Cross-validation folds |
| `kappa_max` | 8.0 | Condition number penalty threshold |
