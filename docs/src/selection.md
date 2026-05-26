# Model Selection

How STMFit chooses the optimal number of lobes (N) from the sweep results.

## Selection Hierarchy

```
Level 1: configured criterion per N (default: GCV)
         score_eff(N) = min(score_circ(N), score_ell_refined(N))
         └─ Circular model is nested within elliptical
         └─ Circular score is a robust fallback if elliptical refinement fails

Level 2: model-specific bests
         └─ N_ell = best valid refined elliptical model
         └─ N_circ = best valid circular model
         └─ N_eff = best effective min(circ, ell) model

Level 3: Final result
         └─ Report N_ell, N_circ, N_eff, residual/support QC, plots
```

## Why GCV by default?

GCV is the default because it is analytical, fast, and does not require
refitting folds. BIC is still computed and remains useful for diagnostics and
legacy comparisons, but it is not the default batch selection criterion.

- `n_eff` (effective sample size) is estimated as `length(zfit) ÷ 9`,
  which may not perfectly capture spatial correlation in STM images.
- BIC assumes all parameters contribute equally, but extra sigma parameters
  can absorb noise without improving predictive accuracy.
- On ambiguous STM images, BIC can marginally prefer over-fit models.

The `selection_criterion` field can be set to `"gcv"`, `"bic"`, `"aicc"`, or
`"cv"`. The chitosan calibration uses `"gcv"` with `cv_method="gcv"`.

## Why min(ell, circ)?

The circular model (σ∥ = σ⟂ per lobe) is a **nested special case** of
the elliptical model. Therefore:

- For any N, there exists an elliptical solution with objective ≤ circular
  objective under the same residual criterion
  (set σ∥ = σ⟂ and you get the circular solution exactly).

If the refined elliptical score is worse than circular at the same N, the
elliptical local refinement did not improve that model. Using `min()` guards
against this failure without introducing any new parameters.

## Effective output columns

- `N_ell`: best valid refined elliptical 2D model. This is the primary model
  count when assessing the elliptical Gaussian chain.
- `N_circ`: best valid circular 2D model.
- `N_eff`: hybrid/effective best using `min(score_circ(N), score_ell(N))` per N.
- `N_1D`: independent 1D slide-profile count used for QC, not to initialize the
  standard 2D circular batch sweep.

## Ambiguity diagnostics

Batch summaries also report close-second-best GCV diagnostics.  They do **not**
change the selected `N`; they only flag cases where the best model and the
second-best distinct `N` are close in relative GCV (default: ΔGCV/GCV ≤ 5%).

- `ambiguous_ell`, `runnerup_N_ell`, `delta_GCV_ell`, `delta_GCV_rel_ell`
  describe the refined elliptical sweep.
- `ambiguous_eff`, `runnerup_N_eff`, `delta_GCV_eff`, `delta_GCV_rel_eff`
  describe the effective `min(ell,circ)` selection.

These columns are intended for QC/visual review of support-sensitive or weakly
identified files, not as an expected-N prior.

For ambiguous files, the best-plot title also includes a warning such as
`ambiguous ell GCV: selected N=5; second best N=6 (ΔGCV=2.5%)`, so close
alternatives remain visible during plot review.

## circ→ell LsqFit Refinement

**All elliptical fitting is now done via circ→ell LsqFit refinement.**
The NLopt global optimizer is intentionally excluded — it always diverges
from the isotropic start in 33D parameter space (see Research Journal §7).

For each N fitted by the circular sweep:
1. Circular solution provides positions, amplitudes, and isotropic σ
2. Params expanded to elliptical format (σ∥ = σ⟂ = σ_circ initially)
3. LsqFit only (skip_global=true, max_iter=50) locally refines sigmas
4. Scores are computed for the configured criterion; the refined elliptical and
   circular results are compared per N.

**Why this works**: The circular solution is near the global elliptical minimum.
LsqFit follows the gradient to the nearest local minimum without exploring
the wider parameter space where NLopt gets lost. The refinement consistently
improves BIC (3-24× better than NLopt elliptical in tests).

## Effective Sample Size

`n_eff = max(10, length(zfit) ÷ 9)`

The divisor ÷9 approximates the spatial correlation area: a 3×3 pixel
block ≈ 1 independent observation. Typical STM images have FWHM ~0.5 nm
≈ 25 pixels, giving a correlation area of ~500 px². The ÷9 factor is
conservative (larger n_eff → larger BIC penalty → favors simpler models).
