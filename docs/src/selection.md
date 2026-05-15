# Model Selection

How STMFit chooses the optimal number of lobes (N) from the sweep results.

## Selection Hierarchy

```
Level 1: min(BIC_circ(N), BIC_ell_refined(N))
         └─ Circular model is nested within elliptical
         └─ circ_BIC is a mathematically valid lower bound
         └─ When ell_BIC > circ_BIC → convergence failure → use circ

Level 2: CV tiebreaker
         └─ Triggered when CV strongly disagrees with BIC
         └─ Rule A: CV_ratio > 2.0 → CV clearly prefers different N
         └─ Rule B: ΔBIC < 100 AND CV prefers simpler N

Level 3: Final result
         └─ Selected N with best (BIC, CV) combination
         └─ Parameters from the winning model variant (circ or ell)
```

## Why Not Just BIC?

BIC is an asymptotic approximation: `BIC = -2·log(L) + k·log(n)`.

- `n_eff` (effective sample size) is estimated as `length(zfit) ÷ 9`,
  which may not perfectly capture spatial correlation in STM images.
- BIC assumes all parameters contribute equally, but extra sigma parameters
  can absorb noise without improving predictive accuracy.
- On ~10% of files, BIC marginally prefers over-fit models (Δ < 100).

Cross-validation (CV) directly estimates out-of-sample prediction error
and is more robust against these effects.

## Why min(ell, circ)?

The circular model (σ∥ = σ⟂ per lobe) is a **nested special case** of
the elliptical model. Therefore:

- For any N, there exists an elliptical solution with BIC ≤ circ_BIC(N)
  (set σ∥ = σ⟂ and you get the circular solution exactly).

If we observe `ell_BIC(N) > circ_BIC(N)`, the elliptical optimizer failed
to find the global minimum. Using `min()` guards against this failure
without introducing any new parameters.

## CV Tiebreaker Details

```julia
if cv_ratio > 2.0
    # CV strongly prefers the simpler model
    # e.g., CV(6)=0.07, CV(8)=0.42 → ratio=5.7 → select N=6
    select N_cv
elseif ΔBIC < 100 && N_cv < N_bic
    # Small BIC margin and CV prefers simpler
    # e.g., BIC(6)=2085, BIC(7)=2078 (Δ=7), CV(6)=0.79, CV(7)=0.97
    select N_cv
end
```

## circ→ell LsqFit Refinement

**All elliptical fitting is now done via circ→ell LsqFit refinement.**
The NLopt global optimizer is intentionally excluded — it always diverges
from the isotropic start in 33D parameter space (see Research Journal §7).

For each N fitted by the circular sweep:
1. Circular solution provides positions, amplitudes, and isotropic σ
2. Params expanded to elliptical format (σ∥ = σ⟂ = σ_circ initially)
3. LsqFit only (skip_global=true, max_iter=50) locally refines sigmas
4. If refined BIC < circ BIC → use refined elliptical result

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
