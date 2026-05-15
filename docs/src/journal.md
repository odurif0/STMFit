# Research Journal

Chronological record of investigations into the 2D elliptical chain model
convergence and model selection problem. Includes both successful and
unsuccessful approaches, with rationale.

---

## Problem Statement (May 2025)

Batch processing of 27 chitosan STM images (240817 dataset). The 2D elliptical
chain model produces N ≠ 6 for 7/27 files despite theoretical expectation of
6 monomers per chain. Three files systematically wrong: 019 (N=8), 026 (N=8),
051 (N=8).

**Goal**: achieve N=6 for all files where the molecule has 6 monomers,
without introducing heuristic/arbitrary parameters.

---

## Investigation Timeline

### 1. Initial Diagnostic (May 14)

**Observed**: 27 files processed with batch v5. Elliptical model selects N=8
for files 019, 026, 051. Circular model agrees on 019/051 but selects N=6
for 026 (PROBLEMATIC class). 1D fit overestimates N vs 2D by ~+1.

**Key data points**:
- 019: ell N=8 BIC=1125, N=6 BIC=1181 (Δ=56). Both ell and circ agree N=8.
- 026: ell N=8 BIC=1563, N=6 BIC=1935. **Circ N=6 BIC=1524** ← circ better!
- 051: N=8 legitimate (support 5.77 nm vs avg 4.13 nm)
- Support mismatch (1D vs 2D): 019=19%, 026=26%

### 2. 026 Deep Dive — Why ell N=6 diverges (May 14-15)

**Finding**: Elliptical N=6 consistently converges to BIC=1935 while circular
N=6 finds BIC=1524. The circular model has FEWER parameters but LOWER BIC
— a mathematical impossibility if both converge to the global minimum.

**Root cause**: The elliptical model's independent σ∥ and σ⟂ parameters create
a loss landscape where the isotropic solution (σ∥=σ⟂) is a **saddle point**.
NLopt+LsqFit slides away from this saddle into a worse basin (BIC=2009).
Even 20 random starts cannot escape this basin.

**Tests performed**:
- Standard elliptical N=6: BIC=2009 ✗
- Warm-start from circular solution: BIC=2009 ✗ (same basin)
- Multistart ×3, ×5, ×10, ×20: all BIC=2009 ✗ (same basin)
- LsqFit-only from circular: BIC=1835 (better, but still worse than circ)

**Bug found**: The `starts` parameter in `_fit_chain_n` was declared but never
used — no multistart loop existed. **Fixed** by implementing the loop
(commit: multistart implementation).

**Solution adopted**: `min(ell_BIC, circ_BIC)` per N. Since the circular model
is nested within elliptical (σ∥=σ⟂ is a special case), circ_BIC is a
mathematically valid lower bound. When ell_BIC > circ_BIC, it indicates
convergence failure. Implemented in `_select_effective_best`.

---

### 3. 019 Investigation — Marginal N=8 preference (May 14)

**Finding**: BIC prefers N=8 by Δ=56. The 1D fit at N=8 has enormous center
errors (2-5 nm with spacing ~0.6 nm), suggesting instability. Edge lobe
amplitudes are 44% of max (0.055 vs 0.134).

**Cross-validation discovery**: The CV score (cross-validation NLL) strongly
prefers N=6 over N=8:
- Circ N=6: CV=0.552, N=8: CV=2.038 (ratio 3.7×)
- 5-fold CV: N=6: CV=0.074, N=8: CV=0.422 (ratio 5.7×)

CV is a more robust overfitting detector than BIC because it estimates
out-of-sample prediction error.

**Solution adopted**: CV tiebreaker in model selection:
- CV ratio > 2.0: CV strongly prefers simpler model → override BIC
- ΔBIC < 100 AND CV prefers simpler N → override BIC

---

### 4. Failed Approach: Endpoint Amplitude Penalty (May 14-15)

**Hypothesis**: Edge lobes with low amplitude are fitting noise, not real
monomers. Penalize them in BIC.

**Implementation**: Added `endpoint_amplitude_ratio` and
`endpoint_amplitude_penalty_weight` to ChainSweepConfig and FitConfig.
Penalty = max(0, threshold - amp/max_amp) × weight × log(N_eff).

**Result**: Rejected.
- Default threshold=0.5, weight=8.0: penalty too weak to flip 019 (Δ=56)
- To flip 019 would need threshold≥0.65, weight≥16 — would also fire on
  clean files (033 has edge ratio 0.58)
- The thresholds are arbitrary heuristics, not derivable from physics

**Lesson**: Edge amplitude ratio alone cannot distinguish real from phantom
lobes — the "right" threshold varies per file. CV tiebreaker is more robust.

**Files affected**: Reverted from types.jl (1D + 2D), core.jl (1D + 2D),
STMMolecularFit.jl.

---

### 5. Failed Approach: Sigma-Ratio Penalty (May 15)

**Hypothesis**: Penalize σ∥/σ⟂ asymmetry in BIC to prevent NLopt from
separating sigmas and diverging.

**Implementation**: Added `sigma_ratio_penalty_weight` to ChainSweepConfig.
Penalty = |σ∥-σ⟂|/σ_mean × weight × log(N_eff) per lobe.

**Result**: Rejected.
- With weight=10, shifts ell N=8 from 1563→1749 on 026, but ell N=6 still
  at 2015 (not fixed)
- Redundant with `min(ell, circ)` which already handles divergence
- Weight is an arbitrary free parameter

**Lesson**: The penalty attacks the symptom (σ∥/σ⟂ asymmetry), not the root
cause (NLopt escaping the isotropic basin). The circular model already
enforces σ∥=σ⟂ structurally — use it directly instead.

**Files affected**: Reverted from types.jl, core.jl.

---

### 6. Failed Approach: Re-parameterization to (σ_iso, Δ) (May 15)

**Hypothesis**: Replace independent σ∥,σ⟂ with coupled (σ_iso, Δ) where
Δ = log(σ∥/σ⟂) is bounded to ±0.25. The bound is derived from physical
constraints (tip convolution isotropy + chain sway amplitude).
Keeps optimizer in the isotropic neighborhood.

**Implementation**: Changed `_decode_chain` and `_pack_chain_initial` to
use σ_iso = √(σ∥×σ⟂) and Δ with sigmoid encoding bounded by
sigma_anisotropy_max.

**Result**: Rejected.
- Created a NEW failure mode: all Δ saturate at +0.25 (the bound), chain
  compresses to 1.24 nm span (should be 3.24 nm), BIC=114065
- The bound creates an attractive basin at the edge of the allowed range
- Same fundamental issue: any free parameter in sigma space creates a
  direction the optimizer exploits

**Lesson**: Any parameterization that allows σ∥≠σ⟂ per lobe creates exploitable
degrees of freedom. Hard bounds just move the divergence to the boundary.
The only invariant solution is the circular model (σ∥=σ⟂ structurally).

**Files affected**: Reverted from types.jl, core.jl.

---

### 7. Successful Approach: circ→ell LsqFit Refinement (May 15)

**Hypothesis**: The circular solution positions are near-optimal. Running
LsqFit-only (local, no global exploration) from the circular solution
in elliptical parameter space should find the true elliptical minimum
without escaping the good basin.

**Test results** (3 files, N=6):

| File | Circ BIC | Ell standalone | LsqFit from circ | Gain |
|------|----------|---------------|------------------|------|
| 033 (clean) | 3268 | 5188 ✗ | **3171** ✓ | 4× better |
| 038 (problem) | 34273 | — | **28360** ✓ | 1.2× better |
| 035 (problem) | 21000 | 15900 | **668** ✓ | 24× better |

**Crucially**: NLopt from circ (even 2s timeout) immediately diverges.
NLopt 8s + LsqFit from circ = same result as standalone ell (5188).
NLopt is **actively harmful** for the elliptical model — any global
exploration escapes the narrow good basin.

**Conclusion**: The optimal elliptical fit is obtained by:
1. Run circular sweep (converges reliably)
2. For each N, warm-start elliptical LsqFit from circular solution
3. Use min(circ_BIC, refined_ell_BIC) for model selection

The elliptical NLopt global optimizer should be **removed entirely** from
the elliptical fitting path — it cannot find the global minimum in 33D
parameter space and always makes the fit worse.

---

### 8. Bug Fix: Multistart Never Executed (May 14)

**Finding**: The `starts` parameter in `_fit_chain_n` was declared
(`starts::Int=ccfg.multistart`) but never used. The function body ran
exactly one optimization regardless of the parameter value.

**Fix**: Implemented a loop over `starts` iterations with random perturbation
of delta and sigma parameters for diversity. Each start runs NLopt+LsqFit
independently; the best result (by RSS) is returned.

**Impact**: Without this fix, all tests of `starts=3,5,10,20` were
meaningless — only 1 start was ever executed. The finding that "20 random
starts give the same result" was actually "1 start gives 1 result, repeated
20 times."

---

### 9. CV Computation Improvement (May 15)

**Finding**: CV computation using `_chain_cv_score` was slow and unreliable:
- Used `starts=3` per fold × 3 folds = 9 NLopt+LsqFit runs per N
- LsqFit Jacobian computation would hang on certain files at N≥7
- Caused timeouts during batch processing

**Fix**: CV fits now use:
- `skip_global=true` (no NLopt, LsqFit only)
- `starts=1` (single start per fold)
- `max_iter=50` (reduced from 300)

This makes CV ~10× faster and eliminates timeout issues.

**Also**: Added CV columns (`cv_nll_mean`, `cv_nll_std`) to score files
for post-hoc analysis.

---

### 10. 5-Fold CV (May 15)

**Finding**: 5-fold CV provides 2× better discrimination than 3-fold:
- 019: CV(6)=0.074, CV(8)=0.422 (ratio 5.7× vs 3.5× for 3-fold)
- Standard deviation is ~3× smaller

**Change**: `cv_folds=5` in batch config (was 3).

---

## Current Pipeline (v6)

```
Step 1: 1D slide profile extraction + peak fitting
        → bootstrap centers and amplitudes for 2D initialization

Step 2: Circular sweep (N = 2..14, adaptive range)
        → reliable convergence, isotropic gaussians

Step 3: Elliptical sweep (N = 2..14)
        → NLopt global + LsqFit local (currently kept for comparison)

Step 4: circ→ell LsqFit refinement at EACH N
        → warm-start from circular solution, local optimization only
        → finds true elliptical minimum without NLopt divergence

Step 5: Model selection = min(BIC_circ, BIC_ell_refined)
        → circular model is nested lower bound
        → refined elliptical when it genuinely improves

Step 6: CV tiebreaker
        → CV ratio > 2.0: override BIC (strong CV preference)
        → ΔBIC < 100 AND CV prefers simpler: override BIC

Step 7: Output best model (N, params, plots, scores)
```

**Selection criteria hierarchy**:
1. `min(BIC_circ(N), BIC_ell_refined(N))` — primary
2. CV tiebreaker — secondary (close calls only)
3. CV ratio > 2.0 — override (strong signal)
4. ΔBIC < 100 + simpler N preferred by CV — override (margin call)

---

## Lessons Learned

1. **Circular model is the anchor**: σ∥=σ⟂ enforced structurally, always
   converges. Use it as the reference in all comparisons.

2. **NLopt global optimizer is harmful for elliptical**: 33D parameter space
   is too large. The isotropic solution is a saddle point that NLopt always
   escapes. LsqFit-only from circular start is optimal.

3. **Min() is more robust than penalties**: Adding penalty terms to BIC
   introduces free parameters. Using `min(ell, circ)` achieves the same
   effect with zero new parameters.

4. **CV detects overfitting that BIC misses**: BIC is an asymptotic
   approximation; CV is a direct estimate of prediction error. Use CV
   as tiebreaker, not as primary.

5. **Re-parameterization doesn't fix optimizer topology**: Changing from
   (σ∥,σ⟂) to (σ_iso,Δ) just moves the divergence point. The fundamental
   issue is that any extra degree of freedom in sigma space can be exploited.

6. **1D over-estimates N**: The 1D fit has more flexibility (no 2D topology
   constraints) and sBIC penalizes less. This is documented but not fixed
   — the 2D selection is what matters for final output.

---

## Open Questions

1. **039 (N=7)** → **RESOLVED (May 15)**: Support detection was too generous
   (th=0.20, pad=0.20). The default support=4.64 nm included ~0.6 nm of edge
   noise, creating room for an artificial 7th lobe. With calibrated support
   (th=0.25, pad=0.05), support drops to 4.21 nm and N=6 wins (Δ=-192).
   See section 12 for calibration details.

2. **n_eff ÷9 factor** → **DEFERRED**: Documented in code as conservative estimate.
   Principled computation from pixel size would give ~÷25 but requires image metadata.
   Low priority — does not affect selection ranking.

3. **1D CV tiebreaker** → **WON'T FIX**: 1D only provides bootstrap initialization for 2D.
   The final N selection is done by the 2D pipeline (min + CV). Adding CV to 1D
   would add complexity without affecting final output.

4. **Remove elliptical NLopt sweep** → **RESOLVED (May 15)**: Replaced by circ→ell LsqFit
   refinement in the batch pipeline. `_refine_circ_to_ell()` runs LsqFit-only elliptical
   warm-started from circular for each N. NLopt is intentionally excluded — proven
   harmful in tests (diverges in 2s even from circular start). See section 7.

---

## 11. Pipeline v6 — NLopt Elliptical Removed (May 15)

**Decision**: The elliptical NLopt global optimizer is removed from the pipeline.
All elliptical fitting is now done via circ→ell LsqFit refinement.

**Rationale**: 
- NLopt diverges in 2s even when warm-started from the circular solution (section 7)
- circ→ell LsqFit gives 3-24× better BIC than NLopt elliptical on all tested files
- The circular model finds good positions; LsqFit locally refines sigmas without escaping

**Implementation**: `_refine_circ_to_ell()` in batch_full.jl:
1. For each N fitted by the circular sweep
2. Expand circ params to elliptical format (duplicate σ for σ∥ and σ⟂)
3. Run LsqFit-only (skip_global=true, max_iter=50) from circ warm-start
4. Finalize BIC and validity checks

**Impact**:
- Pipeline is 2× faster (no separate elliptical NLopt sweep)
- Model selection uses `min(circ_BIC, refined_ell_BIC)` — more accurate
- Zero convergence failures observed across all N on test files

---

## 12. Support Detection Calibration (May 15)

**Finding**: The default support parameters (th=0.20, pad=0.20) detect
supports that are too wide, including edge noise that creates artificial
lobes at higher N. This is the root cause of 039 N=7 and 036 N=7.

**Mechanism**: With generous support detection, the optimizer has room to
fit extra lobes at the edges. For 039, the default support=4.64 nm allowed
7 lobes at 0.74 nm spacing. The actual molecule is ~4.0 nm (6 lobes at
0.67 nm spacing). The extra 0.64 nm of detected "support" was noise.

**Calibration grid search** (039 + 033/036/042):

| th   | pad  | 039 | 036 | 033 | 042 | 019 |
|------|------|-----|-----|-----|-----|-----|
| 0.20 | 0.20 | N=7 | N=7 | N=6 | N=6 | N=8 |
| 0.25 | 0.00 | N=6 | N=7 | N=6 | N=6 | N=6 |
| 0.25 | 0.05 | N=6 | N=6 | N=6 | N=6 | N=6 |
| 0.22 | 0.05 | N=6 | N=6 | N=6 | N=6 | N=6 |

**Selected**: `support_threshold_fraction=0.25`, `support_padding_nm=0.05`.
- Threshold 0.25: requires stronger signal for support (was 0.20)
- Padding 0.05 nm: minimal edge extension (~2.5 STM pixels, was 0.20)
- Supports drop from avg 4.13→3.86 nm on clean files
- BIC values IMPROVE on clean files (033: 1349→1162)

**Impact**: Fixes 039 AND 036 without breaking any clean file. The tighter
support is physically justified — 0.05 nm padding is ~2-3 STM pixels,
which is the minimum needed to include the lobe's sigma tail.

**Absolute BIC tracking**: When calibrating, compare absolute BIC values,
not just which N wins. A setting that makes N=6 win by degrading N=7
(raising both BICs) is worse than a setting that genuinely improves N=6.
For chitosan: 033 BIC 1349→1162 (genuine improvement), 042 BIC 1101→927
(genuine improvement), 036 both BICs rise (trade-off).

**Calibration file**: Created `config/chitosan.toml` with all parameters.
Users can copy and adjust for different molecules/instruments. The batch
script accepts `--config path/to/file.toml`.

---

## File Change Summary (cumulative, v6)

| File | Changes |
|------|---------|
| `GaussianFit2D.jl/src/core.jl` | Multistart loop in `_fit_chain_n`; Lighter CV (`skip_global`, starts=1); n_eff documentation |
| `GaussianFit2D.jl/src/types.jl` | No net changes (penalties added then reverted) |
| `GaussianFit1D.jl/src/core.jl` | No net changes (endpoint penalty added then reverted) |
| `GaussianFit1D.jl/src/types.jl` | No net changes |
| `STMMolecularFit.jl/src/STMMolecularFit.jl` | No net changes |
| `test/batch_full.jl` | `_select_effective_best` (min BIC + CV tiebreaker); `_refine_circ_to_ell` (replaces NLopt ell sweep); `cv_folds=5`; CV in scores; `N_eff`/`eff_source` columns |
| `test/summarize.jl` | Multi-file prefix glob support |
| `test/inspect_one_file.jl` | Single-file diagnostic |
| `docs/` | Full documentation suite (7 files, 700+ lines) |

---

## v7 — Pseudo-Voigt, Covariance, Residual Diagnostics (2026-05-15)

### Features Added

**C — Residual diagnostics** (`STMFitCore.jl`):
- `durbin_watson(residuals)`: tests residual autocorrelation. DW≈2
  = no autocorrelation, DW<1.5 = missed structure.
- `runs_test(residuals)`: Wald-Wolfowitz runs test on residual signs.
  Too few runs = systematic bias, too many = overfitting.
- `ResidualDiagnostics` struct + `compute_residual_diagnostics`: unified API.
- Integrated into `FitResult` (1D) and `ChainModelResult` (2D).
  Exported in results TSV and `print_summary`.

**B — Covariance quantification** (`GaussianFit1D.jl`):
- `FitResult.pcorr`: full parameter correlation matrix from LsqFit
  `estimate_covar`, normalized to unit diagonal.
- `FitResult.center_center_corr`: peak position correlation matrix
  via Jacobian of delta→center transformation.
- `max center-center correlation` shown in `print_summary`.

**B — 2D parameter errors** (`GaussianFit2D.jl`):
- `ChainModelResult.param_perr`: extracted from `estimate_covar` inside
  `_run_one_start`, tracked across multistarts.
- Numeric propagation to axial positions via finite differences in
  `_chain_metrics!`.

**A — Pseudo-Voigt profile** (`GaussianFit1D.jl`):
- `peak_profile` field in `FitConfig`: `:gaussian` (default), `:lorentzian`,
  or `:pseudo_voigt`.
- Global mixing parameter η∈[0,1] shared across all peaks (1 extra param).
- Forwarded through `STMMolecularFit.FitSlideConfig`.
- 2D guarded to `:gaussian` only (`ChainSweepConfig.peak_profile`).

### Bugs Found and Fixed

**Bug 1: NLopt crash on pseudo-Voigt** (`core.jl:324`):
`_make_objective_function` sized `full_buf = 3*n + extra` without
accounting for the η parameter. The NLopt solver passed a vector
including η, causing `DimensionMismatch` on broadcast.
**Fix**: `+ (peak_profile == :pseudo_voigt ? 1 : 0)`.

**Bug 2: Student-t BIC per-fit noise adaptivity** (`core.jl:592-594`):
`compute_metrics` estimated noise as `max(std(resid)*0.1, MAD(resid))`
*per model*. Since better fits have smaller residuals, the noise estimate
shrinks proportionally, keeping `resid/noise` roughly constant. The
Student-t NLL does not decrease with RSS → BIC penalizes better fits.
On 240817_002, N=4 (RSS=0.068, sBIC=511) beat N=6 (RSS=0.003, sBIC=628)
despite R²=0.93 vs 0.997.
**Fix**: Compute noise ONCE from the lowest-RSS model's residuals (MAD),
then recompute sBIC for all models with this fixed reference. Consistent
with 2D (fixed preprocessing noise).

**Bug 3: Incomparable sBIC across profiles** (`core.jl:855-858`):
Each sweep computed its own noise reference → Gaussian and pseudo-Voigt
had different noise estimates (0.00252 vs 0.00288) → sBIC incomparable.
Pseudo-Voigt appeared to beat Gaussian (sBIC 585 vs 618) because its
higher noise reduced NLL.
**Fix**: If `cfg.noise_estimate` is already finite (set from Gaussian
sweep), reuse it. Gaussian sweep auto-stores its noise in
`cfg.noise_estimate` for subsequent profiles.

### Validation on Chitosan Data (240817_002, support 4.59 nm)

| Profile | N | sBIC | R² | DW | η |
|---------|---|------|-----|------|----|
| Gaussian | **6** | 618 | 0.997 | 0.337 | — |
| Pseudo-Voigt | 6 | 699 | 0.996 | 0.292 | 0.000 |
| Lorentzian | 6 | 1868 | 0.973 | 0.044 | — |

- Gaussian wins (Δ=81 over PV) — η=0.000 confirms Gaussian profile.
- DW=0.337 at N=6 vs 0.018 at N=4 — DW detects the underfit that BIC
  alone would miss.
- Lorentzian clearly rejected (Δ=1250).

### Tests

- STMFitCore: 32/32 ✓
- Feature suite (DW logic, covariance, PV, export, params, version):
  53/53 ✓

### Files Modified

| File | Changes |
|------|---------|
| `STMFitCore.jl/src/STMFitCore.jl` | +_norm_cdf, durbin_watson, runs_test, ResidualDiagnostics, compute_residual_diagnostics (+85 lines) |
| `STMFitCore.jl/test/runtests.jl` | +residual diagnostics tests (+40 lines) |
| `GaussianFit1D.jl/src/types.jl` | +peak_profile, pcorr, center_center_corr, residual_diagnostics fields; MGF_VERSION→5.0.0 |
| `GaussianFit1D.jl/src/core.jl` | pseudo-Voigt dispatch in multi_gaussian, global noise fix in run_model_comparison, covariance/correlation computation, residual diags integration, bounds/params for η |
| `GaussianFit1D.jl/src/GaussianFit1D.jl` | +import ResidualDiagnostics, compute_residual_diagnostics |
| `GaussianFit2D.jl/src/types.jl` | +peak_profile, param_perr, residual_diagnostics |
| `GaussianFit2D.jl/src/core.jl` | covariance extraction in _run_one_start, residual diags in _finalize_chain_result! |
| `GaussianFit2D.jl/src/GaussianFit2D.jl` | +import ResidualDiagnostics, compute_residual_diagnostics |
| `STMMolecularFit.jl/src/STMMolecularFit.jl` | +peak_profile forwarding in FitSlideConfig, fit_slide |
| `docs/src/math.md` | New: κ, pseudo-Voigt, Student-t BIC, diagnostics, uncertainty math |
| `docs/src/api.md` | +new functions (DW, runs, ResidualDiagnostics) |
| `docs/src/config.md` | +peak_profile documentation |
| `docs/make.jl` | +math.md page |
| `config/chitosan.toml` | Calibration file with all parameters |
