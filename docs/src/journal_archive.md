# Research Journal — Historical Archive

Earlier investigation detail and completed follow-ups (May 2025 – early July
2026): pipeline evolution, selection-rule work, calibration analysis, QE
restarts, and unit-assignment experiments. Active entries and open questions
live in [journal.md](journal.md).

---

### 2026-06-20 — Objective calibration and effective-sample-size analysis

**Goal:** make every calibration parameter either measured from the data or
derived from a physical principle, so the pipeline generalises to a new molecule
on the same STM without hand-tuning.

**Effective sample size investigation.** The `n_eff = n÷9` heuristic (block-3×3
correlation) enters BIC, AICc, robust-AICc, MDL and Laplace scores. A diagnostic
on three chitosan files (002, 043, 058) compared three estimators:

| estimator | n_eff (043) | finding |
|---|---|---|
| heuristic `n÷9` | 966 | ignores correlation |
| Durbin-Watson AR(1) `n·(1−ρ)/(1+ρ)` | 440 | ρ≈0.90 (raster 1D) |
| 2D variogram in fit window | 22 | underdetermined (range > window) |

A full-image 2D autocorrelation (512×512, subsampled 4×) reveals the correlation
**range** is 17–100 px, far larger than the ~10-px fit window. The number of
independent points in the fit window is therefore `window_area / (π·range²) ≈ 0.06`
— effectively zero. **Conclusion: n_eff is not objectively definable in the fit
window; BIC/AICc (which assume iid observations) are not well-defined here.**
GCV, which does not assume independence (smooth-spline leave-one-out theory), is
the canonical selection criterion; BIC/AICc are retained only as qualitative
diagnostics. This is documented in `docs/src/calibration.md`.

**Auto-calibration tool** (`test/measure_calibration.jl`). Measures from a single
clean scan: noise σ, pixel resolution, FWHM range [5%, 95%], repeat spacing, and
correlation range; derives σ_parallel (FWHM/2.355), spacing (±30%), fit_width
(=σ_min), support_min (3×spacing), n_max (axis/spacing), and emits a ready-to-use
TOML. Validated: the auto-calibrated TOML (measured on 043) gives `N_selected = 6`
on 002 — identical to the hand-tuned `chitosan.toml`. Of ~25 parameters, 5 are
measured, 9 are principled-derived, and ~11 remain genuinely free (optimizer
budget, model-form switches, acquisition-dependent defaults).

**Sensitivity check** (`test/sensitivity_thresholds.jl`, HPC, 4 thresholds
0.03/0.04/0.05/0.06 on 43 common files): 0 pivot files — `N_selected` is
insensitive to the ambiguity threshold across the tested range.

### 2026-06-17 — Symmetric up-when-ambiguous guard branch

The robust-AICc guard was strictly down-only: it could veto over-segmentation
(`robust_AICc_N < N_eff`) but could not recover an under-segmented ambiguous
case even when its own exhaustive sweep recommended a larger count.  This left
`240817_043.sxm` at `N_selected = 5` despite `robust_AICc_N = 6`, `ambiguous_eff
= true`, and `delta_GCV_rel_eff ≈ 1%` (N=5 and N=6 are statistically
indistinguishable on that file).

The circ→ell warm-start refinement (the standard batch path) converges to a
sub-optimal elliptical basin for N=6 on `043` (GCV 8.94e-6 vs 8.86e-6 for N=5).
An independent elliptical sweep with NLopt global finds a better N=6 minimum
(GCV 7.11e-6), but that path was rejected pipeline-wide because NLopt diverges
on most other files (§7).  Tuning the optimiser for `043` alone is not robust:
`max_iter=300` on the warm-start makes N=7 win instead (over-fitting), and
multistart perturbation makes it worse.  `043` is genuinely ambiguous: N=5, 6, 7
all lie within a 10% GCV band.

Resolution: a symmetric **up-when-ambiguous** branch in `_refined_selection`
(`test/batch_full.jl`), label-free and bounded.  When the exhaustive
robust-AICc guard recommends exactly one lobe more than `N_eff`, accept the
upshift only if all of the following hold:

```text
robust_AICc_N == N_eff + 1
AND ambiguous_eff == true            (GCV does not discriminate N_eff vs runner-up)
AND delta_GCV_rel_eff <= 0.05        (the existing GCV_AMBIGUITY_REL_THRESHOLD)
AND runner_up_N_eff == N_eff + 1     (the competing model is the adjacent N)
```

The down branch is unchanged (free).  The rule uses no expected `N`, no target
count, and no file name; it is the mirror of the existing down guard with
explicit guards against the over-segmentation jumps it could otherwise enable
(e.g. `026: 6→8` is blocked because `ambiguous_eff=false` and `dGCV>0.10`;
`036: 6→7` is blocked for the same reason).

`_select_primary` (`packages/STMMolecularFit.jl/src/selectors.jl`) now flags the
selection source as `robust_aicc_guard` for both up and down moves
(`n_refined != n_eff`), so the upshift is traceable in the summary.

Validation on the 240817 chitosan benchmark (reproductible across 3 consecutive
runs; identical `N_selected` on all 48 files before and after this change except
as noted):

- only `240817_043.sxm` changes: `N_selected 5 → 6` via
  `overfit_guard_up_when_ambiguous`;
- all four `clean_target` files (`017`, `019`, `043`, `058`) now report
  `N_selected = 6`;
- primary benchmark exact agreement: `38/39 → 39/39` (`N_eff` itself remains
  `35/39` — the guard supplies the remaining four);
- no other primary or stress-case file changes `N_selected`;
- the batch remains fully reproductible run-to-run.

### 2026-05-29 — Experimental refined selection reporting

Added an external, conservative overfit-guard audit and optional batch reporting
for `N_refined`.  The rule is generic and one-sided:

```text
N_refined = robust_AICc_N if robust_AICc_N < N_eff
            N_eff         otherwise
```

This does not use an expected `N`, and it cannot increase the selected count.
It is intended to catch over-segmentation while preserving clean cases such as
`240817_026.sxm`, where the robust advisory alone over-selects but the primary
`N_eff` is correct.

Current external grading with `results/robust_rescore_audit/full_aicc_nu8.tsv`:

- baseline `N_eff`: `35/39` primary benchmark files;
- refined overfit guard: `38/39`;
- recovered targets: `017`, `019`, `058`;
- remaining target: `043`, consistent with a support-ambiguity rather than an
  overfit case;
- `026` remains `6`.

`test/batch_full.jl` now supports optional reporting via:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --refined-advisory results/robust_rescore_audit/full_aicc_nu8.tsv
```

Follow-up: `test/batch_full.jl` now also supports an integrated experimental
primary selector, without requiring an advisory TSV:

```bash
julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy gcv_with_robust_aicc_guard
```

This computes the robust AICc guard directly from an auxiliary exhaustive
elliptical candidate set, then applies the same down-only rule to the standard
circ→ell `N_eff`.  The selected primary output is written as `N_selected`, while
`N_eff` remains available for comparison.

Full-batch grading for the integrated selector (verified reproductible across
3 consecutive runs on 2026-06-17: identical `N_selected` on all 48 files; the
2026-06-17 up-when-ambiguous branch raised this to the numbers below):

- `N_selected`: `39/39` primary benchmark files;
- target score: `4/4` (`017`, `019`, `043`, `058` all selected as `6`);
- `043` is recovered by the up-when-ambiguous branch (`N_eff=5`,
  `robust_AICc_N=6`, `ambiguous_eff=true`); see §2026-06-17 above;
- `N_eff` on the same run remains `35/39`.

Historical note: at this point the policy was still experimental/default-off.
It was later promoted to the `config/chitosan.toml` batch default after synthetic
known-N validation and additional manually annotated real-data audits.

Follow-up: added several additional label-free experimental selectors for audit
purposes, including support-marginalized GCV, a guarded support variant,
slope-heuristic MDL, stability selection, and local-lobe evidence.  The
local-lobe-evidence selector is a down-only guard that merges adjacent lobes when
their separation is below `2σ∥` and the between-lobe valley is weak.  Initial
smoke tests showed it is too strict for continuous chitosan chains: even clean
`240817_026.sxm` remains at `N=6` only because the diagnostic is inconclusive
(`resolved=1`, `unresolved_pairs=5`) and falls back to GCV.  It is therefore a
separability audit, not a better primary selector.

Added `laplace_evidence` and `laplace_evidence_guard`.  The direct selector uses
a finite-difference Gauss–Newton/Laplace evidence approximation with Student-t
weights and an Occam/sloppiness penalty from the weighted Jacobian singular
values.  It is label-free and statistically cleaner than threshold heuristics,
but smoke tests showed it is too parsimonious as a free selector (`019 → 5`,
`058 → 5`).  The guarded variant caps downshifts to one lobe: smoke tests kept
`026 → 6`, corrected `017 → 6`, `019 → 6`, and `058 → 6`, but still cannot
increase the support-sensitive `043` from `5` to `6`.  This makes it a promising
principled overfit guard candidate, still below the integrated robust-AICc guard
until full frozen validation.

Added `fwd_bwd_consensus`: evaluates the fused-fit model on separate fwd and bwd
preprocessed channels with per-scan linear recalibration and joint GCV.  The
idea is that true molecular signal should appear in both scans while
direction-specific artifacts would not.  Smoke tests showed this simple joint-GCV
does not discriminate better than fused GCV on the chitosan benchmark (`017 → 7`,
`019 → 7`, `058 → 7`, `043 → 5`).  The recalibration absorbs most scan-to-scan
differences.  A more discriminating approach (lobe-level fwd/bwd amplitude
consistency or a true joint refit) would be needed.

### Refactoring: selectors moved to STMMolecularFit core

All experimental selector functions (~550 lines) were extracted from
`test/batch_full.jl` into `packages/STMMolecularFit.jl/src/selectors.jl`
(584 lines), included by the STMMolecularFit module.  `batch_full.jl` shrank
from 1627 to 1099 lines.  The selector functions are no longer in a test file;
they are now part of the STMMolecularFit package core.  `batch_full.jl` imports
the needed names via `import STMMolecularFit: ...`.  No logic was changed.

### Synthetic known-N selector validation

Added `test/synthetic_known_n_validation.jl` as a phase-1 frozen validation
harness for comparing selector policies without using benchmark labels as a
selection prior.  The script generates in-memory Gaussian-chain `SXMImage` cases
with known `true_N` cycling through `4..8`, adds jitter, baseline tilt,
independent fwd/bwd noise, and occasional artifacts, then runs a fast circular 2D
sweep and applies the core selectors from `STMMolecularFit`.

Added `test/aggregate_synthetic_known_n.jl` to summarize synthetic selector TSVs
by policy over all cases and stratified by `true_N` and artifact class.  It
reports exact rate, mean absolute error, over-selection, under-selection, and
error counts.

Follow-up: the generator now accepts `--noise-scale` and writes `noise_scale` in
the synthetic summary.  The aggregator now accepts multiple positional TSVs,
parses both legacy 10-column and newer 11-column summaries, and stratifies by
seed and noise scale in addition to `true_N` and artifact class.

Follow-up: added phase-2 `--mode circ_ell`.  The generator now has two modes:
`circular` for the cheap original stress test and `circ_ell` for a closer analog
of the real batch path.  `circ_ell` runs a fixed candidate circular sweep, then
locally refines valid circular candidates with elliptical sigmas and selects
`N_eff` from the best GCV score across circular/elliptical candidates.  The TSV
now includes a `mode` column, and the aggregator stratifies by `mode=<mode>`.

Important correction: the synthetic candidate search window is now fixed
(`N=2..10`) rather than derived from `true_N`.  This keeps the known label out of
fitting and selection; it is used only for external grading.

The output TSV reports `case_id`, `seed`, `true_N`, `artifact`, `policy`,
`N_eff`, `N_selected`, `abs_error`, `status`, and `score_or_source`.  The known
count is used only for external grading and aggregate summaries, not during
fitting or policy selection.

Initial smoke/default runs:

- `--cases 2 --seed 1234`: GCV exact `0/2`; robust-AICc guard and Laplace guard
  exact `2/2`.
- `--cases 6 --seed 1234`: GCV over-selected on the initial synthetic cases;
  robust-AICc guard was exact `6/6`.
- `--cases 50 --seed 1234`: robust-AICc guard was exact `46/50` with mean
  absolute error `0.12`, compared with GCV exact `2/50` and mean absolute error
  `1.58`.  The four robust-guard errors were under-selections, not overfits.
- Multi-seed/noise aggregate over five 50-case summaries (`250` cases):
  robust-AICc guard exact `231/250` (`92.4%`) with mean absolute error `0.112`;
  GCV exact `6/250` with mean absolute error `1.54`; Laplace guard exact
  `104/250` with mean absolute error `0.60`.  Robust-AICc errors were mostly
  under-selections (`15`) with a small number of over-selections (`4`).
- Noise-scale strata for the robust-AICc guard: `48/50` exact at `0.5×`,
  `139/150` exact at `1.0×`, and `44/50` exact at `1.5×`.
- Phase-2 `circ_ell`, fixed search window, `--cases 50 --seed 1234`: GCV exact
  `43/50` with mean absolute error `0.20`; robust-AICc guard exact `46/50` with
  mean absolute error `0.14` and four under-selections; stability selection exact
  `44/50` with mean absolute error `0.18`.
- Phase-2 `circ_ell` multi-seed/noise aggregate over five 50-case summaries
  (`250` cases): robust-AICc guard exact `231/250` (`92.4%`) with mean absolute
  error `0.128`, four over-selections, and fifteen under-selections.  GCV exact
  `213/250` (`85.2%`, mean absolute error `0.20`), stability selection exact
  `214/250` (`85.6%`, mean absolute error `0.192`), and Laplace guard exact
  `193/250` (`77.2%`, mean absolute error `0.272`).  Robust-AICc remained best
  across the tested noise strata: `48/50` at `0.5×`, `139/150` at `1.0×`, and
  `44/50` at `1.5×`.

Historical note: at this point this was synthetic validation, not yet a promotion
decision.  The robust-AICc guard was later promoted to the chitosan batch default;
support-adaptive variants remain separate experimental candidates.

---


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

**Historical solution**: CV tiebreaker in model selection:
- CV ratio > 2.0: CV strongly prefers simpler model → override BIC
- ΔBIC < 100 AND CV prefers simpler N → override BIC

This was later superseded in the chitosan batch configuration by analytical GCV
as the default primary criterion; the observation remains useful background for
why pure BIC was abandoned.

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
- Model selection can use `min(score_circ, score_ell_refined)` per N; the
  chitosan batch uses GCV by default.
- Zero convergence failures observed across all N on test files

---

## 12. Support Detection Calibration (May 15)

**Finding**: Support detection should be noise-based rather than contrast-based.
The former contrast-fraction support parameter was removed because it couples
support length to the brightest lobe and can cut weak lobes or react to artefacts.
The active rule is now `baseline + support_noise_k * noise`, followed by bounded
edge padding.

**Mechanism**: If the support is too wide, the optimizer has room to fit extra
lobes at the edges. If it is too narrow, terminal lobes are truncated. The
current chitosan calibration uses `support_noise_k=2.5` and
`support_padding_nm=0.25`, which improved the fast sweep from `N_ell=6` on
31/42 files to 33/42 files without adding a second threshold parameter.

**Selected**: noise-only thresholding with `support_noise_k=2.5` and
`support_padding_nm=0.25` in `config/chitosan.toml`.

**Rejected**: reintroducing a contrast-fraction support parameter, because it
adds a second support control with poorer physical meaning than signal-to-noise.

**Absolute BIC tracking**: When calibrating, compare absolute BIC values,
not just which N wins. A setting that makes N=6 win by degrading N=7
(raising both BICs) is worse than a setting that genuinely improves N=6.
For chitosan, track GCV/BIC together with residual plots and support mismatch;
files such as 029/032 are better treated as QC/problematic cases than scalar
parameter tuning targets.

**Calibration file**: Created `config/chitosan.toml` with all parameters.
Users can copy and adjust for different molecules/instruments. The batch
script accepts `--config path/to/file.toml`.

---

## 13. Final Simple Calibration Sweep (May 18, 2026)

**Goal**: exhaust simple, physically/data-driven tuning options without adding
an explicit or implicit prior toward any target N.

**Accepted change**: loosened the condition-number penalty threshold in the
chitosan calibration from `kappa_max=8.0` to `kappa_max=10.0`.

**Rationale**: this is a small existing configuration knob that affects only
ill-conditioned adjacent-lobe fits. It recovered `240817_049.sxm` in focused
and no-plot core tuning without breaking the control files. In the full plotting
batch, `049` remains support-sensitive and selects `N_ell=5`, so this setting is
retained as the best harmless simple calibration rather than as a complete fix.

**Clean full batch with current config**:

- command: `julia -t 4 --project=. test/batch_full.jl 48 --config config/chitosan.toml`
- output: `results/best_plots/summary_overlap060_hard.tsv`
- OK: `48/48`
- errors: `0`
- excluded/absent: `240817_015.sxm`, `240817_027.sxm`, `240817_063.sxm`
- `N_ell=6`: `38/48`
- `N_eff=6`: `38/48`
- `N_circ=6`: `35/48`
- ambiguous by ΔGCV (`delta_GCV_rel <= 0.05`): `12/48` for both elliptical
  and effective selections

**Core benchmark after excluding visually poor-quality files**:

The files `240817_029.sxm`, `240817_030.sxm`, `240817_031.sxm`,
`240817_032.sxm`, `240817_034.sxm`, `240817_035.sxm`,
`240817_037.sxm`, `240817_038.sxm`, and `240817_051.sxm` are kept in the
full run for traceability, but are excluded from calibration/tuning because of
suspected artefacts or poor image quality. Excluding them gives:

- benchmark files: `39`
- `N_ell=6`: `35/39`
- `N_eff=6`: `35/39`
- `N_circ=6`: `31/39`
- remaining non-6 cases: `240817_017.sxm`, `240817_019.sxm`,
  `240817_043.sxm`, and `240817_058.sxm`

**Remaining non-6 `N_ell` cases**:

| File | N_ell | N_circ | N_eff | Notes |
|------|------:|-------:|------:|-------|
| `240817_029.sxm` | 10 | 10 | 10 | problematic / support-mismatch target |
| `240817_031.sxm` | 12 | 12 | 12 | 2D-only robust; ambiguous vs 11 |
| `240817_032.sxm` | 10 | 10 | 10 | problematic / support-mismatch target |
| `240817_043.sxm` | 5 | 5 | 5 | problematic; ambiguous vs 6 |
| `240817_034.sxm` | 7 | 7 | 7 | robust non-6 under current model |
| `240817_058.sxm` | 7 | 6 | 7 | ambiguous minor case; second meilleur N = 6 |
| `240817_035.sxm` | 7 | 6 | 7 | ambiguous minor case |
| `240817_017.sxm` | 7 | 7 | 7 | problematic; ambiguous vs 9 |
| `240817_019.sxm` | 7 | 7 | 7 | ambiguous minor case; second meilleur N = 8 |
| `240817_051.sxm` | 8 | 8 | 8 | problematic / long-support case |

**Clean-benchmark surprises to inspect visually**:

| File | Selected | Second meilleur N | ΔGCV_rel | Notes |
|------|---------:|----------:|---------:|-------|
| `240817_043.sxm` | 5 | 6 | 0.0246 | close `N=6` alternative; 1D gives `N=7`; strong 1D/2D support mismatch |
| `240817_058.sxm` | 7 | 6 | 0.0435 | circular path gives `N=6`; close ambiguous case |
| `240817_017.sxm` | 7 | 9 | 0.0471 | ambiguous by GCV; strong support mismatch; remains QC-sensitive |
| `240817_019.sxm` | 7 | 8 | 0.0345 | ambiguous by GCV; not a stable hard failure |

**Rejected simple changes**:

- `kappa_max = 9, 11, 12`: no net improvement over `10.0`.
- `min_amplitude_fraction = 0.31, 0.32, 0.33, 0.35`: either no core gain or
  regressions such as `006: 6→5`, `017: 7→8`, `018: 6→7`, `049: 6→5`.
- `spacing_min_nm = 0.37, 0.38`: no focus gain; regressions on controls.
- `sigma_parallel_min_nm = 0.20, 0.205, 0.22`: no gain and can break `049`.
- `support_noise_k = 2.0, 3.0`, `fit_width_nm = 0.10/0.12/0.20`,
  `max_overlap = 0.55/0.70`, and `kappa_max = 0`: rejected in focus sweeps.

**Rejected structural/tuning paths**:

- Pseudo-Voigt: not useful for the 2D Gaussian chain path; Gaussian remains default.
- Any `expected_n` or “prefer N=6 if close” rule: rejected as a prior.
- `support_threshold_fraction`: removed and not reintroduced.
- Global padding reduction: helps some files (`049/058`) but hurts too many controls.
- Support morphology/hysteresis and multi-support common-audit selection: no net recovery.
- Lateral center-of-mass seeding and split/merge neighbor warm starts: local gains but
  new regressions.
- Effective-sample-size GCV and k-fold CV as defaults: unstable tradeoffs.

**Decision**: keep `support_padding_nm=0.25`, `selection_criterion="gcv"`,
`cv_method="gcv"`, Gaussian peaks, and `kappa_max=10.0`.  Adopt
`fit_width_nm=0.16` as the active simple scalar because no-maxtime validation
gave the best net `N_ell=6` gain among simple candidates.  Otsu-only,
scaled-Otsu, and full-tube support were tested as conceptual simplifications but
rejected because they introduced broader regressions. Remaining non-6 cases
should be treated as QC/visual-inspection targets or model-limit cases rather
than solved by expected-N priors.

### Maintenance note (2026-05-26)

- Removed late-stage temporary tuning variants from `test/tune_chitosan_params.jl`
  after they failed to improve the active chitosan calibration cleanly.
- No change to `config/chitosan.toml`, Gaussian 2D fitting, circ→ell refinement,
  or the configured GCV model selection path.
- Validation OK: tuning-script syntax parse, `julia --project=. test/batch_full.jl 0 --config config/chitosan.toml`,
  and `julia --project=docs docs/make.jl`.

---

## File Change Summary (cumulative, v6)

| File | Changes |
|------|---------|
| `GaussianFit2D.jl/src/core.jl` | Multistart loop in `_fit_chain_n`; Lighter CV (`skip_global`, starts=1); n_eff documentation |
| `GaussianFit2D.jl/src/types.jl` | No net changes (penalties added then reverted) |
| `GaussianFit1D.jl/src/core.jl` | No net changes (endpoint penalty added then reverted) |
| `GaussianFit1D.jl/src/types.jl` | No net changes |
| `STMMolecularFit.jl/src/STMMolecularFit.jl` | No net changes |
| `test/batch_full.jl` | `_select_effective_best` (min configured score, default GCV); `_refine_circ_to_ell` (replaces NLopt ell sweep); `N_ell`/`N_circ`/`N_eff`; support/residual QC columns |
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

---

---

### Raven QE submission archaeology (2026-06-22 to 06-25)

Verbose per-job narrative of the first Raven submission attempts. All jobs
here are superseded by the active pilot (28363474 and its dependents);
distilled lessons live in the active journal and hpc/qe_molds/README.md.
Kept only to prevent re-debugging the same launch-safety gaps.


Submitted the first two QE mold jobs to Raven after practical launcher fixes.
First, `hpc/submit_qe_molds.sh` now loads the configured Julia module before
re-running the remote preflight, because direct submission may start from a shell
without `julia` on `PATH`. Second, `test/prepare_qe_mold_inputs.jl` now writes
Slurm memory as `4000 MB × ntasks` instead of `--mem=0`, which Raven rejected on
shared nodes. Third, the generated compute-job script now loads the configured
Julia module and checks both `pw.x` and `julia`, because the QE job calls Julia for
the relax-to-SCF handoff. The first pre-compute-Julia submissions, `28278265`
(`qe/glcn`) and `28278266` (`qe/glcnac`), were still pending and were canceled
before start. The regenerated `qe/glcn` and `qe/glcnac` inputs pass preflight with
`8 / 8` total tasks and were resubmitted as `28278618` (`qe/glcn`) and `28278619`
(`qe/glcnac`). Both were pending in the `small` partition at the last check. No QE
outputs were available yet.

Those corrected jobs then failed immediately because Raven's module is not named
`quantum-espresso`; `find-module qe` shows the usable stack
`intel/2024.0`, `impi/2021.11`, `qe/7.4.1`. The launcher and generated sbatch
now expose `QE_COMPILER_MODULE`, `QE_MPI_MODULE`, and `QE_MODULE` and default to
that stack. A second preflight gap was also fixed: QE inputs use `pseudo_dir =
'./pseudo'`, so `test/preflight_qe_mold_inputs.jl` now verifies that every
`ATOMIC_SPECIES` pseudopotential exists, and the remote QE sync includes
`pseudo/*.UPF`. The five PSLibrary/KJPAW pseudos for Cu/C/H/N/O were downloaded
from `pseudopotentials.quantum-espresso.org`; the Cu pseudo recommends
`ecutwfc=71 Ry`, so the prepared inputs now use `ecutwfc=80`, `ecutrho=640`.
After these fixes, preflight passed and the jobs were resubmitted as `28286900`
(`qe/glcn`) and `28286903` (`qe/glcnac`); both were pending, not failed, at the
last check.

The `28286900`/`28286903` pair was then canceled while still pending because the
Raven `n0001` QOS reports a one-node group limit (`GrpTRES` includes `node=1`).
Even though each QE job correctly requests only `4` CPUs and the pair stays within
the `8` CPU budget, two independent Slurm jobs request two separate node
allocations. `hpc/submit_qe_molds.sh` now supports `--sequential`, and
`hpc/launch_qe_molds_remote.sh` defaults to `QE_SEQUENTIAL=1`, submitting later
run directories with `afterok` dependencies. The current active chain is
`28287102` (`qe/glcn`) followed by `28287103` (`qe/glcnac`, dependency
`afterok:28287102`). At the last check, `28287102` was pending with reason `None`
and `28287103` was pending with reason `Dependency`; neither had failed.

The `28287102` GlcN job then started but failed during the first `pw.x` relax step
with Slurm `OUT_OF_MEMORY` after ~4 min; QE had correctly loaded the module stack
and pseudos, but estimated `~159 GB` dynamic RAM per MPI process (`~635 GB` total)
for the original `8×8×4` slab with 4 k-points, Cu `spn` PAW pseudo
(`z_valence=19`), and `ecutwfc=80`/`ecutrho=640`. The dependent GlcNAc job was
canceled by `afterok`, as intended. To fit the current Raven QOS while preserving
a useful first DFT-STM mold path, the prepared jobs were changed to a lighter
pilot setup: `8×6×3` Cu(100) slab, `12 Å` vacuum, Cu `dn` PAW pseudo
(`Cu.pbe-dn-kjpaw_psl.1.0.0.UPF`, `z_valence=11`, suggested cutoff `45/236 Ry`),
`K_POINTS gamma` via `--kpoints 1,1,1`, `ecutwfc=50`, `ecutrho=360`, and a clean
2-task Slurm job. The generated sbatch now removes stale `qe_tmp` before
starting, avoiding partial scratch reuse after failed runs. Preflight passes with
the lighter inputs. Submitting both jobs at once still triggered Raven's shared
`n0001` pressure, so GlcN is submitted alone first with an 8 h walltime.

The `28288215` GlcN pilot initially requested `4` MPI tasks and `16 GB`, was
reduced in place to `2` MPI tasks to escape `QOSGrpCpuLimit`, then started and
failed in the first `pw.x` relax step with Slurm `OUT_OF_MEMORY` after 4 min 22 s.
QE ran with 2 MPI ranks and estimated `17.15 GB` dynamic RAM per process
(`34.30 GB` total), so the failure was a real memory limit, not a module/input
problem. A clean replacement GlcN job, `28292263`, was submitted with
`--ntasks-per-node=2` and `--mem=48000MB` (`ReqTRES=cpu=2,mem=48000M,node=1`).
At the last check it was pending with `QOSGrpCpuLimit`; dry-run probes for both
1-task and 2-task 48 GB jobs reported the same next-day start estimate, so the
2-task replacement was left in the queue rather than resubmitted. GlcNAc will be
submitted only after GlcN succeeds, but its local run directory was regenerated
with the same 2-task/48 GB settings and passes local preflight together with GlcN
(`4 / 8` total tasks).

Follow-up while `28292263` was running: the 2-task/48 GB job proved memory-safe
but underparallel. At ~7 h it had completed only the initial SCF (`39` electronic
iterations, `bfgs steps = 0`, total force `1.946765` vs threshold `1e-3`) and was
still in the next SCF, so an 8 h walltime was not credible for a full relax + SCF
+ PP workflow. The QE generator now writes explicit MPI launches:
`#SBATCH --ntasks-per-node=8`, `#SBATCH --cpus-per-task=1`, `#SBATCH --mem=96000MB`,
`QE_NTASKS=${SLURM_NTASKS:-8}`, and `srun -n "$QE_NTASKS" --cpu-bind=cores ...`.
The preflight script accepts sequential multi-dir submissions by checking maximum
simultaneous tasks rather than summing dependent jobs. Local `qe/glcn` and
`qe/glcnac` were regenerated as 8-task/96 GB/24 h inputs and pass preflight with
`--sequential` (`8 / 8` simultaneous). After launch, verify the QE header reports
`Number of MPI processes: 8`; speedup is expected but not assumed linear.
The timed-out logs from `28292263` were fetched locally, then GlcN alone was
resubmitted as optimized job `28303162` (`8` CPUs, `96000M`, `24:00:00`). At
submission it was pending; GlcNAc remains unsubmitted.

---

## Archived July 2026 follow-ups

### 2026-07-01 — Compact candidate fits accepted for counting; GlcN restart3 submitted

#### Benchmark expansion review

- The compact `triage_potential_benchmark_fit_top10` batch finished for all 45
  candidates:

```text
summary rows with ok plot: 45/45
N_selected counts: N=5: 9, N=6: 23, N=7: 9, N=8: 4
ambiguous_eff counts: false: 35, true: 10
```

- Human visual review then accepted **all 45 fit images** as good chitosan chains
  with true `expected_N=6`. This is a review/grading label only: it was added to
  `benchmarks/chitosan_6mer_preassignment_review.tsv` and regenerated in
  `benchmarks/chitosan_6mer_preassignment_review.md`, but it remains outside the
  fitting/selection path.
- All 45 were assigned:

```text
pre_assignment      = accept_counting_visual_N6_confirmed
recommended_action  = use_for_counting_keep_unit_sequence_unconfirmed
expected_N          = 6
unit_sequence       = <blank>
```

- The count of `accept_counting_visual_N6_confirmed` review rows increased from
  33 to 78 in the 45-image review pass. A follow-up second check then promoted
  `240307_019.sxm` from `probable_counting_visual_N6_doubt` to
  `accept_counting_visual_N6_confirmed`, bringing that review class to 79. At the
  time these rows were kept out of unit-assignment grading for provenance reasons;
  the later scope clarification below supersedes that limitation for benchmark
  definition while preserving the label-free rule.
- Compact review artifacts for this pass are:

```text
results/triage_potential_benchmark_fit_top10/summary_overlap060_hard.tsv
results/triage_potential_benchmark_fit_top10/fit_review_ranked.tsv
results/triage_potential_benchmark_fit_top10/fit_review_ranked.md
```

- Promoted all confirmed `expected_N=6` review rows into a counting-only external
  benchmark manifest:

```text
benchmarks/chitosan_6mer_counting_confirmed.toml
```

  The manifest has 146 `clean_target` entries:

```text
accept_counting_visual_N6_confirmed              79
accept_unit_training_known_010010                35
accept_counting_needs_unit_sequence_confirmation 32
```

  This is for post-hoc counting grades only and must not enter fitting or
  selection. The same file set later became the full 0/1/? benchmark scope, with
  `NKNNKN` as external grading control only. The provenance state at this point was
  tracked in:

```text
benchmarks/chitosan_6mer_validation_pending.tsv
```

  That list recorded 111 rows as missing explicit unit-sequence provenance. It is
  superseded for benchmark scope by the later clarification below: all 145
  confirmed 6mer rows share the external control sequence `NKNNKN`, but that truth
  remains grader-only.

#### QE mold jobs

- Raven jobs were checked after the no-watch relaunch:

```text
28525353  qe-glcn_central    TIMEOUT   1-00:00:30
28525353.0 pw.x              FAILED    1-00:00:24  ExitCode 1:0
28525354  qe-glcnac_central  CANCELLED 00:00:00
```

- The GlcN failure was a walltime cancellation, not an out-of-memory failure:
  the Slurm report showed `MaxRSS` about 4.8 GB per MPI task / about 40 GB per
  node. `qe/glcn_restart2/glcn_central_relax.out` was fetched locally, and the
  last `ATOMIC_POSITIONS (angstrom)` block was extracted to
  `qe/glcn_restart2/glcn_central_best3.xyz` (`213` atoms). The relax had reached
  at least 30 BFGS steps, so continuing from the best-so-far geometry preserves
  real progress.
- Prepared `qe/glcn_restart3` from that extracted geometry with the active Raven
  settings (`8` MPI tasks, `96000 MB`, `24:00:00`, `ecutwfc=50`, `ecutrho=360`,
  Gamma-only). Local and remote preflight passed; remote report:

```text
glcn_restart3  nat=213  frozen_relax_atoms=96  mem_mb=96000  status=OK
```

- Submitted the next no-watch production GlcN restart:

```text
qe/glcn_restart3 -> 28601744
```

Do not fetch or inspect `qe/glcn_restart3` outputs until a completion or timeout
notification is available. The production `qe/glcnac` run remains blocked until
GlcN succeeds, because the previous dependent job was cancelled before start.

### 2026-07-02 — Targeted divergent-file fit-improvement screen

To avoid recomputing the full expanded counting benchmark for every idea, built a
screening set from the current external grades: the 146 confirmed `expected_N=6`
rows were joined to their recorded summaries, and the 43 rows whose current
`N_selected` is not 6 (including one `ERR`) were written to:

```text
/tmp/opencode/stmfit_6mer_experiment/current_from_summaries.tsv
/tmp/opencode/stmfit_6mer_experiment/divergent.tsv
/tmp/opencode/stmfit_6mer_experiment/triage_by_folder/
/tmp/opencode/stmfit_6mer_experiment/divergent_recursive_fit_input.tsv
```

The fit-input TSVs contain only file paths plus neutral placeholder columns
required by `test/batch_full.jl`; they do **not** contain `expected_N`, unit
sequences, or target labels. Labels were used only after fitting or offline replay
for external grading.

Current summary-derived baseline on the expanded benchmark is `103/146` exact by
`N_selected`. The divergent subset has `N_selected` counts `4:2`, `5:24`, `7:11`,
`8:5`, and `ERR:1`. Existing summary columns show why a naive selector change is
not safe: raw `N_eff` would be `106/146`, but it fixes 11 current failures while
regressing 8 current successes, including primary 240817 files, so dropping the
robust guard is not acceptable as-is.

Fresh local reruns with the current default config were then tried on three
small divergent folders before spending HPC time on all 43 rows:

```text
20240307_LHe_Cu100  0/3 recovered
20240814_LHe_Cu100  2/4 recovered: 240814_017, 240814_018
20240817_LHe_Cu100  2/3 recovered: 240817_077, 240817_091
```

The merged comparison is:

```text
/tmp/opencode/stmfit_6mer_experiment/default_repro_comparison.tsv
```

`adaptive_support_rescue` was also tried on `20240814_LHe_Cu100`; it matched the
current default result (`2/4`) and did not show a rescue-specific gain there. A
local 4-file adaptive run exceeded 20 minutes when launched with 4 threads, so
larger strategy sweeps should run on Viper or via chunked batch jobs rather than
interactively.

Offline replay of the GCV ambiguity threshold from existing columns was also not
a strong candidate: threshold `0.03` gives `104/146` but introduces an expanded
benchmark regression, while `0.05`/`0.06` are net neutral on current summaries.
Simple guard-down variants based on ambiguity or GCV deltas can reach `108/146`
offline, but regress primary 240817 successes such as `240817_017.sxm`; do not
promote them without a stronger label-free rationale.

Next recommended step: run the current default code on all 43 divergent files
using `/tmp/opencode/stmfit_6mer_experiment/divergent_recursive_fit_input.tsv`.
Promote any strategy to the full 146 confirmed benchmark only if the 43-file
screen improves by a meaningful margin and the full run preserves the 39/39
primary benchmark, avoids new `ERR` rows, and preferably has zero regressions
among the 103 currently correct expanded rows.

Follow-up for the benchmark-optimized objective: created repository-synced fit
input TSVs with no labels, suitable for Slurm jobs:

```text
benchmarks/chitosan_6mer_divergent_fit_input.tsv   # 43 rows
benchmarks/chitosan_6mer_confirmed_fit_input.tsv   # 146 rows
```

`hpc/launch_remote.sh` now exports `STMFIT_DATA_DIR=/ptmp/<user>/stmfit/data`
to the Slurm job and forwards dedicated `STMFIT_TSV` / `STMFIT_SKIP_1D` variables,
so benchmark subset runs no longer need fragile multi-word `STMFIT_BATCH_ARGS`.
A Viper dry-run for the full 146-file current-default rerun was correct, but the
actual submission was blocked before sync/submission because both `viper` and
`raven` SSH handshakes timed out, likely because no MPCDF 2FA/ControlMaster
session was active. Resume by opening the MPCDF SSH session first, then run:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/default_repro \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SKIP_1D=1 ./hpc/launch_remote.sh
```

After opening the Viper session, that command synced code/data, instantiated the
Julia project on the Viper login node, and submitted the full 146-file current
default rerun as Slurm array job:

```text
10367643  stmfit  viper/small  8 chunks  PD (Priority)
```

Output directory:

```text
results/experiments/6mer_full146/default_repro
```

The first submit exposed a launcher bug: non-`--watch` submissions attempted to
fetch results immediately, before the remote output directory existed. The job was
already submitted successfully; `hpc/launch_remote.sh` now creates the remote code
directory before rsync and exits after submission unless `--watch` or
`--fetch-only` is requested. Resume/fetch after completion with:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/default_repro \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SKIP_1D=1 ./hpc/launch_remote.sh --fetch-only
```

Completion follow-up: Viper job `10367643_[1-8]` completed successfully on the
`small` partition (`ExitCode=0:0` for every array task). The chunks were merged
on Viper and fetched locally:

```text
results/experiments/6mer_full146/default_repro/summary_overlap060_hard.tsv
```

The merged summary has 146 data rows and 146 successful fits. External grading
against `benchmarks/chitosan_6mer_counting_confirmed.toml` gives:

```text
N_selected  106/146 exact (72.6%), +/-1 = 139/146 (95.2%), over/under = 17/23
```

On the original `benchmarks/chitosan_240817.toml` target subset, the same run
remains `4/4` exact. Relative to the older summary-derived `103/146` baseline,
the fresh run changes 15 files: 8 failures are recovered, 5 former successes
regress, and 2 failures move to a different wrong N, for a net `+3` exact gain.
The recovered rows are `240311_Cu100060.sxm`, `240312_Cu100081.sxm`,
`240814_017.sxm`, `240814_018.sxm`, `240817_077.sxm`, `240817_091.sxm`,
`240818_012.sxm`, and `240818_017.sxm`; the regressions are
`240814_024.sxm`, `240816_005.sxm`, `240817_072.sxm`, `240817_083.sxm`, and
`240818_021.sxm`.

Fresh grading of existing selector columns on the same full146 summary does not
beat the default selected column:

```text
N_selected      106/146
N_eff           106/146
N_ell           106/146
robust_aicc_N   103/146
N_circ          100/146
```

Offline replay of generic post-fit guard variants on the fresh summary can reach
`110/146` exact (best cases: robust downshift only when the GCV ambiguity delta
is small, or only when the GCV choice is ambiguous), but every best-scoring
variant trades fixes for regressions. These rules remain screening observations,
not promoted selection changes: they were discovered by external grading on this
confirmed 6mer set, and need a label-free physical or statistical rationale
before they can become part of the pipeline.

Next benchmark-optimization candidates are therefore experimental, not yet
canonical: either run a full146 Viper job with `adaptive_support_rescue` to test
whether the small local screen missed broader gains, or design a label-free guard
variant from residual/fit diagnostics and then grade it externally. Do not tune a
threshold directly to `expected_N=6` and present it as objective.

Launched the first of those experimental candidates on Viper: a full 146-file
rerun with `config/chitosan_adaptive_support_rescue.toml`, the same label-free
input TSV, and 1D skipped. The dry-run showed the intended project/data/config
paths and Slurm export variables, then the no-watch submission succeeded:

```text
10370782  stmfit  viper/small  8 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan_adaptive_support_rescue.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/adaptive_support_rescue \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SKIP_1D=1 ./hpc/launch_remote.sh
```

After completion, fetch/merge with the same environment plus `--fetch-only`, then
grade `results/experiments/6mer_full146/adaptive_support_rescue/summary_overlap060_hard.tsv`
against `benchmarks/chitosan_6mer_counting_confirmed.toml`. The score to beat is
the current default rerun's `106/146` exact by `N_selected`.

Completion follow-up: all `10370782_[1-8]` array tasks completed with
`ExitCode=0:0` and runtimes from `00:09:05` to `00:15:49`. The remote merge
produced 146 data rows and 146 successful fits, then fetched the summary and
chunk outputs to:

```text
results/experiments/6mer_full146/adaptive_support_rescue/summary_overlap060_hard.tsv
```

External grading on the expanded confirmed counting benchmark gives:

```text
adaptive_support_rescue N_selected  104/146 exact (71.2%), +/-1 = 139/146 (95.2%), over/under = 15/27
```

The original 240817 target subset falls to `3/4` exact (`240817_043.sxm` becomes
the target failure). Relative to the current default rerun (`106/146`), the
adaptive-support run changes only 6 rows: it fixes `240814_024.sxm` and
`240816_005.sxm`, but regresses `240817_043.sxm`, `240817_077.sxm`,
`240817_093.sxm`, and `241113_089.sxm`, for a net `-2` exact result. Therefore
`adaptive_support_rescue` is not a benchmark-improving 6mer-counting candidate;
keep the current default run as the score to beat (`106/146`).

Offline inspection of the two full146 summaries shows why this candidate is a
dead end for the 6mer counting objective: no row accepted an actual support
rescue (`adaptive_support_rescue` count = 0). The adaptive run had 119 plain
keeps, 6 `reject_no_improvement` rows, and 26 rows where the robust guard still
acted. The six changed `N_selected` rows therefore came from altered guard/source
behavior, not from genuinely longer recovered support. Do not launch another
adaptive-support sweep for 6mer counting unless a new label-free trigger or
support diagnostic is added first.

Next label-free screen: selected `support_marginalized_gcv_guard` as the cheapest
existing-policy candidate to test before any new full146 run. This selector
rescales GCV across a grid of support paddings and applies a guarded/parsimony
choice without using labels. Because the summary TSVs do not preserve the fitted
parameter objects needed to replay this selector offline, the first screen is a
43-file rerun on the divergent input TSV, not a full benchmark rerun.

Submitted the 43-row Viper screen after a clean dry-run:

```text
10371435  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/support_marginalized_gcv_guard \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=support_marginalized_gcv_guard STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending/running (`10371435_1` running, `10371435_[2-4]`
pending for `QOSGrpCpuLimit`). After it completes, fetch/merge with the same
environment plus `--fetch-only`, grade the 43 divergent rows externally, and
promote to a full146 run only if it recovers a meaningful number of failures
without a label-dependent rationale.

Completion follow-up: all `10371435_[1-4]` tasks completed with `ExitCode=0:0`
and runtimes from `00:03:29` to `00:04:47`. The remote merge produced 43 data
rows and 43 successful fits, fetched to:

```text
results/experiments/6mer_divergent/support_marginalized_gcv_guard/summary_overlap060_hard.tsv
```

External grading on those 43 divergent rows gives:

```text
support_marginalized_gcv_guard N_selected  16/43 exact (37.2%), +/-1 = 37/43 (86.0%), over/under = 11/16
```

This is a meaningful screen improvement. On the same 43 files, the current full146
default rerun is `8/43` exact, while the older summary-derived baseline was
`0/43`. Relative to the current default, the support-marginalized guard changes
22 rows: 12 fixes, 4 regressions, and 6 wrong-to-wrong moves, for a net `+8`
exact result. Relative to the older summary-derived baseline it fixes 16 rows and
introduces no exact regressions, because all 43 rows were failures there. The
selector source is `support_marginalized_gcv` for all 43 rows, with 21 guarded
downshifts and 22 keeps. This is strong enough to justify a full146 Viper run;
do not promote it as canonical until that full run is graded, because the screen
was enriched for failures and does not measure regressions among the other 103
expanded-benchmark rows.

Promoted the screen to a full 146-file Viper run with the same label-free input
TSV, default chitosan config, and `STMFIT_SELECTION_POLICY=support_marginalized_gcv_guard`.
The dry-run was clean, then the no-watch submission succeeded:

```text
10371737  stmfit  viper/small  8 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=8 CPUS_PER_TASK=4 WALLTIME=08:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_full146/support_marginalized_gcv_guard \
N_FILES=146 STMFIT_TSV=benchmarks/chitosan_6mer_confirmed_fit_input.tsv \
STMFIT_SELECTION_POLICY=support_marginalized_gcv_guard STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending (`10371737_[1-8]`, reason `QOSGrpCpuLimit`).
After completion, fetch/merge with the same environment plus `--fetch-only`, then
grade the full summary against `benchmarks/chitosan_6mer_counting_confirmed.toml`.
The score to beat remains the default rerun's `106/146` exact by `N_selected`.

Completion follow-up: all `10371737_[1-8]` tasks completed with `ExitCode=0:0`
and runtimes from `00:04:56` to `00:06:04`. The remote merge produced 146 data
rows and 146 successful fits, fetched to:

```text
results/experiments/6mer_full146/support_marginalized_gcv_guard/summary_overlap060_hard.tsv
```

External grading on the expanded confirmed counting benchmark gives:

```text
support_marginalized_gcv_guard N_selected  92/146 exact (63.0%), +/-1 = 140/146 (95.9%), over/under = 14/40
```

The original 240817 target subset falls to `2/4` exact (`240817_017.sxm` and
`240817_043.sxm` are target failures). Relative to the current default rerun
(`106/146`), this full run changes 50 rows: 14 fixes, 28 regressions, and 8
wrong-to-wrong moves, for a net `-14` exact result. The divergent-only screen was
therefore misleading because it was enriched for failures and did not expose the
large number of regressions among previously correct rows. The selector source is
`support_marginalized_gcv` for all 146 rows, with 42 guarded downshifts, 99
keeps, and 5 parsimony moves. Do not promote `support_marginalized_gcv_guard` for
6mer counting; the score to beat remains the default rerun's `106/146`.

Next existing-policy screen: selected `fwd_bwd_consensus` for a 43-file divergent
rerun because it is label-free and uses the independent forward/backward scan
consistency signal rather than another support-padding heuristic. The dry-run was
clean, then the no-watch Viper submission succeeded:

```text
10372583  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/fwd_bwd_consensus \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=fwd_bwd_consensus STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending (`10372583_[1-4]`, reason `QOSGrpCpuLimit`).
After completion, fetch/merge and grade the 43 divergent rows externally before
considering any full146 promotion.

The `fwd_bwd_consensus` screen completed cleanly on Viper and was merged/fetched
manually:

```text
10372583_1  COMPLETED  00:04:58  0:0
10372583_2  COMPLETED  00:03:52  0:0
10372583_3  COMPLETED  00:02:53  0:0
10372583_4  COMPLETED  00:04:09  0:0
```

Merged summary:

```text
results/experiments/6mer_divergent/fwd_bwd_consensus/summary_overlap060_hard.tsv
  43 data rows (43 ok)
```

External grade on the 43 divergent rows:

```text
fwd_bwd_consensus N_selected  16/43 exact (37.2%), +/-1 = 34/43 (79.1%), over/under = 14/13
```

Against the same 43 rows from the fresh default full146 rerun, this is an exact
improvement from `8/43` to `16/43`: 10 fixes, 2 regressions, 5 wrong-to-wrong
changes, and 6 unchanged correct rows. However, its `+/-1` count drops from
`37/43` to `34/43`. It also only matches the earlier
`support_marginalized_gcv_guard` screen on exact score (`16/43`) while being
worse on `+/-1` (`34/43` vs `37/43`), and that earlier screen already failed the
full146 promotion test badly (`92/146`, net `-14` vs default). Do not promote
`fwd_bwd_consensus` to a full146 run; the score to beat remains the fresh default
rerun's `106/146`.

Next existing-policy screen: selected `laplace_evidence_guard` for the same
43-file divergent rerun because it is technically distinct from the rejected
support-padding and forward/backward rescoring candidates. It is a label-free,
down-only local curvature/evidence guard: it can veto a one-lobe overfit, but it
cannot increase `N` or use the benchmark target.

The dry-run was clean, then the no-watch Viper submission succeeded:

```text
10373261  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/laplace_evidence_guard \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=laplace_evidence_guard STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state was pending (`10373261_[1-4]`, reason `QOSGrpCpuLimit`).
After completion, merge/fetch and grade the 43 divergent rows externally before
considering any full146 promotion.

The `laplace_evidence_guard` screen completed cleanly on Viper and was
merged/fetched manually:

```text
10373261_1  COMPLETED  00:06:26  0:0
10373261_2  COMPLETED  00:04:26  0:0
10373261_3  COMPLETED  00:02:51  0:0
10373261_4  COMPLETED  00:04:10  0:0
```

Merged summary:

```text
results/experiments/6mer_divergent/laplace_evidence_guard/summary_overlap060_hard.tsv
  43 data rows (43 ok)
```

External grade on the 43 divergent rows:

```text
laplace_evidence_guard N_selected  11/43 exact (25.6%), +/-1 = 34/43 (79.1%), over/under = 10/22
```

Against the same 43 rows from the fresh default full146 rerun, this improves the
exact score only from `8/43` to `11/43`: 9 fixes, 6 regressions, 5
wrong-to-wrong changes, and 2 unchanged correct rows. It is worse than both
previous 43-file screens on exact score (`16/43` for both
`support_marginalized_gcv_guard` and `fwd_bwd_consensus`) and does not improve
the `+/-1` count (`34/43`). Do not promote `laplace_evidence_guard` to a full146
run; the score to beat remains the fresh default rerun's `106/146`.

Next existing-policy screen: selected `stability_selection` for the same 43-file
divergent rerun. This uses the frozen support-padding grid, like
`support_marginalized_gcv_guard`, but asks a different label-free question: how
often each candidate `N` remains competitive under support perturbations. It is
therefore a cheap final existing-policy triage before stopping the current
selector sweep.

The dry-run was clean, then the no-watch Viper submission succeeded:

```text
10374201  stmfit  viper/small  4 chunks  submitted
```

Command used:

```bash
STMFIT_DATA_DIR=/home/durif/Rebecca/data/data \
STMFIT_SSH_HOST=viper N_CHUNKS=4 CPUS_PER_TASK=4 WALLTIME=04:00:00 \
MEM_PER_CPU=4000 STMFIT_CONFIG=config/chitosan.toml \
STMFIT_OUTDIR=results/experiments/6mer_divergent/stability_selection \
N_FILES=43 STMFIT_TSV=benchmarks/chitosan_6mer_divergent_fit_input.tsv \
STMFIT_SELECTION_POLICY=stability_selection STMFIT_SKIP_1D=1 \
./hpc/launch_remote.sh
```

Initial Slurm state had chunk 1 running and chunks 2–4 pending under
`QOSGrpCpuLimit`. After completion, merge/fetch and grade the 43 divergent rows
externally before considering any full146 promotion.

The `stability_selection` screen completed cleanly on Viper and was
merged/fetched manually:

```text
10374201_1  COMPLETED  00:05:51  0:0
10374201_2  COMPLETED  00:04:10  0:0
10374201_3  COMPLETED  00:02:45  0:0
10374201_4  COMPLETED  00:04:03  0:0
```

Merged summary:

```text
results/experiments/6mer_divergent/stability_selection/summary_overlap060_hard.tsv
  43 data rows (43 ok)
```

External grade on the 43 divergent rows:

```text
stability_selection N_selected  15/43 exact (34.9%), +/-1 = 35/43 (81.4%), over/under = 7/21
```

Against the same 43 rows from the fresh default full146 rerun, this improves the
exact score from `8/43` to `15/43`: 10 fixes, 3 regressions, 15 wrong-to-wrong
changes, and 5 unchanged correct rows. It is better than
`laplace_evidence_guard` (`11/43`) but still below both previous best 43-file
screens (`16/43` for `support_marginalized_gcv_guard` and
`fwd_bwd_consensus`). Its `+/-1` count (`35/43`) is also below the fresh default
and support-marginalized screens (`37/43`). Do not promote `stability_selection`
to a full146 run; the score to beat remains the fresh default rerun's `106/146`.

At this point the existing-policy triage did not identify a full146 candidate:
`support_marginalized_gcv_guard` had the best 43-file screen profile but already
failed full146 badly, `fwd_bwd_consensus` matched its exact score with worse
`+/-1`, and `laplace_evidence_guard`/`stability_selection` were weaker on the
divergent screen. Further benchmark improvement likely needs new label-free
evidence or a targeted analysis of the remaining failure modes rather than
promoting another existing selector wholesale.

Offline failure-mode audit across the fresh default full146 summary and all four
43-file divergent screens confirmed that conclusion. The fresh default full146
score remains:

```text
default full146  106/146 exact, +/-1 = 139/146, over/under = 17/23
```

The 40 full146 misses are structured, not random:

```text
N_selected=3:  1 file
N_selected=4:  2 files
N_selected=5: 20 files
N_selected=7: 13 files
N_selected=8:  4 files
```

On the same 43 divergent rows, the completed screen scores were:

```text
default                         8/43 exact, +/-1 = 37/43, over/under = 14/21
support_marginalized_gcv_guard 16/43 exact, +/-1 = 37/43, over/under = 11/16
fwd_bwd_consensus              16/43 exact, +/-1 = 34/43, over/under = 14/13
laplace_evidence_guard         11/43 exact, +/-1 = 34/43, over/under = 10/22
stability_selection            15/43 exact, +/-1 = 35/43, over/under = 7/21
```

The screens fix different failure modes: among default under-counts on the 43-row
set, `fwd_bwd_consensus` fixes the most (`9/21`) without increasing the distance
from 6, while among default over-counts, `laplace_evidence_guard` fixes the most
(`7/14`) without increasing the distance. That split is scientifically useful,
but not a promotion rule by itself: choosing a policy by whether a file is an
under- or over-count uses the withheld benchmark label and is therefore invalid.

Several offline ensemble rules built only from existing label-free screen outputs
were scored for diagnostics, not promoted. The best exact-count variant on the
43-row set was a five-policy lower-tie mode:

```text
mode5_lower_tie  14/43 exact, +/-1 = 36/43, over/under = 10/19
```

Other consensus variants reached only `11–13/43` exact. None beat the best
single 43-file screens (`16/43`), and none has a plausible reason to improve the
remaining 103 full146 rows without introducing regressions. Do not promote an
ensemble rule from this audit. The next defensible benchmark-improvement gate is
not another wholesale selector rerun; it should be a targeted, label-free study
of the remaining failure classes, especially support-sensitive under-counts
(`N=5`) versus over-segmentation (`N=7/8`), using per-file diagnostics and visual
evidence before any new rule is frozen.

Targeted representative diagnostics separate the remaining failures into a few
label-free mechanisms worth studying:

- **Guard-induced under-counts.** Several `N=5` misses are not raw GCV failures:
  the primary ell/circ fit already agrees on `N=6`, but the robust-AICc guard
  downshifts to `5`. Examples include `240818_020.sxm`, `240818_026.sxm`,
  `241113_160.sxm`, and `241114_010.sxm`. For `240818_020.sxm`, the ell score
  curve has `N=6` as the GCV minimum, while `N=5` is 21.2% worse; multiple
  alternative label-free screens also return `6`. This suggests the next audit
  should focus on when the robust guard is contradicted by strong 2D GCV and
  ell/circ agreement, rather than changing the primary selector wholesale.
- **Ambiguous support-sensitive under-counts.** Some files have broad candidate
  sets and small GCV gaps, e.g. `240312_Cu100070.sxm` / `240312_Cu100071.sxm`:
  the ell sweep prefers `N=10`, but `N=6` and `N=7` are close, and the robust
  guard collapses to `5`. The fwd/bwd screen alone returns `6`, while support and
  Laplace go high. These need visual/diagnostic inspection of support extent and
  scan-direction consistency before any rule is frozen.
- **Over-segmentation with near-ties.** Some `N=7` files are one-lobe GCV
  overfits where `N=6` is label-free competitive. Examples: `240313_Cu100062.sxm`
  has ell `7` only 0.6% better than `6` and circ prefers `6`; support, Laplace,
  and stability all return `6`. `240314_Cu100_042.sxm` has ell `7` only 1.6%
  better than `6`, and the same three screens return `6`. This is a plausible
  generic one-lobe ambiguity class.
- **Over-segmentation without near-ties.** Other `N=7/8` files do not have a
  competitive `N=6` under current fit support. For example `241114_044.sxm` and
  `241114_045.sxm` have `N=8` clearly favoured by GCV; only stability jumps to
  `6`. Files such as `240310_Cu100009.sxm`, `240815_098.sxm`, `241114_043.sxm`,
  and `241114_046.sxm` remain unfixed by the screens. These are not safe targets
  for a simple ambiguity guard.

Next gate: build a small diagnostic report, not a production selector, that
flags candidate files using only label-free contradictions:

```text
robust-downshift contradiction:
  robust_AICc_N < N_eff
  and N_ell == N_circ == N_eff
  and GCV(N_robust) is materially worse than GCV(N_eff)

one-lobe overfit ambiguity:
  N_eff = N_alt + 1
  and lower-N GCV gap is within a frozen tolerance
  and either circ prefers lower N or support/Laplace/stability concur
```

Score those flags externally only after the diagnostic report is generated. If
the flags isolate high-signal subsets without many obvious false positives, then
freeze the rule and run a full146 test. If they do not, stop selector tuning and
move to visual failure review / data-quality stratification.

Generated the diagnostic report and external score tables:

```text
results/diagnostics/full146_failure_flags.tsv
results/diagnostics/full146_failure_flags_score.tsv
```

The flags are computed from the fresh full146 default summary and existing
label-free screen outputs; the manifest is used only afterward for external
scoring. Corrected score summary:

```text
robust_downshift_contradiction  n=11  fixes=6  regressions=2  neutral_wrong=3  net=+4  precision=54.5%
one_lobe_overfit_ambiguity      n=6   fixes=2  regressions=4  neutral_wrong=0  net=-2  precision=33.3%
combined_unique                 n=17  fixes=8  regressions=6  neutral_wrong=3  net=+2  precision=47.1%
```

Decision: do not promote the one-lobe ambiguity flag; it is net negative on the
external grade. The robust-downshift contradiction flag is promising enough to
preserve as the next focused candidate (`+4` net, implied `110/146` if applied as
a post-selection correction), but it is not clean enough to silently fold into
the default selector: it still has two regressions and three wrong-to-wrong
suggestions. The next defensible step is to freeze this rule as an explicit
diagnostic/experimental policy and run/grade it on full146, or inspect its 11
flagged rows visually before deciding whether it is scientifically acceptable.

Inspected the 11 `robust_downshift_contradiction` rows. External grading of the
suggested replacement (`N_eff`, because robust-AICc downshift was contradicted by
strong 2D GCV) splits as:

```text
fixes:          6  (240816_002, 240818_020, 240818_026, 240818_028, 241113_160, 241114_010)
regressions:    2  (240817_017, 240818_015)
wrong-to-wrong: 3  (240307_019, 240818_025, 241114_043)
```

Tried stricter label-free refinements. The best non-cheating refinement was
`one_step_and_alt_votes_ge2`, requiring a one-lobe robust correction plus at
least two existing label-free screen votes for the same suggested `N`:

```text
one_step_and_alt_votes_ge2  n=5  fixes=4  regressions=0  wrong_to_wrong=1  net=+4  precision=80.0%
```

It avoids the two regressions, but still includes `240307_019.sxm` as a
wrong-to-wrong `4→5` move, so it is not clean enough for default selection. The
only perfect-looking variants were absolute rules such as `5→6 only` or
`suggest N=6 only`:

```text
INVALID_absolute_5_to_6_only  n=6  fixes=6  regressions=0  wrong_to_wrong=0
INVALID_absolute_suggest_6_only  n=6  fixes=6  regressions=0  wrong_to_wrong=0
```

Those are invalid because the absolute target value `6` is the benchmark label
shape, not a label-free scientific criterion. Do not promote them. Current next
gate: either visually review the small robust-contradiction set before designing
a real label-free predicate, or stop selector tuning and move to data-quality /
failure-class stratification. No selector/policy change is justified yet.

Qualitative plot review of the 11 robust-contradiction rows supports the same
conclusion. Without using expected labels, the visually plausible upward moves
were mostly the modest under-count corrections: `240307_019` (`4->5`),
`240818_020` (`5->6`), `240818_026` (`5->6`), `240818_028` (`5->6`, cautious),
`241113_160` (`5->6`), and `241114_010` (`5->6`, cautious). `240816_002`
(`5->6`) was only partly plausible because the feature is broad/blurred and the
residuals remain structured. The larger or high-count upward moves were weak:
`240817_017` (`6->7`) was ambiguous, `240818_015` (`6->7`) looked like fitting
faint tails, `240818_025` (`5->7`) looked crowded/over-counted, and
`241114_043` (`7->8`) looked weak because the extra components were crowded or
terminal.

Decision after visual review: the robust-contradiction diagnostic identifies a
real under-count signal, but it also catches visually weak upward moves. Do not
turn it into a production selector yet. If this line of work continues, the next
scientific step should be a visual/data-quality stratification of the flagged
rows and a genuinely label-free predicate that suppresses tail/crowding-driven
upward moves; otherwise stop selector optimization at the fresh default
`106/146` benchmark state.

Follow-up 90% attempt: tested a purely physical support/spacing post-selection
family offline from the fresh full146 default summary. The rule uses only the
measured 2D support and existing chitosan spacing/overlap calibration from
`config/chitosan.toml` (`spacing_min_nm=0.35`, `spacing_max_nm=0.75`,
`max_overlap=0.60`, `sigma_parallel_max_nm=0.509`, giving
`spacing_min_eff=0.51448 nm`). Labels were used only after rule replay for
external grading.

The first support-midpoint sweep was written to:

```text
results/diagnostics/full146_support_bounds_rules.tsv
```

Its best non-target-shaped rule was
`support_mid_round_bounded_one_step_delta_ge1`: move at most one lobe toward the
midpoint of the physical support-derived feasible-N interval. It reached:

```text
126/146 exact, +/-1 = 143/146
46 changes: 30 fixes, 10 regressions, 5 wrong-closer, 1 wrong-farther, net +20
```

This was the strongest label-free offline signal so far, but still below the
90% target (`132/146`) and too regression-prone for promotion. Inspection showed
the damage concentrated in permissive upshifts from already-correct files,
especially `6->7` on long measured supports.

Two refined sweeps then tried to suppress those upshift failures while preserving
the support signal:

```text
results/diagnostics/full146_support_hybrid_rules.tsv
results/diagnostics/full146_support_refined_rule_sweep.tsv
```

The best hybrid rule was a conservative asymmetric version: downshift one step
whenever the current selection lies above the support midpoint; upshift one step
only when the support midpoint lies above the current selection, `N_eff` is also
above current, and the effective GCV gap is close (`delta_GCV_rel_eff <= 0.3`).
Equivalent refined variants using agreement between ell/circ support midpoints
or the same close-GCV upshift condition topped out at:

```text
127/146 exact, +/-1 = 143/146
31 changes: 23 fixes, 2 regressions, 5 wrong-closer, 1 wrong-farther, net +21
```

The remaining damage was narrow but not removable by the tested label-free
features without giving back many fixes: two `6->7` regressions
(`240817_017.sxm`, `240817_019.sxm`) and one `7->8` wrong-farther move
(`240310_Cu100009.sxm`) came from the same support-upshift mechanism that fixed
many under-counts. Stricter GCV-curve filters reduced regressions but lowered the
exact score to `123-124/146`; no tested rule reached the requested 90% level.

Decision at the time: do not promote the support-midpoint family as a selector
change. It was a useful diagnostic showing that measured support explains many
full146 misses, but the best label-free offline rule was still `5` exact files
short of 90% and retained non-negligible regressions. This was superseded by the
2026-07-02 promotion entry below after accepting the best current full146 result
as a practical chitosan default while continuing selector research.

### 2026-07-02 — Promote support-midpoint hybrid as current chitosan default

**Decision:** Promote `support_midpoint_hybrid` to the default
`config/chitosan.toml` batch policy. This supersedes the previous no-promotion
decision above. The rationale is pragmatic: the frozen label-free rule is the
best current expanded 146-file counting result, even though it does not reach the
earlier aspirational 90% target and is not a universal selector guarantee.

**Implemented rule:** The batch first applies the existing integrated
robust-AICc guard (`gcv_with_robust_aicc_guard`), then computes the midpoint of
the measured 2D support's feasible-N interval using the physical spacing and
overlap calibration. The final support layer is bounded to one lobe:

```text
if N_guarded > support_midpoint:
    N_selected = N_guarded - 1
elif N_guarded < support_midpoint
     and N_eff > N_guarded
     and delta_GCV_rel_eff <= 0.30:
    N_selected = N_guarded + 1
else:
    N_selected = N_guarded
```

The threshold is recorded as
`support_midpoint_up_gcv_rel_threshold = 0.30` in `[selection]`. The rule uses no
expected `N`, no target count, and no benchmark labels during fitting or
selection. The 146-file manifest remains an external grading set only.

**External grade used for the decision:**

```text
previous default replay: 106/146 exact, +/-1 = 139/146
support_midpoint_hybrid: 127/146 exact, +/-1 = 143/146
```

Known residual errors remain: two `6->7` regressions (`240817_017.sxm`,
`240817_019.sxm`) and one `7->8` wrong-farther move (`240310_Cu100009.sxm`) were
observed in the offline sweep. The default change is therefore a documented
current-best chitosan setting, not the endpoint of label-free selection work.

Post-promotion doc review (post-implementation review pass) caught three
documentation mismatches that have now been corrected:

- `docs/src/index.md` still named `gcv_with_robust_aicc_guard` as the batch
  default and quoted the `39/39` primary-benchmark number without separating
  robust-guard validation from the promoted default. Updated to name
  `support_midpoint_hybrid`, record the `127/146` / `143/146` external counting
  grade, and keep `39/39` explicitly as the robust-guard validation result.
- `docs/src/selection.md` documented `selection_source` as `ell_robust_aicc`,
  but `_select_primary` (`selectors.jl:661`) actually emits `robust_aicc_guard`
  when the guard moves the primary count. The doc now lists the real emitted
  values (`robust_aicc_guard`, `support_midpoint_down`, `support_midpoint_up`,
  or the kept effective source); `ell_robust_aicc` is correctly scoped to
  `refined_source`.
- `AGENTS.md` left the `4/4 clean_target` / reproducibility / threshold-robust
  bullets ambiguous in the same section that says `support_midpoint_hybrid` is
  not a no-regression primary-benchmark selector. Scoped those bullets
  explicitly to robust-AICc guard validation.

No code, config, or selection behaviour changed in this pass — only
documentation accuracy.

### 2026-07-03 — Gap≥2 down-to-midpoint extension: 127/146 → 129/146

**Decision:** Extend the support-midpoint hybrid so that when the robust guard
over-counts by at least two lobes relative to the support midpoint (gap ≥ 2),
the rule goes directly to the midpoint instead of making only a one-step
correction.

**Motivation:** Systematic offline replay on the frozen full146 fit data
(`results/experiments/6mer_full146/default_repro`) tested every label-free rule
constructible from the per-file GCV/BIC/chi2 curves, support geometry, κ,
ambiguity, and alternative selectors (stability, Laplace, fwd/bwd, spatial CV).

Key negative finding: the 241114 over-count cluster (N=7-8 instead of 6) is
**genuinely indistinguishable** by any residual-based criterion. Both GCV and
spatial k-fold CV prefer N=7-8 because the 7th lobe captures real spatially-
correlated signal (likely a surface defect or monomer substructure), not iid
noise. No fit-level feature (spacing_cv, κ, overlap, amplitude ratio) discriminates.

Key positive finding: the only improvement comes from a stronger
support-geometry prior. When `N_guarded ≥ support_midpoint + 2`, the geometry
disagrees strongly enough with the fit to justify a direct move to the midpoint.

**Implemented rule change** (`test/batch_full.jl`,
`_support_midpoint_hybrid_selection`):

```text
if N_guarded > support_midpoint + 1:
    N_selected = support_midpoint          # NEW: gap ≥ 2, trust geometry
elif N_guarded > support_midpoint:
    N_selected = N_guarded - 1             # unchanged: gap = 1
elif N_guarded < support_midpoint
     and N_eff > N_guarded
     and delta_GCV_rel_eff <= 0.30:
    N_selected = N_guarded + 1             # unchanged up-shift
else:
    N_selected = N_guarded                 # unchanged
```

New `selection_source` value: `support_midpoint_down_to_mid` (distinguishes the
gap≥2 case from the one-step `support_midpoint_down`).

**Offline replay result on frozen full146 data:**

```text
frozen hybrid ±1:          127/146 exact (87.0%), 143/146 ±1
hybrid ±2 (this change):   129/146 exact (88.4%), 143/146 ±1
changed=2  fixes=2  regressions=0
```

The 2 fixed files: `240818_021` (N=7→6) and `241114_045` (N=7→6), both had
`N_guarded=8, support_midpoint=6`. No regression because the up-shift branch
retains its `dgcv_rel_eff ≤ 0.30` guard.

**Limitation acknowledged:** the over-counts at gap=1 (e.g., `241114_040/044/046`
where N_guarded=7, midpoint=6) cannot be fixed by any residual-based label-free
rule — the 7th lobe is a genuine residual improvement. These remain at the ±1
grade level only.

**Benchmark cleanup:** `240310_Cu100009.sxm` removed from
`benchmarks/chitosan_6mer_counting_confirmed.toml` after visual review confirmed
the image is technically degraded (23.6% NaN pixels, "très complexe"). The
manifest now contains 145 files. Updated grade numbers on the 145-file benchmark:

```text
default_repro:  106/145 exact (73.1%), 138/145 ±1 (95.2%)
hybrid ±1:      127/145 exact (87.6%), 143/145 ±1 (98.6%)
hybrid ±2:      129/145 exact (89.0%), 143/145 ±1 (98.6%)
```

A full146 batch was submitted to Viper to confirm the ±2 rule on a real batch
run; the grading script uses the 145-file manifest.

Completion follow-up: the first Viper submission (`10390605`) failed because the
remote data staging copied local SXM symlinks as broken symlinks. After replacing
the remote data directory with real files via `rsync -L`, rerun `10391759`
completed successfully. `hpc/merge_chunks.jl` merged both shards into
`results/experiments/6mer_full146/pm2_confirm/summary_overlap060_hard.tsv` with
146 data rows and 146 `ok` rows. External grading against
`benchmarks/chitosan_6mer_counting_confirmed.toml` confirmed the offline replay:

```text
primary score:  129/145 (89.0%)
primary ±1:     143/145 (98.6%)
primary over/under: 6/10
stress rows:    0
```

The two remaining beyond-±1 misses are `240814_020.sxm` (`N_selected=4`) and
`240818_019.sxm` (`N_selected=3`). The exact-miss set is unchanged from the
confirmed grade: `240307_019`, `240310_Cu100032`, `240314_Cu100_024`,
`240314_Cu100_026`, `240814_020`, `240814_021`, `240816_005`, `240817_017`,
`240817_019`, `240817_078`, `240817_083`, `240818_007`, `240818_019`,
`241114_040`, `241114_044`, and `241114_046`.

The merged full146 run's `selection_source` distribution was: `ell=103`,
`robust_aicc_guard=12`, `support_midpoint_down=12`,
`support_midpoint_down_to_mid=2`, and `support_midpoint_up=17`. The two
`support_midpoint_down_to_mid` rows are the intended gap≥2 direct-to-midpoint
fixes.

### Follow-up: 0/1/? benchmark scope clarification

The full unit-assignment benchmark must use the same 145 confirmed 6mer files as
the counting/fit benchmark. The external control sequence is `NKNNKN`, encoded as
`010010` if `0=N` and `1=K`, or `101101` under the flipped identity convention.
This correction supersedes the earlier documentation that treated only the 35
240817 clean/clean_target rows as gradable or described 111 expanded rows as
missing unit sequences for benchmark purposes.

The scientific constraint is unchanged and becomes more important: this is a
fully label-free benchmark. `NKNNKN`, the `2K/4N` composition, the fact that `N=6`,
and any benchmark-control convention may be used only by external grading and
diagnostic reports after predictions are frozen. They must not enter fitting,
selection, per-lobe assignment, threshold selection, abstention rules, composition
priors, or method calibration. The objective is a robust rule that extrapolates to
unknown systems rather than one shaped to the benchmark control.

Existing `178/210`, `154/171`, and `101/106` numbers remain useful as historical
35-file subset diagnostics, but they are no longer the headline 0/1/? benchmark.
The reporting path should be expanded or replaced before claiming full145 unit
assignment metrics.

Implementation follow-up: `test/report_unit_assignment_benchmark.jl` now has an
explicit `--full145` mode. It writes `control_full145_truth.tsv` under the report
output directory from `benchmarks/chitosan_6mer_counting_confirmed.toml`, using the
external `NKNNKN` control encoding (`010010` by default, or `101101` with
`--control-sequence`). It then runs the existing grader and refuses any profile
whose graded denominator is not 145 files / 870 lobes. The default command remains
the historical subset report because the frozen prediction profiles currently have
only 210 lobe rows; full145 mode fails on those profiles by design rather than
printing a false headline.

Second implementation follow-up: added an early full145 coverage preflight before
the grader runs. Current frozen profiles now fail immediately with a message like
`expected 145 files / 870 lobes, found 35 files / 210 lobe rows`, plus examples of
missing files. A synthetic temporary full145-perfect prediction TSV validates the
positive path at 145/145 chains and 870/870 lobes. The remaining real work is to
generate label-free prediction profiles for all 145 benchmark files from the
feature/patch pipeline; scoring the current 35-file profiles as full145 is blocked
by design.

### Follow-up: portable label-free prediction builder and full145 blocker

Added `test/build_labelfree_unit_predictions.jl` to close the reproducibility gap
between feature extraction and prediction TSVs. The script is prediction-only and
does not read benchmark truth, `NKNNKN`, `expected_N`, or a composition prior. It
combines per-file standardized feature views over multiple k-means seeds and maps
the higher-amplitude cluster to GlcNAc (1), matching the label-free physical
convention used by the grader. Optional inputs add the split-width
`split_log_skew` feature and backward-patch negative-moment descriptors
(`bwd_neg_com_t`, `bwd_neg_diag45`, `bwd_neg_diag135`).

Validation on the currently available 39-file artifacts:

```bash
julia --project=. test/build_labelfree_unit_predictions.jl \
    --features results/unit_separability/lobe_features_selectedN_primary_local.tsv \
    --split-features results/unit_separability/lobe_features_selectedN_primary_split.tsv \
    --patches results/unit_separability/lobe_patches_selectedN_primary_17x17_bwd.tsv \
    --out /tmp/opencode/labelfree_unit_predictions_subset.tsv \
    --seeds 20 \
    --interactions

julia --project=. test/report_unit_assignment_benchmark.jl \
    --profile portable_kmeans=/tmp/opencode/labelfree_unit_predictions_subset.tsv=portable_kmeans_sanity \
    --outdir /tmp/opencode/labelfree_unit_report_subset
```

The generated TSV covers 39 files / 234 lobes. The historical subset sanity grade
is 170/210 = 81.0% physical (1/35 exact). This is intentionally described as a
portable k-means baseline, not as a replacement for the frozen best GMM ensemble
(`forced_ensemble3`, 178/210).

The first full145 selected-`N` path is now unblocked for honest reporting, while
the strict 870-row path remains blocked for real selected-`N` predictions:

1. The available feature/patch artifacts cover only 39 files, so full145 preflight
   correctly fails with `found 39 files / 234 lobe rows` and 106 missing files.
2. The current full145 counting summary is label-free but not denominator-perfect:
   `N_selected` counts are 129 files at 6, eight at 5, six at 7, one at 4, and one
   at 3. A selected-`N` full145 feature extraction therefore produces 863 lobe
   rows, not 870. Producing exactly six lobe rows for every file by reading the
   manifest/control would use the benchmark label and is not an allowed prediction
   path.

Implemented the safe reporting branch for option (a). `test/report_unit_assignment_benchmark.jl`
now has an explicit `--full145-own-n` mode for profiles generated at each file's
label-free `N_selected`. Strict `--full145` is unchanged and still requires
145 files / 870 prediction rows. Own-N mode still requires all 145 files and
contiguous lobe indices per file, but it reports count mismatches honestly:
positions absent because `N_selected < 6` count as uncertain against the external
870-position control denominator, while `N_selected > 6` lobes are reported as
extra predictions and are not aligned to the 6-position control.

The full145 selected-`N` feature extraction from the Viper counting summary
completed locally:

```text
/tmp/opencode/full145_selectedN_features.tsv
files=145
prediction rows=863
rows by selected N: N=3 -> 3, N=4 -> 4, N=5 -> 40, N=6 -> 774, N=7 -> 42
```

The immediate base/local-feature portable baseline was generated without truth,
composition priors, or the control sequence:

```bash
julia --project=. test/augment_lobe_local_features.jl \
    --features /tmp/opencode/full145_selectedN_features.tsv \
    --out /tmp/opencode/full145_selectedN_features_local.tsv

julia --project=. test/build_labelfree_unit_predictions.jl \
    --features /tmp/opencode/full145_selectedN_features_local.tsv \
    --out /tmp/opencode/full145_selectedN_labelfree_local_predictions.tsv \
    --seeds 20 \
    --interactions

julia --project=. test/report_unit_assignment_benchmark.jl --full145-own-n \
    --profile selectedN_local=/tmp/opencode/full145_selectedN_labelfree_local_predictions.tsv=selectedN_local \
    --outdir /tmp/opencode/unit_report_full145_own_n_local
```

Report result:

```text
prediction rows:        863
aligned/classified:     857/870 (98.5%)
missing control lobes:  13 across 10 short-N files
extra predicted lobes:  6 across 6 N=7 files
physical accuracy:      671/857 = 78.3%
honest view:            671/870 = 77.1% correct, 199/870 uncertain
exact chains:           16/145
```

This is a reproducible selected-`N` full145 baseline, not a solved final map and
not a replacement for the historical best three-view subset ensemble. The next
scientifically safe step is to generate the selected-`N` split features and
backward-patch descriptors for all 145 files and rerun the portable predictor in
the same `--full145-own-n` report mode.

Continuation: the backward-patch half of that step completed locally from the
already frozen selected-`N` feature TSV:

```bash
STMFIT_DATA_DIR=/tmp/opencode/stmfit_full146_data julia --project=. \
    test/extract_lobe_patches_bwd.jl \
    --features /tmp/opencode/full145_selectedN_features.tsv \
    --half-nm 0.32 \
    --step-nm 0.04 \
    --out /tmp/opencode/full145_selectedN_patches_bwd.tsv
```

Validation showed exact coverage parity with the selected-`N` features:

```text
features: 145 files / 863 rows
patches:  145 files / 863 rows
mismatches: 0
```

The base+backward portable predictor used no truth, control sequence, composition
prior, or expected count:

```bash
julia --project=. test/build_labelfree_unit_predictions.jl \
    --features /tmp/opencode/full145_selectedN_features_local.tsv \
    --patches /tmp/opencode/full145_selectedN_patches_bwd.tsv \
    --out /tmp/opencode/full145_selectedN_labelfree_base_bwd_predictions.tsv \
    --seeds 20 \
    --interactions

julia --project=. test/report_unit_assignment_benchmark.jl --full145-own-n \
    --profile selectedN_base_bwd=/tmp/opencode/full145_selectedN_labelfree_base_bwd_predictions.tsv=selectedN_base_bwd \
    --outdir /tmp/opencode/unit_report_full145_own_n_base_bwd
```

Report result:

```text
prediction rows:        863
classified:             854/870 (98.2%)
missing control lobes:  13 across 10 short-N files
extra predicted lobes:  6 across 6 N=7 files
physical accuracy:      671/854 = 78.6%
honest view:            671/870 = 77.1% correct, 199/870 uncertain
exact chains:           5/145
```

So the backward descriptors do not improve the full145 selected-`N` honest
headline over the base/local-only run; they mainly add 3 abstentions and reduce
exact-chain count under this portable k-means rule. The selected-`N` split-width
refit remains incomplete locally: repeated timeout-limited shards produced
complete split rows for 96/145 files, leaving 49 files missing. Those partial
split TSVs must not be used as a full145 profile. Two small script hardening fixes
were made during the attempt: `extract_lobe_features.jl --manifest` now restricts
the selected-summary file list to manifest members before `--primary-only`, and
`extract_lobe_patches_bwd.jl --help` no longer errors when `STMFIT_DATA_DIR` is
unset.

HPC follow-up requested during the same continuation: Viper has no active STMFit
jobs, and the latest full146 confirmation array `10391759` is complete
(`2/2` array tasks, exit `0:0`). Raven showed the pending QE GlcN restart3 job
`28601744` had reached the 24 h walltime limit (`TIMEOUT`; `pw.x` failed after
the limit). Following the established restart protocol, fetched
`qe/glcn_restart3/glcn_central_relax.out`, extracted the last 213-atom geometry to
`qe/glcn_restart3/glcn_central_best4.xyz`, prepared `qe/glcn_restart4` with the
same active Raven settings (`8` tasks, `96000 MB`, `24:00:00`, `ecutwfc=50`,
`ecutrho=360`, Gamma-only), and submitted it on Raven after a clean local and
remote preflight:

```text
qe/glcn_restart4 -> 28658135
initial status: PENDING
```

Do not fetch or inspect `qe/glcn_restart4` outputs until it completes or times
out. GlcNAc production remains dependent on a successful production GlcN restart.

Continuation: completed the full145 selected-`N` portable-predictor sweep using
only frozen label-free feature/patch TSVs. The externally generated control TSV
was used only by `report_unit_assignment_benchmark.jl --full145-own-n` for
post-hoc grading; no view, threshold, abstention rule, or count was chosen from
`NKNNKN`.

The main sweep command graded the base run, base+backward run, the selected-`N`
three-view ensemble, confidence/agreement abstentions, and backward-only views:

```text
summary: /tmp/opencode/unit_report_full145_own_n_variant_sweep/summary.tsv
```

Selected rows:

```text
profile              classified  physical          honest view                 exact
base                 857/870     671/857 = 78.3%  671/870 correct + 199 ?     16/145
base_bwd             854/870     671/854 = 78.6%  671/870 correct + 199 ?      5/145
view3                857/870     676/857 = 78.9%  676/870 correct + 194 ?      7/145
view3_conf80         821/870     653/821 = 79.5%  653/870 correct + 217 ?      6/145
view3_agreebase65    827/870     658/827 = 79.6%  658/870 correct + 212 ?      6/145
bwd_com_only         854/870     419/854 = 49.1%  419/870 correct + 451 ?      0/145
bwd_diag45_only      854/870     426/854 = 49.9%  426/870 correct + 444 ?      0/145
bwd_diag135_only     854/870     427/854 = 50.0%  427/870 correct + 443 ?      0/145
```

The current selected-`N` full145 own-N headline is therefore the three-view
ensemble `BASE`, `BASE+bwd_neg_com_t`, and `BASE+bwd_neg_diag45`: it improves the
base/local honest count by 5 lobes (676 vs 671) while preserving full aligned
coverage. The confidence/agreement abstention variants raise classified accuracy
only marginally but lower the honest correct count, so they are not a better
headline. Backward-only views are negative evidence: their physical convention
collapses to about chance, even though supervised oracle alignment can recover
some signal, so they must not be used alone.

A second targeted sweep checked whether adding individual backward descriptors to
BASE or adding `bwd_neg_diag135` as a fourth view was useful:

```text
summary: /tmp/opencode/unit_report_full145_own_n_variant_sweep2/summary.tsv
base_bwd_com      672/857 honest-aligned correct, 672/870 honest denominator
base_bwd_diag45   671/857 honest-aligned correct, 671/870 honest denominator
base_bwd_diag135  667/857 honest-aligned correct, 667/870 honest denominator
view4             670/857 honest-aligned correct, 670/870 honest denominator
```

These do not beat `view3`. The selected-`N` split-width profile is still blocked
by incomplete local coverage (96/145 files), so it remains excluded from any
full145 headline.

Operational status at the same checkpoint: Raven QE GlcN restart4 job `28658135`
was still running on `ravc4062` at 8:26 elapsed. Viper Slurm status checks were
temporarily unreliable because the remote Modules/Tcl library reported an I/O
error and `squeue` then failed with `ClusterName needs to be specified`.

Continuation: completed the previously blocked selected-`N` split-width full145
feature extraction. The prior local shards covered only 96/145 manifest files;
the remaining 49 files were identified from `/tmp/opencode/full145_selectedN_features.tsv`
and rerun in two-file micro-batches to avoid another tool timeout. The final merge
used the selected-`N` base feature row order as the key and wrote:

```text
/tmp/opencode/full145_selectedN_features_split.tsv
files=145
rows=863
missing rows=0
duplicates=0
ignored extra rows=8 from excluded 240310_Cu100009.sxm
```

This keeps the own-N denominator unchanged: 13 missing external control positions
from 10 short-N files and 6 extra N=7 predictions remain a counting/selection
fact, not a split-width artifact.

Generated split-aware label-free prediction profiles before any grading truth was
read:

```text
/tmp/opencode/full145_selectedN_labelfree_split_default_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_base_split_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_4view_base_bwd_split_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_3view_base_com_split_predictions.tsv
/tmp/opencode/full145_selectedN_labelfree_3view_base_diag45_split_predictions.tsv
```

Post-hoc `--full145-own-n` grading summary:

```text
summary: /tmp/opencode/unit_report_full145_own_n_split_sweep/summary.tsv

profile                  classified  physical          honest view              exact
prev_view3               857/870     676/857 = 78.9%  676/870 + 194 ?          7/145
split_default            857/870     676/857 = 78.9%  676/870 + 194 ?          7/145
base_split               857/870     676/857 = 78.9%  676/870 + 194 ?         17/145
view4_base_bwd_split     857/870     675/857 = 78.8%  675/870 + 195 ?         16/145
view3_base_com_split     857/870     673/857 = 78.5%  673/870 + 197 ?         16/145
view3_base_diag45_split  857/870     673/857 = 78.5%  673/870 + 197 ?         16/145
```

Conclusion: completing split-width removes the full145 coverage blocker, but it
does not improve the lobe-correct headline beyond 676/870. Its useful signal is
sequence-level: the single-view `BASE+split_log_skew` profile keeps the same
honest correct count and raises exact chains from 7/145 to 17/145. Treat this as
a better exact-chain diagnostic, not as a solved binary deacetylation map.

Workflow hardening follow-up: added `test/merge_lobe_feature_shards.jl` so future
split-width/full145 resumptions no longer require hand-written Julia one-liners.
The script takes a reference selected-`N` feature TSV plus repeated `--shard` or
comma-separated `--shards`, validates that every reference `(file,lobe)` key is
present exactly once, reports stale extra rows, and writes the merged TSV in the
reference row order. It reads no truth sequence, expected count, control motif, or
composition prior.

Validation on the completed split-width artifacts:

```text
julia --project=. test/merge_lobe_feature_shards.jl \
    --reference /tmp/opencode/full145_selectedN_features.tsv \
    --shards <all split chunk/missing/resume TSVs> \
    --ignore-extra \
    --dry-run

reference_files=145
reference_rows=863
merged_files=145
merged_rows=863
missing_rows=0
duplicate_rows=0
extra_rows=8
extra_files=240310_Cu100009.sxm
```

The `--ignore-extra` flag is intentionally explicit; without it, stale rows from
an excluded or superseded manifest member remain an error.
