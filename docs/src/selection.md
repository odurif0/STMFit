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
         └─ Chitosan default: robust-AICc guard reports N_selected
         └─ Raw baseline override: --selection-policy gcv reports N_eff
```

## Why GCV as the per-candidate score?

GCV is the default per-candidate sweep score because it is analytical, fast, and does not require
refitting folds. BIC is still computed and remains useful for diagnostics and
legacy comparisons, but it is not the default batch selection criterion.

- `n_eff` (effective sample size) is estimated as `length(zfit) ÷ 9`,
  which may not perfectly capture spatial correlation in STM images.
- BIC assumes all parameters contribute equally, but extra sigma parameters
  can absorb noise without improving predictive accuracy.
- On ambiguous STM images, BIC can marginally prefer over-fit models.

The `selection_criterion` field can be set to `"gcv"`, `"bic"`, `"aicc"`, or
`"cv"`. The chitosan calibration uses `"gcv"` with `cv_method="gcv"`.

## Robust-AICc overfit guard

`test/batch_full.jl` supports an integrated robust primary selection policy.
For `config/chitosan.toml`, this is now the configured default:

```bash
julia -t 4 --project=. test/batch_full.jl 48 --config config/chitosan.toml
```

The raw GCV/effective baseline remains available with `--selection-policy gcv`.

This policy is label-free: it does not use an expected `N`, does not prefer
`N=6`, and does not use benchmark labels during fitting or selection. It is a
conservative overfit guard layered on top of the standard GCV/effective result:

1. Compute the standard circ→ell `N_eff` exactly as in the default pipeline.
2. Fit an auxiliary exhaustive elliptical candidate set.
3. Rescore those candidates with Student-t robust AICc (`nu=8` by default;
   override with `--robust-guard-nu`).
4. Let `robust_aicc_N` be the lowest robust-AICc candidate.
5. Downshift only if robust AICc prefers a simpler model:

```text
N_selected = robust_aicc_N  if robust_aicc_N < N_eff
             N_eff          otherwise
```

Because the rule is down-only, it can veto extra lobes that look like
over-segmentation, but it cannot add lobes or encode a target count.

### Output columns

When the policy is enabled, batch summaries include:

- `N_selected`: the primary result under `--selection-policy`.
- `selection_policy`: usually `gcv` or `gcv_with_robust_aicc_guard`.
- `selection_source`: `N_eff` source when kept, or `robust_aicc_guard` when
  downshifted.
- `N_refined`: same guarded count for compatibility with earlier audit output.
- `refined_policy`: currently `overfit_guard_down_only` when a robust advisory
  was available.
- `refined_source`: `N_eff` when kept, or `ell_robust_aicc` when the integrated
  guard downshifted.
- `robust_aicc_N`: the auxiliary robust-AICc-selected count.

### Current validation status

On the 240817 chitosan clean benchmark, the integrated guard improves exact
agreement from `N_eff = 35/39` to `N_selected = 38/39`, with all four clean
target cases corrected and only `240817_036.sxm` missed.  That file is visually
ambiguous and should not be used to tune another special rule.

This result is validation evidence, not a fitting prior.  Synthetic known-N
validation also favored the guard over raw GCV in the tested `N≈4–8` regime.
Support-adaptive variants remain experimental and should be evaluated as
separate policies rather than folded into the default prematurely.

## Experimental adaptive support rescue

`--selection-policy adaptive_support_rescue` keeps the first-pass support and
GCV selector unless the selected count is at the support-feasibility ceiling.
That trigger is label-free: it compares the detected support length to the
minimum physically usable spacing implied by the calibration.  If triggered, a
second pass uses permissive support parameters and is accepted only when it
increases support length, selects a larger `N`, and keeps circular/elliptical
counts coherent.  Otherwise the standard result is retained.  After any rescue
decision, the validated robust-AICc guard is applied down-only on the active
support. The benchmark and 10–20mer adaptive configs use this same rule; the
10–20mer config differs only by extending `n_max` to allow long chains.
`adaptive_robust_guard_max_drop` exists as a non-default diagnostic option to
cap automatic downshifts, but it is not part of the common benchmark-aligned
workflow.

The reference config is `config/chitosan_adaptive_support_rescue.toml`.  It is
experimental and not the production default.

Operationally, the guard is safe-by-fallback: if the auxiliary robust-AICc sweep
fails, batch processing keeps the standard `N_eff` result and records the guard
failure instead of failing the file.

## Experimental spatial blocked CV selector

`test/batch_full.jl` also supports a more objectivable, but currently diagnostic,
spatial blocked cross-validation selector:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy spatial_blocked_cv \
  --cv-folds 3
```

The rule is pre-registerable and label-free:

1. Compute the standard circ→ell candidate set.
2. Sort fit pixels by their axial coordinate along the molecule.
3. Split them into contiguous spatial/axial blocks.
4. For each candidate `N`, refit on all blocks except one and score the held-out
   block with Student-t negative log-likelihood.
5. Select the `N` with the lowest mean held-out NLL.

This is mathematically attractive because it measures predictive performance on
spatial regions not used for fitting, rather than relying on an analytic GCV
approximation.  It is therefore a good candidate for future frozen validation.

Current status: **diagnostic only**.  Early smoke tests show that the raw blocked
CV selector can be unstable for STM chain counting: for example it downshifted
`240817_002.sxm` from `6` to `5`, selected `10` for `240817_017.sxm` with
3 folds, and kept `240817_043.sxm` at `5`.  This indicates that the idea is
more objectivable, but the current direct selector is not yet robust enough to
replace `N_eff` or the robust-AICc guard.

Use it to study predictive stability, not as the default selector.

## Experimental support-marginalized GCV selector

`test/batch_full.jl` supports a cheap support-sensitivity selector:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy support_marginalized_gcv
```

A conservative guarded variant is also available:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy support_marginalized_gcv_guard
```

This selector does **not** refit models.  It keeps the fitted circ→ell candidate
set and rescores each candidate across a frozen support-padding grid:

```text
support_padding_nm ∈ {0.00, 0.10, 0.20, 0.25, 0.35, 0.50}
```

For each support, it recomputes the fit window, evaluates each candidate's GCV
on that support, converts scores to relative regret against the best candidate
for that support, then selects the `N` with the lowest median regret.  The 75th
percentile regret and lower `N` are used as deterministic tie-breaks.

```text
regret_s(N) = (GCV_s(N) - min_M GCV_s(M)) / |min_M GCV_s(M)|
N_selected = argmin_N (median_s regret_s(N), q75_s regret_s(N), N)
```

The goal is not to prefer a target count, but to ask whether a selected `N` is
stable under plausible support boundaries.  This makes it useful for diagnosing
support-sensitive files such as `240817_043.sxm`.

Current status: **diagnostic/experimental only**.  Smoke tests show sensible
behaviour for some controls (`240817_002.sxm → 6`, `240817_026.sxm → 6`) and
some overfit cases (`240817_019.sxm → 6`), but it does not resolve all target
cases (`240817_017.sxm → 7`, `240817_043.sxm → 5`, `240817_058.sxm → 5` with
the current frozen grid/rescore-only rule).  It is therefore useful as a support
stability diagnostic, not as a default primary selector.

The guarded variant treats the support-marginalized result as bounded evidence
for over-segmentation.  It only downshifts when the support-marginalized median
regret is at least `0.02` better than the GCV `N_eff`, and caps the downshift to
one lobe:

```text
if N_support < N_eff and regret(N_support) + 0.02 < regret(N_eff):
    N_selected = max(N_support, N_eff - 1)
else:
    N_selected = N_eff
```

If there is no clear downshift but `N_eff-1` is within the existing GCV
ambiguity tolerance (`5%`) of `N_eff` under the support-marginalized median and
q75 regrets, the guard may also choose the one-step simpler model by parsimony.
This is a frozen one-standard-error-style rule and never increases `N`.

This avoids aggressive support-only jumps such as `7→5`; in smoke tests it kept
`240817_002.sxm` and `240817_026.sxm` at `6`, changed `240817_019.sxm` to `6`,
and changed `240817_058.sxm` from the raw support choice `5` to the guarded
choice `6`.  The parsimony check did not change `240817_017.sxm` because the
lower model was not stable enough under q75 regret.  It still leaves
`240817_043.sxm` at `5`, so it is a safer diagnostic guard, not a replacement
for the robust-AICc guard.

## Experimental slope-heuristic MDL selector

`test/batch_full.jl` also supports a file-adaptive MDL selector based on the
Birgé–Massart slope heuristic:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy slope_heuristic_mdl
```

For each valid circ/ell candidate, it computes a contrast from the residual sum
of squares and an effective sample size:

```text
C_m = n_eff * log(RSS_m / n_pixels)
D_m = number of free parameters
```

It then estimates the empirical overfit slope from the high-complexity half of
candidates using a robust Theil–Sen slope:

```text
C_m ≈ a - α D_m
```

and selects the model minimizing the doubled-slope MDL score:

```text
score_m = C_m + 2 α D_m
```

This is attractive because the complexity penalty is estimated from the file's
own contrast–dimension curve rather than fixed from a benchmark.  It uses no
expected `N`, no target count, and no benchmark labels.

Current status: **experimental/diagnostic only**.  Smoke tests show that this
selector is conservative on controls (`240817_002.sxm → 6`, `240817_026.sxm →
6`) but does not recover the current target overfit cases (`017 → 7`, `019 →
7`, `058 → 7`) and keeps the support-sensitive `043 → 5`.  It is therefore a
principled diagnostic criterion, not currently a better primary selector.

## Experimental stability-selection selector

`test/batch_full.jl` supports a support-perturbation stability selector:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy stability_selection
```

It uses the same frozen support-padding grid as support-marginalized GCV, but
asks a different question.  Instead of minimizing median regret, it counts how
often each `N` is competitive across support perturbations:

```text
N is competitive on support s if GCV_s(N) ≤ 1.01 * min_M GCV_s(M)
```

The selector ranks candidates by:

```text
1. largest number of competitive supports
2. largest number of feasible supports
3. lower N as parsimony tie-break
```

This is useful as a stability diagnostic because it reports how often a count
remains near-optimal under small support-boundary changes.  It uses no expected
`N`, no target count, and no benchmark labels.

Current status: **experimental/diagnostic only**.  Smoke tests show stable
controls (`240817_002.sxm → 6`, `240817_026.sxm → 6`) and one corrected overfit
case (`240817_019.sxm → 6`), but it leaves `240817_017.sxm → 7`, keeps
`240817_043.sxm → 5`, and is too conservative on `240817_058.sxm → 5`.  It is
therefore not a replacement for the robust-AICc guard.

## Experimental local-lobe-evidence guard

`test/batch_full.jl` also supports a local resolvability diagnostic:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy local_lobe_evidence
```

This policy is a down-only guard on top of `N_eff`.  It decodes the fitted chain,
examines adjacent lobes, and merges a pair as locally unresolved when both of the
following are true:

```text
center separation < 2.0 × mean σ∥
and valley depth is weak: valley SNR < 3.0 or valley fraction < 0.2
```

It then counts connected resolved components.  The highest `N ≤ N_eff` is
accepted only if at most one adjacent pair is unresolved:

```text
N_resolved ≥ N - 1
```

The rule uses no expected `N`, no target count, and no benchmark labels.  It also
never increases `N`; local geometric resolvability can veto redundant lobes but
cannot prove that a missing lobe should be added.

Current status: **diagnostic/inconclusive**.  The criterion is physically
interpretable for isolated peaks, but chitosan lobes are intentionally allowed to
overlap in a continuous molecular chain.  Smoke tests kept controls such as
`240817_026.sxm` at `6`, but the resolved-component count collapsed to one
connected chain (`resolved=1`, `unresolved_pairs=5`), so the guard treated the
local evidence as inconclusive and fell back to GCV.  This makes it useful as an
audit of peak separability, not as a better primary selector.

## Experimental Laplace-evidence selector

`test/batch_full.jl` supports an approximate Bayesian-evidence selector:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy laplace_evidence
```

and a safer one-lobe overfit guard:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy laplace_evidence_guard
```

For each valid circ/ell candidate, it computes a Gauss–Newton/Laplace local
evidence approximation from finite differences of the model prediction with
respect to normalized bounded parameters.  With Student-t residual weights and
singular values `λ` of the weighted Jacobian, the score is:

```text
fit    = 2 * (n_eff / n_pixels) * StudentT_NLL
occam  = Σ log(1 + λ)
d_eff  = Σ λ / (1 + λ)
sloppy = (n_params - d_eff) * log(n_eff)
score  = fit + occam + sloppy
```

The intent is to penalize locally sloppy or unidentified extra parameters using
the fitted model's own curvature, without expected `N`, target counts, or labels.

The direct selector chooses the minimum score.  The guarded variant treats this
as overfit evidence only: if Laplace prefers a lower count than `N_eff`, it
downshifts by at most one lobe; otherwise it keeps `N_eff`.

Current status: **experimental/diagnostic only**.  Smoke tests show that direct
Laplace evidence is too parsimonious for some chitosan files (`019 → 5`, `058 →
5`).  The guarded variant is safer: it keeps `026 → 6`, changes `017 → 6`,
`019 → 6`, and `058 → 6`, but cannot increase the support-sensitive `043` from
`5` to `6`.  It is a principled overfit guard candidate, not yet a replacement
for the integrated robust-AICc guard.

## Experimental fwd/bwd direction-consensus selector

`test/batch_full.jl` supports a forward/backward scan-consensus selector:

```bash
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy fwd_bwd_consensus
```

This selector exploits the physical fact that true molecular signal should be
present in both forward and backward STM scans.  For each candidate `N` fitted on
fused data, it evaluates the model on fwd-only and bwd-only preprocessed pixels,
does a linear recalibration per scan (`z ≈ a × model + b`) to absorb brightness
differences, and selects the `N` with the lowest joint GCV:

```text
joint_GCV(N) = (2·n_pixels) / (2·n_pixels - p_eff)² × (RSS_fwd + RSS_bwd)
p_eff = n_chain_params + 4   (a_fwd, b_fwd, a_bwd, b_bwd)
```

This is label-free, uses only fitted candidates (no refitting), and adds genuine
physical information from scan-direction replication.

Current status: **experimental/diagnostic only**.  Smoke tests show that the
simple joint-GCV with per-scan recalibration does not discriminate better than
fused GCV for the chitosan benchmark overfit cases (`017 → 7`, `019 → 7`,
`058 → 7`), though controls are preserved (`026 → 6`).  The support-sensitive
`043 → 5` is also unchanged.  The recalibration absorbs most scan-to-scan
differences, making the joint GCV essentially equivalent to fused GCV.  A more
discriminating approach (lobe-level fwd/bwd amplitude consistency, or a true
joint refit) would be needed for this information to improve selection.

## Synthetic known-N selector validation

`test/synthetic_known_n_validation.jl` provides a label-free validation harness
for comparing selectors on synthetic STM-like Gaussian chains with known true
counts:

```bash
julia --project=. test/synthetic_known_n_validation.jl \
  --cases 50 \
  --seed 1234 \
  --noise-scale 1.0 \
  --mode circ_ell \
  --out results/synthetic_known_n/summary_50_circ_ell.tsv

julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50_circ_ell.tsv \
  --out results/synthetic_known_n/aggregate_50_circ_ell.tsv

# Fast circular-only mode remains available for cheaper stress tests:
julia --project=. test/synthetic_known_n_validation.jl \
  --cases 50 \
  --seed 1234 \
  --noise-scale 1.0 \
  --mode circular \
  --out results/synthetic_known_n/summary_50.tsv

julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50.tsv \
  --out results/synthetic_known_n/aggregate_50.tsv

julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50.tsv \
  results/synthetic_known_n/summary_50_seed2026.tsv \
  results/synthetic_known_n/summary_50_seed4321.tsv \
  results/synthetic_known_n/summary_50_seed1234_noise05.tsv \
  results/synthetic_known_n/summary_50_seed1234_noise15.tsv \
  --out results/synthetic_known_n/aggregate_multiseed_noise.tsv

julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50_circ_ell.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed2026.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed4321.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed1234_noise05.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed1234_noise15.tsv \
  --out results/synthetic_known_n/aggregate_circ_ell_multiseed_noise.tsv
```

The script generates in-memory `SXMImage` cases whose true `N` cycles through
`4..8`, adds lobe-position/width/amplitude jitter, baseline tilt, independent
fwd/bwd noise, and occasional scan artifacts, then runs the requested synthetic
mode (`circular` or `circ_ell`).
It applies the core selector implementations from `STMMolecularFit` and writes a
TSV with one row per case and policy:

- `case_id`, `seed`, `true_N`, `artifact`, `noise_scale`
- `policy`, `N_eff`, `N_selected`, `abs_error`
- `status`, `score_or_source`

Use `--policies` to restrict the comparison, for example:

```bash
julia --project=. test/synthetic_known_n_validation.jl \
  --policies gcv,gcv_with_robust_aicc_guard,laplace_evidence_guard
```

Current status: **phase-2 validation available**.  The default `--mode circular`
is still a cheap stress test, but `--mode circ_ell` runs the closer batch analog:
fixed-label-free circular candidate sweep followed by local elliptical refinement
of valid circular candidates.  In both modes the candidate search window is fixed
and does not depend on `true_N`; synthetic labels are used only for external
grading (`abs_error` and aggregate exact-match summaries), never inside fitting
or selection.

The companion aggregation script summarizes exact rate, mean absolute error,
over-selection, under-selection, and error counts for all cases and stratified by
`true_N`, artifact class, seed, and noise scale.  It accepts multiple TSV inputs
and can read both legacy 10-column synthetic summaries and newer summaries with
the `noise_scale` column.

In the first 50-case `circ_ell` run (`seed=1234`, `noise_scale=1`), the effective
GCV baseline improved substantially relative to circular-only (`43/50` exact,
mean absolute error `0.20`).  The robust-AICc guard remained best on exact rate
(`46/50`, mean absolute error `0.14`, no over-selections and four
under-selections), while stability selection was close (`44/50`, mean absolute
error `0.18`).

Across the current 250-case phase-2 `circ_ell` grid (three seeds at
`noise_scale=1`, plus `noise_scale=0.5` and `1.5` for seed `1234`), robust-AICc
guard remains best overall: `231/250` exact, mean absolute error `0.128`, with
four over-selections and fifteen under-selections.  GCV is much stronger in
phase 2 than in circular-only mode (`213/250`, mean absolute error `0.20`) but
still trails the robust guard.

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
- `N_selected`: policy-level primary result. With `selection_policy="gcv"` it
  equals `N_eff`; with the chitosan default `gcv_with_robust_aicc_guard` it may
  be a lower robust-AICc guarded count.
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
