# Model selection

How STMFit chooses the number of lobes N of each chain: first per scan, then
across repeated scans of the same molecule. No step reads an expected N, a
benchmark label or a sequence.

## Per-scan hierarchy

```
Level 1  score each candidate N by GCV
         score_eff(N) = min(score_circ(N), score_ell_refined(N))
Level 2  N_circ, N_ell, N_eff = best valid circular / elliptical / effective model
Level 3  N_selected = policy applied to N_eff (chitosan: support_midpoint_hybrid)
Level 4  molecule consensus across repeated scans (test/run_molecule_consensus_chitosan.jl)
```

### Why GCV

GCV (`RSS·n/(n−p)²`) is analytical, needs no refit and needs no effective sample
size. STM residuals are strongly spatially correlated (lag-1 ρ ≈ 0.9–0.95,
correlation range 17–100 px, wider than the ~10 px fit tube), so iid BIC/AICc
values are not calibrated evidence. BIC/AICc are therefore diagnostics or
guards only. Their `n_eff = max(10, n ÷ 9)` is a fixed placeholder; do not
retune or reinterpret it. GCV is still a parameter-count approximation for a
nonlinear constrained fit: it is the canonical practical score, not a
calibrated count uncertainty.

### Why `min(ell, circ)`

The circular model (σ∥ = σ⟂) is nested in the elliptical model. If local
elliptical refinement ends worse than the circular solution at the same N, it
failed to improve that model. Taking the minimum guards against this without a
new parameter.

### Chitosan policy: `support_midpoint_hybrid`

The default in `config/chitosan.toml`. Two label-free layers act on `N_eff`:

**Robust-AICc guard.** An auxiliary exhaustive elliptical candidate set is
rescored with Student-t robust AICc (`robust_guard_nu = 8`); `robust_aicc_N` is
its minimum.

```text
if robust_aicc_N < N_eff:                                   # veto over-segmentation
    N_guarded = robust_aicc_N
elif robust_aicc_N == N_eff + 1 and ambiguous_eff           # recover one missing lobe
     and delta_GCV_rel_eff <= 0.05 and runner_up_N_eff == N_eff + 1:
    N_guarded = robust_aicc_N
else:
    N_guarded = N_eff
```

**Support midpoint.** The measured support length and the physical
spacing/overlap bounds give a feasible interval
`[support_N_min, support_N_max]`; `mid = round((min + max) / 2)`.

```text
if N_guarded > mid + 1:                          N_selected = mid
elif N_guarded > mid:                            N_selected = N_guarded - 1
elif N_guarded < mid and N_eff > N_guarded
     and delta_GCV_rel_eff <= 0.30:              N_selected = N_guarded + 1
else:                                            N_selected = N_guarded
```

The upward branch is deliberately stricter than the downward one. If the
auxiliary guard fit fails, the file keeps `N_eff` and records the failure.

**Provenance.** These rules and thresholds, `fit_width_nm` and
`support_padding_nm` were chosen historically while looking at known-count
benchmark grades. Inference reads no label, but this is not an independent
calibration ([Calibration](calibration.md)).

### Long chains: `adaptive_support_rescue`

Used by `config/chitosan_10_20mer_adaptive_support_rescue.toml` (`n_max = 24`).
If `N_eff` sits at the support-feasibility ceiling, a second pass with
permissive support settings is accepted only when support length and N both
increase and circular/elliptical counts stay coherent. The robust-AICc guard
then acts down-only on the active support. Known limit: the support is a
straight tube, so long curved chains leave it and are undercounted
([journal](journal.md)).

### Other policies

`test/batch_full.jl --selection-policy` also accepts `gcv` (raw `N_eff`) and
`gcv_with_robust_aicc_guard` (guard without the support-midpoint layer).
Diagnostic selectors (spatial blocked CV, support-marginalized GCV,
slope-heuristic MDL, stability selection, local-lobe evidence, Laplace
evidence, fwd/bwd consensus) did not beat the default on the expanded
benchmark or on synthetic known-N data and were removed on 2026-09-30. Their
code and evaluations are at the tag `archive/pre-cleanup-20260929`.

### Output columns (`summary_overlap060_hard.tsv`)

| Column | Meaning |
|---|---|
| `N_selected` | Count under the configured policy (primary result) |
| `N_eff`, `N_ell`, `N_circ` | Best effective / refined elliptical / circular model |
| `selection_policy`, `selection_source` | Policy, and which layer made the final move (`ell`, `circ`, `robust_aicc_guard`, `support_midpoint_down`, `support_midpoint_down_to_mid`, `support_midpoint_up`) |
| `N_refined`, `refined_policy`, `refined_source`, `robust_aicc_N` | Guard stage audit |
| `ambiguous_eff`, `runnerup_N_eff`, `delta_GCV_rel_eff` (and `_ell`) | Close second-best diagnostics (ΔGCV/GCV ≤ 5%); QC only |
| `support_2D_ell_nm`, `support_2D_circ_nm` | Active support length |
| `N_1D` | 1D diagnostic count (`NA` unless `--no-skip-1d`) |

## Repeated-scan molecule consensus

STM sessions often re-image one molecule many times, with drift and different
ranges. `test/build_molecule_consensus.jl` treats these scans as repeated
measurements:

1. **Ordering.** Scans are sorted by acquisition time (`REC_DATE`, `REC_TIME`).
   Only consecutive scans of the same date can be linked.
2. **Registration.** Each scan is resampled in the absolute piezo frame around
   the earlier scan's chain centroid (window ±3 nm, 0.04 nm grid), high-pass
   filtered (σ = 0.32 nm), and registered by a pure translation (≤ 1.6 nm).
   The pair is linked when NCC ≥ 0.5 and the centroids are ≤ 6 nm apart.
3. **Tracks.** Chains of linked consecutive scans form tracks.
4. **Consensus.** A track of ≥ 3 scans with a strict-majority per-scan count
   N* gives N* to a disagreeing scan, provided the reference footprint,
   translated by the accumulated drift, lies inside that scan's frame
   (margin 0.3 nm). Rules: `agrees`, `consensus_applied`, `track_too_short`,
   `no_strict_majority`, `footprint_outside_frame`.
5. **Refit check** (`run_molecule_consensus_chitosan.jl`). The unchanged
   fixed-N fitter must fit each changed scan at N*: first on its own support,
   then once on the registered footprint plus 1.2 nm (`consensus_registered_roi`).
   Otherwise the scan keeps its own count (`consensus_fit_failed`).

Singletons, ties and failed refits keep their per-scan count. Settings are in
`config/molecule_consensus.toml`; their justification is in
[Calibration](calibration.md). On the benchmark the final counts are 137/145
exact (per-scan counting: 123/145), 145/145 within one.
