# Research journal

Dated record of decisions, experiments (including failures) and their reasons.
Add an entry for every experiment, scientific decision, failed approach, bug
fix or parameter change. Keep **Open questions** current.

The complete journal up to the 2026-09-29 consolidation (about 24,000 lines,
with every job, hash and intermediate result) is preserved at the git tag
`archive/pre-cleanup-20260929`:

```bash
git show archive/pre-cleanup-20260929:docs/src/journal.md
git show archive/pre-cleanup-20260929:docs/src/journal_archive.md
```

## Current state (2026-09-29)

**Method** (promoted): per-scan GCV counting (`chitosan.toml`), then
repeated-scan molecule consensus with registered-ROI refits
(`molecule_consensus.toml`), then per-scan assignment learned on corroborated
scans (`unit_assignment_corroborated_training.toml`), then binomial
latent-class fusion. Entry point: `test/run_molecule_consensus_chitosan.jl`.

**Benchmark** (fresh raw run 12025539, `3c414a4`; full145 own-N): 137/145
exact N (145/145 within one), 772/870 correct, 88/145 exact chains, 865/870
classified, 93 errors, against the previous record 129/694/43/855/161. Four
runs give 772, 772, 772 and 769 correct (88 exact chains and 137/145 each);
the spread comes from the time-limited global search.

**Limits:** NKNNKN edge-adjacency confound; exact chains in molecule-sized
blocks; counting parameters with benchmark-informed provenance; molds are
empirical templates; 10–20mer outputs not validated.

## Open questions

1. **Chemical validation.** Does the method detect acetyl groups, or the
   edge-adjacent position that NKNNKN happens to share with GlcNAc? A control
   molecule with a non-palindromic or homo-oligomer sequence is needed.
   Position-neutral features recover edge versus interior, not chemistry.
2. **Long-chain counting.** The axial support is truncated on long chains
   (`251206_013`: 1.5 of 7.5 nm, 4 lobes fresh, 12 manual). The support
   baseline is the 10% quantile of the axial profile, which assumes at least
   10% background bins; on long chains it lands on the molecule. A curved
   centreline does not fix this (2026-09-30). A length-invariant background
   estimate is the candidate fix.
3. **Long-chain assignment consistency.** On the 10–20mers the fusion fit gives
   θ0 = 0.45: per-scan calls disagree across repeat scans of the same molecule.
   What makes long-chain per-scan calls stable?
4. **Independent counting calibration.** The hybrid rule, support padding,
   fit width and κ threshold were chosen with benchmark grades. A label-free
   justification, or validation on an independent labelled set, is still missing.
5. **Physical molds.** The molds keep a legacy cube addressing and read a
   signed LDOS quantity. The fixed-density precision comparison never
   completed. A qualified constant-current observable is open.
6. **Calibrated uncertainty.** Fused posteriors are model posteriors under the
   latent-class model. Unfused confidences are uncalibrated vote margins.
7. **Deterministic fits.** The NLopt global search is time-limited
   (`global_maxtime`), so results depend slightly on node speed (772 or 769
   correct across four runs). An iteration-bounded budget would make runs
   reproducible.

## Lessons learned

**Counting**

1. The circular model is the anchor: σ∥ = σ⟂ always converges. Elliptical fits
   are LsqFit refinements from it. A global elliptical search (NLopt, 33
   parameters) escapes the isotropic saddle and diverges.
2. `min(ell, circ)` per N is more robust than penalty terms and adds no parameter.
3. GCV, not iid BIC/AICc, selects: residuals are spatially correlated over more
   than the fit tube. The 1D fit over-counts and stays diagnostic.
4. Saved results must be reproduced from raw before promotion. The saved hybrid
   counts (129/145) and the record 694/43/855 were not reproduced by fresh
   runs (123/145, 679/29/848).

**Assignment**

5. On a fixed per-scan feature set, dozens of single-change classifier
   variants were neutral or negative. Their changes concentrate on zero-margin
   vote ties and on marginal GlcNAc lobes. With true labels, the same features
   reach at most ~43/137 exact chains per scan (supervised
   leave-one-molecule-out ceiling).
6. Better image agreement is not better recognition. Registration,
   geometric profiling, local frames and local Gaussians all improved fits or
   correlations and lowered the grade.
7. Repeated scans of one molecule are the strongest label-free information.
   Symmetric fusion (mean or majority) locks in each molecule's errors, because
   per-scan errors are asymmetric: GlcNAc is missed, rarely invented. The
   binomial latent-class model captures this asymmetry and keeps a type
   detected in a minority of scans.
8. Whole-cohort learning makes exact chains sensitive to a few scans' counts;
   learning only on corroborated scans reduces it.

**Process**

9. One declared change, one fresh run, one external grade, no retuning after
   the grade. Refactors are verified byte-for-byte on stored outputs.

## History (condensed)

### Counting (May–June 2026)

The 2D elliptical chain model gave N ≠ 6 on 7/27 files of the 240817 set. The
fix was circular sweep → circ→ell LsqFit refinement → GCV. The robust-AICc guard
with an up-when-ambiguous branch then gave 39/39 exact on the primary 240817
set (June 17; `240817_043` recovered by the up branch). The support-midpoint
hybrid became the default after raising the saved full145 grade from 106/145
to 129/145 exact (143/145 within one). That saved summary is not reproduced by
fresh runs (123/145). Diagnostic selectors (blocked CV, support-marginalized
GCV, slope-heuristic MDL, stability, Laplace evidence, fwd/bwd consensus,
local-lobe evidence) did not beat the default on the benchmark or on synthetic
known-N data.

### Unit assignment exploration (June–August 2026)

Per-lobe features (amplitude, widths, split-profile skew, local prominence),
aligned residual patches, k-means and GMM views, self-training, connected
molds (geometric proxies, then DFT LDOS maps), and an empirical Fisher margin
were explored. A generic left/right width asymmetry turned out to reflect
overlap and envelope, not identity, which motivated connected chemical molds.
The best label-free profile was the soft vote of k-means 4-view and GMM 1-view
with constant-current mold and Fisher margins: **677/870 correct, 36/145 exact
chains, 854/870 classified** (August 2, "historical reference"). One of its
original feature producers was lost, so it cannot be rebuilt from raw.

### Native Julia reconstruction (September 16–19)

`test/run_reconstructed_chitosan.jl` rebuilt the pipeline natively (Julia 1.13,
no Python) under a separately named model. Mold and Fisher ports were verified
numerically against the Python reference. Full146 job 11786116 gave 666/849
classified, 24/145 exact chains: below the reference, not promoted.

### Bounded single-change comparisons (September 20–25)

Each row is one declared change, one Viper run and one external full145 grade
against the matched control. Correct /870, exact chains /145 (coverage when it
changed).

| Date | Change | Result | Control | Decision |
|---|---|---|---|---|
| 09-20 | Symmetric fwd/bwd fused-fit correction | 671 / 28 | 666 / 24 | Adopted |
| 09-20 | Physical-u Fisher mirror | 669 / 27 | 671 / 28 | Rejected |
| 09-20 | Matched residual filtering `S(data − model)` | 672 / 29 | 671 / 28 | Adopted |
| 09-20 | Affine half-plane descriptor (first-moment alternative 672/28) | 673 / 33 | 672 / 29 | Adopted |
| 09-21 | Signed CC margins; affine Fisher | 668 / 32; 672 / 32 | 673 / 33 | Rejected |
| 09-21 | Training-mean Fisher centring; final covariance shrinkage | 675 / 33; 665 / 34 | 673 / 33 | Centring adopted |
| 09-21 | Complete-disk patch support (cov. 849 → 852); final Gaussian score | 676 / 34; 666 / 32 | 675 / 33 | Support adopted |
| 09-21 | Complete-patch training (cov. 850) | 674 / 33 | 676 / 34 | Rejected |
| 09-21 | Robust median/IQR GMM normalization | 667 / 10 | 676 / 34 | Rejected |
| 09-21 | Equal-scan GMM weighting | 675 / 36 | 676 / 34 | Trade-off, not adopted |
| 09-21 | Continuous GMM seed vote | 671 / 26 | 676 / 34 | Rejected |
| 09-21 | Whole-scan GMM bagging | 665 / 10 | 676 / 34 | Rejected |
| 09-22 | Tied covariance throughout learning | 666 / 6 | 676 / 34 | Rejected |
| 09-22 | Whole-scan Fisher cross-fitting; within-scan GMM naming | 668 / 34; 676 / 34 | 676 / 34 | No gain |
| 09-22 | Rank-4 factor analyzer; Student-t (df 5) learning | 620 / 23; 631 / 25 | 676 / 34 | Rejected |
| 09-22 | All-update shrinkage; coherent Student density | 637 / 24; 627 / 25 | 676 / 34 | Rejected |
| 09-22 | Fresh raw-GCV counts after symmetric filtering (exact N 106 → 101) | 661 / 32 | 675 / 34 | Rejected |
| 09-22 | Frozen-geometry amplitude/background profiling | 675 / 32 | 676 / 34 | Rejected |
| 09-22 | Geometric variable projection | 632 / 20 | 676 / 34 | Rejected |
| 09-22 | Fixed-geometry acquisition registration (cov. 850) | 671 / 26 | 676 / 34 | Rejected |
| 09-22/23 | Registered native refits; multi-start/subpixel pilot | — | — | Inconclusive (span failures, shift at bound) |
| 09-23 | Local chain-tangent patch frames | 670 / 30 | 676 / 34 | Rejected |
| 09-23 | Locally oriented Gaussian model | 664 / 30 | 674 / 31 | Rejected |
| 09-23 | Cross-view residual subtraction (same / opposite view) | 675 / 30; 673 / 30 | 676 / 34 | Rejected |
| 09-23 | Tangent projection of physical-mold scores | 679 / 33 | 676 / 34 | Per-unit gain only; not adopted |
| 09-23 | Orientation derivative; shared mold state; adjacent-amplitude tangent | 678 / 30; 675 / 32; 671 / 29 | 679 / 33 | Rejected |
| 09-24 | Leave-target-out parity/mirror decoding | 670 / 32 | 679 / 33 | Rejected |
| 09-24 | QE cube-order (last-index-fast) molds | 676 / 27 | 679 / 33 | Not adopted |
| 09-24 | Observed-pixel preprocessing (cov. 848) | 680 / 35 | 678 / 28 | Not a champion |
| 09-24 | Frozen-learning reverse control; whole-scan exclusion; complete-observation training | 680 / 36; 672 / 19; 670 / 16 | 678 / 28 | Diagnostic / rejected |
| 09-25 | Width-conditioned counting (exact N) | 42 / 145 | 101 / 145 | Rejected |
| 09-25 | Exact background conditioning | same RSS, no stationary endpoint | — | Not retained |

**September 24 audits.** Saved hybrid counts combined with fresh assignment
give the record 694/43/855. Two fresh raw runs give 123/145 exact N and
679/29/848, identical to each other, so the record is not reproducible end to
end. The inherited calibration used benchmark grades (provenance limit). The
measurement audit found that the old bootstrap relied on fallbacks for 245/292
width views and 261/292 spacing views.

### DFT observable qualification (September 24–29)

Γ reconstruction, PAW attribution, direct wavefunctions, cutoff and basis
sensitivity, vacuum potentials, surface families, raw-image transfer, and
Davidson/CG precision comparisons. The accepted cubes and the 5e-5 Ry criterion
stand, but no calibrated constant-current observable was established and no
mold was replaced. Both CG cases timed out (12018453, 12023139). Summary in the
[DFT note](dft_calculation_note.md#Qualification-of-the-STM-observable).

### 2026-09-29 — Molecule consensus, latent-class fusion, promotion

**Diagnosis** (labels used only after the fact). The 145 files are repeated
scans of ~36 molecule sites; the largest series has 43 scans. On fresh N = 6
chains, GlcNAc positions 2 and 5 are 55% correct against 90–100% elsewhere.
Usually one of the two GlcNAc is missed, consistently across scans of the
same molecule. 24 of the record's 43 exact chains contain a zero-confidence
0.5 tie. The edge-adjacency confound was quantified: position-neutral
features recover edge versus interior, not chemistry.

**Consensus counting.** Absolute-frame registration (NCC ≥ 0.5; empty gap
0.39–0.55 in the pair distribution), strict majority over ≥ 3 scans, fixed-N
refit check with a registered-ROI retry. Fresh run (12019264): 137/145 exact,
145/145 within one, against 123/143 per scan. Assignment at these counts:
695/30/863, so exact chains did not improve.

**Corroborated training** (`9edcce4`, job 12023870): GMM and k-means learn only
on the 92 scans whose count agrees with the consensus. Result 702/33/863,
161 errors. Seed blocks (0–900) give identical exact chains: the gap follows
the cohort, not seed noise.

**Supervised ceiling** (diagnostic, leave-one-molecule-out with true labels):
at most ~43/137 exact chains per scan from the current features; raw patches
and non-linear models add nothing. Symmetric fusion (majority, mean
probability) raises per-lobe accuracy but lowers exact chains (24).

**Latent-class fusion.** Per-scan errors are asymmetric: missed GlcNAc,
rare false calls. A binomial latent-class model per physical lobe (EM for π,
θ0, θ1) keeps a type detected in a minority of scans. Fresh verification
(12025539, `3c414a4`): π = 0.342, θ0 = 0.035, θ1 = 0.622; **772/88/865/93**.
The independent run 12023870 with the same fusion gives the identical result.
**Promoted.** No measure regresses against the record.

**Application** (job 12023144, 25 targets + 109 collected repeat candidates):
17 counts change. Two long chains fail (4 and 2 lobes): their per-scan
support is truncated, not their consensus. The straight-tube explanation
recorded then was refuted on 2026-09-30 (support baseline, see below).
Fusion θ0 = 0.45. Not chemically interpretable.

### 2026-09-29 — Consolidation of the repository

**Why.** About 30 research branches and ~880 tracked files (rejected variants,
diagnostic campaigns, process apparatus, 24k-line journal) obscured the
promoted method.

**What.**

- Merged `research/labelfree-chain-20260929` and `research/mold-loo-20260924`
  into `consolidate/20260929`. Tagged the pre-cleanup state
  `archive/pre-cleanup-20260929`.
- Removed obsolete experiments and apparatus (879 → ~170 tracked files):
  AuthorityStorage, `integration/`, `openspec/`, the joint-proxy workflow and its
  mold/map templates, QE diagnostic campaigns and launchers, and rejected
  configs, sbatch files and scripts.
- Reduced the kept scripts to the promoted path. Removed acquisition shifts,
  local patch frames, residual features, tangent mold scoring, GMM
  covariance/mixture/bagging variants, training groups and seed-aggregation
  variants, and five config keys that no code read. The GMM builder rejects
  non-promoted policy values.
- Tracked the frozen molds (`templates/chitosan_cc_molds_native_v1.tsv`), their
  provenance and the two DFT frames. The molds rebuild byte-identically from the
  accepted cubes. `hpc/run_molecule_consensus.sbatch` now defaults to the
  promoted configs and tracked molds.
- Added `test/grade_consensus_run.jl` (label-free verification, then external
  grade). It reproduces 772/88/865/93 and 137/145 on run 12025539.
- Rewrote README, AGENTS, the HPC README and all `docs/src` pages around the
  current method. This journal was condensed.

**Verification.** GMM, k-means and patch outputs are byte-identical to stored
runs (fusion, v2 and record caches). All kept test files, package suites and
the synthetic end-to-end pipeline pass. A fresh raw Viper run of the cleaned
code with the new sbatch defaults (job 12028164, source `7c1b424`, COMPLETED
in 1h29, MaxRSS 4.5 GB, 687 output hashes match) passes all label-free checks
of `grade_consensus_run.jl`. It reproduces the promoted run byte for byte: final
predictions, per-scan predictions, final counts and input hashes are identical
to job 12025539 (772/88/865/93, 137/145 exact N, fusion π = 0.342, θ0 = 0.035,
θ1 = 0.622). Later commits change only documentation, tests and five
unread config keys.

### 2026-09-30 — Package pruning

**Why.** The user approved removing experimental options still present in the
packages after the consolidation.

**What.** Removed from `GaussianFit2D`: the `local_sigma_cap` overlap
constraint (width-conditioned counting, rejected 09-25), the `local_tangent`
peak orientation (locally oriented Gaussians, rejected 09-23), the
observed-only fit preprocessing, the optimizer diagnostics hook, and an unused
Current-channel weighting. Removed from `STMMolecularFit/selectors.jl` and
`test/batch_full.jl`: the diagnostic selectors (spatial blocked CV,
support-marginalized GCV, slope-heuristic MDL, stability selection, local-lobe
evidence, Laplace evidence, fwd/bwd consensus), `--cv-folds` and
`--refined-advisory`. Kept policies: `gcv`, `gcv_with_robust_aicc_guard`,
`support_midpoint_hybrid`, `adaptive_support_rescue`.

**Verification.** All test suites pass. On three benchmark scans covering
support-midpoint down/up and the guard's up-when-ambiguous branch
(`240307_015`, `240818_019`, `240817_043`), every count decision is identical
to the unpruned code. **Counting is not bit-reproducible run to run.** The
NLopt global search is time-limited (`global_maxtime = 10 s` per N), so a
rerun of the unpruned code on the same scan changes BIC values in the fifth
significant digit (275.923 → 275.940). The pruning is within that noise.

### 2026-09-30 — Curved centreline: prototype refutes the hypothesis

**Question** (user request). Would a curved centreline, instead of the straight
PCA tube, make long-chain counts trustworthy? Planned test: implement it and
check non-regression on the 6-mers.

**Prototype, before any package change.** Ridge points are intensity centroids
of `u` in 0.175 nm bins of `t`. A polynomial `u = f(t)` is fitted, with the
degree (0–4) chosen by weighted GCV. The support is recomputed in a ±0.16 nm
tube around `f`, with the unchanged `_active_t_support` rule, and compared to
the straight tube. Label-free, no fitting of lobes.

**Result.**

- The axial profiles along the straight and curved tubes are nearly identical
  on `251206_013`, `260220_083` and `240307_015`.
- The curved support is not systematically better on the 25 unknown targets.
  It is longer on some (`260115_040` 6.1 → 7.4 nm), shorter on others
  (`260220_030` 5.3 → 3.9, `260129_026` 3.9 → 3.0), and unchanged on the
  truncated cases (`260220_083` 1.9 → 1.85). `251206_013` only reaches its
  full 7.5 nm through the `no_component_long_enough` fallback.
- The selected degree often hits the maximum. The ridge centroid follows the
  lobe-to-lobe modulation, not a smooth bend.

**Actual cause.** `_active_t_support` sets the baseline to the 10% quantile of
the axial profile bins and the noise to the MAD of the bins below it. This
assumes that at least 10% of the bins are background.

- On a 6-mer, the profile ramps slowly to zero at both ends (`240307_015`:
  baseline 0.03, threshold 0.054, peak 0.16), so the rule works.
- On a long chain that fills about 90% of the profile, the baseline lands on
  the molecule's flanks (`251206_013`: 0.083). The noise, taken from the steep
  end bins, is large, and the threshold reaches 0.149 against a peak of 0.166.
  Only the brightest segments pass.

The same mechanism explains `260215_022` (2.9 of 5.6 nm) and `260221_044`
(5.1 of 7.9 nm).

**Decision.** The curved centreline is not implemented; the evidence does not
support it. The earlier journal and docs statement that "the straight tube
loses curved chains" is withdrawn and corrected. Candidate fix, pending user
approval: a background level and noise estimated independently of chain length
(e.g. from pixels outside the molecule). It would be declared and frozen before
testing, checked for non-regression on the 6-mer benchmark, and compared
descriptively to the manual reads of the long chains. This changes the frozen
counting rule, so it needs a declared experiment.

### 2026-09-30 — Viper verification of the pruned code: run-to-run variability

Job **12030359** (source `ccbd743`, COMPLETED in 1h44, 687 output hashes
match) passes all label-free checks of `grade_consensus_run.jl`. Result:
**137/145 exact N (145 within one), 769/870 correct, 88/145 exact chains,
865/870 classified, 96 errors**. Per-scan calls give 700/33/863/163. Fusion:
π = 0.330, θ0 = 0.035, θ1 = 0.607 (101 calls changed).

Against the promoted runs (772/88/865/93, identical in 12025539, 12023870 and
12028164):

- all 146 per-scan counts and all final counts are identical;
- the fixed-N refits differ on 20 of 146 scans, mostly by less than
  0.003 nm or 0.2% in amplitude, but by up to 0.09 on `240817_050` and
  `241113_089`;
- this changes 5 per-scan calls, and learning on the whole cohort and fusion
  propagate them to −3 correct.

The removed code was unreachable on the promoted path. The difference is the
known non-determinism of the time-limited global search
(`global_maxtime = 10 s` per N). This node was slower (1h44 against 1h29), and a
local rerun of the unpruned code on one scan already changed BIC values in the
fifth digit. **The promoted benchmark numbers are one realization: four runs
give 772, 772, 772 and 769 correct, all with 88 exact chains and 137/145
exact N.** Making the global search iteration-bounded instead of time-bounded
would remove this dependence on node speed. It changes the fitting budget, so
it would be a declared experiment (open question 7).

### 2026-09-30 — Declared experiment: deterministic fits and length-invariant support

**Deliverable** (user request: "fix the support threshold, make the fits
deterministic"). Reproducible fits, and axial supports that are not truncated
on long chains, without regressing the 6-mer benchmark. Declared before any
grade; no label or manual read was used to define the rules.

**Change 1, determinism.** `global_maxtime = 0` (NLopt time limit disabled) in
both counting configs; the global search stays bounded by
`global_maxiter = 10000`. Local evidence: with the time limit, the same scan
gives different BIC values under CPU contention (1 vs 3 parallel fits). Without
it, 1 and 3 threads give identical outputs, with no runtime cost on an
uncontended scan (3:45 vs 3:48).

**Change 2, support rule `half_maximum_cap`.** The legacy threshold
(10% profile quantile + 2.5 × MAD noise) is capped at
`background + 0.5 × (peak − background)`. `background` is the median of the
preprocessed image outside the dilated molecule ROI, and `peak` is the maximum
of the axial profile. Run selection and padding are unchanged.

A pure half-maximum rule was considered first and rejected before any grade,
from label-free profiles only. On `260221_048`, whose profile dips to 0.08
where the chain bends out of the tube, it cut the chain into two segments
(7.2 → 4.6 nm). The cap only ever lowers the threshold, so it keeps the
legacy behaviour wherever that behaviour was not pathological.

Label-free support check, no fitting:

- 6-mers: 52/146 supports get longer (median +0, q95 +0.30 nm). The largest
  gains are `240818_019` (1.8 → 3.7 nm) and `240814_020` (2.4 → 3.5 nm).
- Long chains: 11/25 change, and the truncations disappear: `251206_013`
  1.5 → 7.2, `260220_083` 1.9 → 5.8, `260215_022` 2.9 → 5.3 nm.

**Evaluation, one grade each, no retuning afterwards.**

- Run A (control): change 1 only, legacy support.
- Runs B and B′: changes 1 + 2, two identical runs. B = B′ byte for byte is
  the determinism test.
- Each is the full promoted method from raw on Viper, graded with
  `grade_consensus_run.jl`.
- Acceptance for the benchmark: B does not regress against A on counts (exact
  N, within one), correct, exact chains, classified or errors.
- The 25 long chains are then run with the 10–20mer config and compared
  descriptively with the manual reads. They are not a criterion, and nothing
  is tuned on them.
