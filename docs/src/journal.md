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

**Benchmark** (fresh raw run with deterministic fits, 12042058, `486ddec`;
full145 own-N): 137/145 exact N (145/145 within one), 771/870 correct, 88/145
exact chains, 865/870 classified, 94 errors, against the previous record
129/694/43/855/161. Earlier time-limited runs gave 772 (three runs) and 769.

**Positions:** lobe positions repeat to 0.07 nm along and 0.09–0.11 nm across
the chain over repeated scans (2026-10-04).

**Limits:** NKNNKN edge-adjacency confound. The benchmark GlcNAc signal is
the contrast with the dimmer terminal lobe; without it the image features
separate the classes only weakly (AUC ≤ 0.66, amplitude) or not at all
(shape, DFT molds). Also: exact chains in molecule-sized
blocks; counting parameters with benchmark-informed provenance; molds are
empirical templates; 10–20mer outputs not validated.

## Open questions

1. **Chemical validation.** Does the method detect acetyl groups, or the
   edge-adjacent position that NKNNKN happens to share with GlcNAc? A control
   molecule with a non-palindromic or homo-oligomer sequence is needed.
   Position-neutral features recover edge versus interior, not chemistry.
2. **Long-chain counting vs assignment geometry.** The capped support rule
   (`half_maximum_cap`, 2026-09-30) fixes long-chain supports. Mean
   |N − N_manual| drops from 1.33 to 0.92, and `251206_013` goes from 4 to 11
   lobes (12 read). It also improves 6-mer per-scan counts (123 → 126 exact,
   143 → 145 within one). But it collapses the 6-mer assignment (771/88 →
   702/24), so it is not selected. Counting and assignment must be decoupled,
   or the assignment made robust to end-lobe geometry, before it can be used.
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
7. **Assignment fragility — answered 2026-09-30.** The benchmark signal is
   the contrast with the dimmer terminal lobe. Without it, amplitude separates
   GlcNAc from GlcN at non-terminal positions with AUC ≤ 0.66, and shape, DFT
   molds and Fisher margin do not separate them at all. A robust and accurate
   assignment needs new information (reference oligomers of known
   composition, other imaging conditions), not another classifier variant.

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

### 2026-09-30 — Results: determinism adopted, capped support rejected for the benchmark

Four Viper jobs, all COMPLETED, all hashes match, all label-free checks pass:
A (12042058, `486ddec`, 2h01), B (12042063, `9a9f25d`, 2h31), B′ (12042065,
`9a9f25d`, 1h55) and long chains (12042066, `9a9f25d`, 2h00).

**Determinism: confirmed and adopted.** B and B′, run on different nodes at
different speeds, are identical byte for byte except output paths and logs:
all counts, fits, features, patches, predictions and grades. Run A is the new
reference for the promoted method: **137/145 exact N (145 within one),
771/870 correct, 88/145 exact chains, 865/870 classified, 94 errors**.
Per-scan calls give 702/33/863/161. Fusion: π = 0.342, θ0 = 0.035,
θ1 = 0.622. This is within the earlier time-limited spread (772, 772, 772,
769).

**Capped support (`half_maximum_cap`): rejected for the benchmark** by the
declared criterion.

| Run | Exact N per scan | Within one | Final exact N | Correct | Exact chains | Classified | Errors |
|---|---:|---:|---:|---:|---:|---:|---:|
| A (legacy support) | 123 | 143 | 137 | 771 | 88 | 865 | 94 |
| B (capped support) | 126 | 145 | 137 | 702 | 24 | 865 | 163 |

- Counting improves: per-scan exact N +3, within one +2. The registered-ROI
  refits are no longer needed; 95 scans agree with their consensus instead
  of 92.
- The assignment collapses. Errors at GlcNAc positions 2 and 5 rise from 15
  to 66 and from 43 to 70; positions 1, 3 and 6 slightly improve. Geometry
  changes on 52 of 146 scans (lobe shift q90 0.075 nm), 66 of 876 per-scan
  calls flip (22 of them on scans whose geometry did not change), and fusion
  then propagates the missed GlcNAc: π 0.342 → 0.245, θ0 0.035 → 0.073.
- The unsupervised cohort classifiers are highly sensitive to end-lobe
  geometry. Consistently with the edge-adjacency confound, the benchmark
  GlcNAc signal at positions 2 and 5 depends on how the end lobes are placed.

**Long chains (descriptive, 10–20mer config, capped rule).** Against the 12
manual reads (human judgements, not ground truth), mean |N − N_manual| is:

- per scan: 1.58 → 1.17;
- after consensus: 1.33 → 0.92.

Examples: `251206_013` 4 → 11 (12 read), `260220_083` 2 → 9 (read
rejected), `251206_034` 9 → 10 (9 read), `251206_038` 17 → 16 (18 read).
The fusion false-call rate is 0.12 instead of 0.45, but that comparison is
confounded: the earlier run used the non-corroborated assignment config.

**Decision.** Keep `global_maxtime = 0`; revert both counting configs to
`profile_quantile_noise`. The `half_maximum_cap` rule stays available in
code, inactive. Using it requires decoupling count selection from the
assignment geometry, or an assignment robust to end-lobe placement (open
questions 2 and 7). The previous "no gain from stabilizing counts" lessons
apply: count improvements do not transfer to the assignment through
cohort-wide learning.

### 2026-09-30 — Declared campaign: assignment robust to end-lobe placement

**Deliverable** (user request). GlcN/GlcNAc calls that stay stable when
end-lobe geometry changes and nothing else does. Time box: this session.

**Diagnosis** (label-free: position in the chain is known from the fit).

- In both A and B, the call rates of class 1 by position are 0.01–0.02 at
  chain ends, 0.51–0.57 next to an end, and 0.10 in the interior.
- A↔B agreement is lowest next to an end (252/292, against 286/292 at ends
  and 272/292 in the interior).
- Three of the four local features (`amp_prominence`, `amp_neighbor_ratio`,
  `integrated_prominence`) compare a lobe with the mean of its two
  neighbours. For the lobe next to an end, one of these is the terminal lobe,
  whose fitted amplitude and placement depend on the support.

**Variants, declared before any grade.** Same classifiers, seeds, vote,
corroborated training and fusion as the promoted method; only the feature
views change.

- **V1 interior-referenced.** BASE4 is replaced by `amp_int` and `integ_int`
  (log ratio to the median of the chain's non-terminal lobes) and by
  `amp_prom_int` and `integ_prom_int` (log ratio to the mean of the
  neighbours, where a non-terminal lobe never uses a terminal neighbour).
  Patch features are kept.
- **V2 shape only.** No amplitude features. GMM on descriptor, fwd/bwd molds
  and Fisher; k-means on split skew, backward-patch moments, descriptor and
  molds. Physical naming (brighter cluster = GlcNAc) is unchanged.

**Evaluation.** Each variant runs on the cached features of A (legacy
support) and B (capped support), with a local harness that reruns only the
classifier stages. On V0 it reproduces run A byte for byte.

- **Primary, label-free:** per-scan and fused call agreement between A and B,
  and class-1 rate by position class.
- **Secondary, external:** full145 own-N grade on A and B.
- Adoption requires higher A↔B agreement than V0. The grades are reported
  whatever they are, and nothing is tuned after grading.
- Expected trade-off: on NKNNKN, GlcNAc coincides with "next to an end", so a
  variant that no longer rewards that position may grade lower.

### 2026-09-30 — End-lobe robustness: no robust and accurate variant; the benchmark signal is edge adjacency

All runs use `test/assignment_view_variants.jl` on the cached features of A and B.
V0 reproduces runs A and B byte for byte. Results
(`results/assignment_robustness_20260930/`):

| Variant | Fused A↔B agreement | Class-1 rate A: end / next / interior | Final A (correct/exact/errors) | Final B |
|---|---:|---|---|---|
| V0 current | 784/876 | 0.02 / 0.81 / 0.09 | 771 / 88 / 94 | 702 / 24 / 163 |
| V1 interior-referenced | 859/876 | 0.32 / 0.94 / 0.96 | 487 / 0 / 378 | 485 / 0 / 380 |
| V2 shape only | 474/876 | 0.20 / 0.61 / 0.60 | 526 / 1 / 339 | 510 / 0 / 355 |

- V1 is stable but degenerate. Without neighbour contrasts against terminal
  lobes, the dominant structure is terminal versus non-terminal, and the
  brighter cluster (almost every non-terminal lobe) is named GlcNAc.
- V2 is close to random: patch shape alone gives no consistent clusters.
- Neither is adopted; V0 remains the method.

**Where the signal is** (external diagnostic, labels read for evaluation
only). Among the non-terminal lobes (positions 2–5) of the N = 6 benchmark
chains, the AUC of GlcNAc (2, 5) against GlcN (3, 4) is:

- `amp_prominence` (contrast with neighbours, including the terminal one):
  0.75;
- amplitude referenced to non-terminal lobes only (`amp_int`, `integ_int`,
  `amp_prom_int`, `integ_prom_int`): 0.62–0.66, and this may still include
  position effects;
- widths: 0.54–0.57;
- patch shape: `split_log_skew` 0.45–0.47, descriptor 0.43, `mold_cc_fwd`
  0.49, `mold_cc_bwd` 0.41, `emp_fisher` 0.52.

The same values hold on geometry B.

**Conclusion.** The benchmark accuracy comes almost entirely from the
contrast with the dimmer terminal lobe, i.e. from edge adjacency. The current
image features carry at most a weak position-independent GlcNAc signal
(AUC ≤ 0.66, amplitude only). The DFT molds, the Fisher margin and the patch
asymmetry carry none. An assignment that is robust to end-lobe placement
cannot also be accurate with these features. Chemical validation on the
10–20mers therefore needs new information rather than another classifier
variant, for example:

- reference oligomers of known composition (fully deacetylated or fully
  acetylated chains, or known non-palindromic sequences);
- other imaging conditions (bias, tip);
- chemistry-specific observables.

### 2026-10-04 — 6-mer focus: counting errors and position uncertainty

**Scope** (user decision: concentrate on the 6-mers until reference data
exist). Chemistry cannot be validated on NKNNKN, so the 6-mer work targets
what the benchmark can validate: counts, positions and their uncertainty.

**Counting errors of the deterministic reference** (run A; labels used only
to locate errors). There are 8 final errors.

- **5 under-counts (5 instead of 6), all molecules with fewer than 3 scans**,
  so consensus cannot act. Their visual causes:
  - `240307_019` is S-shaped and `240818_007` is kinked: the straight axis
    cuts across them;
  - `240312_Cu100081`: the robust-AICc guard moves N_eff 6 to 5;
  - `240313_Cu100058` has strong line noise;
  - `240310_Cu100032`: the cause is unclear.
- **3 over-counts: track 22** (`240817_017/018/019`). Two scans give 7, so
  the majority overwrites the correct scan.

Bending does not explain the per-scan errors overall: 22 per-scan errors are
spread evenly over the ridge-offset classes (17%, 15%, 11% of the scans in
the < 0.15, 0.15–0.25 and 0.25–0.40 nm bins). A curved centreline could
address at most 2–3 strongly bent scans. The capped support leaves all final
counts unchanged on the benchmark (identical in runs A and B); it only fixes
per-scan counts that consensus already repaired.

**Position uncertainty (label-free, new script
`test/report_position_repeatability.jl`).** Physical lobes of molecules imaged
several times (20 tracks, 121 lobes in run A) are aligned in the absolute frame
and their scatter measured.

| Alignment | All 20 tracks: along / across | Without track 26: along / across |
|---|---|---|
| Image registration only | 0.126 / 0.206 nm | 0.122 / 0.123 nm |
| + per-scan translation | 0.075 / 0.204 nm | 0.067 / 0.107 nm |
| + per-scan rotation (Procrustes) | 0.069 / 0.106 nm | 0.069 / 0.090 nm |

- Track 26, the 43-scan molecule, rotates by 25–37° on the surface between
  scans `240817_049` and `240817_050`. The rotation stays the same while the
  scan angle varies from +40° to −90°. Elsewhere, scan-angle changes of
  ±40–70° leave residual rotations of 1–3°, so the frame conversion is
  correct.
- **Lobe positions repeat to about 0.07 nm along the chain and 0.09–0.11 nm
  across it (1σ, ends included).** The consecutive-scan image registration
  adds about 0.1 nm.
- The spacing of neighbouring lobes is 0.66 nm on average and repeats to
  0.058 nm.
- This scatter is an empirical upper bound: it includes fit noise, piezo
  distortion between scans and real flexing.

Outputs: `results/support_determinism_20260930/position_repeatability_runA/`.

### 2026-10-04 — Position uncertainty integrated into the pipeline

`run_molecule_consensus_chitosan.jl` now ends with
`report_position_repeatability.jl`, which writes `positions.tsv` and
`position_repeatability/`. Each lobe gets its image and absolute position,
its physical lobe, and a 1σ uncertainty along and across the chain:

- `sd_source = track`: the lobe's own scatter, when its molecule has at least
  3 aligned scans;
- `sd_source = cohort`: otherwise, the pooled value for its position class
  (end or interior);
- `sd_source = none`: the run has no repeated molecule.

`grade_consensus_run.jl` checks the new table.

**Verification.** Viper job 12082827 (source `08c7087`, 1h44, 693 hashes
match) reproduces run A byte for byte on every earlier output (only output
paths and logs differ) and adds the six new files. `positions.tsv` (876
lobes: 651 `track`, 225 `cohort`) is identical to the standalone report on
run A. A 12-scan quickstart run (job 12082847, 16 min) passes all checks.

### 2026-10-04 — Local storage cleanup

The local checkout went from 93 GB to 9.6 GB.

- **`qe/`** keeps only the two accepted states, `glcn_restart5` and `glcnac`
  (SCF wavefunctions, LDOS cubes, logs; 6.4 GB). After the Viper cleanup these
  are the only copy, and they should be backed up. The outputs of the
  September DFT diagnostic campaigns, the preliminary runs and the
  intermediate relaxation steps were deleted (~64 GB).
- **`results/`** keeps the reference evidence:
  - `molecule_consensus_20260929`, `support_determinism_20260930`,
    `assignment_robustness_20260930` and `positions_20261004`;
  - the comparison profiles in `unit_assignment`, `promoted_counts_20260924`
    and `hybrid_reproduction_20260924`;
  - `reconstructed_cc_soft_v1` (benchmark symlinks) and
    `calibration_measurements_20260924`.
- The other 178 result folders, outputs of earlier experiments, were deleted.
  Journal paths to them no longer resolve locally. Their conclusions stay
  here, and their code is at the tag `archive/pre-cleanup-20260929`.
- No human review decision was lost: the review outcomes are in the tracked
  `benchmarks/` files.
