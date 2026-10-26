# Research journal

Dated record of decisions, experiments (including failures) and their reasons.
Add an entry for every experiment, scientific decision, failed approach, bug
fix or parameter change. Keep **Open questions** current.

**Notation (2026-10-22).** In the control NKNNKN, N = GlcNAc and K = GlcN.
Entries before 2026-10-22 call the class of positions 2 and 5 "GlcNAc"; those
units are GlcN. There, read "GlcNAc" as "class 1 (brighter)" and "GlcN" as
"class 0 (dimmer)", except in DFT and mold statements, which are defined
chemically.

The complete journal up to the 2026-09-29 consolidation (about 24,000 lines,
with every job, hash and intermediate result) is preserved at the git tag
`archive/pre-cleanup-20260929`:

```bash
git show archive/pre-cleanup-20260929:docs/src/journal.md
git show archive/pre-cleanup-20260929:docs/src/journal_archive.md
```

## Current state (2026-10-25)

**Method** (promoted): per-scan GCV counting on the capped support, with
fixed-N geometry refits on the legacy support
(`chitosan_decoupled_support.toml`, since 2026-10-25), then repeated-scan
molecule consensus with registered-ROI refits (`molecule_consensus.toml`),
then per-scan assignment learned on corroborated scans
(`unit_assignment_corroborated_training.toml`), then binomial latent-class
fusion. Entry point: `scripts/pipeline/run_molecule_consensus_chitosan.jl`.
The 10–20mer application uses `chitosan_10_20mer_decoupled_support.toml`.

**Benchmark** (fresh raw run, 12380719, `70a9baa`; full145 own-N): 137/145
exact N (145/145 within one), 776/870 correct, 85/145 exact chains, 865/870
classified, 89 errors. The legacy-support run (12042058, `486ddec`) gave
137/771/88/865/94, and the previous record 129/694/43/855/161. The decoupled
supports were adopted by the user after the grade, against the declared rule
(2026-10-25).

**Positions:** lobe positions repeat to 0.07 nm along and 0.09–0.11 nm across
the chain over repeated scans (2026-10-04).

**Data (2026-10-19):** all data live in the untracked `data/` folder.
- One copy of every raw file is in `data/raw/` (`NKN2/<session>`,
  `10_20mer/<session>`, `notes/`).
- Working sets are links in `data/sets/`: `6mer_Cu100_benchmark` (493 =
  145 graded + `240310_Cu100009` + 347 candidates validated by the user on
  2026-10-19), `10_20mer_curated` (19, since 2026-10-22) and
  `10_20mer_targets` (11 former targets in no experimenter list).
- `data/INDEX.tsv` lists 3,101 raw scans with their decision.
- The sets and the index are rebuilt from the tables in `benchmarks/` by
  `scripts/grading/report_data_sets.jl`.
- Discarded scans (656) have no folder, only a motif.
- 10–20mers (2026-10-21): 874 usable Cu(100) scans are in
  `data/sets/10_20mer_candidates` for the user's review (first pass: 279
  kept, 595 discarded), with the decisions in
  `benchmarks/chitosan_10_20mer_cu100_triage.tsv`. The user finds them
  mostly unusable. Since 2026-10-22 the working set is
  `data/sets/10_20mer_curated`: the 19 usable scans the experimenters
  selected (`benchmarks/chitosan_10_20mer_curated.tsv`); 865 candidates
  remain.
- Decoupled supports (2026-10-25, promoted): counting on the capped support
  and geometry on the legacy one remove the truncated 10–20mer supports. The
  counts are closer to the experimenters' length readings (mean |ΔN| 1.47 →
  1.00; run `dec1020_70a9baa`), which are a comparison, not a truth: 10–20mer
  lengths are not known. Against their acetylation readings the promoted run
  still recovers no unit order beyond a composition-preserving shuffle.
  The pattern statistics below are from the earlier run `app1020_5b706ca`.
- 10–20mer test run (2026-10-21, 440 scans): the predicted unit order does
  not match the 15 manual readings beyond chance once composition is
  accounted for. The only robust structure is the edge signature (terminal
  units in the dimmer class 0, alternation near both ends). Against 16
  acetylation readings by the experimenters (2026-10-22, untracked), no
  assignment recovers the order beyond a composition-preserving shuffle, and
  the readings have no end effect: the edge signature is positional.
- 10–20mer provenance (2026-10-22): the files are authentic copies of the
  share folder `Chitosan_10_20mer`, but nothing in the data documents the
  deposited compound. The "manual sequences" encode ring orientation, not
  acetylation. The user's partial triage keeps few scans outside December
  and January. The application is on hold until the sample is
  identified.
- The grading manifest is still the 145-scan
  `chitosan_6mer_counting_confirmed.toml`. The additions are graded in step
  1.5.

**Extended benchmark (2026-10-19, 493 scans, job 12282486):**
- Counting after consensus gives 391/493 exact N (79%) and 489/493 within
  one. The 145 originals are unchanged at 137/145; the 2024 additions reach
  68% and the 2025 additions 87.5%.
- The assignment is cohort-dependent. On the same 145 scans it drops from
  771/870 to 692/870 (28 exact chains) when the run also holds the 348
  added scans, and the fusion prevalence π rises from 0.34 to 0.45 (truth
  0.33). See the 2026-10-19 entry.
- Campaign B (2026-10-20/21) did not generalise to new data. Choosing the
  reference cohort by the benchmark score gave 316/528 on sealed 2025, or
  335 with a stronger search, against 349 for a label-free frozen model and
  361 for the promoted method. The stronger search did gain within 2024
  (+132 units on held-out sessions).
- Campaign A (2026-10-21) is a fully frozen, label-free assignment. Its calls
  are invariant across runs, and it gains on held-out 2024 sessions
  (1979/2430 and 188 exact chains, vs 1847 and 104). It transfers worse to
  2025 (305/528 vs 361). It is a candidate, not promoted.

**Per-unit reading and dataset selection (2026-10-25, campaigns R and D):**
- Per-unit signals repeat across scans of a molecule, but no more in the
  10–20mer than in the fixed-sequence 6-mer (the variance ratio interval is not
  above 1), and amplitude alone reads K against N with a ceiling of about 55%.
  The 6-mer, a position code, cannot validate an N/K reading.
- The 6-mer subset the assignment is learned on moves one 10–20mer call in
  five (noise floor: 81.4% agreement with the full-set model) and the bright
  fraction from 16.5% to 41.3%. No subset is chosen; the full set stays the
  frozen default.
- Counting is the measurable weak point of the 10–20mer: 18 of 65 molecules
  with repeated scans differ by more than 2 units across their scans.
- Elite subsets: the scan-level elite (167 scans) broke the molecules and gave
  44% in its own run (not a bug, diagnosis of 2026-10-25). The molecule-level
  elite v2 (26 molecules, 141 scans), graded on the 493-scan run, gives 601
  correct of 834 classified units (846 possible), 115/141 exact N and 38 exact
  chains. Its own run (Viper job 12387149, classes learned on those scans)
  gives 627 correct of 830 classified (846 possible), 105/141 exact N and 34
  exact chains, with no collapse (fusion π 0.25): +26 correct units, −10 exact
  counts, −4 exact chains against the same scans in the 493-scan run.

**Limits:** NKNNKN edge-adjacency confound (N = GlcNAc, K = GlcN). The
benchmark signal, the bright class on the two GlcN at positions 2 and 5, is
the contrast with the dimmer terminal lobe; without it the image features
separate the classes only weakly (AUC ≤ 0.66, amplitude) or not at all
(shape, DFT molds). The former naming "brighter = GlcNAc" is contradicted by
the control (2026-10-22). Also: exact chains in molecule-sized
blocks; counting parameters with benchmark-informed provenance; molds are
empirical templates; 10–20mer outputs not validated.

## Open questions

1. **Chemical validation.** Does the bright class detect a chemical
   difference, or the edge-adjacent position that NKNNKN happens to share
   with its two GlcN? Either way it falls on the GlcN, opposite to the former
   rule "brighter = GlcNAc" (2026-10-22). A control molecule with a
   non-palindromic or homo-oligomer sequence is needed.
   Campaign R (2026-10-25): amplitude alone cannot read K against N (ceiling
   about 55%). Experiments that can: homo-oligomers (all GlcN, all GlcNAc), a
   bias series, a non-palindromic or interior-K control.
   Position-neutral features recover edge versus interior, not chemistry.
2. **Long-chain counting vs assignment geometry.** The capped support rule
   (`half_maximum_cap`, 2026-09-30) fixes long-chain supports. Mean
   |N − N_manual| drops from 1.33 to 0.92, and `251206_013` goes from 4 to 11
   lobes (12 read). It also improves 6-mer per-scan counts (123 → 126 exact,
   143 → 145 within one). But it collapses the 6-mer assignment (771/88 →
   702/24), so it is not selected. Decoupling them (2026-10-25: count on the
   capped support, refit geometry on the legacy one) keeps every final 6-mer
   count and removes the truncated long-chain supports (mean |N − N_read|
   1.47 → 1.00 against 17 experimenter readings, a comparison only). The 6-mer assignment still moves through cohort
   learning (776 correct but 85 exact chains, against 771/88). The declared
   rule rejected it, and the user then adopted it (2026-10-25). What remains
   is open question 8.
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
8. **Cohort dependence of the assignment (2026-10-19).** The assignment is
   learned without labels on the run's own scans, so adding sessions changes
   the clusters and the fusion prevalence. On the 145-scan benchmark,
   accuracy goes from 88.6% (146-scan run) to 79.5% (493-scan run). Should
   the assignment be trained per session, normalised per session, or frozen
   on a reference cohort? Any of these is a new method campaign. Campaign B
   (2026-10-20) chose the reference cohort by its benchmark score. It failed
   its pre-registered rules: its in-sample gain was within the noise floor,
   and on sealed 2025 it scored 316/528, against 349 for the frozen pool and
   361 for the promoted method. A stronger multistart search (option 1,
   2026-10-21) generalised within 2024 (+132 units held out over the pool,
   33 molecules better and 19 worse) but not to 2025 (335 vs 349). Freezing
   on all corroborated reference scans without labels matched the
   transductive method (1871 vs 1847 held-out in 2024, 349 vs 361 in 2025).
   Campaign A (2026-10-21) freezes everything, Fisher included, on a capped
   label-free reference. It makes calls invariant (142/142 scans across runs,
   against 41/142 for the transductive method) and gains on held-out 2024
   sessions (1979 vs 1847). But it transfers worse to 2025 (305 vs 361).
   The exploratory diagnosis traces the 2025 loss to the per-track cap. The
   cap moves the class boundary, which helps 2024 and hurts 2025, so the
   best boundary differs between campaigns. Open: freeze per campaign, or
   calibrate the boundary per campaign without labels.
9. **A model built on the labelled data, label-free at inference (kept in
   mind at the user's request, 2026-10-21).** Campaign B chose training data
   by the benchmark score. It gained on held-out 2024 sessions but not on
   2025. A supervised model or selection needs more independent labelled
   molecules (or reference oligomers) from several campaigns before its
   transfer to new data can be shown.
7. **Assignment fragility — answered 2026-09-30.** The benchmark signal is
   the contrast with the dimmer terminal lobe. Without it, amplitude is higher
   at GlcN than at GlcNAc at non-terminal positions with AUC ≤ 0.66 (names
   corrected on 2026-10-22), and shape, DFT molds and Fisher margin do not
   separate them at all. The signed mold margin, which the method does not
   use, reaches AUC 0.56 in the chemical direction (2026-10-22). A robust and accurate
   assignment needs new information (reference oligomers of known
   composition, other imaging conditions), not another classifier variant.

10. **Counting reproducibility of the 10–20mer (2026-10-25).** 18 of 65
    molecules with repeated scans differ by more than 2 units across scans,
    while positions repeat. The instability comes from the extent detection of
    heterogeneous chains (a bright end blob). A label-free counting rule that
    reproduces across repeats is the next measurable target.
11. **Which 6-mer subset to freeze (campaign D, 2026-10-25).** The choice moves
    one call in five and the bright fraction between 16.5% and 41.3%; only the
    145 original scans fall beyond the sampling noise. Stability is not
    correctness. Open: report every application unit with its robustness
    across the 20 models, and decide the frozen subset only with new controls.

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

### 2026-10-04 — One plot per molecule in the final outputs

User request: for the shareable pipeline, a plot of every molecule in the
final result, sorted so that a scientist can see quickly what is exact, what
has the right N, and what is problematic.

- **Run (label-free).** The new final stage `test/plot_final_molecules.jl`
  draws one image per scan: the image the fit used (mean of the flattened
  forward and backward Z scans) around the chain, the final fused calls coloured with P(GlcNAc), the 2σ position ellipses from
  `positions.tsv`, and a title with N (and the per-scan N when the consensus
  changed it), the consensus outcome, the call string and the mean
  confidence. Plots go to `plots/all/` and to label-free folders
  (`consensus_agrees`, `consensus_corrected`, `not_checked_by_consensus`,
  `with_uncertain_calls`). This works on unknown data.
- **Grading (labels).** With `--data-dir`, `grade_consensus_run.jl` sorts the
  same plots into `1_all_exact`, `2_N_exact_calls_wrong` (wrong lobes
  circled), `3_N_wrong` and `4_not_in_benchmark`, with the truth in the
  title. The grading step reads label files with its own reader; the
  method's reader keeps refusing label columns.
- **On run A:** 88 / 49 / 8 / 1 plots, consistent with the report (88 exact
  chains, 137 exact N, 145 graded). The 146 plots take about 75 s.
- **Display.** A forward-only display made lobes look about 0.4 nm off the
  ridge on `240816_005`, whose forward and backward scans are offset. The
  plots therefore show the fused image that was fitted.
- **Viper.** The 12-scan quickstart with the integrated stage (job 12084468,
  `a80cad7`) still matches all 8 reference checksums and adds 12 plots.

### 2026-10-05 — 6-mer candidate scans: catalog and central folder

User request (roadmap step 1): look for other usable chitosan 6-mer images for
the benchmark, and gather every real, usable image once in one folder. Sources:
`/home/durif/Rebecca/data/data` (14 sessions, 935 scans) and the institute share
`(NKN)2` (not yet mounted here).

- **Catalog (label-free).** `test/catalog_sxm_candidates.jl` reads headers,
  pixels and the fit's molecule ROI; no fit, no count, no label. Statuses:
  `overview` (frame > 20 nm, the consensus collection limit),
  `object_crosses_frame` (the observed ROI touches two opposite borders),
  `excluded_by_review` (`benchmarks/chitosan_6mer_visual_exclusions.tsv`: the
  34 human rejections of the inventory, mostly non-chitosan) and `duplicate`
  (same bytes or same acquisition). Review flags: `partial_scan`,
  `roi_touches_unscanned`, `roi_touches_border`, `several_objects`,
  `fwd_bwd_mismatch`. `--copy-to` keeps a flat central folder with a
  provenance `MANIFEST.tsv`; a later source only adds the scans it lacks.
- **Rejected criterion.** Requiring at least 95% finite pixels excluded 31 of
  the 145 benchmark scans, which the pipeline handles. Partial scans are a
  flag instead. The frame-crossing status removes 54 scans (2×4 nm zooms on
  part of a chain, substrate steps), none of the 145.
- **Local result.** 935 scans: 679 usable, 195 overview, 54 crossing the frame,
  7 rejected by review (the other 27 rejections were already overview or
  crossing). No duplicate inside the tree; the 146 top-level files of
  `data/chitosan_6mer/` are byte copies of session scans. Central folder:
  `/home/durif/Rebecca/data/chitosan_6mer_usable/` (679 scans, 1.3 GB; moved
  to its `Cu100/` subfolder once the share was added, see below).
- **Candidates (grading side).** `test/report_benchmark_candidates.jl` splits
  the usable scans with the benchmark manifest and the earlier pre-assignment
  review: 145 in the benchmark (all usable), 56 reviewed before and not kept,
  478 never reviewed (194 without a flag, 284 flagged) in eight sessions
  (20240310–14, 20240815, 20241113–14). Many are repeat scans of the same
  molecules. Sheets: `results/benchmark_candidates_20261005/local_review/`.
- **Caveat.** Some earlier pre-exclusions used fit outputs
  (`pre_exclude_likely_nonchain_or_bad_guard_collapse`). Additions must be
  decided by eye, without method outputs: the truth NKNNKN only holds for one
  intact molecule, fully imaged.
- **Institute share.** `hive.fkf.mpg.de/esi/DATA/STM/STM_DATA/(NKN)2`,
  mounted by the user (gvfs), holds the same 14 sessions and 11 more. The
  share copies of the 14 sessions are byte-identical to the local tree (the
  679 central scans come back as duplicates, every other status repeats).
  The new sessions:
  - Cu(100), March–May 2025 (`Cu100LHe*`, 214 scans): 141 usable, 66
    overview, 7 crossing the frame. Added to `chitosan_6mer_usable/Cu100/`
    (820 scans).
  - Ag(111), July 2025 (six `2025*_LHe_Ag111` sessions, 286 scans): 160
    usable, 107 overview, 19 crossing. Kept apart in
    `chitosan_6mer_usable/Ag111/`: the counting calibration is a Cu(100)
    calibration, and the substrate changes the apparent heights and contrast.
  - Bias: 141 of the 145 benchmark scans are at −0.3 V; the new Cu(100) scans
    are mostly at −0.3 V (98 of 141), the Ag(111) scans spread over −0.05 to
    −0.3 V.
  `--root` (repeatable) was added because the gvfs path contains commas.
  Catalogs and sheets: `results/benchmark_candidates_20261005/share_*_catalog/`.
- **Next.** The human decision on the never-reviewed Cu(100) scans (478 from
  2024, 141 from 2025), then a grade of the frozen method on the additions
  before any re-optimisation (step 1.5). Ag(111) is a separate question: a
  new substrate would need its own calibration and its own benchmark.

### 2026-10-05 — Kept/discarded triage of the Cu(100) 6-mer scans (first pass)

User request: two folders, the scans kept for the benchmark and the scans
discarded, with a first triage done by the assistant; the user then checks
both folders and adjusts.

- **Scope.** Every unique (NKN)2 Cu(100) scan of the local tree and the share
  (1,149, deduplicated by sha256). Ag(111) is left out (another substrate, not
  benchmark material).
- **Rules.** Catalog statuses are applied as they are: overview → `vue_large`
  (261), frame crossing → `molecule_coupee` (61), inventory rejections →
  `revue_anterieure_exclu` (7). The 145 manifest scans are kept. The earlier
  human verdicts are kept: 46 ambiguous/excluded scans stay discarded
  (`revue_anterieure_ambigu`) and `240310_Cu100009` (accepted, outside the
  manifest) stays kept. The nine earlier pre-exclusions based on fit outputs
  were judged again by eye.
- **First visual pass.** 628 scans (619 never reviewed plus those nine) were
  judged on fixed-scale review sheets (8 × 8 nm window centred on the
  molecule, 2 nm bar). A scan is kept when it shows one isolated chain, fully
  imaged, with no artefact on the molecule. Motifs: `molecule_coupee` (edge or
  unscanned zone), `plusieurs_objets`, `pas_de_molecule`, `artefact_pointe`
  (line jumps, streaks, doubled image, tip change), `pas_6mere_long`, `pas_6mere_court`, and
  `incertain` for borderline cases. Blur alone is not a reason; the ROI length
  (benchmark scans: 3.6–4.5 nm, 5–95%) only served as an aid for clearly
  long or short objects, and for blurred elongated objects that were not
  clearly the same molecule as a kept series.
  Result: 380 kept, 248 discarded (90 `incertain`).
- **Totals.** 526 kept (145 benchmark + 1 earlier accepted + 380 new) and 623
  discarded: `vue_large` 261, `molecule_coupee` 101, `incertain` 90,
  `revue_anterieure_ambigu` 46, `pas_de_molecule` 36, `artefact_pointe` 35,
  `plusieurs_objets` 22, `pas_6mere_long` 21, `revue_anterieure_exclu` 7,
  `pas_6mere_court` 4.
- **Folders.** `test/report_benchmark_triage.jl` applies the table:
  `/home/durif/Rebecca/data/chitosan_6mer_Cu100/{gardes,ecartes}` (flat, 1.2 GB
  each, checksums verified), `tri.tsv` (decision, motif, source, note,
  source path) and `planches/gardes`, `planches/ecartes/<motif>` for review.
  The flat `chitosan_6mer_usable/Cu100` copy was removed (all 820 scans are
  in the new folders); Ag(111) moved to `chitosan_6mer_Ag111/`.
- **Caveats.** The kept set holds many repeat scans of the same molecules
  (for example about 40 scans of `250430_*`), so the number of independent
  molecules is much smaller than 526. The first pass is an assistant
  judgement on thumbnails: the `incertain` folder is where a human check
  matters most. Benchmark membership changes only after that check; then the
  frozen method is graded on the additions (step 1.5).

### 2026-10-18 — One reference folder for the 6-mer Cu(100) scans

User request: every 6-mer file in use centralised and listed in one place,
very clearly. It was not quite the case: the 146 cohort scans were in
`gardes/` but mixed with the 380 unvalidated candidates; the runs read them
through `full146_raw` symlinks into another copy (`data/chitosan_6mer/`); and
the list was split between the 145-file manifest, the symlinks (the 146th
scan, `240310_Cu100009`) and `tri.tsv`.

- **Layout.** `/home/durif/Rebecca/data/chitosan_6mer_Cu100/` now holds
  `gardes/benchmark/` (146: the 145 graded scans and `240310_Cu100009`, used by
  the runs but not graded), `gardes/candidats/` (380) and `ecartes/<motif>/`
  (623 in ten motif folders). `LISTE.tsv` has one row per scan (folder,
  decision, motif, who decided, note, session, acquisition date and time,
  frame, pixels, origin `<session>/<file>` in the (NKN)2 tree, sha256);
  `LISEZMOI.md` gives the count and meaning of every folder.
- **Single source.** The triage table gained `folder` and `origin` columns;
  `test/report_benchmark_triage.jl` copies every scan to `DEST/<folder>/` and
  writes the list, the readme and per-folder sheets, so the folder is always
  rebuilt from the table. Rebuilt from its own files (no share needed): 1,149
  scans, checksums identical, same set as before.
- **Runs.** Both `full146_raw` symlink sets under
  `results/reconstructed_cc_soft_v1/` now point to `gardes/benchmark/`
  (byte-identical files, same names); `data/chitosan_6mer_Cu100` links to the
  reference folder, and the docs use it.

### 2026-10-18 — All data inside the repository folder

User request: every dataset centralised in `Git/STMFit`, self-contained and
excluded from git.

- **Layout.** `data/` (already ignored by `.gitignore`, except its README and
  `.gitkeep`) now holds `chitosan_6mer_Cu100/` (the reference folder, moved
  from `/home/durif/Rebecca/data`), `chitosan_6mer_Ag111/` (moved) and `raw/`:
  the 2024 lab tree with spectra and annotations, the 10–20mer scans, the 25
  10–20mer targets and the loose files (`SequenceAnalysis.xlsx`, …). File
  counts match the sources (3,416, 6,503 and 25 files).
- **No duplication on disk.** `raw/` files are hard links to the originals in
  `/home/durif/Rebecca/data`, which stay untouched; the 935 reference scans of
  the 2024 sessions are hard links to their raw files (same sha256). The 214
  scans of the 2025 share sessions exist only in the reference folder.
- **Self-contained links.** Every symlink under `results/` is now relative and
  resolves inside the repository: 292 cohort links and 65 links to the
  10–20mer inputs (previously absolute paths into `/home/durif/Rebecca/data`),
  and 1,168 links of eight Viper run-input folders that pointed to `/ptmp`
  (dead locally), repointed by basename to the same 146 cohort scans. No link
  is broken and none leaves the repository.
- **Left as is.** `data/chitosan_6mer/` (1,081 files, the old copy of the lab
  tree and of the cohort) is now redundant: every file is identical to one in
  `raw/` or the reference folder. The absolute paths inside
  `benchmarks/chitosan_6mer_data_inventory.tsv` are historical metadata, not
  read to open files.

### 2026-10-18 — Data layout: raw store, working sets, one index

The first centralisation mixed molecule, substrate and processing stage at
the same level, kept duplicates (the old `data/chitosan_6mer/`, the 25
10–20mer targets copied out of their sessions, the 2024 scans twice), left
the 2025 and Ag(111) scans outside `raw/`, and stored each decision twice
(file location and table). The user agreed to a new layout, with one change:
discarded scans keep no folder, only the reason.

- **`data/raw/`** holds one copy of every raw file, by sample then session.
  - `NKN2/` has the 25 (NKN)2 sessions: 14 Cu(100) 2024 sessions with spectra,
    5 Cu(100) 2025 sessions (214 scans) and 6 Ag(111) sessions (160 scans).
    The 2025 scans moved there from the triage folders under their original
    names.
  - `10_20mer/` has the 37 sessions.
  - `notes/` has the manual annotations, `SequenceAnalysis.xlsx`, the slides and
    the profiles. The `251206_034.sxm` copy found among the loose files is the
    same bytes as the one in its session and was dropped.
- **`data/sets/`** are relative symlinks into `raw/`:
  `6mer_Cu100_benchmark` (146), `6mer_Cu100_candidates` (380) and
  `10_20mer_targets` (25).
- **`data/INDEX.tsv`** has one row per raw scan (2,975): path, sample,
  substrate, session, acquisition date and time, frame, pixels, bias, sha256,
  and status (set, `discarded` 623, `not_triaged` 1,800, `duplicate` 1, the
  `260222_070 - Copy.sxm` of the 10–20mer tree), motif, who decided, note.
- **Decisions** sit only in tracked tables.
  - The triage table now has `set` (empty means discarded, with `motif`) and
    `raw` columns.
  - The new `benchmarks/chitosan_10_20mer_targets.tsv` lists the 25 targets.
  - `test/report_data_sets.jl` checks every row against `raw/` by sha256 and
    rebuilds the sets and the index in about 40 s.
  - `report_benchmark_triage.jl` keeps only the review sheets.
  - `test/test_data_sets.jl` replaces `test/test_benchmark_triage.jl`.
- **Removed after a checksum check** (every `.sxm` identical in `raw/`):
  - the old `data/chitosan_6mer/` (1,081 files);
  - the triage folders `data/chitosan_6mer_Cu100/` and `data/chitosan_6mer_Ag111/`;
  - the copied targets `raw/chitosan_10_20mer_analysis/`.

  The first-pass review sheets moved to
  `results/data_review/6mer_Cu100_first_pass_20261005/`.
- **Run inputs.** The 1,488 symlinks under `results/` now point into `raw/`.
  None is broken and none leaves the repository.
- **Open.** The 2025 share sessions are incomplete in `raw/NKN2`: spectra are
  missing, and so are the 126 Ag(111) scans that were not usable. Copy them
  when the share is mounted.

### 2026-10-18 — `raw/NKN2` completed from the institute share

With the share mounted again (gvfs), `data/raw/NKN2` was compared file by file
with `(NKN)2`.

- **Copied.** 689 files that only the share had (0.36 GB), each verified by
  sha256 after copy and with its modification time kept:
  - the spectra (`.dat`) of the 2025 sessions;
  - the 126 Ag(111) scans that the catalog had not kept;
  - three `.ini` files.

  No common file differed in size.
- **Removed.** The only local-only files were three hidden
  `.sxmviewer_nanonis/` viewer caches (656 exported files, not measurements).
  They were removed from `raw/`; the originals stay in
  `/home/durif/Rebecca/data/data`.
- **Result.** `raw/NKN2` is now an exact mirror: 3,822 files on both sides.
  `data/INDEX.tsv` lists 3,101 raw scans. All 286 Ag(111) 6-mer scans and the
  1,640 10–20mer scans are `not_triaged`. The sets are unchanged (146, 380 and
  25).

### 2026-10-18 — Repository restructured; aggressive cleanup

User request: an optimal layout of the repository, with a fairly aggressive
cleanup. No scientific behavior changes.

- **Layout.**
  - `scripts/` holds the command-line programs by role:
    - `pipeline/` (19): the promoted method and its stages, label-free;
    - `grading/` (7): every script that reads labels;
    - `data/` (2), `calibration/` (2), `tools/` (1).
  - `test/` now holds only tests, with `test/runtests.jl` as the single entry.
    It runs each file in its own process, plus the package tests and the
    end-to-end run (`--fast` skips that run).
  - The former `test/lib/` is the precompiled `STMFit` library in `src/`.
    Its modules are `ScriptUtils`, `ReconstructedUnitAssignment`,
    `UnitAssignmentReport`, `EmpiricalFisherNative`, `CCMoldNative`,
    `PatchPreprocessing`, `MoleculeConsensus` and `MoleculeFusion`; scripts load
    them with `using STMFit.<Module>`.
  - DFT material is grouped in `dft/`: `inputs/` (was `hpc/qe_molds/`),
    `molds/` (was `templates/`) and `runs/` (was `qe/`, untracked). The
    promoted configs are untouched; their comments still name the old
    `test/` paths.
  - Old to new: `test/<script>.jl` → `scripts/<role>/<script>.jl`,
    `test/lib/<x>.jl` → `src/<x>.jl`, `test/lib/cc_mold_builder.py` →
    `test/reference/cc_mold_builder.py`.
- **Removed.**
  - `packages/STMMolecularFitGUI.jl` (unmaintained).
  - The package CLIs and legacy programs: `cli.jl`, `run.jl`, `examples/`, the
    old package READMEs and `.gitignore`s, and the `ArgParse` dependency.
  - The chunked-counting HPC path: `launch_remote.sh`, `batch_array.sbatch`,
    `merge_chunks.jl`, `remote.env.example`, and `summarize.jl`, used only by
    it. The HPC path is now one `sbatch`.
  - Two finished research harnesses: `assignment_view_variants.jl`, and
    `compare_reconstructed_champion.jl` with its test.
  - Root leftovers: agent-tool state (`.agents`, `.codex`, `.codegraph`,
    `.cortexkit`, `.slim`), caches, and the July debug journal (its mold
    conclusion is recorded above).
  - The unused root dependencies `IntervalArithmetic` and `StatsBase`.
- **Environment.** The stdlibs used by the library were added to the root
  project. The Manifest was re-resolved offline with Julia 1.13.2, the user's
  rule being to keep the latest release: package versions are unchanged, and
  only stdlib JLLs moved.
- **Verification.**
  - The pre-refactor suite passed in full: 19 files plus the end-to-end run.
  - After the refactor, `test/runtests.jl` passes all 18 runs (17 test files and package tests, plus the end-to-end run).
  - One test had reached `ScriptUtils` through a script module; it now uses
    `STMFit.ScriptUtils`.
  - Byte-for-byte check on real data, done on Viper with Julia 1.13.2 at
    `b8986a0` (jobs 12272560 and 12272561):
    - The 12-scan quickstart reproduces all 8 reference checksums.
    - The full146 run matches run A (`486ddec`) and the run with positions
      (`08c7087`) on 361 of 367 common files.
    - The remaining differences are only absolute-path columns (`best_plot`,
      `file_dir`, `path`) in the counting summaries and `input_hashes.tsv`.
    - The molecule plots are byte-identical between the pre- and
      post-refactor scripts, and between this machine and Viper. They differ
      from the `a80cad7` quickstart only through the earlier fused-image fix
      (`a977a3a`).
    - The outputs are in `results/{quickstart,full146}_refactor_b8986a0/`.
- **Viper Julia.** Julia 1.13.2 is installed in
  `/viper/u1/<user>/software/julia-1.13.2` (official tarball, checksum
  verified); the 1.13.0 install is no longer used. The user opens the Viper
  connection with `mpcdf-ssh viper`.
- **`release-prep`** (local branch) was rebased onto this layout. Its archive
  list, READMEs, examples and containers use the new paths.

### 2026-10-19 — Review page for the 380 benchmark candidates

User request: help deciding which 6-mer candidates enter the benchmark. A
candidate enters only if it shows one intact (NKN)2 molecule, isolated and
fully imaged, because the benchmark truth is NKNNKN with 6 units. Fragments,
aggregates, cut molecules and non-chitosan objects stay out.

- **Tool.** `scripts/grading/report_candidate_review.jl` builds a local page
  (`index.html`, no upload) with one card per candidate. Each card shows the
  whole frame, an 8 × 8 nm window and the suggestion, with buttons to decide
  and an export of the decisions as TSV.
  - Repeat scans of one molecule are grouped (same session, ≤ 30 min,
    frame centres ≤ 5 nm apart, from the consensus collection settings): 88
    groups for the 380 candidates.
  - One benchmark scan per session is shown for reference.
  - GR plotting is not thread-safe, so the images are drawn serially.
  - The template is `scripts/grading/candidate_review_template.html`.
- **Second visual pass** (assistant suggestions, saved in
  `results/data_review/suggestions_candidats_6mer.tsv`): 347 suggested for the
  benchmark, 33 suggested to discard.
  - Two big lobes only, with no recognisable chain (`240310_Cu100046`–`053`,
    `240313_Cu100031`).
  - Blurred spots: `240312_Cu100006`, `014`–`017`, and `250306_009`–`012`.
  - Folded V shapes, which could be a folded 6-mer or fragments:
    `240815_104`–`110` and `240312_Cu100073`.
  - A lobe detached from one end (`240311_Cu100015`–`018`).
  - A ragged, unstable-tip image (`240815_022`–`023`).
  - A chain against another object (`241113_034`).
  - Blurred scans of a molecule already in the benchmark were kept with it
    (`241114_029`–`030`).
- **Next.** The user's decisions, exported from the page, go into the triage
  table (source `revue_utilisateur`). The data sets are then rebuilt, and the
  frozen method is graded on the additions (step 1.5).

### 2026-10-19 — Candidates decided; HPC account removed from the docs

- **Candidates.** The user accepted the suggestions of the review page: 347
  candidates enter the benchmark set and 33 are discarded (motif `incertain`,
  or `plusieurs_objets` for one), with source `validation_utilisateur`.
  - `data/sets/6mer_Cu100_benchmark` now holds 493 scans; the candidate set is
    empty and was removed.
  - The 145-scan grading manifest is unchanged. Grading the frozen method on
    the 347 additions needs an extended manifest and a Viper run (step 1.5).
    The 2025 sessions should be reported apart, as they are the only new
    molecules; the 2024 additions are mostly repeat scans of benchmark
    molecules.
- **HPC account.** The repository is public. At the user's request, the MPCDF
  account name was replaced by `<user>` in every tracked file (paths
  `/u/<user>/`, `/ptmp/<user>/`, `/viper/u1/<user>/`). No password, token, key
  or OTP material was ever committed: the history was checked. Access needs
  the user's password and OTP, which stay on the local machine. Earlier
  commits still contain the account name; it is not a credential, and
  rewriting the public history was not done.

### 2026-10-19 — Step 1.5 launched: frozen method on the 493-scan benchmark

- **Manifest.** `benchmarks/chitosan_6mer_benchmark_v2.toml` lists the 493
  scans of `data/sets/6mer_Cu100_benchmark`, with a `cohort` field: `full145`
  (145), `added_2024` (260, including `240310_Cu100009`) and `added_2025`
  (88). Keys are the set names, and `raw_path` gives the raw file. The
  145-scan manifest stays the default.
- **Grader.** `grade_consensus_run.jl` gains `--manifest` and `--cohort`.
  Checked on the verified full146 run: the default grade reproduces
  771/88/865/94, and `--manifest v2 --cohort full145` gives identical
  grades.
- **Run.** Viper job 12278834 runs the frozen pipeline (code `b8986a0`,
  identical to HEAD for the pipeline, packages, configs and molds; Julia
  1.13.2) on the 493 scans. The inputs, uploaded to
  `/ptmp/<user>/stmfit/bench493_v1_inputs/raw`, are checksum-identical to the
  local set.
- **Two caveats** for reading the results:
  - The consensus now sees more repeat scans, so the full145 grades may
    differ from run A; that difference is the effect of the larger cohort,
    not a code change.
  - The 88 scans of 2025 are about a dozen molecules, one of them imaged
    45 times, so results are also reported per molecule.

### 2026-10-19 — Step 1.5, first run failed: geometry refit rejected flagged fits

- **Failure.** Viper job 12278834 failed after 3 h 28 in the scan-geometry
  stage, on 12 of the 493 scans (`240313_Cu100023`, `026`, `066`;
  `240314_Cu100_043`; `241113_045`, `097`, `129`, `130`; `250430_022`, `028`,
  `029`, `041`).
  - For these scans the counting step had selected an N (for example 6 for
    `250430_022`, via `support_midpoint_down`) whose fits all carry a quality
    flag (`residual high`).
  - `batch_full.jl` falls back to the best fit in that case, but
    `extract_lobe_features.jl` kept only unflagged fits. It found no fit and
    skipped the scan, and the orchestrator then stopped on the coverage check
    (no partial result).
  - The 146-scan cohort never hits this case.
- **Fix** (`extract_lobe_features.jl`). When the counted N has only flagged
  fits, the geometry describes that same model: the best flagged fit, with
  source `ell_flagged` or `circ_flagged`.
  - The code path runs only where a run used to stop, so every run that
    completed before is unchanged, the verified full146 run included.
  - All 12 scans now give N lobes.
  - A unit test covers the fallback, and all tests pass.
- **Rerun.** The counting summary of the failed run is reused
  (`STMFIT_SELECTED_SUMMARY`; the counting code is unchanged and
  deterministic).

### 2026-10-19 — Step 1.5 results: the frozen method on the 493-scan benchmark

- **Run and grades.**
  - The run is Viper job 12282486 (`a2bcb15`, 1 h 32; counting reused from
    job 12278834). Outputs are in `results/bench493_a2bcb15/`; grades by
    cohort are in `results/bench493_a2bcb15_grade/{all,full145,added_2024,added_2025}/`.
  - Plots sorted by outcome are in `all_with_plots/plots/`.
  - `report_unit_assignment_benchmark.jl` now checks that the graded
    denominator equals the manifest (its files times the control length):
    still 145/870 for the default manifest.
- **Counting** (exact N, and within one, before → after consensus):

  | Cohort | Scans | Per scan | After consensus | Within one |
  |---|---:|---:|---:|---:|
  | full145 | 145 | 123 (84.8%) | 137 (94.5%) | 145 |
  | added_2024 | 260 | 170 (65.4%) | 177 (68.1%) | 256 |
  | added_2025 | 88 | 68 (77.3%) | 77 (87.5%) | 88 |
  | all | 493 | 361 (73.2%) | 391 (79.3%) | 489 |

- **Assignment** (fused calls; correct lobes / possible, exact chains,
  errors):

  | Cohort | Correct | Exact chains | Errors | Per-scan calls (correct) |
  |---|---:|---:|---:|---:|
  | full145 | 692/870 (79.5%) | 28 | 173 | 666 (76.6%) |
  | added_2024 | 1155/1560 (74.0%) | 76 | 343 | 1136 (72.8%) |
  | added_2025 | 361/528 (68.4%) | 17 | 165 | 321 (60.8%) |
  | all | 2208/2958 (74.6%) | 121 | 681 | 2123 (71.8%) |

- **Per molecule** (118 consensus tracks; molecule means):

  | Molecules | Count | Scans | N exact | Assignment |
  |---|---:|---:|---:|---:|
  | Only original scans | 25 | 117 | 88.0% | 78.5% |
  | Original plus 2024 additions | 15 | 127 | 79.2% | 74.8% |
  | Only 2024 additions | 62 | 161 | 63.9% | 67.9% |
  | 2025 (new) | 16 | 88 | 75.0% | 78.3% |

  The 2025 per-scan figures are dominated by `250430`. The consensus splits
  this molecule into a 36-scan track (N always right, assignment 50%) and an
  8-scan track (N wrong).
- **Findings.**
  1. **Counting.** It is stable on the original scans but clearly lower on
     scans that were not in the benchmark. Part of the original 145 was
     chosen with fit outputs (pre-assignment rows used `N_selected`), so its
     94.5% likely overstates accuracy on unselected scans. Some additions may
     also not be intact 6-mers (first-pass triage); the N-wrong plots let the
     user check this before any conclusion on accuracy.
  2. **Assignment.** It is cohort-dependent. The same 145 scans fall from
     771/870 (88 exact chains) to 692/870 (28 exact chains) once the run
     holds the additions. The corroborated k-means/GMM and the fusion EM are
     fitted on the run's own scans; with the added sessions, π goes from
     0.34 to 0.45 (truth 0.33), i.e. too many GlcNAc calls.
  3. **Consequence.** Assignment accuracy is not a property of the method
     alone. This matters for the 10–20mer application, which also mixes
     sessions (open question 8).
- **Headlines.** The 145-scan headline (771/88/865/94) remains the result of
  the 146-scan cohort run and stays as cited. The 493-scan figures are
  recorded here and in Current state.

### 2026-10-20 — Declared campaign B: frozen assignment, reference cohort chosen by its benchmark score

- **Deliverable.** A frozen unit assignment whose reference cohort is the one
  that maximises the absolute number of correct units, with:
  - an honest estimate (nested leave-one-session-out);
  - a noise floor (the same search on random truths);
  - one verdict on sealed sessions.

  Done when the merged report, the sealed verdict and a results entry exist.
  Scope approved by the user on 2026-10-20 ("focus sur B en recherche de
  l'optimum, avec la cohorte optimum, on fera A"), time box three days. A
  (reference fixed without labels) is deferred.
- **Status of B.** The model reads no label at inference. The choice of its
  training scans uses the benchmark labels, so B is a supervised selection of
  training data. Under the scientific invariants it is an experiment and a
  comparison, not a label-free method. Promoting it would need an explicit
  decision by the user.
- **Engine** (`src/frozen_assignment.jl`, module `FrozenAssignment`).
  - It learns the cohort-level parts of the promoted assignment on given
    scans: k-means centroids of the 4 views, GMM of `v_cc`, physical naming.
    It then applies them unchanged to any scan. Features stay standardised
    within each scan, so a lobe's features do not depend on the other scans.
  - Trained on the corroborated scans of the 493-scan run, it reproduces
    that run exactly: `pred_kmeans`, `pred_gmm`, vote and fused calls
    (0 differences out of 2,927 lobes), and fusion π/θ0/θ1 to the printed
    digits.
  - The GMM computes each Cholesky factor once per iteration and uses
    buffers. It is bitwise identical to the production script on random
    subsets, at about 2 s per training instead of minutes.
  - `test/test_frozen_assignment.jl` covers this against the production
    scripts on synthetic data, plus fusion, model I/O and the search.
- **Config** `config/unit_assignment_frozen_reference.toml`. It is the
  corroborated config plus:
  - `assignment_model = "frozen_reference"`;
  - the fusion EM keys;
  - `fusion_emission_source = "reference_sessions_all_scans"`: θ0, θ1 fitted
    on the frozen model's per-scan calls of every scan of the reference
    sessions;
  - `fusion_prior_pi = 0.5`: a neutral prior, so no GlcNAc fraction is
    carried from the reference molecules.
- **Protocol** (`benchmarks/reference_cohort_search_v1.toml`, fixed before any
  result).
  - **Split.**
    - Development: the 13 sessions of 2024 (405 scans).
    - Sealed: the 5 sessions of 2025 (88 scans), read once by `--mode sealed`
      after the final model is written.
  - **Unit.** A physical molecule (consensus track; tracks never cross
    sessions). Its training scans are its corroborated scans. 39 development
    tracks have one; the others cannot change a model.
  - **Objective.** Correct final (fused) units, summed over the held-out
    session of every leave-one-session-out fold, with the
    `grade_unit_assignment.jl` rule. The rule reproduces the official
    2208/2889/121 on the 493 run.
  - **Search.** Floating backward elimination from the full pool. Each step
    removes the unit whose removal gains most, then re-adds removed units
    while that gains. A step needs a gain of at least 1. Ties go to the
    smallest track id.
  - **Honest estimate.** Nested: each development session is held out in
    turn, the whole search runs on the other 12, and the selected cohort and
    the full pool are graded on the held-out session.
  - **Noise floor.** The final search is repeated with per-track random
    truths among the palindromic 2-of-6 patterns (100001, 010010, 001100),
    for 20 seeds.
  - **Secondary criterion** (added the same day, before any result, after a
    review note). Repeat scans of one physical molecule are not independent:
    the 2025 track of `250430` holds 36 of the 88 sealed scans. Every report
    therefore also gives one vote per molecule: the mean over tracks of the
    fraction of correct units, and the number of tracks with a strict
    majority of exact chains.
- **Decision rules** (fixed now).
  - The search counts as a real gain only if three conditions hold:
    1. the nested gain over the pool on held-out sessions is positive;
    2. the in-sample gain of the final search exceeds every null gain;
    3. the per-molecule criterion does not fall.
  - The sealed verdict compares the selected model, the pool model and the
    promoted transductive method on the same 88 scans (361/528 in the
    493-scan run).
  - No parameter of the protocol changes after a result.
- **Limitation.** The empirical Fisher feature is computed once, label-free
  and transductively, on the 493-scan run (lobe-parity CV, 2025 scans
  included). It is not re-learned per cohort; only k-means, GMM, naming and
  emission rates follow the cohort.
- **Disclosure.** A local timing run of the final search was stopped after two
  steps. It showed 1871 correct units for the pool and 1981 after two
  removals (in-sample). The protocol was already written and is unchanged.
- **Execution.** `hpc/run_reference_cohort_search.sbatch` runs the 34
  searches as independent processes on one Viper node (16 × 8 threads),
  then `--mode report`. Results follow in a separate entry.

### 2026-10-20 — Campaign B results: the score-chosen reference cohort does not generalise

- **Runs.**
  - Viper job 12301311 (`e67033a`) ran all 34 searches on one node in 1 h 59
    (the final search alone needed 3,025 trainings). Outputs are in
    `results/cohort_search_e67033a/`.
  - The merged report is in `results/cohort_search_e67033a_report/`, rebuilt
    locally with `e23ef44`, which adds the transductive reference and one
    vote per molecule pooled over held-out sessions.
  - The sealed verdict is in `results/cohort_search_e67033a_sealed/`. It was
    graded once, and `report_unit_assignment_benchmark.jl` gives the same
    numbers.
- **Nested estimate** (the 13 development sessions, each held out once; 405
  scans, 102 molecules):

  | Reference cohort | Correct | Exact chains | Errors | Molecule accuracy | Molecules majority exact |
  |---|---:|---:|---:|---:|---:|
  | Pool (all 39 units, frozen) | 1871/2430 | 118 | 492 | 0.712 | 20 |
  | Selected by the search | 1964/2430 | 176 | 399 | 0.723 | 17 |
  | Promoted method, transductive 493-scan run | 1847/2430 | 104 | 516 | 0.698 | 17 |

  - The held-out gain (+93 units) comes mostly from 20240817 (+43) and
    20240313 (+27).
  - Three sessions lose units (20240815 −11, 20240314 −4, 20240310 −1).
- **Final search** (13 sessions, in-sample).
  - It goes from 1871 to 2010 correct units (+139).
  - It removes 5 tracks (23, 25, 39, 46, 22) and keeps 34.
  - Molecule accuracy goes from 0.712 to 0.748.
- **Noise floor.** On random palindromic truths, the 20 null searches gain 36
  to 185 units in-sample (median 112). The real gain ranks above 13 of the 20.
- **Stability.**
  - Each fold's search removes 2 to 8 units, and 29 of the 39 units are
    removed in at least one fold.
  - Tracks 23 and 22, removed by the final search, are removed in no fold.
- **Sealed verdict** (2025, 88 scans, 16 molecules, read once):

  | Reference cohort | Correct | Exact chains | Errors | Molecule accuracy | Molecules majority exact |
  |---|---:|---:|---:|---:|---:|
  | Pool (frozen) | 349/528 | 10 | 177 | 0.729 | 3 |
  | Selected | 316/528 | 10 | 210 | 0.729 | 3 |
  | Promoted method, transductive | 361/528 | 17 | 165 | 0.771 | 6 |

  The whole loss of the selected cohort is the 36-scan track of `250430`
  (108 → 72 correct). The other 15 molecules gain 3 units, and the molecule
  vote is unchanged.
- **Decision** (rules fixed in the declaration).
  1. Nested gain positive: yes (+93).
  2. In-sample gain above every null gain: **no** (13 of 20).
  3. Per-molecule criterion: mixed. Accuracy is +0.011, but molecules with
     a majority of exact chains fall from 20 to 17.

  The search therefore does not count as a real gain, and the sealed
  sessions agree. The selected model is not used.
- **Paired per-molecule check** (added after the verdict, at the user's
  question; selected minus pool, correct units per physical molecule).
  - **Nested 2024.** Of 102 molecules, 21 are better, 18 worse and 63
    unchanged. The total is +93, with a bootstrap 95% interval over
    molecules of [−1, +214]. Track 58 alone gives +43; without it the gain
    is +50.
  - **Sealed 2025.** Of 16 molecules, 2 are better, 2 worse and 12
    unchanged. The total is −33, interval [−108, +8]. The whole loss is one
    fused physical lobe (`track117_5`, called 1 by the pool and 0 by the
    selection) repeated over 36 scans.
- **Interpretation.**
  - The effects in both directions come from single molecules. At the
    molecule level the selection and the pool are a draw. The correct
    conclusion is "no evidence of a gain", not "a demonstrated loss".
  - The in-sample gain is within what the same search finds on random
    truths, and the choices change from fold to fold. This is consistent
    with the winner's curse discussed before the campaign.
  - One caveat on rule 2: the null searches start from lower accuracy
    (1120–1531 vs 1871 correct), so their gains are not strictly on the same
    scale.
- **For campaign A.** The frozen pool model is learned on every corroborated
  scan of the reference sessions, without any label. It is as good as the
  transductive method:
  - 1871 vs 1847 on held-out 2024 sessions;
  - 349 vs 361 on 2025;
  - molecule accuracy 0.712 vs 0.698 (2024) and 0.729 vs 0.771 (2025).

  Freezing removes the cohort dependence at no clear cost, so it is the
  natural starting point for A.
- **Headlines.** Unchanged; nothing is promoted.

### 2026-10-20 — Declared: campaign B option 1, a stronger cohort search

- **Question** (raised by the user after the v1 results). Is the v1 cohort a
  poor local optimum? Would a stronger search find a cohort that generalises
  better? The user approved testing it ("lance l'option 1"), inside the
  campaign time box.
- **Protocol** (`benchmarks/reference_cohort_search_v2.toml`, fixed before any
  v2 result). Data, split, unit, objective and secondary criterion are those
  of v1. Only the search changes:
  - **Moves.** Each step takes the best flip: remove or add any unit.
  - **Swaps.** When no flip gains at least 1 unit, the search tries swaps:
    one removal plus one addition, among the 5 best removals and the 5 best
    additions of that flip round.
  - **Starts.** Six starts: the full pool, plus five seeded random subsets
    that keep each unit with probability 0.75.
  - **Selection.** The best end of the six (in-sample objective; ties go to
    the lowest start) is the selected cohort, both in the final search and
    inside every nested fold.
  - **Null.** 4 null seeds (random palindromic truths) run the same
    multistart search.
- **Evaluation.**
  - **Evidence.** The nested estimate on the 13 sessions of 2024, compared
    with the pool and with v1 (+93 units held out).
  - **2025.** The 2025 sessions were already read once by v1, so here they
    are a secondary check only.
- **Decision rule** (protocol `[decision]`). The stronger search supports
  "a better optimum generalises" only if all three hold on the held-out
  sessions:
  1. its gain over the pool exceeds +93 units;
  2. the bootstrap 95% interval over molecules excludes 0;
  3. more molecules improve than worsen.
- **Execution.** `hpc/run_reference_cohort_search.sbatch` now takes
  `STMFIT_PROTOCOL`. Each start is its own single-thread process, which
  makes 108 processes on one Viper node. `--mode report` then picks the best
  start of every target and grades the held-out sessions.

### 2026-10-21 — Option 1 results: the stronger search generalises within 2024, not to 2025

- **Run.**
  - Viper job 12302813 (`6461ce4`) ran 108 single-thread processes on one
    node in 36 min. Single-thread processes are about 10 times faster per
    training than the 8-thread processes of v1, presumably because threads
    share garbage collection.
  - Outputs are in `results/cohort_search_v2_6461ce4/` (report in
    `report/`). The 2025 check is in
    `results/cohort_search_v2_6461ce4_check2025/`, confirmed by
    `report_unit_assignment_benchmark.jl`.
- **Nested estimate** (2024 sessions held out, 405 scans, 102 molecules):

  | Reference cohort | Correct | Exact chains | Errors | Molecule accuracy | Molecules majority exact |
  |---|---:|---:|---:|---:|---:|
  | Pool | 1871/2430 | 118 | 492 | 0.712 | 20 |
  | v1 selection (floating) | 1964/2430 | 176 | 399 | 0.723 | 17 |
  | v2 selection (multistart flip/swap) | 2003/2430 | 175 | 360 | 0.742 | 15 |
  | Promoted method, transductive | 1847/2430 | 104 | 516 | 0.698 | 17 |

  - **Paired per molecule** (v2 selection minus pool): 33 better, 19 worse
    and 50 unchanged. The total is +132, with a bootstrap 95% interval of
    [29, 260]. Track 58 gives +43; without it the gain is +89.
  - **Pre-registered rule.** All three conditions hold (+132 > +93; interval
    excludes 0; 33 > 19). Within 2024, a better optimum does generalise to
    held-out sessions.
  - The number of molecules with a majority of exact chains still falls
    (20 → 15). It was not part of the v2 rule.
- **Final search** (13 sessions).
  - The six starts end at 2010, 2026, 2017, 2027, 2038 and 2000; the best is
    start 4.
  - That cohort keeps 24 of 39 units. Among the dropped tracks are three of
    the largest (43, 47, 44).
  - The in-sample gain is +167. The 4 null searches gain 169, 254, 224 and
    75, so the real gain ranks above only one of them.
- **2025 check** (secondary: the sessions were already read by v1):
  - The v2 selection gives 335/528 (9 exact chains, molecule accuracy 0.708).
  - The pool gives 349 (10 exact, 0.729) and the promoted method 361 (17
    exact, 0.771).
  - Per molecule: 2 better, 4 worse and 10 unchanged; total −14, interval
    [−44, +6].
- **Reading.**
  - A stronger search does find cohorts that score better on held-out 2024
    sessions, so the user's suspicion that v1 was a poor optimum is
    confirmed within 2024.
  - The benefit does not reach 2025, a later measurement campaign, where the
    selection is again below the pool and the promoted method.
  - The 10–20mer targets come from December 2025 to March 2026, so transfer
    to new campaigns is what matters for the application.
  - Choosing training data by the score therefore fits structure shared by
    the 2024 sessions, not structure that carries to new data.
  - The in-sample gains remain inside the null range.
- **Decision.** No score-chosen cohort is adopted. The label-free frozen pool
  remains the candidate for campaign A. Headlines unchanged.

### 2026-10-21 — Declared campaign A: fully frozen, label-free assignment

- **Deliverable.** A unit assignment whose cross-scan parts are learned once
  on a reference cohort fixed without labels, then applied unchanged to any
  scan.
  - Done when four things exist: a build/apply path in `scripts/pipeline/`;
    an invariance check; the held-out grade on 2024 sessions; the 2024 →
    2025 transfer grade.
  - Approved by the user on 2026-10-21 ("lançons le modèle A"). They also
    asked to keep in mind a model built on the labelled data but label-free
    at inference (open question 9).
- **What is frozen.** Every part learned across scans:
  - the empirical Fisher feature: one model per lobe parity, fitted on the
    reference lobes, scoring the other parity, as the production lobe-parity
    CV does. With the whole run as reference it reproduces both runs'
    `fisher_cv.tsv` (2927/2927 and 876/876 lobes);
  - the k-means centroids, the GMM components and the physical naming;
  - the fusion emission rates, fitted on every scan of the reference
    sessions.

  All other features are computed scan by scan.
- **Config** `config/unit_assignment_frozen_labelfree.toml`. It is the
  frozen-reference config plus five keys:
  - `reference_scans = "corroborated_counts"`;
  - `reference_max_scans_per_track = 5`, first scans in `file_name` order:
    repeat scans of one molecule are not independent, and v1/v2 showed that
    heavy tracks dominate;
  - `frozen_fisher = "reference_lobe_parity"`;
  - `frozen_gmm_weights = "equal"`: the GMM mixing weights are set to 1/2
    after learning, so no class fraction is carried to new molecules (as
    `fusion_prior_pi = 0.5`).

  These choices are fixed before any grade of A. A is graded once, with no
  variant.
- **Evaluation** (`scripts/grading/grade_frozen_assignment.jl`, on the
  493-scan run).
  1. **Held-out 2024.** Each 2024 session is graded with a model built
     without it: reference rule applied to the other sessions, emission
     rates from all their scans. The comparison is the promoted transductive
     method on the same scans.
  2. **Transfer.** A model built on all 2024 sessions is applied to the 2025
     sessions. These were already read by campaign B, so this is secondary.
  3. **Invariance.** The 2024 model is applied to the 145 original scans as
     they appear in the 146-scan run and in the 493-scan run. Per-scan calls
     must be identical wherever a scan's inputs are identical.
- **Decision rule** (fixed now). A is a candidate to replace the transductive
  assignment if invariance holds, and if on held-out 2024 sessions it gives
  at least as many correct units and at least the same molecule accuracy as
  the transductive method. Promotion would still need the user's decision
  and a full pipeline run.

### 2026-10-21 — Campaign A results: frozen and label-free, better on held-out 2024, worse on 2025

- **Run.** Local, single shot: `grade_frozen_assignment.jl` (`8e6a846`) on
  the 493-scan run, with `--other-run` set to the 146-scan run. Outputs are in
  `results/frozen_labelfree_8e6a846/`.
- **Held-out 2024 sessions** (each graded with a model built without it;
  405 scans, 102 molecules):

  | Assignment | Correct | Exact chains | Errors | Molecule accuracy | Molecules majority exact |
  |---|---:|---:|---:|---:|---:|
  | Frozen, label-free (A) | 1979/2430 | 188 | 384 | 0.736 | 23 |
  | Promoted method, transductive 493-scan run | 1847/2430 | 104 | 516 | 0.698 | 17 |

  - **Paired per molecule** (A minus transductive): 32 better, 14 worse,
    56 unchanged. The total is +132, with a bootstrap 95% interval of
    [41, 251].
  - **Comparison with campaign B.** Without any label, A matches the
    score-chosen cohorts on these sessions (v1 1964, v2 2003) and has the
    most exact chains.
- **Transfer** (model built on all 2024 sessions, applied to 2025; these
  sessions were already read by campaign B):

  | Assignment | Correct | Exact chains | Errors | Molecule accuracy | Molecules majority exact |
  |---|---:|---:|---:|---:|---:|
  | Frozen, label-free (A) | 305/528 | 10 | 221 | 0.740 | 3 |
  | Promoted method, transductive | 361/528 | 17 | 165 | 0.771 | 6 |

  - Per molecule: 2 better, 6 worse, 8 unchanged. The total is −56, with an
    interval of [−139, −3]. Track 117 (36 scans) gives −36 and track 109
    (13 scans) gives −13.
  - The campaign-B pool model (uncapped, learned weights, Fisher from the
    run) gave 349 on the same scans.
- **Invariance** (2024 model applied in both runs).
  - Of the 146 shared scans, 142 have identical inputs. Their per-scan
    frozen calls are identical in 142 of 142.
  - Fused calls are identical in 133, because tracks can differ between
    runs.
  - The transductive calls of the two runs agree on only 41 of the 142
    scans. This is the cohort dependence of the promoted method, now
    measured directly.
- **Decision** (rule fixed in the declaration). Invariance holds, and on
  held-out 2024 sessions A has more correct units (1979 vs 1847) and a higher
  molecule accuracy (0.736 vs 0.698). A is therefore a candidate to replace
  the transductive assignment, but it is not promoted:
  - The transfer to a later campaign is worse than transductive learning on
    that campaign itself (305 vs 361, 6 molecules worse and 2 better). The
    10–20mer targets come from such later campaigns.
  - The cause is not identified. Candidates are the 2024 Fisher model, the
    equal GMM weights, the cap, or a campaign shift in the features. Any
    diagnosis on 2025 is now exploratory.
  - Promotion needs the user's decision and a full pipeline run.

### 2026-10-21 — Campaign A diagnosis (exploratory): the per-track cap moves the GlcNAc boundary, in opposite directions for 2024 and 2025

- **Why exploratory.** Asked by the user after the A results. The 2025
  sessions have now been read several times, so nothing here selects a
  setting.
- **Method.** `grade_frozen_assignment.jl --ablate` undoes the A choices one
  at a time:
  - `cap`: no per-track cap;
  - `weights`: learned GMM weights;
  - `fisher`: the run's transductive emp_fisher instead of a frozen Fisher.

  The report also gives the GlcNAc call rate and the errors by position on
  six-lobe scans. Outputs are in `results/frozen_ablation_20261021/`.
- **Results** (correct units; exact chains; molecule accuracy):

  | Variant | Held-out 2024 (/2430) | 2025 (/528) | GlcNAc call rate 2024 / 2025 |
  |---|---|---|---|
  | A | 1979; 188; 0.736 | 305; 10; 0.740 | 0.316 / 0.376 |
  | A without the cap | 1892; 131; 0.718 | 359; 16; 0.750 | 0.388 / 0.462 |
  | A with learned weights | 1976; 185; 0.735 | 305; 10; 0.740 | 0.312 / 0.376 |
  | A with the run's Fisher | 1978; 185; 0.729 | 301; 10; 0.719 | 0.322 / 0.379 |
  | All three undone (= campaign-B pool) | 1871; 118; 0.712 | 349; 10; 0.729 | 0.359 / 0.477 |
  | Promoted method, transductive | 1847; 104; 0.698 | 361; 17; 0.771 | 0.428 / 0.458 |

  Undoing all three reproduces the campaign-B pool exactly (1871, 349),
  which checks the code paths.
- **Reading.**
  - **Only the cap matters.** The weights and the Fisher change almost
    nothing.
  - **The cap lowers the GlcNAc call rate.** The heavily scanned tracks
    (for example track 58, 39 scans) no longer dominate the clusters.
    - On held-out 2024 this removes false GlcNAc in the middle of the chain
      (171 → 74 at positions 3 and 4).
    - On 2025 it misses true GlcNAc (38 → 88 at positions 2 and 5). There,
      the transductive method and the uncapped model both call about 46%
      GlcNAc.
  - **Per molecule, uncapped minus A.**
    - 2024: 13 better, 21 worse; −87 units, interval [−198, −11].
    - 2025: 6 better, 4 worse; +54 units, interval [+1, +137].
  - **The large swings each come from one long track.** One physical lobe
    is fused differently across all its scans: track 58 accounts for −43,
    tracks 117 and 109 for +36 and +13.
- **Conclusion.**
  - The best GlcNAc decision boundary differs between the 2024 and 2025
    campaigns. No single frozen reference is best for both, and the cap
    cannot be chosen from these data.
  - The transductive method adapts its boundary to each run without labels,
    at the price of cohort dependence.
  - For the 10–20mer campaigns this argues for a reference taken from the
    application's own campaign (frozen per campaign), or for a label-free
    per-campaign calibration of the boundary.

### 2026-10-21 — 10–20mer triage: catalog and first visual pass

- **Deliverable** (roadmap step 2, approved by the user: "go"). The usable
  Cu(100) 10–20mer scans, grouped by molecule and decided by the user,
  recorded in a tracked table. That set is the input of the application run
  (step 3).
- **Catalog** (`catalog_sxm_candidates.jl` on the 33 Cu(100) sessions of
  `data/raw/10_20mer`, outputs in
  `results/data_review/10_20mer_catalog_20261021/`). Of 1553 scans:
  - 898 usable;
  - 582 overviews (frame larger than 20 nm);
  - 70 with the object crossing the frame;
  - 2 unreadable;
  - 1 byte-identical copy.

  The Ag(111) sessions are left out, because the method is calibrated on
  Cu(100).
- **Table** `benchmarks/chitosan_10_20mer_cu100_triage.tsv`.
  - One row per scan except the 25 targets (own table) and the copy, so
    1527 rows: 874 in set `10_20mer_candidates`, the rest discarded by the
    catalog status.
  - Two different scans named `260130_001.sxm` (sessions of 29 and
    30 January) take their session as prefix.
  - `data/sets/10_20mer_candidates/` and `data/INDEX.tsv` were rebuilt.
- **Review page** (`report_candidate_review.jl --profile 10_20mer`).
  - The page now has profiles: candidate and reference sets, decision word
    (`garder`), motifs, texts, and a 16 nm window instead of 8 nm.
  - With the 6-mer profile the rendered template is byte-identical to
    before.
  - The page is at `results/data_review/10_20mer_candidates_review/index.html`:
    372 molecule groups, with one target per session as reference.
- **First visual pass (assistant).** It suggests 279 scans to keep and 595 to
  discard:
  - 234 short fragments;
  - 95 several objects;
  - 94 folded, looped or crossed chains;
  - 51 images with no molecule;
  - 44 objects that do not look like chitosan;
  - 44 uncertain;
  - 25 chains cut by the frame or the end of the scan;
  - 8 tip artifacts.

  Suggestions are in `results/data_review/suggestions_candidats_10_20mer.tsv`.
  - **Length criterion.** The 6-mer control molecules measure about 3 nm on
    these tiles, so a straight chain shorter than about 4 nm has the length
    of a 6–7-mer, not of a 10–20mer, and is suggested as `fragment_court`.
    Several long series (for example `260129_003`–`025`, 23 scans of one
    3 nm molecule) fall here.
  - **Kept chains with features to check.** These carry a note: chains
    crossing a step, chains with a bright end blob (possible adsorbate,
    `260218_021`–`038`), a faint parallel ghost (double tip), and a few
    long irregular chains.
- **English page for colleagues** (same day, at the user's request).
  - The `10_20mer` profile is now in English: interface, motifs
    (`not_chitosan`, `cut_chain`, `folded_or_crossed`, `several_objects`,
    `short_fragment`, `tip_artifact`, `no_molecule`, `uncertain`), decisions
    `keep`/`discard`, tile titles, and a guide to the 4 nm length criterion.
  - A "Your name" field is required for export and written in a `reviewer`
    column, so the files of several reviewers can be merged.
  - The first-pass notes were translated into
    `results/data_review/suggestions_candidates_10_20mer_en.tsv`.
  - The page is in `results/data_review/10_20mer_candidates_review_en/`.
    pngquant compression takes it from 203 to 55 MB; a README for the
    reviewers is included, and `10_20mer_candidates_review_en.zip` is ready
    to share.
  - The 6-mer page still renders byte-identically to the original.
  - The images are unpublished raw data, shared only within the group.
- **Next.** The reviewers export their decisions. These are merged (the
  English motifs map to the table's French codes; disagreements are
  listed) and applied to the table (set `10_20mer_application` for kept
  scans), then the sets are rebuilt. Only after that comes step 3, the
  application run.

### 2026-10-21 — 10–20mer test run and sequence-pattern analysis

- **Request.** "With the selection already made, launch the 10–20mer test;
  try to identify a pattern."
- **Input** (440 Cu(100) scans, staged on Viper as
  `/ptmp/<user>/stmfit/app1020_inputs/raw`):
  - the 279 scans of the assistant's first pass;
  - the 25 targets;
  - 136 whole short chains (7–8 units) that the first pass had set aside
    only for their length.

  The manual readings in `data/raw/notes/SequenceAnalysis.xlsx` contain
  7–8-unit chains (for example `251210_025`, read `0101111`), so the 4 nm
  length rule was too strict. The short chains are reported as their own
  subset. The manual readings are now
  `benchmarks/chitosan_10_20mer_manual_sequences.tsv` (grading side only).
- **Run.**
  - Viper job 12310547 (`5b706ca`), promoted method with the long-chain
    counting config. Counting ran as 96 single-thread shards (new driver
    option `--count-shards`, results independent of the shard count), so
    the whole run took 1 h 34.
  - Outputs are in `results/app1020_5b706ca/`; `grade_consensus_run.jl
    --verify-only` gives ALL PASS.
  - Final count rules: 131 agree, 40 consensus applied, 4 registered ROI,
    88 tracks too short, and 177 without a strict majority. Repeat scans
    of one molecule often disagree on the count.
  - Fusion: π 0.385, θ0 0.191, θ1 0.619 (previous run θ0 0.45; 6-mers
    0.05).
- **Pattern report** (`scripts/grading/report_sequence_patterns.jl`, outputs
  in `results/app1020_5b706ca_patterns/`).
  - There are 133 molecules (tracks), with a median N of 8 (range 2–18);
    99 have N ≥ 7.
  - Three assignments are compared: the run's transductive calls, a frozen
    label-free model per campaign (5 campaigns; the last, with 3 scans, has
    no model, and the second has only 5 reference scans), and the frozen
    6-mer model of campaign A.
  - Every statistic is tested against 10,000 shuffles within each chain
    (composition kept).
- **Results.**
  - **Composition.** GlcNAc (1) fractions are 0.35 (transductive), 0.40
    (frozen per campaign) and 0.36 (6-mer model).
  - **Ends.** All three call the terminal units GlcN more often than chance
    (ends 0.29 vs 0.36 for the transductive calls, p = 0.01; 0.16 vs 0.35
    for the 6-mer model, p < 0.0001). The interior is enriched in GlcNAc.
  - **Alternation.** There is more alternation near the ends than chance
    (longest alternating end segment 2.78 vs 2.53, p = 0.04; 6-mer model
    3.31 vs 2.56, p < 0.0001). With the frozen models the other end
    alternates too (p < 0.0001), and the share of chains alternating at one
    end only is not above chance (p 0.05–0.5).
  - **Manual template.** Agreement with the manual template "1010" + "1"s,
    or with its 0/1 inverse, is at or barely above chance (transductive,
    inverse: 0.660 vs 0.650, p = 0.23).
- **Manual readings** (15 scans read by hand, all in the run).
  - **Length.** The final N equals the manual length for 6 scans and is
    within one for 13. `251206_013` stays at 4 (manual 12).
  - **Convention.** Per-unit agreement is 43–49% with the readings as
    written and 63–70% with 0/1 inverted. (Corrected on 2026-10-22: the
    readings encode ring orientation, 1 = bright "3H-up", not GlcNAc/GlcN;
    see the provenance entry.)
  - **Order.** In both conventions the agreement is at chance level against
    a composition-preserving shuffle of the calls. For example, the
    transductive calls with the inverted convention give 83 observed vs 83
    expected (p = 0.55); the frozen models reach p 0.05–0.09. The agreement
    comes from composition, not from unit order.
- **Reading.**
  - The method does not recover the manual sequence order on these scans.
  - The only robust structure it finds is the 6-mer edge signature:
    terminal units called GlcN, a relatively bright neighbour, alternation
    near both ends. That is the known edge-adjacency confound, not a
    one-sided chemical motif.
  - The repeat scans disagree often, both on counts (177 without a
    majority) and on calls (θ0 0.19). Long-chain counting and per-scan
    assignment consistency (open questions 2 and 3) remain the limits.
- **Not done.** No parameter was changed; this is a test of the promoted
  method plus label-free frozen variants. The meaning of 0/1 in the manual
  sheet is a question for the user.

### 2026-10-22 — 10–20mer provenance check: authentic files, unverified sample

- **Trigger.** After a partial visual triage, the user doubts that many of
  the imaged molecules are chitosan, and questions the provenance of the
  files.
- **User's partial triage** (`/home/durif/tmp/triage_decisions_10_20mer.tsv`,
  reviewer olivier; 362 of 874 cards decided).
  - 67 kept and 295 discarded. The discards include 72 "not chitosan" and
    54 "uncertain", the latter often noted "what is this molecule?".
  - The kept scans come from 6 and 10 December 2025 and 16 and 29 January
    2026. Of the February sessions only 2 scans are kept. The 25 March
    session has 25 "not chitosan".
  - The user agrees with 311 of the assistant's first-pass suggestions and
    disagrees with 51.
- **Chain of custody.**
  - `data/raw/10_20mer` is identical to the institute share folder
    `STM_DATA/Chitosan_10_20mer`: 6503 files, same names and sizes.
  - Every SXM header names that folder as the acquisition location
    (`SCAN_FILE` and Nanonis session path
    `Z:\DATA\STM\STM_DATA\Chitosan_10_20mer\<session>`).
  - The imaging conditions match the (NKN)2 sessions: −0.3 V, 1 pA, same
    instrument.
  - There is no local mix-up: the files are what the instrument saved under
    that project name.
- **Sample identity is not documented in the data.**
  - The header comments are empty in all 1637 SXM files, and the session
    `.ini` files carry no sample note.
  - Only the folder name says "chitosan 10–20mer".
  - The share also holds other chitosan projects (`Chitosan`, `Mucor
    chitosan 1`, `Random_chitosan_polymer_Fabian`) and many sugar or
    glycan projects.
- **The "manual sequences" are not acetylation labels** (correction of the
  2026-10-21 entry).
  - `SequenceAnalysis.xlsx` and the 25 targets come from a third-party
    analysis, `data/raw/notes/_KA_SOLVE.pptx`.
  - Its key reads: bright blob = "3H-up", dimmer blob = "2H-up", which are
    ring orientation states. Units are labelled "Glc". It notes a
    room-temperature deposition with fragmented dimers next to oligomers.
  - So 1/0 there encode an orientation, not GlcNAc/GlcN. The 2026-10-21
    comparison and the "inverted convention" were not meaningful; the
    conclusion that unit order agrees only at chance level stands
    regardless.
  - The table header of `benchmarks/chitosan_10_20mer_manual_sequences.tsv`
    now says so.
- **Physical signature** (median lobe spacing per scan from the
  2026-10-21 run, label-free).
  - 6-mer benchmark: 0.597 nm.
  - Closest: 10 December (0.598) and 29–31 January (0.58–0.61).
  - Larger in most February sessions: 0.63–0.71 nm (18 February 0.683,
    19 February 0.705).
  - The sessions the user keeps are among the closest to the 6-mer
    reference. Much of the February material is physically different
    (another molecule, or another adsorption state).
- **Consequence.**
  - The 10–20mer application is on hold until the deposited sample is
    identified: compound, supplier, DP and DA, deposition method and date,
    from the lab notebook or the experimenter.
  - A GlcN/GlcNAc assignment is meaningless on molecules that are not
    chitosan.
  - The triage table is unchanged, because the user's review is
    incomplete.

### 2026-10-22 — Notation correction: N is GlcNAc in NKNNKN; the classes lose their chemical names

- **Error.** In the control NKNNKN, N = GlcNAc (N-acetyl) and K = GlcN
  (confirmed with the user). The control has four GlcNAc (positions 1, 3, 4,
  6) and two GlcN (positions 2, 5). The truth header
  (`benchmarks/chitosan_240817_unit_sequences.tsv`) defined `010010` with
  1 = GlcNAc, which put the GlcNAc at 2 and 5, and every document repeated it.
- **Unchanged.** The 0/1 strings, every count and every grade. The graders
  compare the method's class codes with the strings, whatever their chemical
  names, so campaigns B, option 1 and A stand as reported.
- **Changed meaning.**
  - The method's class 1 is the brighter cluster (higher mean raw
    amplitude). On the control it falls on the two GlcN. The rule "brighter =
    GlcNAc (the acetyl group is larger)" is contradicted by the benchmark,
    not supported by it. "Correct" in earlier entries means "the bright class
    falls on the GlcN".
  - If the contrast is chemical, the bright units are GlcN; if it is
    positional, they are the neighbours of the chain ends. The control cannot
    tell the two apart (open question 1).
  - The per-feature AUCs of 2026-09-30 used the old names: amplitude is
    higher at the GlcN than at the interior GlcNAc (AUC ≤ 0.66).
  - On the 10–20mers the "edge signature" means terminal units in the dimmer
    class and a bright class enriched in the interior. It says nothing about
    acetylation.
- **DFT molds** (exploratory check, 391 six-lobe scans of
  `results/bench493_a2bcb15`).
  - The molds are defined chemically (GlcN and GlcNAc trimer cubes), so the
    correction does not touch them. The method uses the absolute mold margin.
  - The signed margin, cost(GlcN) − cost(GlcNAc), favours the GlcNAc mold
    slightly more at the interior GlcNAc (3, 4) than at the GlcN (2, 5):
    AUC 0.56 (forward) and 0.57 (backward).
  - The decoded GlcNAc fraction is 0.53–0.60 at the GlcNAc positions and
    0.48–0.53 at the GlcN.
  - That is the chemical direction, but weak, and scans of one molecule are
    correlated. Signed CC margins were rejected on 2026-09-21 when graded
    against the inverted names. Reopening them is a new declared experiment,
    not done here.
- **Changes.** No method behaviour changed and numeric outputs are
  identical.
  - Truth header corrected; strings unchanged.
  - Method code, plots and graders call the classes "class 1 (brighter)" and
    "class 0 (dimmer)". Plot legends read "dimmer (0)" and "brighter (1)".
  - The pattern report columns become `class1_fraction`, `ends_class1`,
    `interior_class1` and, in `profile.tsv`, `class1`.
  - The call profile of `grade_frozen_assignment.jl` names the positions by
    their chemistry.
  - The comments of the two frozen assignment configs now say "no class
    fraction assumed". The edit is comment-only, but their sha256 changes.
  - AGENTS.md (the naming invariant is now physical and chemically neutral),
    README, index, unit assignment, pipeline and calibration pages.
- **10–20mer targets.** The 25 targets are the scans of the former STMFit
  `10_20mer_analysis` run. A copy of those outputs (counts and plots) sits
  in the experimenters' share, so it is not an independent reading.

### 2026-10-22 — 10–20mer: comparison with the experimenters' acetylation readings

- **Source.** At the user's request, the experimenters' 10–20mer analysis
  folder was copied to `data/raw/notes/10_20mer_analysis/` (untracked,
  unpublished).
  - `final_results.md`/`.tsv` and the 25 `*_best.png` are outputs of the
    former STMFit `10_20mer_analysis` run (the 25 targets), not an
    independent reading.
  - The deck `Corrected 10-20mer.pptx` holds manual acetylation readings of
    18 molecules (DP, DA and sequence). Its DA is the fraction of 1s and the
    molar mass from DP and DA matches the mass spectrum, so 1 = GlcNAc.
    These are the first acetylation labels for the 10–20mers.
  - They are in the local table `acetylation_readings.tsv` (grading side
    only). One scan holds several molecules and is not compared, and one is
    not in the run, which leaves 16.
- **Report change.** `report_sequence_patterns.jl` now tests the manual
  agreement against a null that shuffles each scan's calls (composition
  kept). The 2026-10-21 statement "at chance once composition is accounted
  for" is now reproducible; its tables are unchanged.
- **Results** (run `app1020_5b706ca`, 10,000 shuffles; output in
  `data/raw/notes/10_20mer_analysis/comparison_app1020_5b706ca/`).
  - **Counting.** The final N equals the read length in 4 of 16 scans and is
    within one in 14.
  - **Bright class read as GlcN.** 64–67% of units agree, but the null gives
    76–85 of 128 (p 0.08–0.68).
  - **Bright class read as GlcNAc.** 45–49%, below the null (p ≥ 0.88).
  - No assignment (per-scan, transductive, frozen per campaign, frozen 6-mer
    model) recovers the acetylation order. The apparent preference for
    "bright = GlcN" comes from composition: class 1 is the minority class,
    like GlcN in the readings.
  - **Ends.** The readings show no end effect: the GlcNAc fraction at the
    chain ends equals the interior one. The calls put the dim class at the
    ends (class-1 fraction 0.22 at the ends vs 0.30 inside for the
    transductive calls, 0.09 vs 0.45 for the 6-mer model). The edge
    signature is positional, not chemical.
- **Reading.** On these 16 molecules the method gives counts close to the
  manual ones but no acetylation information beyond composition. Better unit
  calls need new information, not another classifier variant (open
  question 7).

### 2026-10-22 — 10–20mer working set: the experimenters' selection (19 scans)

- **Why.** The user's partial review finds the bulk 10–20mer candidates
  mostly unusable (not chitosan, fragments, several objects). The user asked
  for the experimenters' exploitable list instead, even if much shorter.
- **Sources** (searched read-only: the group chat and the experimenters'
  share). No complete list is published. Three partial selections exist:
  - their acetylation readings (corrected slide deck): 18 molecules in 17
    scans, one scan with several molecules, one broken molecule;
  - their January 2026 slides on the December sessions: 8 scans, 4 with
    acetyls read and 4 with lengths only;
  - the third-party list of `benchmarks/chitosan_10_20mer_manual_sequences.tsv`
    (19 scans, 15 read and 4 rejected), which is also the experimenters'
    folder of 19 scans.
  - The chat mentions larger January and February selections whose lists
    are not available.
- **Union: 23 scans**, in `benchmarks/chitosan_10_20mer_curated.tsv`.
  - 19 are usable (`10_20mer_curated`): 17 with an acetylation reading (two
    possibly fragments), 1 with a length only (`251210_029`), and 1 read
    only by the third party (`251206_013`).
  - 4 are set aside: `251206_019` and `260222_010` (several molecules in the
    frame, while the method analyses one chain per scan), `251212_011` (cut
    chain, user's triage) and `260221_044` (broken molecule, per the
    experimenters).
- **Former targets.** 14 of the 25 targets are in the selection. The other
  11 are in no experimenter list (they came from the former STMFit analysis
  folder). They stay in `chitosan_10_20mer_targets.tsv` (set
  `10_20mer_targets`) and are not part of the working set.
- **Tables.** The 23 rows left the targets table (11 rows remain) and the
  triage table (9 candidates left, 865 remain). `data/sets/10_20mer_curated`
  holds 19 links, and `INDEX.tsv` was rebuilt. The review page now takes its
  reference tiles from the curated set; the page itself was not regenerated.
- **Coverage.** 18 of the 19 usable scans were in run `app1020_5b706ca`
  (`260221_056` was not); 16 of them carry acetylation readings (previous
  entry).
- **Next.** Ask the experimenters for their complete January and February
  selections, then run the method on the curated set only.

### 2026-10-25 — Declared campaign: capped support for counting, legacy support for geometry

- **Deliverable** (approved by the user on 2026-10-25, time box until
  2026-10-27). Long-chain counts closer to the experimenters' length
  readings, with no regression on the 6-mer benchmark.
- **Why** (open question 2). The capped support (`half_maximum_cap`)
  removes truncated long-chain supports: mean |N − N_manual| went from 1.33
  to 0.92 on 12 reads (2026-09-30). But it collapses the 6-mer assignment
  (771/88 → 702/24), because the assignment depends on end-lobe geometry.
  So counting and assignment geometry are separated.
- **Change** (one new config key, label-free).
  - `model.geometry_support_threshold_rule` sets the support rule of every
    fixed-N refit after counting: scan geometry, consensus refits and
    assignment features. When absent, it equals `support_threshold_rule`
    and nothing changes.
  - Counting (`batch_full.jl`) selects N on the capped support.
  - Each scan is then refit at that N on the legacy support. If that fit is
    missing or flagged, the scan is refit on the counting support instead.
  - The choice is recorded per scan (column `support_threshold_rule` of
    `geometry_summary.tsv`, carried into the consensus summary). Every later
    refit replays it. Only fit success decides: no label, count or reading.
  - New configs `config/chitosan_decoupled_support.toml` and
    `config/chitosan_10_20mer_decoupled_support.toml` equal the two counting
    configs except `support_threshold_rule = "half_maximum_cap"` and
    `geometry_support_threshold_rule = "profile_quantile_noise"`.
- **Runs** (Viper, one each, no retuning after grading).
  - 6-mer: the 146-scan cohort with `chitosan_decoupled_support.toml`.
    Reference: run A (12042058), reproduced byte for byte at `b8986a0`.
  - 10–20mer: the 440 scans of the test run with
    `chitosan_10_20mer_decoupled_support.toml`. Reference:
    `app1020_5b706ca`. The code path without the new key is unchanged, so
    the reference stands.
- **Criteria, fixed before any run.**
  - 6-mer (external grade, full145 own-N): no regression against A on exact
    N (137), within one (145), correct (771), exact chains (88), classified
    (865) or errors (94).
  - 10–20mer (grading side only): the 17 lengths the experimenters read on
    scans of this cohort (16 acetylation sequences, plus the unit count of
    `251210_029`). The reference run gives mean |N_final − N_read| = 1.47,
    14/17 within one unit and 4/17 exact. Accepted if the mean decreases and
    the within-one count does not.
  - Reported, not criteria:
    - per-scan counts on the same readings;
    - `251206_013` (third-party length reading);
    - the label-free count agreement between repeat scans (reference: 177
      of 440 scans in a track without a strict majority);
    - the number of scans that needed the counting-support fallback.
- **Adoption.** If both criteria hold, the decoupled configs replace the
  promoted counting configs and every cited headline is updated. Otherwise
  the campaign is recorded as rejected.

### 2026-10-25 — Decoupled supports: long-chain counts improve, 6-mer exact chains drop; rejected by the declared rule

- **Runs** (code `70a9baa`, all label-free checks pass, outputs in
  `results/dec6mer_70a9baa/` and `results/dec1020_70a9baa/`):
  - 6-mer, 146 scans, Viper job 12380719 (38 min, 96 counting shards);
  - 10–20mer, the 440 scans of the test run, job 12380720 (1 h 59).
- **Fallback.** Geometry needed the counting support for 1 of 146 6-mer
  scans (`240818_019`, whose legacy support cannot hold its count) and for
  59 of 440 10–20mer scans.
- **6-mer** (external grade, full145 own-N, against run A):

  | Run | Exact N | Within one | Per-scan exact N | Correct | Exact chains | Classified | Errors |
  |---|---:|---:|---:|---:|---:|---:|---:|
  | A (legacy support) | 137 | 145 | 123 | 771 | 88 | 865 | 94 |
  | Decoupled supports | 137 | 145 | 126 | 776 | 85 | 865 | 89 |

  - Every final count equals run A's, on all 146 scans. Eight per-scan
    counts change, and so do their geometries.
  - Those eight scans change the cohort-learned clusters. Final calls change
    on 17 scans in 6 molecules: one 10-scan molecule gains a unit, and three
    chains lose their exactness. Fusion π is 0.330 (A: 0.342; truth 0.333).
  - Exact chains drop from 88 to 85, so the 6-mer criterion fails.
- **10–20mer** (`scripts/grading/report_length_readings.jl`, readings
  untracked):

  | Run | Readings | Exact (final) | Within one | Mean abs ΔN (final) | Mean abs ΔN (per scan) |
  |---|---:|---:|---:|---:|---:|
  | Reference `app1020_5b706ca` | 17 | 4 | 14 | 1.47 | 1.53 |
  | Decoupled supports | 17 | 5 | 14 | 1.00 | 1.06 |

  - The large truncations are gone: `251206_014` goes from 7 to 11 units
    (12 read) and `260220_083` from 2 to 9 (8 read). `251206_013`, read by a
    third party only, goes from 4 to 11 (12 read).
  - Three chains move by one unit: `251210_018` becomes exact, while
    `251206_034` and `251206_038` move one unit away from the reading.
  - The criterion holds: the mean decreases and the within-one count does
    not drop.
  - Repeat scans agree more often: 164 scans sit in a track whose counts
    agree (reference 131), and 145 in a track without a strict majority
    (reference 177).
  - 96 final counts change (70 up, 26 down). Fusion: π 0.48, θ0 0.19,
    θ1 0.71 (reference 0.38 / 0.19 / 0.62).
- **Decision.** By the declared rule, both criteria had to hold. The 6-mer
  one fails on exact chains, so the decoupled configs are not adopted and
  no headline changes. They stay in `config/` as the record of this
  experiment.
- **Reading.** Decoupling achieves its counting goal: long chains are
  counted on full supports, and no 6-mer count is lost. The 6-mer
  assignment no longer moves through geometry, since the scans with
  unchanged N keep run A's geometry. It moves through the cohort-wide
  learning, which a few changed per-scan counts are enough to shift
  (open question 8). Using the decoupled counting only for the 10–20mer
  application config would leave the benchmark untouched. That would be a
  new decision for the user, not part of this campaign's rule.

### 2026-10-25 — Decision: decoupled supports promoted for the 6-mer method and the 10–20mer application

- **Decision** (the user, after the grade): "on peut adopter les supports
  découplés pour 10–20mer ET 6-mer puisque l'ensemble du benchmark est
  équivalent (pour le 6-mer, même si le nombre de chaînes exactes recule, le
  nombre d'erreurs recule aussi)".
- **Against the declared rule.** The campaign rule required no regression
  on any 6-mer figure, and exact chains fell from 88 to 85. This adoption
  overrides that rule after seeing the grade. No parameter was retuned: the
  promoted configs are exactly the declared ones.
- **Grounds.**
  - 6-mer: same final counts (137 exact, 145 within one, on all 146
    scans); correct units 771 → 776; errors 94 → 89; classified unchanged;
    exact chains 88 → 85. Exact chains come in molecule blocks, and here 43
    of 85 come from one 43-scan molecule.
  - 10–20mer: no truncated supports any more, and repeat scans agree more
    often on the count. The counts are closer to 17 experimenter length
    readings (mean |ΔN| 1.47 → 1.00), which are a comparison, not a truth.
- **What changes.**
  - The promoted counting configs are `config/chitosan_decoupled_support.toml`
    (6-mer method) and `config/chitosan_10_20mer_decoupled_support.toml`
    (10–20mer application).
  - `hpc/run_molecule_consensus.sbatch` defaults to the first.
  - `config/chitosan.toml` stays the base config and the default of the
    single-stage scripts `batch_full.jl` and `extract_lobe_features.jl`.
  - Updated headlines: README, docs index, unit-assignment page, runbook,
    selection, config and calibration pages.
- **Not redone.** The 10–20mer pattern report and the comparison with the
  experimenters' acetylation readings are still those of the legacy run
  `app1020_5b706ca`. The 493-scan extended benchmark and the campaigns A and
  B used the legacy support.

### 2026-10-25 — 10–20mer lengths are not a ground truth; N/K comparison on the promoted run

- **User statement.** The experimenters' 10–20mer chain lengths are not
  trusted as results: they are an interesting comparison, never a truth
  criterion. Unlike the 6-mer, whose length is known, the 10–20mer length is
  unknown. The true N/K pattern of the 10–20mer sequence will be provided
  later by the user, and it will be the reference.
- **Consequence.** The length criterion of the decoupled-support campaign
  (mean |ΔN| 1.47 → 1.00 against 17 readings) is a comparison with human
  readings, not evidence of correct counting. The robust part is
  label-free: the truncated supports disappear, and repeat scans agree more
  often on the count (145 scans in a track without a strict majority,
  against 177). No 10–20mer count is validated.
- **N/K assignment against the experimenters' acetylation readings**
  (promoted run `dec1020_70a9baa`, 16 read molecules, 138 units;
  comparison only, readings untracked):

  | Calls | Brighter read as GlcNAc | p | Brighter read as GlcN | Null mean | p |
  |---|---:|---:|---:|---:|---:|
  | per scan | 71 (51%) | 0.96 | 84 (61%) | 84.2 | 0.55 |
  | transductive (promoted, fused) | 73 (53%) | 0.92 | 84 (61%) | 83.8 | 0.53 |
  | frozen per campaign | 65 (47%) | 0.80 | 94 (68%) | 89.1 | 0.14 |
  | frozen 6-mer model | 67 (49%) | 0.97 | 87 (63%) | 86.7 | 0.52 |

  - The null shuffles each scan's calls, keeping its composition.
  - As on the legacy run (61–67% agreement, p 0.08–0.68), the agreement
    comes from composition, not from unit order.
  - The promoted run calls more units brighter: class-1 fraction 0.45,
    against 0.35 on the legacy run. Fusion π is 0.48.
  - The edge signature remains: terminal units are dimmer than chance
    (transductive 0.40 vs 0.45, p = 0.04; frozen 6-mer model 0.15 vs 0.35,
    p < 0.0001). The readings show no such end effect.

### 2026-10-25 — 10–20mer sequences of the promoted run: isolated bright units, alternating ends (exploratory)

- **Data.** Run `dec1020_70a9baa`: 133 molecules (consensus tracks), 65
  imaged at least 3 times, 34 with N ≥ 10. One sequence per molecule (the
  majority of its scans' final calls per lobe); 1 = brighter class.
  Orientation and true length are unknown.
- **Structure against a within-chain shuffle** (composition kept; the
  pattern report, plus orientation-free 3-mer counts on the 95–97 molecules
  with N ≥ 7, 2000 shuffles):
  - Isolated bright units are enriched: `010` 133 vs 105 expected
    (transductive, p < 0.001), 210 vs 115 (frozen 6-mer model, p < 0.001).
    `101` is also enriched. Adjacent bright units are depleted with the
    frozen 6-mer model (`011` 77 vs 114, `111` 13 vs 25).
  - Ends: the terminal unit is dimmer than chance, and the profile oriented
    by the longer alternating end starts 0-1-0 (class-1 rate 0.36 / 0.61 /
    0.55 / 0.32 at positions 1–4, then about 0.45). Chains alternate at one
    end only more often than chance (0.44 vs 0.25, p < 0.001).
  - Class-1 fraction: 0.45 (transductive), 0.35 (frozen 6-mer model).
- **Caveats.** These are image classes, not N/K.
  - Three of the four local features compare a lobe with its neighbours, so
    adjacent calls are anti-correlated by construction. That alone can
    produce isolated bright units and alternation.
  - The frozen 6-mer model learned NKNNKN, whose brighter units are
    isolated and next to the ends.
  - Repeat scans disagree (θ0 0.19), and the experimenters' readings are
    matched no better than chance.
  - The pattern is not interpreted chemically. The user's N/K pattern for
    the 10–20mer will test it.

### 2026-10-25 — Method development strategy (user decision)

- **Optimize on the 6-mer, apply to the 10–20mer.** The 6-mer benchmark,
  whose length and sequence are known, is where the method is optimized and
  selected. The method is then applied unchanged to the 10–20mer, the
  unknown, with no a priori constraint on the result. The user's 10–20mer
  N/K pattern, once given, is a final test only.
- **Transfer.** Optimized parts are frozen before the application,
  assignment included where possible, because transductive re-learning on
  the 10–20mer cohort is not the same method. Optimization uses held-out
  sessions or molecules and block-robust measures (correct units, errors,
  per-molecule scores). In NKNNKN both K are next to a chain end, and no K
  has two interior neighbours. So the benchmark cannot validate K deep
  inside a chain. (Wording corrected the same day: positions 2 and 5 are
  interior units, but each borders a terminal unit.)
- **Dataset.** The optimization set must have impeccable image quality.
  More images is not better when they are unusable: the 493-scan extension
  lowered both counting and assignment (2026-10-19). Scans are selected on
  image quality judged before any method output, and the 10–20mer data are
  held to the same standard.
- Recorded in `AGENTS.md` ("Method development strategy").

### 2026-10-25 — Declared campaign C: optimize the assignment to freeze on the validated 6-mer set, apply it to the curated 10–20mer

- **Deliverable** (user request, 2026-10-25: optimize further on the 6-mer,
  then give the results and sequences of the selected 10–20mer set). Time
  box: 2026-10-25 to 2026-10-26.
  - Strategy: `AGENTS.md`, "Method development strategy". Optimize on the
    6-mer, freeze, apply unchanged to the 10–20mer.
  - Optimization set: the 493 user-validated 6-mer scans
    (`data/sets/6mer_Cu100_benchmark`).
- **Step 1, measurement.** The promoted method (decoupled supports) on the
  493 scans, one Viper run. Graded by cohort (full145, added_2024,
  added_2025) and compared with the legacy-support 493 run (12282486). No
  criterion: this is the reference for steps 2 and 3.
- **Step 2, optimization: which assignment to freeze.** Candidates, all
  label-free at inference:
  - T, the promoted transductive assignment (learned on each run's own
    scans);
  - F(c), the frozen label-free model of campaign A
    (`unit_assignment_frozen_labelfree.toml`), with the reference cap
    `reference_max_scans_per_track` c ∈ {1, 3, 5, 10, none}. Only this key
    varies.

  Each candidate is graded on held-out units only, with
  `grade_frozen_assignment.jl` on the step-1 run:
  - every 2024 session is graded with a model built without it;
  - 2025 is graded with the model built on all 2024 sessions. The 2025
    sessions were already read by campaigns A and B, so they are not
    blind.

  **Rule, fixed now.** The selected candidate has the most correct units
  over all held-out units (2024 held-out plus 2025), with fewer errors as
  the tie-break. Exact chains and paired per-molecule differences are
  reported, not used.
- **Step 3, application.** The promoted counting and geometry
  (`chitosan_10_20mer_decoupled_support.toml`) run on the 19 curated
  10–20mer scans (`data/sets/10_20mer_curated`), one Viper run.
  - If an F(c) wins, its model is rebuilt on all 493 6-mer scans with the
    selected cap and applied unchanged (`apply_frozen_assignment.jl`). If T
    wins, the run's own calls are used.
  - Output: one sequence per scan. The classes are physical (1 = brighter);
    they are shown as N/K under both conventions.
  - Nothing is constrained by length, composition or pattern. The
    experimenters' readings are shown as a comparison only.

### 2026-10-25 — Declared campaign Q: automatic, label-free image-quality decision for the 6-mer optimization set

- **Deliverable** (user request: a blind visual review is hard to make
  objective, so the decision should be automatic). A reproducible per-scan
  quality decision, and with it the 6-mer optimization set. Time box:
  2026-10-25 to 2026-10-26.
- **Measurements.** Image only: no fit, no count, no label, no method
  output. They come from `scripts/data/catalog_sxm_candidates.jl`, which
  now also writes continuous values:
  - technical status: readable, frame at most 20 nm, object not crossing
    the frame, not a duplicate;
  - completeness: observed pixel fraction, and whether the molecule region
    (ROI) touches the border or the unscanned area;
  - isolation: area of the second-largest object above the ROI threshold,
    divided by the largest (`second_object_ratio`);
  - tip stability: forward/backward correlation and NRMSE, line
    discontinuity, stripe periodicity;
  - signal: contrast to noise (99th percentile of the ROI signal minus the
    off-ROI median, divided by the image noise);
  - size: axial length of the ROI.
- **Rule, fixed now** (`config/scan_quality.toml`). A scan is accepted
  when all of the following hold:
  1. its technical status is usable;
  2. the observed pixel fraction is at least 0.99, and the ROI touches
     neither the border nor the unscanned area;
  3. `second_object_ratio` < 0.3;
  4. the forward/backward correlation is at least 0.5;
  5. line discontinuity, stripe periodicity and NRMSE are below
     median + 3 × 1.4826 × MAD;
  6. contrast to noise is above median − 3 × 1.4826 × MAD.

  The values in 2–4 are the catalog's existing flags. The bounds in 5–6
  are computed once, on the technically usable Cu(100) 6-mer scans, and
  written into the config as absolute numbers before any comparison with
  human decisions. The same numbers then apply to any data set, the
  10–20mer included.
- **Gross size check (6-mer only).** The ROI length must lie within
  0.5–1.5 × the median of the accepted scans. One unit is about 17% of a
  6-mer, so 5- or 7-unit counts stay well inside this window; it only
  removes objects of under about 3 or over about 9 units. It applies only
  to the benchmark, whose molecule is known, and never to the 10–20mer,
  whose length is unknown. Results are also reported without it.
- **Validation** (grading side). The automatic decision is compared with
  the earlier human decisions on the same scans
  (`benchmarks/chitosan_6mer_cu100_triage.tsv`, excluding the automatic
  overview rejections):
  - overall agreement;
  - the share of each human rejection reason that the rule catches;
  - the reasons for every human-kept scan the rule rejects.

  No benchmark grade is used, and the rule is not changed after this
  comparison. Any revision would be a new declared step.
- **Optimization set.** Scans kept by the human triage and accepted by the
  rule, size check included, recorded in a tracked table. Campaign C steps
  1–2 are then rerun on it before anything is applied to the 10–20mer.

### 2026-10-25 — Campaign Q, amendment before validation: molecule completeness and log-scale bounds

- **Why.** The catalog of the 1149 Cu(100) 6-mer scans
  (`results/quality_6mer_20261025/catalog/`: 827 usable, 261 overviews, 61
  objects crossing the frame) showed two problems with the declared rule.
  These are label-free distributions; no human decision or grade had been
  compared yet.
  - **Frame completeness.** 345 of the 827 usable scans stopped early
    (observed pixel fraction below 0.99, median 0.74). Their molecule is
    often whole. The criterion for the molecule is "the ROI touches the
    unscanned area" (87 scans), so the frame-level 0.99 bound is dropped.
  - **Bounds.** The values are positive and strongly right-skewed: line
    discontinuity ranges from 5 to 1e11. On the raw scale the bounds
    rejected 27% of scans for line discontinuity and were negative, so
    inactive, for contrast to noise. They are computed on log values
    instead: exp(median ± 3 × 1.4826 × MAD of the log).
- **Fixed bounds** (`config/scan_quality.toml`, derived once from the 827
  usable scans): line discontinuity ≤ 52.6, stripe periodicity ≤ 6012,
  forward/backward NRMSE ≤ 2.21, contrast to noise ≥ 1.66. The other bounds
  are unchanged (second-object ratio < 0.3, forward/backward correlation
  ≥ 0.5, size window 0.5–1.5 for the 6-mer only).

### 2026-10-25 — Campaign Q results: the automatic rule agrees with the human triage on 78.5% of scans; optimization set of 419

- **Decision** (`config/scan_quality.toml`, size check on): of the 1149
  Cu(100) 6-mer scans, 536 are accepted (551 without the size check). The
  median ROI length of accepted scans is 3.94 nm.
- **Against the human triage** (888 scans, the automatic overview
  rejections left out; `results/quality_6mer_20261025/validation/`):

  | Human | Scans | Automatic accept | Automatic reject |
  |---|---:|---:|---:|
  | kept | 493 | 419 | 74 |
  | rejected | 395 | 117 | 278 |

  - Agreement: 697/888 (78.5%), or 684 without the size check.
  - Human rejection reasons the rule catches:
    - cut molecule 101/101;
    - no molecule 35/36;
    - several objects 19/23;
    - earlier-review ambiguous 34/46 and excluded 7/7;
    - not a 6-mer: 14/21 too long, 2/4 too short;
    - uncertain 53/122;
    - tip artefact 13/35.

    Tip artefacts and uncertain scans are what the generic values miss.
  - Of the 74 human-kept scans the rule rejects:
    - reasons: line discontinuity 33, ROI at the edge 24, several objects
      19, forward/backward mismatch 10, size 1;
    - by cohort: full145 30 of 145, added_2024 42 of 260, added_2025 2 of
      88.

    The rejections are therefore not concentrated in the 2024 additions.
- **Optimization set: 419 scans**, kept by the triage and accepted by the
  rule (full145 115, added_2024 218, added_2025 86).
  - Decisions are in `benchmarks/chitosan_6mer_cu100_quality.tsv`, and the
    manifest is `benchmarks/chitosan_6mer_optimization_v1.toml` (entries of
    v2, unchanged).
  - No grade was used. Campaign C steps 1–2 now run on this set. The
    493-scan grading job was cancelled; the 493-scan run stays as a
    measurement.

### 2026-10-25 — Decoupled supports on the 493 scans: counting +26 exact, assignment slightly better (measurement)

- **Run** (code `c2b954a`, Viper job 12384144, 1 h 49; all label-free
  checks pass; `results/dec493_c2b954a/`, grades in
  `results/dec493_c2b954a_grade/`). Geometry needed the counting support for
  17 of 493 scans. Reference: the legacy-support run 12282486 (2026-10-19).
- **Counting** (exact N after consensus, within one):

  | Cohort | Legacy support | Decoupled supports |
  |---|---:|---:|
  | full145 | 137/145 | 138/145 (145) |
  | added_2024 | 177/260 (256) | 201/260 (259) |
  | added_2025 | 77/88 (88) | 78/88 (88) |
  | all | 391/493 (489) | 417/493 (492) |

- **Assignment** (final calls: correct units, exact chains):

  | Cohort | Legacy support | Decoupled supports |
  |---|---:|---:|
  | full145 | 692/870, 28 | 697/870, 34 |
  | added_2024 | 1155/1560, 76 | 1184/1560, 74 |
  | added_2025 | 361/528, 17 | 362/528, 14 |
  | all | 2208/2958, 121 | 2243/2958, 122 (errors 681 → 677) |

- **Reading.** The gain is largest on the 2024 additions, whose counts
  were the weakest. Their previous count errors came partly from truncated
  supports, not only from unusable images. This is a measurement, not a
  selection: it is not the optimization set (campaign C), and the 2025
  cohort was already read by earlier campaigns. The cohort-dependence of
  the assignment remains (full145 697/870 here, against 776/870 in the
  146-scan run).

### 2026-10-25 — Decision: full set of 492 scans; the elite subset is still to be defined

- **Review.** The user viewed the 30 scans of the 145 that the automatic
  rule rejects (`results/quality_6mer_20261025/review_rejected_145/`). Their
  judgment: these scans should not be rejected from the global set, except
  perhaps 240307_019 (several objects), but they are probably not the best
  images, so they do not belong in an elite subset.
- **Decision** (the user).
  - **Full set: 492 scans** = the 493 validated by eye minus 240307_019
    (`benchmarks/chitosan_6mer_full_v1.toml`: 144 + 260 + 88 per cohort).
  - **Elite: a subset to rework.** It is neither the 492 nor the 419. The
    419 scans (`benchmarks/chitosan_6mer_optimization_v1.toml`) are only
    what the current automatic rule accepts, and that rule is rejected as a
    definition of "elite": its thresholds were set by hand.
  - `benchmarks/chitosan_6mer_cu100_quality.tsv` records `in_full`,
    `rule_accepted` and an empty `elite` column.
- **Corrects the earlier wording.** The campaign Q entry called the 419 the
  optimization set; the optimization set will be the elite subset once
  defined. Campaign C steps 1–2 wait for it.
- **Measurement on the 419 run** (job 12384419, `results/opt419_d0d342e/`):
  counting 359/419 exact (419 within one), assignment 1889/2514 correct
  (75.1%), 105 exact chains. The 493-scan run gave 417/493 (84.6%) and
  2243/2958 (75.8%). Filtering by the current rule therefore does not
  improve the assignment: counting gains about one point and the assignment
  loses less than one point. Different cohorts, so this is indicative.

### 2026-10-25 — Declared campaign E: an automatic, label-free definition of the elite subset

- **Deliverable** (user request: an objective, automatic, efficient
  selection; the user does not want to rate images, and may ask colleagues
  later). A frozen image-quality score and an elite subset of the 492-scan
  full set, with the same standard for the 10–20mer. Time box: 2026-10-25
  to 2026-10-26.
- **Score** (image only: no fit, no count, no method output).
  - Seven axes, each oriented so that higher is better: contrast to noise;
    line discontinuity (inverted); stripe periodicity (inverted);
    forward/backward correlation; forward/backward NRMSE (inverted); second-object
    ratio (inverted); molecule-to-border margin in nm (new catalog value
    `roi_border_margin_nm`).
  - Each axis is converted to a percentile rank against a **frozen
    reference**: the 827 technically usable Cu(100) 6-mer scans of
    2026-10-25. The reference deciles are written in
    `config/scan_quality.toml`, so a scan's score does not depend on the
    set it is in, and the 10–20mer data are scored on the same scale.
  - Score = mean of the seven percentile ranks.
- **Validation without new labels.** The earlier human triage is used as
  weak labels: kept (493) versus rejected (395, overview rejections left
  out).
  - Reported: the AUC of the score against the triage, overall and with
    each session held out in turn.
  - The elite threshold is not chosen to maximise any method grade.
- **Threshold, fixed now.** The smallest score above which at least 98% of
  the scans (all scans with a triage decision) are human-kept, found on the
  training sessions only and applied to the held-out sessions. Elite = scans
  of the 492-scan full set above it. The held-out precision and the elite
  size are reported. If the elite has fewer than 80 scans, it is reported as
  too small and the target is not relaxed silently.
- **Then.** The method is graded on the full set and on the elite subset,
  on the existing 493-scan run (same cohort, nothing re-learned), with
  per-molecule results next to per-scan ones. This answers the open
  question on the sub-dataset. It is a measurement: nothing is re-selected
  after it.

### 2026-10-25 — Campaign E results: the declared rule gives an elite of 6 scans; the score is a weak proxy for the triage

- **Score** (`scripts/data/score_scan_quality.jl`, frozen reference deciles in
  `config/scan_quality.toml`, new catalog value `roi_border_margin_nm`).
  Median score 0.51 on the 827 usable scans, as expected for percentile
  ranks.
- **Validation against the earlier human triage**
  (`scripts/grading/report_elite_selection.jl`, outputs in
  `results/quality_6mer_20261025/elite/`):
  - AUC of the score against kept/rejected: **0.661** (827 scans, 493 kept).
  - Threshold at 98% precision on all sessions: 0.748. Leave-one-session-out:
    6 scans predicted elite, 5 of them human-kept (precision 0.833, below
    the 0.98 target).
  - **Elite subset: 6 of the 492 full-set scans.** The declared rule says
    that fewer than 80 scans is reported as too small and the target is not
    relaxed silently. It is reported as such.
- **Reading.**
  - The human triage rejected mostly for molecule integrity (cut molecule,
    no molecule, several objects, not a 6-mer), not for image sharpness.
    The seven axes measure image quality and miss these reasons, so the
    triage is a poor proxy for "best images".
  - With about 60% of the scored scans human-kept, a 98% precision target
    is very demanding for a score with AUC 0.66.
  - Any elite threshold from these data would therefore be a choice, not a
    validated measure. A reference needs human quality grades of the images
    themselves (the closing proposal of the presentation).
- **Track integrity** (`scripts/data/check_track_integrity.jl`, SXM
  headers and fitted lobe positions only; `results/dec493_c2b954a_integrity/`).
  On the 493-scan run: 74 tracks of at least 2 scans (449 scans), none with
  more than one session, no scan of another track interleaved, median span
  5.8 min, median centroid spread 0.6 nm. Four tracks are flagged
  (diagnostic thresholds from the collection limits, nothing selected):
  - track 58 (43 scans, the molecule behind 43 of 85 exact chains): chain
    axis turns by about 40° over the series, in small steps (median 1.2°,
    one step of 27°);
  - track 12 (11 scans): one 5.6 nm jump of the molecule centre between two
    consecutive scans, with the same axis (86–96°) and count;
  - track 25 (17 scans): counts range over 3 units;
  - track 17 (14 scans): a 110 min gap between two scans.
  The 440-scan 10–20mer run is not checked yet.

### 2026-10-25 — Declared: provisional elite subset = best third of the frozen quality score (user: option 1, "I leave the choice of the best to you")

- **Rule, fixed now.** Elite = the scans of the 492-scan full set whose
  frozen quality score (`config/scan_quality.toml`) lies in the top third
  of the full set's scores (about 164 scans). Image values only: no method
  output and no grade enters the choice.
- **Status.** Provisional and arbitrary in its size: the score is a weak
  proxy for the human triage (AUC 0.66), and no human quality grades exist.
  It stands until colleagues grade the images (closing proposal of the
  presentation), and is then redone against that reference.
- **Measurement, not selection.** The method is graded on the full set and
  on the elite, on the existing 493-scan run (classes learned on all 493
  scans, not on the elite alone). Per-scan results sit next to per-molecule
  ones, and nothing is re-selected after the grade.

### 2026-10-25 — Provisional elite results: selecting the best-scored images does not improve the grade

- **Sets.** Full set 492 scans (`benchmarks/chitosan_6mer_full_v1.toml`),
  provisional elite 167 scans (`benchmarks/chitosan_6mer_elite_v1.toml`; the
  top third of the frozen score, 167 rather than 164 because of ties at the
  threshold 0.565): full145 37, added_2024 72, added_2025 58. Graded on the
  existing 493-scan run (job 12384144, classes learned on all 493 scans).
- **Results** (final calls; counts after consensus):

  | Set | Scans | Exact N | Correct units | Exact chains | Errors |
  |---|---:|---:|---:|---:|---:|
  | Full | 492 | 417 (84.8%) | 2240/2952 (75.9%) | 122 | 675 |
  | Elite | 167 | 136 (81.4%) | 675/1002 (67.4%) | 32 | 314 |

- **Pooled, the elite looks worse; per molecule it does not.**
  - Per molecule (tracks), mean fraction of correct units: full 0.728 over
    117 molecules, elite 0.726 over 54.
  - On the 54 molecules in both sets, elite scans minus all scans of the same
    molecule: +0.004, bootstrap 95% interval [−0.005, +0.019].
  - The pooled gap is a composition effect: the elite holds 42 scans of the
    session `Cu100LHe30042025`, which is dominated by the 36-scan molecule
    whose assignment was already at 50% (2026-10-19).
- **Reading.** For the same molecule, the best-scored scans are not graded
  better than the others. With this score, restricting to a sub-dataset
  changes which molecules are in the set, not how well the method reads
  them. This answers the open question on the sub-dataset for this score
  only: it does not say that image quality never matters, only that these
  seven axes do not capture what drives the errors.
- **Not done and not to do.** The elite is not re-chosen on this result
  (that would select on the grade). The elite-only run, with classes
  learned on the elite, was not made.
- **Context for the headline.** The 89% (776/870) on the 145 scans comes
  from a run restricted to 146 scans. The same 145 scans score 697/870 (80%)
  in the 493-scan run, which makes the cohort dependence (open question 8)
  concrete.

### 2026-10-25 — Exploratory: inside the 6-mer, K and N look alike; the bright class mostly marks "not terminal"

- **Why.** The user doubted the edge-adjacency limit (the two GlcN, K, are
  not at the chain ends) and found it odd that K is not clearly brighter.
- **Measurement** (493-scan run `results/dec493_c2b954a/`, 105 molecules
  with 6-lobe scans of the full set, one mean per molecule; labels read
  only here, for grading). Mean amplitude relative to the chain maximum by
  position 1–6 (N K N N K N): 0.60, 0.83, 0.80, 0.82, 0.85, 0.63. Lobe
  widths along the chain: 0.35, 0.42, 0.41, 0.42, 0.42, 0.35 nm.
- **Paired per molecule, K (2, 5) minus interior N (3, 4):**
  - amplitude +0.034 (bootstrap 95% [+0.005, +0.062]; K higher in 70 of 105
    molecules);
  - width along +0.003 nm [−0.009, +0.014]; width across +0.002 nm
    [−0.005, +0.010].
- **K minus terminal N (1, 6):** amplitude +0.225 (K higher in 89 of 105),
  widths +0.07 and +0.05 nm.
- **Reading.**
  - The K are interior units next to the dim, narrow terminal units. Their
    contrast with the ends is large; against the interior N it is 3% in
    amplitude and nil in width.
  - So the bright class on K mostly marks "not terminal". The benchmark
    shows no size difference between K and interior N, which fits the
    earlier AUC ≤ 0.66 for position-independent features.
  - Chemistry predicts that the acetylated unit (N, GlcNAc) is bulkier; the
    data show no such difference inside the chain. The brightness may track
    ring orientation (the experimenters' readings encode "3H-up") more than
    the substituent. Not tested.

- **Follow-up, same day: K against the interior N with three brightness
  measures** (105 molecules, ratio K (2, 5) / interior N (3, 4) per molecule,
  bootstrap over molecules). Earlier figure above: +0.034 on the 0–1
  relative-amplitude scale, about +4%, a difference of means. As a ratio
  per molecule:
  - absolute amplitude +6.3% [+2.6, +10.1], N brighter in 36 of 105 molecules;
  - amplitude over the chain maximum +6.2% [+2.3, +10.0], N brighter in 35 of 105;
  - integrated (amplitude × widths) +9.8% [+3.7, +16.3], N brighter in 39 of 105.
- **By session** (amplitude over the chain maximum, sessions with at least 4
  molecules): K minus interior N goes from −4.3% (20241114, N brighter in 5
  of 8 molecules) and −2.2% (20240310) to +17.7% (20240815). The sign
  changes with the session, so tip or session effects are as large as the
  K–N difference.
- **Reading.** K is slightly brighter than the interior N on average, a few
  percent to about ten, but N is brighter in about a third of the molecules,
  and in two of eleven sessions on average. The brightness contrast between
  the two chemistries is small against the variation between molecules and
  sessions. The weighted position profile (all scans, 43-scan molecule
  included) differs from the per-molecule one and is not used.

### 2026-10-25 — Exploratory: terminal, sub-terminal and core units; a position-only rule beats the method on the 6-mer

- **Question** (user): should terminal N be distinguished from core N?
- **Measurement** (493-scan run, full set, 417 scans counted at 6 units;
  final calls against NKNNKN; labels read only for grading):

  | Units | Method | Rule "ends dim, interior bright" |
  |---|---:|---:|
  | All six | 78.1% | 66.7% |
  | Terminal N (1, 6) | 89.8% | 100% |
  | Interior (2–5) | 72.3% | 50% |
  | K (2, 5) | 79.9% | 100% |
  | Core N (3, 4) | 64.7% | 0% |

  - Per molecule (98 molecules), interior accuracy of the method minus the
    ends-only rule: +18.2 points, bootstrap 95% [+13.0, +23.2].
  - The method calls the K bright in 81% (position 2) and 79% (position 5)
    of the scans, and the core N bright in 41% (3) and 29% (4).
- **The control is a position code.** In NKNNKN a unit is K exactly when it
  is the second from either end. The rule "second from the end is bright"
  therefore scores 100% on every 6-mer, above the method's 78%. Terminal, K
  and core N are three position classes, and each coincides with one
  chemistry.
- **Reading.**
  - The method does more than mark the ends: it separates K from the core
    N, +18 points over the ends-only rule.
  - That separation is likely to come from adjacency to the dim terminal
    lobes (features comparing a lobe with its neighbours; AUC 0.75 with the
    terminal neighbour, 0.62–0.66 without, 2026-09-30), not from the
    chemistry itself. Position-neutral brightness differs by only 6% (above).
  - The 6-mer has no K with two interior neighbours on each side, and only
    two core N. How the method reads a K in the middle of a long chain is
    untested.
- **Consequence for the grades.** Report accuracy by position class
  (terminal, sub-terminal, core) next to the pooled figure, and compare
  with the position-only baseline. The pooled 76–78% includes terminal
  units that a trivial rule gets right.

### 2026-10-25 — Declared campaign R: benchmark-free reliability of per-unit signals, with the fixed-sequence 6-mer as negative control

- **Deliverable** (user, carte blanche on "a good 6-mer set that helps on the
  10–20mer, and objectifying the method without a benchmark"). Criteria that
  need no ground truth on the 10–20mer, and their values for the current
  data. Time box: 2026-10-25 to 2026-10-27.
- **Principle.** A property of a molecule that a method reads must come out
  the same when the same molecule is imaged again. Two label-free facts give
  a benchmark-free test:
  1. *Reproducibility across repeated scans* of each physical lobe bounds
     what any per-unit reading can contain (a reading cannot correlate with
     the truth better than the square root of its reliability).
  2. *The 6-mer is a negative control for sequence.* Every 6-mer molecule has
     the same sequence, and in NKNNKN the sequence is a function of the
     distance d to the nearest chain end. After removing the mean value at
     each d, any unit-level variation that repeats across scans of a 6-mer
     cannot be sequence. If the 10–20mer show no more reproducible unit-level
     variance than that, there is no evidence of sequence information in the
     signal.
- **Exploratory disclosure.** A Python prototype (not in the repository)
  preceded this entry, to check that the design is not degenerate. Its
  magnitudes: reproducibility of the amplitude residual (ICC) 0.3–0.4 in the
  6-mer and 0.2–0.3 in the 10–20mer; between-scan correlation of the
  residual profile 0.8 (0–3 min) to 0.6 (more than 30 min) in the 6-mer, 0.7
  to 0.5 in the 10–20mer; geometric covariates (lateral offset, spacing,
  turning angle, widths) explain 5–7% of the residual variance when
  cross-validated by molecule; no phase-random alternation along the chain.
  The repository script, fixed below, replaces the prototype and its numbers
  are the ones reported.
- **Definitions** (module `STMFit.UnitReliability`, script
  `scripts/data/report_unit_reliability.jl`, label-free).
  - Signal per lobe: log amplitude (primary), log integrated amplitude, log
    width along, log width across (secondary), from `scan_geometry.tsv`.
  - Residual: signal minus its scan median, minus the mean at the lobe's
    distance class d ∈ {1, 2, 3, ≥ 4} (mean of track means within the run).
  - Units: physical lobes seen in at least 3 scans of a track, tracks with
    at least 3 scans, at most 8 scans per track (evenly spaced in
    acquisition order), lobes with d ≥ 2 (d ≥ 3 and d ≥ 1 as sensitivity).
  - Reliability: one-way random-effects variance components over physical
    lobes (unbalanced design): σ_b² between lobes, σ_w² within, ICC =
    σ_b² / (σ_b² + σ_w²). Intervals: bootstrap over tracks, 1000 draws,
    seed 20261025.
  - Time lag: Pearson correlation of the residual profiles of two scans of a
    track (at least 4 shared lobes), by lag 0–3, 3–10, 10–30 and over 30 min,
    averaged by track then over tracks.
  - Count reproducibility: per track, the range of the per-scan counts.
- **Decision rules, fixed now.**
  - *Extra structure in the 10–20mer* only if the 95% interval of
    σ_b²(10–20mer) / σ_b²(6-mer) lies above 1, for the amplitude and for at
    least one other signal. Otherwise the report says that no reproducible
    unit-level information beyond the fixed-sequence level is evidenced.
  - A reading from amplitude is *feasible* only if the contrast it must
    detect (K against core N, measured on the 6-mer, positions 2 and 5
    against 3 and 4) is at least the between-lobe standard deviation σ_b
    (d′ ≥ 1); the accuracy ceiling of a single threshold is Φ(d′/2).
  - No variant is selected on these results: all are reported.

### 2026-10-25 — Elite-only run: learning the assignment on the elite alone gives 44% correct units

- **Run.** `results/elite167_df4a4fa/` (Viper job 12386344, 46 min, code
  `df4a4fa`): the 167 provisional-elite scans through the promoted method
  with the decoupled supports, so counting, consensus, classes and fusion
  are all learned on these scans only. All label-free checks pass; the
  geometry used the counting support for 10 of 167 scans.
- **Grade** (elite manifest, `results/elite167_df4a4fa_grade/`):

  | Run on the elite scans | Exact N | Correct units | Exact chains |
  |---|---:|---:|---:|
  | classes learned on the 493 scans | 136/167 | 675/1002 (67%) | 32 |
  | classes learned on the elite alone | 136/167 | 441/992 (44%) | 0 |

  - Correct units by position 1–6 (N K N N K N): 42, 59, 4, 1, 76, 68%,
    against 68, 62, 51, 74, 79, 86% with the 493-scan classes.
  - Counting is unchanged (136 exact in both): it does not depend on the
    cohort.
- **Mechanism.** The classes learned on 54 molecules put 70% of the units in
  the bright class (44% with the 493-scan classes on the same scans): 96%
  and 99% at the core N (positions 3 and 4) and 58% at the first unit. The
  label-free fusion fit shows it: π 0.66, false-call rate θ0 0.38, detection
  rate θ1 0.96, against 0.42, 0.05 and 0.74 for the 493 scans. The clusters
  no longer separate the ends from the interior.
- **Reading.**
  - This is the dataset dependence of the assignment (open question 8) at
    its extreme: 89% (146-scan run), 80% (the same 145 scans in the 493-scan
    run), 76% (full set), 67% (elite, 493-scan classes) and 44% (elite
    alone), all from one method. "Optimising on a better subset" is not
    supported by this run.
  - A label-free alarm exists: θ0 and π of the fusion fit. θ0 is 0.03–0.05
    in the 6-mer runs of 146 to 493 scans, 0.38 here, and 0.19 on the
    10–20mer run, which therefore sits between the two.

### 2026-10-25 — Campaign R results: per-unit signals repeat, but no more in the 10–20mer than in a fixed-sequence 6-mer; amplitude cannot read K against N

- **Correction of the exploratory disclosure.** The prototype's 10–20mer
  reproducibility (ICC 0.2–0.3) was wrong: it grouped the lobes of scans
  whose count differs from the track's, which have no physical lobe, into one
  pseudo-lobe per track. The repository code skips them. Its numbers
  replace the prototype's (the 6-mer values and the time-lag values were
  unaffected). The prototype's covariate regression (explained variance 5% and
  7%) did not use physical lobes and stands as an exploratory remark.
- **Data.** 6-mer, 493-scan run `dec493_c2b954a`: 53 tracks of at least 3
  scans (407 scans). 10–20mer, 440-scan run `dec1020_70a9baa`: 65 tracks (352
  scans). Outputs in `results/reliability_20261025/` (script
  `scripts/data/report_unit_reliability.jl`) and `results/ceiling_20261025/`
  (`scripts/grading/report_contrast_ceiling.jl`). The module has tests
  (`test/test_unit_reliability.jl`).
- **Reliability of the amplitude residual** (d ≥ 2, at most 8 scans per track;
  between-lobe sd, within-lobe sd, ICC with its 95% interval):

  | Set | Tracks | sd between | sd within | ICC |
  |---|---:|---:|---:|---|
  | 6-mer (fixed sequence) | 52 | 0.145 | 0.168 | 0.43 [0.34, 0.50] |
  | 10–20mer | 50 | 0.130 | 0.151 | 0.43 [0.33, 0.50] |
  | 10–20mer, N 7–9 | 33 | 0.143 | 0.153 | 0.47 [0.35, 0.56] |
  | 10–20mer, N ≥ 10 | 12 | 0.105 | 0.143 | 0.35 [0.22, 0.47] |

  - Every 6-mer has the same sequence and, once the mean at each distance to
    the end is removed, nothing that depends on the sequence is left. The
    persistent unit-level variation of 14.5% in amplitude is therefore not
    sequence. It comes from the molecule and its environment.
- **Rule on extra structure** (fixed in advance): the ratio of between-lobe
  variances 10–20mer over 6-mer is 0.81 [0.52, 1.21] for the amplitude, 0.99
  [0.59, 1.63] for the integrated amplitude, 1.20 [0.63, 2.76] along and 1.20
  [0.61, 2.30] across. No interval lies above 1: **no evidence of
  reproducible unit-level information in the 10–20mer beyond the
  fixed-sequence level.**
- **Persistence.** The correlation of the residual profile between two scans
  of a molecule, amplitude, by lag (0–3, 3–10, 10–30, over 30 min): 6-mer
  0.81, 0.67, 0.73, 0.62; 10–20mer 0.69, 0.61, 0.51, 0.56. It decays slowly, so
  the pattern is a property of the molecule, not of a passing tip state.
- **Order along the chain.** The lag-1 autocorrelation of each molecule's
  persistent profile is −0.267 in the 6-mer and −0.196 in the 10–20mer, against
  −0.243 [−0.329, −0.154] and −0.164 [−0.241, −0.081] when the profiles are
  shuffled along the chain: no alternation beyond chance.
- **Reading ceiling** (control NKNNKN, 105 molecules, `010010`). Contrast of the
  K lobes over the core N lobes (positions 3 and 4), in log amplitude:
  +0.037 [−0.004, +0.074], position effect included, so an upper bound for
  chemistry. Against the between-lobe sd of 0.145, d′ = 0.25 and the best
  single threshold reaches 55% [51, 60] for two equally likely classes
  (53% from one scan). A contrast of 0.24 would be needed for 80%, six times
  the measured one. Width along and across, and the integrated amplitude, give
  54%, 52% and 52%. On the 146-scan run (38 molecules) the amplitude contrast
  is +0.061 [+0.005, +0.124] and the ceiling 59% [51, 69].
- **The 6-mer cannot test an N/K reading.** In NKNNKN a unit is K exactly
  when it is the second from either end, so all the molecules carry the same
  class at a fixed position: there is no variation of the label that is not a
  variation of the position. Any N/K score on this control measures the
  sensitivity to position. A position-aware model (idea of the end effect as
  a nuisance, with per-session normalisation) cannot be graded on it, and is
  not run.
- **Counting is the measurable weak point of the 10–20mer.** Range (max − min)
  of the per-scan counts of one molecule, tracks of at least 3 scans:

  | Run | Tracks | range ≤ 1 | range ≤ 2 | range > 2 |
  |---|---:|---:|---:|---:|
  | 6-mer, decoupled supports | 53 | 50 | 52 | 1 |
  | 10–20mer, legacy support | 65 | 24 | 38 | 27 |
  | 10–20mer, decoupled supports | 65 | 29 | 47 | 18 |

  - The decoupled supports improve the repeat consistency of the 10–20mer
    counts without any reading by the experimenters (27 → 18 tracks with a
    range above 2), which confirms the adoption on an independent criterion.
  - 18 of 65 repeated molecules still differ by more than 2 units between
    scans (up to 9). The scans that deviate from their track carry no more
    quality flags than the others (partial scan 16% against 6%, several
    objects 8% against 1%, nothing else), and most unstable tracks have none.
  - The lobe spacing is constant: 0.607 ± 0.082 nm per gap in the correctly
    counted 6-mers, and the mean spacing of a scan varies by only 0.03 nm.
    In the unstable 10–20mer tracks the spacing stays normal (0.6–0.7 nm)
    while the detected extent of the chain changes between scans of one
    molecule (track 80, 8 scans: 3.5 to 7.4 nm, counts 6 to 12). The
    instability is in the detection of the extent of a heterogeneous chain (a
    bright blob at one end and a dim chain), not in the resolution of the
    units.
  - In the 6-mer the per-scan errors are as often over- as under-counts (57
    and 51), and taking the largest count of a track lowers the exact counts
    from 88.9% (median) to 64.9%: a longest-extent rule is not supported.
  - The image-quality score does not predict the 6-mer count errors either
    (AUC 0.49 for the per-scan count, 0.46 after consensus).
- **Metadata check of the tracks** (`scripts/data/check_track_integrity.jl`).
  10–20mer, 85 tracks of at least 2 scans: none spans two sessions, none has
  a scan of another track interleaved, the median centre spread is 1.1 nm
  (largest 4.3 nm) and one track has an axis spread above 30°. The only
  frequent flag is the count range (20 tracks).
- **What the data support.**
  - Positions repeat (0.07–0.11 nm), the metadata link the scans of a track
    soundly, and the counts of a molecule are stable in about half of the
    repeated 10–20mer.
  - A per-unit N/K reading from these amplitudes is not supported: the
    contrast to detect is below the unit-level variation that is specific to
    each molecule, the 10–20mer show no extra reproducible structure, and
    the 6-mer control cannot validate a reading.
  - Experiments that would change this: homo-oligomers (all GlcN and all
    GlcNAc) give the chemical contrast with no position confound; images at
    several biases test whether the contrast depends on the tunnelling
    conditions; a non-palindromic or interior-K control gives label
    variation at fixed position.

### 2026-10-25 — Declared campaign D: how the 6-mer training subset changes the frozen assignment of the 10–20mer

- **Deliverable** (user idea, 2026-10-25: other groups and sub-data sets of the
  6-mer, and their impact on the 10–20mer). A measurement of the dependence of
  the 10–20mer output on the choice of the 6-mer subset the assignment is
  learned on. Not a selection: no subset is chosen from the result.
- **Models.** The fully frozen, label-free assignment of campaign A
  (`config/unit_assignment_frozen_labelfree.toml`, reference rule: corroborated
  scans, at most 5 per track) learned on a subset of the scans of the 493-scan
  6-mer run, then applied unchanged to the 440-scan 10–20mer run. Features are
  the run's own (standardised per scan, Fisher frozen on the reference), so a
  scan's features do not depend on the subset.
- **Subsets** (of the 492-scan full set unless stated): the full set; the 167
  elite scans; the 419 rule-accepted scans; the 145 original scans (cohort
  `full145`); the 2024 sessions; the 2025 sessions; the full set with a cap of
  1, 3, 10 scans per track and with none; and the noise floor, 5 random
  167-scan subsets (scans drawn at random) and 5 random molecule subsets (random
  tracks until 167 scans are reached), seeds fixed.
- **Metrics.**
  - 10–20mer: the fused call of every lobe; pairwise agreement of the calls of
    two models (lobes classified by both); class-1 fraction by distance class
    to the nearest end (1, 2, 3, ≥ 4).
  - 6-mer (grading side): accuracy by position class (terminal, K, core N) on
    the scans outside the model's training subset, against the position-only
    rule (100%) and the ends-only rule (66.7%); fitted θ0, θ1, π.
- **Reading rules, fixed now.**
  - The noise floor is the pairwise agreement among the random subsets. A
    structured subset changes the 10–20mer output beyond sampling if its
    agreement with the full-set model falls below the 5th percentile of the
    noise floor (agreement of random subsets with the full-set model).
  - A subset is not called better or worse on the 6-mer alone: stability on the
    10–20mer says nothing about correctness, and the 6-mer cannot validate an N/K
    reading (campaign R).

### 2026-10-25 — The elite-only result (44%) is not a bug: the scan-level elite broke the molecules; diagnostics

- **Why this entry.** The user found the 44% suspect: a subset of better
  images should not do worse than the full set. The checks below look for a
  bug before accepting it.
- **No bug, three independent checks.**
  - *Features.* The per-scan features of the 1007 lobes the elite run and the
    493-scan run share are identical for 974 of them (6 scans differ).
  - *Naming.* In the elite-only run the class 1 is the brighter one (mean
    relative amplitude 0.82 against 0.58), as in the 493-scan run (0.88
    against 0.66).
  - *Replay.* Learning the assignment from the elite run's own corroborated
    scans with the frozen-assignment code (no cap, learned weights, the run's own
    Fisher) reproduces the run's calls exactly (100% identical) and its fit
    (θ0 0.382, π 0.66 with fusion). The collapse is therefore reproduced
    offline.
- **Cause.** The elite was chosen scan by scan, and that fragmented the repeat
  structure the method relies on.
  - The 167 elite scans come from 54 molecules (tracks of the 493-scan run);
    38 of them keep fewer than 3 scans (47 scans) and only 16 keep at least 3.
  - In the elite run, 88 scans are corroborated, in 13 molecules; the two
    largest hold 53% (31 and 16 scans). The production assignment has no cap
    per molecule, so the classes are learned on 13 molecules, two of which
    dominate.
  - The same elite-run data with a cap of 5 scans per molecule (frozen design,
    equal weights, frozen Fisher) give 80.1% correct units (terminal 98.5%, K
    64.3%, core N 77.6%), and with a cap of 1 scan 72.7%. Frozen models learned
    on the elite from the 493-scan features give 69–76% on the elite scans and
    about 84% on the other scans, whatever the cap, weights or Fisher.
- **Reading.** The elite-only run measures the fragility of a method that
  re-learns its classes on a small, molecule-poor cohort without a cap. It does
  not show that better images are worse. The cap per molecule of campaign A
  removes the failure. The unit of a dataset selection is the molecule: a
  selection that splits the scans of a molecule below the 3 scans the consensus
  needs destroys the corroboration. A molecule-level elite is declared next.

### 2026-10-25 — Declared: molecule-level elite (v2), after the scan-level elite broke the molecules

- **Rule, fixed now.** The unit of selection is the molecule: a track of the
  493-scan run `dec493_c2b954a` (scans grouped from SXM headers and image
  registration, label-free). Candidates: the tracks of the 492-scan full set
  with at least 3 scans. Molecule score: the mean frozen image-quality score of
  its scans. Elite v2: the better half of the candidate molecules, with all
  their scans, at most 8 per molecule (evenly spaced in file order), so that no
  molecule dominates and every kept molecule can still be corroborated.
- **Status.** Provisional, with the weak score (AUC 0.66 against the triage;
  it does not predict count errors). The scan-level elite (v1) stays recorded
  as a construction that did not work.
- **Measurements, not selection.** (a) The existing 493-scan run graded on
  the elite v2 scans (classes learned on all the scans); (b) an elite-v2-only
  run on Viper (counting, consensus, classes and fusion learned on those scans
  only). Nothing is re-chosen from the results.

### 2026-10-25 — Campaign D results: the 6-mer subset moves one call in five and the bright fraction of the 10–20mer from 17% to 41%

- **Run.** `scripts/grading/report_training_stability.jl` on the 493-scan 6-mer
  run (`dec493_c2b954a`) and the 440-scan 10–20mer run (`dec1020_70a9baa`):
  20 frozen label-free models (the subsets declared above), outputs in
  `results/training_stability_20261025/` (untracked: `models.tsv`,
  `agreement.tsv`, `application_calls.tsv`, `report.md`). The grading side
  reads the control sequence only for the 6-mer accuracies.
- **Noise floor** (rule fixed in the declaration). Agreement with the
  full-set model of the 10 random subsets: mean 81.4%, 5th percentile 78.5%,
  range 78.3–85.7%. Pairwise among random scan subsets 94.9% (93.0–97.8%),
  among random molecule subsets 91.9% (86.6–94.7%). About one lobe in five
  changes call between the full-set model and a model learned on 167 random
  scans.
- **Structured subsets.** Only the 145 original scans fall below the noise
  floor (76.7%, 5th percentile 78.5%). The others are inside it or above:
  cap 1 per molecule 81.0%, cap 3 87.7%, cap 10 96.2%, no cap 95.9%,
  elite 85.2%, rule-accepted 93.4%, sessions 2024 89.3%, sessions 2025 87.3%.
- **Bright fraction of the 10–20mer.** 16.5% to 41.3% over the 20 models
  (full set 37.7%; cap 1: 22.3%, cap 3: 31.7%, cap 10 and no cap: 41.3%;
  random 167-scan subsets 16.5–29.8%). It tends to grow with the size of the
  reference set (63–97 reference scans: 16.5–29.8%; 192: 37.7%; 331: 41.3%).
  A factor 2.5 in the number of units called bright comes from the subset the
  classes were learned on.
- **6-mer, scans outside the training subset.** 76.0–85.6% for every subset
  except sessions 2024: a model learned on the 404 scans of the 2024 sessions
  gives 58.5% on the 88 scans of the 2025 sessions (terminal 76.9%, K 42.9%,
  core N 55.8%), while the reverse (2025 → 2024) gives 85.6%. The held-out
  sets differ between models, so the comparison is indicative; a session
  effect is possible and the pooled benchmark hides it.
- **Per-unit robustness (descriptive, thresholds chosen when reading).** Over
  the 3624 lobes of the 10–20mer run classified by the models, 65% are dim in
  at least 80% of the 20 models, 19% bright in at least 80%, and 15% in
  between. By distance to the nearest end: d = 1: 84% dim, 3% bright, 13% in
  between; d = 2: 42% dim, 40% bright, 18% in between; d = 3: 64 / 21 / 15%;
  d ≥ 4: 70 / 15 / 15%. The stable bright calls sit at d = 2, the position of
  the K of the 6-mer: the 10–20mer output mostly repeats the position code
  of the training molecule (campaign R).
- **Reading, by the rules fixed above.** The output depends on the subset
  beyond sampling only for the 145 original scans, but the sampling noise is
  itself large (one call in five). Stability across subsets says nothing about
  correctness: the 6-mer cannot validate an N/K reading, and a call that
  survives the 20 models is a call the method repeats, not a verified K or N.
  No subset is chosen from this table. For the application, a per-unit
  robustness flag (share of the 20 models that call the unit bright, from
  `application_calls.tsv`) is reported with every sequence, so that unstable
  units are not read as units.

### 2026-10-25 — Molecule-level elite v2, measurement (a): the 493-scan run graded on the 141 elite scans

- **Result.** Classes learned on all 492 scans (run `dec493_c2b954a`), graded
  on the 141 scans of the 26 elite-v2 molecules: 601 correct of 834 classified
  units (846 possible; 72.1%), 115/141 exact N, 38 exact chains.
- **Status.** Measurement only, nothing re-chosen. Measurement (b), the
  elite-v2-only run on Viper (job 12387149, `elite2_241f156`), is pending; the
  72.1% is the baseline it is compared with. Presentation page 4 and the
  fusion/consensus counts of that run are updated once it is graded.

### 2026-10-25 — Molecule-level elite v2, measurement (b): the elite-v2-only run does not collapse and is not clearly better

- **Run.** Viper job 12387149 (35 min, `elite2_241f156`, code `241f156`,
  `config/chitosan_decoupled_support.toml`, the promoted pipeline unchanged):
  the 141 scans of the 26 elite-v2 molecules, counting, consensus, classes and
  fusion learned on those scans only. Graded with
  `benchmarks/chitosan_6mer_elite_v2.toml` (`results/elite2_241f156_grade`,
  untracked).
- **Results** against measurement (a), the 493-scan run on the same 141 scans
  (classes learned on all 492):

  | | exact N | correct units | errors | exact chains |
  |---|---:|---:|---:|---:|
  | (a) 493-scan run | 115/141 | 601/846 (71.0%) | 233 | 38 |
  | (b) elite v2 only | 105/141 | 627/846 (74.1%) | 203 | 34 |

  (b) classifies 830 of 846 units (75.5% correct among them); fusion π 0.25
  (truth 0.33), θ0 0.069, θ1 0.843. Correct units by position 1–6 (scans of
  six units, identity alignment): 93, 71, 88, 79, 41, 98%, against 87, 84, 77,
  63, 61, 83% for (a): the K at position 5 is the weak unit of (b).
- **Reading.** The scan-level collapse (44%, previous entries) came from the
  fragmented molecules, as diagnosed: with whole molecules the elite-only run
  gives 74%, no worse than the full-set classes. It is not better on exact
  counts and exact chains, which the consensus of fewer scans per molecule
  (at most 8) lowers, and the gain in units is +3.1 points of the possible
  units. One run of 141 scans, no repeat: the difference is a measurement, not
  evidence that the elite is better. The quality score that ranked the
  molecules is a weak proxy (AUC 0.66 against the triage), so this does not
  test "better images are better" either. No subset is promoted.
- **Presentation.** Pages 2 and 4 (French and English) show the elite v2 with
  the full-set classes (71%) and learned on the elite (74%); the scan-level
  elite is mentioned only here.

### 2026-10-26 — 10–20mer sequences of the frozen model with per-unit robustness (exploratory, label-free)

- **Deliverable.** The sequences of the curated working set (19 scans, 18 of
  them in the 440-scan run) called by the frozen label-free assignment learned
  on the full 6-mer set, each unit with its robustness over the 20 models of
  campaign D. No subset is chosen and nothing is validated.
- **Tool.** `scripts/data/list_application_sequences.jl` (reads the
  `application_calls.tsv` of campaign D and a table of scans; no label),
  `UnitReliability.call_robustness` (bright if at least 80% of the models call
  the unit bright, dim if at most 20%, unstable between; the 80% is descriptive,
  not tuned), tests `test/test_application_sequences.jl` and
  `test/test_unit_reliability.jl`. Output `results/sequences_20261026/`
  (untracked).
- **Result.** 166 units in 18 scans: 117 robustly dimmer, 35 robustly brighter,
  14 unstable. Of the 35 robust bright units, 20 sit at the second unit from a
  chain end (20 of 36 units at that distance), 5 at the third (of 35) and 10
  further in (of 59); none of the 36 terminal units. The bright calls follow
  the position code of the 6-mer control (K second from either end), as
  campaign R found.
- **Against the earlier listing.** The sequences on the presentation page came
  from the transductive 440-scan run (classes learned on the 10–20mer
  scans). 124 of the 166 units (75%) get the same call with the frozen model,
  and the bright units drop from 76 to 54: about the one call in five that
  campaign D measures between training subsets. The page now shows the frozen
  model, with unstable units drawn as outlines.
- **Documentation.** README, `docs/src/index.md`, the runbook (benchmark-free
  checks and tool table) and `AGENTS.md` carry the conclusions of campaigns R
  and D. The presentation (French and English, `results/presentation_20261025/`,
  untracked) has the elite v2 on pages 2 and 4 and these sequences on page 7.
- **Docs build fix.** `docs/make.jl` failed because the generated `journal.md`
  page (257 KiB) exceeded Documenter's 200 KiB HTML limit. The size check is
  now disabled (`size_threshold = nothing`): the journal is one long page on
  purpose. The build passes (`GKSwstype=100 julia --project=. docs/make.jl
  --build-only`).
- **Tests.** The full suite passed on 2026-10-25 before this change (20 runs,
  end-to-end included); the two test files touched here pass on their own.

### 2026-10-26 — Why the elite v2 scores below the full set: session composition, not image quality (exploratory, grading side)

- **Question (user).** The elite v2 scans score below the full set (71% against
  76% correct units in the 493-scan run). Why would better images do worse?
- **Measurements** (493-scan run, grade of the 492-scan full set, frozen
  quality score of campaign E):
  - Elite v2 scans 71.0% correct units, the other 351 scans 77.8%; per molecule
    (molecules of at least 3 scans) 71.7% against 78.0%, so the gap is not a
    scan-weighting effect. The cap of 8 scans per molecule costs little (the
    elite molecules with all their scans: 73.7%).
  - Within a session, elite and other scans do equally well: the difference,
    session by session and weighted by elite scans, is +0.7 point (11
    sessions hold both).
  - 80 of the 141 elite scans come from the three hardest sessions: 20241113
    (40 scans; 67% correct units for the whole session), 20241114 (24; 54%) and
    30/04/2025 (16; 50%).
  - Spearman correlation of the quality score with per-scan accuracy: −0.27
    over all scans, +0.01 within sessions. The strongest single axis is the
    contrast to noise (−0.34 overall, −0.16 within sessions).
- **Reading.** The elite is not worse because its images are better. The
  score ranks sessions more than scans, and it favours the sessions where the
  method does worst; inside a session it says nothing about how well the
  method reads a scan. Session (tip state, imaging conditions) dominates the
  result, so subsets must be compared within sessions or balanced across them.
  Why high contrast-to-noise sessions read worse is open (a tip state that
  changes the lobe appearance is one hypothesis, not tested). The image score
  is not a valid definition of an elite; a human rating (the proposal of the
  last slide) or an image criterion validated against it is still needed.

### 2026-10-26 — Declared campaign J: joint fit of the repeated scans of one molecule (stage 1: joint count from the per-scan GCV curves)

- **Deliverable** (user, approved 2026-10-26; time box 5 working days, until
  2026-10-31). A molecule-level count, then geometry and unit shape, that uses
  every repeat scan of a molecule in one criterion instead of a per-scan
  choice followed by a majority vote. Target: the 10–20mer counting
  instability (18 of 65 repeated molecules differ by more than 2 units).
  Stages: (1) joint count from the per-scan fits already computed, no new fit;
  (2) a joint fit with shared geometry and an end-lobe type (own amplitude,
  width and spacing for the terminal lobe; not for the second lobe, which is
  always K in the control, so its own parameters would learn the control's
  sequence); (3) unit shape as a diagnostic (repeatability, end against core
  N, which are both N in the control), synthetic long chains for counting.
- **Stage 1 method, fixed now.** Every counting run writes, per scan, the GCV
  of each candidate N for the elliptical and circular models
  (`<scan>/ell_scores.tsv`, `circ_scores.tsv`). The effective score of a scan
  is the minimum of the valid elliptical and circular GCV at each N (the
  Level 1 of the per-scan hierarchy; reconstructed exactly from the tables:
  487/487 scans of the 493-scan run, 429/429 of the 10–20mer run). For a track
  of at least 3 scans (`consensus/consensus.tsv`, as for the consensus), the
  joint score is J(N) = Σ_s [log GCV_s(N) − min_N log GCV_s(N)], over the N
  valid in at least one scan of the track; a scan where N is not valid
  contributes its own largest log ratio. N_joint = argmin J (ties: the smaller
  N); scans without any valid candidate do not contribute. Tracks of fewer
  than 3 scans keep the per-scan count. No parameter is added, no label read.
  The per-scan layers (robust-AICc guard, support midpoint) are not applied to
  J at this stage.
- **Comparisons.** (a) per-scan `N_selected`; (b) the current consensus
  `N_final` (strict majority of `N_selected`); (c) the vote of the Level-1
  counts `N_eff` (strict majority, else median); (d) `N_joint`. On the 6-mer
  (493-scan run `dec493_c2b954a`, 492-scan full set; grading side,
  `grade_chitosan_benchmark.jl`): exact N per scan and per molecule, by
  session. On both runs, label-free: split-half agreement. The scans of each
  track of at least 4 scans are split in time order into odd and even; each
  half gets its own count by each rule (votes: strict majority, else median);
  report the share of tracks where the halves agree and the mean |ΔN|.
- **Decision rule, fixed now.** J becomes the count criterion of stage 2 if,
  on the 6-mer, its exact N per scan is at least that of the current consensus
  and, on the 10–20mer, its split-half agreement is higher than that of the
  current consensus. If J beats the vote of `N_eff` but not the current
  consensus, the per-scan layers carry information that J lacks, and stage 2
  builds them into the joint criterion (pooled support length and guard). If J
  beats neither, stage 2 keeps the vote for counting and works on the shared
  geometry only.
