# STMFit monorepo

Analysis pipeline for STM images of molecular chains (chitosan on Cu(100) and
similar systems). Detects and fits a chain-of-Gaussians model to count the
number of monomer units (lobes) per chain, label-free.

**New measured benchmark reference, September 24: saved promoted counts improve recognition.**
Connecting the saved `support_midpoint_hybrid` counts to unchanged
support assignment gives **694/870 correct units, 43/145 exact chains and
855/870 coverage**. This exceeds historical **677/36/854** on all three measures;
the matched fresh control gives **676/34/852**. Exact counts are **129/145**
versus 106; missing/extra positions fall **16/38 → 13/6**, emitted errors
**176 → 161** against control. Both arms rebuild all geometry and descriptors;
no labels, composition prior or post-grade tuning enter inference.

This is the new **development-benchmark reference**, not a proven optimum,
independent validation or a fresh reproduction of the saved counting policy.
The unknown25 application and its frozen outputs are unchanged. Source
**71de012**, job **11956079**, **0:0 in 35m20s**; all 433 outputs plus the log
are local and verified. A split-kernel allocation regression was repaired before this
matched comparison; timed-fit variability is reported, not assumed absent.
See [the counting-to-recognition result](docs/src/unit_assignment.md#saved-promoted-counts-and-recognition-2026-09-24)
and `results/promoted_counts_20260924/report.md`. Older comparisons below retain
their dated conclusions, not the current reference status.

**Provenance caveat (September 24 audit):** inference reads no labels, but the
historical support/width calibration and hybrid thresholds were selected using
known-count benchmark grades. This is not a wholly label-free calibration or an
independent validation. Two fresh raw-to-prediction repetitions are in progress;
unqualified champion promotion is withheld. See the [audit](docs/src/journal.md#2026-09-24--fresh-hybrid-reproduction-and-historical-label-use-audit).

**Earlier comparison, September 24: correct QE cube addressing does not improve recognition.**
With the remaining surface calibration unchanged, corrected templates give
**676/870 correct, 27/145 exact chains**, versus byte-identical legacy tangent
control **679/33**, both coverage **852/870**. Six scans gain, nine lose; one
exact chain gained, seven lost. Keep the opt-in reader correction, but do not
promote this assignment variant. QE writes the last grid index fastest; old
outputs remain empirical references, not validation of the old coordinates.
N, geometry, patches, Fisher and classifier settings stay fixed; labels enter
only external grading, with no retuning. Historical **677/36/854** is not
surpassed overall. Source **8295965**, job **11955512**, **0:0 in 6m22s**;
all 460 outputs/log fetched and independently checked. See the
[tangent comparisons](docs/src/unit_assignment.md#tangent-physical-mold-scores-2026-09-23-trade-off).

A subsequent read-only audit identifies periodic-substrate contamination,
an inadequate surface-calibration bracket and signed values in the accepted
GlcNAc observable. No new gabarit or grade is claimed; see the
[DFT qualification](docs/src/dft_calculation_note.md#diagnostic-constant-current-observable).

Earlier September 24, leave-target-out parity/mirror selection gives **670/32**
against **679/33**, both coverage **852**: rejected without tuning. Source
**61a56d6**, job **11955091**; `results/mold_loo_20260924/report.md`.

**Previous support reference (saved inputs):** complete-disk support improves
**675 → 676 correct / 870** and **33 → 34 exact chains / 145**. Use the opt-in
`config/unit_assignment_patch_support.toml`. Coverage rises **849 → 852 / 870**;
classified accuracy is **676/852 (79.3%)**. Three newly available predictions
are correct, but two formerly correct decisions regress. This is **not a new
champion**: this support reference trails historical 677/870 and 36/145,
with higher historical coverage (854).
The independent final Gaussian-score arm regresses to **666/870 and 32/145**
and is rejected. All 900 keys and selected counts remain; the support rule imputes no patch pixels,
no variants combined and no settings tuned after grading. Labels are used only
by external grading; this reused development benchmark is not independent
validation. See the [support/score comparison](docs/src/unit_assignment.md#complete-disk-support-and-final-gaussian-score-2026-09-21).

The **September 23 cross-view residual comparison** is negative: same-view
subtraction gives **675/870 correct, 30/145 exact**, opposite-view subtraction
**673/30**, versus exactly replayed support **676/34**, all at **852/870 coverage**.
Only subtraction amplitudes/backgrounds change; N, geometry, main features and
raw patches remain literal. All 292 linear solves converge and all fused means
pass native validity. Retain support, below historical **677/36**; no promotion
or retuning. This is not independent cross-validation.
See the [cross-view result](docs/src/unit_assignment.md#cross-view-residual-subtraction-2026-09-23-negative).

The **September 23 locally oriented Gaussian comparison** is negative:
**664/870 correct, 30/145 exact**, versus matched global refit **674/31** and
byte-identical saved reference **676/34**, all at **852/870 coverage**. The
label-free GCV pool selects local models on 87/146 scans; lower GCV does not
improve recognition. All 584 fits pass native validity but hit their iteration
caps, so they are not converged optima. N, patch sampling and the split cache
stay fixed. Retain support, below historical **677/36**; no promotion or tuning.
See the [local-model result](docs/src/unit_assignment.md#locally-oriented-gaussian-model-2026-09-23-negative).

The earlier **September 23 local patch-frame comparison** is negative: **670/870
correct, 30/145 exact**, versus exactly replayed global support **676/34**,
at unchanged **852/870 coverage**. Two exact chains are gained, six lost.
Only patch sampling rotates along a label-free quadratic centerline; N and
both saved molecular fits remain unchanged. Keep support, below historical
677/36; no angle tuning or promotion. That ablation did not test a locally rotated
molecular fit. See the [local-frame result](docs/src/unit_assignment.md#local-chain-tangent-patch-frames-2026-09-23-negative).

The **original-support registered refit** resolves the three earlier span
failures, but remains **inconclusive**: 145/146 scans complete all four fits;
a different scan fails the unchanged residual guard after registration. All
original supports replay exactly. No partial classifier or grade is produced.
Keep saved support **676/870 correct, 34/145 exact**, below historical **677/36**.
See the [original-support follow-up](docs/src/unit_assignment.md#original-support-registered-refit-2026-09-22-inconclusive).

The **September 23 multi-start/subpixel pilot** adds no fully valid paired fit.
Four deterministic starts improve some circular fits, not eligibility. Residual
2D translation reduces paired RSS **4.93–6.64%** and all eight block-heldout
errors **0.53–13.47%**, but every selected shift saturates **−1 pixel in both
axes**. The bounded pilot is inconclusive for larger shifts, not a proof that
registration cannot help. No new recognition grade or promotion; keep **676/34**,
below **677/36**. See the
[subpixel pilot](docs/src/unit_assignment.md#multi-start-and-subpixel-pilot-2026-09-23-bound-limited).

The preceding **registered native refit** is **inconclusive**, not a new
benchmark result. At saved N, three scans fail the unchanged support-span
constraint (two in the zero-shift control); 143/146 complete all four shape
fits. No scan is excluded and no partial classifier or grade is produced.
The fifteen reference tables replay exactly. Retain support **676/34**, below
historical **677/36**. See the
[registered-refit diagnostic](docs/src/unit_assignment.md#registered-native-refit-2026-09-22-inconclusive).

The earlier fixed-geometry **acquisition registration** comparison is negative:
**671/870 correct, 26/145 exact**, versus **676/34** for both native replay
and the observation-mask control. Coverage falls **852→850/870**. Registration
is accepted on 104/146 scans and improves all 629 comparable patch correlations
(median **0.23→0.99**), but loses eight exact chains and gains none. Keep saved
support; better image agreement is not better recognition here. See the
[acquisition-registration result](docs/src/unit_assignment.md#fixed-geometry-acquisition-registration-2026-09-22-negative).

The matched **geometric variable-projection** comparison is negative:
**632/870 correct, 20/145 exact**, versus **676/34** for both saved support
and the same-start joint optimizer, with unchanged **852/870 coverage**.
The profiled fit lowers RSS on all 146 scans but adds **44 classification
errors**; against the joint control, two exact chains are gained and sixteen
lost. N and classifier settings stay fixed. Retain saved support, below
historical 677/36; no promotion or post-grade tuning. See the
[geometric-profile result](docs/src/unit_assignment.md#matched-geometric-variable-projection-2026-09-22-negative).

Reprofiling Gaussian amplitudes and background at **identical saved N and
geometry** is also negative: **675/870 correct, 32/145 exact**, versus the
byte-identical support control **676/34**, at unchanged **852/870 coverage**.
All 146 linear fits satisfy KKT and reduce RSS (median **0.681%**), but four
scans gain and five lose correct positions; two exact chains are lost and
none gained. Retain saved support, below historical 677/36. That experiment
included no nonlinear refit, selection, threshold adjustment or unknown25 rerun. See the
[frozen-amplitude result](docs/src/unit_assignment.md#frozen-geometry-amplitude-profiling-2026-09-22-negative).

Fresh raw-GCV count selection after symmetric filtering is **negative**:
**661/870 correct, 32/145 exact**, versus **675/34 for a freshly refitted
fixed-N control** and **676/34 for saved support**. Exact counts fall
**106 → 101/145**, coverage **852 → 849/870**, and emitted errors rise
**177 → 188** against the fresh control. The latter differs from saved support
on three decisions; refit variability prevents attributing every change solely
to N. Retain saved support; historical 677/36 remains unexceeded. This tests
the original direct-extractor GCV lineage, not the separately promoted hybrid
counting policy. No tuning or unknown25 rerun follows. See the
[reselection result](docs/src/unit_assignment.md#fresh-raw-gcv-count-selection-2026-09-22-negative).

The complete-patch-training follow-up is **negative**: **674/870 correct,
33/145 exact**, with coverage falling to **850/870**. One scan has no complete
training patch and loses two previously available predictions; another loses
an exact chain. Keep the **676/34 support candidate**, not this alternative.
No normalization fallback or post-grade tuning is added. See the
[complete-training result](docs/src/unit_assignment.md#complete-patch-training-2026-09-21-negative).

The median/IQR GMM-normalization follow-up is also **negative: 667/870 correct,
10/145 exact**, at unchanged **852/870 coverage**. It loses 24 exact chains
without gaining any against the exactly replayed 676/34 control. Reject this
variant and keep `config/unit_assignment_patch_support.toml`; historical 677/36
remains the target. No clipping, rescaling or post-grade tuning follows. See the
[robust-normalization result](docs/src/unit_assignment.md#robust-per-scan-gmm-normalization-2026-09-21-negative).

Equal-scan GMM weighting gives a **tradeoff: 675/870 correct, 36/145 exact**,
at unchanged **852/870 coverage**, versus the replayed 676/34 support control.
It gains two exact chains but loses one correct position; all seven changed
decisions land at the existing zero-margin vote tie. Keep the **676/34 support
candidate** as the primary working reference. The weighted variant is not
promoted and does not exceed historical 677/36. No post-grade tuning follows.
See the [scan-weighting result](docs/src/unit_assignment.md#equal-scan-gmm-weighting-2026-09-21-exact-chain-gain-unit-loss).

Continuous internal GMM voting is **negative: 671/870 correct, 26/145 exact**,
at unchanged **852/870 coverage**, versus the replayed 676/34 control. It gains
one exact chain but loses nine. All 24 final changes leave old zero-margin ties
and switch 1→0; the learned groups, naming and final threshold stay fixed.
Keep the **676/34 support candidate**; no post-grade tie or precision adjustment
follows. Historical 677/36 remains unexceeded. See the
[continuous-vote result](docs/src/unit_assignment.md#continuous-gmm-seed-vote-2026-09-21-negative).

Whole-scan GMM bagging is **negative: 665/870 correct, 10/145 exact**, at
unchanged **852/870 coverage**, versus the exactly replayed 676/34 control.
Twenty bootstrap replicates lose 24 exact chains and gain none. Of 62 final
1→0 changes, 57 leave old GMM=1 / k-means=0 ties under the unchanged fusion
rule. Reject this variant and keep the **676/34 support candidate**; historical
677/36 remains unexceeded. No post-grade adjustment follows. See the
[scan-bagging result](docs/src/unit_assignment.md#whole-scan-gmm-bagging-2026-09-21-negative).

Sharing covariance throughout GMM learning is **negative: 666/870 correct,
6/145 exact**, at unchanged **852/870 coverage**. It loses 28 exact chains
and gains none against the exactly replayed 676/34 support control. Reject
this variant and retain **`config/unit_assignment_patch_support.toml`**;
historical 677/36 remains the target to exceed. No other method or post-grade
adjustment is added. See the
[tied-covariance result](docs/src/unit_assignment.md#tied-gmm-covariance-throughout-learning-2026-09-22-negative).

The two subsequent independent tests bring **no improvement**: whole-scan
Fisher cross-fitting gives **668/870 correct, 34/145 exact**; within-scan GMM
naming gives **676/870, 34/145**, with unchanged predictions. Both retain
**852/870 coverage**. Fisher gains four exact chains but loses four; all 26
changed decisions enter zero-margin vote ties. Relative naming agrees with
the old names for all ten unchanged GMM fits. Keep the **676/34 support
candidate**, below historical 677/36. No combination or post-grade tuning
follows. See the [Fisher/naming comparison](docs/src/unit_assignment.md#whole-scan-fisher-and-relative-gmm-naming-2026-09-22-no-improvement).

The fixed rank-four factor-analyzer and df-five Student learning tests are
**negative: 620/870 correct, 23/145 exact**, and **631/870, 25/145**,
respectively, at unchanged **852/870 coverage**. Keep the **676/34 support
candidate**, below historical 677/36. All ten Student fits converge; the
factor fits reach the fixed 200-update cap, so this is not a verdict on a
fully converged MFA optimum. No iteration/rank/df or tie-rule tuning follows.
See the [mixture comparison](docs/src/unit_assignment.md#factor-analyzer-and-student-learning-2026-09-22-negative).

The subsequent all-update covariance shrinkage and coherent Student-density
tests also regress: **637/870 correct, 24/145 exact**, and **627/870, 25/145**,
respectively, at unchanged **852/870 coverage**. All ten fits per candidate
meet their stopping criterion; control tables replay exactly. Retain the
**676/34 support candidate**, below historical 677/36. No post-grade tuning
or unknown-chain rerun follows. See the
[shrinkage/Student-density comparison](docs/src/unit_assignment.md#all-update-shrinkage-and-coherent-student-decisions-2026-09-22-negative).

**Research branch, September 20:** the symmetric fused-fit correction improves
the regenerated reconstruction from **666 to 671 correct / 870 controls** and
**24 to 28 exact chains / 145**, at identical cached counts and coverage
(**78.4% → 79.0%** classified accuracy). Both arms were refitted; the control
reproduces the archived predictions byte-for-byte. The historical reference
remains better (677 correct / 36 exact), so the correction is **not promoted**.
Count selection was not rerun; the counting numbers below describe the frozen
pre-correction reference.
See the [journal](docs/src/journal.md#2026-09-20--symmetric-fused-fit-correction-and-matched-recognition-comparison).

The subsequent fixed Fisher-mirror experiment is **negative**: physical-u
reflection gives **669/870 correct and 27/145 exact**, versus the reproduced
671/28 symmetric-fusion control at identical coverage. It is not retained;
the target remains exceeding the historical 677/36, not merely the native baseline.

The subsequent matched-residual comparison (`S(data-model)` vs `S(data)-model`)
gains **one correct position and one exact chain**: **672/870, 672/849 (79.2%),
29/145 exact**, at unchanged counts/coverage. Both arms regenerate all patches;
the control is byte-identical to the saved symmetric result. This opt-in variant
became the **working candidate at that stage, not a promoted champion**: history
then led by five correct positions and seven exact chains. No setting was tuned
after grading; this reused benchmark does not establish independent generalization.
See the [matched-residual result](docs/src/unit_assignment.md#matched-residual-filtering-2026-09-20-modest-gain).

The next fixed comparison improves the descriptor to **673/870 correct,
673/849 (79.3%), 33/145 exact chains** by computing half-plane asymmetry after
affine projection of the same patch. The first-moment alternative gives 672/28
and is not retained. `config/unit_assignment_affine_residual.toml` became the
**opt-in working candidate at that stage, not a champion**: history led by four correct
positions and three exact chains. Coverage and counts are unchanged. All sixteen
changed decisions land at the existing vote tie (zero margin), so this modest
development-benchmark gain is not evidence of calibrated chemical confidence.
See the [descriptor comparison](docs/src/unit_assignment.md#two-transverse-descriptor-candidates-2026-09-20-affine-gain-no-promotion).

The September 21 signed-CC / affine-Fisher comparison is **negative for both
independent changes**: signed template margins give **668/870 correct, 32/145
exact**; affine Fisher gives **672/870, 32/145**, against the exactly replayed
673/33 control. Counts and 849/870 coverage are unchanged. Inference uses no
benchmark labels; the separate grade is reused development evidence, not
independent validation. That comparison kept the affine-descriptor candidate at **673/33**;
historical **677/36** remains the target. No combined variant or post-grade
tuning follows. See the [signal comparison](docs/src/unit_assignment.md#signed-cc-and-affine-fisher-candidates-2026-09-21-both-negative).

**Benchmark (6mer):** the robust-AICc guard validates at 39/39 primary 240817
files exact (N=6), reproducible. The frozen pre-correction chitosan reference was promoted for
the expanded 145-file external counting grade: 129/145 exact, 143/145 within one
lobe. The 0/1/? unit-assignment benchmark uses the same 145 files; its external
control sequence is NKNNKN (010010/101101 by convention) for grading only, never
for fitting, selection, thresholding, abstention, or method calibration.
**Historical counting application (10–20mer):** 25/25 files processed, without
known sequences. Visual QC is possible, but it does not establish chemical
assignment accuracy.

**Unit assignment (frozen promoted reference):** label-free soft vote of the
k-means 4-view and the GMM 1-view + adaptive-contour (constant-current) mold
margins + empirical Fisher-discriminant mold margin: 79.3% classified
physical accuracy / 36 exact chains / 677 correct of 854 — the promotion bar
(78.9% / 18 / 677) is met; the Fisher mold generalizes under half-split
cross-validation (66.3% per-lobe, no overfit). Regrading the saved outputs
reproduces these accepted metrics; a full raw-input rebuild of this exact
champion lacks the irrecoverable original `patch_u_asym` producer.

**Approved native reconstruction:** `test/run_reconstructed_chitosan.jl` uses
Julia 1.13 and the separately named `cc_soft_reconstructed_v1`. Its descriptor
is fixed before comparison; the champion is not relabeled as reproduced.
Native component and synthetic end-to-end tests pass. After the adaptive-support
handoff repair, Raven job **30278010** completed the unknown25 application on
September 17 in 6 min 16 s. All **25 scans / 222 selected lobes** are retained;
local integrity checks pass and QC is reproduced exactly. The model outputs
**189 class-0 and 33 class-1 assignments, with no `?`**. However, **19/25 chains
are flagged for review** for low mean confidence. Confidence is an uncalibrated
vote margin; zero abstentions do not establish certainty or chemical accuracy.
The original failed run is preserved. Full146 **Viper job 11786116** completed on
September 18 in 11 min 04 s, retaining all **146 scans / 900 cached GCV lobes**:
**695 class-0 / 198 class-1 / 7 `?`**, with 43 chains flagged for review.
External full145 own-N grading gives **78.4% (666/849 classified positions)**
and **24/145 exact chains**; 16 control positions are missing and 38 predicted
lobes are extra. This is below the frozen reference, so the reconstruction is
**not promoted** or claimed as an exact recovery. The authorized run and external
evaluation are complete; historical-champion recovery and unknown-chain chemical
validation remain unresolved. The 900-lobe cache is not the separate promoted
counting summary (871 lobes). See
[`docs/src/unit_assignment.md`](docs/src/unit_assignment.md#explicit-julia-reconstruction-2026-09-16)
for the command, numerical differences and limitations.

## Packages

| Package | Role |
|---|---|
| `packages/STMFitCore.jl` | Shared physical constraints and scoring (κ penalty, spacing, residual diagnostics). |
| `packages/STMSXMIO.jl` | Shared SXM (Nanonis) I/O: `SXMImage`/`read_sxm` + preprocessing helpers. Used by both fit engines. |
| `packages/GaussianFit1D.jl` | 1D slide-profile fitting engine (now diagnostic-only, off by default). |
| `packages/GaussianFit2D.jl` | 2D Gaussian chain fitting engine (the core). |
| `packages/STMMolecularFit.jl` | Orchestration: slide extraction, selectors, batch comparisons. |
| `packages/STMMolecularFitGUI.jl` | GUI. |

## Quick start

```bash
# Inspect a single file (deep 2D analysis, no batch overhead)
julia --project=. test/inspect_one_file.jl 240817_004.sxm

# Full batch (48 files, default config, 1D diagnostic off for speed)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
    --config config/chitosan.toml

# Summarize a results TSV
julia --project=. test/summarize.jl results/best_plots/summary_overlap060_hard.tsv
```

Runtime outputs → `results/`. Raw SXM data goes in `data/` (see
[`data/README.md`](data/README.md) for where to get the benchmark datasets and
how to organize them).

## Workflows (`test/`, each standalone)

| Script | Purpose |
|---|---|
| `batch_full.jl [N] [--chunk i/n]` | Full 2D batch: fits, plots, enriched summary. `--skip-1d` (default) for speed; `--no-skip-1d` to add the 1D diagnostic. |
| `inspect_one_file.jl <file.sxm>` | Deep 2D ell vs circ on a single file. |
| `measure_calibration.jl <scan.sxm>` | **Auto-calibrate**: derive all objectivable parameters from one clean scan → emits a ready-to-use TOML. |
| `sensitivity_thresholds.jl {generate\|submit\|local\|compare}` | Measure robustness of N_selected to the selection threshold. |
| `diagnose_neff.jl`, `diagnose_fullimg_autocorr.jl` | Effective-sample-size and spatial-correlation diagnostics. |
| `summarize.jl [summary.tsv]` | Print stats from a summary TSV. |
| `run_reconstructed_chitosan.jl` | Native Julia 1.13 reconstructed raw-to-0/1/? pipeline, explicit selected-N caches, QC and maps; no grading. |
| `compare_reconstructed_champion.jl` | Separate keyed comparison with a saved reference; no fitting or parameter selection. |
| `run_unknown_unit_assignment.jl` | Label-free unknown-sequence 0/1/? runner: predictions, summary, manifest, validation logs. |
| `validate_unit_predictions.jl`, `summarize_unknown_unit_qc.jl`, `plot_unit_assignment.py` | Validate, QC, and plot unknown unit-assignment outputs without grading labels. |

## Configs (`config/`)

| Config | Use |
|---|---|
| `chitosan.toml` | Default 6mer chitosan (hand-tuned reference). |
| `chitosan_10_20mer_adaptive_support_rescue.toml` | 10–20mer production (long chains, adaptive support rescue). |
| `chitosan_auto.toml` | Auto-calibrated chitosan (zero hand-tuning, validates the objective method). |
| `template.toml` | Annotated template for calibrating a new molecule. |
| `*_rescue*.toml` | Variants with adaptive support rescue / aggressive settings. |

Each config has `[model]` (physical calibration), `[selection]` (selection
thresholds), and `[preprocessing]` (SXM channel/flatten) sections. See
[**Calibration**](docs/src/calibration.md) for the parameter classification
(measured / principled / free) and the auto-calibration workflow.

## HPC (MPCDF — Raven / Viper)

Large batches run as a Slurm **job array** on the MPCDF cluster: each array task
runs one `--chunk i/n` slice of the file list, so the batch is N× faster with N
parallel tasks. A push-button launcher handles sync + submit + merge + fetch.

```bash
cp hpc/remote.env.example hpc/remote.env && $EDITOR hpc/remote.env
./hpc/launch_remote.sh --dry-run     # preview
./hpc/launch_remote.sh --watch       # sync → submit → wait → merge → fetch
```

See [`hpc/README.md`](hpc/README.md) for setup (SSH/2FA, partitions, Julia
module, where to put code vs data) and tuning.

## Documentation

Full docs in `docs/src/` (built with Documenter):

- [**Pipeline & architecture**](docs/src/pipeline.md) — data flow, component roles.
- [**Selection**](docs/src/selection.md) — the label-free selection rule (GCV + robust-AICc guard + support-midpoint hybrid).
- [**Calibration**](docs/src/calibration.md) — parameter objectivation, auto-calibration, GCV rationale.
- [**Config reference**](docs/src/config.md) — every parameter and flag.
- [**Unit assignment**](docs/src/unit_assignment.md) — GlcNAc/GlcN per-lobe assignment pipeline (label-free, work in progress).
- [**DFT-STM molds**](docs/src/qe_stm_molds.md) — Quantum ESPRESSO LDOS mold workflow for unit assignment.
- [**Math background**](docs/src/math.md) — model form, selection criteria, effective sample size.
- [**API reference**](docs/src/api.md) — package/module index.
- [**Chitosan runbook**](docs/src/chitosan_runbook.md) — reproducible benchmark workflow.
- [**Research journal**](docs/src/journal.md) — dated decision log (the project's memory).
