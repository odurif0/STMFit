# STMFit monorepo

Analysis pipeline for STM images of molecular chains (chitosan on Cu(100) and
similar systems). Detects and fits a chain-of-Gaussians model to count the
number of monomer units (lobes) per chain, label-free.

**Latest working candidate, September 22:** complete-disk support improves
**675 → 676 correct / 870** and **33 → 34 exact chains / 145**. Use the opt-in
`config/unit_assignment_patch_support.toml`. Coverage rises **849 → 852 / 870**;
classified accuracy is **676/852 (79.3%)**. Three newly available predictions
are correct, but two formerly correct decisions regress. This is **not a new
champion**: historical 677/870 and 36/145 still lead, with higher coverage (854).
The independent final Gaussian-score arm regresses to **666/870 and 32/145**
and is rejected. All 900 keys and selected counts remain; no pixels are imputed,
no variants combined and no settings tuned after grading. Labels are used only
by external grading; this reused development benchmark is not independent
validation. See the [support/score comparison](docs/src/unit_assignment.md#complete-disk-support-and-final-gaussian-score-2026-09-21).

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
