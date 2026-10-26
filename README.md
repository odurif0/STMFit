# STMFit

Label-free analysis of STM images of molecular chains (chitosan on Cu(100)):
count the units of each chain, locate them, and assign each unit to one of two
physical classes (0 = dimmer, 1 = brighter) with an explicit uncertainty,
starting from raw Nanonis `.sxm` scans. The chemical identity of the classes
(GlcN or GlcNAc) is not established: on the 6-mer control the brighter class
falls on the two GlcN (journal, 2026-10-22).

## Current method and benchmark result

The promoted method (September 29, 2026; decoupled supports since October 25)
has four stages. None of them reads a benchmark label, expected count,
sequence or class proportion:

1. **Per-scan counting.** A chain of 2D Gaussians is fitted along the molecular
   axis for each candidate N. GCV selects the count, with a robust-AICc guard
   and a support-midpoint adjustment. The count uses the capped axial support;
   every later fixed-N refit uses the legacy support, or the capped one where
   that fit fails (`config/chitosan_decoupled_support.toml`).
2. **Repeated-scan molecule consensus.** Consecutive scans of one molecule are
   registered in the absolute piezo frame (NCC ≥ 0.5). A track of ≥3 scans with
   a strict-majority count refits the disagreeing scans at that count
   (`config/molecule_consensus.toml`).
3. **Per-scan unit assignment.** Local features, backward/forward patches, DFT
   constant-current mold margins and an empirical Fisher margin feed a k-means
   4-view and a GMM 1-view classifier. A soft vote combines them. The
   classifiers are learned only on scans whose count the consensus corroborates
   (`config/unit_assignment_corroborated_training.toml`).
4. **Latent-class fusion.** Per physical lobe of a track, the per-scan calls are
   binomial observations of a hidden type. EM fits the detection rate,
   false-call rate and mixing weight without labels. The posterior is the
   reported per-lobe uncertainty.

Fresh raw-to-prediction run (Viper job 12380719, source `70a9baa`), external
grading on the 145-scan NKNNKN 6-mer benchmark (870 positions; N = GlcNAc,
K = GlcN; "correct" means the brighter class falls on the two GlcN):

| full145 own-N | Exact N /145 | Correct /870 | Exact chains /145 | Classified /870 | Errors |
|---|---:|---:|---:|---:|---:|
| Historical reference (Aug 2026) | 106 | 677 | 36 | 854 | 177 |
| Previous record (saved hybrid counts) | 129 | 694 | 43 | 855 | 161 |
| Legacy support (deterministic, job 12042058) | 137 | 771 | 88 | 865 | 94 |
| **Promoted method** (decoupled supports, job 12380719) | **137** | **776** | **85** | **865** | **89** |

All 145 counts are within one lobe, and every final count equals the
legacy-support run's. The decoupled supports were adopted on 2026-10-25 by the
user's decision after the grade (the declared rule required no regression on
any figure; exact chains went from 88 to 85 while errors went from 94 to 89),
for their long-chain counts. Since 2026-09-30 the fits are deterministic:
identical runs give identical outputs. **Limits:** in the NKNNKN control the two GlcN sit next to the chain
ends, so the brighter class coincides with both GlcN and edge adjacency, and
benchmark accuracy is not chemical validation. The former naming rule
"brighter = GlcNAc" is contradicted by the control. The assignment's benchmark
signal is the contrast with the dimmer terminal lobe. Without it, amplitude is
only weakly higher at GlcN than at interior GlcNAc (AUC ≤ 0.66) and the DFT
molds barely separate them ([journal](docs/src/journal.md)).
Exact chains come in molecule-sized blocks (43 of 85 from one 43-scan molecule).
Only repeatedly imaged molecules are fused. The counting parameters keep their
historical benchmark-informed provenance ([calibration](docs/src/calibration.md)).
On unknown 10–20mers the chain lengths are not known. The decoupled supports
remove the truncated supports, and their counts are closer to 17 lengths read
by the experimenters (mean difference 1.0 unit, against 1.5 with the legacy
support), but those readings are a comparison, not a ground truth. Per-scan
calls stay inconsistent across repeat scans (fitted false-call rate 0.19), and
the unit order does not match the experimenters' acetylation readings beyond
chance. Those outputs are not validated.
Benchmark-free checks (2026-10-25, [journal](docs/src/journal.md)): the per-unit
amplitude repeats across scans of a molecule, but no more in the 10–20mer than
in the fixed-sequence 6-mer, and by itself reads K against N with a ceiling of
about 55%; the 6-mer, a position code, cannot validate an N/K reading. The
6-mer subset the assignment is learned on moves one 10–20mer call in five, and
18 of 65 repeated 10–20mer molecules differ by more than 2 units in count
across their scans. An N/K reading of the 10–20mer is not supported;
homo-oligomers, a bias series or a non-palindromic control are the experiments
that can test it.

## Quick start

Bootstrap once (Julia 1.13; local packages are declared in `Project.toml [sources]`):

```bash
julia --project=. -e '
using Pkg
for p in ["STMFitCore","STMSXMIO","GaussianFit1D","GaussianFit2D","STMMolecularFit"]
    Pkg.develop(PackageSpec(path="packages/$p.jl"))
end
Pkg.instantiate(); Pkg.precompile()'
```

```bash
# Inspect one scan (fits and plots, no batch)
julia --project=. scripts/tools/inspect_one_file.jl data/sets/6mer_Cu100_benchmark/240817_004.sxm

# Per-scan counting only
julia -t 4 --project=. scripts/pipeline/batch_full.jl 146 --config config/chitosan_decoupled_support.toml \
    --data-dir RAW_DIR --outdir results/counting

# Complete method, raw scans -> predictions.tsv (run it on Viper for a cohort)
julia -t 4 --project=. scripts/pipeline/run_molecule_consensus_chitosan.jl \
    --data-dir RAW_DIR --count-config config/chitosan_decoupled_support.toml \
    --config config/unit_assignment_corroborated_training.toml \
    --consensus-config config/molecule_consensus.toml \
    --templates dft/molds/chitosan_cc_molds_native_v1.tsv --outdir NEW_DIR

# External verification and benchmark grade (labels are read only here);
# with --data-dir it also sorts one plot per molecule into
# 1_all_exact / 2_N_exact_calls_wrong / 3_N_wrong
julia --project=. scripts/grading/grade_consensus_run.jl --run NEW_DIR --outdir NEW_GRADE_DIR --data-dir RAW_DIR
```

Raw `.sxm` files are not tracked; see [`data/README.md`](data/README.md).
Outputs go to `results/` (ignored). Cohort runs belong on the MPCDF cluster:
see [`hpc/README.md`](hpc/README.md) and the [runbook](docs/src/chitosan_runbook.md).

## Repository layout

| Path | Content |
|---|---|
| `src/` | `STMFit` library shared by the scripts: molecule consensus, fusion, unit assignment, patch preprocessing, molds, Fisher views. |
| `packages/` | `STMFitCore` (constraints, scoring), `STMSXMIO` (SXM I/O, preprocessing), `GaussianFit1D` (diagnostic 1D fit), `GaussianFit2D` (2D chain fit engine), `STMMolecularFit` (orchestration, selectors). |
| `scripts/pipeline/` | The promoted method, raw scans to predictions, positions and plots. Label-free. |
| `scripts/grading/` | External grading and benchmark reports, the data index and sets: the only code that reads labels. |
| `scripts/data/`, `scripts/calibration/`, `scripts/tools/` | Raw-scan catalog and repeat collection; calibration audit; single-scan inspection. |
| `test/` | Tests only; `test/runtests.jl` runs them all. |
| `config/` | Frozen TOML configs: counting, consensus, assignment (promoted plus the lineage configs used by tests), calibration audit, new-molecule template. |
| `benchmarks/` | External grading manifests, control sequences, data triage tables. Never read by the method. |
| `dft/` | DFT inputs (`inputs/`), frozen constant-current molds with provenance and frames (`molds/`), accepted QE states (`runs/`, untracked). |
| `hpc/` | Slurm script and notes for Raven/Viper. |
| `docs/src/` | Documentation (Documenter): method, configuration, runbook, journal. |
| `data/`, `results/` | All data (raw store, sets, index) and run outputs; untracked. |

## Tests

```bash
julia --project=. test/runtests.jl          # every test, synthetic end-to-end run included
julia --project=. test/runtests.jl --fast   # without the end-to-end run
```

## Documentation

- [Pipeline](docs/src/pipeline.md): data flow, stages and outputs.
- [Model selection](docs/src/selection.md): per-scan count selection and molecule consensus.
- [Unit assignment](docs/src/unit_assignment.md): two-class unit assignment and fusion.
- [Configuration](docs/src/config.md) and [Calibration](docs/src/calibration.md).
- [Runbook](docs/src/chitosan_runbook.md): benchmark, grading and 10–20mer application.
- [DFT molds](docs/src/qe_stm_molds.md) and [DFT calculation note](docs/src/dft_calculation_note.md).
- [Research journal](docs/src/journal.md): decisions, rejected approaches, open questions.

Build: `GKSwstype=100 julia --project=. docs/make.jl --build-only`.
