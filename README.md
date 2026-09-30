# STMFit

Label-free analysis of STM images of molecular chains (chitosan on Cu(100)):
count the units of each chain, locate them, and assign each unit GlcN (0) or
GlcNAc (1) with an explicit uncertainty, starting from raw Nanonis `.sxm` scans.

## Current method and benchmark result

The promoted method (September 29, 2026) has four stages. None of them reads a
benchmark label, expected count, sequence or class proportion:

1. **Per-scan counting.** A chain of 2D Gaussians is fitted along the molecular
   axis for each candidate N. GCV selects the count, with a robust-AICc guard
   and a support-midpoint adjustment (`config/chitosan.toml`).
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

Fresh raw-to-prediction run (Viper job 12025539, source `3c414a4`), external
grading on the 145-scan NKNNKN 6-mer benchmark (870 positions):

| full145 own-N | Exact N /145 | Correct /870 | Exact chains /145 | Classified /870 | Errors |
|---|---:|---:|---:|---:|---:|
| Historical reference (Aug 2026) | 106 | 677 | 36 | 854 | 177 |
| Previous record (saved hybrid counts) | 129 | 694 | 43 | 855 | 161 |
| **Promoted method** | **137** | **772** | **88** | **865** | **93** |

All 145 counts are within one lobe. Fits use a time-limited global search, so
reruns vary slightly: four runs gave 772, 772, 772 and 769 correct, each with
88 exact chains and 137/145 exact N. **Limits:** the NKNNKN control confounds
GlcNAc with edge adjacency, so benchmark accuracy is not chemical validation.
Exact chains come in molecule-sized blocks (43 of 88 from one 43-scan molecule).
Only repeatedly imaged molecules are fused. The counting parameters keep their
historical benchmark-informed provenance ([calibration](docs/src/calibration.md)).
On unknown 10–20mers, some long chains are undercounted (truncated support) and per-scan calls are
inconsistent across repeat scans (fitted false-call rate 0.45). Those outputs
are therefore not chemically validated.

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
julia --project=. test/inspect_one_file.jl data/chitosan_6mer/20240817_LHe_Cu100/240817_004.sxm

# Per-scan counting only
julia -t 4 --project=. test/batch_full.jl 146 --config config/chitosan.toml \
    --data-dir RAW_DIR --outdir results/counting

# Complete method, raw scans -> predictions.tsv (run it on Viper for a cohort)
julia -t 4 --project=. test/run_molecule_consensus_chitosan.jl \
    --data-dir RAW_DIR --count-config config/chitosan.toml \
    --config config/unit_assignment_corroborated_training.toml \
    --consensus-config config/molecule_consensus.toml \
    --templates templates/chitosan_cc_molds_native_v1.tsv --outdir NEW_DIR

# External verification and benchmark grade (labels are read only here)
julia --project=. test/grade_consensus_run.jl --run NEW_DIR --outdir NEW_GRADE_DIR
```

Raw `.sxm` files are not tracked; see [`data/README.md`](data/README.md).
Outputs go to `results/` (ignored). Cohort runs belong on the MPCDF cluster:
see [`hpc/README.md`](hpc/README.md) and the [runbook](docs/src/chitosan_runbook.md).

## Repository layout

| Path | Content |
|---|---|
| `packages/` | `STMFitCore` (constraints, scoring), `STMSXMIO` (SXM I/O, preprocessing), `GaussianFit1D` (diagnostic 1D fit), `GaussianFit2D` (2D chain fit engine), `STMMolecularFit` (orchestration, selectors). `STMMolecularFitGUI` is a separate, unmaintained GUI. |
| `test/` | Command-line workflows (`run_*`, `build_*`, `extract_*`, `grade_*`) and their tests (`test_*.jl`, `molecule_consensus/`). Shared code in `test/lib/`. |
| `config/` | Frozen TOML configs: counting, consensus, assignment (promoted plus the lineage configs used by tests), calibration audit, new-molecule template. |
| `templates/` | Frozen DFT constant-current molds, their provenance and the DFT frames. |
| `benchmarks/` | External grading manifests and control sequences. Never read by the method. |
| `hpc/` | Slurm scripts for Raven/Viper, and the QE input record for the DFT molds. |
| `docs/src/` | Documentation (Documenter): method, configuration, runbook, journal. |

## Tests

```bash
julia --project=. packages/STMFitCore.jl/test/runtests.jl
julia --project=. packages/STMSXMIO.jl/test/runtests.jl
julia --project=. packages/GaussianFit2D.jl/test/runtests.jl
for t in test/test_*.jl test/molecule_consensus/*.jl; do julia --project=. "$t"; done
julia --project=. test/test_reconstructed_pipeline.jl --e2e   # synthetic end to end
```

## Documentation

- [Pipeline](docs/src/pipeline.md): data flow, stages and outputs.
- [Model selection](docs/src/selection.md): per-scan count selection and molecule consensus.
- [Unit assignment](docs/src/unit_assignment.md): GlcN/GlcNAc assignment and fusion.
- [Configuration](docs/src/config.md) and [Calibration](docs/src/calibration.md).
- [Runbook](docs/src/chitosan_runbook.md): benchmark, grading and 10–20mer application.
- [DFT molds](docs/src/qe_stm_molds.md) and [DFT calculation note](docs/src/dft_calculation_note.md).
- [Research journal](docs/src/journal.md): decisions, rejected approaches, open questions.

Build: `GKSwstype=100 julia --project=. docs/make.jl --build-only`.
