# STMFit monorepo

Analysis pipeline for STM images of molecular chains (chitosan on Cu(100) and
similar systems). Detects and fits a chain-of-Gaussians model to count the
number of monomer units (lobes) per chain, label-free.

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
- [**Selection**](docs/src/selection.md) — the label-free selection rule (GCV + robust-AICc guard + up-when-ambiguous).
- [**Calibration**](docs/src/calibration.md) — parameter objectivation, auto-calibration, GCV rationale.
- [**Config reference**](docs/src/config.md) — every parameter and flag.
- [**Chitosan runbook**](docs/src/chitosan_runbook.md) — reproducible benchmark workflow.
- [**Research journal**](docs/src/journal.md) — dated decision log (the project's memory).
