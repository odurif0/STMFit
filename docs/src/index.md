# STMFit — STM Molecular Chain Fitting

Automated pipeline for detecting and fitting 2D Gaussian chain models
to STM images of molecular chains. **Validated** on the 6mer chitosan/Cu(100)
benchmark (39/39 exact); **applied** to the 10–20mer production system (no
ground-truth labels — visual validation is the arbiter). Generalizable to other
chain-like molecules on the same STM via auto-calibration
(see [Calibration](calibration.md)).

## Quick Start

```bash
# Single-file diagnostic
STMFIT_DATA_DIR=/path/to/data julia --project=. test/inspect_one_file.jl 240817_004.sxm

# Full batch (default: chitosan.toml, GCV + robust-AICc guard)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml

# Auto-calibrate for a new molecule from one clean scan
julia --project=. test/measure_calibration.jl path/to/clean_scan.sxm

# Raw GCV baseline (no guard) for comparison
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml --selection-policy gcv

# Summarize results
julia --project=. test/summarize.jl results/best_plots/
```

Experimental selectors (blocked CV, support-marginalized GCV, slope-heuristic
MDL, stability, Laplace, fwd/bwd consensus, local-lobe evidence) and synthetic
validation are documented in [Model Selection](selection.md).

## Calibration

Auto-calibration (`test/measure_calibration.jl`) derives σ, spacing, fit width,
support, and n_max from one clean scan. See
[**Calibration**](calibration.md) for the parameter classification (measured /
principled / free) and why GCV is the canonical criterion (not BIC/AICc).

The default `config/chitosan.toml` is the hand-tuned reference; an auto-derived
equivalent (`config/chitosan_auto.toml`) validates the objective method.

For the chitosan reference set, `benchmarks/chitosan_240817.toml` records
evaluation-only quality classes. It is **not** used by fitting code and must not
become a selection prior.

## Unit Assignment (GlcNAc/GlcN)

The pipeline can assign each fitted lobe a monomer type to produce a
deacetylation map per chain. This is a **work in progress**. See
[Unit Assignment](unit_assignment.md) and [QE STM Molds](qe_stm_molds.md).

## Research Journal

All experimental paths — successful and failed — are documented in the
[Research Journal](journal.md). **Update this journal** whenever you:
- Test a new approach (even if it fails)
- Change the pipeline or model selection logic
- Discover a bug or convergence issue
- Add or remove parameters

## Components

| Component | Role |
|-----------|------|
| `STMFitCore.jl` | Shared utilities: κ penalty, spacing constraints |
| `STMSXMIO.jl` | Shared SXM (Nanonis) I/O: `SXMImage`/`read_sxm` + preprocessing helpers |
| `GaussianFit1D.jl` | 1D slide profile fitting (diagnostic only, off by default) |
| `GaussianFit2D.jl` | 2D chain model: circular + elliptical Gaussian lobes |
| `STMMolecularFit.jl` | Orchestration: SXM I/O, slide extraction, selectors, batch summaries |
