# STMFit — STM Molecular Chain Fitting

Automated pipeline for detecting and fitting 2D Gaussian chain models
to STM images of molecular chains. Historical saved counting runs give 39/39
exact for the robust-AICc guard on the primary benchmark and 129/145 exact
(143/145 within one lobe) for the current support-midpoint hybrid on the expanded
benchmark. Fresh raw-to-prediction reproduction is in progress on September 24;
these saved grades do not establish the current corrected pipeline's score.
Inference reads no labels, but historical calibration used known-count grades;
see [Calibration](calibration.md) for the provenance limit.

The pipeline has also been **applied** to unknown 10–20mer chains. Without
external labels, processing and visual QC do not validate chemical assignment.
Other molecules require their own physical calibration and validation; the
auto-calibration script is a bootstrap, not evidence of generalization.

## Quick Start

```bash
# Single-file diagnostic
STMFIT_DATA_DIR=/path/to/data julia --project=. test/inspect_one_file.jl 240817_004.sxm

# Full batch (default: chitosan.toml, support-midpoint hybrid)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml

# Auto-calibrate for a new molecule from one clean scan
julia --project=. test/measure_calibration.jl path/to/clean_scan.sxm

# Raw GCV baseline (no guard) for comparison
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml --selection-policy gcv

# Original 240817 primary-benchmark validation policy (39/39 exact)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 39 \
  --config config/chitosan.toml --selection-policy gcv_with_robust_aicc_guard

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

The default `config/chitosan.toml` is the historical hand-tuned reference.
`config/chitosan_auto.toml` is an auto-derived comparison, not independent
validation or a certificate of label-free calibration.

For the chitosan reference set, `benchmarks/chitosan_240817.toml` records
evaluation-only quality classes. It is **not** used by fitting code and must not
become a selection prior.

## Unit Assignment (GlcNAc/GlcN)

The pipeline can assign each fitted lobe a monomer type to produce a
deacetylation map per chain. This is a **work in progress**. The 6mer 0/1/?
benchmark uses the same 145 files as the counting benchmark; the control sequence
`NKNNKN` is external grading information only and must not enter the label-free
method. See
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
