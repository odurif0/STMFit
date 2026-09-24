# STMFit — STM Molecular Chain Fitting

Automated pipeline for detecting and fitting 2D Gaussian chain models
to STM images of molecular chains. Two fresh raw-to-prediction repetitions on
September 24 agree: **123/145 exact counts** (143/145 within one lobe),
**679/870 correct units, 29/145 exact chains, 848/870 coverage**. The saved-count
129/145 and 694/43/855 results are not reproduced end to end. Historical
677/36/854 is not surpassed overall; no new champion is promoted.
Inference reads no labels, but historical calibration used known-count grades;
see [Calibration](calibration.md) for the provenance limit.

The pipeline has also been **applied** to unknown 10–20mer chains. Without
external labels, processing and visual QC do not validate chemical assignment.
Other molecules require their own physical calibration and validation; the
measurement diagnostic does not automatically establish a physical calibration.

## Quick Start

```bash
# Single-file diagnostic
STMFIT_DATA_DIR=/path/to/data julia --project=. test/inspect_one_file.jl 240817_004.sxm

# Full batch (default: chitosan.toml, support-midpoint hybrid)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml

# Audit apparent widths/spacings and missing measurements, without fitting
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

`test/measure_calibration.jl` reports apparent widths/spacings in both views and
exposes the old bootstrap's fallback values. It no longer derives physical fit
bounds or writes a production TOML automatically. See [**Calibration**](calibration.md)
for the measurement limits and why GCV remains the canonical criterion.
The completed raw-view audit finds legacy width/spacing fallbacks on 245/292
and 261/292 views; the new diagnostic is not a promoted physical calibration.

The default `config/chitosan.toml` is the historical hand-tuned reference.
`config/chitosan_auto.toml` is a historical bootstrap output, including values
matching the old fallbacks, not independent validation or label-free calibration.

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
