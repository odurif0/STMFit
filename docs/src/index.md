# STMFit — STM Molecular Chain Fitting

Automated pipeline for detecting and fitting 2D Gaussian chain models
to STM images of molecular chains (chitosan on Cu(100)).

## Quick Start

```bash
# Single file diagnostic
julia --project=. test/inspect_one_file.jl data/molecule.sxm

# Batch processing (uses config/chitosan.toml by default)
julia --project=. test/batch_full.jl

# Batch with custom calibration (different molecule/instrument)
julia --project=. test/batch_full.jl --config config/my_system.toml

# Summarize results
julia --project=. test/summarize.jl results/best_plots/
```

## Calibration

System-specific parameters are stored in `config/*.toml` files.
The default is `config/chitosan.toml` (chitosan on Cu(100), Nanomics STM).

To calibrate for a new system:
1. Copy `config/chitosan.toml` → `config/my_system.toml`
2. Adjust the noise-based support parameters (`support_noise_k`,
   `support_padding_nm`), `spacing_min_nm`, `spacing_max_nm`, sigma bounds,
   and the model-selection fields.
3. Run batch with `--config config/my_system.toml`

Track the selected `N_ell`, score curves, residuals, support length, and
absolute score values (not just the winning N) to validate calibration.

For the chitosan reference set, `benchmarks/chitosan_240817.toml` records
evaluation-only quality classes (`clean_target`, `poor_quality`, `excluded`).
It is not used by the fitting code and must not become a selection prior.

## Research Journal

All experimental paths — successful and failed — are documented in the
[Research Journal](journal.md). **Update this journal** whenever you:
- Test a new approach (even if it fails)
- Change the pipeline or model selection logic
- Discover a bug or convergence issue
- Add or remove parameters

The journal is the project's memory. Future you (and future colleagues)
will thank present you.

| Component | Role |
|-----------|------|
| `STMFitCore.jl` | Shared utilities: κ penalty, spacing constraints |
| `GaussianFit1D.jl` | 1D slide profile fitting and QC comparison |
| `GaussianFit2D.jl` | 2D chain model: circular + elliptical Gaussian lobes |
| `STMMolecularFit.jl` | Orchestration: SXM I/O, slide extraction, batch summaries |

## Pipeline (v6)

```
SXM image ─┬→ slide profile (1D) → peak fitting → QC comparison
           │
           └→ 2D circular sweep (N adaptive, independent deterministic init)
              │
         circ→ell LsqFit refinement (∀N, warm-started from circ)
              │
         best N by configured criterion (default GCV)
              │
         outputs: N_ell, N_circ, N_eff, support/residual QC
```
