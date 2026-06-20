# STMFit — STM Molecular Chain Fitting

Automated pipeline for detecting and fitting 2D Gaussian chain models
to STM images of molecular chains. **Validated** on the 6mer chitosan/Cu(100)
benchmark (39/39 exact); **applied** to the 10–20mer production system (no
ground-truth labels — visual validation is the arbiter). Generalizable to other
chain-like molecules on the same STM via auto-calibration
(see [Calibration](calibration.md)).

## Quick Start

```bash
# Single file diagnostic
STMFIT_DATA_DIR=/path/to/data julia --project=. test/inspect_one_file.jl 240817_004.sxm

# Batch processing (default config: chitosan.toml, policy: gcv_with_robust_aicc_guard)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml

# Auto-calibrate for a new molecule from one clean scan
julia --project=. test/measure_calibration.jl path/to/clean_scan.sxm

# Explicit raw GCV/N_eff baseline override
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy gcv

# Diagnostic spatial blocked CV selector (not default)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy spatial_blocked_cv \
  --cv-folds 3

# Diagnostic support-marginalized GCV selector (not default)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy support_marginalized_gcv

# Conservative support-marginalized one-lobe guard (not default)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy support_marginalized_gcv_guard

# File-adaptive slope-heuristic MDL selector (diagnostic)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy slope_heuristic_mdl

# Support-perturbation stability selector (diagnostic)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy stability_selection

# Local lobe-resolvability guard (diagnostic/inconclusive)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy local_lobe_evidence

# Approximate Laplace-evidence guard (diagnostic)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy laplace_evidence_guard

# Fwd/bwd direction-consensus selector (diagnostic)
julia --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml \
  --selection-policy fwd_bwd_consensus

# Synthetic known-N selector validation
julia --project=. test/synthetic_known_n_validation.jl \
  --cases 50 \
  --seed 1234 \
  --noise-scale 1.0 \
  --mode circ_ell \
  --out results/synthetic_known_n/summary_50_circ_ell.tsv

# Aggregate phase-2 synthetic selector validation
julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50_circ_ell.tsv \
  --out results/synthetic_known_n/aggregate_50_circ_ell.tsv

# Fast circular-only synthetic validation
julia --project=. test/synthetic_known_n_validation.jl \
  --cases 50 \
  --seed 1234 \
  --noise-scale 1.0 \
  --mode circular \
  --out results/synthetic_known_n/summary_50.tsv

# Aggregate synthetic selector validation
julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50.tsv \
  --out results/synthetic_known_n/aggregate_50.tsv

# Aggregate several synthetic seeds/noise levels
julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50.tsv \
  results/synthetic_known_n/summary_50_seed2026.tsv \
  results/synthetic_known_n/summary_50_seed4321.tsv \
  results/synthetic_known_n/summary_50_seed1234_noise05.tsv \
  results/synthetic_known_n/summary_50_seed1234_noise15.tsv \
  --out results/synthetic_known_n/aggregate_multiseed_noise.tsv

# Aggregate several phase-2 circ→ell synthetic seeds/noise levels
julia --project=. test/aggregate_synthetic_known_n.jl \
  results/synthetic_known_n/summary_50_circ_ell.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed2026.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed4321.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed1234_noise05.tsv \
  results/synthetic_known_n/summary_50_circ_ell_seed1234_noise15.tsv \
  --out results/synthetic_known_n/aggregate_circ_ell_multiseed_noise.tsv

# Summarize results
julia --project=. test/summarize.jl results/best_plots/
```

## Calibration

System-specific parameters are stored in `config/*.toml` files.
The default is `config/chitosan.toml` (chitosan on Cu(100), Nanonis STM).

To calibrate for a new system, the recommended path is **auto-calibration**:
point `test/measure_calibration.jl` at one clean scan and it derives σ, spacing,
fit width, support length and n_max objectively, emitting a ready-to-use TOML.

```bash
julia --project=. test/measure_calibration.jl path/to/clean_scan.sxm
# → emits <scan>_calibration.toml
```

For the full parameter classification (measured / principled / free), the
effective-sample-size analysis that justifies GCV as the canonical criterion,
and the manual fallback, see **[Calibration](calibration.md)**.

The default `config/chitosan.toml` remains the hand-tuned reference; an
auto-derived equivalent (`config/chitosan_auto.toml`) validates that the
objective method reproduces the benchmark results.

For the chitosan reference set, `benchmarks/chitosan_240817.toml` records
evaluation-only quality classes (`clean_target`, `poor_quality`, `excluded`).
It is not used by the fitting code and must not become a selection prior.

For a reproducible hand-off of the current benchmark-aligned chitosan and
10–20mer workflow, see the [Chitosan Runbook](chitosan_runbook.md).

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
| `STMSXMIO.jl` | Shared SXM (Nanonis) I/O: `SXMImage`/`read_sxm` + preprocessing helpers |
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
               │
          optional experimental robust-AICc guard → N_selected
               │
          optional diagnostic spatial blocked CV → N_selected
               │
          optional support-marginalized GCV diagnostic → N_selected
               │
          optional slope-heuristic MDL diagnostic → N_selected
               │
          optional stability-selection diagnostic → N_selected
               │
          optional local-lobe-evidence diagnostic → N_selected
                │
          optional Laplace-evidence diagnostic → N_selected
                │
          synthetic known-N validation compares selectors externally
```
