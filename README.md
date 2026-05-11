# STMFit monorepo

This repository groups the STM molecular fitting stack while keeping separate packages:

- `packages/STMFitCore.jl`: shared physical constraints and scoring helpers.
- `packages/GaussianFit1D.jl`: 1D Gaussian fitting engine.
- `packages/GaussianFit2D.jl`: 2D Gaussian chain fitting engine.
- `packages/STMMolecularFit.jl`: STM/SXM integration, slide extraction, comparisons, batches.
- `packages/STMMolecularFitGUI.jl`: GUI layer.

The root Julia environment uses local package paths via `[sources]`.

## Tests & workflows

Unit tests:
```bash
julia --project=. test/runtests.jl
```

Integration workflows (`test/scripts/`):
- `compare_circular_elliptical.jl` – single-file 2D ell vs circ vs 1D comparison
- `compare_slide_modes.jl` – compare slide extraction modes
- `batch_triage.jl` – batch 2D chain sweep TSV
- `centralize_best_plots.jl` – full 1D/2D best-plots batch
- `summarize_enriched.jl` – stats from a summary TSV

Helper:
- `common.jl` – shared paths, configs, chunk parsing

Example:
```bash
julia --project=. test/scripts/centralize_best_plots.jl 48 --chunk 1/4 \
  > results/logs/batch.log 2>&1 &
```

Runtime outputs go under `results/`.
