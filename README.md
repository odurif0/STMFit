# STMFit monorepo

This repository groups the STM molecular fitting stack while keeping separate packages:

- `packages/STMFitCore.jl`: shared physical constraints and scoring helpers.
- `packages/GaussianFit1D.jl`: 1D Gaussian fitting engine.
- `packages/GaussianFit2D.jl`: 2D Gaussian chain fitting engine.
- `packages/STMMolecularFit.jl`: STM/SXM integration, slide extraction, comparisons, batches.
- `packages/STMMolecularFitGUI.jl`: GUI layer.

The root Julia environment uses local package paths via `[sources]`.
