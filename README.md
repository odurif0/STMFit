# STMFit monorepo

Packages:
- `packages/STMFitCore.jl` — shared physical constraints and scoring.
- `packages/GaussianFit1D.jl` — 1D Gaussian fitting engine.
- `packages/GaussianFit2D.jl` — 2D Gaussian chain fitting engine.
- `packages/STMMolecularFit.jl` — STM/SXM integration, slide extraction, comparisons.
- `packages/STMMolecularFitGUI.jl` — GUI.

Unit tests:
```bash
julia --project=. test/runtests.jl
```

Workflows (`test/scripts/`), each standalone:

| Script | Purpose |
|---|---|
| `inspect_one_file.jl <file.sxm>` | Deep 1D vs 2D ell vs 2D circ on a single file |
| `batch_triage_2d.jl <dir>` | Fast 2D-only batch, pre-filter valid files |
| `batch_full.jl [N] [--chunk i/n]` | Full 1D + 2D batch: fits, plots, enriched summary |
| `summarize.jl [summary.tsv]` | Print stats from a summary TSV |

Example:
```bash
julia --project=. test/scripts/batch_full.jl 48 --chunk 1/4 \
  > results/logs/batch.log 2>&1 &
```

Runtime outputs → `results/`.
