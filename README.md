# STMFit monorepo

Packages:
- `packages/STMFitCore.jl` — shared physical constraints and scoring.
- `packages/GaussianFit1D.jl` — 1D Gaussian fitting engine.
- `packages/GaussianFit2D.jl` — 2D Gaussian chain fitting engine.
- `packages/STMMolecularFit.jl` — STM/SXM integration, slide extraction, comparisons.
- `packages/STMMolecularFitGUI.jl` — GUI.

Workflows (`test/`), each standalone:

| Script | Purpose |
|---|---|
| `inspect_one_file.jl <file.sxm>` | Deep 1D vs 2D ell vs 2D circ on a single file |
| `batch_full.jl [N] [--chunk i/n]` | Full 1D + 2D batch: fits, plots, enriched summary |
| `summarize.jl [summary.tsv]` | Print stats from a summary TSV |

```bash
julia --project=. test/inspect_one_file.jl 240817_004.sxm
julia --project=. test/batch_full.jl 48 --chunk 1/4 > results/logs/batch.log 2>&1 &
julia --project=. test/summarize.jl results/best_plots/summary_overlap060_hard.tsv
```

Runtime outputs → `results/`.
