# AGENTS.md — guide for AI agents working on STMFit

This file helps a new AI agent (Claude, GPT, etc.) resume work on this project
efficiently. Read it first, then the docs it points to.

## What this project does

STMFit analyzes STM (Scanning Tunneling Microscopy) images of molecular chains
— primarily chitosan on Cu(100) — by fitting a chain-of-Gaussians model to count
the number of monomer units (lobes) per chain. The selection of N (the lobe
count) is **label-free**: it does not use an expected N or benchmark labels.

## Where to look first (read order)

1. **`docs/src/journal.md`** — the dated decision log. This is the project's
   memory: what was tried, what worked, what failed, and *why*. Start here to
   understand the current state and avoid re-treading dead ends.
2. **`docs/src/pipeline.md`** — the data flow and component roles (5 min read).
3. **`docs/src/selection.md`** — the selection rule (GCV + robust-AICc guard +
   up-when-ambiguous). This is the scientific heart.
4. **`docs/src/calibration.md`** — parameter objectivation (which are measured,
   which are free) and why GCV is the canonical criterion (not BIC/AICc).
5. **`docs/src/config.md`** — every parameter, its role, and how it's configured.

## Key conventions

- **Selection is label-free.** Never introduce an expected N, target_N, or
  benchmark label into the selection path. The guard rules (`_refined_selection`,
  `_select_primary` in `selectors.jl`) must stay generic. Tuning against a label
  is explicitly forbidden (see journal entries on 043).
- **GCV is canonical; BIC/AICc are diagnostics only.** The STM residual field is
  so strongly spatially correlated (range 17–100 px, larger than the ~10-px fit
  window) that `n_eff` is effectively undefined. BIC/AICc assume iid — their
  absolute values are not reliable. GCV (valid under correlation) drives
  `N_selected`. See `docs/src/calibration.md`.
- **The 1D fit is off by default.** It never enters `N_selected` (diagnostic
  only). Use `--no-skip-1d` to re-enable it for cross-checking.
- **`config/*.toml` drives everything.** System-specific parameters (σ, spacing,
  support) live in TOML files, not in code defaults. Code defaults are fallbacks.
- **Configs have three sections**: `[model]` (physical), `[selection]`
  (thresholds), `[preprocessing]` (SXM channel/flatten).

## Commands you'll use

```bash
# Single-file inspection (fast, no batch)
julia --project=. test/inspect_one_file.jl <file.sxm>

# Full batch (production)
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
    --config config/chitosan.toml

# Auto-calibrate from one clean scan (for a new molecule)
julia --project=. test/measure_calibration.jl <clean_scan.sxm>

# HPC batch
./hpc/launch_remote.sh --watch    # sync → submit → wait → merge → fetch

# Unit tests for packages
julia --project=packages/STMSXMIO.jl packages/STMSXMIO.jl/test/runtests.jl
julia --project=packages/STMFitCore.jl packages/STMFitCore.jl/test/runtests.jl
```

`batch_full.jl` flags: `--config`, `--data-dir`, `--outdir`, `--chunk i/n`,
`--exclude-from <file>`, `--selection-policy`, `--gcv-ambiguity-rel-threshold`,
`--robust-guard-nu`, `--skip-1d` (default) / `--no-skip-1d`.

## Architecture (5 packages + driver)

```
STMFitCore  ←  STMSXMIO  ←  GaussianFit1D
                       ←  GaussianFit2D  ←  STMMolecularFit (selectors.jl)
                       ←  STMMolecularFit ←  STMMolecularFitGUI
test/batch_full.jl (driver, not a package) orchestrates the batch.
```

- `STMSXMIO.jl` owns the SXM types + reader + shared preprocessing helpers. Both
  engines `using STMSXMIO`. Do **not** redefine SXM types in GF2 or MF.
- `GaussianFit2D.jl/src/core.jl` is the 2D fit engine (~1800 lines, the core).
- `STMMolecularFit.jl/src/selectors.jl` contains the selection logic (the guard,
  up-when-ambiguous rule, ~680 lines).
- `test/batch_full.jl` (~1300 lines) is the production batch driver. It reads the
  TOML, builds configs, runs the sweep, applies selection, writes the summary TSV.

## Known gotchas

- **HPC quota**: user `oldu` has `GrpCPUs=8` on Raven. Use ≤ 8 CPUs total
  (e.g. 2 chunks × 4 CPUs, or 4 chunks × 1 CPU). The cluster is often congested.
- **Wall time**: 10–20mer files with long chains (N up to 25) are slow. Use
  `--time=08:00:00` or more on HPC. The `intelligent_sweep` early-stops, so for
  diagnostic exhaustive sweeps set `intelligent_sweep=false`.
- **NLopt `GN_DIRECT_L`** is deterministic in practice (reproducible run-to-run
  on a given machine). Batch divergences between dates indicate code changes, not
  run-to-run noise.
- **`max_overlap`** (default 0.60) can block high-N fits on dense chains. It's a
  physical prior (Gaussian pair overlap floor), not arbitrary — but verify it
  isn't rejecting good fits if N looks too low on a new molecule.
- **Auto-calibration** (`measure_calibration.jl`) is a bootstrap, not a
  replacement for visual validation. It under-detects on ~4% of files (e.g.
  251206_013) because parameters are coupled. Always spot-check N_selected
  against the visible structure on a new molecule.

## Benchmark validation status (chitosan 6mer)

- `N_selected`: **39/39** primary benchmark exact (N=6).
- 4/4 clean_target files correct.
- Reproducible across 3 consecutive runs (0 files change N_selected).
- Selection threshold robust on [0.03, 0.06] (0 pivot files).
- 10–20mer: 25/25 files processed, N_selected 5–16, control point 260220_083→9.

## What NOT to do

- Do not tune a parameter against a benchmark label and call it objective.
- Do not change `n_eff` (it's undefined in the fit window; the heuristic `n÷9`
  is a placeholder that only affects BIC/AICc diagnostics, not GCV/N_selected).
- Do not re-enable the 1D fit in the selection path (it over-counts).
- Do not hand-edit `Manifest.toml` — regenerate via `Pkg.resolve()` / `Pkg.instantiate()`.
- Do not commit `results/` artifacts or sensitivity test configs (they're in
  `.gitignore`).
