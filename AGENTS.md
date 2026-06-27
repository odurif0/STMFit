# AGENTS.md — guide for AI agents working on STMFit

This file helps a new AI agent (Claude, GPT, etc.) resume work on this project
efficiently. Read it first, then the docs it points to.

## What this project does

STMFit analyzes STM (Scanning Tunneling Microscopy) images of molecular chains
— primarily chitosan on Cu(100) — by fitting a chain-of-Gaussians model to count
the number of monomer units (lobes) per chain. The selection of N (the lobe
count) is **label-free**: it does not use an expected N or benchmark labels.

**Two distinct regimes:**
- **Benchmark (6mer chitosan, 240817):** 39/39 primary files give the correct
  N=6 (100%), reproducible across runs. This is the *validation* set — labels
  exist and the pipeline is graded against them (labels stay outside
  fitting/selection).
- **Application (10–20mer chitosan):** 25/25 files processed (N_selected 5–16).
  **No ground-truth labels** — this is a *real application*, not a benchmark.
  Visual validation is the arbiter here. The pipeline ran successfully and
  produces internally consistent results, but the counts are not "validated" in
  the benchmark sense.

**Unit assignment (GlcNAc/GlcN per lobe):** Beyond counting N, the pipeline can
assign each fitted lobe a type (0 = GlcN, 1 = GlcNAc) to produce a
deacetylation map per chain. This is a **work in progress** (Phases 0–2a
implemented as diagnostics; robust label-free assignment not solved). The same
label-free rule applies: the ground-truth sequence is used for grading only,
never in the fit. See `docs/src/unit_assignment.md`.

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
6. **`docs/src/unit_assignment.md`** — the GlcNAc/GlcN per-lobe assignment
   pipeline (Phases 0–2a implemented as diagnostics). Read this if working on
   unit assignment.

## Key conventions

- **Selection is label-free.** Never introduce an expected N, target_N, or
  benchmark label into the **fitting or selection** path. Using labels for
  **external evaluation/grading only** (counting how many files match the
  expected N) is fine and expected. The guard rules (`_refined_selection`,
  `_select_primary` in `selectors.jl`) must stay generic. Tuning a parameter
  against a benchmark label and presenting it as objective is explicitly
  forbidden (see journal entries on 043).
- **Unit assignment has no composition prior.** Do not assume the number of
  GlcNAc/GlcN units in a chain, even for the 6mer benchmark. Ground-truth
  sequences and composition counts are external grading/diagnostic information
  only. Rules like "top-k lobes are GlcNAc" are not valid label-free assignment.
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

# Unit assignment (GlcNAc/GlcN per lobe) — see docs/src/unit_assignment.md
# and docs/src/qe_stm_molds.md for the full command reference.
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
- **The batch is reproducible run-to-run** (verified: 3 consecutive runs on the
  same machine give identical N_selected on all 48 files). Divergences between a
  past recorded number and a fresh run indicate code changes between the two,
  not run-to-run noise. (NLopt `GN_DIRECT_L` with a fixed `maxtime`/`maxiter`
  budget is deterministic enough in practice on a given machine.)
- **`max_overlap`** (default 0.60) can block high-N fits on dense chains. It's a
  physical prior (Gaussian pair overlap floor), not arbitrary — but verify it
  isn't rejecting good fits if N looks too low on a new molecule.
- **Auto-calibration** (`measure_calibration.jl`) is a bootstrap, not a
  replacement for visual validation. It under-detects on ~4% of files (e.g.
  251206_013) because parameters are coupled. Always spot-check N_selected
  against the visible structure on a new molecule.

## Benchmark vs application status

- **Benchmark (6mer chitosan, 240817):**
  - `N_selected`: **39/39** primary benchmark exact (N=6).
  - 4/4 clean_target files correct.
  - Reproducible across 3 consecutive runs (0 files change N_selected).
  - Selection threshold robust on [0.03, 0.06] (0 pivot files).
  - **Unit assignment (0/1)**: Phases 0–2a implemented as diagnostics.
    Gaussian/local/patch features are bimodal, split-width Gaussians improve
    fixed-N GCV on 36/39 primary files, and manual geometric connected molds
    validate technically, but current label-free assignment does not recover the
    withheld sequence robustly. The refined raw geometric mold reaches 67.9%
    per-lobe and 0/39 exact; non-chemical shape/asymmetry remains the issue.
    `import_stm_mold_maps.jl` is ready for real DFT-STM/LDOS maps in the aligned
    `(t,u)` frame.
- **Application (10–20mer chitosan):**
  - 25/25 files processed (N_selected 5–16).
  - **No ground-truth labels** — this is a real application, not a benchmark.
  - Control point: 260220_083 → N=9 (manual cross-check, not a benchmark label).

## Documentation discipline

**This project's long-term value is in its documentation, not just its code.**
Keep it current as you work — an undocumented change is a change that didn't
happen for the next agent.

**Mandatory updates after any non-trivial change:**
1. **`docs/src/journal.md`** — add a dated entry for every experiment, decision,
   bug fix, or parameter change. Include *why* (not just *what*). Even failed
   approaches must be recorded so they aren't retried.
2. **Benchmark numbers** (39/39, 25/25, etc.) — if a change affects the results,
   re-run the batch and update every doc that cites the old number (AGENTS.md,
   README.md, chitosan_runbook.md, selection.md).
3. **Open Questions** (journal.md §Open Questions) — resolve, defer, or add as
   the work progresses. Do not let this section go stale.
4. **Config docs** (`docs/src/config.md`, `calibration.md`) — if you add or
   rename a parameter, update the reference the same commit.
5. **AGENTS.md itself** — if the architecture, conventions, or gotchas change,
   update this file.

**Rule of thumb:** if a new agent would give a wrong answer because your change
isn't documented, the documentation is broken. Fix it before committing.

## What NOT to do

- Do not tune a parameter against a benchmark label and call it objective.
- Do not introduce the unit sequence (GlcNAc/GlcN labels) into the fitting or
  selection path — grading/external-evaluation only, same rule as for N.
- Do not change `n_eff` (it's undefined in the fit window; the heuristic `n÷9`
  is a placeholder that only affects BIC/AICc diagnostics, not GCV/N_selected).
- Do not re-enable the 1D fit in the selection path (it over-counts).
- Do not hand-edit `Manifest.toml` — regenerate via `Pkg.resolve()` / `Pkg.instantiate()`.
- Do not commit `results/` artifacts or sensitivity test configs (they're in
  `.gitignore`).
