# STMFit agent notes

STMFit is a Julia repo for label-free lobe counting and diagnostic
GlcN/GlcNAc assignment in STM images. Read `docs/src/journal.md` before changing
scientific behavior; it records successful and failed experiments that are not
recoverable from the current code alone.

## Scientific invariants

- Fitting, `N_selected`, unit assignment, calibration, thresholds, and
  abstention must not use expected `N`, `NKNNKN`, class counts, or any benchmark
  label. Labels belong only in external grading/report scripts.
- Unit assignment has no composition prior. “Choose the top k lobes” is not a
  valid label-free method.
- GCV drives base model selection. BIC/AICc are diagnostics or guards because
  spatially correlated residuals make absolute iid information criteria
  unreliable. Do not reinterpret or retune the `n_eff = n ÷ 9` placeholder.
- The 1D fit is diagnostic-only and skipped by default. It must not feed
  `N_selected`; `--no-skip-1d` only enables a cross-check.
- Physical and selection parameters belong in `config/*.toml`, not hidden code
  defaults. Configs have `[model]`, `[selection]`, and `[preprocessing]`.
- Keep benchmark and application claims separate: the 6mer sets have external
  labels; 10–20mer data do not. “Processed” or visually plausible is not
  benchmark validation.
- Relaxed GlcN and GlcNAc DFT cubes now pass the documented production gate at
  the common `5e-5 Ry` acceptance criterion; GlcNAc was validated by byte-identical
  plain/TF LDOS cubes. Read `docs/src/dft_calculation_note.md` before changing the
  criterion or replacing the still-preliminary joint-proxy registry sources.

## Real boundaries and entrypoints

- The root project composes five local packages through `Project.toml [sources]`:
  `STMFitCore`, `STMSXMIO`, `GaussianFit1D`, `GaussianFit2D`, and
  `STMMolecularFit`. `STMMolecularFitGUI` is a separate package and its launcher
  still contains pre-monorepo assumptions; do not treat it as production entry.
- `STMSXMIO` owns `SXMImage`, SXM parsing, channel alignment, and shared
  preprocessing. Do not redefine those types in either fit engine.
- `GaussianFit1D` depends on `STMFitCore`, not `STMSXMIO`.
  `GaussianFit2D` depends on Core, SXM I/O, and the 1D package.
  `STMMolecularFit` orchestrates all four.
- `packages/GaussianFit2D.jl/src/core.jl` is the fit engine;
  `packages/STMMolecularFit.jl/src/selectors.jl` owns selection;
  `test/batch_full.jl` is the production driver and imports some selector
  internals deliberately.
- Root `test/*.jl` files are mostly standalone research/workflow programs, not a
  conventional single test suite. Use their `--help` or header usage before
  assuming arguments.

## Setup and focused verification

CI uses Julia 1.11 and explicitly develops local packages before instantiate.
Use the same bootstrap on a clean depot; do not hand-edit `Manifest.toml`:

```bash
julia --project=. -e '
using Pkg
for p in ["STMFitCore","STMSXMIO","GaussianFit1D","GaussianFit2D","STMMolecularFit"]
    Pkg.develop(PackageSpec(path="packages/$p.jl"))
end
Pkg.instantiate(); Pkg.precompile()'
```

There is no root aggregator for the package tests. Run package suites under the
bootstrapped root environment so local, unregistered dependencies resolve:

```bash
julia --project=. packages/STMFitCore.jl/test/runtests.jl
julia --project=. packages/STMSXMIO.jl/test/runtests.jl
julia --project=. packages/GaussianFit2D.jl/test/runtests.jl
julia --project=. test/joint_proxy/runtests.jl   # joint-proxy workflow suite
```

Build docs from the root environment; `docs/make.jl` activates the root project:

```bash
GKSwstype=100 julia --project=. docs/make.jl
```

## Data and production commands

Raw `.sxm` files are untracked. Point scripts to them with `STMFIT_DATA_DIR` or
`--data-dir`; see `data/README.md` for the expected layout.

```bash
# Fast, focused inspection before a batch
julia --project=. test/inspect_one_file.jl <file.sxm>

# Production batch; the first positional argument limits file count
STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. test/batch_full.jl 48 \
  --config config/chitosan.toml

# Bootstrap a new molecule from one clean scan, then visually check it
julia --project=. test/measure_calibration.jl <clean_scan.sxm>

# Unknown-sequence unit assignment and integrity checks
julia --project=. test/run_unknown_unit_assignment.jl --help
julia --project=. test/validate_unit_predictions.jl --help
julia --project=. test/summarize_unknown_unit_qc.jl --help
```

`measure_calibration.jl` is a bootstrap, not ground truth; parameters are
coupled and known scans can be under-detected. Spot-check visible structure.
If high-N fits disappear, inspect the physical `max_overlap` constraint before
changing selection logic.

## HPC and generated state

- Configure `hpc/remote.env` from `hpc/remote.env.example`; it is personal and
  ignored. Always run `./hpc/launch_remote.sh --dry-run` before `--watch`.
- The launcher order is sync -> instantiate on the login node -> Slurm array ->
  merge chunk TSVs -> fetch. Compute nodes have no internet; never run the STM
  batch or QE compute on login nodes.
- Raven account `oldu` has an observed 8-CPU group quota. The generic launcher
  defaults (`4 chunks x 4 CPUs`) exceed it; configure at most 8 CPUs total.
  `batch_full.jl` caps useful Julia threads at four per task.
- Long 10–20mer chains may need at least 8 hours. `intelligent_sweep` early
  stops; disable it only for intentionally exhaustive diagnostics.
- Run heavyweight diagnostic batches, multi-file fits, and exhaustive whole-ROI
  searches on Viper rather than locally. Prepare a dry-run first and shard them;
  local execution is for focused single-file checks only.
- QE molds use `hpc/launch_qe_molds_remote.sh`, not the STM array launcher, and
  should stay sequential under the one-node Raven QOS limit.
- Do not commit `results/`, raw `data/`, `qe/`, `docs/build/`, generated
  sensitivity configs, mold/map TSVs, or generated QE reports; `.gitignore`
  enumerates them. `/ptmp` is temporary, so fetch HPC results locally.

## Documentation contract

- Add a dated `docs/src/journal.md` entry for every experiment, scientific
  decision, failed approach, bug fix, or parameter change, including why.
- If behavior changes benchmark results, rerun the appropriate grade and update
  every cited headline (`README.md`, runbook, selection/unit-assignment docs),
  rather than copying an old number into this file.
- Keep the journal’s Open Questions current. Parameter additions/renames also
  require `docs/src/config.md` and `docs/src/calibration.md` updates.
- Package READMEs contain legacy standalone paths and 1D/BIC-era descriptions.
  For current methodology and commands, trust root scripts/config, root README,
  and `docs/src/`.
