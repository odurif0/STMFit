# STMFit agent notes

STMFit is a Julia repo for label-free reconstruction of chitosan chains from raw
STM images: unit counting, positions, and GlcN/GlcNAc assignment with explicit
uncertainty. Read `docs/src/journal.md` (current state, open questions,
rejected approaches) before changing scientific behavior. The full pre-2026-09-29
history is at the tag `archive/pre-cleanup-20260929`.

## Scientific invariants

- Fitting, `N_selected`, molecule consensus, unit assignment, fusion,
  calibration, thresholds and abstention must not use expected `N`, `NKNNKN`,
  class counts, benchmark membership or any label. Labels belong only in
  external grading/report scripts (`test/grade_*.jl`, `test/report_*.jl`).
- Unit assignment has no composition prior. "Choose the top k lobes" is not a
  valid label-free method. Cluster naming is physical (brighter = GlcNAc).
- GCV drives base model selection. BIC/AICc are diagnostics or guards because
  spatially correlated residuals make absolute iid information criteria
  unreliable. Do not reinterpret or retune the `n_eff = n ÷ 9` placeholder.
- The 1D fit is diagnostic-only and skipped by default. It must not feed
  `N_selected`; `--no-skip-1d` only enables a cross-check.
- Physical, selection and assignment parameters belong in `config/*.toml`, not
  hidden code defaults. Configs have `[model]`, `[selection]`, `[preprocessing]`.
  Promoted configs are frozen; a change is a new declared experiment.
- Keep benchmark and application claims separate: the 6-mer set has external
  labels; 10–20mer data do not. "Processed" or visually plausible is not
  validation. NKNNKN confounds GlcNAc with edge adjacency, so benchmark
  accuracy is not chemical validation.
- The DFT molds (`templates/chitosan_cc_molds_native_v1.tsv`) are empirical
  templates from the accepted cubes (5e-5 Ry criterion). They keep a legacy
  cube addressing and read a signed LDOS quantity. Read
  `docs/src/dft_calculation_note.md` before changing cubes, criterion or builder;
  new molds get a new versioned file.

## Structure and entrypoints

- The root project composes five local packages through `Project.toml [sources]`:
  `STMFitCore`, `STMSXMIO`, `GaussianFit1D`, `GaussianFit2D`, `STMMolecularFit`.
  `STMMolecularFitGUI` is separate and unmaintained; not a production entry.
- `STMSXMIO` owns `SXMImage`, SXM parsing, channel alignment and shared
  preprocessing. Do not redefine those types in either fit engine.
  `GaussianFit1D` depends on `STMFitCore` only; `GaussianFit2D` on Core, SXM I/O
  and 1D; `STMMolecularFit` orchestrates all four.
- `packages/GaussianFit2D.jl/src/core.jl` is the fit engine;
  `packages/STMMolecularFit.jl/src/selectors.jl` owns selection;
  `test/batch_full.jl` is the counting driver and imports selector internals.
- Promoted method, raw scans to predictions:
  `test/run_molecule_consensus_chitosan.jl` (Slurm: `hpc/run_molecule_consensus.sbatch`).
  It runs `batch_full.jl` → `extract_lobe_features.jl` →
  `build_molecule_consensus.jl` → `run_reconstructed_chitosan.jl` (assignment) →
  `build_molecule_fusion.jl`. Shared code in `test/lib/`.
- Verification and external grading: `test/grade_consensus_run.jl`.
- `test/*.jl` scripts are standalone command-line programs. Use their
  `--help` or header usage before assuming arguments.

## Setup and verification

Julia 1.13 (explicit user decision, 2026-09-16), including tests and
subprocesses; local `julia` resolves to 1.13 via juliaup. The root
`Manifest.toml` is tracked and pins the environment for releases; change it
only through Pkg, never by hand. Bootstrap a clean depot with:

```bash
julia --project=. -e '
using Pkg
for p in ["STMFitCore","STMSXMIO","GaussianFit1D","GaussianFit2D","STMMolecularFit"]
    Pkg.develop(PackageSpec(path="packages/$p.jl"))
end
Pkg.instantiate(); Pkg.precompile()'
```

Tests, run under the root environment:

```bash
julia --project=. packages/STMFitCore.jl/test/runtests.jl
julia --project=. packages/STMSXMIO.jl/test/runtests.jl
julia --project=. packages/GaussianFit2D.jl/test/runtests.jl
for t in test/test_*.jl test/molecule_consensus/*.jl; do julia --project=. "$t"; done
julia --project=. test/test_reconstructed_pipeline.jl --e2e
```

Refactors of kept scripts must be verified byte-for-byte against stored
outputs. Changes to scientific behavior need a fresh raw run and an external grade.

Docs: `GKSwstype=100 julia --project=. docs/make.jl --build-only`.

## Data and commands

Raw `.sxm` files are untracked (`data/README.md`). Pass `--data-dir` (or
`STMFIT_DATA_DIR` for `batch_full.jl`).

```bash
julia --project=. test/inspect_one_file.jl <file.sxm>            # one scan
julia -t 4 --project=. test/batch_full.jl 48 --config config/chitosan.toml --data-dir DIR
julia --project=. test/run_molecule_consensus_chitosan.jl --help  # promoted method
julia --project=. test/grade_consensus_run.jl --run RUN --outdir NEW  # verify + grade
julia --project=. test/measure_calibration.jl <clean_scan.sxm>    # audit, no config written
```

If high-N fits disappear, inspect the physical `max_overlap` constraint before
changing selection logic.

## HPC and generated state

- HPC is deliberately boring: sync committed code (plus the ignored Manifest),
  instantiate on the login node, dry-run, submit one `sbatch`, poll `squeue`,
  fetch outputs locally. No agent orchestration or state on the cluster that
  outlives a job. Never compute on login nodes; compute nodes have no internet.
- Viper sets `SBATCH_EXPORT=NONE`: pass variables with
  `sbatch --export=NONE,VAR=...`; shell variables before `sbatch` are not seen.
- Run cohort fits, multi-file batches and exhaustive searches on Viper, not
  locally. Local runs are for single-file checks and tests.
- Raven account `oldu` has an 8-CPU group quota; the `launch_remote.sh` default
  (4 × 4 CPUs) exceeds it. `batch_full.jl` uses at most 4 threads per task.
- Long 10–20mer counting may need 8 h or more. QE jobs stay sequential under the
  one-node QOS limit.
- Do not commit `results/`, raw `data/`, `qe/`, `docs/build/`, `hpc/remote.env`
  or generated outputs. `/ptmp` is temporary: fetch HPC results.

## Documentation contract

- Add a dated `docs/src/journal.md` entry for every experiment, scientific
  decision, failed approach, bug fix or parameter change, with the reason.
  Keep its Current state and Open questions current.
- If behavior changes benchmark results, rerun the grade and update every
  cited headline (`README.md`, `docs/src/index.md`, runbook, selection and
  unit-assignment pages).
- Parameter additions or renames also require `docs/src/config.md` and
  `docs/src/calibration.md` updates.
- Package READMEs contain legacy standalone descriptions. Trust the root
  README, `docs/src/` and the scripts.

## Process guardrails (2026-09-16, after the T11–T13 apparatus incident)

- Every task starts from a stated scientific deliverable and its definition of
  done. Process artifacts (protocols, ledgers, receipts, review gates, staged
  snapshots) are not deliverables.
- Do not create authorization layers, meta-review cycles, approval cascades or
  "runtime/activation" frameworks. Verification is tests plus one human review.
- New method campaigns need a human-approved scope and a time box. Current
  priority: the unknown 10–20mer application (see the journal's open questions).
- Commit early; never leave validated work uncommitted. Git identity is set
  repo-local (Olivier Durif <o.durif@fkf.mpg.de>). GitHub push goes through the
  local machine (the clusters have no credentials).
- The retired `.omo/` agent state and its September archive were deleted on
  2026-10-04. Never write agent state into the repo outside git.
