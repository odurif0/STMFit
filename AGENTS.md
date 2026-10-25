# STMFit agent notes

STMFit is a Julia repo for label-free reconstruction of chitosan chains from raw
STM images: unit counting, positions, and a two-class unit assignment with
explicit uncertainty (the chemical identity of the classes is not established). Read `docs/src/journal.md` (current state, open questions,
rejected approaches) before changing scientific behavior. The full pre-2026-09-29
history is at the tag `archive/pre-cleanup-20260929`.

## Scientific invariants

- Fitting, `N_selected`, molecule consensus, unit assignment, fusion,
  calibration, thresholds and abstention must not use expected `N`, `NKNNKN`,
  class counts, benchmark membership or any label. Labels belong only in
  external grading/report scripts, all in `scripts/grading/`. Nothing in
  `scripts/pipeline/` or `src/` may read a label file.
- Unit assignment has no composition prior. "Choose the top k lobes" is not a
  valid label-free method. Cluster naming is physical and chemically neutral:
  class 1 = brighter (higher mean raw amplitude), class 0 = dimmer. Method code
  does not call the classes GlcN/GlcNAc; the former rule "brighter = GlcNAc"
  is contradicted by the 6-mer control (2026-10-22 journal entry).
- NKNNKN: N = GlcNAc (N-acetyl), K = GlcN. The control has GlcNAc at
  positions 1, 3, 4, 6 and GlcN at 2, 5; the truth strings `010010` mark the
  GlcN.
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
  validation. The 10–20mer chain lengths are unknown: the experimenters'
  readings (lengths, acetylation) are comparisons, never a truth criterion.
  The reference N/K pattern of the 10–20mer sequence will come from the user.
  In NKNNKN the two GlcN sit next to the chain ends, so the
  bright class coincides with both GlcN and edge adjacency; benchmark accuracy
  is not chemical validation.
- The DFT molds (`dft/molds/chitosan_cc_molds_native_v1.tsv`) are empirical
  templates from the accepted cubes (5e-5 Ry criterion). They keep a legacy
  cube addressing and read a signed LDOS quantity. Read
  `docs/src/dft_calculation_note.md` before changing cubes, criterion or builder;
  new molds get a new versioned file.

## Method development strategy (user, 2026-10-25)

- The 6-mer benchmark (known length and sequence) is where the method is
  optimized and selected. The 10–20mer data are the unknown: the method is
  applied to them unchanged, with no a priori constraint on the result
  (length, composition, pattern). This is what makes the 10–20mer
  measurement trustworthy.
- Freeze what was optimized before applying it, assignment included where
  possible. Re-learning the classes on the 10–20mer cohort is not the same
  method.
- Optimize with held-out sessions or molecules. Prefer correct units, errors
  and per-molecule scores to exact chains, which come in molecule blocks.
  NKNNKN confounds K with edge adjacency and has no interior K, so it cannot
  validate interior assignment.
- The user's 10–20mer N/K pattern, once given, is a final test, never an
  optimization target.
- The optimization dataset decides the result. Its image quality must be
  impeccable, and more images is not better when they are unusable. Select
  scans on image quality judged before any method output (visual review or
  label-free image criteria), never on whether the method succeeds. Apply
  the same quality standard to the 10–20mer data.

## Structure and entrypoints

```
src/STMFit.jl      pipeline library (MoleculeConsensus, MoleculeFusion,
                   ReconstructedUnitAssignment, PatchPreprocessing, CCMoldNative,
                   EmpiricalFisherNative, UnitAssignmentReport, ScriptUtils)
packages/          fit engines and SXM I/O (five local packages)
scripts/pipeline/  the promoted method: label-free
scripts/grading/   external grading and reports: the only code that reads labels
scripts/data/      raw-scan catalog and repeat-scan collection
scripts/calibration/, scripts/tools/   calibration audit, single-scan inspection
test/              tests only (test/runtests.jl runs them all)
config/  benchmarks/  dft/  hpc/  docs/
data/  results/  dft/runs/   untracked
```

- The root project composes five local packages through `Project.toml [sources]`:
  `STMFitCore`, `STMSXMIO`, `GaussianFit1D`, `GaussianFit2D`, `STMMolecularFit`.
- `STMSXMIO` owns `SXMImage`, SXM parsing, channel alignment and shared
  preprocessing. Do not redefine those types in either fit engine.
  `GaussianFit1D` depends on `STMFitCore` only; `GaussianFit2D` on Core, SXM I/O
  and 1D; `STMMolecularFit` orchestrates all four.
- `packages/GaussianFit2D.jl/src/core.jl` is the fit engine;
  `packages/STMMolecularFit.jl/src/selectors.jl` owns selection;
  `scripts/pipeline/batch_full.jl` is the counting driver and imports selector
  internals.
- Promoted method, raw scans to predictions:
  `scripts/pipeline/run_molecule_consensus_chitosan.jl` (Slurm:
  `hpc/run_molecule_consensus.sbatch`). It runs `batch_full.jl` →
  `extract_lobe_features.jl` → `build_molecule_consensus.jl` →
  `run_reconstructed_chitosan.jl` (assignment) → `build_molecule_fusion.jl` →
  `report_position_repeatability.jl` → `plot_final_molecules.jl`. Shared code
  lives in `src/` (`using STMFit.<Module>`).
- Verification and external grading: `scripts/grading/grade_consensus_run.jl`.
- Data: `data/raw/` (one copy of every raw file), `data/sets/` and
  `data/INDEX.tsv`, rebuilt from the tables in `benchmarks/` by
  `scripts/grading/report_data_sets.jl` (`data/README.md`).
- DFT: inputs in `dft/inputs/`, frozen molds in `dft/molds/`, accepted QE
  states in `dft/runs/` (untracked).
- Scripts are standalone command-line programs. Use their `--help` or header
  usage before assuming arguments.

## Setup and verification

The latest Julia release (1.13.2 on 2026-10-18), locally and on Viper, including
tests and subprocesses: regenerate the Manifest when Julia updates. The root
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
julia --project=. test/runtests.jl          # every test, end-to-end run included
julia --project=. test/runtests.jl --fast   # without the end-to-end run
```

Refactors of kept scripts must be verified byte-for-byte against stored
outputs. Changes to scientific behavior need a fresh raw run and an external grade.

Docs: `GKSwstype=100 julia --project=. docs/make.jl --build-only`.

## Data and commands

All data live in the untracked `data/` folder (`data/README.md`). Pass a set or
session folder with `--data-dir` (or `STMFIT_DATA_DIR` for `batch_full.jl`).

```bash
julia --project=. scripts/tools/inspect_one_file.jl <file.sxm>                 # one scan
julia -t 4 --project=. scripts/pipeline/batch_full.jl 48 --config config/chitosan.toml --data-dir DIR
julia --project=. scripts/pipeline/run_molecule_consensus_chitosan.jl --help   # promoted method
julia --project=. scripts/grading/grade_consensus_run.jl --run RUN --outdir NEW  # verify + grade
julia --project=. scripts/calibration/measure_calibration.jl <clean_scan.sxm>  # audit, no config written
```

If high-N fits disappear, inspect the physical `max_overlap` constraint before
changing selection logic.

## HPC and generated state

- HPC is deliberately boring: sync committed code (Manifest included),
  instantiate on the login node, dry-run, submit one `sbatch`, poll `squeue`,
  fetch outputs locally. No agent orchestration or state on the cluster that
  outlives a job. Never compute on login nodes; compute nodes have no internet.
- Viper sets `SBATCH_EXPORT=NONE`: pass variables with
  `sbatch --export=NONE,VAR=...`; shell variables before `sbatch` are not seen.
- Run cohort fits, multi-file batches and exhaustive searches on Viper, not
  locally. Local runs are for single-file checks and tests.
- The Raven group account has an 8-CPU quota. `batch_full.jl` uses at most
  4 threads per task.
- Long 10–20mer counting may need 8 h or more. QE jobs stay sequential under the
  one-node QOS limit.
- Do not commit `results/`, `data/`, `dft/runs/`, `docs/build/` or generated
  outputs. `/ptmp` is temporary: fetch HPC results.

## Documentation contract

- Add a dated `docs/src/journal.md` entry for every experiment, scientific
  decision, failed approach, bug fix or parameter change, with the reason.
  Keep its Current state and Open questions current.
- If behavior changes benchmark results, rerun the grade and update every
  cited headline (`README.md`, `docs/src/index.md`, runbook, selection and
  unit-assignment pages).
- Parameter additions or renames also require `docs/src/config.md` and
  `docs/src/calibration.md` updates.
- Trust the root README, `docs/src/` and the scripts; the packages have no
  README of their own.

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
