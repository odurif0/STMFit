# Runbook

How to run the promoted method, verify a run, grade it on the 6-mer benchmark,
and apply it to unknown 10–20mers. Cohort runs belong on the MPCDF cluster
(Viper). Local execution is for single-file checks and tests.

## Scientific constraints

- No expected N, target N, sequence or class count inside fitting, selection,
  assignment, thresholds, abstention or calibration. Benchmark manifests and
  `NKNNKN` are read only by grading scripts.
- The 1D fit is diagnostic only (`--skip-1d` is the default).
- Configs are frozen for a promoted run. A change is a new, declared experiment
  with a journal entry. It is graded once, without retuning after the grade.
- Keep suspicious 10–20mer images in the cohort and flag them; do not exclude
  them silently.

## Data

Raw `.sxm` files are not tracked ([data/README](https://github.com/odurif0/STMFit/blob/main/data/README.md)).
The 6-mer benchmark cohort is 146 raw scans; 145 of them have a confirmed
count of 6 (`benchmarks/chitosan_6mer_counting_confirmed.toml`), and the
control sequence `NKNNKN` is encoded `010010`
(`benchmarks/chitosan_240817_unit_sequences.tsv`). On Viper, the staged cohort
is `/ptmp/oldu/stmfit/full146_v1_inputs/full146_raw` (symlinks, unique basenames).

## Run the promoted method on Viper

The Slurm script `hpc/run_molecule_consensus.sbatch` defaults to the promoted
configs and the tracked molds. Viper sets `SBATCH_EXPORT=NONE`, so pass every
variable on the `sbatch` command line.

```bash
# local: sync the committed tree (+ the ignored Manifest.toml) to a new directory
C=$(git rev-parse --short HEAD); D=/u/oldu/code/STMFit_$C
git ls-files > /tmp/files.txt && echo Manifest.toml >> /tmp/files.txt
rsync -a --files-from=/tmp/files.txt ./ viper:$D/

# on the login node: instantiate only (compute nodes have no internet)
ssh viper "cd $D && /viper/u1/oldu/software/julia-1.13.0/bin/julia --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'"

# dry-run, then one submission
ssh viper "cd $D && J=/viper/u1/oldu/software/julia-1.13.0/bin/julia \
  IN=/ptmp/oldu/stmfit/full146_v1_inputs OUT=/ptmp/oldu/stmfit/run_$C && \
  STMFIT_PROJECT_DIR=\$PWD STMFIT_INPUT_DIR=\$IN STMFIT_OUTDIR=\$OUT JULIA_BIN=\$J \
    bash hpc/run_molecule_consensus.sbatch --dry-run && \
  sbatch --no-requeue --time=04:00:00 --output=\$OUT-%j.log \
    --export=NONE,STMFIT_PROJECT_DIR=\$PWD,STMFIT_INPUT_DIR=\$IN,STMFIT_OUTDIR=\$OUT,JULIA_BIN=\$J \
    hpc/run_molecule_consensus.sbatch"

# poll, then fetch everything back (/ptmp is temporary)
ssh viper "squeue -u oldu"
rsync -a viper:/ptmp/oldu/stmfit/run_$C results/ && rsync -a "viper:/ptmp/oldu/stmfit/run_$C-*.log" results/
```

Optional variables: `STMFIT_RAW_SUBDIR` (default `full146_raw`),
`STMFIT_COUNT_CONFIG`, `STMFIT_ASSIGNMENT_CONFIG`, `STMFIT_TEMPLATES`, and
`STMFIT_SELECTED_SUMMARY` (reuse per-scan counts; development only). The
full146 run takes about 1.5 h on 4 CPUs (MaxRSS about 5 GB).

## Verify and grade

After fetching the complete run, check integrity and grade outside the method:

```bash
julia --project=. test/grade_consensus_run.jl --run results/run_<C> \
    --outdir results/grade_<C> --data-dir data/chitosan_6mer \
    [--profile name=other_predictions.tsv ...]
```

The script first runs label-free checks (cohort, counts per stage, consensus
rules, fusion integrity, raw hashes). Any failure stops it before grading.
It then filters predictions to the 145 benchmark scans and writes:

- `units/report.md`, `units/summary.tsv`: full145 own-N unit grade for the final
  (`final`) and per-scan (`per_scan`) calls, plus any `--profile`;
- `counts_final.tsv`, `counts_per_scan.tsv`: count grades (exact and ±1);
- with `--data-dir`, `plots/`: one image per molecule sorted by benchmark
  outcome:
  - `1_all_exact/`;
  - `2_N_exact_calls_wrong/` (wrong lobes circled);
  - `3_N_wrong/`;
  - `4_not_in_benchmark/`.

  Each plot carries the truth in its title; `plot_annotations.tsv` holds the
  categories. This sorting reads labels and exists only in the grading step;
  a run's own `plots/` is label-free.

Expected for the promoted method (deterministic fits): final 771/870 correct,
88/145 exact chains, 865/870 classified, 94 errors; counts 137/145 exact,
145/145 within one; per-scan calls 702/33/863/161. Identical runs must be
identical byte for byte, except paths and logs.

A promotion needs a fresh raw run from a committed source, frozen configs,
passing verification, and a complete comparison (counts, correct, exact
chains, errors, abstentions) with no measure regressing.

## Counting only

```bash
julia -t 4 --project=. test/batch_full.jl 146 --config config/chitosan.toml \
    --data-dir RAW_DIR --outdir results/counting_<tag> --skip-1d
julia --project=. test/grade_chitosan_benchmark.jl \
    --manifest benchmarks/chitosan_6mer_counting_confirmed.toml \
    --results results/counting_<tag>/summary_overlap060_hard.tsv \
    --column N_selected --out results/counting_<tag>/grade.tsv
```

Fresh per-scan counting gives 123/145 exact and 143/145 within one. The
historical 129/145 comes from a saved summary that fresh runs do not
reproduce ([journal](journal.md)). For cluster arrays of `batch_full.jl`, see
`hpc/README.md` (`launch_remote.sh`).

## Unknown 10–20mer application

1. Collect repeat-scan candidates of the target molecules (headers only):
   ```bash
   julia --project=. test/collect_repeat_scans.jl --targets targets.txt \
       --raw-root RAW_TREE --config config/molecule_consensus.toml --outdir NEW_DIR
   ```
2. Run the method with the long-chain counting config:
   `STMFIT_COUNT_CONFIG=config/chitosan_10_20mer_adaptive_support_rescue.toml`
   (allow ≥ 8 h). Only registration decides which candidates are the same molecule.
3. Check integrity with `test/grade_consensus_run.jl --run DIR --verify-only`.
   There is no grade, because there are no labels.
4. Read `chain_report.tsv`, `assignment/review_queue.tsv`, the maps in
   `assignment/plots/` and the fusion parameters.

**Status (job 12023144, 25 targets + 109 repeat candidates).** 99 consecutive
pairs are linked and 17 counts change. For the 12 manually read targets, the
mean |N − N_manual| is 1.33. Two long chains fail badly (`251206_013`:
4 lobes, manual 12; `260220_083`: 2, earlier production 9). Their axial
support is truncated because the support baseline (10% quantile of the axial
profile) lands on the molecule when the chain fills most of the profile; a
curved centreline does not change this (prototype, 2026-09-30). The fusion fit
gives θ0 = 0.45, so per-scan calls are inconsistent across repeat scans.
**These outputs are not chemically validated.** A length-invariant support
baseline and consistent per-scan assignment on long chains are prerequisites
([journal](journal.md)).

## Tests

```bash
julia --project=. packages/STMFitCore.jl/test/runtests.jl
julia --project=. packages/STMSXMIO.jl/test/runtests.jl
julia --project=. packages/GaussianFit2D.jl/test/runtests.jl
for t in test/test_*.jl test/molecule_consensus/*.jl; do julia --project=. "$t"; done
julia --project=. test/test_reconstructed_pipeline.jl --e2e
```

`test_cc_mold_native.jl --real-cubes CUBE0 CUBE1 FRAME0 FRAME1 NEW_DIR` also
compares the native mold builder with the Python reference
(`test/lib/cc_mold_builder.py`) on the real DFT cubes.

## Other tools

| Script | Purpose |
|---|---|
| `test/inspect_one_file.jl FILE.sxm` | Circular vs elliptical 2D fits and plots for one scan |
| `test/summarize.jl SUMMARY.tsv` | Count statistics of a batch summary |
| `test/measure_calibration.jl SCAN.sxm` | Apparent width/spacing audit (no config written) |
| `test/validate_unit_predictions.jl`, `test/summarize_unknown_unit_qc.jl` | Prediction integrity and review queue |
| `test/compare_reconstructed_champion.jl` | Keyed comparison of two prediction tables |
| `test/catalog_sxm_candidates.jl --roots DIR[,DIR] ... --outdir NEW [--exclude TSV] [--copy-to CENTRAL]` | Technical triage of raw scans (no fit, no label), contact sheets, deduplicated central folder |
| `test/report_benchmark_candidates.jl --catalog TSV --count-config TOML --outdir NEW` | Grading side: usable scans not yet in the benchmark, with review sheets |
| `test/report_position_repeatability.jl --run RUN --data-dir RAW --outdir NEW` | Label-free lobe-position repeatability across repeated scans ([Pipeline](pipeline.md)) |
| `test/assignment_view_variants.jl RUN VARIANT NEW RAW` | Research harness: rerun the classifier stages of a run with other feature views |
