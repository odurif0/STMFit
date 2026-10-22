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
(`benchmarks/chitosan_240817_unit_sequences.tsv`). All data live in the
untracked `data/` folder ([data/README](https://github.com/odurif0/STMFit/blob/main/data/README.md)):
one copy of every raw file in `data/raw/`, working sets as links in
`data/sets/` (the cohort is `data/sets/6mer_Cu100_benchmark/`), and
`data/INDEX.tsv` with every scan and its decision. The sets and the index are
rebuilt from the tables in `benchmarks/` by `scripts/grading/report_data_sets.jl`. On
Viper, the staged cohort is `/ptmp/<user>/stmfit/full146_v1_inputs/full146_raw`
(symlinks, unique basenames).

## Run the promoted method on Viper

The Slurm script `hpc/run_molecule_consensus.sbatch` defaults to the promoted
configs and the tracked molds. Viper sets `SBATCH_EXPORT=NONE`, so pass every
variable on the `sbatch` command line.

```bash
# local: sync the committed tree (+ the ignored Manifest.toml) to a new directory
C=$(git rev-parse --short HEAD); D=/u/<user>/code/STMFit_$C
git ls-files > /tmp/files.txt && echo Manifest.toml >> /tmp/files.txt
rsync -a --files-from=/tmp/files.txt ./ viper:$D/

# on the login node: instantiate only (compute nodes have no internet)
ssh viper "cd $D && /viper/u1/<user>/software/julia-1.13.2/bin/julia --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'"

# dry-run, then one submission
ssh viper "cd $D && J=/viper/u1/<user>/software/julia-1.13.2/bin/julia \
  IN=/ptmp/<user>/stmfit/full146_v1_inputs OUT=/ptmp/<user>/stmfit/run_$C && \
  STMFIT_PROJECT_DIR=\$PWD STMFIT_INPUT_DIR=\$IN STMFIT_OUTDIR=\$OUT JULIA_BIN=\$J \
    bash hpc/run_molecule_consensus.sbatch --dry-run && \
  sbatch --no-requeue --time=04:00:00 --output=\$OUT-%j.log \
    --export=NONE,STMFIT_PROJECT_DIR=\$PWD,STMFIT_INPUT_DIR=\$IN,STMFIT_OUTDIR=\$OUT,JULIA_BIN=\$J \
    hpc/run_molecule_consensus.sbatch"

# poll, then fetch everything back (/ptmp is temporary)
ssh viper "squeue -u <user>"
rsync -a viper:/ptmp/<user>/stmfit/run_$C results/ && rsync -a "viper:/ptmp/<user>/stmfit/run_$C-*.log" results/
```

Optional variables: `STMFIT_RAW_SUBDIR` (default `full146_raw`),
`STMFIT_COUNT_CONFIG`, `STMFIT_ASSIGNMENT_CONFIG`, `STMFIT_TEMPLATES`, and
`STMFIT_SELECTED_SUMMARY` (reuse per-scan counts; development only). The
full146 run takes about 1.5 h on 4 CPUs (MaxRSS about 5 GB).

## Verify and grade

After fetching the complete run, check integrity and grade outside the method:

```bash
julia --project=. scripts/grading/grade_consensus_run.jl --run results/run_<C> \
    --outdir results/grade_<C> --data-dir data/sets/6mer_Cu100_benchmark \
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
julia -t 4 --project=. scripts/pipeline/batch_full.jl 146 --config config/chitosan.toml \
    --data-dir RAW_DIR --outdir results/counting_<tag> --skip-1d
julia --project=. scripts/grading/grade_chitosan_benchmark.jl \
    --manifest benchmarks/chitosan_6mer_counting_confirmed.toml \
    --results results/counting_<tag>/summary_overlap060_hard.tsv \
    --column N_selected --out results/counting_<tag>/grade.tsv
```

Fresh per-scan counting gives 123/145 exact and 143/145 within one. The
historical 129/145 comes from a saved summary that fresh runs do not
reproduce ([journal](journal.md)).

## Unknown 10–20mer application

0. Triage the raw scans (done for the Cu(100) 10–20mer tree on 2026-10-21):
   `scripts/data/catalog_sxm_candidates.jl` for the technical status, then
   `scripts/grading/report_candidate_review.jl --profile 10_20mer` for the
   review page of the usable scans (decisions go to
   `benchmarks/chitosan_10_20mer_cu100_triage.tsv`). The working set is
   `data/sets/10_20mer_curated` (19 scans the experimenters selected,
   `benchmarks/chitosan_10_20mer_curated.tsv`, 2026-10-22); the bulk
   candidates are mostly unusable.
1. Collect repeat-scan candidates of the target molecules (headers only):
   ```bash
   julia --project=. scripts/data/collect_repeat_scans.jl --targets targets.txt \
       --raw-root RAW_TREE --config config/molecule_consensus.toml --outdir NEW_DIR
   ```
2. Run the method with the long-chain counting config:
   `STMFIT_COUNT_CONFIG=config/chitosan_10_20mer_adaptive_support_rescue.toml`
   (allow ≥ 8 h). Only registration decides which candidates are the same molecule.
3. Check integrity with `scripts/grading/grade_consensus_run.jl --run DIR --verify-only`.
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
julia --project=. test/runtests.jl          # every test, end-to-end run included
julia --project=. test/runtests.jl --fast   # without the end-to-end run
```

`test_cc_mold_native.jl --real-cubes CUBE0 CUBE1 FRAME0 FRAME1 NEW_DIR` also
compares the native mold builder with the Python reference
(`test/reference/cc_mold_builder.py`) on the real DFT cubes.

## Other tools

| Script | Purpose |
|---|---|
| `scripts/tools/inspect_one_file.jl FILE.sxm` | Circular vs elliptical 2D fits and plots for one scan |
| `scripts/calibration/measure_calibration.jl SCAN.sxm` | Apparent width/spacing audit (no config written) |
| `scripts/pipeline/validate_unit_predictions.jl`, `scripts/pipeline/summarize_unknown_unit_qc.jl` | Prediction integrity and review queue |
| `scripts/data/catalog_sxm_candidates.jl --roots DIR[,DIR] ... --outdir NEW [--exclude TSV] [--copy-to CENTRAL]` | Technical triage of raw scans (no fit, no label), contact sheets, deduplicated central folder |
| `scripts/grading/report_benchmark_candidates.jl --catalog TSV --count-config TOML --outdir NEW` | Grading side: usable scans not yet in the benchmark, with review sheets |
| `scripts/grading/report_benchmark_triage.jl --catalog TSV ... --count-config TOML --sheets NEW` | Grading side: fixed-scale review sheets for the visual triage |
| `scripts/grading/report_candidate_review.jl --data data --triage TSV --outdir NEW [--suggestions TSV] [--profile 6mer\|10_20mer]` | Grading side: local review page (cards grouped by molecule, decision export) for 6-mer benchmark candidates or 10–20mer application scans |
| `scripts/grading/report_data_sets.jl --data data --table TSV ...` | Grading side: rebuild `data/sets/` and `data/INDEX.tsv` from the decision tables |
| `scripts/grading/search_reference_cohort.jl --protocol TOML --mode nested\|final\|null\|report\|sealed --outdir NEW` | Grading side, campaign B: reference cohort of a frozen assignment chosen by its benchmark score (`benchmarks/reference_cohort_search_v1.toml`, multistart `_v2.toml`; Viper: `hpc/run_reference_cohort_search.sbatch`) |
| `scripts/pipeline/build_frozen_assignment.jl --run RUN --config TOML --out MODEL.tsv`, `scripts/pipeline/apply_frozen_assignment.jl --run RUN --model MODEL.tsv --config TOML --outdir NEW` | Campaign A: build a fully frozen, label-free assignment model on a reference run, and apply it unchanged to any run ([Pipeline](pipeline.md)) |
| `scripts/grading/grade_frozen_assignment.jl --run RUN --config TOML --outdir NEW [--other-run RUN2]` | Grading side, campaign A: held-out 2024 sessions, 2024 → 2025 transfer, invariance across runs |
| `scripts/pipeline/report_position_repeatability.jl --run RUN --data-dir RAW --outdir NEW` | Label-free lobe-position repeatability across repeated scans ([Pipeline](pipeline.md)) |
