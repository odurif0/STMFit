# Chitosan Benchmark and 10–20mer Runbook

This page is the hand-off document for reproducing the current chitosan workflow
without relying on prior conversation context.  The goal is a label-free fitting
and model-selection pipeline that is externally graded on the benchmark and then
reused unchanged for curated 10–20mer images, except for the allowed `N` range.

## Scientific constraints

- Do not use an expected `N`, target `N`, or preferred `N` inside fitting or
  model selection.
- Benchmark labels are evaluation-only.  They may be used by grading scripts, but
  not by `test/batch_full.jl`, TOML calibration, or model-selection logic.
- Keep the circular-to-elliptical 2D pipeline: circular sweep first, then local
  elliptical refinement from circular solutions.
- Do not globally widen support for all files as a shortcut; global support
  changes were found to regress clean benchmark files.
- Ambiguous/suspicious 10–20mer images should be retained and annotated by QC or
  confidence fields, not silently excluded.

## Current configs

| Purpose | Config | Notes |
|---|---|---|
| Production/default short-chain chitosan | `config/chitosan.toml` | Current default is `support_midpoint_hybrid`, promoted as the best label-free full145 counting rule so far. |
| Benchmark-aligned adaptive workflow | `config/chitosan_adaptive_support_rescue.toml` | Experimental generic support-rescue policy. |
| 10–20mer adaptive workflow | `config/chitosan_10_20mer_adaptive_support_rescue.toml` | Same policy as benchmark adaptive workflow, with `n_max = 24`. |
| Legacy 10–20mer baseline | `config/chitosan_10_20mer.toml` | Raw GCV baseline/comparison. |
| Legacy 10–20mer support-rescue diagnostics | `config/chitosan_10_20mer_rescue*.toml` | Comparison/audit only, not ground truth. |

## Generic adaptive-support policy

`selection_policy = "adaptive_support_rescue"` performs the following steps for
each file:

1. Fit the standard support using the configured support detector.
2. Select `N_eff` from the circular/elliptical 2D sweeps using the configured
   criterion, currently GCV.
3. Trigger a permissive support-rescue pass only if `N_eff` is at the objective
   support-feasibility ceiling.  The ceiling uses the detected support length and
   physical spacing/sigma bounds; it does not use labels.
4. Accept the rescue only if support length increases, selected `N` increases,
   circular/elliptical counts remain coherent, and the selected count is feasible
   on the rescued support.
5. Apply the robust-AICc guard down-only on the active support.  If the guard
   fails, keep the active `N_eff` rather than failing the file.

The optional `adaptive_robust_guard_max_drop` parameter exists only for diagnostic
experiments.  It is not used by the benchmark-aligned workflow.

## Reproduce the short-chain benchmark workflow

Use the benchmark data folder and the robust-AICc guard override when the goal is
to reproduce the original 240817 primary-benchmark validation (`39/39`). The
plain `config/chitosan.toml` default is now optimized for the expanded full145
external counting grade instead.

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 39 \
  --data-dir /home/durif/Rebecca/data/data/20240817_LHe_Cu100 \
  --outdir results/best_plots_240817_robust_guard \
  --tsv /tmp/opencode/chitosan_240817_primary_files.tsv \
  --config config/chitosan.toml \
  --selection-policy gcv_with_robust_aicc_guard
```

Then grade externally against the benchmark manifest.  Labels stay outside the
fit/selection run:

```bash
julia --project=. test/grade_chitosan_benchmark.jl \
  --manifest benchmarks/chitosan_240817.toml \
  --results results/best_plots_240817_robust_guard/summary_overlap060_hard.tsv \
  --out results/benchmark_grades/chitosan_240817_robust_guard_N_selected.tsv \
  --column N_selected
```

Known validation from the current development pass: the robust-AICc guard alone
reaches `N_selected = 39/39` exact on the 240817 primary benchmark, with all four
`clean_target` files (`017`, `019`, `043`, `058`) reporting `N_selected = 6`.
`043` is recovered by the up-when-ambiguous guard branch (its `N_eff = 5`, but
`robust_AICc_N = 6` on an ambiguous file; see the Research Journal §2026-06-17).
The promoted `config/chitosan.toml` default is now `support_midpoint_hybrid`
because it improves the expanded 145-file external counting grade from `106/145`
exact to `129/145` exact (`143/145` within one lobe), confirmed by a full Viper
batch run with 146 `ok` rows and grading against the 145-file manifest.
Re-measure any time with the grade script below.

Reproducibility note: the batch is deterministic run-to-run on a given machine
(verified identical `N_selected` across 3 consecutive runs on 2026-06-17).
Divergences between a past recorded number and a fresh run indicate intervening
code changes, not run-to-run noise.

### Expanded counting benchmark

The canonical 240817 benchmark above remains the reproducible validation set.
For broader external counting grades, use:

```text
benchmarks/chitosan_6mer_counting_confirmed.toml
```

This manifest contains 145 visually/manually confirmed `expected_N=6` files
derived from `benchmarks/chitosan_6mer_preassignment_review.tsv`. It is still an
external grading manifest only and must not be read by fitting or selection code.
For unit assignment, the same 145 files are the benchmark scope and the external
control sequence is `NKNNKN` (`010010` or `101101`, depending on the 0/1 identity
convention). That sequence is grading/control information only: it must not enter
the assignment method, threshold choice, abstention rule, composition prior, or
any calibration intended to extrapolate to unknown systems. Older pending lists
track provenance/curation state and should not be read as a license to shrink the
0/1/? benchmark to the 240817 subset.

## Run the 10–20mer adaptive workflow

Use the same adaptive policy with the 10–20mer config.  `--skip-1d` is
recommended because the 1D panels are not required for this long-chain workflow
and can be expensive.

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 25 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir results/10_20mer_analysis_adaptive_support_rescue \
  --tsv results/10_20mer_analysis_adaptive_support_rescue/triage_unused.tsv \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --skip-1d
```

For a targeted smoke test of the support-rescue behavior:

```bash
julia --project=. test/batch_full.jl 1 \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir /tmp/opencode/stmfit_10_20_adaptive_target \
  --tsv /tmp/opencode/one_260220_083.tsv \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --skip-1d
```

Expected targeted result from the current pass: `260220_083.sxm` accepts support
rescue and reports `N_selected = 9`.

## Restarting and HPC usage

`test/batch_full.jl` appends to `summary_overlap060_hard.tsv` and skips files
already present in the selected `--outdir`.  If a run times out, rerun the same
command with the same output directory to continue from remaining files.

For HPC jobs, keep the same command-line arguments and only change scheduler
details, thread count, and output paths.  Prefer writing to a fresh output
directory per experiment.  Keep the final `summary_overlap060_hard.tsv`, the
per-file folders, and the exact config used.

## Outputs to inspect

The main table is:

```text
<outdir>/summary_overlap060_hard.tsv
```

Important columns:

- `N_eff`: raw effective 2D selection before refined policy/guard.
- `N_selected`: primary reported count after the configured batch policy.
- `selection_policy`: policy requested by config/CLI.
- `selection_source`: source of the primary selection, e.g. `ell` or
  `ell_robust_aicc`, plus `support_midpoint_down` or `support_midpoint_up` when
  the default support-midpoint layer makes the final one-step adjustment.
- `refined_policy`: audit trail for whether the robust-AICc/adaptive support
  stage was kept, accepted, rejected, or guarded before the final selection.
- `robust_aicc_N`: robust-AICc diagnostic/guard count.
- `support_2D_ell_nm`, `support_2D_circ_nm`: active support length after any
  accepted rescue.
- `best_plot`, `file_dir`: locations of visual outputs and per-file artifacts.

The 1D slide fit is **off by default** (it never enters `N_selected`; it is a
diagnostic only). `N_1D` and 1D comparison columns are therefore expected to be
`NA`, and 1D panels are not drawn. To re-enable the 1D diagnostic (e.g. to cross-
check a suspected 2D under-detection), pass `--no-skip-1d`.

## Minimal validation checklist after code/config changes

Run the default config smoke checks:

```bash
julia --project=. test/batch_full.jl 0 --config config/chitosan.toml
julia --project=. test/batch_full.jl 0 --config config/chitosan_10_20mer.toml --skip-1d
git diff --check
```

Run adaptive targeted checks when touching support rescue or selection logic:

```bash
julia --project=. test/batch_full.jl 1 \
  --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --data-dir /home/durif/Rebecca/data/10_20mer_analysis \
  --outdir /tmp/opencode/stmfit_10_20_adaptive_target \
  --tsv /tmp/opencode/one_260220_083.tsv \
  --skip-1d

julia --project=. test/batch_full.jl 2 \
  --config config/chitosan_adaptive_support_rescue.toml \
  --data-dir /home/durif/Rebecca/data/data/20240817_LHe_Cu100 \
  --outdir /tmp/opencode/stmfit_240817_adaptive_guard_targets \
  --tsv /tmp/opencode/two_240817_guard_regressions.tsv
```

Expected targeted results:

- `260220_083.sxm`: `N_selected = 9`.
- `240817_058.sxm`: `N_selected = 6`.
- `240817_019.sxm`: `N_selected = 6`.

## Known caveats

- The support-midpoint hybrid default is the best current label-free full145
  counting rule, not a universal solution; re-validate before using it as a new
  molecule's default.
- The adaptive configs are still experimental and are not the default production
  config until explicitly promoted.
- The legacy 10–20mer finalizer output is a comparison artifact, not ground
  truth.
- Full 10–20mer runs can be slow locally because the robust guard refits a broad
  range.  Use restartable output directories or HPC for complete reruns.
- If a file looks chemically or visually suspect, do not silently exclude it from
  10–20mer analysis; record the QC concern and keep the result available.

## Calibrating a new molecule

The pipeline is molecule-agnostic in its core (chain-of-Gaussians model,
label-free selection); only the calibration constants differ. To analyse a new
chain-like molecule under similar STM conditions:

1. **Copy the template**: `cp config/template.toml config/<molecule>.toml`.
2. **Re-derive the `[model]` values** from a few representative scans — the
   template comments explain each (FWHM → sigma, observed pitch → spacing,
   support length). These are the load-bearing changes.
3. **Treat the `[selection]` defaults** (`gcv_ambiguity_rel_threshold = 0.05`,
   `robust_guard_nu = 8.0`, `support_midpoint_up_gcv_rel_threshold = 0.30`) as
   chitosan-calibrated starting points. If the new molecule's lobe statistics
   differ markedly from chitosan, run
   `test/sensitivity_thresholds.jl` to check whether `N_selected` is robust to
   the threshold; re-calibrate only if it is sensitive.
4. **Exclude non-target files** via `--exclude-from results/<molecule>_exclude.txt`
   (one filename per line) rather than editing the batch code.

The selection path never uses an expected `N` or benchmark label, so the same
logic is label-free. What may need attention is *how often* the
up-when-ambiguous and support-midpoint branches fire for a molecule whose support
geometry or GCV curve has a different shape — hence the sensitivity check.

## Unit assignment workflow (GlcNAc/GlcN)

See [Unit Assignment](unit_assignment.md) for the full pipeline description.

### Unknown chitosan sequence production

<!-- UNKNOWN-CHITOSAN-WORKFLOW:START -->

For a chain with no unit-identity labels, start from the selected-N feature TSVs
produced by the label-free counting pipeline and write production artifacts only:

```bash
julia --project=. test/run_unknown_unit_assignment.jl \
    --features results/unit_assignment/<sample>_features_local.tsv \
    --split-features results/unit_assignment/<sample>_features_split.tsv \
    --patches results/unit_assignment/<sample>_patches_bwd.tsv \
    --profile default \
    --outdir results/unit_assignment/<sample>_unknown

julia --project=. test/validate_unit_predictions.jl \
    --predictions results/unit_assignment/<sample>_unknown/predictions_base_split_log_skew.tsv \
    --features results/unit_assignment/<sample>_features_local.tsv

python3 test/plot_unit_assignment.py \
    --features results/unit_assignment/<sample>_features_local.tsv \
    --predictions results/unit_assignment/<sample>_unknown/predictions_base_split_log_skew.tsv \
    --out-dir results/unit_assignment/<sample>_unknown/plots \
    --mode all

julia --project=. test/summarize_unknown_unit_qc.jl \
    --predictions results/unit_assignment/<sample>_unknown/predictions_base_split_log_skew.tsv \
    --plots-dir results/unit_assignment/<sample>_unknown/plots \
    --out results/unit_assignment/<sample>_unknown/review_queue.tsv
```

The production directory contains fixed-profile prediction TSVs, `summary.tsv`,
`manifest.tsv`, validation logs, plots, and `review_queue.tsv`. Treat `?` as an
explicit abstention requiring review, not as a dropped lobe.

<!-- UNKNOWN-CHITOSAN-WORKFLOW:END -->

### Challenger terminal status (T9 closure)

<!-- T9-TERMINAL-STATUS:START -->

The constant-current T3 lane is terminal `BLOCKED`: accepted GlcNAc cubes have
multiple nonunique isovalue branches, and no branch-selection policy was
predeclared. The hierarchical T5 identifiability gate passed all 13 held-out
dates and the 500-seed scan bootstrap, with a 95% lower bound of
`0.55621530226471905`. T6 integration and leakage checks were confirmed after
the URI-path and partial-view QC corrections.

T7 therefore ended as `NO_ELIGIBLE_CHALLENGER`. The hierarchical lane lacks
durable forward/backward evidence and the required common per-lobe audit table,
so no candidate was frozen. T8 is `SKIPPED_NO_ELIGIBLE_CHALLENGER`; the grader
invocation count was zero, and the one-shot grade budget remains unused.
`config/unit_assignment_candidate.toml` remains
`grade_status = "locked"` with `provenance.status = "pending"` and no frozen
hash. Its unchanged file SHA-256 is
`9dac84437ef0a9c2118b77d4e371efd71a5365595794167a112a338ca3e4a1aa`.

Because no new grade ran, the existing benchmark headline remains historical
and current: `78.9%` classified physical accuracy and `17/145` exact chains.
The classified percentage uses only classified positions; it is not the
fixed-denominator honest view, which keeps missing/abstained benchmark positions
in the denominator. This closure does not establish promotion or new benchmark
validation.

Final review hardened, but did not promote, both diagnostic lanes. Hierarchical
unstable/nonmonotone fits abstain, held-out evaluation rejects them before
scoring, and optional split/backward descriptors now reach fitting through a
records-based pipeline. Prediction provenance binds every consumed primary,
split, and backward artifact by stable role, normalized path, and byte SHA-256,
together with resolved views and model options. Constant-current calibration
scans the full declared interval under the validated
`--isovalue-scan-intervals` policy (default 1024), records
`isovalue_scan_intervals` in provenance, and accepts only one continuous
fixed-support root branch at that declared resolution. Multiple roots or support
discontinuities are rejected rather than resolved by a branch preference. The
accepted GlcNAc cube still violates that uniqueness contract, so this API
behavior does not unblock T3.

The manifest checker rejects multiline TOML strings, requires exactly one real
top-level `[grader_only]` table, and scans forbidden non-grader keys, values, and
paths case-insensitively. Its digest binds every other non-grader byte, including
comments, formatting, and line endings, excluding only real
`candidate.frozen_hash` and lifecycle-only `candidate.grade_status` assignments.
Matching provenance and digest are mandatory in `frozen_once` and `graded`, and
the lifecycle transition does not rehash. This is not stateless historical
proof. The current candidate remains locked and byte-unchanged; these stronger
future-state checks did not freeze it or consume the grade budget.

`hierarchical_equalprior` remains executable for unknown-production diagnostics
but is diagnostic, not frozen or promoted. A reproducible label-free diagnostic
run is:

```bash
julia --project=. test/run_unknown_unit_assignment.jl \
    --features results/unit_assignment/<sample>_features_local.tsv \
    --profile hierarchical_equalprior \
    --outdir results/unit_assignment/<sample>_hierarchical_diagnostic

julia --project=. test/validate_unit_predictions.jl \
    --predictions results/unit_assignment/<sample>_hierarchical_diagnostic/predictions_hierarchical_equalprior.tsv \
    --features results/unit_assignment/<sample>_features_local.tsv
```

<!-- T9-TERMINAL-STATUS:END -->

### Benchmark and post-hoc grading

The workflow for benchmark diagnostics is:

1. **Fill the ground truth** in
   `benchmarks/chitosan_240817_unit_sequences.tsv` (column `sequence`:
   ordered 0/1 along t_nm increasing; 0 = GlcN, 1 = GlcNAc).

2. **Extract per-lobe features** (re-runs the fit, ~10–15 min/file):
   ```bash
   STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. \
       test/extract_lobe_features.jl \
       --config config/chitosan.toml \
       --out results/unit_separability/lobe_features.tsv
   ```

3. **Run separability analysis** (label-free + with-truth):
   ```bash
   julia --project=. test/analyze_unit_separability.jl \
       --features results/unit_separability/lobe_features.tsv \
       --truth benchmarks/chitosan_240817_unit_sequences.tsv \
       --out results/unit_separability
   ```
   Check `results/unit_separability/separability_report.txt`:
   - ΔBIC (k=1−k=2) > 10 → bimodal (two distinct populations)
   - Physical accuracy > 75% → good separability
   - Physical accuracy < 60% → poor; use patch/mold diagnostics or DFT-STM maps

**Decision point:** if the separability report shows physical accuracy > 60%,
continue with label-free assignment diagnostics. If < 60%, prefer the connected
mold / DFT-STM path or a different STM bias/condition.

The ground truth is **never read by the fitter** — only by the separability
analysis (`--with-truth`) and the grading script.
