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
| Production/default short-chain chitosan | `config/chitosan.toml` | Default remains `gcv_with_robust_aicc_guard`. |
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

Use the benchmark data folder and the adaptive config:

```bash
JULIA_NUM_THREADS=4 julia --project=. test/batch_full.jl 39 \
  --data-dir /home/durif/Rebecca/data/data/20240817_LHe_Cu100 \
  --outdir results/best_plots_240817_adaptive_support_rescue \
  --tsv /tmp/opencode/chitosan_240817_primary_files.tsv \
  --config config/chitosan_adaptive_support_rescue.toml
```

Then grade externally against the benchmark manifest.  Labels stay outside the
fit/selection run:

```bash
julia --project=. test/grade_chitosan_benchmark.jl \
  benchmarks/chitosan_240817.toml \
  results/best_plots_240817_adaptive_support_rescue/summary_overlap060_hard.tsv \
  --out results/benchmark_grades/chitosan_240817_adaptive_support_rescue_N_selected.tsv \
  --column N_selected
```

Known validation from the current development pass: the default-config
(`config/chitosan.toml`) workflow reaches `N_selected = 39/39` exact on the
240817 primary benchmark, with all four `clean_target` files (`017`, `019`,
`043`, `058`) reporting `N_selected = 6`.  `043` is recovered by the
up-when-ambiguous guard branch (its `N_eff = 5`, but `robust_AICc_N = 6` on an
ambiguous file; see the Research Journal §2026-06-17).  Guard-sensitive files
`240817_058.sxm` and `240817_019.sxm` should report `N_selected = 6` even though
`N_eff = 7`.  Re-measure any time with the grade script below.

Reproducibility note: the batch is deterministic run-to-run on a given machine
(verified identical `N_selected` across 3 consecutive runs on 2026-06-17).
Divergences between a past recorded number and a fresh run indicate intervening
code changes, not run-to-run noise.

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
  `ell_robust_aicc`.
- `refined_policy`: audit trail for whether adaptive support was kept, accepted,
  rejected, or guarded.
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

- The adaptive configs are currently experimental and are not the default
  production config until explicitly promoted.
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
3. **Keep the `[selection]` defaults** (`gcv_ambiguity_rel_threshold = 0.05`,
   `robust_guard_nu = 8.0`) as a starting point. If the new molecule's lobe
   statistics differ markedly from chitosan, run
   `test/sensitivity_thresholds.jl` to check whether `N_selected` is robust to
   the threshold; re-calibrate only if it is sensitive.
4. **Exclude non-target files** via `--exclude-from results/<molecule>_exclude.txt`
   (one filename per line) rather than editing the batch code.

The selection path never uses an expected `N` or benchmark label, so the same
guard logic applies unchanged. What may need attention is *how often* the
up-when-ambiguous branch fires for a molecule whose GCV curve has a different
shape — hence the sensitivity check.
