# Data directory

Raw Nanonis `.sxm` scans are **not tracked** (size). Place them here, one
subdirectory per dataset, or point scripts elsewhere with `--data-dir` (or
`STMFIT_DATA_DIR` for `batch_full.jl`). Basenames must be unique within a run.

| Dataset | Local | Viper | Content |
|---|---|---|---|
| 6-mer lab tree | `data/chitosan_6mer/` (14 session folders, 935 scans, plus 146 top-level copies of benchmark scans; source `/home/durif/Rebecca/data/data/`) | — | Inventory in `benchmarks/chitosan_6mer_data_inventory.tsv` |
| 6-mer Cu(100) triage | `/home/durif/Rebecca/data/chitosan_6mer_Cu100/` — `gardes/` (526), `ecartes/` (623), `tri.tsv`, `planches/` | — | Every (NKN)2 Cu(100) scan of the local tree and the institute share `hive.fkf.mpg.de/esi/DATA/STM/STM_DATA/(NKN)2`, deduplicated; built by `test/report_benchmark_triage.jl` from `benchmarks/chitosan_6mer_cu100_triage.tsv` (journal 2026-10-05) |
| 6-mer Ag(111) | `/home/durif/Rebecca/data/chitosan_6mer_Ag111/` (160 usable scans, `MANIFEST.tsv`) | — | Share sessions of July 2025; not triaged, not benchmark material |
| 6-mer benchmark cohort | symlinks in `results/reconstructed_cc_soft_v1/full146_raw/` | `/ptmp/oldu/stmfit/full146_v1_inputs/full146_raw/` | 146 scans; 145 graded (`benchmarks/chitosan_6mer_counting_confirmed.toml`) |
| Unknown 10–20mers | `/home/durif/Rebecca/data/10_20mer_analysis/` (25 targets) | `/ptmp/oldu/stmfit/data_10_20mer` | Repeat candidates are collected with `test/collect_repeat_scans.jl` |

The manifests in `benchmarks/` hold external grading metadata only (quality
classes, expected counts, control sequences). No fitting, selection or
assignment code reads them.

A new molecule:

```bash
mkdir data/<molecule>
cp /path/to/scans/*.sxm data/<molecule>/
julia --project=. test/measure_calibration.jl data/<molecule>/<clean_scan>.sxm
```

Then follow `docs/src/calibration.md` (Calibrating a new molecule).
