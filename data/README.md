# Data directory

Raw Nanonis `.sxm` scans are **not tracked** (size). Place them here, one
subdirectory per dataset, or point scripts elsewhere with `--data-dir` (or
`STMFIT_DATA_DIR` for `batch_full.jl`). Basenames must be unique within a run.

| Dataset | Local | Viper | Content |
|---|---|---|---|
| 6-mer lab tree | `data/chitosan_6mer/` (session folders, ~1,080 scans; source `/home/durif/Rebecca/data/data/`) | — | Inventory in `benchmarks/chitosan_6mer_data_inventory.tsv` |
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
