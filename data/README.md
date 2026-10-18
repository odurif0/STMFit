# Data directory

Raw Nanonis `.sxm` scans are **not tracked** (size). Place them here, one
subdirectory per dataset, or point scripts elsewhere with `--data-dir` (or
`STMFIT_DATA_DIR` for `batch_full.jl`). Basenames must be unique within a run.

| Dataset | Local | Viper | Content |
|---|---|---|---|
| **6-mer (NKN)2 on Cu(100): reference** | `/home/durif/Rebecca/data/chitosan_6mer_Cu100/` (link `data/chitosan_6mer_Cu100`) | — | Every Cu(100) scan once (1,149): `gardes/benchmark/` (146, used by the benchmark runs), `gardes/candidats/` (380, to validate), `ecartes/<motif>/` (623). `LISTE.tsv` lists every scan, `LISEZMOI.md` explains the folders. Built by `test/report_benchmark_triage.jl` from `benchmarks/chitosan_6mer_cu100_triage.tsv` |
| 6-mer benchmark cohort (runs) | symlinks in `results/reconstructed_cc_soft_v1/full146_raw/` → `gardes/benchmark/` | `/ptmp/oldu/stmfit/full146_v1_inputs/full146_raw/` | 146 scans; 145 graded (`benchmarks/chitosan_6mer_counting_confirmed.toml`) |
| 6-mer raw source | institute share `hive.fkf.mpg.de/esi/DATA/STM/STM_DATA/(NKN)2`; local copy of the 2024 sessions in `/home/durif/Rebecca/data/data/` (with spectra) | — | Origin of every scan (`origine` column of `LISTE.tsv`); not used directly |
| 6-mer Ag(111) | `/home/durif/Rebecca/data/chitosan_6mer_Ag111/` (160 usable scans, `MANIFEST.tsv`) | — | Share sessions of July 2025; not triaged, not benchmark material |
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
