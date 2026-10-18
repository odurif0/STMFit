# Data directory

All project data live here, inside the repository folder, and are **not
tracked** by git (`.gitignore`: `data/*` except this README and `.gitkeep`).
The folder is self-contained: no data path points outside the repository
(run inputs under `results/` are relative links into `data/`). Pass a folder
with `--data-dir` (or `STMFIT_DATA_DIR` for `batch_full.jl`); basenames must be
unique within a run.

| Folder | Content |
|---|---|
| `chitosan_6mer_Cu100/` | **Reference for the 6-mer (NKN)2 on Cu(100)**: every scan once (1,149, local tree and institute share, deduplicated). `gardes/benchmark/` = the 146 scans the benchmark runs use (145 graded in `benchmarks/chitosan_6mer_counting_confirmed.toml`, plus `240310_Cu100009`); `gardes/candidats/` = 380 kept by the first visual pass, to validate; `ecartes/<motif>/` = 623 discarded. `LISTE.tsv` lists every scan, `LISEZMOI.md` explains the folders. Rebuilt by `test/report_benchmark_triage.jl` from `benchmarks/chitosan_6mer_cu100_triage.tsv` |
| `chitosan_6mer_Ag111/` | 160 usable 6-mer scans on Ag(111) (July 2025 share sessions), `MANIFEST.tsv`; not triaged, not benchmark material |
| `raw/chitosan_6mer_2024_lab_tree/` | Raw 2024 (NKN)2 lab tree: 14 session folders (935 scans, spectra `.dat`/`.3ds`, notes) and `chitosan_manual_annotations.md` |
| `raw/chitosan_10_20mer/` | Raw 10–20mer scans (about 1,670 `.sxm` with their files) |
| `raw/chitosan_10_20mer_analysis/` | The 25 unknown 10–20mer target scans |
| `raw/misc/` | Loose files: `SequenceAnalysis.xlsx` (manual sequence reads), `_KA_SOLVE.pptx`, `251206_034` scan and notes, profiles |

On Viper the staged copies are temporary: `/ptmp/oldu/stmfit/full146_v1_inputs/full146_raw/`
(6-mer cohort) and `/ptmp/oldu/stmfit/data_10_20mer` (10–20mers). The master
copy of the raw (NKN)2 scans is the institute share
`hive.fkf.mpg.de/esi/DATA/STM/STM_DATA/(NKN)2`; the `origine` column of
`LISTE.tsv` gives each scan's `<session>/<file>` there. Files under `raw/` that
also sit in `/home/durif/Rebecca/data` are hard links (one copy on disk).

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
