# Data directory

All project data live here, inside the repository folder, and are **not
tracked** by git (`.gitignore`: `data/*` except this README and `.gitkeep`).
The folder is self-contained: no data path points outside the repository.

```
data/
  INDEX.tsv      one row per raw scan: where, what, what was decided (generated)
  raw/           one copy of every raw file, by sample then session; never edited
    NKN2/<session>/      6-mer (NKN)2: Cu(100) 2024 and 2025, Ag(111) 2025
    10_20mer/<session>/  10-20mer sessions, Cu(100) and Ag(111)
    notes/               SequenceAnalysis.xlsx, manual annotations, slides, profiles
  sets/          working sets: relative symlinks into raw/ (generated)
    6mer_Cu100_benchmark/    benchmark scans (493: the 145-scan graded cohort,
                             240310_Cu100009, 347 added on 2026-10-19)
    10_20mer_targets/        the 25 unknown 10-20mer targets
    10_20mer_candidates/     usable Cu(100) 10-20mer scans awaiting the visual
                             review (2026-10-21)
```

**Rules.**

- `raw/` is never edited.
- Decisions live in tracked tables:
  - `benchmarks/chitosan_6mer_cu100_triage.tsv` holds one row per 6-mer Cu(100)
    scan. A non-empty `set` places the scan in that set; an empty `set` means
    discarded, and `motif` says why. Discarded scans get no folder.
  - `benchmarks/chitosan_10_20mer_targets.tsv` lists the 10-20mer targets.
  - `benchmarks/chitosan_10_20mer_cu100_triage.tsv` holds one row per other
    Cu(100) 10-20mer scan: `10_20mer_candidates` (usable, to review) or
    discarded with a motif (overview, cut by the frame, unreadable).
- To change a decision, edit the table, then rebuild:

  ```bash
  julia -t 4 --project=. scripts/grading/report_data_sets.jl --data data \
      --table benchmarks/chitosan_6mer_cu100_triage.tsv --table benchmarks/chitosan_10_20mer_targets.tsv \
      --table benchmarks/chitosan_10_20mer_cu100_triage.tsv
  ```

  This rewrites `sets/` and `INDEX.tsv` (about 40 s). It refuses to run if a
  table does not match `raw/`.
- `INDEX.tsv` columns:
  - where: `raw_path`;
  - what: `sample`, `substrate`, `session`, `file`, `date`, `time`, `frame_nm`,
    `pixels`, `bias_v`, `sha256`;
  - decision: `status` (set name, `discarded`, `not_triaged` or `duplicate`),
    `motif`, `decided_by`, `note`.
- Run inputs under `results/` are relative links into `raw/`.
- Pass a set or a session folder with `--data-dir` (or `STMFIT_DATA_DIR` for
  `batch_full.jl`). Basenames must be unique within a run.

**Provenance.**

- The master copy of the (NKN)2 scans is the institute share
  `hive.fkf.mpg.de/esi/DATA/STM/STM_DATA/(NKN)2`, and `raw/NKN2` is an exact
  mirror of it (3,822 files, checked 2026-10-18).
- The 2024 sessions and `raw/10_20mer`, `raw/notes` are hard links to the
  originals in `/home/durif/Rebecca/data`, so they take one copy on disk.
- Substrate in `INDEX.tsv` comes from the session name (`Cu100`, `Ag111`); it
  stays empty for `10_20mer/CuLHe03262026_depo2`, whose name gives no facet.
- `raw/10_20mer` keeps the files of the original tree, including derived
  `.pkl`/`.npy`/`.pysxm` and one `- Copy` duplicate, listed as such.
- The review sheets of the first triage pass are in
  `results/data_review/6mer_Cu100_first_pass_20261005/`.
- On Viper, the staged copies are temporary:
  `/ptmp/<user>/stmfit/full146_v1_inputs/full146_raw/` (6-mer cohort) and
  `/ptmp/<user>/stmfit/data_10_20mer` (10-20mers).

The manifests in `benchmarks/` hold external grading metadata only (quality
classes, expected counts, control sequences, triage). No fitting, selection
or assignment code reads them.

A new molecule:

```bash
mkdir data/<molecule>
cp /path/to/scans/*.sxm data/<molecule>/
julia --project=. scripts/calibration/measure_calibration.jl data/<molecule>/<clean_scan>.sxm
```

Then follow `docs/src/calibration.md` (Calibrating a new molecule).
