# Data directory

Raw SXM (Nanonis) scan files live **here**, organized by dataset. They are
**not committed** to git (too large) — obtain them from lab storage or MPCDF
and place them in the subdirectories below.

## Convention

```
data/
├── 240817_chitosan_6mer/      # chitosan/Cu(100), LHe, 6mer benchmark (94 files)
│   └── 240817_001.sxm
│   └── ...
├── 251206_chitosan_10_20mer/  # chitosan/Cu(100), 10–20mer (25 files)
│   └── 260220_083.sxm
│   └── ...
└── <new_molecule>/            # add a subdirectory per new dataset
```

Point the pipeline at a dataset via `STMFIT_DATA_DIR`:

```bash
export STMFIT_DATA_DIR=$(pwd)/data/240817_chitosan_6mer
julia -t 4 --project=. test/batch_full.jl 48 --config config/chitosan.toml
```

## Where to get the data

| Dataset | Source | Files |
|---|---|---|
| `240817_chitosan_6mer` | Lab storage: `/home/durif/Rebecca/data/data/20240817_LHe_Cu100` · MPCDF: `/ptmp/oldu/stmfit/data` | 94 .sxm |
| `251206_chitosan_10_20mer` | Lab storage: `/home/durif/Rebecca/data/10_20mer_analysis` · MPCDF: `/ptmp/oldu/stmfit/data_10_20mer` | 25 .sxm |

The benchmark manifests in `benchmarks/` record which files are clean targets,
poor quality, or excluded — but they contain **no** `.sxm` data, only metadata.

## Adding a new molecule

```bash
mkdir data/<molecule_name>
cp /path/to/scans/*.sxm data/<molecule_name>/
export STMFIT_DATA_DIR=$(pwd)/data/<molecule_name>
julia --project=. test/measure_calibration.jl data/<molecule_name>/<clean_scan>.sxm
```
