# Running STMFit on the MPCDF clusters (Raven, Viper)

HPC use is deliberately simple: sync committed code, instantiate on the login
node, submit one `sbatch`, poll `squeue`, fetch outputs back. No computation on
login nodes (they only run `Pkg.instantiate()` and small merges). Compute nodes
have no internet. `/ptmp` is purged after about 12 weeks, so always fetch results.

| Script | Purpose |
|---|---|
| `run_molecule_consensus.sbatch` | One raw-to-prediction run of the promoted method (4 CPUs, 16 GB, ~1.5 h for 146 scans) |
| `launch_remote.sh` + `batch_array.sbatch` + `merge_chunks.jl` | Push-button `batch_full.jl` counting array: sync → instantiate → array → merge → fetch |
| `remote.env.example` | Template of the personal, ignored `hpc/remote.env` |
| `qe_molds/` | QE structures and input templates of the DFT molds (record only) |

## One-time setup

MPCDF login uses password + OTP (no SSH keys). Add a gateway with connection
sharing to `~/.ssh/config` so the OTP is typed once per session:

```sshconfig
Host gate
    Hostname gate1.mpcdf.mpg.de
    User YOUR_MPCDF_USERNAME
    ControlMaster auto
    ControlPersist 12h
    ControlPath ~/.ssh/master-%C

Host viper raven
    Hostname %h.mpcdf.mpg.de
    User YOUR_MPCDF_USERNAME
    ProxyJump gate
    ControlMaster auto
    ControlPersist 12h
    ControlPath ~/.ssh/master-%C
```

Julia 1.13 on Viper: `/viper/u1/oldu/software/julia-1.13.0/bin/julia`. Code
lives under `/u/<user>/code/`, data and outputs under `/ptmp/<user>/stmfit/`.

## Promoted-method run (`run_molecule_consensus.sbatch`)

Required variables, all absolute paths: `STMFIT_PROJECT_DIR` (committed
checkout), `STMFIT_INPUT_DIR` (contains the raw subdirectory), `STMFIT_OUTDIR`
(must not exist), `JULIA_BIN`. Optional: `STMFIT_RAW_SUBDIR` (default
`full146_raw`), `STMFIT_COUNT_CONFIG`, `STMFIT_ASSIGNMENT_CONFIG`,
`STMFIT_TEMPLATES`, `STMFIT_SELECTED_SUMMARY`. Defaults are the promoted
configs and `templates/chitosan_cc_molds_native_v1.tsv`.

**Viper sets `SBATCH_EXPORT=NONE`.** Pass every variable on the command line:
`sbatch --export=NONE,STMFIT_PROJECT_DIR=...,...`. Variables set in the shell
before `sbatch` are not seen by the job. Always run the script with `--dry-run`
on the login node first (metadata only). The full command sequence is in
the runbook (`docs/src/chitosan_runbook.md`).

## Counting arrays (`launch_remote.sh`)

```bash
cp hpc/remote.env.example hpc/remote.env   # set STMFIT_REMOTE_USER, host, resources
./hpc/launch_remote.sh --dry-run           # always first
./hpc/launch_remote.sh --watch             # sync, instantiate, submit, wait, merge, fetch
./hpc/launch_remote.sh --fetch-only        # fetch a finished job
```

Per-run settings are environment variables: `STMFIT_CONFIG`, `STMFIT_OUTDIR`,
`STMFIT_TSV`, `STMFIT_SKIP_1D`, `N_FILES`, `STMFIT_BATCH_ARGS`. Each array
task runs `batch_full.jl --chunk i/N`. `merge_chunks.jl` concatenates the
shard summaries. Re-submission resumes, because finished files are skipped.

## Limits

- **Raven group quota (account `oldu`): 8 CPUs.** The launcher default of
  4 chunks × 4 CPUs exceeds it; use at most 8 CPUs in total.
- `batch_full.jl` uses at most 4 Julia threads per task; the sbatch scripts set
  `OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = 1`.
- On Viper, do not set `--partition` (the submit filter chooses). On Raven,
  shared jobs need a total `--mem`, not `--mem-per-cpu`.
- Long 10–20mer counting can take 8 h or more.
- QE jobs run sequentially under the one-node QOS limit, never on login nodes.
- Acknowledge MPCDF in publications using these results.
