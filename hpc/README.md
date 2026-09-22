# Running STMFit on the MPCDF HPC cluster

The completed September 22 **registered native refit** uses
`hpc/compare_registered_refit.sbatch`: one two-hour job, four requested CPUs /
16 GB, four one-thread shards, no retry. It reuses the acquisition estimator,
then compares saved support with fresh zero-shift and registered Gaussian/split
fits at saved N. Set the same five absolute-path variables documented below;
run `bash hpc/compare_registered_refit.sbatch --dry-run` before submission and
pass them explicitly through `sbatch --export=ALL,... --no-requeue`.
Sync committed code and the existing Manifest, instantiate on the login node,
compute only in Slurm, then fetch/check all outputs. Source **c7836f9**, job
**11942091**, submitted **20:39:57 CEST**, started **20:42:15** after group CPU quota.
It finishes **21:09:12**, **FAILED 1:0 in 26m57s**: three fixed-N support-span
failures trigger the declared refusal to emit incomplete assignment arms.
All 146 scans are accounted for; 143 complete all four fits. The result is
**inconclusive**, with no candidate grade, not a measured recognition loss.
All **549 outputs and Slurm log** are fetched/checksummed; **22,742** final
checks pass. Fifteen reference tables replay exactly; support **676/34** stays.
Output: `/ptmp/oldu/stmfit/registered_refit_20260922_v1`; do not resubmit it.
Use `test/verify_registered_refit.jl --recorded-failures` for the saved-only
diagnostic, not to authorize partial grading. Full report:
`results/registered_refit_20260922/report.md`.

The completed September 22 acquisition experiment uses
`hpc/compare_acquisition_registration.sbatch`: **one two-hour job, four CPUs /
16 GB, no retry**. Explicitly export `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`
(saved symmetric base/split tables), `STMFIT_INPUT_DIR` (raw scans/templates),
`STMFIT_OUTDIR` (new) and `JULIA_BIN` (1.13) on the `sbatch --export=ALL,...`
command line. Run the script with `--dry-run` first. Up to four one-thread
registration shards precede native replay, observed-mask control and registered
assignment arms; no Gaussian refitting or count selection occurs. Fetch all outputs and
check gates/coverage before external grading. Source **9961973**, job
**11941278**, starts **18:53:18 CEST** on `vipc2227`, after **6m04s** queued on
`QOSGrpCpuLimit`, finishes **19:00:10**, **COMPLETED, 0:0 in 6m52s**. Eight
CPUs allocated, four workers used. All **626 outputs and Slurm log** are fetched
and checksummed. The result is **negative: 671/870 correct, 26/145 exact**,
versus both controls **676/34**, coverage **850 versus 852/870**. Retain saved
support; no retry or tuning. Full report:
`results/acquisition_registration_20260922/report.md`. Do not resubmit
`/ptmp/oldu/stmfit/acquisition_registration_20260922_v1`.

The authorized September 22 matched geometric comparison uses
`hpc/compare_geometry_profile.sbatch`: **one four-hour job, four CPUs / 16 GB**.
Its five explicit exported variables are `STMFIT_PROJECT_DIR`,
`STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`, `STMFIT_OUTDIR` and `JULIA_BIN` (1.13),
with the same cache/raw/template layouts as the amplitude comparison below.
Run `--dry-run` before submitting with explicit `sbatch --export=...`.
Four one-thread shards create common fixed-N native starts and matched
joint/profiled refinements. The job also regenerates the cached support reference
and both unchanged assignment pipelines. Fetch all numerical traces and outputs
before external grading. Source **7c5fb6c**, job **11940220** submitted
**17:24:09 CEST**, starts **17:27:16 CEST** on `vipc2400` after **3m07s**
queued on `QOSGrpCpuLimit`, completes **0:0 in 15m17s**. Eight CPUs allocated;
four workers used. All **782 outputs and Slurm log** are fetched/checksummed.
The comparison is **negative: 632/870 correct, 20/145 exact**, versus same-start
joint and saved support **676/34**, at unchanged **852/870 coverage**. All 146
profiled objectives improve, but 103 searches hit the evaluation limit and
recognition worsens. Retain saved support; no promotion, tuning or retry.
Full report: `results/geometry_profile_20260922/report.md`. The job is closed. Do not resubmit
`/ptmp/oldu/stmfit/geometry_profile_20260922_v1`.

The authorized September 22 fixed-geometry amplitude comparison uses
`hpc/compare_frozen_amplitudes.sbatch`: one two-hour, four-CPU / 16 GB job.
Set `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR` (original symmetric `features.tsv`
and `features_split.tsv`), `STMFIT_INPUT_DIR` (`full146_raw/`, `templates_cc.tsv`),
`STMFIT_OUTDIR` (new), and `JULIA_BIN` (1.13). Run `--dry-run` before submission;
Viper requires an explicit command-line `sbatch --export=...` with all five.
The job regenerates support control, profiles Gaussian amplitudes/background
in four one-thread shards without changing geometry/N, then regenerates the
candidate assignment. Fetch all diagnostics/outputs before separate grading.
Actual **job 11938782** completed **0:0 in 6m17s**, after **5m29s** queued,
source **567cdde**. Four CPUs were requested and eight allocated. All **561
outputs and Slurm log** are fetched/checksummed under
`results/frozen_amplitude_profile_20260922/`. All 146 profiles satisfy KKT;
the control replays saved support exactly. The outcome is negative:
**675/870 correct, 32/145 exact**, versus support **676/34**, with unchanged
**852/870 coverage**, N and geometry. Keep saved support; no promotion or
retry. The job is closed; do not resubmit
`/ptmp/oldu/stmfit/frozen_amplitudes_20260922_v1`.

The authorized September 22 raw-GCV reselection uses
`hpc/compare_gcv_reselection.sbatch`: one four-hour, four-CPU / 16 GB job.
It regenerates the fixed-count support control, performs a four-shard native
GCV sweep, and refits the treatment at its selected counts before assignment.
Set `STMFIT_PROJECT_DIR`, `STMFIT_INPUT_DIR` (with `full146_raw/` and
`templates_cc.tsv`), `STMFIT_SELECTED_SUMMARY` (original cached own-N counts),
`STMFIT_OUTDIR` (new), and `JULIA_BIN` (1.13), all absolute paths. Run
`bash hpc/compare_gcv_reselection.sbatch --dry-run` first, then submit with an
explicit command-line `sbatch --export=...` list of those five variables.
The raw-GCV extractor and its existing early stopping remain unchanged;
neither the batch hybrid guard nor any external label is used. Fetch all
outputs/logs, check coverage, then grade separately.

This comparison completed as **job 11936561**, **0:0 in 1h25m43s**, after
26 seconds queued, source `853ad38`. Viper allocated eight CPUs for four
requested. All **438 outputs and Slurm log** are fetched locally under
`results/gcv_reselection_20260922/`, with matching checksums. The outcome is
negative: **661/870 correct, 32/145 exact**, versus fresh control **675/34**
and saved support **676/34**; exact counts fall **106 → 101/145**. Fresh refits
are not byte-identical even at fixed N; the report quantifies this variability.
Keep saved support, with no promotion or retry. The job is closed; do not
resubmit `/ptmp/oldu/stmfit/gcv_reselection_20260922_v1`.

For the September 20 fixed-count filtering comparison, use
`hpc/compare_fusion.sbatch`. It runs two separate source checkouts sequentially,
regenerating base/split features and all patches in both, then fitting native
assignment on each complete cohort. Its budget is four CPUs, 16 GB, two hours.
Set the six absolute paths documented in that script and run
`bash hpc/compare_fusion.sbatch --dry-run` first. Pass those variables explicitly
in `sbatch --export=...` on Viper (`SBATCH_EXPORT=NONE` otherwise removes them).
Supply the original full146 `selected_from_features.tsv`, not the promoted
counting summary, and fetch both arms and their logs before external grading.
No benchmark labels enter the compute job. Existing completed runs are not
restarted.

This comparison completed as job **11878491** (34m25s, exit 0:0; eight CPUs
allocated by Viper for the four-CPU request). Both complete outputs and the
Slurm log are local under `results/fusion_comparison_20260920/`. The fixed-count
recognition gain is 666 → 671 correct / 870 and 24 → 28 exact chains / 145;
the correction is not promoted. This job is closed; do not resubmit it.

STMFit's batch pipeline (`test/batch_full.jl`) is an *embarrassingly parallel*
sweep over STM image files: each file is fitted independently, and the script
already shards its work list with `--chunk i/n` (round-robin) plus internal
`Threads.@threads` parallelism. That maps directly onto a **Slurm job array** —
one array task per chunk — which is the fastest way to get a large batch done.

This folder adds the three missing pieces to make that painless:

| File | Role |
|---|---|
| `batch_array.sbatch` | Slurm array script (one task = one chunk). Generic for Raven & Viper. |
| `merge_chunks.jl` | Concatenates the per-chunk `summary_*_chunkNNofMM.tsv` into one summary. |
| `launch_remote.sh` | Local push-button launcher: sync code+data → submit → (watch) → merge → fetch. |
| `remote.env` | Your personal config (copied from `remote.env.example`, gitignored). |

> **Why is a merge step needed?** When sharded, `batch_full.jl` writes
> `summary_overlap060_hard_chunk01of04.tsv`, …, `_chunk04of04.tsv` — one per
> chunk. Nothing in the existing code concatenates them. `merge_chunks.jl` does.

---

## 1. Prerequisites (one-time)

### 1a. MPCDF account + 2FA
You need an MPCDF account with **two-factor authentication enabled**. SSH keys
are **not** supported for login — you authenticate with password + OTP. See the
MPCDF [2FA FAQ](https://docs.mpcdf.mpg.de/faq/2fa.html).

### 1b. SSH config (ProxyJump + ControlMaster)
Access goes through a gateway (`gate1.mpcdf.mpg.de` / `gate2.mpcdf.mpg.de`).
Add this to `~/.ssh/config` so you only type password+OTP **once** per session:

```sshconfig
Host gate
    Hostname gate1.mpcdf.mpg.de
    User YOUR_MPCDF_USERNAME
    ServerAliveInterval 120
    ControlMaster auto
    ControlPersist 12h
    ControlPath ~/.ssh/master-%C

Host raven
    Hostname raven.mpcdf.mpg.de
    User YOUR_MPCDF_USERNAME
    ProxyJump gate
    ControlMaster auto
    ControlPersist 12h
    ControlPath ~/.ssh/master-%C

Host viper
    Hostname viper.mpcdf.mpg.de
    User YOUR_MPCDF_USERNAME
    ProxyJump gate
    ControlMaster auto
    ControlPersist 12h
    ControlPath ~/.ssh/master-%C
```

Now `ssh raven` (or `ssh viper`) connects through the gateway and reuses the
master connection for all subsequent `rsync`/`ssh` in the launcher. Source:
[MPCDF connecting guide](https://docs.mpcdf.mpg.de/faq/connecting.html).

### 1c. Pick a system & find Julia
On the login node, check the available Julia module:

```bash
ssh raven          # or viper
find-module julia  # list versions, e.g. julia/1.11.4
module load julia  # or julia/1.11.4 to pin
julia --version
```

- **Raven** — mature Intel CPU cluster. You may set `--partition=general`
  (uncomment the line in `batch_array.sbatch`).
- **Viper** — newer AMD CPU system. **Do not** set `--partition`: a submit
  filter picks the queue automatically from your resource request.
  ([Viper-CPU guide](https://docs.mpcdf.mpg.de/doc/computing/viper-user-guide.html))

### 1d. Where things live (MPCDF convention — paths are derived for you)
| What | Where | Why |
|---|---|---|
| Code (`STMFit/`) | `/u/<user>/code/STMFit` | Permanent home, backed up, quota'd (~1–2.5 TB) |
| `.sxm` data | `/ptmp/<user>/stmfit/data` | Fast scratch for batch I/O |
| Results | `/ptmp/<user>/...` then **fetch back** | `/ptmp` is auto-cleaned after ~12 weeks |

The launcher builds these paths automatically from your **`STMFIT_REMOTE_USER`**
(short MPCDF account name). You don't type full paths anywhere. Override the
bases only in the "Advanced" section of `remote.env` if your layout differs.

> ⚠ `/ptmp` auto-cleans: always `rsync` results back to your machine (the
> launcher does this for you).

---

## 2. Configure

```bash
cp hpc/remote.env.example hpc/remote.env
$EDITOR hpc/remote.env      # the only thing you MUST set is STMFIT_REMOTE_USER
```

Key knobs:

| Variable | Meaning | Default |
|---|---|---|
| `STMFIT_SSH_HOST` | `~/.ssh/config` alias (`raven`/`viper`) | `raven` |
| `STMFIT_REMOTE_USER` | **Short MPCDF account name** for paths (e.g. `yourname`). NOT your SSH `User` (which may be your email for 2FA). — |
| `N_CHUNKS` | Array size / sharding denominator | `4` |
| `CPUS_PER_TASK` | Julia threads per task (code caps at 4) | `4` |
| `WALLTIME` | `HH:MM:SS` | `04:00:00` |
| `MEM_PER_CPU` | MB per CPU | `4000` |
| `JULIA_MODULE_VERSION` | Pin e.g. `1.11.4`; empty = default | — |
| `STMFIT_CONFIG` | TOML config relative to repo | `config/chitosan.toml` |
| `STMFIT_OUTDIR` | Output dir relative to repo | `results/best_plots` |
| `STMFIT_TSV` | Optional triage/input TSV relative to repo | — |
| `STMFIT_SKIP_1D` | `1` to pass `--skip-1d` to `batch_full.jl` | — |
| `STMFIT_BATCH_ARGS` | Extra flags forwarded to `batch_full.jl` | — |
| `N_FILES` | Limit to first N files; empty = all | — |
| `STMFIT_MAIL_USER` | Email for Slurm notifications | — |
| `SSH_CONNECT_TIMEOUT` | SSH connect/banner timeout, useful for slow gateway/password+OTP flows | `180` |
| `SSH_SERVER_ALIVE_INTERVAL` | SSH keepalive interval while commands run | `60` |

`remote.env` is gitignored — it holds your username and paths.

---

## 3. Run (push-button)

```bash
# Preview without changing anything:
./hpc/launch_remote.sh --dry-run

# Sync + instantiate + submit, then exit:
./hpc/launch_remote.sh

# Submit and block until done, then merge + fetch results:
./hpc/launch_remote.sh --watch

# Only fetch results from a job you already submitted:
./hpc/launch_remote.sh --fetch-only
```

The launcher (`--watch`) does, in order:

1. **rsync code** → cluster (`--delete`, excludes `results/`, `.git/`, logs, `Manifest.toml`).
2. **rsync `.sxm`** → `/ptmp` (only `.sxm` files, `--update` so unchanged files skip).
3. **`Pkg.instantiate()`** on the login node (downloads packages — needs internet, only login nodes have it).
4. **`sbatch --array=1-N`** → captures the job id.
5. **polls `squeue`** until the array finishes.
6. **`merge_chunks.jl`** on the login node → one combined summary.
7. **rsync results** back to your local `results/`.

Re-running is **safe and resumable**: `batch_full.jl` skips files already marked
`ok` with an existing plot, so a re-submitted array only finishes the remainder.

The remote batch automatically receives `STMFIT_DATA_DIR=/ptmp/<user>/stmfit/data`
after data sync, so local paths never need to be embedded in Slurm jobs. Prefer
the dedicated `STMFIT_TSV`, `STMFIT_SKIP_1D`, `STMFIT_SELECTION_POLICY`,
`STMFIT_REFINED_ADVISORY`, and `STMFIT_PLOT_MANIFEST` variables over multi-word
`STMFIT_BATCH_ARGS` when possible.

---

## 4. Manual submission (without the launcher)

```bash
ssh raven
cd /u/<user>/code/STMFit
# (first time) julia --project=. -e 'using Pkg; Pkg.instantiate()'

sbatch --export=ALL,N_CHUNKS=8,STMFIT_CONFIG=config/chitosan.toml,STMFIT_OUTDIR=results/best_plots \
       --array=1-8 --cpus-per-task=4 --time=04:00:00 --mem-per-cpu=4000 \
       hpc/batch_array.sbatch

squeue -u $USER          # monitor
scancel <jobid>          # cancel
# after it finishes:
module load julia && julia --project=. hpc/merge_chunks.jl results/best_plots --total 8
```

Then `rsync` the `OUTDIR` back to your machine.

---

## 5. How it fits together

```
                  ┌─ task 1 → batch_full.jl --chunk 1/N ─→ summary_*_chunk01ofN.tsv
sbatch --array=1-N├─ task 2 → batch_full.jl --chunk 2/N ─→ summary_*_chunk02ofN.tsv
   (N parallel)   └─ ...                                      ↓
                                              merge_chunks.jl
                                                      ↓
                                       summary_overlap060_hard.tsv  ← (single, merged)
```

Each task also writes `<file>_best.png` plots and per-file subdirs into `OUTDIR`.
Because sharding is round-robin and per-file outputs are keyed by filename, the
chunks never collide — only the summary TSV needs merging.

---

## 6. Tuning guide

- **How many chunks?** Start at `N_CHUNKS=4–8`. Estimate per-file time locally
  (`julia -t 4 --project=. test/inspect_one_file.jl file.sxm`) and pick `N` so
  `total_files × per_file_time / (N × cpus)` comfortably fits under `WALLTIME`.
- **MPCDF array limit:** you can submit up to ~300 array tasks at once. If you
  need more chunks, split into batches (the launcher is idempotent).
- **CPUs per task:** the code caps Julia threads at `min(4, …)`, so
  `CPUS_PER_TASK=4` is optimal. Don't raise it unless you also edit that cap in
  `test/batch_full.jl:996`.
- **Memory:** `4000 MB/CPU` (16 GB/task at 4 CPUs) is ample for 1024×1024 SXM
  images. Lower it if your images are small; raise for very large scans.
- **BLAS oversubscription:** the sbatch sets `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`.
  Each file is a small least-squares fit; multi-threaded BLAS would fight the
  Julia worker threads.

---

## 7. Troubleshooting

| Symptom | Fix |
|---|---|
| `merge_chunks.jl` reports missing chunks | Some array tasks failed. Re-submit only those: `sbatch --array=2,5 --export=ALL,N_CHUNKS=8 …`, then re-run merge. |
| `Pkg.instantiate()` fails on login node | Manifest drift. Sync is `--exclude=Manifest.toml`; instantiate regenerates it. If it still fails, check your network (login nodes only) or precompile locally with the same Julia version. |
| Job sits in `PD` (pending) forever | Cluster busy / partition full. Smaller `WALLTIME` or fewer CPUs queues faster; or try the other system. |
| `sbatch: error: invalid partition` | You set `--partition` on Viper. Remove it — the submit filter chooses. |
| Password/OTP asked many times | Your `~/.ssh/config` lacks `ControlMaster`/`ProxyJump`. See §1b. |
| Results missing after fetch | Check the cluster: `ssh raven 'ls $STMFIT_OUTDIR'` and the per-task logs in `results/hpc_logs/stmfit_<jobid>_*.err`. |
| `/ptmp` data vanished | Auto-cleaning kicked in (~12 weeks). Re-sync data; results were fetched, so only inputs need restoring. |

---

## 8. QE Mold Jobs

The DFT-STM mold workflow uses separate QE launchers, not the STM image batch
array launcher above.

After preparing `qe/glcn` and `qe/glcnac` locally, run:

```bash
bash hpc/launch_qe_molds_remote.sh --dry-run
bash hpc/launch_qe_molds_remote.sh --watch
```

If the MPCDF gateway is slow or you need more time around password/OTP prompts,
raise the SSH timeout in `hpc/remote.env`:

```bash
SSH_CONNECT_TIMEOUT=300
SSH_SERVER_ALIVE_INTERVAL=60
```

On the cluster itself, from the synced repository root, use:

```bash
bash hpc/submit_qe_molds.sh --watch --sequential
```

Both paths run `test/preflight_qe_mold_inputs.jl` first. `--sequential` submits
the two run directories (`qe/glcn`, `qe/glcnac`) as an `afterok` chain and
enforces the task budget as the maximum simultaneous count (8), so two 8-task
jobs fit within a single-node QOS group limit. QE run directories are
gitignored because they can
contain large scratch and cube files.

---

## 9. Joint-proxy inference array

For label-free joint-proxy inference over a large SXM folder, create a text file
containing one bare `.sxm` basename per line, then submit the eight-way array:

```bash
sbatch --array=1-8 --export=ALL,N_CHUNKS=8,STMFIT_FILES_FROM=hpc/joint_proxy_files_all.txt \
  hpc/joint_proxy_inference_array.sbatch
```

Each task writes and validates a separate shard. After every array task has
finished, submit the merge separately so a pending dependency does not consume
the Raven/Viper group quota:

```bash
sbatch --export=ALL,N_CHUNKS=8 hpc/joint_proxy_inference_merge.sbatch
```

The array copies the configured calibration into each shard under the stable
name `joint_proxy_calibration.toml`; this name is part of the manifest/merge
contract. Scans with no finite count candidate are retained as summary-only
abstentions, and non-finite type evidence becomes `0.5/0.5` with `?` rather than
terminating a shard. This workflow does not read benchmark labels.

## 10. Limits & rules to respect

- **No compute on login nodes** — they're shared and resource-limited. Only run
  `Pkg.instantiate()` and `merge_chunks.jl` (both light) there; never the batch.
- **No internet on compute nodes** — that's why instantiate happens on the login
  node before submission.
- **Acknowledge MPCDF** in publications using these results (see
  [MPCDF help](https://docs.mpcdf.mpg.de/faq/help.html)).

## References
- [Connecting to MPCDF Systems](https://docs.mpcdf.mpg.de/faq/connecting.html)
- [Raven User Guide](https://docs.mpcdf.mpg.de/doc/computing/raven-user-guide.html)
- [Viper-CPU User Guide](https://docs.mpcdf.mpg.de/doc/computing/viper-user-guide.html)
- [Environment Modules](https://docs.mpcdf.mpg.de/doc/computing/software/environment-modules.html)
- [HPC Software FAQ (`find-module`)](https://docs.mpcdf.mpg.de/faq/hpc_software.html)


## Native reconstruction job (Julia 1.13)

`reconstructed_chitosan.sbatch` is the original fresh application/comparison job
for `cc_soft_reconstructed_v1`. After the failure of job 30271312, use the
unknown25 continuation below rather than repeating that fresh run. The original
script runs the 25 unknown chains first, then the saved 146-file geometry cohort. Both calls use the native production CLI; neither
reads benchmark labels nor invokes a grader. At most four single-thread refit
chunks run inside the one allocation. Requested resources are four CPUs,
16,000 MB and 24 hours (Raven's test-only accounting reports eight processors).
Do not submit a second concurrent job under the observed eight-CPU group quota.

Prepare a normal code sync in `~/code/STMFit`, excluding agent state and generated
outputs. The input directory must contain:

- `unknown_raw/`: the 25 raw unknown SXMs, with unique basenames;
- `full146_raw/`: raw SXMs matching all files in the saved geometry table;
- `base_geometry_full146.tsv`: label-free selected-N base geometry, not truth;
- `templates_cc.tsv`: templates from the fixed native CC builder.

Use a verified Julia 1.13 executable, not the cluster's Julia 1.12 module. Only
package setup/precompilation and path/key dry-runs run on the login node. For a
fresh depot with the older root Manifest, load Pkg from the stdlib environment
before activating the project and instantiating dependencies. On Viper, the
Julia 1.12-era dependency graph missed a new stdlib dependency on `Zstd_jll`:
`Pkg.precompile()` could return zero while logging cache errors. A targeted
load of `Zstd_jll`, then `Pkg`, then the production packages under the real root
project prepared those caches. A **new normal four-thread process**, without
those preload commands, then imported `GaussianFit2D`, `STMMolecularFit`,
`Plots` and `Clustering` cleanly, including `Plots.FileIOExt`. Check logs and
fresh imports, not only the bootstrap exit code. Never delete caches broadly,
hand-edit the Manifest or fall back to a different Julia version for this job.

From the remote project directory, set the ordinary job environment:

```bash
export STMFIT_PROJECT_DIR="$HOME/code/STMFit"
export JULIA_BIN="$HOME/software/julia-1.13.0/bin/julia"
export STMFIT_INPUT_DIR="/ptmp/$USER/stmfit/reconstructed_cc_soft_v1_inputs"
export STMFIT_OUTDIR="results/reconstructed_cc_soft_v1/my_new_run"
mkdir -p results/reconstructed_cc_soft_v1
sbatch --test-only --export=ALL hpc/reconstructed_chitosan.sbatch
# After inspecting the dry-run and verifying the environment:
sbatch --export=ALL hpc/reconstructed_chitosan.sbatch
```

`STMFIT_OUTDIR/unknown25` and `STMFIT_OUTDIR/full146` must not already exist.
Monitor with `squeue`/`sacct`, then fetch logs and both output directories. Run
`test/compare_reconstructed_champion.jl` and the existing external grader locally
only after production completes. No cluster-side agent, approval artifact,
or persistent job controller is involved. See
[`docs/src/unit_assignment.md`](../docs/src/unit_assignment.md) for scientific
limits and cache semantics. In particular, a one-chain cohort can be too small
for Fisher PCA10; its explicit `?` output is not a failed count or a valid
chemical benchmark.

### Unknown25 continuation after the support handoff repair

`resume_reconstructed_unknown25.sbatch` completed as **Raven job 30278010** on
September 17, from tested source `aa422a3`: **COMPLETED, exit 0:0**. It ran on
`ravc4043` from 20:55:39 to 21:01:55 CEST (6 min 16 s). Four CPUs, 16,000 MB and
24 hours were requested; accounting records **eight allocated CPUs** and
MaxRSS **3,876,696 KiB (about 3.70 GiB)**. No second allocation was started.
**Do not repeat the submission for this completed run.** Only unknown25 ran;
full146 has since completed separately as job 11786116, including local external
evaluation. See the full146-only section below.

All 25 files / 222 lobe keys are present through the application. Local Julia
1.13 integrity verification passes 179 assertions; prediction validation passes
and saved-component vote/QC replay is byte-identical. All 67 fetched file hashes
match Raven (66 outputs plus Slurm log), with `raw_inputs/` excluded at every
level. The local output root is
`results/reconstructed_cc_soft_v1/raven_20260917_unknown25_resume_v1/`.
The monitor is stopped. The output is not chemical validation: 19/25 chains are
flagged for low mean confidence. See the unit-assignment documentation for the
model-output counts and limitations.

The prepared local input directory is
`results/reconstructed_cc_soft_v1/unknown25_resume_v1_inputs/`:

- `unknown25_selected_summary.tsv`: byte-identical original fresh count summary,
  retaining the adaptive `refined_policy`; an N-only replacement is not valid;
- `unknown25_base_features.tsv`: all 222 base rows for the 25 selected files;
- `unknown25_split_partial.tsv`: nine genuine split rows for each of
  `260215_022.sxm` and `260220_083.sxm`;
- `templates_cc.tsv`: unchanged native CC templates;
- `unknown_raw/`: the same 25 raw SXMs, currently local symlinks.

Sync committed source separately, using tracked files only. Stage the input
directory at a new Raven input path. **Dereference the local raw symlinks** when
copying (`rsync -aL`), so they do not point back to local `/home/...` paths on the
cluster. Inspect an rsync dry-run first. Do not copy agent state, use `--delete`,
overwrite the failed run, or rely on the personal launcher's default host; this
continuation targets **Raven explicitly**.

The procedure used for this run is shown below for reproducibility. Job
30278010 has completed; do not execute these submission commands again. Any
later run needs a separate user decision and a different unused output root.
From the synced project on Raven, after checking its Julia 1.13 environment:

```bash
cd "$HOME/code/STMFit"
export STMFIT_PROJECT_DIR="$PWD"
export JULIA_BIN="$HOME/software/julia-1.13.0/bin/julia"
export STMFIT_INPUT_DIR="/ptmp/$USER/stmfit/unknown25_resume_v1_inputs"
export STMFIT_OUTDIR="results/reconstructed_cc_soft_v1/raven_20260917_unknown25_resume_v1"
# Only the shared log parent may be created; keep STMFIT_OUTDIR absent.
mkdir -p results/reconstructed_cc_soft_v1
bash hpc/resume_reconstructed_unknown25.sbatch --dry-run
sbatch --test-only --export=ALL hpc/resume_reconstructed_unknown25.sbatch
# Historical submission for job 30278010; DO NOT execute again:
sbatch --export=ALL hpc/resume_reconstructed_unknown25.sbatch
```

The script rejects an existing output root, duplicate/gapped/extra keys,
incomplete cached files, mismatched headers and raw-name coverage. It reuses
the existing selected-support parser rather than inventing a policy from N.
The metadata dry-run does not open SXM pixels, validate template numerics or fit
anything; it creates no outputs. Its successful result is not evidence that
remaining scientific stages succeed. Keep scientific execution inside Slurm.

Normal execution excludes the two complete cached files from split extraction.
Only the remaining 23 files / 204 lobes are fitted, using the frozen split config
and original selected summary. At most four one-thread children run in the
allocation. New split coverage must be complete before merging with the 18
cached rows. The native runner then receives the complete base/split tables and
original summary; it computes patches, descriptor, cohort-wide scores,
assignments, QC and maps on all 25 files. No count/base refit, full146 execution
or grader is invoked by this script.

Split tables and their logs remain directly under `STMFIT_OUTDIR`; application
outputs go under `STMFIT_OUTDIR/unknown25/`, with the outer application log in
`STMFIT_OUTDIR/logs/unknown25.log`. Slurm writes
`results/reconstructed_cc_soft_v1/unknown25-resume-<jobid>.log`. For a later
separately approved run, poll its job ID and fetch logs/outputs after it ends.
Use `sacct` independently if `squeue` no longer has the job: an invalid queue ID
is not evidence of failure and must not suppress the accounting query. Coverage
failure stops the job; do not omit files, lower N or change settings to make it
pass. The failed job 30271312 and its artifacts remain unchanged.

Local regression (synthetic fixtures and a fake exporter; no scientific fit):

```bash
julia -t 4 --project=. test/test_resume_reconstructed_unknown25.jl
```

This tests the current script payload and actual shell dry-run/collision paths.
It does not replace the completed-run integrity checks above or establish the
chemical accuracy of any predictions.

### Full146-only fixed-cache reconstruction

The user approved this next application after unknown25 was completed and
reviewed. Use `reconstructed_full146.sbatch`, **not** the earlier two-cohort
script (which would rerun unknown25). The new entrypoint requests one node/task,
four CPUs, 16,000 MB and 24 hours. It does not hardcode a partition; Viper is the
host for this multi-file fit. Authentication and Julia setup are now resolved.
**Job 11786116 was submitted once on September 18 at 00:04:54 CEST**, from
`00528d2`. After waiting for group CPU quota, it completed **`0:0`**, from
**00:36:57 to 00:48:01 CEST (11m04s)** on `vipc2162`, in `small`. Four CPUs,
16,000 MB and 24 hours were requested; accounting confirms **eight allocated
CPUs** and **MaxRSS 3,486,576 KiB (about 3.33 GiB)**. Test-only ID 11786052 was
not a scientific job. The monitor is cancelled. **Do not repeat the submission.**
Outputs, local validation, separate comparison and own-N grading are complete.

The existing `/u/oldu/code/STMFit` tree was left untouched. This run uses:

- project: `/u/oldu/code/STMFit_full146_v1`;
- Julia: `/u/oldu/software/julia-1.13.0/bin/julia`;
- inputs: `/ptmp/oldu/stmfit/full146_v1_inputs`;
- direct output: `results/reconstructed_cc_soft_v1/viper_20260917_full146_v1`
  under that project (the name retains the preparation date);
- Slurm log: `results/reconstructed_cc_soft_v1/full146-11786116.log`.

At staging, all 354 source/lock files and 148 inputs matched the frozen local
snapshot's SHA-256 values; four key Julia runtime files also matched. Raw files are regular, not symlinks.
The actual Julia 1.13 metadata dry-run passed without output creation. Setup
and preflight logs are local under
`results/reconstructed_cc_soft_v1/full146_launch_v1/`.

All **184 outputs (26,642,815 bytes)** plus Slurm log were fetched with
`raw_inputs/` excluded at every level; all remote/local output hashes match.
The remote job source, key runtime files and prepared inputs stayed unchanged. The original direct output
has **900 keys / 146 files**, final **695 zero / 198 one / 7 `?`**, and 43 QC
review chains. Native Julia 1.13 saved-output checks pass **73/73**. No fit,
classifier or additional cluster job was run after completion.

A separate local plot-only correction (`6fa3297`, **267/267** regression checks)
repairs the crowded overview, preserving every standalone PNG byte-for-byte.
The readable 19-page figure set and index are under
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_figures_v1/`;
original job images are retained. This changes only presentation, not assignments.

External-only reports are under
`results/reconstructed_cc_soft_v1/viper_20260917_full146_v1_external/`:
**846/892** labels agree with the frozen reference (zero missing reference keys,
eight extra keys outside the confirmed145 manifest). The full145 own-N grade is
**666/849 (78.4%) / 24 exact chains**, with **16 missing control positions** and
**38 extra predicted lobes**. It is below the frozen reference, **677/854 (79.3%) /
36 exact chains**; no promotion or benchmark-guided retuning follows.

The scientific inputs are fixed at **146 scans / 900 cached GCV lobe keys**.
This is an assignment reconstruction at the cache's own N, not a replay of the
promoted batch count policy. The `pm2_confirm` summary has 871 keys and conflicts
on 28 scans; it must not be supplied alongside this base table. No compatible
full146 split cache is available, so the native runner computes all genuine
split fits in at most four single-thread children. Classifiers then use the
whole cohort. No count/base rerun, unknown25 stage, comparison or grader is
called by this entrypoint.

Stage `results/reconstructed_cc_soft_v1/full146_v1_inputs/` at a new cluster
input directory. It contains exactly:

- `base_geometry_full146.tsv`: unchanged 900-row base geometry;
- `templates_cc.tsv`: unchanged native CC templates;
- `full146_raw/`: all 146 scans, local symlinks that must be dereferenced during
  transfer (`rsync -aL`, after inspecting a dry-run).

Keep source sync limited to tracked runtime files; do not copy agent/archive
state, external truth or reference predictions. Do not use `--delete`, overwrite
previous runs or alter the Manifest. Verify Julia 1.13, source/input hashes and
raw file coverage on the target. No scientific computation runs on login nodes.

`STMFIT_OUTDIR` is the **new direct full146 output directory**, with no implicit
`/full146` suffix. `STMFIT_SELECTED_SUMMARY` is optional, explicit and empty for
this cache; the native runner honestly derives `selected_from_features.tsv`.
An explicitly supplied summary must be compatible; nothing is auto-discovered.
This nonadaptive case does not waive the original-summary requirement for
adaptive-support data.

Invocation reference for a **separately approved new run**, not a replay of
completed job 11786116. After staging and checking the environment:

```bash
# Set absolute STMFIT_PROJECT_DIR, JULIA_BIN, STMFIT_INPUT_DIR and STMFIT_OUTDIR.
# Leave the output directory absent; create only the Slurm log parent.
export STMFIT_SELECTED_SUMMARY=""
bash hpc/reconstructed_full146.sbatch --dry-run
sbatch --test-only --export=ALL hpc/reconstructed_full146.sbatch
# One actual submission only, after successful preflight and scheduler check:
sbatch --export=ALL hpc/reconstructed_full146.sbatch
```

Poll that new job ID, fetch its outputs/logs locally with `raw_inputs/` excluded
at every level, and check all 146/900 keys. Failures retain outputs; do not drop
files, change N or retry with tuned settings. Comparison and the full145 own-N
grade run locally only after integrity validation. Neither grade nor reference
predictions feed back into production.

Focused local checks (synthetic metadata, no scientific fitting):

```bash
julia --project=. test/test_reconstructed_full146_sbatch.jl
julia -t 4 --project=. test/test_reconstructed_pipeline.jl
julia --project=. test/test_compare_reconstructed_champion.jl
```

The entrypoint regression passes 122 assertions (70 shell/capture, 52 native
metadata/boundary checks). The existing pipeline and comparison fixtures pass
19 and nine assertions. Prepared-input metadata verification passes 14 checks;
the actual entrypoint dry-run reports 146 files and creates no output root.
These checks do not validate SXM content, scientific fits or chemical accuracy.


## Bounded label-free diagnostic job (2026-09-18)

This opt-in research pass is separate from the completed reconstruction jobs.
Do not resubmit those jobs. `hpc/label_free_exploration.sbatch` requests four
CPUs, 16 GB and two hours. It runs four independent one-thread acquisition/count
cases, then two sequential fixed native Fisher replays. It does not select a
production N, predict chemical classes, read benchmark labels or grade results.

Use an isolated source tree and exclusive input/output paths. The input layout
is `raw/` (four approved scans), `selected_summary.tsv`, `base_geometry.tsv`,
`candidate_counts.tsv`, and
`fisher/{unknown25,full146}/{patches_fwd17.tsv,fisher_cv.tsv}` plus
`fisher/unit_assignment_reconstructed.toml` copied unchanged. Fisher uses
forward residual patches/prefix `res`, not the backward descriptor patches.
Dereference raw symlinks during transfer; do not stage truth, reference chemical
predictions, other datasets or agent state. Never use `rsync --delete`.

```bash
# After source/input synchronization and hash checks; no fitting on login nodes:
export STMFIT_PROJECT_DIR=/u/oldu/code/STMFit_labelfree_exploration_20260918
export STMFIT_INPUT_DIR=/ptmp/oldu/stmfit/label_free_exploration_20260918_inputs_v1
# Placeholder for a future approved run; do not reuse the failed run's bundle.
export STMFIT_OUTDIR=/ptmp/oldu/stmfit/label_free_exploration_20260918_NEW_RUN
export JULIA_BIN=/u/oldu/software/julia-1.13.0/bin/julia
cd "$STMFIT_PROJECT_DIR"
bash hpc/label_free_exploration.sbatch --dry-run
sbatch --test-only --export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN hpc/label_free_exploration.sbatch
# The first pass already used its one-job allowance (11806180 failed at startup).
# A future submission needs renewed approval and a fresh output path:
sbatch --export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN hpc/label_free_exploration.sbatch
```

Viper currently sets `SBATCH_EXPORT=NONE`. Merely assigning the four variables
before `sbatch` does not export them to the job. The explicit command-line
`--export=...` above overrides that default and passes only the required paths.
A `#SBATCH --export` header alone would lose to `SBATCH_EXPORT`; command-line
options take precedence. Job 11806180 failed at the first required-path check
because its submission omitted this flag. Julia and all scientific stages never
started. The original source/input hashes were unchanged. No automatic retry was submitted.
Shell dry-runs and `sbatch --test-only` do not execute the compute-node environment.
After explicit user approval, corrected job **11812202** used this export list
and completed 0:0 on 2026-09-18 (17:10:22–17:14:12 CEST). Its 73 outputs and Slurm
log were fetched and hash-verified; all 37 source and 12 input hashes stayed
unchanged. Saved-output checks pass: counting 10,365 assertions, Fisher 33,461,
and acquisition 255 grouped checks. All 78 fit rows are retained, including four
invalid native/fixed rows; 25/26 refinements remain evaluation-limited. Both
Fisher score files are byte-identical to the saved originals. Acquisition
registration/noise claims remain blocked; process success is not scientific
validation. See the September 18 journal for measured findings and limits.
The new one-submission allowance is consumed. Do not repeat or requeue the job;
no further scientific work or method promotion is authorized by this pass.

Site account/partition options must match the current Viper allocation. The
metadata-only dry-run reads no SXM pixels, performs no Fisher refit and creates
no result directory. For real execution the native driver requires Slurm and
at least four requested CPUs; useful concurrency stays four even if eight CPUs
are allocated. OpenBLAS/OMP and each Julia process use one thread.

`stages.tsv` records each completed case/stage and exit code immediately.
Per-stage logs, infeasible/failed candidate rows and incomplete outputs remain
available. A failed stage yields a failed job after the other stages finish;
fetch and diagnose it rather than lowering N, changing bounds or submitting a
second job. Fisher export mismatches are scientific output, not grounds to tune
or retry. Use independent `squeue` and `sacct` checks, then fetch the fresh output
directory and Slurm log. No job is submitted merely by installing these scripts.


## Masked-preprocessing four-scan diagnostic

The opt-in pilot uses `hpc/masked_preprocessing.sbatch`, not the previous
exploration or reconstruction launchers. It requests one task, four CPUs,
16000 MB and 30 minutes. Four one-thread Julia 1.13 child commands compare
preprocessing on the same four frozen scans. No count, chemical classifier,
feature-extraction or grading stage runs. Failures are retained in per-case logs
and `stages.tsv`; a failed child makes the batch fail after all four attempts.
Numerical `PARTIAL`/`NONCONVERGED` results remain scientific output, not process
success claims.

Required variables are `STMFIT_PROJECT_DIR`, `STMFIT_INPUT_DIR`, `STMFIT_OUTDIR`
and `JULIA_BIN`. Under Viper's `SBATCH_EXPORT=NONE`, the **actual submit command**
must include:

```text
--export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN
```

Setting variables before `sbatch` or placing an export directive in the script
is not sufficient. Run the script with `--dry-run` first; it prints the four
commands after metadata checks and neither parses raw images nor creates
outputs. Check each per-case CLI's metadata-only dry-run as well. Scheduler
`--test-only` is useful but does not prove that compute-node variables arrive.
Import/precompile checks may run on the login node; real preprocessing must run
inside the one Slurm allocation.

Use a fresh output path with an existing parent and no symlink ancestors. Viper's
`/ptmp` is a symlink, so resolve the intended output parent first: for this pilot
`/ptmp/oldu/stmfit/...` resolves to `/viper/ptmp1/oldu/stmfit/...`. These are the
same scratch location, not two result directories. Keep source/input paths and
hashes, fetch arrays/tables and Slurm logs locally, and use a separate local
folder for saved-only checks and figures. Never overwrite failed or prior runs.
The script does not resubmit itself. The authorized pilot was submitted exactly
once as job **11820414** on 2026-09-18 at 22:25:06 CEST, from scientific
source `cb3b2d9843eccac2612859fe8b19cfe609fdd074`. It completed **0:0** on
`vipc2169`, 23:32:39–23:33:34 CEST (**55 s**). Slurm allocated eight CPUs for
four requested; useful concurrency remained four one-thread commands. All four
stages exited zero. The 37 outputs and Slurm log were fetched once and their
remote/local hashes match; all 33 source and seven input identities stayed
unchanged. The monitor is stopped and the one-job allowance is consumed.
**Never resubmit or requeue this completed pilot.**

Saved-output checking retains **30,811 passes / 15 strict failures**: only the
near-zero stationarity scalar differs after physical-coefficient export. A
separate 118/118 arithmetic check reproduces the export-rounding mechanism and
bounds all 24 view discrepancies; it does **not** recover exact stationarity
replay or exclude every small error. The original checker/tolerance and failed
log remain unchanged. All other plane/row/objective/mask/metric/native-baseline
checks pass. Twenty saved-array figures were generated, with all four common
comparisons inspected. Background spread decreases, but forward/backward RMS
changes are mixed; no observed-pixel loss is added and the incomplete scan stays
`PARTIAL`. Every diagnostic best lag is boundary-limited and ambiguous. No
count/classifier/grade, calibrated noise claim, method selection or production
promotion follows. See the September 18 journal and local report
`results/masked_preprocessing_20260918/real4_v1_report.md` for the measured results
and verification limits.

## One Fisher mirror comparison (September 20)

`hpc/compare_fisher_mirror.sbatch` runs the native assignment pipeline twice on
the same complete saved symmetric-fusion geometry/patches. The only changed
scientific setting is `fisher_layout` (legacy physical-t vs physical-u mirror).
The job requests four CPUs, 16 GB and 30 minutes, without counting or geometry
refits. It does not accept labels or perform grading.

Set absolute `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`
(`full146_raw/`, `templates_cc.tsv`), `STMFIT_OUTDIR` (new) and `JULIA_BIN`.
Run `bash hpc/compare_fisher_mirror.sbatch --dry-run` first, then pass all five
variables explicitly with `sbatch --export=...`. Outputs are `control/` and
`transverse/` plus their logs; failures are retained and the job exits nonzero.
The existing configs/defaults and completed fusion job stay unchanged.

Completed job **11879515** (source `8bb29f4`) exited 0:0 in 2m03s. Its 358
outputs and Slurm log are fetched and checksum-verified under
`results/fisher_mirror_comparison_20260920/`. The transverse candidate regresses
671/28 → 669/27 (correct positions / exact chains); it is not retained. Do not
resubmit the completed job. The report includes every paired gain and loss.

## One matched-residual comparison (September 20)

`hpc/compare_matched_residual.sbatch` requests four CPUs, 16 GB and 30 minutes.
It runs `control/` (`S(data)-model`) and `matched/` (`S(data-model)`) sequentially.
Both arms regenerate all three patch families from the same complete saved
symmetric-fusion geometry; `STMFIT_CACHE_DIR` supplies only `features.tsv` and
`features_split.tsv`. No count or geometry fit, label input or grading runs.

Set absolute `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`
(`full146_raw/`, `templates_cc.tsv`), `STMFIT_OUTDIR` (new, existing parent) and
`JULIA_BIN` (Julia 1.13). First run
`bash hpc/compare_matched_residual.sbatch --dry-run`, then submit with the explicit
command-line option
`--export=STMFIT_PROJECT_DIR,STMFIT_CACHE_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN`.
This is necessary with Viper's `SBATCH_EXPORT=NONE`. Actual execution requires
a Slurm allocation. Logs and failed outputs are preserved; either failed arm
makes the job fail after both have been attempted. Fetch outputs and grade
externally only after checking the regenerated control against its saved output.

Completed job **11883607** (source **dae6395**) exited 0:0 in **5m15s** on
`vipm2063`, with four CPUs requested and eight allocated. All 390 regular output
files and the Slurm log are fetched and checksum-verified under
`results/matched_residual_comparison_20260920/`; raw-input staging links are not
duplicated. The regenerated control reproduces all twelve checked tables
byte-for-byte. Matched filtering gains **671/28 → 672/29** at unchanged coverage;
it is the opt-in working candidate, not a promoted champion. No new tuning or
submission follows this completed pass. **Do not resubmit this job.**

## Two transverse-descriptor candidates (September 20)

`hpc/compare_transverse_descriptors.sbatch` runs the same native pipeline in
`control/`, `moment/` and `affine/` with four CPUs, 16 GB and a 30-minute limit.
The control is the matched-residual candidate; the other configs change only
its descriptor and method name. All patches are regenerated and must match
across arms; count/base/split geometry and templates are reused, not refitted.
The script does no grading, label reading or parameter search.

Use the same five absolute variables as the matched-residual launcher and run
`bash hpc/compare_transverse_descriptors.sbatch --dry-run` first. Submit with
explicit `--export=STMFIT_PROJECT_DIR,STMFIT_CACHE_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN`.
`STMFIT_OUTDIR` must be new. A failed arm is preserved; subsequent arms still
run, and any failure makes the job exit nonzero. Fetch all outputs and verify
the control against the saved 672/29 candidate before external grading.

Completed job **11888011** (source `252edaa`) exits **0:0 in 7m19s** on `vipm2008`;
four CPUs were requested and eight allocated. All 585 regular outputs and the
Slurm log are fetched to `results/transverse_descriptor_comparison_20260920/`,
with matching checksums. All twelve control tables match the saved result.
The first moment gives **672/28** (correct positions / exact chains); affine
residuals give **673/33**, against control **672/29**, all at unchanged coverage.
Retain the explicit affine config as the working candidate, not a champion:
historical 677/36 still leads. All changed affine decisions have zero vote
margin. The report records every gain/loss; no parameter change, hybrid,
counting run or unknown25 application follows. This job is closed; do not
resubmit it as an automatic continuation.

## Independent signed-CC and affine-Fisher comparison (September 20)

`hpc/compare_label_free_signals.sbatch` requests four CPUs, 16 GB and 30 minutes.
It runs `control/` (affine descriptor, 673/33), `signed/` (signed CC costs only)
and `fisher/` (affine Fisher disk only). It does not combine the changes, read
benchmark labels, change counting or search settings. Each arm regenerates
patches; the control must match the saved affine-descriptor result exactly.

Use the five absolute variables `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`,
`STMFIT_INPUT_DIR`, `STMFIT_OUTDIR` and `JULIA_BIN` as in the preceding comparisons.
Run `bash hpc/compare_label_free_signals.sbatch --dry-run` before submitting with
explicit `--export=STMFIT_PROJECT_DIR,STMFIT_CACHE_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN`.
Output must be new. Any failed arm/log is retained, the other arms still run,
and any failure makes the job exit nonzero. Fetch and validate all outputs
before external grading. The raw-patch Fisher attribution diagnostic rejects
the new projected mode rather than misattributing its legacy response identity.
When synchronizing tracked source with `git ls-files`, also copy the existing
root `Manifest.toml` explicitly: it is ignored by Git, and omitting it prevents
dependency loading in a new checkout even when the depot is already populated.
Check its hash; do not resolve new versions to work around a missing copy.

**Completed September 21:** source **3c45c7c**, job **11889205**, exit 0:0 in
**7m05s** on `vipm2008` (four CPUs requested, eight allocated). All 585 outputs
and the Slurm log are fetched and checksum-identical. The twelve control tables
exactly replay 673 correct / 33 exact; all arms keep 900 keys and seven unavailable
assignments. Signed CC gives **668/32**, affine Fisher **672/32**, at the same
849/870 classified coverage. Both are rejected as working replacements; keep
the affine-descriptor 673/33 candidate below historical 677/36. The full result,
all paired gains/losses and replay commands are in
`results/label_free_signals_20260920/report.md`. No follow-up tuning, combined
variant, count refit or unknown25 run is part of this comparison. It is closed;
do not resubmit it automatically.

## Independent Fisher-centering and covariance comparison (September 21)

`hpc/compare_numerical_signals.sbatch` requests four CPUs, 16 GB and 30 minutes.
It runs `control/` (673/33 affine descriptor), `centered/` (training-mean Fisher
origin only), and `shrunk/` (final hard-cluster Ledoit-Wolf covariance only).
Each regenerates patches from the same 146 scans and cached geometry. If all
arms succeed, `test/diagnose_numerical_assignment.jl` verifies fold offsets and
final covariance invariants, and measures GMM sensitivity to two fixed whole-scan
halves. These are numerical/training-perturbation diagnostics, not recognition
grades or chemical confidence. Diagnostic failures preserve their log and fail
the job. No fit runs on the login node.

Set the same five absolute path variables as above; copy the unchanged ignored
`Manifest.toml` explicitly and verify source/input hashes. Run
`bash hpc/compare_numerical_signals.sbatch --dry-run` before the single submission
with explicit
`--export=STMFIT_PROJECT_DIR,STMFIT_CACHE_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN`.
Use a new output directory. Fetch and check all outputs, including exact replay
of the twelve control tables, before the separate local external grade. No
benchmark label, expected N or composition enters this runner. The two changes
are never combined or tuned after grading.

**Completed September 21:** source **0d58b44**, job **11891104**, exit **0:0** in
**9m04s** on `vipl2005` (four CPUs requested, eight allocated). All **590 outputs
and the Slurm log** are fetched and checksum-identical; 76 saved-output checks
pass and the twelve control tables reproduce 673/33 exactly. Centering gives
**675/849 (79.5%), 33/145 exact**, retained as the opt-in working candidate at that stage.
Shrinkage gives **665/849 (78.3%), 34 exact**, not retained as the replacement.
All arms keep the same 900 keys, seven unavailable rows and 849/870 coverage.
Historical 677/36 still leads; there is no champion promotion. The complete
report and all losses are in `results/numerical_signals_20260921/report.md`.
Do not resubmit this completed comparison or add a combined/post-grade arm.

## Independent patch-support and final-score comparison (September 21)

`hpc/compare_patch_support_score.sbatch` requests four CPUs, 16 GB and 30 minutes.
Its three independent arms are `control/` (675/33 centered Fisher), `support/`
(complete central disks, symmetric observed descriptor support), and `volume/`
(final Gaussian-density score on unchanged ridge covariances). It regenerates
patches and assignments for all 146 scans with cached geometry and fixed N.
There is no combined arm, parameter sweep, new normalization or scan weighting.

Use the same five absolute environment variables and explicit `sbatch --export`
list as above. Sync committed source plus the unchanged ignored root Manifest;
run `bash hpc/compare_patch_support_score.sbatch --dry-run` before one submission.
No fit runs on the login node. Fetch and checksum every arm and log, require
byte-identical replay of the twelve saved 675/33 control tables, and retain all
900 keys before external grading. Missing disk pixels remain unavailable; changes
in coverage must be reported alongside all gains and losses. Do not resubmit or
alter a setting in response to a grade.

**Completed September 21:** source **db411a2**, job **11906168**, exit **0:0** in
**7m58s** on `vipc2294` (four CPUs requested, eight allocated). All **585 outputs
and the Slurm log** are fetched and checksum-identical. The twelve control tables
reproduce 675/33 exactly; 76 output checks pass. Support gives **676/852 (79.3%),
34/145 exact**, with coverage rising 849→852/870 and four unavailable keys instead
of seven. It is the latest opt-in candidate, still below historical 677/36.
The Gaussian-score arm gives **666/849 (78.4%), 32 exact** and is rejected.
The full comparison, including every loss, is in
`results/patch_support_score_20260921/report.md`. No automatic resubmission,
combined arm, post-grade tuning, count refit or unknown25 run follows.

## Complete-patch training comparison (September 21)

`hpc/compare_complete_training.sbatch` runs two arms: `control/` replays the
676/34 patch-support candidate; `complete/` fits Fisher/GMM on complete patches
but scores all admissible patches. Four requested CPUs, 16 GB, 30 minutes, one
job. Selected N, geometry, patches, CC, descriptor and k-means stay fixed.
Robust feature scaling, scan weighting and any additional arm are out of scope.

Use the same five absolute environment variables and explicit command-line
`sbatch --export` list. Sync committed source and the unchanged ignored Manifest,
run `bash hpc/compare_complete_training.sbatch --dry-run`, then submit once.
Fetch all outputs/logs, verify exact replay of twelve saved support-control tables,
all 146 files / 900 keys, and both `training_support.tsv` observation counts before
external grading. Report coverage, exact chains and every regression. No fit on
the login node, local real-cohort fit, automatic retry or post-grade tuning.

**Completed September 21:** source **f7a4074**, job **11908779**, exit **0:0**
in **6m42s** on `vipc2262`, after 31m21s queued for resources (four CPUs requested,
eight allocated). All **392 outputs and the Slurm log** are fetched and
checksum-identical. Twelve control tables replay 676/34 exactly; all 900 keys
remain. The complete-training arm gives **674/850 (79.3%), 33/145 exact** and
is rejected. Coverage falls 852→850/870 because one scan lacks any complete
training patch for its normalization. Keep the support candidate at 676/34.
Verification: 1,782 focused assertions, 59 saved-output checks and paired-grade
arithmetic; the initial audit's false per-file completeness assumption is
documented in the report, with no change to the already-frozen abstention rule.
All gains/losses and replay commands: `results/complete_training_20260921/report.md`.
This comparison is complete: no resubmission, normalization fallback, additional
method or post-grade tuning. Historical 677/36 remains the target to exceed.

## Robust per-scan GMM normalization comparison (September 21)

`hpc/compare_robust_normalization.sbatch` runs the 676/34 support `control/` and
one `robust/` candidate using median/Type-7-IQR GMM feature scaling. All admissible
patches remain eligible; upstream signals, k-means and selected N are unchanged.
It requests four CPUs, 16 GB and 30 minutes for one job. No automatic retry or
other method is included. First run `--dry-run` and explicitly pass absolute
`STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`, `STMFIT_OUTDIR` and
`JULIA_BIN` through `sbatch --export=...`; bare environment assignments do not
override Viper's `SBATCH_EXPORT=NONE`. Only compute nodes run scientific stages.
Fetch outputs and logs before external grading.

Source **77b3ce6**, job **11911253**: submitted 13:11:25 CEST, queued 51 seconds,
then **COMPLETED 0:0 in 5m55s**, 13:12:16–13:18:11, on `vipl2005` (eight CPUs
allocated, four Julia threads). All 392 files and the Slurm log are fetched and
checksum-identical. The control exactly reproduces 676/34; median/IQR gives
**667/870 correct, 667/852 (78.3%), 10/145 exact**, with unchanged coverage.
Reject it and retain the support candidate. The comparison is complete; do not
resubmit it or tune a fallback after grading. Full outputs, the corrected
metadata audit and all losses are in
`results/robust_normalization_20260921/report.md`.

## Equal-scan GMM weighting comparison (September 21)

`hpc/compare_scan_weighting.sbatch` runs the 676/34 support `control/` and one
`weighted/` candidate. Only GMM training/naming observation weights change;
mean/std scaling, upstream signals, independent Fisher/k-means heads and N stay
fixed. One job requests four CPUs, 16 GB and 30 minutes; no automatic retry.
Use the five absolute variables `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`,
`STMFIT_INPUT_DIR`, `STMFIT_OUTDIR`, `JULIA_BIN`, run the metadata-only
`--dry-run`, then pass all five explicitly through `sbatch --export=...`.
Sync committed source and the unchanged Manifest, submit once, fetch all outputs
and logs, then check exact control replay and all keys before external grading.
There is no fit on a login node or local real-cohort fit.

Completed job **11915029** (source **5ff68e2**) waited **41m00s** for the shared
CPU quota, then ran **15:36:07–15:42:08 CEST**, **6m01s**, exit **0:0**, on
`vipc2294`. Four CPUs were requested and eight allocated; Julia used four
threads. All **392 output files and the Slurm log** are fetched/checksummed.
The twelve control tables replay exactly; 360 output assertions and paired
grade arithmetic pass after 2,057 pre-run assertions. Both arms retain all
146 scans / 900 keys, unchanged N and upstream/Fisher/k-means tables.

External result: **675/870 correct, 675/852 (79.2%), 36/145 exact**, versus
control **676/852, 34 exact**. It gains two exact chains but loses one correct
position, with all seven final flips ending at zero vote margin. This is a
tradeoff, not a new champion; retain the support candidate as the primary
working reference. Historical 677/36 remains unexceeded. The single comparison
is complete: do not resubmit or tune after grading. Full outputs and all losses
are preserved in `results/scan_weighting_20260921/report.md`.

## Continuous GMM seed-vote comparison (September 21)

`hpc/compare_continuous_vote.sbatch` runs the 676/34 support `control/` and one
`continuous/` candidate. Both regenerate patches with identical geometry and N;
only the aggregation of the same per-seed GMM scores changes. The job requests
four CPUs, 16 GB and 30 minutes, with no automatic retry. Use the same five
absolute variables as the preceding runner: `STMFIT_PROJECT_DIR`,
`STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`, `STMFIT_OUTDIR`, `JULIA_BIN`.
Run `bash hpc/compare_continuous_vote.sbatch --dry-run` first, and pass all five
through an explicit `sbatch --export=...` list. Sync committed source and the
unchanged Manifest, submit once, fetch logs/outputs and verify the control
replay and all identities before grading externally. No real-cohort fitting
runs locally or on a login node.

Scientific source **bef315f**, job **11922389**, submitted at **18:49:32 CEST**
on September 21 with the explicit export list. The source/Manifest checksums and
both 146-scan metadata dry runs match. After **36m23s** under `QOSGrpCpuLimit`,
it runs **19:25:55–19:33:14 CEST**, **7m19s**, and finishes **COMPLETED, 0:0**
on `vipc2034`. Four CPUs were requested, eight allocated; Julia uses four
threads. All **392 output files and the Slurm log** are fetched/checksummed.
All **2,199 pre-run assertions**, **7,240 output checks** and paired-grade
arithmetic pass; the control replays all twelve saved tables byte-for-byte.

External result: **671/870 correct, 671/852 (78.8%), 26/145 exact**, versus
control **676/852, 34 exact** at identical coverage. Reject the continuous
variant and retain the support candidate: one exact chain gained, nine lost.
All 24 final flips are 1→0 from old exact ties. The comparison is complete;
do not resubmit or adjust thresholds/precision after grading. Historical 677/36
remains unexceeded. All outputs, component-score details and paired losses are
preserved in `results/continuous_vote_20260921/report.md`.

## Whole-scan GMM bagging comparison (September 21)

`hpc/compare_scan_bagging.sbatch` runs the 676/34 support `control/` and one
`bagged/` candidate with twenty whole-scan bootstrap replicates. Geometry, N,
patch generation, Fisher, CC and the independent k-means head stay fixed; each
replicate uses the same ten GMM initialization seeds and binary votes.
One job requests four CPUs, 16 GB and **one hour**, with no automatic retry.
Use five absolute variables: `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`,
`STMFIT_INPUT_DIR`, `STMFIT_OUTDIR`, `JULIA_BIN`. Run the metadata-only
`--dry-run` first, then pass all five via an explicit `sbatch --export=...`
list to override Viper's `SBATCH_EXPORT=NONE`. Sync committed source and the
unchanged Manifest; no cohort fits run locally or on a login node.

Fetch outputs and logs, including the candidate's `gmm_scan_bootstrap.tsv`;
check exact control replay, all scan/lobe identities, N, unchanged upstream
outputs and bootstrap multiplicities before external grading. No benchmark
manifest, truth or historical predictions are inference inputs.

Scientific source **da886ce**, job **11925188**: submitted **20:42:12 CEST**,
started **20:42:41** on `vipm2008` after **29 seconds** queued. Four CPUs requested,
eight allocated, four Julia threads, `Requeue=0`. Source and Manifest checksum
checks are empty; both dry runs see 146 scans. The 222 targeted assertions,
189 six-mode pipeline assertions, 108 continuous-vote regression assertions,
20 runner assertions and the documentation build pass. Complementary regressions
also pass: **2,471 assertions in 19 suites** in total, completed early during the
control run. The same job completes **0:0 at 20:49:58**, **7m17s**, with MaxRSS
**1,933,112 K** (about 1.980 GB). No retry or duplicate submission occurs.
All **393 output files and the Slurm log** are fetched locally, excluding
raw-input symlinks, with an empty post-fetch checksum comparison.
**21,927 saved-output checks** pass before grading; all twelve control tables
replay byte-for-byte, all 146 scans / 900 keys and selected N remain, and all
twenty bags match their independently reconstructed draws. All 200 seed fits
are accepted. External grading gives **665/870 correct, 10/145 exact** at
852/870 coverage, versus the replayed **676/34** control. Reject bagging and
retain the support candidate; historical 677/36 remains unexceeded. Seven
component-score checks and paired arithmetic pass; all gains and losses are
saved in `results/scan_bagging_20260921/report.md`. The comparison is complete;
do not resubmit or tune fusion/thresholds after this negative grade.

## Tied-covariance GMM comparison (September 22)

`hpc/compare_tied_covariance.sbatch` runs the 676/34 support `control/` and one
`tied/` candidate. Only covariance sharing throughout GMM learning changes;
both arms regenerate patches at the same cached geometry and selected N.
One job requests four CPUs, 16 GB and **30 minutes**, with `--no-requeue`.
Set five absolute variables: `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`,
`STMFIT_INPUT_DIR`, `STMFIT_OUTDIR`, `JULIA_BIN`. Sync committed code and the
unchanged Manifest, run the metadata-only `--dry-run`, then pass all five in
an explicit `sbatch --export=...` list overriding `SBATCH_EXPORT=NONE`.
No cohort fitting on login nodes or locally, and no duplicate/automatic retry.

Fetch logs and outputs; verify exact control replay, keys, N, unchanged upstream
tables, covariance-mode logs and final vote arithmetic before external grading.
Inference never receives labels, a benchmark manifest or saved reference
predictions.

Completed source **becd0ce**, Viper job **11931106**, **0:0 in 5m53s**, after
15 seconds queued. Slurm allocated eight CPUs for four requested, with
`Requeue=0`, `Restarts=0`; all five variables were explicitly exported. Both
arms finish all 146 scans. All **392 outputs and the Slurm log** are fetched
and checksum-identical; twelve control tables replay the support result exactly.
The **2,716 focused assertions**, **7,338 pre-grade output checks**, five
component-score checks and paired-grade arithmetic pass. All ten tied seed
fits have identical component covariances and freely learned masses.

External grading gives **666/870 correct, 6/145 exact**, versus the replayed
**676/34** control at the same 852/870 coverage. Reject the tied variant;
historical 677/36 remains unexceeded. Complete outputs and losses are saved in
`results/tied_covariance_20260922/report.md`. This comparison is complete: no
retry, other proposed lead, post-grade adjustment, count refit, DFT change or
unknown25 rerun follows. The reused grade is not independent validation.

## Whole-scan Fisher and relative naming comparison (September 22)

`hpc/compare_scan_fisher_naming.sbatch` runs `control/`, `scan_fisher/` and
`relative_naming/` sequentially. The latter two independently change Fisher
grouping and GMM group naming; no combined arm is present. All regenerate
patches at the same cached geometry/N. One job requests four CPUs, 16 GB and
30 minutes, with `--no-requeue`; no automatic retry or login/local cohort fit.

Set `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`, `STMFIT_OUTDIR`
and `JULIA_BIN` to absolute paths. Sync committed code and the unchanged
Manifest, run the metadata-only `--dry-run`, then pass all five explicitly in
`sbatch --export=...` to override Viper's `SBATCH_EXPORT=NONE`. Fetch all results
and logs before cleanup. Verify control replay, all keys/N/upstream tables,
Fisher split/no-overlap logs and naming/unchanged-fit logs before external
grading. Source **43f8993**, job **11931478**, completes **0:0 in 8m30s** after
**3m01s** queued (four CPUs requested, eight allocated, 16 GB requested).
All **588 outputs and the Slurm log** are fetched/checksummed. The 19 suites /
3,035 assertions, 10,987 pre-grade output checks, eight component checks and
paired-grade arithmetic pass. Whole-scan Fisher gives **668/870 correct,
34/145 exact**; relative naming gives **676/870, 34/145**, with unchanged
decisions, at the same **852/870 coverage** as the exactly replayed control.
Keep the **676/34 support candidate**, below historical 677/36. This comparison
is complete; do not automatically resubmit or combine its candidates.
No post-grade tuning, count refit, DFT change or unknown25 rerun follows.
Full gains/losses: `results/scan_fisher_naming_20260922/report.md`.

## Factor-analyzer and Student learning comparison (September 22)

`hpc/compare_factor_student_mixtures.sbatch` runs `control/`,
`factor_analyzer/` and `student_t/` sequentially at the same cached geometry
and selected N. The two learning families are independent candidates with
fixed rank four / df five; no combined arm or parameter search. Patches are
regenerated in each arm. One job requests four CPUs, 16 GB and **one hour**,
with `--no-requeue`; do not compute on login nodes or rerun automatically.

Sync committed source and the unchanged Manifest to
`/u/oldu/code/STMFit_factor_student_mixtures_20260922`. Set the five absolute
variables documented in the runner, use its metadata-only `--dry-run`, and
pass all five explicitly through `sbatch --export=...`. Geometry is
`/ptmp/oldu/stmfit/fusion_comparison_20260920_v1/symmetric`, raw/template input
is `/ptmp/oldu/stmfit/full146_v1_inputs`, Julia is
`/u/oldu/software/julia-1.13.0/bin/julia`, and the new output target is
`/ptmp/oldu/stmfit/factor_student_mixtures_20260922_v1`.

Fetch all outputs/logs locally to `results/factor_student_mixtures_20260922/`.
Before external grading, check exact control replay, all keys/N, unchanged
upstream/Fisher/k-means tables, actual family/seed/convergence logs and final
vote arithmetic. Invalid seeds and unavailable rows must remain visible.
No labels, grade manifest or saved champion predictions enter the job.

Completed source **e2205cb**, job **11935072**, **0:0 in 8m17s**, after
**1m36s** queued; four CPUs requested, eight allocated, 16 GB requested,
MaxRSS 2,238,210 K, no restart/requeue. All **588 outputs plus Slurm log** are
fetched/checksummed. The **3,616 focused assertions**, **11,083 pre-grade
checks**, eight component checks and paired arithmetic pass. Control tables
replay exactly. Fixed factors give **620/870 correct, 23/145 exact**, Student
**631/870, 25/145**, versus **676/34** support, all at **852/870 coverage**.
All Student fits converge; factors reach the fixed 200-update cap. Retain
support, below historical 677/36. This job/comparison is complete: no retry,
cap/rank/df/seed search, combined arm, post-grade fusion tuning or application
rerun follows. Details: `results/factor_student_mixtures_20260922/report.md`.

### All-update shrinkage / coherent Student comparison (2026-09-22)

`hpc/compare_em_shrinkage_student_density.sbatch` runs the unchanged support
control and two independent profiles, `em_shrinkage` and `student_density`.
It reuses saved symmetric geometry and regenerates patches, Fisher and the
two assignment heads for each arm. One job, four requested CPUs, 16 GB,
one hour, no requeue; no parameter sweep, combined arm or count refit.

Set absolute `STMFIT_PROJECT_DIR`, `STMFIT_CACHE_DIR`, `STMFIT_INPUT_DIR`,
`STMFIT_OUTDIR` (new) and `JULIA_BIN` (1.13), then run the script with
`--dry-run` before submission. Supply all five through the actual command's
`sbatch --export=...` list; Viper's `SBATCH_EXPORT=NONE` otherwise drops them.
Use the existing source sync, login-node instantiate, one sbatch, squeue and
fetch workflow. Only compute nodes run the scientific stages. Before grading,
check control byte replay, all 146 scans / 900 keys, unchanged selected N and
upstream tables, seed availability, covariance intensities, convergence and
unchanged vote arithmetic. Nonmonotone ordinary likelihood under shrinkage
is logged without seed selection. Outputs and Slurm logs must be fetched
locally.

Completed source **bb6943b**, job **11936229**, **0:0 in 7m56s**, after
**3s** queued. Four CPUs requested, eight allocated, 16 GB requested,
MaxRSS **2,283,380 K**. All **588 regular outputs and Slurm log** are fetched
and checksummed. **4,277 assertions / 21 suites**, **11,120 pre-grade checks**,
eight component checks and 580 paired scan rows pass. All control tables replay
exactly. Shrinkage gives **637/870 correct, 24/145 exact**; coherent Student
gives **627/870, 25/145**, versus **676/34** support, all at **852/870 coverage**.
Every fit meets the stopping criterion; two shrinkage seeds have reported
ordinary-likelihood decreases, without exclusion. Retain support, below the
historical 677/36 target. No retry, parameter search, combined arm, count refit,
DFT change or unknown25 rerun follows. Full evidence and limitations:
`results/em_shrinkage_student_density_20260922/report.md`.
