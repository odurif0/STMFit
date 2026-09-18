# Running STMFit on the MPCDF HPC cluster

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
