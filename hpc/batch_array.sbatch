#!/bin/bash -l
#SBATCH --job-name=stmfit
#SBATCH --output=results/hpc_logs/stmfit_%A_%a.out
#SBATCH --error=results/hpc_logs/stmfit_%A_%a.err
#
# ── Job array: one Slurm task per chunk of test/batch_full.jl ─────────────────
# Each task runs a disjoint round-robin slice of the file list via --chunk.
# batch_full.jl already skips finished files, so re-submitting is safe & resumes.
#
#   IMPORTANT: the array range / cpus / time below are sensible DEFAULTS for
#   manual use. The launcher (hpc/launch_remote.sh) overrides them on the
#   `sbatch` command line (Slurm: CLI options win over in-script #SBATCH).
#   Slurm parses #SBATCH lines literally — NO shell expansion here.
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=16000MB
#
# ── Partition (optional) ─────────────────────────────────────────────────────
# On Viper the submit filter picks the partition automatically from the
# resources above — do NOT set --partition there. On Raven, leave it unset too:
# the submit filter auto-routes jobs smaller than half a node to the shared
# pool. NOTE: Raven's shared-job filter REQUIRES a total `--mem` (with the MB
# suffix) — `--mem-per-cpu` is NOT accepted and triggers
#   "memory limit must be provided for shared jobs".
##SBATCH --partition=general
#
# ── Notifications ────────────────────────────────────────────────────────────
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Environment (these DO expand — this is bash, run by Slurm on the compute node)
# ─────────────────────────────────────────────────────────────────────────────
PROJECT_DIR="${STMFIT_PROJECT_DIR:-$(pwd)}"
cd "$PROJECT_DIR"
mkdir -p results/hpc_logs

# N_CHUNKS is the sharding DENOMINATOR used in --chunk i/N_CHUNKS. It must be
# stable across partial re-submissions: re-running only task 3 of an 8-way
# shard still needs --chunk 3/8 (NOT 3/1). Export it via the launcher, or with
#   sbatch --export=ALL,N_CHUNKS=8 hpc/batch_array.sbatch
# When unset (pure manual use with the literal array above), infer from the
# array's max index so the default 4-way case just works.
if [[ -z "${N_CHUNKS:-}" ]]; then
    N_CHUNKS="${SLURM_ARRAY_TASK_MAX:-4}"
fi
export N_CHUNKS

# Load Julia. Pin a version after checking `find-module julia` on the login
# node (e.g. julia/1.11.4); bare `module load julia` picks the default.
module purge
module load julia/${JULIA_MODULE_VERSION:-}

# Threading: one Julia thread per allocated CPU. batch_full.jl internally caps
# at min(4, Threads.nthreads(), ntot), so extra threads won't help. Pin BLAS to
# a single thread — each file is a small least-squares problem and the Julia
# threads already parallelise across files.
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

echo "============================================================"
echo " STMFit batch chunk on $(hostname)"
echo " job ${SLURM_JOB_ID} (array task ${SLURM_ARRAY_TASK_ID}/${N_CHUNKS})"
echo " project:  $PROJECT_DIR"
echo " data dir: ${STMFIT_DATA_DIR:-<unset — relies on --data-dir / default>}"
echo " outdir:   ${STMFIT_OUTDIR:-results/best_plots}"
echo " config:   ${STMFIT_CONFIG:-config/chitosan.toml}"
echo " cpus:     ${SLURM_CPUS_PER_TASK:-4}   julia threads: $JULIA_NUM_THREADS"
echo " started:  $(date -Is)"
echo "============================================================"

command -v julia >/dev/null || { echo "ERROR: julia module not loaded"; exit 1; }
julia --version

# ─────────────────────────────────────────────────────────────────────────────
# Assemble the batch_full.jl invocation.
# Chunk index = SLURM_ARRAY_TASK_ID, total = N_CHUNKS.
# ─────────────────────────────────────────────────────────────────────────────
N_FILES_ARG="${N_FILES:-}"
CONFIG="${STMFIT_CONFIG:-config/chitosan.toml}"
OUTDIR="${STMFIT_OUTDIR:-results/best_plots}"
TSP="${SLURM_ARRAY_TASK_ID}/${N_CHUNKS}"

BATCH_ARGS=()
[[ -n "$N_FILES_ARG" ]] && BATCH_ARGS+=("$N_FILES_ARG")
[[ -n "${STMFIT_DATA_DIR:-}" ]] && BATCH_ARGS+=("--data-dir" "$STMFIT_DATA_DIR")
BATCH_ARGS+=("--outdir" "$OUTDIR")
BATCH_ARGS+=("--config" "$CONFIG")
[[ -n "${STMFIT_SELECTION_POLICY:-}" ]] && BATCH_ARGS+=("--selection-policy" "$STMFIT_SELECTION_POLICY")
[[ -n "${STMFIT_REFINED_ADVISORY:-}" ]] && BATCH_ARGS+=("--refined-advisory" "$STMFIT_REFINED_ADVISORY")
[[ -n "${STMFIT_PLOT_MANIFEST:-}" ]] && BATCH_ARGS+=("--plot-manifest" "$STMFIT_PLOT_MANIFEST")
[[ -n "${STMFIT_TSV:-}" ]] && BATCH_ARGS+=("--tsv" "$STMFIT_TSV")
[[ "${STMFIT_SKIP_1D:-0}" == "1" ]] && BATCH_ARGS+=("--skip-1d")

set +e
julia --project=. -t "$JULIA_NUM_THREADS" test/batch_full.jl \
    "${BATCH_ARGS[@]}" \
    ${STMFIT_BATCH_ARGS:-} \
    --chunk "$TSP"
RC=$?
set -e

echo "------------------------------------------------------------"
echo " chunk $TSP finished with exit code $RC at $(date -Is)"
echo "------------------------------------------------------------"
exit $RC
