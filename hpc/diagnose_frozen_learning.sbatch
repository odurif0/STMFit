#!/usr/bin/env bash
#SBATCH --job-name=stmfit-frozen-learning
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16000MB
#SBATCH --no-requeue
set -euo pipefail
# Two independent one-thread repetitions; same first observed-fit inputs.
# Viper: pass every required variable via an explicit sbatch --export list.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set committed source checkout}"
: "${STMFIT_INPUT_DIR:?Set first complete observed-fit output directory}"
: "${STMFIT_OUTDIR:?Set new diagnostic output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 GKSwstype=100
command=("$JULIA_BIN" --startup-file=no -t 1 --project=. test/diagnose_frozen_learning.jl
    --input "$STMFIT_INPUT_DIR" --settings config/frozen_learning.toml)
if [[ $# -eq 1 ]]; then
    exec "${command[@]}" --outdir "$STMFIT_OUTDIR/repeat1" --dry-run
fi
: "${SLURM_JOB_ID:?Run scientific diagnostic on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 2 ]] || { echo "Two CPUs required" >&2; exit 2; }
mkdir "$STMFIT_OUTDIR"
pids=()
for repetition in 1 2; do
    "${command[@]}" --outdir "$STMFIT_OUTDIR/repeat$repetition" > "$STMFIT_OUTDIR/repeat$repetition.log" 2>&1 &
    pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
    if wait "$pid"; then :; else status=1; fi
done
exit "$status"
