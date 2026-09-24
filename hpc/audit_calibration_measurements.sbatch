#!/usr/bin/env bash
#SBATCH --job-name=stmfit-calib-measure
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
# Four ordinary Julia shards. No fit, benchmark, template or calibration output.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set committed source checkout}"
: "${STMFIT_DATA_DIR:?Set raw SXM directory}"
: "${STMFIT_OUTDIR:?Set new output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_DATA_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
if [[ $# -eq 1 ]]; then
    for shard in 1 2 3 4; do
        "$JULIA_BIN" --startup-file=no -t 1 --project=. test/measure_calibration.jl \
            --data-dir "$STMFIT_DATA_DIR" --config config/calibration_measurements.toml \
            --outdir "$STMFIT_OUTDIR/chunk$shard" --chunk "$shard/4" --dry-run
    done
    exit 0
fi
: "${SLURM_JOB_ID:?Run on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
mkdir -p "$STMFIT_OUTDIR"
pids=()
for shard in 1 2 3 4; do
    "$JULIA_BIN" --startup-file=no -t 1 --project=. test/measure_calibration.jl \
        --data-dir "$STMFIT_DATA_DIR" --config config/calibration_measurements.toml \
        --outdir "$STMFIT_OUTDIR/chunk$shard" --chunk "$shard/4" \
        > "$STMFIT_OUTDIR/chunk$shard.log" 2>&1 &
    pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do wait "$pid" || status=1; done
exit "$status"
