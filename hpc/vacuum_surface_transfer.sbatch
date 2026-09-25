#!/usr/bin/env bash
# Four static file shards inside one Julia allocation; no new QE or production fit.
#SBATCH --job-name=vacuum-surface-transfer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=transfer_%j.out
#SBATCH --error=transfer_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 path required}"
: "${STMFIT_TRANSFER_ROOT:?Project path required}"
: "${STMFIT_TRANSFER_DATA:?Raw image path required}"
: "${STMFIT_TRANSFER_OUT:?New output path required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test ! -e "$STMFIT_TRANSFER_OUT"
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_TRANSFER_ROOT" \
    "$STMFIT_TRANSFER_ROOT/test/diagnose_vacuum_surface_transfer.jl" \
    --root "$STMFIT_TRANSFER_ROOT" --data-dir "$STMFIT_TRANSFER_DATA" \
    --config "$STMFIT_TRANSFER_ROOT/config/vacuum_surface_transfer.toml" \
    --outdir "$STMFIT_TRANSFER_OUT"
