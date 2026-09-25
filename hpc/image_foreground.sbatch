#!/usr/bin/env bash
# Two source-only masks, four static file shards, no count/chemical fit or QE.
#SBATCH --job-name=image-foreground
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=foreground_%j.out
#SBATCH --error=foreground_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
: "${STMFIT_FOREGROUND_ROOT:?Project path required}"
: "${STMFIT_FOREGROUND_DATA:?Raw image directory required}"
: "${STMFIT_FOREGROUND_OUT:?New output directory required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test ! -e "$STMFIT_FOREGROUND_OUT"
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_FOREGROUND_ROOT" \
    "$STMFIT_FOREGROUND_ROOT/test/diagnose_image_foreground.jl" \
    --root "$STMFIT_FOREGROUND_ROOT" --data-dir "$STMFIT_FOREGROUND_DATA" \
    --config "$STMFIT_FOREGROUND_ROOT/config/image_foreground.toml" \
    --outdir "$STMFIT_FOREGROUND_OUT"
