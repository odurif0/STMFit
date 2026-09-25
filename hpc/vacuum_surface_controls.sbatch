#!/usr/bin/env bash
# Four static file shards; only frozen observations, no raw count or QE rerun.
#SBATCH --job-name=vacuum-surface-controls
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=controls_%j.out
#SBATCH --error=controls_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 path required}"
: "${STMFIT_CONTROLS_ROOT:?Project path required}"
: "${STMFIT_CONTROLS_OUT:?New output path required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test ! -e "$STMFIT_CONTROLS_OUT"
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_CONTROLS_ROOT" \
    "$STMFIT_CONTROLS_ROOT/test/diagnose_vacuum_surface_controls.jl" \
    --root "$STMFIT_CONTROLS_ROOT" --config "$STMFIT_CONTROLS_ROOT/config/vacuum_surface_controls.toml" \
    --outdir "$STMFIT_CONTROLS_OUT"
