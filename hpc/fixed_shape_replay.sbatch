#!/usr/bin/env bash
# Fixed saved physical candidates, four static file shards, no optimization/QE.
#SBATCH --job-name=fixed-shape-replay
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=fixed_shape_%j.out
#SBATCH --error=fixed_shape_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
: "${STMFIT_REPLAY_ROOT:?Project path required}"
: "${STMFIT_REPLAY_OUT:?New output path required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test ! -e "$STMFIT_REPLAY_OUT"
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_REPLAY_ROOT" \
    "$STMFIT_REPLAY_ROOT/test/diagnose_fixed_shape_replay.jl" \
    --root "$STMFIT_REPLAY_ROOT" --config "$STMFIT_REPLAY_ROOT/config/fixed_shape_replay.toml" \
    --outdir "$STMFIT_REPLAY_OUT"
