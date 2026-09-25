#!/usr/bin/env bash
# One image-envelope diagnostic, four static file shards, no raw/QE recomputation.
#SBATCH --job-name=shared-envelope
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=shared_envelope_%j.out
#SBATCH --error=shared_envelope_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
: "${STMFIT_ENVELOPE_ROOT:?Project path required}"
: "${STMFIT_ENVELOPE_OUT:?New output directory required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test ! -e "$STMFIT_ENVELOPE_OUT"
"$JULIA_BIN" --version
sha256sum "$STMFIT_ENVELOPE_ROOT/config/shared_image_envelope.toml" \
    "$STMFIT_ENVELOPE_ROOT/test/lib/shared_image_envelope.jl" \
    "$STMFIT_ENVELOPE_ROOT/test/diagnose_shared_image_envelope.jl" \
    "$STMFIT_ENVELOPE_ROOT/Project.toml" "$STMFIT_ENVELOPE_ROOT/Manifest.toml"
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_ENVELOPE_ROOT" \
    "$STMFIT_ENVELOPE_ROOT/test/diagnose_shared_image_envelope.jl" \
    --root "$STMFIT_ENVELOPE_ROOT" \
    --config "$STMFIT_ENVELOPE_ROOT/config/shared_image_envelope.toml" \
    --outdir "$STMFIT_ENVELOPE_OUT"
