#!/usr/bin/env bash
# Slurm job: forward residual 9x9 patch extraction for all 145 benchmark files.
# Produces the patch TSV consumed by score_connected_mold_templates.jl against
# the converged DFT-STM v1 templates (same 0.32 nm / 0.08 nm grid).
# Expected env: STMFIT_PROJECT_DIR, STMFIT_DATA_DIR, STMFIT_FEATURES,
# STMFIT_PATCHES_OUT, JULIA_MODULE_VERSION.

#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000MB
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

module purge
module load "julia/${JULIA_MODULE_VERSION}"
cd "${STMFIT_PROJECT_DIR}"
export JULIA_NUM_THREADS=2

echo "=== extract forward patches (half ${STMFIT_HALF_NM:-0.32}, half_u ${STMFIT_HALF_U_NM:-0.32}, step ${STMFIT_STEP_NM:-0.08}) ==="
julia --project=. test/extract_lobe_patches.jl \
    --features "${STMFIT_FEATURES}" \
    --data-dir "${STMFIT_DATA_DIR}" \
    --half-nm "${STMFIT_HALF_NM:-0.32}" --half-u-nm "${STMFIT_HALF_U_NM:-0.32}" --step-nm "${STMFIT_STEP_NM:-0.08}" \
    --out "${STMFIT_PATCHES_OUT}"

echo "=== rows extracted ==="
wc -l "${STMFIT_PATCHES_OUT}"
