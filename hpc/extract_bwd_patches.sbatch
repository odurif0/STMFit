#!/usr/bin/env bash
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
echo "=== extract backward patches (half ${STMFIT_HALF_NM:-0.32}, step ${STMFIT_STEP_NM:-0.04}) ==="
julia --project=. test/extract_lobe_patches_bwd.jl \
    --features "${STMFIT_FEATURES}" \
    --data-dir "${STMFIT_DATA_DIR}" \
    --half-nm "${STMFIT_HALF_NM:-0.32}" --step-nm "${STMFIT_STEP_NM:-0.04}" \
    --out "${STMFIT_PATCHES_OUT}"
wc -l "${STMFIT_PATCHES_OUT}"
