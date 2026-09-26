#!/usr/bin/env bash
# Analysis of saved exports only. No QE executable, SCF, rebuild or requeue.
#SBATCH --job-name=qe-xc-recheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_xc_recheck_%j.out
#SBATCH --error=qe_xc_recheck_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
STMFIT_XCR_RUN=$PWD
STMFIT_XCR_ROOT=$(cd ../.. && pwd)
STMFIT_XCR_SOURCE=$STMFIT_XCR_ROOT/qe/xc_build_control_20260926
test ! -e glcn/controls
test ! -e glcnac/controls
test ! -e glcn/analysis
test ! -e glcnac/analysis
# Hash only the immutable source archive, excluding its copied checkpoints.
(cd "$STMFIT_XCR_SOURCE" && find . -path ./work -prune -o -type f -print0 | sort -z | xargs -0 sha256sum) > source_archive_before.sha256
for molecule in glcn glcnac; do
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_XCR_ROOT" \
        "$STMFIT_XCR_ROOT/test/qe_xc_build_control.jl" report "$STMFIT_XCR_ROOT" . "$molecule"
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_XCR_ROOT" \
    "$STMFIT_XCR_ROOT/test/qe_xc_build_control.jl" check .
for molecule in glcn glcnac; do
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_XCR_ROOT" \
        "$STMFIT_XCR_ROOT/test/qe_xc_build_control.jl" analyze "$STMFIT_XCR_ROOT" . "$molecule"
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_XCR_ROOT" \
    "$STMFIT_XCR_ROOT/test/test_qe_xc_build_control.jl" "$STMFIT_XCR_ROOT" "$STMFIT_XCR_RUN" > independent_verification.out 2>&1
(cd "$STMFIT_XCR_SOURCE" && sha256sum -c --quiet "$STMFIT_XCR_RUN/source_archive_before.sha256")
echo 'XC_SAVED_EXPORT_RECHECK_COMPLETE; no QE or new physical potential'
