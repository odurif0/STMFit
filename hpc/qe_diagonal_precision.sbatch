#!/usr/bin/env bash
# One geometry per allocation; two independent NSCF solves, no retries/new SCF.
#SBATCH --job-name=qe-diag
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_diag_%j.out
#SBATCH --error=qe_diag_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_NTASKS:?} == 8 ]] || exit 1
STMFIT_DIAG_M=$(basename "$PWD")
case "$STMFIT_DIAG_M" in glcn|glcnac) ;; *) exit 1 ;; esac
STMFIT_DIAG_ROOT=$(cd ../../.. && pwd)
STMFIT_DIAG_RUN=$(cd .. && pwd)
STMFIT_DIAG_CASE=$PWD
STMFIT_DIAG_OLD=$STMFIT_DIAG_ROOT/qe/wavefunction_cutoff_20260927/$STMFIT_DIAG_M/work/${STMFIT_DIAG_M}_central.save
sha256sum -c --quiet input_hashes.sha256
test ! -e lo/work
test ! -e hi/work
test ! -e comparison
cmp reference/data-file-schema.xml "$STMFIT_DIAG_OLD/data-file-schema.xml"
find "$STMFIT_DIAG_OLD" -type f -print0 | sort -z | xargs -0 sha256sum > original_before.sha256
trap 'cd "$STMFIT_DIAG_CASE"; sha256sum -c --quiet original_before.sha256' EXIT
module purge
module load intel/2024.0 impi/2021.11 qe/7.4.1 mkl/2024.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sha256sum "$(command -v pw.x)" > executable.sha256
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_DIAG_ROOT" \
    "$STMFIT_DIAG_ROOT/test/qe_diagonal_precision.jl" analyze "$STMFIT_DIAG_RUN" "$STMFIT_DIAG_M" reference "$STMFIT_DIAG_OLD/wfc1.dat"
for STMFIT_DIAG_TRIAL in lo hi; do
    cd "$STMFIT_DIAG_CASE/$STMFIT_DIAG_TRIAL"
    mkdir work
    cp -a "$STMFIT_DIAG_OLD" work/
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_DIAG_ROOT" \
        "$STMFIT_DIAG_ROOT/test/qe_diagonal_precision.jl" check-start "$STMFIT_DIAG_RUN" "$STMFIT_DIAG_M" "$STMFIT_DIAG_TRIAL"
    srun -n 8 --cpu-bind=cores pw.x -in pw_nscf.in > pw_nscf.out
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_DIAG_ROOT" \
        "$STMFIT_DIAG_ROOT/test/qe_diagonal_precision.jl" check-nscf "$STMFIT_DIAG_RUN" "$STMFIT_DIAG_M" "$STMFIT_DIAG_TRIAL"
    find work -type f -print0 | sort -z | xargs -0 sha256sum > converged_before_analysis.sha256
    "$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_DIAG_ROOT" \
        "$STMFIT_DIAG_ROOT/test/qe_diagonal_precision.jl" analyze "$STMFIT_DIAG_RUN" "$STMFIT_DIAG_M" "$STMFIT_DIAG_TRIAL" "work/${STMFIT_DIAG_M}_central.save/wfc1.dat"
    sha256sum -c --quiet converged_before_analysis.sha256
done
cd "$STMFIT_DIAG_CASE"
"$JULIA_BIN" --startup-file=no --project="$STMFIT_DIAG_ROOT" \
    "$STMFIT_DIAG_ROOT/test/qe_diagonal_precision.jl" summarize "$STMFIT_DIAG_RUN" "$STMFIT_DIAG_M"
"$JULIA_BIN" --startup-file=no --project="$STMFIT_DIAG_ROOT" \
    "$STMFIT_DIAG_ROOT/test/test_qe_diagonal_precision.jl" "$STMFIT_DIAG_RUN" "$STMFIT_DIAG_M" > independent_verification.out 2>&1
sha256sum -c --quiet input_hashes.sha256 original_before.sha256
echo 'DIAGONAL_PRECISION_COMPLETE; fixed density; no new SCF, mold, calibration or benchmark adopted'
