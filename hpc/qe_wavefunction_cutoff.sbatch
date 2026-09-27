#!/usr/bin/env bash
# One geometry per job; exactly two sequential jobs, no retry or requeue.
#SBATCH --job-name=qe-wfc60
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_basis_%j.out
#SBATCH --error=qe_basis_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_NTASKS:?} == 8 ]] || exit 1
STMFIT_BASIS_M=$(basename "$PWD")
case "$STMFIT_BASIS_M" in glcn|glcnac) ;; *) exit 1 ;; esac
STMFIT_BASIS_ROOT=$(cd ../../.. && pwd)
STMFIT_BASIS_RUN=$(cd .. && pwd)
STMFIT_BASIS_OLD=$STMFIT_BASIS_ROOT/qe/density_cutoff_20260926/$STMFIT_BASIS_M/work/${STMFIT_BASIS_M}_central.save
sha256sum -c --quiet input_hashes.sha256
test ! -e work
test ! -e comparison
cmp 50/data-file-schema.xml "$STMFIT_BASIS_OLD/data-file-schema.xml"
sha256sum "$STMFIT_BASIS_OLD/data-file-schema.xml" "$STMFIT_BASIS_OLD/wfc1.dat" \
    "$STMFIT_BASIS_OLD/charge-density.dat" > original_before.sha256
trap 'sha256sum -c --quiet original_before.sha256' EXIT
module purge
module load intel/2024.0 impi/2021.11 qe/7.4.1 mkl/2024.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sha256sum "$(command -v pw.x)" > executable.sha256
# Verify the saved 50-Ry reference before spending the SCF allocation.
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_BASIS_ROOT" \
    "$STMFIT_BASIS_ROOT/test/qe_wavefunction_cutoff.jl" analyze "$STMFIT_BASIS_RUN" "$STMFIT_BASIS_M" 50 "$STMFIT_BASIS_OLD/wfc1.dat"
# Fresh atomic potential and atomic+random wavefunctions; no cross-basis file copy.
srun -n 8 --cpu-bind=cores pw.x -in pw_scf.in > pw_scf.out
"$JULIA_BIN" --startup-file=no --project="$STMFIT_BASIS_ROOT" \
    "$STMFIT_BASIS_ROOT/test/qe_wavefunction_cutoff.jl" check-scf "$STMFIT_BASIS_RUN" "$STMFIT_BASIS_M"
find work -type f -print0 | sort -z | xargs -0 sha256sum > converged_before_analysis.sha256
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_BASIS_ROOT" \
    "$STMFIT_BASIS_ROOT/test/qe_wavefunction_cutoff.jl" analyze "$STMFIT_BASIS_RUN" "$STMFIT_BASIS_M" 60 "work/${STMFIT_BASIS_M}_central.save/wfc1.dat"
"$JULIA_BIN" --startup-file=no --project="$STMFIT_BASIS_ROOT" \
    "$STMFIT_BASIS_ROOT/test/qe_wavefunction_cutoff.jl" summarize "$STMFIT_BASIS_RUN" "$STMFIT_BASIS_M"
"$JULIA_BIN" --startup-file=no --project="$STMFIT_BASIS_ROOT" \
    "$STMFIT_BASIS_ROOT/test/test_qe_wavefunction_cutoff.jl" "$STMFIT_BASIS_RUN" "$STMFIT_BASIS_M" > independent_verification.out 2>&1
sha256sum -c --quiet input_hashes.sha256 converged_before_analysis.sha256 original_before.sha256
echo 'WAVEFUNCTION_CUTOFF_SENSITIVITY_COMPLETE; no mold, calibration or benchmark adopted'
