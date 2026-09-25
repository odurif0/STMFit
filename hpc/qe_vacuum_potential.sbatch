#!/usr/bin/env bash
# A single PP-only diagnostic, on complete copies of the accepted states.
# No pw.x, SCF, optimization, density replacement, retry or grade.
#SBATCH --job-name=qe-vacuum-potential
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_potential_%j.out
#SBATCH --error=qe_potential_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_NTASKS:?} == 8 ]] || exit 1
module purge
module load intel/2024.0 impi/2021.11 qe/7.4.1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
STMFIT_POTENTIAL_ROOT=$(cd ../.. && pwd)
STMFIT_POTENTIAL_RUN=$PWD
STMFIT_GLCN_SOURCE=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save
STMFIT_GLCNAC_SOURCE=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save
test ! -e work
test ! -L work
test -s "$STMFIT_GLCN_SOURCE/charge-density.dat"
test -s "$STMFIT_GLCNAC_SOURCE/charge-density.dat"
sha256sum "$STMFIT_GLCN_SOURCE/wfc1.dat" "$STMFIT_GLCNAC_SOURCE/wfc1.dat" \
    "$STMFIT_GLCN_SOURCE/charge-density.dat" "$STMFIT_GLCNAC_SOURCE/charge-density.dat" > original_state_before.sha256
mkdir -p work/glcn work/glcnac
cp -a "$STMFIT_GLCN_SOURCE" work/glcn/
cp -a "$STMFIT_GLCNAC_SOURCE" work/glcnac/
sha256sum -c --quiet accepted_states.sha256
find work -type f -print0 | sort -z | xargs -0 sha256sum > checkpoint_before.sha256
sha256sum "$(command -v pp.x)" > pp_binary.sha256
for molecule in glcn glcnac; do
    cmp "work/$molecule/${molecule}_central.save/data-file-schema.xml" "$molecule/data-file-schema.xml"
    cd "$molecule"
    for observable in total electrostatic repeat_total; do
        srun -n 8 --cpu-bind=cores pp.x -in "$observable.in" > "$observable.out"
        grep -q 'JOB DONE.' "$observable.out"
    done
    cmp total.dat repeat_total.dat
    cmp total.cube repeat_total.cube
    cd "$STMFIT_POTENTIAL_RUN"
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_POTENTIAL_ROOT" \
        "$STMFIT_POTENTIAL_ROOT/test/qe_vacuum_potential.jl" analyze . "$molecule"
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_POTENTIAL_ROOT" \
    "$STMFIT_POTENTIAL_ROOT/test/test_qe_vacuum_potential.jl" . > independent_verification.out 2>&1
sha256sum -c --quiet checkpoint_before.sha256
sha256sum -c --quiet original_state_before.sha256
echo 'VACUUM_POTENTIAL_COMPLETE; no matching plane, wavefunction or reference adopted'
