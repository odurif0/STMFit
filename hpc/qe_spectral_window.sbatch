#!/usr/bin/env bash
# Copied by qe_spectral_window.jl into an isolated ignored qe/ run directory.
#SBATCH --job-name=qe-spectral-window
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_spectral_%j.out
#SBATCH --error=qe_spectral_%j.err
set -euo pipefail
module purge
module load "${QE_COMPILER_MODULE:-intel/2024.0}" "${QE_MPI_MODULE:-impi/2021.11}" "${QE_MODULE:-qe/7.4.1}"
export OMP_NUM_THREADS=1
QE_NTASKS=${SLURM_NTASKS:?Slurm allocation required}
[[ "$QE_NTASKS" == 8 ]] || { echo 'Eight MPI ranks required for the saved distributed checkpoint'; exit 1; }
: "${JULIA_BIN:?Explicit Julia 1.13 binary required}"
STMFIT_SPECTRAL_ROOT=$(cd ../.. && pwd)
STMFIT_SPECTRAL_RUN=$PWD
STMFIT_QE_SOURCE=/u/oldu/code/STMFit/qe
test ! -e work
test ! -L work
test -s "$STMFIT_QE_SOURCE/glcn_restart5/qe_tmp/glcn_central.save/wfc1.dat"
test -s "$STMFIT_QE_SOURCE/glcnac/qe_tmp/glcnac_central.wfc8"
mkdir work
# Full copies: no hard links and no writes to the only surviving old checkpoint.
cp -a "$STMFIT_QE_SOURCE/glcn_restart5/qe_tmp" work/glcn
cp -a "$STMFIT_QE_SOURCE/glcnac/qe_tmp" work/glcnac
find work -type f -print0 | sort -z | xargs -0 sha256sum > checkpoint_before.sha256
srun -n "$QE_NTASKS" --cpu-bind=cores pw.x -in pw_scf.in > restored_glcnac_scf.out
grep -q 'convergence has been achieved' restored_glcnac_scf.out
grep -q 'JOB DONE.' restored_glcnac_scf.out
test -s work/glcnac/glcnac_central.save/wfc1.dat
# Retain the old launcher's literal control invocation for its preflight.
srun -n "$QE_NTASKS" --cpu-bind=cores pp.x -in pp_ldos.in > restored_glcnac_control.out
cmp glcnac_central_ldos_accept_plain.cube glcnac_reference.cube
for molecule in glcn glcnac; do
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_SPECTRAL_ROOT" \
        "$STMFIT_SPECTRAL_ROOT/test/qe_spectral_window.jl" prepare-spectrum \
        --xml "work/$molecule/${molecule}_central.save/data-file-schema.xml" \
        --settings settings.toml --prefix "${molecule}_central" \
        --checkpoint "work/$molecule" --outdir "$molecule"
    cd "$molecule"
    srun -n "$QE_NTASKS" --cpu-bind=cores pp.x -in control.in > control.out
    cmp control.cube "../${molecule}_reference.cube"
    srun -n "$QE_NTASKS" --cpu-bind=cores pp.x -in ildos.in > ildos.out
    grep -q 'JOB DONE.' control.out
    grep -q 'JOB DONE.' ildos.out
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_SPECTRAL_ROOT" \
        "$STMFIT_SPECTRAL_ROOT/test/qe_spectral_window.jl" compare \
        --control control.cube --candidate ildos.cube --reference "../${molecule}_reference.cube" \
        --frame "../${molecule}_frame.tsv" --metadata metadata.toml --settings ../settings.toml --outdir comparison
    cd "$STMFIT_SPECTRAL_RUN"
done
find work -type f -print0 | sort -z | xargs -0 sha256sum > checkpoint_after.sha256
sha256sum settings.toml pw_scf.in pp_ldos.in run_scf_pp.sbatch > run_inputs_after.sha256
echo 'SPECTRAL_DIAGNOSTIC_COMPLETE; no template, prediction or benchmark was produced'
