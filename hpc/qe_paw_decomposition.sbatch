#!/usr/bin/env bash
# One PP-only attribution experiment; no SCF, fitting or density correction.
#SBATCH --job-name=qe-paw-components
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_paw_%j.out
#SBATCH --error=qe_paw_%j.err
set -euo pipefail
module purge
module load "${QE_COMPILER_MODULE:-intel/2024.0}" "${QE_MPI_MODULE:-impi/2021.11}" "${QE_MODULE:-qe/7.4.1}" mkl/2024.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
QE_NTASKS=${SLURM_NTASKS:?Slurm allocation required}
[[ "$QE_NTASKS" == 8 ]] || exit 1
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
STMFIT_PAW_ROOT=$(cd ../.. && pwd)
STMFIT_PAW_RUN=$PWD
STMFIT_QE_BUILD=$STMFIT_PAW_ROOT/qe/gamma_source_20260925
STMFIT_MKL_LIBDIR=${MKLROOT:?MKL build module required}/lib/intel64
STMFIT_MKL_LINK_FLAGS=("-L$STMFIT_MKL_LIBDIR" "-Wl,-rpath,$STMFIT_MKL_LIBDIR" -mkl=sequential)
STMFIT_GLCN_SOURCE=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save
STMFIT_GLCNAC_SOURCE=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save
test ! -e work
test ! -e binaries
test ! -e "$STMFIT_QE_BUILD/make.inc"
test "$(git -C "$STMFIT_QE_BUILD" rev-parse HEAD)" = 500de340b820e1cb8c05f2d8bb8fced102f377c1
git -C "$STMFIT_QE_BUILD" diff --exit-code
git -C "$STMFIT_QE_BUILD" diff --cached --exit-code
test "$(git -C "$STMFIT_QE_BUILD/external/fox" rev-parse HEAD)" = 3453648e6837658b747b895bb7bef4b1ed2eac40
test "$(git -C "$STMFIT_QE_BUILD/external/mbd" rev-parse HEAD)" = 89a3cc199c0a200c9f0f688c3229ef6b9a8d63bd
test "$(git -C "$STMFIT_QE_BUILD/external/devxlib" rev-parse HEAD)" = a6b89ef77b1ceda48e967921f1f5488d2df9226d
for library in mkl_intel_lp64 mkl_sequential mkl_core; do
    test -r "$STMFIT_MKL_LIBDIR/lib$library.so"
done
mkdir -p build_logs
mpiifort "$STMFIT_PAW_ROOT/test/fixtures/qe_mkl_link_probe.f90" \
    "${STMFIT_MKL_LINK_FLAGS[@]}" -o build_logs/mkl_link_probe > build_logs/mkl_link.out 2>&1
build_logs/mkl_link_probe > build_logs/mkl_smoke.out 2>&1
sha256sum "$STMFIT_GLCN_SOURCE/wfc1.dat" "$STMFIT_GLCNAC_SOURCE/wfc1.dat" > original_wfc_before.sha256
mkdir -p work/glcn work/glcnac binaries
cp -a "$STMFIT_GLCN_SOURCE" work/glcn/
cp -a "$STMFIT_GLCNAC_SOURCE" work/glcnac/
sha256sum -c --quiet collected_wfc.sha256
find work -type f -print0 | sort -z | xargs -0 sha256sum > checkpoint_before.sha256
for molecule in glcn glcnac; do
    cmp "work/$molecule/${molecule}_central.save/data-file-schema.xml" "$molecule/data-file-schema.xml"
    cd "$molecule"
    srun -n "$QE_NTASKS" --cpu-bind=cores pp.x -in pp_ldos.in > site_control.out
    grep -q 'JOB DONE.' site_control.out
    cmp site_control.cube reference.cube
    cd "$STMFIT_PAW_RUN"
done
cd "$STMFIT_QE_BUILD"
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=protocol.allow GIT_CONFIG_VALUE_0=never
F90=ifort MPIF90=mpiifort CC=icx BLAS_LIBS="${STMFIT_MKL_LINK_FLAGS[*]}" \
    ./configure --enable-parallel --disable-openmp --with-scalapack=no --with-fox=yes \
    > "$STMFIT_PAW_RUN/build_logs/configure.out" 2>&1
cp install/config.log "$STMFIT_PAW_RUN/build_logs/config.log"
grep -q '^configure: exit 0$' install/config.log
cp make.inc "$STMFIT_PAW_RUN/build_logs/make.inc"
grep -Eq '^LAPACK[[:space:]]*=[[:space:]]*$' make.inc
git apply --check "$STMFIT_PAW_ROOT/hpc/patches/qe-7.4.1-gamma-density.patch"
git apply "$STMFIT_PAW_ROOT/hpc/patches/qe-7.4.1-gamma-density.patch"
make -j 8 pwlibs > "$STMFIT_PAW_RUN/build_logs/libraries.out" 2>&1
make -j 8 -C PP/src pp.x > "$STMFIT_PAW_RUN/build_logs/full_pp.out" 2>&1
cp -p PP/src/pp.x "$STMFIT_PAW_RUN/binaries/pp_full.x"
cp PP/src/local_dos.f90 "$STMFIT_PAW_RUN/build_logs/local_dos.full.f90"
git diff -- PP/src/local_dos.f90 > "$STMFIT_PAW_RUN/build_logs/full.patch"
for component in smooth augmentation; do
    git apply --check "$STMFIT_PAW_ROOT/hpc/patches/qe-7.4.1-ildos-$component.patch"
    git apply "$STMFIT_PAW_ROOT/hpc/patches/qe-7.4.1-ildos-$component.patch"
    make -j 8 -C PP/src pp.x > "$STMFIT_PAW_RUN/build_logs/${component}_pp.out" 2>&1
    cp -p PP/src/pp.x "$STMFIT_PAW_RUN/binaries/pp_$component.x"
    cp PP/src/local_dos.f90 "$STMFIT_PAW_RUN/build_logs/local_dos.$component.f90"
    git diff -- PP/src/local_dos.f90 > "$STMFIT_PAW_RUN/build_logs/$component.patch"
    if [[ "$component" == smooth ]]; then
        git apply -R --check "$STMFIT_PAW_ROOT/hpc/patches/qe-7.4.1-ildos-$component.patch"
        git apply -R "$STMFIT_PAW_ROOT/hpc/patches/qe-7.4.1-ildos-$component.patch"
        cmp PP/src/local_dos.f90 "$STMFIT_PAW_RUN/build_logs/local_dos.full.f90"
    fi
done
cd "$STMFIT_PAW_RUN"
sha256sum binaries/pp_full.x binaries/pp_smooth.x binaries/pp_augmentation.x \
    build_logs/local_dos.full.f90 build_logs/local_dos.smooth.f90 build_logs/local_dos.augmentation.f90 > executable_sources.sha256
for molecule in glcn glcnac; do
    cd "$molecule"
    for observable in full_control full_ildos; do
        srun -n "$QE_NTASKS" --cpu-bind=cores ../binaries/pp_full.x -in "$observable.in" > "$observable.out"
        grep -q 'JOB DONE.' "$observable.out"
    done
    cmp full_control.cube reference.cube
    cmp full_ildos.cube reference_total.cube
    for component in smooth augmentation; do
        for observable in "${component}_control" "${component}_ildos" "repeat_${component}"; do
            srun -n "$QE_NTASKS" --cpu-bind=cores "../binaries/pp_$component.x" -in "$observable.in" > "$observable.out"
            grep -q 'JOB DONE.' "$observable.out"
        done
        cmp "${component}_control.cube" reference.cube
        cmp "${component}_ildos.cube" "repeat_${component}.cube"
    done
    cd "$STMFIT_PAW_RUN"
done
sha256sum -c --quiet checkpoint_before.sha256
sha256sum -c --quiet original_wfc_before.sha256
for molecule in glcn glcnac; do
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_PAW_ROOT" \
        "$STMFIT_PAW_ROOT/test/qe_paw_decomposition.jl" compare "$molecule" settings.toml
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_PAW_ROOT" \
    "$STMFIT_PAW_ROOT/test/verify_qe_paw_decomposition.jl" "$STMFIT_PAW_RUN" \
    > independent_verification.out 2>&1
echo 'PAW_COMPONENT_DIAGNOSTIC_COMPLETE; no mold, calibration, fit or benchmark produced'
