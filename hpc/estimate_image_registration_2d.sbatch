#!/usr/bin/env bash
#SBATCH --job-name=stmfit-image-2d
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set source checkout}"
: "${STMFIT_INPUT_DIR:?Set image input directory}"
: "${STMFIT_OUTDIR:?Set fresh output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
args=(--raw-file "$STMFIT_INPUT_DIR/240817_006.sxm" --base-shifts "$STMFIT_INPUT_DIR/shifts.tsv"
    --config config/chitosan.toml --settings config/image_registration_2d.toml)
if [[ $# -eq 1 ]]; then
    exec "$JULIA_BIN" --startup-file=no -t1 --project=. test/estimate_image_registration_2d.jl \
        "${args[@]}" --outdir "$STMFIT_OUTDIR" --dry-run
fi
: "${SLURM_JOB_ID:?Run image searches on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
mkdir -p "$STMFIT_OUTDIR/logs"
printf 'Job %s on %s; Julia %s; three single-thread image folds\n' "$SLURM_JOB_ID" "${SLURMD_NODENAME:-unknown}" "$JULIA_BIN"
pids=()
for fold in 0 1 2; do
    "$JULIA_BIN" --startup-file=no -t1 --project=. test/estimate_image_registration_2d.jl \
        "${args[@]}" --fold "$fold" --outdir "$STMFIT_OUTDIR/fold$fold" > "$STMFIT_OUTDIR/logs/fold$fold.log" 2>&1 &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do if ! wait "$pid"; then failed=1; fi; done
[[ "$failed" -eq 0 ]] || { echo "Image registration failed; no retry" >&2; exit 1; }
"$JULIA_BIN" --startup-file=no -t1 --project=. test/estimate_image_registration_2d.jl --merge "$STMFIT_OUTDIR"
"$JULIA_BIN" --startup-file=no -t1 --project=. test/verify_image_registration_2d.jl "$STMFIT_OUTDIR"
