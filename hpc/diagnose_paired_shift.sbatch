#!/usr/bin/env bash
#SBATCH --job-name=stmfit-subpixel-pilot
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
: "${STMFIT_INPUT_DIR:?Set saved inputs directory}"
: "${STMFIT_OUTDIR:?Set fresh output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
args=(--native-dir "$STMFIT_INPUT_DIR/native" --reference-dir "$STMFIT_INPUT_DIR/reference"
    --physical-config config/chitosan.toml --assignment-config config/unit_assignment_patch_support.toml
    --model-settings config/paired_acquisition.toml --settings config/paired_convergence.toml
    --solver-settings config/paired_solver.toml --shift-settings config/paired_shift.toml)
if [[ $# -eq 1 ]]; then
    exec "$JULIA_BIN" --startup-file=no -t1 --project=. test/diagnose_paired_shift.jl \
        "${args[@]}" --outdir "$STMFIT_OUTDIR" --dry-run
fi
: "${SLURM_JOB_ID:?Run fits on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
mkdir -p "$STMFIT_OUTDIR/logs"
printf 'Job %s on %s; Julia %s; four single-thread shards\n' "$SLURM_JOB_ID" "${SLURMD_NODENAME:-unknown}" "$JULIA_BIN"
pids=()
for i in 1 2 3 4; do
    "$JULIA_BIN" --startup-file=no -t1 --project=. test/diagnose_paired_shift.jl \
        "${args[@]}" --chunk "$i/4" --outdir "$STMFIT_OUTDIR/chunk$i" > "$STMFIT_OUTDIR/logs/chunk$i.log" 2>&1 &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do if ! wait "$pid"; then failed=1; fi; done
[[ "$failed" -eq 0 ]] || { echo "Pilot failure; outputs retained; no retry" >&2; exit 1; }
"$JULIA_BIN" --startup-file=no -t1 --project=. test/diagnose_paired_shift.jl --merge "$STMFIT_OUTDIR"
