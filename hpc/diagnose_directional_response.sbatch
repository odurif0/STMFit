#!/usr/bin/env bash
#SBATCH --job-name=stmfit-response
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set source checkout}"
: "${STMFIT_INPUT_DIR:?Set input directory}"
: "${STMFIT_OUTDIR:?Set fresh output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
args=(--data-dir "$STMFIT_INPUT_DIR" --base-shifts "$STMFIT_INPUT_DIR/shifts.tsv"
    --config config/chitosan.toml --settings config/directional_response.toml --outdir "$STMFIT_OUTDIR")
if [[ $# -eq 1 ]]; then
    exec "$JULIA_BIN" --startup-file=no -t1 --project=. test/diagnose_directional_response.jl "${args[@]}" --dry-run
fi
: "${SLURM_JOB_ID:?Run real panel searches on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
printf 'Job %s on %s; four Julia threads, one scan per task\n' "$SLURM_JOB_ID" "${SLURMD_NODENAME:-unknown}"
"$JULIA_BIN" --startup-file=no -t4 --project=. test/diagnose_directional_response.jl "${args[@]}"
"$JULIA_BIN" --startup-file=no -t1 --project=. test/verify_directional_response.jl "$STMFIT_OUTDIR"
