#!/usr/bin/env bash
#SBATCH --job-name=stmfit-robust-normalization
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail

# One control and one robust-normalization candidate at fixed geometry and N.
# Both regenerate patches; only per-scan GMM feature normalization changes.
# Set absolute paths via an explicit sbatch --export list on Viper:
# STMFIT_PROJECT_DIR, STMFIT_CACHE_DIR (features.tsv, features_split.tsv),
# STMFIT_INPUT_DIR (full146_raw/, templates_cc.tsv), STMFIT_OUTDIR (new), JULIA_BIN.
# First run --dry-run; actual execution requires a Slurm allocation.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2
    exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set the source checkout}"
: "${STMFIT_CACHE_DIR:?Set the saved symmetric-fusion geometry}"
: "${STMFIT_INPUT_DIR:?Set the raw/template input directory}"
: "${STMFIT_OUTDIR:?Set a new output directory}"
: "${JULIA_BIN:?Set the Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_CACHE_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required: $path" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || {
    echo "Output already exists: $STMFIT_OUTDIR" >&2; exit 2;
}
for file in features.tsv features_split.tsv; do
    [[ -f "$STMFIT_CACHE_DIR/$file" ]] || { echo "Missing geometry: $file" >&2; exit 2; }
done
[[ -d "$STMFIT_INPUT_DIR/full146_raw" && -f "$STMFIT_INPUT_DIR/templates_cc.tsv" ]] || {
    echo "Missing raw/template inputs" >&2; exit 2;
}
cd "$STMFIT_PROJECT_DIR"
if [[ $# -eq 0 ]]; then
    : "${SLURM_JOB_ID:?Scientific comparison must run in a Slurm job}"
    [[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
    mkdir "$STMFIT_OUTDIR"
fi
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=4
status=0
for arm in control robust; do
    case "$arm" in
        control) config=config/unit_assignment_patch_support.toml ;;
        robust) config=config/unit_assignment_robust_normalization.toml ;;
    esac
    command=("$JULIA_BIN" --startup-file=no -t 4 --project=. test/run_reconstructed_chitosan.jl
        --data-dir "$STMFIT_INPUT_DIR/full146_raw" --count-config config/chitosan.toml
        --config "$config" --templates "$STMFIT_INPUT_DIR/templates_cc.tsv"
        --features "$STMFIT_CACHE_DIR/features.tsv" --split-features "$STMFIT_CACHE_DIR/features_split.tsv"
        --outdir "$STMFIT_OUTDIR/$arm")
    if [[ $# -eq 1 ]]; then
        "${command[@]}" --dry-run
    else
        echo "Starting $arm at $(date -Is)"
        if "${command[@]}" > "$STMFIT_OUTDIR/$arm.log" 2>&1; then
            echo "$arm completed at $(date -Is)"
        else
            code=$?
            echo "$arm failed (exit $code) at $(date -Is); see $STMFIT_OUTDIR/$arm.log" >&2
            status=1
        fi
    fi
done
exit "$status"
