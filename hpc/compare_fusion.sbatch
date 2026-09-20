#!/usr/bin/env bash
#SBATCH --job-name=stmfit-fusion
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail

# One paired comparison: two code checkouts, identical cached label-free counts,
# regenerated base/split geometry and patches in both arms, whole-cohort native
# assignment. No external labels or grading in this compute job.
# Required absolute paths: STMFIT_CONTROL_PROJECT, STMFIT_SYMMETRIC_PROJECT,
# STMFIT_INPUT_DIR (full146_raw/ and templates_cc.tsv), STMFIT_SELECTED_SUMMARY,
# STMFIT_OUTDIR (new directory), JULIA_BIN (1.13).
# On Viper, supply all six variables in sbatch's explicit --export list.
# Run this script with --dry-run before submitting. Scientific execution needs
# a Slurm allocation; each arm uses at most four one-thread fit subprocesses.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2
    exit 2
fi
: "${STMFIT_CONTROL_PROJECT:?Set the unchanged control checkout}"
: "${STMFIT_SYMMETRIC_PROJECT:?Set the corrected checkout}"
: "${STMFIT_INPUT_DIR:?Set the common raw/template directory}"
: "${STMFIT_SELECTED_SUMMARY:?Set the original cached-N summary}"
: "${STMFIT_OUTDIR:?Set a new paired output directory}"
: "${JULIA_BIN:?Set the Julia 1.13 executable}"
for path in "$STMFIT_CONTROL_PROJECT" "$STMFIT_SYMMETRIC_PROJECT" "$STMFIT_INPUT_DIR" \
            "$STMFIT_SELECTED_SUMMARY" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required: $path" >&2; exit 2; }
done
[[ "$(realpath "$STMFIT_CONTROL_PROJECT")" != "$(realpath "$STMFIT_SYMMETRIC_PROJECT")" ]] || {
    echo "Control and treatment must use separate source directories" >&2; exit 2;
}
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || {
    echo "Output already exists: $STMFIT_OUTDIR" >&2; exit 2;
}
for file in Project.toml Manifest.toml config/chitosan.toml config/unit_assignment_reconstructed.toml; do
    cmp "$STMFIT_CONTROL_PROJECT/$file" "$STMFIT_SYMMETRIC_PROJECT/$file"
done
[[ -f "$STMFIT_SELECTED_SUMMARY" && -f "$STMFIT_INPUT_DIR/templates_cc.tsv" \
   && -d "$STMFIT_INPUT_DIR/full146_raw" ]] || { echo "Missing common inputs" >&2; exit 2; }
if [[ $# -eq 0 ]]; then
    : "${SLURM_JOB_ID:?Scientific comparison must run in a Slurm job}"
    [[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
    mkdir "$STMFIT_OUTDIR"
fi
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=4
status=0
for arm in control symmetric; do
    if [[ "$arm" == control ]]; then
        project="$STMFIT_CONTROL_PROJECT"
    else
        project="$STMFIT_SYMMETRIC_PROJECT"
    fi
    cd "$project"
    command=("$JULIA_BIN" --startup-file=no -t 4 --project=. test/run_reconstructed_chitosan.jl
        --data-dir "$STMFIT_INPUT_DIR/full146_raw"
        --count-config config/chitosan.toml --config config/unit_assignment_reconstructed.toml
        --selected-summary "$STMFIT_SELECTED_SUMMARY"
        --templates "$STMFIT_INPUT_DIR/templates_cc.tsv" --outdir "$STMFIT_OUTDIR/$arm")
    if [[ $# -eq 1 ]]; then
        "${command[@]}" --dry-run
    else
        echo "Starting $arm in $project at $(date -Is)"
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
