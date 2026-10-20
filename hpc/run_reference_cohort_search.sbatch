#!/usr/bin/env bash
#SBATCH --job-name=stmfit-cohort-search
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=120000MB
set -euo pipefail
# Campaign B (benchmarks/reference_cohort_search_v1.toml): every search of the
# protocol as an independent process of scripts/grading/search_reference_cohort.jl
# (13 nested searches, the final search, the null searches), STMFIT_PARALLEL
# processes at a time with STMFIT_THREADS threads each (default 16 x 8), then
# the merged report. Grading side: reads benchmark labels.
# Required (absolute paths, passed with sbatch --export=NONE,...):
# STMFIT_PROJECT_DIR (committed checkout holding the protocol's run tables
# under results/), STMFIT_OUTDIR (must not exist), JULIA_BIN.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set committed source checkout}"
: "${STMFIT_OUTDIR:?Set new output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
PARALLEL=${STMFIT_PARALLEL:-16}; THREADS=${STMFIT_THREADS:-8}
PROTOCOL=benchmarks/reference_cohort_search_v1.toml
cd "$STMFIT_PROJECT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mapfile -t TASKS < <("$JULIA_BIN" --startup-file=no --project=. -e '
    using TOML; p = TOML.parsefile(ARGS[1]); m = TOML.parsefile(p["data"]["manifest"])["files"]
    run = p["data"]["run"]
    for f in ("assignment/features_predictor.tsv", "assignment/patches_bwd17.tsv",
              "consensus/consensus_summary_final.tsv", "fusion/predictions_fused.tsv")
        isfile(joinpath(run, f)) || error("Missing run table: $(joinpath(run, f))")
    end
    println("final")
    for s in sort(unique(v["session"] for v in values(m) if v["cohort"] in p["split"]["development_cohorts"]))
        println("nested --outer ", s)
    end
    foreach(k -> println("null --seed ", k), p["null"]["seeds"])' "$PROTOCOL")
echo "${#TASKS[@]} searches, $PARALLEL at a time with $THREADS threads; output $STMFIT_OUTDIR"
printf '  %s\n' "${TASKS[@]}"
if [[ $# -eq 1 ]]; then echo "Dry run: nothing computed."; exit 0; fi
: "${SLURM_JOB_ID:?Run on a compute node}"
mkdir -p "$STMFIT_OUTDIR/logs"
export JULIA_BIN PROTOCOL THREADS STMFIT_OUTDIR
printf '%s\n' "${TASKS[@]}" | xargs -P "$PARALLEL" -I{} bash -c '
    set -- {}; mode=$1; shift
    name=$mode${2:+_$2}
    "$JULIA_BIN" --startup-file=no -t "$THREADS" --project=. scripts/grading/search_reference_cohort.jl \
        --protocol "$PROTOCOL" --mode "$mode" "$@" --outdir "$STMFIT_OUTDIR/$name" > "$STMFIT_OUTDIR/logs/$name.log" 2>&1 ||
        { echo "FAILED $name"; exit 1; }
    echo "done $name"'
"$JULIA_BIN" --startup-file=no -t 1 --project=. scripts/grading/search_reference_cohort.jl \
    --protocol "$PROTOCOL" --mode report --campaign "$STMFIT_OUTDIR" --outdir "$STMFIT_OUTDIR/report"
