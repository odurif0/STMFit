#!/usr/bin/env bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000MB
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -euo pipefail
module purge
module load "julia/${JULIA_MODULE_VERSION}"
cd "${STMFIT_PROJECT_DIR}"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
V="${STMFIT_FEATURES_V17}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"
predict() {
    local script="$1" feats="$2" out="$3"; shift 3
    echo "=== $(basename "$out") :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}
BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
CC="patch_u_asym,mold_cc_fwd,mold_cc_bwd"
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v17_ref.tsv" \
    --view "v_cc=${BASE4},${CC},emp_fwd" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v17_avg.tsv" \
    --view "v_cc=${BASE4},${CC},emp_avg" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v17_rawavg.tsv" \
    --view "v_cc=${BASE4},${CC},emp_raw_avg" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v17_all_ni.tsv" \
    --view "v_all=${BASE4},${CC},emp_fwd,emp_bwd,emp_raw_fwd,emp_raw_bwd,emp_diff" \
    --seeds 10 --no-interactions --selftrain 2
echo "V17 SWEEP COMPLETE"
ls -la "${OUT}"
