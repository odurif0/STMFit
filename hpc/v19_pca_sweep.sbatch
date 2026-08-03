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
V="${STMFIT_FEATURES_V19}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"
predict() {
    local script="$1" feats="$2" out="$3"; shift 3
    echo "=== $(basename "$out") :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}
BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
CC="patch_u_asym,mold_cc_fwd,mold_cc_bwd"
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v19_ref.tsv" \
    --view "v_cc=${BASE4},${CC},emp_fwd" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v19_pca.tsv" \
    --view "v_cc=${BASE4},${CC},emp_pca_cv" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v19_both.tsv" \
    --view "v_cc=${BASE4},${CC},emp_fwd,emp_pca_cv" --seeds 10 --interactions --selftrain 2
echo "V19 SWEEP COMPLETE"
ls -la "${OUT}"
