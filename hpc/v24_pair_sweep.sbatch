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
V="${STMFIT_FEATURES_V24}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"
predict() {
    local script="$1" feats="$2" out="$3"; shift 3
    echo "=== $(basename "$out") :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}
BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v24_ref.tsv" \
    --view "v_cc=${BASE4},patch_u_asym,mold_cc_fwd,mold_cc_bwd,emp_fisher" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v24_pair.tsv" \
    --view "v_cc=${BASE4},patch_u_asym,mold_cc_fwd,mold_cc_bwd,emp_fisher,pair_corr_diff" \
    --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v24_pair3.tsv" \
    --view "v_cc=${BASE4},patch_u_asym,mold_cc_fwd,mold_cc_bwd,emp_fisher,pair_corr_next,pair_corr_prev,pair_corr_diff" \
    --seeds 10 --interactions --selftrain 2
echo "V24 SWEEP COMPLETE"
ls -la "${OUT}"
