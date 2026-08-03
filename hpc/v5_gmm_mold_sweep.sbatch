#!/usr/bin/env bash
# Slurm job: GMM 1-view + mold margin refinement sweep (v5).
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V5, STMFIT_OUTPUT_DIR,
# JULIA_MODULE_VERSION.

#SBATCH --time=03:00:00
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

V5="${STMFIT_FEATURES_V5}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"

# baseline: base4 + uasym + mold_margin_avg (best GMM so far: 77.3/35)
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avg.tsv" \
    --view "v_avg=${BASE4},patch_u_asym,mold_margin_avg" --seeds 10 --interactions
# 3-channel mean margin
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avg3.tsv" \
    --view "v_avg3=${BASE4},patch_u_asym,mold_margin_avg3" --seeds 10 --interactions
# 3-channel mean + self-training 2/5 (st2 won +3 exact on the 3-view config)
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avg3_st2.tsv" \
    --view "v_avg3=${BASE4},patch_u_asym,mold_margin_avg3" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avg3_st5.tsv" \
    --view "v_avg3=${BASE4},patch_u_asym,mold_margin_avg3" --seeds 10 --interactions --selftrain 5
# seed scaling on the improved config
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avg3_s20.tsv" \
    --view "v_avg3=${BASE4},patch_u_asym,mold_margin_avg3" --seeds 20 --interactions
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avg3_s40.tsv" \
    --view "v_avg3=${BASE4},patch_u_asym,mold_margin_avg3" --seeds 40 --interactions
# avg + diff margin, no interactions (7 dims)
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v5_avgdiff_ni.tsv" \
    --view "v_avgd=${BASE4},patch_u_asym,mold_margin_avg,mold_margin_diff" --seeds 10 --no-interactions

echo "V5 SWEEP COMPLETE"
ls -la "${OUT}"
