#!/usr/bin/env bash
# Slurm job: k-means seed scaling (corrected: with --interactions) +
# GMM 3-view self-training=2 confirmation at 80 seeds.
# Expected environment: STMFIT_PROJECT_DIR, STMFIT_FEATURES_SPLIT,
# STMFIT_FEATURES_UASYM, STMFIT_OUTPUT_DIR, JULIA_MODULE_VERSION.

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

UASYM="${STMFIT_FEATURES_UASYM}"
SPLIT="${STMFIT_FEATURES_SPLIT}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"

# --- k-means seed scaling, 4-view default config WITH interactions ---
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_i_s20.tsv" --seeds 20 --interactions
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_i_s80.tsv" --seeds 80 --interactions
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_i_s200.tsv" --seeds 200 --interactions
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_i_s500.tsv" --seeds 500 --interactions

# --- GMM 3-view, selftrain=2, 80 seeds (confirm the 29 -> 32 exact effect) ---
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_3view_st2_s80.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --view "base_comt=${BASE4},bwd_neg_com_t" \
    --view "base_diag45=${BASE4},bwd_neg_diag45" --seeds 80 --interactions --selftrain 2

echo "KM SWEEP COMPLETE"
ls -la "${OUT}"
