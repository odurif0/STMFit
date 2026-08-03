#!/usr/bin/env bash
# Slurm job: label-free unit-assignment prediction sweep.
#   - GMM self-training (Mahalanobis, 1-view and 3-view configs)
#   - k-means seed scaling (20/80/200/500 seeds on the 4-view config)
# Each step writes a prediction TSV into STMFIT_OUTPUT_DIR; grading happens
# locally after fetch (post-hoc, label-free until then).
#
# Expected environment (set by hpc/launch_prediction_sweep.sh):
#   STMFIT_PROJECT_DIR       synced repo root on the compute node
#   STMFIT_FEATURES_UASYM    feature TSV incl. patch_u_asym (full145)
#   STMFIT_FEATURES_SPLIT    feature TSV incl. split_log_skew (full145)
#   STMFIT_OUTPUT_DIR        remote output dir (fetched back afterwards)
#   JULIA_MODULE_VERSION     Julia module to load

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

# --- GMM self-training, journal config: single view (base + patch_u_asym), 10 seeds ---
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_1view_st0.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_1view_st5.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions --selftrain 5
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_1view_st10.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions --selftrain 10

# --- GMM self-training, full145 3-view config (20 seeds) ---
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_3view_st0.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --view "base_comt=${BASE4},bwd_neg_com_t" \
    --view "base_diag45=${BASE4},bwd_neg_diag45" --seeds 20 --interactions
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_3view_st2.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --view "base_comt=${BASE4},bwd_neg_com_t" \
    --view "base_diag45=${BASE4},bwd_neg_diag45" --seeds 20 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_3view_st10.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --view "base_comt=${BASE4},bwd_neg_com_t" \
    --view "base_diag45=${BASE4},bwd_neg_diag45" --seeds 20 --interactions --selftrain 10

# --- GMM seed scaling, 3-view config ---
predict build_labelfree_gmm_predictions.jl "${UASYM}" "${OUT}/pred_gmm_3view_s80.tsv" \
    --view "base_uasym=${BASE4},patch_u_asym" --view "base_comt=${BASE4},bwd_neg_com_t" \
    --view "base_diag45=${BASE4},bwd_neg_diag45" --seeds 80 --interactions

# --- k-means seed scaling, default 4-view config (split table) ---
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_s20.tsv" --seeds 20
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_s80.tsv" --seeds 80
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_s200.tsv" --seeds 200
predict build_labelfree_unit_predictions.jl "${SPLIT}" "${OUT}/pred_km_4view_s500.tsv" --seeds 500

echo "SWEEP COMPLETE"
ls -la "${OUT}"
