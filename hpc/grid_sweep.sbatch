#!/usr/bin/env bash
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
F="${STMFIT_FEATURES}"
D="${STMFIT_DATA_DIR}"
O="${STMFIT_OUTPUT_DIR}"
mkdir -p "${O}"
for spec in "16x16:0.30:0.30" "17x18:0.32:0.34" "18x18:0.34:0.34" "16x17:0.30:0.32"; do
  name="${spec%%:*}"; rest="${spec#*:}"
  half="${rest%%:*}"; halfu="${rest##*:}"
  echo "=== $name (t $half, u $halfu) ==="
  julia --project=. test/extract_lobe_patches.jl \
      --features "${F}" --data-dir "${D}" \
      --half-nm "${half}" --half-u-nm "${halfu}" --step-nm 0.04 \
      --out "${O}/lobe_patches_${name}_full145.tsv"
done
echo "GRID SWEEP COMPLETE"
ls -la "${O}"
