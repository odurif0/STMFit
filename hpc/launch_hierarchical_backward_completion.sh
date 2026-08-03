#!/usr/bin/env bash

# Complete the hierarchical unit-assignment lane by adding backward-view
# features that were missing from the original T5 run.
#
# Flow:
#   1. Sync code + base feature TSV to Viper.
#   2. On Viper: extract backward patches (extract_lobe_patches_bwd.jl).
#   3. On Viper: build combined feature TSV (build_combined_hierarchical_features.jl).
#   4. Submit the T5 leave-date-out/bootstrap array with the combined TSV.
#   5. Merge shards + fetch results.
#
# This wrapper reuses the existing hierarchical_unit_assignment_array.sbatch
# and merge_hierarchical_unit_assignment_shards.jl.  The only change is that
# STMFIT_HIERARCHICAL_FEATURES points to the *combined* TSV (with bwd_neg_*
# columns) so that _default_views activates the backward views and
# forward_backward_agreement is no longer NaN.
#
# Usage:
#   STMFIT_HIERARCHICAL_FEATURES=<base_features.tsv> \
#   ./hpc/launch_hierarchical_backward_completion.sh --dry-run

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_FILE="$SCRIPT_DIR/remote.env"

usage() {
    cat <<'EOF'
Usage: launch_hierarchical_backward_completion.sh --dry-run|--submit|--watch

Complete the hierarchical backward-view lane on Viper.

Exactly one mode is mandatory:
  --dry-run  Resolve and print all paths, hashes, commands, and array shape.
  --submit   Sync, extract patches, build combined features, submit array.
  --watch    Sync, extract, build, submit, wait, merge, fetch, validate.

Environment:
  STMFIT_HIERARCHICAL_FEATURES   Local base label-free feature TSV (required).
  STMFIT_HIERARCHICAL_CONFIG     Candidate config [config/unit_assignment_candidate.toml].
  STMFIT_HIERARCHICAL_LOCAL_OUT  Fetch base [results/hierarchical_unit_assignment].
  STMFIT_HIERARCHICAL_SPLIT      Local split-width feature TSV (optional).
  STMFIT_HIERARCHICAL_PATCH_HALF_NM  Backward patch half-size [0.32].
  STMFIT_HIERARCHICAL_PATCH_STEP_NM  Backward patch grid spacing [0.08].
EOF
}

(( $# == 1 )) || { usage >&2; exit 2; }
case "$1" in
    --dry-run) MODE=dry-run ;;
    --submit)  MODE=submit ;;
    --watch)   MODE=watch ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: expected exactly one of --dry-run|--submit|--watch" >&2; usage >&2; exit 2 ;;
esac

[[ -f "$ENV_FILE" ]] || {
    echo "ERROR: $ENV_FILE not found; copy hpc/remote.env.example first" >&2
    exit 1
}
# shellcheck source=/dev/null
source "$ENV_FILE"

: "${STMFIT_REMOTE_USER:?STMFIT_REMOTE_USER is required in hpc/remote.env}"
: "${STMFIT_HIERARCHICAL_FEATURES:?STMFIT_HIERARCHICAL_FEATURES must name the local base feature TSV}"
: "${STMFIT_HIERARCHICAL_CONFIG:=config/unit_assignment_candidate.toml}"
: "${STMFIT_HIERARCHICAL_LOCAL_OUT:=results/hierarchical_unit_assignment}"
: "${STMFIT_HIERARCHICAL_SPLIT:=}"
: "${STMFIT_HIERARCHICAL_PATCH_HALF_NM:=0.32}"
: "${STMFIT_HIERARCHICAL_PATCH_STEP_NM:=0.08}"
: "${JULIA_MODULE_VERSION:=1.12}"
: "${WATCH_POLL:=30}"
: "${SSH_CONNECT_TIMEOUT:=180}"
: "${SSH_SERVER_ALIVE_INTERVAL:=60}"
: "${RSYNC_EXTRA:=}"

SSH_HOST="${STMFIT_HIERARCHICAL_SSH_HOST:-viper}"
[[ "$SSH_HOST" == "viper" ]] || {
    echo "ERROR: the backward-completion gate must target the Viper SSH alias" >&2
    exit 1
}
[[ -f "$STMFIT_HIERARCHICAL_FEATURES" ]] || {
    echo "ERROR: base feature TSV not found: $STMFIT_HIERARCHICAL_FEATURES" >&2
    exit 1
}
CONFIG_PATH="$STMFIT_HIERARCHICAL_CONFIG"
[[ "$CONFIG_PATH" = /* ]] || CONFIG_PATH="$REPO_ROOT/$CONFIG_PATH"
[[ -f "$CONFIG_PATH" ]] || { echo "ERROR: config not found: $CONFIG_PATH" >&2; exit 1; }

sha_file() { sha256sum "$1" | cut -d' ' -f1; }
CONFIG_SHA256=$(sha_file "$CONFIG_PATH")
BASE_INPUT_SHA256=$(sha_file "$STMFIT_HIERARCHICAL_FEATURES")

# Source manifest (same as the existing launcher).
SOURCE_MANIFEST=$(mktemp)
cleanup() { rm -f "$SOURCE_MANIFEST"; }
trap cleanup EXIT
{
    printf 'test/evaluate_hierarchical_unit_assignment.jl\t%s\n' \
        "$(sha_file "$REPO_ROOT/test/evaluate_hierarchical_unit_assignment.jl")"
    printf 'test/lib/hierarchical_unit_assignment.jl\t%s\n' \
        "$(sha_file "$REPO_ROOT/test/lib/hierarchical_unit_assignment.jl")"
    for path in "$REPO_ROOT"/test/lib/hierarchical/*; do
        [[ -f "$path" ]] || continue
        printf 'test/lib/hierarchical/%s\t%s\n' "$(basename "$path")" "$(sha_file "$path")"
    done
    printf 'test/build_combined_hierarchical_features.jl\t%s\n' \
        "$(sha_file "$REPO_ROOT/test/build_combined_hierarchical_features.jl")"
    printf 'test/extract_lobe_patches_bwd.jl\t%s\n' \
        "$(sha_file "$REPO_ROOT/test/extract_lobe_patches_bwd.jl")"
} > "$SOURCE_MANIFEST"
SOURCE_SHA256=$(sha_file "$SOURCE_MANIFEST")

# Resolve fold dates from the base feature TSV.
FOLDS_CSV=$(julia --project="$REPO_ROOT" -e '
include(joinpath(ARGS[1], "test", "evaluate_hierarchical_unit_assignment.jl"))
E = HierarchicalUnitAssignmentEvaluation
paths, _ = E._raw_scan_paths(ARGS[2])
print(join((f.fold_date for f in E.build_date_folds(paths)), ","))
' "$REPO_ROOT" "$STMFIT_HIERARCHICAL_FEATURES")
[[ -n "$FOLDS_CSV" ]] || { echo "ERROR: no leave-date-out folds resolved" >&2; exit 1; }
IFS=',' read -r -a FOLDS <<< "$FOLDS_CSV"
FOLDS_ENV=${FOLDS_CSV//,/:}
N_FOLDS=${#FOLDS[@]}
BLOCKS_PER_FOLD=11
ARRAY_TASKS=$((N_FOLDS * BLOCKS_PER_FOLD))

# Paths.
REMOTE_PROJECT_BASE="${REMOTE_PROJECT_BASE:-/u/$STMFIT_REMOTE_USER/code}"
REMOTE_PROJECT="$REMOTE_PROJECT_BASE/STMFit"
REMOTE_DATA_DIR="${REMOTE_DATA_DIR:-/ptmp/$STMFIT_REMOTE_USER/stmfit/data}"
RUN_ID="bwd-${BASE_INPUT_SHA256:0:12}-${CONFIG_SHA256:0:12}-${SOURCE_SHA256:0:12}"
REMOTE_RUN="/ptmp/$STMFIT_REMOTE_USER/stmfit/hierarchical_backward/$RUN_ID"
REMOTE_FEATURES_BASE="$REMOTE_RUN/inputs/features_base.tsv"
REMOTE_PATCHES="$REMOTE_RUN/inputs/patches_bwd.tsv"
REMOTE_SPLIT="$REMOTE_RUN/inputs/split.tsv"
REMOTE_COMBINED="$REMOTE_RUN/inputs/features_combined.tsv"
REMOTE_OUT="$REMOTE_RUN/output"
LOCAL_RUN="$REPO_ROOT/$STMFIT_HIERARCHICAL_LOCAL_OUT/$RUN_ID"
REMOTE_CONFIG_REL="config/unit_assignment_candidate.toml"

SSH_OPTS=(-o BatchMode=no -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
          -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL")
run_remote() { ssh "${SSH_OPTS[@]}" "$SSH_HOST" "$1"; }

# --- Combined-feature INPUT hash (computed after the build on Viper) ---
# For the dry-run we use the base hash as a placeholder.
COMBINED_INPUT_SHA256="<computed-on-viper>"

# Shard + merge args (same structure as the existing launcher).
SHARD_ARGS=()
for fold in "${FOLDS[@]}"; do
    SHARD_ARGS+=(--shard "$REMOTE_OUT/shards/fold_${fold}_block_unbootstrapped.tsv")
    for start in 0 50 100 150 200 250 300 350 400 450; do
        stop=$((start + 49))
        SHARD_ARGS+=(--shard "$(printf '%s/shards/fold_%s_block_%03d-%03d.tsv' "$REMOTE_OUT" "$fold" "$start" "$stop")")
    done
done
printf -v SHARD_ARGS_Q ' %q' "${SHARD_ARGS[@]}"
MERGED_REMOTE="$REMOTE_OUT/merged.tsv"

# Commands executed on Viper.
EXTRACT_CMD="cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && STMFIT_DATA_DIR=$(printf '%q' "$REMOTE_DATA_DIR") julia --project=. test/extract_lobe_patches_bwd.jl --features $(printf '%q' "$REMOTE_FEATURES_BASE") --out $(printf '%q' "$REMOTE_PATCHES") --data-dir $(printf '%q' "$REMOTE_DATA_DIR") --half-nm $(printf '%q' "$STMFIT_HIERARCHICAL_PATCH_HALF_NM") --step-nm $(printf '%q' "$STMFIT_HIERARCHICAL_PATCH_STEP_NM")"

SPLIT_ARG=""
[[ -z "$STMFIT_HIERARCHICAL_SPLIT" ]] || SPLIT_ARG="--split $(printf '%q' "$REMOTE_SPLIT")"
BUILD_CMD="cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && julia --project=. test/build_combined_hierarchical_features.jl --features $(printf '%q' "$REMOTE_FEATURES_BASE") --patches $(printf '%q' "$REMOTE_PATCHES") $SPLIT_ARG --out $(printf '%q' "$REMOTE_COMBINED") 1>&2 && sha256sum $(printf '%q' "$REMOTE_COMBINED") | cut -d' ' -f1"

MERGE_COMMAND="cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && julia --project=. test/merge_hierarchical_unit_assignment_shards.jl${SHARD_ARGS_Q} --out $(printf '%q' "$MERGED_REMOTE") --expected-folds $(printf '%q' "$FOLDS_CSV") --config-sha256 $CONFIG_SHA256 --input-sha256 \$COMBINED_SHA256 --source-sha256 $SOURCE_SHA256"

EXPORTS="ALL,STMFIT_PROJECT_DIR=$REMOTE_PROJECT,STMFIT_HIERARCHICAL_FEATURES=$REMOTE_COMBINED,STMFIT_HIERARCHICAL_CONFIG=$REMOTE_CONFIG_REL,STMFIT_HIERARCHICAL_OUTBASE=$REMOTE_OUT,STMFIT_HIERARCHICAL_FOLDS=$FOLDS_ENV,STMFIT_HIERARCHICAL_CONFIG_SHA256=$CONFIG_SHA256,STMFIT_HIERARCHICAL_INPUT_SHA256=\$COMBINED_SHA256,STMFIT_HIERARCHICAL_SOURCE_SHA256=$SOURCE_SHA256,JULIA_MODULE_VERSION=$JULIA_MODULE_VERSION"
SBATCH_COMMAND="cd $(printf '%q' "$REMOTE_PROJECT") && COMBINED_SHA256=\$(sha256sum $(printf '%q' "$REMOTE_COMBINED") | cut -d' ' -f1) && sbatch --parsable --array=1-${ARRAY_TASKS}%8 --cpus-per-task=1 --export=$(printf '%q' "$EXPORTS") hpc/hierarchical_unit_assignment_array.sbatch"

cat <<EOF
T5 hierarchical backward-completion gate
  mode:                   $MODE
  host:                   $SSH_HOST
  remote project:         $REMOTE_PROJECT
  remote data dir:        $REMOTE_DATA_DIR
  remote owned run:       $REMOTE_RUN
  remote base features:   $REMOTE_FEATURES_BASE
  remote backward patches:$REMOTE_PATCHES
  remote combined feats:  $REMOTE_COMBINED
  remote output:          $REMOTE_OUT
  local fetch:            $LOCAL_RUN
  folds:                  $FOLDS_CSV ($N_FOLDS)
  bootstrap blocks:       unbootstrapped + 0:49,...,450:499
  array:                  1-${ARRAY_TASKS}%8
  config sha256:          $CONFIG_SHA256
  base input sha256:      $BASE_INPUT_SHA256
  combined input sha256:  $COMBINED_INPUT_SHA256
  source sha256:          $SOURCE_SHA256
  extract command:        $EXTRACT_CMD
  build command:          $BUILD_CMD
  submit command:         $SBATCH_COMMAND
  merge command:          $MERGE_COMMAND
EOF

if [[ "$MODE" == dry-run ]]; then
    echo "DRY RUN PASS: no remote or local state changed"
    exit 0
fi

# --- Execute ---

run_remote "mkdir -p $(printf '%q' "$REMOTE_PROJECT") $(printf '%q' "$REMOTE_RUN/inputs") $(printf '%q' "$REMOTE_OUT/shards") && printf '%s\n' $(printf '%q' "$RUN_ID") > $(printf '%q' "$REMOTE_RUN/.stmfit-bwd-owned")"

# Sync code.
rsync -az --exclude='.git/' --exclude='Manifest.toml' --exclude='results/' \
    --exclude='hpc/remote.env' ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$REPO_ROOT/" "$SSH_HOST:$REMOTE_PROJECT/"

# Sync base feature TSV.
rsync -az ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$STMFIT_HIERARCHICAL_FEATURES" "$SSH_HOST:$REMOTE_FEATURES_BASE"

# Sync optional split TSV.
if [[ -n "$STMFIT_HIERARCHICAL_SPLIT" && -f "$STMFIT_HIERARCHICAL_SPLIT" ]]; then
    rsync -az ${RSYNC_EXTRA} \
        -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
        "$STMFIT_HIERARCHICAL_SPLIT" "$SSH_HOST:$REMOTE_SPLIT"
fi

# Instantiate on login node (internet available).
run_remote "cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && julia --project=. -e 'using Pkg; for p in [\"STMFitCore\",\"STMSXMIO\",\"GaussianFit1D\",\"GaussianFit2D\",\"STMMolecularFit\"]; Pkg.develop(PackageSpec(path=\"packages/\" * p * \".jl\")); end; Pkg.instantiate(); Pkg.precompile()'"

echo "--- extracting backward patches on Viper ---"
run_remote "$EXTRACT_CMD"

echo "--- building combined feature table on Viper ---"
COMBINED_SHA256=$(run_remote "$BUILD_CMD")
[[ "$COMBINED_SHA256" =~ ^[a-f0-9]{64}$ ]] || {
    echo "ERROR: combined feature build did not return a valid SHA-256: $COMBINED_SHA256" >&2
    exit 1
}
echo "combined_sha256=$COMBINED_SHA256"

echo "--- submitting T5 array with combined features ---"
JOB_ID=$(run_remote "COMBINED_SHA256=$COMBINED_SHA256 $SBATCH_COMMAND")
JOB_ID=${JOB_ID%%;*}
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || { echo "ERROR: could not parse Slurm job id: $JOB_ID" >&2; exit 1; }
echo "array_job_id=$JOB_ID"

if [[ "$MODE" == submit ]]; then
    echo "SUBMITTED: rerun --watch only if a new complete run is intended"
    exit 0
fi

# --- Watch: poll, merge, fetch ---
while :; do
    ACTIVE=$(run_remote "squeue -h -j $(printf '%q' "$JOB_ID") 2>/dev/null | wc -l")
    [[ "$ACTIVE" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid squeue response: $ACTIVE" >&2; exit 1; }
    (( ACTIVE == 0 )) && break
    printf 'waiting: job=%s active=%s\n' "$JOB_ID" "$ACTIVE"
    sleep "$WATCH_POLL"
done
SACCT_OUTPUT=$(run_remote "sacct -j $(printf '%q' "$JOB_ID") --format=JobIDRaw,State,ExitCode -n -P")
printf '%s\n' "$SACCT_OUTPUT"

array_total=0; array_failed=0; plain_total=0; plain_failed=0
while IFS='|' read -r job_raw state exit_code _; do
    case "$job_raw" in
        "${JOB_ID}_"*)
            array_total=$((array_total + 1))
            [[ "$state" == COMPLETED && "$exit_code" == 0:0 ]] ||
                array_failed=$((array_failed + 1))
            ;;
        ''|*[!0-9]*) ;;
        *)
            plain_total=$((plain_total + 1))
            [[ "$state" == COMPLETED && "$exit_code" == 0:0 ]] ||
                plain_failed=$((plain_failed + 1))
            ;;
    esac
done <<< "$SACCT_OUTPUT"

if (( array_total > 0 )); then
    accounting_total=$array_total; accounting_failed=$array_failed
else
    accounting_total=$plain_total; accounting_failed=$plain_failed
fi
(( accounting_total == ARRAY_TASKS )) || {
    echo "ERROR: expected $ARRAY_TASKS accounting rows, found $accounting_total" >&2; exit 1; }
(( accounting_failed == 0 )) || {
    echo "ERROR: $accounting_failed of $accounting_total tasks did not exit 0:0" >&2; exit 1; }
echo "accounting_pass=tasks:$accounting_total"

run_remote "COMBINED_SHA256=$COMBINED_SHA256 $MERGE_COMMAND"

mkdir -p "$LOCAL_RUN"
rsync -az --update ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$SSH_HOST:$REMOTE_OUT/" "$LOCAL_RUN/"

LOCAL_SHARDS=()
for remote_path in "${SHARD_ARGS[@]}"; do
    [[ "$remote_path" == --shard ]] && continue
    LOCAL_SHARDS+=(--shard "$LOCAL_RUN/${remote_path#"$REMOTE_OUT/"}")
done
julia --project="$REPO_ROOT" "$REPO_ROOT/test/merge_hierarchical_unit_assignment_shards.jl" \
    "${LOCAL_SHARDS[@]}" --out "$LOCAL_RUN/validated.tsv" \
    --expected-folds "$FOLDS_CSV" --config-sha256 "$CONFIG_SHA256" \
    --input-sha256 "$COMBINED_SHA256" --source-sha256 "$SOURCE_SHA256"

cmp "$LOCAL_RUN/merged.tsv" "$LOCAL_RUN/validated.tsv"
cmp "$LOCAL_RUN/merged.tsv.gate.tsv" "$LOCAL_RUN/validated.tsv.gate.tsv"
grep -q '^PASS' "$LOCAL_RUN/validated.tsv.gate.tsv"

# Check forward_backward_agreement is no longer all-NA.
FB_NA=$(python3 -c "
import csv, sys
rows = list(csv.DictReader(open('$LOCAL_RUN/validated.tsv'), delimiter='\t'))
na = sum(1 for r in rows if r.get('forward_backward_agreement','NA') in ('NA','','NaN'))
print(na)
")
FB_TOTAL=$(python3 -c "
import csv
print(sum(1 for _ in csv.DictReader(open('$LOCAL_RUN/validated.tsv'), delimiter='\t')))
")
echo "forward_backward_agreement: NA=$FB_NA / $FB_TOTAL"
if (( FB_NA == FB_TOTAL )); then
    echo "WARNING: forward_backward_agreement is still all-NA — backward views may not have activated"
elif (( FB_NA > 0 )); then
    echo "NOTE: $FB_NA rows still have NA forward_backward_agreement (some scans may lack backward patches)"
fi

echo "FETCHED ARTIFACT PASS: job=$JOB_ID run=$RUN_ID combined=$COMBINED_SHA256"

run_remote "test \"\$(cat $(printf '%q' "$REMOTE_RUN/.stmfit-bwd-owned"))\" = $(printf '%q' "$RUN_ID") && rm -r -- $(printf '%q' "$REMOTE_RUN")"
echo "remote_owned_run_removed=$REMOTE_RUN"
