#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_FILE="$SCRIPT_DIR/remote.env"

usage() {
    cat <<'EOF'
Usage: launch_hierarchical_unit_assignment.sh --dry-run|--submit|--watch

Exactly one mode is mandatory:
  --dry-run  Resolve and print hashes, Viper paths, array shape, merge, and fetch.
  --submit   Sync, instantiate on the login node, and submit the array.
  --watch    Sync, submit, wait, merge, fetch, and locally validate artifacts.

Environment:
  STMFIT_HIERARCHICAL_FEATURES   Local label-free feature TSV (required).
  STMFIT_HIERARCHICAL_CONFIG     Candidate config [config/unit_assignment_candidate.toml].
  STMFIT_HIERARCHICAL_LOCAL_OUT  Fetch base [results/hierarchical_unit_assignment].
EOF
}

(( $# == 1 )) || { usage >&2; exit 2; }
case "$1" in
    --dry-run) MODE=dry-run ;;
    --submit) MODE=submit ;;
    --watch) MODE=watch ;;
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
: "${STMFIT_HIERARCHICAL_FEATURES:?STMFIT_HIERARCHICAL_FEATURES must name the local label-free feature TSV}"
: "${STMFIT_HIERARCHICAL_CONFIG:=config/unit_assignment_candidate.toml}"
: "${STMFIT_HIERARCHICAL_LOCAL_OUT:=results/hierarchical_unit_assignment}"
: "${JULIA_MODULE_VERSION:=1.12}"
: "${WATCH_POLL:=30}"
: "${SSH_CONNECT_TIMEOUT:=180}"
: "${SSH_SERVER_ALIVE_INTERVAL:=60}"
: "${RSYNC_EXTRA:=}"

SSH_HOST="${STMFIT_HIERARCHICAL_SSH_HOST:-viper}"
[[ "$SSH_HOST" == "viper" ]] || {
    echo "ERROR: the fixed T5 heavy gate must target the Viper SSH alias" >&2
    exit 1
}
[[ -f "$STMFIT_HIERARCHICAL_FEATURES" ]] || {
    echo "ERROR: feature TSV not found: $STMFIT_HIERARCHICAL_FEATURES" >&2
    exit 1
}
CONFIG_PATH="$STMFIT_HIERARCHICAL_CONFIG"
[[ "$CONFIG_PATH" = /* ]] || CONFIG_PATH="$REPO_ROOT/$CONFIG_PATH"
[[ -f "$CONFIG_PATH" ]] || { echo "ERROR: config not found: $CONFIG_PATH" >&2; exit 1; }

sha_file() { sha256sum "$1" | cut -d' ' -f1; }
CONFIG_SHA256=$(sha_file "$CONFIG_PATH")
INPUT_SHA256=$(sha_file "$STMFIT_HIERARCHICAL_FEATURES")
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
} > "$SOURCE_MANIFEST"
SOURCE_SHA256=$(sha_file "$SOURCE_MANIFEST")

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
RUN_ID="${INPUT_SHA256:0:12}-${CONFIG_SHA256:0:12}-${SOURCE_SHA256:0:12}"

REMOTE_PROJECT_BASE="${REMOTE_PROJECT_BASE:-/u/$STMFIT_REMOTE_USER/code}"
REMOTE_PROJECT="$REMOTE_PROJECT_BASE/STMFit"
REMOTE_RUN="/ptmp/$STMFIT_REMOTE_USER/stmfit/hierarchical_unit_assignment/$RUN_ID"
REMOTE_FEATURES="$REMOTE_RUN/inputs/features.tsv"
REMOTE_OUT="$REMOTE_RUN/output"
LOCAL_RUN="$REPO_ROOT/$STMFIT_HIERARCHICAL_LOCAL_OUT/$RUN_ID"
REMOTE_CONFIG_REL="config/unit_assignment_candidate.toml"

SSH_OPTS=(-o BatchMode=no -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
          -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL")
run_remote() { ssh "${SSH_OPTS[@]}" "$SSH_HOST" "$1"; }

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
MERGE_COMMAND="cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && julia --project=. test/merge_hierarchical_unit_assignment_shards.jl${SHARD_ARGS_Q} --out $(printf '%q' "$MERGED_REMOTE") --expected-folds $(printf '%q' "$FOLDS_CSV") --config-sha256 $CONFIG_SHA256 --input-sha256 $INPUT_SHA256 --source-sha256 $SOURCE_SHA256"

EXPORTS="ALL,STMFIT_PROJECT_DIR=$REMOTE_PROJECT,STMFIT_HIERARCHICAL_FEATURES=$REMOTE_FEATURES,STMFIT_HIERARCHICAL_CONFIG=$REMOTE_CONFIG_REL,STMFIT_HIERARCHICAL_OUTBASE=$REMOTE_OUT,STMFIT_HIERARCHICAL_FOLDS=$FOLDS_ENV,STMFIT_HIERARCHICAL_CONFIG_SHA256=$CONFIG_SHA256,STMFIT_HIERARCHICAL_INPUT_SHA256=$INPUT_SHA256,STMFIT_HIERARCHICAL_SOURCE_SHA256=$SOURCE_SHA256,JULIA_MODULE_VERSION=$JULIA_MODULE_VERSION"
SBATCH_COMMAND="cd $(printf '%q' "$REMOTE_PROJECT") && sbatch --parsable --array=1-${ARRAY_TASKS}%8 --cpus-per-task=1 --export=$(printf '%q' "$EXPORTS") hpc/hierarchical_unit_assignment_array.sbatch"

cat <<EOF
T5 hierarchical leave-date-out/bootstrap gate
  mode:                $MODE
  host:                $SSH_HOST
  remote project:      $REMOTE_PROJECT
  remote owned run:    $REMOTE_RUN
  remote features:     $REMOTE_FEATURES
  remote output:       $REMOTE_OUT
  local fetch:         $LOCAL_RUN
  folds:               $FOLDS_CSV ($N_FOLDS)
  bootstrap blocks:    unbootstrapped + 0:49,...,450:499
  array:               1-${ARRAY_TASKS}%8
  CPU/thread per task: 1 / 1
  config sha256:       $CONFIG_SHA256
  input sha256:        $INPUT_SHA256
  source sha256:       $SOURCE_SHA256
  submit command:      $SBATCH_COMMAND
  merge command:       $MERGE_COMMAND
  fetch source:        $SSH_HOST:$REMOTE_OUT/
  fetch destination:   $LOCAL_RUN/
EOF

if [[ "$MODE" == dry-run ]]; then
    echo "DRY RUN PASS: no remote or local state changed"
    exit 0
fi

run_remote "mkdir -p $(printf '%q' "$REMOTE_PROJECT") $(printf '%q' "$REMOTE_RUN/inputs") $(printf '%q' "$REMOTE_OUT/shards") && printf '%s\\n' $(printf '%q' "$RUN_ID") > $(printf '%q' "$REMOTE_RUN/.stmfit-t5-owned")"
rsync -az --exclude='.git/' --exclude='Manifest.toml' --exclude='results/' \
    --exclude='hpc/remote.env' ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$REPO_ROOT/" "$SSH_HOST:$REMOTE_PROJECT/"
rsync -az ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$STMFIT_HIERARCHICAL_FEATURES" "$SSH_HOST:$REMOTE_FEATURES"

run_remote "cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && julia --project=. -e 'using Pkg; for p in [\"STMFitCore\",\"STMSXMIO\",\"GaussianFit1D\",\"GaussianFit2D\",\"STMMolecularFit\"]; Pkg.develop(PackageSpec(path=\"packages/\" * p * \".jl\")); end; Pkg.instantiate(); Pkg.precompile()'"

JOB_ID=$(run_remote "$SBATCH_COMMAND")
JOB_ID=${JOB_ID%%;*}
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || { echo "ERROR: could not parse Slurm job id: $JOB_ID" >&2; exit 1; }
echo "array_job_id=$JOB_ID"

if [[ "$MODE" == submit ]]; then
    echo "SUBMITTED: rerun --watch only if a new complete run is intended; remote shards are independently resumable"
    exit 0
fi

while :; do
    ACTIVE=$(run_remote "squeue -h -j $(printf '%q' "$JOB_ID") 2>/dev/null | wc -l")
    [[ "$ACTIVE" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid squeue response: $ACTIVE" >&2; exit 1; }
    (( ACTIVE == 0 )) && break
    printf 'waiting: job=%s active=%s\n' "$JOB_ID" "$ACTIVE"
    sleep "$WATCH_POLL"
done
SACCT_OUTPUT=$(run_remote "sacct -j $(printf '%q' "$JOB_ID") --format=JobIDRaw,State,ExitCode -n -P")
printf '%s\n' "$SACCT_OUTPUT"

array_total=0
array_failed=0
plain_total=0
plain_failed=0
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
    accounting_total=$array_total
    accounting_failed=$array_failed
else
    # Viper's submit filter materializes array tasks as separate numeric jobs.
    accounting_total=$plain_total
    accounting_failed=$plain_failed
fi
(( accounting_total == ARRAY_TASKS )) || {
    echo "ERROR: expected $ARRAY_TASKS top-level accounting rows, found $accounting_total" >&2
    exit 1
}
(( accounting_failed == 0 )) || {
    echo "ERROR: $accounting_failed of $accounting_total array tasks did not complete with exit 0:0" >&2
    exit 1
}
echo "accounting_pass=tasks:$accounting_total"
run_remote "$MERGE_COMMAND"

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
    --input-sha256 "$INPUT_SHA256" --source-sha256 "$SOURCE_SHA256"

cmp "$LOCAL_RUN/merged.tsv" "$LOCAL_RUN/validated.tsv"
cmp "$LOCAL_RUN/merged.tsv.gate.tsv" "$LOCAL_RUN/validated.tsv.gate.tsv"
grep -q '^PASS' "$LOCAL_RUN/validated.tsv.gate.tsv"
echo "FETCHED ARTIFACT PASS: job=$JOB_ID run=$RUN_ID"

run_remote "test \"\$(cat $(printf '%q' "$REMOTE_RUN/.stmfit-t5-owned"))\" = $(printf '%q' "$RUN_ID") && rm -r -- $(printf '%q' "$REMOTE_RUN")"
echo "remote_owned_run_removed=$REMOTE_RUN"
