#!/usr/bin/env bash

# Launch the label-free unit-assignment prediction sweep on Viper.
#   --dry-run  resolve and print paths, no remote change
#   --watch    sync -> instantiate -> submit -> wait -> fetch
#
# Environment:
#   STMFIT_SWEEP_FEATURES_UASYM  local feature TSV incl. patch_u_asym (required)
#   STMFIT_SWEEP_FEATURES_SPLIT  local feature TSV incl. split_log_skew (required)
#   STMFIT_SWEEP_LOCAL_OUT       fetch base [results/prediction_sweep]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_FILE="$SCRIPT_DIR/remote.env"

usage() {
    cat <<'EOF'
Usage: launch_prediction_sweep.sh --dry-run|--watch

  --dry-run  Resolve and print hashes, Viper paths, and the submit command.
  --watch    Sync, instantiate, submit, wait, fetch.
EOF
}

(( $# == 1 )) || { usage >&2; exit 2; }
case "$1" in
    --dry-run) MODE=dry-run ;;
    --watch) MODE=watch ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: expected exactly one of --dry-run|--watch" >&2; usage >&2; exit 2 ;;
esac

[[ -f "$ENV_FILE" ]] || { echo "ERROR: $ENV_FILE not found" >&2; exit 1; }
# shellcheck source=/dev/null
source "$ENV_FILE"

: "${STMFIT_REMOTE_USER:?STMFIT_REMOTE_USER is required in hpc/remote.env}"
: "${STMFIT_SWEEP_FEATURES_UASYM:?STMFIT_SWEEP_FEATURES_UASYM must name the uasym feature TSV}"
: "${STMFIT_SWEEP_FEATURES_SPLIT:?STMFIT_SWEEP_FEATURES_SPLIT must name the split feature TSV}"
: "${STMFIT_SWEEP_LOCAL_OUT:=results/prediction_sweep}"
: "${JULIA_MODULE_VERSION:=1.12}"
: "${WATCH_POLL:=30}"
: "${SSH_CONNECT_TIMEOUT:=180}"
: "${SSH_SERVER_ALIVE_INTERVAL:=60}"

SSH_HOST="${STMFIT_SSH_HOST:-viper}"
[[ -f "$STMFIT_SWEEP_FEATURES_UASYM" ]] || { echo "ERROR: uasym feature TSV not found" >&2; exit 1; }
[[ -f "$STMFIT_SWEEP_FEATURES_SPLIT" ]] || { echo "ERROR: split feature TSV not found" >&2; exit 1; }

sha_file() { sha256sum "$1" | cut -d' ' -f1; }
UASYM_SHA=$(sha_file "$STMFIT_SWEEP_FEATURES_UASYM")
SPLIT_SHA=$(sha_file "$STMFIT_SWEEP_FEATURES_SPLIT")
SOURCE_MANIFEST=$(mktemp)
trap 'rm -f "$SOURCE_MANIFEST"' EXIT
{
    printf 'test/build_labelfree_gmm_predictions.jl\t%s\n' \
        "$(sha_file "$REPO_ROOT/test/build_labelfree_gmm_predictions.jl")"
    printf 'test/build_labelfree_unit_predictions.jl\t%s\n' \
        "$(sha_file "$REPO_ROOT/test/build_labelfree_unit_predictions.jl")"
} > "$SOURCE_MANIFEST"
SOURCE_SHA=$(sha_file "$SOURCE_MANIFEST")

RUN_ID="${UASYM_SHA:0:12}-${SPLIT_SHA:0:12}-${SOURCE_SHA:0:12}"
REMOTE_PROJECT="/u/$STMFIT_REMOTE_USER/code/STMFit"
REMOTE_RUN="/ptmp/$STMFIT_REMOTE_USER/stmfit/prediction_sweep/$RUN_ID"
REMOTE_UASYM="$REMOTE_RUN/inputs/features_uasym.tsv"
REMOTE_SPLIT="$REMOTE_RUN/inputs/features_split.tsv"
REMOTE_OUT="$REMOTE_RUN/output"
LOCAL_RUN="$REPO_ROOT/$STMFIT_SWEEP_LOCAL_OUT/$RUN_ID"

SSH_OPTS=(-o BatchMode=no -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
          -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL")
run_remote() { ssh "${SSH_OPTS[@]}" "$SSH_HOST" "$1"; }

EXPORTS="ALL,STMFIT_PROJECT_DIR=$REMOTE_PROJECT,STMFIT_FEATURES_UASYM=$REMOTE_UASYM,STMFIT_FEATURES_SPLIT=$REMOTE_SPLIT,STMFIT_OUTPUT_DIR=$REMOTE_OUT,JULIA_MODULE_VERSION=$JULIA_MODULE_VERSION"
SBATCH_COMMAND="cd $(printf '%q' "$REMOTE_PROJECT") && sbatch --parsable --cpus-per-task=2 --export=$(printf '%q' "$EXPORTS") hpc/prediction_sweep.sbatch"

cat <<EOF
Prediction sweep (GMM self-training + seed scaling)
  mode:                $MODE
  host:                $SSH_HOST
  remote project:      $REMOTE_PROJECT
  remote run:          $REMOTE_RUN
  uasym sha256:        $UASYM_SHA
  split sha256:        $SPLIT_SHA
  source sha256:       $SOURCE_SHA
  local fetch:         $LOCAL_RUN
  submit command:      $SBATCH_COMMAND
EOF

if [[ "$MODE" == dry-run ]]; then
    echo "DRY RUN PASS: no remote or local state changed"
    exit 0
fi

run_remote "mkdir -p $(printf '%q' "$REMOTE_PROJECT") $(printf '%q' "$REMOTE_RUN/inputs") $(printf '%q' "$REMOTE_OUT")"
rsync -az --exclude='.git/' --exclude='Manifest.toml' --exclude='results/' \
    --exclude='hpc/remote.env' ${RSYNC_EXTRA:-} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$REPO_ROOT/" "$SSH_HOST:$REMOTE_PROJECT/"
rsync -az ${RSYNC_EXTRA:-} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$STMFIT_SWEEP_FEATURES_UASYM" "$SSH_HOST:$REMOTE_UASYM"
rsync -az ${RSYNC_EXTRA:-} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$STMFIT_SWEEP_FEATURES_SPLIT" "$SSH_HOST:$REMOTE_SPLIT"

run_remote "cd $(printf '%q' "$REMOTE_PROJECT") && module purge && module load julia/$(printf '%q' "$JULIA_MODULE_VERSION") && julia --project=. -e 'using Pkg; for p in [\"STMFitCore\",\"STMSXMIO\",\"GaussianFit1D\",\"GaussianFit2D\",\"STMMolecularFit\"]; Pkg.develop(PackageSpec(path=\"packages/\" * p * \".jl\")); end; Pkg.instantiate(); Pkg.precompile()'"

JOB_ID=$(run_remote "$SBATCH_COMMAND")
JOB_ID=${JOB_ID%%;*}
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || { echo "ERROR: could not parse Slurm job id: $JOB_ID" >&2; exit 1; }
echo "slurm_job_id=$JOB_ID"

while :; do
    ACTIVE=$(run_remote "squeue -h -j $(printf '%q' "$JOB_ID") 2>/dev/null | wc -l")
    [[ "$ACTIVE" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid squeue response" >&2; exit 1; }
    (( ACTIVE == 0 )) && break
    printf 'waiting: job=%s active=%s\n' "$JOB_ID" "$ACTIVE"
    sleep "$WATCH_POLL"
done

SACCT=$(run_remote "sacct -j $(printf '%q' "$JOB_ID") --format=State,ExitCode -n -P")
printf '%s\n' "$SACCT"
grep -q '^COMPLETED|0:0' <<< "$SACCT" || { echo "ERROR: sweep job did not complete cleanly" >&2; exit 1; }

mkdir -p "$LOCAL_RUN"
rsync -az --update ${RSYNC_EXTRA:-} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$SSH_HOST:$REMOTE_OUT/" "$LOCAL_RUN/"
ls "$LOCAL_RUN"/*.tsv
echo "SWEEP FETCHED: job=$JOB_ID run=$RUN_ID -> $LOCAL_RUN"
