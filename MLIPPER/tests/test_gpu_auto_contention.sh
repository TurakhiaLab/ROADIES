#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  tests/test_gpu_auto_contention.sh [options]

Launch multiple concurrent host-side MLIPPER subprocesses with --gpu-auto to
exercise GPU contention and reservation behavior.

Options:
  --workers N              Number of concurrent subprocesses to launch. Default: 2
  --start-delay SEC        Delay before all workers launch together. Default: 3
  --timeout SEC            Per-worker timeout passed to `timeout`. Default: 0 (disabled)
  --gene-dir PATH          Gene input directory. Must contain iter0 ref/query MSA,
                           bestTree, and bestModel. Default:
                           MLIPPER/data/mlipper_placement_required_inputs_v2_split_bundle/gene_638
  --tree-alignment PATH    Explicit tree-alignment path override
  --query-alignment PATH   Explicit query-alignment path override
  --tree PATH              Explicit reference tree path override
  --best-model PATH        Explicit bestModel path override
  --mlipper-bin PATH       MLIPPER executable path. Default: MLIPPER/MLIPPER
  --outdir PATH            Output directory for worker logs and trees. Default:
                           MLIPPER/output/gpu_auto_contention_<timestamp>
  --lock-dir PATH          Shared lock directory used by all test workers. Default:
                           /tmp/mlipper_gpu_contention_<user>_<pid>
  --batch-insert-size N    Forwarded to MLIPPER. Default: 5
  --local-spr-radius N     Forwarded to MLIPPER. Default: 4
  --local-spr-rounds N     Forwarded to MLIPPER. Default: 3
  -h, --help               Show this help

Notes:
  - All workers in one test invocation share the same MLIPPER_GPU_LOCK_DIR.
  - To force visible contention, set --workers greater than the number of visible GPUs.
  - The script invokes the local MLIPPER binary directly and does not use Docker.
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

require_file() {
  local path="$1"
  [[ -e "$path" ]] || die "missing required file: $path"
}

timestamp() {
  date '+%Y-%m-%dT%H:%M:%S%z'
}

summary_gpu_line() {
  local log_path="$1"
  grep -m1 -E 'Auto-reserved CUDA device|Using reserved CUDA device|Using CUDA device|waiting for a free GPU' "$log_path" || true
}

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_GENE_DIR="$REPO_ROOT/data/mlipper_placement_required_inputs_v2_split_bundle/gene_638"
DEFAULT_MLIPPER_BIN="$REPO_ROOT/MLIPPER"
DEFAULT_OUTDIR="$REPO_ROOT/output/gpu_auto_contention_$(date +%Y%m%d_%H%M%S)"

workers=2
start_delay=3
timeout_seconds=0
gene_dir="$DEFAULT_GENE_DIR"
mlipper_bin="$DEFAULT_MLIPPER_BIN"
outdir="$DEFAULT_OUTDIR"
lock_dir=""
tree_alignment_path=""
query_alignment_path=""
tree_path=""
best_model_path=""
batch_insert_size=5
local_spr_radius=4
local_spr_rounds=3

while (($# > 0)); do
  case "$1" in
    --workers)
      workers="$2"
      shift 2
      ;;
    --start-delay)
      start_delay="$2"
      shift 2
      ;;
    --timeout)
      timeout_seconds="$2"
      shift 2
      ;;
    --gene-dir)
      gene_dir="$2"
      shift 2
      ;;
    --tree-alignment)
      tree_alignment_path="$2"
      shift 2
      ;;
    --query-alignment)
      query_alignment_path="$2"
      shift 2
      ;;
    --tree)
      tree_path="$2"
      shift 2
      ;;
    --best-model)
      best_model_path="$2"
      shift 2
      ;;
    --mlipper-bin)
      mlipper_bin="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    --lock-dir)
      lock_dir="$2"
      shift 2
      ;;
    --batch-insert-size)
      batch_insert_size="$2"
      shift 2
      ;;
    --local-spr-radius)
      local_spr_radius="$2"
      shift 2
      ;;
    --local-spr-rounds)
      local_spr_rounds="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[[ "$workers" =~ ^[0-9]+$ ]] || die "--workers must be a non-negative integer"
[[ "$start_delay" =~ ^[0-9]+$ ]] || die "--start-delay must be a non-negative integer"
[[ "$timeout_seconds" =~ ^[0-9]+$ ]] || die "--timeout must be a non-negative integer"
[[ "$batch_insert_size" =~ ^[0-9]+$ ]] || die "--batch-insert-size must be a non-negative integer"
[[ "$local_spr_radius" =~ ^[0-9]+$ ]] || die "--local-spr-radius must be a non-negative integer"
[[ "$local_spr_rounds" =~ ^[0-9]+$ ]] || die "--local-spr-rounds must be a non-negative integer"
(( workers > 0 )) || die "--workers must be > 0"
(( local_spr_rounds > 0 )) || die "--local-spr-rounds must be > 0"

mlipper_bin="$(cd "$(dirname "$mlipper_bin")" && pwd)/$(basename "$mlipper_bin")"
mkdir -p "$outdir"
outdir="$(cd "$outdir" && pwd)"

if [[ -z "$lock_dir" ]]; then
  lock_dir="/tmp/mlipper_gpu_contention_${USER:-user}_$$"
fi
mkdir -p "$lock_dir"
lock_dir="$(cd "$lock_dir" && pwd)"

[[ -x "$mlipper_bin" ]] || die "MLIPPER binary is not executable: $mlipper_bin"

if [[ -n "$tree_alignment_path" || -n "$query_alignment_path" || -n "$tree_path" || -n "$best_model_path" ]]; then
  [[ -n "$tree_alignment_path" ]] || die "--tree-alignment is required when using explicit input overrides"
  [[ -n "$query_alignment_path" ]] || die "--query-alignment is required when using explicit input overrides"
  [[ -n "$tree_path" ]] || die "--tree is required when using explicit input overrides"
  [[ -n "$best_model_path" ]] || die "--best-model is required when using explicit input overrides"
  gene_dir="<explicit-inputs>"
  ref_msa="$tree_alignment_path"
  query_msa="$query_alignment_path"
  best_tree="$tree_path"
  best_model="$best_model_path"
else
  gene_dir="$(cd "$gene_dir" && pwd)"
  ref_msa="$gene_dir/iter0_output_msa_from_ref.fa"
  query_msa="$gene_dir/iter0_output_msa_from_query.fa"
  best_tree="$gene_dir/$(basename "$gene_dir")_filtered.fa.aln.raxml.bestTree"
  best_model="$gene_dir/$(basename "$gene_dir")_filtered.fa.aln.raxml.bestModel"
fi

require_file "$ref_msa"
require_file "$query_msa"
require_file "$best_tree"
require_file "$best_model"

start_epoch=$(( $(date +%s) + start_delay ))

echo "GPU auto contention test"
echo "repo_root=$REPO_ROOT"
echo "mlipper_bin=$mlipper_bin"
echo "gene_dir=$gene_dir"
echo "outdir=$outdir"
echo "lock_dir=$lock_dir"
echo "workers=$workers"
echo "start_delay=$start_delay"
echo "timeout_seconds=$timeout_seconds"
echo "start_epoch=$start_epoch"
echo "batch_insert_size=$batch_insert_size"
echo "local_spr_radius=$local_spr_radius"
echo "local_spr_rounds=$local_spr_rounds"

declare -a pids=()
declare -a worker_dirs=()
declare -a logs=()
declare -a exit_codes=()

cleanup_children() {
  local pid
  for pid in "${pids[@]:-}"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      kill "$pid" 2>/dev/null || true
    fi
  done
}

trap cleanup_children INT TERM

for worker_id in $(seq 1 "$workers"); do
  worker_dir="$outdir/worker_${worker_id}"
  mkdir -p "$worker_dir"
  log_path="$worker_dir/mlipper.log"
  out_tree="$worker_dir/mlipper_gene_tree.nwk"
  worker_dirs+=("$worker_dir")
  logs+=("$log_path")

  cmd=(
    "$mlipper_bin"
    --tree-alignment "$ref_msa"
    --query-alignment "$query_msa"
    --tree "$best_tree"
    --best-model "$best_model"
    --commit-to-tree "$out_tree"
    --no-model-optimization
    --no-global-branch-optimization
    --batch-insert-size "$batch_insert_size"
    --local-spr-radius "$local_spr_radius"
    --local-spr-rounds "$local_spr_rounds"
    --local-spr-dynamic-validation-conflicts
    --gpu-auto
  )

  (
    set -euo pipefail
    {
      echo "worker_id=$worker_id"
      echo "launch_barrier_epoch=$start_epoch"
      echo "lock_dir=$lock_dir"
      echo "launch_ts=$(timestamp)"
      printf 'cmd:'
      printf ' %q' "${cmd[@]}"
      printf '\n'
      while (( $(date +%s) < start_epoch )); do
        sleep 0.05
      done
      echo "exec_ts=$(timestamp)"
      if (( timeout_seconds > 0 )); then
        MLIPPER_GPU_LOCK_DIR="$lock_dir" timeout "$timeout_seconds" "${cmd[@]}"
      else
        MLIPPER_GPU_LOCK_DIR="$lock_dir" "${cmd[@]}"
      fi
    } >"$log_path" 2>&1
  ) &

  current_pid="$!"
  pids+=("$current_pid")
  echo "spawned worker=$worker_id pid=$current_pid log=$log_path"
done

failures=0
for idx in "${!pids[@]}"; do
  pid="${pids[$idx]}"
  if wait "$pid"; then
    status=0
  else
    status=$?
    failures=1
  fi
  exit_codes+=("$status")
  printf '%s\n' "$status" > "${worker_dirs[$idx]}/exit_code.txt"
done

trap - INT TERM

echo
echo "Summary"
for idx in "${!worker_dirs[@]}"; do
  worker_num=$((idx + 1))
  gpu_line="$(summary_gpu_line "${logs[$idx]}")"
  printf 'worker=%d exit=%s\n' "$worker_num" "${exit_codes[$idx]}"
  printf '  log=%s\n' "${logs[$idx]}"
  printf '  tree=%s\n' "${worker_dirs[$idx]}/mlipper_gene_tree.nwk"
  if [[ -n "$gpu_line" ]]; then
    printf '  gpu=%s\n' "$gpu_line"
  fi
done

echo
echo "Lock directory contents:"
find "$lock_dir" -maxdepth 1 -type f | sort || true

if (( failures != 0 )); then
  exit 1
fi
