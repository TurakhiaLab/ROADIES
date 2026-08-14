#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DEFAULT_JOBS="1"

jobs="$DEFAULT_JOBS"
cuda_home="${CUDA_HOME:-}"
pll_inc_dir="/usr/include"
pll_lib_dir=""
skip_apt=0
skip_libpll=0
skip_mlipper=0

usage() {
  cat <<EOF
Usage:
  $(basename "$0") [options]

Set up MLIPPER directly on the host without Docker by following the ROADIES
build path while using the distro-packaged libpll:
1. install apt dependencies
2. install libpll-dev
3. build MLIPPER with USE_DOUBLE=1

This script does not install the CUDA toolkit for you. A working host CUDA
installation with nvcc is required before building MLIPPER.

Options:
  --jobs N              Parallel build jobs
                        Default: $DEFAULT_JOBS
  --cuda-home PATH      Override CUDA_HOME
  --pll-inc-dir PATH    libpll include directory
                        Default: /usr/include
  --pll-lib-dir PATH    libpll library directory
                        Default: auto-detect under /usr/lib/<multiarch>
  --skip-apt            Skip apt dependency installation
  --skip-libpll         Skip apt install of libpll-dev
  --skip-mlipper        Skip MLIPPER build
  --skip-apt-build      Backward-compatible alias of --skip-apt
  --skip-libpll-build   Backward-compatible alias of --skip-libpll
  --skip-mlipper-build  Backward-compatible alias of --skip-mlipper
  -h, --help            Show this message

Examples:
  $(basename "$0")
  $(basename "$0") --jobs 8
  $(basename "$0") --skip-apt --pll-lib-dir /usr/lib/x86_64-linux-gnu
EOF
}

die() {
  echo "$(basename "$0"): $*" >&2
  exit 1
}

run_root() {
  if [[ "$EUID" -eq 0 ]]; then
    "$@"
  elif command -v sudo >/dev/null 2>&1; then
    sudo "$@"
  else
    die "need root privileges for: $* (install sudo or run as root)"
  fi
}

detect_cuda_home() {
  if [[ -n "$cuda_home" ]]; then
    :
  elif [[ -n "${CUDA_HOME:-}" ]]; then
    cuda_home="$CUDA_HOME"
  elif command -v nvcc >/dev/null 2>&1; then
    cuda_home="$(cd "$(dirname "$(command -v nvcc)")/.." && pwd)"
  elif [[ -x /usr/local/cuda/bin/nvcc ]]; then
    cuda_home="/usr/local/cuda"
  elif [[ -x /usr/local/cuda-12/bin/nvcc ]]; then
    cuda_home="/usr/local/cuda-12"
  else
    die "nvcc not found; install CUDA first or pass --cuda-home"
  fi

  [[ -x "$cuda_home/bin/nvcc" ]] || die "nvcc not executable under CUDA_HOME=$cuda_home"
}

detect_pll_paths() {
  if [[ -z "$pll_lib_dir" ]]; then
    local multiarch=""
    if command -v gcc >/dev/null 2>&1; then
      multiarch="$(gcc -print-multiarch 2>/dev/null || true)"
    fi
    if [[ -z "$multiarch" ]] && command -v dpkg-architecture >/dev/null 2>&1; then
      multiarch="$(dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null || true)"
    fi
    if [[ -n "$multiarch" && -d "/usr/lib/$multiarch" ]]; then
      pll_lib_dir="/usr/lib/$multiarch"
    else
      pll_lib_dir="/usr/lib/x86_64-linux-gnu"
    fi
  fi
}

apt_updated=0

apt_update_once() {
  command -v apt-get >/dev/null 2>&1 || die "apt-get not found; this script currently supports apt-based systems only"
  if [[ "$apt_updated" -eq 0 ]]; then
    run_root apt-get update
    apt_updated=1
  fi
}

install_apt_deps() {
  apt_update_once
  run_root apt-get install -y --no-install-recommends \
    autoconf \
    build-essential \
    ca-certificates \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libtbb-dev \
    python3
}

install_libpll() {
  apt_update_once
  run_root apt-get install -y --no-install-recommends libpll-dev
}

build_mlipper() {
  pushd "$REPO_ROOT" >/dev/null
  make clean
  make -j"$jobs" \
    USE_DOUBLE="1" \
    CUDA_HOME="$cuda_home" \
    PLL_INC_DIR="$pll_inc_dir" \
    PLL_LIB_DIR="$pll_lib_dir" \
    MLIPPER
  popd >/dev/null
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)
      jobs="${2:-}"
      shift 2
      ;;
    --cuda-home)
      cuda_home="${2:-}"
      shift 2
      ;;
    --pll-inc-dir)
      pll_inc_dir="${2:-}"
      shift 2
      ;;
    --pll-lib-dir)
      pll_lib_dir="${2:-}"
      shift 2
      ;;
    --skip-apt|--skip-apt-build)
      skip_apt=1
      shift
      ;;
    --skip-libpll|--skip-libpll-build)
      skip_libpll=1
      shift
      ;;
    --skip-mlipper|--skip-mlipper-build)
      skip_mlipper=1
      shift
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

[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
[[ -n "$pll_inc_dir" ]] || die "--pll-inc-dir must not be empty"

if [[ "$skip_mlipper" -eq 0 ]]; then
  detect_cuda_home
fi
if [[ "$skip_mlipper" -eq 0 || "$skip_libpll" -eq 0 ]]; then
  detect_pll_paths
fi

echo "REPO_ROOT=$REPO_ROOT"
echo "CUDA_HOME=$cuda_home"
echo "USE_DOUBLE=1"
echo "JOBS=$jobs"
echo "PLL_INC_DIR=$pll_inc_dir"
echo "PLL_LIB_DIR=$pll_lib_dir"

if [[ "$skip_apt" -eq 0 ]]; then
  install_apt_deps
fi

if [[ "$skip_libpll" -eq 0 ]]; then
  install_libpll
fi

if [[ "$skip_mlipper" -eq 0 ]]; then
  if [[ ! -f "$pll_inc_dir/libpll/pll.h" ]]; then
    die "missing libpll header: $pll_inc_dir/libpll/pll.h"
  fi
  if [[ ! -e "$pll_lib_dir/libpll.so" && ! -e "$pll_lib_dir/libpll.a" && ! -e "$pll_lib_dir/libpll.so.0" ]]; then
    die "could not find libpll library under: $pll_lib_dir"
  fi

  build_mlipper
fi

echo
echo "Done."
echo "Binary: $REPO_ROOT/MLIPPER"
echo "Example:"
echo "  $REPO_ROOT/MLIPPER --help"
