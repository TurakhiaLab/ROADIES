#!/bin/bash
set -euo pipefail

# Debugging: Print current directory and list its contents
echo "Current directory: $(pwd)"
echo "Contents:"
ls -al

mkdir -p $PREFIX/ROADIES

# --- Build libpll (vendored source, no conda package exists) ---
# Installed into $PREFIX so MLIPPER's build below can link against it via
# PLL_INC_DIR/PLL_LIB_DIR, and so the resulting MLIPPER binary's rpath
# ($PREFIX/lib, set by conda's compiler wrappers) resolves libpll.so at runtime.
if [[ ! -f "$PREFIX/lib/libpll.so" ]]; then
    pushd libpll-src
    ./autogen.sh
    ./configure --prefix="$PREFIX"
    make -j"${CPU_COUNT}"
    make install
    popd
fi

# --- Build MLIPPER (vendored source, links against the libpll built above) ---
# IMPORTANT: MLIPPER/MLIPPER is currently a *prebuilt binary checked into git*
# (confirmed via `git ls-files MLIPPER/MLIPPER` - it's tracked, not gitignored).
# That's fine for local dev but wrong for a conda package: it'd ship a binary
# built on one dev machine's glibc/CUDA/GPU-arch combo instead of something
# built reproducibly for the target platform. Rebuild it unconditionally here
# rather than trusting the checked-in one - separately, consider gitignoring
# MLIPPER/MLIPPER upstream so `git status` doesn't show it as dirty after every
# local build either.
# CUDA_HOME/PLL_INC_DIR/PLL_LIB_DIR are all `?=` in MLIPPER/Makefile, i.e.
# overridable - this avoids touching MLIPPER's own Makefile.
# NOTE: verify $BUILD_PREFIX is where cuda-nvcc/cuda-cudart-dev actually place
# nvcc + cudart headers/libs for the pinned cuda-version - this may need
# adjusting to match whatever conda-forge's current cuda packaging layout is
# (check TWILIGHT's build.sh / cuda-nvcc feedstock for the current convention).
pushd MLIPPER
make clean || true
make CUDA_HOME="${BUILD_PREFIX}" \
     PLL_INC_DIR="${PREFIX}/include" \
     PLL_LIB_DIR="${PREFIX}/lib" \
     -j"${CPU_COUNT}"
popd

# Build sampling code
if [[ ! -d "workflow/scripts/sampling/build" ]]; then
    cd workflow/scripts/sampling
    mkdir -p build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX="${PREFIX}"
    make -j"${CPU_COUNT}"
    cd ../../../..
fi

# Debugging: Print current directory and list its contents before copying
echo "Current directory before copying ROADIES: $(pwd)"
echo "Contents before copying ROADIES:"
ls -al

# libpll-src is a build-time-only vendored dependency; don't ship its source
# tree inside the installed ROADIES package.
rm -rf libpll-src

# Copy the entire ROADIES directory to the PREFIX directory
cp -rf * ${PREFIX}/ROADIES

# Debugging: Verify the contents of the PREFIX directory
echo "Contents of PREFIX:"
ls -al $PREFIX
