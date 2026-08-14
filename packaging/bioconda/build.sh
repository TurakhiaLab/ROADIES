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
    # libpll's autotools build has a bison/flex header-generation race under
    # high parallelism (intermittent "parse_utree.h: No such file or
    # directory"); a bare retry succeeds, so do that instead of -j1.
    make -j"${CPU_COUNT}" || make -j"${CPU_COUNT}"
    make install
    popd
fi

# --- Build MLIPPER (vendored source, links against the libpll built above) ---
# MLIPPER/MLIPPER is a prebuilt binary checked into git (see .gitignore/
# README.md) - rebuild it here rather than trust that copy, since a conda
# package needs a binary built for the target platform, not one dev
# machine's glibc/CUDA/GPU-arch.
# CUDA_HOME/PLL_INC_DIR/PLL_LIB_DIR are `?=` in MLIPPER/Makefile (overridable
# without touching it). conda's cuda-nvcc/cuda-cudart-dev put nvcc and its
# headers/libs under $BUILD_PREFIX/targets/x86_64-linux/, not the plain
# include/lib64 paths the Makefile defaults to - point CUDA_HOME there so
# g++ (which compiles the .cpp files, unlike nvcc, has no fallback search
# path) can find cuda_runtime.h; NVCC stays pinned to its standard location.
pushd MLIPPER
make clean || true
make CUDA_HOME="${BUILD_PREFIX}/targets/x86_64-linux" \
     NVCC="${BUILD_PREFIX}/bin/nvcc" \
     PLL_INC_DIR="${PREFIX}/include" \
     PLL_LIB_DIR="${PREFIX}/lib" \
     -j"${CPU_COUNT}"
popd

# --- Vendor-build TWILIGHT from source (mirrors roadies_env.sh's local dev
# setup). Not a conda run: dependency - it's invoked via the hardcoded path
# ${roadies_root}/TWILIGHT/bin/twilight (workflow/scripts/placement.sh), and
# its build vendors its own oneTBB from source, sidestepping the tbb-version
# conflict that makes `twilight`+`kegalign` unsolvable as plain run: deps
# (see meta.yaml).
if [[ ! -f "TWILIGHT/bin/twilight" ]]; then
    if [[ ! -d "TWILIGHT" ]]; then
        git clone --depth 1 https://github.com/TurakhiaLab/TWILIGHT.git
    fi
    # TWILIGHT/CMakeLists.txt hardcodes -march=native, which both breaks
    # nvcc's parsing of newer GCC AVX512BF16/AMX intrinsics headers on
    # capable build hosts and bakes in that host's exact CPU features -
    # wrong for a redistributable package regardless. Portable baseline
    # instead (AVX2/FMA/BMI2, safe since Haswell/2013).
    sed -i 's/-march=native/-march=x86-64-v3/g' TWILIGHT/CMakeLists.txt
    pushd TWILIGHT
    # buildTWILIGHT.sh's cmake step relies on CMake's own check_language(CUDA)
    # auto-detection, which doesn't find nvcc under $BUILD_PREFIX on its own
    # (only an explicit -DCMAKE_CUDA_COMPILER works) - and the script has no
    # set -e, so a failed cmake/make here doesn't stop it or fail the build.
    # Inject the flag via a `cmake` wrapper on PATH instead of patching
    # TWILIGHT's source, since buildTWILIGHT.sh builds its own cmake
    # invocation internally with no env-var hook for extra flags.
    CMAKE_WRAP_DIR="$(mktemp -d)"
    REAL_CMAKE="$(command -v cmake)"
    cat > "${CMAKE_WRAP_DIR}/cmake" <<EOF
#!/bin/bash
if [[ "\$1" == --* ]]; then
    # cmake --build/--install (tool mode) must see that flag first - a
    # prepended -D here would make cmake misparse it as a configure arg.
    exec "${REAL_CMAKE}" "\$@"
else
    # -DCMAKE_CUDA_COMPILER: see above.
    # -DCMAKE_POLICY_VERSION_MINIMUM=3.5: TWILIGHT's vendored oneTBB 2021.9.0
    #   declares cmake_minimum_required(VERSION <3.5), which modern cmake
    #   (>=4) refuses outright; this is cmake's own documented fix.
    # -DTBB_TEST=OFF: skip oneTBB's test suite (not needed - only the
    #   library + TBBConfig.cmake matter here), which otherwise hits a GCC
    #   12.4 false-positive -Wstringop-overflow treated as fatal.
    exec "${REAL_CMAKE}" -DCMAKE_CUDA_COMPILER="${BUILD_PREFIX}/bin/nvcc" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DTBB_TEST=OFF "\$@"
fi
EOF
    chmod +x "${CMAKE_WRAP_DIR}/cmake"
    # buildTWILIGHT.sh's find_tbb_dev() shells out to `dpkg -l | grep
    # libtbb-dev` to decide whether to use the host's system TBB instead of
    # vendoring its own - not sandboxable by conda-build, and host-dependent
    # either way. Force the vendored-oneTBB path deterministically.
    cat > "${CMAKE_WRAP_DIR}/dpkg" <<'EOF'
#!/bin/bash
exit 1
EOF
    chmod +x "${CMAKE_WRAP_DIR}/dpkg"
    PATH="${CMAKE_WRAP_DIR}:${PATH}" bash install/buildTWILIGHT.sh cuda
    rm -rf "${CMAKE_WRAP_DIR}"
    popd
fi
# buildTWILIGHT.sh swallows its own failures (see above) - check for real
# rather than trusting its exit code.
test -f "TWILIGHT/bin/twilight" || { echo "ERROR: TWILIGHT build failed - TWILIGHT/bin/twilight not produced" >&2; exit 1; }

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

# Copy the entire ROADIES directory to the PREFIX directory
cp -rf * ${PREFIX}/ROADIES

# libpll-src is a build-time-only vendored dependency; don't ship its source
# tree inside the installed package. Removed from the PREFIX copy, not from
# $SRC_DIR itself - conda-build's own post-build provenance step still needs
# $SRC_DIR/libpll-src to exist.
rm -rf "${PREFIX}/ROADIES/libpll-src"

# Debugging: Verify the contents of the PREFIX directory
echo "Contents of PREFIX:"
ls -al $PREFIX
