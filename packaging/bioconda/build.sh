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
    # libpll's autotools build has a missing dependency edge between its
    # bison-generated headers (parse_utree.h/parse_rtree.h) and the
    # flex-generated sources that #include them (lex_utree.l/lex_rtree.l) -
    # under high parallelism (-j"${CPU_COUNT}") this occasionally races and
    # fails with "parse_utree.h: No such file or directory". VERIFIED
    # (2026-08-06/07, hit twice this session, both in ad hoc manual testing
    # and in a real conda-build run): a bare retry succeeds cleanly every
    # time, since the header already exists in the source tree from the
    # partial first pass - retry once rather than serializing the whole
    # build to -j1 for a race that's rare in practice.
    make -j"${CPU_COUNT}" || make -j"${CPU_COUNT}"
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
# VERIFIED (2026-08-06, GPU host, manual test): with cuda-nvcc+cuda-cudart-dev
# installed, nvcc/cuda_runtime.h/libcudart/libcuda(stub) all actually live
# under $BUILD_PREFIX/targets/x86_64-linux/{bin,include,lib} - NOT under
# $BUILD_PREFIX/include or $BUILD_PREFIX/lib64 (the latter doesn't even
# exist). A standalone nvcc invocation compiles/links fine with the "wrong"
# -I$(CUDA_HOME)/include/-L$(CUDA_HOME)/lib64 flags because nvcc has its own
# built-in default search under <nvcc-dir>/../targets/<arch>/ that fires
# regardless of explicit -I/-L - CONFIRMED this covers MLIPPER's .cu files
# and its final nvcc-driven link step in the real conda-build sandbox below.
# CORRECTED (2026-08-06/07, real conda-build run): that nvcc-only fallback
# does NOT extend to MLIPPER's .cpp files, which are compiled by plain g++
# (not nvcc) and also need cuda_runtime.h (via src/util/precision.hpp) - g++
# has no such fallback, so real conda-build failed here with
# "cuda_runtime.h: No such file or directory" on every .cpp translation
# unit. Fix: point CUDA_HOME itself at the directory that actually has
# bin/include/lib (targets/x86_64-linux) so the Makefile's own
# -I$(CUDA_HOME)/include resolves correctly for g++ too; separately pin NVCC
# back to the standard-location binary ($BUILD_PREFIX/bin/nvcc, which also
# exists) so its useful default-search fallback stays in play for the final
# link (where CUDA_LIB's -L$(CUDA_HOME)/lib64 is still "wrong" - the real
# dir is lib/, not lib64/ - but harmless given that fallback).
pushd MLIPPER
make clean || true
make CUDA_HOME="${BUILD_PREFIX}/targets/x86_64-linux" \
     NVCC="${BUILD_PREFIX}/bin/nvcc" \
     PLL_INC_DIR="${PREFIX}/include" \
     PLL_LIB_DIR="${PREFIX}/lib" \
     -j"${CPU_COUNT}"
popd

# --- Vendor-build TWILIGHT from source (mirrors roadies_env.sh's local dev
# setup: `git clone TurakhiaLab/TWILIGHT && bash install/buildTWILIGHT.sh
# cuda`). TWILIGHT is NOT a conda package dependency here - it's invoked via
# the hardcoded path ${roadies_root}/TWILIGHT/bin/twilight (see
# workflow/scripts/placement.sh), and its own build script vendors its own
# oneTBB 2021.9.0 from source rather than linking a conda `tbb` package.
# That vendoring is exactly what avoids the tbb version conflict that made
# `twilight`+`kegalign` unsolvable as plain run: deps (see meta.yaml).
if [[ ! -f "TWILIGHT/bin/twilight" ]]; then
    if [[ ! -d "TWILIGHT" ]]; then
        git clone --depth 1 https://github.com/TurakhiaLab/TWILIGHT.git
    fi
    # FOUND (2026-08-06/07, real conda-build run): TWILIGHT/CMakeLists.txt
    # hardcodes -march=native for its C++/CUDA-host-compiler flags. Two
    # separate problems: (1) on this build host's CPU (AVX512-BF16/AMX
    # capable), nvcc's host-compiler pass over gcc 12.4's AVX512BF16/AMX
    # intrinsics headers fails outright ("identifier
    # __builtin_ia32_cvtne2ps2bf16_... is undefined"), a known class of
    # nvcc/newest-glibc-header incompatibility; (2) independent of that
    # error, -march=native is simply wrong for a redistributable conda
    # package regardless of which build host it's compiled on - it bakes in
    # *that specific machine's* exact CPU feature set, so the resulting
    # binary would SIGILL on any end-user CPU lacking the same features.
    # Patch to a portable, well-defined microarchitecture baseline
    # (x86-64-v3: AVX2/FMA/BMI2, ubiquitous on anything from Haswell/2013
    # onward) that also sidesteps the AVX512BF16/AMX headers entirely.
    sed -i 's/-march=native/-march=x86-64-v3/g' TWILIGHT/CMakeLists.txt
    pushd TWILIGHT
    # FOUND (2026-08-06/07, real conda-build run): unlike MLIPPER's own
    # Makefile, buildTWILIGHT.sh's cmake step relies entirely on CMake's own
    # `check_language(CUDA)` compiler auto-detection - it never learned
    # about $BUILD_PREFIX's layout, so it silently fell back to "No CUDA
    # found" and produced no TWILIGHT/bin/twilight at all, with the failure
    # never surfacing: buildTWILIGHT.sh has no `set -e` (no shebang either -
    # it's invoked via `bash install/buildTWILIGHT.sh`), so its internal
    # cmake/make failures don't stop the script, and it exits 0 regardless.
    # CORRECTED further: neither PATH (nvcc does resolve correctly via PATH
    # in this sandbox, confirmed with `command -v nvcc`) nor the documented
    # CUDACXX env var actually make check_language(CUDA) succeed here - only
    # an explicit `-DCMAKE_CUDA_COMPILER=...` flag on the cmake command line
    # works (verified interactively against this exact cmake build). Since
    # buildTWILIGHT.sh builds its own `cmake ... ..` invocation internally
    # with no env-var hook for extra flags, inject the flag via a `cmake`
    # wrapper placed first on PATH rather than patching TWILIGHT's source.
    CMAKE_WRAP_DIR="$(mktemp -d)"
    REAL_CMAKE="$(command -v cmake)"
    # CORRECTED: unconditionally prepending -DCMAKE_CUDA_COMPILER broke
    # `cmake --build ...`/`cmake --install ...` tool-mode calls (used by
    # install_tbb() below) - those flags must be cmake's literal first
    # argument to enter tool mode; prefixing a -D flag first makes cmake
    # instead try (and fail) to treat --build/--parallel/--install as
    # configure-mode arguments. Only inject the flag when the first arg
    # isn't itself a long option (i.e. this is a plain configure call).
    # ALSO FOUND: TWILIGHT's vendored oneTBB 2021.9.0 declares
    # `cmake_minimum_required(VERSION <3.5)` in its own CMakeLists.txt,
    # which this recipe's modern cmake (4.4.2) refuses outright:
    # "Compatibility with CMake < 3.5 has been removed from CMake." CMake's
    # own error message names the fix - add it to the same configure-mode
    # branch of the wrapper.
    # ALSO FOUND: oneTBB's own default CMake target builds its full test
    # suite (we only need the library + TBBConfig.cmake for TWILIGHT to
    # link against), and one of those tests
    # (test/tbb/test_concurrent_monitor.cpp) hits a known GCC 12.4
    # false-positive -Wstringop-overflow in <bits/atomic_base.h> that this
    # build treats as a fatal error. -DTBB_TEST=OFF is oneTBB's own
    # documented option for this and skips the problem entirely; harmless
    # no-op cache var for TWILIGHT's own (unrelated) configure call.
    cat > "${CMAKE_WRAP_DIR}/cmake" <<EOF
#!/bin/bash
if [[ "\$1" == --* ]]; then
    exec "${REAL_CMAKE}" "\$@"
else
    exec "${REAL_CMAKE}" -DCMAKE_CUDA_COMPILER="${BUILD_PREFIX}/bin/nvcc" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DTBB_TEST=OFF "\$@"
fi
EOF
    chmod +x "${CMAKE_WRAP_DIR}/cmake"
    # FOUND (2026-08-06/07, real conda-build run): buildTWILIGHT.sh's
    # find_tbb_dev() shells out to `dpkg -l | grep libtbb-dev` to decide
    # whether to skip its own oneTBB vendor-build and rely on the *host
    # system's* apt-installed TBB instead - a plain `dpkg -l` query that
    # conda-build doesn't (and can't) sandbox. On this host that check
    # matched a stray old `libtbb-dev 2020.1-2` (unrelated apt leftover)
    # which predates oneTBB's CMake package-config support, so
    # `find_package(TBB CONFIG REQUIRED)` then failed to find any
    # TBBConfig.cmake even though the (wrong, old) library was present.
    # Whether or not real bioconda/conda-forge CI images happen to lack
    # libtbb-dev today, depending on that is exactly the kind of
    # host-dependent non-determinism a conda package shouldn't have -
    # force buildTWILIGHT.sh down its own vendored-oneTBB path
    # deterministically by shadowing `dpkg` the same way `cmake` is
    # shadowed above, rather than leaving this to whatever's on the build
    # host.
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
# tree inside the installed package. IMPORTANT: remove it from the PREFIX
# copy, not from $SRC_DIR itself (found via a real conda-build run:
# conda-build's own post-build step calls `source.git_info()` for every
# configured `source:` entry, including this git_url one, to record
# provenance metadata - deleting $SRC_DIR/libpll-src before that step runs
# makes it fail with `AssertionError: assert isdir(src_dir)` and aborts
# packaging entirely).
rm -rf "${PREFIX}/ROADIES/libpll-src"

# Debugging: Verify the contents of the PREFIX directory
echo "Contents of PREFIX:"
ls -al $PREFIX
