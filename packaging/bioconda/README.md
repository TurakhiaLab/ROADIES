# Draft bioconda recipe for ROADIES + ROADIES_XP

This is a staging copy, not a submitted recipe. It's meant to be copied into a
fork of `bioconda/bioconda-recipes` at `recipes/roadies/` once you're ready to
open the version-bump PR. Diffed against the current live recipe
(`recipes/roadies/meta.yaml` @ v0.1.10), the changes are:

1. Added `epa-ng`/`gappa` to `run:` (placement-mode deps; already published
   on bioconda, low-risk) and bumped `python ==3.9` → `python ==3.11` to
   match what `roadies_env.sh`'s own `mamba create` line actually installs.
2. Added a vendored build of MLIPPER (+ its libpll dependency, which has no
   conda package anywhere) and of TWILIGHT (also no run: dependency - see
   below) so `conda install roadies` gets full GPU/placement support in one
   shot, instead of waiting on separate MLIPPER/TWILIGHT recipes.
3. Deliberately did **not** add `twilight` or `kegalign` as `run:`
   dependencies, despite both being published packages - see the
   "twilight/kegalign" finding below for why that would have been wrong.

## What's actually verified vs. what still needs a real build

**Verified this session** (checked against live sources, not assumed):
- `twilight`, `epa-ng`, `gappa` exist on bioconda; `kegalign` exists on
  conda-forge (not bioconda proper, but resolves the same way).
- `libpll` does not exist on bioconda or conda-forge (searched both).
- Ubuntu's `libpll-dev` (0.3.2-2) is Debian source package `libpll`,
  Homepage `http://www.libpll.org/` → upstream is `xflouris/libpll`.
- `xflouris/libpll` tag `0.3.2` (no "v" prefix) exists, and its build system
  is autotools (`Makefile.am`, `configure.ac`, `autogen.sh`).
- `MLIPPER/MLIPPER` in this repo is a **prebuilt binary checked into git**,
  not build output (`git ls-files MLIPPER/MLIPPER` shows it tracked). The
  recipe rebuilds it from source rather than shipping that binary - you
  should also consider gitignoring it upstream.

**Verified 2026-08-06 on a GPU host (peregrine), resolving both prior
open questions:**
- **libpll ABI question — resolved, safe.** The GPU host's existing
  `/usr/local` libpll install turned out to actually be **`xflouris/libpll-2`**
  (Alexey Kozlov's fork used by RAxML-NG), built in single-precision (`float`),
  *not* plain `xflouris/libpll` — a genuinely different repo than what this
  recipe vendors. Its `pll_rnode_t.length` field is `float` there vs. `double`
  in plain `libpll`, a real struct-layout difference. However: MLIPPER's own
  `src/` only calls two libpll functions anywhere
  (`grep -rn "pll_" src/`) — `pll_rtree_parse_newick`/`pll_rtree_destroy`, both
  from the precision-independent Newick-parsing subset of the API — it never
  touches the CLV/likelihood functions where the float/double split actually
  lives. Confirmed empirically: rebuilt MLIPPER from source against a freshly
  built vendored `xflouris/libpll` tag `0.3.2` (autotools build into a scratch
  prefix, exactly as `build.sh` does), ran a real GPU placement
  (`--commit-to-tree`/`--jplace-out` on a synthetic 5-tip reference tree +
  1 query) on an A6000, and compared against the repo's already-tracked
  binary (linked against the host's installed libpll-2): **identical tree
  topology and placement, log-likelihoods agree to 7 significant figures**
  (-95.366047080882 vs -95.366046870634 — the residual is `double` vs `float`
  branch-length parsing precision, not a correctness bug). Conclusion: vendor
  either libpll variant safely, *as long as build and link always use a
  matching header+lib pair* (the Makefile's `PLL_INC_DIR`/`PLL_LIB_DIR`
  pairing already guarantees this) — never mix headers from one variant with
  a runtime `.so` from the other.
- **`CUDA_HOME="${BUILD_PREFIX}"` question — partially wrong as first
  written, corrected after a real `conda-build` run (see below).** Installed
  `cuda-nvcc`+`cuda-cudart-dev` (conda-forge, cuda 12.2) into a scratch env:
  `nvcc`, `cuda_runtime.h`, and `libcudart`/`libcuda` (driver stub) all live
  under `$PREFIX/targets/x86_64-linux/{bin,include,lib}`, **not** under
  `$PREFIX/include` or `$PREFIX/lib64` (the latter doesn't even exist). A
  standalone `nvcc` invocation with the "wrong" `-I$(CUDA_HOME)/include
  -L$(CUDA_HOME)/lib64` flags still compiled/linked fine, because `nvcc` has
  its own built-in default search under `<nvcc-dir>/../targets/<arch>/` that
  fires regardless of explicit `-I`/`-L`. **What that first pass missed:**
  MLIPPER's Makefile also compiles its `.cpp` files with plain `g++` (not
  `nvcc`) - those files include `cuda_runtime.h` too (via
  `src/util/precision.hpp`), and `g++` has no such fallback. A real
  `conda-build` run failed every `.cpp` translation unit with
  `cuda_runtime.h: No such file or directory`, exposing the gap. **Fix**:
  `build.sh` now sets `CUDA_HOME="${BUILD_PREFIX}/targets/x86_64-linux"` (so
  `-I$(CUDA_HOME)/include` resolves for `g++` too) and separately pins
  `NVCC="${BUILD_PREFIX}/bin/nvcc"` back to the standard-location binary, so
  its useful default-search fallback stays in play for the final
  `nvcc`-driven link step (where `CUDA_LIB`'s `-L$(CUDA_HOME)/lib64` is still
  technically wrong - real dir is `lib/`, not `lib64/` - but harmless given
  that fallback, exactly as already verified for the link step alone).

**Found by actually running `conda-build` against this recipe (2026-08-06/07,
GPU host, `-c bioconda -c conda-forge`, local `path:` source override to a
clean `git archive` of the working tree + a `conda_build_config.yaml`
pinning `cuda_compiler_version`/`c_compiler`/`cxx_compiler` since bioconda's
own global pinning isn't available outside `bioconda-utils`) — real
blockers, not hypotheticals, found across several build attempts:**
- **`twilight` and `kegalign` never belonged in `run:` at all - superseding
  the "just bump python/cuda-version" fix attempted first.** The solver
  failure that exposed this: `kegalign` pins `tbb >=2020.2,<2021.0.0a0` on
  *every* published build, while `twilight` pins `tbb >=2021.13.0` on every
  build from 0.1.3 onward - disjoint across each package's entire published
  version history (confirmed via `conda search ... --info` on both), so no
  combination of the two can ever coexist in one environment via published
  artifacts. Chasing *why* upstream apparently never hit this revealed the
  real fix: **neither package is actually meant to live in ROADIES-XP's main
  environment.** `TWILIGHT` is never installed as a conda package at
  all - `roadies_env.sh` git-clones `TurakhiaLab/TWILIGHT` and runs
  `install/buildTWILIGHT.sh [cuda]`, which vendors its own oneTBB 2021.9.0
  build from source (sidestepping the conda `tbb` pin entirely), and is
  invoked via the hardcoded path `${roadies_root}/TWILIGHT/bin/twilight`
  (`workflow/scripts/placement.sh`), never looked up on `$PATH`. `kegalign`
  is invoked exclusively by the `kegalign` Snakemake rule
  (`workflow/rules/pair_align_gpu.smk`) via `conda:
  "../envs/kegalign.yaml"` - Snakemake's `--use-conda` gives that rule its
  own fully isolated environment, so it never needs to coexist with anything
  else. **Fix**: removed both from `run:`; `build.sh` now vendor-builds
  TWILIGHT from source the same way it already does for MLIPPER/libpll
  (`git clone` + `install/buildTWILIGHT.sh cuda`); `kegalign` needs no
  packaging action at all - Snakemake fetches it into its own env
  automatically the first time that rule runs (needs bioconda/conda-forge
  channel access at *that* point, not at `roadies` install time). Also
  reverted the `python ==3.9` → `==3.12` / `cuda-version >=12.9,<13` changes
  from the abandoned kegalign-compat attempt: `python ==3.11` (matching
  `roadies_env.sh`'s actual working dev setup) and the original
  `cuda-version >=12.2,<13` are correct again now that kegalign is gone.
- The `CUDA_HOME`/`g++`/`cuda_runtime.h` finding directly above.
- Both found only because the recipe was actually run through `conda-build`
  in its real isolated `$BUILD_PREFIX`/`$PREFIX` sandbox, not just simulated
  against host paths by hand - worth remembering next time something in this
  recipe "looks right" but hasn't been build-tested for real.

## `conda-build` now succeeds end to end (2026-08-07, GPU host)

After the fixes above, plus five more found only by actually running
`conda-build` repeatedly until it got through the whole pipeline (source
solve → `libpll` → `MLIPPER` → `TWILIGHT`+vendored `oneTBB` → `sampling` →
package → `test:` commands), **a full build+test now completes
successfully**: `roadies-0.2.0-py311pl5321h0f1480a_0.conda` (25.7MB) was
produced and all three `test:` commands passed (`run_roadies.py --help`,
`mashtree --help`, and the MLIPPER binary executable check). The five
additional fixes, all in `build.sh`, all found by reading the actual
`conda-build` sandbox error output rather than guessing:

1. **libpll parallel-build race**: libpll's autotools build has a missing
   make dependency between its bison-generated headers
   (`parse_utree.h`/`parse_rtree.h`) and the flex-generated sources that
   `#include` them - under `-j"${CPU_COUNT}"` this intermittently fails with
   `parse_utree.h: No such file or directory`. Hit twice in this session (once
   in ad hoc manual testing, once in a real `conda-build` run); a bare retry
   succeeded cleanly both times. Fixed with `make -j"${CPU_COUNT}" || make
   -j"${CPU_COUNT}"` rather than serializing the whole build to `-j1`.
2. **`cuda-cudart-dev` needs a version pin matching the compiler**: left
   unpinned, the solver picked the newest available (12.9) while `build:`'s
   `{{ compiler('cuda') }}` was pinned to `cuda-version` 12.2 (via
   `conda_build_config.yaml`) - `nvlink` then refused to link with
   `libcudadevrt.a ... newer than toolkit (129 vs 122)`. Fixed by pinning
   `cuda-cudart-dev >=12.2,<12.3` in `host:`, matching its own self-declared
   `cuda-version` dependency for that build.
3. **TWILIGHT's `check_language(CUDA)` doesn't actually respect PATH or the
   documented `CUDACXX` env var** on this cmake version (4.4.2) - only an
   explicit `-DCMAKE_CUDA_COMPILER=...` flag works (confirmed interactively).
   Since `buildTWILIGHT.sh` builds its own `cmake ... ..` call internally with
   no env-var hook, `build.sh` injects the flag via a `cmake` wrapper placed
   first on `PATH` - careful to only do so for configure-mode calls, since
   naively prepending a `-D` flag in front of `cmake --build`/`cmake
   --install` (used by `install_tbb()`) breaks their tool-mode dispatch
   ("Unknown argument --parallel/--install").
4. **`buildTWILIGHT.sh`'s `find_tbb_dev()` shells out to `dpkg -l | grep
   libtbb-dev`** to decide whether to skip its own oneTBB vendor-build in
   favor of the host system's apt TBB - a check conda-build can't sandbox.
   This host had a stray `libtbb-dev 2020.1-2` (unrelated apt leftover)
   that predates oneTBB's CMake package-config support, so
   `find_package(TBB CONFIG REQUIRED)` then failed outright even though a
   (wrong, old) library was present. Whether or not real bioconda/
   conda-forge CI images have `libtbb-dev` today, depending on that is
   host-dependent non-determinism a conda package shouldn't have - forced
   the vendored-oneTBB path deterministically by shadowing `dpkg` (exits 1)
   in the same wrapper directory.
5. **oneTBB 2021.9.0's own `CMakeLists.txt` requires `cmake_minimum_required
   (VERSION <3.5)`**, which cmake 4.4.2 refuses outright ("Compatibility
   with CMake < 3.5 has been removed"). Fixed via the same wrapper, adding
   `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` (the fix cmake's own error message
   names).
6. **oneTBB's default CMake target builds its full test suite** (not needed
   - only the library + `TBBConfig.cmake` matter for TWILIGHT to link
   against), and one test (`test_concurrent_monitor.cpp`) hits a known GCC
   12.4 false-positive `-Wstringop-overflow` in `<bits/atomic_base.h>`
   treated as fatal. Fixed with oneTBB's own documented `-DTBB_TEST=OFF`
   (harmless no-op for TWILIGHT's own unrelated configure call).
7. **TWILIGHT's `CMakeLists.txt` hardcodes `-march=native`** for its
   CUDA-host-compiler and CPU builds. Two separate problems: on this
   particular build host's CPU (AVX512-BF16/AMX capable), `nvcc`'s
   host-compiler pass over GCC 12's AVX512BF16/AMX intrinsics headers failed
   outright (`identifier "__builtin_ia32_..." is undefined`) - a known class
   of `nvcc`/newest-GCC-header incompatibility, independently confirmed as
   an `nvcc` 12.2-vs-GCC-12.4 compatibility gap (NVIDIA's support matrix for
   CUDA 12.2 tops out around GCC 12.1); separately, and regardless of that
   error, `-march=native` is simply **wrong for a redistributable conda
   package** on any build host - it bakes in that specific machine's exact
   CPU feature set, so the resulting binary would `SIGILL` on any end-user
   CPU lacking the same features. Patched to the portable `-march=x86-64-v3`
   baseline (AVX2/FMA/BMI2, ubiquitous since Haswell/2013) via `sed` on the
   freshly cloned copy - this also happened to sidestep the AVX512BF16/AMX
   headers entirely, but would have needed the GCC pin below regardless.
8. **Belt-and-suspenders**: also pinned `c_compiler_version`/
   `cxx_compiler_version` to `11` in `conda_build_config.yaml` (down from
   `12`) for the local test build, since GCC 12.4 vs. `nvcc` 12.2 is a real,
   independently-documented incompatibility class beyond just this one
   `-march` interaction - worth checking what compiler version bioconda's
   own global pinning actually resolves to for CUDA recipes before assuming
   this is unnecessary in the real submission.

**Known packaging-hygiene follow-up, not a correctness issue:** the built
package is larger than necessary because `build.sh`'s `cp -rf * ${PREFIX}/
ROADIES` copies whole build trees (`MLIPPER/build/`, `TWILIGHT/build/`
including the entire `oneTBB-2021.9.0/build/` tree with every intermediate
`.o` file) rather than just the final binaries - fine for validating that
the recipe *works*, worth trimming before a real submission (e.g. `rm -rf
MLIPPER/build TWILIGHT/build/CMakeFiles TWILIGHT/build/oneTBB-2021.9.0`
after the binaries exist, mirroring the existing `libpll-src` cleanup
pattern).

**Still not verified — the one remaining genuine blocker:**
- `sha256` for the ROADIES source tarball itself is a placeholder
  (`REPLACE_ME`) - fill in once you've cut the actual release tag. Every
  other question this recipe raised has now been checked against a real,
  successful `conda-build` run; this is the only thing that structurally
  cannot be verified before the release exists.

## How to use this

1. Cut the actual ROADIES-XP release tag on `main` (see the branch
   reconciliation note elsewhere - this recipe assumes `main` has the
   ROADIES_XP code, which it doesn't yet as of this draft).
2. Fill in `version` and the `sha256` placeholder in `meta.yaml` (swap the
   `source:` entry back from the local test `path:` override if you copied
   this recipe directly - the version in this repo already points at the
   real `url:`/`sha256: REPLACE_ME` form).
3. Consider trimming the build-tree bloat noted above before submitting.
4. `git clone` your fork of `bioconda-recipes`, copy these two files into
   `recipes/roadies/`, and run bioconda's local build/lint tooling -
   bioconda's own global `conda_build_config.yaml` may pin different
   compiler/CUDA versions than the ad hoc one used for local validation
   here; recheck the GCC-12-vs-nvcc-12.2 finding above against whatever
   bioconda actually resolves.
5. Open the PR.
