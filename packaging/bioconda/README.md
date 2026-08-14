# Draft bioconda recipe for ROADIES + ROADIES_XP

Staging copy, not a submitted recipe. Meant to be copied into a fork of
`bioconda/bioconda-recipes` at `recipes/roadies/` once ready to open the
version-bump PR. Diffed against the live recipe (`recipes/roadies/meta.yaml`
@ v0.1.10):

1. Added `epa-ng`/`gappa` to `run:` (placement-mode deps, already published,
   low-risk) and bumped `python ==3.9` → `==3.11` to match `roadies_env.sh`'s
   actual dev setup.
2. Vendors a build of MLIPPER (+ its libpll dependency, which has no conda
   package anywhere) and of TWILIGHT, so `conda install roadies` gets full
   GPU/placement support in one shot instead of waiting on separate recipes.
3. Deliberately does **not** add `twilight`/`kegalign` to `run:`, despite
   both being published packages - see below.

## What's verified

- `twilight`, `epa-ng`, `gappa` exist on bioconda; `kegalign` exists on
  conda-forge (not bioconda proper, resolves the same way). `libpll` exists
  on neither.
- `libpll` traces to Ubuntu's `libpll-dev` 0.3.2-2 → upstream
  `xflouris/libpll`, tag `0.3.2`, autotools build.
- `MLIPPER/MLIPPER` was a prebuilt binary checked into git, not build
  output - now rebuilt from source in `build.sh` and gitignored upstream.
- **libpll ABI**: the GPU host's existing libpll install is actually
  `libpll-2` (single-precision fork used by RAxML-NG), not plain `libpll` -
  a real struct-layout difference in general. Doesn't matter here: MLIPPER's
  `src/` only calls two libpll functions (`pll_rtree_parse_newick`/
  `pll_rtree_destroy`), both from the precision-independent Newick-parsing
  API. Confirmed by rebuilding MLIPPER against the vendored `libpll` 0.3.2
  and running a real GPU placement - matched the repo's existing
  libpll-2-linked binary on topology, log-likelihoods agreeing to 7
  significant figures (residual is float/double parsing precision, not a
  bug). Safe either way, as long as build and link always use a matching
  header+lib pair (guaranteed by the Makefile's `PLL_INC_DIR`/`PLL_LIB_DIR`).
- **CUDA_HOME layout**: conda's `cuda-nvcc`/`cuda-cudart-dev` put
  `nvcc`/headers/libs under `$BUILD_PREFIX/targets/x86_64-linux/`, not the
  plain `include`/`lib64` MLIPPER's Makefile defaults to. `nvcc` itself has
  a fallback search path that papers over this for `.cu` files, but
  MLIPPER's `.cpp` files are compiled by plain `g++`, which has no such
  fallback - real `conda-build` failed there with `cuda_runtime.h: No such
  file or directory` until `CUDA_HOME` was pointed at the right directory.
- **`twilight`/`kegalign` don't belong in `run:`**: their published builds
  pin disjoint `tbb` version ranges across their entire history (kegalign
  `<2021`, twilight `>=2021.13`) - unsolvable in one environment. Turns out
  moot: TWILIGHT is never installed as a conda package (invoked via a
  hardcoded path, vendors its own oneTBB); kegalign runs exclusively inside
  its own isolated Snakemake `--use-conda` environment. Fixed by removing
  both from `run:` and vendor-building TWILIGHT the same way as MLIPPER.

## `conda-build` succeeds end to end (2026-08-07, GPU host)

A full build+test completed successfully -
`roadies-0.2.0-py311pl5321h0f1480a_0.conda` (25.7MB), all `test:` commands
passing. Fixes along the way, all in `build.sh` unless noted:

1. libpll's autotools build has a bison/flex header-generation race under
   parallel `make` - fixed with a retry instead of forcing `-j1`.
2. `cuda-cudart-dev` needs a version pin matching the build compiler
   (`cuda-version` in `conda_build_config.yaml`), or the solver can pick a
   newer cudart than the toolkit and `nvlink` refuses to link.
3. TWILIGHT's `check_language(CUDA)` doesn't reliably auto-detect `nvcc`
   under `$BUILD_PREFIX` - injected via a `cmake` wrapper on `PATH` rather
   than patching TWILIGHT's source, since `buildTWILIGHT.sh` builds its own
   `cmake` invocation with no env-var hook.
4. `buildTWILIGHT.sh` checks for a host-installed `libtbb-dev` via `dpkg`,
   which isn't sandboxable and is host-dependent either way - forced its
   vendored-oneTBB path deterministically by shadowing `dpkg` too.
5. oneTBB 2021.9.0's own `cmake_minimum_required(<3.5)` is rejected by
   modern cmake (`-DCMAKE_POLICY_VERSION_MINIMUM=3.5` fixes it).
6. oneTBB's default target builds its own test suite, which hits a GCC 12.4
   `-Wstringop-overflow` false positive - skipped with `-DTBB_TEST=OFF`.
7. TWILIGHT's `CMakeLists.txt` hardcodes `-march=native` - both a
   portability bug (SIGILLs on end-user CPUs lacking the build host's exact
   features) and, on this host, an nvcc/GCC-12-AVX512BF16-header conflict.
   Patched to `-march=x86-64-v3` (portable, sidesteps the header issue too).
8. Also pinned `c_compiler_version`/`cxx_compiler_version` to 11 in
   `conda_build_config.yaml` for local testing, since GCC 12 vs. nvcc 12.2
   is an independently-documented incompatibility beyond just `-march`.

**Packaging-hygiene follow-up, not a correctness issue:** the built package
is larger than necessary because `cp -rf * ${PREFIX}/ROADIES` copies whole
build trees (intermediate `.o` files etc.) rather than just final binaries -
worth trimming before a real submission.

**Still not verified:** `sha256` for the ROADIES source tarball is a
placeholder (`REPLACE_ME`) until the release tag exists - the only thing
left that can't be checked before then.

## How to use this

1. Cut the actual ROADIES-XP release tag on `main` (this recipe assumes
   `main` has the ROADIES_XP code, which it doesn't yet as of this draft -
   see the branch-reconciliation note in project memory).
2. Fill in `version` and `sha256` in `meta.yaml`.
3. Consider trimming the build-tree bloat noted above before submitting.
4. Copy into a `bioconda-recipes` fork at `recipes/roadies/` and run their
   local build/lint tooling - bioconda's own global `conda_build_config.yaml`
   may pin different compiler/CUDA versions than the ad hoc one used for
   local validation here; recheck the GCC-12-vs-nvcc-12.2 finding against
   whatever bioconda actually resolves.
5. Open the PR.
