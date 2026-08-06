# Draft bioconda recipe for ROADIES + ROADIES_XP

This is a staging copy, not a submitted recipe. It's meant to be copied into a
fork of `bioconda/bioconda-recipes` at `recipes/roadies/` once you're ready to
open the version-bump PR. Diffed against the current live recipe
(`recipes/roadies/meta.yaml` @ v0.1.10), the changes are:

1. Added `twilight`, `kegalign`, `epa-ng`, `gappa` to `run:` — all four are
   already published (twilight/epa-ng/gappa on bioconda, kegalign on
   conda-forge), so this part is a plain, low-risk addition.
2. Added a vendored build of MLIPPER (+ its libpll dependency, which has no
   conda package anywhere) so `conda install roadies` gets full GPU/placement
   support in one shot, instead of waiting on a separate MLIPPER recipe.

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

**Not verified — needs an actual `conda-build` run before opening the PR**:
- Whether `CUDA_HOME="${BUILD_PREFIX}"` is the right value for where
  `cuda-nvcc`/`cuda-cudart-dev` place `nvcc` and the CUDA headers/libs.
  This convention shifts between conda-forge cuda-toolkit packaging
  generations - cross-check against TWILIGHT's or KegAlign's actual build
  logs (both already build CUDA code successfully on bioconda/conda-forge's
  infra), not just their `meta.yaml`.
- Whether MLIPPER's `src/` was actually developed/tested against
  `xflouris/libpll` 0.3.2 specifically, as opposed to `libpll-2`
  (`ddarriba/pll-modules`) or a different version. The Makefile only checks
  for `libpll/pll.h` by path, not by symbol/version, so a wrong-but-header-
  compatible libpll would fail at link time or (worse) misbehave at runtime.
  Confirm with whoever wrote MLIPPER's `src/` before trusting this pin.
- `sha256` for the ROADIES source tarball itself is a placeholder
  (`REPLACE_ME`) - fill in once you've cut the actual release tag.
- The overall recipe hasn't been run through `conda-build` or bioconda's own
  linter (`bioconda-utils lint`) locally. Do that before opening the PR -
  it'll catch anything wrong with the above faster than reading the YAML.

## How to use this

1. Cut the actual ROADIES-XP release tag on `main` (see the branch
   reconciliation note elsewhere - this recipe assumes `main` has the
   ROADIES_XP code, which it doesn't yet as of this draft).
2. Fill in `version` and the `sha256` placeholder in `meta.yaml`.
3. `git clone` your fork of `bioconda-recipes`, copy these two files into
   `recipes/roadies/`, and run bioconda's local build/lint tooling.
4. Resolve whatever the two "not verified" items above surface - almost
   certainly the `CUDA_HOME` path and possibly the libpll version pin.
5. Open the PR.
