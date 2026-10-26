# Integrated DIPPER Source

This directory contains the DIPPER source used by MLIPPER starting-tree and
divide-and-conquer workflows. It is part of the MLIPPER checkout so both
projects can be changed together without a nested Git repository.

## Provenance

- Upstream: `git@github.com:TurakhiaLab/DIPPER.git`
- Base revision: `09da433943b58dc291153e882de53811ebc0931d`
- Source-only branch: `mlipper-integration`
- Imported source commit: `d3f1614261a923b6b9e45ee788f6c0a7273da7c8`
- License: MIT; see `LICENSE`

The vendored snapshot also contains the MLIPPER integration changes that were
present in the development checkout when it was imported. In particular, the
public `DipperSession` API and D&C synchronization paths differ from the base
revision and are required by MLIPPER.

## Development policy

MLIPPER compiles only the source files listed in its top-level `Makefile`.
Changes required in DIPPER may be made directly under `src/` and committed with
the corresponding MLIPPER change. Keep MLIPPER-only workflow adaptation under
MLIPPER's `src/workflow` when it does not require a change to DIPPER itself.

Do not place a DIPPER `.git` directory here. A nested repository would hide
DIPPER changes from MLIPPER's `git status` and those changes would not be
included in an MLIPPER commit.

The full upstream repository is intentionally not embedded: it contains large
datasets, examples, installation files, and project-site material that are not
part of the MLIPPER build. Upstream synchronization should use a source-only
Git subtree rooted at DIPPER's `src/`, not the complete DIPPER repository.

## Upstream synchronization

The current import is an exact copy of the development checkout's tracked
`src/` files at the time recorded above, including its MLIPPER integration
changes. The subtree relationship is recorded in MLIPPER's Git history.

Pull source updates into a clean MLIPPER worktree with:

```sh
git subtree pull \
  --prefix=third_party/dipper/src \
  git@github.com:TurakhiaLab/DIPPER.git \
  mlipper-integration \
  --squash
```

Publish committed MLIPPER-side DIPPER changes with:

```sh
git subtree push \
  --prefix=third_party/dipper/src \
  git@github.com:TurakhiaLab/DIPPER.git \
  mlipper-integration
```

`mlipper-integration` has the contents of DIPPER's `src/` at its repository
root. Do not open a direct merge from that branch into DIPPER `main`; changes
intended for DIPPER `main` must be reviewed and reapplied below its `src/`
directory.

When updating DIPPER:

1. split or fetch the upstream DIPPER `src/` history;
2. merge it only into `third_party/dipper/src`;
3. review conflicts against the MLIPPER integration changes;
4. update the base revision above;
5. run the CPU tests and both GPU workflows before accepting the update.

`UPSTREAM_README.md` is the README from the imported DIPPER checkout.
