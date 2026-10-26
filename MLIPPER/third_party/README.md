# Third-party source

This directory contains external source components that MLIPPER builds
directly. They are tracked by the MLIPPER repository and may be edited when an
integration change requires it.

## DIPPER

`dipper/` provides MLIPPER's starting-tree implementation. It is vendored so a
clean MLIPPER checkout does not depend on a sibling `DIPPER` repository. See
[`dipper/README.md`](dipper/README.md) for the pinned revision, local integration
changes, license, and update procedure.

Code under this directory follows the upstream project structure and is not
covered by MLIPPER's per-folder README convention.
