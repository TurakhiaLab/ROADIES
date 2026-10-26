# Docker Images

This directory contains reproducible CUDA build/runtime environments. Docker does not change MLIPPER's workflow semantics; it supplies the compiler, libraries, executable, and GPU visibility.

## Files

- `Dockerfile`: the development image used by broader MLIPPER experiments.
- `Dockerfile.runtime`: the smaller multi-stage runtime image. The builder compiles `USE_DOUBLE=1`; the runtime stage contains only MLIPPER and required shared libraries.

Use the runtime image to execute MLIPPER. ROADIES currently calls a local
MLIPPER binary and does not start this image automatically.

## Build

From the repository root:

```sh
docker build -f docker/Dockerfile.runtime -t wenchiehlo/mlipper:runtime .
```

The build context must be the MLIPPER repository root. It includes the pinned DIPPER source under `third_party/dipper`; no sibling checkout is required.

Both Dockerfiles default to CUDA 11.8.0 on Ubuntu 22.04. Docker downloads the
matching NVIDIA CUDA base image automatically, so callers do not need to pass a
CUDA build argument. A different tested CUDA image can be selected explicitly:

```sh
docker build -f docker/Dockerfile.runtime \
  --build-arg CUDA_VERSION=12.3.2 \
  -t wenchiehlo/mlipper:runtime .
```

Published image tags:

- `wenchiehlo/mlipper:runtime` for `docker/Dockerfile.runtime`: MLIPPER and
  its required runtime libraries. This is the recommended image for running
  MLIPPER.
- `wenchiehlo/mlipper:latest` points to the same image as `runtime`.
- `wenchiehlo/mlipper:dev` for `docker/Dockerfile`: the source tree, CUDA
  development toolkit, compilers, EPA-ng, RAxML-ng, and ASTRAL-Pro3. Use it
  to rebuild MLIPPER or run development and downstream analysis tasks.

## Run and GPU Visibility

Manual selection:

```sh
docker run --rm --gpus 'device=0' wenchiehlo/mlipper:runtime --help
```

Automatic admission needs every candidate GPU visible and shared host IPC:

```sh
docker run --rm --gpus all --ipc=host wenchiehlo/mlipper:runtime ... --gpu-auto
```

`--ipc=host` lets separate MLIPPER containers see the same POSIX shared-memory reservation tables. Without it, containers can independently admit themselves to the same GPU and defeat cross-process accounting.

## When Editing a Dockerfile

- Keep CUDA development packages in the builder and runtime libraries in the final stage.
- Keep the production build in double precision.
- Update both images when a required runtime shared library changes.
- Test `MLIPPER --help` and one real GPU workflow from the built image.
