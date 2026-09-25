# Shared Utilities

This directory contains reusable parsing, preprocessing, size-checking, and low-level definitions used across MLIPPER. New workflow orchestration should not be added here.

## Files

- `precision.hpp`: compile-time floating-point aliases and host/device math wrappers.
- `checked_size.hpp`: checked size arithmetic for host and device allocations.
- `model_utils.hpp/.cpp`: model parsing, normalization, and validation helpers.
- `msa_preprocess.hpp/.cpp`: alignment compatibility checks and repeated-site compression.
- `mlipper_util.h`: legacy shared CUDA structures and utility primitives that still have multiple low-level consumers.

## Precision Contract

`USE_DOUBLE=1` defines `MLIPPER_USE_DOUBLE`, making:

- `fp_t = double`
- `fp2_t = double2`
- `fp4_t = double4`

The experimental float build uses the corresponding float types. Use helpers such as `fp_exp`, `fp_log`, `fp_fma`, `fp_dot4`, and `make_fp4` in shared host/device code instead of scattering precision-specific branches.

The production and regression configuration is double precision. A float result must not silently replace a double-precision golden.

## `mlipper_util.h`

This header includes shared operation descriptors and CUDA helper code used by tree traversal, placement, likelihood, and derivative kernels. Before moving or deleting a type, search both `.cpp` and `.cu` consumers.

`NodeOpInfo`, `ClvPool`, and `ClvDir` are the low-level encoding of the [CLV operations explained in the tree guide](../tree/README.md#what-is-a-clv-operation). They live here because both C++ and CUDA code consume them, not because `util/` owns the tree algorithm.

### Shared CUDA reduction

`warp_reduce_sum_double()` and `block_reduce_sum_double()` provide the common sum reduction used by root and placement likelihood kernels. Every thread enters with one `double` partial sum:

```text
32 thread values in each warp
  -> shuffle-down additions
  -> lane 0 writes one warp total to shared memory
  -> synchronize the block
  -> first warp reduces the stored warp totals
  -> thread 0 receives the block total
```

The shared array has 32 entries, supporting up to 32 warps or 1,024 threads in a CUDA block. Threads outside the first warp return zero after contributing their warp totals; only the caller-designated block leader writes the final output. This helper performs no grid-wide reduction. Root likelihood adds a second kernel for that step, while placement keeps one independent block total per candidate. See the [likelihood reduction diagrams](../likelihood/README.md#gpu-parallelism-and-outputs).

`checked_size.hpp` protects allocation arithmetic in two stages:

```text
dimensions -> checked multiplication -> element count
element count * sizeof(T)             -> allocation bytes
```

This prevents integer overflow before allocation. It does not guarantee that enough RAM or VRAM exists; the allocator can still report OOM.

Keep additions narrowly scoped:

- put GPU admission and ownership in `gpu/`;
- put branch/model policy in `optimize/optimization_types.hpp`;
- put tree representations in `tree/tree.hpp`;
- put algorithm-specific helpers in their algorithm directory.

Because this header is widely included, adding heavyweight dependencies or non-inline definitions can increase build time or cause duplicate symbols.
