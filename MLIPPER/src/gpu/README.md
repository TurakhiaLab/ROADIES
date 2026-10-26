# GPU Utilities

This directory manages GPU device selection, memory-demand estimates (for logging/diagnostics), and small pieces of CUDA memory ownership. It does not implement phylogenetic kernels; those live in `tree/`, `likelihood/`, `placement/`, and `pmatrix/`.

For kernel work assignment and outputs, use the component guides: [tree/CLV traversal](../tree/README.md#how-tree-work-is-assigned-on-the-gpu), [likelihood and reductions](../likelihood/README.md#gpu-parallelism-and-outputs), [placement and derivatives](../placement/README.md#gpu-parallelism-during-placement), and [PMAT construction](../pmatrix/README.md#how-gpu-work-is-assigned). GPU selection itself launches no numerical kernel, so it has no block/site mapping.

## Files

- `gpu_admission.hpp/.cpp`: choose a CUDA device and estimate memory demand; prevent two MLIPPER processes from running on the same physical GPU at once.
- `device_buffer.hpp`: minimal move-only RAII wrapper around one `cudaMalloc` allocation.

## GPU Selection

Device selection is a mutual-exclusion lock, not memory-aware admission:

```text
pick a candidate visible device
             ↓
flock() its per-device lock file (non-blocking)
             ↓
got the lock? ── yes ──> use this device, run cudaMalloc later
       │
       no
       ↓
try the next visible device (AutoAdmitAnyVisible),
or wait and retry (AdmitSpecificDevice / all devices busy)
```

`MlipperGpuAcquireMode` has three modes:

- `UseCurrentDevice`: use the current CUDA device with no lock at all.
- `AdmitSpecificDevice`: wait until the requested visible device's lock is free.
- `AutoAdmitAnyVisible`: take the first visible device whose lock is free, waiting and retrying if all are currently held.

The main entrypoint is:

```cpp
DeviceReservation select_device_or_wait_or_throw(
    const MlipperGpuConfig&, int estimated_process_memory_mb);
```

`estimated_process_memory_mb` is accepted for call-site compatibility but unused here — this lock has no notion of a memory budget, just one GPU per concurrent MLIPPER process. The lock file lives under `MLIPPER_GPU_LOCK_DIR` (default `/tmp/mlipper_gpu_locks`), one file per GPU PCI bus ID, created world-writable (umask reset + `fchmod`) so the first user on a shared host to run MLIPPER never blocks any other user from running it afterward.

## Key Functions

- `estimate_mlipper_gpu_process_memory_mb()`: a conservative MiB estimate for a resident small-tip tree, queries, CLVs, PMATs, scalers, and placement scratch buffers — informational only, not consumed by selection.
- `estimate_divide_and_conquer_gpu_process_memory_mb()`: estimate one bounded D&C subtree plus Local SPR workspace.
- `estimate_dipper_and_divide_and_conquer_gpu_process_memory_mb()`: the larger of the known D&C estimate and the conservative DIPPER fallback.
- `select_device_or_wait_or_throw()`: perform the [selection and lock decision](#gpu-selection), wait when every candidate device is currently held, and return an RAII handle.
- `ensure_reservation_capacity_or_wait()`: no-op, kept for call-site compatibility — this lock model has no reservation size to grow.
- `current_device_or_throw()`, `set_device_or_throw()`, and `current_device_properties_or_throw()`: checked CUDA runtime wrappers.

## Reservation Lifetime

`DeviceReservation` is move-only and owns one open file descriptor on its device's lock file. Moving it transfers that ownership. `reset()` or destruction closes the descriptor, which the kernel uses to release the `flock()` automatically — including if the process crashes — so there is no persistent cross-process state to go stale or get stuck.

The `MlipperSession::initializeDivideAndConquerGPU()` overload that accepts a `DeviceReservation` takes it by value and moves ownership into the session. This keeps the lock acquired before DIPPER starting-tree construction held on the same device throughout subsequent MLIPPER initialization and D&C execution.

The lock only prevents two MLIPPER processes from sharing a device; it is not a CUDA-driver memory guarantee. Every real allocation must still check `cudaMalloc` failures, and non-MLIPPER processes on the same GPU are not coordinated with at all.

## `DeviceBuffer<T>`

Use `DeviceBuffer<T>` for a simple resizable scratch allocation. Unlike the [device lock](#gpu-selection), this object owns real CUDA memory:

- `ensureCapacity(count)` keeps an existing large-enough buffer or reallocates.
- `get()` exposes a non-owning raw pointer to kernels.
- `reset()` releases the CUDA allocation.
- copying is disabled; moving transfers ownership.

Do not embed this owner inside `DeviceTree`. Kernel-facing `DeviceTree` must remain trivially copyable and contains only non-owning pointers.

## Runtime Controls

These settings affect the per-device lock used by both `--gpu-id` and `--gpu-auto`.

| Environment variable | Default | Accepted value | Effect |
| --- | ---: | --- | --- |
| `MLIPPER_GPU_LOCK_DIR` | `/tmp/mlipper_gpu_locks` | A directory path | Where per-device `.lock` files are created; override when `/tmp` isn't shared or desired |
| `MLIPPER_GPU_AUTO_POLL_MS` | `1000` | Positive integer, in milliseconds | How often `--gpu-auto` or a busy `--gpu-id` retries for a free device |

Set overrides for one invocation by placing them before the command:

```bash
MLIPPER_GPU_LOCK_DIR=/scratch/me/mlipper_gpu_locks ./MLIPPER --gpu-auto ...
```

For Docker, pass each override with `-e`:

```bash
docker run --gpus all \
  -e MLIPPER_GPU_LOCK_DIR=/scratch/me/mlipper_gpu_locks \
  IMAGE ./MLIPPER --gpu-auto ...
```

When changing this file, update `tests/gpu_admission_test.cpp` and run a concurrent `--gpu-auto` workflow test on a multi-GPU host.
