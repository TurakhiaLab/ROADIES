# GPU Utilities

This directory manages GPU admission and small pieces of CUDA memory ownership. It does not implement phylogenetic kernels; those live in `tree/`, `likelihood/`, `placement/`, and `pmatrix/`.

For kernel work assignment and outputs, use the component guides: [tree/CLV traversal](../tree/README.md#how-tree-work-is-assigned-on-the-gpu), [likelihood and reductions](../likelihood/README.md#gpu-parallelism-and-outputs), [placement and derivatives](../placement/README.md#gpu-parallelism-during-placement), and [PMAT construction](../pmatrix/README.md#how-gpu-work-is-assigned). GPU admission itself launches no numerical kernel, so it has no block/site mapping.

## Files

- `gpu_admission.hpp/.cpp`: choose a CUDA device, estimate memory demand, and coordinate multiple MLIPPER processes through shared reservations.
- `device_buffer.hpp`: minimal move-only RAII wrapper around one `cudaMalloc` allocation.

## GPU Selection

Admission is a scheduling estimate, not an allocation:

```text
estimate this job's peak memory
             ↓
read GPU usage + other MLIPPER reservations
             ↓
enough budget? ── yes ──> reserve a slot ──> run cudaMalloc later
       │
       no
       ↓
wait or reject, depending on the requested mode
```

`MlipperGpuAcquireMode` has three modes:

- `UseCurrentDevice`: use the current CUDA device without shared admission.
- `AdmitSpecificDevice`: wait until a requested visible device is admissible.
- `AutoAdmitAnyVisible`: inspect every visible device and reserve the one with the lowest projected reserved-memory ratio.

The main entrypoint is:

```cpp
DeviceReservation select_device_or_wait_or_throw(
    const MlipperGpuConfig&, int estimated_process_memory_mb);
```

Auto selection rejects a GPU when its SM-utilization threshold is exceeded or when the estimated process memory plus safety margin would exceed its memory budget. Selection and reservation are serialized briefly so concurrent jobs do not all choose the same snapshot.

## Key Functions

- `estimate_mlipper_gpu_process_memory_mb()`: apply the [admission estimate](#gpu-selection) to a resident small-tip tree, queries, CLVs, PMATs, scalers, and placement scratch buffers.
- `estimate_divide_and_conquer_gpu_process_memory_mb()`: estimate one bounded D&C subtree plus Local SPR workspace.
- `estimate_dipper_and_divide_and_conquer_gpu_process_memory_mb()`: reserve the larger of the known D&C estimate and the conservative DIPPER fallback.
- `select_device_or_wait_or_throw()`: perform the [selection and reservation decision](#gpu-selection), wait when capacity is temporarily unavailable, and return an RAII handle.
- `ensure_reservation_capacity_or_wait()`: extend the [reservation lifetime contract](#reservation-lifetime) when a small-tip tree requires more capacity; it never drops the original reservation before the larger one is accepted.
- `current_device_or_throw()`, `set_device_or_throw()`, and `current_device_properties_or_throw()`: checked CUDA runtime wrappers.

## Reservation Lifetime

`DeviceReservation` is move-only. Moving it transfers responsibility for the shared admission entry. `reset()` or destruction removes the entry. The entry records PID and Linux process start time so stale entries and PID reuse can be distinguished.

The `MlipperSession::initializeDivideAndConquerGPU()` overload that accepts a `DeviceReservation` takes it by value and moves ownership into the session. This keeps the reservation acquired before DIPPER starting-tree construction alive on the same device throughout subsequent MLIPPER initialization and D&C execution.

The reservation is accounting, not a CUDA-driver memory guarantee. Every real allocation must still check `cudaMalloc` failures. Programs that do not use MLIPPER admission can consume memory after a reservation is granted.

## `DeviceBuffer<T>`

Use `DeviceBuffer<T>` for a simple resizable scratch allocation. Unlike [admission](#gpu-selection), this object owns real CUDA memory:

- `ensureCapacity(count)` keeps an existing large-enough buffer or reallocates.
- `get()` exposes a non-owning raw pointer to kernels.
- `reset()` releases the CUDA allocation.
- copying is disabled; moving transfers ownership.

Do not embed this owner inside `DeviceTree`. Kernel-facing `DeviceTree` must remain trivially copyable and contains only non-owning pointers.

## Runtime Controls

These advanced settings affect `--gpu-auto` admission. They do not affect runs that select a device with `--gpu-id`.

| Environment variable | Default | Accepted value | Effect of lowering the value |
| --- | ---: | --- | --- |
| `MLIPPER_GPU_AUTO_SM_UTILIZATION_MAX_PCT` | `90` | Integer from `1` to `100` | Rejects GPUs at a lower utilization threshold |
| `MLIPPER_GPU_AUTO_MEMORY_BUDGET_FRACTION` | `0.85` | Decimal in `(0, 1]` | Allows MLIPPER to use a smaller fraction of total VRAM |
| `MLIPPER_GPU_AUTO_MEMORY_SAFETY_MB` | `4096` | Positive integer, in MiB | Leaves less unallocated VRAM as a safety margin |
| `MLIPPER_GPU_AUTO_POLL_MS` | `1000` | Positive integer, in milliseconds | Checks for an admissible GPU more frequently |

Set overrides for one invocation by placing them before the command:

```bash
MLIPPER_GPU_AUTO_MEMORY_BUDGET_FRACTION=0.80 \
MLIPPER_GPU_AUTO_MEMORY_SAFETY_MB=6144 \
./MLIPPER --gpu-auto ...
```

To reuse settings in the current shell, export them. Use `unset` to restore a variable's default:

```bash
export MLIPPER_GPU_AUTO_SM_UTILIZATION_MAX_PCT=75
export MLIPPER_GPU_AUTO_POLL_MS=2000
./MLIPPER --gpu-auto ...

unset MLIPPER_GPU_AUTO_SM_UTILIZATION_MAX_PCT
unset MLIPPER_GPU_AUTO_POLL_MS
```

For Docker, pass each override with `-e`:

```bash
docker run --gpus all \
  -e MLIPPER_GPU_AUTO_MEMORY_BUDGET_FRACTION=0.80 \
  -e MLIPPER_GPU_AUTO_MEMORY_SAFETY_MB=6144 \
  IMAGE ./MLIPPER --gpu-auto ...
```

Lower utilization and memory-budget thresholds make admission more conservative. A larger memory safety margin is also more conservative.

When changing accounting, update `tests/gpu_admission_test.cpp` and run a concurrent `--gpu-auto` workflow test on a multi-GPU host.
