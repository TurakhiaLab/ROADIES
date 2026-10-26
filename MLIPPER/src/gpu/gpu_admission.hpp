#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <sys/types.h>

#include <cuda_runtime.h>

namespace mlipper {

enum class MlipperGpuAcquireMode {
    // Select a device without any cross-process coordination.
    UseCurrentDevice,
    // Wait until gpu_id is not held by another MLIPPER process.
    AdmitSpecificDevice,
    // Wait for and select a visible device not held by another MLIPPER process.
    AutoAdmitAnyVisible
};

struct MlipperGpuConfig {
    MlipperGpuAcquireMode acquire_mode =
        MlipperGpuAcquireMode::UseCurrentDevice;
    int gpu_id = -1;
};

namespace gpu {

// Process-local ownership of an exclusive, advisory lock on one GPU, so two
// MLIPPER processes don't run on the same device at once. The lock is a
// flock() on a per-device file under MLIPPER_GPU_LOCK_DIR (default
// /tmp/mlipper_gpu_locks); the OS releases it automatically when this
// process's file descriptor is closed, including on crash, so there is no
// persistent cross-process state to go stale or get stuck.
struct DeviceReservation {
    DeviceReservation() = default;
    ~DeviceReservation() noexcept;

    DeviceReservation(DeviceReservation&& other) noexcept;
    DeviceReservation& operator=(DeviceReservation&& other) noexcept;

    DeviceReservation(const DeviceReservation&) = delete;
    DeviceReservation& operator=(const DeviceReservation&) = delete;

    void reset() noexcept;

    int device = -1;
    std::string bus_id;
    int fd = -1;
};

/// Returns a conservative MiB estimate for resident tree and placement state.
int estimate_mlipper_gpu_process_memory_mb(
    int node_count,
    int tip_count,
    int query_count,
    size_t site_count,
    int states,
    int rate_cats,
    bool per_rate_scaling,
    bool commit_to_tree);

/// Adds the fixed local-SPR allowance to an existing MiB estimate.
int add_local_spr_gpu_process_memory_mb(int estimated_process_memory_mb);

int estimate_divide_and_conquer_gpu_process_memory_mb(
    int tip_budget,
    size_t site_count,
    int states,
    int rate_cats,
    bool per_rate_scaling);

int estimate_dipper_and_divide_and_conquer_gpu_process_memory_mb(
    int tip_budget,
    size_t site_count,
    int states,
    int rate_cats,
    bool per_rate_scaling);

int current_device_or_throw();
void set_device_or_throw(int device);
cudaDeviceProp current_device_properties_or_throw();

/// Selects and activates a CUDA device, waiting when every candidate device
/// is currently locked by another MLIPPER process. estimated_process_memory_mb
/// is accepted for call-site compatibility but unused: device selection here
/// is a simple per-device mutual-exclusion lock, not memory-aware admission.
/// A returned empty reservation means no cross-process lock was requested
/// (UseCurrentDevice).
DeviceReservation select_device_or_wait_or_throw(
    const MlipperGpuConfig& config,
    int estimated_process_memory_mb);

/// No-op: this simpler lock model has no notion of growing a reservation's
/// memory budget. Kept for call-site compatibility.
void ensure_reservation_capacity_or_wait(
    DeviceReservation& reservation,
    int estimated_process_memory_mb);

} // namespace gpu
} // namespace mlipper
