#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <sys/types.h>

#include <cuda_runtime.h>

namespace mlipper {

enum class MlipperGpuAcquireMode {
    // Select a device without cross-process memory accounting.
    UseCurrentDevice,
    // Wait until gpu_id can accommodate this process's estimated demand.
    AdmitSpecificDevice,
    // Wait for and select the least-loaded admissible visible device.
    AutoAdmitAnyVisible
};

struct MlipperGpuConfig {
    MlipperGpuAcquireMode acquire_mode =
        MlipperGpuAcquireMode::UseCurrentDevice;
    int gpu_id = -1;
};

namespace gpu {

// Process-local ownership of one entry in the shared GPU admission table.
// Destruction removes the reservation on a best-effort basis; stale entries
// are also reclaimed using PID and process-start-time identity checks.
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
    int reserved_memory_mb = 0;
    pid_t pid = 0;
    std::uint64_t process_start_ticks = 0;
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

/// Selects and activates a CUDA device, waiting when admission is temporary.
/// A returned empty reservation means cross-process admission was not requested.
DeviceReservation select_device_or_wait_or_throw(
    const MlipperGpuConfig& config,
    int estimated_process_memory_mb);

/// Grows an existing reservation to estimated_process_memory_mb, if necessary.
void ensure_reservation_capacity_or_wait(
    DeviceReservation& reservation,
    int estimated_process_memory_mb);

} // namespace gpu
} // namespace mlipper
