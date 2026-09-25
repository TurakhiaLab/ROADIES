#include <cassert>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "gpu/gpu_admission.hpp"

int main()
{
    using mlipper::gpu::DeviceReservation;

    static_assert(!std::is_copy_constructible_v<DeviceReservation>);
    static_assert(!std::is_copy_assignable_v<DeviceReservation>);
    static_assert(std::is_nothrow_move_constructible_v<DeviceReservation>);
    static_assert(std::is_nothrow_move_assignable_v<DeviceReservation>);

    const int baseline =
        mlipper::gpu::estimate_mlipper_gpu_process_memory_mb(
            59, 30, 1, 1000, 4, 4, true, false);
    const int more_sites =
        mlipper::gpu::estimate_mlipper_gpu_process_memory_mb(
            59, 30, 1, 2000, 4, 4, true, false);
    assert(baseline > 0);
    assert(more_sites >= baseline);

    const int dnc =
        mlipper::gpu::estimate_divide_and_conquer_gpu_process_memory_mb(
            30, 1000, 4, 4, true);
    assert(dnc == mlipper::gpu::add_local_spr_gpu_process_memory_mb(baseline));
    assert(
        mlipper::gpu::estimate_dipper_and_divide_and_conquer_gpu_process_memory_mb(
            30, 1000, 4, 4, true) >= dnc);

    const int fallback =
        mlipper::gpu::estimate_mlipper_gpu_process_memory_mb(
            0, 0, 0, 0, 0, 0, false, false);
    assert(fallback > 0);

    const int saturated =
        mlipper::gpu::estimate_divide_and_conquer_gpu_process_memory_mb(
            std::numeric_limits<int>::max(),
            std::numeric_limits<size_t>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            true);
    assert(saturated == std::numeric_limits<int>::max());

    DeviceReservation source;
    source.device = 3;
    source.reserved_memory_mb = 1024;
    source.process_start_ticks = 99;
    DeviceReservation destination = std::move(source);
    assert(source.device == -1);
    assert(source.reserved_memory_mb == 0);
    assert(source.process_start_ticks == 0);
    assert(destination.device == 3);
    assert(destination.reserved_memory_mb == 1024);
    assert(destination.process_start_ticks == 99);

    mlipper::MlipperGpuConfig invalid_config;
    invalid_config.acquire_mode =
        mlipper::MlipperGpuAcquireMode::AdmitSpecificDevice;
    bool rejected_invalid_device = false;
    try {
        mlipper::gpu::select_device_or_wait_or_throw(
            invalid_config,
            baseline);
    } catch (const std::runtime_error&) {
        rejected_invalid_device = true;
    }
    assert(rejected_invalid_device);

    return 0;
}
