#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

#include <cuda_runtime.h>

#include "util/checked_size.hpp"

namespace mlipper::gpu {

// Move-only ownership for a cudaMalloc allocation. Kernel-facing structures
// receive raw non-owning pointers via get(); this object remains the owner.
template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    ~DeviceBuffer() { reset(); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& other) noexcept
        : pointer_(std::exchange(other.pointer_, nullptr)),
          capacity_(std::exchange(other.capacity_, 0))
    {}

    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept
    {
        if (this != &other) {
            reset();
            pointer_ = std::exchange(other.pointer_, nullptr);
            capacity_ = std::exchange(other.capacity_, 0);
        }
        return *this;
    }

    void ensureCapacity(size_t count)
    {
        if (count <= capacity_) return;
        const size_t bytes = util::checked_allocation_bytes<T>(
            count, "DeviceBuffer::ensureCapacity");
        T* replacement = nullptr;
        const cudaError_t status = cudaMalloc(&replacement, bytes);
        if (status != cudaSuccess) {
            throw std::runtime_error(
                std::string("cudaMalloc failed: ") +
                cudaGetErrorString(status));
        }

        // Commit only after allocation succeeds so callers retain the previous
        // usable buffer when growth fails.
        reset();
        pointer_ = replacement;
        capacity_ = count;
    }

    void reset() noexcept
    {
        // Destruction cannot report cudaFree failures. A later checked CUDA call
        // remains responsible for surfacing a sticky runtime error.
        if (pointer_) cudaFree(pointer_);
        pointer_ = nullptr;
        capacity_ = 0;
    }

    T* get() noexcept { return pointer_; }
    const T* get() const noexcept { return pointer_; }
    size_t capacity() const noexcept { return capacity_; }
    explicit operator bool() const noexcept { return pointer_ != nullptr; }

private:
    T* pointer_ = nullptr;
    size_t capacity_ = 0;
};

} // namespace mlipper::gpu
