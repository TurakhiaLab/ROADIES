#pragma once
#include <cstdlib>   // abort
#include <csignal>   // raise, SIGTRAP
#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstddef>
#include <sstream>
#include <thread>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#include "tree/tree.hpp"

// Branch length defaults to match epa-ng constants.
constexpr double DEFAULT_BRANCH_LENGTH = 0.10536051565782628; // -log(0.9)
// Match the EPA-ng / PLLMOD branch-length floor for placement optimization.
constexpr double OPT_BRANCH_LEN_MIN = 1.0e-4;
constexpr double OPT_BRANCH_LEN_MAX = 100.0;  // PLLMOD_OPT_MAX_BRANCH_LEN
constexpr double OPT_BRANCH_EPSILON = 1.0e-1;
constexpr double OPT_BRANCH_XTOL = OPT_BRANCH_LEN_MIN / 10.0;

template <typename T>
__host__ __device__ inline T scalar_min(T lhs, T rhs) {
    return lhs < rhs ? lhs : rhs;
}

template <typename T>
__host__ __device__ inline T scalar_max(T lhs, T rhs) {
    return lhs > rhs ? lhs : rhs;
}

template <typename T>
__host__ __device__ inline T clamp_scalar(T value, T lower, T upper) {
    if (upper < lower) upper = lower;
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

namespace mlipper {
namespace env {

inline bool env_flag_enabled(const char* name) {
    const char* value = std::getenv(name);
    return value && value[0] && std::string(value) != "0";
}

inline void set_int_env_if_specified(const char* name, int value) {
    if (value < 0) return;
    setenv(name, std::to_string(value).c_str(), 1);
}

inline void set_double_env_if_specified(const char* name, double value) {
    if (value < 0.0) return;
    setenv(name, std::to_string(value).c_str(), 1);
}

} // namespace env

namespace gpu {

struct DeviceReservation {
    int device = -1;
    int fd = -1;
    std::string bus_id;
};

inline std::string trim_ascii_copy(std::string value);
inline std::string normalize_pci_bus_id_or_throw(const std::string& raw_bus_id);

inline int visible_device_count() {
    int count = 0;
    const cudaError_t err = cudaGetDeviceCount(&count);
    if (err == cudaErrorNoDevice) {
        return 0;
    }
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("cudaGetDeviceCount failed: ") + cudaGetErrorString(err));
    }
    return count;
}

inline int current_device_or_throw() {
    int device = -1;
    const cudaError_t err = cudaGetDevice(&device);
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("cudaGetDevice failed: ") + cudaGetErrorString(err));
    }
    return device;
}

inline void set_device_or_throw(int device) {
    const cudaError_t err = cudaSetDevice(device);
    if (err != cudaSuccess) {
        std::ostringstream oss;
        oss << "cudaSetDevice(" << device << ") failed: "
            << cudaGetErrorString(err);
        throw std::runtime_error(oss.str());
    }
}

inline cudaDeviceProp device_properties_or_throw(int device) {
    cudaDeviceProp props{};
    const cudaError_t err = cudaGetDeviceProperties(&props, device);
    if (err != cudaSuccess) {
        std::ostringstream oss;
        oss << "cudaGetDeviceProperties(" << device << ") failed: "
            << cudaGetErrorString(err);
        throw std::runtime_error(oss.str());
    }
    return props;
}

inline cudaDeviceProp current_device_properties_or_throw() {
    return device_properties_or_throw(current_device_or_throw());
}

inline std::string device_bus_id_or_throw(int device) {
    char bus_id[32] = {};
    const cudaError_t err = cudaDeviceGetPCIBusId(
        bus_id,
        static_cast<int>(sizeof(bus_id)),
        device);
    if (err != cudaSuccess) {
        std::ostringstream oss;
        oss << "cudaDeviceGetPCIBusId(" << device << ") failed: "
            << cudaGetErrorString(err);
        throw std::runtime_error(oss.str());
    }
    return normalize_pci_bus_id_or_throw(bus_id);
}

inline std::string gpu_lock_dir() {
    const char* env_dir = std::getenv("MLIPPER_GPU_LOCK_DIR");
    if (env_dir && env_dir[0]) {
        return std::string(env_dir);
    }
    return "/tmp/mlipper_gpu_locks";
}

inline int gpu_auto_poll_ms() {
    const char* env_ms = std::getenv("MLIPPER_GPU_AUTO_POLL_MS");
    if (env_ms && env_ms[0]) {
        char* end = nullptr;
        const long parsed = std::strtol(env_ms, &end, 10);
        if (end && *end == '\0' && parsed > 0 && parsed <= 60000L) {
            return static_cast<int>(parsed);
        }
    }
    return 1000;
}

inline std::string normalize_pci_bus_id_or_throw(const std::string& raw_bus_id) {
    const std::string value = trim_ascii_copy(raw_bus_id);
    const size_t first_colon = value.find(':');
    const size_t second_colon =
        (first_colon == std::string::npos) ? std::string::npos : value.find(':', first_colon + 1);
    const size_t dot =
        (second_colon == std::string::npos) ? std::string::npos : value.find('.', second_colon + 1);
    if (first_colon == std::string::npos ||
        second_colon == std::string::npos ||
        dot == std::string::npos) {
        throw std::runtime_error("Invalid PCI bus id format: '" + raw_bus_id + "'");
    }

    const unsigned long domain =
        std::stoul(value.substr(0, first_colon), nullptr, 16);
    const unsigned long bus =
        std::stoul(value.substr(first_colon + 1, second_colon - first_colon - 1), nullptr, 16);
    const unsigned long device =
        std::stoul(value.substr(second_colon + 1, dot - second_colon - 1), nullptr, 16);
    const unsigned long function =
        std::stoul(value.substr(dot + 1), nullptr, 16);

    std::ostringstream oss;
    oss << std::hex << std::nouppercase << std::setfill('0')
        << std::setw(8) << domain
        << ":" << std::setw(2) << bus
        << ":" << std::setw(2) << device
        << "." << function;
    return oss.str();
}

inline std::string trim_ascii_copy(std::string value) {
    const auto is_space = [](unsigned char ch) {
        return std::isspace(ch) != 0;
    };
    value.erase(
        value.begin(),
        std::find_if_not(
            value.begin(),
            value.end(),
            is_space));
    value.erase(
        std::find_if_not(
            value.rbegin(),
            value.rend(),
            is_space).base(),
        value.end());
    return value;
}

inline std::string sanitize_lock_token(std::string token) {
    for (char& ch : token) {
        const unsigned char byte = static_cast<unsigned char>(ch);
        if (!((byte >= '0' && byte <= '9') ||
              (byte >= 'A' && byte <= 'Z') ||
              (byte >= 'a' && byte <= 'z'))) {
            ch = '_';
        }
    }
    return token;
}

inline std::string device_lock_path_or_throw(const std::string& bus_id) {
    const std::filesystem::path lock_dir(gpu_lock_dir());
    std::error_code ec;
    std::filesystem::create_directories(lock_dir, ec);
    if (ec) {
        throw std::runtime_error(
            "Failed to create GPU lock directory '" + lock_dir.string() +
            "': " + ec.message());
    }
    return (lock_dir / (sanitize_lock_token(bus_id) + ".lock")).string();
}

inline int open_lock_file_or_throw(const std::string& lock_path) {
    const int fd = ::open(lock_path.c_str(), O_RDWR | O_CREAT, 0666);
    if (fd < 0) {
        std::ostringstream oss;
        oss << "open(" << lock_path << ") failed: " << std::strerror(errno);
        throw std::runtime_error(oss.str());
    }
    return fd;
}

inline bool try_reserve_device(
    int device,
    DeviceReservation* reservation_out) {
    if (reservation_out == nullptr) {
        throw std::runtime_error("try_reserve_device requires a non-null output pointer.");
    }

    DeviceReservation reservation{};
    reservation.device = device;
    reservation.bus_id = device_bus_id_or_throw(device);
    const std::string lock_path = device_lock_path_or_throw(reservation.bus_id);
    reservation.fd = open_lock_file_or_throw(lock_path);

    if (::flock(reservation.fd, LOCK_EX | LOCK_NB) == 0) {
        *reservation_out = reservation;
        return true;
    }

    const int lock_errno = errno;
    ::close(reservation.fd);
    if (lock_errno == EWOULDBLOCK || lock_errno == EAGAIN) {
        return false;
    }

    std::ostringstream oss;
    oss << "flock(" << lock_path << ") failed: "
        << std::strerror(lock_errno);
    throw std::runtime_error(oss.str());
}

inline DeviceReservation reserve_specific_device_or_throw(int device) {
    DeviceReservation reservation{};
    if (try_reserve_device(device, &reservation)) {
        return reservation;
    }

    const cudaDeviceProp props = device_properties_or_throw(device);
    const std::string bus_id = device_bus_id_or_throw(device);
    std::ostringstream oss;
    oss << "CUDA device " << device << " (" << props.name
        << ", PCI " << bus_id
        << ") is already reserved by another MLIPPER process.";
    throw std::runtime_error(oss.str());
}

inline DeviceReservation reserve_any_visible_device_or_wait_or_throw() {
    const int device_count = visible_device_count();
    if (device_count <= 0) {
        throw std::runtime_error("No CUDA devices are visible.");
    }

    const int poll_ms = gpu_auto_poll_ms();
    bool announced_wait = false;

    while (true) {
        for (int device = 0; device < device_count; ++device) {
            DeviceReservation reservation{};
            if (try_reserve_device(device, &reservation)) {
                return reservation;
            }
        }

        if (!announced_wait) {
            std::cerr
                << "All " << device_count
                << " visible CUDA devices are currently reserved by other MLIPPER processes; "
                << "waiting for a free GPU...\n";
            announced_wait = true;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms));
    }
}

inline void release_device_reservation(DeviceReservation* reservation) {
    if (reservation == nullptr || reservation->fd < 0) {
        return;
    }
    ::close(reservation->fd);
    reservation->fd = -1;
}

} // namespace gpu
} // namespace mlipper

__host__ __device__ inline double effective_split_branch_min(
    double total_branch_length,
    double min_branch_length = OPT_BRANCH_LEN_MIN)
{
    if (total_branch_length <= 0.0) return min_branch_length;
    return scalar_min(min_branch_length, 0.5 * total_branch_length);
}

__host__ __device__ inline void normalize_split_branch_lengths(
    double total_branch_length,
    double proposed_proximal_length,
    double min_branch_length,
    double& proximal_length_out,
    double& distal_length_out)
{
    if (total_branch_length <= 0.0) {
        proximal_length_out = min_branch_length;
        distal_length_out = min_branch_length;
        return;
    }

    const double lower_bound = effective_split_branch_min(total_branch_length, min_branch_length);
    const double upper_bound = scalar_max(lower_bound, total_branch_length - lower_bound);
    proximal_length_out = clamp_scalar(proposed_proximal_length, lower_bound, upper_bound);
    distal_length_out = total_branch_length - proximal_length_out;
}

__host__ __device__ inline double sanitize_branch_length(
    double branch_length,
    double min_branch_length = OPT_BRANCH_LEN_MIN,
    double max_branch_length = OPT_BRANCH_LEN_MAX,
    double default_branch_length = DEFAULT_BRANCH_LENGTH)
{
    if (!(branch_length > 0.0)) branch_length = default_branch_length;
    return clamp_scalar(branch_length, min_branch_length, max_branch_length);
}

// Basic tags describing the operation type and CLV buffer selection.
enum NodeOpType : int {
    OP_TIP_TIP = 0,
    OP_TIP_INNER = 1,
    OP_INNER_INNER = 2,
    OP_DOWN_INNER_INNER = 3,
    OP_DOWN_INNER_TIP = 4,
    OP_DOWN_TIP_INNER = 5,
    OP_DOWN_TIP_TIP   = 6
};

enum ClvPool : uint8_t {
    CLV_POOL_UP = 0,
    CLV_POOL_DOWN = 1
};

// Direction tags used by preorder/downward passes.
enum ClvDir : uint8_t {
    CLV_DIR_UNSET      = 0,
    CLV_DIR_UP         = 1, // child -> parent
    CLV_DIR_DOWN_LEFT  = 2, // parent -> left child
    CLV_DIR_DOWN_RIGHT = 3  // parent -> right child
};

struct NodeOpInfo {
    int parent_id = -1;
    int left_id = -1;
    int right_id = -1;
    int left_tip_index = -1;
    int right_tip_index = -1;
    int op_type = OP_TIP_TIP;
    uint8_t clv_pool = static_cast<uint8_t>(CLV_POOL_UP);
    uint8_t dir_tag  = static_cast<uint8_t>(CLV_DIR_UP);
};

// Common CUDA error checking macro used across modules.
#ifndef CHECK_CUDA
#define CHECK_CUDA(call)                                                          \
    do {                                                                          \
        cudaError_t err__ = (call);                                               \
        if (err__ != cudaSuccess) {                                               \
            fprintf(stderr,                                                       \
                "CUDA ERROR: %s (%d)\n"                                            \
                "  at %s:%d\n"                                                    \
                "  call: %s\n",                                                   \
                cudaGetErrorString(err__), (int)err__,                            \
                __FILE__, __LINE__, #call);                                       \
            raise(SIGTRAP);   /* <<< 讓 gdb 停在這裡 */                            \
            abort();                                                             \
        }                                                                         \
    } while (0)
#endif

#define CHECK_CUDA_LAST()                                                        \
    do {                                                                         \
        cudaError_t err__ = cudaGetLastError();                                  \
        if (err__ != cudaSuccess) {                                              \
            fprintf(stderr,                                                      \
                "CUDA KERNEL LAUNCH ERROR: %s (%d)\n"                            \
                "  at %s:%d\n",                                                  \
                cudaGetErrorString(err__), (int)err__,                           \
                __FILE__, __LINE__);                                             \
            raise(SIGTRAP);                                                      \
            abort();                                                            \
        }                                                                        \
    } while (0)

// Fast integer log2 ceil for small unsigned values.
inline unsigned int ceil_log2_u32(unsigned int x) {
    if (x <= 1u) return 0u;
    unsigned int v = x - 1u;
    unsigned int r = 0u;
    while (v) { v >>= 1u; ++r; }
    return r;
}

// ===== Device-side CLV helpers (inlined for reuse across CUDA units) =====
__device__ __forceinline__ size_t per_node_span(const DeviceTree& D) {
    return D.sites * static_cast<size_t>(D.rate_cats) * static_cast<size_t>(D.states);
}

__device__ __forceinline__ size_t scaler_span(const DeviceTree& D) {
    if (D.per_rate_scaling) {
        return D.sites * static_cast<size_t>(D.rate_cats);
    }
    return D.sites;
}

__device__ __forceinline__ size_t scaler_site_offset(
    const DeviceTree& D,
    size_t site)
{
    if (D.per_rate_scaling) {
        return site * static_cast<size_t>(D.rate_cats);
    }
    return site;
}

__device__ __forceinline__ unsigned int* scaler_ptr_for_node(
    unsigned int* base,
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    if (!base) return nullptr;
    if (node_id < 0 || node_id >= D.capacity_N) return nullptr;
    return base + static_cast<size_t>(node_id) * scaler_span(D) + scaler_site_offset(D, site);
}

__device__ __forceinline__ unsigned int* up_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_site_scaler_up, D, node_id, site);
}

__device__ __forceinline__ unsigned int* down_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_site_scaler_down, D, node_id, site);
}

__device__ __forceinline__ unsigned int* mid_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_site_scaler_mid, D, node_id, site);
}

__device__ __forceinline__ unsigned int* mid_base_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_site_scaler_mid_base, D, node_id, site);
}

template <typename T>
__device__ __forceinline__ T* clv_write_pool_base(const DeviceTree& D, const NodeOpInfo& op) {
    return (op.clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN))
        ? reinterpret_cast<T*>(D.d_clv_down)
        : reinterpret_cast<T*>(D.d_clv_up);
}

template <typename T>
__device__ __forceinline__ T* clv_read_pool_base(const DeviceTree& D, const NodeOpInfo& op) {
    return (op.clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN))
        ? reinterpret_cast<T*>(D.d_clv_down)
        : reinterpret_cast<T*>(D.d_clv_up);
}

template <typename T>
__device__ __forceinline__ T* clv_write_ptr_for_node(const DeviceTree& D, const NodeOpInfo& op, int node_id) {
    T* base = clv_write_pool_base<T>(D, op);
    return base ? base + static_cast<size_t>(node_id) * per_node_span(D) : nullptr;
}

template <typename T>
__device__ __forceinline__ T* clv_read_ptr_for_node(const DeviceTree& D, const NodeOpInfo& op, int node_id) {
    T* base = clv_read_pool_base<T>(D, op);
    return base ? base + static_cast<size_t>(node_id) * per_node_span(D) : nullptr;
}

// Variant when the pool is implicitly the "up" pool (used in derivative helpers).
template <typename T>
__device__ __forceinline__ T* clv_read_ptr_for_node(const DeviceTree& D, int node_id) {
    T* base = reinterpret_cast<T*>(D.d_clv_up);
    return base ? base + static_cast<size_t>(node_id) * per_node_span(D) : nullptr;
}

__device__ __forceinline__ unsigned int* site_scaler_ptr_base(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    unsigned int rate_cats)
{
    (void)rate_cats;
    if (op.clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN)) {
        return down_scaler_ptr(D, op.parent_id, site);
    }
    return up_scaler_ptr(D, op.parent_id, site);
}
