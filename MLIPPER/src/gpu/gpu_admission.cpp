#include "gpu/gpu_admission.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/file.h>
#include <sys/stat.h>
#include <thread>
#include <utility>
#include <unistd.h>

#include "util/precision.hpp"
#include "util/checked_size.hpp"

namespace mlipper {
namespace gpu {
namespace {

constexpr int GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB = 8192;
constexpr int GPU_AUTO_ESTIMATE_OVERHEAD_MB = 512;
constexpr int GPU_AUTO_ESTIMATE_ROUND_MB = 128;
constexpr int GPU_AUTO_LOCAL_SPR_EXTRA_MB = 256;
constexpr size_t GPU_AUTO_PLACEMENT_SCRATCH_CHUNK_OPS = 256;
constexpr int GPU_AUTO_BRANCH_LENGTH_BUFFER_COUNT = 5;
constexpr int GPU_AUTO_TREE_CLV_BUFFER_COUNT = 4;
constexpr int GPU_AUTO_TREE_PMAT_BUFFER_COUNT = 4;
constexpr int GPU_AUTO_TREE_SCALER_BUFFER_COUNT = 4;

static std::string trim_ascii_copy(std::string value);
static std::string normalize_pci_bus_id_or_throw(const std::string& raw_bus_id);
static int parse_positive_int_env_or_throw(const char* name, int default_value);
static int gpu_auto_poll_ms_or_default() {
    return parse_positive_int_env_or_throw("MLIPPER_GPU_AUTO_POLL_MS", 1000);
}

static int parse_positive_int_env_or_throw(
    const char* name,
    int default_value) {
    const char* env_value = std::getenv(name);
    if (env_value == nullptr || env_value[0] == '\0') {
        return default_value;
    }

    const std::string value = trim_ascii_copy(env_value);
    if (value.empty()) {
        return default_value;
    }

    size_t parsed_chars = 0;
    long parsed = 0;
    try {
        parsed = std::stol(value, &parsed_chars, 10);
    } catch (const std::exception&) {
        parsed_chars = 0;
    }
    if (parsed_chars != value.size() ||
        parsed <= 0 ||
        parsed > std::numeric_limits<int>::max()) {
        std::ostringstream oss;
        oss << "Invalid " << name << "='" << value
            << "'. Expected a positive integer.";
        throw std::runtime_error(oss.str());
    }
    return static_cast<int>(parsed);
}

static size_t round_up_placement_scratch_ops(size_t required_ops) {
    if (required_ops == 0) return 0;
    const size_t adjusted = util::checked_add_size(
        required_ops,
        GPU_AUTO_PLACEMENT_SCRATCH_CHUNK_OPS - 1,
        "placement scratch rounding");
    return util::checked_mul_size(
        adjusted / GPU_AUTO_PLACEMENT_SCRATCH_CHUNK_OPS,
        GPU_AUTO_PLACEMENT_SCRATCH_CHUNK_OPS,
        "placement scratch rounding");
}

} // namespace

int estimate_mlipper_gpu_process_memory_mb(
    int node_count,
    int tip_count,
    int query_count,
    size_t site_count,
    int states,
    int rate_cats,
    bool per_rate_scaling,
    bool commit_to_tree)
{
    if (node_count <= 0 ||
        tip_count <= 0 ||
        query_count < 0 ||
        site_count == 0 ||
        states <= 0 ||
        rate_cats <= 0) {
        return GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB;
    }

    constexpr long double kMib = 1024.0L * 1024.0L;
    const size_t reserve_inserts =
        commit_to_tree ? static_cast<size_t>(query_count) : 0;
    const size_t capacity_nodes =
        static_cast<size_t>(node_count) + 2 * reserve_inserts;
    const size_t capacity_tips =
        static_cast<size_t>(tip_count) + reserve_inserts;
    const size_t query_capacity = static_cast<size_t>(query_count);
    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    const long double sites = static_cast<long double>(site_count);
    const long double state_count_ld = static_cast<long double>(state_count);
    const long double rate_count_ld = static_cast<long double>(rate_count);
    const long double per_node = sites * rate_count_ld * state_count_ld;
    const long double pmat_per_node =
        rate_count_ld * state_count_ld * state_count_ld;
    const long double scaler_span =
        per_rate_scaling ? (sites * rate_count_ld) : sites;

    long double bytes = 0.0L;
    auto add = [&](long double elements, long double element_bytes) {
        bytes += elements * element_bytes;
    };

    add(rate_count_ld * state_count_ld, sizeof(fp_t));       // d_lambdas
    add(state_count_ld * state_count_ld * 2, sizeof(fp_t));  // d_V, d_Vinv
    add(state_count_ld, sizeof(fp_t));                       // d_frequencies
    add(rate_count_ld, sizeof(fp_t));                        // d_rate_weights

    add(
        static_cast<long double>(capacity_nodes) *
            GPU_AUTO_BRANCH_LENGTH_BUFFER_COUNT,
        sizeof(fp_t));
    add(static_cast<long double>(capacity_tips) * sites, sizeof(uint8_t));
    add(capacity_tips, sizeof(int));

    add(
        static_cast<long double>(capacity_nodes) * per_node *
            GPU_AUTO_TREE_CLV_BUFFER_COUNT,
        sizeof(fp_t));
    add(1, sizeof(double));
    add(
        static_cast<long double>(capacity_nodes) * pmat_per_node *
            GPU_AUTO_TREE_PMAT_BUFFER_COUNT,
        sizeof(fp_t));
    add(sites, sizeof(unsigned));

    if (query_capacity > 0) {
        add(static_cast<long double>(query_capacity) * sites, sizeof(uint8_t));
        add(static_cast<long double>(query_capacity) * per_node, sizeof(fp_t));
    }

    add(
        static_cast<long double>(capacity_nodes) * scaler_span *
            GPU_AUTO_TREE_SCALER_BUFFER_COUNT,
        sizeof(unsigned));
    add((states == 4) ? 16 : (state_count + 1), sizeof(unsigned));

    if (query_capacity > 0) {
        const size_t required_ops =
            commit_to_tree ? capacity_nodes : static_cast<size_t>(node_count);
        const size_t scratch_ops = round_up_placement_scratch_ops(required_ops);
        add(static_cast<long double>(scratch_ops) * per_node, sizeof(fp_t));
        add(scratch_ops, sizeof(fp_t));
        add(static_cast<long double>(scratch_ops) * pmat_per_node, sizeof(fp_t));
        add(static_cast<long double>(scratch_ops), sizeof(fp_t) + sizeof(int));
    }

    bytes += static_cast<long double>(GPU_AUTO_ESTIMATE_OVERHEAD_MB) * kMib;

    const long double raw_mb = std::ceil(bytes / kMib);
    const long double rounded_mb =
        std::ceil(raw_mb / static_cast<long double>(GPU_AUTO_ESTIMATE_ROUND_MB)) *
        static_cast<long double>(GPU_AUTO_ESTIMATE_ROUND_MB);
    if (rounded_mb > static_cast<long double>(std::numeric_limits<int>::max())) {
        return std::numeric_limits<int>::max();
    }
    return std::max(1, static_cast<int>(rounded_mb));
}

int add_local_spr_gpu_process_memory_mb(int estimated_process_memory_mb) {
    if (estimated_process_memory_mb <= 0) {
        estimated_process_memory_mb = GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB;
    }
    if (estimated_process_memory_mb >
        std::numeric_limits<int>::max() - GPU_AUTO_LOCAL_SPR_EXTRA_MB) {
        return std::numeric_limits<int>::max();
    }
    return estimated_process_memory_mb + GPU_AUTO_LOCAL_SPR_EXTRA_MB;
}

int estimate_divide_and_conquer_gpu_process_memory_mb(
    int tip_budget,
    size_t site_count,
    int states,
    int rate_cats,
    bool per_rate_scaling)
{
    const int subtree_tip_budget = std::max(1, tip_budget);
    const long long node_count =
        2LL * static_cast<long long>(subtree_tip_budget) - 1;
    const int estimated_subtree_nodes = static_cast<int>(
        std::min<long long>(
            std::numeric_limits<int>::max(),
            std::max<long long>(3, node_count)));
    return add_local_spr_gpu_process_memory_mb(
        estimate_mlipper_gpu_process_memory_mb(
            estimated_subtree_nodes,
            subtree_tip_budget,
            1,
            site_count,
            states,
            rate_cats,
            per_rate_scaling,
            false));
}

int estimate_dipper_and_divide_and_conquer_gpu_process_memory_mb(
    int tip_budget,
    size_t site_count,
    int states,
    int rate_cats,
    bool per_rate_scaling)
{
    return std::max(
        GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB,
        estimate_divide_and_conquer_gpu_process_memory_mb(
            tip_budget,
            site_count,
            states,
            rate_cats,
            per_rate_scaling));
}

static int visible_device_count() {
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

int current_device_or_throw() {
    int device = -1;
    const cudaError_t err = cudaGetDevice(&device);
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("cudaGetDevice failed: ") + cudaGetErrorString(err));
    }
    return device;
}

void set_device_or_throw(int device) {
    const cudaError_t err = cudaSetDevice(device);
    if (err != cudaSuccess) {
        std::ostringstream oss;
        oss << "cudaSetDevice(" << device << ") failed: "
            << cudaGetErrorString(err);
        throw std::runtime_error(oss.str());
    }
}

namespace {

static cudaDeviceProp device_properties_or_throw(int device) {
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

} // namespace

cudaDeviceProp current_device_properties_or_throw() {
    return device_properties_or_throw(current_device_or_throw());
}

namespace {

static std::string device_bus_id_or_throw(int device) {
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

static std::string normalize_pci_bus_id_or_throw(const std::string& raw_bus_id) {
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

static std::string trim_ascii_copy(std::string value) {
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

static std::string sanitize_lock_token(std::string token) {
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

// Default lock directory for the per-device flock() files below; override
// with MLIPPER_GPU_LOCK_DIR (e.g. when /tmp is not shared or not desired).
static std::string gpu_lock_dir() {
    const char* env_dir = std::getenv("MLIPPER_GPU_LOCK_DIR");
    if (env_dir != nullptr && env_dir[0] != '\0') {
        return std::string(env_dir);
    }
    return "/tmp/mlipper_gpu_locks";
}

static std::string device_lock_path_or_throw(const std::string& bus_id) {
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

// The lock file is a shared resource: whichever user's MLIPPER process runs
// first on this host creates it for a given GPU, and every other user's
// MLIPPER process must still be able to open and flock() it later. Reset
// umask around creation and best-effort fchmod() afterward so it really is
// 0666 regardless of the creating process's umask -- otherwise a lock file
// created by one user silently blocks every other user from ever opening it
// again, which is exactly the bug this replaces.
static int open_world_writable_lock_file_or_throw(const std::string& lock_path) {
    const ::mode_t previous_umask = ::umask(0);
    const int fd = ::open(lock_path.c_str(), O_RDWR | O_CREAT, 0666);
    const int open_errno = errno;
    ::umask(previous_umask);
    if (fd < 0) {
        std::ostringstream oss;
        oss << "open(" << lock_path << ") failed: "
            << std::strerror(open_errno);
        throw std::runtime_error(oss.str());
    }
    ::fchmod(fd, 0666);
    return fd;
}

static bool try_reserve_device(int device, DeviceReservation* reservation_out) {
    if (reservation_out == nullptr) {
        throw std::runtime_error(
            "try_reserve_device requires a non-null output pointer.");
    }

    DeviceReservation reservation{};
    reservation.device = device;
    reservation.bus_id = device_bus_id_or_throw(device);
    const std::string lock_path = device_lock_path_or_throw(reservation.bus_id);
    reservation.fd = open_world_writable_lock_file_or_throw(lock_path);

    if (::flock(reservation.fd, LOCK_EX | LOCK_NB) == 0) {
        *reservation_out = std::move(reservation);
        return true;
    }

    const int lock_errno = errno;
    ::close(reservation.fd);
    reservation.fd = -1;
    if (lock_errno == EWOULDBLOCK || lock_errno == EAGAIN) {
        return false;
    }

    std::ostringstream oss;
    oss << "flock(" << lock_path << ") failed: " << std::strerror(lock_errno);
    throw std::runtime_error(oss.str());
}

static DeviceReservation reserve_specific_device_or_throw(int device) {
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

static DeviceReservation reserve_any_visible_device_or_wait_or_throw() {
    const int device_count = visible_device_count();
    if (device_count <= 0) {
        throw std::runtime_error("No CUDA devices are visible.");
    }

    const int poll_ms = gpu_auto_poll_ms_or_default();
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
                << " visible CUDA devices are currently reserved by other "
                << "MLIPPER processes; waiting for a free GPU...\n";
            announced_wait = true;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms));
    }
}

} // namespace

DeviceReservation select_device_or_wait_or_throw(
    const MlipperGpuConfig& config,
    int /*estimated_process_memory_mb*/)
{
    DeviceReservation reservation;
    switch (config.acquire_mode) {
    case MlipperGpuAcquireMode::UseCurrentDevice:
        if (config.gpu_id >= 0) {
            set_device_or_throw(config.gpu_id);
        }
        break;
    case MlipperGpuAcquireMode::AdmitSpecificDevice:
        if (config.gpu_id < 0) {
            throw std::runtime_error(
                "gpu_id must be >= 0 when admitting a specific CUDA device.");
        }
        reservation = reserve_specific_device_or_throw(config.gpu_id);
        set_device_or_throw(reservation.device);
        break;
    case MlipperGpuAcquireMode::AutoAdmitAnyVisible:
        reservation = reserve_any_visible_device_or_wait_or_throw();
        set_device_or_throw(reservation.device);
        break;
    }
    // Validate UseCurrentDevice as well; this reports a CUDA error immediately
    // when no usable device is available.
    current_device_or_throw();
    return reservation;
}

void ensure_reservation_capacity_or_wait(
    DeviceReservation& /*reservation*/,
    int /*estimated_process_memory_mb*/)
{
    // No-op: this is a binary per-device mutex, not a memory budget, so
    // there is nothing to grow.
}

void DeviceReservation::reset() noexcept {
    if (fd >= 0) {
        ::close(fd);
        fd = -1;
    }
    device = -1;
    bus_id.clear();
}

DeviceReservation::~DeviceReservation() noexcept {
    reset();
}

DeviceReservation::DeviceReservation(DeviceReservation&& other) noexcept
    : device(other.device),
      bus_id(std::move(other.bus_id)),
      fd(other.fd)
{
    other.device = -1;
    other.fd = -1;
}

DeviceReservation& DeviceReservation::operator=(
    DeviceReservation&& other) noexcept
{
    if (this != &other) {
        reset();
        device = other.device;
        bus_id = std::move(other.bus_id);
        fd = other.fd;
        other.device = -1;
        other.fd = -1;
    }
    return *this;
}

} // namespace gpu
} // namespace mlipper
