#include "gpu/gpu_admission.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cctype>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <limits>
#include <mutex>
#include <new>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <utility>
#include <unistd.h>

#include <nvml.h>

#include "util/precision.hpp"
#include "util/checked_size.hpp"

namespace mlipper {
namespace gpu {
namespace {

constexpr int DEFAULT_GPU_AUTO_SM_UTILIZATION_MAX_PCT = 90;
constexpr double DEFAULT_GPU_AUTO_MEMORY_BUDGET_FRACTION = 0.85;
constexpr int DEFAULT_GPU_AUTO_MEMORY_SAFETY_MB = 4096;
constexpr int GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB = 8192;
constexpr int GPU_AUTO_ESTIMATE_OVERHEAD_MB = 512;
constexpr int GPU_AUTO_ESTIMATE_ROUND_MB = 128;
constexpr int GPU_AUTO_LOCAL_SPR_EXTRA_MB = 256;
constexpr size_t GPU_AUTO_PLACEMENT_SCRATCH_CHUNK_OPS = 256;
constexpr size_t GPU_AUTO_SHARED_STATE_CAPACITY = 64;
constexpr int GPU_AUTO_BRANCH_LENGTH_BUFFER_COUNT = 5;
constexpr int GPU_AUTO_TREE_CLV_BUFFER_COUNT = 4;
constexpr int GPU_AUTO_TREE_PMAT_BUFFER_COUNT = 4;
constexpr int GPU_AUTO_TREE_SCALER_BUFFER_COUNT = 4;

struct GpuAutoAdmissionConfig {
    int sm_utilization_max_pct = DEFAULT_GPU_AUTO_SM_UTILIZATION_MAX_PCT;
    double memory_budget_fraction = DEFAULT_GPU_AUTO_MEMORY_BUDGET_FRACTION;
    int memory_safety_mb = DEFAULT_GPU_AUTO_MEMORY_SAFETY_MB;
    int process_memory_mb = GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB;
};

struct GpuRuntimeProbe {
    int sm_utilization_pct = 0;
    int memory_total_mb = 0;
    int memory_used_mb = 0;
};

// Cross-process reservations live in one POSIX shared-memory object per physical
// GPU. CUDA ordinals are process-local and can change under device visibility,
// so object names use the normalized PCI bus ID instead. PID plus Linux process
// start ticks prevents a recycled PID from inheriting a stale reservation.
struct SharedGpuAdmissionEntry {
    pid_t pid = 0;
    int reserved_memory_mb = 0;
    std::uint64_t process_start_ticks = 0;
};

// Fixed-size, versioned POD layout shared by independently started processes.
// Every read or mutation occurs under the corresponding per-device record lock;
// a schema mismatch is reinitialized only while that lock is held.
struct SharedGpuAdmissionState {
    uint32_t magic = 0;
    uint32_t version = 0;
    SharedGpuAdmissionEntry entries[GPU_AUTO_SHARED_STATE_CAPACITY] = {};
};

constexpr uint32_t GPU_AUTO_SHARED_STATE_MAGIC = 0x4d4c4750U; // 'MLGP'
constexpr uint32_t GPU_AUTO_SHARED_STATE_VERSION = 2U;

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

static double parse_fraction_env_or_throw(
    const char* name,
    double default_value) {
    const char* env_value = std::getenv(name);
    if (env_value == nullptr || env_value[0] == '\0') {
        return default_value;
    }

    const std::string value = trim_ascii_copy(env_value);
    if (value.empty()) {
        return default_value;
    }

    size_t parsed_chars = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &parsed_chars);
    } catch (const std::exception&) {
        parsed_chars = 0;
    }
    if (parsed_chars != value.size() ||
        !(parsed > 0.0) ||
        parsed > 1.0) {
        std::ostringstream oss;
        oss << "Invalid " << name << "='" << value
            << "'. Expected a floating-point value in (0, 1].";
        throw std::runtime_error(oss.str());
    }
    return parsed;
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

namespace {

static GpuAutoAdmissionConfig gpu_auto_admission_config_or_throw(
    int process_memory_mb = GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB) {
    GpuAutoAdmissionConfig config{};
    config.sm_utilization_max_pct =
        parse_positive_int_env_or_throw(
            "MLIPPER_GPU_AUTO_SM_UTILIZATION_MAX_PCT",
            DEFAULT_GPU_AUTO_SM_UTILIZATION_MAX_PCT);
    if (config.sm_utilization_max_pct > 100) {
        std::ostringstream oss;
        oss << "Invalid MLIPPER_GPU_AUTO_SM_UTILIZATION_MAX_PCT='"
            << config.sm_utilization_max_pct
            << "'. Expected a value between 1 and 100.";
        throw std::runtime_error(oss.str());
    }
    config.memory_budget_fraction =
        parse_fraction_env_or_throw(
            "MLIPPER_GPU_AUTO_MEMORY_BUDGET_FRACTION",
            DEFAULT_GPU_AUTO_MEMORY_BUDGET_FRACTION);
    config.memory_safety_mb =
        parse_positive_int_env_or_throw(
            "MLIPPER_GPU_AUTO_MEMORY_SAFETY_MB",
            DEFAULT_GPU_AUTO_MEMORY_SAFETY_MB);
    config.process_memory_mb = process_memory_mb > 0
        ? process_memory_mb
        : GPU_AUTO_FALLBACK_PROCESS_MEMORY_MB;
    return config;
}

} // namespace

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

static std::string device_ipc_name_or_throw(
    const char* prefix,
    const std::string& bus_id,
    const char* description) {
    const std::string name = std::string(prefix) + sanitize_lock_token(bus_id);
    if (name.size() >= NAME_MAX) {
        std::ostringstream oss;
        oss << description << " name is too long for '" << bus_id << "'.";
        throw std::runtime_error(oss.str());
    }
    return name;
}

static std::string device_shared_state_name_or_throw(const std::string& bus_id) {
    return device_ipc_name_or_throw(
        "/mlipper_gpu_auto_state_",
        bus_id,
        "POSIX shared-memory");
}

static std::string admission_lock_name_or_throw(const std::string& key) {
    return device_ipc_name_or_throw(
        "/mlipper_gpu_auto_lock_",
        key,
        "GPU admission lock");
}

struct AdmissionLock {
    int fd = -1;

    AdmissionLock() = default;
    ~AdmissionLock() { reset(); }

    void reset() noexcept {
        if (fd >= 0) {
            // POSIX record locks are released when this descriptor is closed.
            ::close(fd);
            fd = -1;
        }
    }

    AdmissionLock(AdmissionLock&& other) noexcept : fd(other.fd) {
        other.fd = -1;
    }
    AdmissionLock& operator=(AdmissionLock&&) = delete;
    AdmissionLock(const AdmissionLock&) = delete;
    AdmissionLock& operator=(const AdmissionLock&) = delete;
};

struct SharedStateMapping {
    int fd = -1;
    SharedGpuAdmissionState* state = nullptr;

    SharedStateMapping() = default;
    ~SharedStateMapping() {
        if (state != nullptr) {
            ::munmap(state, sizeof(SharedGpuAdmissionState));
        }
        if (fd >= 0) {
            ::close(fd);
        }
    }

    SharedStateMapping(SharedStateMapping&& other) noexcept
        : fd(other.fd), state(other.state)
    {
        other.fd = -1;
        other.state = nullptr;
    }
    SharedStateMapping& operator=(SharedStateMapping&&) = delete;
    SharedStateMapping(const SharedStateMapping&) = delete;
    SharedStateMapping& operator=(const SharedStateMapping&) = delete;
};

static AdmissionLock lock_admission_or_throw(const std::string& key) {
    AdmissionLock lock{};
    const std::string lock_name = admission_lock_name_or_throw(key);
    lock.fd = ::shm_open(lock_name.c_str(), O_RDWR | O_CREAT, 0666);
    if (lock.fd < 0) {
        std::ostringstream oss;
        oss << "shm_open(" << lock_name << ") failed: "
            << std::strerror(errno);
        throw std::runtime_error(oss.str());
    }

    struct flock file_lock {};
    file_lock.l_type = F_WRLCK;
    file_lock.l_whence = SEEK_SET;
    while (::fcntl(lock.fd, F_SETLKW, &file_lock) != 0) {
        if (errno == EINTR) {
            continue;
        }
        std::ostringstream oss;
        oss << "fcntl lock(" << lock_name << ") failed: "
            << std::strerror(errno);
        throw std::runtime_error(oss.str());
    }
    return lock;
}

static void nvml_init_or_throw() {
    static std::once_flag init_once;
    static nvmlReturn_t init_status = NVML_SUCCESS;
    std::call_once(
        init_once,
        []() {
            init_status = nvmlInit_v2();
        });
    if (init_status != NVML_SUCCESS) {
        throw std::runtime_error(
            std::string("nvmlInit_v2 failed: ") + nvmlErrorString(init_status));
    }
}

static GpuRuntimeProbe probe_runtime_for_bus_id_or_throw(
    const std::string& bus_id) {
    nvml_init_or_throw();

    nvmlDevice_t device_handle{};
    const nvmlReturn_t handle_status =
        nvmlDeviceGetHandleByPciBusId_v2(bus_id.c_str(), &device_handle);
    if (handle_status != NVML_SUCCESS) {
        std::ostringstream oss;
        oss << "nvmlDeviceGetHandleByPciBusId_v2(" << bus_id << ") failed: "
            << nvmlErrorString(handle_status);
        throw std::runtime_error(oss.str());
    }

    nvmlUtilization_t utilization{};
    const nvmlReturn_t util_status =
        nvmlDeviceGetUtilizationRates(device_handle, &utilization);
    if (util_status != NVML_SUCCESS) {
        std::ostringstream oss;
        oss << "nvmlDeviceGetUtilizationRates(" << bus_id << ") failed: "
            << nvmlErrorString(util_status);
        throw std::runtime_error(oss.str());
    }

    nvmlMemory_t memory{};
    const nvmlReturn_t mem_status =
        nvmlDeviceGetMemoryInfo(device_handle, &memory);
    if (mem_status != NVML_SUCCESS) {
        std::ostringstream oss;
        oss << "nvmlDeviceGetMemoryInfo(" << bus_id << ") failed: "
            << nvmlErrorString(mem_status);
        throw std::runtime_error(oss.str());
    }

    GpuRuntimeProbe probe{};
    probe.sm_utilization_pct = static_cast<int>(utilization.gpu);
    probe.memory_total_mb = static_cast<int>(memory.total / (1024ULL * 1024ULL));
    probe.memory_used_mb = static_cast<int>(memory.used / (1024ULL * 1024ULL));
    return probe;
}

static SharedStateMapping map_shared_state_or_throw(const std::string& bus_id) {
    SharedStateMapping mapping{};
    const std::string state_name = device_shared_state_name_or_throw(bus_id);
    mapping.fd = ::shm_open(state_name.c_str(), O_RDWR | O_CREAT, 0666);
    if (mapping.fd < 0) {
        std::ostringstream oss;
        oss << "shm_open(" << state_name << ") failed: " << std::strerror(errno);
        throw std::runtime_error(oss.str());
    }

    if (::ftruncate(mapping.fd, static_cast<off_t>(sizeof(SharedGpuAdmissionState))) != 0) {
        const int truncate_errno = errno;
        std::ostringstream oss;
        oss << "ftruncate(" << state_name << ") failed: "
            << std::strerror(truncate_errno);
        throw std::runtime_error(oss.str());
    }

    void* mapped = ::mmap(
        nullptr,
        sizeof(SharedGpuAdmissionState),
        PROT_READ | PROT_WRITE,
        MAP_SHARED,
        mapping.fd,
        0);
    if (mapped == MAP_FAILED) {
        const int mmap_errno = errno;
        std::ostringstream oss;
        oss << "mmap(" << state_name << ") failed: "
            << std::strerror(mmap_errno);
        throw std::runtime_error(oss.str());
    }

    mapping.state = static_cast<SharedGpuAdmissionState*>(mapped);
    // The per-device admission lock is held by every caller while mapping and
    // initializing this state, so replacing an old schema cannot race readers.
    if (mapping.state->magic != GPU_AUTO_SHARED_STATE_MAGIC ||
        mapping.state->version != GPU_AUTO_SHARED_STATE_VERSION) {
        ::new (mapping.state) SharedGpuAdmissionState{};
        mapping.state->magic = GPU_AUTO_SHARED_STATE_MAGIC;
        mapping.state->version = GPU_AUTO_SHARED_STATE_VERSION;
    }
    return mapping;
}

static bool pid_is_alive(pid_t pid) {
    if (pid <= 0) {
        return false;
    }
    if (::kill(pid, 0) == 0) {
        return true;
    }
    return errno == EPERM;
}

static std::optional<std::uint64_t> process_start_ticks(pid_t pid) {
    std::ifstream stat_file(
        "/proc/" + std::to_string(static_cast<long long>(pid)) + "/stat");
    std::string stat_line;
    if (!std::getline(stat_file, stat_line)) {
        return std::nullopt;
    }

    // The command name is parenthesized and may contain spaces. Fields after
    // the final ')' begin at field 3; Linux starttime is field 22.
    const size_t command_end = stat_line.rfind(')');
    if (command_end == std::string::npos || command_end + 2 >= stat_line.size()) {
        return std::nullopt;
    }
    std::istringstream fields(stat_line.substr(command_end + 2));
    std::string field;
    for (int field_number = 3; field_number <= 22; ++field_number) {
        if (!(fields >> field)) {
            return std::nullopt;
        }
        if (field_number == 22) {
            try {
                size_t parsed_chars = 0;
                const std::uint64_t ticks =
                    std::stoull(field, &parsed_chars, 10);
                if (parsed_chars == field.size()) {
                    return ticks;
                }
            } catch (const std::exception&) {
            }
            return std::nullopt;
        }
    }
    return std::nullopt;
}

static bool process_identity_is_alive(
    pid_t pid,
    std::uint64_t expected_start_ticks)
{
    if (!pid_is_alive(pid)) {
        return false;
    }
    if (expected_start_ticks == 0) {
        return true;
    }
    const std::optional<std::uint64_t> observed_start_ticks =
        process_start_ticks(pid);
    // PID values are reusable, so start ticks distinguish the reserving process
    // from a later process with the same PID. If /proc is restricted, retain a
    // live PID rather than incorrectly dropping an active reservation.
    return !observed_start_ticks.has_value() ||
        *observed_start_ticks == expected_start_ticks;
}

static void cleanup_stale_shared_admissions(SharedGpuAdmissionState* state) {
    if (state == nullptr) {
        return;
    }
    size_t write_index = 0;
    for (size_t read_index = 0; read_index < GPU_AUTO_SHARED_STATE_CAPACITY; ++read_index) {
        const SharedGpuAdmissionEntry entry = state->entries[read_index];
        if (!process_identity_is_alive(
                entry.pid,
                entry.process_start_ticks) ||
            entry.reserved_memory_mb <= 0) {
            continue;
        }
        if (write_index != read_index) {
            state->entries[write_index] = entry;
        }
        ++write_index;
    }
    for (; write_index < GPU_AUTO_SHARED_STATE_CAPACITY; ++write_index) {
        state->entries[write_index] = SharedGpuAdmissionEntry{};
    }
}

static int total_reserved_shared_memory_mb(const SharedGpuAdmissionState* state) {
    if (state == nullptr) {
        return 0;
    }
    long long total_mb = 0;
    for (const SharedGpuAdmissionEntry& entry : state->entries) {
        if (entry.pid > 0 && entry.reserved_memory_mb > 0) {
            total_mb += entry.reserved_memory_mb;
        }
    }
    if (total_mb > std::numeric_limits<int>::max()) {
        return std::numeric_limits<int>::max();
    }
    return static_cast<int>(total_mb);
}

static int matching_reservation_memory_mb(
    const SharedGpuAdmissionState* state,
    pid_t pid,
    std::uint64_t process_start_ticks_value)
{
    if (state == nullptr) {
        return 0;
    }
    for (const SharedGpuAdmissionEntry& entry : state->entries) {
        if (entry.pid == pid &&
            entry.process_start_ticks == process_start_ticks_value) {
            return std::max(0, entry.reserved_memory_mb);
        }
    }
    return 0;
}

static void add_shared_admission_or_throw(
    SharedGpuAdmissionState* state,
    pid_t pid,
    int reserved_memory_mb,
    std::uint64_t process_start_ticks_value) {
    if (state == nullptr) {
        throw std::runtime_error("Shared GPU admission state is null.");
    }
    for (SharedGpuAdmissionEntry& entry : state->entries) {
        if (entry.pid == pid &&
            entry.process_start_ticks == process_start_ticks_value) {
            entry.reserved_memory_mb = reserved_memory_mb;
            return;
        }
    }
    for (SharedGpuAdmissionEntry& entry : state->entries) {
        if (entry.pid <= 0 || entry.reserved_memory_mb <= 0) {
            entry.pid = pid;
            entry.reserved_memory_mb = reserved_memory_mb;
            entry.process_start_ticks = process_start_ticks_value;
            return;
        }
    }
    throw std::runtime_error(
        "Shared GPU admission table is full. Increase GPU_AUTO_SHARED_STATE_CAPACITY.");
}

static void remove_shared_admission(
    SharedGpuAdmissionState* state,
    pid_t pid,
    std::uint64_t process_start_ticks_value) noexcept {
    if (state == nullptr || pid <= 0) {
        return;
    }
    for (SharedGpuAdmissionEntry& entry : state->entries) {
        if (entry.pid == pid &&
            entry.process_start_ticks == process_start_ticks_value) {
            entry = SharedGpuAdmissionEntry{};
            break;
        }
    }
    cleanup_stale_shared_admissions(state);
}

enum class AdmissionAttempt {
    Admitted,
    TemporarilyBusy,
    InsufficientTotalMemory,
};

struct AdmissionAssessment {
    AdmissionAttempt attempt = AdmissionAttempt::TemporarilyBusy;
    long long projected_memory_mb = 0;
    int memory_budget_mb = 0;
};

static AdmissionAssessment assess_device_admission(
    const SharedGpuAdmissionState* state,
    const GpuRuntimeProbe& probe,
    const GpuAutoAdmissionConfig& config,
    pid_t pid,
    std::uint64_t start_ticks)
{
    const int existing_reservation_mb = matching_reservation_memory_mb(
        state,
        pid,
        start_ticks);
    // NVML usage gradually includes memory from admitted processes, while the
    // reservation table includes their projected peak. Taking the maximum is
    // conservative without double-counting allocations present in both views.
    const int accounted_used_mb = std::max(
        probe.memory_used_mb,
        total_reserved_shared_memory_mb(state));
    const int memory_budget_mb = static_cast<int>(
        static_cast<double>(probe.memory_total_mb) *
        config.memory_budget_fraction);
    const long long minimum_required_memory_mb =
        static_cast<long long>(config.process_memory_mb) +
        config.memory_safety_mb;
    const long long projected_memory_mb =
        static_cast<long long>(accounted_used_mb) +
        std::max(0, config.process_memory_mb - existing_reservation_mb) +
        config.memory_safety_mb;

    AdmissionAssessment assessment;
    assessment.projected_memory_mb = projected_memory_mb;
    assessment.memory_budget_mb = memory_budget_mb;
    if (minimum_required_memory_mb > memory_budget_mb) {
        assessment.attempt = AdmissionAttempt::InsufficientTotalMemory;
    } else if (
        probe.sm_utilization_pct > config.sm_utilization_max_pct ||
        projected_memory_mb > memory_budget_mb) {
        assessment.attempt = AdmissionAttempt::TemporarilyBusy;
    } else {
        assessment.attempt = AdmissionAttempt::Admitted;
    }
    return assessment;
}

static AdmissionAssessment inspect_shared_device(
    int device,
    const GpuAutoAdmissionConfig& config,
    pid_t pid,
    std::uint64_t start_ticks)
{
    const std::string bus_id = device_bus_id_or_throw(device);
    AdmissionLock lock = lock_admission_or_throw(bus_id);
    SharedStateMapping state_mapping = map_shared_state_or_throw(bus_id);
    cleanup_stale_shared_admissions(state_mapping.state);
    const GpuRuntimeProbe probe = probe_runtime_for_bus_id_or_throw(bus_id);
    return assess_device_admission(
        state_mapping.state,
        probe,
        config,
        pid,
        start_ticks);
}

static AdmissionAttempt try_admit_shared_device(
    int device,
    DeviceReservation* reservation_out,
    const GpuAutoAdmissionConfig& config) {
    if (reservation_out == nullptr) {
        throw std::runtime_error("try_admit_shared_device requires a non-null output pointer.");
    }

    const std::string bus_id = device_bus_id_or_throw(device);
    const bool updating_existing =
        reservation_out->device == device &&
        !reservation_out->bus_id.empty() &&
        reservation_out->pid > 0;
    if (updating_existing && reservation_out->bus_id != bus_id) {
        throw std::runtime_error(
            "Existing GPU reservation no longer matches its CUDA device.");
    }
    const pid_t pid = updating_existing
        ? reservation_out->pid
        : ::getpid();
    const std::uint64_t start_ticks = updating_existing
        ? reservation_out->process_start_ticks
        : process_start_ticks(pid).value_or(0);

    AdmissionLock lock = lock_admission_or_throw(bus_id);
    SharedStateMapping state_mapping =
        map_shared_state_or_throw(bus_id);
    cleanup_stale_shared_admissions(state_mapping.state);

    const GpuRuntimeProbe probe =
        probe_runtime_for_bus_id_or_throw(bus_id);
    const AdmissionAssessment assessment = assess_device_admission(
        state_mapping.state,
        probe,
        config,
        pid,
        start_ticks);
    if (assessment.attempt != AdmissionAttempt::Admitted) {
        return assessment.attempt;
    }

    add_shared_admission_or_throw(
        state_mapping.state,
        pid,
        config.process_memory_mb,
        start_ticks);
    reservation_out->device = device;
    reservation_out->bus_id = bus_id;
    reservation_out->reserved_memory_mb = config.process_memory_mb;
    reservation_out->pid = pid;
    reservation_out->process_start_ticks = start_ticks;
    return AdmissionAttempt::Admitted;
}

} // namespace

static DeviceReservation admit_specific_device_or_wait_or_throw(
    int device,
    int estimated_process_memory_mb) {
    const int device_count = visible_device_count();
    if (device_count <= 0) {
        throw std::runtime_error("No CUDA devices are visible.");
    }
    if (device < 0 || device >= device_count) {
        std::ostringstream oss;
        oss << "--gpu-id " << device
            << " is out of range for the current process; "
            << device_count << " CUDA device"
            << (device_count == 1 ? " is" : "s are")
            << " visible.";
        throw std::runtime_error(oss.str());
    }

    const GpuAutoAdmissionConfig config =
        gpu_auto_admission_config_or_throw(estimated_process_memory_mb);
    const int poll_ms = gpu_auto_poll_ms_or_default();
    bool announced_wait = false;
    while (true) {
        DeviceReservation reservation;
        const AdmissionAttempt attempt =
            try_admit_shared_device(device, &reservation, config);
        if (attempt == AdmissionAttempt::Admitted) {
            return reservation;
        }
        if (attempt == AdmissionAttempt::InsufficientTotalMemory) {
            std::ostringstream oss;
            oss << "CUDA device " << device
                << " does not have enough total memory for the "
                << config.process_memory_mb << " MiB estimated process memory "
                << "and " << config.memory_safety_mb
                << " MiB safety margin.";
            throw std::runtime_error(oss.str());
        }
        if (!announced_wait) {
            std::cerr
                << "CUDA device " << device
                << " is currently above the shared admission thresholds; "
                << "waiting to retry...\n";
            announced_wait = true;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms));
    }
}

static DeviceReservation admit_any_visible_device_or_wait_or_throw(
    int estimated_process_memory_mb) {
    const int device_count = visible_device_count();
    if (device_count <= 0) {
        throw std::runtime_error("No CUDA devices are visible.");
    }

    const GpuAutoAdmissionConfig config =
        gpu_auto_admission_config_or_throw(estimated_process_memory_mb);
    const int poll_ms = gpu_auto_poll_ms_or_default();
    const pid_t pid = ::getpid();
    const std::uint64_t start_ticks =
        process_start_ticks(pid).value_or(0);
    bool announced_wait = false;
    while (true) {
        // Lock order is global selection first, then one per-device lock at a
        // time during inspection/admission. Specific-device callers acquire
        // only a device lock, so this order cannot form a reverse dependency.
        // Serialize only the short choose-and-reserve operation to prevent
        // concurrent auto-mode processes from all selecting the same snapshot.
        AdmissionLock selection_lock =
            lock_admission_or_throw("all_visible_devices");

        int best_device = -1;
        long double best_projected_load =
            std::numeric_limits<long double>::infinity();
        int insufficient_device_count = 0;
        for (int device = 0; device < device_count; ++device) {
            const AdmissionAssessment assessment =
                inspect_shared_device(device, config, pid, start_ticks);
            if (assessment.attempt == AdmissionAttempt::InsufficientTotalMemory) {
                ++insufficient_device_count;
                continue;
            }
            if (assessment.attempt != AdmissionAttempt::Admitted) {
                continue;
            }

            const long double projected_load =
                static_cast<long double>(assessment.projected_memory_mb) /
                assessment.memory_budget_mb;
            if (projected_load < best_projected_load) {
                best_device = device;
                best_projected_load = projected_load;
            }
        }

        if (best_device >= 0) {
            DeviceReservation reservation;
            if (try_admit_shared_device(best_device, &reservation, config) ==
                AdmissionAttempt::Admitted) {
                return reservation;
            }
            // A specific-device admission or external GPU workload may have
            // changed the selected device between inspection and reservation.
            continue;
        }

        if (insufficient_device_count == device_count) {
            std::ostringstream oss;
            oss << "No candidate CUDA device has enough total memory for the "
                << config.process_memory_mb << " MiB estimated process memory "
                << "and " << config.memory_safety_mb
                << " MiB safety margin.";
            throw std::runtime_error(oss.str());
        }

        // Do not hold the global selection lock while waiting for capacity.
        selection_lock.reset();
        if (!announced_wait) {
            std::cerr
                << "All " << device_count
                << " visible CUDA devices are currently above the shared "
                << "admission thresholds; waiting to retry...\n";
            announced_wait = true;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms));
    }
}

DeviceReservation select_device_or_wait_or_throw(
    const MlipperGpuConfig& config,
    int estimated_process_memory_mb)
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
        reservation = admit_specific_device_or_wait_or_throw(
            config.gpu_id,
            estimated_process_memory_mb);
        set_device_or_throw(reservation.device);
        break;
    case MlipperGpuAcquireMode::AutoAdmitAnyVisible:
        reservation = admit_any_visible_device_or_wait_or_throw(
            estimated_process_memory_mb);
        set_device_or_throw(reservation.device);
        break;
    }
    // Validate UseCurrentDevice as well; this reports a CUDA error immediately
    // when no usable device is available.
    current_device_or_throw();
    return reservation;
}

void ensure_reservation_capacity_or_wait(
    DeviceReservation& reservation,
    int estimated_process_memory_mb)
{
    if (reservation.device < 0 ||
        estimated_process_memory_mb <= reservation.reserved_memory_mb) {
        return;
    }

    const GpuAutoAdmissionConfig config =
        gpu_auto_admission_config_or_throw(estimated_process_memory_mb);
    const int poll_ms = gpu_auto_poll_ms_or_default();
    bool announced_wait = false;
    while (true) {
        const AdmissionAttempt attempt = try_admit_shared_device(
            reservation.device,
            &reservation,
            config);
        if (attempt == AdmissionAttempt::Admitted) {
            return;
        }
        if (attempt == AdmissionAttempt::InsufficientTotalMemory) {
            std::ostringstream oss;
            oss << "CUDA device " << reservation.device
                << " cannot grow this process reservation to "
                << config.process_memory_mb << " MiB with the configured "
                << config.memory_safety_mb << " MiB safety margin.";
            throw std::runtime_error(oss.str());
        }
        if (!announced_wait) {
            std::cerr
                << "CUDA device " << reservation.device
                << " cannot yet grow this process reservation to "
                << config.process_memory_mb << " MiB; waiting to retry...\n";
            announced_wait = true;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms));
    }
}

void DeviceReservation::reset() noexcept {
    if (!bus_id.empty() && pid > 0) {
        try {
            AdmissionLock lock = lock_admission_or_throw(bus_id);
            SharedStateMapping state_mapping =
                map_shared_state_or_throw(bus_id);
            remove_shared_admission(
                state_mapping.state,
                pid,
                process_start_ticks);
        } catch (...) {
            // reset() and destruction cannot report IPC cleanup failures. The
            // stale-entry sweep will reclaim this reservation after exit.
        }
    }
    device = -1;
    bus_id.clear();
    reserved_memory_mb = 0;
    pid = 0;
    process_start_ticks = 0;
}

DeviceReservation::~DeviceReservation() noexcept {
    reset();
}

DeviceReservation::DeviceReservation(DeviceReservation&& other) noexcept
    : device(other.device),
      bus_id(std::move(other.bus_id)),
      reserved_memory_mb(other.reserved_memory_mb),
      pid(other.pid),
      process_start_ticks(other.process_start_ticks)
{
    other.device = -1;
    other.reserved_memory_mb = 0;
    other.pid = 0;
    other.process_start_ticks = 0;
}

DeviceReservation& DeviceReservation::operator=(
    DeviceReservation&& other) noexcept
{
    if (this != &other) {
        reset();
        device = other.device;
        bus_id = std::move(other.bus_id);
        reserved_memory_mb = other.reserved_memory_mb;
        pid = other.pid;
        process_start_ticks = other.process_start_ticks;
        other.device = -1;
        other.reserved_memory_mb = 0;
        other.pid = 0;
        other.process_start_ticks = 0;
    }
    return *this;
}

} // namespace gpu
} // namespace mlipper
