#include "pmat_gpu.cuh"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "util/checked_size.hpp"

namespace {

constexpr int kMaxStates = 16;

void validate_pmat_inputs(
    int state_count,
    int rate_categories,
    const fp_t* branch_lengths,
    const fp_t* eigenvectors,
    const fp_t* inverse_eigenvectors,
    const fp_t* rate_eigenvalues,
    fp_t* pmats)
{
    if (state_count <= 0 || state_count > kMaxStates ||
        rate_categories <= 0 || !branch_lengths || !eigenvectors ||
        !inverse_eigenvectors || !rate_eigenvalues || !pmats) {
        throw std::invalid_argument("Invalid PMAT rebuild arguments.");
    }
}

__global__ void build_all_branch_pmats_kernel(
    int node_count,
    int state_count,
    int rate_categories,
    const fp_t* branch_lengths,
    const fp_t* eigenvectors,
    const fp_t* inverse_eigenvectors,
    const fp_t* rate_eigenvalues,
    fp_t* pmats)
{
    const int item =
        static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int item_count = node_count * rate_categories;
    if (item >= item_count) return;
    const int node_id = item / rate_categories;
    const int rate_id = item - node_id * rate_categories;
    const size_t matrix_elems =
        static_cast<size_t>(state_count) * state_count;
    pmatrix_from_triple_gpu(
        inverse_eigenvectors,
        eigenvectors,
        rate_eigenvalues +
            static_cast<size_t>(rate_id) * state_count,
        fp_t(1),
        branch_lengths[node_id],
        fp_t(0),
        pmats + static_cast<size_t>(item) * matrix_elems,
        state_count);
}

__global__ void build_single_branch_pmat_kernel(
    int node_id,
    int state_count,
    int rate_categories,
    const fp_t* branch_lengths,
    const fp_t* eigenvectors,
    const fp_t* inverse_eigenvectors,
    const fp_t* rate_eigenvalues,
    fp_t* pmats)
{
    extern __shared__ fp_t diag_expm1[];
    const int tid = static_cast<int>(threadIdx.x);
    const int diagonal_count = rate_categories * state_count;
    const fp_t branch_length = branch_lengths[node_id];
    for (int diagonal = tid; diagonal < diagonal_count;
         diagonal += static_cast<int>(blockDim.x)) {
        diag_expm1[diagonal] = fp_expm1(
            rate_eigenvalues[diagonal] * fp_t(1) * branch_length);
    }
    __syncthreads();

    const size_t matrix_elems =
        static_cast<size_t>(state_count) * state_count;
    const size_t output_base =
        static_cast<size_t>(node_id) * rate_categories * matrix_elems;
    const size_t element_count =
        static_cast<size_t>(rate_categories) * matrix_elems;
    for (size_t element = static_cast<size_t>(tid);
         element < element_count;
         element += blockDim.x) {
        const int rate_id = static_cast<int>(element / matrix_elems);
        const size_t matrix_element = element -
            static_cast<size_t>(rate_id) * matrix_elems;
        const int row = static_cast<int>(matrix_element / state_count);
        const int col = static_cast<int>(matrix_element -
            static_cast<size_t>(row) * state_count);
        fp_t acc = fp_t(0);
        for (int eig_idx = 0; eig_idx < state_count; ++eig_idx) {
            const fp_t v_entry =
                eigenvectors[row * state_count + eig_idx];
            const fp_t vinv_entry =
                inverse_eigenvectors[eig_idx * state_count + col];
            acc += v_entry *
                diag_expm1[rate_id * state_count + eig_idx] *
                vinv_entry;
        }
        if (row == col) acc += fp_t(1);
        if (acc < fp_t(0) && acc > fp_t(-1e-12)) acc = fp_t(0);
        pmats[output_base + element] = acc;
    }
}

} // namespace

__device__ void pmatrix_from_triple_gpu(
    const fp_t* Vinv,
    const fp_t* V,
    const fp_t* rate_eigenvalues,
    fp_t rate_scale,
    fp_t branch_length,
    fp_t pinv,
    fp_t* out_pmat,
    int state_count)
{
    if (!Vinv || !V || !rate_eigenvalues || !out_pmat) return;
    if (state_count <= 0 || state_count > kMaxStates) return;

    fp_t diag_expm1[kMaxStates];
    for (int state_idx = 0; state_idx < state_count; ++state_idx) {
        diag_expm1[state_idx] =
            fp_expm1(rate_eigenvalues[state_idx] * rate_scale * branch_length);
    }

    for (int row = 0; row < state_count; ++row) {
        for (int col = 0; col < state_count; ++col) {
            fp_t acc = fp_t(0);
            for (int eig_idx = 0; eig_idx < state_count; ++eig_idx) {
                const fp_t v_entry = V[row * state_count + eig_idx];
                const fp_t vinv_entry = Vinv[eig_idx * state_count + col];
                acc += v_entry * diag_expm1[eig_idx] * vinv_entry;
            }
            out_pmat[row * state_count + col] = acc;
        }
        out_pmat[row * state_count + row] += fp_t(1.0);
    }

    if (pinv > 0.0) {
        for (int row = 0; row < state_count; ++row) {
            for (int col = 0; col < state_count; ++col) {
                const fp_t identity_entry = (row == col) ? fp_t(1.0) : fp_t(0.0);
                out_pmat[row * state_count + col] =
                    (fp_t(1.0) - pinv) * out_pmat[row * state_count + col] +
                    pinv * identity_entry;
            }
        }
    }

    // Keep PMAT on the raw matrix-exponential path so branch-length
    // derivatives remain consistent with the eigenbasis diagtable used by the
    // Newton kernels. Small floating deviations are preferable to optimizing
    // against a renormalized matrix that the derivative code does not model.
    for (int row = 0; row < state_count; ++row) {
        for (int col = 0; col < state_count; ++col) {
            fp_t& entry = out_pmat[row * state_count + col];
            if (entry < 0.0 && entry > fp_t(-1e-12)) entry = 0.0;
        }
    }
}

void build_all_branch_pmats_gpu(
    int node_count,
    int state_count,
    int rate_categories,
    const fp_t* d_branch_lengths,
    const fp_t* d_eigenvectors,
    const fp_t* d_inverse_eigenvectors,
    const fp_t* d_rate_eigenvalues,
    fp_t* d_pmats,
    cudaStream_t stream)
{
    validate_pmat_inputs(
        state_count, rate_categories, d_branch_lengths, d_eigenvectors,
        d_inverse_eigenvectors, d_rate_eigenvalues, d_pmats);
    if (node_count <= 0) {
        throw std::invalid_argument("PMAT node count must be positive.");
    }
    constexpr int kBlockSize = 256;
    const size_t item_count = mlipper::util::checked_product(
        "all-branch PMAT item count", node_count, rate_categories);
    if (item_count > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("All-branch PMAT item count exceeds kernel indexing.");
    }
    const size_t grid_size =
        (item_count + kBlockSize - 1) / kBlockSize;
    if (grid_size > static_cast<size_t>(std::numeric_limits<unsigned>::max())) {
        throw std::length_error("All-branch PMAT launch dimensions exceed CUDA limits.");
    }
    build_all_branch_pmats_kernel<<<grid_size, kBlockSize, 0, stream>>>(
        node_count, state_count, rate_categories,
        d_branch_lengths, d_eigenvectors, d_inverse_eigenvectors,
        d_rate_eigenvalues, d_pmats);
}

void build_single_branch_pmat_gpu(
    int node_id,
    int state_count,
    int rate_categories,
    const fp_t* d_branch_lengths,
    const fp_t* d_eigenvectors,
    const fp_t* d_inverse_eigenvectors,
    const fp_t* d_rate_eigenvalues,
    fp_t* d_pmats,
    cudaStream_t stream)
{
    validate_pmat_inputs(
        state_count, rate_categories, d_branch_lengths, d_eigenvectors,
        d_inverse_eigenvectors, d_rate_eigenvalues, d_pmats);
    if (node_id < 0) {
        throw std::invalid_argument("PMAT node ID must be non-negative.");
    }
    constexpr int kBlockSize = 128;
    const size_t diagonal_count = mlipper::util::checked_product(
        "single-branch PMAT diagonal count", rate_categories, state_count);
    const size_t shared_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        diagonal_count, "single-branch PMAT shared memory");
    build_single_branch_pmat_kernel<<<1, kBlockSize, shared_bytes, stream>>>(
        node_id, state_count, rate_categories, d_branch_lengths,
        d_eigenvectors, d_inverse_eigenvectors, d_rate_eigenvalues, d_pmats);
}
