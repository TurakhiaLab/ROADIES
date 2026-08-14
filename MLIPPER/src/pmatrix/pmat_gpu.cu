#include "pmat_gpu.cuh"

#include <cmath>

namespace {
constexpr int kMaxStates = 16;
}

__device__ void pmatrix_from_triple_device(
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
