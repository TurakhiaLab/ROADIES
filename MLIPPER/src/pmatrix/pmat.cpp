#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "pmat.h"

#include "util/checked_size.hpp"

extern "C" {
void dsyevd_(char* jobz, char* uplo, int* n,
             double* a, int* lda,
             double* w,
             double* work, int* lwork,
             int* iwork, int* liwork,
             int* info);
}

namespace {

// Remove only roundoff-scale negative probabilities, then repair row sums.
// Materially negative values are retained so upstream numerical defects are
// not silently converted into plausible transition matrices.
void clamp_neg_and_row_norm(double* p, int n) {
    for (int row = 0; row < n; ++row) {
        double sum = 0.0;
        for (int col = 0; col < n; ++col) {
            double& entry = p[row * n + col];
            if (entry < 0.0 && entry > -1e-14) entry = 0.0;
            sum += entry;
        }
        if (sum != 0.0) {
            for (int col = 0; col < n; ++col) p[row * n + col] /= sum;
        }
    }
}

}  // namespace

void pmatrix_from_triple(
    const double* Vinv,
    const double* V,
    const double* lamb,
    double r,
    double t,
    double p,
    double* P,
    int n)
{
    if (!Vinv || !V || !lamb || !P || n <= 0) {
        throw std::invalid_argument("Invalid PMAT construction arguments.");
    }
    const size_t matrix_elems = mlipper::util::checked_product(
        "CPU PMAT matrix elements", n, n);
    std::vector<double> I(matrix_elems, 0.0);
    for (int i = 0; i < n; ++i) I[i * n + i] = 1.0;

    // Compute V*(exp(Lrt)-I)*V^-1 and add I afterward. expm1 preserves small
    // branch changes that would lose precision in exp(x)-1.
    std::vector<double> D(matrix_elems, 0.0);
    for (int j = 0; j < n; ++j) D[j * n + j] = std::expm1(lamb[j] * r * t);

    std::vector<double> T(matrix_elems, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            const double vik = V[i * n + k];
            for (int j = 0; j < n; ++j) T[i * n + j] += vik * D[k * n + j];
        }
    }

    std::fill(P, P + matrix_elems, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            const double tik = T[i * n + k];
            for (int j = 0; j < n; ++j) P[i * n + j] += tik * Vinv[k * n + j];
        }
    }
    for (int i = 0; i < n; ++i) P[i * n + i] += 1.0;

    if (p > 0.0) {
        for (size_t i = 0; i < matrix_elems; ++i) {
            P[i] = (1.0 - p) * P[i] + p * I[i];
        }
    }

    clamp_neg_and_row_norm(P, n);
}

EigResult gtr_eigendecomp_cpu(const double* Q_rowmajor, const double* pi, int n) {
    if (!Q_rowmajor || !pi || n <= 0) {
        throw std::invalid_argument("Invalid eigendecomposition arguments.");
    }
    const size_t matrix_elems = mlipper::util::checked_product(
        "GTR eigendecomposition matrix elements", n, n);

    std::vector<double> sqrtpi(n);
    for (int i = 0; i < n; ++i) {
        if (pi[i] <= 0.0) throw std::invalid_argument("pi must be > 0");
        sqrtpi[i] = std::sqrt(pi[i]);
    }

    // Reversibility makes D^(1/2) Q D^(-1/2) symmetric, allowing the stable
    // symmetric LAPACK solver while retaining a row-major public interface.
    std::vector<double> S(matrix_elems);
    for (int i = 0; i < n; ++i) {
        const double si = sqrtpi[i];
        for (int j = 0; j < n; ++j) S[i * n + j] = si * Q_rowmajor[i * n + j] / sqrtpi[j];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double sym = 0.5 * (S[i * n + j] + S[j * n + i]);
            S[i * n + j] = sym;
            S[j * n + i] = sym;
        }
    }

    std::vector<double> lambdas(n);
    // LAPACK is column-major, so transpose S on input and its eigenvectors on
    // output. U below is row-major with eigenvectors in columns.
    std::vector<double> A(matrix_elems);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) A[j * n + i] = S[i * n + j];
    }

    char jobz = 'V';
    char uplo = 'U';
    int N = n;
    int lda = n;
    int info = 0;
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;

    dsyevd_(&jobz, &uplo, &N,
            A.data(), &lda,
            lambdas.data(),
            &work_query, &lwork,
            &iwork_query, &liwork,
            &info);
    if (info != 0) throw std::runtime_error("dsyevd workspace query failed");

    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);

    dsyevd_(&jobz, &uplo, &N,
            A.data(), &lda,
            lambdas.data(),
            work.data(), &lwork,
            iwork.data(), &liwork,
            &info);
    if (info != 0) throw std::runtime_error("dsyevd failed, info=" + std::to_string(info));

    std::vector<double> U(matrix_elems);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            U[i * n + j] = A[j * n + i];
        }
    }

    EigResult out;
    out.lambdas = lambdas;
    out.V.resize(matrix_elems);
    out.Vinv.resize(matrix_elems);

    for (int i = 0; i < n; ++i) {
        const double invsqrt = 1.0 / sqrtpi[i];
        for (int j = 0; j < n; ++j) out.V[i * n + j] = U[i * n + j] * invsqrt;
    }
    for (int i = 0; i < n; ++i) {
        const double sp = sqrtpi[i];
        for (int j = 0; j < n; ++j) out.Vinv[j * n + i] = U[i * n + j] * sp;
    }

    for (double& lambda : out.lambdas) {
        if (lambda > -1e-15 && lambda < 1e-15) lambda = 0.0;
    }

    return out;
}
