#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "pmat.h"

extern "C" {
void dsyevd_(char* jobz, char* uplo, int* n,
             double* a, int* lda,
             double* w,
             double* work, int* lwork,
             int* iwork, int* liwork,
             int* info);
}

namespace {

void eye(double* matrix, int n) {
    std::fill(matrix, matrix + n * n, 0.0);
    for (int i = 0; i < n; ++i) matrix[i * n + i] = 1.0;
}

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

void pmatrix_from_triple(const double* Vinv, const double* V,
                         const double* lamb, double r, double t, double p,
                         double* P, int n) {
    std::vector<double> I(n * n, 0.0);
    for (int i = 0; i < n; ++i) I[i * n + i] = 1.0;

    std::vector<double> D(n * n, 0.0);
    for (int j = 0; j < n; ++j) D[j * n + j] = std::expm1(lamb[j] * r * t);

    std::vector<double> T(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            const double vik = V[i * n + k];
            for (int j = 0; j < n; ++j) T[i * n + j] += vik * D[k * n + j];
        }
    }

    std::fill(P, P + n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            const double tik = T[i * n + k];
            for (int j = 0; j < n; ++j) P[i * n + j] += tik * Vinv[k * n + j];
        }
    }
    for (int i = 0; i < n; ++i) P[i * n + i] += 1.0;

    if (p > 0.0) {
        for (int i = 0; i < n * n; ++i) P[i] = (1.0 - p) * P[i] + p * I[i];
    }

    clamp_neg_and_row_norm(P, n);
}

EigResult gtr_eigendecomp_cpu(const double* Q_rowmajor, const double* pi, int n) {
    if (n <= 0) throw std::invalid_argument("n must be > 0");

    std::vector<double> sqrtpi(n);
    for (int i = 0; i < n; ++i) {
        if (pi[i] <= 0.0) throw std::invalid_argument("pi must be > 0");
        sqrtpi[i] = std::sqrt(pi[i]);
    }

    std::vector<double> S(n * n);
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
    std::vector<double> A(n * n);
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

    std::vector<double> U(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            U[i * n + j] = A[j * n + i];

    EigResult out;
    out.lambdas = lambdas;
    out.V.resize(n * n);
    out.Vinv.resize(n * n);

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
