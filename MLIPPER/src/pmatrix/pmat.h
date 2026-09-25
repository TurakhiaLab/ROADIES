#pragma once

#include <vector>

struct EigResult {
    // Row-major matrices satisfying Q = V * diag(lambdas) * Vinv.
    std::vector<double> lambdas;
    std::vector<double> V;
    std::vector<double> Vinv;
};

// Eigendecompose a reversible rate matrix using pi for the symmetric transform.
EigResult gtr_eigendecomp_cpu(
    const double* Q_rowmajor,
    const double* pi,
    int n);

// All matrices are row-major n-by-n arrays. V and Vinv must satisfy
// Q = V * diag(lambdas) * Vinv; pi is positive and sums to one.
void pmatrix_from_triple(
    const double* Vinv,
    const double* V,
    const double* lambdas,
    double rate_scale,
    double branch_length,
    double pinv,
    double* output,
    int n);
