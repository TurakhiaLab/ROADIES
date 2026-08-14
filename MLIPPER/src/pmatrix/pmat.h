#pragma once
#include <vector>

struct EigResult {
    std::vector<double> lambdas;
    std::vector<double> V;
    std::vector<double> Vinv;
};

EigResult gtr_eigendecomp_cpu(
    const double* Q_rowmajor,   // Q (row-major, size n*n)
    const double* pi,           // Stationary frequencies (length n, sum=1)
    int n);

void pmatrix_from_triple(const double* Vinv, const double* V,
                                const double* lamb, double r, double t, double p,
                                double* P, int n);
