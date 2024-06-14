#include <cppdlr/cppdlr.hpp>
#include "nda/lapack/geqp3.hpp"
#include "../src/dlr2d.hpp"
#include "nda/layout/policies.hpp"

int main() {

    int n = 1000; // Matrix size
    double eps = 1e-8;

    double alpha = 2;
    int nvec = 100;

    double rate = 1.05;

    // Random matrix with decaying singular values

    // Decaying singular values
    auto svals = nda::matrix<double>(n, n);
    for (int i = 0; i < n; i++) {
        svals(i, i) = pow(rate,-i);
    }

    // Random orthogonal matrices
    auto randmat = nda::array<dcomplex, 2, F_layout>::rand(std::array{n, n});
    auto u = nda::matrix<dcomplex, F_layout>(n, n);
    auto vt = nda::matrix<dcomplex, F_layout>(n, n);
    auto s = nda::vector<double>(n);
    nda::lapack::gesvd(randmat, s, u, vt);

    auto a = nda::matrix<dcomplex, F_layout>(n, n);
    a = u * svals * vt;

    // Pivoted QR decomposition
    auto piv = nda::zeros<int>(n);
    auto tau = nda::vector<dcomplex>(n);
    nda::lapack::geqp3(a, piv, tau);

    // Estimate rank the old-fashioned way
    int rank1 = 0;
    for (int k = 0; k < n; ++k) {
      if (abs(a(k, k)) < eps) {
        rank1 = k;
        break;
      }
    }

    // Estimate rank using randomized algorithm
    int rank2 = estimate_rank(a, eps, alpha, nvec);
    
    PRINT(rank1);
    PRINT(rank2);
    PRINT(log(1/eps)/log(rate));

    // Compute the true spectral norm of R22 for each rank

    auto r22_1 = nda::matrix<dcomplex, F_layout>(n-rank1, n-rank1);
    auto r22_2 = nda::matrix<dcomplex, F_layout>(n-rank2, n-rank2);
    r22_1 = 0;
    r22_2 = 0;
    for (int i = rank1+1; i < n; i++) {
        for (int j = i; j < n; j++) {
            r22_1(i-rank1-1, j-rank1-1) = a(i, j);
        }
    }
    for (int i = rank2+1; i < n; i++) {
        for (int j = i; j < n; j++) {
            r22_2(i-rank2-1, j-rank2-1) = a(i, j);
        }
    }

    auto u1 = nda::matrix<dcomplex, F_layout>(n-rank1, n-rank1);
    auto vt1 = nda::matrix<dcomplex, F_layout>(n-rank1, n-rank1);
    auto s1 = nda::vector<double>(n-rank1);
    nda::lapack::gesvd(r22_1, s1, u1, vt1);

    PRINT(s1(0));

    auto u2 = nda::matrix<dcomplex, F_layout>(n-rank2, n-rank2);
    auto vt2 = nda::matrix<dcomplex, F_layout>(n-rank2, n-rank2);
    auto s2 = nda::vector<double>(n-rank2);
    nda::lapack::gesvd(r22_2, s2, u2, vt2);

    PRINT(s2(0));




}