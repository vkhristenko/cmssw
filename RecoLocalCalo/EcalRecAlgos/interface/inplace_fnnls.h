#ifndef inplace_fnnls_hpp
#define inpalce_fnnls_hpp

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"

namespace ecal { namespace multifit { namespace v1 {

using matrix_t = SampleMatrix;
using vector_t = SampleVector;

__device__
bool
inplace_fnnls(matrix_t const& A,
              vector_t const& b,
              vector_t& x,
              int& npassive,
              const double eps = 1e-11,
              const unsigned int max_iterations = 1000);

__device__
bool inplace_fnnls(matrix_t const& A,
                   vector_t const& b,
                   vector_t& x,
                   int& npassive,
                   const double eps,
                   const unsigned int max_iterations) {
  matrix_t AtA = A.transpose() * A;
  vector_t Atb = A.transpose() *b;
  vector_t s;
  vector_t w;

  Eigen::PermutationMatrix<vector_t::RowsAtCompileTime> permutation;
  permutation.setIdentity();

// main loop
  Eigen::Index w_max_idx_prev = 0;
  float w_max_prev = 0;
  for (auto iter = 0; iter < max_iterations; ++iter) {
    const auto nActive = vector_t::RowsAtCompileTime - npassive;

    if(!nActive)
      break;

    w.tail(nActive) = Atb.tail(nActive) - (AtA * x).tail(nActive);

    // get the index of w that gives the maximum gain
    Eigen::Index w_max_idx;
    const auto max_w = w.tail(nActive).maxCoeff(&w_max_idx);

    // check for convergence
    if (max_w < eps || w_max_idx==w_max_idx_prev && max_w==w_max_prev)
      break;

    w_max_prev = max_w;
    w_max_idx_prev = w_max_idx;

    // need to translate the index into the right part of the vector
    w_max_idx += npassive;

    // swap AtA to avoid copy
    AtA.col(npassive).swap(AtA.col(w_max_idx));
    AtA.row(npassive).swap(AtA.row(w_max_idx));
    // swap Atb to match with AtA
    Eigen::numext::swap(Atb.coeffRef(npassive), Atb.coeffRef(w_max_idx));
    Eigen::numext::swap(x.coeffRef(npassive), x.coeffRef(w_max_idx));
    // swap the permutation matrix to reorder the solution in the end
    Eigen::numext::swap(permutation.indices()[npassive],
                        permutation.indices()[w_max_idx]);

    ++npassive;

// inner loop
    while (true) {
      if (npassive == 0) break;

      s.head(npassive) =
          AtA.topLeftCorner(npassive, npassive).llt().solve(Atb.head(npassive));

      if (s.head(npassive).minCoeff() > 0.) {
        x.head(npassive) = s.head(npassive);
        break;
      }

      auto alpha = std::numeric_limits<double>::max();
      Eigen::Index alpha_idx = 0;

#pragma unroll 10
      for (auto i = 0; i < npassive; ++i) {
        if (s[i] <= 0.) {
          auto const ratio = x[i] / (x[i] - s[i]);
          if (ratio < alpha) {
            alpha = ratio;
            alpha_idx = i;
          }
        }
      }

      if (std::numeric_limits<double>::max() == alpha) {
        x.head(npassive) = s.head(npassive);
        break;
      }

      x.head(npassive) += alpha * (s.head(npassive) - x.head(npassive));
      x[alpha_idx] = 0;
      --npassive;

      AtA.col(npassive).swap(AtA.col(alpha_idx));
      AtA.row(npassive).swap(AtA.row(alpha_idx));
      // swap Atb to match with AtA
      Eigen::numext::swap(Atb.coeffRef(npassive), Atb.coeffRef(alpha_idx));
      Eigen::numext::swap(x.coeffRef(npassive), x.coeffRef(alpha_idx));
      // swap the permutation matrix to reorder the solution in the end
      Eigen::numext::swap(permutation.indices()[npassive],
                          permutation.indices()[alpha_idx]);

    }
  }
  
  x = x.transpose() * permutation.transpose();  
  return true;
}

}}}

#endif
