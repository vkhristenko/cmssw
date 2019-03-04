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
  // Fast NNLS (fnnls) algorithm as per
  // http://users.wfu.edu/plemmons/papers/Chennnonneg.pdf
  // page 8

  // FNNLS memorizes the A^T * A and A^T * b to reduce the computation.
  // The pseudo-inverse obtained has the same numerical problems so
  // I keep the same decomposition utilized for NNLS.

  // pseudoinverse (A^T * A)^-1 * A^T
  // this pseudo-inverse has numerical issues
  // in order to avoid that I substituted the pseudoinverse whit the QR
  // decomposition

  // I'm substituting the vectors P and R that represents active and passive set
  // with a boolean vector: if active_set[i] the i belogns to R else to P

  // bool active_set[VECTOR_SIZE];
  // memset(active_set, true, VECTOR_SIZE * sizeof(bool));

//  x = vector_t<T>::Zero();

  /*
  #ifdef __CUDA_ARCH__
  FixedMatrix AtA;
  matrixMultiplication(A, AtA);
  __syncthreads();
  cudaDeviceSynchronize();
  for(auto i = 0; i < MATRIX_SIZE; ++i){
    for(auto j = 0; j < MATRIX_SIZE; ++j)
      printf("%f ", AtA.data()[j* MATRIX_SIZE + i]);
    printf("\n");
  }
  #else
  FixedMatrix AtA = transpose_multiply(A);
  #endif
  assert(AtA == A.transpose() * A);
  */
//  matrix_t<data_type> AtA = transpose_multiply(A);
#ifdef NNLS_DEBUG
  std::cout << "A = \n" << A << std::endl;
#endif
//  matrix_t<data_type> AtA = transpose_multiply(A);
  matrix_t AtA = A.transpose() * A;
  vector_t Atb = A.transpose() *b;

#ifndef __CUDA_ARCH__
#ifdef FNNLS_DEBUG_CPU
    static int __counter__ = 0;
    if (__counter__ == 113) {
    std::cout << "*** AtA ***" << std::endl;
    std::cout << AtA << std::endl;
    std::cout << "*** Atb ***" << std::endl;
    std::cout << Atb << std::endl;
    }
#endif
#endif

//  matrix_t<data_type> AtA = A.transpose() * A;
  // FixedMatrix AtA = A.transpose() * A;
  vector_t s;
  vector_t w;

  Eigen::PermutationMatrix<vector_t::RowsAtCompileTime> permutation;
  permutation.setIdentity();

#ifdef NNLS_DEBUG
  std::cout << "AtA = \n" << AtA << std::endl;
  std::cout << "Atb = \n" << Atb << std::endl;
#endif

// main loop
  Eigen::Index w_max_idx_prev = 0;
  float w_max_prev = 0;
  for (auto iter = 0; iter < max_iterations; ++iter) {
    const auto nActive = vector_t::RowsAtCompileTime - npassive;

#ifdef NNLS_DEBUG
    std::cout << "***************\n"
        << "iteration = " << iter << std::endl
        << "nactive = " << nActive << std::endl;
    std::cout << "x = \n" << x << std::endl;
#endif

#ifdef DEBUG_FNNLS_CPU
    cout << "iter " << iter << endl;
#endif
    
    if(!nActive)
      break;

#ifdef NNLS_DEBUG
    std::cout << "AtA * x = \n" << AtA*x << std::endl;
#endif

    w.tail(nActive) = Atb.tail(nActive) - (AtA * x).tail(nActive);

#ifdef DEBUG_FNNLS_CPU
    cout << "w" << endl << w.tail(nActive) << endl;
#endif
    // get the index of w that gives the maximum gain
    Eigen::Index w_max_idx;
    const auto max_w = w.tail(nActive).maxCoeff(&w_max_idx);

#ifdef NNLS_DEBUG
    std::cout << "w = \n" << w << std::endl;
    std::cout << "max_w = " << max_w << std::endl;
    std::cout << "w_max_idx = " << w_max_idx << std::endl;
#endif

    // check for convergence
    if (max_w < eps || w_max_idx==w_max_idx_prev && max_w==w_max_prev)
      break;

    w_max_prev = max_w;
    w_max_idx_prev = w_max_idx;

    // cout << "n active " << nActive << endl;
    // cout << "w max idx " << w_max_idx << endl;

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

#ifdef NNLS_DEBUG
    std::cout << "permutation = \n" << permutation.indices() << std::endl;
#endif

    ++npassive;

#ifdef DEBUG_FNNLS_CPU
    cout << "max index " << w_max_idx << endl;
    std::cout << "n_active " << nActive << std::endl;
#endif

// inner loop
    while (true) {
      if (npassive == 0) break;

      s.head(npassive) =
          AtA.topLeftCorner(npassive, npassive).llt().solve(Atb.head(npassive));

#ifdef NNLS_DEBUG
      std::cout << "s = \n" << s << std::endl;
#endif

      if (s.head(npassive).minCoeff() > 0.) {
        x.head(npassive) = s.head(npassive);
        break;
      }

#ifdef DEBUG_FNNLS_CPU
      cout << "s" << endl << s.head(npassive) << endl;
#endif

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

#ifdef DEBUG_FNNLS_CPU

      cout << "alpha " << alpha << endl;

      cout << "x before" << endl << x << endl;

#endif

      x.head(npassive) += alpha * (s.head(npassive) - x.head(npassive));
      x[alpha_idx] = 0;
      --npassive;

#ifdef DEBUG_FNNLS_CPU
      cout << "x after" << endl << x << endl;
#endif
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
