#include <Eigen/Dense>
#include <vector>

#include "RecoLocalCalo/Common/interface/nnls.h"

// #define DEBUG_NNLS_  CPU

#ifdef DEBUG_NNLS_CPU
#include <iostream>
#endif

using namespace std;
using namespace Eigen;

void nnls(const FixedMatrix& A,
                  const FixedVector& b,
                  FixedVector &x,
                  const double eps,
                  const unsigned int max_iterations
                 ) {
  // Fast NNLS (fnnls) algorithm as per
  // http://users.wfu.edu/plemmons/papers/Chennnonneg.pdf
  // page 8

  // pseudoinverse (A^T * A)^-1 * A^T
  // this pseudo-inverse has numerical issues
  // in order to avoid that I substitued the pseudoinvese wiht the QR
  // decomposition

  std::vector<unsigned int> P;
  std::vector<unsigned int> R(VECTOR_SIZE);

// initial set of indexes
  for (unsigned int i = 0; i < VECTOR_SIZE; ++i)
    R[i] = i;
  
  // main loop
  for (unsigned int iter = 0; iter < max_iterations; ++iter) {
#ifdef DEBUG_NNLS_CPU
// cout << "iter " << iter << endl;
#endif

    // NNLS
    // initialize the cost vector
    FixedVector w = A.transpose() * (b - (A * x));

#ifdef DEBUG_NNLS_CPU
// cout << "w" << endl << w << endl;
#endif

    // initialize the value for the while guard
    // max_index will contain the index of the max coeff anf max_w is the max
    // coeff
    unsigned int max_index = R[0];
    unsigned int remove_index = 0;

    for (unsigned int i = 0; i < R.size(); ++i) {
      auto index = R[i];
      if (w[index] > w[max_index]) {
        max_index = index;
        remove_index = i;
      }
    }

#ifdef DEBUG_NNLS_CPU
    // cout << "max index " << max_index << endl;
#endif

    P.emplace_back(max_index);
    R.erase(R.begin() + remove_index);

    // termination condition
    if (R.empty() || w[max_index] < eps)
      break;

#ifdef DEBUG_NNLS_CPU
      // cout << "P " << endl;
      // for (auto elem : P)
      //   cout << elem << " ";
      // cout << endl;
      // cout << "R " << endl;
      // for (auto elem : R)
      //   cout << elem << " ";
      // cout << endl;
#endif

    FixedMatrix A_P = FixedMatrix::Zero();

    for (auto index : P)
      A_P.col(index) = A.col(index);

      // #if DECOMPOSITION == USE_SPARSE_QR
      //     solver.compute(A_P.sparseView());
      // #else
      //     solver.compute(A_P);
      // #endif

#ifdef DEBUG_NNLS_CPU
      // cout << "A_P " << endl << A_P << endl;
#endif

      // FixedVector s = (A_P.transpose()*A_P).inverse() * A_P.transpose() * b;

#if DECOMPOSITION == USE_LLT
    FixedVector s = A_P.llt().matrixL().solve(b);
#elif DECOMPOSITION == USE_LDLT
    FixedVector s = A_P.ldlt().matrixL().solve(b);
#elif DECOMPOSITION == USE_HOUSEHOLDER
    FixedVector s = A_P.colPivHouseholderQr().solve(b);
#endif

    for (auto index : R)
      s[index] = 0;

#ifdef DEBUG_NNLS_CPU
      // cout << "s" << endl << s << endl;
#endif

    // inner loop
    while (true) {
      auto min_s = std::numeric_limits<double>::max();

      for (auto index : P)
        min_s = std::min(s[index], min_s);

#ifdef DEBUG_NNLS_CPU
        // cout << "min_s " << min_s << endl;
#endif

      if (min_s > 0)
        break;

#ifdef DEBUG_NNLS_CPU
      cout << "s" << endl << s << endl;
#endif
      auto alpha = std::numeric_limits<double>::max();
      // auto rm_idx = -1;

      for (unsigned int i = 0; i < P.size(); i++) {
        auto index = P[i];
        if (s[index] <= 0) {
          alpha = std::min(-x[index] / (s[index] - x[index]), alpha);
          // rm_idx = i;
        }
      }
#ifdef DEBUG_NNLS_CPU

      cout << "alpha " << alpha << endl;

      cout << "x before" << endl << x << endl;

#endif

      for (auto index : P) {
#ifdef DEBUG_NNLS_CPU
        cout << "index P " << index << endl;
        cout << "delta " << alpha * (s[index] - x[index]) << endl;
#endif

        x[index] += alpha * (s[index] - x[index]);
      }

#ifdef DEBUG_NNLS_CPU
      cout << "x after" << endl << x << endl;
#endif
      // R.emplace_back(P[rm_idx]);
      // P.erase(P.begin() + rm_idx);
      std::vector<unsigned int> tmp;

#ifdef DEBUG_NNLS_CPU
      cout << "P  before" << endl;
      for (auto elem : P)
        cout << elem << " ";
      cout << endl;
      cout << "R before" << endl;
      for (auto elem : R)
        cout << elem << " ";
      cout << endl;
#endif

      for (int i = P.size() - 1; i >= 0; --i) {
        auto index = P[i];
        if (x[index] == 0) {
          R.emplace_back(index);
          tmp.emplace_back(i);
        }
      }

      for (auto index : tmp)
        P.erase(P.begin() + index);

#ifdef DEBUG_NNLS_CPU
      cout << "P  after" << endl;
      for (auto elem : P)
        cout << elem << " ";
      cout << endl;
      cout << "R after" << endl;
      for (auto elem : R)
        cout << elem << " ";
      cout << endl;
#endif


      A_P.setZero();

      for (auto index : P)
        A_P.col(index) = A.col(index);

#if DECOMPOSITION == USE_LLT
      s = A_P.llt().matrixL().solve(b);
#elif DECOMPOSITION == USE_LDLT
      s = A_P.ldlt().matrixL().solve(b);
#elif DECOMPOSITION == USE_HOUSEHOLDER
      s = A_P.colPivHouseholderQr().solve(b);

#endif

      for (auto index : R)
        s[index] = 0;

#ifdef DEBUG_NNLS_CPU
      cout << "s after " << s << endl;
#endif
    }
    x = s;
  }
}
