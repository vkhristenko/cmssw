#include <Eigen/Dense>
#include "RecoLocalCalo/Common/interface/fnnls.h"

// #define DEBUG_FNNLS_CPU

#ifdef DEBUG_FNNLS_CPU
#include <iostream>
using namespace std;
void print_active_set(bool v[VECTOR_SIZE], std::string s = "") {
  std::cout << s;
  for (size_t i = 0; i < VECTOR_SIZE; i++) {
    std::cout << i << ' ';
  }
  std::cout << std::endl;

  for (size_t i = 0; i < VECTOR_SIZE; i++) {
    std::cout << v[i] << ' ';
  }
  std::cout << std::endl;
}
#endif

using namespace Eigen;

void
fnnls(const FixedMatrix& A,
      const FixedVector& b,
      FixedVector& x,
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

  bool active_set[VECTOR_SIZE];
  memset(active_set, true, VECTOR_SIZE * sizeof(bool));

#ifdef DEBUG_FNNLS_CPU
  print_active_set(active_set);
#endif

  auto AtA = A.transpose() * A;
  auto Atb = A.transpose() * b;

  Eigen::Matrix<double, -1, -1, 0, 10, 10> AtA_P;
  Eigen::Matrix<double, -1, 1, 0, 10, 1> Atb_P;
  Eigen::Matrix<double, -1, 1, 0, 10, 1> tmp_s;

  FixedVector s;
  FixedVector w;

  // main loop
  for (unsigned int iter = 0; iter < max_iterations; ++iter) {
#ifdef DEBUG_FNNLS_CPU
    cout << "iter " << iter << endl;
#endif

    w = Atb - (AtA * x);

#ifdef DEBUG_FNNLS_CPU
    cout << "w" << endl << w << endl;
#endif

    // initialize the value for the while guard
    // max_index will contain the index of the max coeff anf max_w is the max
    // coeff
    auto max_index = 0;
    // start from -1 just because after the next loop one element will be set to
    // inactive
    auto n_active = -1;
    // unsigned int remove_index = 0;

    for (unsigned int i = 0; i < VECTOR_SIZE; ++i) {
      if (active_set[i]) {
        ++n_active;
        if (w[i] > w[max_index])
          max_index = i;
      }
    }

#ifdef DEBUG_FNNLS_CPU
    print_active_set(active_set);
    cout << "max index " << max_index << endl;
#endif

    // setting one element to inactive not decrementing n_active because I
    // already took it into account
    active_set[max_index] = false;

#ifdef DEBUG_FNNLS_CPU
    print_active_set(active_set);
    std::cout << "n_active " << n_active << std::endl;
#endif

    // termination condition
    // if n_active == 0 then exit
    if (!n_active || w[max_index] < eps)
      break;
    //

    // inner loop
    while (true) {
      // cast is needed, otherwise cxx complains signed/unsigned comparisons
      if (n_active >= static_cast<int>(VECTOR_SIZE)) {
        x = s;
        break;
      }
      //
      AtA_P = sub_matrix<FixedMatrix, bool[VECTOR_SIZE]>(
          AtA, active_set, VECTOR_SIZE - n_active);
      Atb_P = sub_vector<FixedVector, bool[VECTOR_SIZE]>(
          Atb, active_set, VECTOR_SIZE - n_active);
      tmp_s = AtA_P.llt().solve(Atb_P);
      //
      for (unsigned int i = 0, idx = 0; i < VECTOR_SIZE; i++) {
        s[i] = active_set[i] ? 0 : tmp_s[idx++];
      }

#ifdef DEBUG_FNNLS_CPU
      cout << "s" << endl << s << endl;
#endif
      auto alpha = std::numeric_limits<double>::max();
      Index alpha_idx = -1;
      for (unsigned int i = 0; i < VECTOR_SIZE; i++) {
        if (!active_set[i] && s[i] <= 0.) {
          auto const tmp = -x[i] / (s[i] - x[i]);
          if (tmp < alpha) {
            alpha = tmp;
            alpha_idx = i;
          }
        }
      }
      if (std::numeric_limits<double>::max() == alpha) {
        x = s;
        break;
      }

#ifdef DEBUG_FNNLS_CPU

      cout << "alpha " << alpha << endl;

      cout << "x before" << endl << x << endl;

#endif

      x += alpha * (s - x);
      x[alpha_idx] = 0;
      active_set[alpha_idx] = true;
      ++n_active;

#ifdef DEBUG_FNNLS_CPU
      cout << "x after" << endl << x << endl;
#endif
    }
  }
}
