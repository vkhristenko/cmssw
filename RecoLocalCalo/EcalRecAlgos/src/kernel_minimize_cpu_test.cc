#include "RecoLocalCalo/EcalRecAlgos/interface/kernel_minimize_cpu_test.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"

#include <iostream>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "cuda.h"

//#define DEBUG

namespace ecal { namespace multifit { namespace v1 {

namespace cpu {

void eigen_solve_submatrix(SampleMatrix& mat, 
                           SampleVector& invec, 
                           SampleVector& outvec, unsigned NP) {
    using namespace Eigen;
    switch( NP ) { // pulse matrix is always square.
    case 10: {   
        Matrix<SampleMatrix::Scalar,10,10> temp = mat.topLeftCorner<10,10>();
        outvec.head<10>() = temp.ldlt().solve(invec.head<10>());
        break;
    }   
    case 9: {
        Matrix<SampleMatrix::Scalar,9,9> temp = mat.topLeftCorner<9,9>();
        outvec.head<9>() = temp.ldlt().solve(invec.head<9>());
        break;
    }   
    case 8: {   
        Matrix<SampleMatrix::Scalar,8,8> temp = mat.topLeftCorner<8,8>();
        outvec.head<8>() = temp.ldlt().solve(invec.head<8>());
        break;
    }   
    case 7: {   
        Matrix<SampleMatrix::Scalar,7,7> temp = mat.topLeftCorner<7,7>();
        outvec.head<7>() = temp.ldlt().solve(invec.head<7>());
        break;
    }
    case 6: {   
        Matrix<SampleMatrix::Scalar,6,6> temp = mat.topLeftCorner<6,6>();
        outvec.head<6>() = temp.ldlt().solve(invec.head<6>());
        break;
    }
    case 5: {   
        Matrix<SampleMatrix::Scalar,5,5> temp = mat.topLeftCorner<5,5>();
        outvec.head<5>() = temp.ldlt().solve(invec.head<5>());
        break;
    }
    case 4: {   
        Matrix<SampleMatrix::Scalar,4,4> temp = mat.topLeftCorner<4,4>();
        outvec.head<4>() = temp.ldlt().solve(invec.head<4>());
        break;
    }
    case 3: {   
        Matrix<SampleMatrix::Scalar,3,3> temp = mat.topLeftCorner<3,3>();
        outvec.head<3>() = temp.ldlt().solve(invec.head<3>());
        break;
    }   
    case 2: {   
        Matrix<SampleMatrix::Scalar,2,2> temp = mat.topLeftCorner<2,2>();
        outvec.head<2>() = temp.ldlt().solve(invec.head<2>());
        break;
    }
    case 1: {   
        Matrix<SampleMatrix::Scalar,1,1> temp = mat.topLeftCorner<1,1>();
        outvec.head<1>() = temp.ldlt().solve(invec.head<1>());
        break;
    }                                                                              
    default:
        return;
    }
}


bool update_covariance(SampleMatrix const& noisecov,
                       FullSampleMatrix const& full_pulse_cov,
                       SampleMatrix& inverse_cov,
                       BXVectorType const& bxs,
                       SampleDecompLLT& covariance_decomposition,
                       SampleVector const& amplitudes) {
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = BXVectorType::RowsAtCompileTime;

    inverse_cov = noisecov;

#ifdef ECAL_RECO_DEBUG 
    std::cout << "*** inverse cov before ***\n";
    std::cout << inverse_cov << std::endl;

    std::cout << "*** amplitudes ***\n";
    std::cout << amplitudes << std::endl;

    std::cout << "*** bxs ***\n";
    std::cout << bxs << std::endl;
#endif

    for (unsigned int ipulse=0; ipulse<npulses; ipulse++) {
        if (amplitudes.coeff(ipulse) == 0) 
            continue;

        int bx = bxs.coeff(ipulse);
        int first_sample_t = std::max(0, bx+3);
        int offset = 7 - 3 - bx;

#ifdef ECAL_RECO_DEBUG

        std::cout << "*** ipulse = " << ipulse << " ***\n";
        std::cout << "*** bx = " << bx << " ***\n";
        std::cout << "*** firstsamplet = " << first_sample_t << " ***\n";
        std::cout << "*** offset = " << offset <<  " ***\n";
#endif

        auto const value = amplitudes.coeff(ipulse);
        auto const value_sq = value*value;

#ifdef ECAL_RECO_DEBUG

        std::cout << "*** ampveccoef = " << value << " ***\n";
        std::cout << "*** ampsq = " << value_sq << " ***\n";
#endif

        unsigned int nsample_pulse = nsamples - first_sample_t;

#ifdef ECAL_RECO_DEBUG

        std::cout << "*** nsamplepulse = " << nsample_pulse << " ***\n";
#endif

        inverse_cov.block(first_sample_t, first_sample_t, 
                          nsample_pulse, nsample_pulse)
            += value_sq * full_pulse_cov.block(first_sample_t + offset,
                                               first_sample_t + offset,
                                               nsample_pulse,
                                               nsample_pulse);
    }

#ifdef ECAL_RECO_DEBUG
    std::cout << "*** noise cov after ***\n";
    std::cout << inverse_cov << std::endl;
#endif

    covariance_decomposition.compute(inverse_cov);
    return true;
}

float compute_chi2(SampleDecompLLT& covariance_decomposition,
                   PulseMatrixType const& pulse_matrix,
                   SampleVector const& amplitudes,
                   SampleVector const& samples) {
    return covariance_decomposition.matrixL()
        .solve(pulse_matrix * amplitudes - samples)
        .squaredNorm();
}

void kernel_minimize(SampleMatrix const* noisecov,
                     FullSampleMatrix const* full_pulse_cov,
                     BXVectorType const* bxs,
                     SampleVector const* samples,
                     SampleVector* amplitudes,
                     float* energies,
                     PulseMatrixType* pulse_matrix, 
                     bool* statuses,
                     float* chi2s,
                     bool const* isSaturated,
                     bool const* hasSwitchToGain6,
                     bool const* hasSwitchToGain1,
                     int nchannels,
                     int max_iterations, 
                     bool gainSwitchUseMaxSample) {
    for (int idx=0; idx<nchannels; idx++) {

#ifdef ECAL_RECO_DEBUG
            std::cout << "*** samples ***\n";
            std::cout << samples[idx] << std::endl;

            std::cout << "*** pulse matrix ***\n";
            std::cout << pulse_matrix[idx] << std::endl;

            std::cout << "*** noise cov ***\n";
            std::cout << noisecov[idx] << std::endl;

            std::cout << "*** pulse cov ***\n";
            std::cout << full_pulse_cov[idx] << std::endl;
#endif

        bool hasGainSwitch = isSaturated[idx] 
            || hasSwitchToGain6[idx]
            || hasSwitchToGain1[idx];
        // TODO: gainSwitchUseMaxSimple depends on eb/ee
        // in principle can be splitted/removed into different kernels
        // for ee non-divergent branch
        if (hasGainSwitch && gainSwitchUseMaxSample)
            return;
        bool status = false;
        int iter = 0;
        SampleDecompLLT covariance_decomposition;
        SampleMatrix inverse_cov;
        int npassive = 0;
        amplitudes[idx] = SampleVector::Zero();
        float chi2 = 0;
        float chi2_now = 0;

        // TODO 
        BXVectorType activeBXs = *bxs;
        permutation_t permutation;
        permutation.setIdentity();

        // loop until ocnverge
        while (true) {
            if (iter >= max_iterations)
                break;

#ifdef ECAL_RECO_DEBUG
            std::cout << "iter = " << iter << std::endl;
#endif

            // TODO
            status = update_covariance(
                noisecov[idx], 
                full_pulse_cov[idx],
                inverse_cov,
                activeBXs,
                covariance_decomposition,
                amplitudes[idx]);
            if (!status) 
                break;

            // TODO
            SampleMatrix A = covariance_decomposition.matrixL()
                .solve(pulse_matrix[idx]);
            SampleVector b = covariance_decomposition.matrixL()
                .solve(samples[idx]);

#ifdef ECAL_RECO_DEBUG
                std::cout << "*** matrix A ***\n";
                std::cout << A.transpose() * A << std::endl;

                std::cout << "*** vector b ***\n";
                std::cout << A.transpose() * b << std::endl;
#endif
            
            status = inplace_fnnls(
                A, b, amplitudes[idx],
                npassive, activeBXs, permutation, pulse_matrix[idx]);

#ifdef ECAL_RECO_DEBUG
                std::cout << "*** solution ***\n";
                std::cout << amplitudes[idx] << std::endl;
#endif
                
            if (!status)
                break;

            // TODO
            chi2_now = compute_chi2(
                covariance_decomposition,
                pulse_matrix[idx],
                amplitudes[idx],
                samples[idx]);
            float deltachi2 = chi2_now - chi2;
            chi2 = chi2_now;
            ++iter;

            if (ecal::abs(deltachi2) < 1e-3)
                break;
        }

        amplitudes[idx] = amplitudes[idx].transpose() * permutation.transpose();
        float energy = amplitudes[idx](5);
        energies[idx] = energy; // according to bxs vector bxs[5] = 0
        statuses[idx] = status;
        chi2s[idx] = chi2;
    }
}

bool inplace_fnnls(matrix_t const& A,
                   vector_t const& b,
                   vector_t& x,
                   int& npassive,
                   BXVectorType& activeBXs,
                   permutation_t& permutation,
                   PulseMatrixType& pulse_matrix,
                   const double eps,
                   const unsigned int max_iterations) {
  matrix_t AtA;
  AtA.noalias() = A.transpose().lazyProduct(A);
  vector_t Atb;
  Atb.noalias() = A.transpose().lazyProduct(b);
  vector_t s;
  vector_t w;

// main loop
  Eigen::Index w_max_idx_prev = 0;
  matrix_t::Scalar w_max_prev = 0;
  double eps_to_use = eps;
  unsigned int iter = 0;
  while (true) {
//  for (unsigned int iter = 0; iter < max_iterations; ++iter) {
    if (iter>0 || npassive==0) {
        const auto nActive = vector_t::RowsAtCompileTime - npassive;

        if(!nActive)
          break;

    //    w.tail(nActive) = Atb.tail(nActive) - (AtA * x).tail(nActive);
        w = Atb - (AtA * x);

        // get the index of w that gives the maximum gain
        Eigen::Index w_max_idx;
        const auto max_w = w.tail(nActive).maxCoeff(&w_max_idx);

        // check for convergence
        if (max_w < eps_to_use || (w_max_idx==w_max_idx_prev && max_w==w_max_prev))
          break;

        // worst case
        if (iter >= 500)
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
        Eigen::numext::swap(activeBXs.coeffRef(npassive), activeBXs.coeffRef(w_max_idx));
        pulse_matrix.col(npassive).swap(pulse_matrix.col(w_max_idx));

        ++npassive;
    }

// inner loop
    while (true) {
      if (npassive == 0) break;

      vector_t s = x;
      eigen_solve_submatrix(AtA, Atb, s, npassive);
      /*
      s.head(npassive) =
          AtA.topLeftCorner(npassive, npassive).llt().solve(Atb.head(npassive));
          */

      // if all coefficients are positive, done for this iteration
      if (s.head(npassive).minCoeff() > 0.) {
        x.head(npassive) = s.head(npassive);
        break;
      }

      auto alpha = std::numeric_limits<matrix_t::Scalar>::max();
      Eigen::Index alpha_idx = 0;

      for (auto i = 0; i < npassive; ++i) {
        if (s[i] <= 0.) {
          auto const ratio = x[i] / (x[i] - s[i]);
          if (ratio < alpha) {
            alpha = ratio;
            alpha_idx = i;
          }
        }
      }

      /*
      if (std::numeric_limits<double>::max() == alpha) {
        x.head(npassive) = s.head(npassive);
        break;
      }*/

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
      Eigen::numext::swap(activeBXs.coeffRef(npassive), 
                          activeBXs.coeffRef(alpha_idx));
      pulse_matrix.col(npassive).swap(pulse_matrix.col(alpha_idx));
    }

    // TODO as in cpu NNLS version
    ++iter;
    if (iter % 16 == 0)
        eps_to_use *= 2;
  }
  
  return true;
}

}

}}}
