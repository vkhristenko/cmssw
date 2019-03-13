#ifndef kernel_minimize_cpu_test_h
#define kernel_minimize_cpu_test_h

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
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"

//#define DEBUG

namespace ecal { namespace multifit { namespace v1 {

namespace cpu {

bool update_covariance(SampleMatrix const& noisecov,
                       FullSampleMatrix const& full_pulse_cov,
                       SampleMatrix& inverse_cov,
                       BXVectorType const& bxs,
                       SampleDecompLLT& covariance_decomposition,
                       SampleVector const& amplitudes);

float compute_chi2(SampleDecompLLT& covariance_decomposition,
                   PulseMatrixType const& pulse_matrix,
                   SampleVector const& amplitudes,
                   SampleVector const& samples);

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
                     bool gainSwitchUseMaxSample);

using matrix_t = SampleMatrix;
using vector_t = SampleVector;
using permutation_t = Eigen::PermutationMatrix<SampleMatrix::RowsAtCompileTime>;

bool
inplace_fnnls(matrix_t const& A,
    vector_t const& b,
    vector_t& x,
    int& npassive,
    BXVectorType& activeBXs,
    permutation_t& permutation,
    PulseMatrixType& pulse_matrix,
    const double eps = 1e-11,
    const unsigned int max_iterations = 500);
}

}}}

#endif
