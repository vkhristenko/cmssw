#ifndef RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2
#define RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2

#include <iostream>
#include <limits>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/DeclsForKernels.h"

#include "cuda.h"

//#define DEBUG

//#define ECAL_RECO_CUDA_DEBUG

namespace ecal { namespace multifit {

__global__
void kernel_update_covariance_matrix(
        SampleMatrix const* noiseCovariance,
        FullSampleMatrix const* fullPulseCovariance,
        SampleVector const* amplitudes,
        BXVectorType const* bxs,
        char const* acState,
        SampleMatrix *updatedNoiseCovariance,
        int nchannels);

__global__
void kernel_matrix_ludecomp(
        SampleMatrix const* covarianceMatrix,
        char const* acState,
        SampleMatrix *Ls,
        int nchannels);

__global__
void kernel_fast_nnls(
        SampleMatrix const* Ls,
        SampleVector const* samples,
        char const* acState,
        PulseMatrixType *pulseMatrix,
        SampleVector *amplitudes,
        int *npassive,
        BXVectorType *activeBXs,
        PermutationMatrix *permutation,
        int nchannels);

__global__
void kernel_compute_chi2_and_propogate_quantities(
        SampleMatrix const* Ls,
        SampleVector const* samples,
        PulseMatrixType const* pulseMatrix,
        PermutationMatrix const* permutation,
        SampleVector* amplitudes,
        float *chi2s,
        char *acState,
        float *energies,
        int nchannels);

__global__
void kernel_reduce_state(
        char const* state,
        char *statePerBlock,
        int nchannels);

__global__
void kernelInitializeBeforeMinimizationProcedure(
        PermutationMatrix *permutation,
        int *npassive,
        char *minimizationStatePerBlock,
        int nchannels, 
        int blocksForStateInitialization);

void minimization_procedure(
        device_data& d_data, 
        host_data& h_data);

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2
