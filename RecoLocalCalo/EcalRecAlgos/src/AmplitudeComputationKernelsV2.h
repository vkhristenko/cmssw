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

#define RUN_UPDATE_COVARIANCE
#ifdef RUN_UPDATE_COVARIANCE
__global__
void kernel_update_covariance_matrix(
        SampleMatrix const* noiseCovariance,
        FullSampleMatrix const* fullPulseCovariance,
        SampleVector const* amplitudes,
        BXVectorType const* bxs,
        char const* acState,
        SampleMatrix *updatedNoiseCovariance,
        int nchannels);
#endif

#define RUN_COVARIANCE_DECOMPOSITION
#ifdef RUN_COVARIANCE_DECOMPOSITION
__global__
void kernel_matrix_ludecomp(
        SampleMatrix const* covarianceMatrix,
        char const* acState,
        SampleMatrix *Ls,
        int nchannels);
#endif

#define RUN_FAST_NNLS
#ifdef RUN_FAST_NNLS
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
#endif

#define RUN_COMPUTE_CHI2
#ifdef RUN_COMPUTE_CHI2
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
#endif

#define RUN_REDUCE_STATE
#ifdef RUN_REDUCE_STATE
__global__
void kernel_reduce_state(
        char const* state,
        char *statePerBlock,
        int nchannels);
#endif

__global__
void kernelInitializeBeforeMinimizationProcedure(
        PermutationMatrix *permutation,
        int *npassive,
        char *minimizationStatePerBlock,
        int nchannels, 
        int blocksForStateInitialization);

#define RUN_MINIMIZATION_PROCEDURE
#ifdef RUN_MINIMIZATION_PROCEDURE
void minimization_procedure(
        device_data& d_data, 
        host_data& h_data);
#endif

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2
