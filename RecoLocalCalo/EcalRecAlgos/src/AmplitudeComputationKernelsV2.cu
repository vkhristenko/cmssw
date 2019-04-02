#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"


#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "inplace_fnnls.h"

#include "AmplitudeComputationKernelsV2.h"

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
        int nchannels) {
    // constants
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = SampleVector::RowsAtCompileTime;

    // indices
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int ch = blockIdx.x;

    auto const state = static_cast<MinimizationState>(acState[ch]);
    if (state == MinimizationState::Finished)
        return;

    // configure shared mem
    //__shared__ SampleVector::Scalar
    
    // load the initial matrix value
    auto matrixValue = noiseCovariance[ch](ty, tx);
    auto const fullMatrixValue1 = fullPulseCovariance[ch](ty, tx);
    auto const fullMatrixValue2 = fullPulseCovariance[ch](ty+9, tx+9);

#pragma unroll
    for (int ipulse=0; ipulse<npulses; ipulse++) {
        auto const amplitude = amplitudes[ch].coeff(ipulse);
        if (amplitude == 0)
            continue;

        int const bx = bxs[ch].coeff(ipulse);
        int const first_sample_t = std::max(0, bx+3);
        int const offset = 7 - 3 - bx;
        int const first_full_sample_t = first_sample_t + offset;
        auto const nsample_pulse = nsamples - first_sample_t;

        // note, by definition tx and ty are both less than the last sample 
        // (nsamples)
        if (!(tx>=first_sample_t && ty>=first_sample_t))
            continue;

        // update the matrix element
        auto const value = amplitude;
        auto const value_sq = value * value;
        matrixValue += value_sq * 
            (first_full_sample_t <= 9 
                ? fullMatrixValue1
                : fullMatrixValue2);
    }

    // store into global
    updatedNoiseCovariance[ch](ty, tx) = matrixValue;

}

__global__
void kernel_matrix_ludecomp(
        SampleMatrix const* covarianceMatrix,
        char const* acState,
        SampleMatrix *Ls,
        int nchannels) {
    int const ch = threadIdx.x + blockIdx.x * blockDim.x;
    if (ch < nchannels) {
        auto const state = static_cast<MinimizationState>(acState[ch]);
        if (state == MinimizationState::Finished)
            return;

        Ls[ch] = covarianceMatrix[ch].llt().matrixL();
    }
}

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
        int nchannels) {
    int const ch = threadIdx.x + blockIdx.x * blockDim.x;
    if (ch < nchannels) {
        auto const state = static_cast<MinimizationState>(acState[ch]);
        if (state == MinimizationState::Finished)
            return;

        auto const L = Eigen::internal::LLT_Traits<SampleMatrix, Eigen::Lower>::getL(Ls[ch]); // this gives just a view - should not generate any memcp
        auto A = L.solve(pulseMatrix[ch]);
        auto b = L.solve(samples[ch]);
        inplace_fnnls(A, b, amplitudes[ch], npassive[ch],
                      activeBXs[ch], permutation[ch], 
                      pulseMatrix[ch]);
    }
}

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
        int nchannels) {
    int const ch = threadIdx.x + blockIdx.x * blockDim.x;

    if (ch < nchannels) {
        auto const state = static_cast<MinimizationState>(acState[ch]);
        if (state == MinimizationState::Finished)
            return;

        // TODO
        // use LLT traits to get triangular view
        auto oldChi2 = chi2s[ch];
        // just a view
        auto const L = Eigen::internal::LLT_Traits<SampleMatrix, Eigen::Lower>::getL(Ls[ch]);
        auto chi2 = L
            .solve(pulseMatrix[ch]*amplitudes[ch] - samples[ch])
            .squaredNorm();
        auto delta = chi2 - oldChi2;
        chi2s[ch] = chi2;

        // if converged, permute amplitudes and set energy value
        // TODO: may be variables propogation could be done more systematically
        // or more efficiently
        if (ecal::abs(delta) < 1e-3) {
            acState[ch] = static_cast<char>(MinimizationState::Finished);
            amplitudes[ch] = amplitudes[ch].transpose() * 
                permutation[ch].transpose();
            energies[ch] = amplitudes[ch](5);
//            statuses[ch] = 
        }
    }
}

__global__
void kernel_reduce_state(
        char const* state,
        char *statePerBlock,
        int nchannels) {
    int const ch = threadIdx.x + blockIdx.x * blockDim.x;

    // configure shared memory and load
    extern __shared__ char smem[];
    char* sState = smem;
    sState[threadIdx.x] = ch<nchannels
        ? state[ch]
        : 1;
    __syncthreads();

    // reduce per block
    for (int s=blockDim.x/2; s>0; s>>=1) {
        if (threadIdx.x < s)
            sState[threadIdx.x] &= sState[threadIdx.x + s];
        __syncthreads();
    }

    // store to global
    if (threadIdx.x == 0)
        statePerBlock[blockIdx.x] = sState[0];
}

__global__
void kernelInitializeBeforeMinimizationProcedure(
        PermutationMatrix *permutation,
        int *npassive,
        char *minimizationStatePerBlock,
        int nchannels, 
        int blocksForStateInitialization) {
    int ch = threadIdx.x + blockIdx.x * blockDim.x;

    if (ch < nchannels) {
        permutation[ch].setIdentity();
        npassive[ch] = 0;

        if (ch < blocksForStateInitialization)
            minimizationStatePerBlock[ch] = 0;
    }
}

void minimization_procedure(
        device_data& d_data, 
        host_data& h_data) {
    int iterations = 0;
    int const maxIterations = 50;

    // TODO: specify all the constants in 1 location
    // initialize all the varibles before starting minimization
    constexpr int blocksForStateInitialization = 1000;
    constexpr int maxChannels = 20000;
    unsigned int threadsInit = 32;
    unsigned int blocksInit = h_data.digis->size()
        ? 1
        : (h_data.digis->size() + threadsInit - 1) / threadsInit;
    kernelInitializeBeforeMinimizationProcedure<<<blocksInit, threadsInit>>>(
        d_data.permutation,
        d_data.npassive,
        d_data.minimizationStatePerBlock,
        h_data.digis->size(),
        blocksForStateInitialization);

    // main loop
    while (true) {
        if (iterations == maxIterations)
            break;
        
        // call kernel to compute update covariance matrix
        dim3 threads_update_cov_matrix{EcalDataFrame::MAXSAMPLES, 
            EcalDataFrame::MAXSAMPLES};
        int blocks_update_cov_matrix = h_data.digis->size();
        kernel_update_covariance_matrix<<<blocks_update_cov_matrix,
                                          threads_update_cov_matrix>>>(
            d_data.noisecov,
            d_data.pulse_covariances,
            d_data.amplitudes,
            d_data.bxs,
            d_data.acState,
            d_data.updatedNoiseCovariance,
            h_data.digis->size());
        cudaDeviceSynchronize();
        ecal::cuda::assert_if_error();

        // call kernel to perform covaraince matrix cholesky decomposition
        unsigned int threadsMatrixLU = 32;
        unsigned int blocksMatrixLU = threadsMatrixLU > h_data.digis->size()
            ? 1
            : (h_data.digis->size() + threadsMatrixLU - 1) / threadsMatrixLU;
        kernel_matrix_ludecomp<<<blocksMatrixLU, threadsMatrixLU>>>(
            d_data.updatedNoiseCovariance,
            d_data.acState,
            d_data.noiseMatrixDecomposition,
            h_data.digis->size());
        cudaDeviceSynchronize();
        ecal::cuda::assert_if_error();

        // call kernel to perform fast nnls
        unsigned int threadsNNLS = 32;
        unsigned int blocksNNLS = threadsNNLS > h_data.digis->size()
            ? 1
            : (h_data.digis->size() + threadsNNLS - 1) / threadsNNLS;
        kernel_fast_nnls<<<blocksNNLS, threadsNNLS>>>(
            d_data.noiseMatrixDecomposition,
            d_data.samples,
            d_data.acState,
            d_data.pulse_matrix,
            d_data.amplitudes,
            d_data.npassive,
            d_data.bxs,
            d_data.permutation,
            h_data.digis->size());
        cudaDeviceSynchronize();
        ecal::cuda::assert_if_error();

#ifdef XXXX
        // call kernel to compute chi2
        unsigned int threadsChi2 = 32;
        unsigned int blocksChi2 = threadsChi2 > h_data.digis->size()
            ? 1
            : (h_data.digis->size() + threadsChi2 - 1) / threadsChi2;
        kernel_compute_chi2_and_propogate_quantities<<<blocksChi2, threadsChi2>>>(
            d_data.noiseMatrixDecomposition,
            d_data.samples,
            d_data.pulse_matrix,
            d_data.permutation,
            d_data.amplitudes,
            d_data.chi2,
            d_data.acState,
            d_data.energies,
            h_data.digis->size());
        cudaDeviceSynchronize();
        ecal::cuda::assert_if_error();

        // call kernel to reduce to generate global state
        unsigned int threadsState = 256;
        unsigned int blocksForStateReduce = threadsState > h_data.digis->size()
            ? 1
            : (h_data.digis->size() + threadsState - 1) / threadsState;
        kernel_reduce_state<<<blocksForStateReduce, threadsState>>>(
            d_data.acState, 
            d_data.minimizationStatePerBlock,
            h_data.digis->size());

        // transfer the reductions per block back
//#define RUN_FINALIZE_MINIMIZATION_PROCEDURE
        cudaMemcpy(d_data.h_minimizationStatesPerBlock.data(),
                   d_data.minimizationStatePerBlock,
                   blocksForStateReduce * sizeof(MinimizationState),
                   cudaMemcpyDeviceToHost);
        // reduce on the host (should be tiny)
        char acc = 1;
        for (unsigned int i=0; i<blocksForStateReduce; i++)
            acc &= d_data.h_minimizationStatesPerBlock[i];
        // global convergence
        if (static_cast<MinimizationState>(acc) == MinimizationState::Finished)
            break;
#endif

        iterations++;
    }
}

}}
