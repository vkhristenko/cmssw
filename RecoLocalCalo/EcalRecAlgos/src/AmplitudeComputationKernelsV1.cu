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

#include "AmplitudeComputationCommonKernels.h"
#include "AmplitudeComputationKernelsV2.h"

//#define DEBUG

//#define ECAL_RECO_CUDA_DEBUG

namespace ecal { namespace multifit {

namespace debug {

__global__
void kernelUpdateCovarianceMatrix(
        SampleMatrix const* noiseCovariance,
        FullSampleMatrix const* fullPulseCovariance,
        SampleVector const* amplitudes,
        BXVectorType const* bxs,
        char const* acState,
        SampleMatrix *updatedNoiseCovariance,
        int const nchannels) {
    // constants
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = SampleVector::RowsAtCompileTime;

    // indices
    int ch = threadIdx.x + blockIdx.x*blockDim.x;

    if (ch < nchannels) {
        auto const state = static_cast<MinimizationState>(acState[ch]);
        if (state != MinimizationState::NotFinished)
            return;

        // start with the default values
        updatedNoiseCovariance[ch] = noiseCovariance[ch];

        for (unsigned int ipulse=0; ipulse<npulses; ipulse++) {
            if (amplitudes[ch].coeff(ipulse) == 0) 
                continue;

            int const bx = bxs[ch].coeff(ipulse);
            int first_sample_t = std::max(0, bx+3);
            int offset = 7 - 3 - bx;

            auto const value = amplitudes[ch].coeff(ipulse);
            auto const value_sq = value * value;

            unsigned int nsample_pulse = nsamples - first_sample_t;
            updatedNoiseCovariance[ch].block(first_sample_t, first_sample_t,
                                         nsample_pulse, nsample_pulse)
                += value_sq * fullPulseCovariance[ch].block(
                    first_sample_t + offset, first_sample_t + offset,
                    nsample_pulse, nsample_pulse);
        }
    }
}

}

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
    if (state != MinimizationState::NotFinished)
        return;

    // configure shared mem
    __shared__ FullSampleMatrix::Scalar smem[
        FullSampleMatrix::RowsAtCompileTime * FullSampleMatrix::ColsAtCompileTime];
    Eigen::Map<FullSampleMatrix> shrFullPulseCovariance{smem};

    for (unsigned int ity=ty; ity<FullSampleMatrix::RowsAtCompileTime; 
         ity+=nsamples)
        for (unsigned int itx=tx; itx<FullSampleMatrix::ColsAtCompileTime;
            itx+=nsamples)
            shrFullPulseCovariance(ity, itx) = fullPulseCovariance[ch](ity, itx);

    // load the initial matrix value
    auto matrixValue = noiseCovariance[ch](ty, tx);

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
        matrixValue += value_sq * shrFullPulseCovariance(ty + offset, tx+offset);
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
        if (state != MinimizationState::NotFinished)
            return;

        Ls[ch] = covarianceMatrix[ch].llt().matrixL();
    }
}

__global__
void kernel_fast_nnls(
        SampleMatrix const* Ls,
        SampleVector const* samples,
        char* acState,
        PulseMatrixType *pulseMatrix,
        SampleVector *amplitudes,
        int *npassive,
        BXVectorType *activeBXs,
        float *chi2s,
        int nchannels) {
    int const ch = threadIdx.x + blockIdx.x * blockDim.x;
    if (ch < nchannels) {
        auto const state = static_cast<MinimizationState>(acState[ch]);
        if (state != MinimizationState::NotFinished)
            return;

        auto const L = Eigen::internal::LLT_Traits<SampleMatrix, Eigen::Lower>::getL(Ls[ch]); // this gives just a view - should not generate any memcp
//        auto L = covarianceMatrix[ch].llt().matrixL();
        SampleMatrix A = L.solve(pulseMatrix[ch]);
        SampleVector b = L.solve(samples[ch]);
        inplace_fnnls(A, b, amplitudes[ch], npassive[ch],
                      activeBXs[ch],
                      pulseMatrix[ch]);

        // ocmpute chi2
        auto const oldChi2 = chi2s[ch];
        auto const chi2 = L.solve(pulseMatrix[ch]*amplitudes[ch] - samples[ch])
            .squaredNorm();
        auto const delta = chi2 - oldChi2;
        if (ecal::abs(delta) < 1e-3)
            acState[ch] = static_cast<char>(MinimizationState::Finished);
        chi2s[ch] = chi2;
    }
}

template<int BLOCKSIZE>
__global__
void kernel_reduce_state(
        char const* acState,
        char *statePerBlock,
        int nchannels) {
    int const ch = threadIdx.x + blockIdx.x * blockDim.x;

    // configure shared memory and load
    __shared__ char sState[BLOCKSIZE];

    sState[threadIdx.x] = ch<nchannels
        ? acState[ch] != 0 ? 1 : 0
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

void minimization_procedure(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream) {
    int iterations = 0;
    int const maxIterations = 50;
    unsigned int const totalChannels = 
        eventInputCPU.ebDigis->size() + eventInputCPU.eeDigis->size();

    // main loop
    while (true) {
        if (iterations == maxIterations)
            break;
 
        dim3 threadsUpdateCov = {EcalDataFrame::MAXSAMPLES,
            EcalDataFrame::MAXSAMPLES};
        unsigned int blocksUpdateCov = totalChannels;
        kernel_update_covariance_matrix<<<blocksUpdateCov, threadsUpdateCov,
                                          0, conf.cuStream>>>(
            scratch.noisecov, 
            scratch.pulse_covariances, // FIXME we can remove additional matrix
            (SampleVector*)eventOutputGPU.amplitudesAll,
            scratch.activeBXs,
            scratch.acState,
            scratch.updatedNoiseCovariance,
            totalChannels);
        cudaCheck( cudaGetLastError() );

        /*

        // call kernel to compute update covariance matrix
        unsigned int threadsUpdateCov = 32;
        unsigned int blocksUpdateCov = threadsUpdateCov > h_data.digis->size()
            ? 1 
            : (h_data.digis->size() + threadsUpdateCov - 1) / threadsUpdateCov;
        kernelUpdateCovarianceMatrix<<<blocksUpdateCov, threadsUpdateCov>>>(
            d_data.noisecov,
            d_data.pulse_covariances,
            d_data.amplitudes,
            d_data.activeBXs,
            d_data.acState,
            d_data.updatedNoiseCovariance,
            h_data.digis->size());
        ecal::cuda::assert_if_error();

        // call kernel to perform covaraince matrix cholesky decomposition
        unsigned int threadsMatrixLU = 32;
        unsigned int blocksMatrixLU = threadsMatrixLU > totalChannels
            ? 1
            : (totalChannels + threadsMatrixLU - 1) / threadsMatrixLU;
        kernel_matrix_ludecomp<<<blocksMatrixLU, threadsMatrixLU,
                                 0, conf.cuStream>>>(
            d_data.updatedNoiseCovariance,
            d_data.acState,
            d_data.noiseMatrixDecomposition,
            totalChannels);
        AssertIfError

        // call kernel to perform fast nnls
        unsigned int threadsNNLS = 32;
        unsigned int blocksNNLS = threadsNNLS > totalChannels
            ? 1
            : (totalChannels + threadsNNLS - 1) / threadsNNLS;
        kernel_fast_nnls<<<blocksNNLS, threadsNNLS,
                           0, conf.cuStream>>>(
            d_data.noiseMatrixDecomposition,
            d_data.samples,
            d_data.acState,
            d_data.pulse_matrix,
            d_data.amplitudes,
            d_data.npassive,
            d_data.activeBXs,
            d_data.chi2, 
            totalChannels);
        AssertIfError

        // call kernel to reduce in order to generate global state
        constexpr unsigned int threadsState = 256;
        unsigned int blocksForStateReduce = threadsState > totalChannels
            ? 1
            : (totalChannels + threadsState - 1) / threadsState;
        kernel_reduce_state<threadsState><<<blocksForStateReduce, threadsState,
                                            0, conf.cuStream>>>(
            d_data.acState, 
            d_data.minimizationStatePerBlock,
            totalChannels);
        AssertIfError

        // transfer the reductions per block back
        cudaMemcpyAsync(d_data.h_minimizationStatesPerBlock.data(),
                   d_data.minimizationStatePerBlock,
                   blocksForStateReduce * sizeof(MinimizationState),
                   cudaMemcpyDeviceToHost,
                   conf.cuStream);
        AssertIfError

        // sync host with the cuda stream before computing global state
        cudaStreamSynchronize(conf.cuStream);

        // reduce on the host (should be tiny)
        bool acc = true;
        for (unsigned int i=0; i<blocksForStateReduce; i++)
            acc = acc && static_cast<bool>(d_data.h_minimizationStatesPerBlock[i]);
        // global convergence
        if (acc)
            break;
        */

        iterations++;
    }

    //
    // permute computed amplitudes
    // and assign the final uncalibared energy value
    //
    unsigned int threadsPermute = 32 * EcalDataFrame::MAXSAMPLES; // 32 * 10
    unsigned int blocksPermute = threadsPermute > 10 * totalChannels
        ? 1
        : (10 * totalChannels + threadsPermute - 1) / threadsPermute;
    int bytesPermute = threadsPermute * sizeof(SampleVector::Scalar);
    kernel_permute_results<<<blocksPermute, threadsPermute, bytesPermute,
                             conf.cuStream>>>(
        eventOutputGPU.amplitudesAll,
        scratch.activeBXs,
        eventOutputGPU.energies,
        scratch.acState,
        totalChannels);
    cudaCheck( cudaGetLastError() );
}

}}
