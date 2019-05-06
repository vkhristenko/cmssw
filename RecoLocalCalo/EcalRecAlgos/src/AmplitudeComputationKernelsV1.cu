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
#include "AmplitudeComputationKernelsV1.h"

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
    if (state != MinimizationState::NotFinished)
        return;

    // configure shared mem
    __shared__ FullSampleMatrix::Scalar __shrFullPulseCovariance[
        FullSampleMatrix::RowsAtCompileTime * FullSampleMatrix::ColsAtCompileTime];
    __shared__ SampleMatrix::Scalar __shrL[nsamples * nsamples];
    __shared__ SampleMatrix::Scalar __shrLSums[nsamples * nsamples];
    __shared__ SampleMatrix::Scalar shrSamples[nsamples],
                                    shrbPrime[nsamples];
    __shared__ SampleMatrix::Scalar __shrPulseMatrix[nsamples*nsamples],
                                    __shrAPrime[nsmaples * nsamples];
    Eigen::Map<SampleMatrix> shrAPrime{__shrAPrime};
    Eigen::Map<SampleMatrix> shrPulseMatrix{__shrPulseMatirx};
    Eigen::Map<FullSampleMatrix> shrFullPulseCovariance{
        __shrFullPulseCovariance};
    Eigen::Map<SampleMatrix> shrL{__shrL};
    Eigen::Map<SampleMatrix> shrLSums{__shrLSums};

    // copy from global mem 'samples' to shared mem
    if (ty==0)
        shrSamples[tx] = samples[ch](tx);

    // shared mem store from global
//    shrPulseMatrix(ty, tx) = pulseMatrix[ch](ty, tx);
    auto pm_i_j = pulseMatrix[ch](ty, tx);

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

    //
    // cholesky decomposition of 10 by 10 matrix 
    // fused with forward substitution in order to use solvers
    //

    // init L sums
    shrLSums(ty, tx) = 0;

    // store only diagonal elements 
    __shared__ SampleMatrix::Scalar shrDiagElements[nsamples];
    if (ty == tx)
        shrDiag[tx] = matrixValue;
    __syncthreads();

    // column 0
    // compute L(ty, 0) for ty >=1
    auto const m_j_j = shrDiag[tx];
    if (tx == 0 && ty >= 1)
        shrL(ty, 0) = matrixValue / std::sqrt(m_j_j);
    else if (tx == 0 && ty == 0)
        shrL(0, 0) = std::sqrt(m_j_j);
    __syncthreads();
    
    // compute the first row for forward substitution
    // compute subs for all elements after the first row
    if (ty == 0) {
        shrAPrime(0, tx) = shrPulseMatrix(0, tx) / shrL(0, 0);
    } else { // ty >= 1
        shrPulseMatrix(ty, tx) -= shrL(ty, 0) * shrPulseMatrix(0, tx) / shrL(0, 0);
    }

    // use the opposite column (threads for that col)
    // to compute values for forward substitution
    if (tx==9 && ty == 0) {
        shrbPrime[0] = shrSamples[0] / shrL(0, 0);
    } else if (tx==9 && ty >= 1) {
        shrSamples[ty] -= shrL(ty, 0) * shrSamples[0] / shrL(0, 0);
    }

    // TODO: verify that the loop is unrolled
    #pragma unroll
    for (int column=1; column<nsamples-1; column++) {
        if (tx==column && ty>=column) {
            // compute L(j, j) = sqrt(M(j, j) - Sum[k](L(j, k) * L(j, k)))
            auto const sumsq = shrLSums(column, column) + 
                shrL(column, column-1) * shrL(column, column-1);
            auto const l_j_j = std::sqrt(m_j_j - sumsq);
            if (ty == column)
                shrL(column, column) = l_j_j;
            else {
                // compute L(i, column) = 
                //   (M(i, column) - Sum[k](L(i, k) * L(column, k))) / L(col, col)
                auto const sumsq_i_j = shrLSums(ty, column) + 
                    shrL(ty, column-1) * shrL(column, column-1);
                auto const tmp = m_i_j - sumsq_i_j;
                auto const l_i_column = tmp / l_j_j;
                shrL(ty, column) = l_i_column;
            }
        }
        if (tx>=column && ty>=column)
            shrLSums(ty, tx) += shrL(ty, column-1) * shrL(tx, column-1);
        __syncthreads();

        // forward substitution computation
        if (ty >= column) {
            if (ty == column) {
                shrAPrime(column, tx) = shrPulseMatrix(column, tx) 
                    / shrL(column, column);
            } else {
                shrPulseMatrix(ty, tx) -= shrL(ty, column) * 
                    shrPulseMatrix(column, tx) / shrL(column, column);
            }
        }

        // use a different group of threads
        // to compute values for forward substitution
        if (tx==column-1 && ty>=column) {
            if (ty == column)
                shrbPrime[column] = shrSamples[column] / shrL(column, column);
            else 
                shrSamples[ty] -= shrL(ty, column) * shrSamples[column] 
                    / shrL(column, column);
        }
    }

    // last iteration for column = 9
    constexpr auto column = nsamples - 1;
    if (tx==column && ty>=column) {
        // compute L(j, j) = sqrt(M(j, j) - Sum[k](L(j, k) * L(j, k)))
        auto const sumsq = shrLSums(column, column) + 
            shrL(column, column-1) * shrL(column, column-1);
        auto const l_j_j = std::sqrt(m_j_j - sumsq);
        if (ty == column)
            shrL(column, column) = l_j_j;
        else {
            // compute L(i, column) = 
            //   (M(i, column) - Sum[k](L(i, k) * L(column, k))) / L(col, col)
            auto const sumsq_i_j = shrLSums(ty, column) + 
                shrL(ty, column-1) * shrL(column, column-1);
            auto const tmp = m_i_j - sumsq_i_j;
            auto const l_i_column = tmp / l_j_j;
            shrL(ty, column) = l_i_column;
        }
    }
    __syncthreads();
    if (ty == 0 )
        shrAPrime(column, tx) = shrPulseMatrix(column, tx) / shrL(column, column);
    // note, we are using different threads
    if (tx == 0 && ty==1)
        shrbPrime[column] = shrSamples[column] / shrL(column, column);

    // store to global memory
    g_L[ch](ty, tx) = shrL(ty, tx);
    __syncthreads();

    // 
    // compute backward substitution
    // reuse memory:
    //   shrPulseMatrix -> A
    //   shrSamples -> b
    //
    auto& shrA = shrPulseMatrix;
    auto& shrb = shrSamples;
    
    // compute row 9
    if (ty==9) {
        shrA(9, tx) = shrAPrime(9, tx) / shrL(9, 9);
        if (tx == 9)
            shrb[9] = shrbPrime[9] / shrL(9, 9);
        else // tx <= 8
            shrbPrime[tx] -= shrL(9, tx) * shrbPrime[9] / shrL(9, 9);
    } else {
        // ty<=8
        shrAPrime(ty, tx) -= shrL(9, ty) * shrAPrime(9, tx) / shrL(9, 9);
    }

    // compute rows 8 - 1
    #pragma unroll
    for (int row=nsamples-2; row>=1; row--) {
        // we need only threads above or equal `row`
        if (ty == row) {
            shrA(row, tx) = shrAPrime(row, tx) / shrL(row, row);
        } else if (ty < row) {
            shrAPrime(ty, tx) -= shrL(row, tx)
                * shrAPrime(row, tx) / shrL(row, row);
        }

        // use different 10 threads than the ones above
        if (ty == nsamples - 1 && tx==row) {
            shrb[row] = shrbPrime[row] / shrL(row, row);
        } else if (ty==nsamples-1 && tx<row) {
            shrbPrime[tx] -= shrL(row, tx) * shrbPrime[row] / shrL(row, row);
        }
        __syncthreads();
    }

    // last iteration
    if (ty == 0)
        shrA(0, tx) = shrAPrime(0, tx) / shr(0, 0);
    if (ty == nsamples-1 && tx==0)
        shrb[0] = shrbPrime[0] / shrL(0, 0);
    __syncthreads();

    //
    // compute AtA and Atb
    //
    auto ata_i_j = 0f;
    for (unsigned int k=0; k<nsamples; k++)
        ata_i_j += shrA(k, ty) * shrA(k, tx);
    g_AtA[ch](ty, tx) = ata_i_j;
    if (ty==0) {
        auto atb_i = 0f;
        for (unsigned int k=0; k<nsamples; k++)
            atb_i += shrA(k, tx) * shrb[k];
        g_Atb[ch](tx) = atb_i;
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
    using DataType = SampleVector::Scalar;

    // indices / constants
    auto const ch = threadIdx.x + blockIdx.x * blockDim.x;
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    
    // remove unneeded threads
    if (ch >= nchannels)
        return;
    
    auto const state = static_cast<MinimizationState>(acState[ch]);
    if (state != MinimizationState::NotFinished)
        return;

    // TODO: add noisecov check
    
    // shared mem conf
    extern __shared__ char* smem;
    DataType* __shrAtA = reinterpret_cast<DataType*>(smem);
    Eigen::Map<SampleMatrix> shrAtA{shrAtARaw};
    DataType* __shrL = reinterpret_cast<DataType*>(
        __shrAtA + nsamples*nsamples*blockDim.x);
    Eigen::Map<SampleMatrix> shrL{__shrL + threadIdx.x*nsamples*nsamples};
    DataType* shrAtb = reinterpret_cast<DataType*>(
        __shrL + nsamples*nsamples*blockDim.x);
    DataType* shrx = reinterpret_cast<DataType*>(
        shrAtb + nsamples*blockDim.x);
    DataType* shrs = reinterpret_cast<DataType*>(
        shrx + nsamples*blockDim.x);
    DataType* shrw = reinterpret_cast<DataType*>(
        shrs + nsamples*blockDim.x);
    char* permutation = reinterpret_cast<char*>(shrw) + nsamples*blockDim.x;

    // load/store to shared/regs
    shrAtA = g_AtA[ch];
    shrAtb = g_Atb[ch];
    auto npassive = g_npassive[ch];

    // 
    Eigen::Index w_max_idx_prev = 0;
    DataType w_max_prev = 0;
    double eps_to_use = eps;
    for (int iteration=0; iteration<maxIterations; iteration++) {
        if (iteration>0 || npassive==0) {
            // we are not infinite
            if (iteration >= 500)
                break;

            auto const nactive = nsamples - npassive;
            if (!nactive)
                break;

            shrw.tail(nactive) = shrAtb.tail(nactive) - (shrAtA*shrx).tail(nactive);

            // get the index of w that gives the max gain
            Eigen::Index w_max_idx;
            auto const w_max = w.tail(nactive).maxCoeff(&w_max_idx);
            if (w_max < eps_to_use || 
                (w_max_idx==w_max_idx_prev && w_max==w_max_prev))
                break;

            w_max_prev = w_max;
            w_max_idx_prev = w_max_idx;
            w_max_idx += npassive;
            Eigen::numext::swap(permutation[npassive], permutation[w_max_idx]);
            ++npassive;
        }
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
        eventInputCPU.ebDigis.size() + eventInputCPU.eeDigis.size();

    // main loop
    while (true) {
        if (iterations == maxIterations)
            break;
 
        dim3 threadsUpdateCov = {EcalDataFrame::MAXSAMPLES,
            EcalDataFrame::MAXSAMPLES};
        unsigned int blocksUpdateCov = totalChannels;
        kernel_update_covariance_matrix<<<blocksUpdateCov, threadsUpdateCov,
                                          0, cudaStream.id()>>>(
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
                             cudaStream.id()>>>(
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.activeBXs,
        eventOutputGPU.amplitude,
        scratch.acState,
        totalChannels);
    cudaCheck( cudaGetLastError() );
}

}}
