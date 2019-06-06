#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "inplace_fnnls.h"
#include "AmplitudeComputationKernelsV2.h"
#include "AmplitudeComputationCommonKernels.h"
#include "KernelHelpers.h"

namespace ecal { namespace multifit { namespace v2 {

#define PRINT_MATRIX_10x10(M)\
            printf("%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n", \
                M(0, 0), M(1, 0), M(2, 0), M(3, 0), M(4, 0), \
                M(5, 0), M(6, 0), M(7, 0), M(8, 0), M(9, 0), \
                M(0, 1), M(1, 1), M(2, 1), M(3, 1), M(4, 1), \
                M(5, 1), M(6, 1), M(7, 1), M(8, 1), M(9, 1), \
                M(0, 2), M(1, 2), M(2, 2), M(3, 2), M(4, 2), \
                M(5, 2), M(6, 2), M(7, 2), M(8, 2), M(9, 2), \
                M(0, 3), M(1, 3), M(2, 3), M(3, 3), M(4, 3), \
                M(5, 3), M(6, 3), M(7, 3), M(8, 3), M(9, 3), \
                M(0, 4), M(1, 4), M(2, 4), M(3, 4), M(4, 4), \
                M(5, 4), M(6, 4), M(7, 4), M(8, 4), M(9, 4), \
                M(0, 5), M(1, 5), M(2, 5), M(3, 5), M(4, 5), \
                M(5, 5), M(6, 5), M(7, 5), M(8, 5), M(9, 5), \
                M(0, 6), M(1, 6), M(2, 6), M(3, 6), M(4, 6), \
                M(5, 6), M(6, 6), M(7, 6), M(8, 6), M(9, 6), \
                M(0, 7), M(1, 7), M(2, 7), M(3, 7), M(4, 7), \
                M(5, 7), M(6, 7), M(7, 7), M(8, 7), M(9, 7), \
                M(0, 8), M(1, 8), M(2, 8), M(3, 8), M(4, 8), \
                M(5, 8), M(6, 8), M(7, 8), M(8, 8), M(9, 8), \
                M(0, 9), M(1, 9), M(2, 9), M(3, 9), M(4, 9), \
                M(5, 9), M(6, 9), M(7, 9), M(8, 9), M(9, 9) \
            )

__device__
__forceinline__
bool update_covariance(SampleMatrix const& noisecov,
                       EcalPulseCovariance const& pulse_cov,
                       SampleMatrix& inverse_cov,
                       BXVectorType const& bxs,
                       SampleDecompLLT& covariance_decomposition,
                       SampleVector const& amplitudes) {
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = BXVectorType::RowsAtCompileTime;

    inverse_cov = noisecov;

    for (unsigned int ipulse=0; ipulse<npulses; ipulse++) {
        if (amplitudes.coeff(ipulse) == 0) 
            continue;

        int bx = bxs.coeff(ipulse);
        int first_sample_t = std::max(0, bx+3);
        int offset = 3 - bx;

        auto const value = amplitudes.coeff(ipulse);
        auto const value_sq = value*value;

        unsigned int nsample_pulse = nsamples - first_sample_t;
        for (int i=first_sample_t; i<first_sample_t+nsample_pulse; ++i)
            for (int j=first_sample_t; j<first_sample_t+nsample_pulse; ++j)
                inverse_cov(i, j) += value_sq * pulse_cov.covval[i+offset][j+offset];

        /*
        inverse_cov.block(first_sample_t, first_sample_t, 
                          nsample_pulse, nsample_pulse)
            += value_sq * full_pulse_cov.block(first_sample_t + offset,
                                               first_sample_t + offset,
                                               nsample_pulse,
                                               nsample_pulse);
                                               */
    }

    return true;
}

__global__
void kernel_update_covariance_compute_cholesky() {
    // constants
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto npulses = nsamples;
    using DataType = SampleVector::Scalar;

    // indices
    auto const tx = threadIdx.x;
    auto const ty = threadIdx.y;
    auto const ch_per_block = threadIdx.z;
    auto const gch = ch_per_block + blockIdx.x*blockDim.z;

    // shared mem configuration
    extern __shared__ char smem[];
    DataType* __shrPulseCovariance = reinterpret_cast<DataType*>(smem);

    // preload some values
    auto noiseValue = noiseCovariance(ty, tx);

    // compute the updated total covariance matrix
    #pragma unroll
    for (int ipulse=0; ipulse<npulses; ipulse++) {
        auto const amplitude = amplitudes(ipulse);
        if (amplitude == 0) 
            continue;

        // note, we do not permute results anymore!
        // therefore straightforward mapping
        auto const bx = ipulse - 5;
        auto const first_sample_t 
    }
}

__global__
void kernel_solve() {

}

__global__
void kernel_mm_mv_mults() {

}

__global__
void kernel_fnnls() {
}

__global__
void kernel_compute_chi2() {

}

__global__
void kernel_minimization_launcher() {
    // we have a single guy in here
    if (threadIdx.x >= 1) return;

    int iter = 0;
    auto channelsLeft = totalChannels;
    while (true) {
        if (iter >= max_iterations)
            break;

        // 
        kernel_update_covariance_compute_cholesky();

        // forward / backward substitution
        kernel_solve();

        // matrix matrix and matrix vector mults
        kernel_mm_mv_mults();

        // fnnls 
        kernel_fnnls();

        // compute chi2
        kernel_compute_chi2();
        cudaDeviceSynchronize();

        // check exit
        channelsLeft = *channelsCounter;
        if (channelsLeft == 0)
            return;
    }
}

__device__
__forceinline__
SampleVector::Scalar compute_chi2(SampleDecompLLT& covariance_decomposition,
                   PulseMatrixType const& pulse_matrix,
                   SampleVector const& amplitudes,
                   SampleVector const& samples) {
    return covariance_decomposition.matrixL()
        .solve(pulse_matrix * amplitudes - samples)
        .squaredNorm();
}

///
/// launch ctx parameters are (nchannels / block, blocks)
/// TODO: trivial impl for now, there must be a way to improve
///
/// Conventions:
///   - amplitudes -> solution vector, what we are fitting for
///   - samples -> raw detector responses
///   - passive constraint - satisfied constraint
///   - active constraint - unsatisfied (yet) constraint
///
__global__
void kernel_minimize(SampleMatrix const* noisecov,
                     EcalPulseCovariance const* pulse_cov,
                     BXVectorType *bxs,
                     SampleVector const* samples,
                     SampleVector* amplitudes,
                     PulseMatrixType* pulse_matrix, 
                     ::ecal::reco::StorageScalarType* chi2s,
                     uint32_t const* dids,
                     char *acState,
                     int nchannels,
                     int max_iterations,
                     unsigned int offsetForHashes) {
    int idx = threadIdx.x + blockDim.x*blockIdx.x;
    if (idx < nchannels) {
        if (static_cast<MinimizationState>(acState[idx]) == 
            MinimizationState::Precomputed)
            return;

        auto const did = DetId{dids[idx]};
        auto const isBarrel = did.subdetId() == EcalBarrel;
        auto const hashedId = isBarrel
            ? hashedIndexEB(did.rawId())
            : offsetForHashes + hashedIndexEE(did.rawId());


        // inits
        bool status = false;
        int iter = 0;
        int npassive = 0;
        
        // FIXME: 
        // need to identify this earrlier and set the state properly
        if (noisecov[idx].isZero(0))
            return;

        // inits
        SampleDecompLLT covariance_decomposition;
        SampleMatrix inverse_cov;
        SampleVector::Scalar chi2 = 0, chi2_now = 0;

#ifdef ECAL_MULTIFIT_KERNEL_MINIMIZE_V1
//    PRINT_MATRIX_10x10(noisecov[idx]);
#endif

        // loop until ocnverge
        while (true) {
            if (iter >= max_iterations)
                break;

            status = update_covariance(
                noisecov[idx], 
                pulse_cov[hashedId],
                inverse_cov,
                bxs[idx],
                covariance_decomposition,
                amplitudes[idx]);

            // compute actual covariance decomposition
            covariance_decomposition.compute(inverse_cov);

            // prepare input matrices for fnnls
            SampleMatrix A = covariance_decomposition.matrixL()
                .solve(pulse_matrix[idx]);
            SampleVector b = covariance_decomposition.matrixL()
                .solve(samples[idx]);
            
            status = inplace_fnnls(
                A, b, amplitudes[idx],
                npassive, bxs[idx], pulse_matrix[idx]);
                
            chi2_now = compute_chi2(
                covariance_decomposition,
                pulse_matrix[idx],
                amplitudes[idx],
                samples[idx]);
            auto deltachi2 = chi2_now - chi2;


#ifdef ECAL_MULTIFIT_KERNEL_MINIMIZE_V1
            if (iter > 10) {
                printf("idx = %d iter = %d chi2 = %f chi2old = %f\n",
                    idx, iter, chi2_now, chi2);
                
                printf("noisecov(0, i): %f %f %f %f %f %f %f %f %f %f\n",
                    noisecov[idx](0, 0),
                    noisecov[idx](0, 1),
                    noisecov[idx](0, 2),
                    noisecov[idx](0, 3),
                    noisecov[idx](0, 4),
                    noisecov[idx](0, 5),
                    noisecov[idx](0, 6),
                    noisecov[idx](0, 7),
                    noisecov[idx](0, 8),
                    noisecov[idx](0, 9));

                printf("ampls: %f %f %f %f %f %f %f %f %f %f\n",
                    amplitudes[idx](0),
                    amplitudes[idx](1),
                    amplitudes[idx](2),
                    amplitudes[idx](3),
                    amplitudes[idx](4),
                    amplitudes[idx](5),
                    amplitudes[idx](6),
                    amplitudes[idx](7),
                    amplitudes[idx](8),
                    amplitudes[idx](9));
            }
#endif

            chi2 = chi2_now;

            if (ecal::abs(deltachi2) < 1e-3)
                break;

            //---- AM: TEST
            //---- it was 3 lines above, now here as in the CPU version
            ++iter;
            
        }

        // the rest will be set later
        chi2s[idx] = chi2;
    }
}

void minimization_procedure(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream,
        unsigned int offsetForHashes) {
    unsigned int totalChannels = eventInputCPU.ebDigis.size() 
        + eventInputCPU.eeDigis.size();
//    unsigned int threads_min = conf.threads.x;
    // TODO: configure from python
    unsigned int threads_min = configParameters.kernelMinimizeThreads[0];
    unsigned int blocks_min = threads_min > totalChannels
        ? 1
        : (totalChannels + threads_min - 1) / threads_min;
    kernel_minimize<<<blocks_min, threads_min, 0, cudaStream.id()>>>(
        scratch.noisecov,
        conditions.pulseCovariances.values,
        scratch.activeBXs,
        scratch.samples,
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.pulse_matrix,
        eventOutputGPU.chi2,
        eventInputGPU.ids,
        scratch.acState,
        totalChannels,
        50,
        offsetForHashes);
    cudaCheck(cudaGetLastError());

    //
    // permute computed amplitudes
    // and assign the final uncalibared energy value
    //
    unsigned int threadsPermute = 32 * EcalDataFrame::MAXSAMPLES; // 32 * 10
    unsigned int blocksPermute = threadsPermute > 10 * totalChannels
        ? 1
        : (10 * totalChannels + threadsPermute - 1) / threadsPermute;
    int bytesPermute = threadsPermute * sizeof(SampleVector::Scalar);
    kernel_permute_results<<<blocksPermute, threadsPermute, 
                             bytesPermute, cudaStream.id()>>>(
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.activeBXs,
        eventOutputGPU.amplitude,
        scratch.acState,
        totalChannels);
    cudaCheck(cudaGetLastError());
}

}

}}
