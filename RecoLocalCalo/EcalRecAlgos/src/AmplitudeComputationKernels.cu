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
#include "KernelHelpers.h"
#include "AmplitudeComputationKernels.h"
#include "AmplitudeComputationCommonKernels.h"

namespace ecal { namespace multifit {

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
bool update_covariance(
        EcalPulseCovariance const& pulse_covariance,
        SampleMatrix& inverse_cov,
        BXVectorType const& bxs,
        SampleVector const& amplitudes) {
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = BXVectorType::RowsAtCompileTime;

    for (unsigned int ipulse=0; ipulse<npulses; ipulse++) {
        auto const amplitude = amplitudes.coeff(ipulse);
        if (amplitude == 0) 
            continue;

        int bx = bxs.coeff(ipulse);
        int first_sample_t = std::max(0, bx+3);
        int offset = -3 - bx;

        auto const value_sq = amplitude * amplitude;

        unsigned int nsample_pulse = nsamples - first_sample_t;

        for (int row=first_sample_t; row<nsamples; row++) {
            for (int col=first_sample_t; col<nsamples; col++) {
                inverse_cov.coeffRef(row, col) += value_sq * 
                    __ldg(&pulse_covariance.covval[row + offset][col + offset]);
            }
        }
    }

    return true;
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
void kernel_minimize(
        uint32_t const* dids_eb,
        uint32_t const* dids_ee,
        SampleMatrix const* noisecov,
        EcalPulseCovariance const* pulse_covariance,
        BXVectorType *bxs,
        SampleVector const* samples,
        SampleVector* amplitudes,
        PulseMatrixType* pulse_matrix, 
        ::ecal::reco::StorageScalarType* chi2s,
        char *acState,
        int nchannels,
        int max_iterations,
        uint32_t const offsetForHashes,
        uint32_t const offsetForInputs) {
    // FIXME: ecal has 10 samples and 10 pulses....
    // but this needs to be properly treated and renamed everywhere
    constexpr auto NSAMPLES = SampleMatrix::RowsAtCompileTime;
    constexpr auto NPULSES = SampleMatrix::RowsAtCompileTime;
    static_assert(NSAMPLES == NPULSES);

    // FIXME: remove eitehr idx or ch -> they are teh same thing
    int idx = threadIdx.x + blockDim.x*blockIdx.x;
    auto const ch = idx;
    if (idx < nchannels) {
        if (static_cast<MinimizationState>(acState[idx]) == 
            MinimizationState::Precomputed)
            return;

        // get the hash
        int const inputCh = ch >= offsetForInputs
            ? ch - offsetForInputs
            : ch;
        auto const* dids = ch >= offsetForInputs
            ? dids_ee
            : dids_eb;
        auto const did = DetId{dids[inputCh]};
        auto const isBarrel = did.subdetId() == EcalBarrel;
        auto const hashedId = isBarrel
            ? hashedIndexEB(did.rawId())
            : offsetForHashes + hashedIndexEE(did.rawId());

        // inits
        int iter = 0;
        int npassive = 0;

        // inits
        SampleDecompLLT covariance_decomposition;
        SampleMatrix inverse_cov;
        SampleVector::Scalar chi2 = 0, chi2_now = 0;

        // loop until ocnverge
        while (true) {
            if (iter >= max_iterations)
                break;

            inverse_cov = noisecov[idx];

            update_covariance(
                pulse_covariance[hashedId],
                inverse_cov,
                bxs[idx],
                amplitudes[idx]);

            // compute actual covariance decomposition
            covariance_decomposition.compute(inverse_cov);
            auto const& matrixL = covariance_decomposition.matrixL();

            //
            // replace eigen
            // prepare input matrices for fnnls
            //SampleMatrix A = covariance_decomposition.matrixL()
            //    .solve(pulse_matrix[idx]);
            //SampleVector b = covariance_decomposition.matrixL()
            //    .solve(samples[idx]);
            SampleMatrix A;
            #pragma unroll
            for (int icol=0; icol<NPULSES; icol++) {
                float reg_b[NSAMPLES];
                float reg_L[NSAMPLES];

                // preload a column and load a column 0 of cholesky
                #pragma unroll
                for (int i=0; i<NSAMPLES; i++) {
                    reg_b[i] = pulse_matrix[idx](i, icol);
                    reg_L[i] = matrixL(i, 0);
                }

                // compute x0 and store it
                auto x_prev = reg_b[0] / reg_L[0];
                A(0, icol) = x_prev;

                // iterate
                #pragma unroll
                for (int iL=1; iL<NSAMPLES; iL++) {
                    // update accum
                    #pragma unroll
                    for (int counter=iL; counter<NSAMPLES; counter++)
                        reg_b[counter] -= x_prev * reg_L[counter];

                    // load the next column of chlesky
                    #pragma unroll
                    for (int counter=iL; counter<NSAMPLES; counter++)
                        reg_L[counter] = matrixL(counter, iL);

                    // ocmpute the next x for M(iL, icol)
                    x_prev = reg_b[iL] / reg_L[iL];

                    // store the result value
                    A(iL, icol) = x_prev;
                }
            }

            float reg_b[NSAMPLES];
            {
                float reg_b_tmp[NSAMPLES];
                float reg_L[NSAMPLES];

                // preload a column and load a column 0 of cholesky
                #pragma unroll
                for (int i=0; i<NSAMPLES; i++) {
                    reg_b_tmp[i] = samples[idx](i);
                    reg_L[i] = matrixL(i, 0);
                }

                // compute x0 and store it
                auto x_prev = reg_b_tmp[0] / reg_L[0];
                reg_b[0] = x_prev;

                // iterate
                #pragma unroll
                for (int iL=1; iL<NSAMPLES; iL++) {
                    // update accum
                    #pragma unroll
                    for (int counter=iL; counter<NSAMPLES; counter++)
                        reg_b_tmp[counter] -= x_prev * reg_L[counter];

                    // load the next column of chlesky
                    #pragma unroll
                    for (int counter=iL; counter<NSAMPLES; counter++)
                        reg_L[counter] = matrixL(counter, iL);

                    // ocmpute the next x for M(iL, icol)
                    x_prev = reg_b_tmp[iL] / reg_L[iL];

                    // store the result value
                    reg_b[iL] = x_prev;
                }
            }

            SampleMatrix AtA;
            SampleVector Atb;
            #pragma unroll
            for (int icol=0; icol<NPULSES; icol++) {
                float reg_ai[NSAMPLES];

                // load column icol
                #pragma unroll
                for (int counter=0; counter<NSAMPLES; counter++)
                    reg_ai[counter] = A(counter, icol);

                // compute diagoanl
                float sum = 0.f;
                #pragma unroll
                for (int counter=0; counter<NSAMPLES; counter++)
                    sum += reg_ai[counter] * reg_ai[counter];

                // store
                AtA(icol, icol) = sum;

                // go thru the other columns
                #pragma unroll
                for (int j=icol+1; j<NPULSES; j++) {
                    // load column j
                    float reg_aj[NSAMPLES];
                    #pragma unroll
                    for (int counter=0; counter<NSAMPLES; counter++)
                        reg_aj[counter] = A(counter, j);

                    // accum
                    float sum = 0.f;
                    #pragma unroll
                    for (int counter=0; counter<NSAMPLES; counter++)
                        sum += reg_aj[counter] * reg_ai[counter];

                    // store
                    AtA(icol, j) = sum;
                    AtA(j, icol) = sum;
                }

                // Atb accum
                float sum_atb = 0.f;
                #pragma unroll
                for (int counter=0; counter<NSAMPLES; counter++)
                    sum_atb += reg_ai[counter] * reg_b[counter];

                // store atb
                Atb(icol) = sum_atb;
            }
            
            inplace_fnnls(
                AtA, Atb, amplitudes[idx],
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

namespace v1 {

void minimization_procedure(
        EventInputDataGPU const& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream) {
    unsigned int totalChannels = eventInputGPU.ebDigis.ndigis
        + eventInputGPU.eeDigis.ndigis;
//    unsigned int threads_min = conf.threads.x;
    // TODO: configure from python
    unsigned int threads_min = configParameters.kernelMinimizeThreads[0];
    unsigned int blocks_min = threads_min > totalChannels
        ? 1
        : (totalChannels + threads_min - 1) / threads_min;
    uint32_t const offsetForHashes = conditions.offsetForHashes;
    uint32_t const offsetForInputs = eventInputGPU.ebDigis.ndigis;
    kernel_minimize<<<blocks_min, threads_min, 0, cudaStream.id()>>>(
        eventInputGPU.ebDigis.ids,
        eventInputGPU.eeDigis.ids,
        scratch.noisecov,
        conditions.pulseCovariances.values,
        scratch.activeBXs,
        scratch.samples,
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.pulse_matrix,
        eventOutputGPU.chi2,
        scratch.acState,
        totalChannels,
        50,
        offsetForHashes,
        offsetForInputs);
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
