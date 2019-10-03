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

/*
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
*/

__device__ __forceinline__
void update_covariance(
        SampleVector const& amplitudes,
        EcalPulseCovariance const& pulse_covariance,
        MapSymM<SampleVector::Scalar, SampleVector::RowsAtCompileTime>& inverse_cov) {
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = BXVectorType::RowsAtCompileTime;

    // FIXME: this cna be made more generic... if needed
    static_assert(npulses == nsamples);

    #pragma unroll
    for (unsigned int ipulse=0; ipulse<npulses; ipulse++) {
        auto const amplitude = amplitudes.coeff(ipulse);
        if (amplitude == 0) 
            continue;

        // in general need soi in here
        int bx = ipulse - 5;
        int first_sample_t = std::max(0, bx+3);
        int offset = -3 - bx;

        auto const value_sq = amplitude * amplitude;
        for (int col=first_sample_t; col<nsamples; col++)
            for (int row=col; row<nsamples; row++)
                inverse_cov(row, col) += value_sq *
                    __ldg(&pulse_covariance.covval[row + offset][col + offset]);
        /*
        for (int row=first_sample_t; row<nsamples; row++) {
            for (int col=row; col<nsamples; col++) {
                inverse_cov.coeffRef(row, col) += value_sq * 
                    __ldg(&pulse_covariance.covval[row + offset][col + offset]);
            }
        }*/
    }
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
        SampleMatrix const* __restrict__ noisecov,
        EcalPulseCovariance const* __restrict__ pulse_covariance,
        SampleVector const* samples,
        SampleVector* amplitudes,
        PulseMatrixType const* __restrict__ pulse_matrix, 
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
    using DataType = SampleMatrix::Scalar;
    static_assert(NSAMPLES == NPULSES);

    // configure shared mem
    extern __shared__ char shrmem[];
    DataType *shrMatrixLFnnlsStorage = 
        reinterpret_cast<float*>(shrmem) + MapSymM<DataType, NPULSES>::total * 
        threadIdx.x;
    DataType *shrAtAStorage = 
        reinterpret_cast<float*>(shrmem) + MapSymM<DataType, NPULSES>::total * (
        threadIdx.x + blockDim.x);

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
        SampleVector::Scalar chi2 = 0, chi2_now = 0;

        ColumnVector<NPULSES, int> pulseOffsets;
        #pragma unroll
        for (int i=0; i<NPULSES; i++)
            pulseOffsets(i) = i;

        // loop until ocnverge
        while (true) {
            if (iter >= max_iterations)
                break;

            printf("idx = %d iter = %d\n", idx, iter);

            // initialize the cov matrix
            DataType covarianceMatrixStorage[MapSymM<DataType, NSAMPLES>::total];
            MapSymM<DataType, NSAMPLES> covarianceMatrix{covarianceMatrixStorage};
            #pragma unroll
            for (int col=0; col<NSAMPLES; col++)
                #pragma unroll
                for (int row=col; row<NSAMPLES; row++)
                    covarianceMatrix(row, col) = 
                        __ldg(&noisecov[idx].coeffRef(row, col));
                //inverse_cov.data()[counter] = __ldg(&noisecov[idx].data()[counter]);

            // update based on the computed amplitudes
            if (iter>0)
                update_covariance(
                    amplitudes[idx],
                    pulse_covariance[hashedId],
                    covarianceMatrix);

            // compute actual covariance decomposition
            //covariance_decomposition.compute(inverse_cov);
            DataType matrixLStorage[MapSymM<DataType, NSAMPLES>::total];
            MapSymM<DataType, NSAMPLES> matrixL{matrixLStorage};
            compute_decomposition_unrolled(matrixL, covarianceMatrix);

            // compute L * A = P
            SampleMatrix A;
            solve_forward_subst_matrix(A, pulse_matrix[idx], matrixL);

            // ocmpute L * b = s
            float reg_b[NSAMPLES];
            solve_forward_subst_vector(reg_b, samples[idx], matrixL);

            // 
            MapSymM<DataType, NPULSES> AtA{shrAtAStorage};
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

            MapSymM<DataType, NPULSES> matrixLForFnnls{shrMatrixLFnnlsStorage};

            fnnls(
                AtA,
                Atb,
                amplitudes[idx],
                npassive,
                pulseOffsets,
                matrixLForFnnls,
                1e-11,
                500);
           
            DataType accum[NSAMPLES];
            {
                DataType results[NPULSES];

                #pragma unroll
                for (int counter=0; counter<NPULSES; counter++)
                    results[counter] = amplitudes[idx].coeff(counter);

                // load accum
                #pragma unroll
                for (int counter=0; counter< NSAMPLES; counter++)
                    accum[counter] = -samples[idx].coeff(counter);

                // iterate
                #pragma unroll
                for (int icol=0; icol<NPULSES; icol++) {
                    DataType pm_col[NSAMPLES];

                    #pragma unroll
                    for (int counter=0; counter<NSAMPLES; counter++)
                        pm_col[counter] = __ldg(
                            &pulse_matrix[idx].coeffRef(counter, icol));

                    // accum
                    #pragma unroll
                    for (int counter=0; counter<NSAMPLES; counter++)
                        accum[counter] += results[icol] * pm_col[counter];
                }
            }

            // solve
            {
                DataType reg_b_tmp[NSAMPLES];
                DataType reg_L[NSAMPLES];
                DataType accumSum = 0;

                // preload
                #pragma unroll
                for (int i=0; i<NSAMPLES; i++) {
                    reg_b_tmp[i] = accum[i];
                    reg_L[i] = matrixL(i, 0);
                }

                // compute x0 and store it
                auto x_prev = reg_b_tmp[0] / reg_L[0];
                accumSum += x_prev * x_prev;
                
                // iterate
                #pragma unroll
                for (int iL=1; iL<NSAMPLES; iL) {
                    // update accum
                    #pragma unroll
                    for (int counter=iL; counter<NSAMPLES; counter++)
                        reg_b_tmp[counter] -= x_prev * reg_L[counter];

                    // load next column of cholesky
                    #pragma unroll
                    for (int counter=iL; counter<NSAMPLES; counter++)
                        reg_L[counter] = matrixL(counter, iL);

                    // ocmpute the next x for M(iL, icol)
                    x_prev = reg_b_tmp[iL] / reg_L[iL];

                    // store teh result value
                    accumSum += x_prev * x_prev;
                }

                chi2_now = accumSum;
            }
              
            auto const deltachi2 = chi2_now - chi2;

            // FIXME add a check for rotation
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
    auto const nbytesShared = 2 * threads_min * 
        MapSymM<SampleVector::Scalar, SampleVector::RowsAtCompileTime>::total * 
        sizeof(SampleVector::Scalar);
    kernel_minimize<<<blocks_min, threads_min, nbytesShared, cudaStream.id()>>>(
        eventInputGPU.ebDigis.ids,
        eventInputGPU.eeDigis.ids,
        scratch.noisecov,
        conditions.pulseCovariances.values,
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
}

}

}}
