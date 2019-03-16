#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo_gpu_new.h"

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

#include "RecoLocalCalo/EcalRecAlgos/interface/inplace_fnnls.h"
//#include "RecoLocalCalo/EcalRecAlgos/interface/kernel_minimize_cpu_test.h"

#include "cuda.h"

//#define DEBUG

//#define ECAL_RECO_CUDA_DEBUG

namespace ecal { namespace multifit { namespace v1 {

///
/// assume kernel launch configuration is 
/// (MAXSAMPLES * nchannels, blocks)
/// TODO: is there a point to split this kernel further to separate reductions
/// 
__global__
void kernel_prep_1d(EcalPulseShape const* shapes_in,
                    FullSampleVector* shapes_out, 
                    uint16_t const* digis_in,
                    SampleVector* amplitudes,
                    SampleGainVector* gainsNoise,
                    SampleGainVector* gainsPedestal,
                    float const* mean_x1,
                    float const* mean_x12,
                    float const* mean_x6,
                    float const* gain6Over1,
                    float const* gain12Over6,
                    bool* hasSwitchToGain6,
                    bool* hasSwitchToGain1,
                    bool* isSaturated,
                    float* energies,
                    char* state_flags,
                    bool gainSwitchUseMaxSample,
                    int nchannels) {
    constexpr bool dynamicPedestal = false;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;
    constexpr int sample_max = 5;
    constexpr int full_pulse_max = 9;
    int tx = threadIdx.x + blockIdx.x*blockDim.x;
    int nchannels_per_block = blockDim.x / nsamples;
    int total_threads = nchannels * nsamples;
    int ch = tx / nsamples;

    if (ch < nchannels) {
        // array of 10 x channels per block
        // TODO: any other way of doing simple reduction
        // assume bool is 1 byte, should be quite safe
        extern __shared__ char shared_mem[];
        bool* shr_hasSwitchToGain6 = reinterpret_cast<bool*>(
            shared_mem);
        bool* shr_hasSwitchToGain1 = shr_hasSwitchToGain6 + 
            nchannels_per_block*nsamples;
        bool* shr_hasSwitchToGain0 = shr_hasSwitchToGain1 + 
            nchannels_per_block*nsamples;
        bool* shr_isSaturated = shr_hasSwitchToGain0 + 
            nchannels_per_block*nsamples;
        bool* shr_hasSwitchToGain0_tmp = shr_isSaturated + 
            nchannels_per_block*nsamples;
        char* shr_counts = reinterpret_cast<char*>(
            shr_hasSwitchToGain0_tmp) + nchannels_per_block*nsamples;
        
        //
        // pulse shape template
        //
        int sample = threadIdx.x % nsamples;
        for (int isample=sample; isample<EcalPulseShape::TEMPLATESAMPLES; 
            isample+=nsamples)
            shapes_out[ch](isample + 7) = shapes_in[ch].pdfval[isample];

        //
        // amplitudes
        //
        int adc = ecalMGPA::adc(digis_in[tx]);
        int gainId = ecalMGPA::gainId(digis_in[tx]);
        SampleVector::Scalar amplitude = 0.;
        SampleVector::Scalar pedestal = 0.;
        SampleVector::Scalar gainratio = 0.;

        shr_hasSwitchToGain6[threadIdx.x] = gainId == EcalMgpaBitwiseGain6;
        shr_hasSwitchToGain1[threadIdx.x] = gainId == EcalMgpaBitwiseGain1;
        shr_hasSwitchToGain0_tmp[threadIdx.x] = gainId == EcalMgpaBitwiseGain0;
        shr_hasSwitchToGain0[threadIdx.x] = shr_hasSwitchToGain0_tmp[threadIdx.x];
        shr_counts[threadIdx.x] = 0;
        __syncthreads();
        
        // non-divergent branch (except for the last 4 threads)
        if (threadIdx.x<=blockDim.x-5) {
            #pragma unroll
            for (int i=0; i<5; i++)
                shr_counts[threadIdx.x] += 
                    shr_hasSwitchToGain0[threadIdx.x+i];
        }
        shr_isSaturated[threadIdx.x] = shr_counts[threadIdx.x] == 5;

        //
        // unrolled reductions
        // TODO
        //
        if (sample < 5) {
            shr_hasSwitchToGain6[threadIdx.x] = 
                shr_hasSwitchToGain6[threadIdx.x] ||
                shr_hasSwitchToGain6[threadIdx.x + 5];
            shr_hasSwitchToGain1[threadIdx.x] =
                shr_hasSwitchToGain1[threadIdx.x] ||
                shr_hasSwitchToGain1[threadIdx.x + 5];
            
            // duplication of hasSwitchToGain0 in order not to
            // introduce another syncthreads
            shr_hasSwitchToGain0_tmp[threadIdx.x] = 
                shr_hasSwitchToGain0_tmp[threadIdx.x] || 
                shr_hasSwitchToGain0_tmp[threadIdx.x+5];
        }
        __syncthreads();
        
        if (sample<2) {
            // note, both threads per channel take value [3] twice to avoid another if
            shr_hasSwitchToGain6[threadIdx.x] = 
                shr_hasSwitchToGain6[threadIdx.x] ||
                shr_hasSwitchToGain6[threadIdx.x+2] || 
                shr_hasSwitchToGain6[threadIdx.x+3];
            shr_hasSwitchToGain1[threadIdx.x] =
                shr_hasSwitchToGain1[threadIdx.x] ||
                shr_hasSwitchToGain1[threadIdx.x+2] || 
                shr_hasSwitchToGain1[threadIdx.x+3];

            shr_hasSwitchToGain0_tmp[threadIdx.x] = 
                shr_hasSwitchToGain0_tmp[threadIdx.x] ||
                shr_hasSwitchToGain0_tmp[threadIdx.x+2] || 
                shr_hasSwitchToGain0_tmp[threadIdx.x+3];

            // sample < 2 -> first 2 threads of each channel will be used here
            // => 0 -> will compare 3 and 4
            // => 1 -> will compare 4 and 5
            shr_isSaturated[threadIdx.x+3] = 
                shr_isSaturated[threadIdx.x+3] || shr_isSaturated[threadIdx.x+4];
        }
        __syncthreads();

        bool check_hasSwitchToGain0 = false;

        if (sample==0) {
            shr_hasSwitchToGain6[threadIdx.x] = 
                shr_hasSwitchToGain6[threadIdx.x] || 
                shr_hasSwitchToGain6[threadIdx.x+1];
            shr_hasSwitchToGain1[threadIdx.x] = 
                shr_hasSwitchToGain1[threadIdx.x] ||
                shr_hasSwitchToGain1[threadIdx.x+1];
            shr_hasSwitchToGain0_tmp[threadIdx.x] =
                shr_hasSwitchToGain0_tmp[threadIdx.x] ||
                shr_hasSwitchToGain0_tmp[threadIdx.x+1];

            hasSwitchToGain6[ch] = shr_hasSwitchToGain6[threadIdx.x];
            hasSwitchToGain1[ch] = shr_hasSwitchToGain1[threadIdx.x];

            // set only for the threadIdx.x corresponding to sample==0
            check_hasSwitchToGain0 = shr_hasSwitchToGain0_tmp[threadIdx.x];

            shr_isSaturated[threadIdx.x+3] = 
                shr_isSaturated[threadIdx.x+3] || 
                shr_isSaturated[threadIdx.x+4];
            isSaturated[ch] = shr_isSaturated[threadIdx.x+3];
        }

        // TODO: divergent branch
        if (gainId==0 || gainId==3) {
            pedestal = mean_x1[ch];
            gainratio = gain6Over1[ch] * gain12Over6[ch];
            gainsNoise[ch](sample) = 2;
            gainsPedestal[ch](sample) = dynamicPedestal ? 2 : -1;
        } else if (gainId==1) {
            pedestal = mean_x12[ch];
            gainratio = 1.;
            gainsNoise[ch](sample) = 0;
            gainsPedestal[ch](sample) = dynamicPedestal ? 0 : -1;
        } else if (gainId==2) {
            pedestal = mean_x6[ch];
            gainratio = gain12Over6[ch];
            gainsNoise[ch](sample)  = 1;
            gainsPedestal[ch](sample) = dynamicPedestal ? 1 : -1;
        }

        // TODO: compile time constant -> branch should be non-divergent
        if (dynamicPedestal)
            amplitude = static_cast<SampleVector::Scalar>(adc) * gainratio;
        else
            amplitude = (static_cast<SampleVector::Scalar>(adc) - pedestal) * gainratio;

        amplitudes[ch][sample] = amplitude;

#ifdef ECAL_RECO_CUDA_DEBUG
        printf("%d %d %d %d %f %f %f\n", tx, ch, sample, adc, amplitude,
            pedestal, gainratio);
        if (adc==0)
            printf("adc is zero\n");
#endif

        if (sample == 0) {

            // this corresponds to cpu branching on lastSampleBeforeSaturation
            // likely false
            if (check_hasSwitchToGain0) {
                // assign for the case some sample having gainId == 0
                energies[ch] = amplitudes[ch][sample_max];

                // check if samples before sample_max have true
                bool saturated_before_max = false;
                #pragma unroll
                for (char ii=0; ii<5; ii++)
                    saturated_before_max = saturated_before_max ||
                        shr_hasSwitchToGain0[threadIdx.x + ii];

                // if saturation is in the max sample and not in the first 5
                if (!saturated_before_max && 
                    shr_hasSwitchToGain0[threadIdx.x + sample_max])
                    energies[ch] = 49140; // 4095 * 12

                // set state flag to terminate further processing of this channel
                state_flags[ch] = 1; 
                return;
            }

            // according to cpu version
            auto max_amplitude = amplitudes[ch][sample_max]; 
            // according to cpu version
            auto shape_value = shapes_out[ch](full_pulse_max); 
            // note, no syncing as the same thread will be accessing here
            bool hasGainSwitch = shr_hasSwitchToGain6[threadIdx.x]
                || shr_hasSwitchToGain1[threadIdx.x]
                || shr_isSaturated[threadIdx.x+3];
            // TODO: divergent branch (groups of 10 will enter)
            if (hasGainSwitch && gainSwitchUseMaxSample) {
                // thread for sample=0 will access the right guys
                energies[ch] = max_amplitude / shape_value;
            }
        }
    }
}

///
/// assume kernel launch configuration is 
/// ([MAXSAMPLES, MAXSAMPLES], nchannels)
///
__global__
void kernel_prep_2d(EcalPulseCovariance const* pulse_cov_in,
                    FullSampleMatrix* pulse_cov_out,
                    SampleGainVector const* gainNoise,
                    SampleMatrixD const* noisecorrs,
                    float const* rms_x12,
                    float const* rms_x6,
                    float const* rms_x1,
                    float const* gain12Over6,
                    float const* gain6Over1,
                    SampleMatrix* noisecov,
                    PulseMatrixType* pulse_matrix,
                    FullSampleVector const* pulse_shape,
                    BXVectorType const* bxs,
                    bool const* hasSwitchToGain6,
                    bool const* hasSwitchToGain1,
                    bool const* isSaturated) {
    int ch = blockIdx.x;
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;
    constexpr float addPedestalUncertainty = 0.f;
    constexpr bool dynamicPedestal = false;
    constexpr bool simplifiedNoiseModelForGainSwitch = true;
    constexpr int template_samples = EcalPulseShape::TEMPLATESAMPLES;

    // only ty == 0 and 1 will go for a second iteration
    for (int iy=ty; iy<template_samples; iy+=nsamples)
        for (int ix=tx; ix<template_samples; ix+=nsamples)
            pulse_cov_out[ch](iy+7, ix+7) = pulse_cov_in[ch].covval[iy][ix];

    /*
    for (int iy=ty, ix=tx; ix<=template_samples && iy<=template_samples; 
        ix+=nsamples, iy+=nsamples)
        pulse_cov_out[ch](iy+7, ix+7) = pulse_cov_in[ch].covval[iy][ix];
        */
    
    bool tmp0 = hasSwitchToGain6[ch];
    bool tmp1 = hasSwitchToGain1[ch];
    bool tmp2 = isSaturated[ch];
    bool hasGainSwitch = tmp0 || tmp1 || tmp2;
    // non-divergent branch for all threads per block
    if (hasGainSwitch) {
        // TODO: did not include simplified noise model
        float noise_value = 0;

        // non-divergent branch - all threads per block
        // TODO: all of these constants indicate that 
        // that these parts could be splitted into completely different 
        // kernels and run one of them only depending on the config
        if (simplifiedNoiseModelForGainSwitch) {
            int isample_max = 5; // according to cpu defs
            int gainidx = gainNoise[ch][isample_max];

            // non-divergent branches
            if (gainidx==0)
                noise_value = rms_x12[ch]*rms_x12[ch]*noisecorrs[0](ty, tx);
            if (gainidx==1) 
                noise_value = gain12Over6[ch]*gain12Over6[ch] * rms_x6[ch]*rms_x6[ch]
                    *noisecorrs[1](ty, tx);
            if (gainidx==2)
                noise_value = gain12Over6[ch]*gain12Over6[ch]
                    * gain6Over1[ch]*gain6Over1[ch] * rms_x1[ch]*rms_x1[ch]
                    * noisecorrs[2](ty, tx);
            if (!dynamicPedestal && addPedestalUncertainty>0.f)
                noise_value += addPedestalUncertainty*addPedestalUncertainty;
        } else {
            // 
            int gainidx=0;
            char mask = gainidx;
            int pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
            noise_value += /* gainratio is 1*/ rms_x12[ch]*rms_x12[ch]
                *pedestal*noisecorrs[0](ty, tx);
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noise_value += /* gainratio is 1 */
                    addPedestalUncertainty*addPedestalUncertainty*pedestal;
            }

            //
            gainidx=1;
            mask = gainidx;
            pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
            noise_value += gain12Over6[ch]*gain12Over6[ch]
                *rms_x6[ch]*rms_x6[ch]*pedestal*noisecorrs[1](ty, tx);
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noise_value += gain12Over6[ch]*gain12Over6[ch]
                    *addPedestalUncertainty*addPedestalUncertainty
                    *pedestal;
            }
            
            //
            gainidx=2;
            mask = gainidx;
            pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
            float tmp = gain6Over1[ch] * gain12Over6[ch];
            noise_value += tmp*tmp * rms_x1[ch]*rms_x1[ch]
                *pedestal*noisecorrs[2](ty, tx);
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noise_value += tmp*tmp * addPedestalUncertainty*addPedestalUncertainty
                    * pedestal;
            }
        }

        noisecov[ch](ty, tx) = noise_value;
    } else {
        auto rms = rms_x12[ch];
        float noise_value = rms*rms * noisecorrs[0](ty, tx);
        if (!dynamicPedestal && addPedestalUncertainty>0.f)
            noise_value += addPedestalUncertainty*addPedestalUncertainty;
        noisecov[ch](ty, tx) = noise_value;
    }

    // pulse matrix
    int bx = (*bxs)(tx);
    int offset = 7 - 3 - bx;
    float value = pulse_shape[ch](offset + ty);
    pulse_matrix[ch](ty, tx) = value;
}

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
bool update_covariance(SampleMatrix const& noisecov,
                       FullSampleMatrix const& full_pulse_cov,
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
        int offset = 7 - 3 - bx;

        auto const value = amplitudes.coeff(ipulse);
        auto const value_sq = value*value;

        unsigned int nsample_pulse = nsamples - first_sample_t;
        inverse_cov.block(first_sample_t, first_sample_t, 
                          nsample_pulse, nsample_pulse)
            += value_sq * full_pulse_cov.block(first_sample_t + offset,
                                               first_sample_t + offset,
                                               nsample_pulse,
                                               nsample_pulse);
    }

    covariance_decomposition.compute(inverse_cov);
    return true;
}

__device__
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
                     float const* rms_x12,
                     char *state_flags,
                     int nchannels,
                     int max_iterations, 
                     bool gainSwitchUseMaxSample) {
    int idx = threadIdx.x + blockDim.x*blockIdx.x;
    if (idx < nchannels) {
        bool hasGainSwitch = isSaturated[idx] 
            || hasSwitchToGain6[idx]
            || hasSwitchToGain1[idx];
        // if should terminate execution
        if (state_flags[idx] == 1)
            return;

        // TODO: gainSwitchUseMaxSimple depends on eb/ee
        // in principle can be splitted/removed into different kernels
        // for ee non-divergent branch
        if (hasGainSwitch && gainSwitchUseMaxSample)
            return;

        // inits
        bool status = false;
        int iter = 0;
        int npassive = 0;
        amplitudes[idx] = SampleVector::Zero();

        //
        // TODO:
        // this happens cause sometimes rms_x12 is 0...
        // needs to be checkec why this is the case
        // general case here is that noisecov is a Zero matrix
        // 
        if (rms_x12[idx] == 0) {
            energies[idx] = 0;
            statuses[idx] = true;
            chi2s[idx] = 0;
            return;
        }
       
        // inits
        SampleDecompLLT covariance_decomposition;
        SampleMatrix inverse_cov;
        SampleVector::Scalar chi2 = 0, chi2_now = 0;


        // TODO 
        BXVectorType activeBXs = *bxs;
        permutation_t permutation;
        permutation.setIdentity();

        // loop until ocnverge
        while (true) {
            if (iter >= max_iterations)
                break;

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
            
            status = inplace_fnnls(
                A, b, amplitudes[idx],
                npassive, activeBXs, permutation, pulse_matrix[idx]);
                
            if (!status)
                break;

            // TODO
            chi2_now = compute_chi2(
                covariance_decomposition,
                pulse_matrix[idx],
                amplitudes[idx],
                samples[idx]);
            auto deltachi2 = chi2_now - chi2;
            chi2 = chi2_now;
            ++iter;

            /*
            if (iter > 10) {
                printf("%d %d %f\n", threadIdx.x, iter, chi2);
            }*/

#ifdef ECAL_RECO_CUDA_DEBUG
            if (iter >= 0 && amplitudes[idx].maxCoeff() == 0) {
        printf("%d %d %f %f %d\n", idx, iter, chi2, chi2_now, 
            hasGainSwitch ? 1 : 0);
        auto matrixL = covariance_decomposition.matrixLLT();
        printf("matrixL max=%f min=%f\n", matrixL.maxCoeff(), matrixL.minCoeff());
        PRINT_MATRIX_10x10(matrixL);
        printf("inverse_cov max=%f min=%f\n", inverse_cov.maxCoeff(), 
                                              inverse_cov.minCoeff());
        PRINT_MATRIX_10x10(inverse_cov);
        printf("noise cov max=%f min=%f\n", noisecov[idx].maxCoeff(),
                                            noisecov[idx].minCoeff());
        PRINT_MATRIX_10x10(noisecov[idx]);
            printf("%d sol amplitudes %f %f %f %f %f %f %f %f %f %f\n", 
                idx,
                amplitudes[idx](0),
                amplitudes[idx](1),
                amplitudes[idx](2),
                amplitudes[idx](3),
                amplitudes[idx](4),
                amplitudes[idx](5),
                amplitudes[idx](6),
                amplitudes[idx](7),
                amplitudes[idx](8),
                amplitudes[idx](9)
            );
            printf("%d samples %f %f %f %f %f %f %f %f %f %f\n", 
                idx,
                samples[idx](0),
                samples[idx](1),
                samples[idx](2),
                samples[idx](3),
                samples[idx](4),
                samples[idx](5),
                samples[idx](6),
                samples[idx](7),
                samples[idx](8),
                samples[idx](9)
            );
            }
#endif

            if (ecal::abs(deltachi2) < 1e-3)
                break;
        }

        amplitudes[idx] = amplitudes[idx].transpose() * permutation.transpose();
        float energy = amplitudes[idx](5);
#ifdef ECAL_RECO_CUDA_DEBUG
        printf("%d %d %f %f %f\n", idx, iter, energy, chi2, chi2_now);
        if (iter==1 && energy==0) {
            printf("%d sol amplitudes %f %f %f %f %f %f %f %f %f %f\n", 
                idx,
                amplitudes[idx](0),
                amplitudes[idx](1),
                amplitudes[idx](2),
                amplitudes[idx](3),
                amplitudes[idx](4),
                amplitudes[idx](5),
                amplitudes[idx](6),
                amplitudes[idx](7),
                amplitudes[idx](8),
                amplitudes[idx](9)
            );
            printf("%d samples %f %f %f %f %f %f %f %f %f %f\n", 
                idx,
                samples[idx](0),
                samples[idx](1),
                samples[idx](2),
                samples[idx](3),
                samples[idx](4),
                samples[idx](5),
                samples[idx](6),
                samples[idx](7),
                samples[idx](8),
                samples[idx](9)
            );
        }
#endif
        energies[idx] = energy; // according to bxs vector bxs[5] = 0
        statuses[idx] = status;
        chi2s[idx] = chi2;
    }
}

///
/// Build an Ecal RecHit.
/// TODO: Use SoA data structures on the host directly
/// the reason for removing this from minimize kernel is to isolate the minimize + 
/// again, building an aos rec hit involves strides... -> bad memory access pattern
///
__global__
void kernel_build_rechit(
    float const* energies,
    float const* chi2s,
    uint32_t* dids,
    EcalUncalibratedRecHit* rechits,
    int nchannels) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < nchannels) {
        rechits[idx] = EcalUncalibratedRecHit{dids[idx], energies[idx],
            0, 0, chi2s[idx], 0};
    }
}

void scatter(host_data& h_data, device_data& d_data, conf_data const& conf) {

/*
void scatter(EcalDigiCollection const& digis,
             EcalUncalibratedRecHitCollection& rechits,
             std::vector<EcalPedestal> const& vpedestals,
             std::vector<EcalMGPAGainRatio> const& vgains,
             std::vector<EcalXtalGroupId> const& vxtals,
             std::vector<EcalPulseShape> const& vpulses,
             std::vector<EcalPulseCovariance> const& vcovariances,
             SampleMatrixGainArray const& noisecors,
             device_data &d_data) {
*/
    auto const& ids = h_data.digis->ids();
    auto const& digis_data = h_data.digis->data();
    using digis_type = std::vector<uint16_t>;
    using dids_type = std::vector<uint32_t>;
    bool barrel = 
        DetId{h_data.digis->begin()->id()}
            .subdetId() == EcalBarrel;
    bool gainSwitchUseMaxSample = barrel; // accodring to the cpu setup

    //
    // TODO: remove per event alloc/dealloc -> do once at the start
    //
    /*
    cudaMalloc((void**)&d_digis_data,
        digis_data.size() * sizeof(digis_type::value_type));
    cudaMalloc((void**)&d_ids,
        ids.size() * sizeof(dids_type::value_type));
    cudaMalloc((void**)&d_pedestals,
        vpedestals.size() * sizeof(EcalPedestal));
    cudaMalloc((void**)&d_gains, 
        vgains.size() * sizeof(EcalMGPAGainRatio));
    cudaMalloc((void**)&d_xtals,
        vxtals.size() * sizeof(EcalXtalGroupId));
    cudaMalloc((void**)&d_shapes,
        vpulses.size() * sizeof(EcalPulseShape));
    cudaMalloc((void**)&d_covariances,
        vcovariances.size() * sizeof(EcalPulseCovariance));
    cudaMalloc((void**)&d_rechits,
        rechits.size() * sizeof(EcalUncalibratedRecHit));
    cudaMalloc((void**)&d_noisecors,
        noisecors.size() * sizeof(SampleMatrix));
    ecal::cuda::assert_if_error();
    */

    //
    // set the state flags
    //
    cudaMemset(d_data.state_flags, 0, h_data.digis->size() * sizeof(char));

    // 
    // copy to the device
    // TODO: can conditions be copied only once when updated?
    //
    cudaMemcpy(d_data.digis_data, digis_data.data(),
        digis_data.size() * sizeof(digis_type::value_type),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.ids, ids.data(),
        ids.size() * sizeof(dids_type::value_type),
        cudaMemcpyHostToDevice);
 
//    cudaMemcpy(d_data.pedestals, h_data.pedestals->data(),
//        h_data.pedestals->size() * sizeof(EcalPedestal),
//        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.mean_x12, h_data.ped_data.mean_x12.data(),
        h_data.ped_data.mean_x12.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.rms_x12, h_data.ped_data.rms_x12.data(),
        h_data.ped_data.rms_x12.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.mean_x6, h_data.ped_data.mean_x6.data(),
        h_data.ped_data.mean_x6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.rms_x6, h_data.ped_data.rms_x6.data(),
        h_data.ped_data.rms_x6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.mean_x1, h_data.ped_data.mean_x1.data(),
        h_data.ped_data.mean_x1.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.rms_x1, h_data.ped_data.rms_x1.data(),
        h_data.ped_data.rms_x1.size() * sizeof(float),
        cudaMemcpyHostToDevice);

//    cudaMemcpy(d_data.gains, h_data.gains->data(),
//        h_data.gains->size() * sizeof(EcalMGPAGainRatio),
//        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.gain12Over6, h_data.gainratio_data.gain12Over6.data(),
        h_data.gainratio_data.gain12Over6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.gain6Over1, h_data.gainratio_data.gain6Over1.data(),
        h_data.gainratio_data.gain6Over1.size() * sizeof(float),
        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.xtals, h_data.xtals->data(),
        h_data.xtals->size() * sizeof(EcalXtalGroupId),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.pulses, h_data.pulse_shapes->data(),
        h_data.pulse_shapes->size() * sizeof(EcalPulseShape),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.covariances, h_data.pulse_covariances->data(),
        h_data.pulse_covariances->size() * sizeof(EcalPulseCovariance),
        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.noisecorrs, h_data.noisecorrs->data(),
        h_data.noisecorrs->size() * sizeof(SampleMatrixD),
        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.sample_mask, h_data.sample_mask,
        sizeof(EcalSampleMask),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EBTimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EBTimeCorrShiftBins, 
        h_data.time_bias_corrections->EBTimeCorrShiftBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EBTimeCorrShiftBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EETimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EETimeCorrAmplitudeBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EETimeCorrAmplitudeBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EETimeCorrShiftBins, 
        h_data.time_bias_corrections->EETimeCorrShiftBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EETimeCorrShiftBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.weights,
        h_data.weights->data(),
        sizeof(EMatrix) * h_data.weights->size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.bxs, h_data.bxs,
        sizeof(BXVectorType),
        cudaMemcpyHostToDevice);
    ecal::cuda::assert_if_error();

    int nthreads_per_block = conf.threads.x;
    int nblocks = (h_data.digis->size() + nthreads_per_block - 1) / nthreads_per_block;

    // 
    // 1d preparation kernel
    //
    unsigned int nchannels_per_block = 32;
    unsigned int threads_1d = 10 * nchannels_per_block;
    unsigned int blocks_1d = threads_1d > 10*h_data.digis->size() 
        ? 1 : (h_data.digis->size()*10 + threads_1d - 1) / threads_1d;
    int shared_bytes = nchannels_per_block * EcalDataFrame::MAXSAMPLES * (
        sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(char)
        + sizeof(bool)
    );
    std::cout << "nchannels = " << h_data.digis->size() << std::endl;
    std::cout << "shared memory per block = " << shared_bytes << "B\n";
    kernel_prep_1d<<<blocks_1d, threads_1d, shared_bytes>>>(
        d_data.pulses, d_data.epulses,
        d_data.digis_data, d_data.samples,
        d_data.gainsNoise,
        d_data.gainsPedestal,
        d_data.mean_x1,
        d_data.mean_x12,
        d_data.mean_x6,
        d_data.gain6Over1,
        d_data.gain12Over6,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.isSaturated,
        d_data.energies,
        d_data.state_flags,
        gainSwitchUseMaxSample,
        h_data.digis->size());
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    std::cout << " after kernel prep 1d\n";

    //
    // 2d preparation kernel
    //
    int blocks_2d = h_data.digis->size();
    dim3 threads_2d{10, 10};
    kernel_prep_2d<<<blocks_2d, threads_2d>>>(
        d_data.covariances, d_data.pulse_covariances,
        d_data.gainsNoise,
        d_data.noisecorrs,
        d_data.rms_x12,
        d_data.rms_x6,
        d_data.rms_x1,
        d_data.gain12Over6,
        d_data.gain6Over1,
        d_data.noisecov,
        d_data.pulse_matrix,
        d_data.epulses,
        d_data.bxs,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.isSaturated);
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    std::cout << "after kernel prep 2d\n";

//#define ECAL_RECO_DEBUG_CPU4GPU
#ifdef ECAL_RECO_DEBUG_CPU4GPU

    // debug quantities before launching minimization
    std::vector<SampleVector> samples(h_data.digis->size());
    std::vector<PulseMatrixType> pulse_matrices(h_data.digis->size());
    std::vector<SampleMatrix> noisecovs(h_data.digis->size());
    std::vector<FullSampleMatrix> pulse_covariances(h_data.digis->size());
    std::vector<SampleVector> amplitudes(h_data.digis->size());
    std::vector<float> energies(h_data.digis->size());
    std::vector<char> statuses(h_data.digis->size());
    std::vector<float> chi2s(h_data.digis->size());
    std::vector<char> isSaturated(h_data.digis->size());
    std::vector<char> hasSwitchToGain6(h_data.digis->size());
    std::vector<char> hasSwitchToGain1(h_data.digis->size());
    cudaMemcpy(samples.data(), d_data.samples,
        h_data.digis->size() * sizeof(SampleVector),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(pulse_matrices.data(), d_data.pulse_matrix,
        pulse_matrices.size() * sizeof(PulseMatrixType),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(noisecovs.data(), d_data.noisecov,
        noisecovs.size() * sizeof(SampleMatrix),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(pulse_covariances.data(), d_data.pulse_covariances,
        pulse_covariances.size() * sizeof(FullSampleMatrix),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(isSaturated.data(), 
        d_data.isSaturated,
        isSaturated.size() * sizeof(bool),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(hasSwitchToGain6.data(), d_data.hasSwitchToGain6,
        hasSwitchToGain6.size() * sizeof(bool),
        cudaMemcpyDeviceToHost);
    cudaMemcpy(hasSwitchToGain1.data(), d_data.hasSwitchToGain1,
        hasSwitchToGain1.size() * sizeof(bool),
        cudaMemcpyDeviceToHost);

    //std::cout << "dumping gpu quantities\n";

    cpu::kernel_minimize(
        noisecovs.data(),
        pulse_covariances.data(),
        h_data.bxs,
        samples.data(),
        amplitudes.data(),
        energies.data(),
        pulse_matrices.data(),
        reinterpret_cast<bool*>(statuses.data()),
        chi2s.data(),
        reinterpret_cast<bool*>(isSaturated.data()),
        reinterpret_cast<bool*>(hasSwitchToGain6.data()),
        reinterpret_cast<bool*>(hasSwitchToGain1.data()),
        h_data.digis->size(),
        50,
        gainSwitchUseMaxSample
    );

    /*
    for (unsigned int i=0; i<samples.size(); i++) {
        auto const& sample_vector = samples[i];
        auto const& pulse_matrix = pulse_matrices[i];
        auto const& noisecov = noisecovs[i];
        auto const& pulse_cov = pulse_covariances[i];

        std::cout << "*** samples ***\n"
            << sample_vector << std::endl;
        std::cout << "*** pulse matrix ***\n"
            << pulse_matrix << std::endl;
        std::cout << "*** noisecov ***\n"
            << noisecov << std::endl;
        std::cout << "*** pulse cov ***\n"
            << pulse_cov << std::endl;
    }
    */

#endif

    cudaEvent_t start_event;
    cudaEvent_t end_event;
    cudaEventCreate(&start_event);
    cudaEventCreate(&end_event);

    unsigned int threads_min = conf.threads.x;
    unsigned int blocks_min = threads_min > h_data.digis->size()
        ? 1 : (h_data.digis->size() + threads_min - 1) / threads_min;
    cudaEventRecord(start_event, 0);
    kernel_minimize<<<blocks_min, threads_min>>>(
        d_data.noisecov,
        d_data.pulse_covariances,
        d_data.bxs,
        d_data.samples,
        d_data.amplitudes,
        d_data.energies,
        d_data.pulse_matrix,
        d_data.statuses,
        d_data.chi2,
        d_data.isSaturated,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.rms_x12,
        d_data.state_flags,
        h_data.digis->size(),
        50,
        gainSwitchUseMaxSample);
    cudaEventRecord(end_event, 0);
    cudaEventSynchronize(end_event);
    float ms;
    cudaEventElapsedTime(&ms, start_event, end_event);
    std::cout << "elapsed time = " << ms << std::endl;
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();
    
/*
    kernel_build_rechit<<<blocks_min, threads_min>>>(
        d_data.energies,
        d_data.chi2,
        d_data.ids,
        d_data.rechits,
        h_data.digis->size());
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();
    */

  /*  kernel_reconstruct<<<nblocks, nthreads_per_block>>>(
        d_data.digis_data,
        d_data.ids,*/
        /* d_rechits, */
/*        d_data.pedestals,
        d_data.gains,
        d_data.xtals,
        d_data.pulses,
        d_data.covariances,
        d_data.rechits,
        d_data.noisecors,
        d_data.sample_mask,
        d_data.EBTimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.size(),
        d_data.EBTimeCorrShiftBins, 
        h_data.time_bias_corrections->EBTimeCorrShiftBins.size(),
        d_data.EETimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EETimeCorrAmplitudeBins.size(),
        d_data.EETimeCorrShiftBins, 
        h_data.time_bias_corrections->EETimeCorrShiftBins.size(),
        d_data.weights,
        h_data.digis->size()
    );*/

    //
    // transfer the results back
    //
//    h_data.rechits_soa.amplitude = std::move(energies);
//    h_data.rechits_soa.chi2 = std::move(chi2s);
    cudaMemcpy(&(*h_data.rechits_soa.amplitude.begin()),
               d_data.energies,
               h_data.rechits_soa.amplitude.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*h_data.rechits_soa.chi2.begin()),
               d_data.chi2,
               h_data.rechits_soa.chi2.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*h_data.rechits_soa.did.begin()),
               d_data.ids,
               h_data.rechits_soa.did.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost);

//    cudaMemcpy(&(*h_data.rechits->begin()), d_data.rechits,
//        h_data.rechits->size() * sizeof(EcalUncalibratedRecHit),
//        cudaMemcpyDeviceToHost);

    // 
    // free all the device ptrs
    // TODO: remove per event dealloc
    //
    /*
    cudaFree(d_digis_data);
    cudaFree(d_ids);
    cudaFree(d_pedestals);
    cudaFree(d_gains);
    cudaFree(d_xtals);
    cudaFree(d_shapes);
    cudaFree(d_covariances);
    cudaFree(d_rechits);
    cudaFree(d_noisecors);
    ecal::cuda::assert_if_error();
    */
}

}}}
