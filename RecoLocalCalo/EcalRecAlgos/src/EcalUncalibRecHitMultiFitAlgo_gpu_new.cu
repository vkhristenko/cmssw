#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo_gpu_new.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"

#include <iostream>
#include <limits>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/inplace_fnnls.h"
//#include "RecoLocalCalo/EcalRecAlgos/interface/kernel_minimize_cpu_test.h"

#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

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
        int adc = ecal::mgpa::adc(digis_in[tx]);
        int gainId = ecal::mgpa::gainId(digis_in[tx]);
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

__device__
bool use_sample(unsigned int sample_mask, unsigned int sample) {
    return sample_mask & (EcalDataFrame::MAXSAMPLES - (sample + 1));
}

#define RUN_NULLHYPOT
#ifdef RUN_NULLHYPOT
__global__
void kernel_time_compute_nullhypot(SampleVector::Scalar const* sample_values,
                                   SampleVector::Scalar const* sample_value_errors,
                                   bool const* useless_sample_values,
                                   SampleVector::Scalar* chi2s,
                                   SampleVector::Scalar* sum0s,
                                   SampleVector::Scalar* sumAAs,
                                   int const nchannels) {
    using ScalarType = SampleVector::Scalar;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int tx = threadIdx.x + blockDim.x*blockIdx.x;
    int ltx = threadIdx.x;
    int ch = tx / nsamples;
    int nchannels_per_block = blockDim.x / nsamples;

    // TODO: make sure that this branch plays nicely with __syncthreads inside
    // can there be a deadlock even if the thread is inactive
    if (ch < nchannels) {
        // 
        int sample = tx % nsamples;

        // shared mem inits
        extern __shared__ char sdata[];
        char* s_sum0 = sdata;
        SampleVector::Scalar* s_sum1 = reinterpret_cast<SampleVector::Scalar*>(
            s_sum0 + nchannels_per_block*nsamples);
        SampleVector::Scalar* s_sumA = s_sum1 + nchannels_per_block*nsamples;
        SampleVector::Scalar* s_sumAA = s_sumA + nchannels_per_block*nsamples;

        // TODO make sure no div by 0
        auto const inv_error = useless_sample_values[tx] 
            ? 0.0 
            : 1.0 / (sample_value_errors[tx] * sample_value_errors[tx]);
        auto const sample_value = sample_values[tx];
        s_sum0[ltx] = useless_sample_values[tx] ? 0 : 1;
        s_sum1[ltx] = inv_error;
        s_sumA[ltx] = sample_value * inv_error;
        s_sumAA[ltx] = sample_value * sample_value * inv_error;
        __syncthreads();

        // 5 threads for [0, 4] samples
        if (sample<5) {
            s_sum0[ltx] += s_sum0[ltx+5];
            s_sum1[ltx] += s_sum1[ltx+5];
            s_sumA[ltx] += s_sumA[ltx+5];
            s_sumAA[ltx] += s_sumAA[ltx+5];
        }
        __syncthreads();

        if (sample<2) {
            // note double counting of sample 3
            s_sum0[ltx] += s_sum0[ltx+2] + s_sum0[ltx+3];
            s_sum1[ltx] += s_sum1[ltx+2] + s_sum1[ltx+3];
            s_sumA[ltx] += s_sumA[ltx+2] + s_sumA[ltx+3];
            s_sumAA[ltx] += s_sumAA[ltx+2] + s_sumAA[ltx+3];
        }
        __syncthreads();

        if (sample == 0) {
            // note, subtract to remove the double counting of sample == 3
            //s_sum0[ltx] += s_sum0[ltx+1] - s_sum0[ltx+3];
            //s_sum1[ltx] += s_sum1[ltx+1] - s_sum1[ltx+3];
            //s_sumA[ltx] += s_sumA[ltx+1] - s_sumA[ltx+3];
            //s_sumAA[ltx] += s_sumAA[ltx+1] - s_sumAA[ltx+3];
            auto const sum0 = s_sum0[ltx] + s_sum0[ltx+1] - s_sum0[ltx+3];
            auto const sum1 = s_sum1[ltx] + s_sum1[ltx+1] - s_sum1[ltx+3];
            auto const sumA = s_sumA[ltx] + s_sumA[ltx+1] - s_sumA[ltx+3];
            auto const sumAA = s_sumAA[ltx] + s_sumAA[ltx+1] - s_sumAA[ltx+3];
            auto const chi2 = sum0>0 
                ? (sumAA - sumA * sumA / sum1) / sum0
                : static_cast<ScalarType>(0);
            chi2s[ch] = chi2;
            sum0s[ch] = sum0;
            sumAAs[ch] = sumAA;
        }
    }
}
#endif

constexpr float fast_expf(float x) { return unsafe_expf<6>(x); }
constexpr float fast_logf(float x) { return unsafe_logf<7>(x); }

#define RUN_MAKERATIO
#ifdef RUN_MAKERATIO
//
// launch ctx parameters are 
// 45 threads per channel, X channels per block, Y blocks
// 45 comes from: 10 samples for i <- 0 to 9 and for j <- i+1 to 9
// TODO: it might be much beter to use 32 threads per channel instead of 45
// to simplify the synchronization
//
__global__
void kernel_time_compute_makeratio(SampleVector::Scalar const* sample_values,
                                   SampleVector::Scalar const* sample_value_errors,
                                   bool const* useless_sample_values,
                                   char const* pedestal_nums,
                                   SampleVector::Scalar const* amplitudeFitParameters,
                                   SampleVector::Scalar const* timeFitParameters,
                                   SampleVector::Scalar const* sumAAsNullHypot,
                                   SampleVector::Scalar const* sum0sNullHypot,
                                   SampleVector::Scalar* tMaxAlphaBetas,
                                   SampleVector::Scalar* tMaxErrorAlphaBetas,
                                   SampleVector::Scalar* g_accTimeMax,
                                   SampleVector::Scalar* g_accTimeWgt,
                                   unsigned int const timeFitParameters_size,
                                   SampleVector::Scalar const timeFitLimits_first,
                                   SampleVector::Scalar const timeFitLimits_second,
                                   int const nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nthreads_per_channel = 45; // n=10, n(n-1)/2
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockDim.x*blockIdx.x;
    int const ch = gtx / nthreads_per_channel;
    int const lch = threadIdx.x / nthreads_per_channel;
    int const ltx = threadIdx.x % nthreads_per_channel;
    int const ch_start = ch*nsamples;
    int const lch_start = lch*nsamples;
    int nchannels_per_block = blockDim.x / nthreads_per_channel;

    extern __shared__ char smem[];
    ScalarType* shr_chi2s = reinterpret_cast<ScalarType*>(smem);
    ScalarType* shr_time_wgt = shr_chi2s + blockDim.x;
    ScalarType* shr_time_max = shr_time_wgt + blockDim.x;
    ScalarType* shrTimeMax = shr_time_max + blockDim.x;
    ScalarType* shrTimeWgt = shrTimeMax + blockDim.x;

    // rmeove inactive threads
    // TODO: need to understand if this is 100% safe in presence of syncthreads
    if (ch >= nchannels) return;

    // map tx -> (sample_i, sample_j)
    int sample_i, sample_j = 0;
    if (ltx>=0 && ltx<=8) {
        sample_i = 0;
        sample_j = 1+ltx;
    } else if (ltx<=16) {
        sample_i = 1;
        sample_j = 2+ltx-9;
    } else if (ltx<=23) {
        sample_i = 2;
        sample_j = 3 + ltx - 17;
    } else if (ltx<=29) {
        sample_i = 3;
        sample_j = 4 + ltx - 24;
    } else if (ltx<=34) {
        sample_i = 4;
        sample_j = 5 + ltx - 30;
    } else if (ltx<=38) {
        sample_i = 5;
        sample_j = 6 + ltx - 35;
    } else if (ltx<=41) {
        sample_i = 6;
        sample_j = 7 + ltx - 39;
    } else if (ltx<=43) {
        sample_i = 7;
        sample_j = 8 + ltx - 42;
    } else if (ltx <= 44) {
        sample_i = 8;
        sample_j = 9;
    } else
        assert(false);

    auto const tx_i = ch_start + sample_i;
    auto const tx_j = ch_start + sample_j;

    //
    // note, given the way we partition the block, with 45 threads per channel
    // we will end up with inactive threads which need to be dragged along
    // through the synching point
    // 
    /*
    bool const condToExit = ch >= nchannels
        ? true
        : useless_sample_values[tx_i] 
          || useless_sample_values[tx_j]
          || sample_values[tx_i]<=1 || sample_values[tx_j]<=1;
          */
    bool const condForUselessSamples = useless_sample_values[tx_i] 
        || useless_sample_values[tx_j]
        || sample_values[tx_i]<=1 || sample_values[tx_j]<=1;

    //
    // see cpu implementation for explanation
    // 
    ScalarType chi2 = std::numeric_limits<ScalarType>::max();
    ScalarType tmax = 0;
    ScalarType tmaxerr = 0;
    if (!condForUselessSamples) {
        auto const rtmp = sample_values[tx_i] / sample_values[tx_j];
        auto const invampl_i = 1.0 / sample_values[tx_i];
        auto const relErr2_i = sample_value_errors[tx_i]*sample_value_errors[tx_i]*
            invampl_i*invampl_i;
        auto const invampl_j = 1.0 / sample_values[tx_j];
        auto const relErr2_j = sample_value_errors[tx_j]*sample_value_errors[tx_j]*
            invampl_j*invampl_j;
        auto const err1 = rtmp * rtmp * (relErr2_i + relErr2_j);
        auto err2 = sample_value_errors[tx_j]*
            (sample_values[tx_i] - sample_values[tx_j])*(invampl_j*invampl_j);
        // TODO non-divergent branch for a block if each block has 1 channel
        // otherwise non-divergent for groups of 45 threads
        // at this point, pedestal_nums[ch] can be either 0, 1 or 2
        if (pedestal_nums[ch]==2)
            err2 *= 0.5;
        auto const err3 = (0.289*0.289) * (invampl_j*invampl_j);
        auto const total_error = std::sqrt(err1 + err2 + err3);

        auto const alpha = amplitudeFitParameters[0];
        auto const beta = amplitudeFitParameters[1];
        auto const alphabeta = alpha * beta;
        auto const invalphabeta = 1.0 / alphabeta;

        // variables instead of a struct
        auto const ratio_index = tx_i;
        auto const ratio_step = tx_j - tx_i;
        auto const ratio_value = rtmp;
        auto const ratio_error = total_error;

        //
        // precompute.
        // in cpu version this was done conditionally
        // however easier to do it here (precompute) and then just filter out
        // if not needed
        // 
        auto const l_timeFitLimits_first = timeFitLimits_first;
        auto const l_timeFitLimits_second = timeFitLimits_second;
        if (ratio_step == 1
            && ratio_value >= l_timeFitLimits_first
            && ratio_value <= l_timeFitLimits_second) {

            auto const time_max_i = static_cast<ScalarType>(ratio_index);
            auto u = timeFitParameters[timeFitParameters_size - 1];
#pragma unroll
            for (int k=timeFitParameters_size-2; k>=0; k--)
                u = u*ratio_value + timeFitParameters[k];

            auto du = (timeFitParameters_size - 1) *
                (timeFitParameters[timeFitParameters_size - 1]);
            for (int k=timeFitParameters_size - 2; k>=1; k--)
                du = du*ratio_value + k*timeFitParameters[k];

            auto const error2 = ratio_error * ratio_error * du * du;
            auto const time_max = error2 > 0
                ? (time_max_i - u) / error2
                : static_cast<ScalarType>(0);
            auto const time_wgt = error2 > 0
                ? 1.0 / error2
                : static_cast<ScalarType>(0);

            // store into shared mem
            // note, this name is essentially identical to the one used 
            // below. 
            shrTimeMax[threadIdx.x] = error2 > 0 ? time_max : 0;
            shrTimeWgt[threadIdx.x] = error2 > 0 ? time_wgt : 0;
        } else {
            shrTimeMax[threadIdx.x] = 0;
            shrTimeWgt[threadIdx.x] = 0;
        }

        // continue with ratios
        auto const stepOverBeta = static_cast<SampleVector::Scalar>(ratio_step) / beta;
        auto const offset = static_cast<SampleVector::Scalar>(ratio_index) + alphabeta;
        auto const rmin = std::max(ratio_value - ratio_error, 0.001);
        auto const rmax = std::min(ratio_value + ratio_error, 
            fast_expf(static_cast<SampleVector::Scalar>(ratio_step) / beta)
            - 0.001);
        auto const time1 = 
            offset - 
            ratio_step / 
                (fast_expf((stepOverBeta - fast_logf(rmin)) / 
                                   alpha) - 1.0);
        auto const time2 = 
            offset - 
            ratio_step /
                (fast_expf((stepOverBeta - fast_logf(rmax)) / 
                                   alpha) - 1.0);

        // set these guys
        auto tmax = 0.5 * (time1 + time2);
        auto tmaxerr = 0.5 * std::sqrt((time1 - time2) * (time1 - time2));

        SampleVector::Scalar sumAf = 0;
        SampleVector::Scalar sumff = 0;
        int const itmin = std::max(-1, static_cast<int>(std::floor(tmax - alphabeta)));
        auto loffset = (static_cast<ScalarType>(itmin) - tmax) * invalphabeta;
        // TODO: data dependence 
        for (int it = itmin+1; it<nsamples; it++) {
            loffset += invalphabeta;
            if (useless_sample_values[ch_start + it])
                continue;
            auto const inverr2 = 1.0 / 
                (sample_value_errors[ch_start + it]*sample_value_errors[ch_start + it]);
            auto const term1 = 1.0 + loffset;
            auto const f = (term1 > 1e-6)
                ? fast_expf(alpha * (fast_logf(term1) - loffset))
                : 0;
            sumAf += sample_values[ch_start+it] * (f * inverr2);
            sumff += f*(f*inverr2);
        }

        auto const sumAA = sumAAsNullHypot[ch];
        auto const sum0 = sum0sNullHypot[ch];
        chi2 = sumAA;
        ScalarType amp = 0;
        // TODO: sum0 can not be 0 below, need to introduce the check upfront
        if (sumff > 0) {
            chi2 = sumAA - sumAf * (sumAf / sumff);
            amp = sumAf / sumff;
        }
        chi2 /= sum0;
    }

    // store into smem
    shr_chi2s[threadIdx.x] = chi2;
    __syncthreads();

    // find min chi2 - quite crude for now
    // TODO validate/check
    char iter = nthreads_per_channel / 2 + nthreads_per_channel % 2;
#pragma unroll
    while (iter>=1) {
        if (ltx < iter)
            // for odd ns, the last guy will just store itself
            // exception is for ltx == 0 and iter==1
            shr_chi2s[threadIdx.x] = iter%2==1 && (ltx==iter-1 && ltx>0)
                ? shr_chi2s[threadIdx.x] 
                : std::min(shr_chi2s[threadIdx.x], shr_chi2s[threadIdx.x+iter]);
        __syncthreads();
        iter = iter==1 ? iter/2 : iter/2 + iter%2;
    }

    // filter out inactive or useless samples threads
    if (!condForUselessSamples) {
        // min chi2, now compute weighted average of tmax measurements
        // see cpu version for more explanation
        auto const chi2min = shr_chi2s[threadIdx.x - ltx];
        auto const chi2Limit = chi2min + 1.0;
        auto const inverseSigmaSquared = 
            chi2 < chi2Limit
                ? 1.0 / (tmaxerr * tmaxerr)
                : 0.0;

        // store into shared mem and run reduction
        // TODO: check if cooperative groups would be better
        // TODO: check if shuffling intrinsics are better
        shr_time_wgt[threadIdx.x] = inverseSigmaSquared;
        shr_time_max[threadIdx.x] = tmax * inverseSigmaSquared;
    } else {
        shr_time_wgt[threadIdx.x] = 0;
        shr_time_max[threadIdx.x] = 0;
    }

    // reduce to compute time_max and time_wgt
    iter = nthreads_per_channel / 2 + nthreads_per_channel % 2;
#pragma unroll
    while (iter>=1) {
        if (ltx < iter) {
            shr_time_wgt[threadIdx.x] = iter%2==1 && (ltx==iter-1 && ltx>0)
                ? shr_time_wgt[threadIdx.x]
                : shr_time_wgt[threadIdx.x] + shr_time_wgt[threadIdx.x+iter];
            shr_time_max[threadIdx.x] = iter%2==1 && (ltx==iter-1 && ltx>0)
                ? shr_time_max[threadIdx.x]
                : shr_time_max[threadIdx.x] + shr_time_max[threadIdx.x+iter];
            shrTimeMax[threadIdx.x] = iter%2==1 && (ltx==iter-1 && ltx>0)
                ? shrTimeMax[threadIdx.x]
                : shrTimeMax[threadIdx.x] + shrTimeMax[threadIdx.x+iter];
            shrTimeWgt[threadIdx.x] = iter%2==1 && (ltx==iter-1 && ltx>0)
                ? shrTimeWgt[threadIdx.x]
                : shrTimeWgt[threadIdx.x] + shrTimeWgt[threadIdx.x+iter];
        }
        
        __syncthreads();
        iter = iter==1 ? iter/2 : iter/2 + iter%2;
    }

    // load from shared memory the 0th guy (will contain accumulated values)
    // compute 
    // store into global mem
    if (ltx == 0) {
        auto const tmp_time_max = shr_time_max[threadIdx.x];
        auto const tmp_time_wgt = shr_time_wgt[threadIdx.x];
        auto const tMaxAlphaBeta = tmp_time_max / tmp_time_wgt;
        auto const tMaxErrorAlphaBeta = 1.0 / std::sqrt(tmp_time_wgt);

        tMaxAlphaBetas[ch] = tMaxAlphaBeta;
        tMaxErrorAlphaBetas[ch] = tMaxErrorAlphaBeta;
        g_accTimeMax[ch] = shrTimeMax[threadIdx.x];
        g_accTimeWgt[ch] = shrTimeWgt[threadIdx.x];
    }
}
#endif

/// launch ctx parameters are 
/// 10 threads per channel, N channels per block, Y blocks
/// TODO: do we need to keep the state around or can be removed?!
#define RUN_FINDAMPLCHI2_AND_FINISH
#ifdef RUN_FINDAMPLCHI2_AND_FINISH
__global__
void kernel_time_compute_findamplchi2_and_finish(
        SampleVector::Scalar const* sample_values,
        SampleVector::Scalar const* sample_value_errors,
        bool const* useless_samples,
        SampleVector::Scalar const* g_tMaxAlphaBeta,
        SampleVector::Scalar const* g_tMaxErrorAlphaBeta,
        SampleVector::Scalar const* g_accTimeMax,
        SampleVector::Scalar const* g_accTimeWgt,
        SampleVector::Scalar const* amplitudeFitParameters,
        SampleVector::Scalar const* sumAAsNullHypot,
        SampleVector::Scalar const* sum0sNullHypot,
        SampleVector::Scalar const* chi2sNullHypot,
        TimeComputationState* g_state,
        SampleVector::Scalar* g_ampMaxAlphaBeta,
        SampleVector::Scalar* g_ampMaxError,
        SampleVector::Scalar* g_timeMax,
        SampleVector::Scalar* g_timeError,
        int const nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants 
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x*blockDim.x;
    int const ch = gtx / nsamples;
    int const sample = threadIdx.x % nsamples;
    int const ch_start = ch * nsamples;

    // configure shared mem
    // per block, we need #threads per block * 2 * sizeof(ScalarType)
    // we run with N channels per block
    extern __shared__ char smem[];
    ScalarType* shr_sumAf = reinterpret_cast<ScalarType*>(smem);
    ScalarType* shr_sumff = shr_sumAf + blockDim.x;

    if (ch >= nchannels) return;

    // TODO is that better than storing into global and launching another kernel
    // for the first 10 threads
    auto const alpha = amplitudeFitParameters[0];
    auto const beta = amplitudeFitParameters[1];
    auto const alphabeta = alpha * beta;
    auto const invalphabeta = 1.0 / alphabeta;
    auto const tMaxAlphaBeta = g_tMaxAlphaBeta[ch];
    auto const tMaxErrorAlphaBeta = g_tMaxErrorAlphaBeta[ch];
    auto const sample_value = sample_values[gtx];
    auto const sample_value_error = sample_value_errors[gtx];
    auto const inverr2 = useless_samples[gtx]
        ? static_cast<ScalarType>(0)
        : 1.0 / (sample_value_error * sample_value_error);
    auto const offset = (static_cast<ScalarType>(sample) - tMaxAlphaBeta) 
        * invalphabeta;
    auto const term1 = 1.0 + offset;
    auto const f = term1 > 1e-6 
        ? fast_expf(alpha * (fast_logf(term1) - offset))
        : static_cast<ScalarType>(0.0);
    auto const sumAf = sample_value * (f * inverr2);
    auto const sumff = f * (f * inverr2);
    auto const sumAA = sumAAsNullHypot[ch];
    auto const sum0 = sum0sNullHypot[ch];
    auto const nullChi2 = chi2sNullHypot[ch];

    // store into shared mem
    shr_sumAf[threadIdx.x] = sumAf;
    shr_sumff[threadIdx.x] = sumff;

    // reduce
    // unroll completely here (but hardcoded)
    if (sample<5) {
        shr_sumAf[threadIdx.x] += shr_sumAf[threadIdx.x+5];
        shr_sumff[threadIdx.x] += shr_sumff[threadIdx.x+5];
    }
    __syncthreads();

    if (sample<2) {
        // will need to subtract for ltx = 3, we double count here
        shr_sumAf[threadIdx.x] += shr_sumAf[threadIdx.x+2] 
            + shr_sumAf[threadIdx.x+3];
        shr_sumff[threadIdx.x] += shr_sumff[threadIdx.x+2] 
            + shr_sumff[threadIdx.x+3];
    }
    __syncthreads();

    if (sample==0) {
        // subtract to avoid double counting
        shr_sumAf[threadIdx.x] += shr_sumAf[threadIdx.x+1] 
            - shr_sumAf[threadIdx.x+3];
        shr_sumff[threadIdx.x] += shr_sumAf[threadIdx.x+1] 
            - shr_sumAf[threadIdx.x+3];

        auto const sumff = shr_sumff[threadIdx.x] 
            + shr_sumff[threadIdx.x+1] 
            - shr_sumff[threadIdx.x+3];
        auto const sumAf = shr_sumAf[threadIdx.x]
            + shr_sumAf[threadIdx.x+1]
            - shr_sumAf[threadIdx.x+3];

        TimeComputationState state = TimeComputationState::NotFinished;
        auto const ampMaxAlphaBeta = sumff>0 ? sumAf / sumff : 0;
        if (sumff > 0) {
            auto const chi2AlphaBeta = (sumAA - sumAf * sumAf / sumff) / sum0;
            if (chi2AlphaBeta > nullChi2) {
                // null hypothesis is better
                state = TimeComputationState::Finished;
            }

            // store to global
            g_ampMaxAlphaBeta[ch] = ampMaxAlphaBeta;
        } else {
            state = TimeComputationState::Finished;
        }

        // store the state to global and finish calcs
        g_state[ch] = state;
        if (state == TimeComputationState::Finished) return;

        auto const ampMaxError = g_ampMaxError[ch];
        auto const test_ratio = ampMaxAlphaBeta / ampMaxError;
        auto const accTimeMax = g_accTimeMax[ch];
        auto const accTimeWgt = g_accTimeWgt[ch];
        // branch to separate large vs small pulses
        // see cpu version for more info
        if (test_ratio > 5.0 && accTimeWgt>0) {
            auto const tMaxRatio = accTimeWgt>0 
                ? accTimeMax / accTimeWgt 
                : static_cast<ScalarType>(0);
            auto const tMaxErrorRatio = accTimeWgt>0 
                ? 1.0 / std::sqrt(accTimeWgt) 
                : static_cast<ScalarType>(0);

            if (test_ratio > 10.0) {
                g_timeMax[ch] = tMaxRatio;
                g_timeError[ch] = tMaxErrorRatio;
            } else {
                auto const timeMax = 
                    (tMaxAlphaBeta * (10.0 - ampMaxAlphaBeta / ampMaxError) + 
                     tMaxRatio * (ampMaxAlphaBeta / ampMaxError - 5.0)) / 5.0;
                auto const timeError = 
                    (tMaxErrorAlphaBeta * (10.0 - ampMaxAlphaBeta / ampMaxError) + 
                     tMaxErrorRatio * (ampMaxAlphaBeta / ampMaxError - 5.0)) / 5.0;
                state = TimeComputationState::Finished;
                g_state[ch] = state;
                g_timeMax[ch] = timeMax;
                g_timeError[ch] = timeError;
            }
        }
        else {
            state = TimeComputationState::Finished;
            g_state[ch] = state;
            g_timeMax[ch] = tMaxAlphaBeta;
            g_timeError[ch] = tMaxErrorAlphaBeta;
        }
    }
}
#endif

__global__
void kernel_time_compute_fixMGPAslew(uint16_t const* digis,
                                     SampleVector::Scalar* sample_values,
                                     SampleVector::Scalar* sample_value_errors,
                                     bool* useless_sample_values,
                                     unsigned int const sample_mask,
                                     int const nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x * blockDim.x;
    int const ch = gtx / nsamples;
    int const sample = threadIdx.x % nsamples;

    // remove thread for sample 0, oversubscribing is easier than ....
    if (ch >= nchannels || sample==0) return;

    if (!use_sample(sample_mask, sample)) return;

    auto const gainIdPrev = ecal::mgpa::gainId(digis[gtx-1]);
    auto const gainIdNext = ecal::mgpa::gainId(digis[gtx]);
    if (gainIdPrev>=1 && gainIdPrev<=3 &&
        gainIdNext>=1 && gainIdNext<=3 && gainIdPrev < gainIdNext) {
        sample_values[gtx-1] = 0;
        sample_value_errors[gtx-1] = 1e+9;
        useless_sample_values[gtx-1] = true;
    }
}

#ifdef RUN_AMPL
__global__
void kernel_time_compute_ampl() {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr ScalarType corr4 = 1.;
    constexpr ScalarType corr6 = 1.;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x * blockDim.x;
    int const ch = gtx / nsamples;
    int const sample = threadIdx.x % nsamples;

    auto const alpha = amplitudeFitParameters[0];
    auto const beta = amplitudeFitParameters[1];
    auto const timeMax = g_timeMax[ch];
    auto const pedestalLimit = timeMax - (alpha * beta) - 1.0;
    auto const sample_value = sample_values[gtx];
    auto const sample_value_error = sample_value_errors[gtx];
    auto const inverr2 = 1. / (sample_value_error * sample_value_error);
    auto const termOne = 1 + (sample - timeMax) / (alpha * beta);
    auto const f = termOne > 1.e-5
        ? fast_expf(alpha * fast_logf(termOne) - 
            (sample - timeMax) / beta)
        : static_cast<ScalarType>(0.); 

    bool const cond = ((sample < pedestalLimit) ||
        (f>0.6*corr6 && sample<=timeMax) ||
        (f>0.4*corr4 && sample>=timeMax)) && !useless_samples[gtx];

    // store into shared mem
    shr_sum1[threadIdx.x] = cond ? inverr2 : static_cast<ScalarType>(0);
    shr_sumA[threadIdx.x] = cond 
        ? sample_value * inverr2 
        : static_cast<ScalarType>(0);
    shr_sumF[threadIdx.x] = cond ? f * inverr2 : static_cast<ScalarType>(0);
    shr_sumAF[threadIdx.x] = cond 
        ? (f*inverr2)*sample_value
        : static_cast<ScalarType>(0);
    shr_sumFF[threadIdx.x] = cond 
        ? f*(f*inverr2)
        : static_cast<ScalarType>(0);

    // reduction
    if (sample <= 4) {
        shr_sum1[threadIdx.x] += shr_sum1[threadIdx.x+5];
        shr_sumA[threadIdx.x] += shr_sumA[threadIdx.x+5];
        shr_sumF[threadIdx.x] += shr_sumF[threadIdx.x+5];
        shr_sumAF[threadIdx.x] += shr_sumAF[threadIdx.x+5];
        shr_sumFF[threadIdx.x] += shr_sumFF[threadIdx.x+5];
    }
    __syncthreads();

    if (sample < 2) {
        // note: we double count sample 3
        shr_sum1[threadIdx.x] += shr_sum1[threadIdx.x+2] + shr_sum1[threadIdx.x+3];
        shr_sumA[threadIdx.x] += shr_sumA[threadIdx.x+2] + shr_sumA[threadIdx.x+3];
        shr_sumF[threadIdx.x] += shr_sumF[threadIdx.x+2] + shr_sumF[threadIdx.x+3];
        shr_sumAF[threadIdx.x] += shr_sumAF[threadIdx.x+2] 
            + shr_sumAF[threadIdx.x+3];
        shr_sumFF[threadIdx.x] += shr_sumFF[threadIdx.x+2] 
            + shr_sumFF[threadIdx.x+3];
    }
    __synchtreads();

    if (sample == 0) {
        auto const sum1 = shr_sum1[threadIdx.x] 
            + shr_sum1[threadIdx.x+1] - shr_sum1[threadIdx.x+3];
        auto const sumA = shr_sumA[threadIdx.x] 
            + shr_sumA[threadIdx.x+1] - shr_sumA[threadIdx.x+3];
        auto const sumF = shr_sumF[threadIdx.x] 
            + shr_sumF[threadIdx.x+1] - shr_sumF[threadIdx.x+3];
        auto const sumAF = shr_sumAF[threadIdx.x] 
            + shr_sumAF[threadIdx.x+1] - shr_sumAF[threadIdx.x+3];
        auto const sumFF = shr_sumFF[threadIdx.x] 
            + shr_sumFF[threadIdx.x+1] - shr_sumFF[threadIdx.x+3];

        auto const denom = sumFF * sum1 - sumF*sumF;
        auto const cond = sum1 > 0 && ecal::abs(denom)>1.e-20;
        auto const amplitudeMax = cond
            ? (sumAF * sum1 - sumA * sumF) / denom
            : static_cast<ScalarType>(0.);

        // store into global mem
        g_amplitudeMax[ch] = amplitudeMax;
    }
}
#endif

#ifdef RUN_TMAXRATIOS
__global__
void kernel_time_compute_tmaxratios(
        g_ratio_step) {
    using ScalarType = SampleVector::Scalar;
    
    // constants
    constexpr int nthreads_per_channel = 45; // n=10, n(n-1)/2
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockDim.x*blockIdx.x;
    int const ch = gtx / nthreads_per_channel;
    int const lch = threadIdx.x / nthreads_per_channel;
    int const lthread = threadIdx.x % nthreads_per_channel;
    int const ch_start = ch*nsamples;
    int const lch_start = lch*nsamples;

    if (ch >= nchannels) return;

    // load from global mem
    auto const ratio_step = g_ratio_step[gtx];
    auto const ratio_index = g_ratio_index[gtx];
    auto const ratio_value = g_ratio_value[gtx];
    auto const ratio_error = g_ratio_error[gtx];

    // load time fit limits
    auto const l_timeFitLimits_first = timeFitLimits_first;
    auto const l_timeFitLimits_second = timeFitLimits_second;

    // filter out threads
    if (!(ratio_step == 1
        && ratio_value >= l_timeFitLimits_first
        && ratio_value <= l_timeFitLimits_second)) return;

    // compute
    ScalarType const time_max_i = static_cast<ScalarType>(ratio_index);

    // polynomial for tmax
    // TODO: these are constants... prolly should go into constant memory
    auto u = timeFitParameters[timeFitParameters_size - 1];
#pragma unroll
    for (int k=timeFitParameters_size-2; k>=0; k--)
        u = u*ratio_value + timeFitParameters[k];

    // calculate derivative for tmax error
    auto du = (timeFitParameters_size - 1) *
        (timeFitParameters[timeFitParameters_size - 1]);
    for (int k=timeFitParameters_size - 2; k>=1; k--)
        du = du * ratio_value + k*timeFitParameters[k];

    // running sums for weighted average
    auto const error2 = ratio_error * ratio_error * du * du;
    auto const time_max = error2>0 
        ? (time_max_i - u) / error2
        : static_cast<ScalarType>(0.0);
    auto const time_wgt = error2>0
        ? 1.0 / error2
        : static_cast<ScalarType>(0.0);

    // store in the shr mem
    shr_time_max[threadIdx.x] = time_max;
    shr_time_wgt[threadIdx.x] = time_wgt;

    // reduce sums
    char iter = nthreads_per_channel/2 + nthreads_per_channel%2;
    while (iter>=1) {
        if (lthread < iter) {
            shr_time_wgt[threadIdx.x] = iter%2==1 && lthread==iter-1
                ? shr_time_wgt[threadIdx.x]
                : shr_time_wgt[threadIdx.x] + shr_time_wgt[threadIdx.x+iter];
            shr_time_max[threadIdx.x] = iter%2==1 && lthread==iter-1
                ? shr_time_max[threadIdx.x]
                : shr_time_max[threadIdx.x] + shr_time_max[threadIdx.x+iter];
        }

        iter = iter/2 + iter%2;
        __synchtreads();
    }
    __synchtreads();

    // shared[0] is the accum
    // TODO: this final step could be refactored into a separate kernel
    // as final adjustments are done on a per channel basis and 
    // memory access are also quite aligned for that
    // TODO: verify that we execute this only when this channel finished 
    // flag is not set
    if (lthread == 0 && !g_finished[ch]) {
        auto const acc_time_max = shr_time_max[threadIdx.x];
        auto const acc_time_wgt = shr_time_wgt[threadIdx.x];
        if (acc_time_wgt > 0) {
            auto const tMaxRatio = acc_time_max / acc_time_wgt;
            auto const tMaxErrorRatio = 1.0 / std::sqrt(time_wgt);

            if (ampMaxAlphaBeta / ampMaxError > 10.0) {
                // use pure ratio method
                g_timeMax[ch] = tMaxRatio;
                g_timeError[ch] = tMaxErrorRatio;
            } else {
                // combine two methods
                auto const timeMax = 
                    (tMaxAlphaBeta * (10.0 - ampMaxAlphaBeta / ampMaxError) + 
                     tMaxRatio * (ampMaxAlphaBeta / ampMaxError - 5.0)) / 5.0;
                auto const timeError = 
                    (tMaxErrorAlphaBeta * (10.0 - ampMaxAlphaBeta / ampMaxError) +
                     tMaxErrorRatio * (ampMaxAlphaBeta / ampMaxError - 5.0)) / 5.0;
                g_timeMax[ch] = timeMax;
                g_timeError[ch] = timeError;
            }
        } else {
            g_timeMax[ch] = tMaxAlphaBeta;
            g_timeError[ch] = tMaxErrorAlphaBeta;
        }
    }
}
#endif

__global__
void kernel_time_computation_init(uint16_t const* digis,
                                  uint32_t const* dids,
                                  float const* rms_x12,
                                  float const* rms_x6,
                                  float const* rms_x1,
                                  float const* mean_x12,
                                  float const* mean_x6,
                                  float const* mean_x1,
                                  float const* gain12Over6,
                                  float const* gain6Over1,
                                  SampleVector::Scalar* sample_values,
                                  SampleVector::Scalar* sample_value_errors,
                                  SampleVector::Scalar* ampMaxError,
                                  bool* useless_sample_values,
                                  char* pedestal_nums,
                                  unsigned int const sample_mask,
                                  int nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int tx = threadIdx.x + blockDim.x*blockIdx.x;
    int ch = tx/nsamples;

    if (ch < nchannels) {
        // indices/inits
        int sample = tx % nsamples;
        int ch_start = ch*nsamples;
        SampleVector::Scalar pedestal = 0.;
        int num = 0;

        // configure shared mem
        extern __shared__ char smem[];
        ScalarType* shrSampleValues = 
            reinterpret_cast<SampleVector::Scalar*>(smem);
        ScalarType* shrSampleValueErrors = shrSampleValues + blockDim.x;

        // 0 and 1 sample values
        auto const adc0 = ecal::mgpa::adc(digis[ch_start]);
        auto const gainId0 = ecal::mgpa::gainId(digis[ch_start]);
        auto const adc1 = ecal::mgpa::adc(digis[ch_start+1]);
        auto const gainId1 = ecal::mgpa::gainId(digis[ch_start+1]);

        // set pedestal
        // TODO this branch is non-divergent for a group of 10 threads
        if (gainId0 == 1 && use_sample(sample_mask, 0)) {
            pedestal = static_cast<SampleVector::Scalar>(adc0);
            num=1;

            auto const diff = adc1 - adc0;
            if (gainId1 == 1 && use_sample(sample_mask, 1)
                && std::abs(diff) < 3*rms_x12[ch]) {
                pedestal = 
                    (pedestal + static_cast<SampleVector::Scalar>(adc1)) / 2.0;
                num=2;
            }
        } else {
            pedestal = mean_x12[ch];
        }

        // ped subtracted and gain-renormalized samples.
        auto const gainId = ecal::mgpa::gainId(digis[tx]);
        auto const adc = ecal::mgpa::adc(digis[tx]);

        bool bad = false;
        SampleVector::Scalar sample_value, sample_value_error;
        // TODO divergent branch
        if (!use_sample(sample_mask, sample)) {
            bad = true;
            sample_value = 0;
            sample_value_error = 0;
        } else if (gainId == 1) {
            sample_value = static_cast<SampleVector::Scalar>(adc) - pedestal;
            sample_value_error = rms_x12[ch];
        } else if (gainId == 2) {
            sample_value =  (static_cast<SampleVector::Scalar>(adc) - mean_x6[ch]) 
                * gain12Over6[ch]; 
            sample_value_error = rms_x6[ch] * gain12Over6[ch];
        } else if (gainId == 3) {
            sample_value = (static_cast<SampleVector::Scalar>(adc) - mean_x1[ch])
                * gain6Over1[ch] * gain12Over6[ch];
            sample_value_error = rms_x1[ch] * gain6Over1[ch] * gain12Over6[ch];
        } else {
            sample_value = 0;
            sample_value_error = 0;
            bad = true;
        }

        // TODO: make sure we save things correctly when sample is useless
        useless_sample_values[tx] = (sample_value_error <= 0) | bad;
        sample_values[tx] = sample_value;
        sample_value_errors[tx] = sample_value_error;

        // store into the shared mem
        shrSampleValues[threadIdx.x] = sample_value_error > 0
            ? sample_value
            : std::numeric_limits<ScalarType>::min();
        shrSampleValueErrors[threadIdx.x] = sample_value_error;
        __syncthreads();

        // perform the reduction with min
        if (sample < 5) {
            // note, if equal -> we keep the value with lower sample as for cpu
            shrSampleValueErrors[threadIdx.x] = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+5] 
                ? shrSampleValueErrors[threadIdx.x+5]
                : shrSampleValueErrors[threadIdx.x];
            shrSampleValues[threadIdx.x] = 
                std::max(shrSampleValues[threadIdx.x], 
                         shrSampleValues[threadIdx.x+5]);
        }
        __syncthreads();

        // a bit of an overkill, but easier than to compare across 3 values
        if (sample<3) {
            shrSampleValueErrors[threadIdx.x] = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+3]
                ? shrSampleValueErrors[threadIdx.x+3]
                : shrSampleValueErrors[threadIdx.x];
            shrSampleValues[threadIdx.x] = 
                std::max(shrSampleValues[threadIdx.x], 
                         shrSampleValues[threadIdx.x+3]);
        }
        __syncthreads();

        if (sample < 2) {
            shrSampleValueErrors[threadIdx.x] = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+2]
                ? shrSampleValueErrors[threadIdx.x+2]
                : shrSampleValueErrors[threadIdx.x];
            shrSampleValues[threadIdx.x] = 
                std::max(shrSampleValues[threadIdx.x], 
                         shrSampleValues[threadIdx.x+2]);
        }
        __syncthreads();
 
        if (sample == 0) {
            // we only needd the max error
            auto const maxSampleValueError = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+1]
                ? shrSampleValueErrors[threadIdx.x+1]
                : shrSampleValueErrors[threadIdx.x];

            // # pedestal samples used
            pedestal_nums[ch] = num;
            // this is used downstream
            ampMaxError[ch] = maxSampleValueError;
        }
    }
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
/*
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
*/

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

    //
    // TODO: this guy can run concurrently with other kernels,
    // there is no dependence on the order of execution
    //
    unsigned int threads_time_init = threads_1d;
    unsigned int blocks_time_init = blocks_1d;
    int sharedBytesInit = 2 * threads_time_init * sizeof(SampleVector::Scalar);
    kernel_time_computation_init<<<blocks_time_init, threads_time_init,
                                   sharedBytesInit>>>(
        d_data.digis_data, d_data.ids,
        d_data.rms_x12,
        d_data.rms_x6,
        d_data.rms_x1,
        d_data.mean_x12,
        d_data.mean_x6,
        d_data.mean_x1,
        d_data.gain12Over6,
        d_data.gain6Over1,
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.ampMaxError,
        d_data.useless_sample_values,
        d_data.pedestal_nums,
        barrel 
            ? h_data.sample_mask.getEcalSampleMaskRecordEB()
            : h_data.sample_mask.getEcalSampleMaskRecordEE(),
        h_data.digis->size()
    );
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    // 
    // TODO: small kernel only for EB. It needs to be checked if 
    /// fusing such small kernels is beneficial in here
    //
    if (barrel) {
        kernel_time_compute_fixMGPAslew<<<blocks_time_init, threads_time_init>>>(
            d_data.digis_data,
            d_data.sample_values,
            d_data.sample_value_errors,
            d_data.useless_sample_values,
            h_data.sample_mask.getEcalSampleMaskRecordEB(),
            h_data.digis->size()
        );
        cudaDeviceSynchronize();
        ecal::cuda::assert_if_error();
    }

    //
    // 
    //
    int sharedBytes = EcalDataFrame::MAXSAMPLES * nchannels_per_block *
        4 * sizeof(SampleVector::Scalar);
    auto const threads_nullhypot = threads_1d;
    auto const blocks_nullhypot = blocks_1d;
    kernel_time_compute_nullhypot<<<blocks_nullhypot, threads_nullhypot, 
                                    sharedBytes>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.useless_sample_values,
        d_data.chi2sNullHypot,
        d_data.sum0sNullHypot,
        d_data.sumAAsNullHypot,
        h_data.digis->size()
    );
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    std::cout << "*** before makeratio kernel ***\n";
    //
    // TODO: configurable parameters for launch context below
    //
    unsigned int nchannels_per_block_makeratio = 10;
    unsigned int threads_makeratio = 45 * nchannels_per_block_makeratio;
    unsigned int blocks_makeratio = threads_makeratio > 45 * h_data.digis->size()
        ? 1
        : (h_data.digis->size() * 45 + threads_makeratio - 1) / threads_makeratio;
    int sharedBytesMakeRatio = 5 * threads_makeratio * sizeof(SampleVector::Scalar);
    kernel_time_compute_makeratio<<<blocks_makeratio, threads_makeratio,
                                    sharedBytesMakeRatio>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.useless_sample_values,
        d_data.pedestal_nums,
        barrel ? d_data.amplitudeFitParametersEB : d_data.amplitudeFitParametersEE,
        barrel ? d_data.timeFitParametersEB : d_data.timeFitParametersEE,
        d_data.sumAAsNullHypot,
        d_data.sum0sNullHypot,
        d_data.tMaxAlphaBetas,
        d_data.tMaxErrorAlphaBetas,
        d_data.accTimeMax,
        d_data.accTimeWgt,
        barrel ? d_data.timeFitParametersSizeEB : d_data.timeFitParametersSizeEE,
        barrel ? d_data.timeFitLimitsFirstEB : d_data.timeFitLimitsFirstEE,
        barrel ? d_data.timeFitLimitsSecondEB : d_data.timeFitLimitsSecondEE,
        h_data.digis->size()
    );
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    //
    //
    //
    auto const threads_findamplchi2 = threads_1d;
    auto const blocks_findamplchi2 = blocks_1d;
    int const sharedBytesFindAmplChi2 = 2 * threads_findamplchi2 * 
        sizeof(SampleVector::Scalar);
    kernel_time_compute_findamplchi2_and_finish<<<blocks_findamplchi2,
                                       threads_findamplchi2,
                                       sharedBytesFindAmplChi2>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.useless_sample_values,
        d_data.tMaxAlphaBetas,
        d_data.tMaxErrorAlphaBetas,
        d_data.accTimeMax,
        d_data.accTimeWgt,
        barrel ? d_data.amplitudeFitParametersEB : d_data.amplitudeFitParametersEE,
        d_data.sumAAsNullHypot,
        d_data.sum0sNullHypot,
        d_data.chi2sNullHypot,
        d_data.tcState,
        d_data.ampMaxAlphaBeta,
        d_data.ampMaxError,
        d_data.timeMax,
        d_data.timeError,
        h_data.digis->size()
    );
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
