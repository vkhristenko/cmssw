#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"

#include "inplace_fnnls.h"
#include "AmplitudeComputationKernelsV1.h"

namespace ecal { namespace multifit {

///
/// assume kernel launch configuration is 
/// (MAXSAMPLES * nchannels, blocks)
/// TODO: is there a point to split this kernel further to separate reductions
/// 
__global__
void kernel_prep_1d_and_initialize(EcalPulseShape const* shapes_in,
                    FullSampleVector* shapes_out, 
                    uint16_t const* digis_in,
                    SampleVector* amplitudes,
                    SampleVector* amplitudesForMinimization,
                    SampleGainVector* gainsNoise,
                    SampleGainVector* gainsPedestal,
                    float const* mean_x1,
                    float const* mean_x12,
                    float const* rms_x12,
                    float const* mean_x6,
                    float const* gain6Over1,
                    float const* gain12Over6,
                    bool* hasSwitchToGain6,
                    bool* hasSwitchToGain1,
                    bool* isSaturated,
                    float* energies,
                    float* chi2,
                    char* acState,
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
        auto const rmsForChecking = rms_x12[ch];
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

        //
        // initialization
        //
        amplitudesForMinimization[ch](sample) = 0;

        if (sample == 0) {
            //
            // initialization
            //
            acState[ch] = static_cast<char>(MinimizationState::NotFinished);
            energies[ch] = 0;
            chi2[ch] = 0;

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
                acState[ch] = static_cast<char>(MinimizationState::Precomputed); 
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
            if (hasGainSwitch && gainSwitchUseMaxSample) {
                // thread for sample=0 will access the right guys
                energies[ch] = max_amplitude / shape_value;
                acState[ch] = static_cast<char>(MinimizationState::Precomputed);
            }
            
            // this happens cause sometimes rms_x12 is 0...
            // needs to be checkec why this is the case
            // general case here is that noisecov is a Zero matrix
            if (rmsForChecking == 0) {
                acState[ch] = static_cast<char>(MinimizationState::Precomputed);
                return;
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

__global__
void kernel_permute_results(
        SampleVector *amplitudes,
        BXVectorType const*activeBXs,
        float *energies,
        char const* acState,
        int const nchannels) {
    // constants
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const tx = threadIdx.x + blockIdx.x * blockDim.x;
    int const ch = tx / nsamples;
    int const iii = tx % nsamples; // this is to address activeBXs

    if (ch >= nchannels) return;

    // configure shared memory and cp into it
    extern __shared__ char smem[];
    SampleVector::Scalar* values = reinterpret_cast<SampleVector::Scalar*>(
        smem);
    values[threadIdx.x] = amplitudes[ch](iii);
    __syncthreads();

    // get the sample for this bx
    auto const sample = static_cast<int>(activeBXs[ch](iii)) + 5;
    auto const state = static_cast<MinimizationState>(acState[ch]);

    // store back to global
    amplitudes[ch](sample) = values[threadIdx.x];

    // store sample 5 separately
    // only for the case when minimization was performed
    // not for cases with precomputed amplitudes
    if (sample == 5 && state != MinimizationState::Precomputed)
        energies[ch] = values[threadIdx.x];
}

///
/// Build an Ecal RecHit.
/// TODO: Use SoA data structures on the host directly
/// the reason for removing this from minimize kernel is to isolate the minimize + 
/// again, building an aos rec hit involves strides... -> bad memory access pattern
///
#ifdef RUN_BUILD_AOS_RECHIT
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
#endif

}}
