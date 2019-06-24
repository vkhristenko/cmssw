#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalSamplesCorrelation.h"

#include "AmplitudeComputationCommonKernels.h"
#include "KernelHelpers.h"

namespace ecal { namespace multifit {

///
/// assume kernel launch configuration is 
/// (MAXSAMPLES * nchannels, blocks)
/// 
// FIXME: add __restrict__
__global__
void kernel_prep_1d_and_initialize(
                    EcalPulseShape const* shapes_in,
                    uint16_t const* digis_in,
                    uint32_t const* dids,
                    SampleVector* amplitudes,
                    SampleVector* amplitudesForMinimization,
                    SampleGainVector* gainsNoise,
                    float const* mean_x1,
                    float const* mean_x12,
                    float const* rms_x12,
                    float const* mean_x6,
                    float const* gain6Over1,
                    float const* gain12Over6,
                    bool* hasSwitchToGain6,
                    bool* hasSwitchToGain1,
                    bool* isSaturated,
                    ::ecal::reco::StorageScalarType* energies,
                    ::ecal::reco::StorageScalarType* chi2,
                    ::ecal::reco::StorageScalarType* g_pedestal,
                    uint32_t *flags,
                    uint32_t *v2rmapping,
                    char *npassive,
                    char *samplesMapping,
                    uint32_t const offsetForHashes,
                    bool const gainSwitchUseMaxSampleEB,
                    bool const gainSwitchUseMaxSampleEE,
                    int const nchannels) {
    constexpr bool dynamicPedestal = false;  //---- default to false, ok
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;
    constexpr int sample_max = 5;
    constexpr int full_pulse_max = 9;
    int const tx = threadIdx.x + blockIdx.x*blockDim.x;
    int const nchannels_per_block = blockDim.x / nsamples;
    int const total_threads = nchannels * nsamples;
    int const ch = tx / nsamples;
    int const sample = threadIdx.x % nsamples;

    if (ch < nchannels) {
        // array of 10 x channels per block
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
        // indices
        //
        auto const did = DetId{dids[ch]};
        auto const isBarrel = did.subdetId() == EcalBarrel;
        // TODO offset for ee, 0 for eb
        auto const hashedId = isBarrel
            ? hashedIndexEB(did.rawId())
            : offsetForHashes + hashedIndexEE(did.rawId());

        // will be used in the future for setting state
        // FIXME: remove these checks
        auto const rmsForChecking = rms_x12[hashedId];

        //
        // amplitudes
        //
        int const adc = ecal::mgpa::adc(digis_in[tx]);
        int const gainId = ecal::mgpa::gainId(digis_in[tx]);
        SampleVector::Scalar amplitude = 0.;
        SampleVector::Scalar pedestal = 0.;
        SampleVector::Scalar gainratio = 0.;

        // store into shared mem for initialization
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
            // => 0 -> will compare 3 and 4 and put into 0
            // => 1 -> will compare 4 and 5 and put into 1
            shr_isSaturated[threadIdx.x] = 
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
                shr_isSaturated[threadIdx.x] || 
                shr_isSaturated[threadIdx.x+1];
            isSaturated[ch] = shr_isSaturated[threadIdx.x+3];
        }

        // TODO: w/o this sync, there is a race
        // if (threadIdx == sample_max) below uses max sample thread, 
        // not for 0 sample                               
        // check if we can remove it
        __syncthreads();

        // TODO: divergent branch
        if (gainId==0 || gainId==3) {
            pedestal = mean_x1[hashedId];
            gainratio = gain6Over1[hashedId] * gain12Over6[hashedId];
            gainsNoise[ch](sample) = 2;
        } else if (gainId==1) {
            pedestal = mean_x12[hashedId];
            gainratio = 1.;
            gainsNoise[ch](sample) = 0;
        } else if (gainId==2) {
            pedestal = mean_x6[hashedId];
            gainratio = gain12Over6[hashedId];
            gainsNoise[ch](sample)  = 1;
        }
        
        // TODO: compile time constant -> branch should be non-divergent
        if (dynamicPedestal)
            amplitude = static_cast<SampleVector::Scalar>(adc) * gainratio;
        else
            amplitude = (static_cast<SampleVector::Scalar>(adc) - pedestal) * gainratio;
        amplitudes[ch][sample] = amplitude;

        //
        // initialization
        //
        amplitudesForMinimization[ch](sample) = 0;
        samplesMapping[ch*nsamples + sample] = sample;

        // select the thread for the max sample 
        //---> hardcoded above to be 5th sample, ok
        if (sample == sample_max) { 
            //
            // initialization
            //
            v2rmapping[ch] = ch;
            energies[ch] = 0;
            chi2[ch] = 0;
            g_pedestal[ch] = 0;
            uint32_t flag = 0;
            npassive[ch] = 0;

            // start of this channel in shared mem
            int const chStart = threadIdx.x - sample_max;
            // thread for the max sample in shared mem
            int const threadMax = threadIdx.x;
            auto const gainSwitchUseMaxSample = isBarrel
                ? gainSwitchUseMaxSampleEB
                : gainSwitchUseMaxSampleEE;
            
            // this flag setting is applied to all of the cases
            if (shr_hasSwitchToGain6[chStart])
                flag |= 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain6;
            if (shr_hasSwitchToGain1[chStart])
                flag |= 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain1;

            // this corresponds to cpu branching on lastSampleBeforeSaturation
            // likely false
            if (check_hasSwitchToGain0) {
                // assign for the case some sample having gainId == 0
                energies[ch] = amplitude;
                amplitudesForMinimization[ch](sample) = amplitude;

                // check if samples before sample_max have true
                bool saturated_before_max = false;
                #pragma unroll
                for (char ii=0; ii<5; ii++)
                    saturated_before_max = saturated_before_max ||
                        shr_hasSwitchToGain0[chStart + ii];

                // if saturation is in the max sample and not in the first 5
                if (!saturated_before_max && 
                    shr_hasSwitchToGain0[threadMax]) {
                    energies[ch] = 49140; // 4095 * 12
                    //---- AM FIXME : no pedestal subtraction???  
                    //It should be "(4095. - pedestal) * gainratio"
                    amplitudesForMinimization[ch](sample) = 49140;
                }

                // set state flag to terminate further processing of this channel
                v2rmapping[ch] = 0xffffffff;
                flag |= 0x1 << EcalUncalibratedRecHit::kSaturated;
                flags[ch] = flag;
                return;
            }

            // according to cpu version
//            auto max_amplitude = amplitudes[ch][sample_max]; 
            auto const max_amplitude = amplitude;
            // according to cpu version
            auto shape_value = shapes_in[hashedId].pdfval[full_pulse_max-7]; 
            // note, no syncing as the same thread will be accessing here
            bool hasGainSwitch = shr_hasSwitchToGain6[chStart]
                || shr_hasSwitchToGain1[chStart]
                || shr_isSaturated[chStart+3];

            // pedestal is final unconditionally
            g_pedestal[ch] = pedestal;
            if (hasGainSwitch && gainSwitchUseMaxSample) {
                // thread for sample=0 will access the right guys
                energies[ch] = max_amplitude / shape_value;
                amplitudesForMinimization[ch](sample) = max_amplitude / shape_value;
                v2rmapping[ch] = 0xffffffff;
                flags[ch] = flag;
                return;
            }
            
            // this happens cause sometimes rms_x12 is 0...
            // needs to be checkec why this is the case
            // general case here is that noisecov is a Zero matrix
            if (rmsForChecking == 0) {
                v2rmapping[ch] = 0xffffffff;
                flags[ch] = flag;
                return;
            }

            // for the case when no shortcuts were taken
            flags[ch] = flag;
        }
    }
}

}}
