#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo_gpu_new.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
#include "CondFormats/EcalObjects/interface/EcalSamplesCorrelation.h"

#include <iostream>
#include <limits>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "cuda.h"

#include "AmplitudeComputationCommonKernels.h"
#include "AmplitudeComputationKernelsV1.h"
#include "AmplitudeComputationKernelsV2.h"
#include "TimeComputationKernels.h"

//#define DEBUG

//#define ECAL_RECO_CUDA_DEBUG

namespace ecal { namespace multifit {
   
void entryPoint(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions, 
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream) {
    using digis_type = std::vector<uint16_t>;
    using dids_type = std::vector<uint32_t>;
    bool const gainSwitchUseMaxSampleEB = true; // accodring to the cpu setup
    bool const gainSwitchUseMaxSampleEE = false;

    uint32_t const offsetForHashes = conditions.offsetForHashes;
    unsigned int totalChannels = eventInputCPU.ebDigis.size() 
        + eventInputCPU.eeDigis.size();
    
    // temporary for recording
    cudaEvent_t start_event;
    cudaEvent_t end_event;
    cudaCheck( cudaEventCreate(&start_event) );
    cudaCheck( cudaEventCreate(&end_event) );

    cudaCheck (cudaEventRecord(start_event, 0) );

    //
    // in what follows we copy eb then ee.
    // offset by size 
    //

    // 
    // copy event data: digis + ids, not really async as vectors have default
    // allocators
    //
    cudaCheck( cudaMemcpyAsync(eventInputGPU.digis, 
               eventInputCPU.ebDigis.data().data(),
               eventInputCPU.eeDigis.data().size() * sizeof(digis_type::value_type),
               cudaMemcpyHostToDevice,
               cudaStream.id()) );
    cudaCheck( cudaMemcpyAsync(eventInputGPU.digis + eventInputCPU.ebDigis.data().size(), 
               eventInputCPU.eeDigis.data().data(),
               eventInputCPU.eeDigis.data().size() * sizeof(digis_type::value_type),
               cudaMemcpyHostToDevice,
               cudaStream.id()));

    cudaCheck( cudaMemcpyAsync(eventInputGPU.ids, 
               eventInputCPU.ebDigis.ids().data(),
               eventInputCPU.ebDigis.ids().size() * sizeof(dids_type::value_type),
               cudaMemcpyHostToDevice,
               cudaStream.id()) );
    cudaCheck (cudaMemcpyAsync(eventInputGPU.ids + eventInputCPU.ebDigis.ids().size(), 
               eventInputCPU.eeDigis.ids().data(),
               eventInputCPU.eeDigis.ids().size() * sizeof(dids_type::value_type),
               cudaMemcpyHostToDevice,
               cudaStream.id()) );

    // 
    // 1d preparation kernel
    //
    unsigned int nchannels_per_block = 32;
    unsigned int threads_1d = 10 * nchannels_per_block;
    unsigned int blocks_1d = threads_1d > 10*totalChannels 
        ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
    int shared_bytes = nchannels_per_block * EcalDataFrame::MAXSAMPLES * (
        sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(char)
        + sizeof(bool)
    );
    std::cout << "nchannels = " << totalChannels << std::endl;
    kernel_prep_1d_and_initialize<<<blocks_1d, threads_1d, 
                                    shared_bytes, cudaStream.id()>>>(
        conditions.pulseShapes.values, 
        scratch.epulses,
        eventInputGPU.digis, 
        eventInputGPU.ids,
        scratch.samples,
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.gainsNoise,
        conditions.pedestals.mean_x1,
        conditions.pedestals.mean_x12,
        conditions.pedestals.rms_x12,
        conditions.pedestals.mean_x6,
        conditions.gainRatios.gain6Over1,
        conditions.gainRatios.gain12Over6,
        scratch.hasSwitchToGain6,
        scratch.hasSwitchToGain1,
        scratch.isSaturated,
        eventOutputGPU.amplitude,
        eventOutputGPU.chi2,
        eventOutputGPU.pedestal,
        eventOutputGPU.flags,
        scratch.acState,
        offsetForHashes,
        gainSwitchUseMaxSampleEB,
        gainSwitchUseMaxSampleEE,
        totalChannels);
    AssertIfError

    //
    // 2d preparation kernel
    //
    int blocks_2d = totalChannels;
    dim3 threads_2d{10, 10};
    kernel_prep_2d<<<blocks_2d, threads_2d, 0, cudaStream.id()>>>(
        conditions.pulseCovariances.values, 
        scratch.pulse_covariances,
        scratch.gainsNoise,
        eventInputGPU.ids,
        conditions.pedestals.rms_x12,
        conditions.pedestals.rms_x6,
        conditions.pedestals.rms_x1,
        conditions.gainRatios.gain12Over6,
        conditions.gainRatios.gain6Over1,
        conditions.samplesCorrelation.EBG12SamplesCorrelation,
        conditions.samplesCorrelation.EBG6SamplesCorrelation,
        conditions.samplesCorrelation.EBG1SamplesCorrelation,
        conditions.samplesCorrelation.EEG12SamplesCorrelation,
        conditions.samplesCorrelation.EEG6SamplesCorrelation,
        conditions.samplesCorrelation.EEG1SamplesCorrelation,
        scratch.noisecov,
        scratch.pulse_matrix,
        scratch.epulses,
        scratch.hasSwitchToGain6,
        scratch.hasSwitchToGain1,
        scratch.isSaturated,
        offsetForHashes);
    AssertIfError
    
    // run minimization kernels
    v1::minimization_procedure(
        eventInputCPU, eventInputGPU, eventOutputGPU,
        scratch, conditions, configParameters, cudaStream);

        /*
    if (conf.runV1)
        v1::minimization_procedure(d_data, h_data, conf);
    else
        v2::minimization_procedure(d_data, h_data, conf);
        */

    /*
    //
    // TODO: this guy can run concurrently with other kernels,
    // there is no dependence on the order of execution
    //
    unsigned int threads_time_init = threads_1d;
    unsigned int blocks_time_init = blocks_1d;
    int sharedBytesInit = 2 * threads_time_init * sizeof(SampleVector::Scalar);
    kernel_time_computation_init<<<blocks_time_init, threads_time_init,
                                   sharedBytesInit, conf.cuStream>>>(
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
        offsetForHashesPlaceholder,
        h_data.sample_mask.getEcalSampleMaskRecordEB(),
        h_data.sample_mask.getEcalSampleMaskRecordEE(),
        totalChannels
    );
    AssertIfError

    // 
    // TODO: small kernel only for EB. It needs to be checked if 
    /// fusing such small kernels is beneficial in here
    //
    // we are running only over EB digis
    // therefore we need to create threads/blocks only for that
    unsigned int const threadsFixMGPA = threads_1d;
    unsigned int const blocksFixMGPA = 
        threadsFixMGPA > 10 * h_data.digisEB->size()
            ? 1
            : (10 * h_data.digisEB->size() + threadsFixMGPA - 1) / threadsFixMGPA;
    kernel_time_compute_fixMGPAslew<<<blocksFixMGPA, threadsFixMGPA, 
                                      0, conf.cuStream>>>(
        d_data.digis_data,
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.useless_sample_values,
        h_data.sample_mask.getEcalSampleMaskRecordEB(),
        totalChannels
    );
    AssertIfError

    //
    // 
    //
    int sharedBytes = EcalDataFrame::MAXSAMPLES * nchannels_per_block *
        4 * sizeof(SampleVector::Scalar);
    auto const threads_nullhypot = threads_1d;
    auto const blocks_nullhypot = blocks_1d;
    kernel_time_compute_nullhypot<<<blocks_nullhypot, threads_nullhypot, 
                                    sharedBytes, conf.cuStream>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.useless_sample_values,
        d_data.chi2sNullHypot,
        d_data.sum0sNullHypot,
        d_data.sumAAsNullHypot,
        totalChannels
    );
    AssertIfError

    unsigned int nchannels_per_block_makeratio = 10;
    unsigned int threads_makeratio = 45 * nchannels_per_block_makeratio;
    unsigned int blocks_makeratio = threads_makeratio > 45 * totalChannels
        ? 1
        : (totalChannels * 45 + threads_makeratio - 1) / threads_makeratio;
    int sharedBytesMakeRatio = 5 * threads_makeratio * sizeof(SampleVector::Scalar);
    kernel_time_compute_makeratio<<<blocks_makeratio, threads_makeratio,
                                    sharedBytesMakeRatio, conf.cuStream>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.ids,
        d_data.useless_sample_values,
        d_data.pedestal_nums,
        d_data.amplitudeFitParametersEB,
        d_data.amplitudeFitParametersEE,
        d_data.timeFitParametersEB,
        d_data.timeFitParametersEE,
        d_data.sumAAsNullHypot,
        d_data.sum0sNullHypot,
        d_data.tMaxAlphaBetas,
        d_data.tMaxErrorAlphaBetas,
        d_data.accTimeMax,
        d_data.accTimeWgt,
        d_data.tcState,
        d_data.timeFitParametersSizeEB, 
        d_data.timeFitParametersSizeEE,
        d_data.timeFitLimitsFirstEB,
        d_data.timeFitLimitsFirstEE,
        d_data.timeFitLimitsSecondEB,
        d_data.timeFitLimitsSecondEE,
        totalChannels
    );
    AssertIfError

    //
    //
    //
    auto const threads_findamplchi2 = threads_1d;
    auto const blocks_findamplchi2 = blocks_1d;
    int const sharedBytesFindAmplChi2 = 2 * threads_findamplchi2 * 
        sizeof(SampleVector::Scalar);
    kernel_time_compute_findamplchi2_and_finish<<<blocks_findamplchi2,
                                       threads_findamplchi2,
                                       sharedBytesFindAmplChi2, conf.cuStream>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.ids,
        d_data.useless_sample_values,
        d_data.tMaxAlphaBetas,
        d_data.tMaxErrorAlphaBetas,
        d_data.accTimeMax,
        d_data.accTimeWgt,
        d_data.amplitudeFitParametersEB,
        d_data.amplitudeFitParametersEE,
        d_data.sumAAsNullHypot,
        d_data.sum0sNullHypot,
        d_data.chi2sNullHypot,
        d_data.tcState,
        d_data.ampMaxAlphaBeta,
        d_data.ampMaxError,
        d_data.timeMax,
        d_data.timeError,
        totalChannels
    );
    AssertIfError
        */

    //
    //
    /*
    auto const threads_ampl = threads_1d;
    auto const blocks_ampl = blocks_1d;
    int const sharedBytesAmpl = 5 * threads_ampl * sizeof(SampleVector::Scalar);
    kernel_time_compute_ampl<<<blocks_ampl, threads_ampl,
                               sharedBytesAmpl, conf.cuStream>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.ids,
        d_data.useless_sample_values,
        d_data.timeMax,
        d_data.amplitudeFitParametersEB,
        d_data.amplitudeFitParametersEE,
        d_data.amplitudeMax,
        totalChannels
    );
    AssertIfError
    */
        /*

    //
    //
    //
    auto const threads_timecorr = 32;
    auto const blocks_timecorr = threads_timecorr > totalChannels
        ? 1 : (totalChannels + threads_timecorr-1) / threads_timecorr;
    kernel_time_correction_and_finalize<<<blocks_timecorr, threads_timecorr,
                                          0, conf.cuStream>>>(
        d_data.energies,
        d_data.digis_data,
        d_data.ids,
        d_data.EBTimeCorrAmplitudeBins,
        d_data.EETimeCorrAmplitudeBins,
        d_data.EBTimeCorrShiftBins,
        d_data.EETimeCorrShiftBins,
        d_data.timeMax,
        d_data.timeError,
        d_data.rms_x12,
        d_data.timeCalibConstants,
        d_data.jitter,
        d_data.jitterError,
        d_data.flags,
        h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.size(),
        h_data.time_bias_corrections->EETimeCorrAmplitudeBins.size(),
        d_data.timeConstantTermEB,
        d_data.timeConstantTermEE,
        d_data.offsetTimeValue,
        d_data.timeNconstEB,
        d_data.timeNconstEE,
        d_data.amplitudeThreshEB,
        d_data.amplitudeThreshEE,
        d_data.outOfTimeThreshG12pEB,
        d_data.outOfTimeThreshG12pEE,
        d_data.outOfTimeThreshG12mEB,
        d_data.outOfTimeThreshG12mEE,
        d_data.outOfTimeThreshG61pEB,
        d_data.outOfTimeThreshG61pEE,
        d_data.outOfTimeThreshG61mEB,
        d_data.outOfTimeThreshG61mEE,
        offsetForHashesPlaceholder,
        totalChannels
    );
    AssertIfError

    //
    // transfer eb then ee
    //

    // amplitude
    cudaMemcpyAsync(h_data.rechits_soa_eb.amplitude.data(),
               d_data.energies,
               h_data.rechits_soa_eb.amplitude.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.amplitude.data(),
               d_data.energies + h_data.rechits_soa_eb.amplitude.size(),
               h_data.rechits_soa_ee.amplitude.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // pedestal
    cudaMemcpyAsync(h_data.rechits_soa_eb.pedestal.data(),
               d_data.pedestal,
               h_data.rechits_soa_eb.pedestal.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.pedestal.data(),
               d_data.pedestal + h_data.rechits_soa_eb.pedestal.size(),
               h_data.rechits_soa_ee.pedestal.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // chi2
    cudaMemcpyAsync(h_data.rechits_soa_eb.chi2.data(),
               d_data.chi2,
               h_data.rechits_soa_eb.chi2.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.chi2.data(),
               d_data.chi2 + h_data.rechits_soa_eb.chi2.size(),
               h_data.rechits_soa_ee.chi2.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // detector ids
    cudaMemcpyAsync(h_data.rechits_soa_eb.did.data(),
               d_data.ids,
               h_data.rechits_soa_eb.did.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.did.data(),
               d_data.ids + h_data.rechits_soa_eb.did.size(),
               h_data.rechits_soa_ee.did.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // flags
    cudaMemcpyAsync(h_data.rechits_soa_eb.flags.data(),
               d_data.flags,
               h_data.rechits_soa_eb.flags.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.flags.data(),
               d_data.flags + h_data.rechits_soa_eb.flags.size(),
               h_data.rechits_soa_ee.flags.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // jitter
    cudaMemcpyAsync(h_data.rechits_soa_eb.jitter.data(),
               d_data.jitter,
               h_data.rechits_soa_eb.jitter.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.jitter.data(),
               d_data.jitter + h_data.rechits_soa_eb.jitter.size(),
               h_data.rechits_soa_ee.jitter.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // jitter error
    cudaMemcpyAsync(h_data.rechits_soa_eb.jitterError.data(),
               d_data.jitterError,
               h_data.rechits_soa_eb.jitterError.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.jitterError.data(),
               d_data.jitterError + h_data.rechits_soa_eb.jitterError.size(),
               h_data.rechits_soa_ee.jitterError.size() * sizeof(float),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError

    // amplitudes per sample
    cudaMemcpyAsync(h_data.rechits_soa_eb.amplitudesAll.data(),
               d_data.amplitudes,
               h_data.rechits_soa_eb.amplitudesAll.size() * 
               sizeof(::ecal::reco::ComputationScalarType),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
    cudaMemcpyAsync(h_data.rechits_soa_ee.amplitudesAll.data(),
               d_data.amplitudes + 
               h_data.rechits_soa_eb.amplitudesAll.size() / EcalDataFrame::MAXSAMPLES,
               h_data.rechits_soa_ee.amplitudesAll.size() * 
               sizeof(::ecal::reco::ComputationScalarType),
               cudaMemcpyDeviceToHost,
               conf.cuStream);
    AssertIfError
        */

    cudaEventRecord(end_event, 0);
    cudaEventSynchronize(end_event);
    float ms;
    cudaEventElapsedTime(&ms, start_event, end_event);
    std::cout << "elapsed time = " << ms << std::endl;
}

}}
