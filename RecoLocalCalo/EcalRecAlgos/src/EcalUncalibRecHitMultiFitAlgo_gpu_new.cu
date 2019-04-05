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
    
void scatter(host_data& h_data, device_data& d_data, conf_data const& conf) {
    using digis_type = std::vector<uint16_t>;
    using dids_type = std::vector<uint32_t>;
    bool barrel = true;
    bool gainSwitchUseMaxSample = barrel; // accodring to the cpu setup

    unsigned int totalChannels = h_data.digisEB->size() + h_data.digisEE->size();
    
    // temporary for recording
    cudaEvent_t start_event;
    cudaEvent_t end_event;
    cudaEventCreate(&start_event);
    cudaEventCreate(&end_event);

    cudaEventRecord(start_event, 0);

    //
    // in what follows we copy eb then ee.
    // offset by size 
    //

    // 
    // copy digis data
    //
    cudaMemcpy(d_data.digis_data, 
               h_data.digisEB->data().data(),
               h_data.digisEB->data().size() * sizeof(digis_type::value_type),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.digis_data + h_data.digisEB->data().size(), 
               h_data.digisEE->data().data(),
               h_data.digisEE->data().size() * sizeof(digis_type::value_type),
               cudaMemcpyHostToDevice);

    // copy ids
    cudaMemcpy(d_data.ids, 
               h_data.digisEB->ids().data(),
               h_data.digisEB->ids().size() * sizeof(dids_type::value_type),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.ids + h_data.digisEB->ids().size(), 
               h_data.digisEE->ids().data(),
               h_data.digisEE->ids().size() * sizeof(dids_type::value_type),
               cudaMemcpyHostToDevice);

    // copy pedestal means/rmss
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

    // copy gains
    cudaMemcpy(d_data.gain12Over6, h_data.gainratio_data.gain12Over6.data(),
        h_data.gainratio_data.gain12Over6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.gain6Over1, h_data.gainratio_data.gain6Over1.data(),
        h_data.gainratio_data.gain6Over1.size() * sizeof(float),
        cudaMemcpyHostToDevice);

    // copy pulse shape information
    cudaMemcpy(d_data.pulses, h_data.pulse_shapes->data(),
        h_data.pulse_shapes->size() * sizeof(EcalPulseShape),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.covariances, h_data.pulse_covariances->data(),
        h_data.pulse_covariances->size() * sizeof(EcalPulseCovariance),
        cudaMemcpyHostToDevice);

    // copy sample correlations
    cudaMemcpy(d_data.G12SamplesCorrelation, 
               barrel
                 ? h_data.noiseCovariances->EBG12SamplesCorrelation.data()
                 : h_data.noiseCovariances->EEG12SamplesCorrelation.data(),
               EcalDataFrame::MAXSAMPLES * sizeof(double),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.G6SamplesCorrelation, 
               barrel
                 ? h_data.noiseCovariances->EBG6SamplesCorrelation.data()
                 : h_data.noiseCovariances->EEG6SamplesCorrelation.data(),
               EcalDataFrame::MAXSAMPLES * sizeof(double),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.G1SamplesCorrelation, 
               barrel
                 ? h_data.noiseCovariances->EBG1SamplesCorrelation.data()
                 : h_data.noiseCovariances->EEG1SamplesCorrelation.data(),
               EcalDataFrame::MAXSAMPLES * sizeof(double),
               cudaMemcpyHostToDevice);

    // copy time bias corrections
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
   
    // copy bxs
    cudaMemcpy(d_data.bxs, h_data.bxs,
        sizeof(BXVectorType),
        cudaMemcpyHostToDevice);
    AssertIfError

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
    kernel_prep_1d_and_initialize<<<blocks_1d, threads_1d, shared_bytes>>>(
        d_data.pulses, d_data.epulses,
        d_data.digis_data, d_data.samples,
        d_data.amplitudes,
        d_data.gainsNoise,
        d_data.gainsPedestal,
        d_data.mean_x1,
        d_data.mean_x12,
        d_data.rms_x12,
        d_data.mean_x6,
        d_data.gain6Over1,
        d_data.gain12Over6,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.isSaturated,
        d_data.energies,
        d_data.chi2,
        d_data.pedestal,
        d_data.flags,
        d_data.acState,
        gainSwitchUseMaxSample,
        totalChannels);
    AssertIfError

    //
    // 2d preparation kernel
    //
    int blocks_2d = totalChannels;
    dim3 threads_2d{10, 10};
    kernel_prep_2d<<<blocks_2d, threads_2d>>>(
        d_data.covariances, d_data.pulse_covariances,
        d_data.gainsNoise,
        d_data.rms_x12,
        d_data.rms_x6,
        d_data.rms_x1,
        d_data.gain12Over6,
        d_data.gain6Over1,
        d_data.G12SamplesCorrelation,
        d_data.G6SamplesCorrelation,
        d_data.G1SamplesCorrelation,
        d_data.noisecov,
        d_data.pulse_matrix,
        d_data.epulses,
        d_data.bxs,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.isSaturated);
    AssertIfError

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

#endif
    
    if (conf.runV1)
        v1::minimization_procedure(d_data, h_data, conf);
    else
        v2::minimization_procedure(d_data, h_data);
    AssertIfError

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
        totalChannels
    );
    AssertIfError

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
            totalChannels
        );
        AssertIfError
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
        d_data.tcState,
        barrel ? d_data.timeFitParametersSizeEB : d_data.timeFitParametersSizeEE,
        barrel ? d_data.timeFitLimitsFirstEB : d_data.timeFitLimitsFirstEE,
        barrel ? d_data.timeFitLimitsSecondEB : d_data.timeFitLimitsSecondEE,
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
        totalChannels
    );
    AssertIfError

    //
    //
    //
    auto const threads_ampl = threads_1d;
    auto const blocks_ampl = blocks_1d;
    int const sharedBytesAmpl = 5 * threads_ampl * sizeof(SampleVector::Scalar);
    kernel_time_compute_ampl<<<blocks_ampl, threads_ampl,
                               sharedBytesAmpl>>>(
        d_data.sample_values,
        d_data.sample_value_errors,
        d_data.useless_sample_values,
        d_data.timeMax,
        barrel ? d_data.amplitudeFitParametersEB : d_data.amplitudeFitParametersEE,
        d_data.amplitudeMax,
        totalChannels
    );
    AssertIfError

    //
    //
    //
    auto const threads_timecorr = 32;
    auto const blocks_timecorr = threads_timecorr > totalChannels
        ? 1 : (totalChannels + threads_timecorr-1) / threads_timecorr;
    kernel_time_correction_and_finalize<<<blocks_timecorr, threads_timecorr>>>(
        d_data.energies,
        d_data.digis_data,
        barrel ? d_data.EBTimeCorrAmplitudeBins : d_data.EETimeCorrAmplitudeBins,
        barrel ? d_data.EBTimeCorrShiftBins : d_data.EETimeCorrShiftBins,
        d_data.timeMax,
        d_data.timeError,
        d_data.rms_x12,
        d_data.timeCalibConstants,
        d_data.jitter,
        d_data.jitterError,
        d_data.flags,
        barrel 
            ? h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.size() 
            : h_data.time_bias_corrections->EETimeCorrAmplitudeBins.size(),
        barrel 
            ? d_data.timeConstantTermEB
            : d_data.timeConstantTermEE,
        d_data.offsetTimeValue,
        barrel 
            ? d_data.timeNconstEB
            : d_data.timeNconstEE,
        barrel 
            ? d_data.amplitudeThreshEB
            : d_data.amplitudeThreshEE,
        barrel
            ? d_data.outOfTimeThreshG12pEB
            : d_data.outOfTimeThreshG12pEE,
        barrel 
            ? d_data.outOfTimeThreshG12mEB
            : d_data.outOfTimeThreshG12mEE,
        barrel
            ? d_data.outOfTimeThreshG61pEB
            : d_data.outOfTimeThreshG61pEE,
        barrel
            ? d_data.outOfTimeThreshG61mEB
            : d_data.outOfTimeThreshG61mEE,
        totalChannels
    );
    AssertIfError

    //
    // transfer eb then ee
    //

    // amplitude
    cudaMemcpy(h_data.rechits_soa_eb.amplitude.data(),
               d_data.energies,
               h_data.rechits_soa_eb.amplitude.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.amplitude.data(),
               d_data.energies + h_data.rechits_soa_eb.amplitude.size(),
               h_data.rechits_soa_ee.amplitude.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // pedestal
    cudaMemcpy(h_data.rechits_soa_eb.pedestal.data(),
               d_data.pedestal,
               h_data.rechits_soa_eb.pedestal.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.pedestal.data(),
               d_data.pedestal + h_data.rechits_soa_eb.pedestal.size(),
               h_data.rechits_soa_ee.pedestal.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // chi2
    cudaMemcpy(h_data.rechits_soa_eb.chi2.data(),
               d_data.chi2,
               h_data.rechits_soa_eb.chi2.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.chi2.data(),
               d_data.chi2 + h_data.rechits_soa_eb.chi2.size(),
               h_data.rechits_soa_ee.chi2.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // detector ids
    cudaMemcpy(h_data.rechits_soa_eb.did.data(),
               d_data.ids,
               h_data.rechits_soa_eb.did.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.did.data(),
               d_data.ids + h_data.rechits_soa_eb.did.size(),
               h_data.rechits_soa_ee.did.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // flags
    cudaMemcpy(h_data.rechits_soa_eb.flags.data(),
               d_data.flags,
               h_data.rechits_soa_eb.flags.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.flags.data(),
               d_data.flags + h_data.rechits_soa_eb.flags.size(),
               h_data.rechits_soa_ee.flags.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // jitter
    cudaMemcpy(h_data.rechits_soa_eb.jitter.data(),
               d_data.jitter,
               h_data.rechits_soa_eb.jitter.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.jitter.data(),
               d_data.jitter + h_data.rechits_soa_eb.jitter.size(),
               h_data.rechits_soa_ee.jitter.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // jitter error
    cudaMemcpy(h_data.rechits_soa_eb.jitterError.data(),
               d_data.jitterError,
               h_data.rechits_soa_eb.jitterError.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.jitterError.data(),
               d_data.jitterError + h_data.rechits_soa_eb.jitterError.size(),
               h_data.rechits_soa_ee.jitterError.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    AssertIfError

    // amplitudes per sample
    cudaMemcpy(h_data.rechits_soa_eb.amplitudesAll.data(),
               d_data.amplitudes,
               h_data.rechits_soa_eb.amplitudesAll.size() * 
               sizeof(::ecal::reco::ComputationScalarType),
               cudaMemcpyDeviceToHost);
    AssertIfError
    cudaMemcpy(h_data.rechits_soa_ee.amplitudesAll.data(),
               d_data.amplitudes + 
               h_data.rechits_soa_eb.amplitudesAll.size() / EcalDataFrame::MAXSAMPLES,
               h_data.rechits_soa_ee.amplitudesAll.size() * 
               sizeof(::ecal::reco::ComputationScalarType),
               cudaMemcpyDeviceToHost);
    AssertIfError

    cudaEventRecord(end_event, 0);
    cudaEventSynchronize(end_event);
    float ms;
    cudaEventElapsedTime(&ms, start_event, end_event);
    std::cout << "elapsed time = " << ms << std::endl;
}

}}
