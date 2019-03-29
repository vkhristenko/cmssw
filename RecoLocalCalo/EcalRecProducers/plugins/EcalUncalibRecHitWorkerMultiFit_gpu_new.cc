#include "RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit_gpu_new.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalWeightXtalGroupsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTBWeightsRcd.h"
#include "CondFormats/DataRecord/interface/EcalSampleMaskRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalSamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseCovariancesRcd.h"

#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include <FWCore/ParameterSet/interface/EmptyGroupDescription.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#define MAX_CHANNELS 50000
#define MAX_TIME_BIAS_CORRECTIONS 100

EcalUncalibRecHitWorkerMultiFitGPUNew::EcalUncalibRecHitWorkerMultiFitGPUNew(const edm::ParameterSet&ps,edm::ConsumesCollector& c) :
  EcalUncalibRecHitWorkerBaseClass(ps,c)
{

  // get the BX for the pulses to be activated
  std::vector<int32_t> activeBXs = ps.getParameter< std::vector<int32_t> >("activeBXs");
  activeBX.resize(activeBXs.size());
  for (unsigned int ibx=0; ibx<activeBXs.size(); ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  }

  // uncertainty calculation (CPU intensive)
  ampErrorCalculation_ = ps.getParameter<bool>("ampErrorCalculation");
  useLumiInfoRunHeader_ = ps.getParameter<bool>("useLumiInfoRunHeader");
  
  if (useLumiInfoRunHeader_) {
    bunchSpacing_ = c.consumes<unsigned int>(edm::InputTag("bunchSpacingProducer"));
    bunchSpacingManual_ = 0;
  } else {
    bunchSpacingManual_ = ps.getParameter<int>("bunchSpacing");
  }

  doPrefitEB_ = ps.getParameter<bool>("doPrefitEB");
  doPrefitEE_ = ps.getParameter<bool>("doPrefitEE");

  prefitMaxChiSqEB_ = ps.getParameter<double>("prefitMaxChiSqEB");
  prefitMaxChiSqEE_ = ps.getParameter<double>("prefitMaxChiSqEE");
  
  dynamicPedestalsEB_ = ps.getParameter<bool>("dynamicPedestalsEB");
  dynamicPedestalsEE_ = ps.getParameter<bool>("dynamicPedestalsEE");
  mitigateBadSamplesEB_ = ps.getParameter<bool>("mitigateBadSamplesEB");
  mitigateBadSamplesEE_ = ps.getParameter<bool>("mitigateBadSamplesEE");
  gainSwitchUseMaxSampleEB_ = ps.getParameter<bool>("gainSwitchUseMaxSampleEB");
  gainSwitchUseMaxSampleEE_ = ps.getParameter<bool>("gainSwitchUseMaxSampleEE");
  selectiveBadSampleCriteriaEB_ = ps.getParameter<bool>("selectiveBadSampleCriteriaEB");
  selectiveBadSampleCriteriaEE_ = ps.getParameter<bool>("selectiveBadSampleCriteriaEE");
  addPedestalUncertaintyEB_ = ps.getParameter<double>("addPedestalUncertaintyEB");
  addPedestalUncertaintyEE_ = ps.getParameter<double>("addPedestalUncertaintyEE");
  simplifiedNoiseModelForGainSwitch_ = ps.getParameter<bool>("simplifiedNoiseModelForGainSwitch");
  
  // algorithm to be used for timing
  auto const & timeAlgoName = ps.getParameter<std::string>("timealgo");
  if(timeAlgoName=="RatioMethod") timealgo_=ratioMethod;
  else if(timeAlgoName=="WeightsMethod") timealgo_=weightsMethod;
  else if(timeAlgoName!="None")
   edm::LogError("EcalUncalibRecHitError") << "No time estimation algorithm defined";

  // ratio method parameters
  EBtimeFitParameters_ = ps.getParameter<std::vector<double> >("EBtimeFitParameters"); 
  EEtimeFitParameters_ = ps.getParameter<std::vector<double> >("EEtimeFitParameters"); 
  EBamplitudeFitParameters_ = ps.getParameter<std::vector<double> >("EBamplitudeFitParameters");
  EEamplitudeFitParameters_ = ps.getParameter<std::vector<double> >("EEamplitudeFitParameters");
  EBtimeFitLimits_.first  = ps.getParameter<double>("EBtimeFitLimits_Lower");
  EBtimeFitLimits_.second = ps.getParameter<double>("EBtimeFitLimits_Upper");
  EEtimeFitLimits_.first  = ps.getParameter<double>("EEtimeFitLimits_Lower");
  EEtimeFitLimits_.second = ps.getParameter<double>("EEtimeFitLimits_Upper");
  EBtimeConstantTerm_=ps.getParameter<double>("EBtimeConstantTerm");
  EEtimeConstantTerm_=ps.getParameter<double>("EEtimeConstantTerm");
  EBtimeNconst_=ps.getParameter<double>("EBtimeNconst");
  EEtimeNconst_=ps.getParameter<double>("EEtimeNconst");
  outOfTimeThreshG12pEB_ = ps.getParameter<double>("outOfTimeThresholdGain12pEB");
  outOfTimeThreshG12mEB_ = ps.getParameter<double>("outOfTimeThresholdGain12mEB");
  outOfTimeThreshG61pEB_ = ps.getParameter<double>("outOfTimeThresholdGain61pEB");
  outOfTimeThreshG61mEB_ = ps.getParameter<double>("outOfTimeThresholdGain61mEB");
  outOfTimeThreshG12pEE_ = ps.getParameter<double>("outOfTimeThresholdGain12pEE");
  outOfTimeThreshG12mEE_ = ps.getParameter<double>("outOfTimeThresholdGain12mEE");
  outOfTimeThreshG61pEE_ = ps.getParameter<double>("outOfTimeThresholdGain61pEE");
  outOfTimeThreshG61mEE_ = ps.getParameter<double>("outOfTimeThresholdGain61mEE");
  amplitudeThreshEB_ = ps.getParameter<double>("amplitudeThresholdEB");
  amplitudeThreshEE_ = ps.getParameter<double>("amplitudeThresholdEE");

  // spike threshold
  ebSpikeThresh_ = ps.getParameter<double>("ebSpikeThreshold");

  ebPulseShape_ = ps.getParameter<std::vector<double> >("ebPulseShape");
  eePulseShape_ = ps.getParameter<std::vector<double> >("eePulseShape");

  // chi2 parameters for flags determination
  kPoorRecoFlagEB_ = ps.getParameter<bool>("kPoorRecoFlagEB");
  kPoorRecoFlagEE_ = ps.getParameter<bool>("kPoorRecoFlagEE");;
  chi2ThreshEB_=ps.getParameter<double>("chi2ThreshEB_");
  chi2ThreshEE_=ps.getParameter<double>("chi2ThreshEE_");

  // threads/blocks conf
  auto vthreads = ps.getParameter<std::vector<int>>("threads");
  conf.threads = {vthreads[0], vthreads[1], vthreads[2]};

  // 
  // TODO
  //
  cudaMalloc((void**)&d_data.digis_data,
    MAX_CHANNELS * EcalDataFrame::MAXSAMPLES * sizeof(uint16_t));
  cudaMalloc((void**)&d_data.ids,
    MAX_CHANNELS * sizeof(uint32_t));
  cudaMalloc((void**)&d_data.amplitudes,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::SampleVector));
  cudaMalloc((void**)&d_data.samples,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::SampleVector));
  cudaMalloc((void**)&d_data.gainsNoise,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::SampleGainVector));
  cudaMalloc((void**)&d_data.gainsPedestal,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::SampleGainVector));

//  cudaMalloc((void**)&d_data.pedestals,
//    MAX_CHANNELS * sizeof(EcalPedestal));

  cudaMalloc((void**)&d_data.mean_x12,
    MAX_CHANNELS * sizeof(float));
  cudaMalloc((void**)&d_data.rms_x12,
    MAX_CHANNELS * sizeof(float));
  cudaMalloc((void**)&d_data.mean_x6,
    MAX_CHANNELS * sizeof(float));
  cudaMalloc((void**)&d_data.rms_x6,
    MAX_CHANNELS * sizeof(float));
  cudaMalloc((void**)&d_data.mean_x1,
    MAX_CHANNELS * sizeof(float));
  cudaMalloc((void**)&d_data.rms_x1,
    MAX_CHANNELS * sizeof(float));

//  cudaMalloc((void**)&d_data.gains, 
//    MAX_CHANNELS * sizeof(EcalMGPAGainRatio));

  cudaMalloc((void**)&d_data.gain12Over6,
    MAX_CHANNELS * sizeof(float));
  cudaMalloc((void**)&d_data.gain6Over1,
    MAX_CHANNELS * sizeof(float));

  cudaMalloc((void**)&d_data.xtals,
    MAX_CHANNELS * sizeof(EcalXtalGroupId));

  cudaMalloc((void**)&d_data.pulses,
    MAX_CHANNELS * sizeof(EcalPulseShape));
  cudaMalloc((void**)&d_data.epulses,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::FullSampleVector));

  cudaMalloc((void**)&d_data.covariances,
    MAX_CHANNELS * sizeof(EcalPulseCovariance));
  cudaMalloc((void**)&d_data.pulse_covariances,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::FullSampleMatrix));

  cudaMalloc((void**)&d_data.rechits,
    MAX_CHANNELS * sizeof(EcalUncalibratedRecHit));

  cudaMalloc((void**)&d_data.noisecorrs,
    3 * sizeof(ecal::multifit::v1::SampleMatrixD)); // size of std::array
  cudaMalloc((void**)&d_data.noisecov,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::SampleMatrix));
  cudaMalloc((void**)&d_data.pulse_matrix,
    MAX_CHANNELS * sizeof(ecal::multifit::v1::PulseMatrixType));
  cudaMalloc((void**)&d_data.bxs,
    sizeof(ecal::multifit::v1::BXVectorType));

  cudaMalloc((void**)&d_data.sample_mask, sizeof(EcalSampleMask));
  cudaMalloc((void**)&d_data.EBTimeCorrAmplitudeBins,
    sizeof(float) * MAX_TIME_BIAS_CORRECTIONS);
  cudaMalloc((void**)&d_data.EBTimeCorrShiftBins,
    sizeof(float) * MAX_TIME_BIAS_CORRECTIONS);
  cudaMalloc((void**)&d_data.EETimeCorrAmplitudeBins,
    sizeof(float) * MAX_TIME_BIAS_CORRECTIONS);
  cudaMalloc((void**)&d_data.EETimeCorrShiftBins,
    sizeof(float) * MAX_TIME_BIAS_CORRECTIONS);
  cudaMalloc((void**)&d_data.weights,
    sizeof(ecal::multifit::v1::EMatrix)*2*MAX_CHANNELS);
  cudaMalloc((void**)&d_data.statuses,
    sizeof(bool) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.chi2,
    sizeof(float)*MAX_CHANNELS);
  cudaMalloc((void**)&d_data.energies,
    sizeof(float)*MAX_CHANNELS);
  cudaMalloc((void**)&d_data.hasSwitchToGain6,
    sizeof(bool) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.hasSwitchToGain1,
    sizeof(bool) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.isSaturated,
    sizeof(bool) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.state_flags,
    sizeof(char) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.sample_values,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS * EcalDataFrame::MAXSAMPLES);
  cudaMalloc((void**)&d_data.sample_value_errors,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS * EcalDataFrame::MAXSAMPLES);
  cudaMalloc((void**)&d_data.useless_sample_values,
    sizeof(bool) * MAX_CHANNELS * EcalDataFrame::MAXSAMPLES);
  cudaMalloc((void**)&d_data.chi2sNullHypot,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.sum0sNullHypot,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.sumAAsNullHypot,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.pedestal_nums,
    sizeof(char) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.amplitudeFitParametersEB,
    sizeof(SampleVector::Scalar) * EBamplitudeFitParameters_.size());
  cudaMalloc((void**)&d_data.amplitudeFitParametersEE,
    sizeof(SampleVector::Scalar) * EEamplitudeFitParameters_.size());
  cudaMalloc((void**)&d_data.timeFitParametersEB,
    sizeof(SampleVector::Scalar) * EBtimeFitParameters_.size());
  cudaMalloc((void**)&d_data.timeFitParametersEE,
    sizeof(SampleVector::Scalar) * EEtimeFitParameters_.size());
  cudaMalloc((void**)&d_data.tMaxAlphaBetas,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.tMaxErrorAlphaBetas,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.accTimeMax,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.accTimeWgt,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.tcState,
    sizeof(ecal::multifit::v1::TimeComputationState) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.ampMaxAlphaBeta,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.ampMaxError,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.timeMax,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.timeError,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.amplitudeMax,
    sizeof(SampleVector::Scalar) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.jitter,
    sizeof(float) * MAX_CHANNELS);
  cudaMalloc((void**)&d_data.jitterError,
    sizeof(float) * MAX_CHANNELS);
  ecal::cuda::assert_if_error();

  //
  // for configuration parameters, transfer once
  //

  // TODO: this copy is done on purpose -> to avoid dealing with the current
  // parameter configuration. replacing double/float ...
  std::vector<SampleVector::Scalar> 
      ebAmplitudeFitParameters(EBamplitudeFitParameters_.size()), 
      eeAmplitudeFitParameters(EEamplitudeFitParameters_.size());
  ebAmplitudeFitParameters[0] = static_cast<SampleVector::Scalar>(
    EBamplitudeFitParameters_[0]);
  ebAmplitudeFitParameters[1] = static_cast<SampleVector::Scalar>(
    EBamplitudeFitParameters_[1]);
  eeAmplitudeFitParameters[0] = static_cast<SampleVector::Scalar>(
    EEamplitudeFitParameters_[0]);
  eeAmplitudeFitParameters[1] = static_cast<SampleVector::Scalar>(
    EEamplitudeFitParameters_[1]);
  cudaMemcpy(d_data.amplitudeFitParametersEB,
             ebAmplitudeFitParameters.data(),
             ebAmplitudeFitParameters.size() * sizeof(SampleVector::Scalar),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_data.amplitudeFitParametersEE,
             eeAmplitudeFitParameters.data(),
             eeAmplitudeFitParameters.size() * sizeof(SampleVector::Scalar),
             cudaMemcpyHostToDevice);

  d_data.timeFitParametersSizeEB = EBtimeFitParameters_.size();
  d_data.timeFitParametersSizeEE = EEtimeFitParameters_.size();
  d_data.timeFitLimitsFirstEB = EBtimeFitLimits_.first;
  d_data.timeFitLimitsSecondEB = EBtimeFitLimits_.second;
  d_data.timeFitLimitsFirstEE = EEtimeFitLimits_.first;
  d_data.timeFitLimitsSecondEE = EEtimeFitLimits_.second;
  std::vector<SampleVector::Scalar> 
      timeFitParametersEB(d_data.timeFitParametersSizeEB), 
      timeFitParametersEE(d_data.timeFitParametersSizeEE);

  // assume the same size for eb/ee
  for (unsigned int i=0; i<d_data.timeFitParametersSizeEB; i++) {
      timeFitParametersEB[i] = static_cast<SampleVector::Scalar>(
        EBtimeFitParameters_[i]);
      timeFitParametersEE[i] = static_cast<SampleVector::Scalar>(
        EEtimeFitParameters_[i]);
  }
  cudaMemcpy(d_data.timeFitParametersEB,
             timeFitParametersEB.data(),
             timeFitParametersEB.size() * sizeof(SampleVector::Scalar),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_data.timeFitParametersEE,
             timeFitParametersEE.data(),
             timeFitParametersEE.size() * sizeof(SampleVector::Scalar),
             cudaMemcpyHostToDevice);

  d_data.timeConstantTermEB = EBtimeConstantTerm_;
  d_data.timeConstantTermEE = EEtimeConstantTerm_;
}

EcalUncalibRecHitWorkerMultiFitGPUNew::~EcalUncalibRecHitWorkerMultiFitGPUNew() {
    //
    // TODO
    //
    if (d_data.digis_data) {
        cudaFree(d_data.digis_data);
        cudaFree(d_data.ids);
        cudaFree(d_data.amplitudes);
        cudaFree(d_data.samples);
        cudaFree(d_data.gainsNoise);
        cudaFree(d_data.gainsPedestal);

//        cudaFree(d_data.pedestals);
        cudaFree(d_data.mean_x12);
        cudaFree(d_data.rms_x12);
        cudaFree(d_data.mean_x6);
        cudaFree(d_data.rms_x6);
        cudaFree(d_data.mean_x1);
        cudaFree(d_data.rms_x1);

//        cudaFree(d_data.gains);
        cudaFree(d_data.gain12Over6);
        cudaFree(d_data.gain6Over1);

        cudaFree(d_data.xtals);
        cudaFree(d_data.pulses);
        cudaFree(d_data.epulses);

        cudaFree(d_data.covariances);
        cudaFree(d_data.pulse_covariances);

        cudaFree(d_data.rechits);
        cudaFree(d_data.noisecorrs);
        cudaFree(d_data.sample_mask);

        cudaFree(d_data.noisecov);
        cudaFree(d_data.pulse_matrix);
        cudaFree(d_data.bxs);

        cudaFree(d_data.EBTimeCorrAmplitudeBins);
        cudaFree(d_data.EBTimeCorrShiftBins);
        cudaFree(d_data.EETimeCorrAmplitudeBins);
        cudaFree(d_data.EETimeCorrShiftBins);
        cudaFree(d_data.weights);
        cudaFree(d_data.statuses);
        cudaFree(d_data.chi2);
        cudaFree(d_data.energies);
        cudaFree(d_data.hasSwitchToGain6);
        cudaFree(d_data.hasSwitchToGain1);
        cudaFree(d_data.isSaturated);
        cudaFree(d_data.state_flags);
        cudaFree(d_data.sample_values);
        cudaFree(d_data.sample_value_errors);
        cudaFree(d_data.useless_sample_values);
        cudaFree(d_data.chi2sNullHypot);
        cudaFree(d_data.sum0sNullHypot);
        cudaFree(d_data.sumAAsNullHypot);
        cudaFree(d_data.pedestal_nums);
        cudaFree(d_data.amplitudeFitParametersEB);
        cudaFree(d_data.amplitudeFitParametersEE);
        cudaFree(d_data.tMaxAlphaBetas);
        cudaFree(d_data.tMaxErrorAlphaBetas);
        cudaFree(d_data.accTimeMax);
        cudaFree(d_data.accTimeWgt);
        cudaFree(d_data.tcState);
        cudaFree(d_data.ampMaxAlphaBeta);
        cudaFree(d_data.ampMaxError);
        cudaFree(d_data.timeMax);
        cudaFree(d_data.timeError);
        cudaFree(d_data.timeFitParametersEB);
        cudaFree(d_data.timeFitParametersEE);
        cudaFree(d_data.amplitudeMax);
        cudaFree(d_data.jitter);
        cudaFree(d_data.jitterError);
        ecal::cuda::assert_if_error();
    }
}

void
EcalUncalibRecHitWorkerMultiFitGPUNew::set(const edm::EventSetup& es)
{

        // common setup
        es.get<EcalGainRatiosRcd>().get(gains);
        es.get<EcalPedestalsRcd>().get(peds);

        // for the multifit method
        if(!ampErrorCalculation_) multiFitMethod_.disableErrorCalculation();
        es.get<EcalSamplesCorrelationRcd>().get(noisecovariances);
        es.get<EcalPulseShapesRcd>().get(pulseshapes);
        es.get<EcalPulseCovariancesRcd>().get(pulsecovariances);

        // weights parameters for the time
        es.get<EcalWeightXtalGroupsRcd>().get(grps);
        es.get<EcalTBWeightsRcd>().get(wgts);

	// which of the samples need be used
	es.get<EcalSampleMaskRcd>().get(sampleMaskHand_);

        // for the ratio method
        es.get<EcalTimeCalibConstantsRcd>().get(itime);
        es.get<EcalTimeOffsetConstantRcd>().get(offtime);

        // for the time correction methods
        es.get<EcalTimeBiasCorrectionsRcd>().get(timeCorrBias_);

        int nnoise = SampleVector::RowsAtCompileTime;
        SampleMatrix &noisecorEBg12 = noisecors_[1][0];
        SampleMatrix &noisecorEBg6 = noisecors_[1][1];
        SampleMatrix &noisecorEBg1 = noisecors_[1][2];
        SampleMatrix &noisecorEEg12 = noisecors_[0][0];
        SampleMatrix &noisecorEEg6 = noisecors_[0][1];
        SampleMatrix &noisecorEEg1 = noisecors_[0][2];
                
        for (int i=0; i<nnoise; ++i) {
          for (int j=0; j<nnoise; ++j) {
            int vidx = std::abs(j-i);
            noisecorEBg12(i,j) = noisecovariances->EBG12SamplesCorrelation[vidx];
            noisecorEEg12(i,j) = noisecovariances->EEG12SamplesCorrelation[vidx];
            noisecorEBg6(i,j)  = noisecovariances->EBG6SamplesCorrelation[vidx];
            noisecorEEg6(i,j)  = noisecovariances->EEG6SamplesCorrelation[vidx];
            noisecorEBg1(i,j)  = noisecovariances->EBG1SamplesCorrelation[vidx];
            noisecorEEg1(i,j)  = noisecovariances->EEG1SamplesCorrelation[vidx];
          }
	}
}

void
EcalUncalibRecHitWorkerMultiFitGPUNew::set(const edm::Event& evt)
{

  unsigned int bunchspacing = 450;

  if (useLumiInfoRunHeader_) {

      edm::Handle<unsigned int> bunchSpacingH;
      evt.getByToken(bunchSpacing_,bunchSpacingH);
      bunchspacing = *bunchSpacingH;
  }
  else {
    bunchspacing = bunchSpacingManual_;
  }

  if (useLumiInfoRunHeader_ || bunchSpacingManual_ > 0){
    if (bunchspacing == 25) {
      activeBX.resize(10);
      activeBX << -5,-4,-3,-2,-1,0,1,2,3,4;
    }
    else {
      //50ns configuration otherwise (also for no pileup)
      activeBX.resize(5);
      activeBX << -4,-2,0,2,4;
    }
  }
 
}

/**
 * Amplitude-dependent time corrections; EE and EB have separate corrections:
 * EXtimeCorrAmplitudes (ADC) and EXtimeCorrShifts (ns) need to have the same number of elements
 * Bins must be ordered in amplitude. First-last bins take care of under-overflows.
 *
 * The algorithm is the same for EE and EB, only the correction vectors are different.
 *
 * @return Jitter (in clock cycles) which will be added to UncalibRechit.setJitter(), 0 if no correction is applied.
 */
double EcalUncalibRecHitWorkerMultiFitGPUNew::timeCorrection(
    float ampli,
	const std::vector<float>& amplitudeBins,
    const std::vector<float>& shiftBins) {

  // computed initially in ns. Than turned in the BX's, as
  // EcalUncalibratedRecHit need be.
  double theCorrection = 0;

  // sanity check for arrays
  if (amplitudeBins.empty()) {
    edm::LogError("EcalRecHitError")
        << "timeCorrAmplitudeBins is empty, forcing no time bias corrections.";

    return 0;
  }

  if (amplitudeBins.size() != shiftBins.size()) {
    edm::LogError("EcalRecHitError")
        << "Size of timeCorrAmplitudeBins different from "
           "timeCorrShiftBins. Forcing no time bias corrections. ";

    return 0;
  }

  // FIXME? what about a binary search?
  int myBin = -1;
  for (int bin = 0; bin < (int) amplitudeBins.size(); bin++) {
    if (ampli > amplitudeBins[bin]) {
      myBin = bin;
    } else {
      break;
	}
  }

  if (myBin == -1) {
    theCorrection = shiftBins[0];
  } else if (myBin == ((int)(amplitudeBins.size() - 1))) {
    theCorrection = shiftBins[myBin];
  } else {
    // interpolate linearly between two assingned points
    theCorrection = (shiftBins[myBin + 1] - shiftBins[myBin]);
    theCorrection *= (((double) ampli) - amplitudeBins[myBin]) /
                     (amplitudeBins[myBin + 1] - amplitudeBins[myBin]);
    theCorrection += shiftBins[myBin];
  }

  // convert ns into clocks
  constexpr double inv25 = 1./25.;
  return theCorrection * inv25;
}

void
EcalUncalibRecHitWorkerMultiFitGPUNew::run( const edm::Event & evt,
                const EcalDigiCollection & digis,
                ecal::SoAUncalibratedRecHitCollection& result )
{
    if (digis.empty())
      return;

    // assume all digis come from the same subdetector (either barrel or endcap)
    DetId detid(digis.begin()->id());
    bool barrel = (detid.subdetId()==EcalBarrel);

    //
    // gather conditions to send to device
    //
//    std::vector<EcalPedestal> vpedestals;
//    std::vector<EcalMGPAGainRatio> vgains;
    ecal::multifit::v1::pedestal_data ped_data;
    ecal::multifit::v1::mgpagain_ratio_data gainratio_data;

    std::vector<EcalXtalGroupId> vxtals;
    std::vector<EcalPulseShape> vpulseshapes;
    std::vector<EcalPulseCovariance> vcovariances;
    std::vector<ecal::multifit::v1::EMatrix> vweights;
    const SampleMatrixGainArray &noisecors = noisecor(barrel);

    // 
    // TODO: employ hashed index on the device directly!
    // need  to resort conditions in the order of digis
    //
//    vpedestals.reserve(digis.size());
//    vgains.reserve(digis.size());
    ped_data.mean_x12.reserve(digis.size());
    ped_data.rms_x12.reserve(digis.size());
    ped_data.mean_x6.reserve(digis.size());
    ped_data.rms_x6.reserve(digis.size());
    ped_data.mean_x1.reserve(digis.size());
    ped_data.rms_x1.reserve(digis.size());

    gainratio_data.gain12Over6.reserve(digis.size());
    gainratio_data.gain6Over1.reserve(digis.size());

    ecal::multifit::v1::BXVectorType bxs;
    bxs << -5, -4, -3, -2, -1, 0, 1, 2, 3, 4;

    vxtals.reserve(digis.size());
    vpulseshapes.reserve(digis.size());
    vcovariances.reserve(digis.size());
    vweights.reserve(2*digis.size());
    for (auto const& digi : digis) {
        DetId detid(digi.id());
        const EcalPedestals::Item * aped = nullptr;
        const EcalMGPAGainRatio * aGain = nullptr;
        const EcalXtalGroupId * gid = nullptr;
        const EcalPulseShapes::Item * aPulse = nullptr;
        const EcalPulseCovariances::Item * aPulseCov = nullptr;
        if (barrel) {
            unsigned int hashedIndex = EBDetId(detid).hashedIndex();
            aped       = &peds->barrel(hashedIndex);
            aGain      = &gains->barrel(hashedIndex);
            gid        = &grps->barrel(hashedIndex);
            aPulse     = &pulseshapes->barrel(hashedIndex);
            aPulseCov  = &pulsecovariances->barrel(hashedIndex);
        } else {
            unsigned int hashedIndex = EEDetId(detid).hashedIndex();
            aped       = &peds->endcap(hashedIndex);
            aGain      = &gains->endcap(hashedIndex);
            gid        = &grps->endcap(hashedIndex);
            aPulse     = &pulseshapes->endcap(hashedIndex);
            aPulseCov  = &pulsecovariances->endcap(hashedIndex);
        }

        EcalTBWeights::EcalTDCId tdcid{1};
        auto const& weightMap = wgts->getMap();
        EcalTBWeights::EcalTBWeightMap::const_iterator wit;
        wit = weightMap.find(std::make_pair(*gid, tdcid));
        if (wit == weightMap.end()) {
            edm::LogError("EcalUncalibRecHitError") 
                << "No weights found for EcalGroupId: "
                << gid->id() << " and  Eca    lTDCId: " << tdcid
                << "\n  skipping digi with     id: " << detid.rawId();
            // TODO: digis array will need to be properly updated if
            // this digi does not need to be sent to the device
            assert(0);
        }
        EcalWeightSet const& wset = wit->second;
        auto const& mat1 = wset.getWeightsBeforeGainSwitch();
        auto const& mat2 = wset.getWeightsAfterGainSwitch();
        ecal::multifit::v1::EMatrix m1,m2;
        for (unsigned int irow=0; 
             irow<EcalWeightSet::EcalWeightMatrix::rep_type::kRows;
             irow++)
            for (unsigned int icol=0;
                 icol<EcalWeightSet::EcalWeightMatrix::rep_type::kCols;
                 icol++) {
                m1(irow, icol) = mat1(irow, icol);
                m2(irow, icol) = mat2(irow, icol);
            }
        vweights.push_back(m1);
        vweights.push_back(m2);

//        vpedestals.push_back(*aped);
//        vgains.push_back(*aGain);
        ped_data.mean_x12.push_back(aped->mean_x12);
        ped_data.rms_x12.push_back(aped->rms_x12);
        ped_data.mean_x6.push_back(aped->mean_x6);
        ped_data.rms_x6.push_back(aped->rms_x6);
        ped_data.mean_x1.push_back(aped->mean_x1);
        ped_data.rms_x1.push_back(aped->rms_x1);

        gainratio_data.gain12Over6.push_back(aGain->gain12Over6());
        gainratio_data.gain6Over1.push_back(aGain->gain6Over1());

        vxtals.push_back(*gid);
        vpulseshapes.push_back(*aPulse);
        vcovariances.push_back(*aPulseCov);
    }
    EcalSampleMask const& sample_mask = *sampleMaskHand_.product();
    
    // 
    // prepare the result
    //
    result.amplitude.resize(digis.size());
    result.chi2.resize(digis.size());
    result.did.resize(digis.size());
    result.jitter.resize(digis.size());
    result.jitterError.resize(digis.size());

    // 
    // launch
    //
    ecal::multifit::v1::host_data h_data{&digis,
                                     ped_data, gainratio_data, sample_mask,
                                     &vxtals, &vpulseshapes, &vcovariances,
                                     &noisecors, &(*timeCorrBias_),
                                     &vweights, &bxs, result};
    ecal::multifit::v1::scatter(h_data, d_data, conf);
}

void
EcalUncalibRecHitWorkerMultiFitGPUNew::run( const edm::Event & evt,
                const EcalDigiCollection & digis,
                EcalUncalibratedRecHitCollection & result )
{
    /*
    if (digis.empty())
      return;

    // assume all digis come from the same subdetector (either barrel or endcap)
    DetId detid(digis.begin()->id());
    bool barrel = (detid.subdetId()==EcalBarrel);

    //
    // gather conditions to send to device
    //
//    std::vector<EcalPedestal> vpedestals;
//    std::vector<EcalMGPAGainRatio> vgains;
    ecal::multifit::v1::pedestal_data ped_data;
    ecal::multifit::v1::mgpagain_ratio_data gainratio_data;

    std::vector<EcalXtalGroupId> vxtals;
    std::vector<EcalPulseShape> vpulseshapes;
    std::vector<EcalPulseCovariance> vcovariances;
    std::vector<ecal::multifit::v1::EMatrix> vweights;
    const SampleMatrixGainArray &noisecors = noisecor(barrel);

    // 
    // TODO: employ hashed index on the device directly!
    // need  to resort conditions in the order of digis
    //
//    vpedestals.reserve(digis.size());
//    vgains.reserve(digis.size());
    ped_data.mean_x12.reserve(digis.size());
    ped_data.rms_x12.reserve(digis.size());
    ped_data.mean_x6.reserve(digis.size());
    ped_data.rms_x6.reserve(digis.size());
    ped_data.mean_x1.reserve(digis.size());
    ped_data.rms_x1.reserve(digis.size());

    gainratio_data.gain12Over6.reserve(digis.size());
    gainratio_data.gain6Over1.reserve(digis.size());

    ecal::multifit::v1::BXVectorType bxs;
    bxs << -5, -4, -3, -2, -1, 0, 1, 2, 3, 4;

    vxtals.reserve(digis.size());
    vpulseshapes.reserve(digis.size());
    vcovariances.reserve(digis.size());
    vweights.reserve(2*digis.size());
    for (auto const& digi : digis) {
        DetId detid(digi.id());
        const EcalPedestals::Item * aped = nullptr;
        const EcalMGPAGainRatio * aGain = nullptr;
        const EcalXtalGroupId * gid = nullptr;
        const EcalPulseShapes::Item * aPulse = nullptr;
        const EcalPulseCovariances::Item * aPulseCov = nullptr;
        if (barrel) {
            unsigned int hashedIndex = EBDetId(detid).hashedIndex();
            aped       = &peds->barrel(hashedIndex);
            aGain      = &gains->barrel(hashedIndex);
            gid        = &grps->barrel(hashedIndex);
            aPulse     = &pulseshapes->barrel(hashedIndex);
            aPulseCov  = &pulsecovariances->barrel(hashedIndex);
        } else {
            unsigned int hashedIndex = EEDetId(detid).hashedIndex();
            aped       = &peds->endcap(hashedIndex);
            aGain      = &gains->endcap(hashedIndex);
            gid        = &grps->endcap(hashedIndex);
            aPulse     = &pulseshapes->endcap(hashedIndex);
            aPulseCov  = &pulsecovariances->endcap(hashedIndex);
        }

        EcalTBWeights::EcalTDCId tdcid{1};
        auto const& weightMap = wgts->getMap();
        EcalTBWeights::EcalTBWeightMap::const_iterator wit;
        wit = weightMap.find(std::make_pair(*gid, tdcid));
        if (wit == weightMap.end()) {
            edm::LogError("EcalUncalibRecHitError") 
                << "No weights found for EcalGroupId: "
                << gid->id() << " and  Eca    lTDCId: " << tdcid
                << "\n  skipping digi with     id: " << detid.rawId();
            // TODO: digis array will need to be properly updated if
            // this digi does not need to be sent to the device
            assert(0);
        }
        EcalWeightSet const& wset = wit->second;
        auto const& mat1 = wset.getWeightsBeforeGainSwitch();
        auto const& mat2 = wset.getWeightsAfterGainSwitch();
        ecal::multifit::v1::EMatrix m1,m2;
        for (unsigned int irow=0; 
             irow<EcalWeightSet::EcalWeightMatrix::rep_type::kRows;
             irow++)
            for (unsigned int icol=0;
                 icol<EcalWeightSet::EcalWeightMatrix::rep_type::kCols;
                 icol++) {
                m1(irow, icol) = mat1(irow, icol);
                m2(irow, icol) = mat2(irow, icol);
            }
        vweights.push_back(m1);
        vweights.push_back(m2);

//        vpedestals.push_back(*aped);
//        vgains.push_back(*aGain);
        ped_data.mean_x12.push_back(aped->mean_x12);
        ped_data.rms_x12.push_back(aped->rms_x12);
        ped_data.mean_x6.push_back(aped->mean_x6);
        ped_data.rms_x6.push_back(aped->rms_x6);
        ped_data.mean_x1.push_back(aped->mean_x1);
        ped_data.rms_x1.push_back(aped->rms_x1);

        gainratio_data.gain12Over6.push_back(aGain->gain12Over6());
        gainratio_data.gain6Over1.push_back(aGain->gain6Over1());

        vxtals.push_back(*gid);
        vpulseshapes.push_back(*aPulse);
        vcovariances.push_back(*aPulseCov);
    }
    EcalSampleMask const& sample_mask = *sampleMaskHand_.product();
    
    // 
    // prepare the result
    //
    result.resize(digis.size());
    EcalUncalibratedRecHitCollection rechits(digis.size());

    // 
    // launch
    //
    ecal::multifit::v1::host_data h_data{&digis, &result,
                                     ped_data, gainratio_data,
                                     &vxtals, &vpulseshapes, &vcovariances,
                                     &noisecors, &sample_mask, &(*timeCorrBias_),
                                     &vweights, &bxs};
    ecal::multifit::v1::scatter(h_data, d_data, conf);
    */
}

edm::ParameterSetDescription 
EcalUncalibRecHitWorkerMultiFitGPUNew::getAlgoDescription() {
  
  edm::ParameterSetDescription psd0;
  psd0.addNode((edm::ParameterDescription<std::vector<double>>("EBPulseShapeTemplate", {1.13979e-02, 7.58151e-01, 1.00000e+00, 8.87744e-01, 6.73548e-01, 4.74332e-01, 3.19561e-01, 2.15144e-01, 1.47464e-01, 1.01087e-01, 6.93181e-02, 4.75044e-02}, true) and
		edm::ParameterDescription<std::vector<double>>("EEPulseShapeTemplate", {1.16442e-01, 7.56246e-01, 1.00000e+00, 8.97182e-01, 6.86831e-01, 4.91506e-01, 3.44111e-01, 2.45731e-01, 1.74115e-01, 1.23361e-01, 8.74288e-02, 6.19570e-02}, true)));
  
  psd0.addNode((edm::ParameterDescription<std::string>("EEdigiCollection", "", true) and
		edm::ParameterDescription<std::string>("EBdigiCollection", "", true) and
		edm::ParameterDescription<std::string>("ESdigiCollection", "", true) and
		edm::ParameterDescription<bool>("UseLCcorrection", true, false) and
		edm::ParameterDescription<std::vector<double>>("EBCorrNoiseMatrixG12", {1.00000, 0.71073, 0.55721, 0.46089, 0.40449, 0.35931, 0.33924, 0.32439, 0.31581, 0.30481 }, true) and
		edm::ParameterDescription<std::vector<double>>("EECorrNoiseMatrixG12", {1.00000, 0.71373, 0.44825, 0.30152, 0.21609, 0.14786, 0.11772, 0.10165, 0.09465, 0.08098 }, true) and
		edm::ParameterDescription<std::vector<double>>("EBCorrNoiseMatrixG06", {1.00000, 0.70946, 0.58021, 0.49846, 0.45006, 0.41366, 0.39699, 0.38478, 0.37847, 0.37055 }, true) and
		edm::ParameterDescription<std::vector<double>>("EECorrNoiseMatrixG06", {1.00000, 0.71217, 0.47464, 0.34056, 0.26282, 0.20287, 0.17734, 0.16256, 0.15618, 0.14443 }, true) and
		edm::ParameterDescription<std::vector<double>>("EBCorrNoiseMatrixG01", {1.00000, 0.73354, 0.64442, 0.58851, 0.55425, 0.53082, 0.51916, 0.51097, 0.50732, 0.50409 }, true) and
		edm::ParameterDescription<std::vector<double>>("EECorrNoiseMatrixG01", {1.00000, 0.72698, 0.62048, 0.55691, 0.51848, 0.49147, 0.47813, 0.47007, 0.46621, 0.46265 }, true) and
		edm::ParameterDescription<bool>("EcalPreMixStage1", false, true) and
		edm::ParameterDescription<bool>("EcalPreMixStage2", false, true)));

  psd0.addOptionalNode((edm::ParameterDescription<std::vector<double>>("EBPulseShapeCovariance", {3.001e-06,  1.233e-05,  0.000e+00, -4.416e-06, -4.571e-06, -3.614e-06, -2.636e-06, -1.286e-06, -8.410e-07, -5.296e-07,  0.000e+00,  0.000e+00, 
 	   1.233e-05,  6.154e-05,  0.000e+00, -2.200e-05, -2.309e-05, -1.838e-05, -1.373e-05, -7.334e-06, -5.088e-06, -3.745e-06, -2.428e-06,  0.000e+00, 
 	   0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 
 	   -4.416e-06, -2.200e-05,  0.000e+00,  8.319e-06,  8.545e-06,  6.792e-06,  5.059e-06,  2.678e-06,  1.816e-06,  1.223e-06,  8.245e-07,  5.589e-07, 
 	   -4.571e-06, -2.309e-05,  0.000e+00,  8.545e-06,  9.182e-06,  7.219e-06,  5.388e-06,  2.853e-06,  1.944e-06,  1.324e-06,  9.083e-07,  6.335e-07, 
 	   -3.614e-06, -1.838e-05,  0.000e+00,  6.792e-06,  7.219e-06,  6.016e-06,  4.437e-06,  2.385e-06,  1.636e-06,  1.118e-06,  7.754e-07,  5.556e-07, 
 	   -2.636e-06, -1.373e-05,  0.000e+00,  5.059e-06,  5.388e-06,  4.437e-06,  3.602e-06,  1.917e-06,  1.322e-06,  9.079e-07,  6.529e-07,  4.752e-07, 
 	   -1.286e-06, -7.334e-06,  0.000e+00,  2.678e-06,  2.853e-06,  2.385e-06,  1.917e-06,  1.375e-06,  9.100e-07,  6.455e-07,  4.693e-07,  3.657e-07, 
 	   -8.410e-07, -5.088e-06,  0.000e+00,  1.816e-06,  1.944e-06,  1.636e-06,  1.322e-06,  9.100e-07,  9.115e-07,  6.062e-07,  4.436e-07,  3.422e-07, 
 	   -5.296e-07, -3.745e-06,  0.000e+00,  1.223e-06,  1.324e-06,  1.118e-06,  9.079e-07,  6.455e-07,  6.062e-07,  7.217e-07,  4.862e-07,  3.768e-07, 
 	   0.000e+00, -2.428e-06,  0.000e+00,  8.245e-07,  9.083e-07,  7.754e-07,  6.529e-07,  4.693e-07,  4.436e-07,  4.862e-07,  6.509e-07,  4.418e-07, 
 	   0.000e+00,  0.000e+00,  0.000e+00,  5.589e-07,  6.335e-07,  5.556e-07,  4.752e-07,  3.657e-07,  3.422e-07,  3.768e-07,  4.418e-07,  6.142e-07}, true) and
     edm::ParameterDescription<std::vector<double>>("EEPulseShapeCovariance", {3.941e-05,  3.333e-05,  0.000e+00, -1.449e-05, -1.661e-05, -1.424e-05, -1.183e-05, -6.842e-06, -4.915e-06, -3.411e-06,  0.000e+00,  0.000e+00, 
 	   3.333e-05,  2.862e-05,  0.000e+00, -1.244e-05, -1.431e-05, -1.233e-05, -1.032e-05, -5.883e-06, -4.154e-06, -2.902e-06, -2.128e-06,  0.000e+00, 
 	   0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 
 	   -1.449e-05, -1.244e-05,  0.000e+00,  5.840e-06,  6.649e-06,  5.720e-06,  4.812e-06,  2.708e-06,  1.869e-06,  1.330e-06,  9.186e-07,  6.446e-07, 
 	   -1.661e-05, -1.431e-05,  0.000e+00,  6.649e-06,  7.966e-06,  6.898e-06,  5.794e-06,  3.157e-06,  2.184e-06,  1.567e-06,  1.084e-06,  7.575e-07, 
 	   -1.424e-05, -1.233e-05,  0.000e+00,  5.720e-06,  6.898e-06,  6.341e-06,  5.347e-06,  2.859e-06,  1.991e-06,  1.431e-06,  9.839e-07,  6.886e-07, 
 	   -1.183e-05, -1.032e-05,  0.000e+00,  4.812e-06,  5.794e-06,  5.347e-06,  4.854e-06,  2.628e-06,  1.809e-06,  1.289e-06,  9.020e-07,  6.146e-07, 
 	   -6.842e-06, -5.883e-06,  0.000e+00,  2.708e-06,  3.157e-06,  2.859e-06,  2.628e-06,  1.863e-06,  1.296e-06,  8.882e-07,  6.108e-07,  4.283e-07, 
 	   -4.915e-06, -4.154e-06,  0.000e+00,  1.869e-06,  2.184e-06,  1.991e-06,  1.809e-06,  1.296e-06,  1.217e-06,  8.669e-07,  5.751e-07,  3.882e-07, 
 	   -3.411e-06, -2.902e-06,  0.000e+00,  1.330e-06,  1.567e-06,  1.431e-06,  1.289e-06,  8.882e-07,  8.669e-07,  9.522e-07,  6.717e-07,  4.293e-07, 
 	   0.000e+00, -2.128e-06,  0.000e+00,  9.186e-07,  1.084e-06,  9.839e-07,  9.020e-07,  6.108e-07,  5.751e-07,  6.717e-07,  7.911e-07,  5.493e-07, 
 	   0.000e+00,  0.000e+00,  0.000e+00,  6.446e-07,  7.575e-07,  6.886e-07,  6.146e-07,  4.283e-07,  3.882e-07,  4.293e-07,  5.493e-07,  7.027e-07}, true)), true);

 edm::ParameterSetDescription psd;
 psd.addNode(edm::ParameterDescription<std::vector<int>>("activeBXs", {-5,-4,-3,-2,-1,0,1,2,3,4}, true) and
	      edm::ParameterDescription<bool>("ampErrorCalculation", true, true) and
	      edm::ParameterDescription<bool>("useLumiInfoRunHeader", true, true) and
	      edm::ParameterDescription<int>("bunchSpacing", 0, true) and
	      edm::ParameterDescription<bool>("doPrefitEB", false, true) and
	      edm::ParameterDescription<bool>("doPrefitEE", false, true) and
	      edm::ParameterDescription<double>("prefitMaxChiSqEB", 25., true) and
	      edm::ParameterDescription<double>("prefitMaxChiSqEE", 10., true) and
	      edm::ParameterDescription<bool>("dynamicPedestalsEB", false, true) and
	      edm::ParameterDescription<bool>("dynamicPedestalsEE", false, true) and
	      edm::ParameterDescription<bool>("mitigateBadSamplesEB", false, true) and
	      edm::ParameterDescription<bool>("mitigateBadSamplesEE", false, true) and
	      edm::ParameterDescription<bool>("gainSwitchUseMaxSampleEB", false, true) and
	      edm::ParameterDescription<bool>("gainSwitchUseMaxSampleEE", false, true) and
	      edm::ParameterDescription<bool>("selectiveBadSampleCriteriaEB", false, true) and
	      edm::ParameterDescription<bool>("selectiveBadSampleCriteriaEE", false, true) and
	      edm::ParameterDescription<double>("addPedestalUncertaintyEB", 0., true) and
	      edm::ParameterDescription<double>("addPedestalUncertaintyEE", 0., true) and
	      edm::ParameterDescription<bool>("simplifiedNoiseModelForGainSwitch", true, true) and
	      edm::ParameterDescription<std::string>("timealgo", "RatioMethod", true) and
	      edm::ParameterDescription<std::vector<double>>("EBtimeFitParameters", {-2.015452e+00, 3.130702e+00, -1.234730e+01, 4.188921e+01, -8.283944e+01, 9.101147e+01, -5.035761e+01, 1.105621e+01}, true) and
	      edm::ParameterDescription<std::vector<double>>("EEtimeFitParameters", {-2.390548e+00, 3.553628e+00, -1.762341e+01, 6.767538e+01, -1.332130e+02, 1.407432e+02, -7.541106e+01, 1.620277e+01}, true) and
	      edm::ParameterDescription<std::vector<double>>("EBamplitudeFitParameters", {1.138,1.652}, true) and
	      edm::ParameterDescription<std::vector<double>>("EEamplitudeFitParameters", {1.890,1.400}, true) and
	      edm::ParameterDescription<double>("EBtimeFitLimits_Lower", 0.2, true) and
	      edm::ParameterDescription<double>("EBtimeFitLimits_Upper", 1.4, true) and
	      edm::ParameterDescription<double>("EEtimeFitLimits_Lower", 0.2, true) and
	      edm::ParameterDescription<double>("EEtimeFitLimits_Upper", 1.4, true) and
	      edm::ParameterDescription<double>("EBtimeConstantTerm", .6, true) and
	      edm::ParameterDescription<double>("EEtimeConstantTerm", 1.0, true) and
	      edm::ParameterDescription<double>("EBtimeNconst", 28.5, true) and
	      edm::ParameterDescription<double>("EEtimeNconst", 31.8, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain12pEB", 5, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain12mEB", 5, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain61pEB", 5, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain61mEB", 5, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain12pEE", 1000, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain12mEE", 1000, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain61pEE", 1000, true) and
	      edm::ParameterDescription<double>("outOfTimeThresholdGain61mEE", 1000, true) and
	      edm::ParameterDescription<double>("amplitudeThresholdEB", 10, true) and
	      edm::ParameterDescription<double>("amplitudeThresholdEE", 10, true) and
	      edm::ParameterDescription<double>("ebSpikeThreshold", 1.042, true) and
	      edm::ParameterDescription<std::vector<double>>("ebPulseShape", {5.2e-05,-5.26e-05 , 6.66e-05, 0.1168, 0.7575, 1.,  0.8876, 0.6732, 0.4741,  0.3194}, true) and
	      edm::ParameterDescription<std::vector<double>>("eePulseShape", {5.2e-05,-5.26e-05 , 6.66e-05, 0.1168, 0.7575, 1.,  0.8876, 0.6732, 0.4741,  0.3194}, true) and
	      edm::ParameterDescription<bool>("kPoorRecoFlagEB", true, true) and
	      edm::ParameterDescription<bool>("kPoorRecoFlagEE", false, true) and
	      edm::ParameterDescription<double>("chi2ThreshEB_", 65.0, true) and
	      edm::ParameterDescription<double>("chi2ThreshEE_", 50.0, true) and
          edm::ParameterDescription<std::vector<int>>("threads", {256, 1, 1}, true)
          and 
	      edm::ParameterDescription<edm::ParameterSetDescription>("EcalPulseShapeParameters", psd0, true));

 return psd;
}



#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalUncalibRecHitWorkerFactory, EcalUncalibRecHitWorkerMultiFitGPUNew, "EcalUncalibRecHitWorkerMultiFitGPUNew" );
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitFillDescriptionWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalUncalibRecHitFillDescriptionWorkerFactory, EcalUncalibRecHitWorkerMultiFitGPUNew, "EcalUncalibRecHitWorkerMultiFitGPUNew" );
