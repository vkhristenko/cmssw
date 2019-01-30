#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo_gpu_cputest.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"

#include <iostream>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS_gpu_cputest.h"

namespace ecal { namespace multifit { namespace test {

class EcalUncalibRecHitMultiFitAlgo
{
      
public:
      
    EcalUncalibRecHitMultiFitAlgo();
      
    EcalUncalibratedRecHit makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, SampleMatrix const* noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX);
     
    void disableErrorCalculation() { _computeErrors = false; }
      
    void setDoPrefit(bool b) { _doPrefit = b; }
      void setPrefitMaxChiSq(double x) { _prefitMaxChiSq = x; }
      void setDynamicPedestals(bool b) { _dynamicPedestals = b; }
      void setMitigateBadSamples(bool b) { _mitigateBadSamples = b; }
      void setSelectiveBadSampleCriteria(bool b) { _selectiveBadSampleCriteria = b; }
      void setAddPedestalUncertainty(double x) { _addPedestalUncertainty = x; }
      void setSimplifiedNoiseModelForGainSwitch(bool b) { _simplifiedNoiseModelForGainSwitch = b; }
      void setGainSwitchUseMaxSample(bool b) { _gainSwitchUseMaxSample = b; }
                                 
private:
    PulseChiSqSNNLS _pulsefunc;
    PulseChiSqSNNLS _pulsefuncSingle;
    bool _computeErrors;
    bool _doPrefit;
    double _prefitMaxChiSq;
    bool _dynamicPedestals;
    bool _mitigateBadSamples;
    bool _selectiveBadSampleCriteria;
    double _addPedestalUncertainty;
    bool _simplifiedNoiseModelForGainSwitch;
    bool _gainSwitchUseMaxSample;
    BXVector _singlebx;
};

  
EcalUncalibRecHitMultiFitAlgo::EcalUncalibRecHitMultiFitAlgo() : 
  _computeErrors(true),
  _doPrefit(false),
  _prefitMaxChiSq(1.0),
  _dynamicPedestals(false),
  _mitigateBadSamples(false),
  _selectiveBadSampleCriteria(false),
  _addPedestalUncertainty(0.),
  _simplifiedNoiseModelForGainSwitch(true),
  _gainSwitchUseMaxSample(false) {
    
  _singlebx.resize(1);
  _singlebx << 0;
  
  /*
  _pulsefuncSingle.disableErrorCalculation();
  _pulsefuncSingle.setMaxIters(1);
  _pulsefuncSingle.setMaxIterWarnings(false);
    */
}

/// compute rechits
  
EcalUncalibratedRecHit 
EcalUncalibRecHitMultiFitAlgo::makeRecHit(const EcalDataFrame& dataFrame, 
                                          const EcalPedestals::Item * aped, 
                                          const EcalMGPAGainRatio * aGain, 
                                          SampleMatrix const* noisecors, 
                                          const FullSampleVector &fullpulse, 
                                          const FullSampleMatrix &fullpulsecov, 
                                          const BXVector &activeBX) {
  uint32_t flags = 0;
  
  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;
  
  double maxamplitude = -std::numeric_limits<double>::max();
  const unsigned int iSampleMax = 5;
  const unsigned int iFullPulseMax = 9;
  
  double pedval = 0.;
    
  SampleVector amplitudes;
  SampleGainVector gainsNoise;
  SampleGainVector gainsPedestal;
  SampleGainVector badSamples = SampleGainVector::Zero();
  bool hasSaturation = dataFrame.isSaturated();
  bool hasGainSwitch = hasSaturation || dataFrame.hasSwitchToGain6() || dataFrame.hasSwitchToGain1();
  
  //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
  bool dynamicPedestal = _dynamicPedestals && !hasGainSwitch;
  
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
        
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    
    double amplitude = 0.;
    int gainId = sample.gainId();
    
    double pedestal = 0.;
    double gainratio = 1.;
    
    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
      gainsNoise[iSample] = 2;
      gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
      gainsNoise[iSample] = 0;
      gainsPedestal[iSample] = dynamicPedestal ? 0 : -1; //-1 for static pedestal
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
      gainsNoise[iSample] = 1;
      gainsPedestal[iSample] = dynamicPedestal ? 1 : -1; //-1 for static pedestals
    }

    if (dynamicPedestal) {
      amplitude = (double)(sample.adc())*gainratio;
    }
    else {
      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    }
    
    if (gainId == 0) {
        // TODO: can we do something better
        printf("EcalUncalibRecHitMultiFitAlgo::Saturation encountered.  Multifit is not intended to be used for saturated channels.\n");
//       edm::LogError("EcalUncalibRecHitMultiFitAlgo")<< "Saturation encountered.  Multifit is not intended to be used for saturated channels.";
      //saturation
      if (dynamicPedestal) {
        amplitude = 4095.*gainratio;
      }
      else {
        amplitude = (4095. - pedestal) * gainratio;
      }
    }
        
    amplitudes[iSample] = amplitude;
    
    if (iSample==iSampleMax) {
      maxamplitude = amplitude;
      pedval = pedestal;
    }
        
  }

  double amplitude, amperr, chisq;
  bool status = false;
    
  //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
  //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
  //option 1: use simple max-sample algorithm
  if (hasGainSwitch && _gainSwitchUseMaxSample) {
    double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
    EcalUncalibratedRecHit rh( dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags );
    rh.setAmplitudeError(0.);
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0) {
        rh.setOutOfTimeAmplitude(bx+5, 0.0);
      }
    }
    return rh;
  }

  //option2: A floating negative single-sample offset is added to the fit
  //such that the affected sample is treated only as a lower limit for the true amplitude
  bool mitigateBadSample = _mitigateBadSamples && hasGainSwitch && iSampleMax>0;
  mitigateBadSample &= (!_selectiveBadSampleCriteria || (gainsNoise.coeff(iSampleMax-1)!=gainsNoise.coeff(iSampleMax)) );
  if (mitigateBadSample) {
    badSamples[iSampleMax-1] = 1;
  }
  
  //compute noise covariance matrix, which depends on the sample gains
  SampleMatrix noisecov;
  if (hasGainSwitch) {
    double pedrmss[3] = {aped->rms_x12, aped->rms_x6, aped->rms_x1};
    double gainratios[3] = { 1., aGain->gain12Over6(), aGain->gain6Over1()*aGain->gain12Over6()};
    if (_simplifiedNoiseModelForGainSwitch) {
      int gainidxmax = gainsNoise[iSampleMax];
      noisecov = gainratios[gainidxmax]*gainratios[gainidxmax]*pedrmss[gainidxmax]*pedrmss[gainidxmax]*noisecors[gainidxmax];
      if (!dynamicPedestal && _addPedestalUncertainty>0.) {
        //add fully correlated component to noise covariance to inflate pedestal uncertainty
        noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
      }
    }
    else {
      noisecov = SampleMatrix::Zero();
      // noisecors is type var[2] is an array of 2
      // size needs to be propogated (on cpu it was std::array)
      for (unsigned int gainidx=0; gainidx<2; ++gainidx) {
        SampleGainVector mask = gainidx*SampleGainVector::Ones();
        SampleVector pedestal = (gainsNoise.array()==mask.array()).cast<SampleVector::value_type>();
        if (pedestal.maxCoeff()>0.) {
          //select out relevant components of each correlation matrix, and assume no correlation between samples with
          //different gain
          noisecov += gainratios[gainidx]*gainratios[gainidx]*pedrmss[gainidx]*pedrmss[gainidx]*pedestal.asDiagonal()*noisecors[gainidx]*pedestal.asDiagonal();
          if (!dynamicPedestal && _addPedestalUncertainty>0.) {
            //add fully correlated component to noise covariance to inflate pedestal uncertainty
            noisecov += gainratios[gainidx]*gainratios[gainidx]*_addPedestalUncertainty*_addPedestalUncertainty*pedestal.asDiagonal()*SampleMatrix::Ones()*pedestal.asDiagonal();
          }
        }
      }
    }
  }
  else {
    noisecov = aped->rms_x12*aped->rms_x12*noisecors[0];
    if (!dynamicPedestal && _addPedestalUncertainty>0.) {
      //add fully correlated component to noise covariance to inflate pedestal uncertainty
      noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
    }
  }
  
  //optimized one-pulse fit for hlt
  bool usePrefit = false;
  if (_doPrefit) {
    status = _pulsefuncSingle.DoFit(amplitudes,noisecov,_singlebx,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    amplitude = status ? _pulsefuncSingle.X()[0] : 0.;
    amperr = status ? _pulsefuncSingle.Errors()[0] : 0.;
    chisq = _pulsefuncSingle.ChiSq();
    
    if (chisq < _prefitMaxChiSq) {
      usePrefit = true;
    }
  }
  
  if (!usePrefit) {
  
    if(!_computeErrors) _pulsefunc.disableErrorCalculation();
    status = _pulsefunc.DoFit(amplitudes,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    chisq = _pulsefunc.ChiSq();
    
    if (!status) {
      // TODO: Can we do better
      printf("EcalUncalibRecHitMultiFitAlgo::makeRecHit failed fit\n");
//      edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
    }

    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      if (_pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }
    
    amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
    amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
  
  }
  
  double jitter = 0.;
  
  EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
  rh.setAmplitudeError(amperr);
  
  if (!usePrefit) {
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0 && ecal::abs(bx)<100) {
        rh.setOutOfTimeAmplitude(bx+5, status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
      else if (bx==(100+gainsPedestal[iSampleMax])) {
        rh.setPedestal(status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
    }
  }
  
  return rh;
}

void kernel_reconstruct(uint16_t const *digis,
                              uint32_t const *ids,
                              EcalPedestal const *pedestals,
                              EcalMGPAGainRatio const *gains,
                              EcalXtalGroupId const *xtals,
                              EcalPulseShape const *shapes,
                              EcalPulseCovariance const *covariances,
                              EcalUncalibratedRecHit *rechits,
                              SampleMatrix const *noisecors,
                              unsigned int size) {
//    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    for (unsigned int idx=0; idx < size; idx++) {
        uint16_t const* p_current_digi = &digis[idx*EcalDataFrame::MAXSAMPLES];
        DetId  current_id{ids[idx]};
        EcalDataFrame edf{edm::DataFrame{ids[idx], p_current_digi, EcalDataFrame::MAXSAMPLES}};
        auto const* aped = &pedestals[idx];
        auto const* aGain = &gains[idx];
        auto const* gid = &xtals[idx];
        auto const* aPulse = &shapes[idx];
        auto const* aPulseCov = &covariances[idx];

        FullSampleVector fullpulse(FullSampleVector::Zero());
        FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());

        double pedVec[3]     = {aped->mean_x12, aped->mean_x6, aped->mean_x1 };
        double pedRMSVec[3]  = {aped->rms_x12,  aped->rms_x6,  aped->rms_x1 };
        double gainRatios[3] = {1., aGain->gain12Over6(), 
                                aGain->gain6Over1()*aGain->gain12Over6()};

        for (int i=0; i<EcalPulseShape::TEMPLATESAMPLES; ++i)
            fullpulse(i+7) = aPulse->pdfval[i];
        for(int i=0; i<EcalPulseShape::TEMPLATESAMPLES;i++)
        for(int j=0; j<EcalPulseShape::TEMPLATESAMPLES;j++)
            fullpulsecov(i+7,j+7) = aPulseCov->covval[i][j];

        int lastSampleBeforeSaturation = -2;
        for(unsigned int iSample = 0; 
            iSample < EcalDataFrame::MAXSAMPLES; iSample++){
            if (edf.sample(iSample).gainId() == 0 ) {
                lastSampleBeforeSaturation=iSample-1;
                break;
            }
        }

        EcalUncalibratedRecHit rh{current_id, 4095*12, 0.f, 0.f, 0.f};

        // === amplitude computation ===
        if ( lastSampleBeforeSaturation == 4 ) { 
            // saturation on the expected max sample
            EcalUncalibratedRecHit tmp{current_id, 4095*12, 0.f, 0.f, 0.f};
            tmp.setFlagBit( EcalUncalibratedRecHit::kSaturated );
            // do not propagate the default chi2 = -1 value 
            // to the calib rechit (mapped to 64), set it to 0 when saturation
            tmp.setChi2(0);
            rechits[idx] = tmp;
        } else if ( lastSampleBeforeSaturation >= -1 ) { 
            // saturation on other samples: cannot extrapolate from the fourth one
            int gainId = edf.sample(5).gainId();
            if (gainId==0) gainId=3;
            auto pedestal = pedVec[gainId-1];
            auto gainratio = gainRatios[gainId-1];
            double amplitude = ((double)(edf.sample(5).adc()) - 
                    pedestal) * gainratio;
            EcalUncalibratedRecHit tmp{current_id, static_cast<float>(amplitude), 0.f, 0.f, 0.f};
            tmp.setFlagBit( EcalUncalibratedRecHit::kSaturated );
            // do not propagate the default chi2 = -1 value to the calib rechit (mapped to 64), set it to 0 when saturation
            tmp.setChi2(0);
            rechits[idx] = tmp;
        } else {
            // multifit
            //result.push_back(multiFitMethod_.makeRecHit(*itdg, aped, aGain, 
            //                 noisecors, fullpulse, fullpulsecov, activeBX));
            BXVector activeBX;
            activeBX.resize(10);
            activeBX << -5,-4,-3,-2,-1,0,1,2,3,4;
            EcalUncalibRecHitMultiFitAlgo algo{};
            rechits[idx] = algo.makeRecHit(edf, aped, aGain, 
                                           noisecors, fullpulse, fullpulsecov, 
                                           activeBX);
        }
        
    }
}

void scatter(EcalDigiCollection const& digis,
             EcalUncalibratedRecHitCollection& vrechits,
             std::vector<EcalPedestal> const& vpedestals,
             std::vector<EcalMGPAGainRatio> const& vgains,
             std::vector<EcalXtalGroupId> const& vxtals,
             std::vector<EcalPulseShape> const& vshapes,
             std::vector<EcalPulseCovariance> const& vcovariances,
             SampleMatrixGainArray const& noisecors) {
    auto const& ids = digis.ids();
    auto const& digis_data = digis.data();
    using digis_type = std::vector<uint16_t>;
    using dids_type = std::vector<uint32_t>;
    digis_type::value_type *d_digis_data;
    dids_type::value_type *d_ids;
    EcalPedestal *d_pedestals;
    EcalMGPAGainRatio *d_gains;
    EcalXtalGroupId *d_xtals;
    EcalPulseShape *d_shapes;
    EcalPulseCovariance *d_covariances;
    EcalUncalibratedRecHit *d_rechits;
    SampleMatrix *d_noisecors;
    
    // 
    // copy to the device
    // TODO: can conditions be copied only once when updated?
    //
    /*
    cudaMemcpy(d_digis_data, digis_data.data(),
        digis_data.size() * sizeof(digis_type::value_type),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_ids, ids.data(),
        ids.size() * sizeof(dids_type::value_type),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_pedestals, vpedestals.data(),
        vpedestals.size() * sizeof(EcalPedestal),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_gains, vgains.data(),
        vgains.size() * sizeof(EcalMGPAGainRatio),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_xtals, vxtals.data(),
        vxtals.size() * sizeof(EcalXtalGroupId),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_shapes, vpulses.data(),
        vpulses.size() * sizeof(EcalPulseShape),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_covariances, vcovariances.data(),
        vcovariances.size() * sizeof(EcalPulseCovariance),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_rechits, &(*rechits.begin()),
        rechits.size() * sizeof(EcalUncalibratedRecHit),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_noisecors, noisecors.data(),
        noisecors.size() * sizeof(SampleMatrix),
        cudaMemcpyHostToDevice);
    ecal::cuda::assert_if_error();
    */
    
    //
    // launch 
    // TODO: further ntreads/nblocks optimizations...
    //
    std::cout << "ecal::multifit::scatter()" << std::endl;
    kernel_reconstruct(
        digis_data.data(),
        ids.data(),
        /* d_rechits, */
        vpedestals.data(),
        vgains.data(),
        vxtals.data(),
        vshapes.data(),
        vcovariances.data(),
        &(*vrechits.begin()),
        noisecors.data(),
        digis.size()
    );
}

}}}

/*
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"

EcalUncalibRecHitMultiFitAlgo::EcalUncalibRecHitMultiFitAlgo() : 
  _computeErrors(true),
  _doPrefit(false),
  _prefitMaxChiSq(1.0),
  _dynamicPedestals(false),
  _mitigateBadSamples(false),
  _selectiveBadSampleCriteria(false),
  _addPedestalUncertainty(0.),
  _simplifiedNoiseModelForGainSwitch(true),
  _gainSwitchUseMaxSample(false){
    
  _singlebx.resize(1);
  _singlebx << 0;
  
  _pulsefuncSingle.disableErrorCalculation();
  _pulsefuncSingle.setMaxIters(1);
  _pulsefuncSingle.setMaxIterWarnings(false);
    
}

/// compute rechits
EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgo::makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX) {

  uint32_t flags = 0;
  
  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;
  
  double maxamplitude = -std::numeric_limits<double>::max();
  const unsigned int iSampleMax = 5;
  const unsigned int iFullPulseMax = 9;
  
  double pedval = 0.;
    
  SampleVector amplitudes;
  SampleGainVector gainsNoise;
  SampleGainVector gainsPedestal;
  SampleGainVector badSamples = SampleGainVector::Zero();
  bool hasSaturation = dataFrame.isSaturated();
  bool hasGainSwitch = hasSaturation || dataFrame.hasSwitchToGain6() || dataFrame.hasSwitchToGain1();
  
  //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
  bool dynamicPedestal = _dynamicPedestals && !hasGainSwitch;
  
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
        
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    
    double amplitude = 0.;
    int gainId = sample.gainId();
    
    double pedestal = 0.;
    double gainratio = 1.;
    
    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
      gainsNoise[iSample] = 2;
      gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
      gainsNoise[iSample] = 0;
      gainsPedestal[iSample] = dynamicPedestal ? 0 : -1; //-1 for static pedestal
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
      gainsNoise[iSample] = 1;
      gainsPedestal[iSample] = dynamicPedestal ? 1 : -1; //-1 for static pedestals
    }

    if (dynamicPedestal) {
      amplitude = (double)(sample.adc())*gainratio;
    }
    else {
      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    }
    
    if (gainId == 0) {
       edm::LogError("EcalUncalibRecHitMultiFitAlgo")<< "Saturation encountered.  Multifit is not intended to be used for saturated channels.";
      //saturation
      if (dynamicPedestal) {
        amplitude = 4095.*gainratio;
      }
      else {
        amplitude = (4095. - pedestal) * gainratio;
      }
    }
        
    amplitudes[iSample] = amplitude;
    
    if (iSample==iSampleMax) {
      maxamplitude = amplitude;
      pedval = pedestal;
    }
        
  }

  double amplitude, amperr, chisq;
  bool status = false;
    
  //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
  //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
  //option 1: use simple max-sample algorithm
  if (hasGainSwitch && _gainSwitchUseMaxSample) {
    double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
    EcalUncalibratedRecHit rh( dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags );
    rh.setAmplitudeError(0.);
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0) {
        rh.setOutOfTimeAmplitude(bx+5, 0.0);
      }
    }
    return rh;
  }

  //option2: A floating negative single-sample offset is added to the fit
  //such that the affected sample is treated only as a lower limit for the true amplitude
  bool mitigateBadSample = _mitigateBadSamples && hasGainSwitch && iSampleMax>0;
  mitigateBadSample &= (!_selectiveBadSampleCriteria || (gainsNoise.coeff(iSampleMax-1)!=gainsNoise.coeff(iSampleMax)) );
  if (mitigateBadSample) {
    badSamples[iSampleMax-1] = 1;
  }
  
  //compute noise covariance matrix, which depends on the sample gains
  SampleMatrix noisecov;
  if (hasGainSwitch) {
    std::array<double,3> pedrmss = {{aped->rms_x12, aped->rms_x6, aped->rms_x1}};
    std::array<double,3> gainratios = {{ 1., aGain->gain12Over6(), aGain->gain6Over1()*aGain->gain12Over6()}};
    if (_simplifiedNoiseModelForGainSwitch) {
      int gainidxmax = gainsNoise[iSampleMax];
      noisecov = gainratios[gainidxmax]*gainratios[gainidxmax]*pedrmss[gainidxmax]*pedrmss[gainidxmax]*noisecors[gainidxmax];
      if (!dynamicPedestal && _addPedestalUncertainty>0.) {
        //add fully correlated component to noise covariance to inflate pedestal uncertainty
        noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
      }
    }
    else {
      noisecov = SampleMatrix::Zero();
      for (unsigned int gainidx=0; gainidx<noisecors.size(); ++gainidx) {
        SampleGainVector mask = gainidx*SampleGainVector::Ones();
        SampleVector pedestal = (gainsNoise.array()==mask.array()).cast<SampleVector::value_type>();
        if (pedestal.maxCoeff()>0.) {
          //select out relevant components of each correlation matrix, and assume no correlation between samples with
          //different gain
          noisecov += gainratios[gainidx]*gainratios[gainidx]*pedrmss[gainidx]*pedrmss[gainidx]*pedestal.asDiagonal()*noisecors[gainidx]*pedestal.asDiagonal();
          if (!dynamicPedestal && _addPedestalUncertainty>0.) {
            //add fully correlated component to noise covariance to inflate pedestal uncertainty
            noisecov += gainratios[gainidx]*gainratios[gainidx]*_addPedestalUncertainty*_addPedestalUncertainty*pedestal.asDiagonal()*SampleMatrix::Ones()*pedestal.asDiagonal();
          }
        }
      }
    }
  }
  else {
    noisecov = aped->rms_x12*aped->rms_x12*noisecors[0];
    if (!dynamicPedestal && _addPedestalUncertainty>0.) {
      //add fully correlated component to noise covariance to inflate pedestal uncertainty
      noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
    }
  }
  
  //optimized one-pulse fit for hlt
  bool usePrefit = false;
  if (_doPrefit) {
    status = _pulsefuncSingle.DoFit(amplitudes,noisecov,_singlebx,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    amplitude = status ? _pulsefuncSingle.X()[0] : 0.;
    amperr = status ? _pulsefuncSingle.Errors()[0] : 0.;
    chisq = _pulsefuncSingle.ChiSq();
    
    if (chisq < _prefitMaxChiSq) {
      usePrefit = true;
    }
  }
  
  if (!usePrefit) {
  
    if(!_computeErrors) _pulsefunc.disableErrorCalculation();
    status = _pulsefunc.DoFit(amplitudes,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    chisq = _pulsefunc.ChiSq();
    
    if (!status) {
      edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
    }

    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      if (_pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }
    
    amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
    amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
  
  }
  
  double jitter = 0.;
  
  EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
  rh.setAmplitudeError(amperr);
  
  if (!usePrefit) {
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0 && std::abs(bx)<100) {
        rh.setOutOfTimeAmplitude(bx+5, status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
      else if (bx==(100+gainsPedestal[iSampleMax])) {
        rh.setPedestal(status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
    }
  }
  
  return rh;
}
*/
