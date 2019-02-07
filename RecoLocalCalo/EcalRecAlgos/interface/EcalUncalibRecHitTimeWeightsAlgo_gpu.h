#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitTimeWeightsAlgo_gpu_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitTimeWeightsAlgo_gpu_HH

/** \class EcalUncalibRecHitTimeWeightsAlgo
  *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
  *  using a weights method
  *
  *  \author J. Bendavid, E. Di Marco
  *  
  *  The chi2 computation with matrix is replaced by the chi2express which is  moved outside the weight algo
  *  (need to clean up the interface in next iteration so that we do not pass-by useless arrays)
  *
  */

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

namespace ecal { namespace multifit {

class EcalUncalibRecHitTimeWeightsAlgo 
{
 public:

  __device__
  EcalUncalibRecHitTimeWeightsAlgo() { };

  /// Compute time
  __device__
  double time(const EcalDataFrame& dataFrame, double amplitudes[10], 
              const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, 
              const FullSampleVector &fullpulse, 
              EMatrix const& weights_0, EMatrix const& weights_1) {
  
    constexpr unsigned int nsample = EcalDataFrame::MAXSAMPLES;
  
    double maxamplitude = -std::numeric_limits<double>::max();
  
    double pulsenorm = 0.;
    int iGainSwitch = 0;

 //   ROOT::Math::SVector<double,nsample> pedSubSamples;
    Eigen::Matrix<double, nsample, 1> pedSubSamples;
    for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    
      const EcalMGPASample &sample = dataFrame.sample(iSample);
    
      double amplitude = 0.;
      int gainId = sample.gainId();
    
      double pedestal = 0.;
      double gainratio = 1.;
        
      if (gainId==0 || gainId==3) {
        pedestal = aped->mean_x1;
        gainratio = aGain->gain6Over1()*aGain->gain12Over6();
        iGainSwitch = 1;
      }
      else if (gainId==1) {
        pedestal = aped->mean_x12;
        gainratio = 1.;
        iGainSwitch = 0;
      }
      else if (gainId==2) {
        pedestal = aped->mean_x6;
        gainratio = aGain->gain12Over6();
        iGainSwitch = 1;
      }

      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    
      if (gainId == 0) {
        //saturation
        amplitude = (4095. - pedestal) * gainratio;
      }
        
      pedSubSamples(iSample) = amplitude;

      if (amplitude>maxamplitude) {
        maxamplitude = amplitude;
      }    
      pulsenorm += fullpulse(iSample);
    }

//    std::vector<double>::const_iterator amplit;
//    for(amplit=amplitudes.begin(); amplit<amplitudes.end(); ++amplit) {
    for (unsigned int i=0; i<10; i++) {
//      int ipulse = std::distance(amplitudes.begin(),amplit);
      int ipulse = i;
      int bx = ipulse - 5;
      int firstsamplet = std::max(0,bx + 3);
      int offset = 7-3-bx;

//      TVectorD pulse;
      Eigen::Matrix<double, nsample, 1> pulse;
//      pulse.ResizeTo(nsample);
      for (unsigned int isample = firstsamplet; isample<nsample; ++isample) {
        pulse(isample) = fullpulse(isample+offset);
        pedSubSamples(isample) = std::max(0., 
            pedSubSamples(isample) - amplitudes[ipulse]*pulse(isample)/pulsenorm);
      }
    }

    // Compute parameters
    double amplitude_(-1.), jitter_(-1.);
    auto param = iGainSwitch==0
        ? weights_0 * pedSubSamples
        : weights_1 * pedSubSamples;
//    ROOT::Math::SVector <double,3> param = (*(weights[iGainSwitch])) * pedSubSamples;
    amplitude_ = param(
        /*EcalUncalibRecHitRecAbsAlgo<EcalDataFrame>::iAmplitude*/ 0);
    if (amplitude_) 
        jitter_ = -param(
            /*EcalUncalibRecHitRecAbsAlgo<EcalDataFrame>::iTime*/ 2) / amplitude_;
    else 
        jitter_ = 0.;
  
    return jitter_;
  }

};

}}
#endif
