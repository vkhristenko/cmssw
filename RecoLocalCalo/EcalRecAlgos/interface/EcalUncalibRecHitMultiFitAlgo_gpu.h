#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgo_gpu_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgo_gpu_HH

/** \class EcalUncalibRecHitMultiFitAlgo
  *  Amplitude reconstucted by the multi-template fit
  *
  *  \author J.Bendavid, E.Di Marco
  */

/*

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"


#include "TMatrixDSym.h"
#include "TVectorD.h"

class EcalUncalibRecHitMultiFitAlgo
{
  
 public:
  
  EcalUncalibRecHitMultiFitAlgo();
  ~EcalUncalibRecHitMultiFitAlgo() { };
  EcalUncalibratedRecHit makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX);
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

*/

#include <vector>

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"


class EcalDigiCollection;
class EcalPedestal;
class EcalMGPAGainRatio;
class EcalXtalGroupId;
class EcalPulseShape;
class EcalPulseCovariance;

namespace ecal { namespace multifit {

class EcalUncalibRecHitMultiFitAlgo
{
  
 public:
  
  __device__ EcalUncalibRecHitMultiFitAlgo();
  __device__ EcalUncalibratedRecHit makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX);
  __device__ void disableErrorCalculation() { _computeErrors = false; }
  __device__ void setDoPrefit(bool b) { _doPrefit = b; }
  __device__ void setPrefitMaxChiSq(double x) { _prefitMaxChiSq = x; }
  __device__ void setDynamicPedestals(bool b) { _dynamicPedestals = b; }
  __device__ void setMitigateBadSamples(bool b) { _mitigateBadSamples = b; }
  __device__ void setSelectiveBadSampleCriteria(bool b) { _selectiveBadSampleCriteria = b; }
  __device__ void setAddPedestalUncertainty(double x) { _addPedestalUncertainty = x; }
  __device__ void setSimplifiedNoiseModelForGainSwitch(bool b) { _simplifiedNoiseModelForGainSwitch = b; }
  __device__ void setGainSwitchUseMaxSample(bool b) { _gainSwitchUseMaxSample = b; }
  
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

void scatter(EcalDigiCollection const&,
             EcalUncalibratedRecHitCollection&,
             std::vector<EcalPedestal> const&,
             std::vector<EcalMGPAGainRatio> const&,
             std::vector<EcalXtalGroupId> const&,
             std::vector<EcalPulseShape> const&,
             std::vector<EcalPulseCovariance> const&,
             SampleMatrixGainArray const&);

}}

#endif
