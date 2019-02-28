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
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"
//#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"

#include <cuda.h>

class EcalDigiCollection;
class EcalPedestal;
class EcalMGPAGainRatio;
class EcalXtalGroupId;
class EcalPulseShape;
class EcalPulseCovariance;
class EcalSampleMask;
class EcalTimeBiasCorrections;

namespace ecal { namespace multifit { namespace v1 {

enum TimeAlgo {noMethod, ratioMethod, weightsMethod};

using EMatrix = Eigen::Matrix<double,
    EcalWeightSet::EcalWeightMatrix::rep_type::kRows,
    EcalWeightSet::EcalWeightMatrix::rep_type::kCols>;

struct device_data {
    uint16_t *digis_data = nullptr;
    uint32_t *ids = nullptr;
    SampleVector* amplitudes = nullptr;
    SampleGainVector* gainsNoise = nullptr;
    SampleGainVector* gainsPedestal = nullptr;
//    EcalPedestal *pedestals = nullptr;
    float* mean_x12 = nullptr;
    float* rms_x12 = nullptr;
    float* mean_x6 = nullptr;
    float* rms_x6 = nullptr;
    float* mean_x1 = nullptr;
    float* rms_x1 = nullptr;
//    EcalMGPAGainRatio *gains = nullptr;
    float* gain12Over6 = nullptr;
    float* gain6Over1 = nullptr;
    EcalXtalGroupId *xtals = nullptr;
    EcalPulseShape *pulses = nullptr;
    FulllSampleVector* epulses = nullptr;
    EcalPulseCovariance *covariances = nullptr;
    FullSampleMatrix* pulse_covariances = nullptr;
    EcalUncalibratedRecHit *rechits = nullptr;
    SampleMatrix *noisecorrs = nullptr; // array of 3
    EcalSampleMask *sample_mask = nullptr;
    SampleMatrix* noisecov = nullptr;
    PulseMatrixType* pulse_matrix = nullptr;
    BXVectorType* bxs = nullptr;
    float *EBTimeCorrAmplitudeBins = nullptr;
    int EBTimeCorrAmplitudeBins_size;
    float *EBTimeCorrShiftBins = nullptr;
    int EBTimeCorrShiftBins_size;
    float *EETimeCorrAmplitudeBins = nullptr;
    int EETimeCorrAmplitudeBins_size;
    float *EETimeCorrShiftBins = nullptr;
    int EETimeCorrShiftBins_size;
    EMatrix *weights = nullptr;
};

struct xyz {
    int x,y,z;
};

struct conf_data {
    xyz threads;
};

struct pedestal_data {
    std::vector<float> mean_x12;
    std::vector<float> rms_x12;
    std::vector<float> mean_x6;
    std::vector<float> rms_x6;
    std::vector<float> mean_x1;
    std::vector<float> rms_x1;
};

struct mgpagain_ratio_data {
    std::vector<float> gain12Over6;
    std::vector<float> gain6Over1;
};

struct host_data {
    EcalDigiCollection const *digis;
    EcalUncalibratedRecHitCollection *rechits;
//    std::vector<EcalPedestal> const *pedestals;
    pedestal_data const& ped_data,
//    std::vector<EcalMGPAGainRatio> const *gains;
    mgpagain_ratio_data const& gainratio_data,
    std::vector<EcalXtalGroupId> const *xtals;
    std::vector<EcalPulseShape> const *pulse_shapes;
    std::vector<EcalPulseCovariance> const *pulse_covariances;
    SampleMatrixGainArray const *noisecorrs;
    EcalSampleMask const *sample_mask;
    EcalTimeBiasCorrections const *time_bias_corrections;
    std::vector<EMatrix> const* weights;
    BXVectorType const* bxs;
};

void scatter(host_data&, device_data&, conf_data const&);

/*
void scatter(EcalDigiCollection const&,
             EcalUncalibratedRecHitCollection&,
             std::vector<EcalPedestal> const&,
             std::vector<EcalMGPAGainRatio> const&,
             std::vector<EcalXtalGroupId> const&,
             std::vector<EcalPulseShape> const&,
             std::vector<EcalPulseCovariance> const&,
             SampleMatrixGainArray const&,
             device_data&);
*/

}}}

#endif
