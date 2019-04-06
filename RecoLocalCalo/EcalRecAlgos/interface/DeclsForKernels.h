#ifndef RecoLocalCalo_EcalRecAlgos_interface_DeclsForKernels_h
#define RecoLocalCalo_EcalRecAlgos_interface_DeclsForKernels_h

#include <vector>

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"
#include "DataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "DataFormats/EcalRecHitSoA/interface/RecoTypes.h"

#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"

class EcalPulseShape;
class EcalSampleMask;
class EcalTimeBiasCorrections;
class EcalPulseCovariance;
class EcalDigiCollection;
class EcalXtalGroupId;
class EcalSamplesCorrelation;
class EBDigiCollection;
class EEDigiCollection;

namespace ecal { namespace multifit {

enum class TimeComputationState : char {
    NotFinished = 0,
    Finished = 1
};
enum class MinimizationState : char {
    NotFinished = 0,
    Finished = 1,
    Precomputed = 2,
};

using EMatrix = Eigen::Matrix<double,
    EcalWeightSet::EcalWeightMatrix::rep_type::kRows,
    EcalWeightSet::EcalWeightMatrix::rep_type::kCols>;

struct device_data {
    uint16_t *digis_data = nullptr;
    uint32_t *ids = nullptr;
    SampleVector* samples = nullptr;
    SampleGainVector* gainsNoise = nullptr;
    SampleGainVector* gainsPedestal = nullptr;
    float* mean_x12 = nullptr;
    float* rms_x12 = nullptr;
    float* mean_x6 = nullptr;
    float* rms_x6 = nullptr;
    float* mean_x1 = nullptr;
    float* rms_x1 = nullptr;
//    EcalMGPAGainRatio *gains = nullptr;
    float* gain12Over6 = nullptr;
    float* gain6Over1 = nullptr;
    float* timeCalibConstants = nullptr;
    EcalPulseShape *pulses = nullptr;
    FullSampleVector* epulses = nullptr;
    EcalPulseCovariance *covariances = nullptr;
    FullSampleMatrix* pulse_covariances = nullptr;
    EcalSampleMask *sample_mask = nullptr;
    SampleMatrix* noisecov = nullptr;
    SampleMatrix* updatedNoiseCovariance = nullptr;
    SampleMatrix* noiseMatrixDecomposition = nullptr;
    PulseMatrixType* pulse_matrix = nullptr;
    BXVectorType* bxs = nullptr;
    BXVectorType* activeBXs = nullptr;
    float *EBTimeCorrAmplitudeBins = nullptr;
    int EBTimeCorrAmplitudeBins_size;
    float *EBTimeCorrShiftBins = nullptr;
    int EBTimeCorrShiftBins_size;
    float *EETimeCorrAmplitudeBins = nullptr;
    int EETimeCorrAmplitudeBins_size;
    float *EETimeCorrShiftBins = nullptr;
    int EETimeCorrShiftBins_size;

    char *acState = nullptr;
    char *minimizationStatePerBlock = nullptr;
    int *npassive = nullptr;
    std::vector<char> h_minimizationStatesPerBlock;

    // rechits
    bool* statuses = nullptr;
    SampleVector* amplitudes = nullptr;
    ::ecal::reco::StorageScalarType* energies = nullptr;
    ::ecal::reco::StorageScalarType* chi2 = nullptr;
    ::ecal::reco::StorageScalarType* pedestal = nullptr;
    uint32_t *flags = nullptr;


    bool* hasSwitchToGain6 = nullptr;
    bool* hasSwitchToGain1 = nullptr;
    bool* isSaturated = nullptr;

    // from timing computation
    SampleVector::Scalar *sample_values, *sample_value_errors;
    bool* useless_sample_values;
    SampleVector::Scalar* chi2sNullHypot;
    SampleVector::Scalar* sum0sNullHypot;
    SampleVector::Scalar* sumAAsNullHypot;
    char* pedestal_nums;

    // TODO: check if we can use __constant__ memory for these guys
    SampleVector::Scalar *amplitudeFitParametersEB, *amplitudeFitParametersEE;

    SampleVector::Scalar *tMaxAlphaBetas, *tMaxErrorAlphaBetas;
    SampleVector::Scalar *accTimeMax, *accTimeWgt;
    SampleVector::Scalar *ampMaxAlphaBeta, *ampMaxError;
    SampleVector::Scalar *timeMax, *timeError;
    float *jitter, *jitterError;
    TimeComputationState *tcState;
    unsigned int timeFitParametersSizeEB, timeFitParametersSizeEE;
    SampleVector::Scalar timeFitLimitsFirstEB, timeFitLimitsSecondEB;
    SampleVector::Scalar timeFitLimitsFirstEE, timeFitLimitsSecondEE;
    // check if cosntant mem is better for these 2
    SampleVector::Scalar *timeFitParametersEB, *timeFitParametersEE;
    SampleVector::Scalar *amplitudeMax;

    double *G12SamplesCorrelationEB, *G6SamplesCorrelationEB, *G1SamplesCorrelationEB;
    double *G12SamplesCorrelationEE, *G6SamplesCorrelationEE, *G1SamplesCorrelationEE;

    // use constant mem?
    SampleVector::Scalar timeConstantTermEB, timeConstantTermEE;

    // time calib constants
    float offsetTimeValue;
    float timeNconstEB, timeNconstEE;
    float amplitudeThreshEE, amplitudeThreshEB;
    float outOfTimeThreshG12pEB, outOfTimeThreshG12mEB;
    float outOfTimeThreshG12pEE, outOfTimeThreshG12mEE;
    float outOfTimeThreshG61pEE, outOfTimeThreshG61mEE;
    float outOfTimeThreshG61pEB, outOfTimeThreshG61mEB;
};

struct xyz {
    int x,y,z;
};

struct conf_data {
    xyz threads;
    bool runV1;
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
    EBDigiCollection const *digisEB;
    EEDigiCollection const *digisEE; 
//    std::vector<EcalPedestal> const *pedestals;
    pedestal_data const& ped_data;
//    std::vector<EcalMGPAGainRatio> const *gains;
    mgpagain_ratio_data const& gainratio_data;
    EcalSampleMask const& sample_mask;
    std::vector<EcalXtalGroupId> const *xtals;
    std::vector<EcalPulseShape> const *pulse_shapes;
    std::vector<EcalPulseCovariance> const *pulse_covariances;
    EcalTimeBiasCorrections const *time_bias_corrections;
    std::vector<float> const& timeCalibConstants;
    BXVectorType const* bxs;
    ecal::UncalibratedRecHit<ecal::Tag::soa>& rechits_soa_eb;
    ecal::UncalibratedRecHit<ecal::Tag::soa>& rechits_soa_ee;
    EcalSamplesCorrelation const* noiseCovariances;
};

}}

#endif
