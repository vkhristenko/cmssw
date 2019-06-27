#ifndef RecoLocalCalo_EcalRecAlgos_interface_DeclsForKernels_h
#define RecoLocalCalo_EcalRecAlgos_interface_DeclsForKernels_h

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/RecoTypes.h"

#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalPedestalsGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalGainRatiosGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalPulseShapesGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalPulseCovariancesGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSamplesCorrelationGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalTimeBiasCorrectionsGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalTimeCalibConstantsGPU.h"

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

// event input data on cpu, just const refs
struct EventInputDataCPU {
    EBDigiCollection const& ebDigis;
    EEDigiCollection const& eeDigis;
};

//
struct EventInputDataGPU {
    uint16_t *digis;
    uint32_t *ids;

    void allocate(uint32_t size) {
        cudaCheck( cudaMalloc((void**)&digis,
            sizeof(uint16_t) * size * EcalDataFrame::MAXSAMPLES) );
        cudaCheck( cudaMalloc((void**)&ids,
            sizeof(uint32_t) * size) );
    }

    void deallocate() {
        cudaCheck( cudaFree(digis) );
        cudaCheck( cudaFree(ids) );
    }
};

// parameters have a fixed type
// Can we go by with single precision
struct ConfigurationParameters {
    using type = double;
    // device ptrs
    type *amplitudeFitParametersEB=nullptr, *amplitudeFitParametersEE=nullptr;

    uint32_t timeFitParametersSizeEB, timeFitParametersSizeEE;
    // device ptrs
    type *timeFitParametersEB=nullptr, *timeFitParametersEE=nullptr;

    type timeFitLimitsFirstEB, timeFitLimitsFirstEE;
    type timeFitLimitsSecondEB, timeFitLimitsSecondEE;

    type timeConstantTermEB, timeConstantTermEE;

    type timeNconstEB, timeNconstEE;

    type amplitudeThreshEE, amplitudeThreshEB;

    type outOfTimeThreshG12pEB, outOfTimeThreshG12mEB;
    type outOfTimeThreshG12pEE, outOfTimeThreshG12mEE;
    type outOfTimeThreshG61pEE, outOfTimeThreshG61mEE;
    type outOfTimeThreshG61pEB, outOfTimeThreshG61mEB;

    std::array<uint32_t, 3> kernelMinimizeThreads;

    bool shouldRunTimingComputation;
};

struct EventOutputDataGPU final : public ::ecal::UncalibratedRecHit<::ecal::Tag::ptr> 
{
    void allocate(ConfigurationParameters const& configParameters, uint32_t size) {
        cudaCheck( cudaMalloc((void**)&amplitudesAll,
            size * sizeof(SampleVector)) );
        cudaCheck( cudaMalloc((void**)&amplitude,
            size * sizeof(::ecal::reco::StorageScalarType)) );
        cudaCheck( cudaMalloc((void**)&chi2,
            size * sizeof(::ecal::reco::StorageScalarType)) );
        cudaCheck( cudaMalloc((void**)&pedestal,
            size * sizeof(::ecal::reco::StorageScalarType)) );

        if (configParameters.shouldRunTimingComputation) {
            cudaCheck( cudaMalloc((void**)&jitter,
                size * sizeof(::ecal::reco::StorageScalarType)) );
            cudaCheck( cudaMalloc((void**)&jitterError,
                size * sizeof(::ecal::reco::StorageScalarType)) );
        }

        cudaCheck( cudaMalloc((void**)&did,
            size * sizeof(uint32_t)) );
        cudaCheck( cudaMalloc((void**)&flags,
            size * sizeof(uint32_t)) );
    }

    void deallocate(ConfigurationParameters const& configParameters) {
        cudaCheck( cudaFree(amplitudesAll) );
        cudaCheck( cudaFree(amplitude) );
        cudaCheck( cudaFree(chi2) );
        cudaCheck( cudaFree(pedestal) );
        if (configParameters.shouldRunTimingComputation) {
            cudaCheck( cudaFree(jitter) );
            cudaCheck( cudaFree(jitterError) );
        }
        cudaCheck( cudaFree(did) );
        cudaCheck( cudaFree(flags) );
    }
};

struct EventDataForScratchGPU {
    SampleVector *samples = nullptr;
    SampleGainVector *gainsNoise = nullptr;

    uint32_t *v2rmapping_1=nullptr;
    uint32_t *v2rmapping_2=nullptr;
    uint32_t *pChannelsCounter = nullptr;

    SampleVector::Scalar *decompMatrixMainLoop = nullptr, *decompMatrixFnnls=nullptr;
    SampleVector::Scalar *AtA=nullptr, *Atb=nullptr;
    char *samplesMapping=nullptr, *npassive=nullptr;
    ::ecal::reco::StorageScalarType *chi2_prev=nullptr;

    bool *hasSwitchToGain6=nullptr,
         *hasSwitchToGain1=nullptr,
         *isSaturated=nullptr;

    SampleVector::Scalar *sample_values, *sample_value_errors;
    bool *useless_sample_values;
    SampleVector::Scalar* chi2sNullHypot;
    SampleVector::Scalar* sum0sNullHypot;
    SampleVector::Scalar* sumAAsNullHypot;
    char* pedestal_nums;
    SampleVector::Scalar *tMaxAlphaBetas, *tMaxErrorAlphaBetas;
    SampleVector::Scalar *accTimeMax, *accTimeWgt;
    SampleVector::Scalar *ampMaxAlphaBeta, *ampMaxError;
    SampleVector::Scalar *timeMax, *timeError;
    TimeComputationState *tcState;

    void allocate(ConfigurationParameters const& configParameters, uint32_t size) {
        cudaCheck( cudaMalloc((void**)&samples,
            size * sizeof(SampleVector)) );
        cudaCheck( cudaMalloc((void**)&gainsNoise,
            size * sizeof(SampleGainVector)) );

        cudaCheck( cudaMalloc((void**)&v2rmapping_1,
            size * sizeof(uint32_t)) );
        cudaCheck( cudaMalloc((void**)&v2rmapping_2,
            size * sizeof(uint32_t)) );
        cudaCheck( cudaMalloc((void**)&pChannelsCounter,
            sizeof(uint32_t)) );
        // FIXME: replace 55 with MapSymM::total
        cudaCheck( cudaMalloc((void**)&decompMatrixMainLoop,
            size * 55 * sizeof(SampleVector::Scalar)) );
        cudaCheck( cudaMalloc((void**)&decompMatrixFnnls,
            size * 55 * sizeof(SampleVector::Scalar)) );
        cudaCheck( cudaMalloc((void**)&AtA,
            size * 55 * sizeof(SampleVector::Scalar)) );
        cudaCheck( cudaMalloc((void**)&Atb,
            size * 10 * sizeof(SampleVector::Scalar)) );
        cudaCheck( cudaMalloc((void**)&samplesMapping,
            size * 10 * sizeof(char)) );
        cudaCheck( cudaMalloc((void**)&npassive,
            size * sizeof(char)) );
        cudaCheck( cudaMalloc((void**)&chi2_prev,
            size * sizeof(::ecal::reco::StorageScalarType)) );

        cudaCheck( cudaMalloc((void**)&hasSwitchToGain6,
            size * sizeof(bool)) );
        cudaCheck( cudaMalloc((void**)&hasSwitchToGain1,
            size * sizeof(bool)) );
        cudaCheck( cudaMalloc((void**)&isSaturated,
            size * sizeof(bool)) );

        if (configParameters.shouldRunTimingComputation) {
            cudaCheck( cudaMalloc((void**)&sample_values,
                size * sizeof(SampleVector)) );
            cudaCheck( cudaMalloc((void**)&sample_value_errors,
                size * sizeof(SampleVector)) );
            cudaCheck( cudaMalloc((void**)&useless_sample_values,
                size * sizeof(bool) * EcalDataFrame::MAXSAMPLES) );
            cudaCheck( cudaMalloc((void**)&chi2sNullHypot,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&sum0sNullHypot,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&sumAAsNullHypot,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&pedestal_nums,
                size * sizeof(char)) );

            cudaCheck( cudaMalloc((void**)&tMaxAlphaBetas,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&tMaxErrorAlphaBetas,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&accTimeMax,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&accTimeWgt,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&ampMaxAlphaBeta,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&ampMaxError,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&timeMax,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&timeError,
                size * sizeof(SampleVector::Scalar)) );
            cudaCheck( cudaMalloc((void**)&tcState,
                size * sizeof(TimeComputationState)) );
        }
    }

    void deallocate(ConfigurationParameters const& configParameters) {
        cudaCheck( cudaFree(samples) );
        cudaCheck( cudaFree(gainsNoise) );

        cudaCheck( cudaFree(v2rmapping_1) );
        cudaCheck( cudaFree(v2rmapping_2) );
        cudaCheck( cudaFree(pChannelsCounter) );

        cudaCheck( cudaFree(decompMatrixMainLoop) );
        cudaCheck( cudaFree(decompMatrixFnnls) );
        cudaCheck( cudaFree(AtA) );
        cudaCheck( cudaFree(Atb) );
        cudaCheck( cudaFree(samplesMapping) );
        cudaCheck( cudaFree(npassive) );
        cudaCheck( cudaFree(chi2_prev) );

        cudaCheck( cudaFree(hasSwitchToGain6) );
        cudaCheck( cudaFree(hasSwitchToGain1) );
        cudaCheck( cudaFree(isSaturated) );

        if (configParameters.shouldRunTimingComputation) {
            cudaCheck( cudaFree(sample_values) );
            cudaCheck( cudaFree(sample_value_errors) );
            cudaCheck( cudaFree(useless_sample_values) );
            cudaCheck( cudaFree(chi2sNullHypot) );
            cudaCheck( cudaFree(sum0sNullHypot) );
            cudaCheck( cudaFree(sumAAsNullHypot) );
            cudaCheck( cudaFree(pedestal_nums) );

            cudaCheck( cudaFree(tMaxAlphaBetas) );
            cudaCheck( cudaFree(tMaxErrorAlphaBetas) );
            cudaCheck( cudaFree(accTimeMax) );
            cudaCheck( cudaFree(accTimeWgt) );
            cudaCheck( cudaFree(ampMaxAlphaBeta) );
            cudaCheck( cudaFree(ampMaxError) );
            cudaCheck( cudaFree(timeMax) );
            cudaCheck( cudaFree(timeError) );
            cudaCheck( cudaFree(tcState) );
        }
    }
};

// const refs products to conditions
struct ConditionsProducts {
    EcalPedestalsGPU::Product const& pedestals;
    EcalGainRatiosGPU::Product const& gainRatios;
    EcalPulseShapesGPU::Product const& pulseShapes;
    EcalPulseCovariancesGPU::Product const& pulseCovariances;
    EcalSamplesCorrelationGPU::Product const& samplesCorrelation;
    EcalTimeBiasCorrectionsGPU::Product const& timeBiasCorrections;
    EcalTimeCalibConstantsGPU::Product const& timeCalibConstants;
    EcalSampleMask const& sampleMask;
    EcalTimeOffsetConstant const& timeOffsetConstant;
    uint32_t offsetForHashes;
};

//*/

struct xyz {
    int x,y,z;
};

struct conf_data {
    xyz threads;
    bool runV1;
    cudaStream_t cuStream;
};

}}

#endif
