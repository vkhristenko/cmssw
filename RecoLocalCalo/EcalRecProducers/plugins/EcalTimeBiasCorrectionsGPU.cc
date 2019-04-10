#include "EcalTimeBiasCorrectionsGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

EcalTimeBiasCorrectionsGPU::EcalTimeBiasCorrectionsGPU(
        EcalTimeBiasCorrections const& values) 
    : EBTimeCorrAmplitudeBins_{values.EBTimeCorrAmplitudeBins}
    , EBTimeCorrShiftBins_{values.EBTimeCorrShiftBins}
    , EETimeCorrAmplitudeBins_{values.EETimeCorrAmplitudeBins}
    , EETimeCorrShiftBins_{values.EETimeCorrShiftBins}
{}

EcalTimeBiasCorrectionsGPU::Product::~Product() {
    // deallocation
    cudaCheck( cudaFree(EBTimeCorrAmplitudeBins) );
    cudaCheck( cudaFree(EBTimeCorrShiftBins) );
    cudaCheck( cudaFree(EETimeCorrAmplitudeBins) );
    cudaCheck( cudaFree(EETimeCorrShiftBins) );
}

EcalTimeBiasCorrectionsGPU::Product const& EcalTimeBiasCorrectionsGPU::getProduct(
        cuda::stream_t<>& cudaStream) const
{
    auto const& product = product_.dataForCurrentDeviceAsync(cudaStream,
        [this](EcalTimeBiasCorrectionsGPU::Product& product, cuda::stream_t<>& cudaStream) {
            // malloc
            cudaCheck( cudaMalloc((void**)&product.EBTimeCorrAmplitudeBins,
                                  this->EBTimeCorrAmplitudeBins_.size() * 
                                  sizeof(float)) );
            cudaCheck( cudaMalloc((void**)&product.EBTimeCorrShiftBins,
                                  this->EBTimeCorrShiftBins_.size() * 
                                  sizeof(float)) );
            cudaCheck( cudaMalloc((void**)&product.EETimeCorrAmplitudeBins,
                                  this->EETimeCorrAmplitudeBins_.size() * 
                                  sizeof(float)) );
            cudaCheck( cudaMalloc((void**)&product.EETimeCorrShiftBins,
                                  this->EETimeCorrShiftBins_.size() * 
                                  sizeof(float)) );
            // transfer 
            cudaCheck( cudaMemcpyAsync(product.EBTimeCorrAmplitudeBins,
                                       this->EBTimeCorrAmplitudeBins_.data(),
                                       this->EBTimeCorrAmplitudeBins_.size() * 
                                       sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
            cudaCheck( cudaMemcpyAsync(product.EBTimeCorrShiftBins,
                                       this->EBTimeCorrShiftBins_.data(),
                                       this->EBTimeCorrShiftBins_.size() * 
                                       sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
            cudaCheck( cudaMemcpyAsync(product.EETimeCorrAmplitudeBins,
                                       this->EETimeCorrAmplitudeBins_.data(),
                                       this->EETimeCorrAmplitudeBins_.size() * 
                                       sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
            cudaCheck( cudaMemcpyAsync(product.EETimeCorrShiftBins,
                                       this->EETimeCorrShiftBins_.data(),
                                       this->EETimeCorrShiftBins_.size() * 
                                       sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
        }
    );

    return product;
}

TYPELOOKUP_DATA_REG(EcalTimeBiasCorrectionsGPU);
