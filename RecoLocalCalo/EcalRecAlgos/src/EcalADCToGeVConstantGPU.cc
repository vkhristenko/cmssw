#include "RecoLocalCalo/EcalRecAlgos/interface/EcalADCToGeVConstantGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

EcalADCToGeVConstantGPU::EcalADCToGeVConstantGPU(EcalADCToGeVConstant const& values) 
: adc2gev_(2)  // size is 2, one form EB and one for EE
{
  adc2gev_[0] = values.getEBValue();
  adc2gev_[1] = values.getEEValue(); 
}

EcalADCToGeVConstantGPU::Product::~Product() {
  // deallocation
  cudaCheck( cudaFree(adc2gev) );
}

EcalADCToGeVConstantGPU::Product const& EcalADCToGeVConstantGPU::getProduct(cuda::stream_t<>& cudaStream) const {
  auto const& product = product_.dataForCurrentDeviceAsync(
    cudaStream,
    [this](EcalADCToGeVConstantGPU::Product& product, cuda::stream_t<>& cudaStream) {
      // malloc
      cudaCheck( cudaMalloc((void**)&product.adc2gev,
                            this->adc2gev_.size() * sizeof(float)) );
      // transfer 
      cudaCheck( cudaMemcpyAsync(product.adc2gev,
                                 this->adc2gev_.data(),
                                 this->adc2gev_.size() * sizeof(float),
                                 cudaMemcpyHostToDevice,
                                 cudaStream.id()) );
     }
  );
  
  return product;
}

TYPELOOKUP_DATA_REG(EcalADCToGeVConstantGPU);
