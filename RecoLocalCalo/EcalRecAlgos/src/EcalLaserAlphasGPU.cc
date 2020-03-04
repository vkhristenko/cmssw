#include "RecoLocalCalo/EcalRecAlgos/interface/EcalLaserAlphasGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

EcalLaserAlphasGPU::EcalLaserAlphasGPU(EcalLaserAlphas const& values) 
    : valuesEB_{values.barrelItems()}
    , valuesEE_{values.endcapItems()}
{}

EcalLaserAlphasGPU::Product::~Product() {
    // deallocation
    cudaCheck( cudaFree(values) );
}

EcalLaserAlphasGPU::Product const& EcalLaserAlphasGPU::getProduct(cuda::stream_t<>& cudaStream) const {
    auto const& product = product_.dataForCurrentDeviceAsync(cudaStream,
        [this](EcalLaserAlphasGPU::Product& product, cuda::stream_t<>& cudaStream) {
            // malloc
            cudaCheck( cudaMalloc((void**)&product.values,
                                  (this->valuesEB_.size() + this->valuesEE_.size()) * 
                                  sizeof(float)) );

            // offset in floats, not bytes
            auto const offset = this->valuesEB_.size();

            // transfer 
            cudaCheck( cudaMemcpyAsync(product.values,
                                       this->valuesEB_.data(),
                                       this->valuesEB_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
            cudaCheck( cudaMemcpyAsync(product.values + offset,
                                       this->valuesEE_.data(),
                                       this->valuesEE_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
        }
    );

    return product;
}

TYPELOOKUP_DATA_REG(EcalLaserAlphasGPU);
