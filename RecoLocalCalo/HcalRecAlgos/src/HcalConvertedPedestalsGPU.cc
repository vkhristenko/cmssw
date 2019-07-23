#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedPedestalsGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

// FIXME: add proper getters to conditions
HcalConvertedPedestalsGPU::HcalConvertedPedestalsGPU(
        HcalPedestals const& pedestals, HcalQIEData const& qieData) 
    : totalChannels_{pedestals.getAllContainers()[0].second.size()
        + pedestals.getAllContainers()[1].second.size()}
    , offsetForHashes_{static_cast<uint32_t>(pedestals.getAllContainers()[0].second.size())}
    , values_(totalChannels_*4)
{
#ifdef HCAL_MAHI_CPUDEBUG
    std::cout << "hello from converted pedestals" << std::endl;
#endif

    // fill in eb
    auto const& barrelValues = pedestals.getAllContainers()[0].second;
    for (uint64_t i=0; i<barrelValues.size(); ++i) {
        values_[i*4] = barrelValues[i].getValue(0);
        values_[i*4 + 1] = barrelValues[i].getValue(1);
        values_[i*4 + 2] = barrelValues[i].getValue(2);
        values_[i*4 + 3] = barrelValues[i].getValue(3);
    }

    // fill in ee
    auto const& endcapValues = pedestals.getAllContainers()[1].second;
    auto const offset = barrelValues.size();
    for (uint64_t i=0; i<endcapValues.size(); ++i) {
        auto const off = offset + i;
        values_[off*4] = endcapValues[i].getValue(0);
        values_[off*4 + 1] = endcapValues[i].getValue(1);
        values_[off*4 + 2] = endcapValues[i].getValue(2);
        values_[off*4 + 3] = endcapValues[i].getValue(3);
    }
}

HcalConvertedPedestalsGPU::Product::~Product() {
    // deallocation
    cudaCheck( cudaFree(values) );    
}

HcalConvertedPedestalsGPU::Product const& HcalConvertedPedestalsGPU::getProduct(
        cuda::stream_t<>& cudaStream) const {
    auto const& product = product_.dataForCurrentDeviceAsync(cudaStream,
        [this](HcalConvertedPedestalsGPU::Product& product, cuda::stream_t<>& cudaStream){
            // malloc
            cudaCheck( cudaMalloc((void**)&product.values, 
                this->values_.size() * sizeof(float)) );
            
            // transfer
            cudaCheck( cudaMemcpyAsync(product.values, 
                                       this->values_.data(),
                                       this->values_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
        });

    return product;
}

TYPELOOKUP_DATA_REG(HcalConvertedPedestalsGPU);
