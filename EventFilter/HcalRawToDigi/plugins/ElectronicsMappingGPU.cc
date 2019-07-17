#include "EventFilter/HcalRawToDigi/plugins/ElectronicsMappingGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"

namespace hcal { namespace raw {

// TODO: 0x3FFFFF * 4B ~= 16MB
// tmp solution for linear mapping of eid -> did
ElectronicsMappingGPU::ElectronicsMappingGPU(HcalElectronicsMap const& mapping) 
    : eid2did_(HcalElectronicsId::maxLinearIndex)
{  
    auto const& eids = mapping.allElectronicsIdPrecision();
    for (uint32_t i=0; i<eids.size(); ++i) {
        auto const& eid = eids[i];
        eid2did_[eid.linearIndex()] = mapping.lookup(eid).rawId();
    }
}

ElectronicsMappingGPU::Product::~Product() {
    // deallocation
    cudaCheck( cudaFree(eid2did) );
}

ElectronicsMappingGPU::Product const& ElectronicsMappingGPU::getProduct(
        cuda::stream_t<>& cudaStream) const
{
    auto const& product = product_.dataForCurrentDeviceAsync(cudaStream,
        [this](ElectronicsMappingGPU::Product& product, cuda::stream_t<>& cudaStream) {
            // malloc
            cudaCheck( cudaMalloc((void**)&product.eid2did,
                                  this->eid2did_.size() * sizeof(uint32_t)) );

            // transfer 
            cudaCheck( cudaMemcpyAsync(product.eid2did,
                                       this->eid2did_.data(),
                                       this->eid2did_.size() * sizeof(uint32_t),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
        }
    );

    return product;
}

}}

TYPELOOKUP_DATA_REG(hcal::raw::ElectronicsMappingGPU);
