#include "RecoLocalCalo/EcalRecAlgos/interface/EcalLinearCorrectionsGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

EcalLinearCorrectionsGPU::EcalLinearCorrectionsGPU(EcalLinearCorrections const& values) 
    : p1_(values.getValueMap().size())
    , p2_(values.getValueMap().size())
{
//     // fill in eb
//     auto const& barrelValues = values.barrelItems();
//     for (unsigned int i=0; i<barrelValues.size(); i++) {
//         p1_[i] = barrelValues[i].p1();
//         p2_[i] = barrelValues[i].p2();
//     }
//     
//     // fill in ee
//     auto const& endcapValues = values.endcapItems();
//     auto const offset = barrelValues.size();
//     for (unsigned int i=0; i<endcapValues.size(); i++) {
//         p1_[offset + i] = endcapValues[i].p1();
//         p2_[offset + i] = endcapValues[i].p2();
//     }
}

EcalLinearCorrectionsGPU::Product::~Product() {
    // deallocation
    cudaCheck( cudaFree(p1) );
    cudaCheck( cudaFree(p2) );
}

EcalLinearCorrectionsGPU::Product const& EcalLinearCorrectionsGPU::getProduct(
        cuda::stream_t<>& cudaStream) const
{
    auto const& product = product_.dataForCurrentDeviceAsync(cudaStream,
        [this](EcalLinearCorrectionsGPU::Product& product, cuda::stream_t<>& cudaStream) {
            // malloc
            cudaCheck( cudaMalloc((void**)&product.p1,
                                  this->p1_.size() * sizeof(float)) );
            cudaCheck( cudaMalloc((void**)&product.p2,
                                  this->p2_.size() * sizeof(float)) );
            // transfer 
            cudaCheck( cudaMemcpyAsync(product.p1,
                                       this->p1_.data(),
                                       this->p1_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
            cudaCheck( cudaMemcpyAsync(product.p2,
                                       this->p2_.data(),
                                       this->p2_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
        }
    );

    return product;
}

TYPELOOKUP_DATA_REG(EcalLinearCorrectionsGPU);
