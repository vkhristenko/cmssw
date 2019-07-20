#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIECodersGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

HcalQIECodersGPU::HcalQIECodersGPU(HcalQIEData const& qiedata) 
    : totalChannels_{qiedata.getAllContainers()[0].second.size()
        + qiedata.getAllContainers()[1].second.size()}
    , offsets_(totalChannels_ * numValuesPerChannel)
    , slopes_(totalChannels_ * numValuesPerChannel)
{
    // fill in hb
    auto const& barrelValues = qiedata.getAllContainers()[0].second;
    for (uint64_t i=0; i<barrelValues.size(); ++i) {
        for (uint32_t k=0; k<4; k++)
            for (uint32_t l=0; l<4; l++) {
                auto const linear = k*4 + l;
                offsets_[i * numValuesPerChannel + linear] = 
                    barrelValues[i].offset(k, l);
                slopes_[i * numValuesPerChannel + linear] = 
                    barrelValues[i].slope(k, l);
            }

        /*
        offset00_[i] = barrelValues[i].offset(0, 0);
        offset01_[i] = barrelValues[i].offset(0, 1);
        offset02_[i] = barrelValues[i].offset(0, 2);
        offset03_[i] = barrelValues[i].offset(0, 3);
        offset10_[i] = barrelValues[i].offset(1, 0);
        offset11_[i] = barrelValues[i].offset(1, 1);
        offset12_[i] = barrelValues[i].offset(1, 2);
        offset13_[i] = barrelValues[i].offset(1, 3);
        offset20_[i] = barrelValues[i].offset(2, 0);
        offset21_[i] = barrelValues[i].offset(2, 1);
        offset22_[i] = barrelValues[i].offset(2, 2);
        offset23_[i] = barrelValues[i].offset(2, 3);
        offset30_[i] = barrelValues[i].offset(3, 0);
        offset31_[i] = barrelValues[i].offset(3, 1);
        offset32_[i] = barrelValues[i].offset(3, 2);
        offset33_[i] = barrelValues[i].offset(3, 3);
        
        slope00_[i] = barrelValues[i].slope(0, 0);
        slope01_[i] = barrelValues[i].slope(0, 1);
        slope02_[i] = barrelValues[i].slope(0, 2);
        slope03_[i] = barrelValues[i].slope(0, 3);
        slope10_[i] = barrelValues[i].slope(1, 0);
        slope11_[i] = barrelValues[i].slope(1, 1);
        slope12_[i] = barrelValues[i].slope(1, 2);
        slope13_[i] = barrelValues[i].slope(1, 3);
        slope20_[i] = barrelValues[i].slope(2, 0);
        slope21_[i] = barrelValues[i].slope(2, 1);
        slope22_[i] = barrelValues[i].slope(2, 2);
        slope23_[i] = barrelValues[i].slope(2, 3);
        slope30_[i] = barrelValues[i].slope(3, 0);
        slope31_[i] = barrelValues[i].slope(3, 1);
        slope32_[i] = barrelValues[i].slope(3, 2);
        slope33_[i] = barrelValues[i].slope(3, 3);
        */
    }

    // fill in he
    auto const& endcapValues = qiedata.getAllContainers()[1].second;
    auto const offset = barrelValues.size();
    for (uint64_t i=0; i<endcapValues.size(); ++i) {
        auto const off = (i + offset) * numValuesPerChannel;
        for (uint32_t k=0; k<4; k++)
            for (uint32_t l=0; l<4; l++) {
                auto const linear = k*4u + l;
                offsets_[off + linear] = endcapValues[i].offset(k, l);
                slopes_[off + linear] = endcapValues[i].slope(k, l);
            }

        /*
        offset00_[i + offset] = endcapValues[i].offset(0, 0);
        offset01_[i + offset] = endcapValues[i].offset(0, 1);
        offset02_[i + offset] = endcapValues[i].offset(0, 2);
        offset03_[i + offset] = endcapValues[i].offset(0, 3);
        offset10_[i + offset] = endcapValues[i].offset(1, 0);
        offset11_[i + offset] = endcapValues[i].offset(1, 1);
        offset12_[i + offset] = endcapValues[i].offset(1, 2);
        offset13_[i + offset] = endcapValues[i].offset(1, 3);
        offset20_[i + offset] = endcapValues[i].offset(2, 0);
        offset21_[i + offset] = endcapValues[i].offset(2, 1);
        offset22_[i + offset] = endcapValues[i].offset(2, 2);
        offset23_[i + offset] = endcapValues[i].offset(2, 3);
        offset30_[i + offset] = endcapValues[i].offset(3, 0);
        offset31_[i + offset] = endcapValues[i].offset(3, 1);
        offset32_[i + offset] = endcapValues[i].offset(3, 2);
        offset33_[i + offset] = endcapValues[i].offset(3, 3);
        
        slope00_[i + offset] = endcapValues[i].slope(0, 0);
        slope01_[i + offset] = endcapValues[i].slope(0, 1);
        slope02_[i + offset] = endcapValues[i].slope(0, 2);
        slope03_[i + offset] = endcapValues[i].slope(0, 3);
        slope10_[i + offset] = endcapValues[i].slope(1, 0);
        slope11_[i + offset] = endcapValues[i].slope(1, 1);
        slope12_[i + offset] = endcapValues[i].slope(1, 2);
        slope13_[i + offset] = endcapValues[i].slope(1, 3);
        slope20_[i + offset] = endcapValues[i].slope(2, 0);
        slope21_[i + offset] = endcapValues[i].slope(2, 1);
        slope22_[i + offset] = endcapValues[i].slope(2, 2);
        slope23_[i + offset] = endcapValues[i].slope(2, 3);
        slope30_[i + offset] = endcapValues[i].slope(3, 0);
        slope31_[i + offset] = endcapValues[i].slope(3, 1);
        slope32_[i + offset] = endcapValues[i].slope(3, 2);
        slope33_[i + offset] = endcapValues[i].slope(3, 3);
        */
    }
}

HcalQIECodersGPU::Product::~Product() {
    // deallocation
    cudaCheck( cudaFree(offsets) );
    cudaCheck( cudaFree(slopes) );
}

HcalQIECodersGPU::Product const& HcalQIECodersGPU::getProduct(
        cuda::stream_t<>& cudaStream) const {
    auto const& product = product_.dataForCurrentDeviceAsync(cudaStream,
        [this](HcalQIECodersGPU::Product& product, cuda::stream_t<>& cudaStream){
            // malloc
            cudaCheck( cudaMalloc((void**)&product.offsets,
                this->offsets_.size() * sizeof(float)) );
            cudaCheck( cudaMalloc((void**)&product.slopes,
                this->slopes_.size() * sizeof(float)) );
            
            // transfer
            // offset
            cudaCheck( cudaMemcpyAsync(product.offsets, 
                                       this->offsets_.data(),
                                       this->offsets_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
            
            // slope
            cudaCheck( cudaMemcpyAsync(product.slopes, 
                                       this->slopes_.data(),
                                       this->slopes_.size() * sizeof(float),
                                       cudaMemcpyHostToDevice,
                                       cudaStream.id()) );
        });

    return product;
}

TYPELOOKUP_DATA_REG(HcalQIECodersGPU);
