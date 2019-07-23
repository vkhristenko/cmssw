#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedPedestalsGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include <cmath>

namespace {
    float convert(
            float const x, int const i,
            HcalQIECoder const& coder, HcalQIEShape const& shape) {
        int const x1 = static_cast<int>(std::floor(x));
        int const x2 = static_cast<int>(std::floor(x+1));
        float const y2 = coder.charge(shape, x2, i);
        float const y1 = coder.charge(shape, x1, i);
        return (y2 - y1) / (x - x1) + y1;
    }
}

// FIXME: add proper getters to conditions
HcalConvertedPedestalsGPU::HcalConvertedPedestalsGPU(
        HcalPedestals const& pedestals, HcalQIEData const& qieData,
        HcalQIETypes const& qieTypes) 
    : totalChannels_{pedestals.getAllContainers()[0].second.size()
        + pedestals.getAllContainers()[1].second.size()}
    , offsetForHashes_{static_cast<uint32_t>(pedestals.getAllContainers()[0].second.size())}
    , values_(totalChannels_*4)
{
#ifdef HCAL_MAHI_CPUDEBUG
    std::cout << "hello from converted pedestals" << std::endl;
#endif

    // have to convert to fc if stored in adc
    auto const unitIsADC = pedestals.isADC();

    // fill in barrel
    auto const& pedestalBarrelValues = pedestals.getAllContainers()[0].second;
    auto const& qieDataBarrelValues = qieData.getAllContainers()[0].second;
    auto const& qieTypesBarrelValues = qieTypes.getAllContainers()[0].second;

#ifdef HCAL_MAHI_CPUDEBUG
    assert(pedestalBarrelValues.size() == qieDataBarrelValues.size());
    assert(pedestalBarrelValues.size() == qieTypesBarrelValues.size());
#endif

    for (uint64_t i=0; i<pedestalBarrelValues.size(); ++i) {
        auto const& qieCoder = qieDataBarrelValues[i];
        auto const qieType = qieTypesBarrelValues[i].getValue();
        auto const& qieShape = qieData.getShape(qieType);

        values_[i*4] = unitIsADC 
            ? convert(pedestalBarrelValues[i].getValue(0), 0, qieCoder, qieShape)
            : pedestalBarrelValues[i].getValue(0);
        values_[i*4 + 1] = unitIsADC 
            ? convert(pedestalBarrelValues[i].getValue(1), 1, qieCoder, qieShape)
            : pedestalBarrelValues[i].getValue(1);
        values_[i*4 + 2] = unitIsADC 
            ? convert(pedestalBarrelValues[i].getValue(2), 2, qieCoder, qieShape)
            : pedestalBarrelValues[i].getValue(2);
        values_[i*4 + 3] = unitIsADC 
            ? convert(pedestalBarrelValues[i].getValue(3), 3, qieCoder, qieShape)
            : pedestalBarrelValues[i].getValue(3);
    }

    // fill in endcap
    auto const& pedestalEndcapValues = pedestals.getAllContainers()[1].second;
    auto const& qieDataEndcapValues = qieData.getAllContainers()[1].second;
    auto const& qieTypesEndcapValues = qieTypes.getAllContainers()[1].second;

#ifdef HCAL_MAHI_CPUDEBUG
    assert(pedestalEndcapValues.size() == qieDataEndcapValues.size());
    assert(pedestalEndcapValues.size() == qieTypesEndcapValues.size());
#endif

    for (uint64_t i=0; i<pedestalEndcapValues.size(); ++i) {
        auto const& qieCoder = qieDataEndcapValues[i];
        auto const qieType = qieTypesEndcapValues[i].getValue();
        auto const& qieShape = qieData.getShape(qieType);

        values_[i*4] = unitIsADC 
            ? convert(pedestalEndcapValues[i].getValue(0), 0, qieCoder, qieShape)
            : pedestalEndcapValues[i].getValue(0);
        values_[i*4 + 1] = unitIsADC 
            ? convert(pedestalEndcapValues[i].getValue(1), 1, qieCoder, qieShape)
            : pedestalEndcapValues[i].getValue(1);
        values_[i*4 + 2] = unitIsADC 
            ? convert(pedestalEndcapValues[i].getValue(2), 2, qieCoder, qieShape)
            : pedestalEndcapValues[i].getValue(2);
        values_[i*4 + 3] = unitIsADC 
            ? convert(pedestalEndcapValues[i].getValue(3), 3, qieCoder, qieShape)
            : pedestalEndcapValues[i].getValue(3);
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
