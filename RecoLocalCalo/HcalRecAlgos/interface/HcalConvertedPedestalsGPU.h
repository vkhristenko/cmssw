#ifndef RecoLocalCalo_HcalRecAlgos_interface_HcalConvertedPedestalsGPU_h
#define RecoLocalCalo_HcalRecAlgos_interface_HcalConvertedPedestalsGPU_h

#include "CondFormats/HcalObjects/interface/HcalPedestals.h"
#include "CondFormats/HcalObjects/interface/HcalQIEData.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

class HcalConvertedPedestalsGPU {
public:
    struct Product {
        ~Product();
        float *values;
    };

#ifndef __CUDACC__
    // FIXME: testing
    HcalConvertedPedestalsGPU() = default;
    // order matters!
    HcalConvertedPedestalsGPU(HcalPedestals const&, HcalQIEData const&);

    // will trigger deallocation of Product thru ~Product
    ~HcalConvertedPedestalsGPU() = default;

    // get device pointers
    Product const& getProduct(cuda::stream_t<>&) const;

    // 
    static std::string name() { return std::string{"hcalConvertedPedestalsGPU"}; }

    uint32_t offsetForHashes() const { return offsetForHashes_; }

private:
    uint64_t totalChannels_;
    uint32_t offsetForHashes_;
    std::vector<float, CUDAHostAllocator<float>> values_;

    CUDAESProduct<Product> product_;
#endif
};

#endif
