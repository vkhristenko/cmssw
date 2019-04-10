#ifndef RecoLocalCalo_EcalRecProducers_src_EcalTimeCalibConstantsGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalTimeCalibConstantsGPU_h

#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"

#include <cuda/api_wrappers.h>

class EcalTimeCalibConstantsGPU {
public:
    struct Product {
        ~Product();
        float *values=nullptr;
    };

    // rearrange pedestals
    EcalTimeCalibConstantsGPU(EcalTimeCalibConstants const&);

    // will call dealloation for Product thru ~Product
    ~EcalTimeCalibConstantsGPU() = default;

    // get device pointers
    Product const& getProduct(cuda::stream_t<>&) const;

    // 
    static std::string name() { return std::string{"ecalTimeCalibConstantsGPU"}; }

private:
    std::vector<float> const& valuesEB_;
    std::vector<float> const& valuesEE_;

    CUDAESProduct<Product> product_;
};


#endif
