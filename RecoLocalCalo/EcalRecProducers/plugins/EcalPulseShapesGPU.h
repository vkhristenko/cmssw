#ifndef RecoLocalCalo_EcalRecProducers_src_EcalPulseShapesGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalPulseShapesGPU_h

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"

#include <cuda/api_wrappers.h>

class EcalPulseShapesGPU {
public:
    struct Product {
        ~Product();
        float *values=nullptr;
    };

    // rearrange pedestals
    EcalPulseShapesGPU(EcalPulseShapes const&);

    // will call dealloation for Product thru ~Product
    ~EcalPulseShapesGPU() = default;

    // get device pointers
    Product const& getProduct(cuda::stream_t<>&) const;

    // 
    static std::string name() { return std::string{"ecalPulseShapesGPU"}; }

private:
    // reuse original vectors (although with default allocator)
    std::vector<EcalPulseShape> const& valuesEB_;
    std::vector<EcalPulseShape> const& valuesEE_;

    CUDAESProduct<Product> product_;
};


#endif
