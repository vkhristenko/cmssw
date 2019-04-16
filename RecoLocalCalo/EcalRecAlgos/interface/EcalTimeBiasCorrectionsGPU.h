#ifndef RecoLocalCalo_EcalRecProducers_src_EcalTimeBiasCorrectionsGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalTimeBiasCorrectionsGPU_h

#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

#include <cuda/api_wrappers.h>

class EcalTimeBiasCorrectionsGPU {
public:
    struct Product {
        ~Product();
        float *EBTimeCorrAmplitudeBins, *EBTimeCorrShiftBins;
        float *EETimeCorrAmplitudeBins, *EETimeCorrShiftBins;
    };

#ifndef __CUDACC__
    // rearrange pedestals
    EcalTimeBiasCorrectionsGPU(EcalTimeBiasCorrections const&);

    // will call dealloation for Product thru ~Product
    ~EcalTimeBiasCorrectionsGPU() = default;

    // get device pointers
    Product const& getProduct(cuda::stream_t<>&) const;

    // 
    static std::string name() { return std::string{"ecalTimeBiasCorrectionsGPU"}; }

private:
    std::vector<float> const& EBTimeCorrAmplitudeBins_;
    std::vector<float> const& EBTimeCorrShiftBins_;
    std::vector<float> const& EETimeCorrAmplitudeBins_;
    std::vector<float> const& EETimeCorrShiftBins_;

    CUDAESProduct<Product> product_;
#endif
};


#endif
