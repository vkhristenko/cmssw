#ifndef RecoLocalCalo_EcalRecAlgos_interfaceEcalMultifitParametersGPU_h
#define RecoLocalCalo_EcalRecAlgos_interfaceEcalMultifitParametersGPU_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/HostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/ESProduct.h"
#endif

class EcalMultifitParametersGPU {
public:
    struct Product {
        ~Product();
        double *amplitudeFitParametersEB, *amplitudeFitParametersEE,
               *timeFitParametersEB, *timeFitParametersEE;
    };

#ifndef __CUDACC__
    EcalMultifitParametersGPU(edm::ParameterSet const&);

    ~EcalMultifitParametersGPU() = default;

    Product const& getProduct(cudaStream_t) const;

private:
    std::vector<double, cms::cuda::HostAllocator<double>> amplitudeFitParametersEB_, amplitudeFitParametersEE_, timeFitParametersEB_, timeFitParametersEE_;

    cms::cuda::ESProduct<Product> product_;
#endif
};

#endif // RecoLocalCalo_EcalRecAlgos_interfaceEcalMultifitParametersGPU_h
