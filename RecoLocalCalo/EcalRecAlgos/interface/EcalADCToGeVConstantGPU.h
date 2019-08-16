#ifndef RecoLocalCalo_EcalRecProducers_src_EcalADCToGeVConstantGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalADCToGeVConstantGPU_h

#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

#include <cuda/api_wrappers.h>

class EcalADCToGeVConstantGPU {
public:
  struct Product {
    ~Product();
    float *adc2gev = nullptr;
  };
  
#ifndef __CUDACC__
  
  // 
  EcalADCToGeVConstantGPU(EcalADCToGeVConstant const&);
  
  // will call dealloation for Product thru ~Product
  ~EcalADCToGeVConstantGPU() = default;
  
  // get device pointers
  Product const& getProduct(cuda::stream_t<>&) const;
  
  // 
  static std::string name() { return std::string{"ecalADCToGeVConstantGPU"}; }
  
private:
  // in the future, we need to arrange so to avoid this copy on the host
  // store eb first then ee
  std::vector<float, CUDAHostAllocator<float>> adc2gev_;
  
  CUDAESProduct<Product> product_;
  
#endif
};


#endif
