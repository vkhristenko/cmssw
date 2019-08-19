#ifndef RecoLocalCalo_EcalRecProducers_src_EcalLaserAlphasGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalLaserAlphasGPU_h

#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

#include <cuda/api_wrappers.h>

class EcalLaserAlphasGPU {
public:
  struct Product {
    ~Product();
    float *values = nullptr;
  };
  
  #ifndef __CUDACC__
  // 
  EcalLaserAlphasGPU(EcalLaserAlphas const&);
  
  // will call dealloation for Product thru ~Product
  ~EcalLaserAlphasGPU() = default;
  
  // get device pointers
  Product const& getProduct(cuda::stream_t<>&) const;
  
  // TODO: do this centrally
  // get offset for hashes. equals number of barrel items
  uint32_t getOffset() const { return valuesEB_.size(); }
  
  // 
  static std::string name() { return std::string{"ecalLaserAlphasGPU"}; }
  
private:
  std::vector<float> const& valuesEB_;
  std::vector<float> const& valuesEE_;
  
  CUDAESProduct<Product> product_;
  #endif
};


#endif
