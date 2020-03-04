#ifndef RecoLocalCalo_EcalRecProducers_src_EcalLaserAPDPNRatiosGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalLaserAPDPNRatiosGPU_h

#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

#include <cuda/api_wrappers.h>

class EcalLaserAPDPNRatiosGPU {
public:
  struct Product {
    ~Product();
    float *p1=nullptr;
    float *p2=nullptr;
    float *p3=nullptr;
    edm::TimeValue_t *t1=nullptr;
    edm::TimeValue_t *t2=nullptr;
    edm::TimeValue_t *t3=nullptr;
  };
  
  #ifndef __CUDACC__
  
  // 
  EcalLaserAPDPNRatiosGPU(EcalLaserAPDPNRatios const&);
  
  // will call dealloation for Product thru ~Product
  ~EcalLaserAPDPNRatiosGPU() = default;
  
  // get device pointers
  Product const& getProduct(cuda::stream_t<>&) const;
  
  // 
  static std::string name() { return std::string{"ecalLaserAPDPNRatiosGPU"}; }
  
private:
  // in the future, we need to arrange so to avoid this copy on the host
  // store eb first then ee
  std::vector<float, CUDAHostAllocator<float> > p1_;
  std::vector<float, CUDAHostAllocator<float> > p2_;
  std::vector<float, CUDAHostAllocator<float> > p3_;
  
  std::vector<edm::TimeValue_t, CUDAHostAllocator<edm::TimeValue_t> > t1_;
  std::vector<edm::TimeValue_t, CUDAHostAllocator<edm::TimeValue_t> > t2_;
  std::vector<edm::TimeValue_t, CUDAHostAllocator<edm::TimeValue_t> > t3_;

  CUDAESProduct<Product> product_;
  
  #endif
};


#endif




    