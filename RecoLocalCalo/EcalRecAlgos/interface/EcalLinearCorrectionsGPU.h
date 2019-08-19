#ifndef RecoLocalCalo_EcalRecProducers_src_EcalLinearCorrectionsGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalLinearCorrectionsGPU_h

#include "CondFormats/EcalObjects/interface/EcalLinearCorrections.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

#include <cuda/api_wrappers.h>

class EcalLinearCorrectionsGPU {
public:
  struct Product {
    ~Product();
    float *p1=nullptr;
    float *p2=nullptr;
    float *p3=nullptr;
    edm::Timestamp *t1=nullptr;
    edm::Timestamp *t2=nullptr;
    edm::Timestamp *t3=nullptr;
  };
  
  #ifndef __CUDACC__
  
  // 
  EcalLinearCorrectionsGPU(EcalLinearCorrections const&);
  
  // will call dealloation for Product thru ~Product
  ~EcalLinearCorrectionsGPU() = default;
  
  // get device pointers
  Product const& getProduct(cuda::stream_t<>&) const;
  
  // 
  static std::string name() { return std::string{"ecalLinearCorrectionsGPU"}; }
  
private:
  // in the future, we need to arrange so to avoid this copy on the host
  // store eb first then ee
  std::vector<float, CUDAHostAllocator<float>> p1_;
  std::vector<float, CUDAHostAllocator<float>> p2_;
  
  CUDAESProduct<Product> product_;
  
  #endif
};


#endif



// public:
//   struct EcalLaserAPDPNpair{
//     EcalLaserAPDPNpair() : p1(0), p2(0), p3(0) {}
//     float p1;
//     float p2;
//     float p3;
//     
//     COND_SERIALIZABLE;
//   };
//   struct EcalLaserTimeStamp{
//     EcalLaserTimeStamp() : t1(), t2(), t3() {}
//     edm::Timestamp t1;
//     edm::Timestamp t2;
//     edm::Timestamp t3;
//     
//     COND_SERIALIZABLE;
//   };
//   
//   typedef EcalCondObjectContainer<EcalLaserAPDPNpair> EcalLinearCorrectionsMap;
//   typedef std::vector<EcalLaserTimeStamp> EcalLaserTimeStampMap;
//   
//   EcalLinearCorrections() : time_map(92) {}; // FIXME
//   ~EcalLinearCorrections() {};
//   
//   void  setValue(uint32_t rawId, const EcalLaserAPDPNpair& value) { laser_map[rawId] = value; };
//   const EcalLinearCorrectionsMap& getLaserMap() const { return laser_map; }
//   
//   void setTime(int hashedIndex, const EcalLaserTimeStamp& value) { time_map[hashedIndex] = value; };
//   const EcalLaserTimeStampMap& getTimeMap() const { return time_map; }  
//   
//   private:
//     EcalLinearCorrectionsMap laser_map;
//     EcalLaserTimeStampMap time_map;
    
    