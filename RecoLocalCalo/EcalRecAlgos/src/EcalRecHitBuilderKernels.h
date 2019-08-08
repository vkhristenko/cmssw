//
// Builder of ECAL RecHits on GPU
//

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/DeclsForKernels.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"


namespace ecal { 
  namespace rechit {
    
    __global__
    void kernel_create_ecal_rehit(
      uint32_t const* dids_eb,
      uint32_t const* dids_ee,
      ::ecal::reco::StorageScalarType const* amplitude_eb,   // in adc counts  
      ::ecal::reco::StorageScalarType const* amplitude_ee,   // in adc counts  
      ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
      int const nchannels,
      uint32_t const offsetForInput
    );
      
    
    // host version, to be called by the plugin
    
    void create_ecal_rehit(
      EventInputDataGPU const& eventInputGPU,
      EventOutputDataGPU&      eventOutputGPU,
      //     eventDataForScratchGPU_,
      //     conditions,
      //     configParameters_,
      uint32_t const offsetForInput, 
      cuda::stream_t<>& cudaStream
    );
    
  }
  
}

