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
      //                              ecal::SoAUncalibratedRecHitCollection const* uncalibRechit_in,
      //                              ecal::SoARecHitCollection* Rechit_out,
      float const* amplitude,  // in adc counts
      type_wrapper<reco::StorageScalarType, Tag::soa>::type * energy,          // in energy
      int nchannels
    ); 
    
    
    // host version, to be called by the plugin
    void create_ecal_rehit(
      float const* amplitude,  // in adc counts
      type_wrapper<reco::StorageScalarType, Tag::soa>::type * energy,          // in energy
      int nchannels
    );
    
  }
  
}

