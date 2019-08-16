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
      // configuration 
      int const* ChannelStatusToBeExcluded,
      uint32_t ChannelStatusToBeExcludedSize,                     
      // conditions
      float const* adc2gev,
      float const* intercalib,
      uint32_t const* status,
      // input
      uint32_t const* did_eb,
      uint32_t const* did_ee,
      ::ecal::reco::StorageScalarType const* amplitude_eb,   // in adc counts  
      ::ecal::reco::StorageScalarType const* amplitude_ee,   // in adc counts  
      ::ecal::reco::StorageScalarType const* time_eb,   
      ::ecal::reco::StorageScalarType const* time_ee,   
      ::ecal::reco::StorageScalarType const* chi2_eb,   
      ::ecal::reco::StorageScalarType const* chi2_ee,   
      // output
      uint32_t *did,
      ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
      ::ecal::reco::StorageScalarType* time,  
      ::ecal::reco::StorageScalarType* chi2,  
      uint32_t* flagBits,
      uint32_t* extra,
      int const nchannels,
      uint32_t const offsetForInput,
      uint32_t const offsetForHashes  
    );
      
    
    // host version, to be called by the plugin
    
    void create_ecal_rehit(
      EventInputDataGPU const& eventInputGPU,
      EventOutputDataGPU&      eventOutputGPU,
      //     eventDataForScratchGPU_,
      ConditionsProducts const& conditions, 
      ConfigurationParameters const& configParameters,
      uint32_t const offsetForInput, 
      cuda::stream_t<>& cudaStream
    );
    
  }
  
}

