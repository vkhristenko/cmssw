#include "cuda.h"

#include "KernelHelpers.h"

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

//
//
#include "EcalRecHitBuilderKernels.h"


namespace ecal {
  namespace rechit {
    
    __global__
    void kernel_create_ecal_rehit(
                    uint32_t const *did_eb,
                    uint32_t const *did_ee,
                    ::ecal::reco::StorageScalarType const* amplitude_eb,   // in adc counts  
                    ::ecal::reco::StorageScalarType const* amplitude_ee,   // in adc counts  
                    ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
                    int const nchannels,
                    uint32_t const offsetForInput
         ) {
      
      
//       
//    NB: energy   "type_wrapper<reco::StorageScalarType, L>::type" most likely std::vector<float>
//       
      
      int ch = threadIdx.x + blockDim.x*blockIdx.x;
      
      int const inputCh = ch >= offsetForInput
                        ? ch - offsetForInput
                        : ch;
      
      uint32_t const * did = ch >= offsetForInput
                        ? did_ee
                        : did_eb;
      
      ::ecal::reco::StorageScalarType const* amplitude = ch >= offsetForInput
                        ? amplitude_ee
                        : amplitude_eb;
   
      
      if (ch < nchannels) {
        
        // simple copy
        energy[ch] = amplitude[inputCh];
      }
      
    }
    
    
    
    // host version, to be called by the plugin
    void create_ecal_rehit(
                  EventInputDataGPU const& eventInputGPU,
                  EventOutputDataGPU&      eventOutputGPU,
                  //     eventDataForScratchGPU_,
                  //     conditions,
                  //     configParameters_,
                  uint32_t const  offsetForInput,
                  cuda::stream_t<>& cudaStream
             ){
    
      int nchannels = 10; // FIXME
//       int offsetForInput = 5; // FIXME
      
      unsigned int totalChannels = 10; //eventInputGPU.ebUncalibRecHits.nchannels +
//       eventInputGPU.eeUncalibRecHits.nchannels;
      
      unsigned int nchannels_per_block = 32;
      unsigned int threads_1d = 10 * nchannels_per_block;
      //   unsigned int blocks_1d = threads_1d > 10*totalChannels  ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
      unsigned int blocks_1d = 2;
      
      // 
      // kernel
      //
      kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
        eventInputGPU.ebUncalibRecHits.did,
        eventInputGPU.eeUncalibRecHits.did,
        eventInputGPU.ebUncalibRecHits.amplitude, 
        eventInputGPU.eeUncalibRecHits.amplitude, 
        eventOutputGPU.energy,
        nchannels,
        offsetForInput
      );
      
      
//       /afs/cern.ch/work/a/amassiro/ECAL/GPU/onGPU/3July2019/CMSSW_10_6_0_Patatrack/src/RecoLocalCalo/EcalRecAlgos/src/EcalRecHitBuilderKernels.cu(66): error: no suitable conversion function from "const std::vector<uint32_t, CUDAHostAllocator<uint32_t, 0U>>" to "const uint32_t *" exists
      
      

      
    }
    
    

    
  }
  
}

