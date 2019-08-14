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
                     uint32_t const offsetForInput
    ) {
      
      
//       
//    NB: energy   "type_wrapper<reco::StorageScalarType, L>::type" most likely std::vector<float>
//       
      
      int ch = threadIdx.x + blockDim.x*blockIdx.x;
      
      int const inputCh = ch >= offsetForInput
                        ? ch - offsetForInput
                        : ch;

//       int const inputCh = ch < offsetForInput
//                         ? ch
//                         : ch - offsetForInput;
                        
//                         
// not used? ... yet ... it will be used to get IC, Laser Correction, ...
//                         
      uint32_t const * didCh = ch >= offsetForInput
                        ? did_ee
                        : did_eb;
      

                        
      // first EB and then EE
                        
      ::ecal::reco::StorageScalarType const* amplitude = ch >= offsetForInput
                        ? amplitude_ee
                        : amplitude_eb;
   
      ::ecal::reco::StorageScalarType const* time_in = ch >= offsetForInput
                        ? time_ee
                        : time_eb;

      ::ecal::reco::StorageScalarType const* chi2_in = ch >= offsetForInput
                        ? chi2_ee
                        : chi2_eb;
      
      if (ch < nchannels) {
        
        // simple copy
        did[ch] = didCh[inputCh];

        energy[ch] = amplitude[inputCh];
        
        // FIXME
        // 
        //  From: https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerSimple.cc 
        //
        // - get ADCToGeVConstant
        // - get IC
        // - get Laser Correction
        //
        // - correct from the "jitter" to "time" properly
        //
        // - what is "extra" ?
        //
        
        // Time is not saved so far, FIXME
//         time[ch] = time_in[inputCh];
        chi2[ch] = chi2_in[inputCh];
        
        // FIXME: calculate the flagBits
        flagBits[ch] = 0;
        extra[ch] = 0;
        
      }
      
    }
    
    
    
    // host version, to be called by the plugin
    void create_ecal_rehit(
                  EventInputDataGPU const& eventInputGPU,
                  EventOutputDataGPU&      eventOutputGPU,
                  //     eventDataForScratchGPU_,
                  ConditionsProducts const& conditions, 
                  //     configParameters_,
                  uint32_t const  offsetForInput,
                  cuda::stream_t<>& cudaStream
             ){
    
      int nchannels = eventInputGPU.ebUncalibRecHits.size + eventInputGPU.eeUncalibRecHits.size ;
      
//       unsigned int totalChannels = 10; //eventInputGPU.ebUncalibRecHits.nchannels +
//       eventInputGPU.eeUncalibRecHits.nchannels;
      
      unsigned int nchannels_per_block = 32;
//       unsigned int threads_1d = 10 * nchannels_per_block;
      unsigned int threads_1d = nchannels_per_block;
      //   unsigned int blocks_1d = threads_1d > 10*totalChannels  ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
      unsigned int blocks_1d = 2;
      
      // 
      // kernel
      //
      kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
// input
        eventInputGPU.ebUncalibRecHits.did,
        eventInputGPU.eeUncalibRecHits.did,
        eventInputGPU.ebUncalibRecHits.amplitude, 
        eventInputGPU.eeUncalibRecHits.amplitude, 
        eventInputGPU.ebUncalibRecHits.jitter, 
        eventInputGPU.eeUncalibRecHits.jitter, 
        eventInputGPU.ebUncalibRecHits.chi2, 
        eventInputGPU.eeUncalibRecHits.chi2, 
// output
        eventOutputGPU.did,
        eventOutputGPU.energy,
        eventOutputGPU.time,
        eventOutputGPU.chi2,
        eventOutputGPU.flagBits,
        eventOutputGPU.extra,
// other
        nchannels,
        offsetForInput
      );
      
      
//       /afs/cern.ch/work/a/amassiro/ECAL/GPU/onGPU/3July2019/CMSSW_10_6_0_Patatrack/src/RecoLocalCalo/EcalRecAlgos/src/EcalRecHitBuilderKernels.cu(66): error: no suitable conversion function from "const std::vector<uint32_t, CUDAHostAllocator<uint32_t, 0U>>" to "const uint32_t *" exists
      
      

      
    }
    
    

    
  }
  
}

