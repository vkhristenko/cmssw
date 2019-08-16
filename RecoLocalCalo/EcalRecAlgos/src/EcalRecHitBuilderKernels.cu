#include "cuda.h"

#include "KernelHelpers.h"

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

//
//
#include "EcalRecHitBuilderKernels.h"


#include "KernelHelpers.h"


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
    ) {
      
      
//       
//    NB: energy   "type_wrapper<reco::StorageScalarType, L>::type" most likely std::vector<float>
//       
      
      int ch = threadIdx.x + blockDim.x*blockIdx.x;

      if (ch < nchannels) {
        
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
        
        // only two values, EB or EE
        float adc2gev_to_use = ch >= offsetForInput
                          ? adc2gev[0]
                          : adc2gev[1];
                          
  
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
        
  
        // simple copy
        did[ch] = didCh[inputCh];
        
        auto const did_to_use = DetId{didCh[inputCh]};
        
        //
        // AM FIXME: move hashedIndexEB outside ecal::multifit
        //
        
        auto const isBarrel = did_to_use.subdetId() == EcalBarrel;
        auto const hashedId = isBarrel
                            ? ecal::multifit::hashedIndexEB(did_to_use.rawId())
                            : offsetForHashes + ecal::multifit::hashedIndexEE(did_to_use.rawId());

        float const intercalib_to_use = intercalib[hashedId];
       
        
        // get laser coefficient
        float lasercalib = 1.;
//         lasercalib = laser->getLaserCorrection(detid, evt.time());
        //
        // AM: ideas
        //
        //    One possibility is to create the map of laser corrections once on CPU
        //    for all crystals and push them on GPU.
        //    Then only if the LS is different, update the laser correction
        //    The variation within a LS is not worth pursuing (<< 0.1% !!)
        //    and below the precision we can claim on the laser corrections (right?).
        //    This will save quite some time (also for the CPU version?)    
        //
        
        
        
        //
        // check for channels to be excluded from reconstruction
        //
//         EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
//         EcalChannelStatusCode::Code dbstatus = chit->getStatusCode();  ------>     Code  getStatusCode() const { return Code(status_&chStatusMask); }
//          static const int chStatusMask      = 0x1F;     
//         // check for channels to be excluded from reconstruction
//         if (!v_chstatus_.empty()) {
//           std::vector<int>::const_iterator res = std::find(v_chstatus_.begin(), v_chstatus_.end(), dbstatus);
//           if (res != v_chstatus_.end())
//             return false;
//         }
        
        static const int chStatusMask      = 0x1F;
        int dbstatus = (status[hashedId]) & chStatusMask;
        if (ChannelStatusToBeExcludedSize != 0) {
          for (int ich_to_check = 0; ich_to_check<ChannelStatusToBeExcludedSize; ich_to_check++) {
            if ( ChannelStatusToBeExcluded[ich_to_check] == dbstatus ) {
              return;
            }
          }
        }
        
        
        
        
        // multiply the adc counts with factors to get GeV
        
        energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib;
        
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
                  ConfigurationParameters const& configParameters,
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
// configuration 
        configParameters.ChannelStatusToBeExcluded,
        configParameters.ChannelStatusToBeExcludedSize,
// conditions
        conditions.ADCToGeV.adc2gev,
        conditions.Intercalib.values,  
        conditions.ChannelStatus.status,  
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
        offsetForInput,
        conditions.offsetForHashes
      );
      
      
//       /afs/cern.ch/work/a/amassiro/ECAL/GPU/onGPU/3July2019/CMSSW_10_6_0_Patatrack/src/RecoLocalCalo/EcalRecAlgos/src/EcalRecHitBuilderKernels.cu(66): error: no suitable conversion function from "const std::vector<uint32_t, CUDAHostAllocator<uint32_t, 0U>>" to "const uint32_t *" exists
      
      

      
    }
    
    

    
  }
  
}

