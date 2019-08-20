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
    void kernel_create_laser_corrections(
      // input
      float const* p1,
      float const* p2,
      float const* p3,
      edm::TimeValue_t const* t1,
      edm::TimeValue_t const* t2,
      edm::TimeValue_t const* t3,
      // output
      float *laser_corrections_out,
      // support
      int const nchannels
    ) {
      
      //
      // Calculate laser corrections:
      //
      
      int ch = threadIdx.x + blockDim.x*blockIdx.x;
      
      if (ch < nchannels) {
        
        laser_corrections_out[ch] = 1.;
                
        //       
        //       https://github.com/cms-sw/cmssw/blob/5d27b6509171a20b6e6a4bbbaf29ca471d612913/CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
        //       ./CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
        //       
        
        
        //       
        //       iLM = MEEBGeom::lmr(ebid.ieta(), ebid.iphi());
        //       iLM = MEEEGeom::lmr(iX, iY, eeid.zside());
        //             
        
        //         int MEEBGeom::lmr(EBGlobalCoord ieta, EBGlobalCoord iphi) {
        //            int idcc = dcc(ieta, iphi);
        //            int ism = idcc - 9;
        //            int iside = side(ieta, iphi);
        //            int ilmr = 1 + 2 * (ism - 1) + iside;
        //            return ilmr;
        //          }
//         
//         
//         int MEEBGeom::dcc(EBGlobalCoord ieta, EBGlobalCoord iphi) {
//             int ism = sm(ieta, iphi);
//             return dccFromSm(ism);
//           }
//           
//           int MEEBGeom::dccFromSm(int ism) {
//                assert(ism >= 1 && ism <= 36);
//                int iz = 1;
//                if (ism > 18)
//                  iz = -1;
//                if (iz == -1)
//                  ism -= 18;
//                assert(ism >= 1 && ism <= 18);
//                int idcc = 9 + ism;
//                if (iz == +1)
//                  idcc += 18;
//                return idcc;
//              }
//         
//         
        //       if (iLM - 1 < (int)laserTimeMap.size()) {
        //         timestamp = laserTimeMap[iLM - 1];
        //       } 
        //       
        //       timestamp = laserTimeMap[iLM - 1];
        //       
        //       
        //       void  setValue(uint32_t rawId, const EcalLaserAPDPNpair& value) { laser_map[rawId] = value; };
        //       const EcalLaserAPDPNRatiosMap& getLaserMap() const { return laser_map; }
        //       
        //       void setTime(int hashedIndex, const EcalLaserTimeStamp& value) { time_map[hashedIndex] = value; };
        //       const EcalLaserTimeStampMap& getTimeMap() const { return time_map; }  
        // 
        
        // interpolation
        
        //
        // FIXME: calculate iLM (i-laser monitoring region)
        //
        int iLM = 1;
       
//         auto const isBarrel = did_to_use.subdetId() == EcalBarrel;
//         auto const hashedId = isBarrel
// 
//         if (isBarrel) {
//           iLM = ...
//         }
//         else {
//           iLM = ...
//         }
//         
//         
        
        //       edm::TimeValue_t t = iTime.value();
        edm::TimeValue_t t = 0.;  // FIXME
        long long t_i = 0, t_f = 0;
        float p_i = 0, p_f = 0;
        long long lt_i = 0, lt_f = 0;
        float lp_i = 0, lp_f = 0;
        
        if (t >= t1[iLM - 1] && t < t2[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[ch];
          p_f = p2[ch];
        } else if (t >= t2[iLM - 1] && t <= t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[ch];
          p_f = p3[ch];
        } else if (t < t1[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[ch];
          p_f = p2[ch];
          
        } else if (t > t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[ch];
          p_f = p3[ch];
        }
        
        float apdpnref = 1.;
        float alpha = 1.;
        
        if (apdpnref != 0 && (t_i - t_f) != 0 && (lt_i - lt_f) != 0) {
          long long tt = t;  // never subtract two unsigned!
          float interpolatedLaserResponse =   p_i / apdpnref + float(tt - t_i)  * (p_f - p_i)   / (apdpnref * float(t_f - t_i));
          float interpolatedLinearResponse = lp_i / apdpnref + float(tt - lt_i) * (lp_f - lp_i) / (apdpnref * float(lt_f - lt_i));  // FIXED BY FC
          
          if (interpolatedLinearResponse > 2.f || interpolatedLinearResponse < 0.1f) {
            interpolatedLinearResponse = 1.f;
          }
          if (interpolatedLaserResponse <= 0.) {
            // how the heck is it possible?
            interpolatedLaserResponse = 0.0001;
          }
          
          float interpolatedTransparencyResponse = interpolatedLaserResponse / interpolatedLinearResponse;
          
          laser_corrections_out[ch] = 1.f / (std::pow(interpolatedTransparencyResponse, alpha) * interpolatedLinearResponse);
        }
       
      
      }
      
    }

    
    __global__
    void kernel_create_ecal_rehit(
                     // configuration 
                     int const* ChannelStatusToBeExcluded,
                     uint32_t ChannelStatusToBeExcludedSize,                     
                     // conditions
                     float const* adc2gev,
                     float const* intercalib,
                     uint32_t const* status,
                     float const* apdpnrefs,
                     float const* alphas,
                     // input for transparency corrections
                     float const* p1,
                     float const* p2,
                     float const* p3,
                     edm::TimeValue_t const* t1,
                     edm::TimeValue_t const* t2,
                     edm::TimeValue_t const* t3,  
                     // input for linear corrections
                     float const* lp1,
                     float const* lp2,
                     float const* lp3,
                     edm::TimeValue_t const* lt1,
                     edm::TimeValue_t const* lt2,
                     edm::TimeValue_t const* lt3,                    
                     // time, used for time dependent corrections
                     edm::TimeValue_t const event_time,
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
                            
  //                         
  // not used? ... yet ... it will be used to get IC, Laser Correction, ...
  //                         
        uint32_t const * didCh = ch >= offsetForInput
                          ? did_ee
                          : did_eb;
        
        // only two values, EB or EE
        // AM : FIXME : why not using "isBarrel" ?    isBarrel ? adc2gev[0] : adc2gev[1]
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
        // AM: move hashedIndexEB outside ecal::multifit   ---> change into reconstruction   [done]
        //
        
        auto const isBarrel = did_to_use.subdetId() == EcalBarrel;
        auto const hashedId = isBarrel
                            ? ecal::reconstruction::hashedIndexEB(did_to_use.rawId())
                            : offsetForHashes + ecal::reconstruction::hashedIndexEE(did_to_use.rawId());

        float const intercalib_to_use = intercalib[hashedId];
       
        
        // get laser coefficient
        float lasercalib = 1.;

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
        //       https://github.com/cms-sw/cmssw/blob/5d27b6509171a20b6e6a4bbbaf29ca471d612913/CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
        //       ./CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
        //       
        
        int iLM = 1;
        
        if (isBarrel) {
          iLM = ecal::reconstruction::laser_monitoring_region_EB (did_to_use.rawId());
        }
        else {
          iLM = ecal::reconstruction::laser_monitoring_region_EE (did_to_use.rawId());
          //           iLM = ...
        }
        
        
        long long t_i = 0, t_f = 0;
        float p_i = 0, p_f = 0;
        long long lt_i = 0, lt_f = 0;
        float lp_i = 0, lp_f = 0;
        
        // laser
        if (event_time >= t1[iLM - 1] && event_time < t2[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[ch];
          p_f = p2[ch];
        } else if (event_time >= t2[iLM - 1] && event_time <= t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[ch];
          p_f = p3[ch];
        } else if (event_time < t1[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[ch];
          p_f = p2[ch];
          
        } else if (event_time > t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[ch];
          p_f = p3[ch];
        }

        
        // linear corrections
        if (event_time >= lt1[iLM - 1] && event_time < lt2[iLM - 1]) {
          lt_i = lt1[iLM - 1];
          lt_f = lt2[iLM - 1];
          lp_i = lp1[ch];
          lp_f = lp2[ch];
        } else if (event_time >= lt2[iLM - 1] && event_time <= lt3[iLM - 1]) {
          lt_i = lt2[iLM - 1];
          lt_f = lt3[iLM - 1];
          lp_i = lp2[ch];
          lp_f = lp3[ch];
        } else if (event_time < lt1[iLM - 1]) {
          lt_i = lt1[iLM - 1];
          lt_f = lt2[iLM - 1];
          lp_i = lp1[ch];
          lp_f = lp2[ch];
          
        } else if (event_time > lt3[iLM - 1]) {
          lt_i = lt2[iLM - 1];
          lt_f = lt3[iLM - 1];
          lp_i = lp2[ch];
          lp_f = lp3[ch];
        }
        
        
        // apdpnref and alpha 
        float apdpnref = apdpnrefs[hashedId];
        float alpha = alphas[hashedId];
        
        // now calculate transparency correction
        if (apdpnref != 0 && (t_i - t_f) != 0 && (lt_i - lt_f) != 0) {
          long long tt = event_time;  // never subtract two unsigned!
          float interpolatedLaserResponse =   p_i / apdpnref + float(tt - t_i)  * (p_f - p_i)   / (apdpnref * float(t_f - t_i));
          float interpolatedLinearResponse = lp_i / apdpnref + float(tt - lt_i) * (lp_f - lp_i) / (apdpnref * float(lt_f - lt_i));  // FIXED BY FC
          
          if (interpolatedLinearResponse > 2.f || interpolatedLinearResponse < 0.1f) {
            interpolatedLinearResponse = 1.f;
          }
          if (interpolatedLaserResponse <= 0.) {
            // AM :  how the heck is it possible?
            interpolatedLaserResponse = 0.0001;
          }
          
          float interpolatedTransparencyResponse = interpolatedLaserResponse / interpolatedLinearResponse;
          
          // ... and now this:
          lasercalib = 1.f / (std::pow(interpolatedTransparencyResponse, alpha) * interpolatedLinearResponse);
        }
        
        
        
        
        //
        // Check for channels to be excluded from reconstruction
        //        
        //
        // Default energy? Not to be updated if "ChannelStatusToBeExcluded"
        // Exploited later by the module "EcalRecHitConvertGPU2CPUFormat"
        //
        energy[ch] = -1;
        
        static const int chStatusMask      = 0x1F;
        int dbstatus = (status[hashedId]) & chStatusMask;
        if (ChannelStatusToBeExcludedSize != 0) {
          for (int ich_to_check = 0; ich_to_check<ChannelStatusToBeExcludedSize; ich_to_check++) {
            if ( ChannelStatusToBeExcluded[ich_to_check] == dbstatus ) {
              return;
            }
          }
        }
        
        //
        // why are these channels killed in a different way?
        // and not via configuration, as the other bits?
        // is it like this from the past development????
        // ----> make it uniform!
        //
//         flagmask_ = 0;
//         flagmask_ |= 0x1 << EcalRecHit::kNeighboursRecovered;
//         flagmask_ |= 0x1 << EcalRecHit::kTowerRecovered;
//         flagmask_ |= 0x1 << EcalRecHit::kDead;
//         flagmask_ |= 0x1 << EcalRecHit::kKilled;
//         flagmask_ |= 0x1 << EcalRecHit::kTPSaturated;
//         flagmask_ |= 0x1 << EcalRecHit::kL1SpikeFlag;
        
        
        
        // multiply the adc counts with factors to get GeV
        
        energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib;
        
        
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
                  edm::TimeValue_t const event_time,
                  cuda::stream_t<>& cudaStream
             ){
    
      int nchannels = eventInputGPU.ebUncalibRecHits.size + eventInputGPU.eeUncalibRecHits.size ;
      
//       unsigned int totalChannels = 10; //eventInputGPU.ebUncalibRecHits.nchannels +
//       eventInputGPU.eeUncalibRecHits.nchannels;
      
      unsigned int nchannels_per_block = 32;
//       unsigned int threads_1d = 10 * nchannels_per_block;
      unsigned int threads_1d = nchannels_per_block;
      //   unsigned int blocks_1d = threads_1d > 10*totalChannels  ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
//       unsigned int blocks_1d = 2;
      unsigned int blocks_1d = (nchannels + threads_1d) / threads_1d; // TEST 
      
      // 
      // kernel update laser corrections
      //
      
//       float* device_laser_corrections;
//       cudaCheck( cudaMalloc((void**)&device_laser_corrections, 
//                             sizeof(float) * nchannels) 
//       );   
      
//       kernel_create_laser_corrections <<< blocks_1d, threads_1d >>> (
//         // input
//         conditions.LaserAPDPNRatios.p1,
//         conditions.LaserAPDPNRatios.p2,
//         conditions.LaserAPDPNRatios.p3,
//         conditions.LaserAPDPNRatios.t1,
//         conditions.LaserAPDPNRatios.t2,
//         conditions.LaserAPDPNRatios.t3,
//         // output
//         device_laser_corrections,
//         // support
//         nchannels       
//         );
        
//       
//       LaserAPDPNRatios   
//       LaserAPDPNRatiosRef
//       LaserAlphas        
//       LinearCorrections  
//       
      
      
      // 
      // kernel create rechit
      //
      kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
// configuration 
        configParameters.ChannelStatusToBeExcluded,
        configParameters.ChannelStatusToBeExcludedSize,
// conditions
        conditions.ADCToGeV.adc2gev,
        conditions.Intercalib.values,  
        conditions.ChannelStatus.status,  
        conditions.LaserAPDPNRatiosRef.values,  
        conditions.LaserAlphas.values,  
// input for transparency corrections
        conditions.LaserAPDPNRatios.p1,
        conditions.LaserAPDPNRatios.p2,
        conditions.LaserAPDPNRatios.p3,
        conditions.LaserAPDPNRatios.t1,
        conditions.LaserAPDPNRatios.t2,
        conditions.LaserAPDPNRatios.t3,
// input for linear corrections
        conditions.LinearCorrections.p1,
        conditions.LinearCorrections.p2,
        conditions.LinearCorrections.p3,
        conditions.LinearCorrections.t1,
        conditions.LinearCorrections.t2,
        conditions.LinearCorrections.t3,
// time, used for time dependent corrections
        event_time,
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

