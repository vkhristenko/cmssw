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
                     bool killDeadChannels,
                     // for flags setting
                     uint32_t const* expanded_v_DB_reco_flags,
                     uint32_t const* expanded_Sizes_v_DB_reco_flags,
                     uint32_t const* expanded_flagbit_v_DB_reco_flags,
                     uint32_t expanded_v_DB_reco_flagsSize,
                     uint32_t flagmask,
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
        
        
        
        // Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
        int iterator_flags = 0;
        bool need_to_exit = false;
        for (unsigned int i = 0; i != expanded_v_DB_reco_flagsSize; ++i) {
          
          for (unsigned int j = 0; j != expanded_Sizes_v_DB_reco_flags[i]; j++) {
            if ( expanded_v_DB_reco_flags[iterator_flags] == dbstatus ) {
              flagBits[ch] =  0x1 << expanded_flagbit_v_DB_reco_flags[i];
              need_to_exit = true;
              break; // from the big loop!!!
            }
            if (need_to_exit) {
              break;
            }
            iterator_flags++;
          }
        }
                    
//         for (unsigned int i = 0; i != v_DB_reco_flags_.size(); ++i) {
//           if (std::find(v_DB_reco_flags_[i].begin(), v_DB_reco_flags_[i].end(), dbstatus) != v_DB_reco_flags_[i].end()) {
//             flagBits[ch] =  0x1 << i;
//             break;
//           }
//         }
        
        
        // AM: FIXME
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
        //       
        //         if (!(flagmask_ & flagBits) || !killDeadChannels_) {
        
        //         bool killDeadChannels = true;
        
        if ( (flagmask & flagBits[ch]) && killDeadChannels ) {
          return;
        }
        
        
        
        
        
        // multiply the adc counts with factors to get GeV
        
        energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib;
        
        
        // Time is not saved so far, FIXME
//         time[ch] = time_in[inputCh];
        chi2[ch] = chi2_in[inputCh];
        
        // FIXME: calculate the flagBits extra
        extra[ch] = 0;
        
        //
        // extra packing ...
        //
        
        uint32_t offset;
        uint32_t width;
        uint32_t value;
        
        float chi2_temp = chi2[ch];
        if (chi2_temp > 64) chi2_temp = 64;
        // use 7 bits
        uint32_t rawChi2 = lround(chi2_temp / 64. * ((1<<7)-1));
        
//         extra_ = setMasked(extra_, rawChi2, 0, 7);
        
        offset = 0;
        width = 7;
        value = 0; // default: https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L65
        
        uint32_t mask = ((1 << width) - 1) << offset;
        value &= ~mask;
        value |= (rawChi2 & ((1U << width) - 1)) << offset;
        
//         extra[ch] = value;
//         
//         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L126-L133
//         
        
        
        uint32_t rawEnergy = 0;
        
        if (energy[ch] > 0.001) {
//           uint16_t exponent = getPower10(energy[ch]);

          static constexpr float p10[] = {1.e-2f,1.e-1f,1.f,1.e1f,1.e2f,1.e3f,1.e4f,1.e5f,1.e6f};
          int b = energy[ch]<p10[4] ? 0 : 5;
          for (;b<9;++b) if (energy[ch]<p10[b]) break;
          
          uint16_t exponent = b;
          
          static constexpr float ip10[] = {1.e5f,1.e4f,1.e3f,1.e2f,1.e1f,1.e0f,1.e-1f,1.e-2f,1.e-3f,1.e-4};
          uint16_t significand = lround( energy[ch] * ip10[exponent]);
          // use 13 bits (3 exponent, 10 significand)
          rawEnergy = exponent << 10 | significand;
          /* here for doc and regression test
           *             u int16_t exponent_old = lround(floor(log10(*energy))) + 3;  
           *             uint16_t significand_old = lround(energy/pow(10, exponent - 5));
           *             std::cout << energy << ' ' << exponent << ' ' << significand 
           *             << ' ' << exponent_old <<     ' ' << significand_old << std::endl;
           *             assert(exponent==exponent_old);
           *             assert(significand==significand_old);
           */
        }
        
//         extra_ = setMasked(extra_, rawEnergy, 8, 13);
        
        offset = 8;
        width = 13;
        // value from last change, ok
        
        mask = ((1 << width) - 1) << offset;
        value &= ~mask;
        value |= (rawEnergy & ((1U << width) - 1)) << offset;
        
        extra[ch] = value;
          
        
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
            
      unsigned int nchannels_per_block = 32;
      unsigned int threads_1d = nchannels_per_block;
      unsigned int blocks_1d = (nchannels + threads_1d) / threads_1d; // TEST : to be optimized (AM)
      
      
      // 
      // kernel create rechit
      //
      kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
// configuration 
        configParameters.ChannelStatusToBeExcluded,
        configParameters.ChannelStatusToBeExcludedSize,
        configParameters.killDeadChannels,
// for flags setting
        configParameters.expanded_v_DB_reco_flags,
        configParameters.expanded_Sizes_v_DB_reco_flags,
        configParameters.expanded_flagbit_v_DB_reco_flags,
        configParameters.expanded_v_DB_reco_flagsSize,
        configParameters.flagmask,
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
      
      
      
    }
    
    

    
  }
  
}

