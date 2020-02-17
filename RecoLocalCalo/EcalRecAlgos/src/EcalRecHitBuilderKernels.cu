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
    
    
    // uncalibrecHit flags
    enum UncalibRecHitFlags {
      kGood=-1,                 // channel is good (mutually exclusive with other states)  setFlagBit(kGood) reset flags_ to zero 
      kPoorReco,                // channel has been badly reconstructed (e.g. bad shape, bad chi2 etc.)
      kSaturated,               // saturated channel
      kOutOfTime,               // channel out of time
      kLeadingEdgeRecovered,    // saturated channel: energy estimated from the leading edge before saturation
      kHasSwitchToGain6,        // at least one data frame is in G6
      kHasSwitchToGain1         // at least one data frame is in G1
    };
    
    
    // recHit flags
    enum RecHitFlags { 
      RecHitFlags_kGood=0,                   // channel ok, the energy and time measurement are reliable
      RecHitFlags_kPoorReco,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
      RecHitFlags_kOutOfTime,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
      RecHitFlags_kFaultyHardware,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
      RecHitFlags_kNoisy,                    // the channel is very noisy
      RecHitFlags_kPoorCalib,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
      RecHitFlags_kSaturated,                // saturated channel (recovery not tried)
      RecHitFlags_kLeadingEdgeRecovered,     // saturated channel: energy estimated from the leading edge before saturation
      RecHitFlags_kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
      RecHitFlags_kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
      RecHitFlags_kDead,                     // channel is dead and any recovery fails
      RecHitFlags_kKilled,                   // MC only flag: the channel is killed in the real detector
      RecHitFlags_kTPSaturated,              // the channel is in a region with saturated TP
      RecHitFlags_kL1SpikeFlag,              // the channel is in a region with TP with sFGVB = 0
      RecHitFlags_kWeird,                    // the signal is believed to originate from an anomalous deposit (spike) 
      RecHitFlags_kDiWeird,                  // the signal is anomalous, and neighbors another anomalous signal  
      RecHitFlags_kHasSwitchToGain6,         // at least one data frame is in G6
      RecHitFlags_kHasSwitchToGain1,         // at least one data frame is in G1
      //
      RecHitFlags_kUnknown                   // to ease the interface with functions returning flags. 
    };
    
    
    // status code
    enum EcalChannelStatusCode_Code {
      kOk=0,
      kDAC,
      kNoLaser,
      kNoisy,
      kNNoisy,
      kNNNoisy,
      kNNNNoisy,
      kNNNNNoisy,
      kFixedG6,
      kFixedG1,
      kFixedG0,
      kNonRespondingIsolated,
      kDeadVFE,
      kDeadFE,
      kNoDataNoTP      
    };
    
    
    
    
    
    __global__
    void kernel_create_ecal_rehit(
      // configuration 
      int const* ChannelStatusToBeExcluded,
      uint32_t ChannelStatusToBeExcludedSize,   
      bool const killDeadChannels,
      bool const recoverEBIsolatedChannels,
      bool const recoverEEIsolatedChannels,
      bool const recoverEBVFE,             
      bool const recoverEEVFE,             
      bool const recoverEBFE,             
      bool const recoverEEFE,              
      float const EBLaserMIN,
      float const EELaserMIN,
      float const EBLaserMAX,
      float const EELaserMAX,
      // for flags setting
      int const* expanded_v_DB_reco_flags,    // FIXME AM: to be checked
      uint32_t const* expanded_Sizes_v_DB_reco_flags,
      uint32_t const* expanded_flagbit_v_DB_reco_flags,
      uint32_t expanded_v_DB_reco_flagsSize,
      uint32_t flagmask,
      // conditions
      float const* adc2gev,
      float const* intercalib,
      uint16_t const* status,
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
      uint32_t const* flags_eb,
      uint32_t const* flags_ee,
      // output
      uint32_t *did,
      ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
      ::ecal::reco::StorageScalarType* time,  
      ::ecal::reco::StorageScalarType* chi2,  
      uint32_t* flagBits,
      uint32_t* extra,
      // other
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
        ? adc2gev[1]  // ee
        : adc2gev[0]; // eb
        
        
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
        
        uint32_t const* flags_in = ch >= offsetForInput
        ? flags_ee
        : flags_eb;
        
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
          p_i = p1[hashedId];
          p_f = p2[hashedId];
        } else if (event_time >= t2[iLM - 1] && event_time <= t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[hashedId];
          p_f = p3[hashedId];
        } else if (event_time < t1[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[hashedId];
          p_f = p2[hashedId];
          
        } else if (event_time > t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[hashedId];
          p_f = p3[hashedId];
        }
        
        
        // linear corrections
        if (event_time >= lt1[iLM - 1] && event_time < lt2[iLM - 1]) {
          lt_i = lt1[iLM - 1];
          lt_f = lt2[iLM - 1];
          lp_i = lp1[hashedId];
          lp_f = lp2[hashedId];
        } else if (event_time >= lt2[iLM - 1] && event_time <= lt3[iLM - 1]) {
          lt_i = lt2[iLM - 1];
          lt_f = lt3[iLM - 1];
          lp_i = lp2[hashedId];
          lp_f = lp3[hashedId];
        } else if (event_time < lt1[iLM - 1]) {
          lt_i = lt1[iLM - 1];
          lt_f = lt2[iLM - 1];
          lp_i = lp1[hashedId];
          lp_f = lp2[hashedId];
          
        } else if (event_time > lt3[iLM - 1]) {
          lt_i = lt2[iLM - 1];
          lt_f = lt3[iLM - 1];
          lp_i = lp2[hashedId];
          lp_f = lp3[hashedId];
        }
        
        
        // apdpnref and alpha 
        float apdpnref = apdpnrefs[hashedId];
        float alpha = alphas[hashedId];
        
        // now calculate transparency correction
        if (apdpnref != 0 && (t_i - t_f) != 0 && (lt_i - lt_f) != 0) {
          long long tt = event_time;  // never subtract two unsigned!
          float interpolatedLaserResponse =   p_i / apdpnref + float(tt - t_i)  * (p_f - p_i)   / (apdpnref * float(t_f - t_i));
          //           float interpolatedLaserResponse =   iLM+10; p_f*=2; p_i*=2;
          //           float interpolatedLaserResponse =   iLM+float(t_f - t_i)+float(tt - t_i); p_f*=2; p_i*=2;
          //           float interpolatedLaserResponse =   iLM+float(t_f - t_i)+float(tt - t_i)+apdpnref+alpha+100; p_f*=2; p_i*=2;
          // float interpolatedLaserResponse =   iLM+float(t_f - t_i)+float(tt - t_i)+apdpnref+alpha+100+p_i; p_f*=2; p_i*=2;
          // float interpolatedLaserResponse =   p_i+apdpnref; p_f*=2; p_i*=2;
          // float interpolatedLaserResponse =   apdpnref; p_f*=2; p_i*=2; // good
          // float interpolatedLaserResponse =   p_i; p_f*=2; p_i*=2; // bad
          // float interpolatedLaserResponse =   iLM; p_f*=2; p_i*=2; // bad
          // float interpolatedLaserResponse =   did_to_use.rawId()/1000000000.; p_f*=2; p_i*=2; // good
          // float interpolatedLaserResponse =   iLM; p_f*=2; p_i*=2; // bad ??  ---> now good??????
          
          // float interpolatedLaserResponse =   tt/6589055619675208704.; p_f*=2; p_i*=2; // good ... but not 100% sure
          // float interpolatedLaserResponse =   1.+float(tt-t_i)*(p_f-p_i)/(apdpnref*float(t_f-t_i)); p_f*=2; p_i*=2; // bad ?
          
          
          float interpolatedLinearResponse = lp_i / apdpnref + float(tt - lt_i) * (lp_f - lp_i) / (apdpnref * float(lt_f - lt_i));  // FIXED BY FC
          
          if (interpolatedLinearResponse > 2.f || interpolatedLinearResponse < 0.1f) {
            interpolatedLinearResponse = 1.f;
          }
          if (interpolatedLaserResponse <= 0.) {
            // AM :  how the heck is it possible?
            //             interpolatedLaserResponse = 0.0001;
            lasercalib = 1.;
            
          }
          else {
            
            float interpolatedTransparencyResponse = interpolatedLaserResponse / interpolatedLinearResponse;
            
            // ... and now this:
            lasercalib = 1.f / ( std::pow(interpolatedTransparencyResponse, alpha) * interpolatedLinearResponse);
            
            //             alpha *= 2;
            //             
            //             lasercalib = interpolatedTransparencyResponse;
            //             
            //             lasercalib = 10.;
            //             
            //             lasercalib = interpolatedLaserResponse;
            
            //             lasercalib = interpolatedLinearResponse; // this is ok
            
          }
        }
        
        
        
        //
        // Check for channels to be excluded from reconstruction
        //        
        //
        // Default energy? Not to be updated if "ChannelStatusToBeExcluded"
        // Exploited later by the module "EcalRecHitConvertGPU2CPUFormat"
        //
        energy[ch] = -1; //---- -1 or -100 FIXME AM
        
        //
        // From EcalRecHitWorkerSimple.cc
        //
        //             std::vector<int> v_chstatus_;
        //
        //     EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
        //     EcalChannelStatusCode::Code  dbstatus = chit->getStatusCode();
        //     
        //     // check for channels to be excluded from reconstruction
        //     if ( !v_chstatus_.empty()) {
        //       
        //       std::vector<int>::const_iterator res = 
        //       std::find( v_chstatus_.begin(), v_chstatus_.end(), dbstatus );
        //       if ( res != v_chstatus_.end() ) return false;
        //       
        //     }
        //
        //         /// return decoded status
        //             Code  getStatusCode() const { return Code(status_&chStatusMask); }
        //         
        //              static const int chStatusMask      = 0x1F;
        //
        //
        static const int chStatusMask      = 0x1F;
        // ChannelStatusToBeExcluded is a "int" then I put "dbstatus" to be the same
        int dbstatus = EcalChannelStatusCode_Code( (status[hashedId]) & chStatusMask );
        //         printf (" ChannelStatusToBeExcludedSize = %d \n", ChannelStatusToBeExcludedSize) ; 
        if (ChannelStatusToBeExcludedSize != 0) {
          for (int ich_to_check = 0; ich_to_check<ChannelStatusToBeExcludedSize; ich_to_check++) {
            //             printf (" ChannelStatusToBeExcluded[%d] = %d  --> dbstatus = %d \n", ich_to_check, ChannelStatusToBeExcluded[ich_to_check], dbstatus) ; 
            if ( ChannelStatusToBeExcluded[ich_to_check] == dbstatus ) {
              //               printf (" killed @ GPU %d -->  ChannelStatusToBeExcluded[%d] = %d  --> dbstatus = %d \n", hashedId, ich_to_check, ChannelStatusToBeExcluded[ich_to_check], dbstatus) ; 
              return; 
            }
          }
        }
        
        
        // Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
        
        //
        // AM: get the smaller "flagbit_counter" with match
        //
        
        uint32_t temporary_flagBits = 0;
        
        int iterator_flags = 0;
        bool need_to_exit = false;
        int flagbit_counter = 0;
        while (!need_to_exit) {
          iterator_flags = 0;
          for (unsigned int i = 0; i != expanded_v_DB_reco_flagsSize; ++i) { 
            // check the correct "flagbit"
            if (expanded_flagbit_v_DB_reco_flags[i] == flagbit_counter) {
              
              for (unsigned int j = 0; j < expanded_Sizes_v_DB_reco_flags[i]; j++) {
                
                if ( expanded_v_DB_reco_flags[iterator_flags] == dbstatus ) {
                  temporary_flagBits =  0x1 << expanded_flagbit_v_DB_reco_flags[i];      
                  need_to_exit = true;
                  break; // also from the big loop!!!
                  
                }
                iterator_flags++;
              }
            }
            else {
              // if not, got to the next bunch directly
              iterator_flags += expanded_Sizes_v_DB_reco_flags[i];
            }
            
            if (need_to_exit) {
              break;
            }
            
          }
          flagbit_counter+=1;
        }
        
        
        //         flagBits[ch] =  dbstatus; // AM: FIXME : just a test
        
        
        //         for (unsigned int i = 0; i != v_DB_reco_flags_.size(); ++i) {
        //           if (std::find(v_DB_reco_flags_[i].begin(), v_DB_reco_flags_[i].end(), dbstatus) != v_DB_reco_flags_[i].end()) {
        //             flagBits[ch] =  0x1 << i;
        //             break;
        //           }
        //         }
        
        
        
        
        //         class EcalRecHit {
        //         public:
        //           typedef DetId key_type;
        //           
        //           // recHit flags
        //           enum Flags { 
        //             kGood=0,                   // channel ok, the energy and time measurement are reliable
        //             kPoorReco,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
        //             kOutOfTime,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
        //             kFaultyHardware,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
        //             kNoisy,                    // the channel is very noisy
        //             kPoorCalib,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
        //             kSaturated,                // saturated channel (recovery not tried)
        //             kLeadingEdgeRecovered,     // saturated channel: energy estimated from the leading edge before saturation
        //             kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
        //             kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
        //             kDead,                     // channel is dead and any recovery fails
        //             kKilled,                   // MC only flag: the channel is killed in the real detector
        //             kTPSaturated,              // the channel is in a region with saturated TP
        //             kL1SpikeFlag,              // the channel is in a region with TP with sFGVB = 0
        //             kWeird,                    // the signal is believed to originate from an anomalous deposit (spike) 
        //             kDiWeird,                  // the signal is anomalous, and neighbors another anomalous signal  
        //             kHasSwitchToGain6,         // at least one data frame is in G6
        //             kHasSwitchToGain1,         // at least one data frame is in G1
        //             //
        //             kUnknown                   // to ease the interface with functions returning flags. 
        //           };
        //           
        //           
        //           
        //           // Associate reco flagbit ( outer vector) to many db status flags (inner vector)
        //           std::vector<std::vector<uint32_t> > v_DB_reco_flags_;
        //           
        //           // Traslate string representation of flagsMapDBReco into enum values 
        //           const edm::ParameterSet & p=ps.getParameter< edm::ParameterSet >("flagsMapDBReco");
        //           std::vector<std::string> recoflagbitsStrings = p.getParameterNames();
        //           v_DB_reco_flags_.resize(32); 
        //         
        //         for (unsigned int i=0;i!=recoflagbitsStrings.size();++i){
        //           EcalRecHit::Flags recoflagbit = (EcalRecHit::Flags)  StringToEnumValue<EcalRecHit::Flags>(recoflagbitsStrings[i]);
        //           std::vector<std::string> dbstatus_s =  p.getParameter<std::vector<std::string> >(recoflagbitsStrings[i]);
        //           std::vector<uint32_t> dbstatuses;
        //           for (unsigned int j=0; j!= dbstatus_s.size(); ++j){
        //             EcalChannelStatusCode::Code  dbstatus  = (EcalChannelStatusCode::Code)
        //             StringToEnumValue<EcalChannelStatusCode::Code>(dbstatus_s[j]);
        //             dbstatuses.push_back(dbstatus);
        //           }
        //           
        //           v_DB_reco_flags_[recoflagbit]=dbstatuses;
        //         }
        //         
        //         uint32_t flagBits = setFlagBits(v_DB_reco_flags_, dbstatus);
        //         
        //         // Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
        //         uint32_t EcalRecHitWorkerSimple::setFlagBits(const std::vector<std::vector<uint32_t> >& map, 
        //                                                      const uint32_t& status  ){
        //           
        //           for (unsigned int i = 0; i!=map.size(); ++i){
        //             if (std::find(map[i].begin(), map[i].end(),status)!= map[i].end()) 
        //               return 0x1 << i;
        //           }
        //           
        //           return 0;
        //                                                      }
        
        
        
        
        
        
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
        //
        //         ----> from EcalRecHitWorkerSimple.cc
        //
        //         -----> I have filled "flagmask" in the EcalRecHitProducerGPU.cc
        //
        //         bool killDeadChannels = true;
        
        //
        // defined here: https://github.com/amassiro/cmssw/blob/amassiro_patatrack_ECALrechit_on_GPU/RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitProducerGPU.cc#L265-L273
        // 
        //          
        //     In EcalRecHitWorkerSimple.cc as:
        // 
        //         if (! (flagmask_ & flagBits ) || !killDeadChannels_) {
        //           
        
        if ( (flagmask & temporary_flagBits) && killDeadChannels ) {
          return;
        }
        
        
        //
        // AM: Seriously, was it useless to set it before???
        //         flagBits[ch] = 0;
        // FIXME
        // AM : no no ... now we store in flagBits[ch] what we had before in temporary_flagBits 
        flagBits[ch] = temporary_flagBits;
        // we don't leave it to 0 !!!
        
        
        
        
        
        
        //
        // multiply the adc counts with factors to get GeV
        //
        
        //         energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use ;
        energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib;
        // AM: lasercalib *was* the wrong one!  ---> now fixed
        
        // Time is not saved so far, FIXME
        //         time[ch] = time_in[inputCh];
        
        
        //      setChi2
        //         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L126-L133
        //                 
        //         
        //         // bound the max value of the chi2
        //         if (chi2 > 64) chi2 = 64;
        //         
        // -> introduce bound also here in order to get cpu = gpu
        //
        //         
        //         chi2[ch] = chi2_in[inputCh];
        //         
        if (chi2_in[inputCh] > 64) chi2[ch] = 64;
        else chi2[ch] = chi2_in[inputCh];
        
        
        
        
        // FIXME: calculate the "flagBits extra"  --> not really "flags", but actually an encoded version of energy uncertainty, time unc., ...
        extra[ch] = 0;  // -1 or 0 as default? ---> it is an unsgned int, it cannot be -1 !
        
        //
        // extra packing ...
        //
        
        //      setChi2
        //         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L126-L133
        //                 
        //         
        //         // bound the max value of the chi2
        //         if (chi2 > 64) chi2 = 64;
        //         
        //         // use 7 bits
        //         uint32_t rawChi2 = lround(chi2 / 64. * ((1<<7)-1));
        //         extra_ = setMasked(extra_, rawChi2, 0, 7);
        //         
        //         
        
        //         
        //         static inline uint32_t setMasked(uint32_t value, uint32_t x, uint32_t offset, uint32_t width) {
        //           const uint32_t mask = ((1 << width) - 1) << offset;
        //           value &= ~mask;
        //           value |= (x & ((1U << width) - 1)) << offset;
        //           return value;
        //         }
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
        
        
        //      setEnergyError
        //         rh.setEnergyError( uncalibRH.amplitudeError()*adcToGeVConstant_*intercalibConstant);
        //         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L141-L163
        //         
        
        // rawEnergy is actually "error" !!!
        uint32_t rawEnergy = 0;
        
        
        // AM: FIXME: this is not propagated currently to the uncalibrecHit collection SOA 
        //            if you want to store this in "extra", we need first to add it to the uncalibrecHit results
        //            then it will be something like the following
        //         amplitudeError[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib
        //         
        //         
        
        float amplitudeError_ch = 0. ; // amplitudeError[ch];
        
        if (amplitudeError_ch > 0.001) {
          //           uint16_t exponent = getPower10(amplitudeError_ch);
          
          static constexpr float p10[] = {1.e-2f,1.e-1f,1.f,1.e1f,1.e2f,1.e3f,1.e4f,1.e5f,1.e6f};
          int b = amplitudeError_ch<p10[4] ? 0 : 5;
          for (;b<9;++b) if (amplitudeError_ch<p10[b]) break;
          
          uint16_t exponent = b;
          
          static constexpr float ip10[] = {1.e5f,1.e4f,1.e3f,1.e2f,1.e1f,1.e0f,1.e-1f,1.e-2f,1.e-3f,1.e-4};
          uint16_t significand = lround( amplitudeError_ch * ip10[exponent]);
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
        
        
        //         // set the energy uncertainty
        //         // (only energy >= 0 will be stored)
        //         void setEnergyError(float energy) {
        //           uint32_t rawEnergy = 0;
        //           if (energy > 0.001) {
        //             uint16_t exponent = getPower10(energy);
        //             static constexpr float ip10[] = {1.e5f,1.e4f,1.e3f,1.e2f,1.e1f,1.e0f,1.e-1f,1.e-2f,1.e-3f,1.e-4};
        //             uint16_t significand = lround(energy*ip10[exponent]);
        //             // use 13 bits (3 exponent, 10 significand)
        //             rawEnergy = exponent << 10 | significand;
        //           }
        //           
        //           extra_ = setMasked(extra_, rawEnergy, 8, 13);
        //         }
        
        
        //         
        //         void setTimeError(uint8_t timeErrBits) {
        //           extra_ = setMasked(extra_, timeErrBits & 0xFF, 24, 8);
        //         }
        //         
        //         uint8_t uncalibRH.jitterErrorBits() ---> set to 0 if not available? FIXME: check default value
        //        
        //         EcalRecHit rh( uncalibRH.id(), energy, time );
        //         rh.setChi2( uncalibRH.chi2() );
        //         rh.setEnergyError( uncalibRH.amplitudeError()*adcToGeVConstant_*intercalibConstant);
        //         rh.setTimeError(uncalibRH.jitterErrorBits());
        //         
        //         
        
        uint32_t jitterErrorBits = 0;
        jitterErrorBits = jitterErrorBits & 0xFF;
        
        
        offset = 24;
        width = 8;
        // value from last change, ok
        
        mask = ((1 << width) - 1) << offset;
        value &= ~mask;
        value |= (jitterErrorBits & ((1U << width) - 1)) << offset;
        
        //
        // now finally set "extra[ch]"
        //
        extra[ch] = value;
        
        
        
        //
        // additional flags setting
        //
        // using correctly the flags as calculated at the UncalibRecHit stage
        //
        // Now fill flags
        
        
        
        
        bool good = true;
        
        //
        //         from: EcalUncalibratedRecHit.cc
        //         
        //         
        //         bool EcalUncalibratedRecHit::checkFlag(EcalUncalibratedRecHit::Flags flag) const {
        //           if(flag == kGood){ if ( ! flags_ ) return true;else return false;} // if all flags are unset, then hit is good
        //           return  flags_ & ( 0x1<<flag);
        //         }
        //         
        //         
        //         
        //         void EcalUncalibratedRecHit::setFlagBit(EcalUncalibratedRecHit::Flags flag){
        //           if  (flag == kGood) {
        //             //then set all bits to zero;
        //             flags_  = 0;
        //             return;
        //           }
        //           // else set the flagbit
        //           flags_|= 0x1 <<  flag;  
        //         }
        //         
        //         
        //  
        //         enum Flags {
        //           kGood=-1,                 // channel is good (mutually exclusive with other states)  setFlagBit(kGood) reset flags_ to zero 
        //           kPoorReco,                // channel has been badly reconstructed (e.g. bad shape, bad chi2 etc.)
        //           kSaturated,               // saturated channel
        //           kOutOfTime,               // channel out of time
        //           kLeadingEdgeRecovered,    // saturated channel: energy estimated from the leading edge before saturation
        //           kHasSwitchToGain6,        // at least one data frame is in G6
        //           kHasSwitchToGain1         // at least one data frame is in G1
        //           
        //         };
        //         
        
        
        if ( flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kLeadingEdgeRecovered) ) ) {
          flagBits[ch]  |=  (0x1 <<  (RecHitFlags::RecHitFlags_kLeadingEdgeRecovered));  
          good = false;          
        }
        
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kSaturated) ) ) {
          // leading edge recovery failed - still keep the information
          // about the saturation and do not flag as dead
          flagBits[ch]  |=  (0x1 <<  (RecHitFlags::RecHitFlags_kSaturated));  
          good = false;
        }
        
        //
        // AM: why do we have two tests one after the other checking almost the same thing??? 
        // Please clean up the code, ... also the original one!
        //        
        // uncalibRH.isSaturated() ---> 
        //         
        //                                   bool EcalUncalibratedRecHit::isSaturated() const {
        //                                     return EcalUncalibratedRecHit::checkFlag(kSaturated);
        //                                   }
        //
        //
        
        //         if (( 0x1 << (UncalibRecHitFlags::kSaturated) ) ) {
        if ( flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kSaturated) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kSaturated));  
          good = false;
        }
        
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kOutOfTime) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kOutOfTime));
          good = false;
        }
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kPoorReco) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kPoorReco));
          good = false;
        }
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kHasSwitchToGain6) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kHasSwitchToGain6));
        }
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kHasSwitchToGain1) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kHasSwitchToGain1));
        }
        
        
        if (good) {
          flagBits[ch] |= (0x1 << (RecHitFlags::RecHitFlags_kGood));
        }
        
        
        
        //         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kLeadingEdgeRecovered)) {
        //           rh.setFlag(EcalRecHit::kLeadingEdgeRecovered);
        //           good = false;
        //         }
        //         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kSaturated)) {
        //           // leading edge recovery failed - still keep the information
        //           // about the saturation and do not flag as dead
        //           rh.setFlag(EcalRecHit::kSaturated);
        //           good = false;
        //         }
        //         if (uncalibRH.isSaturated()) {
        //           rh.setFlag(EcalRecHit::kSaturated);
        //           good = false;
        //         }
        //         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kOutOfTime)) {
        //           rh.setFlag(EcalRecHit::kOutOfTime);
        //           good = false;
        //         }
        //         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kPoorReco)) {
        //           rh.setFlag(EcalRecHit::kPoorReco);
        //           good = false;
        //         }
        //         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kHasSwitchToGain6)) {
        //           rh.setFlag(EcalRecHit::kHasSwitchToGain6);
        //         }
        //         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kHasSwitchToGain1)) {
        //           rh.setFlag(EcalRecHit::kHasSwitchToGain1);
        //         }
        //         
        //         if (good)
        //           rh.setFlag(EcalRecHit::kGood);
        
        
        //        
        //        @ ecalrechit
        //        void setFlag(int flag) {flagBits_|= (0x1 << flag);}
        //        
        
        
        
        
        
        
        
        //         if (detid.subdetId() == EcalBarrel && (lasercalib < EBLaserMIN_ || lasercalib > EBLaserMAX_)) 
        //           myrechit.setFlag(EcalRecHit::kPoorCalib);
        //         if (detid.subdetId() == EcalEndcap && (lasercalib < EELaserMIN_ || lasercalib > EELaserMAX_)) 
        //           myrechit.setFlag(EcalRecHit::kPoorCalib);
        //         
        //         void setFlag(int flag) {flagBits_|= (0x1 << flag);}
        
        
        if (isBarrel  && (lasercalib < EBLaserMIN || lasercalib > EBLaserMAX)) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kPoorCalib));
          
        }
        if (!isBarrel && (lasercalib < EELaserMIN || lasercalib > EELaserMAX)) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kPoorCalib));
        }
        
        
  
  
  
  
        // recover, killing, and other stuff
        
//         if ( recoverEBIsolatedChannels_ || recoverEBFE_ || killDeadChannels_ )
//         {
//           edm::Handle< std::set<EBDetId> > pEBDetId;
//           const std::set<EBDetId> * detIds = nullptr;
//           evt.getByToken( ebDetIdToBeRecoveredToken_, pEBDetId);
//           detIds = pEBDetId.product();
//           
//           
//           if ( detIds ) {
//             edm::ESHandle<EcalChannelStatus> chStatus;
//             es.get<EcalChannelStatusRcd>().get(chStatus);
//             for( std::set<EBDetId>::const_iterator it = detIds->begin(); it != detIds->end(); ++it ) {
//               // get channel status map to treat dead VFE separately
//               EcalChannelStatusMap::const_iterator chit = chStatus->find( *it );
//               EcalChannelStatusCode chStatusCode;
//               if ( chit != chStatus->end() ) {
//                 chStatusCode = *chit;
//               } else {
//                 edm::LogError("EcalRecHitProducerError") << "No channel status found for xtal "
//                 << (*it).rawId()
//                 << "! something wrong with EcalChannelStatus in your DB? ";
//               }
// 
//     AM: 
// EcalUncalibratedRecHit(
//   const DetId& id, float ampl, float ped, float jit, float chi2, uint32_t flags = 0, uint32_t aux = 0);
// 
// 
// 
// 
//               EcalUncalibratedRecHit urh;
//               if ( chStatusCode.getStatusCode()  == EcalChannelStatusCode::kDeadVFE ) { // dead VFE (from DB info)
//                 // uses the EcalUncalibratedRecHit to pass the DetId info
//                 urh = EcalUncalibratedRecHit( *it, 0, 0, 0, 0, EcalRecHitWorkerBaseClass::EB_VFE );
//                 if ( recoverEBVFE_ || killDeadChannels_ ) workerRecover_->run( evt, urh, *ebRecHits );
//               } else {
//                 // uses the EcalUncalibratedRecHit to pass the DetId info
//                 urh = EcalUncalibratedRecHit( *it, 0, 0, 0, 0, EcalRecHitWorkerBaseClass::EB_single );
//                 if ( recoverEBIsolatedChannels_ || killDeadChannels_ ) workerRecover_->run( evt, urh, *ebRecHits );
//               }
//               
//             }
//           }
//         }
  
//   
//   
//   if ( recoverEEIsolatedChannels_ || recoverEEVFE_ || killDeadChannels_ )
//   {
//     edm::Handle< std::set<EEDetId> > pEEDetId;
//     const std::set<EEDetId> * detIds = nullptr;
//     
//     evt.getByToken( eeDetIdToBeRecoveredToken_, pEEDetId);
//     detIds = pEEDetId.product();
//     
//     if ( detIds ) {
//       edm::ESHandle<EcalChannelStatus> chStatus;
//       es.get<EcalChannelStatusRcd>().get(chStatus);
//       for( std::set<EEDetId>::const_iterator it = detIds->begin(); it != detIds->end(); ++it ) {
//         // get channel status map to treat dead VFE separately
//         EcalChannelStatusMap::const_iterator chit = chStatus->find( *it );
//         EcalChannelStatusCode chStatusCode;
//         if ( chit != chStatus->end() ) {
//           chStatusCode = *chit;
//         } else {
//           edm::LogError("EcalRecHitProducerError") << "No channel status found for xtal "
//           << (*it).rawId()
//           << "! something wrong with EcalChannelStatus in your DB? ";
//         }
//         EcalUncalibratedRecHit urh;
//         if ( chStatusCode.getStatusCode()  == EcalChannelStatusCode::kDeadVFE) { // dead VFE (from DB info)
//           // uses the EcalUncalibratedRecHit to pass the DetId info
//           urh = EcalUncalibratedRecHit( *it, 0, 0, 0, 0, EcalRecHitWorkerBaseClass::EE_VFE );
//           if ( recoverEEVFE_ || killDeadChannels_ ) workerRecover_->run( evt, urh, *eeRecHits );
//         } else {
//           // uses the EcalUncalibratedRecHit to pass the DetId info
//           urh = EcalUncalibratedRecHit( *it, 0, 0, 0, 0, EcalRecHitWorkerBaseClass::EE_single );
//           if ( recoverEEIsolatedChannels_ || killDeadChannels_ ) workerRecover_->run( evt, urh, *eeRecHits );
//         }
//       }
//     }
//   }
//   
    
    

    
    //
    // Structure:
    //  EB
    //  EE
    //
    //
    //  - single MVA
    //  - democratic sharing
    //  - kill all the other cases
    //
    
        bool is_Single = false;
        bool is_FE     = false;
        bool is_VFE    = false;
        
        bool is_recoverable = false; // DetIdToBeRecovered
        
        
//         // Integrity errors
//         edm::Handle<EBDetIdCollection> ebIntegrityGainErrors;
//         ev.getByToken(ebIntegrityGainErrorsToken_, ebIntegrityGainErrors);
//         ebDetIdColls.push_back(ebIntegrityGainErrors);
//         
//         edm::Handle<EBDetIdCollection> ebIntegrityGainSwitchErrors;
//         ev.getByToken(ebIntegrityGainSwitchErrorsToken_, ebIntegrityGainSwitchErrors);
//         ebDetIdColls.push_back(ebIntegrityGainSwitchErrors);
//         
//         edm::Handle<EBDetIdCollection> ebIntegrityChIdErrors;
//         ev.getByToken(ebIntegrityChIdErrorsToken_, ebIntegrityChIdErrors);
//         ebDetIdColls.push_back(ebIntegrityChIdErrors);
//         



/* 
 * find isolated dead channels (from DB info)           --> chStatus 10, 11, 12
 * and group of dead channels w/ trigger(from DB info)  --> chStatus 13
 * in interesting regions flagged by SRP
 */
//     const int flag = (*chit).getStatusCode();
//     if (flag >= 10 && flag <= 12) {  // FIXME -- avoid hardcoded values...
//       ebDetIdToRecover->insert(*itId);
//     } else if (flag == 13 || flag == 14) {  // FIXME -- avoid hardcoded values...
//       ebTTDetIdToRecover->insert((*itId).tower());
      
      
        if ( dbstatus == 10 ||  dbstatus == 11 ||  dbstatus == 12 ) {
          is_recoverable = true;
        }
        
        
        if (is_recoverable) {
          if (dbstatus == EcalChannelStatusCode_Code::kDeadVFE) {
            is_VFE = true;
          }
          else if (dbstatus == EcalChannelStatusCode_Code::kDeadVFE) {
            is_FE = true;
          }
          else {
            is_Single = true;
          }
          
          
          // EB
          if (isBarrel) {
            if (is_Single || is_FE || is_VFE) {           
              // single MVA
              if (is_Single && (recoverEBIsolatedChannels || !killDeadChannels) ) {
               
//                 
//                 https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerRecover.cc#L131
//                 
//                 float ebEn = ebDeadChannelCorrector.correct(
//                   detId, result, singleRecoveryMethod_, singleRecoveryThreshold_, sum8RecoveryThreshold_, &AcceptRecHit);
//                 EcalRecHit hit(detId, ebEn, 0., EcalRecHit::kDead);
//                 
                  
              }
              // decmocratic sharing
              else if (is_FE && (recoverEBFE || !killDeadChannels) ) {
               
                
                //                 for (std::vector<DetId>::const_iterator dit = vid.begin(); dit != vid.end(); ++dit) {
                //                   if (alreadyInserted(*dit))
                //                     continue;
                //                   float theta = ebGeom_->getGeometry(*dit)->getPosition().theta();
                //                   float tpEt = ecalScale_.getTPGInGeV(tp->compressedEt(), tp->id());
                //                   if (checkChannelStatus(*dit, dbStatusToBeExcludedEB_)) {
                //                     EcalRecHit hit(*dit, tpEt / ((float)vid.size()) / sin(theta), 0.);
                //                     hit.setFlag(EcalRecHit::kTowerRecovered);
                //                     if (tp->compressedEt() == 0xFF)
                //                       hit.setFlag(EcalRecHit::kTPSaturated);
                //                     if (tp->sFGVB())
                //                       hit.setFlag(EcalRecHit::kL1SpikeFlag);
                //                     insertRecHit(hit, result);
                //                   }
                
                
                
                
              }
              // kill all the other cases
              else {
                energy[ch] = 0.;  // Need to set also the flags ...
              }
            }
          }
          // EE
          else { 
            if (is_Single || is_FE || is_VFE) {           
              // single MVA
              if (is_Single && (recoverEBIsolatedChannels || !killDeadChannels) ) {
                
//                 
//                 https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerRecover.cc#L149
//                 
//                 float eeEn = eeDeadChannelCorrector.correct(
//                   detId, result, singleRecoveryMethod_, singleRecoveryThreshold_, sum8RecoveryThreshold_, &AcceptRecHit);
//                 EcalRecHit hit(detId, eeEn, 0., EcalRecHit::kDead);
//                 
                
              }
              // decmocratic sharing
              else if (is_FE && (recoverEBFE || !killDeadChannels) ) {
  
                
                //                
                //  Code is definitely too long ...              
                //                
                //                
                //                
                //                
                //                
                
              }
              // kill all the other cases
              else {
                energy[ch] = 0.;  // Need to set also the flags ...
              }
            }
          }
          
        }   
    
  
  
//         // EB
//         if (isBarrel) {
//           
//           if ( recoverEBIsolatedChannels || recoverEBFE || killDeadChannels ) {
//             // set energy to 0 by hand
//             //           if ( chStatusCode.getStatusCode()  == EcalChannelStatusCode::kDeadVFE ) { // dead VFE (from DB info)
//             if (dbstatus == EcalChannelStatusCode_Code::kDeadVFE) {
//               // EcalRecHitWorkerBaseClass::EB_VFE
//               energy[ch] = 0.;  // --> not 0!
//             }
//             else {
//               // EcalRecHitWorkerBaseClass::EB_single
//               energy[ch] = 0.;  // --> not 0!
//             }
//           }  
//         }
//         // EE
//         else {
// 
//           if ( recoverEEIsolatedChannels || recoverEBFE || killDeadChannels ) {
//             // set energy to 0 by hand
//             //           if ( chStatusCode.getStatusCode()  == EcalChannelStatusCode::kDeadVFE ) { // dead VFE (from DB info)
//             if (dbstatus == EcalChannelStatusCode_Code::kDeadVFE) {
//               // EcalRecHitWorkerBaseClass::EE_VFE
//               energy[ch] = 0.;  // --> not 0!
//             }      
//             else {
//               // EcalRecHitWorkerBaseClass::EE_single
//               energy[ch] = 0.;  // --> not 0!
//             }
//           }  
//           
//         }
  
  
  
  
  
  
//   if ( recoverEBFE_ || killDeadChannels_ )
//   {
//     edm::Handle< std::set<EcalTrigTowerDetId> > pEBFEId;
//     const std::set<EcalTrigTowerDetId> * ttIds = nullptr;
//     
//     evt.getByToken( ebFEToBeRecoveredToken_, pEBFEId);
//     ttIds = pEBFEId.product();
//     
//     if ( ttIds ) {
//       for( std::set<EcalTrigTowerDetId>::const_iterator it = ttIds->begin(); it != ttIds->end(); ++it ) {
//         // uses the EcalUncalibratedRecHit to pass the DetId info
//         int ieta = (((*it).ietaAbs()-1)*5+1)*(*it).zside(); // from EcalTrigTowerConstituentsMap
//         int iphi = (((*it).iphi()-1)*5+11)%360;             // from EcalTrigTowerConstituentsMap
//         if( iphi <= 0 ) iphi += 360;                        // from EcalTrigTowerConstituentsMap
//         EcalUncalibratedRecHit urh( EBDetId(ieta, iphi, EBDetId::ETAPHIMODE), 0, 0, 0, 0, EcalRecHitWorkerBaseClass::EB_FE );
//         workerRecover_->run( evt, urh, *ebRecHits );
//       }
//     }
//   }
//   
//   if ( recoverEEFE_ || killDeadChannels_ )
//   {
//     edm::Handle< std::set<EcalScDetId> > pEEFEId;
//     const std::set<EcalScDetId> * scIds = nullptr;
//     
//     evt.getByToken( eeFEToBeRecoveredToken_, pEEFEId);
//     scIds = pEEFEId.product();
//     
//     
//     if ( scIds ) {
//       for( std::set<EcalScDetId>::const_iterator it = scIds->begin(); it != scIds->end(); ++it ) {
//         // uses the EcalUncalibratedRecHit to pass the DetId info
//         if (EEDetId::validDetId( ((*it).ix()-1)*5+1, ((*it).iy()-1)*5+1, (*it).zside() )) {
//           EcalUncalibratedRecHit urh( EEDetId( ((*it).ix()-1)*5+1, ((*it).iy()-1)*5+1, (*it).zside() ), 0, 0, 0, 0, EcalRecHitWorkerBaseClass::EE_FE );
//           workerRecover_->run( evt, urh, *eeRecHits );
//         }
//       }
//     }
//   }
//   
  
  
  
//         // EB
//         if (isBarrel) {
//           
//           if ( recoverEBFE || killDeadChannels ) {
//       
//             // set to "0" the amplitudes of all the rechits that were not recovered 
//             // and the seed (really?) of the trigger tower of "ebFEToBeRecoveredToken_" ,
//             // e.g.     ebFEToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:ebFE"),
//             //          eeFEToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:eeFE"),
//             //   --> EcalDetIdToBeRecoveredProducer
//             //            
//             
//             // --> not 0!
//             
//           }
//         }
//         // EE
//         else {
// 
//           if ( recoverEEFE || killDeadChannels ) {
//             
//             // set to "0" the amplitudes of all the rechits that were not recovered 
//             // and the seed (really?) of the trigger tower of "ebFEToBeRecoveredToken_" ,
//             // e.g.     ebFEToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:ebFE"),
//             //          eeFEToBeRecovered = cms.InputTag("ecalDetIdToBeRecovered:eeFE"),
//             //   --> EcalDetIdToBeRecoveredProducer
//             //
//             
//             // --> not 0!
// 
//           }
//           
//         }
        
        
        
        
        
        
        
        
        
        
        
        
        
//         if ( killDeadChannels_ ) {
//           if (    (flags == EcalRecHitWorkerRecover::EB_single && !recoverEBIsolatedChannels_)
//             || (flags == EcalRecHitWorkerRecover::EE_single && !recoverEEIsolatedChannels_)
//             || (flags == EcalRecHitWorkerRecover::EB_VFE && !recoverEBVFE_)
//             || (flags == EcalRecHitWorkerRecover::EE_VFE && !recoverEEVFE_)
//           ) {
//             EcalRecHit hit( detId, 0., 0., EcalRecHit::kDead );
//             hit.setFlag( EcalRecHit::kDead)  ;
//             insertRecHit( hit, result); // insert trivial rechit with kDead flag
//             return true;
//           } 
//           if ( flags == EcalRecHitWorkerRecover::EB_FE && !recoverEBFE_) {
//             EcalTrigTowerDetId ttDetId( ((EBDetId)detId).tower() );
//             std::vector<DetId> vid = ttMap_->constituentsOf( ttDetId );
//             for ( std::vector<DetId>::const_iterator dit = vid.begin(); dit != vid.end(); ++dit ) {
//               EcalRecHit hit( (*dit), 0., 0., EcalRecHit::kDead );
//               hit.setFlag( EcalRecHit::kDead ) ;
//               insertRecHit( hit, result ); // insert trivial rechit with kDead flag
//             }
//             if(logWarningEtThreshold_EB_FE_<0)return true; // if you don't want log warning just return true
//           }
//           if ( flags == EcalRecHitWorkerRecover::EE_FE && !recoverEEFE_) {
//             EEDetId id( detId );
//             EcalScDetId sc( 1+(id.ix()-1)/5, 1+(id.iy()-1)/5, id.zside() );
//             std::vector<DetId> eeC;
//             for(int dx=1; dx<=5; ++dx){
//               for(int dy=1; dy<=5; ++dy){
//                 int ix = (sc.ix()-1)*5 + dx;
//                 int iy = (sc.iy()-1)*5 + dy;
//                 int iz = sc.zside();
//                 if(EEDetId::validDetId(ix, iy, iz)){
//                   eeC.push_back(EEDetId(ix, iy, iz));
//                 }
//               }
//             }
//             for ( size_t i = 0; i < eeC.size(); ++i ) {
//               EcalRecHit hit( eeC[i], 0., 0., EcalRecHit::kDead );
//               hit.setFlag( EcalRecHit::kDead ) ;
//               insertRecHit( hit, result ); // insert trivial rechit with kDead flag
//             }
//             if(logWarningEtThreshold_EE_FE_<0)   return true; // if you don't want log warning just return true
//           }
//         }
//         
//         if ( flags == EcalRecHitWorkerRecover::EB_single ) {
//           // recover as single dead channel
//           ebDeadChannelCorrector.setCaloTopology(caloTopology_.product());
//           
//           // channel recovery. Accepted new RecHit has the flag AcceptRecHit=TRUE
//           bool AcceptRecHit=true;
//           EcalRecHit hit = ebDeadChannelCorrector.correct( detId, result, singleRecoveryMethod_, singleRecoveryThreshold_, &AcceptRecHit);
//           
//           if ( hit.energy() != 0 and AcceptRecHit == true ) {
//             hit.setFlag( EcalRecHit::kNeighboursRecovered ) ;
//           } else {
//             // recovery failed
//             hit.setFlag( EcalRecHit::kDead ) ;
//           }
//           insertRecHit( hit, result );
//           
//         } else if ( flags == EcalRecHitWorkerRecover::EE_single ) {
//           // recover as single dead channel
//           eeDeadChannelCorrector.setCaloTopology(caloTopology_.product());
//           
//           // channel recovery. Accepted new RecHit has the flag AcceptRecHit=TRUE
//           bool AcceptRecHit=true;
//           EcalRecHit hit = eeDeadChannelCorrector.correct( detId, result, singleRecoveryMethod_, singleRecoveryThreshold_, &AcceptRecHit);
//           if ( hit.energy() != 0 and AcceptRecHit == true ) {
//             hit.setFlag( EcalRecHit::kNeighboursRecovered ) ;
//           } else {
//             // recovery failed
//             hit.setFlag( EcalRecHit::kDead ) ;
//           }
//           insertRecHit( hit, result );
//           
//         } else if ( flags == EcalRecHitWorkerRecover::EB_VFE ) {
//           // recover as dead VFE
//           EcalRecHit hit( detId, 0., 0.);
//           hit.setFlag( EcalRecHit::kDead ) ;
//           // recovery not implemented
//           insertRecHit( hit, result );
//         } else if ( flags == EcalRecHitWorkerRecover::EB_FE ) {
//           // recover as dead TT
//           
//           EcalTrigTowerDetId ttDetId( ((EBDetId)detId).tower() );
//           edm::Handle<EcalTrigPrimDigiCollection> pTPDigis;
//           evt.getByToken(tpDigiToken_, pTPDigis);
//           const EcalTrigPrimDigiCollection * tpDigis = nullptr;               
//           tpDigis = pTPDigis.product();
//           
//           EcalTrigPrimDigiCollection::const_iterator tp = tpDigis->find( ttDetId );
//           // recover the whole trigger tower
//           if ( tp != tpDigis->end() ) {
//             //std::vector<DetId> vid = ecalMapping_->dccTowerConstituents( ecalMapping_->DCCid( ttDetId ), ecalMapping_->iTT( ttDetId ) );
//             std::vector<DetId> vid = ttMap_->constituentsOf( ttDetId );
//             float tpEt  = ecalScale_.getTPGInGeV( tp->compressedEt(), tp->id() );
//             float tpEtThreshEB = logWarningEtThreshold_EB_FE_;
//             if(tpEt>tpEtThreshEB){
//               edm::LogWarning("EnergyInDeadEB_FE")<<"TP energy in the dead TT = "<<tpEt<<" at "<<ttDetId;
//             }
//             if ( !killDeadChannels_ || recoverEBFE_ ) {  
//               // democratic energy sharing
//               
//               for ( std::vector<DetId>::const_iterator dit = vid.begin(); dit != vid.end(); ++dit ) {
//                 if (alreadyInserted(*dit)) continue;
//                 float theta = ebGeom_->getGeometry(*dit)->getPosition().theta();
//                 float tpEt  = ecalScale_.getTPGInGeV( tp->compressedEt(), tp->id() );
//                 if(checkChannelStatus(*dit, dbStatusToBeExcludedEB_)){
//                   EcalRecHit hit( *dit, tpEt /((float)vid.size()) / sin(theta), 0.);
//                   hit.setFlag( EcalRecHit::kTowerRecovered ) ;
//                   if ( tp->compressedEt() == 0xFF ) hit.setFlag( EcalRecHit::kTPSaturated );
//                   if ( tp->sFGVB() ) hit.setFlag( EcalRecHit::kL1SpikeFlag );
//                   insertRecHit( hit, result );
//                 }
//               }
//             } else {
//               // tp not found => recovery failed
//               std::vector<DetId> vid = ttMap_->constituentsOf( ttDetId );
//               for ( std::vector<DetId>::const_iterator dit = vid.begin(); dit != vid.end(); ++dit ) {
//                 if (alreadyInserted(*dit)) continue;
//                 EcalRecHit hit( *dit,0., 0. );
//                 hit.setFlag( EcalRecHit::kDead ) ;
//                 insertRecHit( hit, result );
//               }
//             }
//           }
//         } else if ( flags == EcalRecHitWorkerRecover::EE_FE ) {
//           // Structure for recovery:
//           // ** SC --> EEDetId constituents (eeC) --> associated Trigger Towers (aTT) --> EEDetId constituents (aTTC)
//           // ** energy for a SC EEDetId = [ sum_aTT(energy) - sum_aTTC(energy) ] / N_eeC
//           // .. i.e. the total energy of the TTs covering the SC minus 
//           // .. the energy of the recHits in the TTs but not in the SC
//           //std::vector<DetId> vid = ecalMapping_->dccTowerConstituents( ecalMapping_->DCCid( ttDetId ), ecalMapping_->iTT( ttDetId ) );
//           // due to lack of implementation of the EcalTrigTowerDetId ix,iy methods in EE we compute Et recovered energies (in EB we compute E)
//           
//           EEDetId eeId( detId );
//           EcalScDetId sc( (eeId.ix()-1)/5+1, (eeId.iy()-1)/5+1, eeId.zside() );
//           std::set<DetId> eeC;
//           for(int dx=1; dx<=5; ++dx){
//             for(int dy=1; dy<=5; ++dy){
//               int ix = (sc.ix()-1)*5 + dx;
//               int iy = (sc.iy()-1)*5 + dy;
//               int iz = sc.zside();
//               if(EEDetId::validDetId(ix, iy, iz)){
//                 EEDetId id(ix, iy, iz);
//                 if (checkChannelStatus(id,dbStatusToBeExcludedEE_)){
//                   eeC.insert(id);
//                 } // check status
//               }
//             }
//           }
//           
//           edm::Handle<EcalTrigPrimDigiCollection> pTPDigis;
//           evt.getByToken(tpDigiToken_, pTPDigis);
//           const EcalTrigPrimDigiCollection * tpDigis = nullptr;
//           tpDigis = pTPDigis.product();
//           
//           // associated trigger towers
//           std::set<EcalTrigTowerDetId> aTT;
//           for ( std::set<DetId>::const_iterator it = eeC.begin(); it!=eeC.end(); ++it ) {
//             aTT.insert( ttMap_->towerOf( *it ) );
//           }
//           // associated trigger towers: total energy
//           float totE = 0;
//           // associated trigger towers: EEDetId constituents
//           std::set<DetId> aTTC;
//           bool atLeastOneTPSaturated = false;
//           for ( std::set<EcalTrigTowerDetId>::const_iterator it = aTT.begin(); it != aTT.end(); ++it ) {
//             // add the energy of this trigger tower
//             EcalTrigPrimDigiCollection::const_iterator itTP = tpDigis->find( *it );
//             if ( itTP != tpDigis->end() ) {
//               
//               std::vector<DetId> v = ttMap_->constituentsOf( *it );
//               
//               // from the constituents, remove dead channels
//               std::vector<DetId>::iterator ttcons = v.begin();
//               while (ttcons != v.end()){
//                 if (!checkChannelStatus(*ttcons,dbStatusToBeExcludedEE_)){
//                   ttcons=v.erase(ttcons);
//                 } else {
//                   ++ttcons;
//                 }
//               }// while 
//               
//               if ( itTP->compressedEt() == 0xFF ){ // In the case of a saturated trigger tower, a fraction
//                 atLeastOneTPSaturated = true; //of the saturated energy is put in: number of xtals in dead region/total xtals in TT *63.75
//                 
//                 //Alternative recovery algorithm that I will now investigate.
//                 //Estimate energy sums the energy in the working channels, then decides how much energy
//                 //to put here depending on that. Duncan 20101203
//                 
//                 totE += estimateEnergy(itTP->id().ietaAbs(), &result, eeC, v);
//                 
//                 /* 
//                  *                         These commented out lines use
//                  *                         64GeV*fraction of the TT overlapping the dead FE
//                  *                        
//                  *                      int count = 0;
//                  *                      for (std::vector<DetId>::const_iterator idsit = v.begin(); idsit != v.end(); ++ idsit){
//                  *                      std::set<DetId>::const_iterator itFind = eeC.find(*idsit);
//                  *                      if (itFind != eeC.end())
//                  *                      ++count;
//               }
//               //std::cout << count << ", " << v.size() << std::endl;
//               totE+=((float)count/(float)v.size())* ((it->ietaAbs()>26)?2*ecalScale_.getTPGInGeV( itTP->compressedEt(), itTP->id() ):ecalScale_.getTPGInGeV( itTP->compressedEt(), itTP->id() ));*/
//               }
//               else {totE += ((it->ietaAbs()>26)?2:1)*ecalScale_.getTPGInGeV( itTP->compressedEt(), itTP->id() );}
//               
//               
//               // get the trigger tower constituents
//               
//               if (itTP->compressedEt() == 0){ // If there's no energy in TT, the constituents are removed from the recovery.
//                 for (size_t i = 0 ; i < v.size(); ++i)
//                   eeC.erase(v[i]);
//               }
//               else if (itTP->compressedEt()!=0xFF){ //If it's saturated the energy has already been determined, so we do not want to subtract any channels
//                 for ( size_t j = 0; j < v.size(); ++j ) {
//                   aTTC.insert( v[j] );
//                 }
//               }
//               
//             }
//           }
//           // remove crystals of dead SC
//           // (this step is not needed if sure that SC crystals are not 
//           // in the recHit collection)
//           
//           for ( std::set<DetId>::const_iterator it = eeC.begin(); it != eeC.end(); ++it ) {
//             aTTC.erase(*it);
//           }
//           // compute the total energy for the dead SC
//           const EcalRecHitCollection * hits = &result;
//           for ( std::set<DetId>::const_iterator it = aTTC.begin(); it != aTTC.end(); ++it ) {
//             EcalRecHitCollection::const_iterator jt = hits->find( *it );
//             if ( jt != hits->end() ) {
//               float energy = jt->energy(); // Correct conversion to Et
//               float eta = geo_->getPosition(jt->id()).eta();
//               float pf = 1.0/cosh(eta);
//               // use Et instead of E, consistent with the Et estimation of the associated TT
//               totE -= energy*pf;
//             }
//           }
//           
//           
//           float scEt = totE;
//           float scEtThreshEE = logWarningEtThreshold_EE_FE_;
//           if(scEt>scEtThreshEE){
//             edm::LogWarning("EnergyInDeadEE_FE")<<"TP energy in the dead TT = "<<scEt<<" at "<<sc;
//           }
//           
//           // assign the energy to the SC crystals
//           if ( !killDeadChannels_ || recoverEEFE_ ) { // if eeC is empty, i.e. there are no hits 
//             // in the tower, nothing is returned. No negative values from noise.
//             for ( std::set<DetId>::const_iterator it = eeC.begin(); it != eeC.end(); ++it ) {
//               
//               float eta = geo_->getPosition(*it).eta(); //Convert back to E from Et for the recovered hits
//               float pf = 1.0/cosh(eta);
//               EcalRecHit hit( *it, totE / ((float)eeC.size()*pf), 0);
//               
//               if (atLeastOneTPSaturated) hit.setFlag(EcalRecHit::kTPSaturated );                            
//               hit.setFlag(EcalRecHit::kTowerRecovered); 
//               insertRecHit( hit, result );
//               
//             }// for
//           }// if 
//         }
//         return true;
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
  
  
  
  
  
  
  
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
        configParameters.recoverEBIsolatedChannels,
        configParameters.recoverEEIsolatedChannels,
        configParameters.recoverEBVFE,             
        configParameters.recoverEEVFE,             
        configParameters.recoverEBFE,             
        configParameters.recoverEEFE,              
        configParameters.EBLaserMIN,
        configParameters.EELaserMIN,
        configParameters.EBLaserMAX,
        configParameters.EELaserMAX,
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
        eventInputGPU.ebUncalibRecHits.flags, 
        eventInputGPU.eeUncalibRecHits.flags, 
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
