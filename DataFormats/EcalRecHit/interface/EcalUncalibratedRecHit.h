#ifndef DATAFORMATS_ECALUNCALIBRATEDRECHIT
#define DATAFORMATS_ECALUNCALIBRATEDRECHIT

#include <cmath>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class EcalUncalibratedRecHit {

  public:
  
  typedef DetId key_type;

  enum Flags {
          kGood=-1,                 // channel is good (mutually exclusive with other states)  setFlagBit(kGood) reset flags_ to zero 
          kPoorReco,                // channel has been badly reconstructed (e.g. bad shape, bad chi2 etc.)
          kSaturated,               // saturated channel
          kOutOfTime,               // channel out of time
          kLeadingEdgeRecovered,    // saturated channel: energy estimated from the leading edge before saturation
          kHasSwitchToGain6,        // at least one data frame is in G6
          kHasSwitchToGain1         // at least one data frame is in G1
          
  };

    constexpr EcalUncalibratedRecHit() : 
      amplitude_(0.), amplitudeError_(0.), pedestal_(0.), jitter_(0.), 
      chi2_(10000.),
      OOTamplitudes_{0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
      OOTchi2_{0}, flags_{0}, aux_{0}
    {}

    constexpr EcalUncalibratedRecHit(const DetId& id, float ampl, float ped,
                         float jit, float chi2, uint32_t flags = 0, 
                         uint32_t aux = 0)
      : amplitude_(ampl), amplitudeError_(0.), pedestal_(ped), 
        jitter_(jit), chi2_(chi2), 
        OOTamplitudes_{0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
        OOTchi2_{0}, flags_{flags}, aux_{aux}, id_{id} {}


  constexpr float amplitude() const { return amplitude_; }
  constexpr float amplitudeError() const { return amplitudeError_; }
  constexpr float pedestal() const { return pedestal_; }
  constexpr float jitter() const { return jitter_; }
  constexpr float chi2() const { return chi2_; }
  constexpr float outOfTimeAmplitude(int bx) const { return OOTamplitudes_[bx]; }

  constexpr uint32_t flags() const { return flags_; }
    constexpr float jitterError() const
    {
            // stored in ps, but return BXs to match with jitter units
            uint32_t jitterErrorBits = 0xFF & aux_;
            // all bits off --> time reco bailed out (return negative value)
            if( (0xFF & jitterErrorBits) == 0x00)
                    return -1;
            // all bits on  --> time error over 5 ns (return large value)
            if( (0xFF & jitterErrorBits) == 0xFF)
                    return 10000;

            float LSB = 1.26008;
            uint8_t exponent = jitterErrorBits>>5;
            uint8_t significand = jitterErrorBits & ~(0x7<<5);
            return (float)(std::pow(2,exponent)*significand*LSB)/(25.*1000);
    }
    constexpr uint8_t jitterErrorBits() const
    {
            uint8_t jitterErrorBits = 0xFF & aux_;
            return jitterErrorBits;
    }
  constexpr DetId  id() const { return id_; }

  constexpr void setAmplitude( float amplitude ) { amplitude_ = amplitude; }
  constexpr void setAmplitudeError( float amplitudeerror ) { amplitudeError_ = amplitudeerror; }
  constexpr void setPedestal( float pedestal ) { pedestal_ = pedestal; }
  constexpr void setJitter( float jitter ) { jitter_ = jitter; }
  constexpr void setChi2( float chi2 ) { chi2_ = chi2; }
  constexpr void setOutOfTimeAmplitude( int bx, float amplitude ) { OOTamplitudes_[bx] = amplitude; }

    constexpr void setJitterError( float jitterErr )
    {
            // use 8 bits (3 exp, 5 mant) and store in ps
            // has range of 5 ps - 5000 ps
            // expect input in BX units
            // all bits off --> time reco bailed out
            if(jitterErr <= 0)
            {
                    aux_ = (~0xFF & aux_);
                    return;
            }
            // all bits on  --> time error over 5 ns
            if(25*jitterErr >= 5)
            {
                    aux_ = (0xFF | aux_);
                    return;
            }

            float LSB = 1.26008;
            float quantityInLSB = (1000*25*jitterErr)/LSB;
            int log2OfQuantity = (int) (log2( quantityInLSB ));
            int exponentTmp = log2OfQuantity - 4;
            uint8_t exponent=0;
            if (exponentTmp>0) exponent = exponentTmp;
            uint8_t significand = (int) ( std::lround( quantityInLSB / std::pow(2,exponent) )   );
            uint32_t jitterErrorBits = exponent<<5 | significand;
      
            if( (0xFF & jitterErrorBits) == 0xFF)
              jitterErrorBits = 0xFE;
            if( (0xFF & jitterErrorBits) == 0x00)
              jitterErrorBits = 0x01;

            aux_ = (~0xFF & aux_) | (jitterErrorBits & 0xFF);

    }
  constexpr void setFlags( uint32_t flags ) { flags_ = flags; }
  constexpr void setId( DetId id ) { id_ = id; }
  constexpr void setAux( uint32_t aux ) { aux_ = aux; }
    constexpr void setFlagBit(EcalUncalibratedRecHit::Flags flag){
           if  (flag == kGood) {
              //then set all bits to zero;
              flags_  = 0;
              return;
          }
         // else set the flagbit
         flags_|= 0x1 <<  flag;  
    }
    constexpr bool checkFlag(EcalUncalibratedRecHit::Flags flag) const {
           if(flag == kGood){ if ( ! flags_ ) return true;else return false;} // if all flags are unset, then hit is good
           return  flags_ & ( 0x1<<flag);
    }

    constexpr bool isSaturated() const {
      return EcalUncalibratedRecHit::checkFlag(kSaturated);
    }
    constexpr bool isJitterValid() const
    {
            if(jitterError() <= 0)
              return false;
            else
              return true;
    }
    constexpr bool isJitterErrorValid() const
    {
            if(!isJitterValid())
              return false;
            if(jitterError() >= 10000)
              return false;

            return true;
    }

 private:
  float amplitude_;           //< Reconstructed amplitude
  float amplitudeError_;      //< Reconstructed amplitude uncertainty
  float pedestal_;            //< Reconstructed pedestal
  float jitter_;              //< Reconstructed time jitter
  float chi2_;                //< Chi2 of the pulse
  float OOTamplitudes_[EcalDataFrame::MAXSAMPLES];       //< Out-Of-Time reconstructed amplitude, one for each active BX, from readout sample 0 to 9
  float OOTchi2_;             //< Out-Of-Time Chi2 
  uint32_t flags_;            //< flag to be propagated to RecHit
  uint32_t aux_;              //< aux word; first 8 bits contain time (jitter) error
  DetId  id_;                 //< Detector ID
};

#endif
