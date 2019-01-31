#ifndef CondFormats_EcalObjects_EcalSampleMask_H
#define CondFormats_EcalObjects_EcalSampleMask_H
/**
 * Author: Giovanni Franzoni, UMN
 * Created: 09 Apr 2012
 * $Id: EcalSampleMask.h,v 1.1 2012/05/10 08:22:10 argiro Exp $
 **/

#include "CondFormats/Serialization/interface/Serializable.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <vector>
#include <cassert>

using edm::LogError;

class EcalSampleMask {
  public:
    constexpr EcalSampleMask() 
      : sampleMaskEB_{
          static_cast<unsigned int>(std::pow(2, EcalDataFrame::MAXSAMPLES)-1)},
        sampleMaskEE_{
          static_cast<unsigned int>(std::pow(2, EcalDataFrame::MAXSAMPLES)-1)}
    {
      // by default, all samples are set as active
    }

    // construct from pre-organized binary words 
    constexpr EcalSampleMask(const unsigned int ebmask, const unsigned int eemask) 
      : sampleMaskEB_{ebmask},
        sampleMaskEE_{eemask}
    {}
    // constructor from an ordered set of switches, one per sample 
    EcalSampleMask( const std::vector<unsigned int> &ebmask, const std::vector<unsigned int> &eemask);

    constexpr void setEcalSampleMaskRecordEB( const unsigned int mask ) { sampleMaskEB_ = mask; }
    constexpr void setEcalSampleMaskRecordEE( const unsigned int mask ) { sampleMaskEE_ = mask; }
    void setEcalSampleMaskRecordEB( const std::vector<unsigned int> & ebmask );
    void setEcalSampleMaskRecordEE( const std::vector<unsigned int> & eemask );
    
    constexpr float getEcalSampleMaskRecordEB() const { return sampleMaskEB_; }
    constexpr float getEcalSampleMaskRecordEE() const { return sampleMaskEE_; }
    void print(std::ostream& s) const {
      s << "EcalSampleMask: EB " << sampleMaskEB_ << "; EE " << sampleMaskEE_ ;
    }

    bool useSampleEB (const int sampleId) const {
      
      if( sampleId >= EcalDataFrame::MAXSAMPLES ){
        LogError("DataMismatch")<< "in EcalSampleMask::useSampleEB only sampleId up to: "  << EcalDataFrame::MAXSAMPLES 
              << " can be used, while: " << sampleId << " was found. Bailing out." << std::endl;
        assert(0);
      }
      
      // ordering convention:
      // ebmask.at(0)                         refers to the first sample read out and is mapped into the _most_ significant bit of sampleMaskEB_ 
      // ebmask.at(EcalDataFrame::MAXSAMPLES) refers to the last  sample read out and is mapped into the _least_ significant bit of sampleMaskEB_ 
      return ( sampleMaskEB_ & ( 0x1<< (EcalDataFrame::MAXSAMPLES -(sampleId+1) )) );
      
    }
    bool useSampleEE (const int sampleId) const {
      
      if( sampleId >= EcalDataFrame::MAXSAMPLES ){
        LogError("DataMismatch")<< "in EcalSampleMask::useSampleEE only sampleId up to: "  << EcalDataFrame::MAXSAMPLES 
              << " can be used, while: " << sampleId << " was found. Bailing out." << std::endl;
        assert(0);
      }
      
      // ordering convention:
      // ebmask.at(0)                         refers to the first sample read out and is mapped into the _most_ significant bit of sampleMaskEB_ 
      // ebmask.at(EcalDataFrame::MAXSAMPLES) refers to the last  sample read out and is mapped into the _least_ significant bit of sampleMaskEB_ 
      return ( sampleMaskEE_ & ( 0x1<< (EcalDataFrame::MAXSAMPLES -(sampleId+1) )) );
      
    }

    bool useSample  (const int sampleId, DetId &theCrystalId) const {
      
      if( sampleId >= EcalDataFrame::MAXSAMPLES ){
        LogError("DataMismatch")<< "in EcalSampleMask::useSample only sampleId up to: "  << EcalDataFrame::MAXSAMPLES 
              << " can be used, while: " << sampleId << " was found. Bailing out." << std::endl;
        assert(0);
      }
      
      
      if       (theCrystalId.subdetId()==EcalBarrel) {
        return useSampleEB ( sampleId );
      }
      else if  (theCrystalId.subdetId()==EcalEndcap) {
        return useSampleEE ( sampleId );
      }
      else {
        LogError("DataMismatch")<< "EcalSampleMaskuseSample::useSample can only be called for EcalBarrel or EcalEndcap DetID" << std::endl; 
        assert(0);
      }
      
    }

  private:
    unsigned int sampleMaskEB_;
    unsigned int sampleMaskEE_;


  COND_SERIALIZABLE;
};


#endif
