/**
 * Author: Giovanni Franzoni, UMN
 * Created: 08 May 2012
 * $Id: EcalSampleMask.cc,v 1.1 2012/05/10 08:22:09 argiro Exp $
 **/

#include <cassert>
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


using edm::LogError;

EcalSampleMask::EcalSampleMask( const std::vector<unsigned int> &ebmask, const std::vector<unsigned int> &eemask) {
  setEcalSampleMaskRecordEB( ebmask );
  setEcalSampleMaskRecordEE( eemask );
}

void EcalSampleMask::setEcalSampleMaskRecordEB( const std::vector<unsigned int> & ebmask ) {
  
  // check that size of the vector is adequate 
  if( ebmask.size() != static_cast<unsigned int>(EcalDataFrame::MAXSAMPLES) ){
    LogError("DataMismatch")<< " in EcalSampleMask::setEcalSampleMaskRecordEB size of ebmask (" << ebmask.size() << ") need to be: " << EcalDataFrame::MAXSAMPLES 
	      << ". Bailing out."<< std::endl;
    assert(0);
  }

  // check that values of vector are allowed
  for (unsigned int s=0; s<ebmask.size(); s++ ) {
    if    ( ebmask.at(s)==0 || ebmask.at(s)==1  ) {;}
    else {
      LogError("DataMismatch")<< "in EcalSampleMask::setEcalSampleMaskRecordEB ebmask can only have values 0 or 1, while " << ebmask.at(s) << " was found. Bailing out. " << std::endl;
      assert(0);
    }
  }
  
  // ordering of bits:
  // ebmask.at(0)                         refers to the first sample read out and is mapped into the _most_ significant bit of sampleMaskEB_ 
  // ebmask.at(EcalDataFrame::MAXSAMPLES) refers to the last  sample read out and is mapped into the _least_ significant bit of sampleMaskEB_ 
  sampleMaskEB_=0;
  for (unsigned int sampleId=0; sampleId<ebmask.size(); sampleId++ ) {
    sampleMaskEB_ |= (0x1 << (EcalDataFrame::MAXSAMPLES -(sampleId+1) ));
  }

}

void EcalSampleMask::setEcalSampleMaskRecordEE( const std::vector<unsigned int> & eemask ) {

  // check that size of the vector is adequate 
  if( eemask.size() != static_cast<unsigned int>(EcalDataFrame::MAXSAMPLES) ){
      LogError("DataMismatch") << " in EcalSampleMask::setEcalSampleMaskRecordEE size of eemask (" << eemask.size() << ") need to be: " << EcalDataFrame::MAXSAMPLES 
	      << ". Bailing out."<< std::endl;
    assert(0);
  }

  // check that values of vector are allowed
  for (unsigned int s=0; s<eemask.size(); s++ ) {
    if    ( eemask.at(s)==0 || eemask.at(s)==1  ) {;}
    else {
      LogError("DataMismatch") << "in EcalSampleMask::setEcalSampleMaskRecordEE eemask can only have values 0 or 1, while " << eemask.at(s) << " was found. Bailing out. " << std::endl;
      assert(0);
    }
  }
  
  // ordering of bits:
  // eemask.at(0)                         refers to the first sample read out and is mapped into the _most_ significant bit of sampleMaskEE_ 
  // eemask.at(EcalDataFrame::MAXSAMPLES) refers to the last  sample read out and is mapped into the _least_ significant bit of sampleMaskEE_ 
  sampleMaskEE_=0;
  for (unsigned int sampleId=0; sampleId<eemask.size(); sampleId++ ) {
    sampleMaskEE_ |= (0x1 << (EcalDataFrame::MAXSAMPLES -(sampleId+1) ));
  }

}
