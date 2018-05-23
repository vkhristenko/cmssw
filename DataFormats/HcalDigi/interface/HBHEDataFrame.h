#ifndef DIGIHCAL_HBHEDATAFRAME_H
#define DIGIHCAL_HBHEDATAFRAME_H

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "FWCore/Common/interface/CudaDecls.h"
#include <vector>
#include <ostream>

/** \class HBHEDataFrame
      
Precision readout digi for HB and HE.

*/
class HBHEDataFrame {
public:
  typedef HcalDetId key_type; ///< For the sorted collection

  HOST DEVICE HBHEDataFrame() 
      : id_(0), size_(0), hcalPresamples_(0) {
      // for persistence
  }
  HOST DEVICE explicit HBHEDataFrame(const HcalDetId& id);
  
  HOST DEVICE const HcalDetId& id() const { return id_; }
  HOST DEVICE const HcalElectronicsId& elecId() const { return electronicsId_; }
  
  /// total number of samples in the digi 
  HOST DEVICE int size() const { return size_&0xF; }
  /// number of samples before the sample from the triggered beam crossing (according to the hardware)
  HOST DEVICE int presamples() const { return hcalPresamples_&0xF; }
   /// was ZS MarkAndPass?
  HOST DEVICE bool zsMarkAndPass() const { return (hcalPresamples_&0x10); }
  /// was ZS unsuppressed?
  HOST DEVICE bool zsUnsuppressed() const { return (hcalPresamples_&0x20); }
  /// zs crossing mask (which sums considered)
  HOST DEVICE uint32_t zsCrossingMask() const { return (hcalPresamples_&0x3FF000)>>12; }
 
  /// access a sample
  HOST DEVICE const HcalQIESample& operator[](int i) const { return data_[i]; }
  /// access a sample
  HOST DEVICE const HcalQIESample& sample(int i) const { return data_[i]; }

  /// offset of bunch number for this channel relative to nominal set in the unpacker (range is +7->-7.  -1000 indicates the data is invalid/unavailable)
  HOST DEVICE int fiberIdleOffset() const;
  
  /// validate appropriate DV and ER bits as well as capid rotation for the specified samples (default is all)
  HOST DEVICE bool validate(int firstSample=0, int nSamples=100) const;
  
  HOST DEVICE void setSize(int size);
  HOST DEVICE void setPresamples(int ps);
  HOST DEVICE void setZSInfo(bool unsuppressed, bool markAndPass, uint32_t crossingMask=0);
  HOST DEVICE void setSample(int i, const HcalQIESample& sam) { data_[i]=sam; }
  HOST DEVICE void setReadoutIds(const HcalElectronicsId& eid);
  HOST DEVICE void setFiberIdleOffset(int offset);
  
  static const int MAXSAMPLES = 10;
private:
  HcalDetId id_;
  HcalElectronicsId electronicsId_; 
  int size_;
  int hcalPresamples_;
  HcalQIESample data_[MAXSAMPLES];    
};

std::ostream& operator<<(std::ostream&, const HBHEDataFrame&);


#endif
