#ifndef DIGIECAL_EBDATAFRAME_H
#define DIGIECAL_EBDATAFRAME_H

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iosfwd>



/** \class EBDataFrame
      
*/
class EBDataFrame : public EcalDataFrame {
public:
  typedef EBDetId key_type; ///< For the sorted collection
  typedef EcalDataFrame Base;

  constexpr EBDataFrame() {}
  // EBDataFrame(DetId i) :  Base(i) {}
  constexpr EBDataFrame(edm::DataFrame const & base) : Base(base) {}
  constexpr EBDataFrame(EcalDataFrame const & base) : Base(base) {}

  /** estimator for a signal being a spike
   *  based on ratios between 4th, 5th and 6th sample
   */
  constexpr float spikeEstimator() const
  {
        if ( size() != 10 ) {
                edm::LogError("InvalidNumberOfSamples") << "This method only applies to signals sampled 10 times ("
                        << size() << " samples found)";
                return 10.;
        }
        // skip faulty channels
        if ( sample(5).adc() == 0 ) return 10.;
        size_t imax = 0;
        int maxAdc = 0;
        for ( int i = 0; i < size(); ++i ) {
                if ( sample(i).adc() > maxAdc ) {
                        imax = i;
                        maxAdc = sample(i).adc();
                }
        }
        // skip early signals
        if ( imax < 4 ) return 10.;
        float ped = 1./3. * (sample(0).adc() + sample(1).adc() + sample(2).adc());
        return 0.18*(sample(4).adc()-ped)/(sample(5).adc()-ped) + (sample(6).adc()-ped)/(sample(5).adc()-ped);
  }
    
  //~EBDataFrame() override {}

  key_type id() const { return Base::id(); }

};

std::ostream& operator<<(std::ostream&, const EBDataFrame&);

#endif
