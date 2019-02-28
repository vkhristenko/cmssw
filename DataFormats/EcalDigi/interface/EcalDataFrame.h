#ifndef DIGIECAL_ECALDATAFRAME_H
#define DIGIECAL_ECALDATAFRAME_H

#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/DataFrame.h"

#define EcalMgpaBitwiseGain12 1
#define EcalMgpaBitwiseGain6  2
#define EcalMgpaBitwiseGain1  3
#define EcalMgpaBitwiseGain0  0

/** \class EcalDataFrame
      
*/
class EcalDataFrame {
public:
  constexpr EcalDataFrame() {}
  // EcalDataFrame(DetId i) :  m_data(i) {}
  constexpr EcalDataFrame(edm::DataFrame const & iframe) : m_data(iframe){} 

//  virtual ~EcalDataFrame() {} 

  constexpr DetId id() const { return m_data.id();}
    
  constexpr int size() const { return m_data.size();}

  constexpr EcalMGPASample operator[](int i) const { return m_data[i];}
  constexpr EcalMGPASample sample(int i) const { return m_data[i]; }

  // get the leading sample (the first non saturated sample)
  // starting from the fourth sample
  // (it relies on the fact that the unpaker will discard fixed gain0 DataFrame)
  // .. sample numbering: [0, 9]
  // .. return -1 in case of no saturation
  constexpr int lastUnsaturatedSample() const
  {
        int cnt = 0;
        for ( size_t i = 3; i < m_data.size(); ++i ) {
                cnt = 0;
                for ( size_t j = i; j < (i + 5) && j < m_data.size(); ++j ) {
                        if ( ((EcalMGPASample)m_data[j]).gainId() == EcalMgpaBitwiseGain0 ) ++cnt;
                }
                if ( cnt == 5 ) return i-1; // the last unsaturated sample
        }
        return -1; // no saturation found
  }
  // just the boolean method
  constexpr bool isSaturated() const { return ( lastUnsaturatedSample() != -1 ); }
    
  // FIXME (shall we throw??)
  constexpr void setSize(int){}
  // void setPresamples(int ps);
  constexpr void setSample(int i, EcalMGPASample sam) { m_data[i]=sam; }

  constexpr bool hasSwitchToGain6() const
  {
    for(unsigned int u=0; u<m_data.size(); u++) 
    {
      if ( ( static_cast<EcalMGPASample>(m_data[u]) ).gainId() == EcalMgpaBitwiseGain6 ) return true;
    }
    return false;
  }
  constexpr bool hasSwitchToGain1() const
  {
    for(unsigned int u=0; u<m_data.size(); u++) 
    {
      if ( ( static_cast<EcalMGPASample>(m_data[u]) ).gainId() == EcalMgpaBitwiseGain1 ) return true;
    }
    return false;
  }
  
  static constexpr int MAXSAMPLES = 10;

  constexpr edm::DataFrame const & frame() const { return m_data;}
  constexpr edm::DataFrame & frame() { return m_data;}

 private:
 
  edm::DataFrame m_data;
  
};
  
#endif
