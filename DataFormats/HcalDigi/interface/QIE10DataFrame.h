#ifndef DATAFORMATS_HCALDIGI_QIE10DATAFRAME_H
#define DATAFORMATS_HCALDIGI_QIE10DATAFRAME_H

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/Common/interface/DataFrame.h"
#include <ostream>

/** Precision readout digi from QIE10 including TDC information

 */
class QIE10DataFrame {
public:

  static const int WORDS_PER_SAMPLE = 2;
  static const int HEADER_WORDS = 1;
  static const int FLAG_WORDS = 1;

  constexpr QIE10DataFrame() { }
  constexpr QIE10DataFrame(edm::DataFrame const & df) : m_data(df) { }

  class Sample {
  public:
    typedef uint32_t wide_type;

    constexpr Sample(const edm::DataFrame& frame, edm::DataFrame::size_type i) : word1_(frame[i]), word2_(frame[i+1]) { }
    constexpr Sample(const edm::DataFrame::data_type& word1, const edm::DataFrame::data_type& word2) : word1_(word1), word2_(word2) {}
    constexpr explicit Sample(wide_type wide);

    static const int MASK_ADC = 0xFF;
    static const int MASK_LE_TDC = 0x3F;
    static const int MASK_TE_TDC = 0x1F;
    static const int OFFSET_TE_TDC = 6;
    static const int MASK_SOI = 0x2000;
    static const int MASK_OK = 0x1000;
    static const int MASK_CAPID = 0x3;
    static const int OFFSET_CAPID = 12;

    constexpr inline int adc() const { return word1_&MASK_ADC; }
    constexpr inline int le_tdc() const { return word2_&MASK_LE_TDC; }
    constexpr inline int te_tdc() const { return (word2_>>OFFSET_TE_TDC)&MASK_TE_TDC; }
    constexpr inline bool ok() const { return word1_&MASK_OK; }
    constexpr inline bool soi() const { return word1_&MASK_SOI; }
    constexpr inline int capid() const { return (word2_>>OFFSET_CAPID)&MASK_CAPID; }
    constexpr inline edm::DataFrame::data_type raw(edm::DataFrame::size_type i) const
        { return (i > WORDS_PER_SAMPLE) ? 0 : ( (i==1) ? word2_ : word1_ ); }
    constexpr wide_type wideRaw() const;

  private:
    edm::DataFrame::data_type word1_;
    edm::DataFrame::data_type word2_;
  };

  constexpr void copyContent(const QIE10DataFrame& src);

  /// Get the detector id
  constexpr DetId detid() const { return DetId(m_data.id()); }
  constexpr edm::DataFrame::id_type id() const { return m_data.id(); }
  /// more accessors
  constexpr edm::DataFrame::size_type size() const { return m_data.size(); }
  /// iterators
  constexpr edm::DataFrame::iterator begin() { return m_data.begin(); }
  constexpr edm::DataFrame::iterator end() { return m_data.end(); }
  constexpr edm::DataFrame::const_iterator begin() const { return m_data.begin(); }
  constexpr edm::DataFrame::const_iterator end() const { return m_data.end(); }
  /// total number of samples in the digi
  constexpr int samples() const { return (size()-HEADER_WORDS-FLAG_WORDS)/WORDS_PER_SAMPLE; }
  /// for backward compatibility
  constexpr int presamples() const;
  /// get the flavor of the frame
  static const int OFFSET_FLAVOR = 12;
  static const int MASK_FLAVOR = 0x7;
  constexpr int flavor() const { return ((m_data[0]>>OFFSET_FLAVOR)&MASK_FLAVOR); }
  /// was there a link error?
  static const int MASK_LINKERROR = 0x800;
  constexpr bool linkError() const { return m_data[0]&MASK_LINKERROR; }
  /// was this a mark-and-pass ZS event?
  static const int MASK_MARKPASS = 0x100;
  constexpr bool zsMarkAndPass() const {return m_data[0]&MASK_MARKPASS; }
  /// set ZS params
  constexpr void setZSInfo(bool markAndPass);
  /// get the sample
  constexpr inline Sample operator[](edm::DataFrame::size_type i) const { return Sample(m_data,i*WORDS_PER_SAMPLE+HEADER_WORDS); }
  /// set the sample contents
  constexpr void setSample(edm::DataFrame::size_type isample, int adc, int le_tdc, int te_tdc, int capid, bool soi=false, bool ok=true);
  /// get the flag word
  constexpr uint16_t flags() const { return m_data[size()-1]; }
  /// set the flag word
  constexpr void setFlags(uint16_t v);
  
  private:
   edm::DataFrame m_data;

};

std::ostream& operator<<(std::ostream&, const QIE10DataFrame&);


#endif // DATAFORMATS_HCALDIGI_QIE10DATAFRAME_H
