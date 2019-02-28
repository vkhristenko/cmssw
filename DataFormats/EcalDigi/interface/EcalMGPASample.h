#ifndef DIGIECAL_ECALMGPASAMPLE_H
#define DIGIECAL_ECALMGPASAMPLE_H

#include <iosfwd>
#include <boost/cstdint.hpp>

namespace ecalMGPA {
  typedef uint16_t sample_type;
 
  /// get the ADC sample (12 bits)
  constexpr int adc(sample_type sample) { return sample&0xFFF; }
  /// get the gainId (2 bits)
  constexpr int gainId(sample_type sample) { return (sample>>12)&0x3; }
  constexpr sample_type pack(int adc, int gainId) {
    return (adc&0xFFF) | ((gainId&0x3)<<12);
  }
}


/** \class EcalMGPASample
 *  Simple container packer/unpacker for a single sample from teh MGPA electronics
 *
 *
 */
class EcalMGPASample {
 public:
  constexpr EcalMGPASample() : theSample{0} {}
  constexpr EcalMGPASample(uint16_t data) : theSample{data} {}
  constexpr EcalMGPASample(int adc,int gainId) 
      : theSample((adc&0xFFF) | ((gainId&0x3)<<12))
  {}

  /// get the raw word
  constexpr uint16_t raw() const { return theSample; }
  /// get the ADC sample (12 bits)
  constexpr int adc() const { return theSample&0xFFF; }
  /// get the gainId (2 bits)
  constexpr int gainId() const { return (theSample>>12)&0x3; }
  /// for streaming
  constexpr uint16_t operator()() const { return theSample; }
  constexpr operator uint16_t () const { return theSample; }

 private:
  uint16_t theSample;
};

std::ostream& operator<<(std::ostream&, const EcalMGPASample&);
  
#endif
