#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include<iostream>

std::ostream& operator<<(std::ostream& s, const EcalMGPASample& samp) {
  s << "ADC=" << samp.adc() << ", gainId=" << samp.gainId();
  return s;
}
