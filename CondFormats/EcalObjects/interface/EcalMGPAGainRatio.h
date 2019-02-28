#ifndef CondFormats_EcalObjects_EcalMGPAGainRatio_H
#define CondFormats_EcalObjects_EcalMGPAGainRatio_H
/**
 * Author: Shahram Rahatlou, University of Rome & INFN
 * Created: 22 Feb 2006
 * $Id: $
 **/


#include "CondFormats/Serialization/interface/Serializable.h"

#include <iostream>

class EcalMGPAGainRatio {
  public:
    constexpr EcalMGPAGainRatio() : gain12Over6_{2.}, gain6Over1_{6.} {}
    constexpr EcalMGPAGainRatio(const EcalMGPAGainRatio & ratio) 
        : gain12Over6_{ratio.gain12Over6_}, gain6Over1_{ratio.gain6Over1_} {}
        
    constexpr float gain12Over6() const { return gain12Over6_; }
    constexpr float gain6Over1() const { return gain6Over1_; }

    constexpr void setGain12Over6(const float& g) { gain12Over6_ = g; }
    constexpr void setGain6Over1(const float& g)  { gain6Over1_ = g; }

    void print(std::ostream& s) const { s << "gain 12/6: " << gain12Over6_ << " gain 6/1: " << gain6Over1_; }

    constexpr EcalMGPAGainRatio& operator=(const EcalMGPAGainRatio& rhs) {
      gain12Over6_ = rhs.gain12Over6_;
      gain6Over1_ = rhs.gain6Over1_;
      return *this;
    }

  private:
    float gain12Over6_;
    float gain6Over1_;

  COND_SERIALIZABLE;
};
#endif
