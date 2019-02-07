#ifndef EcalPulseShapes_h
#define EcalPulseShapes_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"

struct EcalPulseShape {

public:

  static constexpr int TEMPLATESAMPLES = 12;

  constexpr EcalPulseShape() {}
  
  float pdfval[TEMPLATESAMPLES] = {0,0,0,0,
                                   0,0,0,0,
                                   0,0,0,0};
  
  constexpr float val(int isample) const { return pdfval[isample]; }

  COND_SERIALIZABLE;

};

typedef EcalCondObjectContainer<EcalPulseShape> EcalPulseShapesMap;
typedef EcalPulseShapesMap::const_iterator EcalPulseShapesMapIterator;
typedef EcalPulseShapesMap EcalPulseShapes;

#endif
