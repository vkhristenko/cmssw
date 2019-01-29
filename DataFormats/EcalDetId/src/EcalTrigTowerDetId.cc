#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include <ostream>

std::ostream& operator<<(std::ostream& s,const EcalTrigTowerDetId& id) {
  return s << "(EcalTT subDet " << ((id.subDet()==EcalBarrel)?("Barrel"):("Endcap")) 
	   <<  " iz " << ((id.zside()>0)?("+ "):("- ")) << " ieta " 
	   << id.ietaAbs() << " iphi " << id.iphi() << ')';
}
