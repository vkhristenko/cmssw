#ifndef RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_h
#define RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_h

#include <vector>

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

class HcalRecoParam;
class HcalCalibrations;

namespace hcal { namespace m0 {

// reconstruction
void reco(HBHEChannelInfoCollection&, HBHERecHitCollection&, 
          std::vector<HcalRecoParam> const&, std::vector<HcalCalibrations> const&, bool);

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_h
