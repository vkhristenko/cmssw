#ifndef RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_h
#define RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_h

#include <vector>

class HBHEChannelInfoCollection;
class HBHERecHitCollection;
class HcalRecoParam;
class HcalCalibration;

namespace hcal { namespace m0 {

// reconstruction
void reco(HBHEChannelInfoCollection&, HBHERecHitCollection&, 
          std::vector<HcalRecoParam> const&, std::vector<HcalCalibrations> const&, bool);

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_h
