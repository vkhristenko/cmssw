#ifndef RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_m0_h
#define RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_m0_h

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParam.h"

namespace hcal { namespace m0 {

struct DeviceData {
    HBHEChannelInfo         *vinfos;
    HBHERecHit              *vrechits;
    HcalRecoParam           *vparams;
    HcalCalibrations        *vcalibs;

    void allocate(int size) {
        cudaMalloc((void**)&vinfos, size * sizeof(HBHEChannelInfo));
        cudaMalloc((void**)&vrechits, size * sizeof(HBHERecHit));
        cudaMalloc((void**)&vparams, size * sizeof(HcalRecoParam));
        cudaMalloc((void**)&vcalibs, size* sizeof(HcalCalibrations));
    }
    void free() {
        cudaFree(vinfos);
        cudaFree(vrechits);
        cudaFree(vparams);
        cudaFree(vcalibs);
    }
};

// allocate on the de

// reconstruction
void reco(DeviceData,
          HBHEChannelInfoCollection&, HBHERecHitCollection&, 
          std::vector<HcalRecoParam> const&, std::vector<HcalCalibrations> const&, bool);

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_m0_h
