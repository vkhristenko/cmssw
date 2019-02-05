#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_reco_m0.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_common.h"

#include <iostream>

namespace hcal { namespace m0 {

#define firstSampleShift 0

/*
__device__ float get_energy(HBHEChannelInfo& info, float fc_ampl,
                            double applyContainment, float phasens, int nsamples) {
    
}*/

__device__ float get_time(HBHEChannelInfo& info, float fc_ampl,
                          HcalCalibrations const& calibs, int nsamples_to_add) {
    float time = -9999.f;
    unsigned const nsamples = info.nSamples();
    if (nsamples > 2U) {
        int const soi = info.soi();
        int ibeg = soi + firstSampleShift;
        if (ibeg < 0) ibeg = 0;
        int const iend = ibeg + nsamples_to_add;
        unsigned maxi = info.peakEnergyTS(ibeg, iend);
        if (maxi < HBHEChannelInfo::MAXSAMPLES) {
            if (!maxi) maxi = 1U;
            else if (maxi >= nsamples - 1U) maxi = nsamples - 2U;

            float t0 = info.tsEnergy(maxi - 1U);
            float maxa = info.tsEnergy(maxi);
            float t2 = info.tsEnergy(maxi + 1u);

            float mina = t0;
            if (maxa < mina) mina = maxa;
            if (t2 < mina) mina = t2;
            if (mina < 0.f) { maxa-=mina; t0-=mina; t2-=mina; }
            float wpksamp = t0 + maxa + t2;
            if (wpksamp) wpksamp = (maxa + 2.f*t2) / wpksamp;
            time = (maxi - soi)*25.f;
            //time = (maxi - soi)*25.f + timeshift_ns_hbheho(wpksamp);

            time -= calibs.timecorr();
        }
    }

    return time;
}

/// method 0 kernel
__global__ void kernel_reco(HBHEChannelInfo *vinfos, HBHERecHit *vrechits, 
                            HcalRecoParam *vparams, HcalCalibrations *vcalibs, int size) {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (idx < size) {
        auto info = vinfos[idx];
        auto params = vparams[idx];
        auto calibs = vcalibs[idx];

        int ibeg = static_cast<int>(info.soi()) + firstSampleShift;
        if (ibeg<0) ibeg=0;
        int const nsamples_to_add = params.samplesToAdd();
        double const fc_ampl = info.chargeInWindow(ibeg, ibeg + nsamples_to_add);
        bool const applyContainment = params.correctForPhaseContainment();
        float const phasens = params.correctionPhaseNS();

        // skip containment correction for now...
        float energy = info.energyInWindow(ibeg, ibeg+nsamples_to_add);
        float time = get_time(info, fc_ampl, calibs, nsamples_to_add);

        // set the values for the rec hit
        auto rechit = HBHERecHit(info.id(), energy, time, info.soiRiseTime());
    
        // set the correct rec hit
        vrechits[idx] = rechit;
    }
}

/// method 0 reconstruction
void reco(DeviceData ddata,
          HBHEChannelInfoCollection& vinfos, HBHERecHitCollection& vrechits, 
          std::vector<HcalRecoParam> const& vparams, 
          std::vector<HcalCalibrations> const& vcalibs, bool) {
    // resize the output vector
    vrechits.resize(vinfos.size());

    std::cout << "size = " << vinfos.size() << std::endl;

    // transfer to the device
    cudaMemcpy(ddata.vinfos, &(*vinfos.begin()), vinfos.size() * sizeof(HBHEChannelInfo),
        cudaMemcpyHostToDevice);
    cudaMemcpy(ddata.vparams, &(*vparams.begin()), vinfos.size() * sizeof(HcalRecoParam),
        cudaMemcpyHostToDevice);
    cudaMemcpy(ddata.vcalibs, &(*vcalibs.begin()), vinfos.size() * sizeof(HcalCalibrations),
        cudaMemcpyHostToDevice);

    // call the kernel
    int nthreadsPerBlock = 256;
    int nblocks = (vinfos.size() + nthreadsPerBlock - 1) / nthreadsPerBlock;
    kernel_reco<<<nblocks, nthreadsPerBlock>>>(ddata.vinfos, ddata.vrechits, 
        ddata.vparams, ddata.vcalibs, vinfos.size());
    cudaDeviceSynchronize();
    hcal::cuda::assert_if_error();

    // transfer back
    cudaMemcpy(&(*vrechits.begin()), ddata.vrechits, vinfos.size() * sizeof(HBHERecHit),
        cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    hcal::cuda::assert_if_error();
}

}} 
