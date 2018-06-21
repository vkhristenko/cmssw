#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_reco.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

namespace hcal { namespace m0 {

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
            float maxa = info.tsEnergy(maxI);
            float t2 = info.tsEnergy(maxi + 1u);

            float mina = t0;
            if (maxa < mina) mina = maxa;
            if (t2 < mina) mina = t2;
            if (mina < 0.f) { maxa-=mina; t0-=mina; t2-=mina; }
            float wpksamp = t0 + maxa + t2;
            if (wpksamp) wpksamp = (maxa + 2.f*t2) / wkpsamp;
            time = (maxI - soi)*25.f + timeshift_ns_hbheho(wpksamp);

            time -= calibs.timecorr();
        }
    }

    return time;
}

/// method 0 kernel
__global__ void kernel_reco(HBHEChannelInfo *vinfos, HBHERecHit *vrechits, 
                            HcalRecoParam *vparams, HcalCalibrations *vcalibs, ...) {
    // position in the vector and get the corresponding entries
    int idx = threadIdx.x;
    auto info = vinfos[idx];
    auto rechit = vrechits[idx];
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
    float time = get_time(info, fc_ampl, calibs, nSamplesToAdd);

    // set the values for the rec hit
    rechit = HBHERecHit(info.id(), energy, time, info.soiRiseTime());
}

/// method 0 reconstruction
void reco(HBHEChannelInfoCollection& vinfos, HBHERecHitCollection& vrechits, 
          std::vector<HcalRecoParam> const&, std::vector<HcalCalibrations> const&, bool) {
    // resize the output vector
    vrechits.resize(vinfos.size());

    // allocate memory on the device
    HBHEChannelInfo* d_vinfos;
    HBHERecHit *d_vrechits;
    HcalRecoParam *d_vparams;
    HcalCalibrations *d_vcalibs;
    cudaMalloc((void**)&d_vinfos, vinfos.size() * sizeof(HBHEChannelInfo));
    cudaMalloc((void**)&d_vrechits, vinfos.size() * sizeof(HBHERecHit));
    cudaMalloc((void**)&d_vparams, vinfos.size() * sizeof(HcalRecoParam));
    cudaMalloc((void**)&d_vcalibs, vinfos.size() * sizeof(HcalCalibrations));
    
    // transfer to the device
    cudaMemcpy(d_vinfos, &(*vinfos.begin()), n*sizeof(HBHEChannelInfo),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_vparams, &(*vparams.begin()), n*sizeof(HcalRecoParam),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_vcalibs, &(*vcalibs.begin()), n*sizeof(HcalCalibrations),
        cudaMemcpyHostToDevice);

    // call the kernel
    int nblocks = 1;
    int nthreadsPerBlock = vinfos.size();
    kernel_reco<<<nblocks, nthreadsPerBlock>>>();

    // transfer back
    cudaMemcpy(&(*vrechits.begin()), d_vrechits, n*sizeof(HBHERecHit),
        cudaMemcpyDeviceToHost);
}

}} 
