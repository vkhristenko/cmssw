#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_reco_mahi.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_common.h"

#include <iostream>

#define fullTSofInterest_ 8
#define timeSigmaSiPM_ 2.5
#define timeSigmaHPD_ 5.0
#define ts4Thresh_ 0.0
#define chiSqSwitch_ 15.0

namespace hcal { namespace mahi {

#define firstSampleShift 0

/*
__device__ float get_energy(HBHEChannelInfo& info, float fc_ampl,
                            double applyContainment, float phasens, int nsamples) {
    
}*/

__device__ void do_fit(RecValues &values, int ) {
    // do nothing
    return;
}

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

        RecValues recValues;

        // workspace
        Workspace ws;
        ws.tsSize = info.nSamples();
        ws.tsOffset = info.soi();
        ws.fullTSOffset = fullTSofInterest_ - ws.tsOffset;

        // TODO: understand why we need this:
        bool useTriple = false;

        // 1 sigma time constraint
        if (info.hasTimeInfo()) ws.dt = timeSigmaSiPM_;
        else ws.dt = timeSigmaHPD_;

        // average pedestal width (for covariance matrix constraint)
        float pedVal = 0.25*( info.tsPedestalWidth(0)*info.tsPedestalWidth(0)+
            info.tsPedestalWidth(1)*info.tsPedestalWidth(1)+
            info.tsPedestalWidth(2)*info.tsPedestalWidth(2)+
            info.tsPedestalWidth(3)*info.tsPedestalWidth(3) );

        // 
        ws.pedConstraint = pedVal * SampleMatrix::Ones(ws.tsSize, ws.tsSize);
        ws.amplitudes.resize(ws.tsSize);
        ws.noiseTerms.resize(ws.tsSize);

        // per ts
        double tstot = 0, tstrig = 0;
        for (unsigned int iTS=0; iTS<ws.tsSize; ++iTS) {
            double charge = info.tsRawCharge(iTS);
            double ped = info.tsPedestal(iTS);
            ws.amplitudes.coeffRef(iTS) = charge - ped;

            // adc granularity
            double noiseADC = (1. / sqrt(12)) * info.tsDFcPerADC(iTS);

            // photostatistics
            double noisePhoto = 0;
            auto tmp = charge - ped;
            if (tmp > info.tsPedestalWidth(iTS)) {
                noisePhoto = sqrt(tmp * info.fcByPE());
            }

            // electronic pedestal
            double pedWidth = info.tsPedestalWidth(iTS);

            // total uncertainty from all sources
            ws.noiseTerms.coeffRef(iTS) = noiseADC*noiseADC + noisePhoto*noisePhoto + 
                pedWidth*pedWidth;

            tstot += (charge - ped) * info.tsGain(0);
            if (iTS == ws.tsOffset) 
                tstrig += (charge - ped) * info.tsGain(0);
        }

        //
        if (tstrig >= ts4Thresh_ and tstot>0) {
            useTriple = false;
            
            // only do the prefit with 1 pulse if ichisq threshold is positive
            if (chiSqSwitch_>0) {
                do_fit(recValues, 1);
                if (recValues.chi2 > chiSqSwitch_) {
                    do_fit(recValues, 0); // nbx=0 means to use configured bxs
                    useTriple = true;
                }
            } else {
                do_fit(recValues, 0);
                useTriple = true;
            }
        } else {
            recValues.energy = 0.;
            recValues.time = -9999.;
            recValues.chi2 = -9999.;
        }
        
        // TODO  rewrite these guys
        float energy = recValues.energy * info.tsGain(0);
        float time = recValues.time;
//        float chi2 = recValues.chi2;

        // set the values for the rec hit and put it into the correct array cell
        auto rechit = HBHERecHit(info.id(), energy, time, info.soiRiseTime());
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
