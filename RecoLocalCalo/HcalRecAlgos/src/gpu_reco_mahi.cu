#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_reco_mahi.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_common.h"

#include <iostream>

#define fullTSofInterest_ 8
#define timeSigmaSiPM_ 2.5
#define timeSigmaHPD_ 5.0
#define ts4Thresh_ 0.0
#define chiSqSwitch_ 15.0
#define firstSampleShift 0
#define meanTime_ 0.
#define pedestalBX_ 100
#define nnlsThresh_ 1e-11
#define nMaxItersNNLS_ 500
#define nMaxItersMin_ 500
#define deltaChiSqThresh_ 1e-3
#define timeLimit_ 12.5
#define dynamicPed_ true

__constant__ int const activeBXs_[] = {-1, 0, 1};
static int const bxSizeConf_ = 3;
static int const bxOffsetConf_ = 1;

namespace hcal { namespace mahi {

/*
__device__ float get_energy(HBHEChannelInfo& info, float fc_ampl,
                            double applyContainment, float phasens, int nsamples) {
    
}*/

__device__ void update_pulse_shape(Workspace &ws, FullSampleVector &pshape,
                                   FullSampleVector &pderiv, FullSampleMatrix &pcov,
                                   double itq) {
    //
    float t0 = meanTime_;

    // TODO: to implement
    /*
    if (applyTimeSlew_) {
        if (itq <= 1.0) 
            t0 += tsDelay1GeV_;
        else 
            t0 += hcalTimeSlewDelay_->delay(itq, slewFlavor);
    }*/

    //
    for (int i=0; i<MaxSVSize; ++i) {
        ws.pulseN[i] = 0;
        ws.pulseM[i] = 0;
        ws.pulseP[i] = 0;
    }

    //
    double const xx[4] = {t0, 1.0, 0.0, 3};
    double const xxm[4] = {-ws.dt + t0, 1.0, 0.0, 3};
    double const xxp[4] = {ws.dt + t0, 1.0, 0.0, 3};

    // 
    ws.functor.singlePulseShapeFunc(xx);
    ws.functor.getPulseShape(ws.pulseN);
    ws.functor.singlePulseShapeFunc(xxm);
    ws.functor.getPulseShape(ws.pulseM);
    ws.functor.singlePulseShapeFunc(xxp);
    ws.functor.getPulseShape(ws.pulseP);

    //
    int delta = ws.tsOffset == 3 ? 1 : 0;
    for (unsigned int its=ws.fullTSOffset; its<ws.fullTSOffset + ws.tsSize; ++its) {
        pshape.coeffRef(its) = ws.pulseN[its - ws.fullTSOffset + delta];
        pderiv.coeffRef(its) = 0.5*(ws.pulseM[its - ws.fullTSOffset+delta] + 
            ws.pulseP[its - ws.fullTSOffset + delta])/(2 * ws.dt);
        ws.pulseM[its - ws.fullTSOffset] -= ws.pulseN[its - ws.fullTSOffset];
        ws.pulseP[its - ws.fullTSOffset] -= ws.pulseN[its - ws.fullTSOffset];
    }

    //
    for (unsigned int its=ws.fullTSOffset; its<ws.fullTSOffset+ws.tsSize; ++its) {
        for (unsigned int jts=ws.fullTSOffset; jts<its+1; ++jts) {
            double tmp = 0.5*(ws.pulseP[its - ws.fullTSOffset + delta]*
                ws.pulseP[jts - ws.fullTSOffset + delta] + 
                ws.pulseM[its - ws.fullTSOffset + delta] * 
                ws.pulseM[jts - ws.fullTSOffset + delta]);

            pcov(its, jts) += tmp;
            pcov(jts, its) += tmp;
        }
    }
}

__device__ void update_covariance(Workspace &ws) {
    //
    ws.invCovMat.resize(ws.tsSize, ws.tsSize);
    ws.invCovMat = ws.noiseTerms.asDiagonal();
    ws.invCovMat = ws.pedConstraint;

    // 
    for (unsigned int ibx=0; ibx<ws.nPulseTot; ++ibx) {
        if (ws.ampVec.coeff(ibx) == 0) continue;
        int offset = ws.bxs.coeff(ibx);
        if (offset == pedestalBX_) continue;
        else {
            ws.invCovMat += ws.ampVec.coeff(ibx)*ws.ampVec.coeff(ibx)*
                ws.pulseCovArray[offset + ws.bxOffset].block(
                    ws.fullTSOffset - offset, ws.fullTSOffset - offset,
                    ws.tsSize, ws.tsSize);
        }
    }

    // recompute Cholesky decomposition
    // TODO
    ws.covDecomp.compute(ws.invCovMat);
}

__device__ void nnls_unconstrain_parameter(Workspace &ws, Eigen::Index idxp) {
    ws.aTaMat.col(ws.nP).swap(ws.aTaMat.col(idxp));
    ws.aTaMat.row(ws.nP).swap(ws.aTaMat.row(idxp));
    ws.pulseMat.col(ws.nP).swap(ws.pulseMat.col(idxp));
    
    // 
    Eigen::numext::swap(ws.aTbVec.coeffRef(ws.nP), ws.aTbVec.coeffRef(idxp));
    Eigen::numext::swap(ws.ampVec.coeffRef(ws.nP), ws.ampVec.coeffRef(idxp));
    Eigen::numext::swap(ws.bxs.coeffRef(ws.nP), ws.bxs.coeffRef(idxp));
    ++ws.nP;
}

__device__ void nnls_constrain_parameter(Workspace &ws, Eigen::Index minratioidx) {
    ws.aTaMat.col(ws.nP-1).swap(ws.aTaMat.col(minratioidx));
    ws.aTaMat.row(ws.nP-1).swap(ws.aTaMat.row(minratioidx));
    ws.pulseMat.col(ws.nP-1).swap(ws.pulseMat.col(minratioidx));

    //
    Eigen::numext::swap(ws.aTbVec.coeffRef(ws.nP-1), ws.aTbVec.coeffRef(minratioidx));
    Eigen::numext::swap(ws.ampVec.coeffRef(ws.nP-1), ws.ampVec.coeffRef(minratioidx));
    Eigen::numext::swap(ws.bxs.coeffRef(ws.nP-1), ws.bxs.coeffRef(minratioidx));
    --ws.nP;
}

/*
__device__ void solve_submatrix(Workspace &ws) {
    switch (ws.nP) {
    case 10:
        {
            auto tmp = ws.aTaMat;

            aTbVec
                ampvecpermtest
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    case 9:
        {
        }
        break;
    default:
        return;
    }
}
*/

__device__ void nnls(Workspace &ws) {
    unsigned int const npulse = ws.nPulseTot;
    for (unsigned int ibx=0; ibx<npulse; ++ibx) {
        int offset = ws.bxs.coeff(ibx);
        if (offset == pedestalBX_)
            ws.pulseMat.col(ibx) = SampleVector::Ones(ws.tsSize);
        else 
            ws.pulseMat.col(ibx) = ws.pulseShapeArray[offset + ws.bxOffset].segment(
                ws.fullTSOffset - offset, ws.tsSize);
    }

    // TODO
    ws.invcovp = ws.covDecomp.matrixL().solve(ws.pulseMat);
    ws.aTaMat = ws.invcovp.transpose().lazyProduct(ws.invcovp);
    auto tmp = ws.covDecomp.matrixL().solve(ws.amplitudes);
    //auto tmp = SampleVector::Ones(ws.tsSize);
    ws.aTbVec = ws.invcovp.transpose().lazyProduct(tmp);

    //
    int iter = 0;
    double wmax = 0.0;
    double threshold = nnlsThresh_;
    ws.nP = 0;
    Eigen::Index idxwmax = 0;
    while (true) {
        if (iter>0 or ws.nP==0) {
            if (ws.nP == std::min(npulse, ws.tsSize)) break;

            unsigned int const nactive = npulse - ws.nP;
            ws.updateWork = ws.aTbVec - ws.aTaMat*ws.ampVec;
            Eigen::Index idxwmaxprev = idxwmax;
            double wmaxprev = wmax;
            wmax = ws.updateWork.tail(nactive).maxCoeff(&idxwmax);
            if (wmax < threshold or (idxwmax == idxwmaxprev and wmax == wmaxprev))
                break;
            if (iter >= nMaxItersNNLS_)
                break;

            //
            Eigen::Index idxp = ws.nP + idxwmax;
            nnls_unconstrain_parameter(ws, idxp);
        }

        //
        while (true) {
            if (ws.nP == 0) break;
            ws.ampvecpermtest = ws.ampVec;
    
            // solveSubmatrix()
            // TODO
            //solve_submatrix(ws.aTaMat, ws.aTbVec, ws.ampvecpermtest, ws.nP);
            ws.ampvecpermtest.head(ws.nP) = ws.aTaMat.topLeftCorner(ws.nP, ws.nP).ldlt().solve(ws.aTbVec.head(ws.nP));

            // check solution
            bool positive = true;
            for (unsigned int i=0; i<ws.nP; i++)
                positive &= ws.ampvecpermtest(i) > 0;
            if (positive) {
                ws.ampVec.head(ws.nP) = ws.ampvecpermtest.head(ws.nP);
                break;
            }

            // update parameter vector
            Eigen::Index minratioidx = 0;
            double minratio = std::numeric_limits<double>::max();
            for (unsigned int ipulse=0; ipulse<ws.nP; ++ipulse) {
                if (ws.ampvecpermtest.coeff(ipulse) <= 0.) {
                    double const c_ampvec = ws.ampVec.coeff(ipulse);
                    double const ratio = c_ampvec / 
                        (c_ampvec - ws.ampvecpermtest.coeff(ipulse));
                    if (ratio < minratio) {
                        minratio = ratio;
                        minratioidx = ipulse;
                    }
                }
            }

            //
            ws.ampVec.head(ws.nP) += minratio * (ws.ampvecpermtest.head(ws.nP) - 
                ws.ampVec.head(ws.nP));

            ws.ampVec.coeffRef(minratioidx) = 0.;
            nnls_constrain_parameter(ws, minratioidx);
        }

        ++iter;
        // adaptive convergence threshold to avoid infinite loops but still
        // ensure best value is used
        if (iter % 10 == 0)
            threshold *= 10.;
    }
}

__device__ void one_pulse_minimize(Workspace &ws) {
    //ws.invcovp = ws.covDecomp.matrixL().solve(ws.pulseMat);

    // TODO 
    auto aTamatval = ws.invcovp.transpose()*ws.invcovp;
    auto tmp = ws.covDecomp.matrixL().solve(ws.amplitudes);
    auto aTbvecval = ws.invcovp.transpose()*tmp;
    //SingleVector aTbvecval = SingleVector::Ones(ws.tsSize);
    ws.ampVec.coeffRef(0) = std::max(0., aTbvecval.coeff(0)/aTamatval.coeff(0));
}

// TODO
__device__ double calculate_chi2(Workspace &ws) {
    return ws.covDecomp.matrixL().solve(ws.pulseMat*ws.ampVec - ws.amplitudes).squaredNorm();
}

__device__ double minimize(Workspace &ws) {
    double oldchi2 = 9999;
    double chi2 = oldchi2;

    for (int it=1; it<nMaxItersMin_; ++it) {
        update_covariance(ws);

        // 
        if (ws.nPulseTot > 1)
            nnls(ws);
        else 
            one_pulse_minimize(ws);

        // 
        double newchi2 = calculate_chi2(ws);
        double deltachi2 = newchi2 - chi2;
        if (newchi2 == oldchi2 and newchi2<chi2)
            break;
        oldchi2 = chi2;
        chi2 = newchi2;
        if (std::abs(deltachi2) < deltaChiSqThresh_)
            break;
    }

    return chi2;
}

__device__ double calc_arrival_time(Workspace &ws) {
    ws.pulseDerivMat.resize(ws.tsSize, ws.nPulseTot);
    int itIndex = 0;
    for (unsigned int ibx=0; ibx<ws.nPulseTot; ++ibx) {
        int offset = ws.bxs.coeff(ibx);
        if (offset == 0) itIndex =ibx;
        if (offset == pedestalBX_) 
            ws.pulseDerivMat.col(ibx) = SampleVector::Zero(ws.tsSize);
        else 
            ws.pulseDerivMat.col(ibx) = ws.pulseDerivArray[
                offset + ws.bxOffset].segment(ws.fullTSOffset - offset, 
                ws.tsSize);
    }

    // TODO: this needs to be implemented!
    PulseVector solution = ws.pulseDerivMat.colPivHouseholderQr().solve(ws.residuals);
    //PulseVector solution = PulseVector::Ones(ws.tsSize);
    float t = solution.coeff(itIndex) / ws.ampVec.coeff(itIndex);
    t = (t > timeLimit_) ? timeLimit_ :
        ((t < -timeLimit_) ? -timeLimit_ : t);

    return t;
}

__device__ void do_fit(Workspace &ws, RecValues &values, int nbx) {
    unsigned int bxSize = 1;

    //
    if (nbx==1) {
         ws.bxOffset = 0;
        ws.bxs.resize(bxSize);
        ws.bxs.coeffRef(0) = 0;
    } else {
        bxSize = bxSizeConf_;
        ws.bxOffset = bxOffsetConf_;

        ws.bxs.resize(bxSize);
        for (unsigned int ibx=0; ibx<bxSize; ++ibx) 
            ws.bxs.coeffRef(ibx) = activeBXs_[ibx];
    }

    //
    ws.nPulseTot = bxSize;
    if (dynamicPed_) {
        ws.nPulseTot++;
        ws.bxs.resize(ws.nPulseTot);
        ws.bxs[ws.nPulseTot-1] = pedestalBX_;
    }

    //
    ws.pulseMat.resize(ws.tsSize, ws.nPulseTot);
    ws.ampVec = PulseVector::Zero(ws.nPulseTot);
    ws.errVec = PulseVector::Zero(ws.nPulseTot);
    int offset = 0;

    //
    for (unsigned int ibx=0; ibx<ws.nPulseTot; ++ibx) {
        offset = ws.bxs.coeff(ibx);

        //
        ws.pulseShapeArray[ibx] = FullSampleVector::Zero(MaxFSVSize);
        ws.pulseDerivArray[ibx] = FullSampleVector::Zero(MaxFSVSize);
        ws.pulseCovArray[ibx] = FullSampleMatrix::Constant(0);

        if (offset == pedestalBX_) 
            ws.ampVec.coeffRef(ibx) = 0;
        else {
            update_pulse_shape(ws, 
                               ws.pulseShapeArray[ibx],
                               ws.pulseDerivArray[ibx],
                               ws.pulseCovArray[ibx],
                               ws.amplitudes.coeff(ws.tsOffset + offset));
            ws.ampVec.coeffRef(ibx) = 0;
            ws.pulseMat.col(ibx) = ws.pulseShapeArray[ibx].segment(
                ws.fullTSOffset - offset, ws.tsSize);
        }
    }

    // 
    ws.pulseMat.col(ws.nPulseTot - 1) = SampleVector::Ones(ws.tsSize);
    ws.aTaMat.resize(ws.nPulseTot, ws.nPulseTot);
    ws.aTbVec.resize(ws.nPulseTot);

    //
    // minimization
    //
    double chi2 = minimize(ws);

    //
    // compute residuals
    //
    ws.residuals = ws.pulseMat*ws.ampVec - ws.amplitudes;

    //
    bool foundintime = false;
    unsigned int ipulseintime = 0;
    for (unsigned int ibx=0; ibx<ws.nPulseTot; ++ibx) {
        if (ws.bxs.coeff(ibx)==0) {
            ipulseintime = ibx;
            foundintime = true;
        }
    }

    //
    if (foundintime) {
        values.energy = ws.ampVec.coeff(ipulseintime);
        if (values.energy != 0) {
            double arrival_time = calc_arrival_time(ws);
            values.energy = arrival_time;
        } else 
            values.time = -9999;
        values.chi2 = chi2;
    }
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
                            HcalRecoParam *vparams, HcalCalibrations *vcalibs, 
                            int* hashes, float* psdata, int size) {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (idx < size) {
        auto info = vinfos[idx];
        auto params = vparams[idx];
        auto calibs = vcalibs[idx];

        // reconstructed values
        RecValues recValues;

        // workspace
        auto *pshape = psdata + hashes[info.recoShape()];
        Workspace ws;
        ws.functor.assign(pshape, false, false, false, 1, 0, 0, 10);
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
                do_fit(ws, recValues, 1);
                if (recValues.chi2 > chiSqSwitch_) {
                    do_fit(ws, recValues, 0); // nbx=0 means to use configured bxs
                    useTriple = true;
                }
            } else {
                do_fit(ws, recValues, 0);
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
        float chi2 = recValues.chi2;

        // set the values for the rec hit and put it into the correct array cell
        auto rechit = HBHERecHit(info.id(), energy, time, info.soiRiseTime());
        vrechits[idx] = rechit;
    }
}

/// method 0 reconstruction
void reco(DeviceData ddata,
          HBHEChannelInfoCollection& vinfos, HBHERecHitCollection& vrechits, 
          std::vector<HcalRecoParam> const& vparams, 
          std::vector<HcalCalibrations> const& vcalibs, 
          PulseShapeData &psdata, bool) {
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
        ddata.vparams, ddata.vcalibs, psdata.hashes, psdata.data, vinfos.size());
    cudaDeviceSynchronize();
    hcal::cuda::assert_if_error();

    // transfer back
    cudaMemcpy(&(*vrechits.begin()), ddata.vrechits, vinfos.size() * sizeof(HBHERecHit),
        cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    hcal::cuda::assert_if_error();
}

}} 
