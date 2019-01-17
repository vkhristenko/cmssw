/*#include <algorithm>

#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCorrectionFunctions.h"

#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHitAuxSetter.h"
#include "DataFormats/METReco/interface/HcalPhase1FlagLabels.h"
#include "CondFormats/DataRecord/interface/HcalTimeSlewRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Maximum fractional error for calculating Method 0
// pulse containment correction
constexpr float PulseContainmentFractionalError = 0.002f;

SimpleHBHEPhase1Algo::SimpleHBHEPhase1Algo(
    const int firstSampleShift,
    const int samplesToAdd,
    const float phaseNS,
    const float timeShift,
    const bool correctForPhaseContainment,
    std::unique_ptr<PulseShapeFitOOTPileupCorrection> m2,
    std::unique_ptr<HcalDeterministicFit> detFit,
    std::unique_ptr<MahiFit> mahi)
    : pulseCorr_(PulseContainmentFractionalError),
      firstSampleShift_(firstSampleShift),
      samplesToAdd_(samplesToAdd),
      phaseNS_(phaseNS),
      timeShift_(timeShift),
      runnum_(0),
      corrFPC_(correctForPhaseContainment),
      psFitOOTpuCorr_(std::move(m2)),
      hltOOTpuCorr_(std::move(detFit)),
      mahiOOTpuCorr_(std::move(mahi))
{
  hcalTimeSlew_delay_ = nullptr;
}

void SimpleHBHEPhase1Algo::beginRun(const edm::Run& r,
                                    const edm::EventSetup& es)
{
    edm::ESHandle<HcalTimeSlew> delay;
    es.get<HcalTimeSlewRecord>().get("HBHE", delay);
    hcalTimeSlew_delay_ = &*delay;
  
    runnum_ = r.run();
    pulseCorr_.beginRun(es);
}

void SimpleHBHEPhase1Algo::endRun()
{
    runnum_ = 0;
}

HBHERecHit SimpleHBHEPhase1Algo::reconstruct(const HBHEChannelInfo& info,
                                             const HcalRecoParam* params,
                                             const HcalCalibrations& calibs,
                                             const bool isData)
{
    HBHERecHit rh;

    const HcalDetId channelId(info.id());

    // Calculate "Method 0" quantities
    float m0t = 0.f, m0E = 0.f;
    {
        int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int nSamplesToAdd = params ? params->samplesToAdd() : samplesToAdd_;
        const double fc_ampl = info.chargeInWindow(ibeg, ibeg + nSamplesToAdd);
        const bool applyContainment = params ? params->correctForPhaseContainment() : corrFPC_;
        const float phasens = params ? params->correctionPhaseNS() : phaseNS_;
        m0E = m0Energy(info, fc_ampl, applyContainment, phasens, nSamplesToAdd);
        m0E *= hbminusCorrectionFactor(channelId, m0E, isData);
        m0t = m0Time(info, fc_ampl, calibs, nSamplesToAdd);
    }

    // Run "Method 2"
    float m2t = 0.f, m2E = 0.f, chi2 = -1.f;
    bool useTriple = false;
    const PulseShapeFitOOTPileupCorrection* method2 = psFitOOTpuCorr_.get();
    if (method2)
    {
        psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(info.recoShape()),
                                               !info.hasTimeInfo(),info.nSamples(),hcalTimeSlew_delay_);
        // "phase1Apply" call below sets m2E, m2t, useTriple, and chi2.
        // These parameters are pased by non-const reference.
        method2->phase1Apply(info, m2E, m2t, useTriple, chi2);
        m2E *= hbminusCorrectionFactor(channelId, m2E, isData);
    }

    // Run "Method 3"
    float m3t = 0.f, m3E = 0.f;
    const HcalDeterministicFit* method3 = hltOOTpuCorr_.get();
    if (method3)
    {
        // "phase1Apply" sets m3E and m3t (pased by non-const reference)
        method3->phase1Apply(info, m3E, m3t, hcalTimeSlew_delay_);
        m3E *= hbminusCorrectionFactor(channelId, m3E, isData);
    }

    // Run Mahi
    float m4E = 0.f, m4chi2 = -1.f;
    float m4T = 0.f;
    bool m4UseTriple=false;

    const MahiFit* mahi = mahiOOTpuCorr_.get();

    if (mahi) {
      mahiOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(info.recoShape()),hcalTimeSlew_delay_);
      mahi->phase1Apply(info,m4E,m4T,m4UseTriple,m4chi2);
      m4E *= hbminusCorrectionFactor(channelId, m4E, isData);
    }

    // Finally, construct the rechit
    float rhE = m0E;
    float rht = m0t;
    float rhX = -1.f;
    if (mahi) 
    {
      rhE = m4E;
      rht = m4T;
      rhX = m4chi2;
    }
    else if (method2)
    {
        rhE = m2E;
        rht = m2t;
	rhX = chi2;
    }
    else if (method3)
    {
        rhE = m3E;
        rht = m3t;
    }
    float tdcTime = info.soiRiseTime();
    if (!HcalSpecialTimes::isSpecial(tdcTime))
        tdcTime += timeShift_;
    rh = HBHERecHit(channelId, rhE, rht, tdcTime);
    rh.setRawEnergy(m0E);
    rh.setAuxEnergy(m3E);
    rh.setChiSquared(rhX);

    // Set rechit aux words
    HBHERecHitAuxSetter::setAux(info, &rh);

    // Set some rechit flags (here, for Method 2/Mahi)
    if (useTriple || m4UseTriple)
       rh.setFlagField(1, HcalPhase1FlagLabels::HBHEPulseFitBit);

    return rh;
}

float SimpleHBHEPhase1Algo::hbminusCorrectionFactor(const HcalDetId& cell,
                                                    const float energy,
                                                    const bool isRealData) const
{
    float corr = 1.f;
    if (isRealData && runnum_ > 0)
        if (cell.subdet() == HcalBarrel)
        {
            const int ieta = cell.ieta();
            const int iphi = cell.iphi();
            corr = hbminus_special_ecorr(ieta, iphi, energy, runnum_);
        }
    return corr;
}

float SimpleHBHEPhase1Algo::m0Energy(const HBHEChannelInfo& info,
                                     const double fc_ampl,
                                     const bool applyContainmentCorrection,
                                     const double phaseNs,
                                     const int nSamplesToAdd)
{
    int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
    if (ibeg < 0)
        ibeg = 0;
    double e = info.energyInWindow(ibeg, ibeg + nSamplesToAdd);

    // Pulse containment correction
    {    
        double corrFactor = 1.0;
        if (applyContainmentCorrection)
            corrFactor = pulseCorr_.get(info.id(), nSamplesToAdd, phaseNs)->getCorrection(fc_ampl);
        e *= corrFactor;
    }

    return e;
}

float SimpleHBHEPhase1Algo::m0Time(const HBHEChannelInfo& info,
                                   const double fc_ampl,
                                   const HcalCalibrations& calibs,
                                   const int nSamplesToExamine) const
{
    float time = -9999.f; // historic value

    const unsigned nSamples = info.nSamples();
    if (nSamples > 2U)
    {
        const int soi = info.soi();
        int ibeg = soi + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int iend = ibeg + nSamplesToExamine;
        unsigned maxI = info.peakEnergyTS(ibeg, iend);
        if (maxI < HBHEChannelInfo::MAXSAMPLES)
        {
            if (!maxI)
                maxI = 1U;
            else if (maxI >= nSamples - 1U)
                maxI = nSamples - 2U;

            // The remaining code in this scope emulates
            // the historic algorithm
            float t0 = info.tsEnergy(maxI - 1U);
            float maxA = info.tsEnergy(maxI);
            float t2 = info.tsEnergy(maxI + 1U);

            // Handle negative excursions by moving "zero"
            float minA = t0;
            if (maxA < minA) minA = maxA;
            if (t2 < minA)   minA=t2;
            if (minA < 0.f) { maxA-=minA; t0-=minA; t2-=minA; }
            float wpksamp = (t0 + maxA + t2);
            if (wpksamp) wpksamp = (maxA + 2.f*t2) / wpksamp;
            time = (maxI - soi)*25.f + timeshift_ns_hbheho(wpksamp);

            // Legacy QIE8 timing correction
            time -= hcalTimeSlew_delay_->delay(std::max(1.0, fc_ampl), HcalTimeSlew::Medium);
            // Time calibration
            time -= calibs.timecorr();
        }
    }
    return time;
}

}}
*/

#include <iostream>

#include "DataFormats/HcalRecHit/interface/HBHERecHitAuxSetter.h"
#include "DataFormats/METReco/interface/HcalPhase1FlagLabels.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/MahiFit_gpu.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCorrectionFunctions.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/SimpleHBHEPhase1Algo_gpu.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_common.h"

namespace hcal { namespace mahi {

class SimpleHBHEPhase1Algo
{
public:
    // Constructor arguments:
    //
    //   firstSampleShift -- first TS w.r.t. SOI to use for "Method 0"
    //                       reconstruction.
    //
    //   samplesToAdd     -- default number of samples to add for "Method 0"
    //                       reconstruction. If, let say, SOI = 4,
    //                       firstSampleShift = -1, and samplesToAdd = 3
    //                       then the code will add time slices 3, 4, and 5.
    //
    //   phaseNS          -- default "phase" parameter for the pulse
    //                       containment correction
    //
    //   timeShift        -- time shift for QIE11 TDC times
    //
    //   correctForPhaseContainment -- default switch for applying pulse
    //                                 containment correction for "Method 0"
    //
    //   m2               -- "Method 2" object
    //
    //   detFit           -- "Method 3" (a.k.a. "deterministic fit") object
    //
    __device__
    SimpleHBHEPhase1Algo(int firstSampleShift,
                         int samplesToAdd,
                         float phaseNS,
                         float timeShift,
                         bool correctForPhaseContainment);
    __device__
    SimpleHBHEPhase1Algo();

    __device__
    inline bool isConfigurable() const {return false;}

    __device__
    HBHERecHit reconstruct(const HBHEChannelInfo& info,
                           const HcalRecoParam* params,
                           const HcalCalibrations& calibs,
                           bool isRealData,
                           float const* pshape /* this is added */);

    // Basic accessors
    __device__
    inline int getFirstSampleShift() const {return firstSampleShift_;}
    __device__
    inline int getSamplesToAdd() const {return samplesToAdd_;}
    __device__
    inline float getPhaseNS() const {return phaseNS_;}
    __device__
    inline float getTimeShift() const {return timeShift_;}
    __device__
    inline bool isCorrectingForPhaseContainment() const {return corrFPC_;}

    // TODO: resolve this
//    const HcalTimeSlew* hcalTimeSlew_delay_;

protected:
    // Special HB- correction
    __device__
    float hbminusCorrectionFactor(const HcalDetId& cell,
                                  float energy, bool isRealData) const;

    // "Method 0" rechit energy. Calls a non-const member of
    // HcalPulseContainmentManager, so no const qualifier here.
    // HB- correction is not applied inside this function.
    __device__
    float m0Energy(const HBHEChannelInfo& info,
                   double reconstructedCharge,
                   bool applyContainmentCorrection,
                   double phaseNS, int nSamplesToAdd);

    // "Method 0" rechit timing (original low-pileup QIE8 algorithm)
    __device__
    float m0Time(const HBHEChannelInfo& info,
                 double reconstructedCharge,
                 const HcalCalibrations& calibs,
                 int nSamplesToExamine) const;
private:
    // TODO:resolve this
 //   HcalPulseContainmentManager pulseCorr_;

    int firstSampleShift_;
    int samplesToAdd_;
    float phaseNS_;
    float timeShift_;
    bool corrFPC_;
    int runnum_{999999};
};

constexpr float PulseContainmentFractionalError = 0.002f;

__device__
SimpleHBHEPhase1Algo::SimpleHBHEPhase1Algo(
    const int firstSampleShift,
    const int samplesToAdd,
    const float phaseNS,
    const float timeShift,
    const bool correctForPhaseContainment)
    : firstSampleShift_(firstSampleShift),
      samplesToAdd_(samplesToAdd),
      phaseNS_(phaseNS),
      timeShift_(timeShift),
      corrFPC_(correctForPhaseContainment)
{}

__device__
SimpleHBHEPhase1Algo::SimpleHBHEPhase1Algo() 
    : firstSampleShift_{0}, samplesToAdd_{2},
      phaseNS_{6.0},
      timeShift_{0.0},
      corrFPC_{true}
{}

__device__
HBHERecHit SimpleHBHEPhase1Algo::reconstruct(const HBHEChannelInfo& info,
                                             const HcalRecoParam* params,
                                             const HcalCalibrations& calibs,
                                             const bool isData,
                                             float const* pshape) {
    HBHERecHit rh;

    const HcalDetId channelId(info.id());

    // Calculate "Method 0" quantities
    float m0t = 0.f, m0E = 0.f;
    {
        int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int nSamplesToAdd = params ? params->samplesToAdd() : samplesToAdd_;
        const double fc_ampl = info.chargeInWindow(ibeg, ibeg + nSamplesToAdd);
        const bool applyContainment = params ? params->correctForPhaseContainment() : corrFPC_;
        const float phasens = params ? params->correctionPhaseNS() : phaseNS_;
        m0E = m0Energy(info, fc_ampl, applyContainment, phasens, nSamplesToAdd);
        m0E *= hbminusCorrectionFactor(channelId, m0E, isData);
        m0t = m0Time(info, fc_ampl, calibs, nSamplesToAdd);
    }

    // Run "Method 2"
    float m2t = 0.f, m2E = 0.f, chi2 = -1.f;
    bool useTriple = false;
    /*
    const PulseShapeFitOOTPileupCorrection* method2 = psFitOOTpuCorr_.get();
    if (method2)
    {
        psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(info.recoShape()),
                                               !info.hasTimeInfo(),info.nSamples(),hcalTimeSlew_delay_);
        // "phase1Apply" call below sets m2E, m2t, useTriple, and chi2.
        // These parameters are pased by non-const reference.
        method2->phase1Apply(info, m2E, m2t, useTriple, chi2);
        m2E *= hbminusCorrectionFactor(channelId, m2E, isData);
    }*/

    // TODO: resolve this
    // Run "Method 3"
    float m3t = 0.f, m3E = 0.f;
    /*
    const HcalDeterministicFit* method3 = hltOOTpuCorr_.get();
    if (method3)
    {
        // "phase1Apply" sets m3E and m3t (pased by non-const reference)
        method3->phase1Apply(info, m3E, m3t, hcalTimeSlew_delay_);
        m3E *= hbminusCorrectionFactor(channelId, m3E, isData);
    }*/

    // Run Mahi
    float m4E = 0.f, m4chi2 = -1.f;
    float m4T = 0.f;
    bool m4UseTriple=false;

    MahiFit mfit{pshape};
    mfit.phase1Apply(info, m4E, m4T, m4UseTriple, m4chi2);
    m4E *= hbminusCorrectionFactor(channelId, m4E, isData);

    /*
    const MahiFit* mahi = mahiOOTpuCorr_.get();

    if (mahi) {
      mahiOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(info.recoShape()),hcalTimeSlew_delay_);
      mahi->phase1Apply(info,m4E,m4T,m4UseTriple,m4chi2);
      m4E *= hbminusCorrectionFactor(channelId, m4E, isData);
    }
    */

    // Finally, construct the rechit
    float rhE = m0E;
    float rht = m0t;
    float rhX = -1.f;
      rhE = m4E;
      rht = m4T;
      rhX = m4chi2;
    float tdcTime = info.soiRiseTime();
    if (!HcalSpecialTimes::isSpecial(tdcTime))
        tdcTime += timeShift_;
    rh = HBHERecHit(channelId, rhE, rht, tdcTime);
    rh.setRawEnergy(m0E);
    rh.setAuxEnergy(m3E);
    rh.setChiSquared(rhX);

    // Set rechit aux words
    HBHERecHitAuxSetter::setAux(info, &rh);

    // Set some rechit flags (here, for Method 2/Mahi)
    if (useTriple || m4UseTriple)
       rh.setFlagField(1, HcalPhase1FlagLabels::HBHEPulseFitBit);

    return rh;
}

__device__
float SimpleHBHEPhase1Algo::hbminusCorrectionFactor(const HcalDetId& cell,
                                                    const float energy,
                                                    const bool isRealData) const
{
    float corr = 1.f;
    if (isRealData && runnum_ > 0)
        if (cell.subdet() == HcalBarrel)
        {
            const int ieta = cell.ieta();
            const int iphi = cell.iphi();
            corr = hbminus_special_ecorr(ieta, iphi, energy, runnum_);
        }
    return corr;
}

__device__
float SimpleHBHEPhase1Algo::m0Energy(const HBHEChannelInfo& info,
                                     const double fc_ampl,
                                     const bool applyContainmentCorrection,
                                     const double phaseNs,
                                     const int nSamplesToAdd)
{
    int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
    if (ibeg < 0)
        ibeg = 0;
    double e = info.energyInWindow(ibeg, ibeg + nSamplesToAdd);

    /* TODO: resolve this correction
    // Pulse containment correction
    {    
        double corrFactor = 1.0;
        if (applyContainmentCorrection)
            corrFactor = pulseCorr_.get(info.id(), nSamplesToAdd, phaseNs)->getCorrection(fc_ampl);
        e *= corrFactor;
    }*/

    return e;
}

__device__
float SimpleHBHEPhase1Algo::m0Time(const HBHEChannelInfo& info,
                                   const double fc_ampl,
                                   const HcalCalibrations& calibs,
                                   const int nSamplesToExamine) const
{
    float time = -9999.f; // historic value

    const unsigned nSamples = info.nSamples();
    if (nSamples > 2U)
    {
        const int soi = info.soi();
        int ibeg = soi + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int iend = ibeg + nSamplesToExamine;
        unsigned maxI = info.peakEnergyTS(ibeg, iend);
        if (maxI < HBHEChannelInfo::MAXSAMPLES)
        {
            if (!maxI)
                maxI = 1U;
            else if (maxI >= nSamples - 1U)
                maxI = nSamples - 2U;

            // The remaining code in this scope emulates
            // the historic algorithm
            float t0 = info.tsEnergy(maxI - 1U);
            float maxA = info.tsEnergy(maxI);
            float t2 = info.tsEnergy(maxI + 1U);

            // Handle negative excursions by moving "zero"
            float minA = t0;
            if (maxA < minA) minA = maxA;
            if (t2 < minA)   minA=t2;
            if (minA < 0.f) { maxA-=minA; t0-=minA; t2-=minA; }
            float wpksamp = (t0 + maxA + t2);
            if (wpksamp) wpksamp = (maxA + 2.f*t2) / wpksamp;
            time = (maxI - soi)*25.f + timeshift_ns_hbheho(wpksamp);

            // Legacy QIE8 timing correction
            // TODO: resolve this correction
            //time -= hcalTimeSlew_delay_->delay(std::max(1.0, fc_ampl), HcalTimeSlew::Medium);
            // Time calibration
            time -= calibs.timecorr();
        }
    }

    return time;
}

/// mahi/multifit kernel
__global__ 
void kernel_reconstruct(HBHEChannelInfo *vinfos, HBHERecHit *vrechits,
                            HcalRecoParam *vparams, HcalCalibrations *vcalibs,
                            int* hashes, float* psdata, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx < size) {
        auto info = vinfos[idx];
        auto params = vparams[idx];
        auto calibs = vcalibs[idx];
        auto *pshape = psdata + hashes[info.recoShape()] * 256;

        SimpleHBHEPhase1Algo algo{};
        auto rh = algo.reconstruct(info, &params, calibs,  
            /* TODO: drag this boolean from the host */ true,
            pshape);

        vrechits[idx] = rh;
    }
}

void reconstruct(DeviceData ddata,
          HBHEChannelInfoCollection& vinfos, HBHERecHitCollection& vrechits, 
          std::vector<HcalRecoParam> const& vparams, 
          std::vector<HcalCalibrations> const& vcalibs, 
          PulseShapeData &psdata, bool) {
    // resize the output vector
    vrechits.resize(vinfos.size());

    std::cout << "size = " << vinfos.size() << std::endl;

    // transfer to the device
    cudaMemcpy(ddata.vinfos, &(*vinfos.begin()), 
        vinfos.size() * sizeof(HBHEChannelInfo),
        cudaMemcpyHostToDevice);
    cudaMemcpy(ddata.vparams, &(*vparams.begin()), 
        vinfos.size() * sizeof(HcalRecoParam),
        cudaMemcpyHostToDevice);
    cudaMemcpy(ddata.vcalibs, &(*vcalibs.begin()), 
        vinfos.size() * sizeof(HcalCalibrations),
        cudaMemcpyHostToDevice);

    // call the kernel
    int nthreadsPerBlock = 256;
    int nblocks = (vinfos.size() + nthreadsPerBlock - 1) / nthreadsPerBlock;
    kernel_reconstruct<<<nblocks, nthreadsPerBlock>>>(ddata.vinfos, ddata.vrechits,
        ddata.vparams, ddata.vcalibs, psdata.hashes, psdata.data, vinfos.size());
        cudaDeviceSynchronize();
        hcal::cuda::assert_if_error();

    // transfer back
    cudaMemcpy(&(*vrechits.begin()), ddata.vrechits, 
        vinfos.size() * sizeof(HBHERecHit),
        cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        hcal::cuda::assert_if_error();
}

}}
