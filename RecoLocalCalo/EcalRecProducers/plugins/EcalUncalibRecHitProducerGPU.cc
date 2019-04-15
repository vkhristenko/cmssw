// framework
#include "FWCore/Framework/interface/stream/EDProducer.h"
//#include "HeterogeneousCore/Producer/interface/HeterogeneousEDProducer.h"
//#include "HeterogeneousCore/Producer/interface/HeterogeneousEvent.h"

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

// algorithm specific
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include <iostream>

#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseCovariancesRcd.h"
#include "CondFormats/DataRecord/interface/EcalSamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"

#include "RecoLocalCalo/EcalRecProducers/interface/EcalPedestalsGPU.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalGainRatiosGPU.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalPulseShapesGPU.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalPulseCovariancesGPU.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalSamplesCorrelationGPU.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalTimeBiasCorrectionsGPU.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalTimeCalibConstantsGPU.h"

class EcalUncalibRecHitProducerGPU
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit EcalUncalibRecHitProducerGPU(edm::ParameterSet const& ps);
    ~EcalUncalibRecHitProducerGPU() override;

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<EBDigiCollection> digisTokenEB_;
    edm::EDGetTokenT<EEDigiCollection> digisTokenEE_;

    std::string recHitsLabelEB_, recHitsLabelEE_;
    
    // conditions handles
    edm::ESHandle<EcalPedestalsGPU> pedestalsHandle_;
    edm::ESHandle<EcalGainRatiosGPU> gainRatiosHandle_;
    edm::ESHandle<EcalPulseShapesGPU> pulseShapesHandle_;
    edm::ESHandle<EcalPulseCovariancesGPU> pulseCovariancesHandle_;
    edm::ESHandle<EcalSamplesCorrelationGPU> samplesCorrelationHandle_;
    edm::ESHandle<EcalTimeBiasCorrectionsGPU> timeBiasCorrectionsHandle_;
    edm::ESHandle<EcalTimeCalibConstantsGPU> timeCalibConstantsHandle_;

    // configuration parameters

    CUDAContextToken ctxToken_;
};

EcalUncalibRecHitProducerGPU::EcalUncalibRecHitProducerGPU(
        const edm::ParameterSet& ps) 
{
    digisTokenEB_ = consumes<EBDigiCollection>(
        ps.getUntrackedParameter<edm::InputTag>("digisLabelEB"));
    digisTokenEE_ = consumes<EEDigiCollection>(
        ps.getUntrackedParameter<edm::InputTag>("digisLabelEE"));

    recHitsLabelEB_ = ps.getUntrackedParameter<std::string>("recHitsLabelEB");
    recHitsLabelEE_ = ps.getUntrackedParameter<std::string>("recHitsLabelEE");

    EBamplitudeFitParameters = ps.getParameter<std::vector<float>>(
        "EBamplitudeFitParameters");
    EEamplitudeFitParameters = ps.getParameter<std::vector<float>>(
        "EEamplitudeFitParameters");

    produces<ecal::SoAUncalibratedRecHitCollection>(recHitsLabelEB_);
    produces<ecal::SoAUncalibratedRecHitCollection>(recHitsLabelEE_);

    //
    // configuration and physics parameters: done once
    // assume there is a single device
    // use sync copying
    //

    // amplitude fit parameters copying
    cudaMemcpy(configParameters.amplitudeFitParametersEB,
        EBamplitudeFitParameters.data(),
        EBamplitudeFitParameters.size() * sizeof(SampleVector::Scalar),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.amplitudeFitParametersEE,
        EEamplitudeFitParameters.data(),
        EEamplitudeFitParameters.size() * sizeof(SampleVector::Scalar),
        cudaMemcpyHostToDevice);

      d_data.timeFitParametersSizeEB = EBtimeFitParameters_.size();
        d_data.timeFitParametersSizeEE = EEtimeFitParameters_.size();
          d_data.timeFitLimitsFirstEB = EBtimeFitLimits_.first;
            d_data.timeFitLimitsSecondEB = EBtimeFitLimits_.second;
              d_data.timeFitLimitsFirstEE = EEtimeFitLimits_.first;
                d_data.timeFitLimitsSecondEE = EEtimeFitLimits_.second;
}

EcalUncalibRecHitProducerGPU::~EcalUncalibRecHitProducerGPU() {
    //
}

void EcalUncalibRecHitProducerGPU::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder holder) 
{
    DurationMeasurer<std::chrono::milliseconds> timer{std::string{"acquire duration"}};

    // raii
    CUDAScopedContext ctx{event.streamID(), std::move(holder)};

    setup.get<EcalPedestalsRcd>().get(pedestalsHandle_);
    setup.get<EcalGainRatiosRcd>().get(gainRatiosHandle_);
    setup.get<EcalPulseShapesRcd>().get(pulseShapesHandle_);
    setup.get<EcalPulseCovariancesRcd>().get(pulseCovariancesHandle_);
    setup.get<EcalSamplesCorrelationRcd>().get(samplesCorrelationHandle_);
    setup.get<EcalTimeBiasCorrectionsRcd>().get(timeBiasCorrectionsHandle_);
    setup.get<EcalTimeCalibConstantsRcd>().get(timeCalibConstantsHandle_);

    auto const& pedProduct = pedestalsHandle_->getProduct(ctx.stream());
    auto const& gainsProduct = gainRatiosHandle_->getProduct(ctx.stream());
    auto const& pulseShapesProduct = pulseShapesHandle_->getProduct(ctx.stream());
    auto const& pulseCovariancesProduct = pulseCovariancesHandle_->getProduct(ctx.stream());
    auto const& samplesCorrelationProduct = samplesCorrelationHandle_->getProduct(ctx.stream());
    auto const& timeBiasCorrectionsProduct = timeBiasCorrectionsHandle_->getProduct(ctx.stream());
    auto const& timeCalibConstantsProduct = timeCalibConstantsHandle_->getProduct(ctx.stream());
    
    //
    // retrieve collections
    //
    Handle<EBDigiCollection> ebDigis;
    Hanlde<EEDigicollection> eeDigis;
    evt.getByToken(digisTokenEB_, ebDigis);
    evt.getByToken(digisTokenEE_, eeDigis);

    //
    // run algorithms
    //
    ecal::multifit::entryPoint();

    // preserve token
    ctxToken_ = ctx.toToken();
}

void EcalUncalibRecHitProducerGPU::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    DurationMeasurer<std::chrono::milliseconds> timer{std::string{"produce duration"}};
    CUDAScopedContext ctx{std::move(ctxToken_)};

    // rec hits
    auto ebRecHits = std::make_unique<ecal::UncalibratedRecHit<ecal::Tag::soa>>();
    auto eeRecHits = std::make_unique<ecal::UncalibratedRecHit<ecal::Tag::soa>>();
    ebRecHits.resize(ebDigisSize);
    eeRecHits.resize(eeDigisSize);
   
    // transfer results back
    ecal::multifit::transferToHost();

    evt.put(std::move(ebRecHits), recHitsLabelEB_);
    evt.put(std::move(eeRecHits), recHitsLabelEE_);
}

DEFINE_FWK_MODULE(EcalUncalibRecHitProducerGPU);
