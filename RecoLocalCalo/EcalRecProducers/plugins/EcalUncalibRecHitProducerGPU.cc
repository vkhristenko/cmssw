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

    produces<ecal::SoAUncalibratedRecHitCollection>(recHitsLabelEB_);
    produces<ecal::SoAUncalibratedRecHitCollection>(recHitsLabelEE_);
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
    std::cout << "acquire\n";

    ctxToken_ = ctx.toToken();
}

void EcalUncalibRecHitProducerGPU::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    DurationMeasurer<std::chrono::milliseconds> timer{std::string{"produce duration"}};
    CUDAScopedContext ctx{std::move(ctxToken_)};

    std::cout << "produce\n";
}

/*
void EcalUncalibRecHitProducerGPU::produce(edm::Event& e, const edm::EventSetup&) {
    // input products
    edm::Handle<EBDigiCollection> digisEB;
    edm::Handle<EEDigiCollection> digisEE;
    e.getByToken(digisTokenEB_, digisEB);
    e.getByToken(digisTokenEE_, digisEE);

    // output products
    auto recHitsEB =
        std::make_unique<ecal::UncalibratedRecHit<ecal::Tag::soa>>();
    auto recHitsEE =
        std::make_unique<ecal::UncalibratedRecHit<ecal::Tag::soa>>();

    e.put(std::move(recHitsEB), recHitsLabelEB_);
    e.put(std::move(recHitsEE), recHitsLabelEE_);
}
*/

DEFINE_FWK_MODULE(EcalUncalibRecHitProducerGPU);
