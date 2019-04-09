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

#include <iostream>

#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "EcalPedestalsGPU.h"

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
    // raii
    CUDAScopedContext ctx{event.streamID(), std::move(holder)};

    // retrieve device ptrs to conditions
    edm::ESHandle<EcalPedestalsGPU> pedestals;
    setup.get<EcalPedestalsRcd>().get(pedestals);

    std::cout << "acquire\n";

    ctxToken_ = ctx.toToken();
}

void EcalUncalibRecHitProducerGPU::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
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
