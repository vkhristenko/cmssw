// framework
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

// algorithm specific
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"

class EcalUncalibRecHitProducerGPU : public edm::stream::EDProducer<> {
public:
    explicit EcalUncalibRecHitProducerGPU(edm::ParameterSet const& ps);
    ~EcalUncalibRecHitProducerGPU() override;

    // main entry point
    void produce(edm::Event&, const edm::EventSetup& es) override;
    void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
        override;

private:
    edm::EDGetTokenT<EBDigiCollection> digisTokenEB_;
    edm::EDGetTokenT<EEDigiCollection> digisTokenEE_;

    std::string recHitsLabelEB_, recHitsLabelEE_;

    struct ConditionsPerDetector {
    //    std::array<SampleMatrixGainArray, 2> noiseCorrections;
    };

    struct ConditionsPerChannel {
        
    };

    ConditionsPerDetector condsPerDetector_;
    ConditionsPerChannel condsPerChannel_;
};

EcalUncalibRecHitProducerGPU::EcalUncalibRecHitProducerGPU(
        const edm::ParameterSet& ps) {
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

void EcalUncalibRecHitProducerGPU::beginLuminosityBlock(
        edm::LuminosityBlock const&,
        edm::EventSetup const&) {
    //
    // per detector conditions
    //
    

    // 
    // per channel conditions
    //
}

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

DEFINE_FWK_MODULE(EcalUncalibRecHitProducerGPU);
