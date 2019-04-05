#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

class EcalUncalibRecHitProducerGPU : public edm::stream::EDProducer<> {
public:
    explicit EcalUncalibRecHitProducerGPU(const edm::ParameterSet& ps);
    ~EcalUncalibRecHitProducerGPU() override;

    // main entry point
    void produce(edm::Event&, const edm::EventSetup& es) override;

private:
};

EcalUncalibRecHitProducerGPU::EcalUncalibRecHitProducerGPU(
        const edm::ParameterSet& ps) {
    // 
}

EcalUncalibRecHitProducerGPU::~EcalUncalibRecHitProducerGPU() {
    //
}

void EcalUncalibRecHitProducerGPU::produce(edm::Event&, const edm::EventSetup&) {
    //
}

DEFINE_FWK_MODULE(EcalUncalibRecHitProducerGPU);
