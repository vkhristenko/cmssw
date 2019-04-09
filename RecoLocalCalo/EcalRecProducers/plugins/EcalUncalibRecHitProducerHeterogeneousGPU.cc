// framework
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "HeterogeneousCore/Producer/interface/HeterogeneousEDProducer.h"
#include "HeterogeneousCore/Producer/interface/HeterogeneousEvent.h"
#include "HeterogeneousCore/CUDACore/interface/GPUCuda.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

// algorithm specific
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"

class EcalUncalibRecHitProducerHeterogeneousGPU
    : public HeterogeneousEDProducer<
        heterogeneous::HeterogeneousDevices<heterogeneous::GPUCuda,
                                            heterogeneous::CPU>
      > 
{
public:
    explicit EcalUncalibRecHitProducerHeterogeneousGPU(edm::ParameterSet const& ps);
    ~EcalUncalibRecHitProducerHeterogeneousGPU() override;

private:
    void acquireGPUCuda(edm::HeterogeneousEvent const&,
                        edm::EventSetup const&,
                        cuda::stream_t<>& cudaStream) override;
    void produceCPU(edm::HeterogeneousEvent&, edm::EventSetup const&) override;
    void produceGPUCuda(edm::HeterogeneousEvent&, 
                        edm::EventSetup const&,
                        cuda::stream_t<>& cudaStream) override;

private:
    edm::EDGetTokenT<EBDigiCollection> digisTokenEB_;
    edm::EDGetTokenT<EEDigiCollection> digisTokenEE_;

    std::string recHitsLabelEB_, recHitsLabelEE_;
};

EcalUncalibRecHitProducerHeterogeneousGPU::EcalUncalibRecHitProducerHeterogeneousGPU(
        const edm::ParameterSet& ps) 
    : HeterogeneousEDProducer{ps}
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

EcalUncalibRecHitProducerHeterogeneousGPU::~EcalUncalibRecHitProducerHeterogeneousGPU() {
    //
}

void EcalUncalibRecHitProducerHeterogeneousGPU::acquireGPUCuda(
        edm::HeterogeneousEvent const& event,
        edm::EventSetup const& setup,
        cuda::stream_t<>& cudaStream) 
{

}

void EcalUncalibRecHitProducerHeterogeneousGPU::produceCPU(
        edm::HeterogeneousEvent& event, 
        edm::EventSetup const& setup) 
{
    throw cms::Exception{"NotImplmented"} << "CPU version is not supported";
}

void EcalUncalibRecHitProducerHeterogeneousGPU::produceGPUCuda(
        edm::HeterogeneousEvent& event,
        edm::EventSetup const& setup,
        cuda::stream_t<>& cudaStream) 
{
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

DEFINE_FWK_MODULE(EcalUncalibRecHitProducerHeterogeneousGPU);
