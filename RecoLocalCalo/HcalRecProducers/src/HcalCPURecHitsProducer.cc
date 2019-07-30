#include <iostream>

// framework
#include "FWCore/Framework/interface/stream/EDProducer.h"
//#include "HeterogeneousCore/Producer/interface/HeterogeneousEDProducer.h"
//#include "HeterogeneousCore/Producer/interface/HeterogeneousEvent.h"

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "CUDADataFormats/HcalRecHitSoA/interface/RecHitCollection.h"

class HcalCPURecHitsProducer
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit HcalCPURecHitsProducer(edm::ParameterSet const& ps);
    ~HcalCPURecHitsProducer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<CUDAProduct<hcal::RecHitCollection<hcal::Tag::ptr>>> 
        recHitsM0TokenIn_;
    edm::EDPutTokenT<hcal::RecHitCollection<hcal::Tag::soa>> 
        recHitsM0TokenOut_;

    // to pass from acquire to produce
    hcal::RecHitCollection<hcal::Tag::soa> tmpRecHits_;
};

void HcalCPURecHitsProducer::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("recHitsM0LabelIn", 
        edm::InputTag{"hbheRecHitProducerGPU", "recHitsM0HBHE"});
    desc.add<std::string>("recHitsM0LabelOut", "recHitsM0HBHE");

    std::string label = "hcalCPURecHitsProducer";
    confDesc.add(label, desc);
}

HcalCPURecHitsProducer::HcalCPURecHitsProducer(
        const edm::ParameterSet& ps) 
    : recHitsM0TokenIn_{
        consumes<CUDAProduct<hcal::RecHitCollection<hcal::Tag::ptr>>>(
            ps.getParameter<edm::InputTag>("recHitsM0LabelIn"))}
    , recHitsM0TokenOut_{
        produces<hcal::RecHitCollection<hcal::Tag::soa>>("recHitsM0LabelOut")}
{}

HcalCPURecHitsProducer::~HcalCPURecHitsProducer() {}

void HcalCPURecHitsProducer::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder taskHolder) 
{
    // retrieve data/ctx
    auto const& recHitsProduct = event.get(recHitsM0TokenIn_);
    CUDAScopedContextAcquire ctx{recHitsProduct, std::move(taskHolder)};
    auto const& recHits = ctx.get(recHitsProduct);

    // resize tmp buffers
    tmpRecHits_.resize(recHits.size);

    auto lambdaToTransfer = [&ctx](auto& dest, auto* src) {
        using vector_type = typename std::remove_reference<decltype(dest)>::type;
        using type = typename vector_type::value_type;
        cudaCheck( cudaMemcpyAsync(dest.data(),
                                   src,
                                   dest.size() * sizeof(type),
                                   cudaMemcpyDeviceToHost,
                                   ctx.stream().id()) );
    };

    lambdaToTransfer(tmpRecHits_.energy, recHits.energy);
    lambdaToTransfer(tmpRecHits_.chi2, recHits.chi2);
    lambdaToTransfer(tmpRecHits_.energyM0, recHits.energyM0);
    lambdaToTransfer(tmpRecHits_.timeM0, recHits.timeM0);
    lambdaToTransfer(tmpRecHits_.did, recHits.did);
}

void HcalCPURecHitsProducer::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    event.put(recHitsM0TokenOut_, 
        std::make_unique<hcal::RecHitCollection<hcal::Tag::soa>>(
            std::move(tmpRecHits_)));
}

DEFINE_FWK_MODULE(HcalCPURecHitsProducer);
