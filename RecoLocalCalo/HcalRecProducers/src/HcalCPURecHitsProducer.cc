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
    hcal::RecHitCollection<hcal::Tag::soa> tmpRecHits;
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
    auto const& recHitsM0Product = event.get(recHitsM0TokenIn_);
    CUDAScopedContextAcquire ctx{recHitsM0Product, std::move(taskHolder)};
    auto const& recHitsM0 = ctx.get(recHitsM0Product);

    // resize tmp buffers
    tmpRecHits.resize(recHitsM0.size);

    auto lambdaToTransfer = [&ctx](auto& dest, auto* src) {
        using vector_type = typename std::remove_reference<decltype(dest)>::type;
        using type = typename vector_type::value_type;
        cudaCheck( cudaMemcpyAsync(dest.data(),
                                   src,
                                   dest.size() * sizeof(type),
                                   cudaMemcpyDeviceToHost,
                                   ctx.stream().id()) );
    };

    lambdaToTransfer(tmpRecHits.energy, recHitsM0.energy);
    lambdaToTransfer(tmpRecHits.time, recHitsM0.time);
    lambdaToTransfer(tmpRecHits.did, recHitsM0.did);

    // enqeue transfers
    /*
    cudaCheck( cudaMemcpyAsync(dataf01he.data(),
                               f01HERecHits.data,
                               dataf01he.size() * sizeof(uint16_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(dataf5hb.data(),
                               f5HBRecHits.data,
                               dataf5hb.size() * sizeof(uint16_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(idsf01he.data(),
                               f01HERecHits.ids,
                               idsf01he.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(idsf5hb.data(),
                               f5HBRecHits.ids,
                               idsf5hb.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
                               */
}

void HcalCPURecHitsProducer::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    event.put(recHitsM0TokenOut_, 
        std::make_unique<hcal::RecHitCollection<hcal::Tag::soa>>(
            std::move(tmpRecHits)));

    // output collections
    /*
    auto f01he = std::make_unique<edm::DataFrameContainer>(
        stridef01he, HcalEndcap, idsf01he.size());
    auto f5hb = std::make_unique<edm::DataFrameContainer>(
        stridef5hb, HcalBarrel, idsf5hb.size());
    
    // cast constness away
    // use pointers to buffers instead of move operator= semantics (or swap)
    // cause we have different allocators in there...
    auto *dataf01hetmp = const_cast<uint16_t*>(f01he->data().data());
    auto *dataf5hbtmp = const_cast<uint16_t*>(f5hb->data().data());

    auto *idsf01hetmp = const_cast<uint32_t*>(f01he->ids().data());
    auto idsf5hbtmp = const_cast<uint32_t*>(f5hb->ids().data());

    // copy data
    std::memcpy(dataf01hetmp, dataf01he.data(), dataf01he.size() * sizeof(uint16_t));
    std::memcpy(dataf5hbtmp, dataf5hb.data(), dataf5hb.size() * sizeof(uint16_t));
    std::memcpy(idsf01hetmp, idsf01he.data(), idsf01he.size() * sizeof(uint32_t));
    std::memcpy(idsf5hbtmp, idsf5hb.data(), idsf5hb.size() * sizeof(uint32_t));

    event.put(digisF01HETokenOut_, std::move(f01he));
    event.put(digisF5HBTokenOut_, std::move(f5hb));
    */
}

DEFINE_FWK_MODULE(HcalCPURecHitsProducer);
