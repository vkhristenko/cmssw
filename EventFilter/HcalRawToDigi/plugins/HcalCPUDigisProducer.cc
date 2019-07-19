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

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CUDADataFormats/HcalDigi/interface/DigiCollection.h"

class HcalCPUDigisProducer
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit HcalCPUDigisProducer(edm::ParameterSet const& ps);
    ~HcalCPUDigisProducer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<CUDAProduct<hcal::DigiCollection<hcal::Flavor01>>> 
        digisF01HETokenIn_;
    edm::EDGetTokenT<CUDAProduct<hcal::DigiCollection<hcal::Flavor5>>> 
        digisF5HBTokenIn_;

    edm::EDPutTokenT<edm::DataFrameContainer> digisF01HETokenOut_,
        digisF5HBTokenOut_;

    // FIXME better way to pass pointers from acquire to produce?
    std::vector<uint32_t, CUDAHostAllocator<uint32_t>> idsf01he, idsf5hb;
    std::vector<uint16_t, CUDAHostAllocator<uint16_t>> dataf01he, dataf5hb;
    uint32_t stridef01he, stridef5hb;
};

void HcalCPUDigisProducer::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("digisLabelF01HEIn", 
        edm::InputTag{"hcalRawToDigiGPU", "f01HEDigisGPU"});
    desc.add<edm::InputTag>("digisLabelF5HBIn", 
        edm::InputTag{"hcalRawToDigiGPU", "f5HBDigisGPU"});
    desc.add<std::string>("digisLabelF01HEOut", "f01HEDigis");
    desc.add<std::string>("digisLabelF5HBOut", "f5HBDigis");

    std::string label = "hcalCPUDigisProducer";
    confDesc.add(label, desc);
}

HcalCPUDigisProducer::HcalCPUDigisProducer(
        const edm::ParameterSet& ps) 
    : digisF01HETokenIn_{consumes<CUDAProduct<hcal::DigiCollection<hcal::Flavor01>>>(
        ps.getParameter<edm::InputTag>("digisLabelF01HEIn"))}
    , digisF5HBTokenIn_{consumes<CUDAProduct<hcal::DigiCollection<hcal::Flavor5>>>(
        ps.getParameter<edm::InputTag>("digisLabelF5HBIn"))}
    , digisF01HETokenOut_{produces<edm::DataFrameContainer>(
        ps.getParameter<std::string>("digisLabelF01HEOut"))}
    , digisF5HBTokenOut_{produces<edm::DataFrameContainer>(
        ps.getParameter<std::string>("digisLabelF5HBOut"))}
{}

HcalCPUDigisProducer::~HcalCPUDigisProducer() {}

void HcalCPUDigisProducer::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder taskHolder) 
{
    // retrieve data/ctx
    auto const& f01HEProduct = event.get(digisF01HETokenIn_);
    auto const& f5HBProduct = event.get(digisF5HBTokenIn_);
    CUDAScopedContextAcquire ctx{f01HEProduct, std::move(taskHolder)};
    auto const& f01HEDigis = ctx.get(f01HEProduct);
    auto const& f5HBDigis = ctx.get(f5HBProduct);

    // resize out tmp buffers
    idsf01he.resize(f01HEDigis.ndigis);
    dataf01he.resize(f01HEDigis.ndigis * f01HEDigis.stride);
    idsf5hb.resize(f5HBDigis.ndigis);
    dataf5hb.resize(f5HBDigis.ndigis * f5HBDigis.stride);
    stridef01he = f01HEDigis.stride;
    stridef5hb = f5HBDigis.stride;

    // enqeue transfers
    cudaCheck( cudaMemcpyAsync(dataf01he.data(),
                               f01HEDigis.data,
                               dataf01he.size() * sizeof(uint16_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(dataf5hb.data(),
                               f5HBDigis.data,
                               dataf5hb.size() * sizeof(uint16_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(idsf01he.data(),
                               f01HEDigis.ids,
                               idsf01he.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(idsf5hb.data(),
                               f5HBDigis.ids,
                               idsf5hb.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
}

void HcalCPUDigisProducer::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    // output collections
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
}

DEFINE_FWK_MODULE(HcalCPUDigisProducer);
