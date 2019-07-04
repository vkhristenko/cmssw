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

// algorithm specific

#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>

class EcalRawToDigiGPU
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit EcalRawToDigiGPU(edm::ParameterSet const& ps);
    ~EcalRawToDigiGPU() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<FEDRawDataCollection> rawDataToken_;

    CUDAContextState cudaState_;

    std::vector<int> fedsToUnpack_;
};

void EcalRawToDigiGPU::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("InputLabel", edm::InputTag("rawDataCollector"));
    std::vector<int> feds(54);
    for (uint32_t i=0; i<54; ++i)
        feds[i] = i+601;
    desc.add<std::vector<int>>("FEDs", feds);

    std::string label = "ecalRawToDigiGPU";
    confDesc.add(label, desc);
}

EcalRawToDigiGPU::EcalRawToDigiGPU(
        const edm::ParameterSet& ps) 
    : rawDataToken_{consumes<FEDRawDataCollection>(ps.getParameter<edm::InputTag>(
        "InputLabel"))}
    , fedsToUnpack_{ps.getParameter<std::vector<int>>("FEDs")}
{}

EcalRawToDigiGPU::~EcalRawToDigiGPU() {
}

void EcalRawToDigiGPU::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder holder) 
{
    // FIXME: remove debugging
    auto start = std::chrono::high_resolution_clock::now();

    //DurationMeasurer<std::chrono::milliseconds> timer{std::string{"acquire duration"}};

    // raii
    CUDAScopedContextAcquire ctx{event.streamID(), std::move(holder), cudaState_};

    // event data
    edm::Handle<FEDRawDataCollection> rawDataHandle;
    event.getByToken(rawDataToken_, rawDataHandle);

    // iterate over feds
    for (auto const& fed : fedsToUnpack_) {
        //std::cout << "fed: " << fed << std::endl;
        auto const& data = rawDataHandle->FEDData(fed);
        auto const nbytes = data.size();

        printf("fed = %d nbytes = %lu\n", fed, nbytes);
    }

    // FIXME: remove debugging
    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)    .count(); 
    std::cout << "acquire  duration = " << duration << std::endl;
}

void EcalRawToDigiGPU::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    //DurationMeasurer<std::chrono::milliseconds> timer{std::string{"produce duration"}};
    CUDAScopedContextProduce ctx{cudaState_};
}

DEFINE_FWK_MODULE(EcalRawToDigiGPU);
