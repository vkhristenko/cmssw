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

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
//#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
//#include "CUDADataFormats/HcalDigi/interface/DigisCollection.h"

//#include "CondFormats/DataRecord/interface/HcalMappingElectronicsRcd.h"
//#include "EventFilter/HcalRawToDigi/interface/ElectronicsMappingGPU.h"

#include "EventFilter/HcalRawToDigi/plugins/DeclsForKernels.h"
#include "EventFilter/HcalRawToDigi/plugins/DecodeGPU.h"
#include "EventFilter/HcalRawToDigi/plugins/ElectronicsMappingGPU.h"
#include "CondFormats/DataRecord/interface/HcalElectronicsMapRcd.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

class HcalRawToDigiGPU
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit HcalRawToDigiGPU(edm::ParameterSet const& ps);
    ~HcalRawToDigiGPU() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<FEDRawDataCollection> rawDataToken_;
    //edm::EDPutTokenT<CUDAProduct<hcal::DigisCollection>> digisEBToken_, 
    //    digisEEToken_;

    CUDAContextState cudaState_;

    std::vector<int> fedsToUnpack_;

    hcal::raw::ConfigurationParameters config_;
    // FIXME move this to use raii
    hcal::raw::InputDataCPU inputCPU_;
    hcal::raw::InputDataGPU inputGPU_;
    //hcal::raw::OutputDataGPU outputGPU_;
    //hcal::raw::ScratchDataGPU scratchGPU_;
    //hcal::raw::OutputDataCPU outputCPU_;
};

void HcalRawToDigiGPU::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("InputLabel", edm::InputTag("rawDataCollector"));
    auto nFeds = FEDNumbering::MAXHCALuTCAFEDID - FEDNumbering::MINHCALuTCAFEDID + 1;
    std::vector<int> feds(nFeds);
    for (int i=0; i<nFeds; ++i)
        feds[i] = i + FEDNumbering::MINHCALuTCAFEDID;
    desc.add<std::vector<int>>("FEDs", feds);
    desc.add<uint32_t>("maxChannels", 20000);
//    desc.add<std::string>("digisLabelEB", "ebDigisGPU");
//    desc.add<std::string>("digisLabelEE", "eeDigisGPU");

    std::string label = "hcalRawToDigiGPU";
    confDesc.add(label, desc);
}

HcalRawToDigiGPU::HcalRawToDigiGPU(
        const edm::ParameterSet& ps) 
    : rawDataToken_{consumes<FEDRawDataCollection>(ps.getParameter<edm::InputTag>(
        "InputLabel"))}
    //, digisEBToken_{produces<CUDAProduct<hcal::DigisCollection>>(
    //    ps.getParameter<std::string>("digisLabelEB"))}
    //, digisEEToken_{produces<CUDAProduct<hcal::DigisCollection>>(
    //    ps.getParameter<std::string>("digisLabelEE"))}
    , fedsToUnpack_{ps.getParameter<std::vector<int>>("FEDs")}
{
    config_.maxChannels = ps.getParameter<uint32_t>("maxChannels");

    inputCPU_.allocate();
    inputGPU_.allocate();
    /*
    outputGPU_.allocate(config_);
    scratchGPU_.allocate(config_);
    outputCPU_.allocate();
    */
}

HcalRawToDigiGPU::~HcalRawToDigiGPU() {
    inputGPU_.deallocate();
    /*
    outputGPU_.deallocate(config_);
    scratchGPU_.deallocate(config_);
    */
}

void HcalRawToDigiGPU::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder holder) 
{
    // raii
    CUDAScopedContextAcquire ctx{event.streamID(), std::move(holder), cudaState_};

    // conditions
    edm::ESHandle<hcal::raw::ElectronicsMappingGPU> eMappingHandle;
    setup.get<HcalElectronicsMapRcd>().get(eMappingHandle);
    auto const& eMappingProduct = eMappingHandle->getProduct(ctx.stream());

    // bundle up conditions
    hcal::raw::ConditionsProducts conditions{eMappingProduct};

    // event data
    edm::Handle<FEDRawDataCollection> rawDataHandle;
    event.getByToken(rawDataToken_, rawDataHandle);

    // iterate over feds
    // TODO: another idea
    //   - loop over all feds to unpack and enqueue cuda memcpy 
    //   - accumulate the sizes
    //   - after the loop launch cuda memcpy for sizes
    //   - enqueue the kernel
    uint32_t currentCummOffset = 0;
    uint32_t counter = 0;
    for (auto const& fed : fedsToUnpack_) {
        auto const& data = rawDataHandle->FEDData(fed);
        auto const nbytes = data.size();

        // skip empty feds
        if (nbytes < hcal::raw::empty_event_size)
            continue;
        
        printf("fed = %d nbytes = %lu\n", fed, nbytes);

        // copy raw data into plain buffer
        std::memcpy(inputCPU_.data.data() + currentCummOffset, data.data(), nbytes);
        // set the offset in bytes from the start
        inputCPU_.offsets[counter] = currentCummOffset;
        inputCPU_.feds[counter] = fed;

        // this is the current offset into the vector
        currentCummOffset += nbytes;
        ++counter;
    }

    hcal::raw::entryPoint(
        inputCPU_, inputGPU_, /*outputGPU_, scratchGPU_, outputCPU_,*/
        /*conditions,*/ ctx.stream(), counter, currentCummOffset);
}

void HcalRawToDigiGPU::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    CUDAScopedContextProduce ctx{cudaState_};

    /*
    // get the number of channels 
    auto const nchannelsEB = outputCPU_.nchannels[0];
    auto const nchannelsEE = outputCPU_.nchannels[1];
    
    hcal::DigisCollection digisEB{outputGPU_.idsEB, 
        outputGPU_.samplesEB, nchannelsEB};
    hcal::DigisCollection digisEE{outputGPU_.idsEE,
        outputGPU_.samplesEE, nchannelsEE};

    ctx.emplace(event, digisEBToken_, std::move(digisEB));
    ctx.emplace(event, digisEEToken_, std::move(digisEE));
    */
}

DEFINE_FWK_MODULE(HcalRawToDigiGPU);
