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

#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CondFormats/DataRecord/interface/HcalRecoParamsRcd.h"
#include "CondFormats/DataRecord/interface/HcalGainWidthsRcd.h"
#include "CondFormats/DataRecord/interface/HcalGainsRcd.h"
#include "CondFormats/DataRecord/interface/HcalLUTCorrsRcd.h"
#include "CondFormats/DataRecord/interface/HcalQIEDataRcd.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"
#include "CondFormats/DataRecord/interface/HcalTimeCorrsRcd.h"
#include "CondFormats/DataRecord/interface/HcalQIETypesRcd.h"
#include "CondFormats/DataRecord/interface/HcalSiPMParametersRcd.h"
#include "CondFormats/DataRecord/interface/HcalSiPMCharacteristicsRcd.h"

//#include "CondFormats/DataRecord/interface/HcalPedestalsRcd.h"
#include "RecoLocalCalo/HcalRecProducers/src/HcalCombinedRecordsGPU.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRecoParamsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalLUTCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalPedestalWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIECodersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRespCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalTimeCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIETypesGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSiPMParametersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSiPMCharacteristicsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedPedestalsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedEffectivePedestalsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedPedestalWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedEffectivePedestalWidthsGPU.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/DeclsForKernels.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/MahiGPU.h"

class HBHERecHitProducerGPU : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit HBHERecHitProducerGPU(edm::ParameterSet const&);
    ~HBHERecHitProducerGPU() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&,
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

    edm::EDGetTokenT<CUDAProduct<hcal::DigiCollection<hcal::Flavor01>>> 
        digisTokenF01HE_;
    edm::EDGetTokenT<CUDAProduct<hcal::DigiCollection<hcal::Flavor5>>> 
        digisTokenF5HB_;
    edm::EDPutTokenT<CUDAProduct<hcal::RecHitCollection<hcal::Tag::ptr>>> rechitsM0Token_;

    hcal::mahi::ConfigParameters configParameters_;
    hcal::mahi::OutputDataGPU outputGPU_;
    CUDAContextState cudaState_;
};

HBHERecHitProducerGPU::HBHERecHitProducerGPU(edm::ParameterSet const& ps) 
    : digisTokenF01HE_{
        consumes<CUDAProduct<hcal::DigiCollection<hcal::Flavor01>>>(
            ps.getParameter<edm::InputTag>("digisLabelF01HE"))}
    , digisTokenF5HB_{
        consumes<CUDAProduct<hcal::DigiCollection<hcal::Flavor5>>>(
            ps.getParameter<edm::InputTag>("digisLabelF5HB"))}
    , rechitsM0Token_{
        produces<CUDAProduct<hcal::RecHitCollection<hcal::Tag::ptr>>>(
            ps.getParameter<std::string>("recHitsLabelM0HBHE"))}
{
    configParameters_.maxChannels = ps.getParameter<uint32_t>("maxChannels");
    configParameters_.kprep1dChannelsPerBlock = ps.getParameter<uint32_t>(
        "kprep1dChannelsPerBlock");
    configParameters_.sipmQTSShift = ps.getParameter<int>("sipmQTSShift");
    configParameters_.sipmQNTStoSum = ps.getParameter<int>("sipmQNTStoSum");
    configParameters_.firstSampleShift = ps.getParameter<int>("firstSampleShift");
    configParameters_.useEffectivePedestals = ps.getParameter<bool>("useEffectivePedestals");

    outputGPU_.allocate(configParameters_);
}

HBHERecHitProducerGPU::~HBHERecHitProducerGPU() {
    outputGPU_.deallocate(configParameters_);
}

void HBHERecHitProducerGPU::fillDescriptions(edm::ConfigurationDescriptions& cdesc) {
    edm::ParameterSetDescription desc;
    desc.add<uint32_t>("maxChannels", 10000u);
    desc.add<uint32_t>("kprep1dChannelsPerBlock", 32);
    desc.add<edm::InputTag>("digisLabelF01HE", 
        edm::InputTag{"hcalRawToDigiGPU", "f01HEDigisGPU"});
    desc.add<edm::InputTag>("digisLabelF5HB", 
        edm::InputTag{"hcalRawToDigiGPU", "f5HBDigisGPU"});
    desc.add<std::string>("recHitsLabelM0HBHE", "recHitsM0HBHE");
    desc.add<int>("sipmQTSShift", 0);
    desc.add<int>("sipmQNTStoSum", 3);
    desc.add<int>("firstSampleShift", 0);
    desc.add<bool>("useEffectivePedestals", true);

    std::string label = "hbheRecHitProducerGPU";
    cdesc.add(label, desc);
}

void HBHERecHitProducerGPU::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder holder) 
{
#ifdef HCAL_MAHI_CPUDEBUG
    auto start = std::chrono::high_resolution_clock::now();
#endif

    // input + raii
    auto const& f01HEProduct = event.get(digisTokenF01HE_);
    auto const& f5HBProduct = event.get(digisTokenF5HB_);
    CUDAScopedContextAcquire ctx{f01HEProduct, std::move(holder), cudaState_};
    auto const& f01HEDigis = ctx.get(f01HEProduct);
    auto const& f5HBDigis = ctx.get(f5HBProduct);

    hcal::mahi::InputDataGPU inputGPU{f01HEDigis, f5HBDigis};

    // conditions
    edm::ESHandle<HcalRecoParamsGPU> recoParamsHandle;
    setup.get<HcalRecoParamsRcd>().get(recoParamsHandle);
    auto const& recoParamsProduct = recoParamsHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalGainWidthsGPU> gainWidthsHandle;
    setup.get<HcalGainWidthsRcd>().get(gainWidthsHandle);
    auto const& gainWidthsProduct = gainWidthsHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalGainsGPU> gainsHandle;
    setup.get<HcalGainsRcd>().get(gainsHandle);
    auto const& gainsProduct = gainsHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalLUTCorrsGPU> lutCorrsHandle;
    setup.get<HcalLUTCorrsRcd>().get(lutCorrsHandle);
    auto const& lutCorrsProduct = lutCorrsHandle->getProduct(ctx.stream());
   
    // use only 1 depending on useEffectivePedestals
    edm::ESHandle<HcalConvertedPedestalWidthsGPU> pedestalWidthsHandle;
    edm::ESHandle<HcalConvertedEffectivePedestalWidthsGPU> 
        effectivePedestalWidthsHandle;
    if (configParameters_.useEffectivePedestals)
        setup.get<HcalConvertedEffectivePedestalWidthsRcd>()
            .get(effectivePedestalWidthsHandle);
    else
        setup.get<HcalConvertedPedestalWidthsRcd>()
            .get(pedestalWidthsHandle);
    auto const& pedestalWidthsProduct = 
        configParameters_.useEffectivePedestals
            ? effectivePedestalWidthsHandle->getProduct(ctx.stream())
            : pedestalWidthsHandle->getProduct(ctx.stream());

    edm::ESHandle<HcalConvertedPedestalsGPU> pedestalsHandle;
    setup.get<HcalConvertedPedestalsRcd>().get(pedestalsHandle);
    auto const& pedestalsProduct = pedestalsHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalConvertedEffectivePedestalsGPU> effectivePedestalsHandle;
    if (configParameters_.useEffectivePedestals)
        setup.get<HcalConvertedEffectivePedestalsRcd>().get(effectivePedestalsHandle);
    auto const* effectivePedestalsProduct = configParameters_.useEffectivePedestals
        ? &effectivePedestalsHandle->getProduct(ctx.stream())
        : nullptr;
    
    edm::ESHandle<HcalQIECodersGPU> qieCodersHandle;
    setup.get<HcalQIEDataRcd>().get(qieCodersHandle);
    auto const& qieCodersProduct = qieCodersHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalRespCorrsGPU> respCorrsHandle;
    setup.get<HcalRespCorrsRcd>().get(respCorrsHandle);
    auto const& respCorrsProduct = respCorrsHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalTimeCorrsGPU> timeCorrsHandle;
    setup.get<HcalTimeCorrsRcd>().get(timeCorrsHandle);
    auto const& timeCorrsProduct = timeCorrsHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalQIETypesGPU> qieTypesHandle;
    setup.get<HcalQIETypesRcd>().get(qieTypesHandle);
    auto const& qieTypesProduct = qieTypesHandle->getProduct(ctx.stream());

    edm::ESHandle<HcalTopology> topologyHandle;
    setup.get<HcalRecNumberingRecord>().get(topologyHandle);
    edm::ESHandle<HcalDDDRecConstants> recConstantsHandle;
    setup.get<HcalRecNumberingRecord>().get(recConstantsHandle);

    edm::ESHandle<HcalSiPMParametersGPU> sipmParametersHandle;
    setup.get<HcalSiPMParametersRcd>().get(sipmParametersHandle);
    auto const& sipmParametersProduct = sipmParametersHandle->getProduct(ctx.stream());
    
    edm::ESHandle<HcalSiPMCharacteristicsGPU> sipmCharacteristicsHandle;
    setup.get<HcalSiPMCharacteristicsRcd>().get(sipmCharacteristicsHandle);
    auto const& sipmCharacteristicsProduct = sipmCharacteristicsHandle->getProduct(ctx.stream());

    // bundle up conditions
    hcal::mahi::ConditionsProducts conditions{
        gainWidthsProduct, gainsProduct, lutCorrsProduct,
        pedestalWidthsProduct, pedestalsProduct,
        qieCodersProduct, recoParamsProduct,
        respCorrsProduct, timeCorrsProduct,
        qieTypesProduct, 
        sipmParametersProduct, sipmCharacteristicsProduct,
        effectivePedestalsProduct,
        topologyHandle.product(), recConstantsHandle.product(),
        pedestalsHandle->offsetForHashes()
    };

    hcal::mahi::entryPoint(inputGPU, outputGPU_, conditions,
        configParameters_, ctx.stream());

#ifdef HCAL_MAHI_CPUDEBUG
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "acquire  duration = " << duration << std::endl;
#endif
}

void HBHERecHitProducerGPU::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    CUDAScopedContextProduce ctx{cudaState_};

    // copy construct a new guy/view and place cuda product into the event
    hcal::RecHitCollection<hcal::Tag::ptr> recHits = outputGPU_.recHits;
    ctx.emplace(event, rechitsM0Token_, std::move(recHits));
}

DEFINE_FWK_MODULE(HBHERecHitProducerGPU);
