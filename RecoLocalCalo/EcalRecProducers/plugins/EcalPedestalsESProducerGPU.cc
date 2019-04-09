#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "EcalPedestalsGPU.h"

#include <iostream>

class EcalPedestalsESProducerGPU : public edm::ESProducer {
public:
    explicit EcalPedestalsESProducerGPU(edm::ParameterSet const&);
    
    std::unique_ptr<EcalPedestalsGPU> produce(EcalPedestalsRcd const&);

    static void fillDescriptions(edm::ConfigurationDescriptions& desc);

private:
    std::string pedestalsLabel_;
};

EcalPedestalsESProducerGPU::EcalPedestalsESProducerGPU(edm::ParameterSet const& ps) 
    : pedestalsLabel_{ps.getParameter<std::string>("pedestalsLabel")}
{
    std::string name = ps.getParameter<std::string>("ComponentName");
    setWhatProduced(this, name);
}

std::unique_ptr<EcalPedestalsGPU> EcalPedestalsESProducerGPU::produce(
        EcalPedestalsRcd const& record)
{
    // retrieve conditions in old format 
    edm::ESTransientHandle<EcalPedestals> pedestals;
    record.get(pedestalsLabel_, pedestals);

    std::cout << __PRETTY_FUNCTION__ << std::endl;

    return std::make_unique<EcalPedestalsGPU>(*pedestals);
}

void EcalPedestalsESProducerGPU::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    std::string label = "ecalPedestalsESProducerGPU";
    desc.add<std::string>("ComponentName", "");
    desc.add<std::string>("pedestalsLabel", "")->setComment("Pedestals Label");
    confDesc.add(label, desc);
}

DEFINE_FWK_EVENTSETUP_MODULE(EcalPedestalsESProducerGPU);
