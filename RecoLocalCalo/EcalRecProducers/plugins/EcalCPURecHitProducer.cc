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

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

class EcalCPURecHitProducer
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit EcalCPURecHitProducer(edm::ParameterSet const& ps);
    ~EcalCPURecHitProducer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<CUDAProduct<ecal::RecHit<ecal::Tag::ptr>>> recHitsInEBToken_, recHitsInEEToken_;
    edm::EDPutTokenT<ecal::RecHit<ecal::Tag::soa>> recHitsOutEBToken_, recHitsOutEEToken_;

    ecal::RecHit<ecal::Tag::soa> recHitsEB_, recHitsEE_;
    bool containsTimingInformation_;
};

void EcalCPURecHitProducer::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("recHitsInLabelEB", edm::InputTag{"ecalRecHitProducerGPU", "EcalRecHitsEB"});
    desc.add<edm::InputTag>("recHitsInLabelEE", edm::InputTag{"ecalRecHitProducerGPU", "EcalRecHitsEE"});
    desc.add<std::string>("recHitsOutLabelEB", "EcalRecHitsEB");
    desc.add<std::string>("recHitsOutLabelEE", "EcalRecHitsEE");
    desc.add<bool>("containsTimingInformation", false);

    std::string label = "ecalCPURecHitProducer";
    confDesc.add(label, desc);
}

EcalCPURecHitProducer::EcalCPURecHitProducer(
        const edm::ParameterSet& ps) 
    : recHitsInEBToken_{consumes<CUDAProduct<ecal::RecHit<ecal::Tag::ptr>>>(ps.getParameter<edm::InputTag>("recHitsInLabelEB"))}
    , recHitsInEEToken_{consumes<CUDAProduct<ecal::RecHit<ecal::Tag::ptr>>>(ps.getParameter<edm::InputTag>("recHitsInLabelEE"))}
    , recHitsOutEBToken_{produces<ecal::RecHit<ecal::Tag::soa>>(ps.getParameter<std::string>("recHitsOutLabelEB"))}
    , recHitsOutEEToken_{produces<ecal::RecHit<ecal::Tag::soa>>(ps.getParameter<std::string>("recHitsOutLabelEE"))}
    , containsTimingInformation_{ps.getParameter<bool>("containsTimingInformation")}
{}

EcalCPURecHitProducer::~EcalCPURecHitProducer() {}

void EcalCPURecHitProducer::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder taskHolder) 
{
    // retrieve data/ctx
    auto const& ebRecHitsProduct = event.get(recHitsInEBToken_);
    auto const& eeRecHitsProduct = event.get(recHitsInEEToken_);
    CUDAScopedContextAcquire ctx{ebRecHitsProduct, std::move(taskHolder)};
    auto const& ebRecHits = ctx.get(ebRecHitsProduct);
    auto const& eeRecHits = ctx.get(eeRecHitsProduct);

    // resize the output buffers
    recHitsEB_.resize(ebRecHits.size);
    recHitsEE_.resize(eeRecHits.size);

//     std::cout << " [EcalCPURecHitProducer::acquire] ebRecHits.size = " << ebRecHits.size << std::endl;
//     std::cout << " [EcalCPURecHitProducer::acquire] eeRecHits.size = " << eeRecHits.size << std::endl;
    
   
    // AM: FIXME : why all "uint32_t" and not "float" where needed?
    
    
    // enqeue transfers
    cudaCheck( cudaMemcpyAsync(recHitsEB_.did.data(),
                               ebRecHits.did,
                               recHitsEB_.did.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.did.data(),
                               eeRecHits.did,
                               recHitsEE_.did.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
// 
//     ./CUDADataFormats/EcalRecHitSoA/interface/RecoTypes.h:using StorageScalarType = float;
// 
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.energy.data(),
                               ebRecHits.energy,
                               recHitsEB_.energy.size() * sizeof(::ecal::reco::StorageScalarType),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.energy.data(),
                               eeRecHits.energy,
                               recHitsEE_.energy.size() * sizeof(::ecal::reco::StorageScalarType),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.chi2.data(),
                               ebRecHits.chi2,
                               recHitsEB_.chi2.size() * sizeof(::ecal::reco::StorageScalarType),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.chi2.data(),
                               eeRecHits.chi2,
                               recHitsEE_.chi2.size() * sizeof(::ecal::reco::StorageScalarType),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );

    cudaCheck( cudaMemcpyAsync(recHitsEB_.extra.data(),
                               ebRecHits.extra,
                               recHitsEB_.extra.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.extra.data(),
                               eeRecHits.extra,
                               recHitsEE_.extra.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.flagBits.data(),
                               ebRecHits.flagBits,
                               recHitsEB_.flagBits.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.flagBits.data(),
                               eeRecHits.flagBits,
                               recHitsEE_.flagBits.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
 
    
    
    
//     for (unsigned int ieb = 0; ieb <  ebRecHits.size ; ieb++) {
//       if (recHitsEB_.extra[ieb] != 0 ) std::cout << " [ " << ieb << " :: " << ebRecHits.size << " ] [ " << recHitsEB_.did[ieb] << " ] eb extra = " << recHitsEB_.extra[ieb] << std::endl;
//     }
    
//     
//     for (unsigned int ieb = 0; ieb <  ebRecHits.size ; ieb++) {
//       if (recHitsEB_.energy[ieb] != 0 ) std::cout << " [ " << ieb << " :: " << ebRecHits.size << " ] [ " << recHitsEB_.did[ieb] << " ] eb energy = " << recHitsEB_.energy[ieb] << std::endl;
//     }
//     
//     for (unsigned int iee = 0; iee <  eeRecHits.size ; iee++) {
//       if (recHitsEE_.energy[iee] != 0 ) std::cout << " [ " << iee << " :: " << eeRecHits.size << " ] [ " << recHitsEE_.did[iee] << " ] ee energy = " << recHitsEE_.energy[iee] << std::endl;
//     }
//     
    
    
    
    
}

void EcalCPURecHitProducer::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    // tmp vectors
    auto recHitsOutEB = std::make_unique<ecal::RecHit<ecal::Tag::soa>>(std::move(recHitsEB_));
    auto recHitsOutEE = std::make_unique<ecal::RecHit<ecal::Tag::soa>>(std::move(recHitsEE_));

    // put into event
    event.put(recHitsOutEBToken_, std::move(recHitsOutEB));
    event.put(recHitsOutEEToken_, std::move(recHitsOutEE));
}

DEFINE_FWK_MODULE(EcalCPURecHitProducer);


