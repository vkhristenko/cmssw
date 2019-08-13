// framework
#include "FWCore/Framework/interface/stream/EDProducer.h"


#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

// 
// 
// 

// format
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/RecoTypes.h"


// the kernels
#include "RecoLocalCalo/EcalRecAlgos/src/EcalRecHitBuilderKernels.h"


class EcalRecHitProducerGPU: public edm::stream::EDProducer<edm::ExternalWork> {

public:
  explicit EcalRecHitProducerGPU(edm::ParameterSet const& ps);
  ~EcalRecHitProducerGPU() override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  
private:
  
  using RecHitType = ecal::RecHit<ecal::Tag::soa>;
  void acquire(edm::Event const&, 
               edm::EventSetup const&,
               edm::WaitingTaskWithArenaHolder) override;
  void produce(edm::Event&, edm::EventSetup const&) override;
               
  void transferToHost(RecHitType& ebRecHits, RecHitType& eeRecHits, cuda::stream_t<>& cudaStream);
  
  
private:
  
//   std::string recHitsLabelEB_, recHitsLabelEE_;
//   
//   edm::EDGetTokenT<ecal::SoAUncalibratedRecHitCollection> uncalibrechitToken_; // EB and EE input together
  
  // gpu input  
  edm::EDGetTokenT<CUDAProduct< ecal::UncalibratedRecHit<ecal::Tag::ptr> > > uncalibRecHitsInEBToken_;
  edm::EDGetTokenT<CUDAProduct< ecal::UncalibratedRecHit<ecal::Tag::ptr> > > uncalibRecHitsInEEToken_;
   
  
    
  // event data
  ecal::rechit::EventOutputDataGPU eventOutputDataGPU_;
  bool shouldTransferToHost_{true};
  
  CUDAContextState cudaState_;
  
  // gpu output
  edm::EDPutTokenT<CUDAProduct<ecal::RecHit<ecal::Tag::ptr>>>  recHitsTokenEB_, recHitsTokenEE_;
  
  
  // configuration parameters
  
  uint32_t maxNumberHits_;
  uint32_t neb_, nee_; // extremely important, in particular neb_
  
};



void EcalRecHitProducerGPU::fillDescriptions(
  edm::ConfigurationDescriptions& confDesc) 
{
  
  edm::ParameterSetDescription desc;
  
  //     desc.add<edm::InputTag>("digisLabelEB", edm::InputTag("ecalDigis", "ebDigis"));
  
}


EcalRecHitProducerGPU::EcalRecHitProducerGPU(const edm::ParameterSet& ps)   {
  
  //---- input
  uncalibRecHitsInEBToken_ = consumes<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>>(ps.getParameter<edm::InputTag>("uncalibrecHitsInLabelEB"));
  uncalibRecHitsInEEToken_ = consumes<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>>(ps.getParameter<edm::InputTag>("uncalibrecHitsInLabelEE"));
        
  //---- output
  recHitsTokenEB_ = produces<CUDAProduct<ecal::RecHit<ecal::Tag::ptr>>>( ps.getParameter<std::string>("recHitsLabelEB") );
  recHitsTokenEE_ = produces<CUDAProduct<ecal::RecHit<ecal::Tag::ptr>>>( ps.getParameter<std::string>("recHitsLabelEE") );
  
  
  
  
  // max number of digis to allocate for
  maxNumberHits_ = ps.getParameter<uint32_t>("maxNumberHits");
  
  // allocate event output data
//   eventOutputDataGPU_.allocate(configParameters_, maxNumberHits_);
  eventOutputDataGPU_.allocate(maxNumberHits_);
  
  
  
}


EcalRecHitProducerGPU::~EcalRecHitProducerGPU() {
  
  // free event ouput data 
//   eventOutputDataGPU_.deallocate(configParameters_);
  eventOutputDataGPU_.deallocate();
  
}


void EcalRecHitProducerGPU::acquire(
  edm::Event const& event,
  edm::EventSetup const& setup,
  edm::WaitingTaskWithArenaHolder holder) 
{
  // cuda products
  auto const& ebUncalibRecHitsProduct = event.get(uncalibRecHitsInEBToken_);
  auto const& eeUncalibRecHitsProduct = event.get(uncalibRecHitsInEEToken_);
  // raii
  CUDAScopedContextAcquire ctx{ebUncalibRecHitsProduct, std::move(holder), cudaState_};
  // get actual object
  auto const& ebUncalibRecHits = ctx.get(ebUncalibRecHitsProduct);
  auto const& eeUncalibRecHits = ctx.get(eeUncalibRecHitsProduct);
  
  ecal::rechit::EventInputDataGPU inputDataGPU{ebUncalibRecHits, eeUncalibRecHits};
  
  
  
  neb_ = ebUncalibRecHits.size;
  nee_ = eeUncalibRecHits.size;
  
//   std::cout << " [EcalRecHitProducerGPU::acquire]  neb_:nee_ = " << neb_ << " : " << nee_ << std::endl;
  
  
  
  int nchannelsEB = ebUncalibRecHits.size;
  int offsetForInput = nchannelsEB;  // first EB and then EE
  
//   int totalChannels = 10000; // FIXME
  
  // 
  // kernel
  //
  unsigned int nchannels_per_block = 32;
  unsigned int threads_1d = 10 * nchannels_per_block;
//   unsigned int blocks_1d = threads_1d > 10*totalChannels  ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
  unsigned int blocks_1d = 2;
  
  //
  // schedule algorithms
  //
  ecal::rechit::create_ecal_rehit(
    inputDataGPU,
    eventOutputDataGPU_,
//     eventDataForScratchGPU_,
//     conditions,
//     configParameters_,
    offsetForInput,
    ctx.stream()
  );
  
  
  
  
  cudaCheck(cudaGetLastError());
  
   
}

void EcalRecHitProducerGPU::produce(
  edm::Event& event, 
  edm::EventSetup const& setup) 
{
  //DurationMeasurer<std::chrono::milliseconds> timer{std::string{"produce duration"}};
  CUDAScopedContextProduce ctx{cudaState_};
  
  // copy construct output collections
  // note, output collections do not own device memory!
  ecal::RecHit<ecal::Tag::ptr> ebRecHits{eventOutputDataGPU_};
  ecal::RecHit<ecal::Tag::ptr> eeRecHits{eventOutputDataGPU_};
  
  
  
  // set the size of eb and ee
  ebRecHits.size = neb_;
  eeRecHits.size = nee_;
  
  // shift ptrs for ee
  eeRecHits.energy += neb_;
  eeRecHits.chi2 += neb_;
  eeRecHits.did += neb_;
   
  // put into the event
  ctx.emplace(event, recHitsTokenEB_, std::move(ebRecHits));
  ctx.emplace(event, recHitsTokenEE_, std::move(eeRecHits));
  
}




void EcalRecHitProducerGPU::transferToHost(
  RecHitType& ebRecHits, RecHitType& eeRecHits,
  cuda::stream_t<>& cudaStream) {
  
  // copy from eventOutputDataGPU_ to ebRecHits/eeRecHits
//   cudaCheck( cudaMemcpyAsync(ebRecHits.energy.data(),
//                              eventOutputDataGPU_.energy,
//                              ebRecHits.energy.size() * sizeof(ecal::reco::StorageScalarType),
//                              cudaMemcpyDeviceToHost,
//                              cudaStream.id()) );
//   cudaCheck( cudaMemcpyAsync(eeRecHits.energy.data(),
//                              eventOutputDataGPU_.energy + ebRecHits.energy.size(),
//                              eeRecHits.energy.size() * sizeof(ecal::reco::StorageScalarType),
//                              cudaMemcpyDeviceToHost,
//                              cudaStream.id()) );
  
    
}
  
  

    
DEFINE_FWK_MODULE(EcalRecHitProducerGPU);

