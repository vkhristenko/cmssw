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
  
  std::string recHitsLabelEB_, recHitsLabelEE_;
//   
//   edm::EDGetTokenT<ecal::SoAUncalibratedRecHitCollection> uncalibrechitToken_; // EB and EE input together
  
  // gpu input  
  edm::EDGetTokenT<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr> > > uncalibRecHitsInEBToken_;
  edm::EDGetTokenT<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr> > > uncalibRecHitsInEEToken_;
  
  
    
  // event data
//   ecal::rechit::EventInputDataGPU eventInputDataGPU_;
  ecal::rechit::EventOutputDataGPU eventOutputDataGPU_;
//   ecal::rechit::EventDataForScratchGPU eventDataForScratchGPU_;
  bool shouldTransferToHost_{true};
  
  CUDAContextState cudaState_;
  
  // output
  std::unique_ptr< ecal::RecHit<ecal::Tag::soa> > ebRecHits_{nullptr};
  std::unique_ptr< ecal::RecHit<ecal::Tag::soa> > eeRecHits_{nullptr};
  
  // configuration parameters
//   ecal::multifit::ConfigurationParameters configParameters_;
  
  uint32_t maxNumberHits_;
  
};



void EcalRecHitProducerGPU::fillDescriptions(
  edm::ConfigurationDescriptions& confDesc) 
{
  
  edm::ParameterSetDescription desc;
  
  //     desc.add<edm::InputTag>("digisLabelEB", edm::InputTag("ecalDigis", "ebDigis"));
  
}


EcalRecHitProducerGPU::EcalRecHitProducerGPU(const edm::ParameterSet& ps)   {
  
  //---- input
  uncalibRecHitsInEBToken_ = consumes<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>>(ps.getParameter<edm::InputTag>("recHitsInLabelEB"));
  uncalibRecHitsInEEToken_ = consumes<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>>(ps.getParameter<edm::InputTag>("recHitsInLabelEE"));
        
  //---- output
  recHitsLabelEB_ = ps.getParameter<std::string>("recHitsLabelEB");
  recHitsLabelEE_ = ps.getParameter<std::string>("recHitsLabelEE");
  
  produces<ecal::SoARecHitCollection>(recHitsLabelEB_);
  produces<ecal::SoARecHitCollection>(recHitsLabelEE_);

  
  
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
  // retrieve data/ctx
  auto const& ebUncalibRecHitsProduct = event.get(uncalibRecHitsInEBToken_);
  auto const& eeUncalibRecHitsProduct = event.get(uncalibRecHitsInEEToken_);
  CUDAScopedContextAcquire ctx{ebUncalibRecHitsProduct, std::move(holder)};
  auto const& ebUncalibRecHits = ctx.get(ebUncalibRecHitsProduct);
  auto const& eeUncalibRecHits = ctx.get(eeUncalibRecHitsProduct);
  
   
//   int nchannelsEB = ebUncalibRecHits.size;
  int nchannelsEB = 10;
  
  
  int totalChannels = 10000;
  // 
  // kernel
  //
  unsigned int nchannels_per_block = 32;
  unsigned int threads_1d = 10 * nchannels_per_block;
//   unsigned int blocks_1d = threads_1d > 10*totalChannels  ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
  unsigned int blocks_1d = 2;
  
  
//   edm::Handle<ecal::SoAUncalibratedRecHitCollection> hRecHitsGPUEB, hRecHitsGPUEE;
//   event.getByToken(recHitsGPUEB_, hRecHitsGPUEB);
//   event.getByToken(recHitsGPUEE_, hRecHitsGPUEE);
//   
//   auto recHitsCPUEB = std::make_unique<EBUncalibratedRecHitCollection>();
//   auto recHitsCPUEE = std::make_unique<EEUncalibratedRecHitCollection>();
//   recHitsCPUEB->reserve(hRecHitsGPUEB->amplitude.size());
//   recHitsCPUEE->reserve(hRecHitsGPUEE->amplitude.size());
//   
//   for (uint32_t i=0; i<hRecHitsGPUEB->amplitude.size(); ++i) {
//     recHitsCPUEB->emplace_back(
//       DetId{hRecHitsGPUEB->did[i]},
//       hRecHitsGPUEB->amplitude[i],
//       hRecHitsGPUEB->pedestal[i],
//       hRecHitsGPUEB->jitter[i],
//       hRecHitsGPUEB->chi2[i],
//       hRecHitsGPUEB->flags[i]
//     );
//     
    
    
    
//   ecal::rechit::create_ecal_rehit (
//                      ebUncalibRecHits.amplitude,
//                      ebRecHits_->energy,
//                      nchannelsEB
//   );
  
//   ecal::rechit::kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
//                          ebUncalibRecHits.amplitudes,
//                          ebRecHits_.energy,
//                          nchannelsEB
//   );
  
  cudaCheck(cudaGetLastError());
  
  
  
  
  // output
  ebRecHits_ = std::move( std::make_unique<ecal::RecHit<ecal::Tag::soa>>() );
  eeRecHits_ = std::move( std::make_unique<ecal::RecHit<ecal::Tag::soa>>() );
  
  
  
  
  
}

void EcalRecHitProducerGPU::produce(
  edm::Event& event, 
  edm::EventSetup const& setup) 
{
  //DurationMeasurer<std::chrono::milliseconds> timer{std::string{"produce duration"}};
  CUDAScopedContextProduce ctx{cudaState_};
  
  if (shouldTransferToHost_) {
    // rec hits objects were not originally member variables
    transferToHost(*ebRecHits_, *eeRecHits_, ctx.stream());
    
    // TODO
//     for now just sync on the host when transferring back products
        cudaStreamSynchronize(ctx.stream().id());
  }
  
  
  
  
  event.put(std::move(ebRecHits_), recHitsLabelEB_);
  event.put(std::move(eeRecHits_), recHitsLabelEE_);
  
}




void EcalRecHitProducerGPU::transferToHost(
  RecHitType& ebRecHits, RecHitType& eeRecHits,
  cuda::stream_t<>& cudaStream) {
  
  // copy from eventOutputDataGPU_ to ebRecHits/eeRecHits
  cudaCheck( cudaMemcpyAsync(ebRecHits.energy.data(),
                             eventOutputDataGPU_.energy,
                             ebRecHits.energy.size() * sizeof(ecal::reco::StorageScalarType),
                             cudaMemcpyDeviceToHost,
                             cudaStream.id()) );
  cudaCheck( cudaMemcpyAsync(eeRecHits.energy.data(),
                             eventOutputDataGPU_.energy + ebRecHits.energy.size(),
                             eeRecHits.energy.size() * sizeof(ecal::reco::StorageScalarType),
                             cudaMemcpyDeviceToHost,
                             cudaStream.id()) );
  
    
}
  
  

    
DEFINE_FWK_MODULE(EcalRecHitProducerGPU);

