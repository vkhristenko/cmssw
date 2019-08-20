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

// conditions cpu
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/DataRecord/interface/EcalLinearCorrectionsRcd.h"


// conditions gpu
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalADCToGeVConstantGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalIntercalibConstantsGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalChannelStatusGPU.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalLaserAPDPNRatiosGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalLaserAPDPNRatiosRefGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalLaserAlphasGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalLinearCorrectionsGPU.h"


// configuration
#include "CommonTools/Utils/interface/StringToEnumValue.h"


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
  
private:
 
  // data
  uint32_t neb_, nee_; // extremely important, in particular neb_
  
  
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
  ecal::rechit::ConfigurationParameters configParameters_;
  uint32_t maxNumberHits_;
  
  
  // conditions handles
  edm::ESHandle<EcalADCToGeVConstantGPU>    ADCToGeVConstantHandle_;
  edm::ESHandle<EcalIntercalibConstantsGPU> IntercalibConstantsHandle_;
  edm::ESHandle<EcalChannelStatusGPU>       ChannelStatusHandle_;
  
  edm::ESHandle<EcalLaserAPDPNRatiosGPU>    LaserAPDPNRatiosHandle_;
  edm::ESHandle<EcalLaserAPDPNRatiosRefGPU> LaserAPDPNRatiosRefHandle_;
  edm::ESHandle<EcalLaserAlphasGPU>         LaserAlphasHandle_;
  edm::ESHandle<EcalLinearCorrectionsGPU>   LinearCorrectionsHandle_;

  // configuration
  std::vector<int> v_chstatus_;
  
  
  
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
  
  
  //---- db statuses to be exluded from reconstruction
  v_chstatus_ = StringToEnumValue<EcalChannelStatusCode::Code>( ps.getParameter<std::vector<std::string> >("ChannelStatusToBeExcluded"));
  
  
  
  
  // max number of digis to allocate for
  maxNumberHits_ = ps.getParameter<uint32_t>("maxNumberHits");
  
  // allocate event output data
  eventOutputDataGPU_.allocate(configParameters_, maxNumberHits_);
//   eventOutputDataGPU_.allocate(maxNumberHits_);
  
  
//   configParameters_.timeFitLimitsFirstEE = EEtimeFitLimits.first;
//   configParameters_.timeFitLimitsSecondEE = EEtimeFitLimits.second;
  
  
  configParameters_.ChannelStatusToBeExcludedSize = v_chstatus_.size();
  
  cudaCheck( cudaMalloc((void**)&configParameters_.ChannelStatusToBeExcluded, 
                        sizeof(int) * v_chstatus_.size()) 
           );
  cudaCheck( cudaMemcpy(configParameters_.ChannelStatusToBeExcluded, 
                        v_chstatus_.data(),
                        v_chstatus_.size() * sizeof(int),
                        cudaMemcpyHostToDevice) );
  
  
}


EcalRecHitProducerGPU::~EcalRecHitProducerGPU() {
  
  // free event ouput data 
  eventOutputDataGPU_.deallocate(configParameters_);
//   eventOutputDataGPU_.deallocate();
  
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

  // conditions
  // - laser correction 
  // - IC
  // - adt2gev
  

//   
  setup.get<EcalADCToGeVConstantRcd>()   .get(ADCToGeVConstantHandle_);
  setup.get<EcalIntercalibConstantsRcd>().get(IntercalibConstantsHandle_);
  setup.get<EcalChannelStatusRcd>()      .get(ChannelStatusHandle_);
  
  setup.get<EcalLaserAPDPNRatiosRcd>()     .get(LaserAPDPNRatiosHandle_);
  setup.get<EcalLaserAPDPNRatiosRefRcd>()  .get(LaserAPDPNRatiosRefHandle_);
  setup.get<EcalLaserAlphasRcd>()          .get(LaserAlphasHandle_);
  setup.get<EcalLinearCorrectionsRcd>()    .get(LinearCorrectionsHandle_);
  
//   

  auto const& ADCToGeVConstantProduct    = ADCToGeVConstantHandle_    -> getProduct(ctx.stream());
  auto const& IntercalibConstantsProduct = IntercalibConstantsHandle_ -> getProduct(ctx.stream());
  auto const& ChannelStatusProduct       = ChannelStatusHandle_       -> getProduct(ctx.stream());
  
  auto const& LaserAPDPNRatiosProduct     = LaserAPDPNRatiosHandle_     -> getProduct(ctx.stream());
  auto const& LaserAPDPNRatiosRefProduct  = LaserAPDPNRatiosRefHandle_  -> getProduct(ctx.stream());
  auto const& LaserAlphasProduct          = LaserAlphasHandle_          -> getProduct(ctx.stream());
  auto const& LinearCorrectionsProduct    = LinearCorrectionsHandle_    -> getProduct(ctx.stream());
  
  
  // bundle up conditions
  ecal::rechit::ConditionsProducts conditions {
      ADCToGeVConstantProduct,
      IntercalibConstantsProduct,
      ChannelStatusProduct,
//       
      LaserAPDPNRatiosProduct,   
      LaserAPDPNRatiosRefProduct,
      LaserAlphasProduct,        
      LinearCorrectionsProduct,  
//       
      IntercalibConstantsHandle_->getOffset()
  };
  
  
  //
  // schedule algorithms
  //
  
  edm::TimeValue_t event_time = event.time().value();
  
  
  ecal::rechit::create_ecal_rehit(
    inputDataGPU,
    eventOutputDataGPU_,
//     eventDataForScratchGPU_,
    conditions,  
    configParameters_,
    offsetForInput,
    event_time,
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
  eeRecHits.energy   += neb_;
  eeRecHits.chi2     += neb_;
  eeRecHits.did      += neb_;
  eeRecHits.time     += neb_;
  eeRecHits.flagBits += neb_;
  eeRecHits.extra    += neb_;
 
  // put into the event
  ctx.emplace(event, recHitsTokenEB_, std::move(ebRecHits));
  ctx.emplace(event, recHitsTokenEE_, std::move(eeRecHits));
  
}



    
DEFINE_FWK_MODULE(EcalRecHitProducerGPU);

