#include "RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitProducer.h"

#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalUncalibRecHitWorkerFactory.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

HGCalUncalibRecHitProducer::HGCalUncalibRecHitProducer(const edm::ParameterSet& ps) :
  eeDigiCollection_( consumes<HGCalDigiCollection>( ps.getParameter<edm::InputTag>("HGCEEdigiCollection") ) ),
  hefDigiCollection_( consumes<HGCalDigiCollection>( ps.getParameter<edm::InputTag>("HGCHEFdigiCollection") ) ),
  hebDigiCollection_( consumes<HGCalDigiCollection>( ps.getParameter<edm::InputTag>("HGCHEBdigiCollection") ) ),
  eeHitCollection_( ps.getParameter<std::string>("HGCEEhitCollection") ),
  hefHitCollection_( ps.getParameter<std::string>("HGCHEFhitCollection") ),
  hebHitCollection_( ps.getParameter<std::string>("HGCHEBhitCollection") ),
  worker_{ HGCalUncalibRecHitWorkerFactory::get()->create(ps.getParameter<std::string>("algo"), ps) }
{

  produces< HGCeeUncalibratedRecHitCollection >(eeHitCollection_);
  produces< HGChefUncalibratedRecHitCollection >(hefHitCollection_);
  produces< HGChebUncalibratedRecHitCollection >(hebHitCollection_);

}

HGCalUncalibRecHitProducer::~HGCalUncalibRecHitProducer() {
}

void
HGCalUncalibRecHitProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  using namespace edm;
  
  
 
  // tranparently get things from event setup
  worker_->set(es);
  
  // prepare output
  auto eeUncalibRechits = std::make_unique<HGCeeUncalibratedRecHitCollection>();
  auto hefUncalibRechits = std::make_unique<HGChefUncalibratedRecHitCollection>();
  auto hebUncalibRechits = std::make_unique<HGChebUncalibratedRecHitCollection>();
  
  // loop over HGCEE digis
  edm::Handle< HGCalDigiCollection > pHGCEEDigis;
  evt.getByToken( eeDigiCollection_, pHGCEEDigis);
  const HGCalDigiCollection* eeDigis = 
    pHGCEEDigis.product(); // get a ptr to the product
  eeUncalibRechits->reserve(eeDigis->size());
  for(auto itdg = eeDigis->begin(); itdg != eeDigis->end(); ++itdg) {
    worker_->run1(evt, itdg, *eeUncalibRechits);
  }
  
  edm::Handle< HGCalDigiCollection > pHGCHEFDigis;
  evt.getByToken( hefDigiCollection_, pHGCHEFDigis);
  const HGCalDigiCollection* hefDigis = 
    pHGCHEFDigis.product(); // get a ptr to the product
  hefUncalibRechits->reserve(hefDigis->size());
  for(auto itdg = hefDigis->begin(); itdg != hefDigis->end(); ++itdg) {
    worker_->run2(evt, itdg, *hefUncalibRechits);
  }

  edm::Handle< HGCalDigiCollection > pHGCHEBDigis;
  evt.getByToken( hebDigiCollection_, pHGCHEBDigis);
  const HGCalDigiCollection* hebDigis = 
    pHGCHEBDigis.product(); // get a ptr to the product
  hebUncalibRechits->reserve(hebDigis->size());
  for(auto itdg = hebDigis->begin(); itdg != hebDigis->end(); ++itdg) {
    worker_->run3(evt, itdg, *hebUncalibRechits);
  }
  
  // put the collection of recunstructed hits in the event
  evt.put(std::move(eeUncalibRechits), eeHitCollection_);
  evt.put(std::move(hefUncalibRechits), hefHitCollection_);
  evt.put(std::move(hebUncalibRechits), hebHitCollection_);
}

#include "FWCore/Framework/interface/MakerMacros.h"                                                                                                            
DEFINE_FWK_MODULE( HGCalUncalibRecHitProducer );
