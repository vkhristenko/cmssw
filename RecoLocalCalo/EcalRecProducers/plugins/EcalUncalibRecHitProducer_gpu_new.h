#ifndef RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitProducerGPUNew_hh
#define RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitProducerGPUNew_hh

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"

#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerBaseClass.h"


class EBDigiCollection;
class EEDigiCollection;

class EcalUncalibRecHitProducerGPUNew : public edm::stream::EDProducer<> {

        public:
                explicit EcalUncalibRecHitProducerGPUNew(const edm::ParameterSet& ps);
                ~EcalUncalibRecHitProducerGPUNew() override;
                void produce(edm::Event& evt, const edm::EventSetup& es) override;
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        private:

		edm::EDGetTokenT<EBDigiCollection>  ebDigiCollectionToken_; 
                edm::EDGetTokenT<EEDigiCollection>  eeDigiCollectionToken_; 

                std::string ebHitCollection_; 
                std::string eeHitCollection_; 
		
                EcalUncalibRecHitWorkerBaseClass * worker_;
};
#endif
