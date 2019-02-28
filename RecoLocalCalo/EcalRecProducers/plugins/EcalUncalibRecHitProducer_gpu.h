#ifndef RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitProducerGPU_hh
#define RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitProducerGPU_hh

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"

#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerBaseClass.h"


class EBDigiCollection;
class EEDigiCollection;

class EcalUncalibRecHitProducerGPU : public edm::stream::EDProducer<> {

        public:
                explicit EcalUncalibRecHitProducerGPU(const edm::ParameterSet& ps);
                ~EcalUncalibRecHitProducerGPU() override;
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
