// -*- C++ -*-
//
// Package:    ComparisonPlots/cmp_ecal_uncalibrechits
// Class:      cmp_ecal_uncalibrechits
//
/**\class cmp_ecal_uncalibrechits cmp_ecal_uncalibrechits.cc ComparisonPlots/cmp_ecal_uncalibrechits/plugins/cmp_ecal_uncalibrechits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mariarosaria D'Alfonso
//         Created:  Mon, 17 Dec 2018 16:22:58 GMT
//
//


// system include files                                                                                                                                                                                                                                               
#include <memory>
#include <string>
#include <map>
#include <iostream>
using namespace std;


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"


//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/Framework/interface/Event.h"

//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "TH2F.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

#include "common.h"

class cmp_ecal_uncalibrechits : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit cmp_ecal_uncalibrechits(const edm::ParameterSet&);
  ~cmp_ecal_uncalibrechits();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //  void ClearVariables();
  
    //
    comparator<100, 0, 100> cmp_amplitude;
    comparator<100, 0, 50> cmp_amplitude_error;
    comparator<100, 0, 50> cmp_pedestal;
    comparator<100, -20, 20> cmp_jitter;
    comparator<100, 0, 100> cmp_chi2;
    comparator<100, 0, 100> cmp_ootamplitudes[EcalDataFrame::MAXSAMPLES];
    comparator<100, 0, 100> cmp_ootchi2;


  // create the output file
  edm::Service<TFileService> fs;
  // create the token to retrieve hit information
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> tok_cpu;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> tok_gpu;
  
  // crap for trouble-shooting, create a TCanvas here to print out pulse shapes of problem channels
  //  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  int nProblems = 0;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
cmp_ecal_uncalibrechits::cmp_ecal_uncalibrechits(const edm::ParameterSet& iConfig)
// :
  //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
  usesResource("TFileService");  

    cmp_amplitude.create("amplitude", fs);
    cmp_amplitude_error.create("amplitude_error", fs);
    cmp_pedestal.create("pedestal", fs);
    cmp_jitter.create("jitter", fs);
    cmp_chi2.create("chi2", fs);
    for (unsigned int i=0; i<EcalDataFrame::MAXSAMPLES; i++) 
        cmp_ootamplitudes[i].create("ootamplitude_" + std::to_string(i), fs);
    cmp_ootchi2.create("ootchi2", fs);

  tok_cpu = consumes<EcalUncalibratedRecHitCollection>(
    iConfig.getParameter<edm::InputTag>("tok_cpu"));
  tok_gpu = consumes<EcalUncalibratedRecHitCollection>(
    iConfig.getParameter<edm::InputTag>("tok_gpu"));
}


cmp_ecal_uncalibrechits::~cmp_ecal_uncalibrechits()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
cmp_ecal_uncalibrechits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Read events
   Handle<EcalUncalibratedRecHitCollection> col_cpu; // create handle
   iEvent.getByToken(tok_cpu, col_cpu); // get events based on token

   Handle<EcalUncalibratedRecHitCollection> col_gpu;
   iEvent.getByToken(tok_gpu, col_gpu);

   for (unsigned int i=0; i<col_cpu->size(); i++) {
       auto const& rh = (*col_cpu)[i];
        cmp_amplitude.v1->Fill(rh.amplitude());
        cmp_amplitude_error.v1->Fill(rh.amplitudeError());
        cmp_pedestal.v1->Fill(rh.pedestal());
        cmp_jitter.v1->Fill(rh.jitter());
        cmp_chi2.v1->Fill(rh.chi2());
        for (unsigned int bx=0; bx<EcalDataFrame::MAXSAMPLES; bx++)
            cmp_ootamplitudes[bx].v1->Fill(rh.outOfTimeAmplitude(bx));
//        cmp_ootchi2.v1->Fill(rh.)
   }

   for (unsigned int i=0; i<col_gpu->size(); i++) {
       auto const& rh = (*col_gpu)[i];
        cmp_amplitude.v2->Fill(rh.amplitude());
        cmp_amplitude_error.v2->Fill(rh.amplitudeError());
        cmp_pedestal.v2->Fill(rh.pedestal());
        cmp_jitter.v2->Fill(rh.jitter());
        cmp_chi2.v2->Fill(rh.chi2());
        for (unsigned int bx=0; bx<EcalDataFrame::MAXSAMPLES; bx++)
            cmp_ootamplitudes[bx].v2->Fill(rh.outOfTimeAmplitude(bx));
   }

   for (unsigned int i=0; i<col_cpu->size(); i++) {
       auto const& hit_cpu = (*col_cpu)[i];
       for (unsigned int j=0; j<col_gpu->size(); j++) {
           auto const& hit_gpu = (*col_gpu)[j];
           if (hit_cpu.id() == hit_gpu.id()) {
               cmp_amplitude.correlation->Fill(
                    hit_cpu.amplitude(), hit_gpu.amplitude());
               cmp_amplitude.ratio->Fill(
                    hit_cpu.amplitude() / hit_gpu.amplitude());
               cmp_amplitude_error.correlation->Fill(
                    hit_cpu.amplitudeError(), hit_gpu.amplitudeError());
               cmp_amplitude_error.ratio->Fill(
                    hit_cpu.amplitudeError() / hit_gpu.amplitudeError());

               cmp_pedestal.correlation->Fill(
                    hit_cpu.pedestal(), hit_gpu.pedestal());
               cmp_pedestal.ratio->Fill(
                    hit_cpu.pedestal() / hit_gpu.pedestal());
               cmp_jitter.correlation->Fill(
                    hit_cpu.jitter(), hit_gpu.jitter());
               cmp_jitter.ratio->Fill(
                    hit_cpu.jitter() / hit_gpu.jitter());
               cmp_chi2.correlation->Fill(
                    hit_cpu.chi2(), hit_gpu.chi2());
               cmp_chi2.ratio->Fill(
                    hit_cpu.chi2() / hit_gpu.chi2());
           }
       }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void
cmp_ecal_uncalibrechits::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
cmp_ecal_uncalibrechits::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
cmp_ecal_uncalibrechits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


/*
void 
cmp_ecal_uncalibrechits::ClearVariables(){
  RecHitEnergy = 0;
  depth=0;
  iEta = 0;
  iPhi = 0;
  RecHitTime = 0;
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(cmp_ecal_uncalibrechits);
