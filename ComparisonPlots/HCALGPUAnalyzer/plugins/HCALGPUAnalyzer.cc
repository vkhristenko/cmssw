// -*- C++ -*-
//
// Package:    ComparisonPlots/HCALGPUAnalyzer
// Class:      HCALGPUAnalyzer
//
/**\class HCALGPUAnalyzer HCALGPUAnalyzer.cc ComparisonPlots/HCALGPUAnalyzer/plugins/HCALGPUAnalyzer.cc

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

#include "TH2F.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

class HCALGPUAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HCALGPUAnalyzer(const edm::ParameterSet&);
  ~HCALGPUAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //  void ClearVariables();
  
  // some variables for storing information
  double Method0Energy, Method0EnergyGPU;
  double RecHitEnergy, RecHitEnergyGPU;
  double RecHitTime, RecHitTimeGPU;
  double iEta, iEtaGPU;
  double iPhi, iPhiGPU; 
  int depth, depthGPU;
  
  TH2F *hCheckEnergy_2dMahi;
  TH2F *hCheckEnergy_2dM0;
  TH2F *hCheckTime_2dMahi;

  TH2F *Unmatched;
  TH2F *Matched;
  TH1F *hCheckEnergy_cpu;
  TH1F *hCheckEnergy_gpu;
  TH1F *hCheckEnergyM0_cpu;
  TH1F *hCheckEnergyM0_gpu;
  TH1F *hCheckTime_cpu;
  TH1F *hCheckTime_gpu;   


  // create the output file
  edm::Service<TFileService> FileService;
  // create the token to retrieve hit information
  edm::EDGetTokenT<HBHERecHitCollection>    hRhToken;
  edm::EDGetTokenT<HBHERecHitCollection> hRhTokenGPU;
  
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
HCALGPUAnalyzer::HCALGPUAnalyzer(const edm::ParameterSet& iConfig)
// :
  //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{

  usesResource("TFileService");  

  hRhToken = consumes<HBHERecHitCollection >(iConfig.getUntrackedParameter<string>("HBHERecHits","hbheprereco"));
  hRhTokenGPU = consumes<HBHERecHitCollection >(iConfig.getUntrackedParameter<string>("HBHERecHits","hbheprerecogpu"));

  hCheckEnergy_2dM0 = FileService->make<TH2F>("hCheckEnergy_2dM0","hCheckEnergy_2dM0", 1000, 0., 100. , 1000, 0., 100.);
  hCheckEnergy_2dM0->GetXaxis()->SetTitle("Cpu M0 Energy");
  hCheckEnergy_2dM0->GetYaxis()->SetTitle("GPU M0 Energy");

  hCheckEnergy_2dMahi = FileService->make<TH2F>("hCheckEnergy_2dMahi","hCheckEnergy_2dMahi", 1000, 0., 100. , 1000, 0., 100.);
  hCheckEnergy_2dMahi->GetXaxis()->SetTitle("CPU Energy");
  hCheckEnergy_2dMahi->GetYaxis()->SetTitle("GPU Energy");

  hCheckTime_2dMahi = FileService->make<TH2F>("hCheckTime_2dMahi","hCheckTime_2dMahi", 250, -12.5, 12.5 , 250, -12.5, 12.5);
  hCheckTime_2dMahi->GetXaxis()->SetTitle("Mahi Time CPU");
  hCheckTime_2dMahi->GetYaxis()->SetTitle("Mahi Time GPU");

  //

  hCheckEnergyM0_cpu = FileService->make<TH1F>("hCheckEnergyM0_cpu","hCheckEnergyM0_cpu", 100, 0., 100.);
  hCheckEnergyM0_cpu->GetXaxis()->SetTitle("CPU Energy");

  hCheckEnergy_cpu = FileService->make<TH1F>("hCheckEnergy_cpu","hCheckEnergy_cpu", 50, 0., 50.);
  hCheckEnergy_cpu->GetXaxis()->SetTitle("CPU Energy");

  hCheckEnergy_gpu = FileService->make<TH1F>("hCheckEnergy_gpu","hCheckEnergy_gpu", 50, 0., 50.);
  hCheckEnergy_gpu->GetXaxis()->SetTitle("GPU Energy");

  //

  hCheckTime_cpu = FileService->make<TH1F>("hCheckTime_cpu","hCheckTime_cpu", 50, -25., 25.);
  //  hCheckTime_cpu = FileService->make<TH1F>("hCheckTime_cpu","hCheckTime_cpu", 50, 0., 25.);
  hCheckTime_cpu->GetXaxis()->SetTitle("CPU Time");

  hCheckTime_gpu = FileService->make<TH1F>("hCheckTime_gpu","hCheckTime_gpu", 50, -25., 25.);
  //  hCheckTime_gpu = FileService->make<TH1F>("hCheckTime_gpu","hCheckTime_gpu", 50, 0., 25.);
  hCheckTime_gpu->GetXaxis()->SetTitle("GPU Time");

  Unmatched= FileService->make<TH2F>("Unmatched","Unmatched (eta,phi)",100,-50.,50.,85,0.,85.);
  Matched= FileService->make<TH2F>("Matched","Matched (eta,phi)",100,-50.,50.,85,0.,85.);

  //now do what ever initialization is needed

}


HCALGPUAnalyzer::~HCALGPUAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HCALGPUAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   // Read events
   Handle<HBHERecHitCollection> hRecHits; // create handle
   iEvent.getByToken(hRhToken, hRecHits); // get events based on token

   Handle<HBHERecHitCollection> hRecHitsGPU;
   iEvent.getByToken(hRhTokenGPU, hRecHitsGPU);
  
   // Loop over all rechits in one event
   for(int i = 0; i < (int)hRecHits->size(); i++) {
     //     if (i==0) cout << "CPU RecHits ---------------------";
     //     ClearVariables(); // sets a bunch of stuff to zero
    
     // get ID information for the reconstructed hit
     HcalDetId detID_rh = (*hRecHits)[i].id().rawId();
    
     // ID information can get us detector coordinates
     depth = (*hRecHits)[i].id().depth();
     iEta = detID_rh.ieta();
     iPhi = detID_rh.iphi();
    
     // get some variables
     Method0Energy = (*hRecHits)[i].eraw();
     RecHitEnergy = (*hRecHits)[i].energy();
     RecHitTime = (*hRecHits)[i].time();

     hCheckEnergy_cpu->Fill(RecHitEnergy);
     hCheckTime_cpu->Fill(RecHitTime);

     /*
     cout << "Run " << i << ": ";
     cout << "Method0Energy: " << Method0Energy;
     cout << "RecHitEnergy: " << RecHitEnergy;
     cout << "depth: " << depth;
     cout << "iEta: " << iEta;
     cout << "iPhi: " << iPhi;
     cout << "RecHitTime" << RecHitTime;
     */

   }

   for(int i = 0; i < (int)hRecHitsGPU->size(); i++) {
     //     if (i==0) cout << "GPU RecHits ----------------------"
     //		 ClearVariables(); // sets a bunch of stuff to zero
    
     // get ID information for the reconstructed hit
     HcalDetId detID_rh = (*hRecHitsGPU)[i].id().rawId();
    
     // ID information can get us detector coordinates
     depthGPU = (*hRecHitsGPU)[i].id().depth();
     iEtaGPU = detID_rh.ieta();
     iPhiGPU = detID_rh.iphi();
    
     // get some variables
     Method0EnergyGPU = (*hRecHitsGPU)[i].eraw();
     RecHitEnergyGPU = (*hRecHitsGPU)[i].energy();
     RecHitTimeGPU = (*hRecHitsGPU)[i].time();

     hCheckEnergy_gpu->Fill(RecHitEnergyGPU);
     hCheckTime_gpu->Fill(RecHitTimeGPU);

     /*
     cout << "Run " << i << ": ";
     cout << "Method0Energy: " << Method0EnergyGPU;
     cout << "RecHitEnergy: " << RecHitEnergyGPU;
     cout << "depth: " << depthGPU;
     cout << "iEta: " << iEtaGPU;
     cout << "iPhi: " << iPhiGPU;
     cout << "RecHitTime" << RecHitTimeGPU;
     */

   }


   // Loop over all rechits in one event
   for(int i = 0; i < (int)hRecHits->size(); i++) {

     HcalDetId detID_rh = (*hRecHits)[i].id().rawId();

     bool unmatched=true;
     //     cout << "--------------------------------------------------------" << endl; 

     for(int j = 0; j < (int)hRecHitsGPU->size(); j++) {

       HcalDetId detID_gpu = (*hRecHitsGPU)[j].id().rawId();

      
       if ( (detID_rh == detID_gpu)) {

	 /*
	 cout << "Mtime(cpu)" << (*hRecHits)[i].time() << endl; 
	 cout << "     Mtime(gpu)" << (*hRecHitsGPU)[j].time() << endl;

	 cout << "ME(cpu)" << (*hRecHits)[i].energy() << endl; 
	 cout << "     ME(gpu)" << (*hRecHitsGPU)[j].energy() << endl;

	 cout << "M0E(cpu)" << (*hRecHits)[i].eraw() << endl; 
	 cout << "     M0E(gpu)" << (*hRecHitsGPU)[j].eraw() << endl;
	 */

	 hCheckEnergy_2dM0->Fill((*hRecHits)[i].eraw(),(*hRecHitsGPU)[j].eraw());
	 hCheckEnergy_2dMahi->Fill((*hRecHits)[i].energy(),(*hRecHitsGPU)[j].energy());
	 hCheckTime_2dMahi->Fill((*hRecHits)[i].time(),(*hRecHitsGPU)[j].time());


	 Matched->Fill(detID_rh.ieta(),detID_rh.iphi());

	 unmatched=false;

       }

     }

     ///

     if(unmatched) { 
       Unmatched->Fill(detID_rh.ieta(),detID_rh.iphi());
       //       cout << "   recHit not matched ="  << detID_rh << "  E(raw)=" << (*hRecHits)[i].eraw() << " E=" << (*hRecHits)[i].energy() << endl;
     }

   }



}


// ------------ method called once each job just before starting event loop  ------------
void
HCALGPUAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HCALGPUAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HCALGPUAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
HCALGPUAnalyzer::ClearVariables(){
  RecHitEnergy = 0;
  depth=0;
  iEta = 0;
  iPhi = 0;
  RecHitTime = 0;
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(HCALGPUAnalyzer);
