// -*- C++ -*-
//
// Package:    CUDA/PixelGPU
// Class:      PixelGPU
// 
/**\class PixelGPU PixelGPU.cc CUDA/PixelGPU/plugins/PixelGPU.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Viktor Khristenko
//         Created:  Thu, 18 Jan 2018 09:20:14 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"


//
// class declaration
//

class PixelGPU : public edm::stream::EDProducer<> {
   public:
      explicit PixelGPU(const edm::ParameterSet&);
      ~PixelGPU();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      virtual void produce(edm::Event&, const edm::EventSetup&) override;
   private:
      edm::EDGetTokenT<FEDRawDataCollection> m_tRawCollection;
      std::vector<unsigned int> m_fedIds;
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
PixelGPU::PixelGPU(const edm::ParameterSet& iConfig)
{
    // get the input label for the raw colletion
    m_tRawCollection = consumes<FEDRawDataCollection>(
        iConfig.getParameter<edm::InputTag>("InputLabel"));

    // initialize for now...
    for (int i = FEDNumbering::MINSiPixelFEDID; 
         i <= FEDNumbering::MAXSiPixelFEDID; i++)
        m_fedIds.push_back(i);
    for (int i = FEDNumbering::MINSiPixeluTCAFEDID; 
         i <= FEDNumbering::MAXSiPixeluTCAFEDID; i++)
        m_fedIds.push_back(i);
    for (int i = FEDNumbering::MINSiPixel2nduTCAFEDID; 
         i <= FEDNumbering::MAXSiPixel2nduTCAFEDID; i++)
        m_fedIds.push_back(i);

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


PixelGPU::~PixelGPU()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PixelGPU::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // extract the fed raw collection from the event
   edm::Handle<FEDRawDataCollection> hRawCollection;
   iEvent.getByToken(m_tRawCollection, hRawCollection);

   for (auto fid : m_fedIds) {
       // get the RAW Data for the given FED
       FEDRawData const& rawData = hRawCollection->FEDData(fid);

       std::cout << "fed = " << fid << "  "
           << "buffer length = " << rawData.size() << std::endl;
       printf("dump the buffer:\n");
       unsigned char const *data = rawData.data();
       for (size_t ic=0; ic<rawData.size(); ic++) {
           if (ic % 20 == 0)
               printf("\n");
           printf("%#x ", data[ic]);
       }
       printf("\n");
   }
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*void
PixelGPU::beginStream(edm::StreamID)
{
}*/

// ------------ method called once each stream after processing all runs, lumis and events  ------------
/*void
PixelGPU::endStream() {
}*/

// ------------ method called when starting to processes a run  ------------
/*
void
PixelGPU::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PixelGPU::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PixelGPU::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PixelGPU::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PixelGPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelGPU);
