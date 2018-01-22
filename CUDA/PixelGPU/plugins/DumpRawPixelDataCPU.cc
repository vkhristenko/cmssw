// -*- C++ -*-
//
// Package:    CUDA/DumpRawPixelDataCPU
// Class:      DumpRawPixelDataCPU
// 
/**\class DumpRawPixelDataCPU DumpRawPixelDataCPU.cc CUDA/DumpRawPixelDataCPU/plugins/DumpRawPixelDataCPU.cc

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
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/src/fed_header.h"

#include "CUDA/DataFormats/interface/DigiFrame.h"

#define DONOTPRINT

// some constants
constexpr int BITS_LINK = 6;
constexpr int BITS_ROC = 5;
constexpr int BITS_DCOL = 5;
constexpr int BITS_PIXEL = 8;
constexpr int BITS_ADC = 8;

constexpr int SHIFT_LINK = BITS_ROC + BITS_DCOL + BITS_PIXEL + BITS_ADC;
constexpr int SHIFT_ROC = BITS_DCOL + BITS_PIXEL + BITS_ADC;
constexpr int SHIFT_DCOL = BITS_PIXEL + BITS_ADC;
constexpr int SHIFT_PIXEL = BITS_ADC;

constexpr int MASK_LINK = 0x0000003F;
constexpr int MASK_ROC = 0x0000001F;
constexpr int MASK_DCOL = MASK_ROC;
constexpr int MASK_PIXEL = 0x000000FF;
constexpr int MASK_ADC = MASK_PIXEL;

//
// DUmmy Producer to dump the RAW Data for the Pixel FEDs
//
class DumpRawPixelDataCPU : public edm::stream::EDProducer<> {
   public:
      // some type aliasing
      using Word64 = unsigned long;
      using Word32 = unsigned int;
      using DataWord = Word32;
      using Product = std::vector<testpixel::DigiFrame>;

      explicit DumpRawPixelDataCPU(const edm::ParameterSet&);
      ~DumpRawPixelDataCPU();

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
DumpRawPixelDataCPU::DumpRawPixelDataCPU(const edm::ParameterSet& iConfig)
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
    produces<std::vector<testpixel::DigiFrame> >("PixelDigisCPU") ; 
}


DumpRawPixelDataCPU::~DumpRawPixelDataCPU()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DumpRawPixelDataCPU::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // extract the fed raw collection from the event
   edm::Handle<FEDRawDataCollection> hRawCollection;
   iEvent.getByToken(m_tRawCollection, hRawCollection);

   // initialize the collection to be put into iEvent
   Product digis;

   for (auto fid : m_fedIds) {
       // get the RAW Data for the given FED
       FEDRawData const& rawData = hRawCollection->FEDData(fid);

       // skip FEDs without any raw buffers
       if (rawData.size()==0) continue;
       
       int nWord64s = rawData.size() / sizeof(Word64);
#ifndef DONOTPRINT
       printf("----------------------------fed: %d byte length: %lu Word64 length: %d--------------------------\n",
              fid, reinterpret_cast<long unsigned int>(rawData.size()), nWord64s);
#endif

       // dump the whole buffer without any interpretation
       unsigned char const *data = rawData.data();
#ifndef DONOTPRINT
       printf("dump the whole buffer:\n");
       for (size_t ic=0; ic<rawData.size(); ic++) {
           if (ic % 20 == 0)
               printf("\n");
           printf("%#x ", data[ic]);
       }
       printf("\n\n");
       
#endif
       // let's try to interpret the buffer
#ifndef DONOTPRINT
       printf("interpreting the raw byte buffer:\n");
#endif

       // header
#ifndef DONOTPRINT
       printf("header:\n");
#endif
       Word64 const *header = reinterpret_cast<Word64 const*>(data);
       FEDHeader fedHeader(reinterpret_cast<unsigned char const *>(header));
#ifndef DONOTPRINT
       printf("\tsourceId = %d\n", fedHeader.sourceID());
#endif
       fedh_t const * test = reinterpret_cast<fedh_t const*>(data);
#ifndef DONOTPRINT
       printf("\ttest sourceId = %d\n", test->sourceid >> 8);
       printf("\tmoreHeaders = %s\n", fedHeader.moreHeaders() ? "yes" : "no");
#endif

       // trailer words
#ifndef DONOTPRINT
       printf("trailer:\n");
#endif
       Word64 const *trailerWord = reinterpret_cast<Word64 const*>(data) + 
           nWord64s - 1;
       FEDTrailer fedTrailer(reinterpret_cast<unsigned char const*>(trailerWord));
#ifndef DONOTPRINT
       printf("\tfragmentLength = %d\n", fedTrailer.fragmentLength());
       printf("\tmoreTrailers = %s\n", fedTrailer.moreTrailers() ? "yes" : "no");
#endif

       // data
#ifndef DONOTPRINT
       printf("data:\n");
#endif
       int numData32Words = 2*(nWord64s - 2); // - 4 32-bit (2 64-bit) h/t words
       DataWord const* firstWord = (DataWord const*)(header+1);
       DataWord const* lastWord = ((DataWord const*)(trailerWord-1))+1;
       // the last word might be 0 in case of number of pixels being odd
//       if (*lastWord == 0) 
//           lastWord--;
       int iWord = 0;
       for (auto dataWord = firstWord; dataWord <= lastWord; dataWord++) {
#ifndef DONOTPRINT
           printf("\tword %d raw data: %#x\n", iWord, *dataWord);
#endif
 //          if (*dataWord == 0) continue;

           // unpack 1 pixel data/coordinates
           int link = (*dataWord >> SHIFT_LINK) & MASK_LINK;
           int roc = (*dataWord >> SHIFT_ROC) & MASK_ROC;
           int dcol = (*dataWord >> SHIFT_DCOL) & MASK_DCOL;
           int pixel = (*dataWord >> SHIFT_PIXEL) & MASK_PIXEL;
           int adc = (*dataWord) & MASK_ADC;
//           DigiFrame frame {link, roc, dcol, pixel, adc};
           digis.emplace_back(link, roc, dcol, pixel, adc);
#ifndef DONOTPRINT
           printf("\tunpacked data:\n");
           printf("\t\tlink = %d roc = %d dcol = %d pixel = %d adc = %d\n", 
                  link, roc, dcol, pixel, adc);
#endif

           // just to keep track
           iWord++;
       }
   }

   //
   // put the unpacked digis into the collection
   //
   iEvent.put(std::make_unique<Product>(digis), "PixelDigisCPU");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*void
DumpRawPixelDataCPU::beginStream(edm::StreamID)
{
}*/

// ------------ method called once each stream after processing all runs, lumis and events  ------------
/*void
DumpRawPixelDataCPU::endStream() {
}*/

// ------------ method called when starting to processes a run  ------------
/*
void
DumpRawPixelDataCPU::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DumpRawPixelDataCPU::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DumpRawPixelDataCPU::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DumpRawPixelDataCPU::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DumpRawPixelDataCPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DumpRawPixelDataCPU);
