// -*- C++ -*-
//
// Package:    CUDA/DumpRawPixelDataGPUExtWork
// Class:      DumpRawPixelDataGPUExtWork
// 
/**\class DumpRawPixelDataGPUExtWork DumpRawPixelDataGPUExtWork.cc CUDA/DumpRawPixelDataGPUExtWork/plugins/DumpRawPixelDataGPUExtWork.cc

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
#include <thread>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Concurrency/interface/WaitingTaskWithArenaHolder.h"

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
#include "CUDA/PixelGPU/interface/kernels.h"

#include <cuda.h>
#include <cuda_runtime.h>

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

// just a safe observation
constexpr int NUM_PIXELS_PER_FED = 2000;

//
// DUmmy Producer to dump the RAW Data for the Pixel FEDs
//
class DumpRawPixelDataGPUExtWork : public edm::stream::EDProducer<> {
   public:
      // some type aliasing
      using Word64 = unsigned long;
      using Word32 = unsigned int;
      using DataWord = Word32;
      using DataType = testpixel::DigiFrame;
      using Product = std::vector<DataType>;

      explicit DumpRawPixelDataGPUExtWork(const edm::ParameterSet&);
      ~DumpRawPixelDataGPUExtWork();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      virtual void produce(edm::Event&, const edm::EventSetup&) override;
 //     virtual void acquire(edm::Event const&, edm::EventSetup const&, 
 //                          edm::WaitingTaskWithArenaHolder) override; 
   private:
      edm::EDGetTokenT<FEDRawDataCollection> m_tRawCollection;
      std::vector<unsigned int> m_fedIds;

      cudaStream_t m_stream;
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
DumpRawPixelDataGPUExtWork::DumpRawPixelDataGPUExtWork(const edm::ParameterSet& iConfig)
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

    // number of feds
    int size = FEDNumbering::MAXSiPixelFEDID - FEDNumbering::MINSiPixelFEDID + 1 + 
        FEDNumbering::MAXSiPixeluTCAFEDID - FEDNumbering::MINSiPixeluTCAFEDID + 1 + 
        FEDNumbering::MAXSiPixel2nduTCAFEDID - FEDNumbering::MINSiPixel2nduTCAFEDID + 1;

    // create a stream
    cudaStreamCreate(&m_stream);

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
    produces<std::vector<testpixel::DigiFrame> >("PixelDigisGPUExtWork") ; 
}


DumpRawPixelDataGPUExtWork::~DumpRawPixelDataGPUExtWork()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// produce
//
//void DumpRawPixelDataGPUExtWork::produce(edm::Event& iEvent,
//                                  edm::EventSetup const& iSetup) {
// iEvent.put(std::make_unique<Product>(), "PixelDigisGPUExtWork");
//}

//
// acquire
//
void
DumpRawPixelDataGPUExtWork::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
 //                            edm::WaitingTaskWithArenaHolder holder) {
    // get the collection
    edm::Handle<FEDRawDataCollection> hRawCollection;
    iEvent.getByToken(m_tRawCollection, hRawCollection);

    // launch the thread
//    std::thread(
//        [holder, &hRawCollection, this] {

            // some initialization
            int totalNumWords = 0;
            std::vector<std::tuple<DataWord*, DataType*, int>> dataRecords;

            // enqueue a kernel and input/output data transfers for each fed
            for (auto fid : m_fedIds) {
                // get he RAW Data for the given FED
                FEDRawData const& rawData = hRawCollection->FEDData(fid);


                // skip if no data
                if (rawData.size() == 0) continue;

                // get the number of 64 bit words and a pointer to the beginning
                int nWords64 = rawData.size() / sizeof(Word64);
                unsigned char const *data = rawData.data();

                // get the header
                Word64 const* header = reinterpret_cast<Word64 const*>(data);
                FEDHeader fedHeader(data);
                
                // get the trailer 
                Word64 const *trailer = reinterpret_cast<Word64 const*>(data) + nWords64 - 1;
                FEDTrailer fedTrailer(reinterpret_cast<unsigned char const*>(trailer));

                // get data buffer and size
                int numData = 2*(nWords64 - 1);
                DataWord const* firstWord = (DataWord const*)(header+1);

                // allocate memory on the GPUExtWork for input/output
                DataWord *d_data;
                DataType *d_digis;
                cudaMalloc(&d_data, numData * sizeof(DataWord));
                cudaMalloc(&d_digis, numData * sizeof(DataType));

                // transfer input data
                cudaMemcpyAsync(d_data, firstWord, numData * sizeof(DataWord),
                    cudaMemcpyHostToDevice, m_stream);

                // enqueue the kernel
                testpixel::wrap_raw2digi_simple(m_stream, d_data, d_digis, numData);

                // remember pointers
                dataRecords.push_back(std::make_tuple(d_data, d_digis, numData));

                totalNumWords += numData;
            }

            // synch that everyone is finished
            cudaStreamSynchronize(m_stream);

            // transfer the results back and free defice memory
            Product digis(totalNumWords);
            Product::iterator currentDigi = digis.begin();
            for (auto const& t : dataRecords) {
                // extract the saved pointers
                DataWord *d_data; 
                DataType *d_digis;
                int size;
                std::tie(d_data, d_digis, size) = t;

                // transfer the data from device back to the host
                cudaMemcpyAsync(&(*currentDigi), d_digis, sizeof(DataType) * size,
                        cudaMemcpyDeviceToHost, m_stream);
                
                // move the pointer further
                currentDigi += size;

                // release memory on the device
                cudaFree(d_data);
                cudaFree(d_digis);
            }

            // sync that the copy is finished
            cudaStreamSynchronize(m_stream);

    iEvent.put(std::make_unique<Product>(digis), "PixelDigisGPUExtWork");

            // signal to tbb that we are done with offloading
//            edm::WaitingTaskWithArenaHolder newh = std::move(holder);
//            std::exception_ptr exc;
//            newh.doneWaiting(exc);
 //       }
//    ).detach();
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*void
DumpRawPixelDataGPUExtWork::beginStream(edm::StreamID

{
}*/

// ------------ method called once each stream after processing all runs, lumis and events  ------------
/*void
DumpRawPixelDataGPUExtWork::endStream() {
}*/

// ------------ method called when starting to processes a run  ------------
/*
void
DumpRawPixelDataGPUExtWork::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DumpRawPixelDataGPUExtWork::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DumpRawPixelDataGPUExtWork::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DumpRawPixelDataGPUExtWork::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DumpRawPixelDataGPUExtWork::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DumpRawPixelDataGPUExtWork);
