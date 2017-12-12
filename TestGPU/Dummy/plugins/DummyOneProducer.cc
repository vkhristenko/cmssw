// -*- C++ -*-
//
// Package:    TestGPU/Dummy
// Class:      DummyOneProducer
// 
/**\class DummyOneProducer DummyOneProducer.cc TestGPU/Dummy/plugins/DummyOneProducer.cc

 Description: Dummy EDM One Producer with GPU offload

*/
//
// Original Author:  Viktor Khristenko
//         Created:  Tue, 14 Nov 2017 10:33:24 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// use the new product
#include "TestGPU/Dummy/interface/Vector.h"
#include "TestGPU/Dummy/interface/gpu_kernels.h"

#include <cuda.h>
#include <cuda_runtime.h>

//
// class declaration
//
class DummyOneProducer : public edm::one::EDProducer<> {
public:
    // some type aliasing
    using DataType = int;

    // ctor and dtor
    explicit DummyOneProducer(const edm::ParameterSet&);
    ~DummyOneProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void beginJob() override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    cudaStream_t m_stream;
    cudaEvent_t m_estart, m_estop;

    int m_size;
    DataType *m_ha, *m_hb, *m_hc;
    DataType *m_da, *m_db, *m_dc;

};

//
// constructors and destructor
//
DummyOneProducer::DummyOneProducer(const edm::ParameterSet& iConfig)
{
    //
    // get the size of vectors to be used
    //
    m_size = iConfig.getUntrackedParameter<int>("size", 1000);

    //
    // Initialize the start/stop  Events
    //
    cudaEventCreate(&m_estart);
    cudaEventCreate(&m_estop);
    
    //
    // Initialize a cuda stream
    //
    cudaStreamCreate(&m_stream);

    //
    // Perform the memory allocs only once per single edm::producer!
    //
    cudaMalloc(&m_da, m_size * sizeof(DataType));
    cudaMalloc(&m_db, m_size * sizeof(DataType));
    cudaMalloc(&m_dc, m_size * sizeof(DataType));

    //
    // Let the framework know that we are going to put Vector into the event
    //
    produces<testgpu::Vector<DataType> >("VectorForGPU");

/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


DummyOneProducer::~DummyOneProducer()
{
    //
    // free memory on the host side
    //
    delete [] m_ha;
    delete [] m_hb;
    delete [] m_hc;

    //
    // free just once at destruction
    //
    cudaFree(m_da);
    cudaFree(m_db);
    cudaFree(m_dc);

    // 
    // destroy Events and Ones
    //
    cudaEventDestroy(m_estart);
    cudaEventDestroy(m_estop);
    cudaStreamDestroy(m_stream);
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DummyOneProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //
    // initialize the values in the vector
    //
    for (auto i=0; i<m_size; i++) {
        m_ha[i] = i;
        m_hb[i] = i+i;
    }

    //
    // record the start of event
    //
    cudaEventRecord(m_estart, 0);

    // 
    // Perform memcpy
    //
    cudaMemcpyAsync(m_da, m_ha, m_size * sizeof(DataType),
                    cudaMemcpyHostToDevice, m_stream);
    cudaMemcpyAsync(m_db, m_hb, m_size * sizeof(DataType),
                    cudaMemcpyHostToDevice, m_stream);

    // 
    // call the kernel wrapper
    //
    testgpu::wrapperVectorAdd<DataType>(m_da, m_db, m_dc, m_stream, m_size);

    // 
    // copy the results
    //
    cudaMemcpyAsync(m_hc, m_dc, m_size * sizeof(DataType),
                    cudaMemcpyDeviceToHost, m_stream);

    // 
    // synch with the stream
    //
    cudaStreamSynchronize(m_stream);

    //
    // Record the stop
    //
    cudaEventRecord(m_estop, 0);
    cudaEventSynchronize(m_estop);
    float elapsedTime {0.0};
    cudaEventElapsedTime(&elapsedTime, m_estart, m_estop);
    printf("Time Taken: %3.1f ms\n", elapsedTime);

    //
    // put the computed vector into the event
    //
    testgpu::Vector<int> v;
    v.m_values = std::vector<int>(m_hc, m_hc + m_size);
    iEvent.put(std::make_unique<testgpu::Vector<int> >(v), "VectorForGPU");
}

void
DummyOneProducer::beginJob()
{
}

void
DummyOneProducer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DummyOneProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyOneProducer);
