// -*- C++ -*-
//
// Package:    TestGPU/Dummy
// Class:      DummyStreamProducer
// 
/**\class DummyStreamProducer DummyStreamProducer.cc TestGPU/Dummy/plugins/DummyStreamProducer.cc

 Description: A simple Dummy EDM Stream Producer with GPU offload

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
#include "FWCore/Framework/interface/stream/EDProducer.h"

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
class DummyStreamProducer : public edm::stream::EDProducer<> {
public:
    // some type aliasing
    using DataType = int;

    // ctor and dtor
    explicit DummyStreamProducer(const edm::ParameterSet&);
    ~DummyStreamProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    cudaStream_t m_stream;
    cudaEvent_t m_estart, m_estop;

    bool m_isPinned;

    int m_size;
    DataType *m_ha, *m_hb, *m_hc;
    DataType *m_da, *m_db, *m_dc;

};

//
// constructors and destructor
//
DummyStreamProducer::DummyStreamProducer(const edm::ParameterSet& iConfig)
{
    //
    // get the size of vectors to be used
    //
    m_size = iConfig.getUntrackedParameter<int>("size");

    //
    // should we use pinned memory
    //
    m_isPinned = iConfig.getUntrackedParameter<bool>("isPinned");

    // 
    // allocate stuff on the host's side
    //
    if (m_isPinned) {
        cudaHostAlloc((void**)&m_ha, m_size * sizeof(DataType), cudaHostAllocDefault);
        cudaHostAlloc((void**)&m_hb, m_size * sizeof(DataType), cudaHostAllocDefault);
        cudaHostAlloc((void**)&m_hc, m_size * sizeof(DataType), cudaHostAllocDefault);
    } else {
        m_ha = new DataType[m_size];
        m_hb = new DataType[m_size];
        m_hc = new DataType[m_size];
    }

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
}


DummyStreamProducer::~DummyStreamProducer()
{
    //
    // free memory on the host side
    //
    if (m_isPinned) {
        cudaFreeHost(m_ha);
        cudaFreeHost(m_hb);
        cudaFreeHost(m_hc);
    } else {
        delete [] m_ha;
        delete [] m_hb;
        delete [] m_hc;
    }

    //
    // free just once at destruction
    //
    cudaFree(m_da);
    cudaFree(m_db);
    cudaFree(m_dc);

    // 
    // destroy Events and Streams
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
DummyStreamProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    cudaEventRecord(m_estart, m_stream);

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
    cudaEventRecord(m_estop, m_stream);
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
DummyStreamProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

void
DummyStreamProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DummyStreamProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyStreamProducer);
