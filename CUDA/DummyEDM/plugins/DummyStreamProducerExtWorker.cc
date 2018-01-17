// -*- C++ -*-
//
// Package:    TestGPU/Dummy
// Class:      DummyStreamProducerExtWorker
// 
/**\class DummyStreamProducerExtWorker DummyStreamProducerExtWorker.cc TestGPU/Dummy/plugins/DummyStreamProducerExtWorker.cc

 Description: A simple Dummy EDM Stream Producer with GPU offload

*/
//
// Original Author:  Viktor Khristenko
//         Created:  Tue, 14 Nov 2017 10:33:24 GMT
//
//


// system include files
#include <memory>
#include <chrono>
#include <thread>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Concurrency/interface/WaitingTaskWithArenaHolder.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// use the new product
#include "CUDA/DummyEDM/interface/Vector.h"
#include "CUDA/DummyEDM/interface/gpu_kernels.h"

#include <cuda.h>
#include <cuda_runtime.h>

//
// class declaration
//
class DummyStreamProducerExtWorker : public edm::stream::EDProducer<edm::ExternalWork> {
public:
    // some type aliasing
    using DataType = int;

    // ctor and dtor
    explicit DummyStreamProducerExtWorker(const edm::ParameterSet&);
    ~DummyStreamProducerExtWorker();

    void acquire(edm::Event const&, edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;
private:

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
DummyStreamProducerExtWorker::DummyStreamProducerExtWorker(const edm::ParameterSet& iConfig)
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
    produces<testgpu::Vector<DataType> >("VectorForGPUWithExtWorker");
}


DummyStreamProducerExtWorker::~DummyStreamProducerExtWorker()
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

// acquire
void DummyStreamProducerExtWorker::acquire(edm::Event const& iEvent,
                                           edm::EventSetup const& iSetup,
                                           edm::WaitingTaskWithArenaHolder holder) {
    std::thread(
        // capture list
          [holder, this]{
            // do the work in this thread
            std::cout << "starting the work in external thread" << std::endl;

            // initialize the values in the input vector
            for (auto i=0; i<m_size; i++) {
                m_ha[i] = i;
                m_hb[i] = i+i;
            }

            // record the start of the event
            cudaEventRecord(m_estart, m_stream);

            // perform the copy to device
            auto startCopy = std::chrono::high_resolution_clock::now();
            cudaMemcpyAsync(m_da, m_ha, m_size * sizeof(DataType),
                            cudaMemcpyHostToDevice, m_stream);
            auto finishCopy = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = finishCopy - startCopy;
            printf("Duration of cudaMemcpyAsync function call: %f s\n", diff.count());
            startCopy = std::chrono::high_resolution_clock::now();
            cudaMemcpyAsync(m_db, m_hb, m_size * sizeof(DataType),
                            cudaMemcpyHostToDevice, m_stream);
            finishCopy = std::chrono::high_resolution_clock::now();
            diff = finishCopy - startCopy;
            printf("Duration of cudaMemcpyAsync function call: %f s\n", diff.count());

            // call the kernel
            testgpu::wrapperVectorAdd<DataType>(m_da, m_db, m_dc, m_stream, m_size);

            // copy the results back
            startCopy = std::chrono::high_resolution_clock::now();
            cudaMemcpyAsync(m_hc, m_dc, m_size * sizeof(DataType),
                            cudaMemcpyDeviceToHost, m_stream);
            finishCopy = std::chrono::high_resolution_clock::now();
            diff = finishCopy - startCopy;
            printf("Duration of cudaMemcpyAsync function call: %f s\n", diff.count());

            // sync the stream
            cudaStreamSynchronize(m_stream);

            // record the stop
            cudaEventRecord(m_estop, m_stream);
            cudaEventSynchronize(m_estop);
            float elapsedTime {0.0};
            cudaEventElapsedTime(&elapsedTime, m_estart, m_estop);
            printf("Time taken: %3.1f ms\n", elapsedTime);

            // move the holder and call done
            edm::WaitingTaskWithArenaHolder newh = std::move(holder);
            std::exception_ptr exc;
            newh.doneWaiting(exc);
         }
    ).detach();
}

// ------------ method called to produce the data  ------------
void
DummyStreamProducerExtWorker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // all the work has been done in acquire
    // put the vector into the edm::Event
    testgpu::Vector<int> v;
    v.m_values = std::vector<int>(m_hc, m_hc + m_size);
    iEvent.put(std::make_unique<testgpu::Vector<int>>(v), "VectorForGPUWithExtWorker");
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyStreamProducerExtWorker);
