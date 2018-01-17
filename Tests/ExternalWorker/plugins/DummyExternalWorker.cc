// -*- C++ -*-
//
// Package:    TestGPU/Dummy
// Class:      DummyExternalWorker
// 
/**\class DummyExternalWorker DummyExternalWorker.cc TestGPU/Dummy/plugins/DummyExternalWorker.cc

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
#include "FWCore/Concurrency/interface/WaitingTaskWithArenaHolder.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//
class DummyExternalWorker : public edm::stream::EDProducer<edm::ExternalWork> {
public:
    // some type aliasing
    using DataType = int;

    // ctor and dtor
    explicit DummyExternalWorker(const edm::ParameterSet&);
    ~DummyExternalWorker();

    void acquire(edm::Event const&, edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    // ----------member data ---------------------------
};

//
// constructors and destructor
//
DummyExternalWorker::DummyExternalWorker(const edm::ParameterSet& iConfig)
{
}


DummyExternalWorker::~DummyExternalWorker()
{
}

// acquire
void DummyExternalWorker::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup,
                                  edm::WaitingTaskWithArenaHolder)

// ------------ method called to produce the data  ------------
void DummyExternalWorker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyExternalWorker);
