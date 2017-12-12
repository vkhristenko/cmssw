import FWCore.ParameterSet.Config as cms
import os

#
# define a new process
#
process = cms.Process("TestGPU")

#
# Message Logger 
#
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring(
        'cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG')
    ),
    debugModules = cms.untracked.vstring("*")
)

#
# 10 events max to process
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

#
# Nothing to read -> empty source
#
process.source = cms.Source("EmptySource")

#
# Declare an edm::stream producer to put into the pipeline
#
process.testGPU = cms.EDProducer('DummyStreamProducer',
    size = cms.untracked.int32(10000),
    # allocates the data as Page-Locked
    isPinned = cms.untracked.bool(True)
)

#
# the pipeline
#
process.p = cms.Path(process.testGPU)

#
# Preserve the edm Event content into the ROOT file
#
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("test_streamproducer.root")
)

#
# ROOT file output should be run at the very end
#
process.finalize = cms.EndPath(process.out)

#
# Set the number of threads and CMSSW streams
#
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(4),
    numberOfStreams = cms.untracked.uint32(4)
)
