import FWCore.ParameterSet.Config as cms

process = cms.Process("TestModules")

# msg logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring(
        "cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("DEBUG")),
    debugModules = cms.untracked.vstring("*")
)

# source cfg
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# modules cfg
process.prod1 = cms.EDProducer('TestProducer1')
process.prod2 = cms.EDProducer('TestProducer2')
process.prod3 = cms.EDProducer('TestProducer3')

process.filter1 = cms.EDFilter("TestFilter1")
process.filter2 = cms.EDFilter("TestFilter2")
process.filter3 = cms.EDFilter("TestFilter3")

process.anaalyzer1 = cms.EDAnalyzer("TestAnalyzer1")
process.anaalyzer2 = cms.EDAnalyzer("TestAnalyzer2")
process.anaalyzer3 = cms.EDAnalyzer("TestAnalyzer3")

# output cfg
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

# path config
process.p = cms.Path(process.myProducerLabel)
process.e = cms.EndPath(process.out)
