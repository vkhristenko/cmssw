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

process.DependencyGraph = cms.Service(
    "DependencyGraph",
    fileName = cms.untracked.string("test_graph.svg"),
    showPathDependencies = cms.untracked.bool(True)
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

process.analyzer1 = cms.EDAnalyzer("TestAnalyzer1")
process.analyzer2 = cms.EDAnalyzer("TestAnalyzer2")
process.analyzer3 = cms.EDAnalyzer("TestAnalyzer3")

# output cfg
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_modules_out.root')
)

# path config
process.p1 = cms.Path(process.prod1 + process.filter1 + process.analyzer1)
process.p2 = cms.Path(process.prod2 + process.filter2 + process.analyzer3)
process.p3 = cms.Path(process.prod3 + process.filter3 + process.analyzer3)
process.p4 = cms.Path(process.prod1 + process.prod2 + process.prod3 + process.analyzer3)
process.e = cms.EndPath(process.out)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(8),
    numberOfStreams = cms.untracked.uint32(8)
)
