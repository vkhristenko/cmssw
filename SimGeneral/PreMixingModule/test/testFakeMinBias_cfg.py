import FWCore.ParameterSet.Config as cms

process = cms.Process("MinBias")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(10)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('testFakeMinBias.root'),
)

process.e = cms.EndPath(process.out)
