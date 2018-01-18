import FWCore.ParameterSet.Config as cms

#
# Define the CMSSW Process
#
process = cms.Process("DumpRawPixelData")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

#
# Source 
#
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        "file:step2_DIGIPREMIX_S2_DATAMIX_L1_DIGI2RAW_HLT.root"
#        "root://eoscms.cern.ch//eos/cms/store/data/Run2017B/ZeroBias/RAW/v1/000/298/991/00000/70D06660-4768-E711-A6FE-02163E019DE0.root"
        "file:0049088C-8DF6-E711-A790-0CC47A7C340C.root"
#        "root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_0_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/100X_upgrade2018_realistic_v4_HS-v1/20000/0049088C-8DF6-E711-A790-0CC47A7C340C.root"
    )
)

#
# Here is our GPU Producer
#
process.dump = cms.EDProducer(
    'DumpRawPixelDataCPU',
    InputLabel = cms.InputTag("rawDataCollector")
)

#
# OutputModule
#
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.dump)

process.e = cms.EndPath(process.out)

#
# Set the threads/streams options
#
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1)
)
