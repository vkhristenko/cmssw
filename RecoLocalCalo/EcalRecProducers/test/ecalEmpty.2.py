import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
#from Configuration.ProcessModifiers.gpu_cff import gpu

process = cms.Process('RECO', eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeterogeneousCore.CUDAServices.CUDAService_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_HLT_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_dataRun3_HLT_v1', '')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# load data using the DAQ source
import sys, os, inspect
sys.path.append(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
process.load('sourceFromRawCmggpu_cff')

#-----------------------------------------
# CMSSW/Hcal non-DQM Related Module import
#-----------------------------------------
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")

process.load("RecoLocalCalo.EcalRecProducers.ecalChannelStatusGPUESProducer_cfi")
    
    
    
#    
# addition
#


process.load("RecoLocalCalo.EcalRecProducers.ecalPedestalsGPUESProducer_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalGainRatiosGPUESProducer_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalPulseShapesGPUESProducer_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalPulseCovariancesGPUESProducer_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalSamplesCorrelationGPUESProducer_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalTimeBiasCorrectionsGPUESProducer_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalTimeCalibConstantsGPUESProducer_cfi")

#
# so far it works
#


process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")

process.load("EventFilter.EcalRawToDigi.ecalRawToDigiGPU_cfi")
process.load("EventFilter.EcalRawToDigi.ecalElectronicsMappingGPUESProducer_cfi")

    
   
#
# so far it works
#



process.ecalDigis = process.ecalEBunpacker.clone()
process.ecalDigis.InputLabel = cms.InputTag('rawDataCollector')



process.digiPath = cms.Path(
    process.ecalDigis
    *process.ecalRawToDigiGPU
)
 
    
    
    
    
    
    
    
    
    
    
    
    

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("test.root")
)

process.finalize = cms.EndPath(process.out)

process.bunchSpacing = cms.Path(
    process.bunchSpacingProducer
)


process.schedule = cms.Schedule(
    process.bunchSpacing,
    process.digiPath,    # new
    process.finalize
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(4),
    numberOfStreams = cms.untracked.uint32(4),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# report CUDAService messages
process.MessageLogger.categories.append("CUDAService")
