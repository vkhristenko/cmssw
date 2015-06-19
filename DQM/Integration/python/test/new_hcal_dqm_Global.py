#-------------------------------------
#	Hcal DQM Application using New DQM Sources/Clients
#-------------------------------------

#-------------------------------------
#	Standard Python Imports
#-------------------------------------
import os, sys, socket, string

#-------------------------------------
#	Standard CMSSW Imports/Definitions
#-------------------------------------
import FWCore.ParameterSet.Config as cms
process			= cms.Process('HCALDQM')
subsystem		= 'Hcal'
cmssw			= os.getenv("CMSSW_VERSION").split("_")
debugstr		= "### HcalDQM::cfg::DEBUG: "
warnstr			= "### HcalDQM::cfg::WARN: "
errorstr		= "### HcalDQM::cfg::ERROR:"
useOfflineGT	= False
useFileInput	= False
useMap		= False

#-------------------------------------
#	Central DQM Stuff imports
#-------------------------------------
from DQM.Integration.test.online_customizations_cfi import *
if useOfflineGT:
	process.load('DQM.Integration.test.FrontierCondition_GT_Offline_cfi')
	process.GlobalTag.globaltag = 'GR_P_V56::All'
else:
	process.load('DQM.Integration.test.FrontierCondition_GT_cfi')
if useFileInput:
	process.load("DQM.Integration.test.fileinputsource_cfi")
else:
	process.load('DQM.Integration.test.inputsource_cfi')
process.load('DQMServices.Components.DQMEnvironment_cfi')
process.load('DQM.Integration.test.environment_cfi')

#-------------------------------------
#	Central DQM Customization
#-------------------------------------
process.dqmEnv.subSystemFolder = subsystem
referenceFileName = '/dqmdata/dqm/reference/hcal_reference.root'
process.DQMStore.referenceFileName = referenceFileName
process = customise(process)

#	Note, runType is obtained after importing DQM-related modules
#	=> DQM-dependent
runType			= process.runType.getRunType()
print debugstr, "Running with run type= ", runType

#-------------------------------------
#	CMSSW/Hcal non-DQM Related Module import
#-------------------------------------
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("L1Trigger.Configuration.L1DummyConfig_cff")
process.load("EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

#-------------------------------------
#	CMSSW/Hcal non-DQM Related Module Settings
#	-> runType
#	-> Generic Input tag for the Raw Collection
#	-> cmssw version
#	-> Turn off default blocking of dead channels from rechit collection
#	-> Drop Channel Status Bits (had benn 'HcalCellOff', "HcalCellDead")
#	-> For Trigger Primitives Emulation
#	-> L1 GT setting
#	-> Rename the hbheprereco to hbhereco
#-------------------------------------
runType			= process.runType.getRunType()
cmssw			= os.getenv("CMSSW_VERSION").split("_")
rawTag			= cms.InputTag("rawDataCollector")
process.essourceSev = cms.ESSource(
		"EmptyESSource",
		recordName		= cms.string("HcalSeverityLevelComputerRcd"),
		firstValid		= cms.vuint32(1),
		iovIsRunNotTime	= cms.bool(True)
)
process.hcalRecAlgos.DropChannelStatusBits = cms.vstring('')
process.emulTPDigis = \
		process.simHcalTriggerPrimitiveDigis.clone()
process.emulTPDigis.inputLabel = \
		cms.VInputTag("hcalDigis", 'hcalDigis')
process.emulTPDigis.FrontEndFormatError = \
		cms.bool(True)
process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)
process.emulTPDigis.FG_threshold = cms.uint32(2)
process.emulTPDigis.InputTagFEDRaw = rawTag
process.l1GtUnpack.DaqGtInputTag = rawTag
process.hbhereco = process.hbheprereco.clone()

#-------------------------------------
#	Hcal DQM Tasks and Clients import
#-------------------------------------
process.load("DQM.HcalTasks.HcalDigiTask")
process.load("DQM.HcalTasks.HcalDeadCellTask")
process.load("DQM.HcalTasks.HcalHotCellTask")
process.load("DQM.HcalTasks.HcalLEDTask")
process.load("DQM.HcalTasks.HcalLaserTask")
process.load("DQM.HcalTasks.HcalNoiseTask")
process.load("DQM.HcalTasks.HcalPedestalTask")
process.load("DQM.HcalTasks.HcalRawTask")
process.load("DQM.HcalTasks.HcalRecHitTask")
process.load("DQM.HcalTasks.HcalTPTask")
process.load("DQM.HcalTasks.HcalTimingTask")
process.load("DQM.HcalTasks.HcaluTCATask")
process.load("DQM.HcalTasks.HcalPhaseScanTask")
#-------------------------------------
#	To force using uTCA
#	Will not be here for Online DQM
#-------------------------------------
if useMap:
	process.es_pool = cms.ESSource("PoolDBESSource",
			process.CondDBSetup,
			timetype = cms.string('runnumber'),
			toGet = cms.VPSet(
				cms.PSet(
					record = cms.string(
						"HcalElectronicsMapRcd"
					),
					tag = cms.string(
						"HcalElectronicsMap_v7.05_hlt"
					)
				)
			),
			connect = cms.string(
				'frontier://FrontierProd/CMS_CONDITIONS'),
			authenticationMethod = cms.untracked.uint32(0)
	)	
	process.es_prefer_es_pool = cms.ESPrefer('PoolDBESSource', 'es_pool')

#-------------------------------------
#	For Debugginb
#-------------------------------------
#process.hcalTPTask.moduleParameters.debug = 0

#-------------------------------------
#	Some Settings before Finishing up
#-------------------------------------
process.hcalDigiTask.moduleParameters.subsystem = cms.untracked.string(subsystem)
process.hcalDeadCellTask.moduleParameters.subsystem = cms.untracked.string(
		subsystem)
process.hcalHotCellTask.moduleParameters.subsystem = cms.untracked.string(
		subsystem)
process.hcalLEDTask.moduleParameters.subsystem = cms.untracked.string(subsystem)
process.hcalLaserTask.moduleParameters.subsystem = cms.untracked.string(subsystem)
process.hcalNoiseTask.moduleParameters.subsystem = cms.untracked.string(subsystem)
process.hcalPedestalTask.moduleParameters.subsystem = cms.untracked.string(
		subsystem)
process.hcalRawTask.moduleParameters.subsystem = cms.untracked.string(subsystem)
process.hcalRecHitTask.moduleParameters.subsystem = cms.untracked.string(
		subsystem)
process.hcalTPTask.moduleParameters.subsystem = cms.untracked.string(subsystem)
process.hcalTimingTask.moduleParameters.subsystem = cms.untracked.string(
		subsystem)
process.hcalPhaseScanTask.moduleParameters.subsystem = cms.untracked.string(
		subsystem)

#-------------------------------------
#	Hcal DQM Tasks Sequence Definition
#-------------------------------------
process.tasksSequence = cms.Sequence(
		process.hcalDigiTask
		*process.hcalDeadCellTask
		*process.hcalHotCellTask
		*process.hcalNoiseTask
		*process.hcalRawTask
		*process.hcalRecHitTask
		*process.hcalTPTask
		*process.hcalTimingTask
#		*process.hcaluTCATask
		*process.hcalPhaseScanTask
)

#-------------------------------------
#	Execution Sequence Definition
#-------------------------------------
process.p = cms.Path(
					process.hcalDigis
                    *process.emulTPDigis
                    *process.l1GtUnpack
                    *process.horeco
                    *process.hfreco
                    *process.hbhereco
					*process.tasksSequence
                    #*process.qTester
                    *process.dqmEnv
                    *process.dqmSaver)

#-------------------------------------
#	Scheduling
#-------------------------------------
process.options = cms.untracked.PSet(
		Rethrow = cms.untracked.vstring(
			"ProductNotFound",
			"TooManyProducts",
			"TooFewProducts"
		)
)
