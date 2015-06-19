import FWCore.ParameterSet.Config as cms 

import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcaluTCATask"
#	Modify whatever is in StandardSet importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.Standard2DMap.path			= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.Standard2DMap.desc			= cms.untracked.string(
	"Some uTCA Task 2D Map")

hfVMEFEDs	= [718, 719]
hfuTCAFEDs	= [1118, 1119, 1120]

#	Main Task Description
hcaluTCATask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,

		#	Plots to compare uTCA vs VME
		HF_ADC			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("ADC Comparison"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(130),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(129.5),
				title	= cms.untracked.string("VME FEDs (%s)" % hfVMEFEDs)
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(130),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(129.5),
				title	= cms.untracked.string("uTCA FEDs (%s)" % hfuTCAFEDs)
			)
		)
	)
)
