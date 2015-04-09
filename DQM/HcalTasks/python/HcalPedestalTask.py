import FWCore.ParameterSet.Config as cms 

import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalPedestalTask"
#	Modify whatever is in StandardSet importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.Standard2DMap.path			= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.Standard2DMap.desc			= cms.untracked.string(
	"Some Pedestal Task 2D Map")

#	Main Task Description
hcalPedestalTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,
		
		HEPedestalShape				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE Pedestal Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HFPedestalShape				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF Pedestal Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HOPedestalShape				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO Pedestal Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)	
		),
		PedestalSizeCheck			= StandardSet.Standard2DMap 
#		me4			= cms.untracked.PSet(
#			path	= cms.untracked.string("Hcal/%s/" % moduleName),
#			kind	= cms.untracked.string("PROF"),
#			desc	= cms.untracked.string("Example ME4"),
#			xaxis	= cms.untracked.PSet(
#				edges	= cms.untracked.bool(False),
#				nbins	= cms.untracked.int32(200),
#				min		= cms.untracked.double(-100),
#				max		= cms.untracked.double(100),
#				title	= cms.untracked.string("me4-X")
#			),
#			yaxis	= cms.untracked.PSet(
#				wnbins	= cms.untracked.bool(True),
#				nbins	= cms.untracked.int32(100),
#				min		= cms.untracked.double(-50),
#				max		= cms.untracked.double(50),
#				title	= cms.untracked.string("me4-Y")
#			)
#		)
	)
)
