import FWCore.ParameterSet.Config as cms 

#	import standard cfg and clone the parameters
import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalBeamTask"
#	Modify whatever is in standard importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.Standard2DMap.path			= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.Standard2DMap.desc			= cms.untracked.string(
	"Some Beam Task 2D Map")

#	Main Task Description
hcalBeamTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,
		
		HEBeamShape				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE Beam Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HFBeamShape				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF Beam Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HOBeamShape				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO Beam Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)	
		),
		BeamSizeCheck			= StandardSet.Standard2DMap 
	)
)
