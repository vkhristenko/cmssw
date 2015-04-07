import FWCore.ParameterSet.Config as cms 
from DQM.HcalCommon.HcalDQSourceEx import hcalDQSourceEx

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalBeamTask"
hcalBeamTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= cms.untracked.PSet(
		name		= cms.untracked.string(moduleName),
		debug		= cms.untracked.int32(10),
		calibTypes	= cms.untracked.vint32(0),
		runType		= cms.untracked.string("TEST"),
		mtype		= cms.untracked.string("SOURCE"),
		Labels		= hcalDQSourceEx.moduleParameters.Labels
	),
	MEs					= cms.untracked.PSet(
		HEBeamShape				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
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
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
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
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
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
		BeamSizeCheck			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("Beam Size Check for SubSystems"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5),
				title	= cms.untracked.string("Subsystem")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(20),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(20),
				title	= cms.untracked.string("Beam Size")
			)
		),
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
