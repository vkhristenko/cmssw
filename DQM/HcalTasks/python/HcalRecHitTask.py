import FWCore.ParameterSet.Config as cms 

import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalRecHitTask"
#	Modify whatever is in StandardSet importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"Hcal/%s/" % moduleName)

folderSum = "HcalSum2DMaps"

#	Creating 2D Maps
HcalMap = [StandardSet.Standard2DMap.clone() for x in range(9)]
for i in range(3):
	HcalMap[i].path			= cms.untracked.string("Hcal/%s/%s" % (moduleName,
		folderSum)
	)
	HcalMap[i+3].path		= cms.untracked.string("Hcal/%s/%s" % (moduleName,
		folderSum)
	)
	HcalMap[i+6].path		= cms.untracked.string("Hcal/%s/" % moduleName)
	HcalMap[i].desc			= cms.untracked.string(
			"HB HE HF Depth%d Energy Sum (GeV)" % (i+1))
	HcalMap[i+3].desc		= cms.untracked.string(
			"HB HE HF Depth%d Time Sum (ns)" % (i+1))
	HcalMap[i+6].desc		= cms.untracked.string(
			"HB HE HF Depth%d Occupancy" % (i+1))

#	Extract Phi/Eta profiles and define axis for that
enaxis = cms.untracked.PSet(
	wnbins								= cms.untracked.bool(True),
	nbins								= cms.untracked.int32(200),
	min									= cms.untracked.double(0),
	max									= cms.untracked.double(400),
	title								= cms.untracked.string("Energy (GeV)")
)
StandardSet.StandardPhiProf.yaxis			= enaxis
StandardSet.StandardPhiProf.path			= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.StandardPhiProf.desc			= cms.untracked.string(
	"Hcal Energy vs iphi Profile")
StandardSet.StandardEtaProf.yaxis			= enaxis
StandardSet.StandardEtaProf.path			= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.StandardEtaProf.desc			= cms.untracked.string(
	"Hcal Energy vs ieta Profile")

#	Main Task Description
hcalRecHitTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,
		HBRecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB RecHit Energy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HERecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE RecHit Energy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HORecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO RecHit Energy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HFRecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF RecHit Energy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HBRecHitOcc					= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB RecHit Occupancy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HERecHitOcc					= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE RecHit Occupancy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HORecHitOcc					= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO RecHit Occupancy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HFRecHitOcc					= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF RecHit Occupancy"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HBRecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB RecHit Time"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HERecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE RecHit Time"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HORecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO RecHit Time"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HFRecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF RecHit Time"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(210),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HcalEnergyEta				= StandardSet.StandardEtaProf,
		HcalEnergyPhi				= StandardSet.StandardPhiProf,
		HBHEHFEnergySumMapD1		= HcalMap[0],
		HBHEHFEnergySumMapD2		= HcalMap[1],
		HBHEHFEnergySumMapD3		= HcalMap[2],
		HBHEHFTimeSumMapD1			= HcalMap[3],
		HBHEHFTimeSumMapD2			= HcalMap[4],
		HBHEHFTimeSumMapD3			= HcalMap[5],
		HBHEHFOccupancyMapD1		= HcalMap[6],
		HBHEHFOccupancyMapD2		= HcalMap[7],
		HBHEHFOccupancyMapD3		= HcalMap[8],
		RecHitSizeCheck			= StandardSet.Standard2DMap 
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
