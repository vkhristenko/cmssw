import FWCore.ParameterSet.Config as cms 

import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalRecHitTask"
#	Modify whatever is in StandardSet importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"%s/" % moduleName)

#	Defs for Reset	
resetEv = cms.untracked.string("EVENT")
resetLS = cms.untracked.string("LS")

#	Define some Useful strings
cutstr_no		= " No Cuts Applied. "
cutstr_eg2		= " DQM Cut (Energy > 2GeV) "
cutstr_eg0		= " DQM Cut (Energy > 0GeV) "
sumFolder		= "/HcalSum2DMaps"

#	Axis Defs
#	Energy, Time, Fraction
enProfaxis = cms.untracked.PSet(
	wnbins								= cms.untracked.bool(True),
	nbins								= cms.untracked.int32(100),
	min									= cms.untracked.double(-10),
	max									= cms.untracked.double(300),
	title								= cms.untracked.string("Energy (GeV)")
)
enTHaxis = cms.untracked.PSet(
		edges	= cms.untracked.bool(False),
		nbins	= cms.untracked.int32(100),
		min		= cms.untracked.double(-10.),
		max		= cms.untracked.double(300),
		title	= cms.untracked.string("Energy (GeV)")
)
timeProfaxis = cms.untracked.PSet(
	wnbins								= cms.untracked.bool(True),
	nbins								= cms.untracked.int32(200),
	min									= cms.untracked.double(-50),
	max									= cms.untracked.double(50),
	title								= cms.untracked.string("Time (ns)")
)
timeTHaxis = cms.untracked.PSet(
	edges								= cms.untracked.bool(False),
	nbins								= cms.untracked.int32(200),
	min									= cms.untracked.double(-50),
	max									= cms.untracked.double(50),
	title								= cms.untracked.string("Time (ns)")
)
fractionProfaxis = cms.untracked.PSet(
	wnbins								= cms.untracked.bool(False),
	min									= cms.untracked.double(0),
	max									= cms.untracked.double(1.05),
	title								= cms.untracked.string("Fraction ()")
)

#	Main Task Description
hcalRecHitTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,

		#	RecHit Energy Plots
		HB_RecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB RecHit Energy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(300),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(500),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HE_RecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE RecHit Energy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(300),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(500),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HO_RecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO RecHit Energy" + 
				cutstr_eg0),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(300),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(500),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HF_RecHitEnergy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF RecHit Energy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(300),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(500),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),

		HB_RecHitEnergyVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HB RecHit Energy vs ieta" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= enProfaxis
		),
		HE_RecHitEnergyVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HE RecHit Energy vs ieta" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= enProfaxis
		),
		HO_RecHitEnergyVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HO RecHit Energy vs ieta" + 
				cutstr_eg0),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= enProfaxis
		),
		HF_RecHitEnergyVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HF RecHit Energy vs ieta" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= enProfaxis
		),
		HBHEHFD1_EnergyMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HBHEHF D1 Average Energy Map"
					+ cutstr_eg2),
				xaxis		= StandardSet.ietaAxis.clone(),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(0),
					max			= cms.untracked.double(300)
				)
		),
		HBHEHFD2_EnergyMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HBHEHF D2 Average Energy Map"
					+ cutstr_eg2),
				xaxis		= StandardSet.ietaAxis.clone(),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(0),
					max			= cms.untracked.double(300)
				)
		),
		HBHEHFD3_EnergyMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HBHEHF D3 Average Energy Map"
					+ cutstr_eg2),
				xaxis		= StandardSet.ietaAxis.clone(),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(0),
					max			= cms.untracked.double(300)
				)
		),
		HOD4_EnergyMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HO D4 Average Energy Map"
					+ cutstr_eg0),
				xaxis		= cms.untracked.PSet(
					edges	= cms.untracked.bool(False),
					nbins	= cms.untracked.int32(33),
					min		= cms.untracked.double(-16.5),
					max		= cms.untracked.double(16.5),
					title	= cms.untracked.string("ieta")
				),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(0),
					max			= cms.untracked.double(100)
				)
		),
		
		#	Rec Hit Occupancies Plots
		HB_RecHitOccupancy					= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB RecHit Occupancy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HE_RecHitOccupancy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE RecHit Occupancy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HO_RecHitOccupancy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO RecHit Occupancy" + 
				cutstr_eg0),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HF_RecHitOccupancy				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF RecHit Occupancy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(500),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(2500),
				title	= cms.untracked.string("Occupancy")
			)
		),
		HB_RecHitOccupancyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HB RecHit Occupancy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(False),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5000),
				title	= cms.untracked.string("# RecHits")
			)
		),
		HE_RecHitOccupancyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HE RecHit Occupancy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(False),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5000),
				title	= cms.untracked.string("# RecHits")
			)
		),
		HO_RecHitOccupancyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HO RecHit Occupancy" + 
				cutstr_eg0),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(False),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5000),
				title	= cms.untracked.string("# RecHits")
			)
		),
		HF_RecHitOccupancyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HF RecHit Occupancy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(False),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5000),
				title	= cms.untracked.string("# RecHits")
			)
		),

		#	Plots took from FSQ meeting. Energy vs LS
		HB_RecHitEnergyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HB RecHit Energy vs LS" + 
				cutstr_no),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(250),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(250),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HE_RecHitEnergyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HE RecHit Energy vs LS" + 
				cutstr_no),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(250),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(250),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HO_RecHitEnergyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HO RecHit Energy vs LS" + 
				cutstr_no),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(250),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(250),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),
		HF_RecHitEnergyVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HF RecHit Energy vs LS" + 
				cutstr_no),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(250),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(250),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(200),
				title	= cms.untracked.string("Energy (GeV)")
			)
		),

		HBHEHFD1_RecHitOccupancy			= cms.untracked.PSet(
				path	= cms.untracked.string("%s" % moduleName),
				kind	= cms.untracked.string("TH2D"),
				desc	= cms.untracked.string("HBHEHF D1 RecHit Occupancy" + 
					cutstr_no),
				xaxis	= StandardSet.ietaAxis.clone(),
				yaxis	= StandardSet.iphiAxis.clone()
		),
		HBHEHFD2_RecHitOccupancy			= cms.untracked.PSet(
				path	= cms.untracked.string("%s" % moduleName),
				kind	= cms.untracked.string("TH2D"),
				desc	= cms.untracked.string("HBHEHF D2 RecHit Occupancy" + 
					cutstr_no),
				xaxis	= StandardSet.ietaAxis.clone(),
				yaxis	= StandardSet.iphiAxis.clone()
		),
		HBHEHFD3_RecHitOccupancy			= cms.untracked.PSet(
				path	= cms.untracked.string("%s" % moduleName),
				kind	= cms.untracked.string("TH2D"),
				desc	= cms.untracked.string("HBHEHF D3 RecHit Occupancy" + 
					cutstr_no),
				xaxis	= StandardSet.ietaAxis.clone(),
				yaxis	= StandardSet.iphiAxis.clone()
		),
		HOD4_RecHitOccupancy			= cms.untracked.PSet(
				path	= cms.untracked.string("%s" % moduleName),
				kind	= cms.untracked.string("TH2D"),
				desc	= cms.untracked.string("HO D4 RecHit Occupancy" + 
					cutstr_no),
				xaxis	= cms.untracked.PSet(
					edges	= cms.untracked.bool(False),
					nbins	= cms.untracked.int32(33),
					min		= cms.untracked.double(-16.5),
					max		= cms.untracked.double(16.5),
					title	= cms.untracked.string("ieta")
				),
				yaxis	= StandardSet.iphiAxis.clone()
		),

		#	Rec Hit Timing Plots
		HB_RecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB RecHit Time" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(-100.),
				max		= cms.untracked.double(100),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HE_RecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE RecHit Time" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(-100.),
				max		= cms.untracked.double(100),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HO_RecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO RecHit Time" + 
				cutstr_eg0),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(-100.),
				max		= cms.untracked.double(100),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HF_RecHitTime				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF RecHit Time" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(-100.),
				max		= cms.untracked.double(100),
				title	= cms.untracked.string("Time (ns)")
			)
		),
		HB_RecHitTimeVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HB RecHit Time vs ieta" + 
				cutstr_eg2),
			xaxis	= StandardSet.ietaAxis.clone(),
			yaxis	= enProfaxis
		),
		HE_RecHitTimeVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HE RecHit Time vs ieta" + 
				cutstr_eg2),
			xaxis	= StandardSet.ietaAxis.clone(),
			yaxis	= enProfaxis
		),
		HO_RecHitTimeVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HO RecHit Time vs ieta" + 
				cutstr_eg0),
			xaxis	= StandardSet.ietaAxis.clone(),
			yaxis	= enProfaxis
		),
		HF_RecHitTimeVSieta		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string("HF RecHit Time vs ieta" + 
				cutstr_eg2),
			xaxis	= StandardSet.ietaAxis.clone(),
			yaxis	= enProfaxis
		),
		HBHEHFD1_TimingMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HBHEHF D1 Average Timing Map"
					+ cutstr_eg2),
				xaxis		= StandardSet.ietaAxis.clone(),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(-100),
					max			= cms.untracked.double(100)
				)
		),
		HBHEHFD2_TimingMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HBHEHF D2 Average Timing Map"
					+ cutstr_eg2),
				xaxis		= StandardSet.ietaAxis.clone(),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(-100),
					max			= cms.untracked.double(100)
				)
		),
		HBHEHFD3_TimingMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HBHEHF D3 Average Timing Map"
					+ cutstr_eg2),
				xaxis		= StandardSet.ietaAxis.clone(),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(-100),
					max			= cms.untracked.double(100)
				)
		),
		HOD4_TimingMap			= cms.untracked.PSet(
				path		= cms.untracked.string("%s" % moduleName),
				kind		= cms.untracked.string("PROF2D"),
				desc		= cms.untracked.string("HO D4 Average Timing Map"
					+ cutstr_eg0),
				xaxis		= cms.untracked.PSet(
					edges	= cms.untracked.bool(False),
					nbins	= cms.untracked.int32(33),
					min		= cms.untracked.double(-16.5),
					max		= cms.untracked.double(16.4),
					title	= cms.untracked.string("ieta")
				),
				yaxis		= StandardSet.iphiAxis.clone(),
				zaxis	= cms.untracked.PSet(
					wnbins		= cms.untracked.bool(False),
					min			= cms.untracked.double(-100),
					max			= cms.untracked.double(100)
				)
		),

		#	Time vs Energy Plots
		HB_RecHitTimeVSenergy		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HB" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HB RecHit Time vs Energy" + 
				cutstr_eg2),
			xaxis	= enTHaxis,
			yaxis	= timeTHaxis
		),
		HE_RecHitTimeVSenergy		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HE" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HE RecHit Time vs Energy" + 
				cutstr_eg2),
			xaxis	= enTHaxis,
			yaxis	= timeTHaxis
		),
		HO_RecHitTimeVSenergy		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HO" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HO RecHit Time vs Energy" + 
				cutstr_eg0),
			xaxis	= enTHaxis,
			yaxis	= timeTHaxis
		),
		HF_RecHitTimeVSenergy		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string("HF RecHit Time vs Energy" + 
				cutstr_eg2),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(50),
				min		= cms.untracked.double(-10.),
				max		= cms.untracked.double(100),
				title	= cms.untracked.string("Energy (GeV)")
			),
			yaxis	= timeTHaxis
		)
	)
)











