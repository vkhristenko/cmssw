import FWCore.ParameterSet.Config as cms 

#	import standard cfg and clone the parameters
import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalDigiTask"
#	Modify whatever is in standard importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"Hcal/%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"Hcal/%s/" % moduleName)



HcalMap = [StandardSet.Standard2DMap.clone() for x in range(3)]
for i in range(3):
	HcalMap[i].path						= cms.untracked.string("Hcal/%s/" %
			moduleName)
	HcalMap[i].desc						= cms.untracked.string(
	"HB HE HF Depth%d Occupancy" % (i+1))



#	Main Task Description
hcalDigiTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,	
		HBNumberDigis			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB/Diagnostics" 
				% moduleName),
			kind	= cms.untracked.string("INT")
		),
		HENumberDigis			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE/Diagnostics" 
				% moduleName),
			kind	= cms.untracked.string("INT")
		),
		HONumberDigis			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO/Diagnostics" 
				% moduleName),
			kind	= cms.untracked.string("INT")
		),
		HFNumberDigis			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF/Diagnostics" 
				% moduleName),
			kind	= cms.untracked.string("INT")
		),
		HBDigiShape			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB Digi Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HEDigiShape			= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE Digi Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HFDigiShape				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF Digi Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)
		),
		HODigiShape				= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO Digi Shape"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(10),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(10.),
				title	= cms.untracked.string("TS")
			)	
		),
		HBADCCountPerTS		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB ADC Counts per 1TS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(130),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(130),
				title	= cms.untracked.string("TS")
			)	
		),
		HEADCCountPerTS		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE ADC Counts per 1TS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(130),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(130),
				title	= cms.untracked.string("TS")
			)	
		),
		HOADCCountPerTS		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO ADC Counts per 1TS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(130),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(130),
				title	= cms.untracked.string("TS")
			)	
		),
		HFADCCountPerTS		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF ADC Counts per 1TS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(130),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(130),
				title	= cms.untracked.string("TS")
			)	
		),
		HBPresamples		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB Number of Presamples"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(20),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(20),
				title	= cms.untracked.string("# Presamples")
			)	
		),
		HEPresamples		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE Number of Presamples"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(20),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(20),
				title	= cms.untracked.string("# Presamples")
			)	
		),
		HOPresamples		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO Number of Presamples"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(20),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(20),
				title	= cms.untracked.string("# Presamples")
			)	
		),
		HFPresamples		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF Number of Presamples"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(20),
				min		= cms.untracked.double(0.),
				max		= cms.untracked.double(20),
				title	= cms.untracked.string("# Presamples")
			)	
		),
		HBCapId		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB Cap ID"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string("Cap ID")
			)	
		),
		HECapId		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE Cap ID"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string("Cap ID")
			)	
		),
		HOCapId		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO Cap ID"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string("Cap ID")
			)	
		),
		HFCapId		= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF Cap ID"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string("Cap ID")
			)	
		),
		HBbcnOffset	= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HB" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HB Fiber Idle BCN Offset"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(16),
				min		= cms.untracked.double(-8),
				max		= cms.untracked.double(8),
				title	= cms.untracked.string("BCN Offset")
			)	
		),
		HEbcnOffset	= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HE" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HE Fiber Idle BCN Offset"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(16),
				min		= cms.untracked.double(-8),
				max		= cms.untracked.double(8),
				title	= cms.untracked.string("BCN Offset")
			)	
		),
		HObcnOffset	= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HO" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HO Fiber Idle BCN Offset"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(16),
				min		= cms.untracked.double(-8),
				max		= cms.untracked.double(8),
				title	= cms.untracked.string("BCN Offset")
			)	
		),
		HFbcnOffset	= cms.untracked.PSet(
			path	= cms.untracked.string("Hcal/%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("HF Fiber Idle BCN Offset"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(16),
				min		= cms.untracked.double(-8),
				max		= cms.untracked.double(8),
				title	= cms.untracked.string("BCN Offset")
			)	
		),
		HBHEHFOccupancyMapD1	= HcalMap[0],
		HBHEHFOccupancyMapD2	= HcalMap[1],
		HBHEHFOccupancyMapD3	= HcalMap[2],
		DigiSizeCheck			= StandardSet.Standard2DMap 
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
