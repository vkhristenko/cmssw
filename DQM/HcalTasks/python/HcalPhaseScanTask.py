import FWCore.ParameterSet.Config as cms 

#	import standard cfg and clone the parameters
import DQM.HcalCommon.HcalDQStandard as standard
StandardSet = standard.StandardSet.clone()

#	List of FEDs
lFEDs = [x+700 for x in range(32)] + [929, 1118, 1120, 1122]

moduleName = "HcalPhaseScanTask"
#	Modify whatever is in standard importing
StandardSet.moduleParameters.name		= cms.untracked.string(moduleName)
StandardSet.EventsProcessed.path		= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.EventsProcessedPerLS.path	= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.Standard2DMap.path			= cms.untracked.string(
	"%s/" % moduleName)
StandardSet.Standard2DMap.desc			= cms.untracked.string(
	"Some PhaseScan Task 2D Map")

#	Main Task Description
hcalPhaseScanTask = cms.EDAnalyzer(
	moduleName,
	moduleParameters	= StandardSet.moduleParameters,
	MEs					= cms.untracked.PSet(
		#	A Must!
		EventsProcessed			= StandardSet.EventsProcessed,
		EventsProcessedPerLS	= StandardSet.EventsProcessedPerLS,
		
		#	Source-specific
		HFP_Shape				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFP Signal Shape. Nominal fC are on the Y-axis"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string('TS')
			)
		),
		HFM_Shape				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFM Signal Shape. Nominal fC are on the Y-axis"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string('TS')
			)
		),
		HFP_Shape_3TSQg20		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFP Signal Shape (Using Sum NomfC of 3TS>%d). Nominal fC are on the Y-axis" % 20),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string('TS')
			)
		),
		HFM_Shape_3TSQg20		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFM Signal Shape (Using Sum NomfC of 3TS>%d). Nominal fC are on the Y-axis" % 20),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(5),
				min		= cms.untracked.double(-0.5),
				max		= cms.untracked.double(4.5),
				title	= cms.untracked.string('TS')
			)
		),
		HFP_Timing				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFP Timing (Nominal fC-weighted time average). After the Cut"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0.5),
				max		= cms.untracked.double(3.5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HFM_Timing				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFM Timing (Nominal fC-weighted time average). After the Cut"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HF_TimingVSieta			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HF Timing (Nominal fC-weighted time average). After the Cut vs ieta"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HF_TimingVSieta2D		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HF Timing (Nominal fC-weighted time average). After the Cut vs ieta"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0.5),
				max		= cms.untracked.double(3.5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HFP_TimingVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFP Timing (Nominal fC-weighted time average). After the Cut vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HFM_TimingVSls			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFM Timing (Nominal fC-weighted time average). After the Cut vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HFP_TimingVSls2D			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HFP Timing (Nominal fC-weighted time average). After the Cut vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0.5),
				max		= cms.untracked.double(3.5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		HFM_TimingVSls2D			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HFM Timing (Nominal fC-weighted time average). After the Cut vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0.5),
				max		= cms.untracked.double(3.5),
				title	= cms.untracked.string("TS(ns-like)")
			)
		),
		SumQ_3TS				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string("3TS nominal fC Sum"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(30),
				title	= cms.untracked.string("Lin. ADC")
			)
		),
		HFP_QTS2QTS12				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFP Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2))")
			)
		),
		HFM_QTS2QTS12				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFM Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2))")
			)
		),
		HFP_QTS2QTS23				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFP Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3))")
			)
		),
		HFM_QTS2QTS23				= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH1D"),
			desc	= cms.untracked.string(
				"HFM Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3))")
			)
		),
		HFP_QTS2QTS12vsLS			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFP Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2)")
			)
		),
		HFM_QTS2QTS12vsLS			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFM Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2)")
			)
		),
		HFP_QTS2QTS23vsLS			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFP Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3)")
			)
		),
		HFM_QTS2QTS23vsLS			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFM Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3)")
			)
		),
		HFMiphi43_QTS2QTS23vsLS			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFM iphi43 Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3)")
			)
		),
		HFMiphi43_QTS2QTS12vsLS			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HFM iphi43 Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2)")
			)
		),
		HFP_QTS2QTS12vsLS2D			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HFP Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2)")
			)
		),
		HFM_QTS2QTS12vsLS2D			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HFM Ratio Q(TS=2)/sum(Q(TS=1-2)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=1-2))")
			)
		),
		HFP_QTS2QTS23vsLS2D			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HFP Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3))")
			)
		),
		HFM_QTS2QTS23vsLS2D			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(
				"HFM Ratio Q(TS=2)/sum(Q(TS=2-3)), SumQ > 20 linADC vs LS"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(1000),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("Q(TS=2)/sum(Q(TS=2-3)")
			)
		),
		HFM_OccupancyietavsLS		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(	"HFM Occupancy [ieta vs LS]"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(13),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(-28.5),
				title	= cms.untracked.string("ieta")
			)
		),
		HFP_OccupancyietavsLS		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(	"HFP Occupancy [ieta vs LS]"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(100),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1000),
				title	= cms.untracked.string("LS")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(13),
				min		= cms.untracked.double(28.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			)
		),
		HF_OccupancyD1		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(	"HF Occupancy D1"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(72),
				min		= cms.untracked.double(0.5),
				max		= cms.untracked.double(72.5),
				title	= cms.untracked.string("iphi")
			)
		),
		HF_OccupancyD2		= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("TH2D"),
			desc	= cms.untracked.string(	"HF Occupancy D2"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(72),
				min		= cms.untracked.double(0.5),
				max		= cms.untracked.double(72.5),
				title	= cms.untracked.string("iphi")
			)
		),
		HF_OccupancyVSieta			= cms.untracked.PSet(
			path	= cms.untracked.string("%s/HF" % moduleName),
			kind	= cms.untracked.string("PROF"),
			desc	= cms.untracked.string(
				"HF Occupancy VS ieta"),
			xaxis	= cms.untracked.PSet(
				edges	= cms.untracked.bool(False),
				nbins	= cms.untracked.int32(83),
				min		= cms.untracked.double(-41.5),
				max		= cms.untracked.double(41.5),
				title	= cms.untracked.string("ieta")
			),
			yaxis	= cms.untracked.PSet(
				wnbins	= cms.untracked.bool(True),
				nbins	= cms.untracked.int32(200),
				min		= cms.untracked.double(0),
				max		= cms.untracked.double(1.05),
				title	= cms.untracked.string("# Hits")
			)
		)
	)
)










