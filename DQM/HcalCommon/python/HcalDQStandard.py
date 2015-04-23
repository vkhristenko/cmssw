import FWCore.ParameterSet.Config as cms

#	Generate a list of FEDs
lFEDs = [x+700 for x in range(32)] + [1118, 1120, 1122]

StandardSet		= cms.untracked.PSet(
	moduleParameters	= cms.untracked.PSet(
			name			= cms.untracked.string("HcalDQStandard"),
		debug			= cms.untracked.int32(10),
		calibTypes		= cms.untracked.vint32(0),
		runType			= cms.untracked.string("TEST"),
		mtype			= cms.untracked.string("SOURCE"),
		FEDs			= cms.untracked.vint32(lFEDs)
	),
	
	Labels				= cms.untracked.PSet(
		HBHEDigi		= cms.untracked.InputTag("hcalDigis"),
		HFDigi			= cms.untracked.InputTag("hcalDigis"),
		HODigi			= cms.untracked.InputTag("hcalDigis"),
		RAW				= cms.untracked.InputTag("rawDataCollector"),
		HBHERecHit		= cms.untracked.InputTag("hbhereco"),
		HFRecHit		= cms.untracked.InputTag("hfreco"),
		HORecHit		= cms.untracked.InputTag("horeco"),
		L1GT			= cms.untracked.InputTag("l1GtUnpack"),
		HLTResults		= cms.untracked.InputTag("TriggerResults"),
		DCS				= cms.untracked.InputTag("scalersRawToDigi"),
		UnpackerReport	= cms.untracked.InputTag("hcalDigis")
	),

	EventsProcessed		= cms.untracked.PSet(
		path			= cms.untracked.string("Hcal/HcalDQStandard/"),
		kind			= cms.untracked.string("INT"),
	#	desc			= cms.untracked.string("Processed Events Total"),
	),

	EventsProcessedPerLS = cms.untracked.PSet(
		path			= cms.untracked.string("Hcal/HcalDQStandard/"),
		kind			= cms.untracked.string("INT"),
	#	desc			= cms.untracked.
	),

	Standard2DMap		= cms.untracked.PSet(
		path			= cms.untracked.string("Hcal/HcalDQStandard/"),
		kind			= cms.untracked.string("TH2D"),
		desc			= cms.untracked.string("Standard 2D Map"),
		xaxis			= cms.untracked.PSet(
			edges		= cms.untracked.bool(False),
			nbins		= cms.untracked.int32(84),
			min			= cms.untracked.double(-42),
			max			= cms.untracked.double(42),
			title		= cms.untracked.string("ieta")
		),
		yaxis			= cms.untracked.PSet(
			edges		= cms.untracked.bool(False),
			nbins		= cms.untracked.int32(72),
			min			= cms.untracked.double(0),
			max			= cms.untracked.double(72),
			title		= cms.untracked.string("iphi")
		)
	),

	Standard2DProf		= cms.untracked.PSet(
		path			= cms.untracked.string("Hcal/HcalDQStandard/"),
		kind			= cms.untracked.string("PROF2D"),
		desc			= cms.untracked.string("Standard 2D Profile"),
		xaxis			= cms.untracked.PSet(
			edges		= cms.untracked.bool(False),
			nbins		= cms.untracked.int32(84),
			min			= cms.untracked.double(-42),
			max			= cms.untracked.double(42),
			title		= cms.untracked.string("ieta")
		),
		yaxis			= cms.untracked.PSet(
			edges		= cms.untracked.bool(False),
			nbins		= cms.untracked.int32(72),
			min			= cms.untracked.double(0),
			max			= cms.untracked.double(72),
			title		= cms.untracked.string("iphi")
		)
	),

	StandardPhiProf		= cms.untracked.PSet(
		path			= cms.untracked.string("Hcal/HcalDQStandard/"),
		kind			= cms.untracked.string("PROF"),
		desc			= cms.untracked.string("Standard Phi Profile"),
		xaxis			= cms.untracked.PSet(
			edges		= cms.untracked.bool(False),
			nbins		= cms.untracked.int32(72),
			min			= cms.untracked.double(0),
			max			= cms.untracked.double(72),
			title		= cms.untracked.string("iphi")
		)
	),

	StandardEtaProf		= cms.untracked.PSet(
		path			= cms.untracked.string("Hcal/HcalDQStandard/"),
		kind			= cms.untracked.string("PROF"),
		desc			= cms.untracked.string("Standard Eta Profile"),
		xaxis			= cms.untracked.PSet(
			edges		= cms.untracked.bool(False),
			nbins		= cms.untracked.int32(84),
			min			= cms.untracked.double(-42),
			max			= cms.untracked.double(42),
			title		= cms.untracked.string("ieta")
		)
	)
)

StandardSet.moduleParameters.Labels = StandardSet.Labels

