
#include "DQM/HcalTasks/interface/LEDTask.h"

using namespace hcaldqm;
LEDTask::LEDTask(edm::ParameterSet const& ps):
	DQTask(ps)
{
	//	Containers
	_cSignalMean_Subdet.initialize(_name, "SignalMean",
		hashfunctions::fSubdet, 
		new quantity::ValueQuantity(quantity::ffC_3000),
		new quantity::ValueQuantity(quantity::fN, true));
	_cSignalRMS_Subdet.initialize(_name, "SignalRMS",
		hashfunctions::fSubdet, 
		new quantity::ValueQuantity(quantity::ffC_1000),
		new quantity::ValueQuantity(quantity::fN, true));
	_cTimingMean_Subdet.initialize(_name, "TimingMean",
		hashfunctions::fSubdet, 
		new quantity::ValueQuantity(quantity::fTiming_TS200),
		new quantity::ValueQuantity(quantity::fN, true));
	_cTimingRMS_Subdet.initialize(_name, "TimingRMS",
		hashfunctions::fSubdet, 
		new quantity::ValueQuantity(quantity::fTiming_TS200), 
		new quantity::ValueQuantity(quantity::fN, true));

	_cShapeCut_FEDSlot.initialize(_name, "Shape", 
		hashfunctions::fFEDSlot,
		new quantity::ValueQuantity(quantity::fTiming_TS),
		new quantity::ValueQuantity(quantity::ffC_3000));

	_cSignalMean_depth.initialize(_name, "SignalMean",
		hashfunctions::fdepth, 
		new quantity::DetectorQuantity(quantity::fieta), 
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::ffC_3000));
	_cSignalRMS_depth.initialize(_name, "SignalRMS",
		hashfunctions::fdepth, 
		new quantity::DetectorQuantity(quantity::fieta), 
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::ffC_1000));
	_cTimingMean_depth.initialize(_name, "TimingMean",
		hashfunctions::fdepth, 
		new quantity::DetectorQuantity(quantity::fieta), 
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fTiming_TS200));
	_cTimingRMS_depth.initialize(_name, "TimingRMS",
		hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta), 
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fTiming_TS200));

	_cTimingVME.initialize(_name, "TimingMean",
		new quantity::ElectronicsQuantity(quantity::fFEDVME),
		new quantity::ElectronicsQuantity(quantity::fSpigot),
		new quantity::ValueQuantity(quantity::fTiming_TS200));
	_cSignalVME.initialize(_name, "SignalMean",
		new quantity::ElectronicsQuantity(quantity::fFEDVME),
		new quantity::ElectronicsQuantity(quantity::fSpigot),
		new quantity::ValueQuantity(quantity::ffC_3000));
	_cTiminguTCA.initialize(_name, "TimingMean",
		new quantity::ElectronicsQuantity(quantity::fFEDuTCA),
		new quantity::ElectronicsQuantity(quantity::fSlotuTCA),
		new quantity::ValueQuantity(quantity::fTiming_TS200));
	_cSignaluTCA.initialize(_name, "SignalMean",
		new quantity::ElectronicsQuantity(quantity::fFEDuTCA),
		new quantity::ElectronicsQuantity(quantity::fSlotuTCA),
		new quantity::ValueQuantity(quantity::ffC_3000));

	_cOccupancyVME.initialize(_name, "Occupancy",
		new quantity::ElectronicsQuantity(quantity::fFEDVME),
		new quantity::ElectronicsQuantity(quantity::fSpigot),
		new quantity::ValueQuantity(quantity::fN));
	_cOccupancyuTCA.initialize(_name, "Occupancy",
		new quantity::ElectronicsQuantity(quantity::fFEDuTCA),
		new quantity::ElectronicsQuantity(quantity::fSlotuTCA),
		new quantity::ValueQuantity(quantity::fN));

	_cMissing_depth.initialize(_name, "Missing",
		hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta),
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fN));

	
	//	initialize compact containers
	_cSignals_DChannel.initialize(hashfunctions::fDChannel);
	_cTiming_DChannel.initialize(hashfunctions::fDChannel);

	//	tags
	_tagHBHE = ps.getUntrackedParameter<edm::InputTag>("tagHBHE",
		edm::InputTag("hcalDigis"));
	_tagHO = ps.getUntrackedParameter<edm::InputTag>("tagHO",
		edm::InputTag("hcalDigis"));
	_tagHF = ps.getUntrackedParameter<edm::InputTag>("tagHF",
		edm::InputTag("hcalDigis"));
	_tagTrigger = ps.getUntrackedParameter<edm::InputTag>("tagTrigger",
		edm::InputTag("tbunpacker"));
	_tokHBHE = consumes<HBHEDigiCollection>(_tagHBHE);
	_tokHO = consumes<HODigiCollection>(_tagHO);
	_tokHF = consumes<HFDigiCollection>(_tagHF);
	_tokTrigger = consumes<HcalTBTriggerData>(_tagTrigger);

	//	constants
	_lowHBHE = ps.getUntrackedParameter<double>("lowHBHE",
		20);
	_lowHO = ps.getUntrackedParameter<double>("lowHO",
		20);
	_lowHF = ps.getUntrackedParameter<double>("lowHF",
		20);
}

/* virtual */ void LEDTask::bookHistograms(DQMStore::IBooker &ib,
	edm::Run const& r, edm::EventSetup const& es)
{
	if (_ptype==fLocal)
		if (r.runAuxiliary().run()==1)
			return;

	char cutstr[20];
	sprintf(cutstr, "sumQHBHE%dHO%dHF%d", int(_lowHBHE),
		int(_lowHO), int(_lowHF));

	DQTask::bookHistograms(ib, r, es);

	edm::ESHandle<HcalDbService> dbService;
	es.get<HcalDbRecord>().get(dbService);
	_emap = dbService->getHcalMapping();

	_cSignalMean_Subdet.book(ib, _emap);
	_cSignalRMS_Subdet.book(ib, _emap);
	_cTimingMean_Subdet.book(ib, _emap);
	_cTimingRMS_Subdet.book(ib, _emap);

	_cSignalMean_depth.book(ib, _emap);
	_cSignalRMS_depth.book(ib, _emap);
	_cTimingMean_depth.book(ib, _emap);
	_cTimingRMS_depth.book(ib, _emap);

	_cShapeCut_FEDSlot.book(ib, _emap);

	_cTimingVME.book(ib, _subsystem, std::string("VME"));
	_cSignalVME.book(ib, _subsystem, std::string("VME"));
	_cTiminguTCA.book(ib, _subsystem, std::string("uTCA"));
	_cSignaluTCA.book(ib, _subsystem, std::string("uTCA"));
	_cOccupancyVME.book(ib, _subsystem, std::string("VME"));
	_cOccupancyuTCA.book(ib, _subsystem, std::string("uTCA"));

	_cMissing_depth.book(ib, _emap);

	//	book compact containers
	_cSignals_DChannel.book(_emap);
	_cTiming_DChannel.book(_emap);
}

/* virtual */ void LEDTask::_resetMonitors(UpdateFreq uf)
{
	DQTask::_resetMonitors(uf);
}

/* virtual */ void LEDTask::_dump()
{
	_cSignalMean_Subdet.reset();
	_cSignalRMS_Subdet.reset();
	_cTimingMean_Subdet.reset();
	_cTimingRMS_Subdet.reset();
	_cSignalMean_depth.reset();
	_cSignalRMS_depth.reset();
	_cTimingMean_depth.reset();
	_cTimingRMS_depth.reset();

	_cSignals_DChannel.dump(&_cSignalMean_Subdet, &_cMissing_depth, true);
	_cSignals_DChannel.dump(&_cSignalRMS_Subdet, &_cMissing_depth,false);
	_cTiming_DChannel.dump(&_cTimingMean_Subdet, &_cMissing_depth,true);
	_cTiming_DChannel.dump(&_cTimingRMS_Subdet, &_cMissing_depth,false);
	_cSignals_DChannel.dump(&_cSignalMean_depth, &_cMissing_depth,true);
	_cSignals_DChannel.dump(&_cSignalRMS_depth, &_cMissing_depth,false);
	_cTiming_DChannel.dump(&_cTimingMean_depth,&_cMissing_depth, true);
	_cTiming_DChannel.dump(&_cTimingRMS_depth,&_cMissing_depth, false);
}

/* virtual */ void LEDTask::_process(edm::Event const& e,
	edm::EventSetup const& es)
{
	edm::Handle<HBHEDigiCollection>		chbhe;
	edm::Handle<HODigiCollection>		cho;
	edm::Handle<HFDigiCollection>		chf;

	if (!e.getByToken(_tokHBHE, chbhe))
		_logger.dqmthrow("Collection HBHEDigiCollection isn't available "
			+ _tagHBHE.label() + " " + _tagHBHE.instance());
	if (!e.getByToken(_tokHO, cho))
		_logger.dqmthrow("Collection HODigiCollection isn't available "
			+ _tagHO.label() + " " + _tagHO.instance());
	if (!e.getByToken(_tokHF, chf))
		_logger.dqmthrow("Collection HFDigiCollection isn't available "
			+ _tagHF.label() + " " + _tagHF.instance());

//	int currentEvent = e.eventAuxiliary().id().event();

	for (HBHEDigiCollection::const_iterator it=chbhe->begin();
		it!=chbhe->end(); ++it)
	{
		const HBHEDataFrame digi = (const HBHEDataFrame)(*it);
		double sumQ = utilities::sumQ<HBHEDataFrame>(digi, 2.5, 0, 
			digi.size()-1);
		if (sumQ<_lowHBHE)
			continue;
		HcalDetId did = digi.id();
		HcalElectronicsId eid = digi.elecId();

		double aveTS = utilities::aveTS<HBHEDataFrame>(digi, 2.5, 0,
			digi.size()-1);
		_cSignals_DChannel.fill(did, sumQ>0 ? sumQ : GARBAGE_VALUE);
		_cTiming_DChannel.fill(did, sumQ>0 ? aveTS : GARBAGE_VALUE);

		if (eid.isVMEid())
		{
			_cTimingVME.fill(eid, aveTS);
			_cSignalVME.fill(eid, sumQ);
			_cOccupancyVME.fill(eid);
		}
		else
		{
			_cTiminguTCA.fill(eid, aveTS);
			_cSignaluTCA.fill(eid, sumQ);
			_cOccupancyuTCA.fill(eid);
		}

		for (int i=0; i<digi.size(); i++)
			_cShapeCut_FEDSlot.fill(eid, i, 
				digi.sample(i).nominal_fC()-2.5);
	}
	for (HODigiCollection::const_iterator it=cho->begin();
		it!=cho->end(); ++it)
	{
		const HODataFrame digi = (const HODataFrame)(*it);
		double sumQ = utilities::sumQ<HODataFrame>(digi, 8.5, 0, 
			digi.size()-1);
		if (sumQ<_lowHO)
			continue;
		HcalDetId did = digi.id();
		HcalElectronicsId eid = digi.elecId();

		double aveTS = utilities::aveTS<HODataFrame>(digi, 8.5, 0,
			digi.size()-1);
		_cSignals_DChannel.fill(did, sumQ>0 ? sumQ : GARBAGE_VALUE);
		_cTiming_DChannel.fill(did, sumQ>0 ? aveTS : GARBAGE_VALUE);

		if (eid.isVMEid())
		{
			_cTimingVME.fill(eid, aveTS);
			_cSignalVME.fill(eid, sumQ);
			_cOccupancyVME.fill(eid);
		}
		else
		{
			_cTiminguTCA.fill(eid, aveTS);
			_cSignaluTCA.fill(eid, sumQ);
			_cOccupancyuTCA.fill(eid);
		}

		for (int i=0; i<digi.size(); i++)
			_cShapeCut_FEDSlot.fill(eid, i, 
				digi.sample(i).nominal_fC()-8.5);
	}
	for (HFDigiCollection::const_iterator it=chf->begin();
		it!=chf->end(); ++it)
	{
		const HFDataFrame digi = (const HFDataFrame)(*it);
		double sumQ = utilities::sumQ<HFDataFrame>(digi, 2.5, 0, 
			digi.size()-1);
		if (sumQ<_lowHF)
			continue;
		HcalDetId did = digi.id();
		HcalElectronicsId eid = digi.elecId();

		double aveTS = utilities::aveTS<HFDataFrame>(digi, 2.5, 0,
			digi.size()-1);
		_cSignals_DChannel.fill(did, sumQ>0 ? sumQ : GARBAGE_VALUE);
		_cTiming_DChannel.fill(did, sumQ>0 ? aveTS : GARBAGE_VALUE);

		if (eid.isVMEid())
		{
			_cTimingVME.fill(eid, aveTS);
			_cSignalVME.fill(eid, sumQ);
			_cOccupancyVME.fill(eid);
		}
		else
		{
			_cTiminguTCA.fill(eid, aveTS);
			_cSignaluTCA.fill(eid, sumQ);
			_cOccupancyuTCA.fill(eid);
		}

		for (int i=0; i<digi.size(); i++)
			_cShapeCut_FEDSlot.fill(eid, i, 
				digi.sample(i).nominal_fC()-2.5);
	}

	if (_ptype==fOnline && _evsTotal>0 &&
		_evsTotal%constants::CALIBEVENTS_MIN==0)
		this->_dump();
}

/* virtual */ bool LEDTask::_isApplicable(edm::Event const& e)
{
	if (_ptype!=fOnline)
	{
		//	local
		edm::Handle<HcalTBTriggerData> ctrigger;
		if (!e.getByToken(_tokTrigger, ctrigger))
			_logger.dqmthrow("Collection HcalTBTriggerData isn't available "
				+ _tagTrigger.label() + " " + _tagTrigger.instance());
		return ctrigger->wasLEDTrigger();
	}

	return false;
}

DEFINE_FWK_MODULE(LEDTask);


