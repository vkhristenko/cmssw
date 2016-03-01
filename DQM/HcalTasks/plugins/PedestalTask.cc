
#include "DQM/HcalTasks/interface/PedestalTask.h"

using namespace hcaldqm;
PedestalTask::PedestalTask(edm::ParameterSet const& ps):
	DQTask(ps)
{
	//	Containers
	_cMean_Subdet.initialize(_name, "Mean",hashfunctions::fSubdet, 
		new quantity::ValueQuantity(quantity::fADC_15),
		new quantity::ValueQuantity(quantity::fN, true));
	_cRMS_Subdet.initialize(_name, "RMS", hashfunctions::fSubdet, 
		new quantity::ValueQuantity(quantity::fADC_5),
		new quantity::ValueQuantity(quantity::fN, true));
	_cMean_depth.initialize(_name, "Mean", hashfunctions::fdepth, 
		new quantity::DetectorQuantity(quantity::fieta), 
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fADC_15));
	_cRMS_depth.initialize(_name, "RMS", hashfunctions::fdepth, 
		new quantity::DetectorQuantity(quantity::fieta), 
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fADC_5));
	_cPed.initialize(hashfunctions::fDChannel);
	_cPedRef.initialize(hashfunctions::fDChannel, compact::fDoubleValued);
	_cMeanDBRef_Subdet.initialize(_name, "MeanDBRef", hashfunctions::fSubdet,
		new quantity::ValueQuantity(quantity::fAroundZero),
		new quantity::ValueQuantity(quantity::fN, true));
	_cRMSDBRef_Subdet.initialize(_name, "RMSDBRef", hashfunctions::fSubdet,
		new quantity::ValueQuantity(quantity::fAroundZero),
		new quantity::ValueQuantity(quantity::fN, true));
	_cMeanDBRef_depth.initialize(_name, "MeanDBRef", hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta),
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fAroundZero));
	_cRMSDBRef_depth.initialize(_name, "RMSDBRef", hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta),
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fAroundZero));
	_cMissing_depth.initialize(_name, "Missing", hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta),
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fN));
	_cMeanBad_depth.initialize(_name, "MeanBad", hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta),
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fN));
	_cRMSBad_depth.initialize(_name, "RMSBad", hashfunctions::fdepth,
		new quantity::DetectorQuantity(quantity::fieta),
		new quantity::DetectorQuantity(quantity::fiphi),
		new quantity::ValueQuantity(quantity::fN));


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
}

/* virtual */ void PedestalTask::bookHistograms(DQMStore::IBooker &ib,
	edm::Run const& r, edm::EventSetup const& es)
{
	if (_ptype==fLocal)
		if (r.runAuxiliary().run()==1)
			return;
	DQTask::bookHistograms(ib, r, es);

	edm::ESHandle<HcalDbService> dbService;
	es.get<HcalDbRecord>().get(dbService);
	_emap = dbService->getHcalMapping();

	//	book plots
	_cMean_Subdet.book(ib, _emap);
	_cRMS_Subdet.book(ib, _emap);
	_cMean_depth.book(ib, _emap);
	_cRMS_depth.book(ib, _emap);
	_cMeanDBRef_Subdet.book(ib, _emap);
	_cRMSDBRef_Subdet.book(ib, _emap);
	_cMeanDBRef_depth.book(ib, _emap);
	_cRMSDBRef_depth.book(ib, _emap);
	_cMissing_depth.book(ib, _emap);
	_cMeanBad_depth.book(ib, _emap);
	_cRMSBad_depth.book(ib, _emap);
	
	//	book compact containers and load conditions pedestals
	_cPed.book(_emap);
	_cPedRef.book(_emap);
	_cPedRef.loadPedestals(dbService);
}

/* virtual */ void PedestalTask::_resetMonitors(UpdateFreq uf)
{
	DQTask::_resetMonitors(uf);
}

/* virtual */ void PedestalTask::endRun(edm::Run const& r, 
	edm::EventSetup const&)
{	
	if (_ptype==fLocal)
		if (r.runAuxiliary().run()==1)
			return;

	this->_dump();
/*
	/////////
//	DQMStore *store = edm::Service<DQMStore>().operator->();
//	store->showDirStructure();
//	std::cout <<"CURRENT DIRECTORY: "<< store->pwd() << std::endl;
//	store->open("DQM_V0001_Hcal_R000262665.root",
//		false, "", "HCAL", DQMStore::KeepRunDirs);
//	store->load("DQM_V0001_Hcal_R000262665.root");
//	store->setCurrentFolder("./Run 262665");
//	store->showDirStructure();
//	store->cd("Run ");
//	std::cout <<"11111111 CURRENT DIRECTORY: "<< store->pwd() << std::endl;
//	store->cd("../");
//	std::cout <<"2222222 CURRENT DIRECTORY: "<< store->pwd() << std::endl;
//	store->goUp();
//	std::cout <<"3333333 CURRENT DIRECTORY: "<< store->pwd() << std::endl;

//	sprintf(cutstr, "_sumQHBHE%dHO%dHF%d", int(20),
//		int(20), int(20));
//	Container2D	cMeans;
//	Container2D cRMSs;
//	cMeans.initialize(_name, "Mean", hashfunctions::fdepth, 
//		new quantity::DetectorQuantity(quantity::fieta), 
//		new quantity::DetectorQuantity(quantity::fiphi),
//		new quantity::ValueQuantity(quantity::fADC_15), 10);
//	cRMSs.initialize(_name, "RMS", hashfunctions::fdepth, 
//		new quantity::DetectorQuantity(quantity::fieta), 
//		new quantity::DetectorQuantity(quantity::fiphi),
//		new quantity::ValueQuantity(quantity::fADC_5));
//	cMeans.load(store, _emap, _subsystem, std::string(""),
//		std::string("HCAL/Run 262665"), 
//		DQMStore::KeepRunDirs);
//	cRMSs.load(store, _emap, _subsystem, std::string(""),
//		std::string("HCAL/Run 262665"),
//		DQMStore::KeepRunDirs);
//	cMeans.print();
//	cRMSs.print();
//
//
	std::cout << "SIZE=" << _cPedRef.size() << std::endl;
	_cPedRef.print();
//	_cPedRef.load(&cMeans, &cRMSs);
	_cPedRef.print();
	hcaldqm::QualityWrapper qwrap(&hcaldqm::pedestal_quality2, 0.1);
	_cPedRef.compare(_cPeds, &_cMeanRef_Subdet, &_cAbsent_depth,
		&_cBad_depth, &hcaldqm::diff, qwrap);
	_cPedRef.compare(_cPeds, &_cMeanRef_depth, &_cAbsent_depth,
		&_cBad_depth, &hcaldqm::diff, qwrap);
	_cPedRef.compare(_cPeds, &_cRMSRef_Subdet, &_cAbsent_depth,
		&_cBad_depth, &hcaldqm::diff, qwrap, false);
	_cPedRef.compare(_cPeds, &_cRMSRef_depth, &_cAbsent_depth,
		&_cBad_depth, &hcaldqm::diff, qwrap, false);
//	_cPeds.compare(_cPedRef, &_cMeanRef_Subdet);
//	_cPeds.compare(_cPedRef, &_cMeanRef_depth);
	_cPedRef.print();
	store->cd();
	std::cout <<"44444444 CURRENT DIRECTORY: "<< store->pwd() << std::endl;
	store->cd("Hcal");
	std::cout <<"55555555 CURRENT DIRECTORY: "<< store->pwd() << std::endl;
	store->cd("");
	store->cd("HCAL/Run 262665");
	std::cout <<"666666666 CURRENT DIRECTORY: "<< store->pwd() << std::endl;
	store->rmdir(store->pwd());
	
	store->showDirStructure();
	*/
}

/* virtual */ void PedestalTask::_dump()
{
	//	dump the pedestals first
	_cMean_Subdet.reset();
	_cRMS_Subdet.reset();
	_cMean_depth.reset();
	_cRMS_depth.reset();
	_cPed.dump(&_cMean_Subdet, true);
	_cPed.dump(&_cRMS_Subdet, false);
	_cPed.dump(&_cMean_depth, true);
	_cPed.dump(&_cRMS_depth, false);

	//	compare with Conditions Pedestals and dump
	_cMeanDBRef_Subdet.reset();
	_cMeanDBRef_depth.reset();
	_cRMSDBRef_Subdet.reset();
	_cRMSDBRef_depth.reset();
	_cMissing_depth.reset();
	_cMeanBad_depth.reset();
	_cRMSBad_depth.reset();
	hcaldqm::QualityWrapper qwrap(&hcaldqm::pedestal_quality2, 0.2);
	std::vector<Container1D*> vMean;
	std::vector<Container1D*> vRMS;
	vMean.push_back(&_cMeanDBRef_Subdet);
	vMean.push_back(&_cMeanDBRef_depth);
	vRMS.push_back(&_cRMSDBRef_Subdet);
	vRMS.push_back(&_cRMSDBRef_depth);
	_cPedRef.compare(_cPed, vMean, &_cMissing_depth,
		&_cMeanBad_depth, &hcaldqm::diff, qwrap); //	means
	_cPedRef.compare(_cPed, vRMS, &_cMissing_depth,
		&_cRMSBad_depth, &hcaldqm::diff, qwrap, false); //	rmss
}

/* virtual */ void PedestalTask::_process(edm::Event const& e,
	edm::EventSetup const& es)
{
	edm::Handle<HBHEDigiCollection>		chbhe;
	edm::Handle<HODigiCollection>		cho;
	edm::Handle<HFDigiCollection>		chf;

	if (!e.getByToken(_tokHBHE, chbhe))
		_logger.dqmthrow("Collection HBHEDigiCollection isn't available"
			+ _tagHBHE.label() + " " + _tagHBHE.instance());
	if (!e.getByToken(_tokHO, cho))
		_logger.dqmthrow("Collection HODigiCollection isn't available"
			+ _tagHO.label() + " " + _tagHO.instance());
	if (!e.getByToken(_tokHF, chf))
		_logger.dqmthrow("Collection HFDigiCollection isn't available"
			+ _tagHF.label() + " " + _tagHF.instance());

	for (HBHEDigiCollection::const_iterator it=chbhe->begin();
		it!=chbhe->end(); ++it)
	{
		const HBHEDataFrame digi = (const HBHEDataFrame)(*it);
		HcalDetId did = digi.id();
		int digiSizeToUse = floor(digi.size()/constants::CAPS_NUM)*
			constants::CAPS_NUM-1;
		for (int i=0; i<digiSizeToUse; i++)
		{
			_cPed.fill(did, it->sample(i).adc());
		}
	}
	for (HODigiCollection::const_iterator it=cho->begin();
		it!=cho->end(); ++it)
	{
		const HODataFrame digi = (const HODataFrame)(*it);
		HcalDetId did = digi.id();
		int digiSizeToUse = floor(digi.size()/constants::CAPS_NUM)*
			constants::CAPS_NUM-1;
		for (int i=0; i<digiSizeToUse; i++)
		{
			_cPed.fill(did, it->sample(i).adc());
		}
	}
	for (HFDigiCollection::const_iterator it=chf->begin();
		it!=chf->end(); ++it)
	{
		const HFDataFrame digi = (const HFDataFrame)(*it);
		HcalDetId did = digi.id();
		int digiSizeToUse = floor(digi.size()/constants::CAPS_NUM)*
			constants::CAPS_NUM-1;
		for (int i=0; i<digiSizeToUse; i++)
		{
			_cPed.fill(did, it->sample(i).adc());
		}
	}

	if (_ptype==fOnline && _evsTotal>0 && 
		_evsTotal%constants::CALIBEVENTS_MIN==0)
		this->_dump();
}

/* virtual */ bool PedestalTask::_isApplicable(edm::Event const& e)
{
	if (_ptype==fOnline)
	{
		//	online-global
		return this->_getCalibType(e)==hc_Pedestal;
	}
	else 
	{
		//	local
		edm::Handle<HcalTBTriggerData>	ctrigger;
		if (!e.getByToken(_tokTrigger, ctrigger))
			_logger.dqmthrow("Collection HcalTBTriggerData isn't available"
				+ _tagTrigger.label() + " " + _tagTrigger.instance());
		return ctrigger->wasSpillIgnorantPedestalTrigger();
	}

	return false;
}

DEFINE_FWK_MODULE(PedestalTask);


