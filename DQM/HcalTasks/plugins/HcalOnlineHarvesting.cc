#include "DQM/HcalTasks/interface/HcalOnlineHarvesting.h"

HcalOnlineHarvesting::HcalOnlineHarvesting(edm::ParameterSet const& ps) :
	DQHarvester(ps), _reportSummaryMap(NULL)
{
	_vsumgen.resize(nSummary);
	_vnames.resize(nSummary);
	_vmarks.resize(nSummary);
	_vsumgen[fRaw] = new RawRunSummary("RawRunSummary",
		"RawTask", ps);
	_vsumgen[fDigi] = new DigiRunSummary("DigiRunSummary", 
		"DigiTask",ps);
	_vsumgen[fReco] = new RecoRunSummary("RecoRunSummary",
		"RecHitTask", ps);
	_vsumgen[fTP] = new TPRunSummary("TPRunSummary",
		"TPTask", ps);
	_vnames[fRaw] = "RawTask";
	_vnames[fDigi] = "DigiTask";
	_vnames[fReco] = "RecHitTask";
	_vnames[fTP] = "TPTask";

	for (uint32_t i=0; i<_vsumgen.size(); i++)
		_vmarks.push_back(false);
}

/* virtual */ void HcalOnlineHarvesting::beginRun(
	edm::Run const& r, edm::EventSetup const& es)
{
	DQHarvester::beginRun(r,es);
	for (std::vector<DQClient*>::const_iterator it=_vsumgen.begin();
		it!=_vsumgen.end(); ++it)
		(*it)->beginRun(r,es);
}

/* virtual */ void HcalOnlineHarvesting::_dqmEndLuminosityBlock(
	DQMStore::IBooker& ib,
	DQMStore::IGetter& ig, edm::LuminosityBlock const&, 
	edm::EventSetup const&)
{
	//	DETERMINE WHICH MODULES ARE PRESENT IN DATA
	if (ig.get(_subsystem+"/RawTask/EventsTotal")!=NULL)
		_vmarks[fRaw]=true;
	if (ig.get(_subsystem+"/DigiTask/EventsTotal")!=NULL)
		_vmarks[fDigi]=true;
	if (ig.get(_subsystem+"/TPTask/EventsTotal")!=NULL)
		_vmarks[fTP]=true;
	if (ig.get(_subsystem+"/RecHitTask/EventsTotal")!=NULL)
		_vmarks[fReco]=true;

	//	CREATE SUMMARY REPORT MAP FED vs LS and LOAD MODULE'S SUMMARIES
	//	NOTE: THIS STATEMENTS WILL BE EXECUTED ONLY ONCE!
	if (!_reportSummaryMap)
	{
		ig.setCurrentFolder(_subsystem+"/EventInfo");
		_reportSummaryMap = ib.book2D("reportSummaryMap", "reportSummaryMap",
			_maxLS, 1, _maxLS+1, _vFEDs.size(), 0, _vFEDs.size());
		for (uint32_t i=0; i<_vFEDs.size(); i++)
		{
			char name[5];
			sprintf(name, "%d", _vFEDs[i]);
			_reportSummaryMap->setBinLabel(i+1, name, 2);
		}
		//	set LS bit to mark Xaxis as LS axis
		_reportSummaryMap->getTH1()->SetBit(BIT(BIT_OFFSET+BIT_AXIS_LS));

		// INITIALIZE ALL THE MODULES
		for (uint32_t i=0; i<_vnames.size(); i++)
			_vcSummaryvsLS.push_back(ContainerSingle2D(_vnames[i],
				"SummaryvsLS",
				new quantity::LumiSection(_maxLS),
				new quantity::FEDQuantity(_vFEDs),
				new quantity::ValueQuantity(quantity::fState)));

		//	LOAD ONLY THOSE MODULES THAT ARE PRESENT IN DATA
		for (uint32_t i=0; i<_vmarks.size(); i++)
		{
			if (_vmarks[i])
				_vcSummaryvsLS[i].load(ig, _subsystem);
		}
	}
/*	THIS DOES WORK
 *
 *	if (!_me)
	{
		ib.setCurrentFolder("Hcal/Folder1");
		_me = ib.book1D("MonitorElement", "MonitorElement",
			10, 0, 10);
		for (int i=0; i<10; i++)
			_me->Fill(i);
		std::cout << "11111111111" << std::endl;
	}
*/	

	int ifed=0;
	for (std::vector<uint32_t>::const_iterator it=_vhashFEDs.begin();
		it!=_vhashFEDs.end(); ++it)
	{
		HcalElectronicsId eid(*it);
		flag::Flag fSum("Status", flag::fNCDAQ);
		for (uint32_t im=0; im<_vmarks.size(); im++)
			if (_vmarks[im])
			{
				std::cout << "FED="<< _vFEDs[ifed] << std::endl;
				std::cout << _vnames[im] << std::endl;
				std::cout << eid << std::endl;
				int x = _vcSummaryvsLS[im].getBinContent(eid, _currentLS);
				std::cout << "x=" << x << std::endl;
				flag::Flag flag("Status", (flag::State)x);
				std::cout << flag._name << "  " <<  flag._state << std::endl;
				fSum+=flag;
			}
		_reportSummaryMap->setBinContent(_currentLS, ifed+1, int(fSum._state));
		ifed++;
	}
}

/* virtual */ void HcalOnlineHarvesting::_dqmEndJob(DQMStore::IBooker& ib,
	DQMStore::IGetter& ig)
{
	/*
	 *	THIS CODE DOESN"T WORK!
	 */
	ib.setCurrentFolder("Hcal/Folder1");
	_me = ib.book1D("MonitorElement", "MonitorElement",
		10, 0, 10);
	for (int i=0; i<10; i++)
		_me->Fill(i);
		
	std::cout << "11111111111" << std::endl;
	//	iterate over Run Summary Clients and generate Run Summary
//	for (int ii=fRaw; ii<nSummary; ii++)
//		if (_vmarks[ii])
//			_vsumgen[ii]->endJob(ib, ig);
}

DEFINE_FWK_MODULE(HcalOnlineHarvesting);
