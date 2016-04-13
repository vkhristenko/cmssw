#include "DQM/HcalTasks/interface/HcalOfflineHarvesting.h"

HcalOfflineHarvesting::HcalOfflineHarvesting(edm::ParameterSet const& ps) :
	DQHarvester(ps), _reportSummaryMap(NULL)
{
	_vsumgen.push_back(new RawRunSummary("RawRunSummary",
		"RawTask", ps));
	_vsumgen.push_back(new DigiRunSummary("DigiRunSummary",
		"DigiTask", ps));
	_vsumgen.push_back(new RecoRunSummary("RecoRunSummary",
		"RecHitTask", ps));
	_vsumgen.push_back(new TPRunSummary("TPRunSummary",
		"TPTask", ps));

	for (uint32_t i=0; i<_vsumgen.size(); i++)
		_vmarks.push_back(false);
}

/* virtual */ void HcalOfflineHarvesting::beginRun(
	edm::Run const& r, edm::EventSetup const& es)
{
	DQHarvester::beginRun(r,es);

	for (std::vector<DQClient*>::const_iterator it=_vsumgen.begin();
		it!=_vsumgen.end(); ++it)
		(*it)->beginRun(r,es);
}

//
//	For OFFLINE there is no per LS evaluation
//
/* virtual */ void HcalOfflineHarvesting::_dqmEndLuminosityBlock(
	DQMStore::IBooker& ib,
	DQMStore::IGetter& ig, edm::LuminosityBlock const& lb, 
	edm::EventSetup const& es)
{
	int ii=0;
	for (std::vector<DQClient*>::const_iterator it=_vsumgen.begin();
		it!=_vsumgen.end(); ++it)
	{	
		//	run only if have to
		if (_vmarks[ii])
			(*it)->endLuminosityBlock(ib,ig,lb,es);
		ii++;
	}
}

//
//	Evaluate and Generate Run Summary
//
/* virtual */ void HcalOfflineHarvesting::_dqmEndJob(DQMStore::IBooker& ib,
	DQMStore::IGetter& ig)
{
	//	OBTAIN/SET WHICH MODULES ARE PRESENT
	int num=0; std::map<std::string, uint32_t> datatiers;
	if (ig.get(_subsystem+"/RawTask/EventsTotal")!=NULL)
	{
		_vmarks[fRaw]=true;num++;
		datatiers.insert(std::pair<std::string, uint32_t>("RAW",num));
	}
	if (ig.get(_subsystem+"/DigiTask/EventsTotal")!=NULL)
	{
		_vmarks[fDigi]=true;num++;
		datatiers.insert(std::pair<std::string, uint32_t>("DIGI",num));
	}
	if (ig.get(_subsystem+"/TPTask/EventsTotal")!=NULL)
	{
		_vmarks[fTP]=true;num++;
		datatiers.insert(std::pair<std::string, uint32_t>("TP",num));
	}
	if (ig.get(_subsystem+"/RecHitTask/EventsTotal")!=NULL)
	{
		_vmarks[fReco]=true;num++;
		datatiers.insert(std::pair<std::string, uint32_t>("RECO",num));
	}

	//	CREATE THE REPORT SUMMARY MAP
	if (!_reportSummaryMap)
	{
		ib.setCurrentFolder(_subsystem+"/EventInfo");
		_reportSummaryMap = ib.book2D("reportSummaryMap", "reportSummaryMap",
			_vFEDs.size(), 0, _vFEDs.size(), num,0,num);
		//	x axis labels
		for (uint32_t i=0; _vFEDs.size(); i++)
		{
			char name[5];
			sprintf(name, "%d", _vFEDs[i]);
			_reportSummaryMap->setBinLabel(i+1, name, 1);
		}
		//	y axis lables
		for (std::map<std::string, uint32_t>::const_iterator
			it=datatiers.begin(); it!=datatiers.end(); ++it)
		{
			_reportSummaryMap->setBinLabel(it->second, it->first, 2);
		}
	}

	//	iterate over all summary generators and get the flags
	int ii=0;
	for (std::vector<DQClient*>::const_iterator it=_vsumgen.begin();
		it!=_vsumgen.end(); ++it)
	{
		if (!_vmarks[ii])
		{ii++;continue;}

		//	do only if have to
		std::vector<flag::Flag> flags = (*it)->endJob(ib,ig);
		for (uint32_t ifed=0; ifed<_vFEDs.size(); ifed++)
		{
			_reportSummaryMap->setBinContent(ifed+1, 
				datatiers[flags[ifed]._name], (int)flags[ifed]._state);
		}
		ii++;
	}
}

DEFINE_FWK_MODULE(HcalOfflineHarvesting);
