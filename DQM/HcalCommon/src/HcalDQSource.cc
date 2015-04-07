

#include "DQM/HcalCommon/interface/HcalDQSource.h"

namespace hcaldqm
{
	HcalDQSource::HcalDQSource(edm::ParameterSet const& ps) :
		HcalDQMonitor(ps.getUntrackedParameterSet("moduleParameters")), 
		_mes(ps.getUntrackedParameterSet("MEs"))
	{
		_si.currentCalibType = -1;
	}

	/* virtual */HcalDQSource::~HcalDQSource() 
	{
		this->warn_("Calling Destructor...");
	}

	//	Functions to be reimplemented from DQMEDAnalyzer
	/* virtual */ void HcalDQSource::analyze(edm::Event const &e,
			edm::EventSetup const& es)
	{
		this->extractCalibType(e);
		if (this->isAllowedCalibType()==false)
			return;
		this->doWork(e, es);
	}

	/* virtual */ void HcalDQSource::bookHistograms(DQMStore::IBooker &ib,
			edm::Run const& r, edm::EventSetup const& es)
	{
		_mes.book(ib);
	}

	/* virtual */ void HcalDQSource::dqmBeginRun(edm::Run const& r,
			edm::EventSetup const& es)
	{}

	//	extract Event Calibration Type from FEDs
	void HcalDQSource::extractCalibType(edm::Event const&e)
	{
		edm::Handle<FEDRawDataCollection> craw;
		INITCOLL(_labels["RAW"], craw);

		//	for now
		_si.currentCalibType = 0;
	}

	//	Check if calibration type set is among allowed
	bool HcalDQSource::isAllowedCalibType()
	{
		for (std::vector<int>::const_iterator it=_mi.calibTypesAllowed.begin();
				it!=_mi.calibTypesAllowed.end(); ++it)
			if (_si.currentCalibType==*it)
				return true;

		return false;
	}
}









