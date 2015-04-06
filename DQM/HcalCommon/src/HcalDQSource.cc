

#include "DQM/HcalCommon/interface/HcalDQSource.h"

namespace hcaldqm
{
	HcalDQSource::HcalDQSource(edm::ParameterSet const& ps) :
		HcalDQMonitor(ps.getUntrackedParameterSet("moduleParameters")), 
		_mes(ps.getUntrackedParameterSet("MEs"))
	{}

	HcalDQSource::~HcalDQSource() {}

	//	Functions to be reimplemented from DQMEDAnalyzer
	/* virtual */ void HcalDQSource::analyze(edm::Event const &e,
			edm::EventSetup const& es)
	{
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
}









