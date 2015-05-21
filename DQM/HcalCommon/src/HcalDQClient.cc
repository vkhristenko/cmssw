#include "DQM/HcalCommon/interface/HcalDQClient.h"

namespace hcaldqm
{
	HcalDQClient::HcalDQClient(edm::ParameterSet const& ps)
		: HcalDQMonitor(ps.getUntrackedParameterSet("moduleParameters")),
		_bmes(ps.getUntrackedParameterSet("bookMEs"), _mi.debug),
		_rmes(ps.getUntrackedParameterSet("retrieveMEs"), _mi.debug)
	{}

	/* virtual */HcalDQClient::~HcalDQClient()
	{
		this->debug_("Calling Destructor");
	}

	//	Function to be reimplemented from DQMEDAnalyzer
	//	Executed at the end of the job
	/* virtual */ void HcalDQClient::dqmEndJob(DQMStore::IBooker& ib,
			DQMStore::IGetter& ig)
	{
		_bmes.book(ib);	
		doWork(ib, ig);
	}

	//	Function to be reimplemented from the DQMEDAnalyzer
	//	Executed at the edn of LS
	/* virtual */ void HcalDQClient::dqmEndLuminosityBlock(DQMStore::IGetter& ig,
			edm::LuminosityBlock const& ls, edm::EventSetup const& es)
	{
		_rmes.retrieve(ig);
		doWork(ig, ls, es);
	}

	//	reset
	/* virtual */ void HcalDQClient::reset(int const periodflag)
	{
		//	Collection Class determines itself who needs a reset and when
		//	Do it only for Monitor Modules which have been booked in this client
		_bmes.reset(periodflag);

		if (periodflag==0)
		{
			//	each event 
		}
		else if (periodflag==1)
		{
			//	each LS
			_mi.evsPerLS = 0;
		}
	}
}












