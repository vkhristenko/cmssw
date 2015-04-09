#include "DQM/HcalCommon/interface/HcalDQMonitor.h"

namespace hcaldqm
{
	//	Constructor
	HcalDQMonitor::HcalDQMonitor(edm::ParameterSet const& ps) : 
		_labels(ps.getUntrackedParameterSet("Labels"))
	{
		_mi.type		= ps.getUntrackedParameter<std::string>("mtype");
		_mi.runType		= ps.getUntrackedParameter<std::string>("runType");
		_mi.calibTypesAllowed	= ps.getUntrackedParameter<std::vector<int> >(
				"calibTypes");
		_mi.feds		= ps.getUntrackedParameter<std::vector<int> >(
				"FEDs");
		_mi.name		= ps.getUntrackedParameter<std::string>("name");
		_mi.debug		= ps.getUntrackedParameter<int>("debug");
	}

	//	Destructor
	/* virtual */ HcalDQMonitor::~HcalDQMonitor()
	{}
}
