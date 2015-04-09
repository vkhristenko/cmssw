//	cmssw includes
#include "DQM/HcalTasks/interface/HcalLaserTask.h"

//	system includes
#include <iostream>
#include <string>

HcalLaserTask::HcalLaserTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalLaserTask::~HcalLaserTask()
{
}

/* virtual */ void HcalLaserTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalLaserTask);



