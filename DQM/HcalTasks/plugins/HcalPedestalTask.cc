//	cmssw includes
#include "DQM/HcalTasks/interface/HcalPedestalTask.h"

//	system includes
#include <iostream>
#include <string>

HcalPedestalTask::HcalPedestalTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalPedestalTask::~HcalPedestalTask()
{
}

/* virtual */ void HcalPedestalTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalPedestalTask);



