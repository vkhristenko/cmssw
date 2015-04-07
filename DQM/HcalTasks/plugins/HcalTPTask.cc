//	cmssw includes
#include "DQM/HcalTasks/interface/HcalTPTask.h"

//	system includes
#include <iostream>
#include <string>

HcalTPTask::HcalTPTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalTPTask::~HcalTPTask()
{
}

/* virtual */ void HcalTPTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	this->info_("Plugged and Running...");
}

DEFINE_FWK_MODULE(HcalTPTask);



