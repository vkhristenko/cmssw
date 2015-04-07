//	cmssw includes
#include "DQM/HcalTasks/interface/HcalRawTask.h"

//	system includes
#include <iostream>
#include <string>

HcalRawTask::HcalRawTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalRawTask::~HcalRawTask()
{
}

/* virtual */ void HcalRawTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	this->info_("Plugged and Running...");
}

DEFINE_FWK_MODULE(HcalRawTask);



