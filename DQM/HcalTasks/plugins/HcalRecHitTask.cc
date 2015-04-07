//	cmssw includes
#include "DQM/HcalTasks/interface/HcalRecHitTask.h"

//	system includes
#include <iostream>
#include <string>

HcalRecHitTask::HcalRecHitTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalRecHitTask::~HcalRecHitTask()
{
}

/* virtual */ void HcalRecHitTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	this->info_("Plugged and Running...");
}

DEFINE_FWK_MODULE(HcalRecHitTask);



