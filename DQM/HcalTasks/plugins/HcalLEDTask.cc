//	cmssw includes
#include "DQM/HcalTasks/interface/HcalLEDTask.h"

//	system includes
#include <iostream>
#include <string>

HcalLEDTask::HcalLEDTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalLEDTask::~HcalLEDTask()
{
}

/* virtual */ void HcalLEDTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalLEDTask);



