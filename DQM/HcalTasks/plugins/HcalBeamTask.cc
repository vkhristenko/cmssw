//	cmssw includes
#include "DQM/HcalTasks/interface/HcalBeamTask.h"

//	system includes
#include <iostream>
#include <string>

HcalBeamTask::HcalBeamTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalBeamTask::~HcalBeamTask()
{
}

/* virtual */ void HcalBeamTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalBeamTask);



