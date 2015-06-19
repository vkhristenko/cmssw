//	cmssw includes
#include "DQM/HcalTasks/interface/HcalTimingTask.h"

//	system includes
#include <iostream>
#include <string>

HcalTimingTask::HcalTimingTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalTimingTask::~HcalTimingTask()
{
}

/* virtual */ void HcalTimingTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalTimingTask);



