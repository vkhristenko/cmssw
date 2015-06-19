//	cmssw includes
#include "DQM/HcalTasks/interface/HcalNoiseTask.h"

//	system includes
#include <iostream>
#include <string>

HcalNoiseTask::HcalNoiseTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalNoiseTask::~HcalNoiseTask()
{
}

/* virtual */ void HcalNoiseTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalNoiseTask);



