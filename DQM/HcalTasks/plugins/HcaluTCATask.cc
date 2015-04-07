//	cmssw includes
#include "DQM/HcalTasks/interface/HcaluTCATask.h"

//	system includes
#include <iostream>
#include <string>

HcaluTCATask::HcaluTCATask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcaluTCATask::~HcaluTCATask()
{
}

/* virtual */ void HcaluTCATask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	this->info_("Plugged and Running...");
}

DEFINE_FWK_MODULE(HcaluTCATask);



