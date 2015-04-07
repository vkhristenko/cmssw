//	cmssw includes
#include "DQM/HcalTasks/interface/HcalDeadCellTask.h"

//	system includes
#include <iostream>
#include <string>

HcalDeadCellTask::HcalDeadCellTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalDeadCellTask::~HcalDeadCellTask()
{
}

/* virtual */ void HcalDeadCellTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalDeadCellTask);



