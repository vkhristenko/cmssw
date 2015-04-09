//	cmssw includes
#include "DQM/HcalTasks/interface/HcalHotCellTask.h"

//	system includes
#include <iostream>
#include <string>

HcalHotCellTask::HcalHotCellTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalHotCellTask::~HcalHotCellTask()
{
}

/* virtual */ void HcalHotCellTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
}

DEFINE_FWK_MODULE(HcalHotCellTask);



