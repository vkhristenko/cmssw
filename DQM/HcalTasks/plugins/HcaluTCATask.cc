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

/* virtual */void HcaluTCATask::beginLuminosityBlock(
		edm::LuminosityBlock const& lb, edm::EventSetup const& es)
{
	HcalDQSource::beginLuminosityBlock(lb, es);
}

/* virtual */void HcaluTCATask::endLuminosityBlock(
		edm::LuminosityBlock const& lb, edm::EventSetup const& es)
{
	HcalDQSource::endLuminosityBlock(lb, es);
}

/* virtual */ void HcaluTCATask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	edm::Handle<HFDigiCollection>		chf_vme;
	edm::Handle<HFDigiCollection>		chf_utca;

	INITCOLL(_labels["HFDigi"], chf_utca);
	INITCOLL(_labels["HFDigi_VME"], chf_vme);
	this->process(*chf_vme, *chf_utca, std::string("HF"));
}

//	specializer
template<typename Hit>
void HcaluTCATask::specialize(Hit const& hit1, Hit const& hit2, 
		std::string const& nameRes)
{

}

//	Reset
/* virtual */void HcaluTCATask::reset(int const periodflag)
{
	HcalDQSource::reset(periodflag);
	if (periodflag==0)
	{
		//	each event
	}
	else if (periodflag==1)
	{
		//	Each LS
	}
}

DEFINE_FWK_MODULE(HcaluTCATask);



