//	cmssw includes
#include "DQM/HcalTasks/interface/HcalDigiTask.h"

//	system includes
#include <iostream>
#include <string>

HcalDigiTask::HcalDigiTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalDigiTask::~HcalDigiTask()
{
}

/* virtual */ void HcalDigiTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	edm::Handle<HBHEDigiCollection>			chbhe;
	edm::Handle<HODigiCollection>			cho;
	edm::Handle<HFDigiCollection>			chf;
	
	INITCOLL(_labels["HBHEDigi"], chbhe);
	INITCOLL(_labels["HODigi"], cho);
	INITCOLL(_labels["HFDigi"], chf);

	this->process(*chbhe, std::string("HBHE"));
	this->process(*cho, std::string("HO"));
	this->process(*chf, std::string("HF"));

	//process<HBHEDigiCollection, HBHEDataFrame>(chbhe);

	/*
	for (HBHEDigiCollection::const_iterator it=chbhe.begin();
			it!=chbhe.end(); ++it)
	{
		const HBHEDataFrame digi = (const HBHEDataFrame)(*it);
		for (int i=0; i<digi.size(); i++)
			_mes["HBHEDigiShape"].Fill(i+1, digi.sample(i).adc());
	}*/
}

//	specializer
template<typename Hit>
void HcalDigiTask::specialize(Hit const& hit, std::string const& nameRes)
{
	for (int i=0; i<hit.size(); i++)
		_mes[nameRes+"DigiShape"].Fill(i+1, hit.sample(i).adc());
}

DEFINE_FWK_MODULE(HcalDigiTask);



