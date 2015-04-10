//	cmssw includes
#include "DQM/HcalTasks/interface/HcalDigiTask.h"

//	system includes
#include <iostream>
#include <string>

HcalDigiTask::HcalDigiTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{
	_ornMsgTime = ps.getUntrackedParameter<int>("OrnMsgTime", 3559);
}

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

	this->process(*chbhe, std::string("HB"));
	this->process(*chbhe, std::string("HE"));
	this->process(*cho, std::string("HO"));
	this->process(*chf, std::string("HF"));

}

//	specializer
template<typename Hit>
void HcalDigiTask::specialize(Hit const& hit, std::string const& nameRes)
{
	//	offset of BCN for this channel(digi) relative to the nominal 
	//	in the unpacker. range(-7, +7), -1000 indicates invalid data
	int offset = hit.fiberIdleOffset();
	if (offset!=INVALID_FIBER_IDLEOFFSET && _ornMsgTime>-1)
		_mes[nameRes+"bcnOffset"].Fill(offset);
	
	//	
	for (int i=0; i<hit.size(); i++)
	{
		_mes[nameRes+"DigiShape"].Fill(i, hit.sample(i).adc());
		_mes[nameRes+"ADCCountPerTS"].Fill(hit.sample(i).adc());
		_mes[nameRes+"Presamples"].Fill(hit.presamples());
		_mes[nameRes+"CapId"].Fill(hit.sample(i).capid());
	}
}

DEFINE_FWK_MODULE(HcalDigiTask);



