//	cmssw includes
#include "DQM/HcalTasks/interface/HcalRecHitTask.h"

//	system includes
#include <iostream>
#include <string>

HcalRecHitTask::HcalRecHitTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{
}

/* virtual */ HcalRecHitTask::~HcalRecHitTask()
{
}

/* virtual */ void HcalRecHitTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	edm::Handle<HBHERecHitCollection>	chbhe;
	edm::Handle<HORecHitCollection>		cho;
	edm::Handle<HFRecHitCollection>		chf;

	INITCOLL(_labels["HBHERecHit"], chbhe);
	INITCOLL(_labels["HORecHit"], cho);
	INITCOLL(_labels["HFRecHit"], chf);

	this->process(*chbhe, std::string("HB"));
	this->process(*chbhe, std::string("HE"));
	this->process(*cho, std::string("HO"));
	this->process(*chf, std::string("HF"));
}

//	specializer
template<typename Hit>
void HcalRecHitTask::specialize(Hit const& hit, std::string const& nameRes)
{
	float en	= hit.energy();
	float time	= hit.time();
	int ieta	= hit.id().ieta();
	int iphi	= hit.id().iphi();
	int depth	= hit.id().depth();

	//	fill Hcal-generic plots
	_mes["HcalEnergyPhi"].Fill(iphi, en);
	_mes["HcalEnergyEta"].Fill(ieta, en);
	if (nameRes!="HO")
	{
		_mes["HBHEHFEnergySumMapD" + 
			boost::lexical_cast<std::string>(depth)].Fill(
			ieta, iphi, en);
		_mes["HBHEHFTimeSumMapD" + 
			boost::lexical_cast<std::string>(depth)].Fill(
			ieta, iphi, en);
		_mes["HBHEHFOccupancyMapD" + 
			boost::lexical_cast<std::string>(depth)].Fill(
			ieta, iphi);
	}

	//	subsystem specific
	_mes[nameRes+"RecHitEnergy"].Fill(en);
	_mes[nameRes+"RecHitTime"].Fill(time);
}

DEFINE_FWK_MODULE(HcalRecHitTask);



