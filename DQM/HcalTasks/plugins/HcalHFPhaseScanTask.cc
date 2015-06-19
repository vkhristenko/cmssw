//	cmssw includes
#include "DQM/HcalTasks/interface/HcalHFPhaseScanTask.h"

//	system includes
#include <iostream>
#include <string>

HcalHFPhaseScanTask::HcalHFPhaseScanTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalHFPhaseScanTask::~HcalHFPhaseScanTask()
{
}

/* virtual */ void HcalHFPhaseScanTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	edm::Handle<HBHEDigiCollection>			cdhbhe;
	edm::Handle<HODigiCollection>			cdho;
	edm::Handle<HFDigiCollection>			cdhf;
	edm::Handle<HBHERecHitCollection>		crhbhe;
	edm::Handle<HORecHitCollection>			crho;
	edm::Handle<HFRecHitCollection>			crhf;

	INITCOLL(_labels["HBHEDigi"], cdhbhe);
	INITCOLL(_labels["HODigi"], cdho);
	INITCOLL(_labels["HFDigi"], cdhf);
	INITCOLL(_labels["HBHERecHit"], crhbhe);
	INITCOLL(_labels["HORecHit"], crho);
	INITCOLL(_labels["HFRecHit"], crhf);
/*
	this->process(*cdhbhe, std::string("HB"));
	this->process(*cdhbhe, std::string("HE"));
	this->process(*cdho, std::string("HO"));
	this->process(*cdhf, std::string("HF"));
	this->process(*crhbhe, std::string("HB"));
	this->process(*crhbhe, std::string("HE"));
	this->process(*crho, std::string("HO"));
	this->process(*crhf, std::string("HF"));
	*/
}

//	specializer
//template<typename Hit>
//void HcalDigiTask::specialize(Hit const& hit, std::string cosnt& nameRes)
//{
	/*
	int ieta = hit.id().ieta();
	int iphi = hit.id().iphi();
	int depth = hit.id().depth();

	//	the only module that uses both Digis and RecHits
	if (std::typeid(hit)==typeid(HBHEDataFrame) || 
			std::typeid(hit)==typeid(HFDataFrame) ||
			std::typeid(hit)==typeid(HODataFrame))
	{

	}
	else if (std::typeid(hit)==typeid(HBHERecHit) ||
			std::typeid(hit)==typeid(HFRecHit) ||
			std::typeid(hit)==typeid(HORecHit))
	{

	}
	*/
//}

DEFINE_FWK_MODULE(HcalHFPhaseScanTask);



